!========================================================================
!
!  BVECTORLIB.f90
!
!  vector routines for BIGSTICK
!
!  routines to carry out vector algebra and storage
!  including for parallel implementations  
!  initiated June 2013 @ SDSU
!
module bvectorlib_mod
contains
!=======================================================================
!  SUBROUTINES IN THIS FILE
!
!

!  MODIFIED FROM BIGSTICK ROUTINE setup_localvectors
!
!  ALLOCATE "LOCAL" FRAGMENTS OF VECTORS
!
!  here basestart(f), basestop(f) are start and stop
!  of the basis states for fragment f
!
!  Default for 1 processor is frag1 = frag2 = 1 (only one "fragment")
!
! Also builds fragment level communicators
!
! CALLED BY
!    lanczos_p
!    exactdiag_p
!    density1b_from_oldwfn
!    overlap
!    applicator_h   ( apply scalar 2-body opeator)
!    applicator1b   ( apply nonscalar 1-body operator)
!    particle_occupation_p
!    particle_occupation_p_orig
!
!  SUBROUTINES CALLED:
!   setup_mpi_hist
!   memreport
!
subroutine setup_localvectors_both
   use flagger
   use localvectors
   use io
   use nodeinfo
   use fragments
!   use lanczos_info
   use basis
!   use tribution
   use mod_reorthog
   use bmpi_mod
   use butil_mod
   implicit none
   character :: whichfile

   integer :: key, ierr
   integer :: aerr
   integer :: ii
   integer :: num_threads
   integer(kind=basis_prec) :: vi
   integer :: tid
   
   !integer(kind=basis_prec),pointer :: dimbasis_local

   ! OpenMP functions
   integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
   
!   if(whichfile/='I' .and. whichfile/='F')then
!	   print*, ' problem with which file ',whichfile
!	   stop
!   end if
   
   frag1 = nodal(iproc)%ifragment
   frag2 = nodal(iproc)%ffragment
   if(frag1 < 1 .or. frag1 > nfragments) then ! This nfragments reference is ok
      print *, "setup_localvectors: ", iproc, " out of range, frag1=", frag1
      !write(logfile,*) "setup_localvectors: ", iproc, " out of range, frag1=", frag1
!      flush(logfile)
      stop 1
   end if
   if(frag2 < 1 .or. frag2 > nfragments) then ! This nfragments reference is ok
      print *, "setup_localvectors: ", iproc, " out of range, frag2=", frag2
      !write(logfile,*) "setup_localvectors: ", iproc, " out of range, frag2=", frag2
!      flush(logfile)
      stop 1
   end if
   ! use to access slice regardless of how vec1/vec2 are allocated
   v1s = basestart_i(frag1)
   v1e = basestop_i(frag1)
   v2s = basestart_f(frag2)
   v2e = basestop_f(frag2)
   !print*, "setuplocalvectors: ", iproc, ", frag1=", frag1, ", frag2=", frag2
   !print*, "setuplocalvectors: ", iproc, ", v1s=", v1s, ", v1e=", v1e, ", v2s=", v2s, ", v2e=", v2e
   !print*, "setuplocalvectors: ", iproc, "
   if(v1e < v1s .or. v2e < v2s) then
      print *, "setuplocalvectors: ", iproc, ", range error ve < vs, stopping"
	   print*,v1e,v1s,v2e,v2s
      stop 1
   end if
   
      allocate( vec1(v1s:v1e), stat=aerr )       ! holds previous vector
      if(aerr /= 0) call memerror("setup_localvectors 3")
      allocate( vec2(v2s:v2e), stat=aerr )       ! holds new vector
      if(aerr /= 0) call memerror("setup_localvectors 4")

   ! allocate memory for each thread to write to in OpenMP mode
!$omp parallel
   num_threads = omp_get_num_threads();
!$omp end parallel
   ompNumThreads = num_threads
   useVec2Thread = useNewReorthog .and. (ompNumThreads > 1) .and. wantUseVec2Thread
   if(useVec2Thread) then
      if(iproc == 0) print *, "Allocating vec2thread, threadcount=", ompNumThreads
      !! allocate(vec2thread(v2s:v2e, 0:ompNumThreads-1), stat=aerr)
      vec2threadchunk = v2e - v2s + 1          ! size of chunk
      vec2threadchunkm1 = vec2threadchunk - 1  ! to find ent of slice
      vec2threadend = ompNumThreads * (v2e - v2s + 1) - 1  ! end of buffer
      allocate(vec2threadflat(0 : vec2threadend), stat=aerr)
      if(aerr /= 0) call memerror("setup_localvectors vec2thread - set wantUseVec2Thread=.false. to save memory");
      !! vec2thread(:, :) = 0.0
      vec2threadflat(:) = 0.0
   end if

   if(storevectorsenabled .and. nproc>=niter)then
      writetodisk = .false.
   else
      writetodisk = .true.
   end if

   ! Seems like a good place to set up history mechanism
   ! HOWEVER MAY WAIT ON THIS
!   if(iproc == 0) print *, "KSM:  Calling setup_mpi_hist(", niter, ")"
!   call setup_mpi_hist(niter)
!   if(noisy0) print *, "Returned from setup_mpi_hist"
!   call memreport     ! generates a report on storage ADDED 7.3.8
!   if(noisy0) print *, "Generated memory report"

   ! KSM - 15Aug2014
   ! Build communicators for allreduce and for broadcast from root node of each fragment
   ! will work fine even if only one fragment
   
   !print*,"iproc",iproc,"isfragroot",isfragroot
   if(isfragroot) then
       key = 0  ! root has to get rank 0
   else
      key = iproc
   endif

   !if(iproc == 0) print*,"right before BMPI_COMM."
   ! These communicators are split by the fragment they serve
   !print*,"icomm",icomm,"frag1",frag1,"frag2",frag2,"key",key
   call BMPI_COMM_SPLIT(icomm, frag2, key, fcomm2, ierr)  
   call BMPI_COMM_SPLIT(icomm, frag1, key, fcomm1, ierr)
   !print*,"fcomm1",fcomm1
   !print*,"ierr",ierr
   !if(noisy0) print *, "Done making fcomm1, and fcomm2"
   !if(noisy0) print *, "Done making fcomm1, and fcomm2"
   if(iproc == 0) print*,"finished setting up local vectors..."
   return
end subroutine setup_localvectors_both
!
!=======================================================================
!  double-precision normalization
!
!  -- for 'new' parallelization
!
! vchar :  information on vector (normal or reverse ordering
! ipoint:  point to 'i' (initial) or 'f' (final) vector
!
! dnorm : double-precision norm of dvec
! smallflag: logical flag = .true. if norm is smaller than dtol
!
! CALLED BY:
!    initialize_lanczos_vector
!    lanczos_p
!    random_restart_p
!    expectator_p
!==========================================================
subroutine dnormvec_p(vchar,ipoint,dnorm,smallflag)
  use precisions
  use nodeinfo
!  use timing
  use localvectors
  use fragments
  use bmpi_mod
  implicit none

  integer(4) :: ierr
  character(1) :: vchar,ipoint
  real(kind=lanc_prec), pointer :: dvec(:)
  real(kind=8) :: dnorm
  logical smallflag
!------------------ INTERMEDIATE -------------------------------------------
  real(kind=8) :: d
  real(kind=8) :: dtol
  real(kind=8) :: tdnorm

  integer(4)   :: i

  integer(kind=basis_prec) :: il,vstart,vstop

 ! call clocker('dot','sta')
  dtol = 1.d-8
  dnorm = 0.d0
  smallflag = .false.

!..... choose which vector to normalize
!      complicated because during lanczos 
!      we switch between vec1 and vec2 (stored in module localvectors
!      also can normalize initial or final vector

  if( (vchar == 'n' .and. ipoint == 'i' ) .or. & 
        (vchar == 'r' .and. ipoint == 'f') )then
        dvec => vec1
        vstart = v1s
        vstop  = v1e
  else
        dvec => vec2
        vstart = v2s
        vstop  = v2e
  end if
!$omp parallel do private(il,d), shared(dvec, vstart,vstop), reduction(+ : dnorm) 
   do il= vstart,vstop
        d = real(dvec(il),kind=8)
        dnorm = dnorm+d*d
   end do
!$omp end parallel do

!..... IF USING MPI with fragments, MUST COMBINE HERE
! KSM:  We only want to provide a result from the nodes
!       labled isfragroot.  This selects one representative
!       for each fragment.  Otherwise we would be double counting.
! This routine should only be used sparingly because it wastes
! many nodes.   On the other hand, the cost of distributing the
! vector slices using breorthog.f90 can be high.
!
  if(.not. isfragroot .and. nproc > 1) dnorm = 0.d0
  
  call BMPI_ALLREDUCE(dnorm, tdnorm, 1,  MPI_SUM, icomm, ierr)
  dnorm = tdnorm

  if ( dnorm < dtol ) then
     smallflag = .true.
     if ( iproc == 0 )write(6,*)' zero vector ',dnorm
     return
  end if

  dnorm = dsqrt(dnorm)
  d = 1.d0 / dnorm

!$omp parallel do private(il), shared(dvec, d, vstart,vstop)
     do il = vstart,vstop
        dvec(il) = real(dvec(il)*d,kind=lanc_prec)
     end do
!$omp end parallel do

!  call clocker('dot','end')
  return
end subroutine dnormvec_p

!=================================================================
!  double-precision projection
!
!
! n: dimension of vector
! dvec1: double-precision vector
! dvec2: double-preciscion vector
!
! dsclrprod  = dvec1*dvec2
!
! dvec1 -> dvec1 - dvec2*dsclrprod
! returns \alpha - the overlap with the previous vector
!
! CALLED BY: reorthogonalize_a
!====================================================================
subroutine dvecproj_p(vchar,dsclrprod)
  use precisions
  use nodeinfo
!  use timing
!  use lanczos_info
  use localvectors
  use fragments
  use bmpi_mod
  implicit none

  character(1) :: vchar
  integer(4)               :: ierr
  real(kind=8)             :: dsclrprod
  real(kind=lanc_prec), pointer :: dvec1(:),dvec2(:)

!------------------ INTERMEDIATE -------------------------------------------
  real(kind=8)             :: d1,d2
  integer(4)               :: i

  integer(kind=basis_prec) :: il,vstart,vstop

  if(nfragments > 1) then   ! This nfragments reference is ok, not about new reorthog
     if(iproc == 0) print *, "dvecproj_p: not supported with new reorthog"
     stop 1 ! has msg
  end if

  if(.not.storelanczosincoreMPI .and. iproc > 0)return
!...............................

!  call clocker('pro','sta')
  dsclrprod = 0.d0

       ! This code doesn't work when frag1 /= frag2
       vstart = v1s  !    basestart(frag1)
       vstop  = v1e  !    basestop (frag1)

  select case (vchar)
          case ('n')
            dvec1 => vec1
            dvec2 => vec2

          case ('r')
            dvec1 => vec2
            dvec2 => vec1

        case default
           stop 1   ! prevent uninitialized compiler msg
  end select

!$omp parallel do private(il, d1, d2), shared(vstart,vstop, dvec1, dvec2), reduction(+ : dsclrprod)
  do il =vstart,vstop
        d1 = real(dvec1(il),kind=8)
        d2 = real(dvec2(il),kind=8)
        dsclrprod = dsclrprod+d1*d2
  end do
!$omp end parallel do

!$omp parallel do private(il, d1, d2), shared(vstart,vstop, dvec1, dvec2,dsclrprod)
  do il = vstart,vstop
        d1 = real(dvec1(il),kind=8)
        d2 = real(dvec2(il),kind=8)
        dvec1(il) = real(d1 - d2*dsclrprod,kind=lanc_prec)
  end do
!$omp end parallel do
!  call clocker('pro','end')
  return
end subroutine dvecproj_p

!===================================================
      subroutine dvecdot_p(vchar,dsclrprod)
!
!  double-precision projection
!
!  vchar: if = 'n' vec1 = initial, vec2 = final
!         if = 'r' vec2 = initial, vec1 = final
!
! dsclrprod  = dvec1*dvec2
!
!
      use precisions
	use nodeinfo
      use localvectors
      use fragments
      implicit none
      character(1) :: vchar
      real(kind=lanc_prec), pointer :: dvec1(:),dvec2(:)
      real(kind=8) :: dsclrprod

!------------------ INTERMEDIATE -------------
      real(kind=8) :: d1,d2
      integer      :: i

      integer(kind=basis_prec) :: il,vstart,vstop

   if(nfragments > 1) then  ! this nfragments reference is ok
      if(iproc == 0) print *, "dvecdot_p: not supported with nfragments>1"
      stop 1
   end if

	! This code didn't work before when frag1 /= frag2
	vstart = v1s  !    basestart(frag1)
	vstop  = v1e  !    basestop (frag1)

   select case (vchar)
    case ('n')
      dvec1 => vec1
      dvec2 => vec2

    case ('r')
      dvec1 => vec2
      dvec2 => vec1

    case default
      stop 1   ! prevent uninitialized compiler msg

   end select
   dsclrprod = 0.d0
!$OMP PARALLEL SHARED(dvec1,dvec2,vstart,vstop), PRIVATE(i,d1,d2), REDUCTION(+:dsclrprod)
!$OMP DO SCHEDULE(STATIC)
   do il = vstart,vstop
         d1 = dvec1(il)
         d2 = dvec2(il)
         dsclrprod = dsclrprod+d1*d2
   enddo
!$OMP END DO
!$OMP END PARALLEL

   return
end subroutine dvecdot_p

!================================================================
      subroutine doverlapvec(dsclrprod)
!
!  double-precision overlap
!  added in 7.7.0
!
! dsclrprod  = dvec1*dvec2
!
! CALLED BY overlap
!
      use precisions
	  use nodeinfo
      use localvectors
      use fragments
	  use bmpi_mod
      implicit none
      real(kind=8) :: dsclrprod,tmpprod

!------------------ INTERMEDIATE -------------
      real(kind=8) :: d1,d2
      integer      :: i

      integer(kind=basis_prec) :: il,vstart,vstop
	  integer :: ierr


	vstart = v1s  !    basestart(frag1)
	vstop  = v1e  !    basestop (frag1)

   dsclrprod = 0.d0
   if(isfragroot .or. nproc==1)then
!$OMP PARALLEL SHARED(vec1,vec2,vstart,vstop), PRIVATE(i,d1,d2), REDUCTION(+:dsclrprod)
!$OMP DO SCHEDULE(STATIC)
   do il = vstart,vstop
         d1 = vec1(il)
         d2 = vec2(il)
         dsclrprod = dsclrprod+d1*d2
   enddo
!$OMP END DO
!$OMP END PARALLEL

end if

call BMPI_REDUCE(dsclrprod,1,MPI_SUM,0,icomm,ierr)
!call BMPI_BARRIER(icomm,ierr)
!print*,iproc,vstart,vstop,dsclrprod

!dsclprod=tmpprod
   return
end subroutine doverlapvec


!================================================================
      subroutine doverlapvec_new(dsclrprod)
!
!  double-precision overlap
!  added in 7.7.0
!
! dsclrprod  = dvec1*dvec2
!
! CALLED BY overlap
!
      use precisions
	  use nodeinfo
      use localvectors
      use fragments
	  use bmpi_mod
	  use basis_map
	  use sectors
      implicit none
      real(kind=8) :: dsclrprod,tmpprod

!------------------ INTERMEDIATE -------------
      real(kind=8) :: d1,d2
      integer      :: i,csi
	  
	  integer    :: isp,isn

      integer(kind=basis_prec) :: il,vstart,vstop,isdp,isdn,fsdp,fsdn
	  integer(kind=basis_prec) :: ipstate,fpstate,instate,fnstate
	  integer(2) :: pphase,nphase
	  integer :: ierr


	vstart = v1s  !    basestart(frag1)
	vstop  = v1e  !    basestop (frag1)

   dsclrprod = 0.d0
   if(isfragroot .or. nproc==1)then
!!$OMP PARALLEL SHARED(vec1,vec2,vstart,vstop), PRIVATE(i,d1,d2), REDUCTION(+:dsclrprod)
!!$OMP DO SCHEDULE(STATIC)
!   do il = vstart,vstop
!         d1 = vec1(il)
!         d2 = vec2(il)
!         dsclrprod = dsclrprod+d1*d2
!   enddo
!!$OMP END DO
!!$OMP END PARALLEL

!......... NEW ALGORITHM ALLOWS FOR MAPPING......
!          BECAUSE I AM COMPARING 2 BASES IT LOOKS COMPLICATED

do isp = 1,nsectors_i(1)  ! loop over proton sectors
	do isdp = xsd_i(1)%sector(isp)%xsdstart,xsd_i(1)%sector(isp)%xsdend	  ! LOOP OVER SDs in that initial sector
		fsdp = psdmap(isdp)    ! MAP initial proton SD to final proton SD
		if(fsdp==-1)cycle      ! if no mapping, cycle
		ipstate = pstart_i(isdp)    ! find index of initial proton state
		fpstate = pstart_f(fsdp)     ! find index of final proton state
		pphase  = psdmapphase(isdp)   ! find any phase from mapping
		
		do csi= 1,xsd_i(1)%sector(isp)%ncsectors     ! loop over conjugate neutron sectors
			isn = xsd_i(1)%sector(isp)%csector(csi)   ! find the index of the conjugate neutron sector
			do isdn = xsd_i(2)%sector(isn)%xsdstart,xsd_i(2)%sector(isn)%xsdend	
				fsdn = nsdmap(isdn)
				if(fsdn==-1)cycle
				nphase = nsdmapphase(isdn)
				instate = nstart_i(isdn)
				fnstate = nstart_f(fsdn)
!				print*,ipstate,instate,fpstate,fnstate
				
				d1=vec1(ipstate+instate)
				d2 = vec2(fpstate+fnstate)
!				print*,isdn,fsdn,instate,fnstate
!print*,d1,d2*pphase*nphase,ipstate,fpstate,instate,fnstate
				dsclrprod = dsclrprod+d1*d2*pphase*nphase
				
				
				
			end do
			
			
		end do !csi
		
		
	end do	
	
	
	
	
end do ! isp

end if

call BMPI_REDUCE(dsclrprod,1,MPI_SUM,0,icomm,ierr)
!call BMPI_BARRIER(icomm,ierr)
!print*,iproc,vstart,vstop,dsclrprod

!dsclprod=tmpprod
   return
end subroutine doverlapvec_new



!================================================================
      subroutine projectvec
! modified by -RMZ 3/6/22 RZ added MPI functionality
! 
! CALLED BY overlap
!
     use precisions
	  use nodeinfo
     use localvectors
     use fragments
	  use bmpi_mod
	  use basis_map
	  use sectors


      implicit none
      real(kind=8) :: dsclrprod,tmpprod

!------------------ INTERMEDIATE -------------
      real(kind=8) :: d1,d2
      integer      :: i,csi
	  
	  integer    :: isp,isn

   integer(kind=basis_prec) :: il,vstart,vstop,isdp,isdn,fsdp,fsdn
	  integer(kind=basis_prec) :: ipstate,fpstate,instate,fnstate
	  integer(2) :: pphase,nphase
     integer(8) st, et
	  integer :: ierr
     integer :: i_start, i_end

   !usage of proc_range subroutine in rhfragments.f90
   

	vstart = v1s  !    basestart(frag1)
	vstop  = v1e  !    basestop (frag1)

   !print*,"v2s---->", v2s, "v2e---->",v2e
   !break up projection algorithm over proton sectors given to each fragment

   if(nproc .gt. 1)then
      i_start = fragmentlist_i(frag1)%ssectorstart
      i_end = fragmentlist_i(frag1)%ssectorend

      st = 1
      et = nsectors_i(1)
      !if mpi enabled and vec1 fits on each frag then just divide up istart and stop over mpi ranks
      if(nproc .gt. 1) then
         if((i_start .eq. 1) .and. ((i_end .eq. nsectors_i(1))) ) then
            call proc_range(st, et, nproc, iproc, i_start, i_end) 
         endif 
      endif

   else ! serial
      i_start = 1
      i_end = nsectors_i(1)
   endif

   
!!$OMP PARALLEL SHARED(vec1,vec2,vstart,vstop), PRIVATE(i,d1,d2), REDUCTION(+:dsclrprod)
!!$OMP DO SCHEDULE(STATIC)
!   do il = vstart,vstop
!         d1 = vec1(il)
!         d2 = vec2(il)
!         dsclrprod = dsclrprod+d1*d2
!   enddo
!!$OMP END DO
!!$OMP END PARALLEL

!......... NEW ALGORITHM ALLOWS FOR MAPPING......
!          BECAUSE I AM COMPARING 2 BASES IT LOOKS COMPLICATED
!print*,"iproc",iproc,"i_start", i_start, "i_end",i_end
do isp = i_start,i_end
	do isdp = xsd_i(1)%sector(isp)%xsdstart,xsd_i(1)%sector(isp)%xsdend	  ! LOOP OVER SDs in that initial sector	
		fsdp = psdmap(isdp)
		if(fsdp==-1)cycle
		ipstate = pstart_i(isdp)
		fpstate = pstart_f(fsdp)
		pphase  = psdmapphase(isdp)
		
		do csi= 1,xsd_i(1)%sector(isp)%ncsectors
			isn = xsd_i(1)%sector(isp)%csector(csi)
			do isdn = xsd_i(2)%sector(isn)%xsdstart,xsd_i(2)%sector(isn)%xsdend	
				fsdn = nsdmap(isdn)
				if(fsdn==-1)cycle
				nphase = nsdmapphase(isdn)
				instate = nstart_i(isdn)
				fnstate = nstart_f(fsdn)
				
            !if(iproc==0) print*,"iproc",iproc,"vec2",fpstate+fnstate,"vec1",ipstate+instate
            !if(iproc==0) print*,"isp",isp,"fpstate",fpstate,"fnstate",fnstate,"ipstate",ipstate,"instate",instate
				vec2(fpstate+fnstate) = vec1(ipstate+instate)*pphase*nphase
				
				
			end do
			
			
		end do !csi
		
		
	end do	
	
	
	
!print*,"iproc",iproc,"vec2",fpstate+fnstate,"vec1",ipstate+instate	
end do ! isp
!print*,"------------------------------------------------------------"

!call BMPI_REDUCE(vec2,1,MPI_SUM,0,icomm,ierr)
!call BMPI_BARRIER(icomm,ierr)
   return
end subroutine projectvec

!================================================
subroutine initialize_final(vchar)
   use nodeinfo
   use fragments
   use localvectors
   use precisions
   implicit none
   character(1) :: vchar

   integer(kind=basis_prec) :: i
   integer :: tid
   real(kind=lanc_prec), pointer :: v2p(:)
   ! OpenMP functions
   integer(kind=4) :: omp_get_thread_num, omp_get_num_threads

   if(useVec2Thread) then
!$omp parallel private(i, tid, v2p)  &
!$omp       shared(v2s, v2e)         &
!$omp       shared(vec2threadchunk, vec2threadchunkm1, vec2threadflat)
      tid = omp_get_thread_num();
      i = tid * vec2threadchunk
      v2p(v2s:v2e) => vec2threadflat(i: i + vec2threadchunkm1) 
      v2p(:) = 0.0
!$omp end parallel
      ! fall through and initialize vec2 as well
      ! needed anyway for final reduce
   end if

   select case(vchar)
      case('n')

      do i = basestart(frag2), basestop(frag2)
         vec2(i) = 0.0e0_lanc_prec
      end do

      case('r')

      do i = basestart(frag1), basestop(frag1)
         vec1(i) = 0.0e0_lanc_prec
      end do
   end select

   return
end subroutine initialize_final
!  

end module bvectorlib_mod
