
!=================================================================================
!
!  bfragments.f90
!
!  routines to break up Lanczos vectors into "fragments"  
!
!============================================================
!
!  The next module is used for 
!  breaking up the Lanczos vectors 
!
!=========================================================
!
!  FRAGMENTS
!
!   contiguous (sub) sectors of Lanczos vectors grouped together
!   
! started 8/2011 by CWJ
!

module fragments
   use precisions
!   use operation_stats
   implicit none

   logical :: fragment_vectors   ! flag to break up vectors into fragments
   integer :: nfragments         ! number of fragments
   logical :: useNewReorthog1 = .true.  ! Use new reorthog even for single fragment
   logical :: useNewReorthog  ! set to (nfragments > 1) .or. useNewReorthog1

   type frag
        integer :: nsectors  ! # of subsectors
        integer(8) :: ssectorstart,ssectorend  ! start, end of sectors; these are defined by "proton" sectors
		  integer(8) :: csectorstart,csectorend  ! start, end of conjugate sectors, defined by "neutron" sectors, must check contiguity
        integer(8):: basisstart,basisend  ! where the basis starts and stops for this sector
        integer(8):: localdim      ! local dimension = basisend-basisstart+1
        real(8)    :: totops    ! tot# of operations to/from this fragment

   end type frag

   type (frag), allocatable,target :: fragmentlist_i(:),fragmentlist_f(:)

!   integer, allocatable :: sectormap2frag(:)  ! maps a sector to a fragment  SCHEDULED FOR OBSOLESCENCE

!--------- THIS ALLOWS US TO STORE JUST FRAGMENTS OF A VECTOR --

   integer (kind=basis_prec), allocatable :: basestart(:),basestop(:)  !NEED TO GET RID OF LATER

   integer (kind=basis_prec), allocatable,target :: basestart_i(:), basestop_i(:), basestart_f(:), basestop_f(:)
   
   logical :: test_fragments = .false.  ! added in 7.4.6, used to test new fragment scheme on a single processor
   
!------------- JUMP STORAGE ACROSS FRAGMENTS (f2f matvec blocks)-------------- added in 7.6.2 -----------------

   type f2fjump
	   integer(8) :: nX1bjump(2),nXX2bjump(2),nXXX3bjump(2)
	   real(8)    :: totjumpstorage
	   integer(4) :: minprocs     ! minimal assigned MPI procs/nodes
	   integer(4) :: nprocs       ! # of actual assigned procs in the end
   end type f2fjump   
   
   type (f2fjump), allocatable,target :: f2fjumpstat(:,:)

!------- INFORMATION ON WHAT FRAGMENTS ON WHAT MPI PROCESSES
!        previously in module tribution, moved in 7.6.7
	   
	type whatsonanode
	      integer :: ifragment,ffragment   ! which fragments are here
	      logical :: ifirst ! first node in group servicing ifragment
		  
!	      integer :: sfragment  ! which fragment is stored here; later might want to store more than one
!	      real :: frac_ops   ! what fraction of possible ops are stored here      
	end type whatsonanode

	type (whatsonanode), allocatable :: nodal(:)	   

contains
!======================================================================
!
!  SUBROUTINES CONTAINED IN THIS FILE
!    fragmenter:  master routine to break up Lanczos vectors into fragments
!        which calls:
!    vectorstorageceiling
!    defaultfragments:  routine to set fragment information when only 1 fragment
!    count_create_fragments_blocks: updated fragmentation routine which breaks fragments along haiku blocks
!    reconstruct_sector: when sectors are subdivided, reconstruct proton sectors
!    fragmentvectorsize:  sets base start/stop for a given fragment
!       also in this file (called by distro_opbundles_over_fragments )
!    divide_procs_over_fragments: routine to divide up available nodes over fragments
!
!=======================================================================
!
!  Notes 7.4.6 (March 2015) regarding "new" fragmentation based upon proton haiku blocks
!  To test/debug fragmentation, run on a single processor;
!  in order to do this, we must modify routines which divide up work over MPI processes
!   
!  key routines include: 
!          setnodaltribution  (in bparallel_lib1.f90)  assigns initial, final fragments to MPI procs/compute nodes
!          setup_localvectors (in blanczoslib1.f90)  allocates memory for (fragmented) lanczos vectors
!  If one simulates a fragmented case on a single processor, the above have to be modified
!

!======================================================================
!
!  proc_range
!   divides up list(sectors) based on nproc and current mpi rank.
!   returns start and stop of sectors for each rank
!
! CALLED BY:
!  proj_boss rhprojlib.f90
!
subroutine proc_range(lv, hv, nprocs, irank, istart, iend)
   implicit none
   integer(8) :: lv, hv ! Low and high iteration variables.
   Integer(4) :: nprocs, irank, istart, iend
   integer(4) ::  q, r ! quotient and remainder
   
   q = (hv-lv + 1)/nprocs
   r = MOD(hv-lv+1,nprocs)
   
   istart = lv + irank*q + MIN(irank,r) ! starting iter
   iend = istart + q - 1 ! final iter
   
   if (r .gt. irank) then
      iend = iend + 1
   endif
   return
end subroutine proc_range

!======================================================================
!
!  fragmentvectorsize
!     sets base start/stop for a given fragment
!
! CALLED BY:
!  set_nfragments
!
subroutine fragmentvectorsize
   !   use fragments
      use basis
      use nodeinfo
      use sectors
      use bmpi_mod
      use localvectors
      use flagger
      implicit none
      integer :: fg
      integer :: aerr

   !assign intial and final basis start/stop to fragments based on fragment sectors. 
      if(break_vectors) then
      do fg = 1, nfragments
         basestart_i(fg)=xsd_i(1)%sector(fragmentlist_i(fg)%ssectorstart)%basisstart
         basestop_i(fg)= xsd_i(1)%sector(fragmentlist_i(fg)%ssectorend)%basisend

         basestart_f(fg)=xsd_f(1)%sector(fragmentlist_f(fg)%ssectorstart)%basisstart
         basestop_f(fg)= xsd_f(1)%sector(fragmentlist_f(fg)%ssectorend)%basisend
       enddo
      else ! 1 fragment
         basestart_i(1)=1
         basestop_i(1)=dimbasis_i
         basestart_f(1)=1
         basestop_f(1)=dimbasis_f
      endif
      return
   end subroutine fragmentvectorsize
   
!======================================================================
   


!======================================================================
!  subroutine set_nfragments
! Sets nfragments and side effects of assignment
subroutine set_nfragments(nf)
   !  use fragments
     use basis
     use nodeinfo
     use sectors
     use bmpi_mod
     use localvectors
     use flagger

     implicit none
     integer :: jproc, nf
     integer :: i ! mpiS
     integer :: is, ie
     integer :: ierr, aerr
     integer :: basis_istart, basis_iend
     integer :: sectstart, sectend
     integer(8) :: st, fg
     st = 1
     is = 0
     ie = 0
  
     nfragments = nf
     useNewReorthog = useNewReorthog1 .or. nfragments > 1
     if(iproc == 0) print *, "useNewReorthog=", useNewReorthog
     allocate(nodal(0:nproc-1))
     
     
     !if (iproc == 0)then
     ! print*,"printing intitial space sector information..."
     !   do i = 1, nsectors_i(1)
     !      print*, "i", i, "xsd_i(1)%sector(i)%jzX", xsd_i(1)%sector(i)%jzX,"xsd_i(1)%sector(i)%parx",xsd_i(1)%sector(i)%parx
     !   end do
     !   print*,"finished"
     !end if 

     
    
     do jproc = 0,nproc-1
        nodal(jproc)%ifragment= 1 + jproc
        nodal(jproc)%ffragment= 1 + jproc
     end do
     
     allocate(basestart_i(nf), basestop_i(nf))
     allocate(basestart_f(nf), basestop_f(nf))
   
     allocate(fragmentlist_i(nf), stat=aerr)
     allocate(fragmentlist_f(nf), stat=aerr)
     
     
     !if(nfragments == 1) then
     !   basestart_i(1)=1
     !   basestop_i(1)=dimbasis_i
     !   basestart_f(1)=1
     !   basestop_f(1)=dimbasis_f
  
     !else if(nfragments .gt. 1) then

   ! used temporarily to store fragment start stop sector info
   !   allocate(fragmentlist_i(nf), stat=aerr)
   !   allocate(fragmentlist_f(nf), stat=aerr)

      !assign sectors from the intital and final basis to each fragment
      !if (iproc == 0 )then
      !  print*,"printing final sector information..."
      !     do i = 1, nsectors_f(1)
      !        print*, "i", i, "xsd_f(1)%sector(i)%jzX", xsd_f(1)%sector(i)%jzX,"xsd_f(1)%sector(i)%parx",xsd_f(1)%sector(i)%parx
      !     end do
      !endif

      !we divide the final basis into nproc fragments by sector
      ! each fragment has some range [istart,istop] of contigous sectors assigned to it
      !call create_fragments

      !set basis start/stop using fragment sector info
      !call fragmentvectorsize

     !if (iproc == 0 )then
     ! print*,"printing final sector information..."
     !   do i = 1, nsectors_f(1)
     !      print*, "i", i, "xsd_f(1)%sector(i)%jzX", xsd_f(1)%sector(i)%jzX,"xsd_f(1)%sector(i)%parx",xsd_f(1)%sector(i)%parx
     !   end do
   ! endif
     !end if 
    
     !call BMPI_BARRIER(icomm,ierr)
     return
  
  end subroutine set_nfragments
  
  subroutine assign_fragments(tranf)
   !  use fragments
     use basis
     use nodeinfo
     use sectors
     use bmpi_mod
     use localvectors
     use flagger
     implicit none
     integer :: is, ic, fg, n, localdim, dim_sum
     integer :: basesplit_ind, jztemp_scan, partemp_scan 
     integer :: jztemp, partemp 
     integer :: jz_i, jz_e, par_i, par_e
     integer :: s_shift, i, remaindr, fragnum
     
     integer (kind=basis_prec), pointer :: basestart_x(:), basestop_x(:)
     type (frag), pointer :: fragmentlist_x(:), fragmentlist_y(:)
     integer,pointer:: nsectors_x(:), nsectors_y(:)
     integer(8):: dimbasis_x
     integer :: nbreaks
     type(mastersect),pointer :: xsd_x(:), xsd_y(:)
     character,intent(in):: tranf
     integer :: ierr
     logical :: trackf
     integer,allocatable :: break_index(:) !tracks sector indices denoting fragment breaks
     integer,allocatable :: sec_tag(:) ! used for dividing target basis sectors
     allocate(break_index(nfragments-1))
     

     select case (tranf)
      !'forward' projection vec1 < vec2
      !projecting from smaller basis into larger
       case('F') 
       nsectors_x => nsectors_f
       nsectors_y => nsectors_i
       basestart_x => basestart_f
       basestop_x  => basestop_f
       fragmentlist_x => fragmentlist_f
       fragmentlist_y => fragmentlist_i
       xsd_x => xsd_f
       xsd_y => xsd_i
       dimbasis_x = dimbasis_f
  
       case('R')
      !'Reverse' projection vec1 > vec2
      !projecting from larger basis to smaller
        nsectors_x => nsectors_i
        nsectors_y => nsectors_f
        basestart_x => basestart_i
        basestop_x  => basestop_i
        fragmentlist_x => fragmentlist_i
        fragmentlist_y => fragmentlist_f
        xsd_x => xsd_i
        xsd_y => xsd_f
        dimbasis_x = dimbasis_i
       case default
       
       print*,' WRONG CHOICE OF TRANSFORMATION ',tranf
       stop
       
       end select
       n = nfragments
       nbreaks = 0
       remaindr = 0
       dim_sum = 0
       localdim = 0
       basesplit_ind = 0

       fg = 0
       do is = 1, nsectors_x(1)
         !try to keep fragments around the size of user described maxfragmentsize
         if(n .eq. 0 ) n=1
         basesplit_ind = abs((dimbasis_x - remaindr)/(n-nbreaks))! approximate fragment break
         !if(maxfragmentsize .le. basesplit_ind) basesplit_ind = maxfragmentsize
         localdim = xsd_x(1)%sector(is)%basisend - xsd_x(1)%sector(is)%basisstart + 1
         !print*,"localdim", localdim
            dim_sum = dim_sum + localdim
            !record fragment sector start stop
            if (dim_sum .ge. (remaindr+basesplit_ind)) then ! remaindr holds dim_sum from previous iteration
               nbreaks = nbreaks + 1
               break_index(nbreaks) = is
               remaindr = dim_sum
		if(iproc==0)then
                  print*,"---------------------------------------------------------------------"
                  print*,"basesplit_ind",basesplit_ind,"localdim",localdim
                  print*,"dim_sum",dim_sum,"r+basesplit_ind",remaindr+basesplit_ind
                  print*,"---------------------------------------------------------------------"
               endif 
            endif
            if(nbreaks .ge. nfragments-1) exit
         enddo
	if(iproc == 0) then
         	print*,"fragment sector breaks() on indices: "
         	do i = 1,nfragments-1
            		print*,break_index(i)
         	enddo 
         endif
         !check to make sure sectors with different w are all on the same fragment.
         do i = 1, nbreaks
            s_shift = 0
            if(break_index(i)+1 .le. nsectors_x(1)-1 ) then
               do ic = break_index(i)+1, nsectors_x(1)-1
                  jztemp = xsd_x(1)%sector(ic)%jzX 
                  partemp = xsd_x(1)%sector(ic)%parX
                  if(iproc == 0) then 
                  print*,"jztemp",jztemp,"partemp",partemp,"xsd_x(1)%sector(i)%jzX",xsd_x(1)%sector(i)%jzX,&
                  &"xsd_x(1)%sector(i)%parX",xsd_x(1)%sector(i)%parX,"i",i,"ic",ic
                  endif
                  if((jztemp .eq. xsd_x(1)%sector(break_index(i))%jzX).and.(partemp .eq. xsd_x(1)%sector(break_index(i))%parX))then
                     s_shift = s_shift + 1
                  endif
               enddo
               break_index(i) = break_index(i) + s_shift
            else
               exit
            endif
         enddo
	if(iproc == 0) then
         print*,"fragment sector breaks() after shift on indices: "
         do i = 1,nbreaks
            print*,break_index(i)
         enddo 
	endif
	
	!now set fragment sector start/stops based on larger basis break indices
	do fg = 1, nfragments
		if( fg .eq. 1) then 
		   fragmentlist_x(fg)%ssectorstart = 1
         fragmentlist_x(fg)%ssectorend  = break_index(fg)
		else if(fg .eq. nfragments) then
		   fragmentlist_x(fg)%ssectorstart = fragmentlist_x(fg-1)%ssectorend + 1
         fragmentlist_x(fg)%ssectorend  = nsectors_x(1)
		else
		   fragmentlist_x(fg)%ssectorstart = fragmentlist_x(fg-1)%ssectorend + 1
         fragmentlist_x(fg)%ssectorend  = break_index(fg)
		endif 
	enddo

   !now divide the target basis sectors onto fragments based on the par, Jz 
   allocate(sec_tag(nsectors_y(1)))
   break_index= 0
   nbreaks = 0
   do fg = 1, nfragments
      !range of par and Jz allowed on fragment
      jz_i = xsd_x(1)%sector(fragmentlist_x(fg)%ssectorstart)%jzX
      par_i = xsd_x(1)%sector(fragmentlist_x(fg)%ssectorstart)%parX
      jz_e = xsd_x(1)%sector(fragmentlist_x(fg)%ssectorend)%jzX
      par_e = xsd_x(1)%sector(fragmentlist_x(fg)%ssectorend)%parX
      
      trackf = .False.
      do ic = 1, nsectors_y(1)
         ! if basis is small enough then just assign all sectors to each fragment.
         if(.not. trackf) then
            if(ic .eq. nsectors_y(1)) then
            !   
               do n = 1, nfragments - 1
                  nbreaks = nbreaks + 1
                  break_index(n) = nsectors_y(1)
               enddo
               exit
            endif
         endif
         !check fragment Jz and par of each sector in the target space
         ! against allowed Jz and par of fragment. label each sector with its respective
         ! fragment number.
         jztemp = xsd_y(1)%sector(ic)%jzX
         partemp = xsd_y(1)%sector(ic)%parX
         if(iproc==0)print*,"ic",ic,"jz_e",jz_e,"par_e",par_e,"jztemp",jztemp,"partemp",partemp,"trackf",trackf

         if(trackf)then
            if((jztemp .eq. jz_e) .and.(partemp .eq. par_e))then
               cycle
            endif
         endif

         if(.not. trackf)then
            if((jztemp .eq. jz_e) .and.(partemp .eq. par_e))then
               trackf = .True.
               cycle
            endif
         endif
          
         if(trackf)then
            if((jztemp .gt. jz_e) .and.(partemp .ge. par_e))then
               trackf = .False.
               nbreaks = nbreaks + 1
               break_index(nbreaks) = ic-1
               exit
            endif
         endif
  
      enddo

      if(iproc==0)print*,"--------------------New fragment-----------------------------"
      if(nbreaks .ge. nfragments-1) exit
	enddo

   !determine fragment start stops from frag number in sec_tag
   !break_index = 0
   
   !do fg = 1, nfragments
   !   n = 0 
   !   do i = 1, nsectors_y(1)
   !      if(sec_tag(i) .eq. fg) then
   !         n = n + 1
   !      endif
   !   enddo
   !   if( fg .eq. 1) then 
   !      fragmentlist_y(fg)%ssectorstart = 1
   !      fragmentlist_y(fg)%ssectorend  = n
   !  else if(fg .eq. nfragments) then
   !     fragmentlist_y(fg)%ssectorstart = fragmentlist_y(fg-1)%ssectorend + 1
   !      fragmentlist_y(fg)%ssectorend  = nsectors_y(1)
   !  else
   !      fragmentlist_y(fg)%ssectorstart = fragmentlist_y(fg-1)%ssectorend + 1
   !      fragmentlist_y(fg)%ssectorend  = n
   !  endif 
   !enddo 
   if(iproc == 0) then
      print*,"fragment sector breaks() 2ND basis: "
      do i = 1,nfragments-1
         print*,"i",i,"break_index",break_index(i)
      enddo 
endif

   !now set fragment sector start/stops based on small basi sector break indices
	do fg = 1, nfragments
		if( fg .eq. 1) then 
		   fragmentlist_y(fg)%ssectorstart = 1
         fragmentlist_y(fg)%ssectorend  = break_index(fg)
		else if(fg .eq. nfragments) then
		   fragmentlist_y(fg)%ssectorstart = fragmentlist_y(fg-1)%ssectorend + 1
         fragmentlist_y(fg)%ssectorend  = nsectors_y(1)
		else
		   fragmentlist_y(fg)%ssectorstart = fragmentlist_y(fg-1)%ssectorend + 1
         fragmentlist_y(fg)%ssectorend  = break_index(fg)
		endif 
	enddo



   if(iproc == 0) print*,"finished assigning fragments..."
	call BMPI_BARRIER(icomm,ierr)
	return
   end subroutine assign_fragments




  ! subroutine create_fragments. This subroutine calculates the number of 
  !basis elements to be given to the set of remaining mpi ranks.
  ! Makes sure more even distribution so the nproc-1 fragment doesn't have too few basis elements.
  subroutine create_fragments(tranf)
   !  use fragments
     use basis
     use nodeinfo
     use sectors
     use bmpi_mod
     use localvectors
     use flagger
     implicit none
     integer :: is, ic, fg, localdim, dim_sum
     integer :: basesplit_ind, jztemp_i, partemp_i 
     integer :: jztemp, partemp, jztemp_f2, partemp_f2
     integer :: s_shift, i
     integer :: s, r, n !used to more evenly distribute elements
     !pointers 
     integer (kind=basis_prec), pointer :: basestart_x(:), basestop_x(:)
     type (frag), pointer :: fragmentlist_x(:), fragmentlist_y(:)
     integer,pointer:: nsectors_x(:)
     integer(8):: dimbasis_x
     type(mastersect),pointer :: xsd_x(:), xsd_y(:)
     character,intent(in):: tranf



     select case (tranf)
    !'Reverse' projection vec1 > vec2
    !projecting from larger basis to smaller

     case('R') 
     nsectors_x => nsectors_i
     basestart_x => basestart_i
     basestop_x  => basestop_i
     fragmentlist_x => fragmentlist_i
     xsd_x => xsd_i
     dimbasis_x = dimbasis_i

    !'forward' projection vec1 < vec2
    !projecting from smaller basis into larger
     case('F')
      nsectors_x => nsectors_f
      basestart_x => basestart_f
      basestop_x  => basestop_f
      fragmentlist_x => fragmentlist_f
      xsd_x => xsd_f
      dimbasis_x = dimbasis_f
     case default
     
     print*,' WRONG CHOICE OF TRANSFORMATION ',tranf
     stop
     
     end select

     dim_sum = 0
     s_shift = 0
     r = 0 ! running sum of sector dimensions
     s = 0 ! ranks with fragments already
     fg = 1
     n = 1

     do is = 1, nsectors_x(1)
      !try to keep fragments around the size of user described maxfragmentsize
      
      if(s .ne. nproc ) then
         basesplit_ind = abs((dimbasis_x - r)/(nproc - s))! approximate fragment break
         !if(maxfragmentsize .le. basesplit_ind) basesplit_ind = maxfragmentsize
      endif
      localdim = xsd_x(1)%sector(is)%basisend - xsd_x(1)%sector(is)%basisstart + 1
      !print*,"localdim", localdim
         dim_sum = dim_sum + localdim
         if(iproc==0)then
            !print*,"---------------------------------------------------------------------"
            !print*,"is",is,"s",s,"r",r,"basesplit_ind",basesplit_ind,"localdim",localdim
            !print*,"dim_sum",dim_sum,"r+basesplit_ind",r+basesplit_ind
            !print*,"---------------------------------------------------------------------"
         endif
         !record fragment sector start stop
         if (dim_sum .ge. (r+basesplit_ind)) then ! r holds dim_sum from previous iteration
            s_shift = 0
            !check next sectors to make sure we get all sectors of the same Jz and parity on the same fragment
            do ic = is+1, nsectors_x(1)-1
               jztemp_f2 = xsd_x(1)%sector(ic)%jzX 
               partemp_f2 = xsd_x(1)%sector(ic)%parX
               if((jztemp_f2 .eq. xsd_x(1)%sector(is)%jzX).and.(partemp_f2 .eq. xsd_x(1)%sector(is)%parX) ) then
                  s_shift = s_shift + 1
               end if
            enddo

            if(s .eq. 0) then ! iproc
               fragmentlist_x(fg)%ssectorstart = 1 
               fragmentlist_x(fg)%ssectorend  = is + s_shift
              

               print*,"here a ---"
               print*,"is:",is,"running sum:", dim_sum,"basesplit_ind",basesplit_ind,"shift",s_shift 
               print*,"sec start:",fragmentlist_x(fg)%ssectorstart,"sec end:",fragmentlist_x(fg)%ssectorend
               
               !increment
               s = s + 1
               n = n + 1
               fg = fg + 1
               r = dim_sum
               
               cycle
            else
               if(fg .eq. nproc)then
                  if(is .le. nsectors_x(1)) then
                     fragmentlist_x(fg)%ssectorstart = fragmentlist_x(fg-1)%ssectorend + 1
                     fragmentlist_x(fg)%ssectorend  = nsectors_x(1)
                     print*,"here b ---"
                     print*,"is:",is,"running sum:", dim_sum,"basesplit_ind",basesplit_ind,"shift",s_shift
                     print*,"sec start:",fragmentlist_x(fg)%ssectorstart,"sec end:",fragmentlist_x(fg)%ssectorend
                     !increment
                     s = s + 1
                     n = n + 1
                     fg = fg + 1
                     r = dim_sum
                     
                     cycle
                  endif
               else
                  s_shift = 0
                  do ic = is+1, nsectors_x(1)-1
                     jztemp_f2 = xsd_x(1)%sector(ic)%jzX 
                     partemp_f2 = xsd_x(1)%sector(ic)%parX
                     if((jztemp_f2 .eq. xsd_x(1)%sector(is)%jzX).and.(partemp_f2 .eq. xsd_x(1)%sector(is)%parX) ) then
                        s_shift = s_shift + 1
                     end if
                  enddo
                  if((fragmentlist_x(fg-1)%ssectorend + 1) .ge. is ) cycle !in case at boundary of Jz/M par and shift < end sector for frag
                     fragmentlist_x(fg)%ssectorstart = fragmentlist_x(fg-1)%ssectorend + 1
                     fragmentlist_x(fg)%ssectorend  = is + s_shift
                     print*,"here c ---"
                     print*,"is:",is,"running sum:", dim_sum,"basesplit_ind",basesplit_ind,"shift",s_shift
                     print*,"sec start:",fragmentlist_x(fg)%ssectorstart,"sec end:",fragmentlist_x(fg)%ssectorend


                  !increment
                  s = s + 1
                  n = n + 1
                  fg = fg + 1
                  r = dim_sum
                  
                  cycle
               endif
               
   
            end if 
            
            !n = n + 1
            !fg = fg + 1
            
            
         end if
      enddo

      !now switch pointers to sort the other basis sectors to the right fragment.
      select case (tranf)
          case('R') 
          nsectors_x => nsectors_f
          basestart_x => basestart_f
          basestop_x  => basestop_f
          fragmentlist_x => fragmentlist_f
          fragmentlist_y => fragmentlist_i
          xsd_x => xsd_f
          xsd_y => xsd_i
         !'forward' projection vec1 < vec2
         !projecting from smaller basis into larger
          case('F')
           nsectors_x => nsectors_i
           basestart_x => basestart_i
           basestop_x  => basestop_i
           fragmentlist_x => fragmentlist_i
           xsd_x => xsd_i
           xsd_y => xsd_f
          case default
          
          print*,'Problem in fragementer sorting 2nd basis sectors',tranf
          stop
          
          end select
      do fg = 1,nfragments !iterate over each fragment
         do is = 1, nsectors_x(1)
            !organization...
            jztemp_i = xsd_x(1)%sector(is)%jzX   
            partemp_i = xsd_x(1)%sector(is)%parX
            !jztemp_f1 = xsd_f(1)%sector(fragmentlist_f(fg)%ssectorstart)%jzX 
            !partemp_f1 = xsd_f(1)%sector(fragmentlist_f(fg)%ssectorstart)%parX
            jztemp_f2 = xsd_y(1)%sector(fragmentlist_y(fg)%ssectorend)%jzX 
            partemp_f2 = xsd_y(1)%sector(fragmentlist_y(fg)%ssectorend)%parX
            print*,"--------------------------------------------------------"
            print*,"fg:",fg,"is:",is
            print*,"jztemp_i",jztemp_i,"partemp_i",partemp_i,"jztemp_f2", jztemp_f2,"partemp_f2",partemp_f2
            print*,"--------------------------------------------------------"

            if(fg .eq. 1) then ! frag 1 setup
               
               !if((partemp_i .eq. partemp_f2) .and. (jztemp_i .eq. jztemp_f2 ) .and. (is .ne. nsectors_x(1)) ) then
               if((partemp_i .eq. partemp_f2) .and. (jztemp_i .eq. jztemp_f2 ) .and. (is .ne. nsectors_x(1)) ) then
                  !check next sectors to make sure we get all sectors of the same Jz and parity on the same fragment
                  if(is .ne. nsectors_x(1))then
                     s_shift = 0
                     do ic = is+1, nsectors_x(1)
                        jztemp = xsd_x(1)%sector(ic)%jzX 
                        partemp = xsd_x(1)%sector(ic)%parX
                        if((jztemp .eq. jztemp_f2).and.(partemp .eq. partemp_f2) ) then
                           s_shift = s_shift + 1
                        end if
                     enddo
                  endif
                  fragmentlist_x(fg)%ssectorstart = 1
                  fragmentlist_x(fg)%ssectorend = is + s_shift
                  print*,"a"
                  print*,"sec start:",fragmentlist_x(fg)%ssectorstart,"sec end:",fragmentlist_x(fg)%ssectorend
                  exit
               else if((partemp_i .eq. partemp_f2) .and. (jztemp_i .gt. jztemp_f2) ) then
                 
                  fragmentlist_x(fg)%ssectorstart = 1
                  fragmentlist_x(fg)%ssectorend = is-1
                  print*,"b"
                  print*,"sec start:",fragmentlist_x(fg)%ssectorstart,"sec end:",fragmentlist_x(fg)%ssectorend
                  exit
               else if( is .eq. nsectors_x(1)) then
                  fragmentlist_x(fg)%ssectorstart = 1
                  fragmentlist_x(fg)%ssectorend = nsectors_x(1)
                  print*,"c"
                  print*,"sec start:",fragmentlist_x(fg)%ssectorstart,"sec end:",fragmentlist_x(fg)%ssectorend
                  exit
               endif
            
            else
               if((partemp_i .eq. partemp_f2) .and. (jztemp_i .eq. jztemp_f2) ) then
                  if(is .ne. nsectors_x(1))then
                     s_shift = 0
                     do ic = is+1, nsectors_x(1)
                        jztemp = xsd_x(1)%sector(ic)%jzX 
                        partemp = xsd_x(1)%sector(ic)%parX
                        if((jztemp .eq. jztemp_f2).and.(partemp .eq. partemp_f2) ) then
                           s_shift = s_shift + 1
                        end if
                     enddo
                  endif
                  !in case at boundary of Jz/M par and shift < end sector for frag
                  if((fragmentlist_x(fg-1)%ssectorend + 1) .ge. is ) cycle
                  if((is + s_shift) .ge. nsectors_x(1)) then
                     fragmentlist_x(fg)%ssectorstart = fragmentlist_x(fg-1)%ssectorend + 1
                     fragmentlist_x(fg)%ssectorend = is
                  else
                     fragmentlist_x(fg)%ssectorstart = fragmentlist_x(fg-1)%ssectorend + 1
                     fragmentlist_x(fg)%ssectorend = is + s_shift
                  end if
                  print*,"d"
                  print*,"shift",s_shift, "is", is
                  print*,"sec start:",fragmentlist_x(fg)%ssectorstart,"sec end:",fragmentlist_x(fg)%ssectorend
                  exit
               else if((partemp_i .eq. partemp_f2) .and. (jztemp_i .gt. jztemp_f2) .and. (is .ne. nsectors_x(1))) then 
                  fragmentlist_x(fg)%ssectorstart = fragmentlist_x(fg-1)%ssectorend + 1
                  fragmentlist_x(fg)%ssectorend = is-1
                  print*,"e"
                  print*,"sec start:",fragmentlist_x(fg)%ssectorstart,"sec end:",fragmentlist_x(fg)%ssectorend
                  exit
               else if( is .eq. nsectors_x(1)) then

                  if((fragmentlist_x(fg-1)%ssectorend + 1) .gt. nsectors_x(1)) then
                     !just put the entire vec1 on the node
                     fragmentlist_x(fg)%ssectorstart = 1
                     fragmentlist_x(fg)%ssectorend = nsectors_x(1)
                  else
                     fragmentlist_x(fg)%ssectorstart = fragmentlist_x(fg-1)%ssectorend + 1
                     fragmentlist_x(fg)%ssectorend = is
                  end if
                  print*,"f"
                  print*,"sec start:",fragmentlist_x(fg)%ssectorstart,"sec end:",fragmentlist_x(fg)%ssectorend
                  exit
               endif


            endif

         enddo
      enddo

      
     return
  
  end subroutine create_fragments
! 
  subroutine check_dist
   !  use fragments
     use basis
     use nodeinfo
     use sectors
     use bmpi_mod
     use localvectors
     use flagger
     implicit none
     integer :: fragdim_est, ierr
   if(nproc .gt. 1)then
      fragdim_est = dimbasis_i/nproc
      if(iproc == 0)then
         print*,"checking mpi distribution..."
         
         if(dimbasis_i .gt. dimbasis_f) then
            print*,"largest basis dimension(initital): ", dimbasis_i
            fragdim_est = dimbasis_i/nproc
         else 
            print*,"largest basis dimension(final): ", dimbasis_f
            fragdim_est = dimbasis_f/nproc
         endif


         print*,"number of requested mpi ranks is: ",nproc
         print*,"crudely distributing largest basis over all requested mpi ranks (dimbasis/nproc)"
         print*,"would require a fragment size of approximately:  ",fragdim_est
         print*,"the default min/max fragment sizes are: ", minpiece, maxfragmentsize_default
      endif 

         !recommend distribution
         if((fragdim_est .ge. minpiece) .and. (fragdim_est .le. maxfragmentsize_default)) then
            if(iproc == 0) print*," will try to create fragments using requested mpi ranks..." 
               break_vectors = .true.
         else if(fragdim_est .le. minpiece) then
            if(iproc == 0) print*,"the fragment dimension is too small for mpi, consider running in serial." 
            if(iproc == 0) print*,"STOPPING RUN"
               break_vectors = .false.
               call BMPI_ABORT(icomm,101,ierr)
               stop  
         else
            if(iproc == 0) print*,"the fragment dimension is large for the number mpi ranks, consider running with more processes."
            if(iproc == 0) print*,"STOPPING RUN"   
               break_vectors = .false.
               call BMPI_ABORT(icomm,101,ierr)
               stop  
         endif
   else
      print*,"using a single fragment..."
      
      if(nproc == 1) then
         print*,"you requested only one mpi process, use serial version."
         print*,"STOPPING RUN"
         call BMPI_ABORT(icomm,101,ierr)
         stop
      endif
      break_vectors = .false.
   endif
         
         
  end subroutine check_dist

  subroutine printfragsectorinfo
   !  use fragments
     use basis
     use nodeinfo
     use sectors
     use bmpi_mod
     use localvectors
     implicit none
     integer :: fg
     !real :: basesplit_ind
     print*,"total number of fragments created:",nfragments
     print*,"printing fragment sector information..."
     print*,"====================================================================================================================="
     do fg = 1, nfragments
         print*,"fragment:",fg,"int. sec. start:",fragmentlist_i(fg)%ssectorstart,"int. sec. end:",fragmentlist_i(fg)%ssectorend
         print*,"fragment:",fg,"final. sec. start:",fragmentlist_f(fg)%ssectorstart,"final. sec. end:",fragmentlist_f(fg)%ssectorend
      enddo
      print*,"===================================================================================================================="
   return
  end subroutine printfragsectorinfo

  subroutine printsectorinfo
   !  use fragments
     use basis
     use nodeinfo
     use sectors
     use bmpi_mod
     use localvectors
     implicit none
     integer :: fg, i
     !real :: basesplit_ind
     print*,"printing fragment sector information..."

     
     if (iproc == 0)then
      print*,"====================================================================================================================="
      print*,"printing intitial space sector information..."
        do i = 1, nsectors_i(1)
           print*, "i", i, "jzX", xsd_i(1)%sector(i)%jzX,"parx",xsd_i(1)%sector(i)%parx,"wX",xsd_i(1)%sector(i)%wX
        end do
      print*,"====================================================================================================================="

      print*,"printing final sector information..."
        do i = 1, nsectors_f(1)
         print*, "i", i, "jzX", xsd_f(1)%sector(i)%jzX,"parx",xsd_f(1)%sector(i)%parx,"wX",xsd_f(1)%sector(i)%wX
        end do
  
      print*,"====================================================================================================================="
   endif
   return
  end subroutine printsectorinfo


! This subroutine is not used for the mpi projection routine.
  subroutine vectorstorageceiling
   use nodeinfo
   use system_parameters
   use basis
   use menu_choices
   use flagger
   use nodeinfo
   use bmpi_mod
!   use fragments
   use io
   implicit none
   integer(8) :: minmaxfragsize  ! largest block found by splitting along haikus
   integer :: ierr
   character :: dummchar
   logical   :: interrflag  ! error in reading an integer

!...... parameters ask2break_vectors and maxfragmentsize_default found in module flagger..
!... MODIFIED IN 7.6.5 to ask for dummy fragmentsize even if not in MPI...

   if(ask2break_vectors)then  ! should always default to this...
      if(nproc > 1)then
         if(iproc==0)then
            print*,' '
            print*,' Enter desired limit on fragment size for breaking Lanczos vectors'
            print*,' To keep the memory on a single mpi node less than 1Gb keep max maxfragmentsize < 200000000'
            print*,' the maximum'
            print*,' Enter 0 to use default, maxfragsize =', maxfragmentsize_default
            print*,'The default minimum fragment size is, minpiece =', minpiece
            read(5,*,err=1111)maxfragmentsize
            if(maxfragmentsize .gt. maxfragmentsize_default ) then
               print*,"Fragment size too large, may run out of memory..."
               print*,"Setting maxfragment size to the default value..."
               maxfragmentsize = maxfragmentsize_default
            else if((maxfragmentsize .le. maxfragmentsize_default) .and. (maxfragmentsize .ge. minpiece )) then
               print*,"Deviating from default fragment size..."
               print*,"maxfragmentsize = ", maxfragmentsize
              
            else
               print*,"stopping run, check fragment size input..."
               print*,"maxfragmentsize = ", maxfragmentsize
               goto 1111
            endif
            
            
         
         else
            maxfragmentsize=maxfragmentsize_default
         endif

      else         
         maxfragmentsize=maxfragmentsize_default ! this isn't used for now in the non-mpi case
   end if
   call BMPI_BARRIER(icomm,ierr)
   call BMPI_Bcast(maxfragmentsize,1,0,icomm,ierr)
   !print*,"break_vectors_enabled", break_vectors_enabled
      if(break_vectors_enabled .and. dimbasis_f > maxfragmentsize)then
         if(nproc ==1) then
   !		  if(.not.test_fragments)then
               print*,' Only one process; cannot break vector '
               print*,' You may run out of memory '
               break_vectors=.false.
         else
            if(iproc==0)then
               print*,' '
               print*,' Basis dimension greater than limit of ',maxfragmentsize
               print*,' Breaking up into fragments '
               print*,' (NOTE: You can set parameter maxfragmentsize in file rhmodule_flags.f90 ) '
               print*,' '
               print*, "setting break_vectors True"
            end if
            break_vectors = .true.

         end if
      else if (break_vectors_enabled .and. dimbasis_f < maxfragmentsize) then
            break_vectors = .false.
            if(iproc == 0) print*, "setting break_vectors FALSE"
            if(iproc == 0) print*, "the fragment size is larger than the final basis"
      else
            break_vectors = .false.
      end if
   endif
   if(iproc == 0) print*,"max fragmentsize is", maxfragmentsize
   return
1111 continue
   if(iproc==0)then
	   print*,' ERROR you did not set the fragmentsize ERROR'
	   !write(resultfile,*)' ERROR you did not set the fragmentsize ERROR'
	   !write(logfile,*)' ERROR you did not set the fragmentsize ERROR'
	   !close(resultfile)
	   !close(logfile)
   endif
   call BMPI_ABORT(icomm,101,ierr)   
   stop
end subroutine vectorstorageceiling

!  subroutine fragmenter
!
!  master routine to break up Lanczos vectors into "fragments"
!  when necessary it breaks up sectors into "subsectors"
!  by breaking along the proton basis Slater Determinants
!
!  The "non-optimized" version simply breaks lanczos vectors by sectors,
!  because that is easiest to implement.
!
!  The "block" version breaks up sectors by proton haiku blocks, which is a natural division
!  and preserves most of the downstream routines intact; the divided sectors must be then redefined
!
!  first count up the number of break points
!
!  CALLED BY:
!    main routine
!
!  SUBROUTINES CALLED:
!    vectorstorageceiling   checks to see if vectors get broken up, and how
!    count_create_fragments_noopt
!        [this creates non-optimized fragments by simply bundling together sectors;
!         later will use optimized fragments with subroutine count_create_fragments ]
!       OR
!    count_create_fragments_blocks
!    reconstruct_sector
!    [map_sectors2noopt_fragments] 
!    defaultfragments      only one fragment
!    fragmentvectorsize    sets start and stop for basis in fragments
!    fragmentconjugatesectors  checks which neutron sectors belong (added 7.6.2)
!
subroutine fragmenter
   
   use flagger
   use basis
   use nodeinfo
!   use sectors
!   use menu_choices

   implicit none

   
   !call vectorstorageceiling  !  check to see if vectors get broken up into fragments
   
   call check_dist ! set break_vectors and recommend mpi dist
   print*,"break_vectors:", break_vectors
   if(iproc==0) print*,"finished checking mpi distribution..."
   if(break_vectors)then
	print*,"made it to here"
      call set_nfragments(nproc)
      if(iproc==0) print*,"finished setting nfragments..."
   !we divide the final basis into nproc fragments by sector
   !each fragment has some range [istart,istop] of contigous sectors assigned to it
      if(dimbasis_f .gt. dimbasis_i ) then
         call assign_fragments('F')
         !call create_fragments('F')
	 !call assign_fragments('F')
      else
          call assign_fragments('R')
         !call create_fragments('R')
      endif
      if(iproc==0) print*,"finished creating fragments..." 
      !print info
      if(iproc==0) call printsectorinfo
      if(iproc==0) call printfragsectorinfo
      
   else
      call set_nfragments(1)
   end if

   !set basis start/stop using fragment sector info
   call fragmentvectorsize
   if(iproc==0) print*,"finished assigning basis start stop for fragments..."
   return
end subroutine fragmenter

!===============================================================
!
!  VECTORSTORAGECEILING
!
!  computes a ceiling on storing pieces of vectors, if any
!  stops run if number of mpi procs requested is too large for the fragment size
!
!  initiated 8/2011 by CWJ @ SDSU
!  modified 2/1/22 by RMZ @ SDSU
!
!  CALLED BY:
!     fragmenter
!
!  CALLS:
!    fragmentsurveysays





end module fragments

