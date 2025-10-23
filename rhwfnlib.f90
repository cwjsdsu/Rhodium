! NOTE: reorganized in 7.6.6 by CWJ to group like subroutines together
!
!  modules with global data and subroutines to open/close/read/write wfn files
!

module wfn_mod
   use io
   use nodeinfo
   use fragments
   use bmpi_mod
   implicit none

   logical :: wfn_mpiio

   integer, parameter :: wfn_magic = 31415926
   integer, parameter :: wfn_vecmagic = 27182818
   integer, parameter :: wfn_version = 1
   integer, parameter :: wfn_after_header_word = 3
   ! hack - kind's do not have to be the number of 
   ! bytes.  They are symbolic labels for options
   ! could use bit_size to test.
   integer, parameter :: wfn_integer_kind = 4
   integer :: offset_i = 0
   integer :: offset_f = 0 ! store offset from initial file, and compare to final before wfn readins
   ! this is to attempt to get the correct wfn pos for readin when the intitial/final basis files
   ! have different sps. 
   logical:: diffPrincipleQn = .false. ! logical to track if principle Qn of basis files are different.
   ! remember locations in file for later update
   integer(kind=MPI_OFFSET_KIND) :: wfn_after_header_pos  ! location just after header
   
  logical :: wfnverbose = .false.  ! added in 0.8.7
   
contains
	
! CONTAINS THE FOLLOWING SUBROUTINES
!
!
!
!	

!-------------------- WAVEFUNCTION OPEN/CLOSE/READ/WRITE SUBROUTINES---------------
!
! Check if we should be doing MPI IO
subroutine wfn_checkmpiio()
   use localvectors
   implicit none
   if(iproc == 0) return ! ok
   if(isfragroot) return ! ok
   print *, "iproc=", iproc, "Should not be doing MPI IO"
   stop 1 ! has msg
end subroutine wfn_checkmpiio

!============== MASTER ROUTINES FOR WRITING WFN FILE(S) =================
!-------------------------------------------------------------
!   subroutine wrfn_writeeigenvec (formerly: writeeigenvec_p)
!
!  writes eigenvector to disk
!
!  INPUT
!       filenumber = which file to write to
!       fg = fragment held in this processor, frag1 or frag2
!      v(:) either (v1s:v1e) or (v2s:v2e)
!       i  = index of eigenvector
!       e  = energy of eigenvector
!       xj = j (as real number) of eigenvector
!       xt2 = T(T+1) [as real number] of eigenvector
!
!   CALLS:
!     wfn_seekvecfragpos
!
subroutine wfn_writeeigenvec(filenumber, fg, v, i,e,xj,xt2)  ! formerly writeeigenvec_p
!   use lanczos_info
   use basis
   use precisions
   use io
   use nodeinfo
   use localvectors
   use fragments
   implicit none

   character(1) :: vchar
   real(kind=lanc_prec), intent(in) :: v(:) ! either (v1s:v1e) or (v2s:v2e)
   real :: e,xj,xt2
   integer(4):: i,j
   integer :: filenumber
   integer :: fg
   integer (kind=MPI_OFFSET_KIND) :: offset
   integer :: status(MPI_STATUS_SIZE)
   integer:: ierr
   integer(kind=basis_prec) :: nw
   integer(kind=basis_prec) :: jl

   !print*,"writeout",writeout
   if(.not.writeout)return   
   if(wfn_mpiio) then
      if(.not. isfragroot) return
!---------------- move to start of vec/fragment
      !if(iproc==0)print*, " made it to a"
      call wfn_seekvecfragpos(filenumber, i, fg)
      !if(iproc==0)print*, " made it to b"
      call BMPI_FILE_get_position(filenumber, offset, ierr)
      !if(iproc==0)print*, " made it to c"
    	! start writing
      if(fg == 1) then
         ! first fragment writes the overhead
         call BMPI_FILE_WRITE(filenumber, wfn_vecmagic, 1, status, ierr)
         call BMPI_FILE_WRITE(filenumber, i, 1, status, ierr)
         call BMPI_FILE_WRITE(filenumber, e, 1, status, ierr)
         call BMPI_FILE_WRITE(filenumber, xj, 1, status, ierr)
         call BMPI_FILE_WRITE(filenumber, xt2, 1,  status, ierr)
      end if
      nw = basestop_f(fg)- basestart_f(fg) + 1;
      call BMPI_FILE_WRITE(filenumber, v, nw, status, ierr)
      call BMPI_FILE_GET_POSITION(filenumber, offset, ierr)
      !if(iproc==0)print*, " made it to d"
   else 
      if(iproc/=0) return
      inquire(filenumber, pos=offset)
      write(filenumber) wfn_vecmagic,i,e,xj,xt2
      write(filenumber) (v(jl),jl=1,dimbasis_f)
      inquire(filenumber, pos=offset)
   end if
   !print*,"iproc finished writing wfn=", iproc
   return
end subroutine wfn_writeeigenvec  ! formerly:  writeeigenvec_p

!=========================================================
! ====  MASTER ROUTINE(S) FOR READING WFN FILE(S) ========
!=========================================================
!   subroutine wfn_readeigenvec
!
!  reads eigenvector from disk
!
!  always reads into a vector formatted like vec1.
!
!  INPUT:
!       filenumber = which unit to fetch vector from
!       fg = fragment held in this processor, frag1 or frag2
!       fc = communicator for update of other nodes belonging to fragment fg
!            either fcomm1 or fcomm2
!       i  = index of eigenvector - checks that we get the right one
!  OUTPUT
!    v  ! Either v1s:v1e, or v2s:v2e
!       e  = energy of eigenvector
!       xj = j (as real number) of eigenvector
!
!  SUBROUTINES CALLED:
!   wfn_seekvecfragpos
!
!
! MODIFIED FROM BIGSTICK ROUTINE
subroutine wfn_readeigenvec_select(whichfile,filenumber, fg, fc, v, i, e, xj, xt2)
!	use lanczos_info
	use basis
	use precisions
	use io
	use nodeinfo
	use localvectors
   use mod_reorthog
   use bmpi_mod
	implicit none
	! arguments
	character,intent(in):: whichfile
	integer, intent(in) :: filenumber
	integer, intent(in) :: fc ! communicator for vec1/2 slicing
	integer, intent(in) :: fg ! fragment this node works on
	integer, intent(in) :: i
	real, intent(out) :: e,xj,xt2
   ! we have an interface, so should pick up bounds of vector v
   real (kind=lanc_prec) :: v(:)  ! Either v1s:v1e, or v2s:v2e
	! locals
	real(kind=4) :: ejtbuf(3) ! read/bcast buffer for e, xj, xt2
	integer(4) :: j
	integer(kind=basis_prec) :: jl,dimbasis_local
	integer :: ierr
	integer :: status(MPI_STATUS_SIZE)
   integer(4) :: vmagic, vnum
   integer :: rootproc
   character(len=*), parameter :: msg_wrongvec = "wfn_readeigenvec_select: read wrong vector"

   
   
   if(nproc .gt. 1 ) then
      select case(whichfile)
      case('I')
         dimbasis_local = basestop_i(fg)-basestart_i(fg) + 1
      !dimbasis_local = dimbasis_i
         
      case('F')

      dimbasis_local = basestop_f(fg)-basestart_f(fg) + 1 
      
      case default
      
      print*,' error in using whichfile ',whichfile
      stop
      end select
   else
      select case(whichfile)
      case('I')

      dimbasis_local = dimbasis_i
         
      case('F')

      dimbasis_local = dimbasis_f

      case default
      
      print*,' error in using whichfile ',whichfile
      stop
   end select

   endif

	if(wfn_mpiio) then
      isfragroot = .true.
      if(isfragroot) then
			! move to start of vec/fragment
			call wfn_seekvecfragpos_select(whichfile,filenumber, i, fg)
         
			! start reading
         if(fg == 1) then
            !fragment 1 owns the vector overhead
            call BMPI_FILE_READ(filenumber, vmagic, 1, status, ierr)
           
            call BMPI_FILE_READ(filenumber, vnum, 1, status, ierr)
    
            if(i /= vnum) then
				   call BMPI_ERR(msg_wrongvec)
              
			   end if 
            
            call BMPI_FILE_READ(filenumber, ejtbuf, 3, status, ierr)
            
         end if
         ! read the vector

         call BMPI_FILE_READ(filenumber, v, size(v,1,basis_prec), status, ierr)
		end if
      
      ! All procs have to use the fg=1, isfragroot process
      ! as rootproc for BCAST:   br_rank_of_frag1_root
      rootproc = br_rank_of_frag1_root
      !print*,"iproc-->",iproc, "rootproc", rootproc, 'v', v
      ! note size with dflt kind here is ok, we are only sending fragment
   
      call BMPI_BCAST(v, size(v), rootproc, fc, ierr)
 
   else
		if(iproc==0)then
		    ! no seeking needed, we read the whole thing
         call wfn_seekvecfragpos_select(whichfile,filenumber, i, fg)
         
		    ! read i, e, xj, xt2 from file
         
			read(filenumber) vmagic, vnum, (ejtbuf(j), j = 1,3)
         
         if(i /= vnum) then
            if(nproc > 1) call BMPI_ERR(msg_wrongvec)
            print *, msg_wrongvec
            print*,"i-->",i, "vnum",vnum
            stop 1 ! has msg
         end if
			read(filenumber)(v(jl),jl=1,dimbasis_local)
		end if
      ! BMPI_BCAST does chunking
      rootproc = 0
      call BMPI_BCAST(v, size(v), rootproc, fc, ierr)
	end if
	! update all nodes with vec properties . This works  with or without fragments.
   if(iproc == rootproc) then
      if(wfnverbose)print*, "vmagic", "wfn_vecmagic", vmagic, wfn_vecmagic
      if(vmagic /= wfn_vecmagic) then
         print *, "vector", i, " is missing magic number"
         stop 1 ! has msg
      end if
   end if
   !print*, "made it here 9"
	call BMPI_BCAST(ejtbuf, size(ejtbuf), rootproc, icomm, ierr)
   !print*, "made it here 10"
	e = ejtbuf(1)
	xj = ejtbuf(2)
	xt2 = ejtbuf(3)
   !print*, "made it here 11"
	return
end subroutine wfn_readeigenvec_select

! ============================================================
!   HEADER WRITE/READ ROUTINES
! ============================================================
!
!  writes header to .wfn file with information about the basis
!
!  ONLY WRITES TO 'FINAL' BASIS 
!
subroutine write_wfn_header(filenumber)
  use system_parameters
  use sporbit
  use spstate
  use w_info
  use basis
  implicit none
  integer :: ierr
  integer :: filenumber
  integer :: it,n, dummy, v
  integer (kind=MPI_OFFSET_KIND) :: offset, dimbasispos
  type (spst),pointer ::xspsqn(:)
	  

  if(.not.writeout) return
  dimbasispos = 0
  ! KSM:  should really use root process of frag1 with mpi, which may not be iproc == 0
  if(iproc==0) then
     !
     ! Write file type magic number and
     ! version number
     ! print *, "writing from iproc=0"

     call wfn_write_int4(filenumber, wfn_magic)   ! ids file

     call wfn_write_int4(filenumber, wfn_version) ! version number
     ! save position of offset to after header
     call wfn_write_int4(filenumber, 0) ! write over later
     !
     ! Write some zeroed space at the end for offsets
     ! The first offset will be to after the header
     ! write offset past header
     dummy = 0
     do n = 1, 5
       call wfn_write_int4(filenumber, dummy)
     end do
 
     !------------- WRITE INFO ON VALENCE PARTICLES
     !              PAY ATTENTION TO P-H CONJUGATION
     !        
     if(phconj_f(1))np_f(1) = -np_f(1)
     if(phconj_f(2))np_f(2) = -np_f(2)
     ! write(filenumber)np(1),np(2)
     call wfn_write_int4(filenumber, np_f(1))
     call wfn_write_int4(filenumber, np_f(2))
     if(phconj_f(1))np_f(1) = -np_f(1)
     if(phconj_f(2))np_f(2) = -np_f(2)
 
     !------------ WRITE OUT INFORMATION ON S.P. ORBITS 
     !  write(filenumber)isoflag
     !  write(filenumber)numorb(1),numorb(2)
     call wfn_write_logical(filenumber, isoflag)
     call wfn_write_int4(filenumber, numorb_f(1))
     call wfn_write_int4(filenumber, numorb_f(2))
     do it = 1,2
       do n = 1,numorb_f(it)
       	 ! all fields are integers
         !  write(filenumber)orbqn(it,n)%nr,orbqn(it,n)%j,orbqn(it,n)%l,  & 
         !            orbqn(it,n)%par, orbqn(it,n)%w
         call wfn_write_int4(filenumber, orbqn_f(it,n)%nr)
         call wfn_write_int4(filenumber, orbqn_f(it,n)%j)
         call wfn_write_int4(filenumber, orbqn_f(it,n)%l)
         call wfn_write_int4(filenumber, orbqn_f(it,n)%par)
         call wfn_write_int4(filenumber, orbqn_f(it,n)%w)
       enddo
     enddo
 
     !-------------- WRITE INFO ON S.P. STATES   -- CWJ 11/09  
     !               fixes a subtle bug when reading in without w-cuts
     !  write(filenumber)nsps(1),nsps(2)
     call wfn_write_int4(filenumber, nsps_f(1))
     call wfn_write_int4(filenumber, nsps_f(2))
      do it = 1,2
        do n = 1,nsps_f(it)
          call wfn_write_int4(filenumber, spsqn_f(it,n)%nr)
          call wfn_write_int4(filenumber, spsqn_f(it,n)%j)
          call wfn_write_int4(filenumber, spsqn_f(it,n)%m)
          call wfn_write_int4(filenumber, spsqn_f(it,n)%l)
          call wfn_write_int4(filenumber, spsqn_f(it,n)%w)
          call wfn_write_int4(filenumber, spsqn_f(it,n)%par)
          call wfn_write_int4(filenumber, spsqn_f(it,n)%orb)
          call wfn_write_int4(filenumber, spsqn_f(it,n)%group)
        enddo
     enddo
 
     !------------- WRITE OUT INFORMATION ON JZ, PARITY, W
     call wfn_write_int4(filenumber, jz_f)
     v = ICHAR(cparity_f)
     call wfn_write_int4(filenumber, v)
     call wfn_write_int4(filenumber, maxWtot_f)
     call wfn_write_logical(filenumber, allsameparity_f)
     call wfn_write_logical(filenumber, allsameW_f)
     call wfn_write_logical(filenumber, spinless)
 
     !------------ WRITE OUT BASIS DIMENSION

     if(nfragments == 1) inquire(filenumber, POS=dimbasispos)
     call wfn_write_bkint(filenumber, dimbasis_f)
  end if
  ! All processes reach here.  Broadcast offset if we are using
  ! MPI
  if(wfn_mpiio) then
    ! get offset
    offset = 13
    if(iproc==0) then
       call BMPI_FILE_get_position(filenumber, offset, ierr)
       offset = offset + 4  ! allocate space for nkeep between header and vectors
       dummy = int(offset,4)  ! use 4 bytes for offset.
       ! KSM - FIXME Something wrong here
       !print*,"dummy1", dummy
       call wfn_write_int4_at(filenumber, wfn_after_header_word*4, dummy)
       !print*,"dummy2", dummy
       call BMPI_FILE_seek(filenumber, offset, MPI_SEEK_SET, ierr)
       !print*,"offset 2", offset
       call wfn_write_int4(filenumber, dummy)
    end if
    ! sync all processes
    ! if(iproc == 0) print *, "dummy=", dummy
    call BMPI_BCAST(dummy, 1, 0, icomm, ierr)
    offset = dummy
    ! save position after header
    !wfn_after_header_pos = offset; ! don't update wfn_header_pos,
      if(wfnverbose) print *, "iproc=", iproc, ", position after header=", offset
    ! Set position for writing of nkeep
    call BMPI_FILE_SEEK(filenumber, offset-4, MPI_SEEK_SET, ierr)
  else if(iproc == 0) then
    ! no fragments, iproc==0 does all the work
    inquire(UNIT=filenumber, POS=offset)
    offset = offset + 4  ! allocate space for nkeep between header and vectors
    dummy = int(offset-1,4)  ! positions start at 1, MPI_positions start at 0
    v = wfn_after_header_word*4+1
    write(filenumber, pos=v) dummy
      if(wfnverbose) print *, "Wrote vec offset at ", v, ", value=", dummy
    ! move back to end
    write(filenumber, pos=dimbasispos) dimbasis_f
    ! write(filenumber) 1234
    !wfn_after_header_pos = offset;
    if(wfnverbose)print *, "no-frag, position after header=", offset
  end if
  return
end subroutine write_wfn_header
!===========================================================================
!
!
!  reads header from .wfn file with information about the basis
!
! CALLED BY :
!    menu
!
subroutine read_wfn_header_select(filenumber,verbose,compare,whichfile)
  use system_parameters
  use sporbit
  use spstate
  use w_info
  use basis
  use butil_mod
  use reporter
  implicit none
  integer:: filenumber
  logical:: verbose,dummyread
  logical :: compare  ! whether to compare with old results or not
  logical :: okay    ! returns TRUE if all in agreement
  character :: whichfile
  integer:: it,n, v
  integer:: ierr
  integer :: fpos
  integer :: aerr
  integer :: nprottmp,nneuttmp,jztmp,iparitytmp,maxWtotX
  logical :: allsameparityX,allsameWX
  character(len=1) :: cparitytmp
  integer :: numorbx(2),nspsx(2),npmaxx
  type (orb),pointer :: orbqnx(:,:),orbqny(:,:)
  type (spst),pointer :: spsqnx(:,:),spsqny(:,:)
  logical :: bothparitiestmp
  


  if(iproc == 0) then
     ! print *, "reading wfn header, filenumber=", filenumber
     ! read initial 8 word table
     call wfn_read_int4(filenumber, v)
	 if(wfnverbose)print*,' testing magic number ',v,wfn_magic
      if(wfnverbose .and. iproc == 0) print *, "read magic=", v
     if(v /= wfn_magic) then
        if(iproc == 0) print *,filenumber,v,wfn_magic, "WFN file is missing magic number"
        stop 1 ! has msg
     end if
     call wfn_read_int4(filenumber, v) ! version number
     ! print *, "read version=", v
     if(v /= wfn_version) then
        if(iproc == 0) print *, "WFN file is wrong version"
        stop 1 ! has msg
     end if
     call wfn_read_int4(filenumber, v) ! word not used yet
     call wfn_read_int4(filenumber, v) ! word 3 (starting with word 0)
      if(wfnverbose)print *, "read vec pos=", v
     if(wfnverbose)print*,"wfn_mpiio:", wfn_mpiio
     if(.not. wfn_mpiio) v = v + 1 ! NON MPI read/write starts at 1
     wfn_after_header_pos = v ! diff size
     do n=1, 4
       call wfn_read_int4(filenumber, v)
     end do
  !------------- READ INFO ON VALENCE PARTICLES
     ! read(filenumber)np(1),np(2)
     call wfn_read_int4(filenumber, nprottmp)
     call wfn_read_int4(filenumber, nneuttmp)
!		 if(nprottmp< 0 .or.nneuttmp < 0)then
!			 print*,' NEED TO FIX READING IN p-h conjugation '
!			 print*,nprottmp,nneuttmp
!		 end if
	  

     !------------ READ OUT INFORMATION ON S.P. ORBITS 
     ! read(filenumber)isoflag
     call wfn_read_logical(filenumber, isoflag)
     ! read(filenumber)numorb(1),numorb(2)
     call wfn_read_int4(filenumber, numorbx(1))   
     call wfn_read_int4(filenumber, numorbx(2))
  endif


  call BMPI_BCAST(wfn_after_header_pos, 1, 0, icomm, ierr)
  call BMPI_BCAST(nprottmp, 1, 0, icomm, ierr)
  call BMPI_BCAST(nneuttmp, 1, 0, icomm, ierr)
  call BMPI_BCAST(numorbx(1), 1, 0, icomm, ierr)
  call BMPI_BCAST(numorbx(2),1,0, icomm, ierr)
  call BMPI_BCAST(isoflag, 1, 0, icomm, ierr)
  
  if(compare)then
  
	  if(whichfile=='I')then
		  if(numorb_i(1)/=numorbx(1) .or. numorb_i(2)/=numorbx(2) )then
			  print*,' Mismatch in initial orbitals '
			  print*,' Expected from .bas file: ',numorb_i
			  print*,' Found in .wfn file: ',numorbx
			  stop
		  end if

	  else
		  if(numorb_f(1)/=numorbx(1) .or. numorb_f(2)/=numorbx(2) )then
			  print*,' Mismatch in final orbitals '
			  print*,' Expected from .bas file: ',numorb_f
			  print*,' Found in .wfn file: ',numorbx
			  stop
		  end if
	  end if  
else    ! don't compare, read in
  	if(whichfile=='I')then
		  numorb_i(1)=numorbx(1)
	  	numorb_i(2)=numorbx(2)
	else
	  	numorb_f(1)=numorbx(1)
	  	numorb_f(2)=numorbx(2)
  	end if
end if

 numorbmax = bmax(numorbx(1),numorbx(2))
  
!------------------ALLOCATE MEMORY------------------------------------------

if(compare)then
	
	if(associated(orbqnx))nullify(orbqnx)
	if(associated(orbqny))nullify(orbqny)

	allocate(orbqnx(2,numorbmax),orbqny(2,numorbmax))
	if(whichfile=='I')then
      orbqny => orbqn_i
	else
    	orbqny => orbqn_f
	end if

else

	if(whichfile=='I')then
  	  if(.not.allocated(orbqn_i)) then
     	 allocate ( orbqn_i(2,numorbmax), stat=aerr )
     	if(aerr /= 0) call memerror("read_wfn_header 1")
  	  end if
      orbqnx => orbqn_i
	else
    	if(.not.allocated(orbqn_f)) then
      	  allocate ( orbqn_f(2,numorbmax), stat=aerr )
       	  if(aerr /= 0) call memerror("read_wfn_header 2")
    	end if
    	orbqnx => orbqn_f
	end if
end if
  do it = 1,2
      do n = 1,numorbx(it)
        if(iproc==0) then
           ! read(filenumber)orbqn(it,n)%nr,orbqn(it,n)%j,orbqn(it,n)%l,  & 
           !           orbqn(it,n)%par, orbqn(it,n)%w 
    	     call wfn_read_int4(filenumber, orbqnx(it, n)%nr)
    	     call wfn_read_int4(filenumber, orbqnx(it, n)%j)
    	     call wfn_read_int4(filenumber, orbqnx(it, n)%l)
    	     call wfn_read_int4(filenumber, orbqnx(it, n)%par)
    	     call wfn_read_int4(filenumber, orbqnx(it, n)%w)
        end if  
    	! update all other nodes
        call BMPI_BCAST(orbqnx(it,n)%nr,1,0,icomm,ierr)
        call BMPI_BCAST(orbqnx(it,n)%j ,1,0,icomm,ierr)
        call BMPI_BCAST(orbqnx(it,n)%l ,1,0,icomm,ierr)
        call BMPI_BCAST(orbqnx(it,n)%par,1,0,icomm,ierr)
        call BMPI_BCAST(orbqnx(it,n)%w ,1,0,icomm,ierr)
      enddo

  enddo

  if(compare)then
	  do it = 1,2
	      do n = 1,numorbx(it)
			  if((orbqnx(it,n)%nr /= orbqny(it,n)%nr) .or.     &
			     (orbqnx(it,n)%j /= orbqny(it,n)%j)  .or.      &
			     (orbqnx(it,n)%l /= orbqny(it,n)%l ) .or.      &
			     (orbqnx(it,n)%par /= orbqny(it,n)%par ) .or.  &
			     (orbqnx(it,n)%w /= orbqny(it,n)%w ) )then
				 print*,' Mismatch in orbits ',it,whichfile,n
				 print*,' expected: ',orbqny(it,n)%nr,orbqny(it,n)%j,orbqny(it,n)%l,orbqny(it,n)%par,orbqny(it,n)%w
				 print*,' Found: ',orbqnx(it,n)%nr,orbqnx(it,n)%j,orbqnx(it,n)%l,orbqnx(it,n)%par,orbqnx(it,n)%w
				 stop
				 
			 end if
			  
		  end do ! n
		  
	  end do ! it
	  
  end if
  
  
  if(iproc==0)then
     ! read(filenumber)nsps(1),nsps(2)
     call wfn_read_int4(filenumber, nspsx(1))
     call wfn_read_int4(filenumber, nspsx(2))
  end if
  call BMPI_BCAST(nspsx,2,0,icomm,ierr)
  
  if(compare)then
	  if(whichfile=='I')then
		  if(nsps_i(1)/=nspsx(1).or. nsps_i(2)/=nspsx(2))then
			  print*,' Mismatch in initial nsps'
			  print*,' Expected from .bas file: ',nsps_i
			  print*,' Found in .wfn file: ',nspsx
			  stop 2	  
		  end if
		  
	  else
		  if(nsps_f(1)/=nspsx(1).or. nsps_f(2)/=nspsx(2))then
			  print*,' Mismatch in final nsps'
			  print*,' Expected from .bas file: ',nsps_f
			  print*,' Found in .wfn file: ',nspsx
			  stop 3		  
		  end if
		  
	  end if
	  
  else
  	if(whichfile=='I')then
	    nsps_i(1)=nspsx(1)
	  	nsps_i(2)=nspsx(2)
  	else
	    nsps_f(1)=nspsx(1)
	    nsps_f(2)=nspsx(2)
    end if
  end if
  
  if(compare)then
	  if(associated(spsqnx))nullify(spsqnx)
	  allocate(spsqnx(2,bmax(nspsx(1),nspsx(2))), stat=aerr)
    	if(whichfile=='I')then
  		    spsqny => spsqn_i
  	    else

  		    spsqny => spsqn_f
  	    end if

  else
	  
  	if(whichfile=='I')then
  		if(.not.allocated(spsqn_i)) then
    		 allocate(spsqn_i(2,bmax(nspsx(1),nspsx(2))), stat=aerr)
     		if(aerr /= 0) call memerror("read_wfn_header 2a")
  		end if
		spsqnx => spsqn_i
	else
  		if(.not.allocated(spsqn_f)) then
    		 allocate(spsqn_f(2,bmax(nspsx(1),nspsx(2))), stat=aerr)
     		 if(aerr /= 0) call memerror("read_wfn_header 2b")
  		 end if
		 spsqnx => spsqn_f
	 end if

end if
  
  do it = 1,2
     do n = 1,nspsx(it)
       if(iproc==0) then
    	  call wfn_read_int4(filenumber, spsqnx(it, n)%nr)
    	  call wfn_read_int4(filenumber, spsqnx(it, n)%j)
    	  call wfn_read_int4(filenumber, spsqnx(it, n)%m)
    	  call wfn_read_int4(filenumber, spsqnx(it, n)%l)
    	  call wfn_read_int4(filenumber, spsqnx(it, n)%w)
    	  call wfn_read_int4(filenumber, spsqnx(it, n)%par)
    	  call wfn_read_int4(filenumber, spsqnx(it, n)%orb)
    	  call wfn_read_int4(filenumber, spsqnx(it, n)%group)
       end if
       ! update other nodes
       call BMPI_BCAST(spsqnx(it,n)%nr,1,0,icomm,ierr)
       call BMPI_BCAST(spsqnx(it,n)%j ,1,0,icomm,ierr)
       call BMPI_BCAST(spsqnx(it,n)%m ,1,0,icomm,ierr)
       call BMPI_BCAST(spsqnx(it,n)%l ,1,0,icomm,ierr)
       call BMPI_BCAST(spsqnx(it,n)%w,1,0,icomm,ierr)
       call BMPI_BCAST(spsqnx(it,n)%par,1,0,icomm,ierr)
       call BMPI_BCAST(spsqnx(it,n)%orb,1,0,icomm,ierr)
       call BMPI_BCAST(spsqnx(it,n)%group,1,0,icomm,ierr)
    enddo

  enddo
  
  if(compare)then
	  do it = 1,2
	     do n = 1,nspsx(it)
			 if(spsqnx(it,n)%nr /= spsqny(it,n)%nr .or. & 
  			    spsqnx(it,n)%j /= spsqny(it,n)%j.or. & 
  			    spsqnx(it,n)%m /= spsqny(it,n)%m .or. & 
  			    spsqnx(it,n)%l /= spsqny(it,n)%l .or. & 
  			    spsqnx(it,n)%par /= spsqny(it,n)%par .or. & 
  			    spsqnx(it,n)%w /= spsqny(it,n)%w)then
			    print*,' Mismatch in s.p. states ',it,whichfile,n
			    print*,' expected: ',spsqny(it,n)%nr,spsqny(it,n)%j,spsqny(it,n)%m,spsqny(it,n)%par,spsqny(it,n)%w
			    print*,' Found: ',spsqnx(it,n)%nr,spsqnx(it,n)%j,spsqnx(it,n)%m,spsqnx(it,n)%par,spsqnx(it,n)%w
			    stop 4
			end if						 
		 end do
	 end do 
  end if  

!------- ALLOCATE MEMORY  -- CWJ 11/09
!------------- READ IN INFORMATION ON JZ, PARITY, W
  if(iproc==0) then
     ! read(filenumber)jz,cparity,maxWtot
     call wfn_read_int4(filenumber, jztmp)
     call wfn_read_int4(filenumber, n)
     cparitytmp = ACHAR(n) ! convert integer back to character
     call wfn_read_int4(filenumber, maxWtotX)
  end if
  call BMPI_BCAST(jztmp,1,0,icomm,ierr)
  call BMPI_BCAST(maxWtotX,1,0,icomm,ierr)
  call BMPI_BCAST(cparitytmp,1,0,icomm,ierr)
  select case (cparitytmp)
    case ('+')
	iparitytmp = 1
	bothparitiestmp = .false.
    case ('0')
	iparitytmp = 1
	bothparitiestmp = .true.
	   
    case ('-')
       iparitytmp = 2
   	bothparitiestmp = .false.
	   
    case default
       write(6,*)' Problem with parity ',cparitytmp
       stop
  end select
  
  
  if(compare)then
	  if(whichfile =='I')then
		  if(iparity_i/=iparitytmp .or. cparity_i/= cparitytmp)then
			  print*,' Mismatch in initial parity '
			  print*,' Expect ',iparity_i,cparity_i
			  print*,' Found ',iparitytmp,cparitytmp
		  end if
		  if(maxWtot_i /= maxWtotX)then
			  print*,' Mismatch in maxWtot '
			  print*,' Expect ',maxWtot_i
			  print*,' Found ',maxWtotX
		  end if

	  else
		  if(iparity_f/=iparitytmp .or. cparity_f/= cparitytmp)then
			  print*,' Mismatch in final parity '
			  print*,' Expect ',iparity_f,cparity_f
			  print*,' Found ',iparitytmp,cparitytmp
		  end if
		  if(maxWtot_f /= maxWtotX)then
			  print*,' Mismatch in maxWtot '
			  print*,' Expect ',maxWtot_f
			  print*,' Found ',maxWtotX
		  end if
		  
	  end if
	  
  else
	  
  	if(whichfile =='I')then
	  iparity_i = iparitytmp
	  cparity_i = cparitytmp
	  maxWtot_i = maxWtotX
	  bothparities_i = bothparitiestmp
  	else
	  iparity_f = iparitytmp
	  cparity_f = cparitytmp
	  maxWtot_f = maxWtotX
	  bothparities_f = bothparitiestmp
	  
  	end if
  
end if

  if(iproc==0) then 
     ! read(filenumber)allsameparity, allsameW,spinless
     call wfn_read_logical(filenumber, allsameparityX)
     call wfn_read_logical(filenumber, allsameWX)
     call wfn_read_logical(filenumber, spinless)
  end if
  call BMPI_BCAST(allsameparityX,1,0,icomm,ierr)
  call BMPI_BCAST(allsameWX,1,0,icomm,ierr)
  call BMPI_BCAST(spinless,1,0,icomm,ierr)
  
  if(compare)then
	  if(whichfile =='I')then
		  if(allsameparityX .neqv. allsameparity_i)then
			  print*,' Mismatch in initial allsameparity ',allsameparity_i,allsameparityX			  
		  end if
		  if(allsameWX .neqv. allsameW_i)then
			  print*,' Mismatch in initial allsameW ',allsameW_i,allsameWX
			  
		  end if
	  else
		  if(allsameWX .neqv. allsameW_f)then
			  print*,' Mismatch in final allsameW ',allsameW_f,allsameWX
			  
		  end if		  
		  if(allsameparityX .neqv. allsameparity_f )then
			  print*,' Mismatch in final allsameparity ',allsameparity_f,allsameparityX			  
		  end if
	  end if
	  
  else
	  if(whichfile =='I')then
		  allsameparity_i = allsameparityX
		  allsameW_i      = allsameWX
		  
	  else
		  allsameparity_f = allsameparityX
		  allsameW_f      = allsameWX		  
		  
	  end if
	  
  end if
!------------ WRITE OUT BASIS DIMENSION

  if(iproc==0)then
     ! read(filenumber)dimbasischeck
     call wfn_read_bkint(filenumber, dimbasischeck)
     print *, "dimbasischeck=", dimbasischeck
  end if
  call BMPI_BCAST(dimbasischeck,1,0,icomm,ierr)

  if(compare)then
	  
	  do it = 1,2 ! check for p-h conjugation
		  if(whichfile=='I')then
			  if(phconj_i(it))then
				  if(it==1 .and. nprottmp > 0)then
					  print*,' some problem with initial proton p-h conjugation '
					  print*,np_i(it),nprottmp
					  stop
				  end if
				  
				  if(it==1)nprottmp = -nprottmp
				  if(it==2 .and. nneuttmp > 0)then
					  print*,' some problem with initial neutron p-h conjugation '
					  print*,phconj_i(2),np_i(it),nneuttmp
					  stop
				  end if
				  if(it==2)nneuttmp = -nneuttmp
				  
			  end if
		  
		  else
			  if(phconj_f(it))then
				  if(it==1 .and. nprottmp > 0)then
					  print*,' some problem with final proton p-h conjugation '
					  print*,np_f(it),nprottmp
					  stop
				  end if
				  
				  if(it==1)nprottmp = -nprottmp
				  if(it==2 .and. nneuttmp > 0)then
					  print*,' some problem with final neutron p-h conjugation '
					  print*,phconj_f(2),np_f(it),nneuttmp
					  stop
				  end if
				  if(it==2)nneuttmp = -nneuttmp
				  
			  end if			  
			  
		  end if
		  
		  
	  end do
	  
	  if(whichfile=='I')then
		  if(np_i(1)/=nprottmp .or. np_i(2)/=nneuttmp)then
			  print*,' Mismatch in initial Z,N '
			  print*,' Expect : ',np_i
			  print*,' Found : ',nprottmp,nneuttmp
			  stop
		  end if
		  if(jztmp /= jz_i)then
			  print*,' Mismatch in initial M '
			  print*,' Expect : ',jz_i
			  print*,' Found : ',jztmp
			  stop
			  
		  endif 
		  
	  else
		  if(np_f(1)/=nprottmp .or. np_f(2)/=nneuttmp)then
			  print*,' Mismatch in final Z,N '
			  print*,' Expect : ',np_f
			  print*,' Found : ',nprottmp,nneuttmp
			  stop
		  end if	
		  if(jztmp /= jz_f)then
			  print*,' Mismatch in initial M '
			  print*,' Expect : ',jz_f
			  print*,' Found : ',jztmp
			  stop
			  
		  endif 	  
		  
	  end if
	  
  else
	  
  select case (whichfile)
  case('I')
  dimbasis_i = dimbasischeck
  np_i(1)=nprottmp
  np_i(2)=nneuttmp
  jz_i   =jztmp
  case('F')
  dimbasis_f = dimbasischeck
  
  np_f(1)=nprottmp
  np_f(2)=nneuttmp
  jz_f   =jztmp
  case default
  print*,' Error in whichfile ',whichfile
  stop
  end select


!--- ADDED 0.9.1 (June 2025 by CWJ) -- add in p-h conjugation 

do it = 1,2
    npmaxx = 0
  do n = 1,numorbx(it)
	  if(orbqnx(it,n)%w < 99)npmaxx=npmaxx + orbqnx(it,n)%j+1
  end do ! n
  
  if(whichfile == 'I')then
	  npmax_i(it)=npmaxx
	  phconj_i(it) = .false.
	  npeff_i(it) = np_i(it)
  else
	  npmax_f(it)=npmaxx
	  phconj_f(it) = .false.
	  npeff_f(it) = np_f(it)

  end if
  
  print*,' temp ',nprottmp,nneuttmp
  
  if(it == 1 .and. nprottmp < 0 .or. it==2 .and. nneuttmp < 0)then
	  if(whichfile == 'I')then
		  npeff_i(it)=np_i(it)+ npmaxx
		  np_i(it)=-np_i(it)
		  phconj_i(it)=.true.
		  
          if(iproc==0 .and. it==1)print*,np_i(it),' proton holes = ',npeff_i(it),' protons ' 
          if(iproc==0 .and. it==2)print*,np_i(it),' neutron holes = ',npeff_i(it),' neutrons ' 
		  
	  else
		  npeff_f(it)=np_f(it)+ npmaxx
		  np_f(it)=-np_f(it)
		  phconj_f(it)=.true.
		  
          if(iproc==0 .and. it==1)print*,np_f(it),' proton holes = ',npeff_f(it),' protons ' 
          if(iproc==0 .and. it==2)print*,np_f(it),' neutron holes = ',npeff_f(it),' neutrons ' 
		  
	  end if
	  
	  
  end if
  
end do ! it  
  

end if
!-------------------------------
  if(verbose .and. iproc==0)then
     print*,' Valence Z, N = ',nprottmp,nneuttmp
     print*,' Single particle space : '
     if(isoflag)then
         print*,'   N    L  2xJ '
         do n = 1,numorbx(1)
            write(6,'(3i5)')orbqnx(1,n)%nr,orbqnx(1,n)%l,orbqnx(1,n)%j
         end do
         print*,' '
         print*,' Total # of orbits = ',numorbx(1)
     else
         print*,' Oops not yet set up '
     endif
     if(allsameparityX)then
        print*,' 2 x Jz = ',jztmp
     else
        if(iparitytmp == 1)then
           print*,' 2 x Jz, parity = ',Jztmp,'+'
        else
           print*,' 2 x Jz, parity = ',Jztmp,'-'
        endif
     endif
     it = 1
     if(writeout)then
        if(isoflag)then
            write(resultfile,*) ' N       L       2xJ '
            do n = 1,numorbx(1)
               write(resultfile,*) orbqnx(it,n)%nr,orbqnx(it,n)%l,orbqnx(it,n)%j
            end do
        else
            print*,' Oops not yet set up '
        endif
     endif
  endif
  
  if(compare .and. iproc==0)print*,' .bas and .wfn headers agree '
  return
end subroutine read_wfn_header_select

!=====================================================================
!   subroutine (filenumber, i, fg)
!
!  seeks to start of fragment fg in vector i in the wfn file
!
!  note:  The first fragment carries the vector overhead.  This way
!         the file is independent of fragmentation choices
!
!  INPUT
!       filenumber = file id to operate on
!       i = vector number, starting with vector 1
!       fragment = piece of vector - see module fragments
! 
!  CALLED BY:
!   wfn_readeigenvec
!   wfn_writeeigenvec
!
!=====================================================================
subroutine wfn_seekvecfragpos(filenumber, i, fg)
   use fragments
	use localvectors
	implicit none
	integer, intent(in) :: filenumber, i, fg
	integer :: fi
	integer(kind=MPI_OFFSET_KIND) :: vs, pos, fragpos, veclen, nw, offset
	integer:: ierr
   integer :: overhead
   integer(kind=4) :: dummy4

   overhead = 5 * 4 ! (i, e, xj, xt2) + vec data 
   if(lanc_prec == 4) then
      vs = 4
   else
      vs = 8
   end if
	! should do this part once, could call from end of header read
   pos = 0
   fragpos = 0
   do fi = 1, nfragments
      if(fi == fg) fragpos = pos;
      nw = basestop_f(fi)- basestart_f(fi) + 1;
      pos = pos + nw*vs ! (i, e, xj, xt2) + vec data
      !print*,"fi",fi,"pos--", pos, "nw--", nw,"bs_f",basestart_f(fi),"be_f", basestop_f(fi) 
   end do
   veclen = pos + overhead  ! total vector size, all fragments
   ! compute offset to what we want to write
   offset = wfn_after_header_pos + (i-1)*veclen + fragpos
   if(fg /= 1) offset = offset + overhead
   !print*, "---nfragments",nfragments,"nw",nw,"pos",pos,"veclen",veclen,"offset",offset

   if(wfn_mpiio) then
   !if(nproc .gt. 1) then
      call BMPI_FILE_seek(filenumber, offset, MPI_SEEK_SET, ierr)
   else
      ! seek to start of vector

	  ! rewind necessary to clear state of filenumber before read when using ifort
	  ! otherwise, read will hang the second time around. - KSM
      rewind(filenumber)  

      read(filenumber, pos=offset-4) dummy4
   end if
	return
end subroutine wfn_seekvecfragpos
!=====================================================================
!   subroutine wfn_seekvecfragpos(filenumber, i, fg)
!
!  seeks to start of fragment fg in vector i in the wfn file
!
!  note:  The first fragment carries the vector overhead.  This way
!         the file is independent of fragmentation choices
!
!  INPUT
!       filenumber = file id to operate on
!       i = vector number, starting with vector 1
!       fragment = piece of vector - see module fragments
! 
!  CALLED BY:
!   wfn_readeigenvec
!   wfn_writeeigenvec
!
!=====================================================================
subroutine wfn_seekvecfragpos_select(whichfile,filenumber, i, fg)
   use fragments
	use localvectors
	implicit none
	character,intent(in):: whichfile
	integer, intent(in) :: filenumber, i, fg
	integer :: fi
	integer(kind=MPI_OFFSET_KIND) :: vs, pos, fragpos, veclen, nw, offset
	integer:: ierr
   integer :: overhead
   integer(kind=4) :: dummy4
   integer (kind=basis_prec), pointer :: basestart_x(:), basestop_x(:)
   
   select case (whichfile)
   
   case('I')
   
   basestart_x => basestart_i
   basestop_x  => basestop_i
   
   case('F')
   
   basestart_x => basestart_f
   basestop_x  => basestop_f
   
   case default
   
   print*,' WRONG CHOICE OF WHICHFILE ',whichfile
   stop
   
   end select

   overhead = 5 * 4 ! (i, e, xj, xt2) + vec data 
   if(lanc_prec == 4) then
      vs = 4
   else
      vs = 8
   end if
	! should do this part once, could call from end of header read
   pos = 0
   fragpos = 0
   do fi = 1, nfragments
      if(fi == fg) fragpos = pos;
      nw = basestop_x(fi)- basestart_x(fi) + 1;
      pos = pos + nw*vs ! (i, e, xj, xt2) + vec data
      if(wfnverbose)print*,"fi",fi,"pos--", pos, "nw--", nw,"bs_x",basestart_x(fi),"be_x", basestop_x(fi)
   end do
   veclen = pos + overhead  ! total vector size, all fragments
   ! compute offset to what we want to write
   offset = wfn_after_header_pos + (i-1)*veclen + fragpos
   if(fg /= 1) offset = offset + overhead
   if(wfnverbose)print*, "fragpos",fragpos, "wfn_after_header_pos",wfn_after_header_pos
   if(wfnverbose)print*, "offset",offset, "veclen", veclen, "overhead", overhead, "nw", nw

   if(wfn_mpiio) then

      call BMPI_FILE_seek(filenumber, offset, MPI_SEEK_SET, ierr)
  
   else
      ! seek to start of vector

	  ! rewind necessary to clear state of filenumber before read when using ifort
	  ! otherwise, read will hang the second time around. - KSM
      rewind(filenumber)  

      read(filenumber, pos=offset-4) dummy4
   end if
	return
end subroutine wfn_seekvecfragpos_select


!========================================================

! Wrapper for rewinding to the front of the file that
! handles the difference between MPI-IO and normal IO
subroutine wfn_rewind(fn)
   implicit none
   integer, intent(in) :: fn
   integer:: ierr

   if(wfn_mpiio) then
      call wfn_checkmpiio()
      ! pos 0 is beginning for MPI_FILE_seek
      call BMPI_FILE_seek(fn, int(0, MPI_OFFSET_KIND), MPI_SEEK_SET, ierr)
   else
      ! pos 1 is beginning for inquire/fseek/... on normal files
      ! all file I/O happens in iproc==0 for normal files
      if(iproc == 0) rewind(fn)
   end if
end subroutine wfn_rewind

!===========================================================
!  MISC READ/WRITE ROUTINES
!===========================================================
! sometimes the file is read manually and nkeep is
! discarded without overwriting the global
subroutine wfn_write_nkeep(nkeep)
   implicit none
   integer, intent(in):: nkeep
   integer :: ierr
   integer(kind=4) :: nkeep4
   if(writeout .and. write_wfn .and. iproc == 0) then
      nkeep4 = nkeep
      call wfn_write_int4(wfnfile_f, nkeep4)
   end if
end subroutine wfn_write_nkeep

subroutine wfn_read_nkeep(fn, nkeep)
   implicit none
   integer, intent(in) :: fn
   integer, intent(out) :: nkeep
   integer :: ierr
   integer(kind=4) :: nkeep4
   
   if(iproc == 0) then
      call wfn_read_int4(fn, nkeep4)
      nkeep = nkeep4
   end if
   if(nproc > 1) call BMPI_BCAST(nkeep,1, 0, icomm, ierr)
end subroutine wfn_read_nkeep

!===========================================================================

!
! Write one integer to the wfn file. 
! 
subroutine wfn_write_int4(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer, intent(in) :: v
   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
   integer(kind=4) :: data4

   data4 = v
   if(wfn_mpiio) then
      call wfn_checkmpiio()
      call BMPI_FILE_WRITE(fn, data4, 1, status, ierr)
   else
      write(fn) data4
   end if 
end subroutine wfn_write_int4
! write one integer8 to the file
subroutine wfn_write_int8(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(in) :: v
   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr

   if(wfn_mpiio) then
      call wfn_checkmpiio()
      call BMPI_FILE_WRITE(fn, v, 1, status, ierr)
   else
      write(fn) v
   end if 
end subroutine wfn_write_int8
! read one integer8 from the file
subroutine wfn_read_int8(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(out) :: v
   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr

   if(wfn_mpiio) then
      call wfn_checkmpiio()
      call BMPI_FILE_READ(fn, v, 1, status, ierr)
   else
      read(fn, err=911, iostat=ierr) v
      return
911   print *, "read error: ", ierr
!      call printstacktrace
      stop 1 ! has msg
   end if
end subroutine wfn_read_int8
!========================================================
! This writes a "basis kind" integer
! We just convert to 8 byte
subroutine wfn_write_bkint(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer (kind=basis_prec) :: v
   integer (kind=8) :: v8

   v8 = v
   call wfn_write_int8(fn, v8)
end subroutine wfn_write_bkint
!========================================================
! read one integer from the file
! just read as 8 byte integer and assign
subroutine wfn_read_bkint(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer (kind = basis_prec), intent(out) :: v
   integer (kind=8) :: v8

   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr

   call wfn_read_int8(fn, v8)
   v = v8
end subroutine wfn_read_bkint
!========================================================
! write integer at specified position.   Used to update header of
! file
subroutine wfn_write_int4_at(fn, fpos, v)
   implicit none
   integer, intent(in) :: fn, fpos
   integer(kind=4), intent(in) :: v
   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
   integer (kind=MPI_OFFSET_KIND) :: offset

   if(wfn_mpiio) then
      call wfn_checkmpiio()
      offset = fpos;
      call BMPI_FILE_seek(fn, offset, MPI_SEEK_SET, ierr)
      call BMPI_FILE_WRITE(fn, v, 1, status, ierr)
   else
      write(fn, pos=fpos) v
   end if 
end subroutine wfn_write_int4_at
!========================================================
! read one integer from the file
subroutine wfn_read_int4(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer, intent(out) :: v
   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
   integer(kind=4) :: data4

   if(wfn_mpiio) then
      call wfn_checkmpiio()
      !print*,"made it line 1206"
      call BMPI_FILE_READ(fn, data4, 1, status, ierr)
      v = data4
   else
      read(fn, err=911, iostat=ierr) data4
      v = data4
      return
911   print *, "read error: ", ierr
!      call printstacktrace
      stop 1 ! has msg
   end if
end subroutine wfn_read_int4
!========================================================
! read one integer at a position
subroutine wfn_read_int4_at(fn, fpos, v)
   implicit none
   integer, intent(in) :: fn, fpos
   integer, intent(out) :: v
   integer(kind=4) :: data4
   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
   integer (kind=MPI_OFFSET_KIND) :: offset

   if(wfn_mpiio) then
      call wfn_checkmpiio()
      offset = fpos;
      call BMPI_FILE_SEEK(fn, offset, MPI_SEEK_SET, ierr)
      call BMPI_FILE_READ(fn, data4, 1, status, ierr)
      v = data4
   else
      read(fn, pos=fpos, err=911, iostat=ierr) data4
      v = data4
      return
911   print *, "read error: ", ierr
!      call printstacktrace
      stop 1 ! has msg
   end if
end subroutine wfn_read_int4_at
!========================================================
! write a boolean into a file.  Uses same space as integer
subroutine wfn_write_logical(fn, lg)
   implicit none
   logical, intent(in) :: lg
   integer, intent(in) :: fn
   integer(kind=4) :: v
   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr

   v = 0
   if(lg) v = 1
   if(wfn_mpiio) then
      call wfn_checkmpiio()
      call BMPI_FILE_WRITE(fn, v, 1, status, ierr)
   else
      write(fn) v
   endif 
end subroutine wfn_write_logical
!========================================================
! read boolean from file
subroutine wfn_read_logical(fn, lg)
   implicit none
   integer, intent(in) :: fn
   logical, intent(out) :: lg
   integer(kind=4) :: v
   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr

   if(wfn_mpiio) then
      call wfn_checkmpiio()
      call BMPI_FILE_READ(fn, v, 1,  status, ierr)
      lg = v .ne. 0
   else
      read(fn, err=911, iostat=ierr) lg
      return
911   print *, "read error: ", ierr
!      call printstacktrace
      stop 1 ! has msg
   end if
end subroutine wfn_read_logical

!----------------------------------------------------
! SUBROUTINES TO OPEN FILES
!========================================================
! fragment mode:  Used to open an MPI file for writing out the
! eigenvectors (wavefunctions) for the first N states
! Only iproc==0 has afilename correct at entry so we
! have to broadcast it.
! Then we open the file.
!
subroutine wfn_wopen_file_frag(filenumber, afilename)
   implicit none
   integer, intent(inout) :: filenumber
   integer:: ierr
   character(len=*), intent(in) :: afilename
   character(len=1024) :: wfnfilename
   integer(kind=MPI_OFFSET_KIND), parameter :: zsize = 0

   ! distribute the file name
   wfnfilename = afilename
   call BMPI_BCAST(wfnfilename, LEN(wfnfilename) , 0, icomm, ierr)
   ! print *, "iproc=", iproc, ", wfnfilename=", trim(wfnfilename)

   call BMPI_FILE_OPEN(icomm, TRIM(wfnfilename), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, filenumber, ierr)
   if(ierr .ne. 0) then
      if(iproc == 0) print *, "Can't wopen .wfn file: ", TRIM(wfnfilename),ierr
      stop 1 ! has msg
   end if
   ! truncate the file
   call BMPI_FILE_SET_SIZE(filenumber, zsize, ierr)
end subroutine wfn_wopen_file_frag
!========================================================
! fragment mode:  Used to open the MPI file for the WFN file
subroutine wfn_ropen_file_frag(filenumber, afilename)
   implicit none
   integer, intent(inout) :: filenumber
   integer:: ierr
   character(len=*), intent(in) :: afilename
   character(len=1024) :: wfnfilename

   ! distribute the file name
   wfnfilename = afilename
   call BMPI_BCAST(wfnfilename, LEN(wfnfilename), 0, icomm, ierr)
   call BMPI_FILE_open(icomm, wfnfilename, MPI_MODE_RDONLY, MPI_INFO_NULL, filenumber, ierr)
   if(ierr .ne. 0) then
      if(iproc == 0) print *, "Can't ropen .wfn file: ", TRIM(wfnfilename),ierr
	  print*,MPI_MODE_RDONLY, MPI_INFO_NULL,filenumber
      stop 1 ! has msg
   end if

end subroutine wfn_ropen_file_frag

!=========================================================================
!   Open the wfn file for reading only.
!
!  modification of BIGSTICK routine wfn_ropen_file
!
subroutine wfn_ropen_file_select(filenumber,whichfile)
	use reporter
	use io
   implicit none
   integer, intent(inout) :: filenumber
   character :: whichfile  ! = 'I' or 'F' for initial, final
   integer :: ilast
   logical :: found
   character (len=1024) :: wfnfilename
   integer :: ierr
   character(55) :: infile

   found = .false.
   if(whichfile /='I' .and. whichfile/= 'F')then      ! ....... ERROR TRAP .....
	   print*,' Error in wfn_ropen_file_select ',whichfile
	   stop
   end if  ! ................END ERROR TRAP...........
   if(iproc == 0) then
      do while(.not.found)
         print*,' Enter input name of .wfn file '
         if(auto_input)then
            read(autoinputfile,'(a)')infile
         else
            read(5,'(a)')infile
            write(autoinputfile,'(a)')infile
         end if
         ilast = index(infile,' ')-1
         wfnfilename = infile(1:ilast)//'.wfn'
		 wfn_in_filename = infile(1:ilast)//'.wfn'
!...... MODIFIED in 7.6.5 to INQUIRE statement		 
         inquire(file=wfnfilename, &
             exist=found)
	     if(.not.found)then
				 print*,infile(1:ilast)//'.wfn does not exist '
		elseif(auto_input)then
		    print*,' reading from ',infile(1:ilast)//'.wfn  '
					 
		 end if

      end do
   end if
   
   select case(whichfile)
   
   case('I')
   
   wfnfilename_i = wfnfilename
   
   case('F')
   
   wfnfilename_f = wfnfilename
   
   case default
   
   print*,' Some error on whichfile ',whichfile
   stop
   
   end select
   ! all procs get here
   ! KSM - open on all mpi nodes
   ! The open above is just a trial run.
   wfn_mpiio = .false.
   if(nproc == 1) then
       open(unit=filenumber, file=wfnfilename, &
   	      action='READ', access='stream', form='unformatted', &
           status = 'old')	   
   else
      wfn_mpiio = .true.
      ! We successfully opened it on the root.  Now we close,
      ! broadcast the filename, and open on all processes.
      ! with MPI IO
      call wfn_ropen_file_frag(filenumber, wfnfilename)
   end if
   return
end subroutine wfn_ropen_file_select

!=====================================================================
!=========================================================================
!   Open the TEMP.wfn file for reading only.
!
!  modification of BIGSTICK routine wfn_ropen_file
!
subroutine wfn_ropen_file_TEMP(filenumber,whichfile)
	use reporter
	use io
   implicit none
   integer, intent(inout) :: filenumber
   character :: whichfile  ! = 'I' or 'F' for initial, final
   integer :: ilast
   logical :: found
   character (len=1024) :: wfnfilename
   integer :: ierr
   character(55) :: infile

   found = .false.
   if(whichfile /='I' .and. whichfile/= 'F')then      ! ....... ERROR TRAP .....
	   print*,' Error in wfn_ropen_file_select ',whichfile
	   stop
   end if  ! ................END ERROR TRAP...........
   if(iproc == 0) then
	   
         wfnfilename = 'TEMP.wfn'
		 wfn_in_filename = 'TEMP.wfn'
!...... MODIFIED in 7.6.5 to INQUIRE statement		 
         inquire(file=wfnfilename, &
             exist=found)
	     if(.not.found)then
				 print*,infile(1:ilast)//'.wfn does not exist '
				 stop
					 
		 end if

   end if
   
   select case(whichfile)
   
   case('I')
   
   wfnfilename_i = wfnfilename
   
   case('F')
   
   wfnfilename_f = wfnfilename
   
   case default
   
   print*,' Some error on whichfile ',whichfile
   stop
   
   end select
   ! all procs get here
   ! KSM - open on all mpi nodes
   ! The open above is just a trial run.
   wfn_mpiio = .false.
   if(nproc == 1) then
       open(unit=filenumber, file=wfnfilename, &
   	      action='READ', access='stream', form='unformatted', &
           status = 'old')	   
   else
      wfn_mpiio = .true.
      ! We successfully opened it on the root.  Now we close,
      ! broadcast the filename, and open on all processes.
      ! with MPI IO
      call wfn_ropen_file_frag(filenumber, wfnfilename)
   end if
   return
end subroutine wfn_ropen_file_TEMP

!=====================================================================
! Open wfn file for writing.
! If we have fragments, open with MPI call
! We still do trial open on iproc==0 to make
! sure file is writeable.
!
subroutine wfn_wopen_file(filenumber,basisfile)
   implicit none
   integer, intent(inout) :: filenumber
   logical, intent(in) :: basisfile   ! added in 7.7.2; this way we can reuse the same routines
   integer :: ilast
   character (len=1024) :: wfnfilename
   character (len=100) :: basewfnfilename
   logical :: found
   integer :: stat
   logical :: dolink
   
   ilast = index(outfile,' ')-1
   dolink = .false.
   if(basisfile)then
      basewfnfilename = outfile(1:ilast)//'.bas'
   else
      basewfnfilename = outfile(1:ilast)//'.wfn'
  end if

!   if(scratch_dir == '.') then
      wfnfilename = basewfnfilename
!   else
!      wfnfilename = TRIM(scratch_dir) // '/' // TRIM(basewfnfilename)
!      dolink = .true.
!   end if
   found = .false.
   if ( iproc == 0 .and. (write_wfn .or. basisfile)) then
      open(unit=filenumber, iostat=stat, file=wfnfilename, status='old')
      if(stat == 0) close(filenumber, status='delete')
!	  print*,basewfnfilename,index(basewfnfilename,' ')
!	  print*,wfnfilename,index(wfnfilename,' ')
      open(unit=filenumber, file=wfnfilename, &
           access='stream', form='unformatted', action='WRITE', & 
           status = 'unknown', err=101)
      found = .true.
      ! make link from working directory to file in scratch directory
      ! so the file is where we expect it
      if(dolink) then
         call system("ln -s " // TRIM(wfnfilename) // " .")
      end if
101   continue
      if(.not.found)then
    	 print *, "Can't open .wfn file for writing: ", wfnfilename
         print*," (in wfn_wopen_file) "
    	 stop 1 ! has msg
      end if
   end if
   ! all procs get here
   ! KSM - open on all mpi nodes
   wfn_mpiio = .false.
   if(nproc > 1) then
      wfn_mpiio = .true.
      ! We successfully opened it on the root.  Now we close
      ! and delete it, broadcast the filename, and open on all proceses.
      ! with MPI IO
      if(iproc == 0) close(unit = filenumber, status='delete')
      call wfn_wopen_file_frag(filenumber, wfnfilename)
   end if
   return
end subroutine wfn_wopen_file

!-------------------------------------
! routine to open a TEMP wfn file
subroutine wfn_wopen_file_TEMP(filenumber)
   implicit none
   integer, intent(inout) :: filenumber
   integer :: ilast
   character (len=1024) :: wfnfilename
   character (len=100) :: basewfnfilename
   logical :: found
   integer :: stat
   logical :: dolink
   
   ilast = index(outfile,' ')-1
   dolink = .false.

      basewfnfilename = 'TEMP.wfn'


!   if(scratch_dir == '.') then
      wfnfilename = basewfnfilename
!   else
!      wfnfilename = TRIM(scratch_dir) // '/' // TRIM(basewfnfilename)
!      dolink = .true.
!   end if
   found = .false.
   if ( iproc == 0 .and. (write_wfn)) then
      open(unit=filenumber, iostat=stat, file=wfnfilename, status='old')
      if(stat == 0) close(filenumber, status='delete')
!	  print*,basewfnfilename,index(basewfnfilename,' ')
!	  print*,wfnfilename,index(wfnfilename,' ')
      open(unit=filenumber, file=wfnfilename, &
           access='stream', form='unformatted', action='WRITE', & 
           status = 'unknown', err=101)
      found = .true.
      ! make link from working directory to file in scratch directory
      ! so the file is where we expect it
      if(dolink) then
         call system("ln -s " // TRIM(wfnfilename) // " .")
      end if
101   continue
      if(.not.found)then
    	 print *, "Can't open .wfn file for writing: ", wfnfilename
         print*," (in wfn_wopen_file) "
    	 stop 1 ! has msg
      end if
   end if
   ! all procs get here
   ! KSM - open on all mpi nodes
   wfn_mpiio = .false.
   if(nproc > 1) then
      wfn_mpiio = .true.
      ! We successfully opened it on the root.  Now we close
      ! and delete it, broadcast the filename, and open on all proceses.
      ! with MPI IO
      if(iproc == 0) close(unit = filenumber, status='delete')
      call wfn_wopen_file_frag(filenumber, wfnfilename)
   end if
   return
end subroutine wfn_wopen_file_TEMP

!-------------------------------------
! SUBROUTINE TO CLOSE FILE(S)
!========================================================
subroutine wfn_close_file(fn)
   implicit none
   integer, intent(inout) :: fn
   integer :: ierr

   if(fn <= 0) then
      print *, "wfn_close_file: wfnfile is not open"
      stop 1 ! has msg
   end if
   if(wfn_mpiio) then
      call BMPI_FILE_close(fn, ierr)
   else
      close(fn)
   end if
   fn = 0
end subroutine wfn_close_file


  
!=====================================================================
!============== SLATED FOR OBSOLESCENCE===============================
!=====================================================================
!--------------------------------------
!
! CALLS: wfn_get_position
!
subroutine wfn_wpos(fn, msg)
   implicit none
   integer :: fn
   character (len=*) :: msg
   integer :: fpos

   if(iproc /= 0) return
   call wfn_get_position(fn, fpos)
   write(6, "(A,A,O8)") msg, ", pos=", fpos
end subroutine wfn_wpos
!========================================================
!
! CALLED BY: wfn_wpos
!
subroutine wfn_get_position(filenumber, fpos)
   implicit none
   integer :: ierr
   integer:: filenumber, fpos
   integer (kind=MPI_OFFSET_KIND) :: offset

   if(wfn_mpiio) then
      call wfn_checkmpiio()
      call BMPI_FILE_GET_POSITION(filenumber, offset, ierr)
      fpos = int(offset,4) ! size mismatch
   else
      inquire(UNIT=filenumber, POS=fpos)
   end if
end subroutine wfn_get_position
!================================

end module wfn_mod



