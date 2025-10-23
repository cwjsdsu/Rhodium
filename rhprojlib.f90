
!
! routines for projecting wave functions into a different basis
!

module proj_mod
	
	use basis
	use basis_map
	implicit none

contains
! BOSS ROUTINE FOR PROJECTION FROM ONE SPACE INTO ANOTHER
!

	subroutine proj_boss
		use wfn_mod
		use bvectorlib_mod
		use fragments
		use localvectors
		use io
		use mod_reorthog
		
		implicit none
		integer(4) :: nkeep_i,nkeep_f
	    integer ikeep, nfr
	    integer i,j,n, tidx
	    real :: xe,xj,xt2
		real(8) ::dovlp
		integer :: ilast
		integer :: ierr, pi
		
!		!set wfn mpi flag
!		!if(nproc .ge. 1) then
!		!	wfn_mpiio = .true.
!		!endif
!		if(iproc == 0) print*, "Flag wfn_mpiio = ", wfn_mpiio
!	

		! set fragroot for setup local vectors.
		if(iproc == 0) isfragroot = .true.


		call phasechecker
		call basismapper(1)
		call basismapper(2)
		! added RMZ 2/1/22 -SDSU
		call fragmenter
		call BMPI_BARRIER(icomm,ierr)
		!Find the root of the first fragment.
		br_rank_of_frag1_root = -1
		do pi = 0, nproc-1
		   if(nodal(pi)%ifragment == 1 .and. nodal(pi)%ffragment == 1) then
			  br_rank_of_frag1_root = pi
			  exit
		   end if
		end do
		if(isfragroot .and. frag1 == 1) then
		   if(br_rank_of_frag1_root /= iproc) then
			  print *, "picked wrong root for read/write"
			  stop 23
		   end if
		end if

	    frag1= nodal(iproc)%ifragment
	    frag2= nodal(iproc)%ffragment
		!call BMPI_BARRIER(icomm,ierr)
		!if(iproc == 0)then 
		!	print*,' running setup local vectors...'
		!endif

		call setup_localvectors_both

		!if(iproc == 0)then 
		!	print*,'finished running setup local vectors...'
		!endif
		
!............ OPEN INITIAL .WFN FILE AND READ HEADER.....

		if(iproc == 0)then 
			print*,' Initial state wave function file'
		endif

		call wfn_ropen_file_select(wfnfile_i,'I') ! wfn_mpiio set here
		call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')	
		call wfn_read_nkeep(wfnfile_i,nkeep_i)
		if(iproc == 0)then
			print*,' There are ',nkeep_i,' initial wave functions '		
		endif
!........... OPEN FINAL .WFN FILE AND WRITE HEADER.....

       if( iproc ==0) then
			write(6,*)' Enter name of output .wfn file ' 
	   
	   		read(5,'(a)')outfile
	   endif
	  
	   call BMPI_BCAST(outfile,55,0,icomm,ierr)

       call wfn_wopen_file(wfnfile_f,.false.)
       call write_wfn_header(wfnfile_f)
	   call wfn_write_nkeep(nkeep_i) ! write number of vectors to wfn file
	  
	   
!........... NOW PROJECT..........................
       
	   if(iproc == 0) print*, "nkeep = ",nkeep_i
       do i = 1,nkeep_i
           ! new interface - we say which vec to read, it checks
           ! KSM:  This will be very slow, only need to read head part of each vector
		   !call wfn_readeigenvec_select('I',wfnfile_i,frag1,fcomm1,vec1,i,xe,xj,xt2)
		call wfn_readeigenvec_select('I',wfnfile_i,frag1,fcomm1,vec1,i,xe,xj,xt2)    
           if(iproc==0)print*,i,xe,xj,xt2	
		   call projectvec
	
		   !call BMPI_BARRIER(icomm,ierr)
           !call br_grab_vec2()
		   !print*,"iproc",iproc,"vec2",vec2
		   !call BMPI_BARRIER(icomm,ierr)
		   call wfn_writeeigenvec(wfnfile_f, frag2, vec2, i,xe,xj,xt2)
		   !call BMPI_BARRIER(icomm,ierr)
	   end do
	   !tidy up

	end subroutine proj_boss
	
end module proj_mod