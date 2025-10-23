

module dot_mod
	use basis
	use basis_map
	implicit none
	
	
contains

!
!  master routine for dot products
!
!  Store FINAL wave functions in br_history
!	
	subroutine dot_master
		use wfn_mod
		use bvectorlib_mod
		use fragments
		use localvectors
		
		implicit none
		
		integer(4) :: nkeep_i,nkeep_f
	    integer ikeep
	    integer i,j,n, tidx
	    real :: xe,xj,xt2
		real(8) ::dovlp
		
		call phasechecker
		call basismapper(1)
		call basismapper(2)
		
		call set_nfragments(1)
	    frag1= nodal(iproc)%ifragment
	    frag2= nodal(iproc)%ffragment
!	    call overlaptribution  ! DO I NEED TO ADD THIS?
		
		call setup_localvectors_both
		
!............ OPEN INITIAL .WFN FILE AND READ HEADER.....

		print*,' Initial state wave function file'

		call wfn_ropen_file_select(wfnfile_i,'I')
		call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')	
		call wfn_read_nkeep(wfnfile_i,nkeep_i)
		print*,' There are ',nkeep_i,' initial wave functions '


!............ OPEN FINAL .WFN FILE AND READ HEADER.....
			
		print*,' '
		print*,' Final state wave function file'
		print*,' Important: cannot be the same file as initial '

		call wfn_ropen_file_select(wfnfile_f,'F')
		call read_wfn_header_select(wfnfile_f,.false.,.true.,'F')	
				
		call wfn_read_nkeep(wfnfile_f,nkeep_f)
		print*,' There are ',nkeep_f,' final wave functions '
		
!.................. START READING IN WAVEFUNCTIONS AND TAKE DOT PRODUCTS...........
!
! look at subroutine overlap in bdenslib1.f90
! save vectors to history
!  
	    do i = 1,nkeep_i
	       ! new interface - we say which vec to read, it checks
	       ! KSM:  This will be very slow, only need to read head part of each vector
	       call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
	       if(iproc==0)then
			   write(6,*)' Initial Energy   J     T^2'
			   write(6,'(i6,f10.4,f6.1,f7.2)')i,xe,xj,xt2
		   end if
	       if(iproc==0)then
			   write(6,*)' Final  Energy    J     T^2    < f | i >'
		   end if
		   do j = 1,nkeep_f
		       call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,j,xe,xj,xt2)  
			   call doverlapvec_new(dovlp)
		       if(iproc==0)then
				   write(6,'(i6,f10.4,f6.1,f7.2,f10.6)')j,xe,xj,xt2,dovlp
			   end if
			   
			   
		   end do
		   
	    enddo
		
		return
		
	end subroutine dot_master
	
end module dot_mod
	