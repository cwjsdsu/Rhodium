!
!
!

module entropy1body
	use precisions
	use dens1body
	implicit none

	
contains

!
!  SUBROUTINES CALLED
!
!  CALLED BY: main_menu
!	
subroutine entropy1b_boss
   	    use basis
	    use basis_map
		use wfn_mod
		use bvectorlib_mod
		use fragments
		use localvectors
		use io
		use mod_reorthog
		use spstate
		use sporbit
		use system_parameters
		use parr
		
		implicit none
		integer(4) :: nkeep_i,nkeep_f
	    integer ikeep
	    integer i,j,n, tidx
	    real :: xe,xj,xt2,ef,ei
		real :: xt
		real(8) ::dovlp
		integer :: ilast
		integer :: ierr
	    integer(4) :: ji,ti,jf,tf,jt,tt,jmin,jmax
	    integer(kind=basis_prec) :: k
	    logical :: evenAJ,evenAT    ! needed to force "correct" J, T
		logical :: zeroflag
	    real, allocatable :: denmat(:,:,:)
		integer :: numorbdim
		integer :: a,b
		integer :: start_i, stop_i,start_f,stop_f
		real(8) :: Stot,Sx
		integer :: mypar
		
		if(np_i(1)/=np_f(1) .or. np_i(2)/=np_f(2))then
			print*,' Not number conserving '
			print*,np_i, np_f
			stop
		end if
		

	    numorbdim=max(numorb_i(1),numorb_i(2) )
		numorbdim=max(numorbdim,numorb_f(1))
		numorbdim=max(numorbdim,numorb_f(2))

	       allocate( denmat( numorbdim, numorbdim, 0:1 ) )
	    if( mod(np_i(1)+np_i(2),2) == 1 )then


	         evenAJ = .false.
	       evenAT = .false.
	    else
	       evenAJ = .true.
	       evenAT = .true.
	    end if
		call parity_boss
		call phasechecker
		call basismapper(1)
		call basismapper(2)
		
		call set_nfragments(1)
	    frag1= nodal(iproc)%ifragment
	    frag2= nodal(iproc)%ffragment		
		call setup_localvectors_both
		
!............ OPEN INITIAL .WFN FILE AND READ HEADER.....

		print*,' Initial state wave function file'

		call wfn_ropen_file_select(wfnfile_i,'I')
		call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')	
		call wfn_read_nkeep(wfnfile_i,nkeep_i)
		print*,' There are ',nkeep_i,' initial wave functions '		
		print*,' Enter start, stop to compute density matrices '
		print*,' (Enter 0,0 to take all )'
		read(5,*)start_i,stop_i
		if(start_i < 1)start_i = 1
		if(stop_i < 1 .or. stop_i > nkeep_i)stop_i = nkeep_i
		if(start_i > stop_i)then
			print*,' start, stop INITIAL ', start_i,stop_i
			stop
		end if
		
        call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,1,ei,xj,xt2)  
		
		write(resultfile,*)'# Wavefunctions -- initial '
		do i = start_i,stop_i
            call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
            ji = closest2J(evenAJ,xj)
	        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			call assign_parity('i',mypar)
			
			call write_ejt(i, xe, ei,xj, xt,mypar)
			
		end do
		
!		rewind(wfnfile_i)
!		call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')	
		
!........... OPEN FINAL .WFN FILE AND WRITE HEADER.....

!       write(6,*)' Enter name of output .wfn file ' 
!	   read(5,'(a)')outfile	   
!       call wfn_wopen_file(wfnfile_f,.false.)
!       call write_wfn_header(wfnfile_f)
!	   call wfn_write_nkeep(nkeep_i) ! write number of vectors to wfn file

!        print*,' '
!        print*,' Final state wave function file'
!        print*,' Important: cannot be the same file as initial '

!        call wfn_ropen_file_select(wfnfile_f,'F')
!        call read_wfn_header_select(wfnfile_f,.false.,.true.,'F')	
		
!        call wfn_read_nkeep(wfnfile_f,nkeep_f)
!        print*,' There are ',nkeep_f,' final wave functions '


		
!		rewind(wfnfile_f)
!		call read_wfn_header_select(wfnfile_f,.false.,.true.,'F')	
		
		
		allocate(p1bopme(nrhsps_i(1),nrhsps_f(1)))
		allocate(n1bopme(nrhsps_i(2),nrhsps_f(2)))

!........... NOW PROJECT..........................

       !
       do i = start_i,stop_i
           ! new interface - we say which vec to read, it checks
           ! KSM:  This will be very slow, only need to read head part of each vector
           call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
           ji = closest2J(evenAJ,xj)
		   
		   ei =xe
!           if(iproc==0)print*,i,xe,xj,xt2	

!		   do j = start_i,stop_i
!		       call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,i,xe,xj,xt2)  
			   vec2= vec1
!	           jf = closest2J(evenAJ,xj)
!			   ef = xe
!	           jmax = (jf + ji)/2
!	           jmin = abs( jf -ji)/2
			   
!		       if(iproc==0)then
!				   write(6,'(i6,f10.4,f6.1,f7.2,f10.6)')j,xe,xj,xt2,dovlp
!			   end if
			   call p1bdensity
			   call n1bdensity

	            Stot = 0.d0
				call entropy1b(1,Sx)
				Stot= Stot+Sx
				call entropy1b(2,Sx)
				Stot= Stot+Sx
				print*,Stot,' Entropy ',ei,ji
				write(resultfile,*)ei,Stot,ji


!		   end do
		   

!		   call wfn_writeeigenvec(wfnfile_f, frag2, vec2, i,xe,xj,xt2) 

		
	   end do
!
if(iproc == 0)then
   write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'


   write(6,*)' '
   write(6,*)'All done with entropy '
   write(6,*)' '
end if
return
end subroutine entropy1b_boss	

!==================================================================================

subroutine entropy1b(it,s)
	use spstate
	use system_parameters
	implicit none
	integer :: it
	real(8) :: S
    real(kind=obs_prec),pointer   :: x1bopme(:,:)
	real(kind=8),allocatable :: mat(:,:),eval(:),work(:)
	integer maxdim
	integer :: i,j
	real :: sum
	
	integer :: info
	
	if(it==1)then
		x1bopme => p1bopme
	else
		x1bopme => n1bopme
	end if
	
	
	S = 0.d0
	
	maxdim = max(nrhsps_i(1),nrhsps_f(1))
	maxdim = max(nrhsps_i(2),maxdim)
	maxdim = max(nrhsps_f(2),maxdim)
	if(.not.allocated(mat))allocate(mat(maxdim,maxdim))
	if(.not.allocated(eval))allocate(eval(maxdim))
	if(.not.allocated(work))allocate(work(3*maxdim))
	
	mat = 0.d0
	do i = 1,nrhsps_i(it)
		do j = 1,nrhsps_i(it)
			mat(i,j)=x1bopme(i,j)
		end do
! 		mat(i,i)=mat(i,i) - real(np_i(it)-1,kind=8)/real(nrhsps_i(it),8)
	end do
	
	call DSYEV( 'N','U',nrhsps_i(it) , mat, maxdim, eval, WORK, 3*maxdim, INFO )
!	print*,eval(1:nrhsps_i(it))
    sum = 0.d0
	do i = 1,nrhsps_i(it)
!		if(eval(i)< 0.0)print*,' neg value ',eval(i)
!		if(eval(i)> 1.0)print*,' too much value ',eval(i)
!		print*,eval(i),1.d0-eval(i)
		eval(i)=eval(i)/real(np_i(it),kind=8)
		if(eval(i)> 0.d0 .and. eval(i) < 1.d0)S = S -eval(i)*log(eval(i)) - (1.d0 - eval(i))*log(1.d0 - eval(i))
!		print*,S,' so far ', eval(i), -eval(i)*log(eval(i))
		sum = sum + eval(i)
	end do
!	print*,'sum rule ',sum
	
	
	
	return
end subroutine entropy1b
!==================================================================================
	
end module entropy1body
	
	
	