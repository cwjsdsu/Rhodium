
!
! routines for projecting wave functions into a different basis
!

module lincom_mod
	
	use basis
	use basis_map
	implicit none

	logical :: norm_lincom = .true. ! normalize each eigenstate v_i*|i> 


contains
! BOSS ROUTINE FOR CONSTRUCTING LINEAR COMBINATION OF WAVEFUNCTIONS


	subroutine lincom_boss
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
		real(8) ::dovlp, wfn_amp
		real(kind=8) :: dnorm
		integer :: ilast
		integer :: ierr, pi
		integer :: nsource, which_wfn
		logical :: smallflag



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

					!open output file and bcast filename
		!........... OPEN FINAL .WFN FILE AND WRITE HEADER.....

		if( iproc ==0) then
			write(6,*)' Enter name of output .wfn file ' 
	   
	   		read(5,'(a)')outfile
	   endif
	  
	   call BMPI_BCAST(outfile,55,0,icomm,ierr)
       call wfn_wopen_file(wfnfile_f,.false.)
       call write_wfn_header(wfnfile_f)

		if(iproc==0) then
			write(6,*)' Enter the number of quantity of source files'
			read(5,*) nsource
		endif


		call BMPI_BARRIER(icomm,ierr)
	  	call BMPI_BCAST(nsource,1,0,icomm,ierr)
		if(iproc ==0) print*, "number of vectors to be written to outfile: ", nsource
		call wfn_write_nkeep(nsource) ! write number of vectors to wfn file
		
		if(iproc == 0)then 
			print*,'starting lincom main loop...'
		endif
		!now begin main event loop
			do i = 1, nsource
			 call wfn_ropen_file_select(wfnfile_i,'I') ! wfn_mpiio set here
			 call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')	
			 call wfn_read_nkeep(wfnfile_i,nkeep_i)
			 if(iproc == 0)then
			 	 print*,' There are ',nkeep_i,' initial wave functions '		
			 endif

			 if(iproc==0) then
				write(6,*)' Enter which wave function (state number) to be used in linear combination'
				read(5,*) which_wfn
			endif
			!if(iproc ==0) print*, "state number selected from input wfn: ", which_wfn
			call BMPI_BARRIER(icomm,ierr)
			call BMPI_BCAST(which_wfn,1,0,icomm,ierr)

			if(iproc==0) then
				write(6,*)'Enter amplitude'
				read(5,*) wfn_amp
			endif
			call BMPI_BARRIER(icomm,ierr)
			call BMPI_BCAST(wfn_amp,1,0,icomm,ierr)

			 !if(iproc ==0) print*, "which_wfn ", which_wfn
			 if(iproc ==0) print*, " input wfn_amp: ", wfn_amp
			 call wfn_readeigenvec_select('I',wfnfile_i,frag1,fcomm1,vec1,which_wfn,xe,xj,xt2)
				if(iproc==0)print*,"state(iter,ex,J,T):",i,xe,xj,xt2	
				!call lincomvec
			!multiply eigenvector by amplitude
			call dvecdot_s(wfn_amp)
			!normalize result after amplitude is applied.
			if(norm_lincom) then
				if(iproc ==0) print*, "normalizing eigenstate..."
				 call dnormvec_p('n','f',dnorm,smallflag)
			endif 

			call wfn_writeeigenvec(wfnfile_f, frag2, vec2,i,xe,xj,xt2)
			close(wfnfile_i)
			!if(iproc ==0) print*,"finished writing writing wfn to file -iter:", i
		 end do

	   
	   if(iproc ==0) print*,"finished lincom_boss"

	end subroutine lincom_boss

	!===================================================
	subroutine dvecdot_s(sclr)
		!
		! double-precision scalr*vec
		! dvec2 = dvec1*scalar
		!
		!
			  use precisions
			  use nodeinfo
			  use localvectors
			  use fragments
			  implicit none
			  !character(1) :: vchar

			  real(kind=lanc_prec), pointer :: dvec1(:), dvec2(:)
			  real(kind=8) :: sclr
		
		!------------------ INTERMEDIATE -------------
			  real(kind=8) :: d1, d2
			  integer      :: iter, i
		
			  integer(kind=basis_prec) :: il,vstart,vstop

			! This code didn't work before when frag1 /= frag2
			vstart = v1s  !    basestart(frag1)
			vstop  = v1e  !    basestop (frag1)

			dvec1 => vec1
			dvec2 => vec2

		!$OMP PARALLEL SHARED(dvec1,dvec2,vstart,vstop,sclr), PRIVATE(il,d1)
		!$OMP DO SCHEDULE(STATIC)
		   do il = vstart,vstop
				 d1 = dvec1(il)
				 !d2 = dvec2(il)
				 dvec2(il) = d1*sclr
		   enddo
		!$OMP END DO
		!$OMP END PARALLEL
		
		   return
		end subroutine dvecdot_s
		
	
end module lincom_mod