!
!
!

module dens1body
	use precisions
	use dot_mod
	implicit none
    logical :: pndensities    ! flag to denote printing out densities in proton-neutron formalism
	
    real(kind=obs_prec), allocatable,target   :: p1bopme(:,:), n1bopme(:,:)
    
	type densitypack
   	   integer :: Jmin,Jmax
   	   logical,allocatable :: zeroflag(:)
	   logical :: empty ! used when some matrix elements simply aren't used
   	   real(kind=obs_prec),allocatable :: denmat(:,:,:,:) ! denmat(J,a,b,T)
    end type densitypack
 
    logical :: writeunformatted1bden = .false.   ! if TRUE then write out as unformatted TO BE ADDED
  
    type (densitypack), allocatable,target :: densitybag(:,:)   ! densitybag(istate,fstate)
		
	integer :: start_i_default, stop_i_default ! used for passing start, stop 
	
	real, allocatable :: elistall(:)
	integer, allocatable :: j2listall(:)
	real, allocatable :: tlistall(:)
	integer, allocatable :: parlistall(:)
 
    logical :: usedensitybag = .false.     ! IF TRUE then use new density master routines
    										! automatically uses symmetry and writes out
  	
contains

!
!  SUBROUTINES CALLED
!
!  CALLED BY: main_menu
!	this deals with general formats
!
subroutine density1b_boss(missin,samelevels)
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
		logical :: missin ! looking for missing densities
		logical :: samelevels ! if TRUE, only read in one set of levels
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
		integer :: mypar
		integer :: flevel_shift ! used to help order final states
		
		if(np_i(1)/=np_f(1) .or. np_i(2)/=np_f(2))then
			print*,' Not number conserving '
			print*,np_i, np_f
			stop
		end if
		call parity_boss

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
		if(.not.missin .or. .not.samelevels)then
		   call phasechecker
		   call basismapper(1)
		   call basismapper(2)
		
		   call set_nfragments(1)
	       frag1= nodal(iproc)%ifragment
	       frag2= nodal(iproc)%ffragment		
		   call setup_localvectors_both
	    end if
		
!............ OPEN INITIAL .WFN FILE AND READ HEADER.....
        if(missin .and. samelevels)then
			start_i = start_i_default
			stop_i = stop_i_default
			call wfn_rewind(wfnfile_i)
 		    call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')				
 		    call wfn_read_nkeep(wfnfile_i,nkeep_i)						
		else

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
	     end if
		
        call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,1,ei,xj,xt2)  
		
		write(resultfile,"('# Wavefunctions -- initial ',2i4)")np_i(1),np_i(2)
		write(resultfile,"('  State       E          Ex           J        T     par ')")
		do i = start_i,stop_i
            call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
            ji = closest2J(evenAJ,xj)
			if(.not.missin .or. ji > 0)then
				call assign_parity('i',mypar)
			else
				mypar = 1   ! for missin=T and jf = 0, will skip anyway
			end if
			
	        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			
			call write_ejt(i, xe,ei, xj, xt,mypar)
			
		end do
        print*,' '
        print*,' Final state wave function file'
		   if(missin .and. samelevels)then
			   		 call wfn_ropen_file_TEMP(wfnfile_f,'F') ! A KLUGE; I HAVE TO CLOSE AND REOPEN, dont know why
		  else
              print*,' Important: cannot be the same file as initial '
             call wfn_ropen_file_select(wfnfile_f,'F')
	      endif
		
        call read_wfn_header_select(wfnfile_f,.false.,.true.,'F')	
		
        call wfn_read_nkeep(wfnfile_f,nkeep_f)
        if(missin .and. samelevels)then
			start_f = start_i
			stop_f = stop_i
		else
           print*,' There are ',nkeep_f,' final wave functions '
		   print*,' Enter start, stop to compute density matrices '
		   print*,' (Enter 0,0 to take all )'
		   read(5,*)start_f,stop_f
	    endif
		
		flevel_shift = 0
		if(missin .and. .not.samelevels)flevel_shift = stop_i 
		
		if(start_f < 1)start_f = 1
		if(stop_f < 1 .or. stop_f > nkeep_f)stop_f= nkeep_f
		
		if(start_f > stop_f)then
			print*,' start, stop FINAL ', start_f,stop_f
			stop
		end if
        call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,1,ei,xj,xt2)  
		
		if(.not.missin)then
		   write(resultfile,*)'# '
		   write(resultfile,"('# Wavefunctions -- final   ',2i4)")np_f(1),np_f(2)
	    endif
		do i = start_f,stop_f
            call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,i,xe,xj,xt2)  
			
            jf = closest2J(evenAJ,xj)
			if(.not.missin .or. .not.samelevels)then
				call assign_parity('f',mypar)
						call write_ejt(i+flevel_shift, xe,ei, xj, xt,mypar)
			end if
			
		end do
		
!		rewind(wfnfile_f)
!		call read_wfn_header_select(wfnfile_f,.false.,.true.,'F')	
		
        call pocc_write_orbits_alt(resultfile)
		if(allocated(p1bopme))deallocate(p1bopme)
		if(allocated(n1bopme))deallocate(n1bopme)
		
		allocate(p1bopme(nrhsps_f(1),nrhsps_i(1)))   !originally had f, i swapped incorrectly
		allocate(n1bopme(nrhsps_f(2),nrhsps_i(2)))

!........... NOW PROJECT..........................

       !
       do i = start_i,stop_i
           ! new interface - we say which vec to read, it checks
           ! KSM:  This will be very slow, only need to read head part of each vector
           call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
           ji = closest2J(evenAJ,xj)
		   if(missin .and. ji==0)cycle
		   
		   ei =xe
           if(isoflag)then
              ti = closest2J(evenAT,-0.5 + sqrt(xt2+0.25)) ! nint(-1 + sqrt(4*xtt+1))
           endif
		   do j = start_f,stop_f
		       call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,j,xe,xj,xt2)  
	           jf = closest2J(evenAJ,xj)
			   if(missin .and. jf==0)cycle
			   ef = xe
	           jmax = (jf + ji)/2
	           jmin = abs( jf -ji)/2			   
			   jmin = max (jmin,abs(Jz_i-Jz_f)/2)
	           if(isoflag)then
	              tf = closest2J(evenAT,-0.5 + sqrt(xt2+0.25)) ! nint(-1 + sqrt(4*xtt+1))
	           endif
!		       if(iproc==0)then
!				   write(6,'(i6,f10.4,f6.1,f7.2,f10.6)')j,xe,xj,xt2,dovlp
!			   end if
			   call p1bdensity
			   call n1bdensity
	           if ( iproc == 0 .and. isoflag) then
	                 write(resultfile,*)' '
	                 write(resultfile,333)i,ei,Ji,Ti 
	                 write(resultfile,334)j+flevel_shift,ef,Jf,Tf
	           end if
	           if ( iproc == 0 .and. .not.isoflag) then
	                 write(resultfile,*)' '
	                 write(resultfile,433)i,ei,Ji
	                 write(resultfile,434)j+flevel_shift,ef,Jf
	           end if
	   333     format(' Initial state #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
	   334     format(' Final state   #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
	   433     format(' Initial state #',i5,' E = ',f10.5,' 2xJ   = ',i4) 
	   434     format(' Final state   #',i5,' E = ',f10.5,' 2xJ   = ',i4) 

			   do Jt = jmin,jmax
				   
				   if(missin .and. (-1)**(Jt + (jf+ji)/2)==1)cycle
				   denmat = 0.0
   			        call density1b_coupled(Jt,Ji,Ti,Jf,Tf,Jz_i,Jz_f,np_i(1)-np_i(2),zeroflag,numorbdim,denmat)
	                if(zeroflag)cycle
                    if(pndensities)then
						write(resultfile,32)Jt
					else
						write(resultfile,322)Jt
					end if
   32               format(' Jt = ',i3,', proton      neutron ')
   322               format(' Jt = ',i3,', Tt = 0        1 ')
   
              do a = 1,numorbdim
                 do b = 1,numorbdim
                    if ( (abs(denmat(a,b,0)) + abs ( denmat(a,b,1))> 1.e-6) )then
                       write(resultfile,36)a,b,(denmat(a,b,tt),tt=0,1)
                    end if
36                  format(2i5, 2f12.7)
                 end do ! b
              end do   ! a
					
					
				end do
		   end do
		
	   end do
!
if(iproc == 0)then
   write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
   write(resultfile,*)' '
   write(resultfile,*)' Definition of density matrices : '
   write(resultfile,*)' rho_K(a^+b) =   (Jf || (a^+ b)_K || Ji) / sqrt(2K+1) '
   write(resultfile,*)'  where reduced matrix element is convention of Edmonds '
   write(resultfile,*)' (note: if isospin is good symmetry, then '
   write(resultfile,*)'  doubly reduced/ divided by sqrt(2T+1) as well'
   write(resultfile,*)' '
   write(resultfile,*)' Note time-reversal symmetry relation: '
   write(resultfile,*)' rho_K(a^+b, Ji->Jf) = (-1)^(ja-jb + Ji-Jf) rho_K(b^+a, Jf->Ji) '
   write(resultfile,*)' For isospin, add in factor (-1)^(Ti-Tf)'
   write(resultfile,*)' '
   write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
   write(resultfile,*)' '

   write(6,*)' '
   write(6,*)'All done with density matrices'
   write(6,*)' '
end if
return
end subroutine density1b_boss	

!==================================================================================
!
!  this uses density_bag and time-reversal (TR) to fill in missing densities
!  TASKS:
!   set up master list of energies and Js, Ts, especially if not samelevels
!   Note: another common application is when initial states are one parity and final are the other
!   for this, one can have either missin = TRUE or FALSE
!
!  INPUTS:
!    missin: if TRUE then only compute 'missing' (due to vanishing Clebsch-Gordan coeff) density matrix elements
!    samelevels: if TRUE, then the initial and final levels are the same; this often occurs when missin = TRUE
!
subroutine density1b_bag_boss(missin,samelevels,samebasis)
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
		logical :: missin ! looking for missing densities
		logical :: samelevels ! if TRUE, only read in one set of levels
		logical :: samebasis
		integer(4) :: nkeep_i,nkeep_f
	    integer ikeep
	    integer i,j,n, tidx
	    real :: xe,xj,xt2,ef,ei
		real :: xt
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
		integer :: mypar
		integer :: flevel_shift ! used to help order final states
		
		real, allocatable :: elist1(:),elist2(:)
		integer, allocatable :: j2list1(:),j2list2(:)
		real, allocatable :: tlist1(:),tlist2(:)
		integer, allocatable :: parlist1(:),parlist2(:)
		logical :: testsum,sameif
		integer :: phase,ja,jb
		real(kind=8) :: dovlp
		
		if(np_i(1)/=np_f(1) .or. np_i(2)/=np_f(2))then
			print*,' Not number conserving '
			print*,np_i, np_f
			stop
		end if
		call parity_boss

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
		if(.not.missin .or. .not.samelevels)then
		   call phasechecker
		   call basismapper(1)
		   call basismapper(2)
		
		   call set_nfragments(1)
	       frag1= nodal(iproc)%ifragment
	       frag2= nodal(iproc)%ffragment		
		   call setup_localvectors_both
	    end if
		
!............ OPEN INITIAL .WFN FILE AND READ HEADER.....
        if(missin .and. samelevels)then
			start_i = start_i_default ! this is set in apply1b_boss when writing to a temp file
			stop_i = stop_i_default   ! and specifically when 'raising' using J+
			call wfn_rewind(wfnfile_i)
 		    call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')				
 		    call wfn_read_nkeep(wfnfile_i,nkeep_i)						
		else   ! read in start, stop for initial states

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
	     end if
		 allocate(elist1(start_i:stop_i),j2list1(start_i:stop_i)) ! this is used for writing out
		 allocate(tlist1(start_i:stop_i),parlist1(start_i:stop_i))

        call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,1,ei,xj,xt2)  
		
!		write(resultfile,"('# Wavefunctions ',2i4)")np_i(1),np_i(2)
		write(resultfile,"('  State       E          Ex           J        T     par ')")
		do i = start_i,stop_i
            call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
			elist1(i)=xe ! make a list of energies, j, t, parity for thse initial states
            ji = closest2J(evenAJ,xj)
			j2list1(i)=ji
			
			if(.not.missin .or. ji > 0)then
				call assign_parity('i',mypar)
			else
				mypar = 1   ! for missin=T and jf = 0, will skip anyway
			end if
	        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			tlist1(i)=xt
			parlist1(i)=mypar
			call write_ejt(i-start_i+1, xe,ei, xj, xt,mypar)
			
		end do
        print*,' '
        print*,' Final state wave function file'
		   if(missin .and. samelevels)then
			   		 call wfn_ropen_file_TEMP(wfnfile_f,'F') ! A KLUGE; I HAVE TO CLOSE AND REOPEN, dont know why
		  else
              print*,' Important: cannot be the same file as initial '
             call wfn_ropen_file_select(wfnfile_f,'F')
	      endif
		
        call read_wfn_header_select(wfnfile_f,.false.,.true.,'F')	
		
        call wfn_read_nkeep(wfnfile_f,nkeep_f)
        if(missin .and. samelevels)then ! automatically set start, stop to be the same
			start_f = start_i
			stop_f = stop_i
		else
           print*,' There are ',nkeep_f,' final wave functions '
		   print*,' Enter start, stop to compute density matrices '
		   print*,' (Enter 0,0 to take all )'
		   read(5,*)start_f,stop_f
	    endif
		flevel_shift = 0
		if(missin .and. .not.samelevels)flevel_shift = stop_i 
		if(.not.missin .and. .not.samelevels)flevel_shift = stop_i +1 - start_i
		if(start_f < 1)start_f = 1
		if(stop_f < 1 .or. stop_f > nkeep_f)stop_f= nkeep_f
		
		if(start_f > stop_f)then
			print*,' start, stop FINAL ', start_f,stop_f
			stop
		end if
		if(.not.samelevels)then
   		 allocate(elist2(start_f:stop_f),j2list2(start_f:stop_f))
   		 allocate(tlist2(start_f:stop_f),parlist2(start_f:stop_f))
		end if
        call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,1,ei,xj,xt2)  
		
		do i = start_f,stop_f
            call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,i,xe,xj,xt2)  
            jf = closest2J(evenAJ,xj)
			if(.not.samelevels)then
				if(.not.missin .or. jf > 0)then
					call assign_parity('f',mypar)
				else
					mypar = 1   ! for missin=T and jf = 0, will skip anyway
				end if
						call write_ejt(i - start_f +1 +(stop_i-start_i+1), xe,ei, xj, xt,mypar)
			end if
			if(.not.samelevels)then
				elist2(i)=xe
				j2list2(i)=jf
				tlist2(i)=xt2
				parlist2(i)=mypar
			end if
			
		end do
		
!		rewind(wfnfile_f)
!		call read_wfn_header_select(wfnfile_f,.false.,.true.,'F')	
		
        call pocc_write_orbits_alt(resultfile)
		if(allocated(p1bopme))deallocate(p1bopme)
		if(allocated(n1bopme))deallocate(n1bopme)
		
		allocate(p1bopme(nrhsps_f(1),nrhsps_i(1)))   !originally had f, i swapped incorrectly
		allocate(n1bopme(nrhsps_f(2),nrhsps_i(2)))
!...... create masterlists of energies, j, t, parity....
! I assume that either the levels are the same, or no overlaps
! for more general case, do not call this subroutine!
!
if(samelevels)then
	allocate(elistall(stop_i+1-start_i))
	allocate(j2listall(stop_i+1-start_i))
	allocate(tlistall(stop_i+1-start_i))
	allocate(parlistall(stop_i+1-start_i))
	do i = start_i,stop_i
		elistall(i+1-start_i)=elist1(i)		
		j2listall(i+1-start_i)=j2list1(i)		
		tlistall(i+1-start_i)=tlist1(i)		
		parlistall(i+1-start_i)=parlist1(i)		

	end do
	allocate(densitybag(stop_i+1-start_i,stop_i+1-start_i))
	
else
	allocate(elistall(stop_i+1-start_i + stop_f +1 -start_f))
	allocate(j2listall(stop_i+1-start_i + stop_f +1 -start_f))
	allocate(tlistall(stop_i+1-start_i + stop_f +1 -start_f))
	allocate(parlistall(stop_i+1-start_i + stop_f +1 -start_f))
	allocate(densitybag(stop_i+1-start_i+ stop_f +1 -start_f,stop_i+1-start_i+ stop_f +1 -start_f))
	
	do i = start_i,stop_i
		elistall(i+1-start_i)=elist1(i)		
		j2listall(i+1-start_i)=j2list1(i)		
		tlistall(i+1-start_i)=tlist1(i)		
		parlistall(i+1-start_i)=parlist1(i)		

	end do
	do i = start_f,stop_f
		elistall(i+1-start_f+flevel_shift)=elist2(i)		
		j2listall(i+1-start_f+flevel_shift)=j2list2(i)		
		tlistall(i+1-start_f+flevel_shift)=tlist2(i)		
		parlistall(i+1-start_f+flevel_shift)=parlist2(i)		

	end do	
	
	
end if
		
densitybag(:,:)%empty=.true.  ! default
!........... NOW PROJECT..........................
       do i = start_i,stop_i
           ! new interface - we say which vec to read, it checks
           ! KSM:  This will be very slow, only need to read head part of each vector
           call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
           ji = closest2J(evenAJ,xj)
		   if(missin .and. ji==0)cycle
		   
		   ei =xe
           if(isoflag)then
              ti = closest2J(evenAT,-0.5 + sqrt(xt2+0.25)) ! nint(-1 + sqrt(4*xtt+1))
           endif
		   do j = start_f,stop_f
			   if(samelevels .and. j < i )cycle
			   if(samelevels .and. i==j)then
				   sameif = .true.
			   else
				   sameif = .false.
			   end if
		       call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,j,xe,xj,xt2)  
	           jf = closest2J(evenAJ,xj)
			   if(missin .and. jf==0)cycle
			   ef = xe
	           jmax = (jf + ji)/2
	           jmin = abs( jf -ji)/2
			   jmin = max (jmin,abs(Jz_i-Jz_f)/2)
			   
	           if(isoflag)then
	              tf = closest2J(evenAT,-0.5 + sqrt(xt2+0.25)) ! nint(-1 + sqrt(4*xtt+1))
	           endif
			   
			   if(samebasis .and. (phconj_i(1) .or.phconj_i(2) ))then
				   call doverlapvec_new(dovlp)
			   else
				   dovlp = 0.d0
			   end if


			   call p1bdensity
			   if(phconj_i(1))then  ! added to handle hole formalism
	              call phconjdensity(1,dovlp)
			   end if
			   call n1bdensity
			   if(phconj_i(2))then ! added to handle hold formalism
 	              call phconjdensity(2,dovlp)
	
			   end if
			   
               allocate(densitybag(i+1-start_i,j+1-start_f+flevel_shift)%zeroflag(jmin:jmax))
			   densitybag(i+1-start_i,j+1-start_f+flevel_shift)%zeroflag(:) = .false.
			   
               if(i+1-start_i /= j+1-start_f+flevel_shift)then
				   allocate(densitybag(j+1-start_f+flevel_shift,i+1-start_i)%zeroflag(jmin:jmax))
				   densitybag(j+1-start_f+flevel_shift,i+1-start_i)%zeroflag(:) = .false.
				   
			   end if
			   
               densitybag(i+1-start_i,j+1-start_f+flevel_shift)%jmin = jmin
               densitybag(i+1-start_i,j+1-start_f+flevel_shift)%jmax = jmax
               densitybag(j+1-start_f+flevel_shift,i+1-start_i)%jmin = jmin
               densitybag(j+1-start_f+flevel_shift,i+1-start_i)%jmax = jmax
               allocate(densitybag(i+1-start_i,j+1-start_f+flevel_shift)%denmat(jmin:jmax,numorbmax,numorbmax,0:1))
               densitybag(i+1-start_i,j+1-start_f+flevel_shift)%denmat= 0.0
               if(i+1-start_i /= j+1-start_f+flevel_shift)then
                  allocate(densitybag(j+1-start_f+flevel_shift,i+1-start_i)%denmat(jmin:jmax,numorbmax,numorbmax,0:1))
                  densitybag(j+1-start_f+flevel_shift,i+1-start_i)%denmat= 0.0
		       end if
			   
			   do Jt = jmin,jmax
				   if(missin .and. (-1)**(Jt + (jf+ji)/2)==1)then ! don't bother with cases we don't care about
					   densitybag(i+1-start_i,j+1-start_f+flevel_shift)%zeroflag(jt) = .true.
					   densitybag(j+1-start_f+flevel_shift,i+1-start_i)%zeroflag(jt) = .true.

					   cycle
				   end if
				   denmat = 0.0
   			        call density1b_coupled(Jt,Ji,Ti,Jf,Tf,Jz_i,Jz_f,np_i(1)-np_i(2),zeroflag,numorbdim,denmat)
	                if(zeroflag)cycle

              testsum = .true.
              do a = 1,numorbdim
				  ja = orbqn_i(1,a)%j
                 do b = 1,numorbdim
   				    jb = orbqn_f(1,b)%j
					 
					 if(abs(denmat(a,b,0))> 1.e-6 .or. abs(denmat(a,b,1))> 1.e-6)testsum = .false.
					 densitybag(i+1-start_i,j+1-start_f+flevel_shift)%denmat(jt,a,b,0)=denmat(a,b,0)
					 densitybag(i+1-start_i,j+1-start_f+flevel_shift)%denmat(jt,a,b,1)=denmat(a,b,1)
!					 print*,denmat
! NOW PUT IN TIME REVERSE
                     phase=(-1)**((Ji-Jf +ja-jb)/2)
                     if(isoflag.and..not.pndensities) phase=phase*(-1)**((Ti-Tf)/2)
					 if(i+1-start_i .ne. j+1-start_f+flevel_shift)then
						 densitybag(j+1-start_f+flevel_shift,i+1-start_i)%denmat(jt,b,a,0)=denmat(a,b,0)*phase
						 densitybag(j+1-start_f+flevel_shift,i+1-start_i)%denmat(jt,b,a,1)=denmat(a,b,1)*phase
						 
					 end if
                 end do ! b
              end do   ! a
			  densitybag(i+1-start_i,j+1-start_f+flevel_shift)%zeroflag(jt)=testsum
			  densitybag(j+1-start_f+flevel_shift,i+1-start_i)%zeroflag(jt)=testsum
			  if(.not.testsum)then !not empty
				  densitybag(i+1-start_i,j+1-start_f+flevel_shift)%empty = .false.
				  densitybag(j+1-start_f+flevel_shift,i+1-start_i)%empty = .false.
				  
			  end if
				end do ! jt
		   end do !j
		
	   end do ! i
	   
!
if(iproc == 0)then
	if(samelevels)then
	   call density1b_writeout_bag(stop_i - start_i +1)
   else
	   call density1b_writeout_bag(stop_i -start_i + 1+stop_f - start_f +1)
   end if
end if
return
end subroutine density1b_bag_boss	

!==================================================================================

subroutine density1b_coupled(Jt,Ji,Ti,Jf,Tf,Jzi,Jzf,Tz,zeroflag, ndim,denmat)
    use system_parameters
	
	use spstate
	use sporbit
	use precisions
	implicit none
	
    integer Jt, Tt      ! J, T of transition operator
    integer Ji,Ti       ! initial wfn J, T 
    integer Jf,Tf       ! final wfn J, T
    integer Jzi,Jzf          ! of wfns
    integer Tz          ! = (Z -N )/2
    integer it
    logical zeroflag     ! flag to indicate no matrix elements possible
    integer:: ndim        ! dimensions of coupled density matrices
    real denmat(ndim,ndim,0:1)
    real cleb
    real cg
    integer a,b,ja,ma,jb,mb
    integer ia,ib
    integer asps,ath,bsps,bth,asgn,bsgn
    integer iphase,tsign
    real cgt
    real, parameter :: cgtol = 5.0e-6
    logical altflag
    real(kind=obs_prec), pointer :: x1bopme(:,:)
    integer :: ierr

    denmat(:,:,:) = 0.0

    cg = cleb(Ji,Jzi,2*Jt,Jzf-Jzi,Jf,Jzf)
    cg = cg*sqrt(float(2*Jt+1)) /sqrt(float( (jf+1)))
    if(abs(cg) < cgtol)then
       zeroflag = .true.
       return
    else
      zeroflag = .false.
    endif
    do it = 1,2
       if(np_i(it) ==0 .or. np_f(it)==0)cycle
       if(it==1)then
          x1bopme => p1bopme
           tsign = 1
       else
          x1bopme => n1bopme
          tsign =-1
       end if
       do ia = 1, nrhsps_f(it)  ! important: note here that label a -> final state 
          a = rhspsqn_f(it,ia)%orb
          ja= rhspsqn_f(it,ia)%j
          ma= rhspsqn_f(it,ia)%m
          do ib = 1, nrhsps_i(it)  !...and label b --> initial state
 
             b = rhspsqn_i(it,ib)%orb
             jb= rhspsqn_i(it,ib)%j
             mb= rhspsqn_i(it,ib)%m
             if( ma /= mb+ (Jzf-Jzi)) cycle
             if(Jt > (ja + jb)/2 )cycle
             if(Jt < abs(ja - jb)/2 ) cycle
             iphase = (-1)**( (jb -mb)/2)
             if(isoflag .and. .not. pndensities)then
                 do tt = 0,1
                    cgt = cleb(Ti,Tz,2*Tt,0,Tf,Tz)*sqrt(float(2*Tt+1)/float(Tf+1))
                    if(abs(cgt) < cgtol)cycle
                denmat(a,b,tt) = denmat(a,b,tt) + cleb(ja,ma,jb,-mb,2*Jt,ma-mb)*iphase & 
              * tsign* x1bopme(ia,ib)*cleb(1,tsign,1,-tsign,2*Tt,0) /(cg*cgt)
!			  if(x1bopme(ia,ib)/=0.0)print*,'(T)',cleb(ja,ma,jb,-mb,2*Jt,ma-mb),x1bopme(ia,ib),ia,ib
			  
                 end do
             else				 
                 denmat(a,b,it-1) = denmat(a,b,it-1)+ cleb(ja,ma,jb,-mb,2*Jt,ma-mb)*iphase & 
              *  x1bopme(ia,ib)/cg
!			  print*,it,x1bopme(ia,ib),ia,ib
!			  if(it==2)print*,a,b,x1bopme(ia,ib),cleb(ja,ma,jb,-mb,2*Jt,ma-mb)
             endif

         end do  ! bsps
      enddo  ! asps
    end do
	
	return
end subroutine density1b_coupled

! write out orbit quantum numbers with label so that
! density matrix report is understandable
! written in alternate form
subroutine pocc_write_orbits_alt(fn)
   use sporbit
   implicit none
   integer :: fn ! file to write to
   integer :: i
   
   write(fn, *) " "  ! blank line
   write(fn,*) ' Single particle state quantum numbers'
   write(fn,"(A)", advance='yes') 'ORBIT      N     L   2 x J  '
   do i=1,numorb_f(1)
      write(fn, 101, advance='yes') i,orbqn_f(1,i)%nr,orbqn_f(1,i)%l, orbqn_f(1,i)%j
   end do

   write(fn, *) " "  ! newline at end
   write(fn, *) " "  ! blank line
101 format(4i6)
end subroutine pocc_write_orbits_alt

!==================================================================================
!
! base routine to compute proton 1-body density matrices
! one issue to watch is truncation
!
subroutine p1bdensity
	use jumps1body
	use sectors
	use basis_map
	use localvectors
	use w_info
	use spstate ! for testing only
	use sporbit,only:sameweightscheme
	implicit none
	
	
	integer :: ipsj    ! label of proton sector jump
	integer :: is,fs,ics,fcs      ! label of initial/final proton sector, conjugate neutron sector(s)
	integer :: ic
	integer :: ncsi,ncsf     ! # of initial, final neutron sectors
	integer :: ipjmp   ! label of proton jumps
	
	integer(kind=basis_prec) :: isdp, fsdp   ! initial, final proton sds
	integer(kind=basis_prec) :: isdn, fsdn   ! initial, final neutron sds
	                                           ! in principle insdi=insdf, but may be different
											   ! if bases are different
	integer :: pphase,nphase
	integer :: iopp       ! label of one-body proton operator
    integer(kind=basis_prec) :: ipstate,fpstate,instate,fnstate
	real(8) :: xme
	integer :: iop1body
	integer :: ac,bd    ! proton creation/annihilation operators
	integer(4) :: Wpf,Wnf
	
	p1bopme=0.d0
!
!$OMP PARALLEL DO PRIVATE(is,fs,ncsi,iop1body,wnf,ipjmp,isdn,fsdn,instate,fnstate,nphase,iopp,bd,ac) &
!$OMP PRIVATE(ic,ics,wpf,isdp,fsdp,pphase,ipstate,fpstate,xme), REDUCTION(+:p1bopme)	
	do ipsj = 1,nsectorjumps(1)  ! loop over proton sector jumps
		is = p1bsectjump(ipsj)%isector
		fs = p1bsectjump(ipsj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
										! NOTE ADDED MAY 2024: but if cross-parity, cannot have is = fs
		ncsi = 	xsd_i(1)%sector(is)%ncsectors				
		iop1body = 			p1bsectjump(ipsj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',ipsj,iop1body
			cycle
		end if
        wpf = p1bsectjump(ipsj)%wXf
!						print*,' map ',nsdmap				
		do ipjmp = 1,p1bsectjump(ipsj)%njumps
			isdp = p1bsectjump(ipsj)%isd(ipjmp)
			fsdp = p1bsectjump(ipsj)%fsd(ipjmp)
			if(fsdp==-1)cycle
			ipstate = pstart_i(isdp)
			fpstate = pstart_f(fsdp)
			pphase= p1bsectjump(ipsj)%phase(ipjmp)  ! fetch phase from applying the operator
			iopp   = p1bsectjump(ipsj)%iXop(ipjmp)  ! find index of one-body operator
			bd = pops1body(iop1body)%Xop(iopp,1)
			ac = pops1body(iop1body)%Xop(iopp,2)
!			print*,isdp,fsdp,ac,bd
!			print*,'(A)',ipjmp,isdp,fsdp,ac,bd,ncsi
!............ LOOP OVER CONJUGATE NEUTRON SDS....................			
			
			do ic = 1,ncsi   ! loop over conjugate neutron sectors
				ics = xsd_i(1)%sector(is)%csector(ic)  ! find index of conjugate neutron sector
				                                       ! HOWEVER must be careful when truncating--
				                                       ! not all FINAL neutron SD will be allowed
				                                       ! if the sum of final Ws goes too high.
													  
				if(sameweightscheme)then              ! can check this issue directly
					wnf = xsd_i(2)%sector(ics)%wx
					
					if(wpf + wnf > maxWtot_f)cycle									   					
				end if
				
				do isdn = xsd_i(2)%sector(ics)%xsdstart,xsd_i(2)%sector(ics)%xsdend   !loop over neutron SDs 
					fsdn= nsdmap(isdn)  ! need to map this initial, unchanged neutron SD
					                      ! to a final neutron SD in the final basis
										  
					if(.not.sameweightscheme)then
							wnf=nsdmap_w(isdn)		
							if(wpf + wnf > maxWtot_f)cycle									   					
									  
					end if
					if(fsdn ==-1)cycle
					nphase = nsdmapphase(isdn)   ! also may get a phase
					instate = nstart_i(isdn)
					fnstate = nstart_f(fsdn)
!					print*,ipstate,instate,fpstate,fnstate
!print*,'(B)',isdn,fsdn,instate,ipstate,fnstate,fpstate

					xme = vec1(ipstate+instate)*vec2(fpstate+fnstate)*pphase*nphase
!					print*,'(b)',ipstate+instate,fpstate+fnstate,xme
!print*,ac,bd,ipstate+instate,fpstate+fnstate,isdp,fsdp,isdn,fsdn
					p1bopme(ac,bd)=p1bopme(ac,bd)+xme   ! note here: ac acts on final state, bd acts on initial state
					
				end do
				
			end do
			
		end do ! ipjmp
		
		
		
	end do ! ipsj

	
	return
end subroutine p1bdensity
!==================================================================================
!
! base routine to compute neutron 1-body density matrices
! one issue to watch is truncation
!
subroutine n1bdensity
	use jumps1body
	use sectors
	use basis_map
	use localvectors
	use spstate ! for testing only
	use sporbit,only:sameweightscheme
	use W_info,only:maxWtot_f
	implicit none
	
	
	integer :: insj    ! label of proton sector jump
	integer :: is,fs,ics,fcs      ! label of initial/final proton sector, conjugate neutron sector(s)
	integer :: ic
	integer :: ncsi,ncsf     ! # of initial, final neutron sectors
	integer :: injmp   ! label of proton jumps
	
	integer(kind=basis_prec) :: isdp, fsdp   ! initial, final proton sds
	integer(kind=basis_prec) :: isdn, fsdn   ! initial, final neutron sds
	                                           ! in principle insdi=insdf, but may be different
											   ! if bases are different
	integer :: pphase,nphase
	integer :: iopn      ! label of one-body proton operator
    integer(kind=basis_prec) :: ipstate,fpstate,instate,fnstate
	real(8) :: xme
	integer :: iop1body
	integer :: ac,bd    ! proton creation/annihilation operators
	integer(4) :: Wpf,Wnf
	
	
	n1bopme=0.d0
!$OMP PARALLEL DO PRIVATE(is,fs,ncsi,iop1body,wnf,injmp,isdn,fsdn,instate,fnstate,nphase,iopn,bd,ac) &
!$OMP PRIVATE(ic,ics,wpf,isdp,fsdp,pphase,ipstate,fpstate,xme), REDUCTION(+:n1bopme)
	do insj = 1,nsectorjumps(2)  ! loop over proton sector jumps
		is = n1bsectjump(insj)%isector
		fs = n1bsectjump(insj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
!										print*,'sectors ',is,fs,p1bsectjump(ipsj)%dJz,p1bsectjump(ipsj)%dW
		ncsi = 	xsd_i(2)%sector(is)%ncsectors				
		iop1body = 			n1bsectjump(insj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',insj,iop1body
			cycle
!		else
!			print*,' non empty sector ',ipsj,iop1body
		end if
!		ncsf = 	xsd_f(1)%sector(fs)%ncsectors							
        wnf = n1bsectjump(insj)%wXf
										
		do injmp = 1,n1bsectjump(insj)%njumps
			isdn = n1bsectjump(insj)%isd(injmp)
			fsdn = n1bsectjump(insj)%fsd(injmp)
			if(fsdn==-1)cycle
			instate = nstart_i(isdn)
			fnstate = nstart_f(fsdn)
			nphase= n1bsectjump(insj)%phase(injmp)  ! fetch phase from applying the operator
			iopn   = n1bsectjump(insj)%iXop(injmp)  ! find index of one-body operator
			bd = nops1body(iop1body)%Xop(iopn,1)
			ac = nops1body(iop1body)%Xop(iopn,2)
			
!............ LOOP OVER CONJUGATE NEUTRON SDS....................			
			
			do ic = 1,ncsi   ! loop over conjugate neutron sectors
				ics = xsd_i(2)%sector(is)%csector(ic)  ! find index of conjugate proton sector
				                                       ! HOWEVER must be careful when truncating--
				                                       ! not all FINAL neutron SD will be allowed
				                                       ! if the sum of final Ws goes too high.
				if(sameweightscheme)then              ! can check this issue directly
					wpf = xsd_i(1)%sector(ics)%wx
					
					if(wpf + wnf > maxWtot_f)cycle	
					
				else
					
					
				end if									   
				do isdp = xsd_i(1)%sector(ics)%xsdstart,xsd_i(1)%sector(ics)%xsdend   !loop over proton SDs 
					fsdp= psdmap(isdp)  ! need to map this initial, unchanged proton SD
					                      ! to a final proton SD in the final basis
					if(fsdp ==-1)cycle
					
					if(.not.sameweightscheme)then
							wpf=psdmap_w(isdp)		
							if(wpf + wnf > maxWtot_f)cycle									   					
									  
					end if
					pphase = psdmapphase(isdp)   ! also may get a phase
					ipstate = pstart_i(isdp)
					fpstate = pstart_f(fsdp)
!					print*,ipstate,instate,fpstate,fnstate
					xme = vec1(ipstate+instate)*vec2(fpstate+fnstate)*pphase*nphase
!					print*,ipstate+instate,fpstate+fnstate,xme
					n1bopme(ac,bd)=n1bopme(ac,bd)+xme   ! note here: ac acts on final state, bd acts on initial state
					
				end do
				
			end do
			
		end do ! ipjmp
		
		
		
	end do ! ipsj
	return
end subroutine n1bdensity
!=========================================================================
!
!  write out one-body densities at the end
!  NEW VERSION added to RHODIUM 0.9.0
!     
! intended when having a consolidated list
!
! CALLED BY:
!
!

subroutine density1b_writeout_bag(nlevels)
  use system_parameters
  use sporbit
  use spstate
  use localvectors
  use nodeinfo
  use io
  use basis
  use obs
!  use lanczos_info
  use localvectors
!  use densities
!  use haiku_info
  use precisions
!  use fragments
!  use mod_reorthog
  use wfn_mod
  use butil_mod
  use bvectorlib_mod
!  use para_main_mod
!  use jump_mod
!  use lanczos_util
!  use pocc_mod
  
  implicit none
  integer(4),intent(in) :: nlevels !istart,istop,fstart,fstop ! start, stop for initial, final states

  real(4) :: xj,xt,ei,ef,xjj !,xtt,xjj
  integer :: pari,parf
  integer(4) :: i,j,m,n
  integer(4) :: ji,ti,jf,tf,jt,tt,jmin,jmax
  logical :: evenAJ,evenAT    ! needed to force "correct" J, T

  integer :: ierr
!  integer :: aerr
!--------------------- CONVERGENCE -----------------------------

  logical smallflag,zeroflag
  
  real(8), allocatable :: denmat(:,:,:)
  real :: dentol = 1.0e-6
  real :: nparticles
  logical numberflag  ! flag to check on total number of particles; used for debugging
  integer inode
 ! character(1) :: vchar
  real(4) :: nxprot,nxneut
  real(4) :: xtz

  if(iproc/=0)return
    
!  call clocker('obw','sta')

!........ DETERMINE IF EVEN OR ODD SYSTEM............

  if( mod(np_i(1)+np_i(2),2) == 1 )then

! if "spinless" then J is always integer (because it is really L)
! however T can be half-integer (if "spinless" then this is really S)
     if(spinless)then
       evenAJ = .true.
     else
       evenAJ = .false.
     end if
     evenAT = .false.
  else
     evenAJ = .true.
     evenAT = .true.
  end if
  do i = 1, nlevels !istop-istart+1 !istart,istop
        ei = elistall(i)
!        xj = 0.5*j2listall(i)  
        ji = j2listall(i)  

!        xtt= 0.5*(tlistall(i)   
        if(isoflag)ti= tlistall(i)   
		pari = parlistall(i)
!  	   if(only_odd .and. (ji/4)*4==ji)cycle		

     do j = 1,nlevels 
		 if(densitybag(i,j)%empty)cycle

           ef = elistall(j)

           jf = j2listall(j)
		   
		    !         
           if(isoflag)then
              tf = closest2J(evenAT,tlistall(j)) ! 
           endif
!	  	   if(only_odd .and. (jf/4)*4==jf)cycle
		   
        jmax = (jf + ji)/2
        jmin = abs( jf -ji)/2
	   jmin = max (jmin,abs(Jz_i-Jz_f)/2)
		
		parf = parlistall(j)
        if (  isoflag) then
              write(resultfile,*)' '
              write(resultfile,333)i,ei,Ji,Ti 
              write(resultfile,334)j,ef,Jf,Tf
        end if
        if ( .not.isoflag) then
              write(resultfile,*)' '
              write(resultfile,433)i,ei,Ji
              write(resultfile,434)j,ef,Jf 
        end if
333     format(' Initial state #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
334     format(' Final state   #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
335     format(' State #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
336     format(' State   #',i5,' E = ',f10.5,' 2xJ = ',i4) 

433     format(' Initial state #',i5,' E = ',f10.5,' 2xJ   = ',i4) 
434     format(' Final state   #',i5,' E = ',f10.5,' 2xJ   = ',i4) 


!---------------------------------------------------------------

!---------------------- WRITE OUT DENSITY MATRIX ELEMENTS
             do jt = jmin,jmax
                 if(densitybag(i,j)%zeroflag(jt))cycle

              if(isoflag .and. .not.pndensities)then
                 write(resultfile,31)Jt
31               format(' Jt = ',i3,', Tt = 0        1 ' )
              else
                 write(resultfile,32)Jt
32               format(' Jt = ',i3,', proton      neutron ')

              endif
              do m = 1,numorbmax
                 do n = 1,numorbmax
					 
					if(orbqn_i(1,m)%j + orbqn_i(1,n)%j < 2*Jt .or. abs(orbqn_i(1,m)%j - orbqn_i(1,n)%j) > 2*Jt)cycle
                    if ( abs(densitybag(i,j)%denmat(jt,m,n,0)) > dentol .or.  & 
					     abs(densitybag(i,j)%denmat(jt,m,n,1)) > dentol)then		
                       write(resultfile,36)m,n,(densitybag(i,j)%denmat(jt,m,n,tt),tt=0,1)
                    end if
36                  format(2i5, 2f12.7)
                 end do ! b
				 if(i==j .and. jt ==0)then
                     xjj = 0.5* orbqn_i(1,m)%j
					 
					 if(isoflag .and. .not. pndensities)then
						 xtz = 0.5*(np_i(1)-np_i(2))
						 nxprot = densitybag(i,j)%denmat(jt,m,m,0) 
						 if(ti >0 )nxprot = nxprot + sqrt(3.)*xtz/sqrt(0.5*Ti*(0.5*Ti+1.0))* densitybag(i,j)%denmat(jt,m,m,1)
						 nxneut = densitybag(i,j)%denmat(jt,m,m,0)
						 if(ti> 0)nxneut=nxneut - sqrt(3.)*xtz/sqrt(0.5*Ti*(0.5*Ti+1.0))* densitybag(i,j)%denmat(jt,m,m,1)
						 
						 nxprot = nxprot* sqrt(2.)*sqrt(2*xjj+1)*0.5/sqrt((Ji+1.0)*(Ti+1.0))
						 nxneut = nxneut* sqrt(2.)*sqrt(2*xjj+1)*0.5/sqrt((Ji+1.0)*(Ti+1.0))
!						 write(occresultfile,701)m,nxprot,nxneut
					 else
					 
!					 write(occresultfile,701)m,(densitybag(i,j)%denmat(jt,m,m,tt)*sqrt(2*xjj+1.0)/sqrt(jf+1.0),tt=0,1)
				     end if
				 end if
701 format(i4,2f10.4)
				 
				 
              end do   ! a
           end do ! jt

        end do   !j

  end do  ! i
     
  if(iproc == 0)then
     write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
     write(resultfile,*)' '
     write(resultfile,*)' Definition of density matrices : '
     write(resultfile,*)' rho_K(a^+b) =   (Jf || (a^+ b)_K || Ji) / sqrt(2K+1) '
     write(resultfile,*)'  where reduced matrix element is convention of Edmonds '
     write(resultfile,*)' (note: if isospin is good symmetry, then '
     write(resultfile,*)'  doubly reduced/ divided by sqrt(2T+1) as well'
     write(resultfile,*)' '
     write(resultfile,*)' Note time-reversal symmetry relation: '
     write(resultfile,*)' rho_K(a^+b, Ji->Jf) = (-1)^(ja-jb + Ji-Jf) rho_K(b^+a, Jf->Ji) '
     write(resultfile,*)' For isospin, add in factor (-1)^(Ti-Tf)'
     write(resultfile,*)' '
     write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
     write(resultfile,*)' '
!	 write(occresultfile,*)' '

  end if
!  call clocker('obw','end')

end subroutine density1b_writeout_bag

!==================================================================================

subroutine write_ejt(i, xe, e0, xj, xt,xpar)
   use nodeinfo
   use io
   implicit none
   integer :: i
   real :: xe
   real    :: xj, xt,e0
   integer :: xpar
   if(iproc == 0) then
      write(6, 100) i, xe, xe-e0, xj, xt,xpar
      if(writeout) write(resultfile, 100) i, xe,xe-e0, xj, xt,xpar
   end if
100 format(i5,3x,2f12.5,2x,2f8.3,2x,i2)
end subroutine write_ejt
!================================================
!
!  function to force conversion of unconverged xJ to integer J
!  that is, odd 2 x J for odd A, and even for even A
!
  function closest2J(evenA,xj)

  implicit none
  integer closest2J
  real xj
  logical evenA

  if(evenA)then
     closest2J = 2*nint(xj)
     if(closest2J < 0)closest2J = 0
  else
     closest2J = 2*nint(xj-0.5)+1
     if(closest2J < 1)closest2J = 1
  end if

  return
  end function closest2J
!
!======================================================
! ORIGINAL CODE TAKEN FROM BIGSTICK file bdenslib1.f90

!
!  swaps indices (+phase) if there is particle-hole conjugation
!  (Originally added 7.2.5 by CWJ @ SDSU )
!
!  NOTES added in version 7.5.6 by CWJ 9/2015
!  a hole creation operator of angular momentum j,m is equivalent to 
!  a particle annihilation operator j,-m with phase (-1)^(j+m)
!
!  The arrays x1bopme(asps,bsps) is empty unless m_asps = m_bsps
!  To properly map this, we must find the conjugates for asps and bsps
!  (and the associated phase)
!  Note: the time-reverse state is already stored in hspsqn(it,i)%tr
!
!  NOTE: MIGHT NOT WORK FOR "SPINLESS" SYSTEMS
!
subroutine phconjdensity(it,dovlp)
   use system_parameters
   use spstate
   use sporbit
!   use densities
   use precisions
!   use haiku_info
   implicit none
!........ INPUT...............
   integer :: it   ! species 
   real(kind=8) :: dovlp
!....... INTERNAL.......
   real(kind=obs_prec), pointer :: x1bopme(:,:)
   integer a,b,ja,ma,jb,mb,la,lb
   integer ia,ib
!   integer asps,ath,bsps,bth,asgn,bsgn
   integer :: iatr, ibtr       ! time-reversed states
   real(kind=obs_prec) :: diag
   real(kind=obs_prec) :: tmpden
   integer :: tsign
   integer :: iphase
   logical, allocatable :: checkoff(:,:)

   if(.not. phconj_i(it) .or. npeff_i(it) < 1 ) return
   
   
   if(it==1)then
      x1bopme=> p1bopme
      tsign = 1
   else
      x1bopme=> n1bopme
      tsign = -1 
   end if
   allocate(checkoff(nrhsps_f(it),nrhsps_i(it)))
   
   checkoff =.false.
       do ia = 1, nrhsps_f(it)  ! important: note here that label a -> final state 
          a = rhspsqn_f(it,ia)%orb
		  if(orbqn_f(it,a)%w==99)cycle
          ja= rhspsqn_f(it,ia)%j
          la= rhspsqn_f(it,ia)%l

          ma= rhspsqn_f(it,ia)%m
		  iatr = rhspsqn_f(it,ia)%tr
		  
          do ib = 1, nrhsps_i(it)  !...and label b --> initial state
			  
 
             b = rhspsqn_i(it,ib)%orb
   		     if(orbqn_i(it,b)%w==99)cycle
			 
             jb= rhspsqn_i(it,ib)%j
             lb= rhspsqn_i(it,ib)%l
			 
             mb= rhspsqn_i(it,ib)%m
             if( ma /= mb+ (Jz_f-Jz_i)) cycle

             iphase = (-1)**( (ja +ma +jb+mb)/2)
   		     ibtr = rhspsqn_i(it,ib)%tr
			 
			 if(checkoff(ia,ib) .or. checkoff(ibtr,iatr))cycle
			 checkoff(ia,ib) =.true.  ! to prevent redoing
			 checkoff(ibtr,iatr)=.true.
			 if(ia==ib .and. ma==mb .and. a == b .and. ja==jb .and. la==lb)then
				 diag=dovlp
			 else
				 diag = 0.0
			 end if
			 tmpden = x1bopme(ia,ib)
			 x1bopme(ia,ib) = diag-x1bopme(ibtr,iatr)*iphase
			 x1bopme(ibtr,iatr) = diag- tmpden*iphase

         end do  ! bsps
      enddo  ! asps
	  deallocate(checkoff)
   return
end subroutine phconjdensity
!======================================================
  
  
  
  
	
end module dens1body
	
	
	