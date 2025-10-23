!
! charge-changing densities
!

module dens1bodycc
	use precisions
	use dens1body
	use hops
	implicit none
	
    real(kind=obs_prec), allocatable,target   :: cc1bopme(:,:)
    
	
contains

!
!  SUBROUTINES CALLED
!
!  CALLED BY: main_menu
!	
subroutine ccdensity1b_boss
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
		integer(4) :: jn,jp
	    integer(kind=basis_prec) :: k
	    logical :: evenAJ,evenAT    ! needed to force "correct" J, T
		logical :: zeroflag
	    real, allocatable :: denmat(:,:,:)
		integer :: numorbdim
		integer :: a,b
		integer :: start_i, stop_i,start_f,stop_f
		logical :: p2n
		
		integer :: maxjt
		integer :: mypar
		
		
		if(np_f(1)-np_i(1)==1)then
			p2n=.false.
		
			if(np_f(2)-np_i(2)/=-1)then  ! error trap
				print*,' something wrong in ccdensity1b_coupled (1)'
				stop
			end if
		else
			p2n=.true.
		
			if(np_f(2)-np_i(2)/=1)then  ! error trap
				print*,' something wrong in ccdensity1b_coupled (2)'
				stop
			end if
			if(np_f(1)-np_i(1)/=-1)then  ! error trap
				print*,' something wrong in ccdensity1b_coupled (3)'
				stop
			end if
		
		end if
		call parity_boss

	    numorbdim=max(numorb_i(1),numorb_i(2) )
		numorbdim=max(numorbdim,numorb_f(1))
		numorbdim=max(numorbdim,numorb_f(2))

	       allocate( denmat( numorbdim, numorbdim, 1 ) )
	    if( mod(np_i(1)+np_i(2),2) == 1 )then
	         evenAJ = .false.
	       evenAT = .false.
	    else
	       evenAJ = .true.
	       evenAT = .true.
	    end if
		
		call phasechecker
!		call basismapper(1)
!		call basismapper(2)
		
		call set_nfragments(1)
	    frag1= nodal(iproc)%ifragment
	    frag2= nodal(iproc)%ffragment		
		call setup_localvectors_both
		
!....... find maxjt..........

        maxjt = 0
			jp = 0
			jn = 0
			do i = 1,numorb_i(1)
				jp = max(jp,orbqn_i(1,i)%j)				
			end do
			do i = 1,numorb_f(1)
				jp = max(jp,orbqn_f(1,i)%j)				
			end do
			do i = 1,numorb_i(2)
				jn = max(jn,orbqn_i(2,i)%j)				
			end do
			do i = 1,numorb_f(1)
				jn = max(jn,orbqn_f(2,i)%j)				
			end do
			maxjt = (jn+jp)/2
		
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
		
		write(resultfile,"('#                  valence    Z    N')")
		write(resultfile,"('# Wavefunctions -- initial ',2i4)")np_i(1),np_i(2)
		do i = start_i,stop_i
            call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
            ji = closest2J(evenAJ,xj)
	        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			call assign_parity('i',mypar)
			
			call write_ejt(i, xe, ei,xj, xt,mypar)
			
		end do

        print*,' '
        print*,' Final state wave function file'
        print*,' Important: cannot be the same file as initial '

        call wfn_ropen_file_select(wfnfile_f,'F')
        call read_wfn_header_select(wfnfile_f,.false.,.true.,'F')	
		
        call wfn_read_nkeep(wfnfile_f,nkeep_f)
        print*,' There are ',nkeep_f,' final wave functions '
		print*,' Enter start, stop to compute density matrices '
		print*,' (Enter 0,0 to take all )'
		read(5,*)start_f,stop_f
		if(start_f < 1)start_f = 1
		if(stop_f < 1 .or. stop_f > nkeep_f)stop_f= nkeep_f
		
		if(start_f > stop_f)then
			print*,' start, stop FINAL ', start_f,stop_f
			stop
		end if
        call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,1,ei,xj,xt2)  
		
		write(resultfile,'("#")')
		write(resultfile,"('# Wavefunctions -- final   ',2i4)")np_f(1),np_f(2)
		do i = start_f,stop_f
            call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,i,xe,xj,xt2)  
            jf = closest2J(evenAJ,xj)
	        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			call assign_parity('f',mypar)
			
			call write_ejt(i, xe, ei,xj, xt,mypar)
			
		end do
!		rewind(wfnfile_f)
!		call read_wfn_header_select(wfnfile_f,.false.,.true.,'F')	
        call pocc_write_orbits_alt(resultfile)

		
		allocate(cc1bopme(max(nrhsps_f(1),nrhsps_i(1)),max(nrhsps_f(2),nrhsps_i(2))))   !originally had f, i swapped incorrectly

!........... NOW PROJECT..........................

       !
       do i = start_i,stop_i
           ! new interface - we say which vec to read, it checks
           ! KSM:  This will be very slow, only need to read head part of each vector
           call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
           ji = closest2J(evenAJ,xj)
           xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
		   ti = closest2J(evenAJ,xt)
		   ei =xe
!           if(iproc==0)print*,i,xe,xj,xt2	

		   do j = start_f,stop_f
			   
		       call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,j,xe,xj,xt2)  
	           jf = closest2J(evenAJ,xj)
	           xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			   tf = closest2J(evenAJ,xt)
			   ef = xe
	           jmax = (jf + ji)/2
	           jmin = abs( jf -ji)/2
			   
!		       if(iproc==0)then
!				   write(6,'(i6,f10.4,f6.1,f7.2,f10.6)')j,xe,xj,xt2,dovlp
!			   end if
			   call x2ydensity
	           if ( iproc == 0 .and. isoflag) then
	                 write(resultfile,*)' '
	                 write(resultfile,333)i,ei,Ji,Ti 
	                 write(resultfile,334)j,ef,Jf,Tf
	           end if
	           if ( iproc == 0 .and. .not.isoflag) then
	                 write(resultfile,*)' '
	                 write(resultfile,433)i,ei,Ji
	                 write(resultfile,434)j,ef,Jf
	           end if
	   333     format(' Initial state #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
	   334     format(' Final state   #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
	   433     format(' Initial state #',i5,' E = ',f10.5,' 2xJ   = ',i4) 
	   434     format(' Final state   #',i5,' E = ',f10.5,' 2xJ   = ',i4) 

			   do Jt = jmin,min(maxjt,jmax)
				   denmat = 0.0
   			        call ccdensity1b_coupled(Jt,Ji,Ti,Jf,Tf,Jz_i,Jz_f,zeroflag,numorbdim,denmat)
	                if(zeroflag)cycle
                    if(p2n)then
						write(resultfile,32)Jt
   32               format(' Jt = ',i3,', neutron <- proton ')
				    else
						write(resultfile,320)Jt
   320               format(' Jt = ',i3,', proton <- neutron ')						
						
					end if


			if(p2n)then
                do b = 1,numorbdim
                   do a = 1,numorbdim
                      if ( (denmat(a,b,1) /=  0.0) )write(resultfile,36)b,a,(denmat(a,b,1))

                   end do ! b
                end do   ! a
			else
              do a = 1,numorbdim
                 do b = 1,numorbdim
                    if ( (denmat(a,b,1) /=  0.0) )write(resultfile,36)a,b,(denmat(a,b,1))

                 end do ! b
              end do   ! a


			 
		    end if
36                  format(2i5, 2f11.6)
		  
					
				end do
		   end do
		   

!		   call wfn_writeeigenvec(wfnfile_f, frag2, vec2, i,xe,xj,xt2) 

		
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
end subroutine ccdensity1b_boss	

!==================================================================================
!
! modified for charge-changing
!
subroutine ccdensity1b_coupled(Jt,Ji,Ti,Jf,Tf,Jzi,Jzf,zeroflag, ndim,denmat)
    use system_parameters
	
	use spstate
	use sporbit
	use precisions
	implicit none
	
    integer Jt, Tt      ! J, T of transition operator
    integer Ji,Ti       ! initial wfn J, T 
    integer Jf,Tf       ! final wfn J, T
    integer Jzi,Jzf          ! of wfns
    integer Tzi,Tzf          ! = (Z -N )/2
    integer it
    logical zeroflag     ! flag to indicate no matrix elements possible
    integer:: ndim        ! dimensions of coupled density matrices
    real denmat(ndim,ndim,1)
    real cleb
    real cg
    integer a,b,ja,ma,jb,mb
    integer ia,ib
    integer asps,ath,bsps,bth,asgn,bsgn
    integer iphase,tsign
    real cgt
    real, parameter :: cgtol = 5.0e-6
    logical altflag
    integer :: ierr
	
	integer(4) :: nprotsps,nneutsps  ! depends on direction
	
	logical :: p2n

    denmat(:,:,:) = 0.0
	
	Tzi = np_i(1)-np_i(2)
	Tzf = np_f(1)-np_f(2)

    cg = cleb(Ji,Jzi,2*Jt,Jzf-Jzi,Jf,Jzf)
    cg = cg*sqrt(float(2*Jt+1)) /sqrt(float( (jf+1)))
    if(abs(cg) < cgtol)then
       zeroflag = .true.
       return
    else
      zeroflag = .false.
    endif
	
	if(np_f(1)-np_i(1)==1)then
		p2n=.false.
		
		if(np_f(2)-np_i(2)/=-1)then  ! error trap
			print*,' something wrong in ccdensity1b_coupled (1)'
			stop
		end if
	else
		p2n=.true.
		
		if(np_f(2)-np_i(2)/=1)then  ! error trap
			print*,' something wrong in ccdensity1b_coupled (2)'
			stop
		end if
		if(np_f(1)-np_i(1)/=-1)then  ! error trap
			print*,' something wrong in ccdensity1b_coupled (3)'
			stop
		end if
		
	end if
	
	if(p2n)then
		nprotsps = nrhsps_i(1)
		nneutsps = nrhsps_f(2)
	else
		nprotsps = nrhsps_f(1)
		nneutsps = nrhsps_i(2)
		
	end if
       do ia = 1, nprotsps  ! important: note here that label a -> proton state
		   
		  if(p2n)then
              a = rhspsqn_i(1,ia)%orb
              ja= rhspsqn_i(1,ia)%j
              ma= rhspsqn_i(1,ia)%m
		  else
              a = rhspsqn_f(1,ia)%orb
              ja= rhspsqn_f(1,ia)%j
              ma= rhspsqn_f(1,ia)%m
			  
		  endif
		  
		  
          do ib = 1, nneutsps !...and label b --> neutron state
   		     if(p2n)then
                 b = rhspsqn_f(2,ib)%orb
                 jb= rhspsqn_f(2,ib)%j
                 mb= rhspsqn_f(2,ib)%m
	             iphase = (-1)**( (ja -ma)/2)

   		     else
                 b = rhspsqn_i(2,ib)%orb
                 jb= rhspsqn_i(2,ib)%j
                 mb= rhspsqn_i(2,ib)%m 
	             iphase = (-1)**( (jb -mb)/2)
				 
   		     endif

             if( ma /= mb+ (Jzf-Jzi)) cycle !...? is this correct if we have switched order of protons and neutrons?
			 
             if(Jt > (ja + jb)/2 )cycle
             if(Jt < abs(ja - jb)/2 ) cycle

			 
             if(isoflag .and. .not. pndensities)then
                    cgt = cleb(Ti,Tzi,2,Tzf-Tzi,Tf,Tzf)*sqrt(float(2*Tt+1)/float(Tf+1))
!					print*,' CGT ',cgt
                    if(abs(cgt) < cgtol)cycle
					
					if(p2n)then
						cgt = cgt/cleb(1,1,2,-2,1,-1)
						
					else
						cgt = cgt/cleb(1,-1,2,2,1,1)
						
						
					end if
					
                denmat(a,b,1) = denmat(a,b,1) + cleb(ja,ma,jb,-mb,2*Jt,ma-mb)*iphase & 
              * cc1bopme(ia,ib) /(cg*cgt)
			  
             else
				 
				 if(p2n)then
                     denmat(a,b,1) = denmat(a,b,1)+ cleb(jb,mb,ja,-ma,2*Jt,mb-ma)*iphase & 
              *  cc1bopme(ia,ib)/cg
			  
		        else
                    denmat(a,b,1) = denmat(a,b,1)+ cleb(ja,ma,jb,-mb,2*Jt,ma-mb)*iphase & 
             *  cc1bopme(ia,ib)/cg
		        end if
			  
             endif

         end do  ! bsps
      enddo  ! asps
!	  print*,denmat
	return
end subroutine ccdensity1b_coupled

!==================================================================================
!
! base routine to compute charge-changing 1-body density matrices
! one issue to watch is truncation
!
! changes a proton to a neutron or vice verse
!
subroutine x2ydensity
	use jumps1body
	use sectors
	use basis_map
	use localvectors
	use w_info
	use spstate ! for testing only
	use sporbit
	use system_parameters
	
	implicit none
	
	
	integer :: ipsj ,insj   ! label of proton sector jump
	integer :: isp,fsp,ics,fcs      ! label of initial/final proton sector, conjugate neutron sector(s)
	integer :: isn,fsn  ! label of initial/final neutron sectors
	integer :: ic
	integer :: ncsi,ncsf     ! # of initial, final neutron sectors
	integer :: ipjmp,injmp   ! label of proton jumps
	
	integer(kind=basis_prec) :: isdp, fsdp   ! initial, final proton sds
	integer(kind=basis_prec) :: isdn, fsdn   ! initial, final neutron sds
	                                           ! in principle insdi=insdf, but may be different
											   ! if bases are different
	integer :: pphase,nphase
	integer :: iopp,iopn       ! label of one-body proton operator
    integer(kind=basis_prec) :: ipstate,fpstate,instate,fnstate
	real(8) :: xme
	integer :: iop1bodyp,iop1bodyn
	integer :: ac,bd    ! proton creation/annihilation operators
	integer(4) :: Wpi,Wpf,Wni,Wnf,dJzp,dJzn,parp,parn
	integer(4) :: dJz,dparity
	integer(4) :: xopp,xopn ! labels
	
!........ SET UP CHANGE OF QUANTUM NUMBERS AS CHECK....	
    dParity = parmult(iparity_i,iparity_f)
	dJz = jz_f - jz_i  !check if this is correct, otherwise switch
	cc1bopme=0.d0
	do ipsj = 1,nsectorjumps(1)  ! loop over proton sector hops
		isp = psecthop(ipsj)%isector
		fsp = psecthop(ipsj)%fsector  
		iop1bodyp = psecthop(ipsj)%iop
		if(iop1bodyp==0)then
			print*,'  Empty proton hop sector ',ipsj,iop1bodyp
			cycle
		end if
		wpi = xsd_i(1)%sector(isp)%wx
		wpf = xsd_f(1)%sector(fsp)%wx
		dJzp = psecthop(ipsj)%dJz
		parp= psecthop(ipsj)%dparity
		do insj = 1,nsectorjumps(2)  ! loop over neutron sector hops
				isn = nsecthop(insj)%isector
				fsn = nsecthop(insj)%fsector 
!.................. check if properly conjugate
                if(.not.isconjugate(.true.,isp,isn))cycle
                if(.not.isconjugate(.false.,fsp,fsn))cycle
				
				iop1bodyn = nsecthop(insj)%iop
				if(iop1bodyn==0)then
					print*,'  Empty neutron hop sector ',insj,iop1bodyn
					cycle
				end if
!.................. check quantum numbers........				
				wni = xsd_i(2)%sector(isn)%wx
				wpf = xsd_f(1)%sector(fsn)%wx
				dJzn = nsecthop(insj)%dJz
				parn= nsecthop(insj)%dparity
				
!				print*,insj,ipsj, dJzp+dJzn,dJz
				if(wpi + wni > maxWtot_i)cycle									   									
				if(wpf + wnf > maxWtot_f)cycle									   								   					
				if(parn*parp /= dparity)cycle
				if(dJzp+dJzn /=dJz)cycle
				
				do ipjmp = 1,psecthop(ipsj)%njumps ! loop over proton hops
					isdp = psecthop(ipsj)%isd(ipjmp)
					fsdp = psecthop(ipsj)%fsd(ipjmp)
					if(fsdp==-1)cycle
					ipstate = pstart_i(isdp)
					fpstate = pstart_f(fsdp)
					pphase= psecthop(ipsj)%phase(ipjmp)  ! fetch phase from applying the operator
					iopp   = psecthop(ipsj)%iXop(ipjmp)  ! find index of one-body operator
					bd = phopop(iop1bodyp)%Xop(iopp,1)
					ac = phopop(iop1bodyp)%Xop(iopp,2)
					xopp = ac+bd  ! oneof these is zero but to reduce switching I use both

					
					do injmp = 1,nsecthop(insj)%njumps  ! loop over neutron hops
						isdn = nsecthop(insj)%isd(injmp)
						fsdn = nsecthop(insj)%fsd(injmp)
						if(fsdn==-1)cycle

						instate = nstart_i(isdn)
						fnstate = nstart_f(fsdn)
						nphase= nsecthop(insj)%phase(injmp)  ! fetch phase from applying the operator
						iopn   = nsecthop(insj)%iXop(injmp)  ! find index of one-body operator
						bd = nhopop(iop1bodyn)%Xop(iopn,1)
						ac = nhopop(iop1bodyn)%Xop(iopn,2)
						xopn= ac+bd		! oneof these is zero but to reduce switching I use both
						xme = vec1(ipstate+instate)*vec2(fpstate+fnstate)*pphase*nphase

						cc1bopme(xopp,xopn)=cc1bopme(xopp,xopn)+xme  
					end do ! injmp
			end do ! ipjmp
		end do ! insj
	end do ! ipsj
	
	return
end subroutine x2ydensity
!==================================================================================
!  check if jsn is a conjugate sector of jsp
!
logical function isconjugate(init_flag,jsp,jsn)
use sectors
logical :: init_flag

integer :: jsp,jsn
integer :: ics
type (mastersect), pointer :: xsd(:)
	logical success

if(init_flag)then
	xsd => xsd_i
else
	xsd => xsd_f
end if

isconjugate = .false.

do ics = 1,xsd(1)%sector(jsp)%ncsectors
	if(jsn == xsd(1)%sector(jsp)%csector(ics))then
		isconjugate=.true.
		exit
	end if
	
end do

return 
end function isconjugate
	
end module dens1bodycc
	
	
	