!================================================================2
!
!
!  routines for  applying one-body, charge-changing operators
!=================================================================
!
! these routines read in a one-body operator and apply it to a vector
! operators are read in as reduced matrix elements, using the conventions
! from Edmonds "Angular momentum in quantum mechanics," Eq. (5.4.1)
!
!(Jf || O_K || Ji) = [Jf] ( Jf Mf, KM | Ji Mi)^-1 ( Jf Mf | O_KM | Ji Mi)
!
!Special case: number operator:  One can easily see that 
!rho(a^+a) = sqrt(2 j_a + 1)/ sqrt( 2 Jf +1 ) n(a)
!If one includes isospin then
!rho(a^+a) = sqrt(2*(2 j_a +1 ))/ [ sqrt(2 Jf +1) sqrt(2 Tf +1)] * n(a)
!
!
!=====================================================

module apply1bodycc
	use dens1bodycc
	use apply1body
	
contains


!=====================================================

!
!
! CALLED BY main routine
!
!

  subroutine decouple1bodyop_cc

  use opmatrixelements
  use onebodypot
  use sporbit
  use spstate
  use precisions
  use nodeinfo
  use system_parameters
  use bmpi_mod

  implicit none

  integer a,b
  integer ja,jb
  integer ma,mb
  integer ia,ib
  integer asps, bsps
  integer asgn, bsgn, ath, bth
  integer iphase
  real   cleb
  integer :: aerr
  
  real cg

  integer :: ierr

integer(4) :: nprotsps,nneutsps  ! depends on direction

logical :: p2n
  
!------------- ALLOCATE UNCOUPLED ARRAYS ------------------

	if(np_f(1)-np_i(1)==1)then
		p2n=.false.
		
		if(np_f(2)-np_i(2)/=-1)then  ! error trap
			print*,' something wrong in decouple1bodyop_cc (1)'
			stop
		end if
	else
		p2n=.true.
		
		if(np_f(2)-np_i(2)/=1)then  ! error trap
			print*,' something wrong in decouple1bodyop_cc (2)'
			stop
		end if
		if(np_f(1)-np_i(1)/=-1)then  ! error trap
			print*,' something wrong in decouple1bodyop_cc (3)'
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
	allocate(cc1bopme(nprotsps,nneutsps))
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

             if( p2n .and.  mb /= ma+ (Jz_f-Jz_i)) cycle   ! not sure about this for swapping protons, neutrons, need to check
             if( .not.p2n .and.  ma/= mb+ (Jz_f-Jz_i)) cycle   ! not sure about this for swapping protons, neutrons, need to check
			 
             if(Jop > (ja + jb)/2 )cycle
             if(Jop < abs(ja - jb)/2 ) cycle

				 
				 if(p2n)then
		             cc1bopme(ia,ib) = iphase*cleb(jb,mb,ja,-ma,2*Jop,mb-ma)* & 
		                    op1bod(a,b)/sqrt(3.)/sqrt(2.*jop+1.)*sqrt(0.5)
		        else
	             cc1bopme(ia,ib) = iphase*cleb(ja,ma,jb,-mb,2*jop,ma-mb)* & 
	                    op1bod(a,b)/sqrt(3.)/sqrt(2.*jop+1.)*sqrt(0.5)
		        end if
			  

         end do  ! bsps
      enddo  ! asps

!	  print*, ' cc1bopme ',cc1bopme

  return
  
  end subroutine decouple1bodyop_cc

!=====================================================
!==================================================================================
!
! base routine to compute proton 1-body density matrices
! one issue to watch is truncation
!
subroutine apply_x2y
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
	!
!........ SET UP CHANGE OF QUANTUM NUMBERS AS CHECK....	
    dParity = parmult(iparity_i,iparity_f)
	dJz = jz_f - jz_i  !check if this is correct, otherwise switch
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
						xme = cc1bopme(xopp,xopn)*pphase*nphase
						vec2(fpstate+fnstate)=vec2(fpstate+fnstate)+ vec1(ipstate+instate)*xme

					end do ! injmp
			end do ! ipjmp
		end do ! insj
	end do ! ipsj
	
	return
end subroutine apply_x2y
!==================================================================================

subroutine apply1b_cc_boss
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
		use nodeinfo
		use bmpi_mod
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
		
		logical :: p2n
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
		
	    if( mod(np_i(1)+np_i(2),2) == 1 )then
			
	         evenAJ = .false.
	       evenAT = .false.
	    else
	       evenAJ = .true.
	       evenAT = .true.
	    end if
		call parity_boss
		call phasechecker
!		call basismapper(1)  ! this is only called when one species or another can go unchanged
!		call basismapper(2)
		
		call set_nfragments(1)
	    frag1= nodal(iproc)%ifragment
	    frag2= nodal(iproc)%ffragment		
		call setup_localvectors_both
		
		call readin1bodyop
		call decouple1bodyop_cc
!............ OPEN INITIAL .WFN FILE AND READ HEADER.....

		print*,' Initial state wave function file'

		call wfn_ropen_file_select(wfnfile_i,'I')
		call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')	
		call wfn_read_nkeep(wfnfile_i,nkeep_i)
		print*,' There are ',nkeep_i,' initial wave functions '		
		print*,' Enter start, stop to apply 1-body operator '
		print*,' (Enter 0,0 to take all )'
		read(5,*)start_i,stop_i
		if(start_i < 1)start_i = 1
		if(stop_i < 1 .or. stop_i > nkeep_i)stop_i = nkeep_i
		if(start_i > stop_i)then
			print*,' start, stop INITIAL ', start_i,stop_i
			stop
		end if
		
		nkeep_f = stop_i-start_i+1
		
		if(nkeep_f < 1)then
			
			print*, ' wrong number to keep ',nkeep_f
			stop
		end if
		
        call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,1,ei,xj,xt2)  
		
		
		write(resultfile,"('# Wavefunctions -- initial ',2i4)")np_i(1),np_i(2)
		do i = start_i,stop_i
            call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
            ji = closest2J(evenAJ,xj)
	        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			call assign_parity('i',mypar)
			
			call write_ejt(i, xe, ei,xj, xt,mypar)
			
		end do
		
!........... NOW PROJECT..........................
!........... OPEN FINAL .WFN FILE AND WRITE HEADER.....

       write(6,*)' Enter name of output .wfn file ' 
	   read(5,'(a)')outfile
	   
       call wfn_wopen_file(wfnfile_f,.false.)
       call write_wfn_header(wfnfile_f)
	   
	   call wfn_write_nkeep(nkeep_f) ! write number of vectors to wfn file
	   

!........... NOW PROJECT..........................
       do i = start_i,stop_i
           ! new interface - we say which vec to read, it checks
           ! KSM:  This will be very slow, only need to read head part of each vector
           call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
           if(iproc==0)print*,i,xe,xj,xt2	
		   vec2= 0.0
		   call apply_x2y

		   call wfn_writeeigenvec(wfnfile_f, frag2, vec2, i,xe,xj,xt2) 
		
	   end do
	   
	   print*,' closing... '
	   call wfn_close_file(wfnfile_f)
	   print*,' closed '

!

return
end subroutine apply1b_cc_boss	

!==================================================================================
!==================================================================================

end module apply1bodycc
