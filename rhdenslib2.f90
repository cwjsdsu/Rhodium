!================================================================2
!
! modified from BIGSTICK file BDENSLIB1.f90
!
! master routines for computing density matrices 
!                     and applying one-body operators
!  (reorganized in 7.6.7)
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

module apply1body
	use dens1body
	
contains

!============================================================
!
!  CALLED BY main routine 
!

  subroutine readin1bodyop
  use opmatrixelements
  use sporbit
  use nodeinfo
  use bmpi_mod
  use butil_mod
  
  use system_parameters
  implicit none

  logical success
  character*15 :: appfile
  character*40 :: title
  integer ilastap
  integer i,i1,i2
  integer a,b
  real xj
  real :: xxx,yyy
 
!-------------- DUMMIES FOR CONSISTENCY CHECK --------------
  integer n,l,j
  integer :: numorbdim
  integer :: aerr,ierr
  logical :: xpn_flag  
  
!------------- READ IN OPERATOR FILE AND DECOUPLE

if(iproc == 0)then
  success = .false.
  do while(.not.success)  
     write(6,*)' Enter name of .opme  file '
     read(5,'(a)')appfile
     ilastap = index(appfile,' ') -1
     open(unit = 23,file=appfile(1:ilastap)//'.opme',status = 'old',err=301)
     success = .true.
     cycle
301  continue
     print*,'File ',appfile(1:ilastap),'.opme does not exist '
   enddo
!------- READ IN HEADER???
   read(23,'(a)')title
   print*,title
   
   select case (title(1:3))
      case ('iso')
	  pnoperators = .false.
	  xpn_flag =.false.
	  print*,' Reading operator file in isospin format '
  
	  case ('pns')
	  pnoperators = .true.
	  xpn_flag = .false.
	  print*,' Reading operator in two-column pns format'
	  
	  case('xpn')
	  pnoperators = .true.
	  xpn_flag = .true.
	  print*,' Reading operator in one-column xpn format '
	
      case default
	  print*,'  Something wrong with first line of .opme file '
	  print*,title(1:3)
	  print*,title
	  stop
   end select
end if

!.......... SET UP OPERATOR ARRAYS.................

numorbdim=max(numorb_i(1),numorb_i(2) )
numorbdim=max(numorbdim,numorb_f(1))
numorbdim=max(numorbdim,numorb_f(2))

call BMPI_BCAST(pnoperators,1,0,icomm,ierr)
if(pnoperators)then
  if(np_i(1)> 0)then
	  
     allocate ( pop1bod( numorbdim,numorbdim ), stat=aerr)
     if(aerr /= 0) call memerror("readin1bodyop -- p")
     pop1bod = 0.0
  end if
  if(np_i(2)> 0)then
	  
     allocate ( nop1bod( numorbdim,numorbdim ), stat=aerr)
     if(aerr /= 0) call memerror("readin1bodyop -- n")
     nop1bod = 0.0
  end if   	
	
else
  allocate ( op1bod( numorbdim,numorbdim ), stat=aerr)
  if(aerr /= 0) call memerror("readin1bodyop -- iso")
  op1bod = 0.0	
end if

if(iproc==0)then
   
   
!------------- CHECK FOR CONSISTENCY----------------
!              CHECK THEY MATCH INITIAL SPACE; may modify later
!
   read(23,*)n
   if(n /= numorb_i(1))then
     print*,' Mismatch in single particle orbits ',n,numorb_i(1) 
     stop
   endif

   do i = 1, numorb_i(1)
	   if(xpn_flag)then
		   read(23,*),i1,i2,n,l,xj
		   if(i2 /= i1+numorb_i(1))then
			   print*,' Some problem in interpreting proton, neutron orbit labels in xpn'
			   print*,i1,i2,numorb_i(1)
			   stop
		   end if
	   else
         read(23,*)i1,n,l,xj
	   end if
     if(i1 /= i)then
         print*,' ooopsie ',i1,i
         stop
     endif
     j = nint(2*xj)
     if( n/= orbqn_i(1,i)%nr .or. l /= orbqn_i(1,i)%l .or. j /= orbqn_i(1,i)%j)then
        print*,i,' Mismatch in (initial) orbits ',n,l,j
        print*, orbqn_i(1,i)%nr, orbqn_i(1,i)%l, orbqn_i(1,i)%j
        stop
     endif

   enddo
!---------- READ IN J,T for operator----------------------
  if(.not.pnoperators)then
    read(23,*)jop,top      ! J and T of operator
    write(6,*)' Operator has J, T = ',jop,top
  else
      read(23,*)jop     ! J and T of operator
      write(6,*)' Operator has J ',jop
	  top = 2
  end if

!--- MKGK option - only for r^2 Y00 ----
!--- for strength only to excited states and not back to the gs
  if (jop == 0 .and. subtract_enabled) then 
     subtract = .true.
  else
     subtract = .false.
  end if

!------------ READ IN REDUCED MATRIX ELEMENTS ----------------
  do i = 1, 1000
	  if(pnoperators)then
		  if(xpn_flag)then
			 read(23,*,end=348)a,b,xxx
			 yyy = 0.0
			 if(a > numorb_i(1) .and. b > numorb_i(1))then ! this is a neutron matrix element
				 a = a-numorb_i(1)
				 
				 b = b-numorb_i(1)
				 yyy = xxx
				 xxx = 0.0
				 
			 end if
			 if( (a > numorb_i(1) .and. b .le. numorb_i(1)) .or. & 
			  (a < numorb_i(1)+1 .and. b > numorb_i(1)))then
			      print*,' some mismatch in having proton-neutron matrix elements '
			      print*,a,b,numorb_i(1)
			      stop
		      endif
			 
		  else
			  
		     read(23,*,end=348)a,b,xxx,yyy
		   end if
	  else
        read(23,*,end =348)a,b,xxx
	  end if
!------------- ERROR TRAP---------------------
    if( orbqn_i(1,a)%j+orbqn_i(1,b)%j < 2*jop .or. abs( orbqn_i(1,a)%j-orbqn_i(1,b)%j) > 2*jop) then
       if(iproc == 0)then
         write(6,*)' Mismatch in angular momentum of operator '
         write(6,*)' one-body states : ',a,b
         write(6,*)orbqn_i(1,a)%j/2.,orbqn_i(1,b)%j/2. , jop

       endif
       stop
    endif
	if(pnoperators)then
		if(np_i(1)>0)pop1bod(a,b)=xxx
		if(np_i(2)>0)nop1bod(a,b)=yyy
	else
        op1bod(a,b) = xxx
    end if
  enddo

348      continue

  close(unit=23)
end if
!--------------- NOW BROADCAST -------------------------------
!    a bit of a kludge
call BMPI_BCAST(jop,1,0,icomm,ierr)
call BMPI_BCAST(top,1,0,icomm,ierr)
do a = 1,numorb_i(1)
	do b = 1,numorb_i(1)
		if(pnoperators)then
            if(np_i(1)> 0)call BMPI_BCAST(pop1bod(a,b),1,0,icomm,ierr)
            if(np_i(2)> 0)call BMPI_BCAST(nop1bod(a,b),1,0,icomm,ierr)
			
			
		else
           call BMPI_BCAST(op1bod(a,b),1,0,icomm,ierr)
	   
       end if
    end do
end do
!print*,iproc,' op1bod mes ',op1bod
  return
  end subroutine readin1bodyop
!=====================================================
!  automatically constructs matrix elements for J+ raising operator

subroutine raising_op
    use opmatrixelements
    use sporbit
    use nodeinfo
    use bmpi_mod
    use butil_mod
  
    use system_parameters
	implicit none
	integer :: aerr
	integer :: numorbdim
	integer :: a
	integer :: jj
	real :: jx
	pnoperators = .true.
	jop = 1
	numorbdim=max(numorb_i(1),numorb_i(2) )
	numorbdim=max(numorbdim,numorb_f(1))
	numorbdim=max(numorbdim,numorb_f(2))

	  
	     allocate ( pop1bod( numorbdim,numorbdim ), stat=aerr)
	     if(aerr /= 0) call memerror("raising_op -- p")
	     pop1bod = 0.0
   	  if(np_i(1)> 0)then
   	     do a = 1,numorb_i(1)
   		     jj = orbqn_i(1,a)%j
   		     pop1bod(a,a)= -0.5*sqrt(float(jj*(jj+1)*(jj+2)))
   	     end do
	 endif
	  
	     allocate ( nop1bod( numorbdim,numorbdim ), stat=aerr)
	     if(aerr /= 0) call memerror("raising_op -- n")
	     nop1bod = 0.0
   	  if(np_i(2)> 0)then
   	     do a = 1,numorb_i(2)
   		     jj = orbqn_i(2,a)%j
   		     nop1bod(a,a)= -0.5*sqrt(float(jj*(jj+1)*(jj+2)))
   	     end do
	  end if   
	
	  return
	  
end subroutine raising_op

!=====================================================

!
!
! CALLED BY main routine
!
!

  subroutine decouple1bodyop

  use opmatrixelements
  use onebodypot
  use sporbit
  use spstate
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

!------------- ALLOCATE UNCOUPLED ARRAYS ------------------


allocate(p1bopme(nrhsps_i(1),nrhsps_f(1)))
allocate(n1bopme(nrhsps_i(2),nrhsps_f(2)))
!  allocate(pPOT_h(nrhsps_f(1), nrhsps_i(1)), stat=aerr)
!  if(aerr /= 0) call memerror("decouple1bodyop 1")
!  allocate(nPOT_h(nrhsps_f(2), nrhsps_i(2)), stat=aerr)
!  if(aerr /= 0) call memerror("decouple1bodyop 2")

 ! pPOT_h = 0.0
 ! nPOT_h = 0.0
 p1bopme = 0.0
 n1bopme = 0.0
 
 if(np_i(1)/=np_f(1) .or. np_i(2) /= np_f(2))then
	 if(iproc==0)print*,' Mismatch in particle numbers ',np_i,np_f
	 call BMPI_Abort(icomm,101,aerr)
	 stop
	 
 end if
!--------------- PROTONS ------------------------
  do ia = 1, nrhsps_f(1)
     a = rhspsqn_f(1,ia)%orb
     ja= rhspsqn_f(1,ia)%j
     ma= rhspsqn_f(1,ia)%m


     do ib = 1, nrhsps_i(1)
        b = rhspsqn_i(1,ib)%orb
        jb= rhspsqn_i(1,ib)%j
        mb= rhspsqn_i(1,ib)%m
        if( ma /= mb+ (Jz_f-Jz_i)) cycle

        if(pnoperators )then			
			if(np_i(1)==0)cycle
            if(pop1bod(a,b) == 0.0)cycle
			
 !           if( ma /= mb) cycle
            if(jop > (ja + jb)/2 )cycle
            if(jop < abs(ja - jb)/2 ) cycle
            iphase = (-1)**( (jb -mb)/2)

            p1bopme(ia,ib) = iphase*cleb(ja,ma,jb,-mb,2*jop,ma-mb)* & 
                   pop1bod(a,b)/sqrt(2.*jop+1.)			
		else			
           if(op1bod(a,b) == 0.0)cycle
!           if( ma /= mb) cycle
           if(jop > (ja + jb)/2 )cycle
           if(jop < abs(ja - jb)/2 ) cycle
           iphase = (-1)**( (jb -mb)/2)

           p1bopme(ia,ib) = iphase*cleb(ja,ma,jb,-mb,2*jop,ma-mb)* & 
                  op1bod(a,b)/sqrt(2.*top+1)/sqrt(2.*jop+1.)*cleb(1,1,1,-1,2*top,0)				  
				  
	    end if
     enddo  ! bsps
  enddo   !asps
 
!-------------------- NEUTRONS------------------
  do ia = 1, nrhsps_f(2)


     a = rhspsqn_f(2,ia)%orb
     ja= rhspsqn_f(2,ia)%j
     ma= rhspsqn_f(2,ia)%m


     do ib = 1, nrhsps_i(2)
        b = rhspsqn_i(2,ib)%orb
        jb= rhspsqn_i(2,ib)%j
        mb= rhspsqn_i(2,ib)%m
        if( ma /= mb+ (Jz_f-Jz_i)) cycle  ! seems to be crucial
		
        if(pnoperators)then
			if(np_i(2)==0)cycle
            if(nop1bod(a,b) == 0.0)cycle
			
!            if( ma /= mb) cycle ! this is wrong
            if(jop > (ja + jb)/2 )cycle
            if(jop < abs(ja - jb)/2 ) cycle
            iphase = (-1)**( (jb -mb)/2)

            n1bopme(ia,ib) = iphase*cleb(ja,ma,jb,-mb,2*jop,ma-mb)* & 
                   nop1bod(a,b)/sqrt(2.*jop+1.)
			
		else
           if(op1bod(a,b) == 0.0)cycle
   !        if( ma /= mb) cycle
           if(jop > (ja + jb)/2 )cycle
           if(jop < abs(ja - jb)/2 ) cycle
!           iphase = (-1)**( (jb -mb)/2)
           iphase = (-1)**( (jb -mb)/2)*(-1) ! for neutrons

           n1bopme(ia,ib) = iphase*cleb(ja,ma,jb,-mb,2*jop,ma-mb)* & 
                 op1bod(a,b)/sqrt(2.*top+1)/sqrt(2.*jop+1.)*cleb(1,-1,1,1,2*top,0)				  
	    end if


     enddo  ! bsps
  enddo   !asps

  return
  end subroutine decouple1bodyop

!=====================================================
!==================================================================================
!
! base routine to compute proton 1-body density matrices
! one issue to watch is truncation
!
subroutine apply_p1b
	use jumps1body
	use sectors
	use basis_map
	use localvectors
	use w_info
	use spstate ! for testing only
	use sporbit,only:sameweightscheme
	use system_parameters
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
	
	do ipsj = 1,nsectorjumps(1)  ! loop over proton sector jumps
		is = p1bsectjump(ipsj)%isector
		fs = p1bsectjump(ipsj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
!										print*,'sectors ',is,fs,p1bsectjump(ipsj)%dJz,p1bsectjump(ipsj)%dW
		ncsi = 	xsd_i(1)%sector(is)%ncsectors				
		iop1body = 			p1bsectjump(ipsj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',ipsj,iop1body
			cycle
!		else
!			print*,' non empty sector ',ipsj,iop1body
		end if
!		ncsf = 	xsd_f(1)%sector(fs)%ncsectors							
        wpf = p1bsectjump(ipsj)%wXf
										
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

!					xme = vec1(ipstate+instate)*vec2(fpstate+fnstate)*pphase*nphase
!					print*,'(b)',ipstate+instate,fpstate+fnstate,xme
!					p1bopme(ac,bd)=p1bopme(ac,bd)+xme
                    vec2(fpstate+fnstate)=vec2(fpstate+fnstate)+ vec1(ipstate+instate)*pphase*nphase*p1bopme(ac,bd)
					
					
				end do
				
			end do
			
		end do ! ipjmp
		
		
		
	end do ! ipsj
!	print*,p1bopme
!xme = 0.0
!do ic = 1,nrhsps_i(1)
!	xme = xme + p1bopme(ic,ic)
!end do
!print*,' trace = ',xme
	
	return
end subroutine apply_p1b
!==================================================================================
!
! base routine to compute neutron 1-body density matrices
! one issue to watch is truncation
!
subroutine apply_n1b
	use jumps1body
	use sectors
	use basis_map
	use localvectors
	use spstate ! for testing only
	use sporbit,only:sameweightscheme
	use W_info,only:maxWtot_f
	use system_parameters
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
	
	
	do insj = 1,nsectorjumps(2)  ! loop over proton sector jumps
		is = n1bsectjump(insj)%isector
		fs = n1bsectjump(insj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
		ncsi = 	xsd_i(2)%sector(is)%ncsectors				
		iop1body = 			n1bsectjump(insj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',insj,iop1body
			cycle

		end if
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
                     vec2(fpstate+fnstate)=vec2(fpstate+fnstate)+ vec1(ipstate+instate)*pphase*nphase*n1bopme(ac,bd)			
				end do
				
			end do
			
		end do ! ipjmp
		
		
		
	end do ! ipsj

	
	return
end subroutine apply_n1b
!===========================================
subroutine apply1b_boss(raisin,tempest)
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
		logical,intent(in) :: raisin ! if TRUE then apply raising operator
		logical,intent(in) :: tempest ! if TRUE then write to temp file
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
		logical :: smallflag
		integer :: mypar
		
		if(np_i(1)/=np_f(1) .or. np_i(2)/=np_f(2))then
			if(iproc==0)then
				print*,' Not number conserving '
			    print*,np_i, np_f
			end if
			call BMPI_Abort(icomm,101,ierr)
			stop
		end if
		call parity_boss
		
!............ CHECK THAT Mf = Mi + 1
		
	    if( mod(np_i(1)+np_i(2),2) == 1 )then
			
	         evenAJ = .false.
	       evenAT = .false.
	    else
	       evenAJ = .true.
	       evenAT = .true.
	    end if
		
		call phasechecker
		call basismapper(1)
		call basismapper(2)
		
		call set_nfragments(1)
	    frag1= nodal(iproc)%ifragment
	    frag2= nodal(iproc)%ffragment		
		call setup_localvectors_both
		
		if(raisin)then
			call raising_op
		else
		    call readin1bodyop
	    end if
		call decouple1bodyop
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
		
		if(tempest)then
			start_i_default = start_i
			stop_i_default = stop_i	
		end if
		
        call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,1,ei,xj,xt2)  
		
		write(resultfile,"('# Wavefunctions -- initial ',2i4)")np_i(1),np_i(2)
		do i = start_i,stop_i
            call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
			call assign_parity('i',mypar)
			
            ji = closest2J(evenAJ,xj)
	        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			
			call write_ejt(i, xe, ei,xj, xt,mypar)
			
		end do

		
!........... NOW PROJECT..........................
!........... OPEN FINAL .WFN FILE AND WRITE HEADER.....
       if(tempest)then
		   call wfn_wopen_file_TEMP(wfnfile_f)
		   
	   else
		   
          write(6,*)' Enter name of output .wfn file ' 
	      read(5,'(a)')outfile
	   
          call wfn_wopen_file(wfnfile_f,.false.)
       end if
	   
       call write_wfn_header(wfnfile_f)
	   
	   call wfn_write_nkeep(nkeep_f) ! write number of vectors to wfn file
	   

!........... NOW PROJECT..........................
if(iproc==0)then
	if(raisin)then
		print*,' Zero vector means a vector with J = 0, cannot be raised '
	else
		print*,' Zero vector means a vector with zero normalization after application '
	end if
end if
       do i = start_i,stop_i
           ! new interface - we say which vec to read, it checks
           ! KSM:  This will be very slow, only need to read head part of each vector
           call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
		   vec2= 0.0
		   if(np_i(1) > 0)call apply_p1b
		   if(np_i(2)> 0 )call apply_n1b
		   
		   if(raisin)then ! must normalize...and maybe multiply by -1?
			   call dnormvec_p('n','f',dovlp,smallflag)
	!		   print*,' Normalization is ',dovlp,', expect ',sqrt( xj*(xj+1.) -0.5*Jz_i*(0.5*Jz_i+1.0) )
		   end if
		   call wfn_writeeigenvec(wfnfile_f, frag2, vec2, i,xe,xj,xt2) 
		
	   end do
	   
!	   if(.not.tempest)then
	     print*,' closing... '
	     call wfn_close_file(wfnfile_f)
	     print*,' closed '
!	 else
!		 print*,' wfn ',wfnfile_f, ' not closed'
 !      end if

!

return
end subroutine apply1b_boss	

!==================================================================================
!==================================================================================

end module apply1body
