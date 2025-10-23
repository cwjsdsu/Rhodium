!
!  data and routines for 1-body spectroscopic factors
!
!  based upon rhdenslib1.f90
!
!  added 0.8.1  March 2020
!
module spfac
	use precisions
	use dens1body
	use hops
	implicit none
	
    real(kind=obs_prec), allocatable,target   :: pspopme(:), nspopme(:) ! M-scheme
    real(kind=obs_prec), allocatable,target   :: psamp(:),nsamp(:) ! J-scheme amplitudes
	
	logical, allocatable :: xsp_used(:)
	
contains

!
!  SUBROUTINES CALLED
!
!  CALLED BY: main_menu
!	
subroutine spfac_boss
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
	    integer(4) :: ji,ti,jf,tf,jt2,tt,j2min,j2max
	    integer(kind=basis_prec) :: k
	    logical :: evenAJi,evenATi    ! needed to force "correct" J, T
	    logical :: evenAJf,evenATf    ! needed to force "correct" J, T

		logical :: zeroflag
	    real, allocatable :: spfacmat(:)
		integer :: numorbdim
		integer :: a,b
		integer :: start_i, stop_i,start_f,stop_f
		logical :: addflag
		integer :: it,mypar
!--------- CAN TEST BY SUM RULES----------
!  if ADDING a particle, then 
! sum_i  |A_a^if|^2 = < f | n_a  | f >
! sum_f (2Jf+1)/(2Ji+1) |A_a^if|^2 = < i | n_a | i> 
!
! if REMOVING a particle, then reversed:
! sum_f  |A_a^if|^2 = < i | n_a  | i >
! sum_i (2Ji+1)/(2Jf+1) |A_a^if|^2 = < f | n_a | f > 
!		
		real(4), allocatable :: sumrule_i(:,:),sumrule_f(:,:)
		real(4) :: numpart

        print*,' Need to check number '		

	    numorbdim=max(numorb_i(1),numorb_i(2) )
		numorbdim=max(numorbdim,numorb_f(1))
		numorbdim=max(numorbdim,numorb_f(2))

	    allocate( spfacmat( numorbdim ) )
	    if( mod(np_i(1)+np_i(2),2) == 1 )then
	       evenAJi = .false.
	       evenATi = .false.
	    else
	       evenAJi = .true.
	       evenATi = .true.
	    end if
		evenAJf = .not. evenAJi
		evenATf = .not. evenATi
		call parity_boss
		call phasechecker
		if(dnp(1)==0)call basismapper(1)
		if(dnp(2)==0)call basismapper(2)
		
		call set_nfragments(1)
	    frag1= nodal(iproc)%ifragment
	    frag2= nodal(iproc)%ffragment		
		call setup_localvectors_both
		
		if(dnp(1) /= 0)then
			
			it= 1

		else
			it =2 
		end if
		if(dnp(it)==1)then
			addflag = .true.
		else
			addflag = .false.
		end if
			
!............ OPEN INITIAL .WFN FILE AND READ HEADER.....

		print*,' Initial state wave function file'

		call wfn_ropen_file_select(wfnfile_i,'I')
		call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')	
		call wfn_read_nkeep(wfnfile_i,nkeep_i)
		print*,' There are ',nkeep_i,' initial wave functions '		
		print*,' Enter start, stop to compute spectroscopic factors '
		print*,' (Enter 0,0 to take all )'
		read(5,*)start_i,stop_i
		if(start_i < 1)start_i = 1
		if(stop_i < 1 .or. stop_i > nkeep_i)stop_i = nkeep_i
		if(start_i > stop_i)then
			print*,' start, stop INITIAL ', start_i,stop_i
			stop
		end if
		
        call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,1,ei,xj,xt2)  
		
		write(resultfile,"('# Wavefunctions -- initial ',2i4)")np_i(1),np_i(2)
		do i = start_i,stop_i
            call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
            ji = closest2J(evenAJi,xj)
	        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			call assign_parity('i',mypar)
			call write_ejt(i, xe, ei,xj, xt, mypar)
			
		end do

        print*,' '
        print*,' Final state wave function file'
        print*,' Important: cannot be the same file as initial '

        call wfn_ropen_file_select(wfnfile_f,'F')
        call read_wfn_header_select(wfnfile_f,.false.,.true.,'F')	
		
        call wfn_read_nkeep(wfnfile_f,nkeep_f)
        print*,' There are ',nkeep_f,' final wave functions '
		print*,' Enter start, stop to compute spectroscopic factors '
		print*,' (Enter 0,0 to take all )'
		read(5,*)start_f,stop_f
		if(start_f < 1)start_f = 1
		if(stop_f < 1 .or. stop_f > nkeep_f)stop_f= nkeep_f
		
		if(start_f > stop_f)then
			print*,' start, stop FINAL ', start_f,stop_f
			stop
		end if
        call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,1,ei,xj,xt2)  
		
		write(resultfile,*)'# '
		write(resultfile,"('# Wavefunctions -- final   ',2i4)")np_f(1),np_f(2)
		do i = start_f,stop_f
            call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,i,xe,xj,xt2)  
            jf = closest2J(evenAJf,xj)
			call assign_parity('f',mypar)
	        xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			
			call write_ejt(i, xe, ei,xj, xt,mypar)
			
		end do
		call pocc_write_orbits_alt(resultfile)
!----------------- ALLOCATE STORAGE----------------------------		
	    allocate(pspopme(max(nrhsps_f(1),nrhsps_i(1))))
		allocate(nspopme(max(nrhsps_i(2),nrhsps_f(2))))
		
		allocate(sumrule_i(start_i:stop_i,numorbdim))
		sumrule_i=0.0
		allocate(sumrule_f(start_f:stop_f,numorbdim))
		sumrule_f = 0.0

!........... NOW PROJECT..........................

       !
       do i = start_i,stop_i
           ! new interface - we say which vec to read, it checks
           ! KSM:  This will be very slow, only need to read head part of each vector
           call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
           ji = closest2J(evenAJi,xj)
		   
		   ei =xe
!           if(iproc==0)print*,i,xe,xj,xt2	

		   do j = start_f,stop_f
		       call wfn_readeigenvec_select('F',wfnfile_f,frag2, fcomm2, vec2,j,xe,xj,xt2)  
	           jf = closest2J(evenAJf,xj)
			   ef = xe
!	           j2max = (jf + ji)    !still x 2
!	           j2min = abs( jf -ji)
			   
			   if(dnp(1)/=0)call applyspfac_p
			   if(dnp(2)/=0)call applyspfac_n
!.................. DO NOT NEED TO WRITE OUT ISOSPIN AGAIN.....			   
!	           if ( iproc == 0 .and. isoflag) then
!	                 write(resultfile,*)' '
!	                 write(resultfile,333)i,ei,Ji,Ti 
!	                 write(resultfile,334)j,ef,Jf,Tf
!	           end if
!	           if ( iproc == 0 .and. .not.isoflag) then
		           if ( iproc == 0 ) then

	                 write(resultfile,*)' '
	                 write(resultfile,433)i,ei,Ji
	                 write(resultfile,434)j,ef,Jf
	           end if
	   333     format(' Initial state #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
	   334     format(' Final state   #',i5,' E = ',f10.5,' 2xJ, 2xT = ',2i4) 
	   433     format(' Initial state #',i5,' E = ',f10.5,' 2xJ   = ',i4) 
	   434     format(' Final state   #',i5,' E = ',f10.5,' 2xJ   = ',i4) 
	   500     format( ' orbit  amp ')

			   spfacmat = 0.0
   			   call spfac_coupled(it,addflag,Ji,Jf,Jz_i,Jz_f,zeroflag,numorbdim,spfacmat)
	           if(zeroflag)cycle
               write(resultfile,500)
               do a = 1,numorbdim
                    if ( (spfacmat(a) /=  0.0) )then
!						if(.not.addflag)spfacmat(a)= spfacma
                       write(resultfile,36)a,spfacmat(a)
					   if(addflag)then
						   sumrule_f(j,a)=sumrule_f(j,a)+ spfacmat(a)**2
						   sumrule_i(i,a)=sumrule_i(i,a)+ spfacmat(a)**2*(Jf+1.)/(Ji+1.)
					   else
						   sumrule_f(j,a)=sumrule_f(j,a)+ spfacmat(a)**2	
						   sumrule_i(i,a)=sumrule_i(i,a)+ spfacmat(a)**2*(Jf+1.)/(Ji+1.)						   
					   end if
					   
                    end if
36                  format(i5, f11.6)
              end do   ! a
					
					
		   end do  !j
		   

!		   call wfn_writeeigenvec(wfnfile_f, frag2, vec2, i,xe,xj,xt2) 

		
	   end do  ! i
	   
	   
if(iproc == 0)then
	write(resultfile,*)' '
   write(resultfile,*)'# Initial state sum rules ############################################### '
   do i = start_i,stop_i
	   write(resultfile,*)i,' initial state '
	   write(resultfile,*)'# orbit   < i | n(a) | i > '
	   numpart = 0.0
	   do a = 1,numorbdim
		   
		   if(addflag)then
			   sumrule_i(i,a)= orbqn_i(it,a)%j+1-sumrule_i(i,a)
		   end if
		   numpart =numpart+ sumrule_i(i,a)
		   write(resultfile,603)a,sumrule_i(i,a)
	   end do
	   write(resultfile,*)numpart,' particles total '
   end do
write(resultfile,*)' '
   
   write(resultfile,*)'# Final state sum rules ############################################### '
   do i = start_f,stop_f
	   write(resultfile,*)i,' final state '
	   write(resultfile,*)'# orbit   < f | n(a) | f > '
	   numpart = 0.0

	   do a = 1,numorbdim
		   if(.not.addflag)then
			      sumrule_f(i,a)= orbqn_f(it,a)%j+1-sumrule_f(i,a)
		   end if
		   numpart =numpart+ sumrule_f(i,a)
		   
		   write(resultfile,603)a,sumrule_f(i,a)
	   end do
	   write(resultfile,*)numpart,' particles total '
	   
   end do
603 format(i6,f11.6)   
   write(resultfile,*)' '
   write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
   write(resultfile,*)' '
   write(resultfile,*)' Definition of spectroscopic amplitudes : '
   if(addflag)then
      write(resultfile,*)' A_j^fi =   -(Jf || a^+_j || Ji) / sqrt(2Jf+1) '
   else
       write(resultfile,*)' A_j^fi =   -(Jf || a_j || Ji) / sqrt(2Jf+1) '
	   
   end if
   write(resultfile,*)'  where reduced matrix element is convention of Edmonds '
   write(resultfile,*)' '
   write(resultfile,*)'+++++++++++++++++++++++++++++++++++++++++++++'
   write(resultfile,*)' '

   write(6,*)' '
   write(6,*)'All done with spectroscopic factors'
   write(6,*)' '
end if
return
end subroutine spfac_boss	

!==================================================================================

subroutine spfac_coupled(it,createflag,Ji,Jf,Jzi,Jzf,zeroflag, ndim,spfacmat)
    use system_parameters
	
	use spstate
	use sporbit
	use precisions
	implicit none
    integer it	
	logical :: createflag  !if creations or destruction
    integer Jt, Tt      ! J, T of transition operator
    integer Ji,Ti       ! initial wfn J, T 
    integer Jf,Tf       ! final wfn J, T
    integer Jzi,Jzf          ! of wfns
    integer Tz          ! = (Z -N )/2
    logical zeroflag     ! flag to indicate no matrix elements possible
    integer:: ndim        ! dimensions of coupled density matrices
    real spfacmat(ndim)
    real cleb
    real cg
    integer a,b,ja,ma,jb,mb
    integer ia,ib
    integer asps,ath,bsps,bth,asgn,bsgn
    integer iphase,tsign
    real cgt
    real, parameter :: cgtol = 5.0e-6
    logical altflag
    real(kind=obs_prec), pointer :: xspopme(:)
    integer :: ierr


    zeroflag = .true.
    spfacmat(:) = 0.0

       if(it==1)then
          xspopme => pspopme
           tsign = 1
       else
          xspopme => nspopme
          tsign =-1
       end if
	   
!--- THIS (MAY) NEED TO BE REVISED ---	   
       do ia = 1, nrhsps_f(it)
          a = rhspsqn_f(it,ia)%orb
          ja= rhspsqn_f(it,ia)%j
          ma= rhspsqn_f(it,ia)%m
		  
		  if(createflag .and. ma /= dnp(it)*(Jzf-Jzi))cycle

	      cg = cleb(Ji,Jzi,ja,ma*dnp(it),Jf,Jzf)

	      if(abs(cg) > cgtol)then
	         zeroflag = .false.  ! check to make sure there is something to print ou t
		 else
			 cycle

	      endif
		  
		  spfacmat(a)= spfacmat(a) - xspopme(ia)/cg

      enddo  ! asps
	
	return
end subroutine spfac_coupled

!==================================================================================
!
! base routine to compute proton 1-body spectroscopic amplitudes
! one issue to watch is truncation
!
subroutine applyspfac_p
	use jumps1body
	use sectors
	use basis_map
	use localvectors
	use w_info
	use spstate ! for testing only
	use sporbit,only:sameweightscheme
	use hops
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
	integer :: xop ! = ac+bd
	integer(4) :: Wpf,Wnf
	
	pspopme=0.d0
	do ipsj = 1,nsectorjumps(1)  ! loop over proton sector jumps
		is = psecthop(ipsj)%isector
		fs = psecthop(ipsj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
		ncsi = 	xsd_i(1)%sector(is)%ncsectors				
		iop1body = 			psecthop(ipsj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',ipsj,iop1body
			cycle
		end if
        wpf = psecthop(ipsj)%wXf
										
		do ipjmp = 1,psecthop(ipsj)%njumps
			isdp = psecthop(ipsj)%isd(ipjmp)
			fsdp = psecthop(ipsj)%fsd(ipjmp)
			if(fsdp==-1)cycle
			ipstate = pstart_i(isdp)
			fpstate = pstart_f(fsdp)
			pphase= psecthop(ipsj)%phase(ipjmp)  ! fetch phase from applying the operator
			iopp   = psecthop(ipsj)%iXop(ipjmp)  ! find index of one-body operator
			bd = phopop(iop1body)%Xop(iopp,1)
			ac = phopop(iop1body)%Xop(iopp,2)
			xop = ac+bd  ! oneof these is zero but to reduce switching I use both
			
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


					xme = vec1(ipstate+instate)*vec2(fpstate+fnstate)*pphase*nphase
					pspopme(xop)=pspopme(xop)+xme
!					if(abs(xme)> 1.e-4 .and. xop==9)			print*,' huh ',xme,pphase,nphase

				end do
				
			end do
			
		end do ! ipjmp
		
		
		
	end do ! ipsj
	
	return
end subroutine applyspfac_p
!==================================================================================
!
! base routine to compute neutron 1-body spectroscopic amplitudes
! one issue to watch is truncation
!
subroutine applyspfac_n
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
	integer :: ac,bd ,xop   ! proton creation/annihilation operators
	integer(4) :: Wpf,Wnf
	
	
	nspopme=0.d0
	do insj = 1,nsectorjumps(2)  ! loop over neutron sector jumps
		is = nsecthop(insj)%isector
		fs = nsecthop(insj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
		ncsi = 	xsd_i(2)%sector(is)%ncsectors				
		iop1body = 			nsecthop(insj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',insj,iop1body
			cycle

		end if
        wnf = nsecthop(insj)%wXf
										
		do injmp = 1,nsecthop(insj)%njumps
			isdn = nsecthop(insj)%isd(injmp)
			fsdn = nsecthop(insj)%fsd(injmp)
			if(fsdn==-1)cycle
			instate = nstart_i(isdn)
			fnstate = nstart_f(fsdn)
			nphase= nsecthop(insj)%phase(injmp)  ! fetch phase from applying the operator
			iopn   = nsecthop(insj)%iXop(injmp)  ! find index of one-body operator
			bd = nhopop(iop1body)%Xop(iopn,1)
			ac = nhopop(iop1body)%Xop(iopn,2)
			xop= ac+bd
!............ LOOP OVER CONJUGATE NEUTRON SDS....................			
			
			do ic = 1,ncsi   ! loop over conjugate proton sectors
				ics = xsd_i(2)%sector(is)%csector(ic)  ! find index of conjugate proton sector
				                                       ! HOWEVER must be careful when truncating--
				                                       ! not all FINAL neutron SD will be allowed
				                                       ! if the sum of final Ws goes too high.
				if(sameweightscheme)then              ! can check this issue directly
					wpf = xsd_i(1)%sector(ics)%wx
					
					if(wpf + wnf > maxWtot_f)cycle	
					
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
					xme = vec1(ipstate+instate)*vec2(fpstate+fnstate)*pphase*nphase
					nspopme(xop)=nspopme(xop)+xme
					
				end do
				
			end do
			
		end do ! ipjmp
		

	end do ! ipsj
	
	return
end subroutine applyspfac_n
!==================================================================================

subroutine appsp_boss
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
		integer :: it
		logical :: normalize
		character :: ychar
		real(8) :: dnorm
		logical :: smallflag
		integer :: mypar
		
		it = -1
		if(abs(np_i(1)-np_f(1))==1 .and. np_i(2)==np_f(2))it = 1
		if(abs(np_i(2)-np_f(2))==1 .and. np_i(1)==np_f(1))it = 2
		
		
		if(it==-1)then
			if(iproc==0)then
				print*,' Not number changing'
			    print*,np_i, np_f
			end if
			call BMPI_Abort(icomm,101,ierr)
			stop
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
		if(dnp(1)==0)call basismapper(1)
		if(dnp(2)==0)call basismapper(2)
		
		call set_nfragments(1)
	    frag1= nodal(iproc)%ifragment
	    frag2= nodal(iproc)%ffragment		
		call setup_localvectors_both
		
		
		call readin1spamp
!.... ASK IF NORMALIZE...........
        normalize = .false.
		print*,' Do you want to normalize final wfns? (y/n)'
		read(5,*)ychar
		if(ychar=='y' .or. ychar=='Y')then
			normalize = .true.
			print*,' Normalizing final wave functions '
		else
			print*,' Not normalizing final wave functions '
			
		end if
!............ OPEN INITIAL .WFN FILE AND READ HEADER.....

		print*,' Initial state wave function file'

		call wfn_ropen_file_select(wfnfile_i,'I')
		call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')	
		call wfn_read_nkeep(wfnfile_i,nkeep_i)
		print*,' There are ',nkeep_i,' initial wave functions '		
		print*,' Enter start, stop to apply one particle creation/removal'
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
		   if(it==1)call applysphop_p
		   if(it==2)call applysphop_n
		   if(normalize)then
			   call dnormvec_p('n','f',dnorm,smallflag)
			   if(iproc==0)print*,' Normalization ',dnorm
		   end if
		   call wfn_writeeigenvec(wfnfile_f, frag2, vec2, i-start_i+1,xe,xj,xt2) ! need to make sure we start at 1
		
	   end do
	   
	   print*,' closing... '
	   call wfn_close_file(wfnfile_f)
	   print*,' closed '	

return
end subroutine appsp_boss


!==================================================================================
! routine to take a single, fixed initial wfn, and then apply/ remove a single particle
! from each of a set of designated orbitals, writing to output
!
subroutine allsp_boss
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
	    integer i,j,n, tidx,f
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
		integer :: ichoose
		integer :: it
		logical :: normalize
		character :: ychar
		real(8) :: dnorm
		logical :: smallflag
		integer :: mypar
		integer :: iorb
		
		
		it = -1
		if(abs(np_i(1)-np_f(1))==1 .and. np_i(2)==np_f(2))it = 1
		if(abs(np_i(2)-np_f(2))==1 .and. np_i(1)==np_f(1))it = 2
		
		
		if(it==-1)then
			if(iproc==0)then
				print*,' Not number changing'
			    print*,np_i, np_f
			end if
			call BMPI_Abort(icomm,101,ierr)
			stop
		end if
		
		if(it==1)then
	        allocate(xsp_used(max(nrhsps_f(1),nrhsps_i(1))))
	        allocate(pspopme(max(nrhsps_f(1),nrhsps_i(1))))
			
		else
			print*,' Sizes  ',nrhsps_f(2),nrhsps_i(2)
	        allocate(xsp_used(max(nrhsps_f(2),nrhsps_i(2))))
	        allocate(nspopme(max(nrhsps_f(2),nrhsps_i(2))))
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
		if(dnp(1)==0)call basismapper(1)
		if(dnp(2)==0)call basismapper(2)
		
		call set_nfragments(1)
	    frag1= nodal(iproc)%ifragment
	    frag2= nodal(iproc)%ffragment		
		call setup_localvectors_both
				
!.... ASK IF NORMALIZE...........
        normalize = .false.
!		print*,' Do you want to normalize final wfns? (y/n)'
!		read(5,*)ychar
!		if(ychar=='y' .or. ychar=='Y')then
!			normalize = .true.
!			print*,' Normalizing final wave functions '
!		else
!			print*,' Not normalizing final wave functions '			
!		end if
!............ OPEN INITIAL .WFN FILE AND READ HEADER.....

		print*,' Initial state wave function file'

		call wfn_ropen_file_select(wfnfile_i,'I')
		call read_wfn_header_select(wfnfile_i,.false.,.true.,'I')	
		call wfn_read_nkeep(wfnfile_i,nkeep_i)
		print*,' There are ',nkeep_i,' initial wave functions '		
		
		if(nkeep_i > 1)then
 	       if(iproc==0 )then
 			   write(6,*)' Initial Energy   J     T^2'
		   end if
	    do i = 1,nkeep_i
	       ! new interface - we say which vec to read, it checks
	       ! KSM:  This will be very slow, only need to read head part of each vector
	       call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
	       if(iproc==0 )then
	!		   write(6,*)' Initial Energy   J     T^2'
			   write(6,'(i6,f10.4,f6.1,f7.2)')i,xe,xj,xt2
		   end if
	    end do ! i
		print*,' Choose which one to apply one particle creation/removal'
!		print*,' (Enter 0,0 to take all )'
		read(5,*)ichoose
       	else
	       	ichoose = 1
	    end if

		if(ichoose < 1)ichoose = 1
		if(ichoose > nkeep_i)ichoose = nkeep_i
		
        call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,ichoose,ei,xj,xt2)  
		
		write(resultfile,"('# Wavefunctions -- initial ',2i4)")np_i(1),np_i(2)
        call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,ichoose,xe,xj,xt2)  
        ji = closest2J(evenAJ,xj)
	    xt = real( -0.5 + sqrt(xt2 + 0.25), kind=4)
			call assign_parity('i',mypar)
			
			call write_ejt(i, xe, ei,xj, xt,mypar)
			
			if(it==1)call survey_sphop_p
			if(it==2)call survey_sphop_n
			
	        call pocc_write_orbits_alt(resultfile)
			
!
!............ COUNT UP # OF POSSIBLE FINAL STATES		
		nkeep_f =0 

		do i = 1,size(xsp_used)
			if(xsp_used(i))nkeep_f = nkeep_f+1

		end do		

		
		if(iproc==0)print*,' There are ',nkeep_f,' final states possible '					
!........... NOW PROJECT..........................
!........... OPEN FINAL .WFN FILE AND WRITE HEADER.....

       write(6,*)' Enter name of output .wfn file ' 
	   read(5,'(a)')outfile
	   
       call wfn_wopen_file(wfnfile_f,.false.)
       call write_wfn_header(wfnfile_f)
	   
	   call wfn_write_nkeep(nkeep_f) ! write number of vectors to wfn file
	   

!........... NOW PROJECT..........................
       write(resultfile,*)' '
       write(resultfile,*)' FINAL STATE    ORBIT ADDED/REMOVED'
       do f = 1,nkeep_f
		   call update_spop(it,f,iorb)
		   
           ! new interface - we say which vec to read, it checks
           ! KSM:  This will be very slow, only need to read head part of each vector
 !          call wfn_readeigenvec_select('I',wfnfile_i,frag1, fcomm1, vec1,i,xe,xj,xt2)  
!           if(iproc==0)print*,i,xe,xj,xt2	
		   vec2= 0.0
!........ modify the vector		   
		   
		   if(it==1)call applysphop_p
		   if(it==2)call applysphop_n
!		   if(normalize)then
!			   call dnormvec_p('n','f',dnorm,smallflag)
!			   if(iproc==0)print*,' Normalization ',dnorm
!		   end if
		   call wfn_writeeigenvec(wfnfile_f, frag2, vec2, f,xe,xj,xt2) ! need to make sure we start at 1
		   write(resultfile,*)f,iorb
	   end do
	   write(resultfile,*)' '
	   close(resultfile)
	   print*,' closing... '
	   call wfn_close_file(wfnfile_f)
	   print*,' closed '	

return
end subroutine allsp_boss

!==================================================================================
!routine to read in amplitudes for applying/removing a single particle
!==================================================================================

  subroutine readin1spamp
!  use opmatrixelements
  use sporbit
  use spstate
  use nodeinfo
  use bmpi_mod
  use butil_mod
  
  use system_parameters
  implicit none

  logical success
  character*15 :: appfile
  character*40 :: title
  integer ilastap
  integer i,i1
  integer a,b
  real xj
  real :: xxx,yyy
  integer :: it
 
!-------------- DUMMIES FOR CONSISTENCY CHECK --------------
  integer n,l,j
  integer :: numorbdim
  integer :: aerr,ierr
  integer :: ia,ja,ma
  
!------------- READ IN OPERATOR FILE AND DECOUPLE

if(iproc == 0)then
  success = .false.
  do while(.not.success)  
     write(6,*)' Enter name of .spa file '
     read(5,'(a)')appfile
     ilastap = index(appfile,' ') -1
     open(unit = 23,file=appfile(1:ilastap)//'.spa',status = 'old',err=301)
     success = .true.
     cycle
301  continue
     print*,'File ',appfile(1:ilastap),'.spa does not exist '
   enddo
!------- READ IN HEADER???
   read(23,'(a)')title
   print*,title

end if

!.......... SET UP OPERATOR ARRAYS.................

numorbdim=max(numorb_i(1),numorb_i(2) )
numorbdim=max(numorbdim,numorb_f(1))
numorbdim=max(numorbdim,numorb_f(2))

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
     read(23,*)i1,n,l,xj
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

   read(23,*)it  ! species
   
   if(dnp(it) == 0)then
	   print*,' '
	   print*,' ERROR -- need to change species in the .spa file '
	   print*,' '
	   stop
   end if
   if(it==1)then
	   print*,' '
	   print*,' ** ADDING/REMOVING A PROTON ** '
	   print*,' '
   else
	   print*,' '
	   print*,' ** ADDING/REMOVING A NEUTRON ** '
	   print*,' '
   end if
	   
end if

call BMPI_BCAST(it,1,0,icomm,ierr)
  if(it==1)then
	  
     allocate ( psamp( numorbdim ), stat=aerr)
     if(aerr /= 0) call memerror("readinspamp -- p")
     allocate(pspopme(max(nrhsps_f(1),nrhsps_i(1))))
	 psamp = 0.0
     pspopme = 0.0
  end if
  if(it==2)then
	  
     allocate ( nsamp( numorbdim ), stat=aerr)
     if(aerr /= 0) call memerror("readinspamp -- n")
	 allocate(nspopme(max(nrhsps_i(2),nrhsps_f(2))))
	 nsamp = 0.0
     nspopme = 0.0
  end if   		

!------------ READ IN REDUCED MATRIX ELEMENTS ----------------
if(iproc==0)then
  do i = 1, 1000

        read(23,*,end =348)a,xxx

		if(it==1)then
			psamp(a)=xxx
		else
			nsamp(a)=xxx
			
		end if

  enddo

348      continue

  close(unit=23)
end if
!--------------- NOW BROADCAST -------------------------------

do a = 1,numorb_i(1)
            if(it==1)call BMPI_BCAST(psamp(a),1,0,icomm,ierr)
            if(it==2)call BMPI_BCAST(nsamp(a),1,0,icomm,ierr)

end do
call MPI_Barrier(icomm,ierr)

do ia = 1, nrhsps_f(it)
   a = rhspsqn_f(it,ia)%orb
   ja= rhspsqn_f(it,ia)%j
   ma= rhspsqn_f(it,ia)%m
  
  if(ma /= dnp(it)*(Jz_f-Jz_i))cycle
  if(it==1)pspopme(ia)=psamp(a)
  if(it==2)nspopme(ia)=nsamp(a)
  
end do

!print*,' testing p ',pspopme
!print*,' testing n ',nspopme


  return
  end subroutine readin1spamp
!============================================================
!
! base routine to apply proton 1-body addition/removal
! one issue to watch is truncation
!
subroutine applysphop_p
	use jumps1body
	use sectors
	use basis_map
	use localvectors
	use w_info
	use spstate ! for testing only
	use sporbit,only:sameweightscheme
	use hops
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
	integer :: xop ! = ac+bd
	integer(4) :: Wpf,Wnf
	
	do ipsj = 1,nsectorjumps(1)  ! loop over proton sector jumps
		is = psecthop(ipsj)%isector
		fs = psecthop(ipsj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
		ncsi = 	xsd_i(1)%sector(is)%ncsectors				
		iop1body = 			psecthop(ipsj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',ipsj,iop1body
			cycle
		end if
        wpf = psecthop(ipsj)%wXf
										
		do ipjmp = 1,psecthop(ipsj)%njumps
			isdp = psecthop(ipsj)%isd(ipjmp)
			fsdp = psecthop(ipsj)%fsd(ipjmp)
			if(fsdp==-1)cycle
			ipstate = pstart_i(isdp)
			fpstate = pstart_f(fsdp)
			pphase= psecthop(ipsj)%phase(ipjmp)  ! fetch phase from applying the operator
			iopp   = psecthop(ipsj)%iXop(ipjmp)  ! find index of one-body operator
			bd = phopop(iop1body)%Xop(iopp,1)
			ac = phopop(iop1body)%Xop(iopp,2)
			xop = ac+bd  ! oneof these is zero but to reduce switching I use both
!			print*,'check 1 ', ac,bd,xop,pspopme(xop)
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
					xme = pspopme(xop)*pphase*nphase
					vec2(fpstate+fnstate)= vec2(fpstate+fnstate)+xme*vec1(ipstate+instate)

				end do
			end do
			
		end do ! ipjmp
	end do ! ipsj
	
	return
end subroutine applysphop_p
!==================================================================================
!
! base routine to compute neutron 1-body addition/removal
! one issue to watch is truncation
!
subroutine applysphop_n
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
	integer :: ac,bd ,xop   ! proton creation/annihilation operators
	integer(4) :: Wpf,Wnf
	
	
	do insj = 1,nsectorjumps(2)  ! loop over neutron sector jumps
		is = nsecthop(insj)%isector
		fs = nsecthop(insj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
		ncsi = 	xsd_i(2)%sector(is)%ncsectors				
		iop1body = 			nsecthop(insj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',insj,iop1body
			cycle

		end if
        wnf = nsecthop(insj)%wXf
										
		do injmp = 1,nsecthop(insj)%njumps
			isdn = nsecthop(insj)%isd(injmp)
			fsdn = nsecthop(insj)%fsd(injmp)
			if(fsdn==-1)cycle
			instate = nstart_i(isdn)
			fnstate = nstart_f(fsdn)
			nphase= nsecthop(insj)%phase(injmp)  ! fetch phase from applying the operator
			iopn   = nsecthop(insj)%iXop(injmp)  ! find index of one-body operator
			bd = nhopop(iop1body)%Xop(iopn,1)
			ac = nhopop(iop1body)%Xop(iopn,2)
			xop= ac+bd
!............ LOOP OVER CONJUGATE NEUTRON SDS....................			
			
			do ic = 1,ncsi   ! loop over conjugate proton sectors
				ics = xsd_i(2)%sector(is)%csector(ic)  ! find index of conjugate proton sector
				                                       ! HOWEVER must be careful when truncating--
				                                       ! not all FINAL neutron SD will be allowed
				                                       ! if the sum of final Ws goes too high.
				if(sameweightscheme)then              ! can check this issue directly
					wpf = xsd_i(1)%sector(ics)%wx
					
					if(wpf + wnf > maxWtot_f)cycle	
					
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
					xme = nspopme(xop)*pphase*nphase
					vec2(fpstate+fnstate)= vec2(fpstate+fnstate)+xme*vec1(ipstate+instate)
!					print*,' hey ', ipstate+instate,fpstate+fnstate,xop,xme
					
				end do
				
			end do
			
		end do ! ipjmp
		

	end do ! ipsj
!	print*,vec2
	return
end subroutine applysphop_n
!==================================================================================	
!
! routines to determine the available one-particle hops
!
subroutine survey_sphop_p
	use jumps1body
	use sectors
!	use basis_map
!	use localvectors
!	use w_info
!	use spstate ! for testing only
!	use sporbit,only:sameweightscheme
	use hops
	implicit none
	
	integer :: ipsj    ! label of proton sector jump
!	integer :: is,fs,ics,fcs      ! label of initial/final proton sector, conjugate neutron sector(s)
!	integer :: ic
!	integer :: ncsi,ncsf     ! # of initial, final neutron sectors
	integer :: ipjmp   ! label of proton jumps
	
	integer(kind=basis_prec) :: isdp, fsdp   ! initial, final proton sds
!	integer(kind=basis_prec) :: isdn, fsdn   ! initial, final neutron sds
	                                           ! in principle insdi=insdf, but may be different
											   ! if bases are different
!	integer :: pphase,nphase
	integer :: iopp       ! label of one-body proton operator
!    integer(kind=basis_prec) :: ipstate,fpstate,instate,fnstate
!	real(8) :: xme
	integer :: iop1body
	integer :: ac,bd    ! proton creation/annihilation operators
	integer :: xop ! = ac+bd
	integer(4) :: Wpf,Wnf
		
	xsp_used=.false.
	do ipsj = 1,nsectorjumps(1)  ! loop over proton sector jumps
!		is = psecthop(ipsj)%isector
!		fs = psecthop(ipsj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
!		ncsi = 	xsd_i(1)%sector(is)%ncsectors				
		iop1body = 			psecthop(ipsj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',ipsj,iop1body
			cycle
		end if
!        wpf = psecthop(ipsj)%wXf
										
		do ipjmp = 1,psecthop(ipsj)%njumps
!			isdp = psecthop(ipsj)%isd(ipjmp)
			fsdp = psecthop(ipsj)%fsd(ipjmp)
			if(fsdp==-1)cycle
!			ipstate = pstart_i(isdp)
!			fpstate = pstart_f(fsdp)
!			pphase= psecthop(ipsj)%phase(ipjmp)  ! fetch phase from applying the operator
			iopp   = psecthop(ipsj)%iXop(ipjmp)  ! find index of one-body operator
			bd = phopop(iop1body)%Xop(iopp,1)
			ac = phopop(iop1body)%Xop(iopp,2)
			xop = ac+bd  ! oneof these is zero but to reduce switching I use both
			xsp_used(xop)=.true.
			
		end do ! ipjmp
	end do ! ipsj
		
	return
	
end subroutine survey_sphop_p
!==================================================================================	
subroutine survey_sphop_n_new
	use jumps1body
	use sectors
!	use basis_map
!	use localvectors
	use w_info
!	use spstate ! for testing only
	use sporbit,only:sameweightscheme
	use hops
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
!	integer :: pphase,nphase
	integer :: iopn       ! label of one-body proton operator
!    integer(kind=basis_prec) :: ipstate,fpstate,instate,fnstate
!	real(8) :: xme
	integer :: iop1body
	integer :: ac,bd    ! proton creation/annihilation operators
	integer :: xop ! = ac+bd
	integer(4) :: Wpf,Wnf
		
	xsp_used=.false.
	
	do insj = 1,nsectorjumps(2)  ! loop over neutron sector jumps
!		is = nsecthop(insj)%isector
!		fs = nsecthop(insj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
!		ncsi = 	xsd_i(2)%sector(is)%ncsectors				
		iop1body = 			nsecthop(insj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',insj,iop1body
			cycle

		end if
 !       wnf = nsecthop(insj)%wXf
										
		do injmp = 1,nsecthop(insj)%njumps
!			isdn = nsecthop(insj)%isd(injmp)
!			fsdn = nsecthop(insj)%fsd(injmp)
			if(fsdn==-1)cycle
!			instate = nstart_i(isdn)
!			fnstate = nstart_f(fsdn)
!			nphase= nsecthop(insj)%phase(injmp)  ! fetch phase from applying the operator
			iopn   = nsecthop(insj)%iXop(injmp)  ! find index of one-body operator
			bd = nhopop(iop1body)%Xop(iopn,1)
			ac = nhopop(iop1body)%Xop(iopn,2)
			xop= ac+bd
!
!............ LOOP OVER CONJUGATE NEUTRON SDS....................			
			
			do ic = 1,ncsi   ! loop over conjugate proton sectors
				ics = xsd_i(2)%sector(is)%csector(ic)  ! find index of conjugate proton sector
				                                       ! HOWEVER must be careful when truncating--
				                                       ! not all FINAL neutron SD will be allowed
				                                       ! if the sum of final Ws goes too high.
				if(sameweightscheme)then              ! can check this issue directly
					wpf = xsd_i(1)%sector(ics)%wx
					
					if(wpf + wnf > maxWtot_f)cycle	
					
				end if									   
				do isdp = xsd_i(1)%sector(ics)%xsdstart,xsd_i(1)%sector(ics)%xsdend   !loop over proton SDs 
					fsdp= psdmap(isdp)  ! need to map this initial, unchanged proton SD
					                      ! to a final proton SD in the final basis
					if(fsdp ==-1)cycle
					
					if(.not.sameweightscheme)then
							wpf=psdmap_w(isdp)		
							if(wpf + wnf > maxWtot_f)cycle									   					
									  
					end if			
			        xsp_used(xop)=.true.
					go to 1        ! if true, can 
				end do
				
			end do
1           continue
			
		end do
	end do
!............ LOOP OVER CONJUGATE NEUTRON SDS....................			
	
	return
	
	
end subroutine survey_sphop_n_new
!==================================================================================	
!==================================================================================	
subroutine survey_sphop_n
	use jumps1body
	use sectors
!	use basis_map
!	use localvectors
!	use w_info
!	use spstate ! for testing only
!	use sporbit,only:sameweightscheme
	use hops
	implicit none
	
	integer :: insj    ! label of proton sector jump
!	integer :: is,fs,ics,fcs      ! label of initial/final proton sector, conjugate neutron sector(s)
!	integer :: ic
!	integer :: ncsi,ncsf     ! # of initial, final neutron sectors
	integer :: injmp   ! label of proton jumps
	
	integer(kind=basis_prec) :: isdp, fsdn   ! initial, final proton sds
!	integer(kind=basis_prec) :: isdn, fsdn   ! initial, final neutron sds
	                                           ! in principle insdi=insdf, but may be different
											   ! if bases are different
!	integer :: pphase,nphase
	integer :: iopn       ! label of one-body proton operator
!    integer(kind=basis_prec) :: ipstate,fpstate,instate,fnstate
!	real(8) :: xme
	integer :: iop1body
	integer :: ac,bd    ! proton creation/annihilation operators
	integer :: xop ! = ac+bd
	integer(4) :: Wpf,Wnf
		
	xsp_used=.false.
	
	do insj = 1,nsectorjumps(2)  ! loop over neutron sector jumps
!		is = nsecthop(insj)%isector
!		fs = nsecthop(insj)%fsector  ! in principle, is = fs
                                        ! note however value may be different
										! if bases are different
!		ncsi = 	xsd_i(2)%sector(is)%ncsectors				
		iop1body = 			nsecthop(insj)%iop
		if(iop1body==0)then
			print*,'  Empty sector ',insj,iop1body
			cycle

		end if
 !       wnf = nsecthop(insj)%wXf
										
		do injmp = 1,nsecthop(insj)%njumps
!			isdn = nsecthop(insj)%isd(injmp)
!			fsdn = nsecthop(insj)%fsd(injmp)
			if(fsdn==-1)cycle
!			instate = nstart_i(isdn)
!			fnstate = nstart_f(fsdn)
!			nphase= nsecthop(insj)%phase(injmp)  ! fetch phase from applying the operator
			iopn   = nsecthop(insj)%iXop(injmp)  ! find index of one-body operator
			bd = nhopop(iop1body)%Xop(iopn,1)
			ac = nhopop(iop1body)%Xop(iopn,2)
			xop= ac+bd
			xsp_used(xop)=.true.
		end do
	end do
!............ LOOP OVER CONJUGATE NEUTRON SDS....................			
	
	return
	
	
end subroutine survey_sphop_n
!==================================================================================	
! this just updates the array p/nspopme so that a single s.p. is used
!
subroutine update_spop(it,f,iorb)
	use spstate
	implicit none
	integer :: it
	integer :: f,ff ! these count up how many orbitals are available
	integer :: j  ! this is the actual index of the orbital
	integer :: iorb  ! final orbit

	if(it==1)then
		pspopme =0.0
		ff = 0
		do j = 1,size(xsp_used)
			if(xsp_used(j))ff= ff+1
			if(ff==f)exit
		end do
		if(ff==f)pspopme(j)=1.0
	else
		nspopme =0.0
		ff = 0
		do j = 1,size(xsp_used)
			if(xsp_used(j))ff= ff+1
			if(ff==f)exit
		end do
		if(ff==f)nspopme(j)=1.0		
		
	end if
	
	iorb = rhspsqn_f(it,j)%orb
	
	if(iorb == 0)then
!		print*, xsp_used
		print*,f,' zero orbit ? ',j,nspopme(j)
		print*,rhspsqn_f(it,j)%orb
	else
		print*,' orbital occupied ',j,iorb,f,rhspsqn_f(it,j)%j, rhspsqn_f(it,j)%m,rhspsqn_f(it,j)%nr
	end if
	
!	print*,' update ',f,iorb,rhspsqn_f(it,j)%orb

	
	return
end subroutine update_spop
!==================================================================================	
end module spfac
	
	
	
