!
!  in RHODIUM, jumps are organized into OPTEAMS. OPTEAMS are similar to bundles and sectorjumps in BIGSTICK.
! 
!  types of teams:
!     dn1b      -- one-body density, non-charge-changing
!     cr1p      -- create one particle
!     an1h      -- anihilate to get one hole
!     d1cc      -- one-body charge-changing
!
!
module jumpteam_def
	use precisions
	implicit none
	
	type team
		character(6) :: teamtype    ! type of team
		logical :: skipme           ! if skipme, then cannot reach
		integer(kind = 4) :: ipsector,fpsector,insector,fnsector  ! initial/final proton/neutron sectors
		integer(kind=8) :: njumpsx(2)
		integer(kind=8), allocatable :: ipsd(:),fpsd(:),insd(:),fnsd(:)
		integer,allocatable :: pphase(:),nphase(:)
		integer :: dWx(2),dparx(2),dMx2(2)      !change in quantum numbers
		
		
	end type team
	
	integer(kind=4)          :: nopteams
	type (team), allocatable :: opteam(:)
	
	
	
end module jumpteam_def
!===========================================================================
!
! SUBROUTINES CALLED:
!  count_tot_sectors
!
module jumpteam_setup
	use jumpteam_def
	implicit none
	
contains
	
! MAJOR ROUTINE FOR ALLOCATING JUMP TEAMS
!	
	subroutine setup_jumpteams
		use sectors
		implicit none
		
		integer :: nsectorstot_i,nsectorstot_f
		integer :: isector,fsector,csector_i,csector_f,ics,fcs
		integer :: nteam
		integer :: Wpi,Wni,Wpf,Wnf,Jzpi,Jzni,Jzpf,Jznf,parpi,parni,parpf,parnf
		
		call count_tot_sectors('INI',nsectorstot_i)
		call count_tot_sectors('FIN',nsectorstot_f)
		
		nopteams = nsectorstot_i*nsectorstot_f
		
		allocate(opteam(nopteams))
		
		nteam = 0
		do isector = 1,nsectors_i(1)
			Wpi = xsd_i(1)%sector(isector)%wX
			Jzpi = xsd_i(1)%sector(isector)%jzx
			parpi = xsd_i(1)%sector(isector)%parx

			do csector_i = 1,xsd_i(1)%sector(isector)%ncsectors
				ics = xsd_i(1)%sector(isector)%csector(csector_i)
				Wni = xsd_i(2)%sector(ics)%wX
				Jzni = xsd_i(2)%sector(ics)%jzx
				parni = xsd_i(2)%sector(ics)%parx
				
				do fsector=1,nsectors_f(1)
					Wpf = xsd_f(1)%sector(fsector)%wX
					Jzpf = xsd_f(1)%sector(fsector)%jzx
					parpf = xsd_f(1)%sector(fsector)%parx
					do csector_f = 1,xsd_f(1)%sector(fsector)%ncsectors
						fcs = xsd_f(1)%sector(fsector)%csector(csector_f)
						Wnf= xsd_f(2)%sector(fcs)%wX
						Jznf = xsd_f(2)%sector(fcs)%jzx
						parnf = xsd_f(2)%sector(fcs)%parx
						nteam = nteam+1
						opteam(nteam)%ipsector = isector
						opteam(nteam)%fpsector = fsector
						opteam(nteam)%insector = ics
						opteam(nteam)%fnsector = fcs
						
						opteam(nteam)%dWx(1)= Wpf-Wpi
						opteam(nteam)%dMx2(1)= jzpf-jzpi
						opteam(nteam)%dparx(1)= parpi*parpf
						opteam(nteam)%dWx(2)= Wnf-Wni
						opteam(nteam)%dMx2(2)= jznf-jzni
						opteam(nteam)%dparx(2)= parni*parnf
						
					end do ! csector_f
				end do  ! fsector
			end do      ! csector_i
		end do         ! isector
		return
	
	end subroutine setup_jumpteams

!---------------------------------------------------
!
!	routine to count up all allowed combinations of proton-neutron sectors
!   can also assign these sectors to jumpteams
!
	subroutine count_tot_sectors(whichbasis,nsectorstot)
		use sectors
		implicit none
		character(3) :: whichbasis  ! either 'INI' or 'FIN' for initial or final
		integer      :: nsectorstot
		type (mastersect), pointer :: xsd(:)
		integer, pointer :: nsectors(:)	
		integer      :: isector	
			
			
		select case (whichbasis)
		
		case ('INI','ini')
		
		xsd => xsd_i
		nsectors => nsectors_i
		
		case ('FIN','fin')
		
		xsd => xsd_f
		nsectors => nsectors_f
		
		case default
		
		print*,' Should not be here ',whichbasis
		stop
		
		
  	    end select	
			
		nsectorstot = 0
		
		do isector = 1,nsectors(1)
			nsectorstot = nsectorstot + xsd(1)%sector(isector)%ncsectors
		end do			
		
		return
	end subroutine count_tot_sectors 	 
	
end module jumpteam_setup