module parr
	use basis
	use basis_map
	!... ADDED 0.8.9: parity information
	  integer(2), allocatable, target :: pparlist_i(:),nparlist_i(:),  pparlist_f(:),nparlist_f(:)
contains
	!call this routine to set up parity information
	subroutine parity_boss
		implicit none
	
		call parity_search(1,'i')
		call parity_search(1,'f')
		call parity_search(2,'i')
		call parity_search(2,'f')
		return
	end subroutine parity_boss
	!called by subroutine parity_boss
	    subroutine parity_search(it,which)
			use sectors
			use system_parameters
			implicit none
			integer :: it
			character :: which
			integer :: nXsdx
			integer(2),pointer :: xparlist_x(:)
			integer :: isector, nsectorx
			type (mastersect),pointer :: xsdx 
			integer :: parx
			integer(kind=8) :: isd
			integer :: npx
			integer, allocatable :: occx(:)
			integer(kind=8),pointer:: xsdlist_x(:,:)
			

			select case(which)
		
			case('i')
			nXsdx = nXsd_i(it)
			nsectorx = nsectors_i(it)
			xsdx => xsd_i(it)
			npx = np_i(it)
			if(it==1)then
				if(.not.allocated(pparlist_i))allocate(pparlist_i(nXsd_i(it)))
				xparlist_x=> pparlist_i
				xsdlist_x=> psdlist_i
				
			else
				if(.not.allocated(nparlist_i))allocate(nparlist_i(nXsd_i(it)))
				xparlist_x=> nparlist_i
				xsdlist_x=> nsdlist_i
				
			end if
	
	
			case('f')
			nXsdx = nXsd_f(it)
			nsectorx = nsectors_f(it)
			xsdx => xsd_f(it)	
			npx = np_f(it)
			
			if(it==1)then
				if(.not.allocated(pparlist_f))allocate(pparlist_f(nXsd_f(it)))
				xparlist_x=> pparlist_f
				xsdlist_x=> psdlist_f
				
			else
				if(.not.allocated(nparlist_f))allocate(nparlist_f(nXsd_f(it)))
				xparlist_x=> nparlist_f
				xsdlist_x=> nsdlist_f
				
			
			end if
	
			case default	
			print*,' problem in parity_of_occ '
			stop

	     	end select
		
			xparlist_x = 0
			if(npx ==0 )then
				xparlist_x = 1
				return
			end if
			allocate(occx(npx))
			do isd = 1,nXsdx
				if(which == 'i')then
				  call convert_bitrepsd2occ(it,'INI',npx,xsdlist_x(isd,:),occx)
			    else
					
  				  call convert_bitrepsd2occ(it,'FIN',npx,xsdlist_x(isd,:),occx)
					
				end if
				parx = parity_of_occ(it,which,npx,occx)
				xparlist_x(isd)=parx
			end do
		
		
		
			!search over the sectors
		
!			do isector = 1,nsectorx
!				parx = xsdx%sector(isector)%parx
!				print*,isector,parx,' sector parity'
!				do isd= xsdx%sector(isector)%xsdstart,xsdx%sector(isector)%xsdend
!					xparlist_x(isd)=parx
!				end do			
!			end do
            if(allocated(occx))deallocate(occx)
			return
		end subroutine parity_search
		! routine to find the parity of a given vector
		!  which = i or f meaning vectors 1, or 2
		!  returns mypar
		subroutine assign_parity(which,mypar)
			use localvectors
			use precisions
			use basis
			use sectors
			implicit none
			character :: which
			integer :: mypar
			real(kind=lanc_prec), pointer :: vecx(:)
			integer(kind=basis_prec) :: dimbasis_x,i,ii
			real :: tol
			integer,pointer :: nsectors_x(:)
			type (mastersect),pointer :: xsdx(:)
			integer :: isp,isn,isdp,isdn,csi
			integer(kind=basis_prec),pointer :: pstart_x(:),nstart_x(:)
			integer(kind=basis_prec) :: ibasis,ipstate,instate
			integer(2),pointer :: pparlist_x(:),nparlist_x(:)
			
			
			select case (which)
			
			case('i')
			
			vecx => vec1
			nsectors_x => nsectors_i
			xsdx => xsd_i
			pstart_x => pstart_i
			nstart_x => nstart_i

			dimbasis_x = dimbasis_i
			pparlist_x => pparlist_i
			nparlist_x => nparlist_i
			
			case('f')
			
			vecx => vec2
			nsectors_x => nsectors_f
			xsdx => xsd_f
			pstart_x => pstart_f
			nstart_x => nstart_f
			dimbasis_x = dimbasis_f
			pparlist_x => pparlist_f
			nparlist_x => nparlist_f
			
			case default
			
			print*,' that choice not allowed ',which
			stop
			
   		    end select
			
			tol = 1.d0/sqrt(real(dimbasis_x+1,8))
			ii = 0
			do i = 1,dimbasis_x ! search for a nontrivial component
				if(abs(vecx(i))> tol)then
					ii = i
					exit
				end if
			end do
			if(ii == 0)then
				print*,' problem in assign_parity ', which,dimbasis_x,tol
				stop
			end if	
			do isp = 1,nsectors_x(1)  ! loop over proton sectors
				do isdp = xsdx(1)%sector(isp)%xsdstart,xsdx(1)%sector(isp)%xsdend	  ! LOOP OVER SDs in that initial sector
					
					ipstate = pstart_x(isdp)    ! find index of initial proton state
					do csi= 1,xsdx(1)%sector(isp)%ncsectors     ! loop over conjugate neutron sectors
						isn = xsdx(1)%sector(isp)%csector(csi)   ! find the index of the conjugate neutron sector
						do isdn = xsdx(2)%sector(isn)%xsdstart,xsdx(2)%sector(isn)%xsdend	
							instate = nstart_x(isdn)
							ibasis = ipstate+instate
!							write(99,*)ibasis
							if(ibasis==ii)then
								mypar = pparlist_x(isdp)*nparlist_x(isdn)
						!		print*,ii,isdp,isdn,pparlist_x(isdp),nparlist_x(isdn)
								return
							end if
						end do
					end do
				end do
			end do
			
			print*,' never found parity '
			print*,ii,dimbasis_x,vecx(ii)
!			print*,pparlist_x
!			print*,nparlist_x
					
			stop 
		end subroutine assign_parity
	
	
end module parr