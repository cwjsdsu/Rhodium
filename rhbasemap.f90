!
! data and routines to map basis states
!
! only works when same M, par

module basis_map
	use basis
	use precisions
	use spstate
	implicit none
	
	integer(kind=basis_prec), allocatable,target :: psdmap(:),nsdmap(:)
	integer(4), allocatable,target :: psdmap_w(:),nsdmap_w(:)

	integer(kind=2), allocatable, target :: psdmapphase(:),nsdmapphase(:)
	
	logical :: needphase
	
contains

!....... ROUTINE TO READ IN SLATER DETERMINANT IN THE INITIAL BASIS
!        AND CREATE ITS EQUIVALENT IN THE FINAL BASIS
!	CALLED BY:
!    basismapper
!
!  SUBROUTINES CALLED:
!   convert_bitrepsd2occ
!   sortoccwphase
!   convert_occ2bitrepsd
!
	subroutine sdmapper(it,sdin,occin,sdout,success,occout,phase)
		use system_parameters
		implicit none
		integer :: it
		integer(8) :: sdin(Nword_i(it)),sdout(Nword_f(it))
		integer   :: occin(np_i(it)),occout(np_f(it))
		logical   :: success
		integer :: ip
		integer :: isp,fsp
		integer(2) :: phase
		
		if(np_i(it)/= np_f(it))then
			print*,' mismatch in # of particles (sdmapper)',it,np_i(it),np_f(it)
!			stop
		end if
			
		call convert_bitrepsd2occ(it,'INI',np_i(it),sdin,occin)
		
		occout=0
		sdout = 0
		success = .false.
		do ip = 1,np_i(it)
			isp = occin(ip)
			fsp = rhspsqn_i(it,isp)%ifmap
			if(fsp==0)return  ! could not map
			occout(ip)=fsp 
		end do
!......... HERE IS WHERE PHASES WOULD COME IN	

        if(needphase)then
			call sortoccwphase(np_f(it),occout,phase)
		else
			phase=1
			
		end if	

        success = .true.
        call convert_occ2bitrepsd(it,'FIN',np_f(it),occout,sdout)
		
		return
	end subroutine sdmapper
!=======================================================================
!....... ROUTINE TO READ IN SLATER DETERMINANT IN THE INITIAL BASIS
!        AND CREATE ITS EQUIVALENT IN THE FINAL BASIS
! ..... AUGMENTED TO ALLOW TO ADD/SUBTRACT PARTICLES
!
!  ... NEED TO FIX IF DESTROY MORE PARTICLES THAN CREATED... TO BE DONE!
!
!   CALLED BY:
!   count_fill_1body_sectorjumps
!	
!  SUBROUTINES CALLED:
!   convert_bitrepsd2occ
!   sortoccwphase
!   convert_occ2bitrepsd
!
	subroutine sdmapperplus(it,sdin,occin,ndestroy,destroy,ncreate,create,sdout,success,occout,phase)
		use system_parameters
		implicit none
		integer :: it
		integer(8) :: sdin(Nword_i(it)),sdout(Nword_f(it))
		integer   :: occin(np_i(it)),occout(np_f(it))
		integer(4),intent(in) :: ndestroy,ncreate
		integer(4),intent(in) :: destroy(ndestroy),create(ncreate)
		integer :: idestroy,icreate
		logical   :: success,found
		integer :: ip,jp
		integer :: isp,fsp
		integer(2) :: phase
		
		
		if(np_f(it)/= np_i(it) + ncreate-ndestroy)then
			print*,' mismatch in # of particles (sdmapperplus)'
			print*,it,np_i(it),np_f(it),ncreate,ndestroy
			stop
		end if
		call convert_bitrepsd2occ(it,'INI',np_i(it),sdin,occin)
		
		occout=0
		sdout = 0
		success = .false.
		
!......... REMOVE DESTROYED OPERATORS......replace them with 9999
   if(ndestroy> 0 )then
        do idestroy = 1,ndestroy
			isp = destroy(idestroy)
			found=.false.
            do ip = 1,np_i(it)
				if(isp==occin(ip))then
					occin(ip)=9999
				    found=.true.
				    exit
				end if
		    end do
			if(.not.found)return
	    end do
	end if
	
		icreate = 0
		
!.......... IF JUST DESTROYING, NEED TO SORT FIRST.........		
		
        if(ndestroy > ncreate)then
			if (needphase)then
			   call sortoccwphase(np_i(it),occin,phase)
		    else
			    phase=1
			end if
			
		end if
		
		jp =0
		
		do ip = 1,np_i(it) !min(np_i(it),np_f(it))  !don't go beyond either occout or occin
			isp = occin(ip)
			if(isp==9999)then   !this state was destroyed, replace with a 'created' one
				icreate = icreate+1
				if(icreate <= ncreate)then  !  keeps from creating too many
					jp = jp+1
				     occout(jp)=create(icreate)
			     end if
			 else
			   fsp = rhspsqn_i(it,isp)%ifmap
			   if(fsp==0)return  ! could not map
			   jp = jp+1
			   occout(jp)=fsp 			
		    end if
		end do
		if(ndestroy > ncreate)then
			success=.true.
	        call convert_occ2bitrepsd(it,'FIN',np_f(it),occout,sdout)
			
			return
		end if
		
!.......... ADD ANY ADDITIONAL CREATION not yet covered
        do while(icreate < ncreate)
		   icreate = icreate+1
		   jp = jp+1
		   occout(jp)=create(icreate)	
		end do ! while

!......... HERE IS WHERE PHASES WOULD COME IN	

        if(needphase)then
			call sortoccwphase(np_f(it),occout,phase)
		else
			phase=1
			
		end if	

        success = .true.
		
        call convert_occ2bitrepsd(it,'FIN',np_f(it),occout,sdout)
		
		return
	end subroutine sdmapperplus
!================================================================
    subroutine phasechecker
		use nodeinfo ! mpi
		implicit none
		integer :: isps
		integer :: it
		integer :: fsps,fsps1
		
		needphase = .false.
		do it = 1,2
		  do isps = 1,nrhsps_i(it)-1
			  fsps = rhspsqn_i(it,isps)%ifmap
			  fsps1 = rhspsqn_i(it,isps+1)%ifmap

			  if(fsps> fsps1)then
				  needphase=.true.
				  if(iproc==0)print*,' Need to include phases in basis mapping '
				  return
			  end if
			  
		
		  end do
	    end do
		if(iproc==0)print*,' No need for phases in basis mapping '
		return
	end subroutine phasechecker	
!=======================================================================
	
	subroutine basismapper(it)
		use precisions
		use system_parameters
		use sectors
		use nodeinfo ! mpi
		use sporbit,only:sameweightscheme
		
		implicit none
		integer :: it
		integer(8),pointer :: xsdmap(:)
		integer(2),pointer :: xsdmapphase(:)
		integer(kind=basis_prec) :: isd,fsd
		integer,allocatable :: occin(:),occout(:)
		integer,pointer:: xsdmap_w(:)
		integer(kind=8),pointer:: xsdlist_i(:,:)
		integer(kind=8),allocatable :: sdout(:)
		logical :: success,foundit
		integer :: is,fs
		integer :: Mxi,parxi,Wxf
		integer :: notfound
		integer(2) :: phase
		
		if(np_i(it)/= np_f(it))then
			print*,' Note mismatch in # of particles (basismapper)'
!			stop
		end if
		
		if(jz_i /= jz_f)then
			print*,' Note mismatch in M  '
!			stop
		end if
		
		if(iparity_i /= iparity_f)then
			print*,' Note mismatch in parity '
!			stop
		end if
		
		if(it==1)then
			allocate(psdmap(nxsd_i(1)))
			xsdmap => psdmap
			allocate(psdmapphase(nxsd_i(1)))
			xsdmapphase => psdmapphase
			if(.not.sameweightscheme)then
				allocate(psdmap_w(nxsd_i(1)))
				xsdmap_w => psdmap_w
			end if
		else
			allocate(nsdmap(nxsd_i(2)))
			xsdmap => nsdmap
			allocate(nsdmapphase(nxsd_i(2)))
			xsdmapphase => nsdmapphase
			
			if(.not.sameweightscheme)then
				allocate(nsdmap_w(nxsd_i(2)))
				xsdmap_w => nsdmap_w
			end if
		end if		
		xsdmap = -1
		xsdmapphase=1
		
		if(np_i(it)==0)then
			xsdmap(1)=1
			return
		end if
		
		allocate(occin(np_i(it)),occout(np_f(it)))
		
		if(it==1)then
			xsdlist_i=> psdlist_i
		else
			xsdlist_i=> nsdlist_i
		end if
		allocate(sdout(Nword_f(it)))
		
!........... TO SPEED UP SEARCHING, DIVVY UP BY SECTORS AND MATCH Jz,PARITY
        notfound = 0
        do is = 1,nsectors_i(it)  ! LOOP OVER INITIAL SECTORS
			Mxi = xsd_i(it)%sector(is)%jzx
!			print*,it,Mxi, xsd_i(it)%sector(is)%xsdstart,xsd_i(it)%sector(is)%xsdend
			parxi=xsd_i(it)%sector(is)%parx
			do isd = xsd_i(it)%sector(is)%xsdstart,xsd_i(it)%sector(is)%xsdend	  ! LOOP OVER SDs in that initial sector		
			   call sdmapper(it,xsdlist_i(isd,:),occin,sdout,success,occout,phase) ! TRANSLATE TO FINAL SPACE			   
			   if(.not.success)then
				   xsdmap(isd)=-1
				   notfound = notfound+1
				   !print*,Mxi,' Jz not found (1)'
			   else        ! NOW SEARCH
				   Wxf = w_of_occ_f(it,np_f(it),occout)
				   call searcher(it,sdout,Mxi,parxi,Wxf,foundit,fsd)
				   if(foundit)then
					   xsdmap(isd)=fsd
					   xsdmapphase(isd)=phase
					   if(.not.sameweightscheme)xsdmap_w(isd)=Wxf
				   else
					   xsdmap(isd)=-1
					   notfound = notfound+1
!					   print*,Mxi,' Jz not found (2)',isd
					   
				   end if
				   
			   end if
			   
		   end do
		end do
		if(iproc==0) then
			if(notfound > 0)print*,notfound,' SDs not mapped (probably due to conflicting quantum numbers)'
		endif
		return
	end subroutine basismapper
!=======================================================================
!
! routine to (unordered) search for suitable final basis states
! restrict to search within Jz/M, parity, and W of final Slater determinant
! note, because more general, much slower than BIGSTICK routines!

!------------------------------
subroutine searcher(it,sdout,Mxf,parxf,Wxf,foundit,fsd)
	use system_parameters
	use sectors
	implicit none
	integer :: it
	integer(8) :: sdout(Nword_f(it))
	integer :: Mxf,parxf,Wxf
	logical :: foundit,success
	integer(kind=basis_prec) :: fsd,myfsd,targetfsd
	integer :: iword
	integer(kind=8),pointer:: xsdlist_f(:,:)
	integer :: fs
	integer :: omp_get_num_threads
	integer :: nthreads
	
!$OMP PARALLEL
  nthreads=omp_get_num_threads() ! to see if one uses threaded search or not
!$OMP END PARALLEL	
	fsd =-1
	foundit=.false.
	if(it==1)then
		xsdlist_f=> psdlist_f
	else
		xsdlist_f=> nsdlist_f
	end if	
	do fs = 1,nsectors_f(it)
		if(Mxf==xsd_f(it)%sector(fs)%jzx .and. parxf == xsd_f(it)%sector(fs)%parx .and. Wxf==xsd_f(it)%sector(fs)%Wx )then
			if(nthreads==1 .or. Nword_f(it) < 2)then
						
			do fsd= xsd_f(it)%sector(fs)%xsdstart,xsd_f(it)%sector(fs)%xsdend
				foundit = .true.
				
				do iword = 1,Nword_f(it)
				    if(sdout(iword)/= xsdlist_f(fsd,iword))then
					    foundit=.false.
					    exit
				    end if
			    end do
				if(foundit)exit

			end do
			if(foundit)then
				targetfsd = fsd
				return
			else
				targetfsd = 0
			end if

			else   ! OPENMP THREADED SEARCH
				myfsd = 0
!$OMP PARALLEL DO PRIVATE(iword,foundit)  SHARED(myfsd)  SCHEDULE(static)	
				do fsd= xsd_f(it)%sector(fs)%xsdstart,xsd_f(it)%sector(fs)%xsdend
					foundit = .true.
					
					    do iword = 1,Nword_f(it)
					       if(sdout(iword)/= xsdlist_f(fsd,iword))then
						      foundit=.false.
						      exit
					       end if
				        end do
						if(foundit)then
							myfsd = fsd
						end if

				end do 			
!$OMP END PARALLEL	DO	
                fsd = myfsd
				if(myfsd /=0)foundit=.true.


           end if		
			if(foundit)return
			
		end if
	end do
	
	return
	
end subroutine searcher
!------------------------------
! subroutine to order an occupation, with a phase
!
! this assumes mostly in order, tries to minimizes sorting
! probably could be made faster, not sure how
!
! INPUT:
!   npx = number of particles and dimension of occ
!   occ(1:npx) = array of occupied single particle states, needs to be in order
! OUTPUT:
!   phase


subroutine sortoccwphase(npx,occ,phase)
	implicit none
	integer :: npx
	integer :: occ(npx)
	integer(2) :: phase
	integer :: i,j,k,m,ok
	logical :: mustsort
	
	phase = 1
    mustsort=.false.
	do i = 1,npx-1
		if(occ(i+1)< occ(i))then
			mustsort = .true.
			exit

		end if
	end do

	if(.not.mustsort)return
		
	do i =1,npx-1
		k=i
		ok = occ(i)
		do j =i+1,npx
			if(occ(j)< ok)then
				ok = occ(j)
				k = j
			end if
		end do
		if(k >i)then  ! swap
			occ(k)=occ(i)
			occ(i)=ok
			phase = -phase
			
		end if

	end do
	
	
	return
	
end subroutine sortoccwphase
!------------------------------
!
! function to extract W value of final state

function w_of_occ_f(it,npx,occ)
	implicit none
	integer :: it
	integer :: npx
	integer :: occ(npx)
	integer :: w_of_occ_f
	integer :: ip,isps
	
	w_of_occ_f = 0
	if(npx==0)return
	
	do ip = 1,npx
		isps = occ(ip)
		w_of_occ_f = w_of_occ_f + rhspsqn_f(it,isps)%w
	end do 
	
	return
end function w_of_occ_f
!...... ROUTINES FOR HANDLING PARITY....................
!this function may not be needed
function parity_of_occ(it,which,npx,occ)
	implicit none
	integer :: it
	character :: which
	integer :: npx
	integer :: occ(npx)
	type (spst),pointer :: rhspsqn_x(:,:)
	integer :: ip,isps
	
	integer(2) :: parity_of_occ
	
	select case(which)
	
	case('i')
	
	rhspsqn_x => rhspsqn_i
	
	case('f')
	rhspsqn_x => rhspsqn_f
	
	case default
	
	print*,' problem in parity_of_occ '
	stop
	
    end select
	
	parity_of_occ = 1
	if(npx==0)return
	
	do ip = 1,npx
		isps=occ(ip)
		parity_of_occ = parity_of_occ * (-1)**rhspsqn_x(it,isps)%l
	end do
		
	return
end function parity_of_occ

end module basis_map
	
	
	