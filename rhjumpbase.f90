!
!  RHODIUM 1-body jumps
!  started late July 2018 by CWJ @ SDSU
!
! NOTE Aug 2018: one potential problem is change in truncation wt W
! Nov 2018: If this is the case, must search through all Ws
! THiS PROBLEM WAS FIXED

module jumpbase
	use precisions

	implicit none
	
	
	type basejump
		integer(4) :: isector,fsector   ! initial, final sector
		integer(4) :: dJz,dParity    ! change in quantum numbers
		integer(4) :: dW !,dWmin,dWmax   ! account for change in conjugate W
		integer(4) :: wXf   ! final W
!		integer(4) :: nop1body  ! # of distinct quantum numbers
		integer(4) :: nop  ! # of distinct quantum numbers
!		integer(4) :: iop1body !(:)       ! which of Xops1body contain information on operators
		integer(4) :: iop !(:)       ! which of Xops contain information on operators
		integer(kind=basis_prec) :: njumps
		integer(kind=8), allocatable :: isd(:),fsd(:)
		integer(kind=4),allocatable :: iXop(:)
		integer(kind=2),allocatable :: phase(:)
		
	end type basejump
	
	integer nsectorjumps(2)

   
    type baseoperator
		integer(4) :: d2M,dW,dParity
		integer(4) :: nops
		integer,allocatable :: Xop(:,:)
    end type baseoperator
	
	
	
contains

!==================================================================
!	
!  generic routine to count up jumps between sectors
!
! MODIFIED MAY 2024: option to impose that initial and final neutron sectors overlap
!
!  INPUT: 
!   it = species; 1 = proton, 2 = neutron
!   fixedq : if true, then do not change quantum numbers
!   fill : if true, then start to fill in 
!
! CALLED BY:
!  boss1bjumps
!
	subroutine countsectorjumps(it,xsectjump,fixedJz,fixedpar,fill)
		use sectors
		use system_parameters
		use sporbit
		implicit none
		integer, intent(in) :: it  ! species, 1=proton,2=neutron
		logical,intent(in)  :: fixedpar,fixedJz,fill
		integer :: is,fs            ! initial, final 
		integer :: dJz,dParity,dW
		type (basejump) :: xsectjump(:)

		integer :: itc
		integer :: ics,cs
				
		itc = 3-it
				
		dJz = jz_f - jz_i
		if(fixedJz)print*,' change in Jz = ',dJz
		dParity = parmult(iparity_i,iparity_f)
		if(fixedpar)print*,' change in parity = ',3-2*dParity

		nsectorjumps(it) = 0
				
		do is = 1,nsectors_i(it)
			
			do fs = 1,nsectors_f(it)
				if( xsd_f(it)%sector(fs)%jzX -xsd_i(it)%sector(is)%jzX/= dJz .and. fixedJz)then
					cycle
				else
					dJz = xsd_f(it)%sector(fs)%jzX -xsd_i(it)%sector(is)%jzX
				end if
!				print*, 'testing parity ', xsd_i(it)%sector(is)%parX, xsd_f(it)%sector(fs)%parX
				if( parmult(xsd_i(it)%sector(is)%parX, xsd_f(it)%sector(fs)%parX)/= dParity .and. fixedpar)then
!					print*,'cycling '
					cycle
				else
					 dParity= parmult(xsd_i(it)%sector(is)%parX, xsd_f(it)%sector(fs)%parX)
				end if
				nsectorjumps(it)=nsectorjumps(it) + 1
				
				if(fill)then
					xsectjump(nsectorjumps(it))%isector = is
					xsectjump(nsectorjumps(it))%fsector = fs
					xsectjump(nsectorjumps(it))%dJz = dJz
					xsectjump(nsectorjumps(it))%dParity = dParity
                    dW = xsd_f(it)%sector(fs)%Wx -xsd_i(it)%sector(is)%Wx
					xsectjump(nsectorjumps(it))%dW = dW 	
					xsectjump(nsectorjumps(it))%wXf = xsd_f(it)%sector(fs)%Wx 			
				end if
				
			end do ! fs
			
		end do
		if(.not.fill)then
			print*,nsectorjumps(it),' sector jumps for species ',it
		endif
		
		return
	end subroutine countsectorjumps
	
!........... SUBROUTINES TO DEFINE AND SET UP ONE-BODY OPERATORS................	
!            FIRST COUNT UP # OF DISTINCT SETS OF QUANTUM NUMBERS	

subroutine countup_quantum_numbers(it,nsjumps,xsjumps,nqns,fill,xops)
	use sporbit,only:sameweightscheme
	implicit none
	integer,intent(in) :: it
	integer,intent(in) :: nsjumps
	type (basejump),intent(inout) :: xsjumps(nsjumps)
	integer(4),intent(out) :: nqns
	logical,intent(in) :: fill
	type (baseoperator),intent(inout),pointer :: xops(:)
	
	integer :: isj,jsj
	logical :: foundcopy	
	integer :: ndW,idW,dW	
	integer :: iqns,n
	

	if(fill)then  ! NEED TO ALLOCATE ARRAY OF ALLOWED OPERATOR QUANTUM NUMBERS
		do isj = 1,nsjumps
			n = xsjumps(isj)%nop
			xsjumps(isj)%iop=0

		end do
	end if
    xsjumps(:)%nop = 0		
	nqns = 0

	do isj = 1,nsjumps

		    foundcopy = .false.
			dW = xsjumps(isj)%dW

!			print*,isj,' quantum number (X)',xsjumps(isj)%dJz

            if(isj==1 .or. .not. sameweightscheme)then
				foundcopy = .false.
			else		! look to see if found previously
			   do jsj = 1,isj-1
				  if(xsjumps(isj)%dJz /= xsjumps(jsj)%dJz .or. & 
			   	   xsjumps(isj)%dParity /= xsjumps(jsj)%dParity)cycle
				   if( ( sameweightscheme .and. dW == xsjumps(jsj)%dW))then 
			          foundcopy = .true.
			          if(fill)then
				          xsjumps(isj)%iop = xsjumps(jsj)%iop	 
			          end if
			   
			          exit
		           end if
					
		       end do  ! jsj
     		end if
		    xsjumps(isj)%nop = xsjumps(isj)%nop+1
		    if(.not.foundcopy)then
			    nqns = nqns+1
			    if(fill)then
				   xops(nqns)%d2M = xsjumps(isj)%dJz
				   xops(nqns)%dParity = xsjumps(isj)%dParity
				   xops(nqns)%dW = dW
				   xsjumps(isj)%iop = nqns
			    end if		

				
		    end if
	 end do  ! isj
	
	return
end subroutine countup_quantum_numbers

!......................................................................
!
!  COUNT UP (AND FILL) ONE-BODY JUMPS IN A SECTOR JUMP
!
!      apply ac^+ bd  (ac= creation, bd = destruction)
!  to count up:
!     loop over SDs in initial sector
!         loop over allowed 1-body operators
!           check site bd is occupied, site ac unoccupied or = bd
!
! ALTERNATELY for 'hops' just apply ac^+  or bd (create/destroy one particle)


subroutine count_fill_1body_sectorjumps(it,xops1body,fill,xsjump)
	use basis
	use sectors
	use spstate
	use system_parameters
	use basis_map
	use sporbit
	implicit none
	integer,intent(in)  :: it
	logical,intent(in)  :: fill
	type (basejump),intent(inout) :: xsjump
	
	type (baseoperator) :: xops1body(:)
	integer :: isector,fsector
		
	integer :: iop1body	
	integer(8) :: i,f
	integer :: iop
	integer :: ac,bd
	integer(8) :: njumps
	integer(8),pointer :: isd(:,:),fsd(:,:)
    type (sectorinfo), pointer :: sector
	integer(4)  :: icount
	integer  :: create(1),destroy(1)
	integer(8) :: sdin(Nword_i(it)),sdout(Nword_f(it))
	integer   :: occin(np_i(it)),occout(np_f(it))
	integer(2) :: phase
	integer(4) :: Mxf, Wxf,parxf
	logical :: success,foundit,search4it
	
	integer :: ndestroy, ncreate  
	
	if(np_i(it)==np_f(it))then
		ndestroy = 1
		ncreate =  1
	elseif(np_i(it) -1 == np_f(it))then
		ndestroy = 1
		ncreate = 0
	else
		ndestroy=0
		ncreate = 1
	end if
	
	needphase=fill
	isector  = xsjump%isector
	fsector  = xsjump%fsector
	njumps = 0
	

    iop1body = xsjump%iop
    if(iop1body==0)then
		print*,' missing operators '
		return
	end if
    do icount = 1,xsjump%nop

	   if(it==1)then
		  isd=>psdlist_i
		  fsd=>psdlist_f
	   else
		  isd=>nsdlist_i
		  fsd=>nsdlist_f
	   end if

	   sector => Xsd_i(it)%sector(isector)
	   Mxf = xsd_i(it)%sector(isector)%jzX + xops1body(iop1body)%d2M
	   parxf = parmult(xsd_i(it)%sector(isector)%parX, xops1body(iop1body)%dparity)
	
	   do i = sector%xsdstart,sector%xsdend
		  do iop = 1,xops1body(iop1body)%nops
				bd = xops1body(iop1body)%Xop(iop,1)     ! destruction operator label, in initial space
				ac = xops1body(iop1body)%Xop(iop,2) 	! creation operator label, in final space		
!............... GENERAlIZE TO ACCOUNT FOR HOPS where either ac or bd = 0...................				

                search4it = .false.
				if(ac==0 .and. bd==0)then
					print*,' some problem, too many zero operators '
					stop
				end if
				if( ac /= 0 .and. bd /= 0)then   ! one-body operator; check that BD is occupied and AC is not occupied in ISD
					search4it=	(check_bitocc('INI',isd(i,:),bd) .and. & 
				                (ac==rhspsqn_i(it,bd)%ifmap .or. .not. check_bitocc('INI',isd(i,:),rhspsqn_f(it,ac)%ifmap)) )
				else
					if(ac == 0)then						
						search4it=check_bitocc('INI',isd(i,:),bd)   ! check that BD is occupied in ISD
					else						
						search4it=.not. check_bitocc('INI',isd(i,:),rhspsqn_f(it,ac)%ifmap)   ! check that AC is not occupied in ISD
					end if
					 
			    end if
!				print*,' searching? ',search4it,fill,sameweightscheme

                if(search4it)then
				    if(.not.sameweightscheme .or. fill)then  ! replace and search
						create(1)=ac
						destroy(1)=bd
						call sdmapperplus(it,isd(i,:),occin,ndestroy,destroy,ncreate,create,sdout,success,occout,phase)				
						if(success)then
		 				    Wxf = w_of_occ_f(it,np_f(it),occout)
							
							! NEED TO CHECK THIS IS IN THE CORRECT FINAL SECTOR
							if(Wxf /= Xsd_f(it)%sector(fsector)%Wx)cycle
							
							call searcher(it,sdout,Mxf,parxf,Wxf,foundit,f)
							if(.not.foundit)cycle
							njumps = njumps + 1
							
							if(fill)then
								xsjump%phase(njumps)=phase
								xsjump%isd(njumps)=i
								xsjump%fsd(njumps)=f
								xsjump%iXop(njumps)=iop
																											
							end if
						end if
					else   ! we automatically know it exists
						njumps = njumps + 1						
						
				    end if			
				end if
				
		   end do ! iop		
	   end do ! i
    end do ! icount
	if(fill .and. njumps /=xsjump%njumps)then
		print*,' search failed to get correct number of jumps ',njumps,xsjump%njumps
		stop
	end if
	xsjump%njumps = njumps
	if(.not.fill)then  ! allocate
         allocate(xsjump%phase(njumps))
         allocate(xsjump%isd(njumps))
         allocate(xsjump%fsd(njumps))
         allocate(xsjump%iXop(njumps))
		
	end if
			
	return
end subroutine count_fill_1body_sectorjumps

!==================================================================
!
! LOOP OVER INITIAL ORBITS AND SEE IF THEY HAVE THE SAME WEIGHT SCHEME AS FINAL
!
subroutine checkorbitweightscheme
	use sporbit
	implicit none
	
	integer :: it
	integer :: i,f
	
	sameweightscheme = .true.
	
	do it = 1,2
		do i = 1,numorb_i(it)
			f = orbqn_i(it,i)%ifmap
			if(f==0)cycle
			if(orbqn_i(it,i)%w/=orbqn_f(it,f)%w)then
				sameweightscheme=.false.
				return
			end if
		end do		
		
	end do
	
	return
	
end subroutine checkorbitweightscheme	
!==================================================================

	
end module jumpbase
