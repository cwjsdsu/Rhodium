!
!  data and routines for mapping initial -> final:
!    single particle states;
!    Slater determinants
!
!   Note: biggest issue is out-of-order of single particle states, which
!   can occur with different truncations
!
!  There are a number of maps which are needed. There are maps between initial and final basis data.
!  In addition, because of the way BIGSTICK groups data, we often need to map information.
!  For example, the order in which s.p. state occupations are defined internally
!  follow the so-called haiku ordering, which is given by the hspsqn_x information
!  Thus we have to map between hspsqn and spsqn.
!    
!
module spmap_mod
	use spstate
	implicit none
		
contains
!-----------------------------------------------------------------------------
!  this routine sets up the rhspsqn_x data 
!  these are the actual single particle states used and in the order
!  used to create the Slater determinants
!-----------------------------------------------------------------------------
subroutine setup_rhsps(whichbasis)
	use bitstuff
	use basis,only:Nword_i,Nword_f
	use sporbit
	use nodeinfo !mpi
	implicit none
	
	character(3) :: whichbasis
	integer :: it,ith
	integer :: isps,jsps,hsgn
	integer, pointer :: nspsx(:),nhspsx(:),nrhspsx(:)
	type (spst),pointer :: hspsqnx(:,:),rhspsqnx(:,:)
    type (orb),pointer :: myorbqnx(:,:)
	integer,pointer :: numorb_x(:)
	integer :: nunused  ! any unused s.p. states?
	integer :: nrhspsmax
	integer :: Nwordx
	integer :: iorb

	if(whichbasis=='INI')then
		nspsx=> nsps_i
		nhspsx=>nhsps_i
		nrhsps_i(1)=nhsps_i(1)+nhsps_i(-1)
		nrhsps_i(2)=nhsps_i(2)+nhsps_i(-2)
		nrhspsmax = max(nrhsps_i(1),nrhsps_i(2))
		nrhspsx=>nrhsps_i
		allocate(rhspsqn_i(2,nrhspsmax))
		rhspsqnx=>rhspsqn_i
		hspsqnx=>hspsqn_i	
		myorbqnx=> orbqn_i
		numorb_x=> numorb_i
	else
		nspsx=> nsps_f
		nhspsx=>nhsps_f
		nrhsps_f(1)=nhsps_f(1)+nhsps_f(-1)
		nrhsps_f(2)=nhsps_f(2)+nhsps_f(-2)
		nrhspsmax = max(nrhsps_f(1),nrhsps_f(2))
		nrhspsx=>nrhsps_f
		allocate(rhspsqn_f(2,nrhspsmax))
		
		rhspsqnx=>rhspsqn_f
		hspsqnx=>hspsqn_f		
		myorbqnx=> orbqn_f
		numorb_x=> numorb_f
			
	end if
	
	nunused = nspsx(1)-nhspsx(1)-nhspsx(-1)
	if(nunused /= 0)then
		print*,whichbasis,nunused,' unused proton s.p. states '
	end if
	nunused = nspsx(2)-nhspsx(2)-nhspsx(-2)
	if(nunused /= 0)then
		print*,whichbasis,nunused,' unused neutron s.p. states '
	end if
	rhspsqnx(:,:)%ifmap = 0
	
	
	do it = 1,2
		jsps = 0
		do hsgn = -1,1,2
			ith =it*hsgn
			do isps = 1,nhspsx(ith)
				jsps = jsps+1
				rhspsqnx(it,jsps)%nr = hspsqnx(ith,isps)%nr
				rhspsqnx(it,jsps)%j  = hspsqnx(ith,isps)%j
				rhspsqnx(it,jsps)%m  = hspsqnx(ith,isps)%m
				rhspsqnx(it,jsps)%w  = hspsqnx(ith,isps)%w
				rhspsqnx(it,jsps)%l  = hspsqnx(ith,isps)%l
				rhspsqnx(it,jsps)%par= (3-hspsqnx(ith,isps)%par)/2   ! CONVERT
!				print*,' TESTING PARITY ', hspsqnx(ith,isps)%par
				rhspsqnx(it,jsps)%ifmap= 0
				
!..................... MAP TO AN ORBIT..............	
                rhspsqnx(it,jsps)%orb=0			
                do iorb = 1,numorb_x(it)
					if(rhspsqnx(it,jsps)%nr /=myorbqnx(it,iorb)%nr)cycle
					if(rhspsqnx(it,jsps)%j /=myorbqnx(it,iorb)%j)cycle
					if(rhspsqnx(it,jsps)%l /=myorbqnx(it,iorb)%l)cycle
					rhspsqnx(it,jsps)%orb =iorb
					exit
				end do
				if(rhspsqnx(it,jsps)%orb==0)then
					print*,' could not find orbit '
					stop
				end if
				
			end do  ! isps
			
		end do  ! hsgn
		
	end do ! it
	
! .... NOW SET UP # OF WORDS NEEDED FOR BIT REPRESENTATION...

do it =1 ,2
     Nwordx = nrhspsx(it)/max_bit_word+1
	 
!	 print*,nrhspsx(it),max_bit_word,nrhspsx(it)/max_bit_word, Nwordx
	 if(whichbasis=='INI')then
		 Nword_i(it)=Nwordx
	 else
		 Nword_f(it)=Nwordx
	 end if
	 print*,whichbasis,it,nwordx,' Nword '
	 
 end do  	

	return			
end subroutine setup_rhsps

!=====================================================	
!
!  MAPS INITIAL single-particle states to final single particle states
!
!  One thing important to note: the spsqn arrays are defined by the initial sps definition
!  while the hspsqn states are defined by the states actually used
!  So for example in truncations not all s.p. states may be used.
!
!  So we use the data bundle rhspsqn_x (x=i,f) with only the used single-particle states
! 
!   THIS ROUTINE maps rhspsqn_x initial <-> final via rhspsqn_x%ifmap
!
	subroutine map_initial2final_rhsps
		use nodeinfo !mpi
		use system_parameters,only:phconj_i
		implicit none
		integer :: isps,fsps
		integer :: it
		
		integer :: unused  ! # of unused states
		
!......... COUNT UP UNUSED STATES..................	
unused = nrhsps_i(1)-nrhsps_f(1)
if(iproc == 0)print*,' Initial proton s.p. spaces uses ',unused,' more states than final s.p. space'

unused = nrhsps_i(2)-nrhsps_f(2)
if(iproc == 0)print*,' Initial neutron s.p. spaces uses ',unused,' more states than final s.p. space'
!.....................................	
		
		do it = 1,2
			do isps = 1,nrhsps_i(it)
				do fsps=1,nrhsps_f(it)
					if(rhspsqn_i(it,isps)%nr == rhspsqn_f(it,fsps)%nr .and. & 
					   rhspsqn_i(it,isps)%j == rhspsqn_f(it,fsps)%j .and. &
					   rhspsqn_i(it,isps)%m == rhspsqn_f(it,fsps)%m .and. &
					   rhspsqn_i(it,isps)%l == rhspsqn_f(it,fsps)%l .and. &
					   rhspsqn_i(it,isps)%par == rhspsqn_f(it,fsps)%par )then
					   
				      if(rhspsqn_i(it,isps)%ifmap/=0 .or. rhspsqn_f(it,fsps)%ifmap/=0)then
						  print*,' already mapped, huh? '
						  stop
					  end if
				      rhspsqn_i(it,isps)%ifmap = fsps
				      rhspsqn_f(it,fsps)%ifmap = isps
					   
					   

						  
						  exit 
					  end if
				end do
			end do			
			! NEW! ADDED in 0.9.1   FIND TIME-REVERSE if needed for hole-formalism		
			rhspsqn_i(it,:)%tr = -1
			rhspsqn_f(it,:)%tr = -1
			if(.not.phconj_i(it))cycle
			do isps = 1,nrhsps_i(it)
				do fsps=1,nrhsps_f(it)
					if(rhspsqn_i(it,isps)%nr == rhspsqn_f(it,fsps)%nr .and. & 
					   rhspsqn_i(it,isps)%j == rhspsqn_f(it,fsps)%j .and. &
					   rhspsqn_i(it,isps)%m == -rhspsqn_f(it,fsps)%m .and. &
					   rhspsqn_i(it,isps)%l == rhspsqn_f(it,fsps)%l .and. &
					   rhspsqn_i(it,isps)%par == rhspsqn_f(it,fsps)%par )then
					   					      rhspsqn_i(it,isps)%tr = fsps
					   					      rhspsqn_f(it,fsps)%tr = isps
						  
						  exit 
					  end if
				end do
			end do			

		end do
		
		return

	end subroutine map_initial2final_rhsps
	
end module spmap_mod
