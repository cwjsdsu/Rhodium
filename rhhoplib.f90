!
!  RHODIUM hops--single creation/annihilation
!  started late July 2018 by CWJ @ SDSU
!
! NOTE Aug 2018: one potential problem is change in truncation wt W
! Nov 2018: If this is the case, must search through all Ws
! 

module hops
	use precisions
	use jumpbase
	implicit none
	
	type (basejump), allocatable, target :: psecthop(:),nsecthop(:)
			
	integer :: dnp(2)  ! change in proton, neutron numbers	
		
!........... NEED TO DEFINE ONE-BODY OPERATORS..........
!            first define change in quantum numbers
	
	type (baseoperator), pointer :: phopop(:),nhopop(:)
	integer(4),target :: nphopops,nnhopops   ! # of unique quantum numbers
	
	integer :: dWhopmin(2),dWhopmax(2)   ! limits of how much W can change
	
	logical,parameter ::  verbosehops=.false.
	
contains


!===================================================================
!........ SET UP HOPS ..................................
!  of the form a^dagger(ac) a(bd)
!               a_create   b_destroy
!
subroutine setup_hopops(it,xhopop,fill)
	use spstate
	use sporbit
	
	implicit none
	integer,intent(in) :: it ! species: it = 1 protons, 2 neutrons
	type (baseoperator) :: xhopop
	logical,intent(in) :: fill  
	
	integer :: d2M,dW,dParity
	integer :: nhops
	integer :: ac, bd
	
	
	nhops = 0
	
	d2M = xhopop%d2M
	dParity=xhopop%dParity
	dW = xhopop%dW
	
   if(dnp(it)< 0)then   ! destroy a particle
!.......... LOOP OVER (INITIAL) DESTRUCTION OPERATOR
        ac = 0
		if(verbosehops)print*,' destroying a particle '
        do bd = 1,nrhsps_i(it)


			 if( - rhspsqn_i(it,bd)%m /= d2M)cycle
			 if( sameweightscheme .and. - rhspsqn_i(it,bd)%W /= dW)cycle

			 if(  rhspsqn_i(it,bd)%par /= dParity)cycle
			 
			 nhops = nhops + 1
			 if(verbosehops)print*,' counting....', nhops
			 if(fill)then
				 xhopop%xop(nhops,1)=bd
				 xhopop%xop(nhops,2)=ac	 ! = 0
			 end if			 
        end do   ! bd
	else
		bd = 0
		if(verbosehops)print*,' creating a particle '

!.......... LOOP OVER (FINAL) CREATION OPERATOR	

         do ac = 1,nrhsps_f(it)
			 if( rhspsqn_f(it,ac)%m /= d2M)cycle
			 if( sameweightscheme .and. rhspsqn_f(it,ac)%W  /= dW)cycle

			 if( rhspsqn_f(it,ac)%par  /= dParity)cycle
			 
			 nhops = nhops + 1
			 if(verbosehops)print*,' counting....', nhops
			 
			 if(fill)then
				 xhopop%xop(nhops,1)=bd  ! = 0
				 xhopop%xop(nhops,2)=ac	 
			 end if			 
		 end do
	end if
	if(.not.fill)then
!		print*,it,' there are ',nhops,' 1-body operators ',d2M,dParity,dW
		xhopop%nops = nhops
		allocate(xhopop%Xop(nhops,2))
	end if
	
	return
	
end subroutine setup_hopops



!==================================================================
! find min, max of change in W
!
subroutine find_hops_dWminmax
	use sporbit
	implicit none
	integer :: it
	integer :: i,f,Wi,Wf
	
	do it = 1,2
	    dWhopmin(it) = 0
	    dWhopmax(it) = 0
		if(dnp(it)==0)cycle
	    dWhopmin(it) = 100000
	    dWhopmax(it) = -100000
		
		if(dnp(it)> 0)then  ! adding a particle
		   do f = 1,numorb_f(it)
				Wf = orbqn_f(it,f)%W
				
				dWhopmax(it) = max(dWhopmax(it),Wf)
				dWhopmin(it) = min(dWhopmin(it),Wf)			
		   end do
	    else			! subtracting a particle
		   do i = 1,numorb_i(it)
			  Wi = orbqn_i(it,i)%W
			dWhopmax(it) = max(dWhopmax(it),Wi)
			dWhopmin(it) = min(dWhopmin(it),Wi)		
	
		   end do
	   end if
	end do
	return
	
end subroutine find_hops_dWminmax

!==================================================================
!
!.......... BOSS ROUTINE to set up hops...............................
!
!
!  SUBROUTINES CALLED
!   find_hops_dWminmax
!   checkorbitweightscheme
!   countectorjumps
!   countup_quantum_numbers
!   setup_1b_operators
!   count_fill_1body_sectorjumps
!
!  CALLED BY
!     main_menu
!
subroutine bosshops(fixedJz,fixedpar)
	use system_parameters
	use sporbit,only:sameweightscheme
	
	implicit none
	logical :: fixedJz,fixedpar
	integer :: it
	
	integer :: iqs,isj
	integer(8) :: tothops
	

    do it = 1,2
		
		dnp(it) = np_f(it)-np_i(it)
	
	end do 
	
	nsectorjumps(:)= 0
	call find_hops_dWminmax
	call checkorbitweightscheme
	if(sameweightscheme)then
		print*,' orbits have same weights '
	else
		print*,' orbits have different weights '
		
	end if

	do it = 1,2
		if(dnp(it)==0)cycle
		
		print*,' hops for species ',it
		tothops = 0
!		if(np_i(it)==0)cycle   ! THIS NEEDS TO BE CHANGED TO NPEFF
		if(it==1)then
			call countsectorjumps(it,psecthop,fixedJz,fixedpar,.false.)
			
			allocate(psecthop(nsectorjumps(it)))
			psecthop%nop = 0
			call countsectorjumps(it,psecthop,fixedJz,fixedpar,.true.)
			
		else
			call countsectorjumps(it,nsecthop,fixedJz,fixedpar,.false.)
			allocate(nsecthop(nsectorjumps(it)))
			nsecthop%nop = 0
			call countsectorjumps(it,nsecthop,fixedJz,fixedpar,.true.)
			
			
		end if
		if(it==1)then			
		   call countup_quantum_numbers(it,nsectorjumps(it),psecthop,nphopops,.false.,phopop)
		   print*,nphopops,' proton operators'
		   allocate(phopop(nphopops))
		   call countup_quantum_numbers(it,nsectorjumps(it),psecthop,nphopops,.true.,phopop)		  
		   print*,nphopops,' quantum numbers for species ',it
			
		   do iqs = 1,nphopops
		      call setup_hopops(it,phopop(iqs),.false.)
		      call setup_hopops(it,phopop(iqs),.true.)

 		   end do		   
! WOULD LIKE TO OPENMP HERE...BUT CAUSES PROBLEMS... in routine convert_occ2bitrepsd  (not sure of problem)
		   do isj = 1,nsectorjumps(it)
			   call count_fill_1body_sectorjumps(it,phopop,.false.,psecthop(isj))
			   call count_fill_1body_sectorjumps(it,phopop,.true.,psecthop(isj))
		   end do 
		   do isj = 1,nsectorjumps(it)  !COUNT UP
			   tothops = tothops+psecthop(isj)%njumps
		   end do 
		   
	   else
		   
		   call countup_quantum_numbers(it,nsectorjumps(it),nsecthop,nnhopops,.false.,nhopop)
		   print*,nnhopops,' neutron operators'
		   
		   allocate(nhopop(nnhopops))
		   call countup_quantum_numbers(it,nsectorjumps(it),nsecthop,nnhopops,.true.,nhopop)
		   print*,nnhopops,' quantum numbers for species ',it
		   do iqs = 1,nnhopops
 		      call setup_hopops(it,nhopop(iqs),.false.)
 		      call setup_hopops(it,nhopop(iqs),.true.)
 		   end do
! WOULD LIKE TO OPENMP HERE...BUT CAUSES PROBLEMS
		   do isj = 1,nsectorjumps(it)	  
			   call count_fill_1body_sectorjumps(it,nhopop,.false.,nsecthop(isj))
			   call count_fill_1body_sectorjumps(it,nhopop,.true.,nsecthop(isj))
		   end do
		   
		   do isj = 1,nsectorjumps(it)	  
			   tothops = tothops+nsecthop(isj)%njumps
		   end do
		   
	   end if	   
	   print*,' total hops = ',tothops
	end do
	
	print*,' finished hopping '

	return
	
end subroutine bosshops
	
	
end module hops
