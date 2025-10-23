!
!  RHODIUM 1-body jumps
!  started late July 2018 by CWJ @ SDSU
!
! NOTE Aug 2018: one potential problem is change in truncation wt W
! Nov 2018: If this is the case, must search through all Ws
! 

module jumps1body
	use precisions
	use jumpbase

	implicit none
	
	type (basejump), allocatable, target :: p1bsectjump(:),n1bsectjump(:)
		
		
!........... NEED TO DEFINE ONE-BODY OPERATORS..........
!            first define change in quantum numbers
	
	type (baseoperator), pointer :: pops1body(:),nops1body(:)
	integer(4),target :: npops1body,nnops1body   ! # of unique quantum numbers
	
	integer :: dW1bmin(2),dW1bmax(2)   ! limits of how much W can change
	
	
contains


!===================================================================
!........ SET UP ONE-BODY OPERATORS..................................
!  of the form a^dagger(ac) a(bd)
!               a_create   b_destroy
!
subroutine setup_1b_operators(it,xops1body,fill)
	
	use spstate
	use sporbit
	
	implicit none
	integer,intent(in) :: it ! species: it = 1 protons, 2 neutrons
	type (baseoperator) :: xops1body
	logical,intent(in) :: fill  
	
	integer :: d2M,dW,dParity
	integer :: n1bops
	integer :: ac, bd
	
	n1bops = 0
	
	d2M = xops1body%d2M
	dParity=xops1body%dParity
	dW = xops1body%dW	
	
!.......... LOOP OVER (INITIAL) DESTRUCTION OPERATOR

    do bd = 1,nrhsps_i(it)

!.......... LOOP OVER (FINAL) CREATION OPERATOR	

         do ac = 1,nrhsps_f(it)
			 if( rhspsqn_f(it,ac)%m - rhspsqn_i(it,bd)%m /= d2M)cycle
			 if( sameweightscheme .and. rhspsqn_f(it,ac)%W - rhspsqn_i(it,bd)%W /= dW)cycle

			 if( parmult(rhspsqn_f(it,ac)%par , rhspsqn_i(it,bd)%par) /= dParity)cycle
			 
			 n1bops = n1bops + 1
			 
			 if(fill)then
				 xops1body%xop(n1bops,1)=bd
				 xops1body%xop(n1bops,2)=ac	 
			 end if			 
		 end do
    end do   ! bd
	
	if(.not.fill)then
		xops1body%nops = n1bops
		allocate(xops1body%Xop(n1bops,2))
	end if
	
	return
	
end subroutine setup_1b_operators



!==================================================================
! find min, max of change in W
!
subroutine find_1body_dWminmax
	use sporbit
	implicit none
	integer :: it
	integer :: i,f,Wi,Wf
	
	do it = 1,2
	dW1bmin(it) = 100000
	dW1bmax(it) = -100000
	
		do i = 1,numorb_i(it)
			Wi = orbqn_i(it,i)%W
			
			do f = 1,numorb_f(it)
				Wf = orbqn_f(it,f)%W
				
				dW1bmax(it) = max(dW1bmax(it),Wf-Wi)
				dW1bmin(it) = min(dW1bmin(it),Wf-Wi)			
			end do			
		end do
	end do
	return
	
end subroutine find_1body_dWminmax

!==================================================================
!
!.......... BOSS ROUTINE to set up 1-body jumps................................
!
!
!  SUBROUTINES CALLED
!   find_1body_dWminmax
!   checkorbitweightscheme
!   count1bsectorjumps
!   countup_quantum_numbers
!   setup_1b_operators
!   count_fill_1body_sectorjumps
!
!  CALLED BY
!     main_menu
!
subroutine boss1bjumps
	use system_parameters
	use sporbit !,only:sameweightscheme
	
	implicit none
	integer :: it
	
	integer :: iqs,isj
	integer(8) :: totjumps,alljumps
	logical :: fixedJz,fixedpar
	
	
	fixedJz = .true.
	fixedpar = .true.
	nsectorjumps(:)= 0
	call find_1body_dWminmax
	call checkorbitweightscheme
	if(sameweightscheme)then
		print*,' orbits have same weights '
	else
		print*,' orbits have different weights '
		
	end if

    alljumps = 0
	do it = 1,2
		totjumps = 0
		if(np_i(it)==0)cycle   ! 
		if(it==1)then
			call countsectorjumps(it,p1bsectjump,fixedJz,fixedpar,.false.)
			
			allocate(p1bsectjump(nsectorjumps(it)))
			p1bsectjump%nop = 0
			call countsectorjumps(it,p1bsectjump,fixedJz,fixedpar,.true.)
			
		else
			call countsectorjumps(it,n1bsectjump,fixedJz,fixedpar,.false.)
			
			allocate(n1bsectjump(nsectorjumps(it)))
			n1bsectjump%nop = 0
			call countsectorjumps(it,n1bsectjump,fixedJz,fixedpar,.true.)
			
		end if
		if(it==1)then			
		   call countup_quantum_numbers(it,nsectorjumps(it),p1bsectjump,npops1body,.false.,pops1body)
!		   print*,' There are ',npops1body,' unique sets of proton quantum numbers'
		   allocate(pops1body(npops1body))
		   call countup_quantum_numbers(it,nsectorjumps(it),p1bsectjump,npops1body,.true.,pops1body)		  
		   
		   print*,npops1body,' operators for protons '
		   if(npops1body==0)then
			   print*,' There are NO one-body proton operators that can connect states '
			   print*,' Most likely the change in Jz is too large '
			   
		   end if
		    
		   do iqs = 1,npops1body
		      call setup_1b_operators(it,pops1body(iqs),.false.)
		      call setup_1b_operators(it,pops1body(iqs),.true.)
 		   end do		   
		   
		   do isj = 1,nsectorjumps(it)
			   
			   call count_fill_1body_sectorjumps(it,pops1body,.false.,p1bsectjump(isj))
			   totjumps = totjumps+p1bsectjump(isj)%njumps
			   call count_fill_1body_sectorjumps(it,pops1body,.true.,p1bsectjump(isj))
			   
		   end do 
		   
	   else
		   
		   call countup_quantum_numbers(it,nsectorjumps(it),n1bsectjump,nnops1body,.false.,nops1body)
		   allocate(nops1body(nnops1body))
		   call countup_quantum_numbers(it,nsectorjumps(it),n1bsectjump,nnops1body,.true.,nops1body)
		   
		   if(nnops1body==0)then
			   print*,' There are NO one-body neutron operators that can connect states '
			   print*,' Most likely the change in Jz is too large '
			   cycle
			   
		   end if
		   
!$OMP PARALLEL DO private(iqs)		   

		   do iqs = 1,nnops1body
		      call setup_1b_operators(it,nops1body(iqs),.false.)
		      call setup_1b_operators(it,nops1body(iqs),.true.)
 		   end do
!$OMP END PARALLEL DO		   

		   do isj = 1,nsectorjumps(it)
			  
			   call count_fill_1body_sectorjumps(it,nops1body,.false.,n1bsectjump(isj))
			   totjumps = totjumps+n1bsectjump(isj)%njumps
			   call count_fill_1body_sectorjumps(it,nops1body,.true.,n1bsectjump(isj))
			   
		   end do

	   end if	  
	   

	   print*,' total jumps = ',totjumps
	   alljumps = alljumps + totjumps
	   	   
	end do
    if(npops1body==0 .and. nnops1body==0)then
	   print*,' No one-body operators can connect the initial, final basis '
	   print*,' Most likely the change in Jz is too large '
	   print*,' Stopping now '
	   stop
    end if
    if(alljumps==0)then
		print*,' '
		print*,' * * * * * * * * * * * * * * * * *'
		print*,' '

	   print*,' No one-body jumps can connect the initial, final basis '
	   print*,' Most likely the change in Jz is too large '
	   print*,' Stopping now '
	print*,' '
	print*,' * * * * * * * * * * * * * * * * *'
	print*,' '
	   stop
    end if

	return
	
end subroutine boss1bjumps	
	
	
end module jumps1body
