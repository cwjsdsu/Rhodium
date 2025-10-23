!
!  MENU OPTIONS: 6 letter code
!  see manual

module menu_mod
	implicit none
	
	character(6) :: menu_choice
	
contains
	
	
	subroutine main_menu
		use dot_mod
		use proj_mod
		use io
		use jumps1body
		use dens1body
		use dens1bodycc
		use nodeinfo ! mpi
		use entropy1body
		use hops
		use spfac
		use system_parameters
		use apply1body
		use apply1bodycc
		use bmpi_mod
		use sporbit,only:isoflag
		use lincom_mod !for wfn linear combination
		implicit none
		logical :: success
		logical :: change1,changepn
		logical :: combinedlist, samebasis
		integer :: ierr
		character :: ychar
		
		change1= .false.
		changepn = .false.
		if(np_i(1)/=np_f(1) .and. np_i(2)==np_f(2))change1=.true.
		if(np_i(2)/=np_f(2).and. np_i(1)==np_f(1))change1=.true.
		
		if(np_i(1)==np_f(1)+1 .and. np_i(2)==np_f(2)-1)changepn=.true.
		if(np_i(2)==np_f(2)+1 .and. np_i(1)==np_f(1)-1)changepn=.true.
		
		if(np_i(1)==np_f(1) .and. np_i(2)==np_f(2)  & 
		.and. Jz_i == Jz_f .and. iparity_i == iparity_f)then
		    samebasis = .true.
	    else
			samebasis = .false.
  	    end if
		
		if(iproc == 0) then 
		write(6,*)' '
		write(6,*)' ~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~ '
		write(6,*)' '
		end if
		success = .false.

		do while(.not.success)
			if(iproc == 0) then
				write(6,*)' Choose one of the following menu options '
				write(6,*)'  '
			endif
			if(change1)then
			    if(iproc == 0) then
198 continue
 				     write(6,*)' Enter the six-character code for your choice::'
				     write(6,*)' (spamp1) 1-particle spectroscopic amplitude (change in number of particles)'
				     write(6,*)' (appsp1) Add/remove a single particle (linear combination)'
				!	 write(6,*)' (allsp1) Add/remove a single particle in different orbitals from initial wfn (see manual)'
	 				 read(5,'(a)')menu_choice
	                 if(menu_choice /= 'spamp1' .and. menu_choice /='appsp1' .and. menu_choice /='allsp1')then
						 print*,' Try again '
						 goto 198
					 end if
				endif


				elseif(changepn)then
				if(iproc == 0)	then
	!				write(6,*)' [CALCULATE CHARGE-CHANGING MATRIX ELEMENTS BETWEEN STATES]'
					write(6,*)' (dn1bcc) Charge-changing 1-body densities (pn)'	
				endif
				if(isoflag)then
					if(iproc == 0)	then
						write(6,*)' (d1iscc) Charge-changing 1-body densities (iso)'	

						
					endif
					
						
					
				end if
				write(6,*)' (app1cc) Apply charge-changing 1-body operator to wave functions and write out '	
				
				write(6,*)'  '
			
				write(6,*)' Enter the six-character code for your choice::'
				read(5,'(a)')menu_choice
					
				else
			if(iproc == 0) then
			   write(6,*)' [CALCULATE MATRIX ELEMENTS BETWEEN STATES]'
			   write(6,*)' (dotwfn) Dot product between wave functions in different bases (different truncations)'
			   write(6,*)' (dn1iso) 1-body density matrices, non-charge changing, good isospin '
			   write(6,*)' (dn1xpn) 1-body density matrices, non-charge changing, explicit proton-neutron '			
			   write(6,*)' (ent1bd) 1-body entropy '
			   write(6,*)' (app1bd) Apply non-charge-changing 1-body operator to wave functions and write out'
			   write(6,*)' (projct) Project a state into a basis with different truncation '
			   write(6,*)' (lincom) construct linear combinations of wave functions '
			   if(Jz_f == Jz_i +2)then
				   write(6,*)' (jraise) Apply J+ raising operator to wave functions and write out '
				   write(6,*)' (dn1jrm) Compute missing density matrix elements by applying J+ raising'
				   write(6,*)' (dn1mis) Compute missing density matrix elements '
				   
			   end if

			endif
		   end if
						
			write(6,*)'  '
			
			if(.not.change1 .and. .not. changepn)then
				if(iproc==0) then
					write(6,*)' Enter the six-character code for your choice::'
				
					read(5,'(a)')menu_choice
				endif
		    endif
		
			call BMPI_BARRIER(icomm,ierr)
  			call BMPI_BCAST(menu_choice,6,0,icomm,ierr)
			
			select case (menu_choice)
			
			case('dotwfn', 'DOTWFN') !***************************************************
			if(iproc==0) then
				print*,' You are calculating the dot (inner) product between wave functions, in bases with different truncations '
			endif
			success = .true.
			
			call dot_master
			
			case('dn1iso','DN1ISO')	  !***************************************************
			if(iproc==0) then		
				print*,' You are calculating non-charge-changing one-body density matrices with good isospin '
				print*,' Do you want to combine initial/final levels into a single list? (y/n)'
				read(5,'(a)')ychar
			endif
			call BMPI_BCAST(ychar,1,0,icomm,ierr)
			if(ychar=='y' .or. ychar=='Y')then
				combinedlist=.true.
				if(iproc==0)print*,' Combining initial/final levels into a single list'
			else
				combinedlist = .false.
				if(iproc==0)print*,' Keeping initial/final levels in separate lists'
				
			end if
			
			success = .true.
			call getoutputfile
			
			call boss1bjumps
			pndensities=.false.
			if(combinedlist)then
				call density1b_bag_boss(.false.,.false.,samebasis)
			else
   			  call density1b_boss(.false.,samebasis)
		    end if
			


			case('dn1xpn','DN1XPN') !***************************************************
			success = .true.
			if(iproc==0) then
				print*,' You are calculating non-charge-changing one-body density matrices, proton-neutron formalism'	
				print*,' Do you want to combine initial/final levels into a single list? (y/n)'
				read(5,'(a)')ychar
            endif
			call BMPI_BCAST(ychar,1,0,icomm,ierr)
			if(ychar=='y' .or. ychar=='Y')then
				combinedlist=.true.
				if(iproc==0)print*,' Combining initial/final levels into a single list'
			else
				combinedlist = .false.
				if(iproc==0)print*,' Keeping initial/final levels in separate lists'
				
			end if
			call getoutputfile
			call boss1bjumps
			pndensities=.true.
			if(combinedlist)then
			    call density1b_bag_boss(.false.,.false.,samebasis)					
			else
			    call density1b_boss(.false.,samebasis)		
			end if
			
			case('jraise','JRAISE') !***************************************************
			
			if(Jz_f /= Jz_i + 2)then
				print*,'Wrong Jz values for bases to raise; try again ',Jz_i,Jz_f
				success = .false.
				cycle
			end if
			success =.true.
			writeout=.true.
			
			write_wfn= .true.

			if(iproc==0)then
				print*,' You are applying J+ raising operator to change Jz (M) and write out '
				print*,' (Will scale by appropriate factor to keep normalized)'
			end if
			
			call boss1bjumps
			call apply1b_boss(.true.,.false.)

			case('dn1jrm','DN1JRM') !***************************************************
			
			if(Jz_f /= Jz_i + 2)then
				print*,'Wrong Jz values for bases to raise; try again ',Jz_i,Jz_f
				success = .false.
				cycle
			end if
			
			if( (Jz_f/2)*2 /= Jz_f)then
				print*,' Only need to do this for even number of particles '
				success = .false.
				cycle
			end if
			success =.true.
			writeout=.true.
			
			write_wfn= .true.

			if(iproc==0)then
				print*,' You are computing missing 1-body density-matrix elements '
				print*,' by applying J+ raising operator to change Jz (M) '
				print*,' (raised wave functions written temporarily to TEMP.wfn)'
			end if
			
			call boss1bjumps
			call apply1b_boss(.true.,.true.)	
			pndensities=.true.
			call getoutputfile
			
			call density1b_bag_boss(.true.,.true.,.false.)			
			
			case('dn1mis','DN1MIS') !***************************************************
			
			if(Jz_f /= Jz_i + 2)then
				print*,'Wrong Jz values for bases to raise; try again ',Jz_i,Jz_f
				success = .false.
				cycle
			end if
			
			if( (Jz_f/2)*2 /= Jz_f)then
				print*,' Only need to do this for even number of particles '
				success = .false.
				cycle
			end if
			
			if(iproc==0) then
				print*,' You are calculating non-charge-changing one-body density matrices, proton-neutron formalism'	
				print*,' Do you want to combine initial/final levels into a single list? (y/n)'
				read(5,'(a)')ychar
            endif
			call BMPI_BCAST(ychar,1,0,icomm,ierr)
			if(ychar=='y' .or. ychar=='Y')then
				combinedlist=.true.
				if(iproc==0)print*,' Combining initial/final levels into a single list'
			else
				combinedlist = .false.
				if(iproc==0)print*,' Keeping initial/final levels in separate lists'
				
			end if
			
			success =.true.
			writeout=.true.
			
			write_wfn= .true.

			if(iproc==0)then
				print*,' You are computing missing 1-body density-matrix elements '
				print*,' Both initial and final wave functions must already be computed;'
				print*,' They must have M values differing by 1'
			end if
			
			call boss1bjumps
			pndensities=.true.
			call getoutputfile
			
			if(combinedlist)then
			    call density1b_bag_boss(.true.,.false.,.false.)					
			else
			    call density1b_boss(.true.,.false.)		
			end if
			
			case('lincom','LINCOM')	!***************************************************
			if(iproc==0) then		
				print*,' You are constructing a linear combinations of wave functions'
			endif
			success = .true.
			writeout=.true.
			write_wfn= .true.
			
			call lincom_boss

			case('ent1bd','ENT1BD') !***************************************************
			success = .true.
			if(iproc==0) then
			print*,' You are calculating one-body entropy (from m-scheme one-body density matrix )'	
			endif
            call getoutputfile
			call boss1bjumps
			pndensities=.true.
			call entropy1b_boss		
			
			case('spamp1','SPAMP1') !***************************************************
			success = .true.
			if(iproc==0) then
				print*,' You are calculating spectroscopic amplitudes with changes in the number of particles'	
            endif
			call getoutputfile
			call bosshops(.true.,.true.)  ! needs to be fixed
!			pndensities=.true.
			call spfac_boss		

			case('appsp1','APPSP1') !***************************************************
			success =.true.
			writeout=.true.
			
			write_wfn= .true.
			if(iproc==0) then
				print*,' You are adding/removing a (linear combination) of single-particles '	
            endif
!			call getoutputfile
			call bosshops(.true.,.true.)  
!			pndensities=.true.
			call appsp_boss	
			
			!
			case('allsp1','ALLSP1') !***************************************************
			success =.true.
			writeout=.true.
			
			write_wfn= .true.
			if(iproc==0) then
				print*,' THIS OPTION CURRENTLY HAS A BUG AND NEEDS TO BE FIXED'
				print*,' You are adding/removing a single-particle from a fixed initial wfn, '
				print*,' but from a designated set of orbitals (multiple output wfns)'	
            endif
			call getoutputfile
			call bosshops(.true.,.true.)  
!			pndensities=.true.
			call allsp_boss				
				
			

			case( 'app1bd','APP1BD') !***************************************************
			success =.true.
			writeout=.true.
			
			write_wfn= .true.
			
			if(iproc==0)then 
				print*,' You are applying a non-charge-changing 1-body operator and writing out '
			end if
  !          call getoutputfile
			call boss1bjumps
			call apply1b_boss(.false.,.false.)
			
			
			case('projct','PROJCT') !***************************************************
			if(iproc==0) then
				print*,' You are projecting a state into a basis with different truncation '
			end if
			success = .true.	
			writeout=.true.
			write_wfn= .true.
			

			call proj_boss
			
			case ('dn1bcc') !***************************************************
			success=.true.
			if(iproc==0) then
				print*,' You are computing the charge-changing, one-body densities '
				print*,' in p-n formalism '
			endif
			pndensities=.true.
			
			call getoutputfile
			call bosshops(.false.,.false.) ! needs to be fixed
			
			if(iproc==0) then
				print*,' all hops set up, now computing densities...'
			endif
			call ccdensity1b_boss
			
			case ('d1iscc') !***************************************************
			success=.true.
			if(iproc==0) then
				print*,' You are computing the charge-changing, one-body densities '
				print*,' in iso formalism '
			endif

			call getoutputfile
			call bosshops(.false.,.false.)
			pndensities=.false.
			if(iproc==0) then
			 print*,' all hops set up, now computing densities...'
			endif
			call ccdensity1b_boss
			
			
			case ('app1cc') !***************************************************
			success =.true.
			writeout=.true.
			
			write_wfn= .true.
			
			if(iproc==0)then 
				print*,' You are applying a charge-changing 1-body operator and writing out '
			end if
  
            call bosshops(.false.,.false.) ! needs to be fixed
  			call apply1b_cc_boss			
			
			
		
			case default
			if(iproc==0) then 
				 print*,' That choice not currently allowed '
			endif
			success = .false.
			
	   	    end select
			
		end do ! while not success

		return
		
	end subroutine main_menu
	
!
!=====================================================================
!
!  MAIN ROUTINE TO OPEN OUTPUT FILE
!
subroutine getoutputfile
  use menu_choices
  use program_info
  use io
  use nodeinfo
  use wfn_mod
  use system_parameters
  implicit none

  integer(4) :: ierr
  integer(4) :: ilast
  character(len=8) :: res_suf

1 continue
  if(iproc==0)then
     if(auto_input)then
        read(autoinputfile,'(a)')outfile
        write(6,*)outfile
     else 
        if(densityflag .or. trdensout)then
           print*,' Enter output name (required for your chosen option)!'
        else
           print*,' Enter output name (enter "none" if none)'
        endif
        read(5,'(a)')outfile
        write(autoinputfile,'(a)')outfile
	end if
	
  end if
  call BMPI_BARRIER(icomm,ierr)
  call BMPI_BCAST(outfile,20,0,icomm,ierr)
  ilast = index(outfile,' ')-1  
  if(densityflag .and. outfile(1:ilast)=='none')goto 1
  
  if(outfile(1:ilast)=='none')then
     writeout = .false.
  else
     writeout = .true.
     if ( iproc == 0 ) then
        res_suf = ".res"
        if(menu_choice == 'dn1xpn') res_suf = ".dres"
        if(menu_choice == 'dn1iso') res_suf = ".dres"
        if(menu_choice == 'dn1mis') res_suf = ".dres"
        if(menu_choice == 'dn1jrm') res_suf = ".dres"
		
		if(menu_choice=='spamp1')res_suf=".spres"
        open(unit=resultfile,file=outfile(1:ilast)// TRIM(res_suf),status = 'unknown')
        write(resultfile,'(a18,a6,a9)')"  RHODIUM version ",version,lastmodified
		
		if(menu_choice =='dn1xpn' .or. menu_choice == 'dn1iso' .or. menu_choice == 'dn1mis' .or. menu_choice == 'dn1jrm')then
   		    write(resultfile,'(" ")')
   		    write(resultfile,*)np_i(1),np_f(2)
		end if
		
   	    select case(menu_char)
       	 case('dx')
   		 write(6,*)' Density matrices written to :', outfile(1:ilast)//".dres" 
   		 write(logfile,*)' Density matrices written to :', outfile(1:ilast)//".dres" 

       	 case('t')
   		 write(6,*)' Output written to :', outfile(1:ilast)//".res" 
   		 write(6,*)' Wfn info written to :', outfile(1:ilast)//".trwfn" 
   		 write(logfile,*)' Output written to :', outfile(1:ilast)//".res" 
   		 write(logfile,*)' Wfn info written to :', outfile(1:ilast)//".trwfn" 
		 
   		 case('d')
   		 write(6,*)' Output written to :', outfile(1:ilast)//".res" 
   		 write(logfile,*)' Output written to :', outfile(1:ilast)//".res" 
		 
		 
   		 case default

        end select
		
     endif

     if(base_scratch_dir /= '.') then
        scratch_dir = TRIM(base_scratch_dir) // '/' // outfile(1:ilast)
        if(iproc == 0) print *, "set scratch_dir=", TRIM(scratch_dir)
        call system('mkdir -p ' // TRIM(scratch_dir))
     endif
  endif
  return
end subroutine getoutputfile

!===========================================================================	
	

end module menu_mod
