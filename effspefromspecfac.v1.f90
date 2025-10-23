!
! code for computing effective s.p.e.s from spectroscopic factors
! 

module filestuff
	implicit none
	
	integer :: unitin = 57
	integer :: unitout= 58
	
contains
	
	subroutine open_inputfile
		implicit none
		
		character*80 :: filenamein
		integer :: ilast
		
111     continue		
		print*,' Enter name of .spres file (do not include extension)'
		read(5,'(a)')filenamein
		ilast = index(filenamein,' ')-1
		open(unit=unitin,file=filenamein(1:ilast)//".spres",status = 'old',err = 112)
		
		return
112		continue
        print*,' That file does not exist '
		go to 111

	end subroutine open_inputfile	
	
end module filestuff

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module sporbitals
	implicit none
	
	integer :: numorb
	integer, allocatable :: jorb(:),lorb(:),nrorb(:)
	
contains
	
	subroutine readinorbitals(iunit)
		implicit none
		integer,intent(in) :: iunit
		character*80 :: line
		logical :: found
		integer :: iorb
		integer :: n,i1,i2,i3
		found = .false.
		do while (.not. found)
			read(iunit,'(a)')line
			if(line(1:3)=="ORB")found = .true.
		end do ! while not found
		
		numorb = 0
		do iorb = 1,10000
			read(iunit,*,err=101)n,i1,i2,i3
			numorb = numorb +1
			if(numorb /=n)then
				print*,' Some problem reading in ',numorb,iorb,n
				stop
			end if
		end do
101     continue		
		print*,numorb,' orbitals '
		rewind(iunit)
		found =.false.
		do while (.not. found)
			read(iunit,'(a)')line
			if(line(1:3)=="ORB")found = .true.
		end do ! while not found
		allocate(nrorb(numorb),jorb(numorb),lorb(numorb))
		
		do iorb = 1,numorb
			read(iunit,*)n,i1,i2,i3
			if (iorb /= n)stop 9999
			nrorb(iorb)= i1
			lorb(iorb) = i2
			jorb(iorb) = i3
		end do		
		
		return
		
	end subroutine readinorbitals
	
end module sporbitals
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module levels
	implicit none
	
	integer,target :: nlevelsi,nlevelsf
	real,allocatable,target :: energyi(:),energyf(:),xji(:),xjf(:)
	
contains
	
	subroutine readlevels(iunit,whichlevel)
		
		implicit none
		integer,intent(in) :: iunit
		character,intent(in) :: whichlevel ! can be 'i' for initial or 'f' for final
		
		real,pointer :: energyx(:),xjx(:)
		character*3 :: target_phrase,test_phrase
		logical :: success
		integer,pointer :: nlevelsx
		integer :: ilevel, n
		real :: e,ex,xj
		
		select case (whichlevel)
		
		case ('i')
		
		target_phrase = 'ini'
		nlevelsx => nlevelsi
		
		case ('f')
		
		target_phrase = 'fin'
		nlevelsx => nlevelsf
		
		case default
		
		print*,' error in subroutine readlevesl ',whichlevel
		stop
		
	    end select
		
		success = .false.
		
		do while(.not.success)
			read(iunit,'(19x,a3)')test_phrase
			if(test_phrase==target_phrase)then
				success = .true.
			end if
		end do ! while not success
! count up # of levels		
		do ilevel = 1,100000
			read(iunit,*,err=333)n,e,ex,xj
			nlevelsx = n
			if(n /= ilevel)then
				print*,' miscount ',n,ilevel
				stop
			end if
		end do
		
333     continue		
        rewind(iunit)
		success = .false.
		
		do while(.not.success)
			read(iunit,'(19x,a3)')test_phrase
			if(test_phrase==target_phrase)then
				success = .true.
			end if
		end do ! while not success
		select case (whichlevel)
		
		case ('i')		
		
		allocate(energyi(nlevelsi),xji(nlevelsi))
		print*,' There are ',nlevelsx,' initial levels '
		energyx => energyi
		xjx     => xji
		
		case ('f')
		
		allocate(energyf(nlevelsf),xjf(nlevelsf))
		print*,' There are ',nlevelsx,' final levels '
		energyx => energyf
		xjx     => xjf		
		
     	end select
		
		do ilevel = 1,nlevelsx
			read(iunit,*)n,e,ex,xj
			energyx(n)= e
			xjx(n)    =xj
		end do
		
		return
	end subroutine readlevels
	
end module levels
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module spamps
	use levels
	use sporbitals
	implicit none
	
	integer :: itarget ! which target level to choose for the initial state
	real, allocatable :: spamp(:,:)  ! spamp(orbital,final_level)
	
	
contains
	
	subroutine get_spamps(iunit)
		implicit none
		
		integer,intent(in) :: iunit 
		integer :: flevel
		character*80 :: line
		logical :: success
		integer :: i,f,iorb,ii
		real :: x
		
		
		if(nlevelsi == 1)then
			itarget = 1
			print*,' Only one initial level; using that'
		else
1           continue			
			print*,' Which initial level to use? Up to  ',nlevelsi
			read(5,*)itarget
			if(itarget > nlevelsi)goto 1
		end if
		allocate(spamp(numorb,nlevelsf))
		spamp = 0.0
		backspace(iunit)
		success = .false.
	    do while (.not.success)
		    read(iunit,'(a)')line
		    if(line(2:4)=='Ini')then
			    backspace(iunit)
			    read(iunit,'(18x,i4)')i
			    print*,i
			    if(i==itarget)success=.true.
		    end if
	    end do
		print*,nlevelsf
		do f = 1,nlevelsf
			success = .false.
		    do while (.not.success)
			    read(iunit,'(a)')line
!				if(f > 997)write(17,*)line
!				print*,line
			    if(line(2:4)=='Fin')then
				    backspace(iunit)
				    read(iunit,'(17x,i5)')i
				    if(i==f)success=.true.
			    end if
		    end do
			read(iunit,'(a)')line
			do iorb = 1,numorb
				read(iunit,*,err=99)ii,x
				spamp(ii,f)=x
!			    print*,f,ii,x	
			end do
99          continue
            backspace(iunit)
			
	    end do
		close(iunit)
		
		open(unit=77,file="spfac.dat",status='unknown')
		do f = 1,nlevelsf
			write(77,*)energyf(f)-energyi(itarget),(spamp(iorb,f)**2,iorb = 1,numorb)
		end do
		close(77)
		
		
	end subroutine get_spamps
	
	
	subroutine spes
		implicit none
		real, allocatable :: sum(:),esum(:),spe(:)
		
		integer :: ilevel,iorb
		
		allocate(sum(numorb),esum(numorb),spe(numorb))
		sum = 0.0
		esum = 0.0
		spe = 0.0
		
		do ilevel = 1,nlevelsf
			do iorb = 1,numorb
				sum(iorb)=sum(iorb)+ spamp(iorb,ilevel)**2
				esum(iorb)=esum(iorb)+ spamp(iorb,ilevel)**2*(energyf(ilevel)-energyi(itarget))
				
			end do
			
		end do
		print*,' sums ',sum
		do iorb = 1,numorb
			print*,iorb,esum(iorb)/sum(iorb)
		end do
		return
	end subroutine spes
	
end module spamps
! main program


use levels
use sporbitals
use filestuff
use spamps


call open_inputfile
call readlevels(unitin,'i')
call readlevels(unitin,'f')
call readinorbitals(unitin)
call get_spamps(unitin)
call spes

end