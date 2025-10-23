!
!
!

module basisfilehandler
	use io
	use nodeinfo
	use bmpi_mod
implicit none

contains

!----- ROUTINE openbasisfile -----
!      opens .bas files created by BIGSTICK 
!  INPUT:
!    character(3) : whichfile = INI or FIN
!

subroutine openbasisfile(whichfile)

implicit none

!....INPUT VARIABLES....
character(3) :: whichfile
!....INTERNAL VARIABLES...
integer ilast
integer :: ierr
logical :: success

if(iproc==0)then

select case (whichfile)

   case('INI')
   
    success=.false.
    do while(.not.success)
        write(6,*)' Enter name of INITIAL .bas file '
        read(5,'(a)')basfilename_i
        ilast = index(basfilename_i,' ')-1
        inquire(file=basfilename_i(1:ilast)//".bas",exist=success)
        
        if(.not.success)then
            write(6,*)' file ',basfilename_i(1:ilast),'.bas does not exist '
        end if
    end do ! while
    open(unit=basfile_i,file=basfilename_i(1:ilast)//".bas",status='old',access='stream',form='unformatted',err=1111)   
   
   case('FIN')
   
    success=.false.
    do while(.not.success)
        write(6,*)' Enter name of FINAL .bas file '
        read(5,'(a)')basfilename_f
        ilast = index(basfilename_f,' ')-1
        inquire(file=basfilename_f(1:ilast)//".bas",exist=success)
        
        if(.not.success)then
            write(6,*)' file ',basfilename_f(1:ilast),'.bas does not exist '
        end if
    end do ! while
    open(unit=basfile_f,file=basfilename_f(1:ilast)//".bas",status='old',access='stream',  form='unformatted',err=1111)
    
   case default
   
   write(6,*)' Wrong input to routine openbasisfile ',whichfile
   stop

end select

write(6,*)' Basis file successfully opened '
end if

call BMPI_Barrier(icomm,ierr)
return
1111 continue
       print*,' File exists but did not open successfully '
	   call BMPI_Abort(icomm,101,ierr)
       stop
end subroutine openbasisfile


!=======================================================
!
! routine for reading in basis file
!
!  based upon routine basis_out4postprocessing in BIGSTICK file boutputlib.f90
!

subroutine readbasisfile(whichfile)
		
	use precisions
	use spstate
	use sporbit
	use io
	use sectors
	use basis
	use system_parameters
	use wfn_mod
	use spmap_mod

	implicit none

!....INPUT VARIABLES....
	character(3) :: whichfile
!.....POINTERS.........

   integer,pointer :: nhsps(:)	
   integer, pointer :: basisfile
   type (spst),pointer :: hspsqn(:,:)
   integer,pointer:: nsectors(:)
   integer(8),pointer :: nXsd(:)
   type(mastersect),pointer :: xsd(:)	  
   integer(4),allocatable :: occ(:) , occtest(:)
   integer(4),pointer :: xocc(:,:),np(:)

!...... OTHER INTERNAL VARIABLES....

   integer :: ith   ,i,isd,ps,ns,xs
!   integer :: jzZ,parX,wX
   integer :: ncsectors
!   integer(8) :: xsdstart,xsdend,nxsd
   integer :: nstate
   integer(kind=basis_prec)::ip
   integer(kind=basis_prec) :: px,nx

   integer(kind=basis_prec),pointer :: pstart(:),nstart(:)
   integer(8), pointer :: xsdlist(:,:)
   integer :: occoffset
   integer :: hsps_size
   integer :: maxWtotX
   logical :: sameparity
   integer :: ierr
	
   select case (whichfile)
   
   case('INI')
   
   nhsps=> nhsps_i
   basisfile=>basfile_i
   nsectors=> nsectors_i
   nXsd=> nXsd_i
   np=>np_i     ! CHECK ON P-H CONJUGATION
   if(iproc==0)rewind(basfile_i)
   call read_wfn_header_select(basfile_i,.true.,.false.,'I')
   
   sameparity = allsameparity_i
   case('FIN')
   
   nhsps => nhsps_f
   basisfile=>basfile_f
   nsectors=> nsectors_f
   nXsd=> nXsd_f 
   np=>np_f     ! CHECK ON P-H CONJUGATION
   call read_wfn_header_select(basfile_f,.true.,.false.,'F')
   sameparity = allsameparity_f
   
   
   case default
   
   if(iproc==0)print*,' some error in which file ',whichfile
   call BMPI_abort(icomm,101,ierr)
   
   end select
   call BMPI_Barrier(icomm,ierr)
   if(iproc==0)then
     print*,' '
     print*,' Header read'
     print*,' '
   end if
   nstate = 0

!...... ALLOCATE SPACE FOR HSPS STATES....



!............. READ IN HSPS STATES............

ith = -1
if(iproc==0)read(basisfile)nhsps
call BMPI_Bcast(nhsps,5,0,icomm,ierr)
hsps_size = max(nhsps(-1),nhsps(-2))
hsps_size = max(nhsps(1),hsps_size)
hsps_size = max(nhsps(2),hsps_size)

if(whichfile=='INI')then
   allocate(hspsqn_i(-2:2,hsps_size))
   hspsqn=> hspsqn_i
else
   allocate(hspsqn_f(-2:2,hsps_size))
   hspsqn=> hspsqn_f	
end if   

do i = 1,nhsps(ith)  ! left states
  nstate = nstate + 1
  
  if(iproc==0)read(basisfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
  hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par
  call BMPI_Bcast(nstate,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%nr,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%l,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%j,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%m,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%w,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%par,1,0,icomm,ierr)

enddo  ! i

ith = 1
if(nhsps(ith)> hsps_size)then
	if(iproc==0)print*,' bad allocation on # of h.sps states ',ith,nhsps(ith)
	call BMPI_Abort(icomm,111,ierr)
	stop
end if

do i = 1,nhsps(ith)  ! left states
  nstate = nstate + 1
  
  if(iproc==0)read(basisfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
  hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par
  
  call BMPI_Bcast(nstate,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%nr,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%l,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%j,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%m,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%w,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%par,1,0,icomm,ierr)
enddo  ! i	
ith = -2

if(nhsps(ith)> hsps_size)then
	if(iproc==0)print*,' bad allocation on # of h.sps states ',ith,nhsps(ith),hsps_size
    call BMPI_Abort(icomm,101,ierr)
	stop
end if
do i = 1,nhsps(ith)  ! left states
  nstate = nstate + 1
  if(iproc==0)read(basisfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
  hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par
  call BMPI_Bcast(nstate,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%nr,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%l,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%j,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%m,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%w,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%par,1,0,icomm,ierr)
enddo  ! i
ith = 2
if(nhsps(ith)> hsps_size)then
	if(iproc==0)print*,' bad allocation on # of h.sps states ',ith,nhsps(ith),hsps_size
    call BMPI_Abort(icomm,101,ierr)
	stop
end if
do i = 1,nhsps(ith)  ! left states
  nstate = nstate + 1
  if(iproc==0)read(basisfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
  hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par
  call BMPI_Bcast(nstate,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%nr,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%l,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%j,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%m,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%w,1,0,icomm,ierr)
  call BMPI_Bcast(hspsqn(ith,i)%par,1,0,icomm,ierr)
enddo  ! i		

call BMPI_Barrier(icomm,ierr)
!......... NOW THAT S.P. INFORMATION READ IN, TAKE CARE OF SOME BOOK KEEPING

call setup_rhsps(whichfile)

!..... READ IN PROTON SECTORS, SDs..................................
if(iproc==0)then
	read(basisfile)nsectors(1),nxSD(1)
	print*,nsectors(1),' proton sectors ',nxsd(1),' SDs '
end if
call BMPI_Bcast(nsectors(1),1,0,icomm,ierr)
call BMPI_Bcast(nxSD(1),1,0,icomm,ierr)

if(whichfile=='INI')then
	
	allocate(xsd_i(1)%sector(nsectors(1)))
	xsd=>xsd_i
	if(storebasisasocc)then
		allocate(pocc_i(nxsd(1),np(1)))
		xocc => pocc_i
	else
		if(nrhsps_i(1)> 0)then
		   allocate(psdlist_i(nxsd(1),nrhsps_i(1)))
	   else
		   allocate(psdlist_i(nxsd(1),1))
	   end if
		xsdlist => psdlist_i
	end if
	
	allocate(pstart_i(nxsd(1)))
	pstart => pstart_i
	
else
	allocate(xsd_f(1)%sector(nsectors(1)))
	xsd=>xsd_f
	if(storebasisasocc)then
		allocate(pocc_f(nxsd(1),np(1)))
		xocc => pocc_f
	else
		if(nrhsps_f(1)> 0)then
		
		    allocate(psdlist_f(nxsd(1),nrhsps_f(1)))
		else
		    allocate(psdlist_f(nxsd(1),1))

		end if
		xsdlist => psdlist_f
	end if
	allocate(pstart_f(nxsd(1)))
	pstart => pstart_f
end if
xsdlist = 0

if(np(1)>0)then
	allocate(occ(np(1)))
	allocate(occtest(np(1)))
endif

do ps = 1,nsectors(1)
	if(iproc==0)then
      read(basisfile)xs,xsd(1)%sector(ps)%jzX,xsd(1)%sector(ps)%parX,xsd(1)%sector(ps)%Wx
	  read(basisfile)xsd(1)%sector(ps)%xsdstart,xsd(1)%sector(ps)%xsdend,xsd(1)%sector(ps)%nxsd
	  read(basisfile)xsd(1)%sector(ps)%ncsectors
    end if
	call BMPI_Bcast(xs,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(1)%sector(ps)%jzX,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(1)%sector(ps)%parX,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(1)%sector(ps)%Wx,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(1)%sector(ps)%xsdstart,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(1)%sector(ps)%xsdend,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(1)%sector(ps)%nxsd,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(1)%sector(ps)%ncsectors,1,0,icomm,ierr)
	
	if(whichfile=='INI')then
		allocate(xsd_i(1)%sector(ps)%csector(xsd(1)%sector(ps)%ncsectors ))
	else
		allocate(xsd_f(1)%sector(ps)%csector(xsd(1)%sector(ps)%ncsectors ))	
		
	end if
	
	if(iproc==0)then
	  read(basisfile)(xsd(1)%sector(ps)%csector(i),i=1,xsd(1)%sector(ps)%ncsectors)
	  read(basisfile)xsd(1)%sector(ps)%basisstart,xsd(1)%sector(ps)%basisend	
    end if
	call BMPI_Bcast(xsd(1)%sector(ps)%csector,xsd(1)%sector(ps)%ncsectors,0,icomm,ierr)
	call BMPI_Bcast(xsd(1)%sector(ps)%basisstart,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(1)%sector(ps)%basisend,1,0,icomm,ierr)
	
!................ NOW READ IN SLATER DETERMINANTS..............

    do isd = 1,xsd(1)%sector(ps)%nxsd
		if(np(1)> 0)then
			
		   if(iproc==0)read(basisfile)ip,px,(occ(i),i=1,np(1))
		   call BMPI_Bcast(ip,1,0,icomm,ierr)
		   call BMPI_Bcast(px,1,0,icomm,ierr)
		   call BMPI_Bcast(occ,np(1),0,icomm,ierr)
		   
	   else
		   if(iproc==0)read(basisfile)ip,px
		   call BMPI_Bcast(ip,1,0,icomm,ierr)
		   call BMPI_Bcast(px,1,0,icomm,ierr)
	   end if
		pstart(ip)=px
		if(storebasisasocc)then
			do i = 1,np(1)
				xocc(ip,i)=occ(i)  ! check if this is best order
			end do
		else			
			call convert_occ2bitrepsd(1,whichfile,np(1),occ,xsdlist(ip,:))						
			call convert_bitrepsd2occ(1,whichfile,np(1),xsdlist(ip,:),occtest)						
			do i = 1,np(1)
				if(occ(i)/=occtest(i))then
					if(iproc==0)then
					  print*,'original ',occ
					  print*,'conversion', occtest
				    end if
					call BMPI_Abort(icomm,101,ierr)
					stop
				end if
			end do
			
		end if		
	end do		
	
end do ! ps
if(np(1)>0)then
	deallocate(occ)
	deallocate(occtest)
end if

call BMPI_Barrier(icomm,ierr)
!..... READ IN NEUTRON SECTORS, SDs..................................
if(iproc==0)then
  read(basisfile)nsectors(2),nxSD(2)
  print*,nsectors(2),' neutron sectors ',nxsd(2),' SDs '
end if
call BMPI_Bcast(nsectors(2),1,0,icomm,ierr)
call BMPI_Bcast(nxSD(2),1,0,icomm,ierr)

if(whichfile=='INI')then
	
	allocate(xsd_i(2)%sector(nsectors(2)))
	xsd=>xsd_i
	
	if(storebasisasocc)then
		allocate(nocc_i(nxsd(2),np(2)))
		xocc => nocc_i
	else
		if(nrhsps_i(2)> 0)then
		   allocate(nsdlist_i(nxsd(2),nrhsps_i(2)))
	    else
 		   allocate(nsdlist_i(nxsd(2),1))
		end if
		xsdlist => nsdlist_i
	end if
	
	allocate(nstart_i(nxsd(2)))
	nstart => nstart_i
	
else
	allocate(xsd_f(2)%sector(nsectors(2)))
	xsd=>xsd_f
	
	if(storebasisasocc)then
		allocate(nocc_f(nxsd(2),np(2)))
		xocc => nocc_f
	else
		if(nrhsps_f(2)> 0)then
		
		   allocate(nsdlist_f(nxsd(2),nrhsps_f(2)))
	    else
 		   allocate(nsdlist_f(nxsd(2),1))
	    end if
		xsdlist => nsdlist_f
	end if

	allocate(nstart_f(nxsd(2)))
	nstart => nstart_f
end if
xsdlist = 0

if(np(2)>0)allocate(occ(np(2)),occtest(np(2)))

do ns = 1,nsectors(2)
	if(iproc==0)then
      read(basisfile)xs,xsd(2)%sector(ns)%jzX,xsd(2)%sector(ns)%parX,xsd(2)%sector(ns)%Wx
	  read(basisfile)xsd(2)%sector(ns)%xsdstart,xsd(2)%sector(ns)%xsdend,xsd(2)%sector(ns)%nxsd
	  read(basisfile)xsd(2)%sector(ns)%ncsectors
    end if
	call BMPI_Bcast(xs,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(2)%sector(ns)%jzX,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(2)%sector(ns)%parX,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(2)%sector(ns)%Wx,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(2)%sector(ns)%xsdstart,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(2)%sector(ns)%xsdend,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(2)%sector(ns)%nxsd,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(2)%sector(ns)%ncsectors,1,0,icomm,ierr)
	if(whichfile=='INI')then
		allocate(xsd_i(2)%sector(ns)%csector(xsd(2)%sector(ns)%ncsectors ))
	else
		allocate(xsd_f(2)%sector(ns)%csector(xsd(2)%sector(ns)%ncsectors ))	
		
	end if
	if(iproc==0)then
	  read(basisfile)(xsd(2)%sector(ns)%csector(i),i=1,xsd(2)%sector(ns)%ncsectors)
	  read(basisfile)xsd(2)%sector(ns)%basisstart,xsd(2)%sector(ns)%basisend	
    end if
	call BMPI_Bcast(xsd(2)%sector(ns)%csector,xsd(2)%sector(ns)%ncsectors,0,icomm,ierr)
	call BMPI_Bcast(xsd(2)%sector(ns)%basisstart,1,0,icomm,ierr)
	call BMPI_Bcast(xsd(2)%sector(ns)%basisend,1,0,icomm,ierr)
!................ NOW READ IN SLATER DETERMINANTS..............
    do isd = 1,xsd(2)%sector(ns)%nxsd
		if(np(2)> 0)then
		   if(iproc==0)read(basisfile)ip,nx,(occ(i),i=1,np(2))
		   call BMPI_Bcast(ip,1,0,icomm,ierr)
		   call BMPI_Bcast(nx,1,0,icomm,ierr)
		   call BMPI_Bcast(occ,np(2),0,icomm,ierr)
	   else
		   if(iproc==0)read(basisfile)ip,nx
		   call BMPI_Bcast(ip,1,0,icomm,ierr)
		   call BMPI_Bcast(nx,1,0,icomm,ierr)
	   end if
	   		   
!		print*,ip,nx,occ
		if(whichfile=='INI')then
			occoffset = nrhsps_i(1)
		else
			occoffset = nrhsps_f(1)			
		end if
!        occoffset = nhsps
		do i=1,np(2)
			occ(i)=occ(i)-occoffset
		end do 
		nstart(ip)=nx
		if(storebasisasocc)then
			do i = 1,np(2)
				xocc(ip,i)=occ(i)  ! check if this is best order
			end do
		else
!			print*,ip,occ
			
			call convert_occ2bitrepsd(2,whichfile,np(2),occ,xsdlist(ip,:))
			call convert_bitrepsd2occ(2,whichfile,np(2),xsdlist(ip,:),occtest)		
			do i = 1,np(2)
				if(occ(i)/=occtest(i))then
					if(iproc==0)then
					   print*,'original ',occ
					   print*,'conversion', occtest
				    end if
					call BMPI_Abort(icomm,101,ierr)
					stop
				end if
			end do

		end if
		
	end do	
	
end do ! ns
if(np(2)>0)then
	deallocate(occ)
	deallocate(occtest)
end if

if(iproc==0)then
  print*,' DONE READING IN '
  close(basisfile)
end if

call BMPI_Barrier(icomm,ierr)
return
end subroutine readbasisfile



end module basisfilehandler