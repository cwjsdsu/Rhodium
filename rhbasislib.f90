!...........................................................................
!
!  basis dimension
!  and information on how to find basis label from proton, neutron labels
!
!  if proton Slater determinant is ip and neutron Slater determinant is in
!  basis = pstart(ip)+nstart(in)
!
!  IMPORTANT:  ip and in must come from "conjugate" sectors
!  also note: neutrons form inner loop
!
module basis
  use precisions

  implicit none
  integer,target  :: Nword_i(2),Nword_f(2)     ! # of words needed to construct a Slater determinant
  
  integer(kind=basis_prec),target :: dimbasis_i,dimbasis_f 		    ! dimension of combined basis
  integer(kind=basis_prec),target :: dimbasis_proc_i,dimbasis_proc_f         ! local dimension of combined basis
  integer(kind=basis_prec),target :: nXsd_i(2),nXsd_f(2)         ! dimensions of proton, neutrons SDs
  integer(kind=basis_prec),allocatable,target :: pstart_i(:),nstart_i(:),pstart_f(:),nstart_f(:)   ! used for finding basis
!---- STORAGE OF BASIS; TWO VERSIONS-----
!     VERSION 1: LIST OF OCCUPIED STATES  
  integer(4),allocatable,target :: pocc_i(:,:),nocc_i(:,:),pocc_f(:,:),nocc_f(:,:)
  logical    :: storebasisasocc ! if TRUE, store as occupation
!     VERSION 2: STORED AS BITS
  integer(8),allocatable,target :: psdlist_i(:,:),nsdlist_i(:,:),psdlist_f(:,:),nsdlist_f(:,:)

!----------------------------------------  
contains
	
!
      function bit_occ(hsd,isps)
      use bitstuff
      implicit none
      logical bit_occ
      integer,pointer :: hsd(:)
      integer isps

      integer iword,i
      integer imask
      
      iword = (isps-1)/(max_bit_word) + 1
      i = isps - (iword -1)*(max_bit_word)
      imask = ishft(1,i-1)
	  bit_occ = .not.(iand(imask,hsd(iword)) == 0)   ! should be faster than an if statement

      return
      end function bit_occ

!
!
!
!  converts a   (binary word) to an array ASSUMING fixed numpart
!  NEW VERSION March 2013 to speed up usage for certain cases
!  reduces division which is slow
!  note this has a different definition of occ from convert_haiku_array
!
!  INPUT:
!   ith = species/handness
!   nh = # of paricles
!   hsd(:) = pointer to binary words that make up the haiku
!
!  OUTPUT:
!    occ(:)  array = 0. or 1 designating occupied s.p. states
!
!
      subroutine convert_bitrepsd2occ(it,whichbasis,numpart,hsd,occ)
      use bitstuff
	  use spstate
      implicit none


      integer :: it
	  character(3) :: whichbasis
      integer :: numpart
      integer(8) :: hsd(:)
      integer :: occ(:)

      integer isps
      integer n
      integer i, iword, imask,ilast,jhsd
	  integer,pointer :: nrhspsx(:),Nwordx(:)
	  
	  if(numpart==0)return
	  
	  if(whichbasis=='INI')then
		  nrhspsx=>nrhsps_i
		  Nwordx => Nword_i
	  else
		  nrhspsx=>nrhsps_f
		  Nwordx => Nword_f
	  end if
	  
      occ = 0
	  if(numpart==0)return
      n = 0

      isps = 0
      do iword = 1,Nwordx(it )
         jhsd=hsd(iword)
!............. SET UP SIZE OF INNER DO LOOP .....
         if(iword == Nwordx( it ))then
            ilast = nrhspsx(it)-(Nwordx(it) -1)*(max_bit_word)
         else
            ilast = max_bit_word
         end if
         imask = 1
         do i = 1,ilast
            isps = isps + 1
!			print*,jhsd,imask,iand(imask,jhsd)
            if(iand(imask,jhsd) /= 0)then
               n = n + 1
               occ(n) = isps
            endif
            imask = ishft(imask,1)  ! shift by 1
         end do
         if(n==numpart)return

      end do ! iword

      if(n < numpart)then
        print*,' problem converting ',n,numpart
		print*,it,whichbasis,Nwordx(it)
		print*,hsd(1:Nwordx(it))
 !       stop
      endif
      return
      end subroutine convert_bitrepsd2occ

!
!
!  check if a bit is occupied or unoccupied 
! 
!  note this has a different definition of occ from convert_haiku_array
!
!  INPUT:
!   it = species/handness
!   hsd(:) = pointer to binary words that make up the haiku
!   isps = occupied state
!
!  OUTPUT:
!    if TRUE then occupied
!
      logical function check_bitocc(whichbasis,hsd,isps)
      use bitstuff
	  use spstate
      implicit none

	  character(3) :: whichbasis
      integer :: numpart
      integer(8) :: hsd(:)
      integer :: isps,iword
	  integer(8) :: imask,ibit

      integer n
!      integer i, iword, imask,ilast,jhsd
	  integer,pointer :: nrhspsx(:)
	  
	  if(isps==0)then
		  check_bitocc=.false.
		  return
	  end if
	  	  
	  if(whichbasis=='INI')then
		  nrhspsx=>nrhsps_i
	  else
		  nrhspsx=>nrhsps_f
	  end if

	  iword = 1
	  ibit = isps-(iword-1)*max_bit_word
	  do while(ibit > max_bit_word)
		  ibit =ibit - max_bit_word
		  iword =iword+1
	  end do

	  imask = 1
	  imask = ishft(imask,ibit-1)

      check_bitocc = (iand(imask,hsd(iword))/=0)

      return
      end function check_bitocc	  	
!
!  converts a   (binary word) to an array ASSUMING fixed numpart
!  NEW VERSION March 2013 to speed up usage for certain cases
!  reduces division which is slow
!  note this has a different definition of occ from convert_haiku_array
!
!  INPUT:
!   ith = species/handness
!   nh = # of paricles
!    occ(:)  array = 0. or 1 designating occupied s.p. states
!
!  OUTPUT:
!   hsd(:) = pointer to binary words that make up the haiku
!
!
      subroutine convert_occ2bitrepsd(it,whichbasis,numpart,occ,hsd)
      use bitstuff
	  use spstate
      implicit none


      integer :: it
	  character(3) :: whichbasis
      integer :: numpart
      integer(8) :: hsd(:)
      integer :: occ(:)

      integer isps
      integer n
      integer i, iword, ilast,jhsd
	  integer(8) :: imask,ibit
	  integer,pointer :: nrhspsx(:),Nwordx(:)
	  
	  
	  if(whichbasis=='INI')then
		  nrhspsx=>nrhsps_i
		  Nwordx => Nword_i
	  else
		  nrhspsx=>nrhsps_f
		  Nwordx => Nword_f
	  end if
	  
      hsd =0
	  if(numpart==0)return
      n = 0

      isps = 0
	  iword = 1
	  do i = 1,numpart
		  ibit = occ(i)-(iword-1)*max_bit_word
	!	  if(i==5)print*,ibit,occ(i), iword
		  do while(ibit > max_bit_word)
			  ibit =ibit - max_bit_word
			  iword =iword+1
		  end do
		  if(ibit==0)then
			  print*,' some problem with ibit ',i,occ(i),iword,ibit
			  print*,whichbasis,it,numpart,occ(:),iword
			  stop
		  end if
		  if(iword > Nwordx(it))then
			  print*,' some trouble with iword ',i,occ(i),iword,ibit,Nwordx(it)
			  print*,nrhspsx(it)
			  stop
		  end if
		  imask = 1
		  imask = ishft(imask,ibit-1)
		  if(imask < 1)then
			  print*,it,' some problem with mask ',whichbasis,occ
			  print*,i,numpart,ibit,iword,imask
			  print*,bit_size(imask)
			  print*,i,occ(i),occ(i)-(iword-1)*max_bit_word
			  stop
		  end if
		  hsd(iword)=ior(imask,hsd(iword))
	  end do

      return
      end subroutine convert_occ2bitrepsd
	
	  
end module basis