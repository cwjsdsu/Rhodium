!
!  "GROUPS" are sets of single-particle states which have the same abelian quantum numbers, i.e.,
!  Jz,parity, and W (but not J and not explicitly L)
!

module groups
	use spstate
	implicit none
	
	
	
contains
!
!===========================================================================
!  Calls routines to divide hsp spaces into "groups", by W,jz,par
!===========================================================================
!
! CALLED BY:
! SUBROUTINES CALLED
!	group_hsp
!
subroutine grouper

  implicit none
  
  call group_sp(1,'INI')
  call group_sp(1,'FIN')
  call group_sp(2,'INI')
  call group_sp(2,'FIN')

  return
end subroutine grouper
      
!===========================================================================
!  Divides up hsp of species it by W, jz, par into groups
!
! called by: grouper
!
!===========================================================================
subroutine group_sp(it,whichbasis)

  implicit none
  
  integer :: it         ! species
  character(3) :: whichbasis

  integer :: i
  integer :: n
  type(grouplist), pointer :: groupnow
  integer ::  ngr 
  
  if(whichbasis=='INI')then
	  groupnow => group_i(it)
  else
	  groupnow => group_f(it)
  end if
  
!------------------COUNT UP # OF SECTIONS-----------------------------------
  ngr = 1
  do i = 2, nhsps(ith)
     if ( hspsqn(ith,i)%w /= hspsqn(ith,i-1)%w .or. hspsqn(ith,i)%m /= hspsqn(ith,i-1)%m & 
         .or. hspsqn(ith,i)%par /= hspsqn(ith,i-1)%par ) then
        ngroups(ith) = ngroups(ith) + 1
     end if
  end do ! it
!      print*,' there are ',nsections(it),' sections ',nhsps(it)
!------------------ALLOCATE-------------------------------------------------

  allocate ( group(ith)%jz(ngroups(ith)),  group(ith)%par(ngroups(ith)), group(ith)%w(ngroups(ith)) )
  allocate ( group(ith)%start(ngroups(ith)), group(ith)%fin(ngroups(ith)) )

!------------------Find START, FIN of each section--------------------------
  n = 1
  group(ith)%w(1) = hspsqn(ith,1)%w
  group(ith)%jz(1) = hspsqn(ith,1)%m
  group(ith)%par(1) = hspsqn(ith,1)%par
  hspsqn(ith,1)%group = n
  group(ith)%start(1) = 1
  do i = 2, nhsps(ith)
     if ( hspsqn(ith,i)%w > hspsqn(ith,i-1)%w .or.  hspsqn(ith,i)%m /= hspsqn(ith,i-1)%m & 
         .or. hspsqn(ith,i)%par /= hspsqn(ith,i-1)%par ) then
        n = n + 1
        group(ith)%w(n) = hspsqn(ith,i)%w
        group(ith)%par(n) = hspsqn(ith,i)%par
        group(ith)%jz(n) = hspsqn(ith,i)%m

        group(ith)%start(n) = i
        group(ith)%fin(n-1) = i - 1
     end if
     hspsqn(ith,i)%group = n

  end do
  group(ith)%fin(n) = nhsps(ith)

  return
end subroutine group_sp
!===========================================================================
!	
	
	
	
end module groups
	
	
	
	
	