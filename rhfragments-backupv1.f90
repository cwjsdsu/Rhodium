
!=================================================================================
!
!  bfragments.f90
!
!  routines to break up Lanczos vectors into "fragments"  
!
!============================================================
!
!  The next module is used for 
!  breaking up the Lanczos vectors 
!
!=========================================================
!
!  FRAGMENTS
!
!   contiguous (sub) sectors of Lanczos vectors grouped together
!   
! started 8/2011 by CWJ
!

module fragments
   use precisions
!   use operation_stats
   implicit none

   logical :: fragment_vectors   ! flag to break up vectors into fragments
   integer :: nfragments         ! number of fragments
   logical :: useNewReorthog1 = .true.  ! Use new reorthog even for single fragment
   logical :: useNewReorthog  ! set to (nfragments > 1) .or. useNewReorthog1

   type frag
        integer :: nsectors  ! # of subsectors
        integer(8) :: ssectorstart,ssectorend  ! start, end of sectors; these are defined by "proton" sectors
		  integer(8) :: csectorstart,csectorend  ! start, end of conjugate sectors, defined by "neutron" sectors, must check contiguity
        integer(8):: basisstart,basisend  ! where the basis starts and stops for this sector
        integer(8):: localdim      ! local dimension = basisend-basisstart+1
        real(8)    :: totops    ! tot# of operations to/from this fragment

   end type frag

   type (frag), allocatable,target :: fragmentlist_i(:),fragmentlist_f(:)

!   integer, allocatable :: sectormap2frag(:)  ! maps a sector to a fragment  SCHEDULED FOR OBSOLESCENCE

!--------- THIS ALLOWS US TO STORE JUST FRAGMENTS OF A VECTOR --

   integer (kind=basis_prec), allocatable :: basestart(:),basestop(:)  !NEED TO GET RID OF LATER

   integer (kind=basis_prec), allocatable,target :: basestart_i(:), basestop_i(:), basestart_f(:), basestop_f(:)
   
   logical :: test_fragments = .false.  ! added in 7.4.6, used to test new fragment scheme on a single processor
   
!------------- JUMP STORAGE ACROSS FRAGMENTS (f2f matvec blocks)-------------- added in 7.6.2 -----------------

   type f2fjump
	   integer(8) :: nX1bjump(2),nXX2bjump(2),nXXX3bjump(2)
	   real(8)    :: totjumpstorage
	   integer(4) :: minprocs     ! minimal assigned MPI procs/nodes
	   integer(4) :: nprocs       ! # of actual assigned procs in the end
   end type f2fjump   
   
   type (f2fjump), allocatable,target :: f2fjumpstat(:,:)

!------- INFORMATION ON WHAT FRAGMENTS ON WHAT MPI PROCESSES
!        previously in module tribution, moved in 7.6.7
	   
	type whatsonanode
	      integer :: ifragment,ffragment   ! which fragments are here
	      logical :: ifirst ! first node in group servicing ifragment
		  
!	      integer :: sfragment  ! which fragment is stored here; later might want to store more than one
!	      real :: frac_ops   ! what fraction of possible ops are stored here      
	end type whatsonanode

	type (whatsonanode), allocatable :: nodal(:)	   

contains
!======================================================================
!
!  SUBROUTINES CONTAINED IN THIS FILE
!    fragmenter:  master routine to break up Lanczos vectors into fragments
!        which calls:
!    vectorstorageceiling
!    defaultfragments:  routine to set fragment information when only 1 fragment
!    count_create_fragments_blocks: updated fragmentation routine which breaks fragments along haiku blocks
!    reconstruct_sector: when sectors are subdivided, reconstruct proton sectors
!    fragmentvectorsize:  sets base start/stop for a given fragment
!       also in this file (called by distro_opbundles_over_fragments )
!    divide_procs_over_fragments: routine to divide up available nodes over fragments
!
!=======================================================================
!
!  Notes 7.4.6 (March 2015) regarding "new" fragmentation based upon proton haiku blocks
!  To test/debug fragmentation, run on a single processor;
!  in order to do this, we must modify routines which divide up work over MPI processes
!   
!  key routines include: 
!          setnodaltribution  (in bparallel_lib1.f90)  assigns initial, final fragments to MPI procs/compute nodes
!          setup_localvectors (in blanczoslib1.f90)  allocates memory for (fragmented) lanczos vectors
!  If one simulates a fragmented case on a single processor, the above have to be modified
!

!======================================================================
!
!  proc_range
!   divides up list(sectors) based on nproc and current mpi rank.
!   returns start and stop of sectors for each rank
!
! CALLED BY:
!  proj_boss rhprojlib.f90
!
subroutine proc_range(lv, hv, nprocs, irank, istart, iend)
   implicit none
   integer(8) :: lv, hv ! Low and high iteration variables.
   Integer(4) :: nprocs, irank, istart, iend
   integer(4) ::  q, r ! quotient and remainder
   
   q = (hv-lv + 1)/nprocs
   r = MOD(hv-lv+1,nprocs)
   
   istart = lv + irank*q + MIN(irank,r) ! starting iter
   iend = istart + q - 1 ! final iter
   
   if (r .gt. irank) then
      iend = iend + 1
   endif
   return
end subroutine proc_range

!======================================================================
!
!  fragmentvectorsize
!     sets base start/stop for a given fragment
!
! CALLED BY:
!  set_nfragments
!
subroutine fragmentvectorsize
   !   use fragments
      use sectors
      implicit none
      integer :: fg
      integer :: aerr

   !assign intial and final basis start/stop to fragments based on fragment sectors. 
      do fg = 1, nfragments
         basestart_i(fg)=xsd_i(1)%sector(fragmentlist_i(fg)%ssectorstart)%basisstart
         basestop_i(fg)= xsd_i(1)%sector(fragmentlist_i(fg)%ssectorend)%basisend
 
         basestart_f(fg)=xsd_f(1)%sector(fragmentlist_f(fg)%ssectorstart)%basisstart
         basestop_f(fg)= xsd_f(1)%sector(fragmentlist_f(fg)%ssectorend)%basisend
       enddo
      return
   end subroutine fragmentvectorsize
   
!======================================================================
   


!======================================================================
!  subroutine set_nfragments
! Sets nfragments and side effects of assignment
subroutine set_nfragments(nf)
   !  use fragments
     use basis
     use nodeinfo
     use sectors
     use bmpi_mod
     use localvectors
     implicit none
     integer :: jproc, nf
     integer :: i ! mpiS
     integer :: is, ie
     integer :: ierr, aerr
     integer :: basis_istart, basis_iend
     integer :: sectstart, sectend
     integer(8) :: st, fg
     st = 1
     is = 0
     ie = 0
  
     nfragments = nf
     useNewReorthog = useNewReorthog1 .or. nfragments > 1
     if(iproc == 0) print *, "useNewReorthog=", useNewReorthog
     allocate(nodal(0:nproc-1))
     
     
     if (iproc == 0 )then
      print*,"printing intitial space sector information..."
        do i = 1, nsectors_i(1)
           print*, "i", i, "xsd_i(1)%sector(i)%jzX", xsd_i(1)%sector(i)%jzX,"xsd_i(1)%sector(i)%parx",xsd_i(1)%sector(i)%parx
        end do
     end if 

     print*,"finished"
    
     do jproc = 0,nproc-1
        nodal(jproc)%ifragment= 1 + jproc
        nodal(jproc)%ffragment= 1 + jproc
     end do
     
     allocate(basestart_i(nf), basestop_i(nf))
     allocate(basestart_f(nf), basestop_f(nf))
     
     if(nfragments == 1) then
        basestart_i(1)=1
        basestop_i(1)=dimbasis_i
        basestart_f(1)=1
        basestop_f(1)=dimbasis_f
  
     else if(nfragments .gt. 1) then

   ! used temporarily to store fragment start stop sector info
      allocate(fragmentlist_i(nf), stat=aerr)
      allocate(fragmentlist_f(nf), stat=aerr)

      !assign sectors from the intital and final basis to each fragment
      if (iproc == 0 )then
         print*,"printing final sector information..."
           do i = 1, nsectors_f(1)
              print*, "i", i, "xsd_f(1)%sector(i)%jzX", xsd_f(1)%sector(i)%jzX,"xsd_f(1)%sector(i)%parx",xsd_f(1)%sector(i)%parx
           end do
         endif

      !we divide the final basis into nproc fragments by sector
      ! each fragment has some range [istart,istop] of contigous sectors assigned to it
      call assign_fragments

      !set basis start/stop using fragment sector info
      call fragmentvectorsize


      !assign intial and final basis start/stop to fragments based on fragment sectors. 
      !do fg = 1, nfragments
      !  basestart_i(fg)=xsd_i(1)%sector(fragmentlist_i(fg)%ssectorstart)%basisstart
      !  basestop_i(fg)= xsd_i(1)%sector(fragmentlist_i(fg)%ssectorend)%basisend
      !
      !  basestart_f(fg)=xsd_f(1)%sector(fragmentlist_f(fg)%ssectorstart)%basisstart
      !  basestop_f(fg)= xsd_f(1)%sector(fragmentlist_f(fg)%ssectorend)%basisend
      !enddo

      
     if (iproc == 0 )then
      print*,"printing final sector information..."
        do i = 1, nsectors_f(1)
           print*, "i", i, "xsd_f(1)%sector(i)%jzX", xsd_f(1)%sector(i)%jzX,"xsd_f(1)%sector(i)%parx",xsd_f(1)%sector(i)%parx
        end do
      endif
     end if 
    
     call BMPI_BARRIER(icomm,ierr)
     return
  
  end subroutine set_nfragments
  
  ! subroutine distbasis. This subroutine calculates the number of 
  !basis elements to be given to the set of remaining mpi ranks.
  ! Makes sure more even distribution so the nproc-1 fragment doesn't have too few basis elements.
  subroutine assign_fragments
   !  use fragments
     use basis
     use nodeinfo
     use sectors
     use bmpi_mod
     use localvectors
     implicit none
     integer :: is, ic, fg, localdim, dim_sum
     integer :: basesplit_ind, jztemp_i, partemp_i 
     integer :: jztemp, partemp, jztemp_f2, partemp_f2
     integer :: s_shift, i
     integer :: s, r, n !used to more evenly distribute elements
     logical :: startfound
     

     startfound = .false.
     dim_sum = 0
     s_shift = 0
     r = 0 ! running sum of sector dimensions
     s = 0 ! ranks with fragments already
     fg = 1
     n = 1

     
 
     do is = 1, nsectors_f(1)
      if(s .ne. nproc ) basesplit_ind = abs((dimbasis_f - r)/(nproc - s))! approximate fragment break
      localdim = xsd_f(1)%sector(is)%basisend - xsd_f(1)%sector(is)%basisstart + 1
      !print*,"localdim", localdim
         dim_sum = dim_sum + localdim
         print*,"-----------------------------------------------------------"
         print*,"is",is,"s",s,"r",r,"basesplit_ind",basesplit_ind,"localdim",localdim
         print*,"dim_sum",dim_sum,"r+basesplit_ind",r+basesplit_ind
         print*,"-----------------------------------------------------------"
         !record fragment sector start stop
         if (dim_sum .ge. (r+basesplit_ind)) then ! r holds dim_sum from previous iteration
            s_shift = 0
            !check next sectors to make sure we get all sectors of the same Jz and parity on the same fragment
            do ic = is+1, nsectors_f(1)-1
               jztemp_f2 = xsd_f(1)%sector(ic)%jzX 
               partemp_f2 = xsd_f(1)%sector(ic)%parX
               if((jztemp_f2 .eq. xsd_f(1)%sector(is)%jzX).and.(partemp_f2 .eq. xsd_f(1)%sector(is)%parX) ) then
                  s_shift = s_shift + 1
               end if
            enddo

            if(s .eq. 0) then ! iproc
               fragmentlist_f(fg)%ssectorstart = 1 
               fragmentlist_f(fg)%ssectorend  = is + s_shift
              

               print*,"here a ---"
               print*,"is:",is,"running sum:", dim_sum,"basesplit_ind",basesplit_ind,"shift",s_shift 
               print*,"sec start:",fragmentlist_f(fg)%ssectorstart,"sec end:",fragmentlist_f(fg)%ssectorend
               
               !increment
               s = s + 1
               n = n + 1
               fg = fg + 1
               r = dim_sum
               
               cycle
            else
               if(fg .eq. nproc)then
                  if(is .le. nsectors_f(1)) then
                     fragmentlist_f(fg)%ssectorstart = fragmentlist_f(fg-1)%ssectorend + 1
                     fragmentlist_f(fg)%ssectorend  = nsectors_f(1)
                     print*,"here b ---"
                     print*,"is:",is,"running sum:", dim_sum,"basesplit_ind",basesplit_ind,"shift",s_shift
                     print*,"sec start:",fragmentlist_f(fg)%ssectorstart,"sec end:",fragmentlist_f(fg)%ssectorend
                     !increment
                     s = s + 1
                     n = n + 1
                     fg = fg + 1
                     r = dim_sum
                     
                     cycle
                  endif
               else
                  s_shift = 0
                  do ic = is+1, nsectors_f(1)-1
                     jztemp_f2 = xsd_f(1)%sector(ic)%jzX 
                     partemp_f2 = xsd_f(1)%sector(ic)%parX
                     if((jztemp_f2 .eq. xsd_f(1)%sector(is)%jzX).and.(partemp_f2 .eq. xsd_f(1)%sector(is)%parX) ) then
                        s_shift = s_shift + 1
                     end if
                  enddo
                  if((fragmentlist_f(fg-1)%ssectorend + 1) .ge. is ) cycle !incase at boundary of Jz/M par and shift = 0
                     fragmentlist_f(fg)%ssectorstart = fragmentlist_f(fg-1)%ssectorend + 1
                     fragmentlist_f(fg)%ssectorend  = is + s_shift
                     print*,"here c ---"
                     print*,"is:",is,"running sum:", dim_sum,"basesplit_ind",basesplit_ind,"shift",s_shift
                     print*,"sec start:",fragmentlist_f(fg)%ssectorstart,"sec end:",fragmentlist_f(fg)%ssectorend


                  !increment
                  s = s + 1
                  n = n + 1
                  fg = fg + 1
                  r = dim_sum
                  
                  cycle
               endif
               
   
            end if 
            
            !n = n + 1
            !fg = fg + 1
            
            
         end if
      enddo



      ! now match Jz and par with the initial basis

      do fg = 1,nfragments !iterate over each fragment
         do is = 1, nsectors_i(1)
            !organization...
            jztemp_i = xsd_i(1)%sector(is)%jzX   
            partemp_i = xsd_i(1)%sector(is)%parX
            !jztemp_f1 = xsd_f(1)%sector(fragmentlist_f(fg)%ssectorstart)%jzX 
            !partemp_f1 = xsd_f(1)%sector(fragmentlist_f(fg)%ssectorstart)%parX
            jztemp_f2 = xsd_f(1)%sector(fragmentlist_f(fg)%ssectorend)%jzX 
            partemp_f2 = xsd_f(1)%sector(fragmentlist_f(fg)%ssectorend)%parX
            !print*,"--------------------------------------------------------"
            !print*,"fg:",fg,"is:",is
            !print*,"jztemp_i",jztemp_i,"partemp_i",partemp_i,"jztemp_f2", jztemp_f2,"partemp_f2",partemp_f2
            !print*,"--------------------------------------------------------"

            if(fg .eq. 1) then ! frag 1 setup
               
               if((partemp_i .eq. partemp_f2) .and. (jztemp_i .eq. jztemp_f2 ) .and. (is .ne. nsectors_i(1)) ) then
                  !check next sectors to make sure we get all sectors of the same Jz and parity on the same fragment
                  if(is .ne. nsectors_i(1))then
                     s_shift = 0
                     do ic = is+1, nsectors_i(1)
                        jztemp = xsd_i(1)%sector(ic)%jzX 
                        partemp = xsd_i(1)%sector(ic)%parX
                        if((jztemp .eq. jztemp_f2).and.(partemp .eq. partemp_f2) ) then
                           s_shift = s_shift + 1
                        end if
                     enddo
                  endif
                  fragmentlist_i(fg)%ssectorstart = 1
                  fragmentlist_i(fg)%ssectorend = is + s_shift
                 ! print*,"a"
                 ! print*,"sec start:",fragmentlist_i(fg)%ssectorstart,"sec end:",fragmentlist_i(fg)%ssectorend
                  exit
               else if((partemp_i .eq. partemp_f2) .and. (jztemp_i .gt. jztemp_f2) ) then
                 
                  fragmentlist_i(fg)%ssectorstart = 1
                  fragmentlist_i(fg)%ssectorend = is-1
                 ! print*,"b"
                 ! print*,"sec start:",fragmentlist_i(fg)%ssectorstart,"sec end:",fragmentlist_i(fg)%ssectorend
                  exit
               else if( is .eq. nsectors_i(1)) then
                  fragmentlist_i(fg)%ssectorstart = 1
                  fragmentlist_i(fg)%ssectorend = nsectors_i(1)
                  !print*,"c"
                  !print*,"sec start:",fragmentlist_i(fg)%ssectorstart,"sec end:",fragmentlist_i(fg)%ssectorend
                  exit
               endif
            
            else
               if((partemp_i .eq. partemp_f2) .and. (jztemp_i .eq. jztemp_f2 ) ) then
                  if(is .ne. nsectors_i(1))then
                     s_shift = 0
                     do ic = is+1, nsectors_i(1)
                        jztemp = xsd_i(1)%sector(ic)%jzX 
                        partemp = xsd_i(1)%sector(ic)%parX
                        if((jztemp .eq. jztemp_f2).and.(partemp .eq. partemp_f2) ) then
                           s_shift = s_shift + 1
                        end if
                     enddo
                  endif

                  if((is + s_shift) .ge. nsectors_i(1)) then
                     fragmentlist_i(fg)%ssectorstart = fragmentlist_i(fg-1)%ssectorend + 1
                  fragmentlist_i(fg)%ssectorend = is
                  else
                     fragmentlist_i(fg)%ssectorstart = fragmentlist_i(fg-1)%ssectorend + 1
                     fragmentlist_i(fg)%ssectorend = is + s_shift
                  end if
                  !print*,"d"
                  !print*,"shift",s_shift, "is", is
                  !print*,"sec start:",fragmentlist_i(fg)%ssectorstart,"sec end:",fragmentlist_i(fg)%ssectorend
                  exit
               else if((partemp_i .eq. partemp_f2) .and. (jztemp_i .gt. jztemp_f2) .and. (is .ne. nsectors_i(1))) then 
                  fragmentlist_i(fg)%ssectorstart = fragmentlist_i(fg-1)%ssectorend + 1
                  fragmentlist_i(fg)%ssectorend = is-1
                  !print*,"e"
                  !print*,"sec start:",fragmentlist_i(fg)%ssectorstart,"sec end:",fragmentlist_i(fg)%ssectorend
                  exit
               else if( is .eq. nsectors_i(1)) then

                  if((fragmentlist_i(fg-1)%ssectorend + 1) .gt. nsectors_i(1)) then
                     !just put the entire vec1 on the node
                     fragmentlist_i(fg)%ssectorstart = 1
                     fragmentlist_i(fg)%ssectorend = nsectors_i(1)
                  else
                     fragmentlist_i(fg)%ssectorstart = fragmentlist_i(fg-1)%ssectorend + 1
                     fragmentlist_i(fg)%ssectorend = is
                  end if
                  !print*,"f"
                  !print*,"sec start:",fragmentlist_i(fg)%ssectorstart,"sec end:",fragmentlist_i(fg)%ssectorend
                  exit
               endif


            endif

         enddo
      enddo

      
     return
  
  end subroutine assign_fragments

  subroutine printfragsectorinfo
   !  use fragments
     use basis
     use nodeinfo
     use sectors
     use bmpi_mod
     use localvectors
     implicit none
     integer :: is, fg, localdim, dim_sum
     integer :: basesplit_ind, jztemp_i, partemp_i 
     integer :: jztemp_f1, partemp_f1, jztemp_f2, partemp_f2
     !real :: basesplit_ind
     integer :: s, r, n !used to more evenly distribute elements
     logical :: startfound
     print*,"printing fragment sector information..."
     print*,"======================================================================"
     do fg = 1, nfragments
         print*," initial sectors on fragment:",fg, "start sector:",fragmentlist_i(fg)%ssectorstart
         print*,"initial sectors on fragment:",fg, "stop sector:",fragmentlist_i(fg)%ssectorend
         print*," final sectors on fragment:",fg, "start sector:",fragmentlist_f(fg)%ssectorstart
         print*,"final sectors on fragment:",fg, "stop sector:",fragmentlist_f(fg)%ssectorend
      enddo
     print*,"======================================================================"
   return
  end subroutine printfragsectorinfo


!======================================================================
!============================Archived==================================
!======================================================================
!  subroutine fragmenter
!
!  master routine to break up Lanczos vectors into "fragments"
!  when necessary it breaks up sectors into "subsectors"
!  by breaking along the proton basis Slater Determinants
!
!  The "non-optimized" version simply breaks lanczos vectors by sectors,
!  because that is easiest to implement.
!
!  The "block" version breaks up sectors by proton haiku blocks, which is a natural division
!  and preserves most of the downstream routines intact; the divided sectors must be then redefined
!
!  first count up the number of break points
!
!  CALLED BY:
!    main routine
!
!  SUBROUTINES CALLED:
!    vectorstorageceiling   checks to see if vectors get broken up, and how
!    count_create_fragments_noopt
!        [this creates non-optimized fragments by simply bundling together sectors;
!         later will use optimized fragments with subroutine count_create_fragments ]
!       OR
!    count_create_fragments_blocks
!    reconstruct_sector
!    [map_sectors2noopt_fragments] 
!    defaultfragments      only one fragment
!    fragmentvectorsize    sets start and stop for basis in fragments
!    fragmentconjugatesectors  checks which neutron sectors belong (added 7.6.2)
!
subroutine fragmenter
   
   use flagger
   use basis
   use nodeinfo
!   use sectors
!   use menu_choices

   implicit none

   
   !call vectorstorageceiling  !  check to see if vectors get broken up into fragments
!   if(break_vectors .or. test_fragments)then
!	   call test_sector_breaks(83)  ! only for testing
!       call count_create_fragments_blocks(.false.)
!       call count_create_fragments_blocks(.true.)
!       call reconstruct_sectors   
!   call test_sector_breaks(84)  ! only for testing
!   else
!      call defaultfragments
!   endif
!   call fragmentvectorsize
!   call fragmentconjugatesectors
   return
end subroutine fragmenter


end module fragments

