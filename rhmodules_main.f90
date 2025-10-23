!===========================================================================
! MODULES FOR RHODIUM
!===========================================================================



!-------------- OVERALL QUANTUM NUMBERS OF THE SYSTEMS----------------------
module system_parameters
  implicit none
  integer,target :: np_i(2),np_f(2) ! Number of particles  in initial, final systems
  integer,target :: jz_i, jz_f  ! M value for initial, final systems
  integer :: iparity_i,iparity_f  ! parity for initial, final systems
  character(len=1) :: cparity_i,cparity_f
  logical :: bothparities_i,bothparities_f
  integer :: Acore    ! a of the core
  
!...... INFORMATION FOR PARTICLE-HOLE CONJUGATION........
  logical :: phconj(2)
  logical, target:: phconj_i(2),phconj_f(2)  ! THESE NEED TO BE THE SAME
  integer :: npmax_i(2),npmax_f(2)
  integer :: npeff_i(2),npeff_f(2)  ! effective number of particles regarding particle-hole conjugation;
                       ! added in 7.5.6

end module system_parameters

!...........................................................................
module menu_choices
  implicit none
  character(len=3)  :: menu_char
  logical :: menu_dx_omathematica ! density in mathematica format
end module menu_choices

!------------------------------------------
! MODULE PROGRAM_INFO
!
! contains information about code
!

  module program_info
     implicit none

     character(6) :: version = '0.9.1 '
     character(9) :: lastmodified = 'June 2025'

  end module program_info


!-------------------------------------------------------------------
!  MODULE VERBOSITY
!-------------------------------------------------------------------
!  Various flags for printing out information.
!  This allows one to easily print out or suppress information.
!
!  NOTE: verbose_X  means write out information to terminal
!        chatty_X  means write out information as being computed (rare)
!        print_X    means write to file
!        X_file     is unit number for file
!-----------------------------------------------------------------
module verbosity
  implicit none

  logical :: print4modelinfo       ! prints out information needed for modeling
  integer :: modelinfo_file    =  44
  logical :: verbose_winnow    = .false.
  logical :: print_groups      = .false.
  integer :: group_file        = 37
  logical :: print_hsps        = .false.
  integer :: hsps_file         = 37
  logical :: verbose_templates =  .false. 
  logical :: verbose_haikus    = .false.
  logical :: print_haikus      = .false.  
  integer :: haiku_file        =  37
  logical :: verbose_blocks    = .false.
  logical :: verbose_hops      = .false.
  logical :: print_hops        = .false.
  integer :: hop_file          =  38 
  logical :: verbose_sectors   = .false.
  logical :: print_sectors     = .false.
  integer :: sector_file       = 37  
  logical :: chatty_descent    = .false.
  logical :: chatty_links      = .false.

  logical :: print_genes       = .false.
  integer :: gene_file         = 38
  logical :: print_basis       = .false.
  integer :: basis_file        = 31

  integer :: xsd_file          = 32  ! added 4/2011 CWJ
  logical :: verbose_uncouple  = .false.
  
  logical :: verbose_jumps  = .false.

  logical :: print_matrix      = .false.
  integer :: matrix_file       = 55
  logical :: verbose_orthog    = .false.

  logical :: print_block_ops   = .false.
  logical :: print_node_load   = .false.

  logical :: print_sectorjumps = .false.    ! added 5/2010 CWJ
  integer :: sectorjump_file   = 39

  logical :: print_jumps = .false.   ! added 4/2011 CWJ
  integer :: jump_file = 46

  logical :: verbose_3body_readin = .false.
  logical :: verbose_3bodyjumps = .false.

  logical :: verbose_distrib = .false.  ! added 8/2011 CWJ
  logical :: verbose_bundles  = .true.  ! added 10/2011 CWJ

  logical :: verbose_fragments = .false.
  logical :: verbose_distro    =.false.

end module verbosity

!...........................................................................
module bitstuff
  implicit none
  integer max_bit_word         ! the max. # of usable bits in a
  parameter (max_bit_word=30)  ! given word to describe an SD
end module bitstuff

!---------------------------------------------------------------------------
!  declaring size of variables
!---------------------------------------------------------------------------
module precisions
  integer,parameter :: lanc_prec =    4
  integer,parameter :: basis_prec   = 8
  integer,parameter :: obs_prec     = 8
  integer,parameter :: egv_prec     = 8
  integer(4) :: MPI_lanc
end module precisions


!...........................................................................
!
!  GATHERS MISCELLANEOUS INFORMATION FOR END-OF-RUN REPORT
!  
module reporter
   implicit none
   character(len=40) :: spfilename
   character(len=50) :: tbmefilename(50)    ! assumes no more than 50 input files
   real              :: tbmescale(50),spescale(50)      ! scaling of 2-body matrix elements,spes
   character(1)      :: tbmetype(50)        ! type of TBME
   integer           :: ntbmefiles          ! # of tbme files; needs to be initialized
   character(len=50) :: wfn_in_filename

end module reporter
!...........................................................................
module io
  use precisions
  implicit none
  
  character(55),target :: wfnfilename_i,wfnfilename_f   ! initial, final wave function file names
  character(55),target :: basfilename_i,basfilename_f   ! initial, final basis function file names
  integer :: wfnfile_i = 26   
  integer :: wfnfile_f = 27
  integer,target :: basfile_i = 28
  integer,target :: basfile_f = 29
  
  character(55) :: outfile              ! filename for outputs
  logical       :: writeout
  integer       :: autoinputfile   = 8
  integer       :: resultfile      = 12
  integer       :: wfnfile         = 13
  integer       :: oldwfnfile      = 14
  integer       :: mdensfile       = 15 ! mathematica density file
  integer       :: logfile         = 36 ! human-readable log file
  integer       :: binlogfile      = 22 ! binary log file, used for speeding up creations
  integer       :: fragfile        = 46 ! file for fragment info
  logical       :: auto_readin         ! flag to read in information from .wfn file
  logical       :: auto_input          ! flag to read from autoinput.bigstick file
  logical       :: ham_readin          ! flag to read in hamiltonian files
  logical       :: op_readin           ! flag to read in (nonscalar) one-body operator
  logical       :: strengthflag        ! flag to compute strength functions starting from a pivot
  logical       :: densityflag         ! to compute density matrices inline
  logical       :: trdensout           ! a flag to write to file suitable for Navratil's TRDENS code
  logical       :: restart_enabled = .true.  ! if set = .false., makes for easier restart
                                       ! else lanczos vectors automatically erased
  logical       :: write_wfn           ! flag to write wavefunctions to file (.wfn) (added by P.G.K)
  logical       :: get_JT               ! flag to enable calculation of J, T; added 7.6.8
  

  integer(kind=basis_prec) :: dimbasischeck

  ! KSM:  Added to support writing large files to scratch directories
  !       on HPC machines
  character (len=64) :: nersc_host
  character (len=1024) :: base_scratch_dir  ! like /tmp or $SCRATCH
  character (len=1024) :: scratch_dir       ! add on resultfile
end module io

!--------------------------------------------------------------------------
!  single-particle ORBIT information
!  key derived variable: orbqn
!---------------------------------------------------------------------------
module sporbit
   implicit none
   integer,target :: numorb_i(2),numorb_f(2)  ! # of s.p. orbits
   integer :: numorbmax	        ! max of numorb
   logical :: isoflag           ! in isospin formalism
   logical :: pnwtflag          ! flag to allow p,n to have different weights even if in isospin formalism
                                ! added in version 7.2.8 by CWJ
   logical :: spinless         ! for "spinless" fermions
   logical :: allsameparity_i,allsameparity_f     ! flag for all sp orbits same parity
   logical :: allsameW_i,allsameW_f         ! flag for all sp orbits same W
!------------ CREATE A DEFINED TYPE-----------------------------------------
   type orb
      integer :: nr            ! radial quantum number
      integer :: j             ! 2 x j
      integer :: l             ! L
      integer :: par           ! parity = +/- 1
      integer :: w             ! excitation 
	  integer :: ifmap         ! map initial <--> final orbits
   end type orb
   type (orb),allocatable,target :: orbqn_i(:,:),orbqn_f(:,:) ! orbqn(species,iorb)

   character*200 :: sps_path    ! environment variable for searching for directory containing single particle orbits
   integer :: length_sps_path
   integer :: maxorblabel       ! added 7.4.1; used for error traps in reading in matrix element files     
logical :: sameweightscheme  ! IF .TRUE., then the orbits all have the same W, and we can map in a straightforward fashion
                             ! OTHERWISE, HAVE TO DO MORE SEARCHING 
contains
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! "multiplies" parities where 1 = + and 2 = -
!
!  INPUT: par1,par2
!
      integer function parmult(par1,par2)

      implicit none
      integer par1,par2  ! input
!..........................................
      integer i1,i2,i3

      if(par1 < 1 .or. par2 > 2)then
        print*,' wrong parity ',i1,i2
        stop
      endif

      if(par2 < 1 .or. par2 > 2)then
        print*,' wrong parity (2)',i1,i2
        stop
      endif

      i1 = 3-2*par1
      i2 = 3-2*par2
      i3 = i1*i2
      parmult = (3-i3)/2

      return
      end function parmult
	  
!........... ROUTINE TO MAP ORBITS BACK AND FORM

      subroutine maporbits
		  implicit none
		  integer :: it
		  integer :: i,f
		  
		  do it = 1,2
			  orbqn_f(it,:)%ifmap = 0
			  do i = 1,numorb_i(it)
				  do f = 1,numorb_f(it)
					  if(orbqn_f(it,f)%ifmap/=0)cycle
					  if( orbqn_i(it,i)%j/=orbqn_f(it,f)%j)cycle
					  if( orbqn_i(it,i)%nr/=orbqn_f(it,f)%nr)cycle
					  if( orbqn_i(it,i)%l/=orbqn_f(it,f)%l)cycle
					  orbqn_i(it,i)%ifmap = f
					  orbqn_f(it,f)%ifmap = i
					  exit
					  
				  end do			  
			  end do
	  
		  end do
		  return
		  
	  end subroutine maporbits	  
	  
end module sporbit

!---------------------------------------------------------------------------
!  single-particle STATE information
!  key dervied variable: spsqn
!---------------------------------------------------------------------------
module spstate
  implicit none
  integer, target :: nsps_i(2),nsps_f(2) ! # of s.p. states
  integer :: nspsmax
!------------ CREATE A DEFINED TYPE-----------------------------------------
  type spst
     integer :: nr       ! radial quantum number
     integer :: j        ! 2 x j
     integer :: m        ! 2 x jz
     integer :: l        ! L
     integer :: w        ! excitation 
     integer :: par      ! parity
     integer :: orb      ! orbit label
     integer :: group    ! which group -- labels by jz, w, par
	 integer :: ifmap      ! ADDED FOR RHODIUM: maps initial s.p. states to final s.p. states AND VICE VERSA
	 integer :: tr      ! ADDED FOR HOLE FORMALISM -- map of conjugate single particle state
  end type spst
  type (spst),allocatable,target :: spsqn_i(:,:), spsqn_f(:,:)  ! spsqn_x(it,istate)

!--------- INFO ON HAIKU S.P. STATES:---------------------------------------
!  haikus are divided up into "left" built from s.p. states with jz < 0
!                         and "right" built from s.p. states with jz >= 0
!------------------ CONVENTION for hspsqn (and elsewhere)-------------------
!                   if ispecies < 0, then "left" haiku  jz < 0
!                               >0        "right" haiku  jz >= 0
!...........................................................................
  integer :: nhsps0(-2:2)  
  integer :: nhspsmax
  integer,target :: nhsps_i(-2:2),nhsps_f(-2:2)
  type (spst),allocatable,target :: hspsqn_i(:,:),hspsqn_f(:,:) ! hspsqn(ispecies,istate)

!---------------------- GROUPS -- group together hsps by par, jz, par-------
  integer, target :: ngroups_i(2),ngroups_f(2)

  type grouplist
     integer, pointer :: jz(:)
     integer, pointer :: w(:)
     integer, pointer :: par(:)
     integer, pointer :: start(:)
     integer, pointer :: fin(:)
     logical, pointer :: active(:)   
  end type grouplist
  
  type (grouplist), target :: group_i(2),group_f(2)
	  
!----------- FOR USE IN RHODIUM ONLY------
!            SINGLE PARTICLE STATES ACTUALLY USED AND 
!            IN THE ORDER USED FOR CREATING SLATER DETS
!   see file rhmaplib.f90 for useful routines in creating these

   integer,target :: nrhsps_i(2),nrhsps_f(2)
   type (spst),allocatable,target :: rhspsqn_i(:,:),rhspsqn_f(:,:) ! hspsqn(ispecies,istate)
	   
   logical :: i2fspoutoforder ! if TRUE then initial -> final s.p. states out of order

end module spstate

!---------------------------------------------------------------------------
!  Information on W (excitation weighting) 
! (equivalent to Nhw in REDSTICK, t in ANTOINE)
!---------------------------------------------------------------------------
module W_info
  implicit none
  integer :: maxW(2)      ! max excitation allowed for p, n
  integer :: minW(2)      ! min excitation allowed for p, n 
  integer :: maxWtot_i,maxWtot_f,minWtot     ! max,min total excitation
  integer :: wceiling(2)  ! ceiling to limit creating too many haikus
end module W_info

! NOTE: in 7.6.9 module spsections MOVED to file bsections.f90



!---------------------------------------------------------------------------
!  sectors of the basis are divided by (first) jz, then parity, finally W
!  sectors are constructed from blocks of haikus
!---------------------------------------------------------------------------
module sectors

  implicit none
  integer,target     :: nsectors_i(2),nsectors_f(2)                 ! # of sections
  type sectorinfo
     integer(8)  :: xsdstart                          ! where, in the list of X-species sds, does this start
     integer(8)  :: xsdend
     integer(8)  :: nxsd                              ! # of SDs in this sector
	 integer(8)  :: ncxsd
     integer  :: jzX
     integer  :: parX
     integer  :: wX
     integer  :: ncsectors                         ! # of opposite-species sections can couple to
     integer, pointer :: csector(:)                ! list of opposite-species sections
     integer  :: nhblocks                          ! # of pairs of haiku blocks used to construct this section
     integer, pointer :: rhblock(:),lhblock(:)     ! list of right- and left-haiku blocks 
     integer(8), pointer :: blockstart(:),blockend(:) ! where, for a given block, this starts and ends
     integer,pointer :: cut(:)                     ! ADDED 7.4.1 to be used to cut proton sectors for fragments
     integer         :: ncuts
     integer(kind=8) :: basisstart,basisend        ! ADDED 7.4.9: for protons, where basis starts and stops
     logical :: ibelong                            ! ADDED 7.5.9: denotes if a sector belongs on a given MPI process
                                                   ! especially when vectors are broken into fragments
  end type sectorinfo
  type mastersect
     type (sectorinfo), pointer :: sector(:)
  end type mastersect
  
  type (mastersect),target :: xsd_i(2),xsd_f(2)
end module sectors



!===========================================================================
! Observables, densities, applied operators, etc.
!===========================================================================
module obs
  use precisions
  implicit none
  real(kind=obs_prec) :: xj2,xt2
  real(kind=4),allocatable :: energy(:),xjlist(:), xtlist(:)
  real(kind=obs_prec) :: xj2tmp,xt2tmp ! for global sum /allreduce/ (added by P.G.K.)
  logical :: twoobsflag
  
end module obs
   
!=======================================================================
!module densities moved to bdensities.f90
!
!===================================================================
!
! coupled matrix elements for transition operators   ADDED 7.5.2
! FOR application to a wavefunction
! (corresponding uncoupled matrix elements found in module onebodypot)

module opmatrixelements
	implicit none
	logical :: pnoperators    ! flag to signal if operators are in pn format or isospin format
	  real, allocatable :: op1bod(:,:)  ! reduced matrix elements for one-body operartor; to be deleted
	
	real,allocatable,target :: pop1bod(:,:),nop1bod(:,:)  ! singly-reduced matrix elements for proton, neutron 1-body operators
	integer :: jop                         ! J of transition operator
	integer :: top                         ! (optional) T of transition operator
	
    logical :: subtract_enabled=.false.       ! for MKGK application of r^2 Y00
    logical :: subtract
	
end module opmatrixelements

!===================================================================
!
! application of one-body (uncoupled) potentials, including self-consistent mean-field for starting
! the original coupled matrix elements found in module opmatrixelements
!
      module onebodypot

      implicit none

      real, pointer :: pPOT(:,:),nPOT(:,:)

      real, allocatable, target :: ppot_obs(:,:),npot_obs(:,:)
      real, allocatable, target :: ppot_h(:,:),npot_h(:,:)
	  
	  logical :: meanie   ! to do self-consistent mean-field

      end module onebodypot

