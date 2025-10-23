!
!  MPI master notes:
!    module nodeinfo :
!     nproc = # of MPI proceses
!     iproc = rank of each MPI process (from 0, which is master root node, to nproc-1)
!    (note: for simulations, nprocs = simulated # of MPI process, and iproc_dummy simulated rank)
!     icomm = communicator group for all processes
!
!    module localvectors
!      on each MPI process:
!      frag1, frag2 = assigned fragments for vec1, vec2; these are set in setup_localvectors
!      v1s, v1e are start, end for basis indices on vec1, so vec1(v1s:v1e)
!           note: same as basestart(frag1), basestop(frag1)
!      v2s, v2e are start, end for basis indices on vec2
! 
!      fcomm1, fcomm1 = communicators for reducing vec1 and vec2; set in BMPI_COMM_SPLIT (called by setup_localvectors)
!      isfragroot: logical flag to denote a "root" process assigned a fragment; set in setup_mpi_hist_comm
!
!=============================================================
!
!  for "new" parallelization, store vectors (or fragments of vectors) 
!  in this module; this reduces the # of information for subroutine calls
!
!  by making the vectors targets, we can use pointers to switch between them
!
!  when operating in parallel, can swap between vec1 and vec2
!  vecstore is used when storing the lanczos vectors all in RAM
!
!  With fragments, we only compute one direction
!
module localvectors
   use precisions
   implicit none
   ! useVec2Thread gives each thread an independent buffer to write to.
   logical :: wantUseVec2Thread = .false.   ! combine with numthreads > 1 to set useVec2Thread
   logical :: useHZSomp = .false.    ! flag for new OMP merged in 7.7.3 from HZS 
   double precision :: mpiStartTime
   
   logical :: useVec2Thread
   integer :: ompNumThreads   ! top level number of OpenMP threads

   real(kind=lanc_prec), allocatable, target :: vec1(:), vec2(:)
   ! KSM  Add thread local vec2 clones for performance testing
   !! real(kind=lanc_prec), allocatable, target :: vec2thread(:, :) ! (v2s:v2e, 0:num_threads-1)
   ! make new vec2thread flat to avoid bug in old fortran compilers
   integer(kind=basis_prec) :: vec2threadchunk, vec2threadchunkm1, vec2threadend
   real(kind=lanc_prec), allocatable, target :: vec2threadflat(:) ! (0 : num_threads * (v2e - v2s + 1))
   ! KSM  {v1s,v1e,v2s,v2e}={basestart(frag1), basestop(frag1), basestart(frag2), basestop(frag2)}
   integer(kind=basis_prec) :: v1s, v1e, v2s, v2e   ! slice start and end for vec1 and vec2
   integer, target :: frag1,frag2  ! which fragments each vector is
   ! note:   frag1 = nodal(iproc)%ifragment,  see setup_localvectors

   ! KSM - communicators for broadcast and reduction across
   ! one fragment.  We need one for each direction of vchar
   integer :: fcomm1, fcomm2  ! communicators to use for bcast/reduce on frag1/frag2
   logical :: isfragroot  ! node to scatterv/gatherv to slices, bcast back to fragment matvec nodes

!........... FOR DISTRIBUTED STORAGE AND REORTHOGONALIZATION....
!            the following arrays are used for controlling
!            distributed storage of Lanczos vectors
!            and their reorthgonalization
!
!  We do not store an entire Lanczos vector on a node
!  instead we store a "piece" (a fraction of a vector) 
!  and on a given node we store that same gene for all 
!  Lanczos vectors 
! 
   logical :: storelanczosincoreMPI     ! store lanczos vectors across MPI cores
   logical :: storelanczosincore1       ! store lanczos vectors on a sing
   integer   :: npiece   ! # of pieces of the vector
   integer(8):: Lpiece  ! size of a "piece" of the vector
   integer(8) :: Ldim, Lstart
   integer(8), allocatable :: Ldim_r(:), Lstart_r(:)
   integer(kind=8) :: buff_max
   
   integer :: niter  ! temporary
  
end module localvectors
!=============================================================

!
module nodeinfo
  implicit none
  integer(4) :: iproc,nproc,icomm

  integer(4) :: nprocs    ! used as a dummy
  integer(4) :: iproc_dummy   ! used for simulations

  logical :: distributeMPI =  .true.    
  logical :: simulateMPIdefault   =  .false.
  logical :: simulateMPI  

!...........................................................................
end module nodeinfo

