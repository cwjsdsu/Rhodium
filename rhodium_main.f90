


use io
use wfn_mod
use nodeinfo
use program_info
use basisfilehandler
use menu_mod
use spmap_mod
use sporbit
use bmpi_mod

implicit none
integer :: ierr
integer :: nthreads
integer :: omp_get_num_threads

call BMPI_INIT(ierr)
icomm = MPI_COMM_WORLD
call BMPI_COMM_RANK(icomm,iproc,ierr)
call BMPI_COMM_SIZE(icomm,nproc,ierr)

if(iproc==0)then
write(6,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*)'~                                                     ~'
write(6,*)'~  Welcome to RHODIUM                                 ~'
write(6,*)'~  A post-processing companion to BIGSTICK            ~'
write(6,*)'~ ', version,lastmodified ,'                      '
write(6,*)'~                                                     ~'
write(6,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*)' '
!$OMP PARALLEL
   nthreads = omp_get_num_threads()
!$OMP END PARALLEL
print*,' Using ',nthreads,' OpenMP threads'
write(6,*)' '
endif

print*,' '


!............. READ FIRST .BAS FILE..............

call openbasisfile('INI')

if(iproc==0)print*,' about to read in basis'
call readbasisfile('INI')


!............. READ FINAL .BAS FILE..............
call openbasisfile('FIN')
call readbasisfile('FIN')
call maporbits
call map_initial2final_rhsps

call main_menu
!call BMPI_BARRIER(icomm,ierr)
call BMPI_FINALIZE(ierr)

end



