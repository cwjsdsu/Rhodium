###########################################################################
#  NOTE: FOR LIST OF OPTIONS 
#  "make help"
#                                                                          #
#                     Make file for BIGSTICK                               #
#                                                                          # 
#   1) Correctly insert the names of your compilers into the               #
#      enviornment variables at the top. FXX abd SERIAXX_COMPILER          #
#      should be the same. MPIXX_COMPILER is the name of you f90 and f77   #
#      MPI compiler.                                                       #
#                                                                          #
#   2) For OpenMP complilation, specify the proper compiler and linker     #
#      directives in the variables OPENMP_COMPILE and OPENMP_LINK, e.g.,   #
#      for ifort, this would be -openmp.                                   #
#                                                                          #
#   The target builds are:                                                 #
#                                                                          #
#   1) serial: 'make serial' builds the serial version with executable     #
#      name bigstick.x                                                     #
#                                                                          #
#   2) openmp: 'make openmp' builds the openmp version with executable     #
#      name bigstick-openmp.x                                              #
#                                                                          #
#   3) mpi: 'make mpi' builds the mpi version with executable              #
#      name bigstick-mpi.x                                                 #
#                                                                          #
#   4) mpi+openmp: 'make openmp-mpi'  builds the mpi+openmp version with   #
#      executable name bigstick-mpi-omp.x                                   #
#                                                                          #
#   In general, a simple 'make' will not work to build the serial version  #
#   as the mpif.h file may not be present, which is copied from the        #
#   file mpif.h.template. Like dummy_MPI_libs.f90, this is a set of        #
#   "dummy" set up information and routines in order to compile without    #
#   mpi activated.                                                         #
#                                                                          #
#   The current version has the ability to select the precision for        #
#   Lanczos vectors. This is done in the module Lanczos_precision          #
#   located in bmodules.f90. BIGSTICK is delivered with single-precision   #
#   as the default. To change to double precision, change the comments     #
#   in the module.                                                         #
#                                                                          #
############################################################################
F90=ifort                   #  default FORTRAN 90 compiler name
F77=ifort                   #  default FORTRAN 77 compiler name
PGI_COMPILER=pgfortran                  # Portland Group compiler
SERIAL90_COMPILER=ifort     #  serial FORTRAN 90 compiler name
SERIAL77_COMPILER=ifort     #  serial FORTRAN 77 compiler name
GFORTRAN_COMPILER =  gfortran # GNU fortran compiler
G95_COMPILER = g95 # # alternate gfortran compiler
MPI90_COMPILER=mpif90  # mpiifort     #  MPI FORTRAN 90 compiler name
MPI77_COMPILER=mpif90 #mpiifort     #  MPI FORTRAN 77 compiler name
MPI90_COMPILER_SIERRA= mpif90 # mpiifort     #  MPI FORTRAN 90 compiler name for SIERRA
MPI77_COMPILER_SIERRA= mpif90 # mpiifort     #  MPI FORTRAN 77 compiler name for SIERRA
MPI90_COMPILER_VULCAN=  mpixlf77_r  # mpigfortran  MPI FORTRAN 90 compiler name for VULCAN
MPI77_COMPILER_VULCAN=  mpixlf77_r     # mpigfortran # MPI FORTRAN 77 compiler name for VULCAN
MPI90_COMPILER_MIRA=  mpixlf90_r  # mpigfortran  MPI FORTRAN 90 compiler name for VULCAN
MPI77_COMPILER_MIRA=  mpixlf77_r     # mpigfortran # MPI FORTRAN 77 compiler name for VULCAN
EDISONF90_COMPILER=ftn     #  MPI FORTRAN 90 compiler name for EDISON
EDISONF77_COMPILER=ftn     #  MPI FORTRAN 77 compiler name for EDISON
MPI90_COMPILER_EDISON= ftn # mpiifort     #  MPI FORTRAN 90 compiler name for EDISON
MPI77_COMPILER_EDISON= ftn # mpiifort     #  MPI FORTRAN 77 compiler name for EDISON
OPT_edison="-c -fast -no-ipo"
FGCFLAGS       = "-c -O2 -std=legacy #-fbounds-check"
F77CFLAGS      =  -c -O2  #-fallow-argument-mismatch #-fbounds-check      #  Fortran 77 compiler directives
F90CFLAGS      =  -c  -O2 -std=legacy  #-fallow-argument-mismatch #-fbounds-check    #  Fortran 90 compiler directives
#INTELCFLAGS    =  -c -O2    # default Intel compiler flags
PGCFLAGS       ="-c -fast  " # PGI Fortran compiler directives
LFLAGS_all     =-llapack # -fbounds-check       #  Linker directives
LFLAGS_ifort   =          # -Wl
LFLAGS_VULCAN   =          # -Wl
LFLAGS_MIRA  =          # -Wl
OPENMP_COMPILE=-openmp      #  OpenMP compiler directive
OPENMP_LINK=-openmp         #  OpenMP linker directive
OPENMP_PGI=-mp              # OpenMP directive for Portland Group
OPENMP_G95=-openmp              # OpenMP directive for GNU compilers
OPENMP_VULCAN_COMPILE= -qsmp=omp      # -fopenmp # OpenMP compiler directive
OPENMP_VULCAN_LINK=-qsmp=omp         #  -fopenmp # OpenMP linker directive
OPENMP_MIRA_COMPILE=" -qsmp=omp -qnostrict -O3 -q64 -qtune=qp -qarch=qp -qsimd=auto "     # -fopenmp # OpenMP compiler directive
OPENMP_MIRA_LINK=-qsmp=omp         #  -fopenmp # OpenMP linker directive
OPENMP_EDISON_COMPILE=-openmp # OpenMP directive for Edison
OPENMP_CORI_COMPILE=-qopenmp # OpenMP directive for Edison
OPENMP_EDISON_LINK="-openmp -lifcoremt" # tmp directive to account for error in ftn wrapper
JAGFLAGS       =   -c -fast # compiler directives for Jaguar
GFORTRAN_MISMATCH_FLAG = -fallow-argument-mismatch #to deal with mpi integer type error with legacy code

XDENSE_OBJECTS =     \
rhmodules_main.o      \
rhmodules_parallel.o \
rhmodules_flags.o \
bmpi.o               \
rhutil.o              \
rhbasislib.o      \
rhspmap.o \
rhbasemap.o \
rhfragments.o \
rhlinalglib.o \
rhvectorlib.o \
rhwfnlib.o \
rhbasefilelib.o \
rhteamlib.o \
rhdotlib.o \
rhparity.o \
rhprojlib.o \
rhlincomlib.o \
rhjumpbase.o \
rhjumplib1b.o \
rhhoplib.o \
rhdenslib1.o \
rhdenslib2.o \
rhdenslib3.o \
rhdenslib4.o \
rhentropylib1.o \
rhspfaclib1.o \
rhmenu.o \
rhodium_main.o      \
lapackeig.o \
dlapackeig.o \
dlamc3.o \
libra.o              
 
   

%.o : %.f90
	$(F90) $(F90CFLAGS) $(MP_FLAGS) $(<F)  

%.o : %.f
	$(F77) $(F77CFLAGS) $(MP_FLAGS) $(<F)  

default: $(XDENSE_OBJECTS) $(OPENMP_OBJECTS) $(MPI_OBJECTS)
	$(F90) -o ${EX}.x $(MP_LINK) $(XDENSE_OBJECTS) $(OPENMP_OBJECTS) $(MPI_OBJECTS) $(LFLAGS) $(LIBS) $(GFORT_FLAG)

serial :
	cp mpif.h.template mpif.h
	$(MAKE) F90=$(SERIAL90_COMPILER) F77=$(SERIAL77_COMPILER) LFLAGS=$(LFLAGS_ifort) OPENMP_OBJECTS=dummy_OpenMP_libs.o MPI_OBJECTS=dummy_MPI_libs.o EX=rhodium
	rm mpif.h

pgserial :
	cp mpif.h.template mpif.h
	$(MAKE) F90=$(PGI_COMPILER) F77=$(PGI_COMPILER) F90CFLAGS=$(PGCFLAGS) F77CFLAGS=$(PGCFLAGS)  OPENMP_OBJECTS=dummy_OpenMP_libs.o MPI_OBJECTS=dummy_MPI_libs.o EX=bigstick
	rm mpif.h

gfortran :
	cp mpif.h.template mpif.h
#	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER)  MP_FLAGS=-fopenmp  MP_LINK=-fopenmp MPI_OBJECTS=dummy_MPI_libs.o EX=gbigstick
	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) LFLAGS=-llapack MPI_OBJECTS=dummy_MPI_libs.o OPENMP_OBJECTS=dummy_OpenMP_libs.o EX=rhodium
	rm mpif.h

profile :
	cp mpif.h.template mpif.h
#	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER)  MP_FLAGS=-fopenmp  MP_LINK=-fopenmp MPI_OBJECTS=dummy_MPI_libs.o EX=gbigstick
	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) LFLAGS=-llapack MPI_OBJECTS=dummy_MPI_libs.o OPENMP_OBJECTS=dummy_OpenMP_libs.o MP_FLAGS="-pg" MP_LINK="-pg" EX=rhodium-prof
	rm mpif.h
		
rhodium-mpi :
	$(MAKE) F90=$(MPI90_COMPILER) F77=$(MPI77_COMPILER) LFLAGS=-llapack OPENMP_OBJECTS=dummy_OpenMP_libs.o EX=rhodium-mpi

g95 :
	cp mpif.h.template mpif.h
#	$(MAKE) F90=$(G95_COMPILER) F77=$(G95_COMPILER) MP_FLAGS=-fopenmp MP_LINK=-fopenmp MPI_OBJECTS=dummy_MPI_libs.o EX=gbigstick
	$(MAKE) F90=$(G95_COMPILER) F77=$(G95_COMPILER) MPI_OBJECTS=dummy_MPI_libs.o OPENMP_OBJECTS=dummy_OpenMP_libs.o EX=bigstick
	rm mpif.h
sierra-serial :
	cp mpif.h.template mpif.h
	$(MAKE) F90=$(SERIAL90_COMPILER) F77=$(SERIAL77_COMPILER) LFLAGS=$(LFLAGS_ifort) OPENMP_OBJECTS=dummy_OpenMP_libs.o MPI_OBJECTS=dummy_MPI_libs.o EX=bigstick
	rm mpif.h

openmp :
	cp mpif.h.template mpif.h
	$(MAKE) F90=$(SERIAL90_COMPILER) F77=$(SERIAL77_COMPILER) LFLAGS=$(LFLAGS_ifort) MPI_OBJECTS=dummy_MPI_libs.o \
MP_FLAGS=$(OPENMP_COMPILE) MP_LINK=$(OPENMP_LINK) EX=rhodium-openmp
	rm mpif.h

pgopenmp :
	cp mpif.h.template mpif.h
	$(MAKE) F90=$(PGI_COMPILER) F77=$(PGI_COMPILER) MPI_OBJECTS=dummy_MPI_libs.o \
MP_FLAGS=$(OPENMP_PGI) MP_LINK=$(OPENMP_PGI) EX=bigstick-openmp
	rm mpif.h

g95openmp :
	cp mpif.h.template mpif.h
#	$(MAKE) F90=$(G95_COMPILER) F77=$(G95_COMPILER) MP_FLAGS=-fopenmp MP_LINK=-fopenmp MPI_OBJECTS=dummy_MPI_libs.o EX=gbigstick
	$(MAKE) F90=$(G95_COMPILER) F77=$(G95_COMPILER) MPI_OBJECTS=dummy_MPI_libs.o  MP_FLAGS=$(OPENMP_G95) EX=bigstick-openmp.x
	rm mpif.h

gfortran-openmp :
	cp mpif.h.template mpif.h
#	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) MP_FLAGS=-fopenmp MP_LINK=-fopenmp MPI_OBJECTS=dummy_MPI_libs.o EX=gbigstick
	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) MPI_OBJECTS=dummy_MPI_libs.o  MP_FLAGS=-fopenmp MP_LINK=-fopenmp EX=rhodium-openmp
	rm mpif.h

debug :
	cp mpif.h.template mpif.h
#	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) MP_FLAGS=-fopenmp MP_LINK=-fopenmp MPI_OBJECTS=dummy_MPI_libs.o EX=gbigstick
	$(MAKE) F90=$(GFORTRAN_COMPILER) F77=$(GFORTRAN_COMPILER) MPI_OBJECTS=dummy_MPI_libs.o  MP_FLAGS="-fopenmp -fbounds-check" MP_LINK="-fopenmp -fbounds-check" EX=rhodium-dbg
	rm mpif.h
	
sierra-openmp :
	cp mpif.h.template mpif.h
	$(MAKE) F90=$(SERIAL90_COMPILER) LFLAGS=$(LFLAGS_ifort) F77=$(SERIAL77_COMPILER) MPI_OBJECTS=dummy_MPI_libs.o \
MP_FLAGS=$(OPENMP_COMPILE) MP_LINK=$(OPENMP_LINK) EX=bigstick-openmp
	rm mpif.h

mpi :
	$(MAKE) F90=$(MPI90_COMPILER) F77=$(MPI77_COMPILER) OPENMP_OBJECTS=dummy_OpenMP_libs.o  EX=rhodium-mpi

sierra-mpi :
	@echo ------------------------
	@echo  Before beginning be sure to load correct tools
	@echo /usr/global/tools/dotkit/init.sh
	@echo use mvapich2-intel-1.7
	@echo -------------------------
	$(MAKE) F90=$(MPI90_COMPILER_SIERRA) F77=$(MPI77_COMPILER_SIERRA) LFLAGS=$(LFLAGS_ifort) OPENMP_OBJECTS=dummy_OpenMP_libs.o  EX=bigstick-mpi
	@echo  ----------------------------
	@echo  Note when submitting MPI jobs 
	@echo  be sure to include these lines: 
	@echo /usr/global/tools/dotkit/init.sh
	@echo use mvapich2-intel-1.7
	@echo export LD_LIBRARY_PATH=:/opt/intel-13.0/compiler/lib/intel64/:/usr/local/tools/mvapich2-intel-1.7/lib/
	@echo 	--------------------

openmp-mpi :
	$(MAKE) F90=$(MPI90_COMPILER) F77=$(MPI77_COMPILER)LFLAGS=$(LFLAGS_ifort) MP_FLAGS=$(OPENMP_COMPILE) \
MP_LINK=$(OPENMP_LINK) EX=bigstick-mpi-omp

sierra-openmp-mpi :
	$(MAKE) F90=$(MPI90_COMPILER_SIERRA) F77=$(MPI77_COMPILER_SIERRA) LFLAGS=$(LFLAGS_ifort) MP_FLAGS=$(OPENMP_COMPILE) \
MP_LINK=$(OPENMP_LINK) EX=bigstick-mpi-omp

edison-mpi :
	$(MAKE) F90=$(MPI90_COMPILER_EDISON) F77=$(MPI77_COMPILER_EDISON) F90CFLAGS=$(OPT_edison) F77CFLAGS=$(OPT_edison) OPENMP_OBJECTS=dummy_OpenMP_libs.o  EX=bigstick-mpi

edison-openmp-mpi :
	$(MAKE) F90=$(MPI90_COMPILER_EDISON) F77=$(MPI77_COMPILER_EDISON) F90CFLAGS=$(OPT_edison) F77CFLAGS=$(OPT_edison) MP_FLAGS=$(OPENMP_EDISON_COMPILE) \
MP_LINK=$(OPENMP_EDISON_LINK) EX=bigstick-mpi-omp

cori-mpi :
	$(MAKE) F90=$(MPI90_COMPILER_EDISON) F77=$(MPI77_COMPILER_EDISON) F90CFLAGS=$(OPT_edison) F77CFLAGS=$(OPT_edison) OPENMP_OBJECTS=dummy_OpenMP_libs.o  EX=bigstick-mpi

cori-openmp-mpi :
	$(MAKE) F90=$(MPI90_COMPILER_EDISON) F77=$(MPI77_COMPILER_EDISON) F90CFLAGS=$(OPT_edison) F77CFLAGS=$(OPT_edison) MP_FLAGS=$(OPENMP_CORI_COMPILE) \
MP_LINK=$(OPENMP_CORI_LINK) EX=bigstick-mpi-omp

jaguar :
	$(MAKE) F90=$(JAGF90_COMPILER) F77=$(JAGF77_COMPILER) \
OPENMP_OBJECTS=dummy_OpenMP_libs.o EX=bigstick-jaguar

vulcan-openmp-mpi :
	@echo ------------------------
	@echo  Before beginning be sure to load correct tools
	@echo  use bggcc-4.7.2
	$(MAKE) F90=$(MPI90_COMPILER_VULCAN) F77=$(MPI77_COMPILER_VULCAN) LFLAGS=$(LFLAGS_VULCAN) MP_FLAGS=$(OPENMP_VULCAN_COMPILE) \
MP_LINK=$(OPENMP_VULCAN_LINK) EX=bigstick-mpi-openmp-vulcan

mira-openmp-mpi :
	$(MAKE) F90=$(MPI90_COMPILER_MIRA) F77=$(MPI77_COMPILER_MIRA) LFLAGS=$(LFLAGS_MIRA) MP_FLAGS=$(OPENMP_MIRA_COMPILE) \
MP_LINK=$(OPENMP_MIRA_LINK) EX=bigstick-mpi-openmp

clean : 
	rm *.o *.mod *rhodium*.x mpif.h

options  :
	@echo To see a list of options, make help
help  :
	@echo
	@echo Here are some compile options:
	@echo Note default compiler is intel ifort compiler
	@echo ----------------------------------------------
	@echo       default compiler is intel ifort 
	@echo make serial -- default serial, rhodium.x
	@echo make openmp -- OpenMP parallelism,  rhodium-openmp.x
	@echo make mpi    -- MPI parallelism,  rhodium-mpi.x
	@echo make openmp-mpi -- hybrid OpenMP+MPI,  rhodium-mpi-omp.x
	@echo -----------------------------------------------
	@echo make gfortran -- serial with GNU compiler gfortran, rhodium.x
	@echo make gfortran-openmp -- OpenMP with GNU compiler gfortran, rhodium-openmp.x
	@echo make debug -- OpenM and gfortran with bounds checking, rhodium-dbg.x
#	@echo make g95      -- serial with GNU compiler g95, bigstick.x
#	@echo make g95openmp -- OpenMP with GNU compiler g95, bigstick.x NOT WORKING
#	@echo -----------------------------------------------
#	@echo make sierra -- serial with intel compiler on Sierra, bigstick.x
#	@echo make sierra-mpi -- MPI with intel compiler on Sierra, bigstick-mpi.x
#	@echo make sierra-openmp -- OpenMP with intel compiler on Sierra, bigstick-omp.x
#	@echo make sierra-openmp-mpi  -- hybrid OpenMP+MPI on Sierra,  bigstick-mpi-omp.x
#	@echo -----------------------------------------------
#	@echo make edison-mpi -- MPI with ftn compiler on Edison, bigstick-mpi.x
#	@echo make edison-openmp-mpi  -- hybrid OpenMP+MPI on Edison,  bigstick-mpi-omp.x
#	@echo -----------------------------------------------
#	@echo make cori-mpi -- MPI with ftn compiler on Cori, bigstick-mpi.x
#	@echo make cori-openmp-mpi  -- hybrid OpenMP+MPI on Cori,  bigstick-mpi-omp.x
#	@echo -----------------------------------------------
#	@echo make vulcan-openmp-mpi -- hybrid OpenMP+MPI on Vulcan,  bigstick-mpi-omp.x
#	@echo -----------------------------------------------
#	@echo make mira-openmp-mpi -- hybrid OpenMP+MPI on Mira,  bigstick-mpi-omp.x
