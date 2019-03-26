# suitable for BlueGene/Q containing Fortran XL Compiler and MPI library for BlueGene Series
# MPIF90 is a MPI wrapper for Fortran and is used to compile reaxFF code in src directory
# FC is a fortran compiler and is used to compite geninit.F90 code in init directory so as to create the input configuration for reaxFF simulation
# MPIF90_FLAGS is the flag for MPIF90
# FC_FLAGS is the  flag for FC
# LINKER is linker (optional)
# LIBS are extra needed libraries (optional) 

MPIF90 = mpif90
FC = xlf

# xl fortran
MPIF90_FLAGS = -O3 -qhot -g
#MPIF90_FLAGS = -O3 -qhot -qsmp=omp
#MPIF90_FLAGS = -C -g

FC_FLAGS = -O3 
#FC_FLAGS = -O3 -g

# LINKER = 
LIBS = /soft/perftools/hpctw/lib/libmpihpm_r.a /bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a /bgsys/drivers/ppcfloor/spi/lib/libSPI_upci_cnk.a
