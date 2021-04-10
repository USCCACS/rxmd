include ../make.inc

# suffix rules
.SUFFIXES: .o .F90 .cpp
.F90.o:
	$(MPIF90) $(MPIF90_FLAGS) -c $<

%.o: %.mod

SRCS := utils.F90 memory.F90 velocity.F90 pot.F90 fileio.F90 comm.F90 init.F90 bo.F90 qeq.F90 reaxff.F90 main.F90 pqeq.F90 cmdline.F90 fnn.F90 lists.F90 msd.F90 ensembles.F90 constant.F90 remd_pot_ANN_original.F90 shortrep.F90 cg.F90
OBJS := $(SRCS:.F90=.o)

EXE := ../rxmd

$(EXE): $(OBJS)
	$(MPIF90) $(MPIF90_FLAGS) -o $(EXE) $(OBJS) module.o $(LIBS)

cg.o : fileio.o comm.o

shortrep.o : module.o fileio.o utils.o

remd_pot_ANN_original.o : constant.o

module.o : utils.o memory.o remd_pot_ANN_original.o

fileio.o : module.o
comm.o : module.o
velocity.o : module.o

lists.o : module.o utils.o
fnn.o : module.o comm.o lists.o fileio.o shortrep.o

cmdline.o : utils.o module.o memory.o
msd.o : utils.o module.o
reaxff.o : cmdline.o utils.o fileio.o comm.o velocity.o msd.o cg.o
ensembles.o : module.o

bo.o : reaxff.o
pot.o : reaxff.o bo.o comm.o pqeq.o qeq.o
init.o : reaxff.o fileio.o pot.o fnn.o shortrep.o
qeq.o : reaxff.o comm.o lists.o
pqeq.o : reaxff.o comm.o lists.o

main.o : init.o module.o


nompi: nompi.f90 
	$(FC) $(FC_FLAGS) -DNOMPI -c nompi.f90 utils.F90 memory.F90
	make MPIF90=$(FC) MPIF90_FLAGS="$(FC_FLAGS) -DNOMPI" OBJS="$(OBJS) nompi.o" $(EXE) 

librxmd: 
	make nompi
	$(FC) $(FC_FLAGS) -DLIBRXMD -DNOMPI -c main.F90
	ar crv librxmd.a *.o

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h *.a ../rxmd

rebuild:
	make clean
	make -j
