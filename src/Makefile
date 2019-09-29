include ../make.inc

# suffix rules
.SUFFIXES: .o .F90 .cpp
.F90.o:
	$(MPIF90) $(MPIF90_FLAGS) -c $<

%.o: %.mod

SRCS := utils.F90 memory.F90 velocity.F90 pot.F90 fileio.F90 comm.F90 init.F90 bo.F90 qeq.F90 reaxff.F90 main.F90 pqeq.F90 cmdline.F90 fnn.F90 lists.F90 msd.F90 stats.F90
OBJS := $(SRCS:.F90=.o)

EXE := ../rxmd

$(EXE): $(OBJS)
	$(MPIF90) $(MPIF90_FLAGS) -o $(EXE) $(OBJS) module.o $(LIBS)

#cg.o : fileio.o
#main.F90 : cg.o

module.o : utils.o memory.o

fileio.o : module.o
comm.o : module.o
velocity.o : module.o

lists.o : module.o utils.o
fnn.o : module.o comm.o lists.o fileio.o

cmdline.o : utils.o module.o memory.o
msd.o : utils.o module.o
stats.o: utils.o lists.o module.o
reaxff.o : cmdline.o utils.o fileio.o comm.o velocity.o msd.o stats.o

bo.o : reaxff.o
pot.o : reaxff.o bo.o comm.o pqeq.o qeq.o
init.o : reaxff.o fileio.o pot.o fnn.o
qeq.o : reaxff.o comm.o lists.o
pqeq.o : reaxff.o comm.o lists.o

main.o : init.o


nompi: nompi.f90
	$(FC) $(FC_FLAGS) -DNOMPI -c nompi.f90 utils.F90
	make MPIF90=$(FC) MPIF90_FLAGS="$(FC_FLAGS) -DNOMPI" OBJS="$(OBJS) nompi.o" $(EXE) 

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h ../rxmd

rebuild:
	make clean
	make -j
