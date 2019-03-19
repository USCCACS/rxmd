include ../make.inc

# suffix rules
.SUFFIXES: .o .F90 .cpp
.F90.o:
	$(MPIF90) $(MPIF90_FLAGS) -c $<

%.o: %.mod

SRCS := cg.F90 pot.F90 fileio.F90 comm.F90 init.F90 bo.F90 qeq.F90 param.F90 main.F90 pqeq.F90 cmdline.F90
OBJS := $(SRCS:.F90=.o)

EXE := rxmd

TRACE := /soft/perftools/hpctw/lib/libmpihpm_r.a /bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a /bgsys/drivers/ppcfloor/spi/lib/libSPI_upci_cnk.a

../$(EXE): $(OBJS)
	$(MPIF90) $(MPIF90_FLAGS) -o ../$(EXE) $(OBJS) module.o $(LIBS)  $(TRACE)

$(OBJS) : base.mod

base.mod: module.F90
	$(MPIF90) $(MPIF90_FLAGS) -c module.F90

main.F90 : cg.mod

cg.mod : cg.o

cg.o : fileio.o

fileio.o : param.o

cmdline.o : base.mod

param.o : cmdline.o

init.o : fileio.o

main.o : init.o

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h ../rxmd rxmd

rebuild:
	make clean
	make -j
