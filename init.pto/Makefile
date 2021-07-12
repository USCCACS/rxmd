include ../make.inc

# suffix rules
.SUFFIXES: .o .f90 .F90
.f90.o:
	$(FC)  $(FC_FLAGS) -c $<

.F90.o:
	$(FC)  $(FC_FLAGS) -c $<


TARGETS = geninit
all: $(TARGETS)

# explicit rules
geninit: geninit.o
	$(FC) $(FC_FLAGS) -o $@ $^
	./geninit
	cp -v rxff.bin ../DAT

nocopy: geninit.o
	$(FC) $(FC_FLAGS) -o $@ $^
	./geninit

rebuild:
	make clean
	make

clean:
	rm -f all.bin rxff.bin ../DAT/rxff.bin geninit.xyz real.xyz norm.xyz *.o *.mod *.MOD $(TARGETS)
