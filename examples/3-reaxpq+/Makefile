run:
	make set
	mpirun -np 1 ../../rxmd 

set:
	../../init/geninit -i input.xyz -f ffield -o DAT/ -mc 2 3 5 

clean:
	rm -fv geninit* rxmd DAT/*
