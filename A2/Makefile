PROCESSES = 2
N = 6000

build: mpib.c mpinb.c
	mpicc mpib.c -o mpib
	mpicc mpinb.c -o mpinb
	mpicc mpic.c -o mpic
mpib: build
	mpirun -n $(PROCESSES) ./mpib $(N) 32 $(N)
mpinb: build
	mpirun -n $(PROCESSES) ./mpinb $(N) 32 $(N)
mpic: build
	mpirun -n $(PROCESSES) ./mpic $(N) 32 $(N)
