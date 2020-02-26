PROCESSES = 4
N = 4000

build: mpib.c mpinb.c
	mpicc mpib.c -o mpib
	mpicc mpinb.c -o mpinb
mpib: build
	mpirun -n $(PROCESSES) ./mpib $(N) 32 $(N)
mpinb: build
	mpirun -n $(PROCESSES) ./mpinb $(N) 32 $(N)
