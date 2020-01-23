N = 8000
THREADS = 1

build: lud.c
	gcc -o lud lud.c -lm
	gcc -o ludpt ludpt.c -lm -lpthread
	gcc -o ludmp ludmp.c -lm -fopenmp
sq: build
	./lud $(N)
pt: build
	./ludpt $(N) $(THREADS)
mp: build
	export OMP_NUM_THREADS=$(THREADS); ./ludmp $(N)
mp_default: build
	unset OMP_NUM_THREADS; ./ludmp $(N)