N = 100
THREADS = 1
INFILE = B_100.txt

build: lud.c
	gcc -o lud lud.c -lm
	gcc -o ludpt ludpt.c -g -lm -lpthread
	gcc -o ludmp ludmp.c -lm -fopenmp
	gcc -o ludmp1d ludmp1d.c -lm -fopenmp
sq: build
	./lud $(N) $(INFILE)
pt: build
	./ludpt $(N) $(THREADS) $(INFILE)
mp: build
	export OMP_NUM_THREADS=$(THREADS); ./ludmp $(N) $(INFILE)
mp1d: build
	export OMP_NUM_THREADS=$(THREADS); ./ludmp1d $(N) $(INFILE)
