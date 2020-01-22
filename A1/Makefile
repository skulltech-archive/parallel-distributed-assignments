N = 2000
THREADS = 2

build: lud.c
	gcc -o lud lud.c -lm
	gcc -o ludpt ludpt.c -lm -lpthread
sq: build
	./lud $(N)
pt: build
	./ludpt $(N) $(THREADS)
