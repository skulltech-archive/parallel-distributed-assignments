#include "helpers.h"
#include <pthread.h>
#include <semaphore.h>

// Tolerance value for detecting if a float is close enough to zero
#define TOLERANCE 0.000001
// Macro for getting the min of two numbers
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))


// The struct to be used as argument for the parallel threads
struct pllop1_arg {
	int n, i, j1, j2, t;
    double **A;
    sem_t sem;
};

// Barrier for synchronising the threads
pthread_barrier_t barrier;

// The function for the threads
void *pllop1(void *arguments) {
	struct pllop1_arg *args = arguments;

	while (1) {
		// Wait for the signal from the main thread
		sem_wait(&args->sem); 
		int n = args->n, i = args->i, j1 = args->j1, j2 = args->j2;
		double **A = args->A;
		for (int j = j1; j < j2; ++j) {
			A[j][i] /= A[i][i];
			for (int k = i + 1; k < n; ++k) {
				A[j][k] -= A[j][i] * A[i][k];
			}
		}
		// Wait for all the threads to have finished
		pthread_barrier_wait(&barrier); 
		// If the function is called by the main threads, exit
		if (!args->t) { break; }
	}
}

// The LU decompose algorithm
int LUDecompose(int n, double **A, int *Pi, int threads) {
	pthread_t tids[threads];
	struct pllop1_arg argses[threads];
	pthread_barrier_init(&barrier, NULL, threads);

	for (int t = 0; t < threads; ++t) {
		argses[t].n = n;
		argses[t].A = A;
		argses[t].t = t;
		sem_init(&argses[t].sem, 0, 0);
		// The argses[0] is for the main thread, create a thread for every other args
		if (t > 0) {
			if (pthread_create(&tids[t], NULL, pllop1, &argses[t])) {
				fprintf(stderr, "Error creating thread\n");
			}
		}
	}

	for (int i = 0; i < n; ++i) {
		Pi[i] = i;
	}

	for (int i = 0; i < n; ++i) {
		double abs, max = 0.0;
		int id = i;
		
		for (int j = i; j < n; ++j) {
			abs = fabs(A[j][i]);
			if (abs > max) {
				max = abs;
				id = j;
			}
		}

		if (max < TOLERANCE) {
			return 0;
		}

		if (id != i) {
			int temp = Pi[i];
			Pi[i] = Pi[id];
			Pi[id] = temp;

			double *ptr = A[i];
			A[i] = A[id];
			A[id] = ptr;
		}

		// Set the arguments for the different threads, and then post the semaphore signal
		for (int t = 0; t < threads; ++t) {
			argses[t].i = i;
			int per_thread = ((n - i - 1) / threads) + 1;
			argses[t].j1 = MIN((t * per_thread) + (i + 1), n);
			argses[t].j2 = MIN(((t + 1) * per_thread) + (i + 1), n);
			sem_post(&argses[t].sem);
		}
		// The main thread doesn't sit idle, it also does a slice of the work, using argses[0]
		pllop1(&argses[0]);
	}

	// Kill every child thread
	for (int t = 1; t < threads; ++t) {
		pthread_cancel(tids[t]);
	}
	return 1;
}


int main(int argc, char const *argv[]) {
	struct timespec start, finish;
	double elapsed;
	// We need system time, not CPU time
	clock_gettime(CLOCK_MONOTONIC, &start);
	int n = atoi(argv[1]), threads = atoi(argv[2]);
	bool in = false;
	const char *infile;
	if (argc > 3) {
		in = true;
		infile = argv[3];
	}

	double **A = malloc(n * sizeof(double *));
	double **Ainit = malloc(n * sizeof(double *));
	double **U = malloc(n * sizeof(double *));
	double **L = malloc(n * sizeof(double *));
	for (int i = 0; i < n; ++i) {
		A[i] = malloc(n * sizeof(double));
		Ainit[i] = malloc(n * sizeof(double));
		U[i] = malloc(n * sizeof(double));
		L[i] = malloc(n * sizeof(double));
	}
	int *Pi = malloc(sizeof(int[n]));
	
	if (in) {
		readarr(infile, n, A, Ainit);
	} else {
		randarr(n, A, Ainit);
	}
	LUDecompose(n, A, Pi, threads);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j < i) {
				U[i][j] = 0;
				L[i][j] = A[i][j];
			} else { 
				U[i][j] = A[i][j];
				L[i][j] = 0;
			}
		}
		L[i][i] = 1;
	}

	clock_gettime(CLOCK_MONOTONIC, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("Time taken for LU decomposition of a %i X %i matrix: %f\n", n, n, elapsed);

	// Write the matrices, including the resultant ones, to a file
	write(n, Ainit, U, L, Pi);
	return 0;
}
