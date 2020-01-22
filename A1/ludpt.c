#include "helpers.h"
#include <pthread.h>

#define TOLERANCE 0.000001
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))


struct pllop1_arg {
	int n, i, j1, j2;
    double **A;
    bool execute, terminate, ready;
};


void *pllop1(void *arguments) {
	struct pllop1_arg *args = arguments;
	
	while (!args->terminate) {
		if (args->execute) {
			int n = args->n, i = args->i, j1 = args->j1, j2 = args->j2;
			double **A = args->A;

			// printf("i: %i j1: %i, j2: %i\n", i, j1, j2);
			for (int j = j1; j < j2; ++j) {
				A[j][i] /= A[i][i];
				for (int k = i + 1; k < n; ++k) {
					// printf("%i\t%i\t%i\n", i, j, k);
					A[j][k] -= A[j][i] * A[i][k];
				}
			}
			args->execute = false;
			args->ready = true;
		}
	}
}


int LUDecompose(int n, double **A, int *Pi, int threads) {
	pthread_t tids[threads];
	struct pllop1_arg argses[threads];

	for (int t = 0; t < threads; ++t) {
		argses[t].n = n;
		argses[t].A = A;
		argses[t].execute = false;
		argses[t].terminate = false;
		argses[t].ready = true;
		if (pthread_create(&tids[t], NULL, pllop1, &argses[t])) {
			fprintf(stderr, "Error creating thread\n");
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

		for (int t = 0; t < threads; ++t) {
			argses[t].i = i;
			int per_thread = ((n - i - 1) / threads) + 1;
			argses[t].j1 = MIN((t * per_thread) + (i + 1), n);
			argses[t].j2 = MIN(((t + 1) * per_thread) + (i + 1), n);
			argses[t].execute = true;
			argses[t].ready = false;
		}
		
		while (1) {
			bool allready = true;
			for (int t = 0; t < threads; ++t) {
				allready = allready && argses[t].ready;
			}
			if (allready) { break; }
		}
	}
	for (int t = 0; t < threads; ++t) {
		argses[t].terminate = true;
	}
	for (int t = 0; t < threads; ++t) {
		pthread_join(tids[t], NULL);
	}
	return 1;
}


int main(int argc, char const *argv[]) {
	clock_t start = clock();
	int n = atoi(argv[1]), threads = atoi(argv[2]);

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
	
	randarr(n, A, Ainit);
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

	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("Time taken for LU decomposition of a %i X %i matrix: %f\n", n, n, seconds);

	FILE *ofile = fopen("matrices", "w");
	write(n, Ainit, U, L, Pi, ofile);
	fclose(ofile);
	return 0;
}
