#include "helpers.h"
#define TOLERANCE 0.000001


int LUDecompose(int n, double **A, int *Pi) {
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

		for (int j = i + 1; j < n; ++j) {
			A[j][i] /= A[i][i];

			for (int k = i + 1; k < n; ++k) {
				A[j][k] -= A[j][i] * A[i][k];
			}
		}
	}
	return 1;
}


int main(int argc, char const *argv[]) {
	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_MONOTONIC, &start);
	int n = atoi(argv[1]);
	bool in = false;
	const char *infile;
	if (argc > 2) {
		in = true;
		infile = argv[2];
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
	LUDecompose(n, A, Pi);

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
