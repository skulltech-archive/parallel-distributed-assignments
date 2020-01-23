#include "helpers.h"
#include <omp.h>

#define TOLERANCE 0.000001


int LUDecompose(int n, double *A, int *Pi) {
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; ++i) {
		Pi[i] = i;
	}

	for (int i = 0; i < n; ++i) {
		double abs, max = 0.0;
		int id = i;

		for (int j = i; j < n; ++j) {
			abs = fabs(A[j*n + i]);
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

			double tmp;
			for (int j = 0; j < n; ++j) {
				double tmp = A[i*n + j];
				A[i*n + j] = A[id*n + j];
				A[id*n + j] = tmp;
			}
		}

		#pragma omp parallel for schedule(static)
		for (int j = i + 1; j < n; ++j) {
			A[j*n + i] /= A[i*n + i];

			for (int k = i + 1; k < n; ++k) {
				A[j*n + k] -= A[j*n + i] * A[i*n + k];
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

	double *A = malloc(sizeof(double[n][n]));
	double *Ainit = malloc(sizeof(double[n][n]));
	double *U = malloc(sizeof(double[n][n]));
	double *L = malloc(sizeof(double[n][n]));
	int *Pi = malloc(sizeof(int[n]));
	
	randarromp1d(n, A, Ainit);
	LUDecompose(n, A, Pi);

	#pragma omp parallel for schedule(static) collapse(2)
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			L[i*n + i] = 1;
			if (j < i) {
				U[i*n + j] = 0;
				L[i*n + j] = A[i*n + j];
			} else {
				U[i*n + j] = A[i*n + j];
				L[i*n + j] = 0;
			}
		}
	}

	clock_gettime(CLOCK_MONOTONIC, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("Time taken for LU decomposition of a %i X %i matrix: %f\n", n, n, elapsed);

	FILE *ofile = fopen("matrices", "w");
	write1d(n, Ainit, U, L, Pi, ofile);
	fclose(ofile);
	return 0;
}
