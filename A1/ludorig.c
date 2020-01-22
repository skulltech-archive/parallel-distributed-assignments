#include "helpers.h"
#define TOLERANCE 0.000001


int LUDecompose(int n, double **A, double **L, double **U, int *Pi) {
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

			for (int j = 0; j < i - 1; ++j) {
				temp = L[i][j];
				L[i][j] = L[id][j];
				L[id][j] = temp;
			}
		}

		U[i][i] = A[i][i];
		for (int j = i + 1; j < n; ++j) {
			L[j][i] = A[j][i] / U[i][i];
			U[i][j] = A[i][j];

		}
		for (int j = i + 1; j < n; ++j) {
			for (int k = i + 1; k < n; ++k) {
				A[j][k] -= L[j][i] * U[i][k];
			}
		}
	}
	return 1;
}


int main(int argc, char const *argv[]) {
	clock_t start = clock();
	int n = 2000;

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
	
	initarrs(n, A, Ainit, L, U, Pi);
	LUDecompose(n, A, L, U, Pi);

	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("Time taken for LU decomposition of a %i X %i matrix: %f\n", n, n, seconds);

	// FILE *ofile = fopen("matrices", "w");
	// write(n, Ainit, U, L, Pi, ofile);
	// fclose(ofile);
	// return 0;
}
