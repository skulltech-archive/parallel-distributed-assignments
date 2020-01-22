#include "helpers.h"
#define TOLERANCE 0.000001


int LUDecompose(int n, double **A, int *Pi) {
	for (int i = 0; i < n ; ++i) {
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
	int n = 4;

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
	int *P = malloc(sizeof(int[n]));
	
	randarr(n, A);
	for (int i = 0; i < n; ++i) {
		memcpy(Ainit[i], A[i], n * sizeof(double));
	}

	LUDecompose(n, A, P);

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

	FILE *ofile = fopen("matrices", "w");
	write(n, Ainit, U, L, P, ofile);
	fclose(ofile);
	return 0;
}
