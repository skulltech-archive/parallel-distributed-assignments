#include "helpers.h"




int LUPDecompose(double **A, int N, double Tol, int *P) {
	int i, j, k, imax;
	double maxA, *ptr, absA;

	for (i = 0; i <=N ; ++i) {
		P[i] = i;
	}
	for (i = 0; i < N; i++) {
		maxA = 0.0;
		imax = i;

		for (k = i; k < N; k++) {
			if ((absA = fabs(A[k][i])) > maxA) {
				maxA = absA;
				imax = k;
			}
		}

		if (maxA < Tol) {
			return 0;
		}

		if (imax != i) {
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			ptr = A[i];
			A[i] = A[imax];
			A[imax] = ptr;

			P[N]++;
		}

		for (j = i + 1; j < N; j++) {
			A[j][i] /= A[i][i];

			for (k = i +1; k < N; k++) {
				A[j][k] -= A[j][i] * A[i][k];
			}
		}
	}
	return 1;
}

int main(int argc, char const *argv[])
{
	int n = 4;

	double **A = malloc(n * sizeof(double *));
	double **Ainit = malloc(n * sizeof(double *));
	// double *A = malloc(sizeof(double[n][n]));
	// double *Ainit = malloc(sizeof(double[n][n]));
	double **U = malloc(n * sizeof(double *));
	double **L = malloc(n * sizeof(double *));
	for (int i = 0; i < n; ++i) {
		A[i] = malloc(n * sizeof(double));
		Ainit[i] = malloc(n * sizeof(double));
		U[i] = malloc(n * sizeof(double));
		L[i] = malloc(n * sizeof(double));
	}
	int *P = malloc(sizeof(int[n+1]));
	double Tol = 0.000001;
	
	randarr(n, A);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			Ainit[i][j] = A[i][j];
		}
	}
	printf("here\n");
	// randarrw(n, A);
	// memcpy(Ainit, A, sizeof(double[n][n]));

	LUPDecompose(A, n, Tol, P);
	// Doolittle_LU_Decomposition_with_Pivoting(A, P, n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j < i) {
				U[i][j] = 0;
				// L[i][j] = A[i*n + j];
				L[i][j] = A[i][j];
			} else { 
				// U[i][j] = A[i* n + j];
				U[i][j] = A[i][j];
				L[i][j] = 0;
			}
		}
		L[i][i] = 1;
	}


	FILE *ofile = fopen("matrices", "w");
	writep(n, Ainit, U, L, P, ofile);
	// writepw(n, Ainit, U, L, P, ofile);
	fclose(ofile);
	return 0;
}
