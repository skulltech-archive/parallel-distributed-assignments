#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>


void matmul(int n, double **a, double **b, double **result);
void matdif(int n, double **a, double **b, double **result);
double l21norm(int n, double **a);
void write(int n, double **a, double **u, double **l, int *pi, FILE *ofile);
void write1d(int n, double *a, double *u, double *l, int *pi, FILE *ofile);
void randarr(int n, double **A, double **Ainit);
void randarromp(int n, double **A, double **Ainit);
void randarromp1d(int n, double *A, double *Ainit);


// Function for creating a random array and a copy thereof
void randarr(int n, double **A, double **Ainit) {
	srand48(time(NULL));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[i][j] = drand48();
		}
		memcpy(Ainit[i], A[i], n * sizeof(double));
	}
}

// Parallelised version of randarr
void randarromp(int n, double **A, double **Ainit) {
	srand48(time(NULL));
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[i][j] = drand48();
		}
		memcpy(Ainit[i], A[i], n * sizeof(double));
	}
}

// Parallelised version of randarr that works with 1d arrays
void randarromp1d(int n, double *A, double *Ainit) {
	srand48(time(NULL));
	#pragma omp parallel for schedule(static) collapse(2)
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[i*n + j] = drand48();
		}
	}
	memcpy(Ainit, A, sizeof(double[n][n]));
}

// Unused function that initializes A, Ainit, L, U and Pi, all at once
void initarrs(int n, double **A, double **Ainit, double **L, double **U, int *Pi) {
	srand48(time(NULL));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[i][j] = drand48();
			if (j < i) {
				U[i][j] = 0;
				L[i][j] = drand48();
			} else { 
				U[i][j] = drand48();
				L[i][j] = 0;
			}
		}
		L[i][i] = 1;
		Pi[i] = i;
		memcpy(Ainit[i], A[i], n * sizeof(double));
	}
}

// Function for multiplying two matrices
void matmul(int n, double **a, double **b, double **result) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i][j] = 0;
			for (int k = 0; k < n; ++k) {
				result[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}

// Calculates the difference of two matrices
void matdif(int n, double **a, double **b, double **result) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i][j] = a[i][j] - b[i][j];
		}
	}
}

// Calculates the l21 norm of a matrix
double l21norm(int n, double **a) {
	double sum = 0;
	for (int j = 0; j < n; ++j) {
		double sqrsum = 0;
		for (int i = 0; i < n; ++i) {
			sqrsum = sqrsum + (a[i][j] * a[i][j]);
		}
		sum = sum + sqrt(sqrsum);
	}
	return sum;
}

// Writes the matrices A, U, L and Pi to a file
void write(int n, double **a, double **u, double **l, int *pi, FILE *ofile) {
	fprintf(ofile, "%i\n", n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", a[i][j]);
		}
		fprintf(ofile, "\n");
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", u[i][j]);
		}
		fprintf(ofile, "\n");
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", l[i][j]);
		}
		fprintf(ofile, "\n");
	}
	for (int i = 0; i < n; ++i) {
		fprintf(ofile, "%i\t", pi[i]);
	}
	fprintf(ofile, "\n");
	
}

// A version of the function write that works with 1d matrices
void write1d(int n, double *a, double *u, double *l, int *pi, FILE *ofile) {
	fprintf(ofile, "%i\n", n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", a[i*n + j]);
		}
		fprintf(ofile, "\n");
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", u[i*n + j]);
		}
		fprintf(ofile, "\n");
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", l[i*n + j]);
		}
		fprintf(ofile, "\n");
	}
	for (int i = 0; i < n; ++i) {
		fprintf(ofile, "%i\t", pi[i]);
	}
	fprintf(ofile, "\n");
	
}

// Unused code for verifying the result of the decomposition.
/*int verify() {
	double (*p)[n] = malloc(sizeof(double[n][n]));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			p[i][j] = 0;
		}
		p[i][pi[i]] = 1;
	}
	double (*pa)[n] = malloc(sizeof(double[n][n]));
	double (*lu)[n] = malloc(sizeof(double[n][n]));
	double (*residual)[n] = malloc(sizeof(double[n][n]));
	matmul(n, p, ainit, pa);
	matmul(n, l, u, lu);
	matdif(n, pa, lu, residual);
	double norm = l21norm(n, residual);
	printf("L21 norm of the residual: %f\n", norm);

}*/
