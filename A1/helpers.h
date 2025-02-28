#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>


void write(int n, double **a, double **u, double **l, int *pi);
void write1d(int n, double *a, double *u, double *l, int *pi);
void randarr(int n, double **A, double **Ainit);
void randarromp(int n, double **A, double **Ainit);
void randarromp1d(int n, double *A, double *Ainit);
void readarr(const char filename[], int n, double **A, double **Ainit);
void readarr1d(const char filename[], int n, double *A, double *Ainit);


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

// Function for creating an array and a copy thereof, by reading from a file
void readarr(const char filename[], int n, double **A, double **Ainit) {
	FILE *infile;
	infile = fopen(filename, "r");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fscanf(infile, "%lf ", &A[i][j]);
		}
		fscanf(infile, "\n");
		memcpy(Ainit[i], A[i], n * sizeof(double));
	}

	fclose(infile);
}

// Parallelised version of readarr that works with 1d arrays
void readarr1d(const char filename[], int n, double *A, double *Ainit) {
	FILE *infile;
	infile = fopen(filename, "r");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fscanf(infile, "%lf ", &A[i*n + j]);
		}
		fscanf(infile, "\n");
	}
	memcpy(Ainit, A, sizeof(double[n][n]));
}

// Writes the matrices A, U, L and Pi to a file
void write(int n, double **a, double **u, double **l, int *pi) {
	FILE *permutation = fopen("Permutation", "w");
	FILE *lowertri = fopen("LowerTri", "w");
	FILE *uppertri = fopen("UpperTri", "w");
	FILE *matrix = fopen("Matrix", "w");
	for (int i = 0; i < n; ++i) {
		fprintf(permutation, "%i ", pi[i]);
		for (int j = 0; j < n; ++j) {
			fprintf(lowertri, "%lf ", l[i][j]);
			fprintf(uppertri, "%lf ", u[i][j]);
			fprintf(matrix, "%lf ", a[i][j]);
		}
		fprintf(lowertri, "\n");
		fprintf(uppertri, "\n");
		fprintf(matrix, "\n");
	}
	fclose(permutation);
	fclose(lowertri);
	fclose(uppertri);
}

// A version of the function write that works with 1d matrices
void write1d(int n, double *a, double *u, double *l, int *pi) {
	FILE *permutation = fopen("Permutation", "w");
	FILE *lowertri = fopen("LowerTri", "w");
	FILE *uppertri = fopen("UpperTri", "w");
	FILE *matrix = fopen("Matrix", "w");
	for (int i = 0; i < n; ++i) {
		fprintf(permutation, "%i ", pi[i]);
		for (int j = 0; j < n; ++j) {
			fprintf(lowertri, "%lf ", l[i*n + j]);
			fprintf(uppertri, "%lf ", u[i*n + j]);
			fprintf(matrix, "%lf ", u[i*n + j]);
		}
		fprintf(lowertri, "\n");
		fprintf(uppertri, "\n");
		fprintf(matrix, "\n");
	}
	fclose(permutation);
	fclose(lowertri);
	fclose(uppertri);
	fclose(matrix);	
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
