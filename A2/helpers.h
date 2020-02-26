#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define TOLERANCE 0.000001


void Multiply_serial(float *A, float *B, float *C, int m, int n, int p);
int IsEqual(float *A, float *B, int m, int n);
void randarr(int n, float *A);


// Function for serial multiplication of two matrices
void Multiply_serial(float *A, float *B, float *C, int m, int n, int p) {
	for (int i = 0; i < m; ++i){
		for (int j = 0; j < p; ++j){
			C[i*p + j] = 0.0;
			for (int k = 0; k < n; ++k)
				C[i*p + j] += A[i*n + k] * B[k*p + j];
		}
	}
}

// Function that checks if two matrices are equal, up to the defined tolerance
int IsEqual(float *A, float *B, int m, int n) {
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			if (fabs(A[i*m + j] - B[i*m + j]) > TOLERANCE) {
				return 0;
			}
		}
	}
	return 1;
}

// Function for creating a random array
void randarr(int n, float *A) {
	srand(time(NULL));
	for (int i = 0; i < n; ++i) {
		A[i] = (float) rand() / (float) RAND_MAX;
	}
}
