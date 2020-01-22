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
void randarr(int n, double **A);


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

void randarr(int n, double **A) {
	srand48(time(NULL));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[i][j] = drand48();
		}
	}
}

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

void matdif(int n, double **a, double **b, double **result) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i][j] = a[i][j] - b[i][j];
		}
	}
}

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
