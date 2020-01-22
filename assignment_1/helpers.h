#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>


void swapsin(int *a, int *b);
void swaparr(double *a, double *b, int n);
void matmul(int n, double (*a)[n], double (*b)[n], double (*result)[n]);
void matdif(int n, double (*a)[n], double (*b)[n], double (*result)[n]);
double l21norm(int n, double (*a)[n]);
void write(int n, double (*a)[n], double (*u)[n], double (*l)[n], int *pi, FILE *outfile);
void randarr(int n, double **A);


void randarr(int n, double **A) {
	srand48(time(NULL));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[i][j] = drand48();
		}
	}
}

void randarrw(int n, double *A) {
	srand48(time(NULL));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[(i * n) + j] = drand48();
		}
	}
}

void swapsin(int *a, int *b) {
	int temp = *a;
	*a = *b;
	*b = temp;
}

void swaparr(double *a, double *b, int n) {
	double *temp = malloc(sizeof(double[n]));
	for (int i = 0; i < n; ++i) {
		temp[i] = a[i];
		a[i] = b[i];
		b[i] = temp[i];
	}
}

void matmul(int n, double (*a)[n], double (*b)[n], double (*result)[n]) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i][j] = 0;
			for (int k = 0; k < n; ++k) {
				result[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}

void matdif(int n, double (*a)[n], double (*b)[n], double (*result)[n]) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i][j] = a[i][j] - b[i][j];
		}
	}
}

double l21norm(int n, double (*a)[n]) {
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

void write(int n, double (*a)[n], double (*u)[n], double (*l)[n], int *pi, FILE *ofile) {
	fprintf(ofile, "%i\n", n);
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", a[i][j]);
		}
		fprintf(ofile, "\n");
	}
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", u[i][j]);
		}
		fprintf(ofile, "\n");
	}
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", l[i][j]);
		}
		fprintf(ofile, "\n");
	}
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		fprintf(ofile, "%i\t", pi[i]);
	}
	fprintf(ofile, "\n");
	
}

void writep(int n, double **a, double **u, double **l, int *pi, FILE *ofile) {
	fprintf(ofile, "%i\n", n);
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", a[i][j]);
		}
		fprintf(ofile, "\n");
	}
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", u[i][j]);
		}
		fprintf(ofile, "\n");
	}
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", l[i][j]);
		}
		fprintf(ofile, "\n");
	}
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		fprintf(ofile, "%i\t", pi[i]);
	}
	fprintf(ofile, "\n");
	
}

void writepw(int n, double *a, double **u, double **l, int *pi, FILE *ofile) {
	fprintf(ofile, "%i\n", n);
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", a[i *n + j]);
		}
		fprintf(ofile, "\n");
	}
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", u[i][j]);
		}
		fprintf(ofile, "\n");
	}
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(ofile, "%f\t", l[i][j]);
		}
		fprintf(ofile, "\n");
	}
	// fprintf(ofile, "\n");
	for (int i = 0; i < n; ++i) {
		fprintf(ofile, "%i\t", pi[i]);
	}
	fprintf(ofile, "\n");
	
}

