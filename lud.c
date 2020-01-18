#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


void swaparr(double *a, double *b, int n);
void swapsin(int *a, int *b);
void status(int n, double (*a)[n], double (*u)[n], double (*l)[n], int *pi);
double l21norm(int n, double (*a)[n]);


int lud_sequential(int n) {
	double (*a)[n] = malloc(sizeof(double[n][n]));
	double (*u)[n] = malloc(sizeof(double[n][n]));
	double (*l)[n] = malloc(sizeof(double[n][n]));
	int *pi = malloc(sizeof(int[n]));
	srand48(time(NULL));

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			a[i][j] = drand48();
			if (j < i) {
				u[i][j] = 0;
				l[i][j] = drand48();
			} else { 
				u[i][j] = drand48();
				l[i][j] = 0;
			}
		}
		l[i][i] = 1;
		pi[i] = i;
	}
	status(n, a, u, l, pi);
	printf("l21 norm: %f\n", l21norm(n, a));

	double max;
	int kd;
	for (int k = 0; k < n; ++k) {
		max = 0;
		for (int i = k; i < n; ++i) {
			if (max < fabs(a[i][k])) {
				max = fabs(a[i][k]);
				kd = i;
			}
		}
		if (max == 0) {
			printf("error, singular matrix\n");
			return 1;
		}
		swapsin(&pi[k], &pi[kd]);
		swaparr(a[k], a[kd], n);
		swaparr(&l[k][0], &l[kd][0], k-1);
		u[k][k] = a[k][k];

		for (int i = k+1; i < n; ++i) {
			l[i][k] = a[i][k] / u[k][k];
			u[k][i] = a[k][i];
		}
		for (int i = k+1; i < n; ++i) {
			for (int j = k+1; i < n; ++i) {
				a[i][j] = a[i][j] - l[i][k] * u[k][j];
			}
		}
	}
	printf("\n");
	status(n, a, u, l, pi);
	printf("l21 norm: %f\n", l21norm(n, a));
	return 1;
}


int main(int argc, char const *argv[]) {
	int r = lud_sequential(4);
	return 0;
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

void status(int n, double (*a)[n], double (*u)[n], double (*l)[n], int *pi) {
	printf("a:\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			printf("%f\t", a[i][j]);
		}
		printf("\n");
	}
	printf("u:\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			printf("%f\t", u[i][j]);
		}
		printf("\n");
	}
	printf("l:\n");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			printf("%f\t", l[i][j]);
		}
		printf("\n");
	}
	printf("pi:\n");
	for (int i = 0; i < n; ++i) {
		printf("%i\t", pi[i]);
	}
	printf("\n");
}