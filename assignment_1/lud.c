#include "helpers.h"


void lud(int n, bool output, bool verify) {
	clock_t start = clock();
	const char outfile[] = "matrices";
	double (*a)[n] = malloc(sizeof(double[n][n]));
	double (*u)[n] = malloc(sizeof(double[n][n]));
	double (*l)[n] = malloc(sizeof(double[n][n]));
	int *pi = malloc(sizeof(int[n]));
	srand48(time(NULL));
	FILE *ofile = fopen(outfile, "w");

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
	double (*ainit)[n] = malloc(sizeof(double[n][n]));
	memcpy(ainit, a, sizeof(double[n][n]));

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
			printf("[*] Error! Singular matrix.\n");
			return;
		}
		swapsin(&pi[k], &pi[kd]);
		swaparr(a[k], a[kd], n);
		// algo 1
		swaparr(&l[k][0], &l[kd][0], k-1);
		u[k][k] = a[k][k];

		for (int i = k+1; i < n; ++i) {
			l[i][k] = a[i][k] / u[k][k];
			u[k][i] = a[k][i];
		}
		for (int i = k+1; i < n; ++i) {
			for (int j = k+1; i < n; ++i) {
				a[i][j] -= l[i][k] * u[k][j];
			}
		}
		// algo 2
		/*for (int i = k+1; i < n; ++i) {
			a[i][k] = a[i][k] / a[k][k];
			for (int j = k+1; i < n; ++i) {
				a[i][j] -= a[i][k] * a[k][j];
			}
		}*/
		// algo 3
		/*for (int i = k; i < n; ++i) {
			for (int j = 0; j < k; ++j) {
				a[k][i] -= a[k][j] * a[j][i];
			}
		}
		for (int i = k+1; i < n; ++i) {
			for (int j = 0; j < k; ++j) {
				a[i][k] -= a[i][j] * a[j][k];
			}
			a[i][k] = a[i][k] / a[k][k];
		}*/
	}
	// algo 2 and 3
	/*for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j < i) {
				u[i][j] = 0;
				l[i][j] = a[i][j];
			} else { 
				u[i][j] = a[i][j];
				l[i][j] = 0;
			}
		}
		l[i][i] = 1;
	}*/

	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("Time taken for LU decomposition of a %i X %i matrix: %f\n", n, n, seconds);

	// Write matrices to file
	if (output) {
		write(n, ainit, u, l, pi, ofile);
	}

	// Verification part
	if (verify) {
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
	}

	fclose(ofile);
}


int main(int argc, char const *argv[]) {
	if (argc < 2) {
		printf("The first argument must be n.\n");
		exit(1);
	} else if (argc == 2) {
		lud(atoi(argv[1]), false, false);
	} else {
		bool output = (argv[2][0] == 'o' || argv[2][1] == 'o');
		bool verify = (argv[2][0] == 'v' || argv[2][1] == 'v');
		lud(atoi(argv[1]), output, verify);
	}
	return 0;
}
