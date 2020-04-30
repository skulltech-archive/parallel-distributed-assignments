#include "helpers.h"
#include <mpi.h>


int main(int argc, char const *argv[]) {
	int m = atoi(argv[1]);
	int n = atoi(argv[2]);
	int p = atoi(argv[3]);
	double start, elapsed;

	MPI_Init(NULL, NULL);
	int world_size, world_rank, name_len;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Get_processor_name(processor_name, &name_len);
	printf("Hello world from processor %s, rank %d out of %d processors\n", processor_name, world_rank, world_size);
	
	float *A = malloc(m * n * sizeof(float));
	float *B = malloc(n * p * sizeof(float));
	float *C = malloc(m * p * sizeof(float));
	float *C_serial = malloc(m * p * sizeof(float));

	if (world_rank == 0) {
		randarr(m*n, A);
	} else if (world_rank == 1) {
		randarr(n*p, B);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Serial multiplication
	if (!world_rank) {
		start = MPI_Wtime();
		Multiply_serial(A, B, C_serial, m, n, p);
		elapsed = MPI_Wtime() - start;
		printf("[*] Serial multiplication: %f seconds\n", elapsed);
	}

	// Parallel multiplication
	int mpart = m / world_size;
	float *Apart = malloc(mpart * n * sizeof(float));
	float *Cpart = &C[world_rank * mpart * p];

	start = MPI_Wtime();
	MPI_Scatter(A, mpart * n, MPI_FLOAT, Apart, mpart * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// Muliplication part
	for (int i = 0; i < mpart; ++i) {
		for (int j = 0; j < p; ++j) {
			Cpart[i*p + j] = 0.0;
			for (int k = 0; k < n; ++k) {
				C[i*p + j] += + Apart[i*n + k] * B[k*p + j];
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	// Master node
	if (!world_rank) {
		elapsed = MPI_Wtime() - start;
		printf("[*] Parallel multiplication: %f seconds\n", elapsed);

		int correct = IsEqual(C, C_serial, m, p);
		printf("[*] Parallel correctness: %d\n", correct);
	}

	MPI_Finalize();
	return 0;
}
