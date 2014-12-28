#include "mpi.h"
#include <stdio.h>
#include <math.h>

#define N 128

int rank, np;
int length, begin, end;
MPI_Status status;
float f[N], g[N];

void init() {
	length = N / np;
	begin = rank * length;
	end = begin + length - 1;
}

void set_f() {
	int i;

	for(i=0; i<N; i++) {
		f[i] = sin(i * (1.0/N));
	}
}

void calc_g() {
	int i;

	for(i=begin; i<=end; i++) {
		g[i] = 2.0 * f[i];
	}
}

void show_g() {
	int i;

	printf("g[]:\n");
	for(i=0; i<N; i++) {
		printf(" %.2f", g[i]);
	}
	printf("/n");
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	init();
	if(rank==0) {
		set_f();
	}
	MPI_Scatter(&f[begin], length, MPI_FLOAT, &f[begin], length, MPI_FLOAT, 0, MPI_COMM_WORLD);
	calc_g();
	MPI_Gather(&f[begin], length, MPI_FLOAT, &f[begin], length, MPI_FLOAT, 0, MPI_COMM_WORLD);
	if(rank==0) {
		show_g();
	}

	MPI_Finalize();
	return 0;
}
