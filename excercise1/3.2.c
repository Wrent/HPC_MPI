#include "mpi.h"
#include <stdio.h>
#include <math.h>

int np, rank;

int main(int argc, char **argv) {
  	float f, g;
  	MPI_Status status;

  	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &np);
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  	if(rank==0) {   
		printf("Enter a float:\n");
		scanf("%f", &f);
	}
  	
	MPI_Bcast(&f, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	f = pow(f, rank);
	MPI_Reduce(&f, &g, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	if(rank==0) {
		printf("Reduced: %f.\n", g);
	}
	
  	MPI_Finalize();
  	return 0;
}
