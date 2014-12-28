#include "mpi.h"
#include <stdio.h>

int np, rank;

int main(int argc, char **argv) {
  	float f;
  	MPI_Status status;

  	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &np);
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  	if(rank==0) {   
		f = (float) rank + 0.2;
	}
  	
	MPI_Bcast(&f, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

	printf("Node %d received %f.\n", rank, f);
	
  	MPI_Finalize();
  	return 0;
}
