#include "mpi.h"
#include <stdio.h>

int np, rank;

int main(int argc, char **argv)
{
  int source;
  float f;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if(rank==0) {
    source = 1;
    while(source < np) {    
	MPI_Recv(&f, 1, MPI_FLOAT, source, 11, MPI_COMM_WORLD, &status);
    	printf("Node %d sends float %f\n", source++, f);
    }
  } else {
    f = (float) rank + 0.2;
    MPI_Send(&f, 1, MPI_FLOAT, 0, 11, MPI_COMM_WORLD);
  }
  MPI_Finalize();
  return 0;
}
