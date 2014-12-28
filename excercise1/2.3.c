#include "mpi.h"
#include <stdio.h>

int np, rank;

int main(int argc, char **argv)
{
  int dest;
  float f;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if(rank==0) {
    dest = 1;
    while(dest < np) {    
	f = (float) rank + dest + 0.2;
	MPI_Send(&f, 1, MPI_FLOAT, dest++, 11, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(&f, 1, MPI_FLOAT, 0, 11, MPI_COMM_WORLD, &status);
    printf("Node 0 sends float %f\n", f);
  }
  MPI_Finalize();
  return 0;
}
