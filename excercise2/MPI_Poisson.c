/*
 * MPI_Poisson.c
 * 2D Poison equation solver
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DEBUG 1

#define max(a,b) ((a)>(b)?a:b)

enum
{
  X_DIR, Y_DIR
};

/* global variables */
int gridsize[2];
double precision_goal;		/* precision_goal of solution */
int max_iter;			/* maximum number of iterations alowed */
int proc_rank;

int offset[2];

int proc_coord[2];
int proc_top, proc_right, proc_bottom, proc_left;

int P;
int P_grid[2];
MPI_Comm grid_comm;
MPI_Status status;

MPI_Datatype border_type[2];

/* benchmark related variables */
clock_t ticks;			/* number of systemticks */
int timer_on = 0;		/* is timer running? */
double wtime;

/* local grid related variables */
double **phi;			/* grid */
int **source;			/* TRUE if subgrid element is a source */
int dim[2];			/* grid dimensions */

void Setup_Grid();
double Do_Step(int parity);
void Solve();
void Write_Grid();
void Clean_Up();
void Debug(char *mesg, int terminate);
void start_timer();
void resume_timer();
void stop_timer();
void print_timer();
void Setup_Proc_Grid(int argc, char **argv);
void Setup_MPI_Datatypes();
void Exchange_Borders();

void start_timer()
{
  if (!timer_on)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    ticks = clock();
    wtime = MPI_Wtime();
    timer_on = 1;
  }
}

void resume_timer()
{
  if (!timer_on)
  {
    ticks = clock() - ticks;
    wtime = MPI_Wtime() - wtime;
    timer_on = 1;
  }
}

void stop_timer()
{
  if (timer_on)
  {
    ticks = clock() - ticks;
    wtime = MPI_Wtime() - wtime;
    timer_on = 0;
  }
}

void print_timer()
{
  if (timer_on)
  {
    stop_timer();
    printf("(%i) Elapsed Wtime: %14.6f s (%5.1f%% CPU)\n", proc_rank, wtime, 100.0*ticks*(1.0/CLOCKS_PER_SEC) / wtime);
    resume_timer();
  }
  else
    printf("(%i) Elapsed Wtime: %14.6f s (%5.1f%% CPU)\n", proc_rank, wtime, 100.0*ticks*(1.0/CLOCKS_PER_SEC) / wtime);
}

void Debug(char *mesg, int terminate)
{
  if (DEBUG || terminate)
    printf("(%i) %s\n", proc_rank, mesg);
  if (terminate)
    exit(1);
}

void Setup_Grid()
{
  int x, y, s;
  double source_x, source_y, source_val;
  FILE *f;
  int upper_offset[2];

  Debug("Setup_Subgrid", 0);

  if(proc_rank==0) {
    f = fopen("input.dat", "r");
    if (f == NULL)
      Debug("Error opening input.dat", 1);
    fscanf(f, "nx: %i\n", &gridsize[X_DIR]);
    fscanf(f, "ny: %i\n", &gridsize[Y_DIR]);
    fscanf(f, "precision goal: %lf\n", &precision_goal);
    fscanf(f, "max iterations: %i\n", &max_iter);
  }

  MPI_Bcast(&gridsize, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&precision_goal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&max_iter, 1, MPI_INT, 0, MPI_COMM_WORLD);


  offset[X_DIR] = gridsize[X_DIR] * proc_coord[X_DIR] / P_grid[X_DIR];
  offset[Y_DIR] = gridsize[Y_DIR] * proc_coord[Y_DIR] / P_grid[Y_DIR];
  upper_offset[X_DIR] = gridsize[X_DIR] * (proc_coord[X_DIR] + 1) / P_grid[X_DIR];
  upper_offset[Y_DIR] = gridsize[Y_DIR] * (proc_coord[Y_DIR] + 1) / P_grid[Y_DIR];
  /* Calculate dimensions of local subgrid */
  dim[Y_DIR] = upper_offset[Y_DIR] - offset[Y_DIR];
  dim[X_DIR] = upper_offset[X_DIR] - offset[X_DIR];

  dim[Y_DIR] += 2;
  dim[X_DIR] += 2;

  if(DEBUG)
    printf("(%i) dim X: %i, dim Y: %i\n", proc_rank, dim[X_DIR], dim[Y_DIR]);



  /* allocate memory */
  if ((phi = malloc(dim[X_DIR] * sizeof(*phi))) == NULL)
    Debug("Setup_Subgrid : malloc(phi) failed", 1);
  if ((source = malloc(dim[X_DIR] * sizeof(*source))) == NULL)
    Debug("Setup_Subgrid : malloc(source) failed", 1);
  if ((phi[0] = malloc(dim[Y_DIR] * dim[X_DIR] * sizeof(**phi))) == NULL)
    Debug("Setup_Subgrid : malloc(*phi) failed", 1);
  if ((source[0] = malloc(dim[Y_DIR] * dim[X_DIR] * sizeof(**source))) == NULL)
    Debug("Setup_Subgrid : malloc(*source) failed", 1);
  for (x = 1; x < dim[X_DIR]; x++)
  {
    phi[x] = phi[0] + x * dim[Y_DIR];
    source[x] = source[0] + x * dim[Y_DIR];
  }

  /* set all values to '0' */
  for (x = 0; x < dim[X_DIR]; x++)
    for (y = 0; y < dim[Y_DIR]; y++)
    {
      phi[x][y] = 0.0;
      source[x][y] = 0;
    }

  /* put sources in field */
  do
  {
    if(proc_rank == 0)
      s = fscanf(f, "source: %lf %lf %lf\n", &source_x, &source_y, &source_val);
    
    MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (s==3)
    {
      MPI_Bcast(&source_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&source_y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&source_val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      x = source_x * gridsize[X_DIR];
      y = source_y * gridsize[Y_DIR];
      x += 1;
      y += 1;
      x = x - offset[X_DIR];
      y = y - offset[Y_DIR];
      if(x>0 && x < dim[X_DIR] - 1 && y>0 && y < dim[Y_DIR] - 1) {
        phi[x][y] = source_val;
        source[x][y] = 1;
      }
    }
  }
  while (s==3);

  if(proc_rank == 0) 
    fclose(f);
}

double Do_Step(int parity)
{
  int x, y;
  double old_phi;
  double max_err = 0.0;

  /* calculate interior of grid */
  for (x = 1; x < dim[X_DIR] - 1; x++)
    for (y = 1; y < dim[Y_DIR] - 1; y++)
      if ((x + y) % 2 == parity && source[x][y] != 1)
      {
	old_phi = phi[x][y];
	phi[x][y] = (phi[x + 1][y] + phi[x - 1][y] +
		     phi[x][y + 1] + phi[x][y - 1]) * 0.25;
	if (max_err < fabs(old_phi - phi[x][y]))
	  max_err = fabs(old_phi - phi[x][y]);
      }

  return max_err;
}

void Solve()
{
  int count = 0;
  double delta;
  double delta1, delta2;

  Debug("Solve", 0);

  /* give global_delta a higher value then precision_goal */
  delta = 2 * precision_goal;

  while (delta > precision_goal && count < max_iter)
  {
    Debug("Do_Step 0", 0);
    delta1 = Do_Step(0);
    Exchange_Borders();
    Debug("Do_Step 1", 0);
    delta2 = Do_Step(1);
    Exchange_Borders();
    delta = max(delta1, delta2);
    printf("(%i) %f\n", delta);
    count++;
  }

  printf("(%i) Number of iterations : %i\n", proc_rank, count);
}

void Write_Grid()
{
  int x, y;
  FILE *f;

  char filename[40];

  sprintf(filename, "output%i.dat", proc_rank);

  if ((f = fopen(filename, "w")) == NULL)
    Debug("Write_Grid : fopen failed", 1);

  Debug("Write_Grid", 0);

  for (x = 1; x < dim[X_DIR] - 1; x++)
    for (y = 1; y < dim[Y_DIR] - 1; y++)
      fprintf(f, "%i %i %f\n", x, y, phi[x][y]);

  fclose(f);
}

void Clean_Up()
{
  Debug("Clean_Up", 0);

  free(phi[0]);
  free(phi);
  free(source[0]);
  free(source);
}

void Setup_Proc_Grid(int argc, char **argv) {
  int wrap_around[2];
  int reorder;

  Debug("My_MPI_Init", 0);

  MPI_Comm_size(MPI_COMM_WORLD, &P);

  if(argc > 2) {
    P_grid[X_DIR] = atoi(argv[1]);
    P_grid[Y_DIR] = atoi(argv[2]);
    if(P_grid[X_DIR] * P_grid[Y_DIR] != P)
      Debug("ERROR: Process grid dimensions do not match with P ", 1);
  }
  else
    Debug("ERROR: Wrong parameter input", 1);

  wrap_around[X_DIR] = 0;
  wrap_around[Y_DIR] = 0;
  reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, 2, P_grid, wrap_around, reorder, &grid_comm);

  MPI_Comm_rank(grid_comm, &proc_rank);
  MPI_Cart_coords(grid_comm, proc_rank, 2, proc_coord);

  printf("(%i) (x,y)=(%i,%i)\n", proc_rank, proc_coord[X_DIR], proc_coord[Y_DIR]);

  MPI_Cart_shift(grid_comm, Y_DIR, 1, &proc_bottom, &proc_top);
  MPI_Cart_shift(grid_comm, X_DIR, 1, &proc_left, &proc_right);

  if(DEBUG)
    printf("(%i) top %i, right %i, bottom %i, left %i\n", proc_rank, proc_top, proc_right, proc_bottom, proc_left);
}

void Setup_MPI_Datatypes() {
  Debug("Setup_MPI_Datatypes", 0);

  printf("(%i) dim X: %i, dim Y: %i\n", proc_rank, dim[X_DIR], dim[Y_DIR]);
  //vertical
  MPI_Type_vector(dim[X_DIR] - 2, 1, dim[Y_DIR], MPI_DOUBLE, &border_type[Y_DIR]);
  MPI_Type_commit(&border_type[Y_DIR]);

  //horizontal
  MPI_Type_vector(dim[Y_DIR] - 2, 1, dim[X_DIR], MPI_DOUBLE, &border_type[X_DIR]);
  MPI_Type_commit(&border_type[X_DIR]);
}

void Exchange_Borders() {
  Debug("Exchange_Borders", 0);

/*
  //traffic in top direction
  MPI_Sendrecv(&phi[1][], 1, border_type[Y_DIR], proc_top, 0, &phi[1][], 1, border_type[Y_DIR], proc_bottom, 0, grid_comm, &status);
  //traffic in bottom direction
  MPI_Sendrecv(&phi[1][], 1, border_type[Y_DIR], proc_top, 0, &phi[1][], 1, border_type[Y_DIR], proc_bottom, 0, grid_comm, &status);
  //traffic in left direction
  MPI_Sendrecv(&phi[1][], 1, border_type[Y_DIR], proc_top, 0, &phi[1][], 1, border_type[Y_DIR], proc_bottom, 0, grid_comm, &status);
  //traffic in right direction
  MPI_Sendrecv(&phi[1][], 1, border_type[Y_DIR], proc_top, 0, &phi[1][], 1, border_type[Y_DIR], proc_bottom, 0, grid_comm, &status);
  */
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  Setup_Proc_Grid(argc, argv);
  Setup_MPI_Datatypes();

  start_timer();

  Setup_Grid();

  Solve();

  Write_Grid();

  print_timer();

  Clean_Up();
  MPI_Finalize();
  return 0;
}
