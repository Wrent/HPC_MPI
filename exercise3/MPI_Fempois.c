/*
 * MPI_Fempois.c
 * 2D Poisson equation solver with MPI and FEM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

#define DEBUG 0

#define TYPE_GHOST 1
#define TYPE_SOURCE 2

#define MAXCOL 20

typedef struct
{
  int type;
  double x, y;
}
Vertex;

typedef int Element[3];

typedef struct
{
  int Ncol;
  int *col;
  double *val;
}
Matrixrow;

/* global variables */
double precision_goal;		/* precision_goal of solution */
int max_iter;			/* maximum number of iterations alowed */
int P;				/* total number of processes */
int P_grid[2];			/* processgrid dimensions */
MPI_Comm grid_comm;		/* grid COMMUNICATOR */
MPI_Status status;

/* benchmark related variables */
clock_t ticks;			/* number of systemticks */
double wtime;			/* wallclock time */
int timer_on = 0;		/* is timer running? */
double computation_time = 0;
double neighbours_time = 0;
double global_time = 0;
double idle_time = 0;

/* local process related variables */
int proc_rank;			/* rank of current process */
int proc_coord[2];		/* coordinates of current procces in processgrid */
int N_neighb;			/* Number of neighbouring processes */
int *proc_neighb;		/* ranks of neighbouring processes */
MPI_Datatype *send_type;	/* MPI Datatypes for sending */
MPI_Datatype *recv_type;	/* MPI Datatypes for receiving */

/* local grid related variables */
Vertex *vert;			/* vertices */
double *phi;			/* vertex values */
int N_vert;			/* number of vertices */
Matrixrow *A;			/* matrix A */

void Setup_Proc_Grid();
void Setup_Grid();
void Build_ElMatrix(Element el);
void Sort_MPI_Datatypes();
void Setup_MPI_Datatypes(FILE *f);
void Exchange_Borders(double *vect);
void Solve();
void Write_Grid();
void Clean_Up();
void Debug(char *mesg, int terminate);
void start_timer();
void resume_timer();
void stop_timer();
void print_timer();

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
    ticks = clock()-ticks;
    wtime = MPI_Wtime()-wtime;
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
  int timer_was_on = 0;
  if (timer_on)
  {
    timer_was_on = 1;
    stop_timer();
  }
    printf("(%i) Elapsed Wtime: %14.6f s (%5.1f%% CPU)\n",
     proc_rank, wtime, 100.0 * ticks * (1.0 / CLOCKS_PER_SEC) / wtime);

    printf("(%i) Computation Wtime: %14.6f s\n", proc_rank, computation_time);
    printf("(%i) Neighbours Wtime: %14.6f s\n", proc_rank, neighbours_time);
    printf("(%i) Global Wtime: %14.6f s\n", proc_rank, global_time);
    printf("(%i) Idle Wtime: %14.6f s\n", proc_rank, idle_time);
    double global_comp;
    double global_nei;
    double global_comm;
    MPI_Allreduce(&computation_time, &global_comp, 1, MPI_DOUBLE, MPI_SUM, grid_comm);
    MPI_Allreduce(&neighbours_time, &global_nei, 1, MPI_DOUBLE, MPI_SUM, grid_comm);
    MPI_Allreduce(&global_time, &global_comm, 1, MPI_DOUBLE, MPI_SUM, grid_comm);
    printf("(%i) cmp vs comm: %14.6f %14.6f s\n", proc_rank, global_comp, global_nei + global_comm);
    if (timer_was_on)
  {
    resume_timer();
  }
}

void Debug(char *mesg, int terminate)
{
  if (DEBUG || terminate)
    printf("(%i) %s\n", proc_rank, mesg);
  if (terminate)
  {
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }
}

void Setup_Proc_Grid()
{
  FILE *f = NULL;
  char filename[25];
  int i;
  int N_nodes = 0, N_edges = 0;
  int *index, *edges, reorder;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
  Debug("My_MPI_Init", 0);

  /* Retrieve the number of processes and current process rank */
  MPI_Comm_size(MPI_COMM_WORLD, &P);

  /* Create process topology (Graph) */
  if (proc_rank == 0)
  {
    sprintf(filename, "mapping%i.dat", P);
    if ((f = fopen(filename, "r")) == NULL)
      Debug("My_MPI_Init : Can't open mapping inputfile", 1);

    /* after reading N_nodes, a line is skipped */
    fscanf(f, "N_proc : %i\n%*[^\n]\n", &N_nodes);
    if (N_nodes != P)
      Debug("My_MPI_Init : Mismatch of number of processes in mapping inputfile", 1);
  }
  else
    N_nodes = P;

  if ((index = malloc(N_nodes * sizeof(int))) == NULL)
      Debug("My_MPI_Init : malloc(index) failed", 1);

  if (proc_rank == 0)
  {
    for (i = 0; i < N_nodes; i++)
      fscanf(f, "%i\n", &index[i]);
  }

  MPI_Bcast(index, N_nodes, MPI_INT, 0, MPI_COMM_WORLD);

  N_edges = index[N_nodes - 1];
  if (N_edges>0)
  {
    if ((edges = malloc(N_edges * sizeof(int))) == NULL)
      Debug("My_MPI_Init : malloc(edges) failed", 1);
  }
  else
    edges = index; /* this is actually nonsense,
                      but 'edges' needs to be a non-null pointer */

  if (proc_rank == 0)
  {
    fscanf(f, "%*[^\n]\n");		/* skip a line of the file */
    for (i = 0; i < N_edges; i++)
      fscanf(f, "%i\n", &edges[i]);

    fclose(f);
  }

  MPI_Bcast(edges, N_edges, MPI_INT, 0, MPI_COMM_WORLD);

  reorder = 1;
  MPI_Graph_create(MPI_COMM_WORLD, N_nodes, index, edges, reorder, &grid_comm);

  /* Retrieve new rank of this process */
  MPI_Comm_rank(grid_comm, &proc_rank);

  if (N_edges>0)
    free(edges);
  free(index);
}

void Setup_Grid()
{
  int i, j, v;
  Element element;
  int N_elm;
  char filename[25];
  FILE *f;

  Debug("Setup_Grid", 0);

  /* read general parameters (precision/max_iter) */
  if (proc_rank==0)
  {
    if ((f = fopen("input.dat", "r")) == NULL)
      Debug("Setup_Grid : Can't open input.dat", 1);
    fscanf(f, "precision goal: %lf\n", &precision_goal);
    fscanf(f, "max iterations: %i", &max_iter);
    fclose(f);
  }
  MPI_Bcast(&precision_goal, 1, MPI_DOUBLE, 0, grid_comm);
  MPI_Bcast(&max_iter, 1, MPI_INT, 0, grid_comm);

  /* read process specific data */
  sprintf(filename, "input%i-%i.dat", P, proc_rank);
  if ((f = fopen(filename, "r")) == NULL)
    Debug("Setup_Grid : Can't open data inputfile", 1);
  fscanf(f, "N_vert: %i\n%*[^\n]\n", &N_vert);

  /* allocate memory for phi and A */
  if ((vert = malloc(N_vert * sizeof(Vertex))) == NULL)
    Debug("Setup_Grid : malloc(vert) failed", 1);
  if ((phi = malloc(N_vert * sizeof(double))) == NULL)
    Debug("Setup_Grid : malloc(phi) failed", 1);

  if ((A = malloc(N_vert * sizeof(*A))) == NULL)
    Debug("Setup_Grid : malloc(*A) failed", 1);
  for (i=0; i<N_vert; i++)
  {
    if ((A[i].col=malloc(MAXCOL*sizeof(int)))==NULL)
      Debug("Setup_Grid : malloc(A.col) failed", 1);
    if ((A[i].val=malloc(MAXCOL*sizeof(double)))==NULL)
      Debug("Setup_Grid : malloc(A.val) failed", 1);
  }

  /* init matrix rows of A */
  for (i = 0; i < N_vert; i++)
      A[i].Ncol = 0;

  /* Read all values */
  for (i = 0; i < N_vert; i++)
  {
    fscanf(f, "%i", &v);
    fscanf(f, "%lf %lf %i %lf\n", &vert[v].x, &vert[v].y,
	   &vert[v].type, &phi[v]);
  }

  /* build matrix from elements */
  fscanf(f, "N_elm: %i\n%*[^\n]\n", &N_elm);
  for (i = 0; i < N_elm; i++)
  {
    fscanf(f, "%*i");  /* we are not interested in the element-id */
    for (j = 0; j < 3; j++)
    {
      fscanf(f, "%i", &v);
      element[j] = v;
    }
    fscanf(f, "\n");
    Build_ElMatrix(element);
  }

  Setup_MPI_Datatypes(f);

  fclose(f);
}

void Add_To_Matrix(int i, int j, double a)
{
  int k;
  k=0;
  
  while ( (k<A[i].Ncol) && (A[i].col[k]!=j) )
    k++;
  if (k<A[i].Ncol)
    A[i].val[k]+=a;
  else
  {
    if (A[i].Ncol>=MAXCOL)
      Debug("Add_To_Matrix : MAXCOL exceeded", 1);
    A[i].val[A[i].Ncol]=a;
    A[i].col[A[i].Ncol]=j;
    A[i].Ncol++;
  }
}

void Build_ElMatrix(Element el)
{
  int i, j;
  double e[3][2];
  double s[3][3];
  double det;

  e[0][0] = vert[el[1]].y - vert[el[2]].y;	/* y1-y2 */
  e[1][0] = vert[el[2]].y - vert[el[0]].y;	/* y2-y0 */
  e[2][0] = vert[el[0]].y - vert[el[1]].y;	/* y0-y1 */
  e[0][1] = vert[el[2]].x - vert[el[1]].x;	/* x2-x1 */
  e[1][1] = vert[el[0]].x - vert[el[2]].x;	/* x0-x2 */
  e[2][1] = vert[el[1]].x - vert[el[0]].x;	/* x1-x0 */

  det = e[2][0] * e[0][1] - e[2][1] * e[0][0];
  if (det == 0.0)
    Debug("One of the elements has a zero surface", 1);

  det = fabs(2 * det);

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      s[i][j] = (e[i][0] * e[j][0] + e[i][1] * e[j][1]) / det;

  for (i = 0; i < 3; i++)
    if (!((vert[el[i]].type & TYPE_GHOST) |
	  (vert[el[i]].type & TYPE_SOURCE)))
      for (j = 0; j < 3; j++)
        Add_To_Matrix(el[i],el[j],s[i][j]);
}

void Sort_MPI_Datatypes()
{
  int i, j;
  MPI_Datatype data2;
  int proc2;

  for (i=0;i<N_neighb-1;i++)
    for (j=i+1;j<N_neighb;j++)
      if (proc_neighb[j]<proc_neighb[i])
      {
        proc2 = proc_neighb[i];
        proc_neighb[i] = proc_neighb[j]; 
        proc_neighb[j] = proc2;
        data2 = send_type[i];
        send_type[i] = send_type[j];
        send_type[j] = data2;
        data2 = recv_type[i];
        recv_type[i] = recv_type[j];
        recv_type[j] = data2;
      }
}

int sum_count = 0;

void Setup_MPI_Datatypes(FILE * f)
{
  int i, s;
  int count;
  int *indices;
  int *blocklens;

  Debug("Setup_MPI_Datatypes", 0);

  fscanf(f, "Neighbours: %i\n", &N_neighb);

  /* allocate memory */

  if (N_neighb>0)
  {
    if ((proc_neighb = malloc(N_neighb * sizeof(int))) == NULL)
        Debug("Setup_MPI_Datatypes: malloc(proc_neighb) failed", 1);
    if ((send_type = malloc(N_neighb * sizeof(MPI_Datatype))) == NULL)
      Debug("Setup_MPI_Datatypes: malloc(send_type) failed", 1);
    if ((recv_type = malloc(N_neighb * sizeof(MPI_Datatype))) == NULL)
      Debug("Setup_MPI_Datatypes: malloc(recv_type) failed", 1);
  }
  else
  {
    proc_neighb = NULL;
    send_type = NULL;
    recv_type = NULL;
  }

  if ((indices = malloc(N_vert * sizeof(int))) == NULL)
      Debug("Setup_MPI_Datatypes: malloc(indices) failed", 1);
  if ((blocklens = malloc(N_vert * sizeof(int))) == NULL)
      Debug("Setup_MPI_Datatypes: malloc(blocklens) failed", 1);

  for (i = 0; i < N_vert; i++)
    blocklens[i] = 1;

  /* read vertices per neighbour */
  for (i = 0; i < N_neighb; i++)
  {
    fscanf(f, "from %i :", &proc_neighb[i]);
    s = 1;
    count = 0;
    while (s == 1)
    {
      s = fscanf(f, "%i", &indices[count]);
      if ((s == 1) && !(vert[indices[count]].type & TYPE_SOURCE))
	count++;
    }
    fscanf(f, "\n");
    MPI_Type_indexed(count, blocklens, indices, MPI_DOUBLE, &recv_type[i]);
    MPI_Type_commit(&recv_type[i]);

    fscanf(f, "to %i :", &proc_neighb[i]);
    s = 1;
    sum_count += count;
    count = 0;
    while (s == 1)
    {
      s = fscanf(f, "%i", &indices[count]);
      if ((s == 1) && !(vert[indices[count]].type & TYPE_SOURCE))
	count++;
    }
    fscanf(f, "\n");
    MPI_Type_indexed(count, blocklens, indices, MPI_DOUBLE, &send_type[i]);
    MPI_Type_commit(&send_type[i]);
    sum_count += count;
  }

  Sort_MPI_Datatypes();

  free(blocklens);
  free(indices);
}

int data_sent = 0;

void Exchange_Borders(double *vect)
{
  int i;
  data_sent += N_neighb * sizeof(MPI_Datatype);
  Debug("Exchange_Borders", 0);

  for (i = 0; i < N_neighb; i++)
    MPI_Sendrecv(vect, 1, send_type[i], proc_neighb[i], 0,
		 vect, 1, recv_type[i], proc_neighb[i], 0,
		 grid_comm, &status);
}


void Solve()
{
  int count = 0;
  int i, j;
  double *r, *p, *q;
  double a, b, r1, r2 = 1;

  double sub;

  Debug("Solve", 0);

  double start_time = wtime;


  if ((r = malloc(N_vert * sizeof(double))) == NULL)
      Debug("Solve : malloc(r) failed", 1);
  if ((p = malloc(N_vert * sizeof(double))) == NULL)
      Debug("Solve : malloc(p) failed", 1);
  if ((q = malloc(N_vert * sizeof(double))) == NULL)
      Debug("Solve : malloc(q) failed", 1);

  /* Implementation of the CG algorithm : */
  //neighbours phase
  stop_timer();
  start_time = wtime;
  resume_timer();
  MPI_Barrier(grid_comm);
  stop_timer();
  idle_time += wtime - start_time;
  start_time = wtime;
  resume_timer();
  Exchange_Borders(phi);
  stop_timer();
  neighbours_time += wtime - start_time;
  start_time = wtime;
  resume_timer();

  //computation phase
  /* r = b-Ax */
  for (i = 0; i < N_vert; i++)
  {
    r[i] = 0.0;
    for (j = 0; j < A[i].Ncol; j++)
      r[i] -= A[i].val[j] * phi[A[i].col[j]];
  }

  r1 = 2 * precision_goal;
  while ((count < max_iter) && (r1 > precision_goal))
    {
      /* r1 = r' * r */
      sub = 0.0;
      for (i = 0; i < N_vert; i++)
        if (!(vert[i].type & TYPE_GHOST))
  	sub += r[i] * r[i];

        stop_timer();
        computation_time += wtime - start_time;
        start_time = wtime;
        resume_timer();
  MPI_Barrier(grid_comm);
  stop_timer();
  idle_time += wtime - start_time;
  start_time = wtime;
  resume_timer();
  //global communicication
      MPI_Allreduce(&sub, &r1, 1, MPI_DOUBLE, MPI_SUM, grid_comm);
      stop_timer();
      global_time += wtime - start_time;
      start_time = wtime;
      resume_timer();
  //computation
      if (count == 0)
      {
        /* p = r */
        for (i = 0; i < N_vert; i++)
  	p[i] = r[i];
      }
      else
      {
        b = r1 / r2;

        /* p = r + b*p */
        for (i = 0; i < N_vert; i++)
  	p[i] = r[i] + b * p[i];
      }
      stop_timer();
        computation_time += wtime - start_time;
        start_time = wtime;
        resume_timer();
      //neihgbours
          MPI_Barrier(grid_comm);
  stop_timer();
  idle_time += wtime - start_time;
  start_time = wtime;
  resume_timer();
      Exchange_Borders(p);
      stop_timer();
        neighbours_time += wtime - start_time;
        start_time = wtime;
        resume_timer();
  //computation
      /* q = A * p */
      for (i = 0; i < N_vert; i++)
      {
        q[i] = 0;
        for (j = 0; j < A[i].Ncol; j++)
          q[i] += A[i].val[j] * p[A[i].col[j]];
      }

      /* a = r1 / (p' * q) */
      sub = 0.0;
      for (i = 0; i < N_vert; i++)
        if (!(vert[i].type & TYPE_GHOST))
  	sub += p[i] * q[i];
  //global
     stop_timer();
        computation_time += wtime - start_time;
        start_time = wtime;
        resume_timer();
          MPI_Barrier(grid_comm);
  stop_timer();
  idle_time += wtime - start_time;
  start_time = wtime;
  resume_timer();
      MPI_Allreduce(&sub, &a, 1, MPI_DOUBLE, MPI_SUM, grid_comm);
      stop_timer();
        global_time += wtime - start_time;
        start_time = wtime;
        resume_timer();
      a = r1 / a;

  //computation
      /* x = x + a*p */
      for (i = 0; i < N_vert; i++)
        phi[i] += a * p[i];

      /* r = r - a*q */
      for (i = 0; i < N_vert; i++)
        r[i] -= a * q[i];

      r2 = r1;

      count++;
    }
  free(q);
  free(p);
  free(r);

  if (proc_rank == 0)
    printf("Number of iterations : %i\n", count);
}

void Write_Grid()
{
  int i;
  char filename[25];
  FILE *f;

  Debug("Write_Grid", 0);

  sprintf(filename, "output%i-%i.dat", P, proc_rank);
  if ((f = fopen(filename, "w")) == NULL)
    Debug("Write_Grid : Can't open data outputfile", 1);

  for (i = 0; i < N_vert; i++)
    if (!(vert[i].type & TYPE_GHOST))
      fprintf(f, "%f %f %f\n", vert[i].x, vert[i].y, phi[i]);

  fclose(f);
}

void Clean_Up()
{
  int i;
  Debug("Clean_Up", 0);

  if (N_neighb>0)
  {
    free(recv_type);
    free(send_type);
    free(proc_neighb);
  }

  for (i=0;i<N_vert;i++)
  {
    free(A[i].col);
    free(A[i].val);
  }
  free(A);
  free(vert);
  free(phi);
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  start_timer();
  stop_timer();
  Setup_Proc_Grid();

  Setup_Grid();

  Solve();

  Write_Grid();

  Clean_Up();

  print_timer();

  Debug("MPI_Finalize", 0);
  printf("(%i) %i\n", proc_rank, sum_count);
  printf("(%i) %i\n", proc_rank, data_sent);

  MPI_Finalize();

  return 0;
}

