/***
 * Parallel code for NBODY problem using the cell-algoritm and MPI
 *
 * Note: First look for `TODO' comments if the code doesn't work.
 *
 * Version 1.7, by Sander Flobbe and Fons van Hees, august 19, 1999.
 ***/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

#define DEBUG 0
#define max(a,b) ((a)>(b)?a:b)
#define min(a,b) ((a)<(b)?a:b)
#define RANMAX 2147483648


/* Main parameters                                                          */
double sizeX, sizeY;      /* Size of the global domain                      */
int nPartX, nPartY;       /* Particle number in x, y direction              */
int nPart;                /* Total number of particles                      */
int nCellX, nCellY;       /* (Global) number of cells in x, y direction     */
int steps;                /* Number of timesteps                            */
double dt;                /* Stepsize for timesteps                         */
int logs;                 /* Whether or not we want to keep logfiles        */


/* Cells, global and local domain                                           */
double beginx,  beginy;    /* Begin of domain in x, y direction             */
double endx,    endy;      /* End of domain in x, y direction               */
int    startcx, startcy;   /* Cell-begin in x, y direction                  */
int    nextcx,  nextcy;    /* End of cell-domain in x, y direction          */
int    procNcx, procNcy;   /* Number of internal cells on this process      */


/* Intermediate results for calculations                                    */
double sizeCellX, sizeCellY; /* length of a cell in x, y direction          */
double dt2;


/* MPI related variables                                                    */
MPI_Comm grid_comm;       /* Communicator for process-grid                  */
MPI_Status status;        /* Status                                         */
int nProc;                /* Number of processes                            */
int nProcX, nProcY;       /* Number of processes in x, y direction          */
int px, py;               /* Processor coordinates (in grid)                */
int rank;                 /* Rank of this process                           */


/* Particle numbers                                                         */
int nPartProc = 0;        /* Number of particles on this process            */
int maxNpart;             /* Maximum number of particles on this process    */


/* Benchmarking variables                                                   */
clock_t ticks;            /* Number of clockticks (cpu-time)                */
double wtime;             /* Wall-clock time                                */
int timer_on = 0;         /* Whether or not the timer is running            */


/* Logfile names                                                            */
const char *energyLog;    /* Track time, eLJ, eKin, eTotal                  */
const char *trajectLog;   /* Track one particle (x,y)                       */
const char *particleDump; /* Track all particles (x,y)                      */


typedef struct
{
  double x, y;
  double oldx, oldy;
  double fx, fy;
  int cx, cy;             /* Cell-indices                                   */
  int localindex;         /* Local index of particle in a cell's list       */
} particleType;


typedef struct
{
  int nParticles;
  particleType **particles;
} cellType;


typedef struct
{
  int rank;
  int dx, dy;
  int nCells;
  int sendStartCx, sendStartCy;
  int recvStartCx, recvStartCy;
  int *sendNumList, *recvNumList;
  int nPartSend, nPartReceived;
  particleType **sendParticles, **recvParticles;
  double (*sendList)[2], (*recvList)[2];
} neighbourType;


particleType    *particle;
cellType        **cell;
neighbourType   neighbour[8];


/* Print a debug line and terminate the program if necessary                */
void Debug(char *mesg, int terminate)
{
  if(DEBUG || terminate)
    fprintf(stderr, "(%i) %s\n", rank, mesg);
  if(terminate)
  {
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }
}


void StartTimer()
{
  if(!timer_on)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    ticks = clock();
    wtime = MPI_Wtime();
    timer_on = 1;
  }
}


void ResumeTimer()
{
  if(!timer_on)
  {
    ticks = clock() - ticks;
    wtime = MPI_Wtime() - wtime;
    timer_on = 1;
  }
}


void StopTimer()
{
  if(timer_on)
  {
    ticks = clock() - ticks;
    wtime = MPI_Wtime() - wtime;
    timer_on = 0;
  }
}


void PrintTimer()
{
  if(timer_on)
  {
    StopTimer();
    printf("(%i) Elapsed time: %14.6f s (%5.1f%% CPU)\n",
      rank, wtime, 100.0 * ticks * (1.0 / CLOCKS_PER_SEC) / wtime);
    ResumeTimer();
  }
  else
  {
    printf("(%i) Elapsed time: %14.6f s (%5.1f%% CPU)\n",
      rank, wtime, 100.0 * ticks * (1.0 / CLOCKS_PER_SEC) / wtime);
  }
}


/* Return a random double in the range 0..1.                                */
double RandomDouble()
{
  return(random() * 1.0 / RANMAX);
}


/* Return a random integer in the range 0..range.                           */
double RandomInteger(int range)
{
  return(random() % range);
}


/* Lennard-Jones potential;  input: r^2, output: v_{LJ}(r)                  */
double vLJ(double r2)
{
  double r2inv, r6inv;

  r2inv = 1.0 / r2;
  r6inv = r2inv * r2inv * r2inv;
  return(r6inv*r6inv-r6inv);
}


/* Lennard-Jones force divided by r, input: r^2, output f_{LJ}(r)/r          */

double fLJr(double r2)
{
  double r2inv, r6inv;

  r2inv = 1.0 / r2;
  r6inv = r2inv * r2inv * r2inv;
  return((-12.0 *  r6inv + 6.0) * r6inv * r2inv);
}


/* Unassign a particle from its cell                                        */
void CellUnassign(particleType *p)
{
  cellType *c;
  int idx, n;

  c   = &cell[p->cx][p->cy]; /* Current cell this particle is assigned to   */
  idx = p->localindex;       /* Local index of particle in the cell's list  */
  n   = c->nParticles;       /* Number of particles in the cell             */

  /* Move the last particle to the vacant list-entry:                       */
  c->particles[idx] = c->particles[n-1];
  c->particles[idx]->localindex = idx;

  /* We just removed a particle from the array:                             */
  c->nParticles--;
}


/* Remove a particle from this process                                      */
void RemoveParticle(particleType *p)
{
  CellUnassign(p);
  if(p != &particle[nPartProc-1])
  {
    /* Copy last particle into address of the removed particle              */
    memcpy(p, &particle[nPartProc-1], sizeof(particleType));
    /* Update the list of the cell the last particle is assigned to         */
    cell[p->cx][p->cy].particles[p->localindex] = p;
  }
  nPartProc--;
}


/* Move particle on this process to a specific cell                         */
void MoveParticle(particleType *p, int cx, int cy)
{
  int n;
  cellType *c;

  CellUnassign(p);
  c               = &cell[cx][cy];      
  n               = c->nParticles;
  c->particles[n] = p;
  p->localindex   = n;
  p->cx           = cx;
  p->cy           = cy;
  c->nParticles++;
}


/* Add particles to this process and assign them to a specific cell         */
void AddParticle(int cx, int cy, int nAdd)
{
  particleType *p;
  cellType *c;
  int i, n;

  /* printf("[%i] Add %i particles in (%i,%i).\n", rank, nAdd, cx, cy);     */
  if(nPartProc + nAdd > maxNpart)
    Debug("ERROR in AddParticle(): limit (maxNpart) reached!", 1);

  for (i=0; i<nAdd; i++)
  {  
    p = &particle[nPartProc];
    c = &cell[cx][cy];  
    n = c->nParticles;
    c->particles[n] = p;  
    p->localindex = n;
    p->cx = cx;
    p->cy = cy;
    c->nParticles++;
    nPartProc++;
  }
}


/* Do initial work to setup process-grid                                    */
void SetupProcGrid(int argc, char *argv[])
{
  int dims[2], periods[2], reorder, coord[2];

  if(argc!=3)
    Debug("ERROR: Need 2 parameters!", 1);
  
  nProcX = atoi(argv[1]);
  nProcY = atoi(argv[2]);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);

  if(nProc != (nProcX*nProcY))
    Debug("ERROR: Incorrect number of processes!", 1);

  dims[0] = nProcX;
  dims[1] = nProcY;
  periods[0] = 0; /* Non-periodic grid                                      */
  periods[1] = 0;
  reorder = 1;    /* Yes, we want MPI to reorder the ranks                  */

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &grid_comm);
  MPI_Comm_rank(grid_comm, &rank);
  MPI_Cart_coords(grid_comm, rank, 2, coord);

  px = coord[0];
  py = coord[1];
}


/* Read input from a file                                                   */
void ReadInput(const char *fname)
{
  FILE *fp;
  char c;

  Debug("ReadInput", 0);
  if(rank == 0)
  {
    fp = fopen(fname, "r");
    if(!fp) Debug("Cannot open input file", 1);
    if(fscanf(fp, "sizeX: %lf\n", &sizeX) != 1) Debug("sizeX?",  1);
    if(fscanf(fp, "sizeY: %lf\n", &sizeY) != 1) Debug("sizeY?",  1);
    if(fscanf(fp, "nPartX:%i\n", &nPartX) != 1) Debug("nPartX?", 1);
    if(fscanf(fp, "nPartY:%i\n", &nPartY) != 1) Debug("nPartY?", 1);
    if(fscanf(fp, "nCellX:%i\n", &nCellX) != 1) Debug("nCellX?", 1);
    if(fscanf(fp, "nCellY:%i\n", &nCellY) != 1) Debug("nCellY?", 1);
    if(fscanf(fp, "steps: %li\n", &steps) != 1) Debug("steps?",  1);
    if(fscanf(fp, "dt:    %lf\n", &dt)    != 1) Debug("dt?",     1);
    if(fscanf(fp, "logs:  %c\n",  &c)     != 1) Debug("logs?",   1);
    logs = (c == 'y');
    fclose(fp);
  }


  MPI_Bcast(&sizeX, 1, MPI_DOUBLE, 0, grid_comm);  
  MPI_Bcast(&sizeY, 1, MPI_DOUBLE, 0, grid_comm);
  MPI_Bcast(&nPartX,1, MPI_INT,    0, grid_comm);  
  MPI_Bcast(&nPartY,1, MPI_INT,    0, grid_comm);
  MPI_Bcast(&nCellX,1, MPI_INT,    0, grid_comm);
  MPI_Bcast(&nCellY,1, MPI_INT,    0, grid_comm);
  MPI_Bcast(&steps, 1, MPI_INT,    0, grid_comm);
  MPI_Bcast(&dt,    1, MPI_DOUBLE, 0, grid_comm);
  MPI_Bcast(&logs,  1, MPI_INT,    0, grid_comm);
  nPart = nPartX * nPartY;
  dt2 = dt * dt;
}


/* Init all cells                                                           */
void InitCells()
{
  int cx, cy;
  int maxCell;          /* Maximum number of particles in a cell            */

  Debug("InitCells", 0);

  sizeCellX = sizeX / nCellX;     
  sizeCellY = sizeY / nCellY;

  startcx = (px*nCellX) / nProcX;  /* index of my first cell                */
  startcy = (py*nCellY) / nProcY;

  nextcx = ((px+1)*nCellX) / nProcX; 
  nextcy = ((py+1)*nCellY) / nProcY;

  procNcx = nextcx - startcx;
  procNcy = nextcy - startcy;

  beginx = startcx * sizeCellX;
  beginy = startcy * sizeCellY;

  endx = beginx + sizeCellX * procNcx;
  endy = beginy + sizeCellY * procNcy;

  /* TODO: Chose the maximum number of particles in a cell                  */
  maxCell = nPart / (nCellX * nCellY) + 10;
  
  /* Allocate memory for internal cells plus ghostcells                     */
  cell = (cellType**) malloc((procNcx+2)*sizeof(cellType*));
  for(cx=0; cx<procNcx+2; cx++)
  {
    cell[cx] = (cellType*) malloc((procNcy+2)*sizeof(cellType));
    for(cy=0; cy<procNcy+2; cy++)
    {
      cell[cx][cy].nParticles = 0;
      cell[cx][cy].particles =
        (particleType**) malloc(maxCell * sizeof(particleType*));
    }
  }
}


/* Set initial positions and velocities for all particles.                  */
void InitParticles()
{
  int i, j, idx;
  int cx, cy;
  int nx, ny, n;
  double dx, dy, x, y;
  double vmax = 2.; /* Maximum initial velocity component                   */

  particleType *p;
  cellType *c;

  Debug("InitParticles", 0);

/* TODO: Chose maximum number of particles on a process                     */
  maxNpart = nPart;
  particle = (particleType*) malloc(maxNpart*sizeof(particleType));

  dx = sizeX / nPartX;
  dy = sizeY / nPartY;

  for(i=0; i<nPartX; i++)
  {
    x = (i+0.5) * dx;
    if (x >= beginx && x < endx)
    { 
      for(j=0; j<nPartY; j++)
      {
        y = (j+0.5) * dy;
        if (y >= beginy && y < endy)
        {
          p = &particle[nPartProc];
          cx = floor(x / sizeCellX) - startcx + 1;
          cy = floor(y / sizeCellY) - startcy + 1;
          AddParticle(cx, cy, 1);
          p->x = x;
          p->y = y;
          p->oldx = x + dt * vmax * (2. * RandomDouble() - 1.);
          p->oldy = y + dt * vmax * (2. * RandomDouble() - 1.);
        }
      }
    }
  }
}


/* Setup communication variables for neighbour cells.                       */
void SetupNeighbours()
{
  int dx, dy;
  int idx;
  int coords[2];
  int maxList; /* Maximum number of particles in a send/recv-list           */
  neighbourType *nb;

  Debug("SetupNeighbours", 0);

  /* idx runs from 0..7 to indicate        5     6     7                    */
  /* all neighbours in this way:           3  -  +  -  4                    */
  /*                                       0     1     2                    */
  idx = 0;

  for(dy=-1; dy<=1; dy++)
  {
    for(dx=-1; dx<=1; dx++)
    {
      if(dx!=0 || dy!=0)
      {
        nb = &neighbour[idx];
        coords[0] = px + dx;
        coords[1] = py + dy;

        /* Get the rank of this neighbour according to its coordinates      */
        if(coords[0]>=0 && coords[0]<nProcX && 
           coords[1]>=0 && coords[1]<nProcY)
        { 
          MPI_Cart_rank(grid_comm, coords, &nb->rank);
          nb->dx = dx;
          nb->dy = dy;

          /* Set the number of cells to communicate for send/recv           */
          if(dx==0) 
            nb->nCells = procNcx;
          else if(dy==0) 
            nb->nCells = procNcy;
          else 
            nb->nCells = 1;

          /* Set the starting cell-coordinates for send/receive operations  */
          nb->sendStartCx = (dx==1)?procNcx:1;
          nb->sendStartCy = (dy==1)?procNcy:1;
          nb->recvStartCx = nb->sendStartCx + dx;
          nb->recvStartCy = nb->sendStartCy + dy;
          
          nb->sendNumList = (int*) malloc(nb->nCells*sizeof(int));
          nb->recvNumList = (int*) malloc(nb->nCells*sizeof(int));

          /* TODO: Chose the maximum number of particles to send to or      */
          /*       receive from a neighbour                                 */
          maxList = nPart / min(nCellX, nCellY) + 100;

          /* Pointers to the particles that are communicated                */
          nb->sendParticles = (particleType**) 
            malloc(maxList*sizeof(particleType*));
          nb->recvParticles = (particleType**)
            malloc(maxList*sizeof(particleType*));

          /* Buffer for sending and receiving of x and y coordinates        */
          nb->sendList = malloc(2*maxList*sizeof(double));
          nb->recvList = malloc(2*maxList*sizeof(double));
        }
        else
        {
          nb->rank = MPI_PROC_NULL;
        }
        idx++;
      }
    }
  }
}


/* Calculate the Lennard-Jones and kinetic energy  using the cell algorithm */
void CalcEnergy(double *eLJ, double *eKin)
{
  particleType *p1, *p2;
  double dx, dy, rr;
  int cx1, cy1, cx2, cy2, j1, j2;

  *eLJ  = 0.0;
  *eKin = 0.0;

  /* Run over internal cells                                                */
  for(cx1=1; cx1<=procNcx; cx1++)
  {
    for(cy1=1; cy1<=procNcy; cy1++)
    {
      for(j1=0; j1<cell[cx1][cy1].nParticles; j1++)
      {
        p1  = cell[cx1][cy1].particles[j1];
        dx = p1->x - p1->oldx;
        dy = p1->y - p1->oldy;
        *eKin += 0.5 * (dx*dx + dy*dy) / dt2;
 
        /* Run over neighbouring cells                                      */
        for(cx2=cx1-1; cx2<=cx1+1; cx2++)
        {
          for(cy2=cy1-1; cy2<=cy1+1; cy2++)
          {
            for(j2=0; j2<cell[cx2][cy2].nParticles; j2++)
            {
              p2  = cell[cx2][cy2].particles[j2];
              if(p1 != p2)
              {
                dx = p2->x - p1->x;
                dy = p2->y - p1->y;
                rr = dx*dx + dy*dy;  /* the squared distance !!!            */
                *eLJ += 0.5 * vLJ(rr);
              }
            }
          }
        }
      }
    }
  }
}


/* Calculate the forces on all particles using the cell algoritm.           */
void CalcForces()
{
  particleType *p1, *p2;
  double dx, dy, rr, f;
  int cx1, cy1, cx2, cy2, j1, j2;

  /* Run over internal cells                                                */
  for(cx1=1; cx1<=procNcx; cx1++)
  {
    for(cy1=1; cy1<=procNcy; cy1++)
    {
      for(j1=0; j1<cell[cx1][cy1].nParticles; j1++)
      {
        p1  = cell[cx1][cy1].particles[j1];
        p1->fx = 0.0;
        p1->fy = 0.0;

        /* Run over neighbouring cells                                      */
        for(cx2=cx1-1; cx2<=cx1+1; cx2++)
        {
          for(cy2=cy1-1; cy2<=cy1+1; cy2++)
          {
            for(j2=0; j2<cell[cx2][cy2].nParticles; j2++)
            {
              p2  = cell[cx2][cy2].particles[j2];
              if(p1 != p2)
              {
                dx  = p2->x - p1->x;
                dy  = p2->y - p1->y;
                rr  = dx*dx + dy*dy;   /* the squared distance !!!          */
                f   = fLJr(rr);        /* the force divided by r !!!        */
                p1->fx += dx * f;
                p1->fy += dy * f;
              }
            }
          }
        }
      }
    }
  }
}


/* Use the Verlet algoritm to propagate particles one timestep              */
void Verlet()
{
  particleType *p;
  double oldx, oldy;
  int i;

  for(i=0; i<nPartProc; i++)
  {
    p = &particle[i];
    oldx = p->x;
    oldy = p->y;
    p->x += p->x - p->oldx + p->fx * dt2;
    p->y += p->y - p->oldy + p->fy * dt2;
    p->oldx = oldx;
    p->oldy = oldy;
  }
}


/* Setup the logfiles                                                       */
void SetupLogs()
{
  FILE *fp;

  energyLog    = "energy.dat";
  particleDump = "pdump.dat";
  trajectLog   = "traject.dat";

  if(rank==0) /* Process 0 clears up all files */
  {
    fp = fopen(energyLog, "w");
    close(fp);
    fp = fopen(particleDump, "w");
    close(fp);
    fp = fopen(trajectLog, "w");
    close(fp);
  }
}


/* Update the logfiles                                                      */
void UpdateLogs(int iter)
{
  FILE *fp;
  double time;
  double eLJ, eKin;
  double totLJ, totKin;
  int cx, cy, i, dummy;

  time = iter * dt;

  /* Logging of the energy.                                                 */
  CalcEnergy(&eLJ, &eKin);

  MPI_Reduce(&eLJ, &totLJ, 1, MPI_DOUBLE, MPI_SUM, 0, grid_comm);
  MPI_Reduce(&eKin, &totKin, 1, MPI_DOUBLE, MPI_SUM, 0, grid_comm);
  if(rank==0)
  {
    fp = fopen(energyLog, "a"); /* Append to file                           */
    fprintf(fp, "%10.3e %12.6e %12.6e %12.6e\n",
            time, totLJ, totKin, totLJ + totKin);
    fclose(fp);
  }

  /* Write particle positions at the last timestep.                         */
  /* Each process writes to the dump file waiting for its turn.             */
  if(iter == steps)
  {
    if(rank>0) /* If we're not the first process                            */
    {
      /* Wait until previous process has finished writing                   */
      MPI_Recv(&dummy, 1, MPI_INT, rank-1, 19, grid_comm, &status);
    }
    fp = fopen(particleDump, "a");
    for(cx=1; cx<procNcx+1; cx++) /* Run over all (internal) cells          */
    {
      for(cy=1; cy<procNcy+1; cy++)
      {
        for(i=0; i<cell[cx][cy].nParticles; i++)
        {
          fprintf(fp,"%10.4e %10.4e\n",
            cell[cx][cy].particles[i]->x, cell[cx][cy].particles[i]->y);
        }
      }
    }
    fclose(fp);
    if(rank<nProc-1) /* If we're not the last process                       */
    {
      /* Next process may now append its part to the logfile                */
      MPI_Send(&dummy, 1, MPI_INT, rank+1, 19, grid_comm);
    }
  }

  /* Log the trajectory of the first particle in the first process          */
  i = 1;
  if(rank == 0)
  {
    fp = fopen(trajectLog, "a");
    fprintf(fp, "%10.4e %10.4e\n", particle[i].x, particle[i].y);
    fclose(fp);
  }
}


/* Show the cell occupations on all processes                               */
void ShowCells()
{
  int cx, cy, dummy;

  MPI_Barrier(grid_comm); /* Hopefully, flush all output buffers            */

  if(rank>0) MPI_Recv(&dummy, 1, MPI_INT, rank-1, 11, grid_comm, &status);

  printf("Process %i reporting in...\n", rank);

  for(cy=procNcy+1; cy>=0; cy--)
  {
    printf("(process %i)  ", rank);
    for(cx=0; cx<=procNcx+1; cx++)
    {
      if(cx==0 || cx==procNcx+1 || cy==0 || cy==procNcy+1)
        printf("[%2i]", cell[cx][cy].nParticles);
      else
        printf(" %2i ", cell[cx][cy].nParticles);
    }
    printf("\n");
  }
  printf("\n");

  if(rank<nProc-1) MPI_Send(&dummy, 1, MPI_INT, rank+1, 11, grid_comm);

  MPI_Barrier(grid_comm);
}


/* Reassign particles in a process to the appropriate cell                  */
void ReassignParticles()
{
  particleType *p;
  int i;
  int newgcx, newgcy;         /* New global cell coordinates                */
  int cx, cy, newcx, newcy;   /* Local cell coordinates                     */

  for (i=nPartProc-1; i>=0; i--)
  {
    p = &particle[i];
    /* The current cell indices for this particle                           */
    cx = p->cx;
    cy = p->cy;

    /* Obtain new global cell indices                                       */
    newgcx = floor(p->x / sizeCellX);
    newgcy = floor(p->y / sizeCellY);
    if(newgcx<0)       newgcx = 0;
    if(newgcx>=nCellX) newgcx = nCellX-1;
    if(newgcy<0)       newgcy = 0;
    if(newgcy>=nCellY) newgcy = nCellY-1;

    /* ... and corresponding local indices                                  */
    newcx = newgcx - startcx + 1;
    newcy = newgcy - startcy + 1;

    if (newcx != cx || newcy != cy)
    {
      if (abs(newcx-cx) > 1 || abs(newcy-cy) > 1)
      {
        Debug("ReassignParticles: Reduce speed or reassign more often", 1);
      }
      else if (newgcx<startcx-1 || newgcx>=nextcx+1 ||
               newgcy<startcy-1 || newgcy>=nextcy+1)
      { /* New indices are outside the range for this process               */
        RemoveParticle(p);
      }
      else  
      { /*  Move particle within process to a neighbouring cell             */
        MoveParticle(p, newcx, newcy);  
      }
    }
  }
}


/* Setup or reset the communication with neighbours                         */
void SetupBorderExchange()
{
  int cx, cy;
  int idx, i, k, cnt, n, nAdd;
  neighbourType *nbS, *nbR;
  particleType *p;

  Debug("SetupBorderExchange", 0);

  for(idx=0; idx<8; idx++)
  {
    nbS = &neighbour[idx];  /* Neighbour to send to      */
    nbR = &neighbour[idx];  /* Neighbour to receive from */

    if (nbS->rank != MPI_PROC_NULL)
    {
      /* Prepare sendNumList[] and sendParticles[]                         */
      cx = nbS->sendStartCx;
      cy = nbS->sendStartCy;
      cnt = 0;

      for(k=0; k<nbS->nCells; k++)
      {
        n = cell[cx][cy].nParticles;
        nbS->sendNumList[k] = n;
        for(i=0; i<n; i++)
        {
          nbS->sendParticles[cnt+i] = cell[cx][cy].particles[i];
        }
        if (nbS->dx != 0) cy++;
        if (nbS->dy != 0) cx++;
        cnt += n;
      }
      nbS->nPartSend = cnt;
      printf("[%i] sends %i particles to [%i] \n", rank, cnt, nbS->rank);
    }

    MPI_Sendrecv(nbS->sendNumList, nbS->nCells, MPI_INT, nbS->rank, 14, 
                 nbR->recvNumList, nbR->nCells, MPI_INT, nbR->rank, 14, 
                 grid_comm, &status);

    if (nbR->rank != MPI_PROC_NULL)
    {
      /* Prepare recvNumList[] and recvparticles[]                          */
      cx = nbR->recvStartCx;
      cy = nbR->recvStartCy;
      cnt = 0;

      for(k=0; k<nbR->nCells; k++)
      {
        n = nbR->recvNumList[k];

        /* Particles entered ghostcells, so add as many as necessary        */
        nAdd = n-cell[cx][cy].nParticles;
        if (nAdd>0) AddParticle(cx, cy, nAdd);
        if (nAdd<0) Debug("SetupBorderExchange: nAdd<0", 1); 

        for(i=0; i<n; i++)
        {
          nbR->recvParticles[cnt+i] = cell[cx][cy].particles[i];
        }
        if (nbR->dx != 0) cy++;
        if (nbR->dy != 0) cx++;
        cnt += n;
      }
      nbR->nPartReceived = cnt;
    }
  }
}


/* Exchange the (x,y) coordinates of bordercell/ghostcell particles         */
void ExchangeBorders()
{
  int idx, j, nSend, nRecv;
  neighbourType *nbS, *nbR;
  double (*sendList)[2], (*recvList)[2];
  particleType **sendParticles, **recvParticles;

  Debug("ExchangeBorders", 0);

  for(idx=0; idx<8; idx++)
  {
    nbS            = &neighbour[idx]; /* neighbour to send to       */
    nbR            = &neighbour[idx]; /* neighbour to receive from  */

    sendList      = nbS->sendList;
    sendParticles = nbS->sendParticles;
    nSend         = nbS->nPartSend;

    recvList      = nbR->recvList;
    recvParticles = nbR->recvParticles;
    nRecv         = nbR->nPartReceived;

    /* Prepare the sendLists for all neighbours                             */
    if (nbS->rank != MPI_PROC_NULL)
      for(j=0; j<nSend; j++)
      {
        sendList[j][0] = sendParticles[j]->x;
        sendList[j][1] = sendParticles[j]->y;
      }

    MPI_Sendrecv(sendList, 2*nSend, MPI_DOUBLE, nbS->rank, 13,
                 recvList, 2*nRecv, MPI_DOUBLE, nbR->rank, 13, 
                 grid_comm, &status);
  
    /* Update particles in ghostcells upon received recvLists               */
    if(nbR->rank != MPI_PROC_NULL)
      for(j=0; j<nRecv; j++)
      {
        recvParticles[j]->x = recvList[j][0];
        recvParticles[j]->y = recvList[j][1];
      }
  }
}


/* Run                                                                      */
void Run(int numberOfSteps)
{
  int step;

  Debug("Run", 0);

  SetupBorderExchange();
  ExchangeBorders();

  for(step=1; step<=numberOfSteps; step++)
  {
    CalcForces();
    Verlet(); 
    ExchangeBorders();

    if(logs==1) {UpdateLogs(step);}

    if(step%10==0)  /* After each 10 steps we update the cell-assignments   */
    {
      ReassignParticles();
      SetupBorderExchange();
      ExchangeBorders();
      printf("[%i] After step %i I have %i part.\n", rank, step, nPartProc);
    }
  }
}


/* Main program (with command-line parameters)                              */
int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  StartTimer();
  SetupProcGrid(argc, argv);
  ReadInput("input.dat");
  InitCells();
  InitParticles();
  SetupNeighbours();
  SetupLogs();
  Run(steps);
  ShowCells();

  PrintTimer();

  MPI_Finalize();
  return 0;
}

