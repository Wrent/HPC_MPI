#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG 0

#define TYPE_GHOST 1
#define TYPE_SOURCE 2

enum
{
  X_DIR, Y_DIR
};

int gridsize[2];
int P_grid[2];			/* processgrid dimensions */
int N_sources;			/* number of sources */
int *source;			/* global vertex id of source */
double *source_val;		/* value of sources */
int do_adapt;			/* perfrom grid adaptation */

void Debug(char *mesg, int terminate);
void Setup_Grid(int argc, char **argv);
void Write_Grid();
void Write_Datafiles();
void Write_GraphMap();

#include "grid.c"

void Debug(char *mesg, int terminate)
{
  if (DEBUG || terminate)
    printf("%s\n", mesg);
  if (terminate)
    exit(1);
}

void Setup_Grid(int argc, char **argv)
{
  int i, x, y, N = 0;
  int wrong_param = 0;
  double source_x, source_y;
  FILE *f;

  Debug("Setup_Grid", 0);

  gridsize[X_DIR] = gridsize[Y_DIR] = 0;

  do_adapt = 0;

  if ( (argc != 5) && (argc != 6) )
    wrong_param = 1;
  else
  {
    P_grid[X_DIR] = atoi(argv[1]);
    P_grid[Y_DIR] = atoi(argv[2]);
    N = P_grid[X_DIR] * P_grid[Y_DIR];
    gridsize[X_DIR] = atoi(argv[3]);
    gridsize[Y_DIR] = atoi(argv[4]);
    if ((N == 0) || (gridsize[X_DIR] * gridsize[Y_DIR] == 0))
      wrong_param = 1;
    if (argc==6)
    {
      if (strcmp(argv[5],"adapt"))
        wrong_param = 1;
      else
        do_adapt = 1;
    }
  }
  if (wrong_param)
    Debug("Wrong number of parameters.\nUse : GridDist <Px> <Py> <dim_x> <dim_y> [adapt]", 1);

/****/
  nx = gridsize[X_DIR];
  ny = gridsize[Y_DIR];
  nm = 50;
/****/

  if ((f = fopen("sources.dat", "r")) == NULL)
    Debug("Can't open sources.dat", 1);

  fscanf(f, "%i\n", &N_sources);
  if (N_sources > 0)
  {
    if ((source = malloc(N_sources * sizeof(int))) == NULL)
        Debug("Error: malloc 'source'", 1);
    if ((source_val = malloc(N_sources * sizeof(double))) == NULL)
        Debug("Error: malloc 'source_val'", 1);

/****/
    nsource = N_sources;
    if ((source2=malloc(nsource*sizeof(sourcepoint)))==NULL)
      Debug("Out of memory", 1);
/****/

    for (i = 0; i < N_sources; i++)
    {
      fscanf(f, "source: %lf %lf %lf\n", &source_x, &source_y, &source_val[i]);
      x = floor(0.5 + source_x * (gridsize[X_DIR] - 1));
      y = floor(0.5 + source_y * (gridsize[Y_DIR] - 1));
      source[i] = y * gridsize[X_DIR] + x;
      source2[i].xpos = source_x; /****/
      source2[i].ypos = source_y; /****/
      source2[i].value = source_val[i]; /****/
/*
      printf("(%f,%f) --> (%f,%f) --> (%i,%i) --> %i\n", source_x, source_y,
	     (double) x / (gridsize[X_DIR] - 1),
	     (double) y / (gridsize[Y_DIR] - 1), x, y, source[i]);
*/
    }
  }
/*
  for (i = 0; i < N_sources; i++)
    printf("source %i : %f\n", source[i], source_val[i]);
*/
}

void Write_Grid()
{
  int x, y;
  char filename[25];
  FILE *f;

  Debug("Write_Grid", 0);
  printf("Writing Grid\n");
  sprintf(filename, "grid%i.dat", P_grid[X_DIR] * P_grid[Y_DIR]);

  if ((f=fopen(filename,"w")) == NULL)
    Debug("Write_Grid: Could not open data outputfile", 1);

  for (y=0;y<gridsize[Y_DIR];y++)
    for (x=0;x<gridsize[X_DIR];x++)
      fprintf(f,"%f %f\n", grid[x+y*gridsize[X_DIR]+1].xpos,
                           grid[x+y*gridsize[Y_DIR]+1].ypos);
}


void Write_Datafiles()
{
  int i, x, y, t, v;
  int px, py;
  int x_off, y_off, x_dim, y_dim;
  int N_vert, N_elm;
  int top, left, right, bottom;
  int start, end;
  double s_val = 0;
  char filename[25];
  FILE *f;

  Debug("Write_Datafiles", 0);

  printf("Writing file");
  for (py = 0; py < P_grid[Y_DIR]; py++)
    for (px = 0; px < P_grid[X_DIR]; px++)
    {
      printf(" %i", py * P_grid[X_DIR] + px);
      fflush(stdout);
      sprintf(filename, "input%i-%i.dat", P_grid[X_DIR] * P_grid[Y_DIR],
	      py * P_grid[X_DIR] + px);
      if ((f = fopen(filename, "w")) == NULL)
	Debug("Write_Datafiles: Could not open data outputfile", 1);
      x_off = gridsize[X_DIR] * px / P_grid[X_DIR];
      y_off = gridsize[Y_DIR] * py / P_grid[Y_DIR];
      x_dim = gridsize[X_DIR] * (px + 1) / P_grid[X_DIR] - x_off;
      y_dim = gridsize[Y_DIR] * (py + 1) / P_grid[Y_DIR] - y_off;
      top = bottom = right = left = start = end = 0;
      if (py != 0)
	top = 1;
      if (py != P_grid[Y_DIR] - 1)
	bottom = 1;
      if (px != 0)
	left = 1;
      if (px != P_grid[X_DIR] - 1)
	right = 1;
      if (top && left)
	start = -1;
      if (bottom && right)
	end = -1;
      x_dim += left + right;
      y_dim += top + bottom;
      x_off -= left;
      y_off -= top;

      /* print vertices */
      N_vert = x_dim * y_dim + start;
      if (bottom && right)
	N_vert--;
      fprintf(f, "N_vert: %i\n", N_vert);
      fprintf(f, "id x y type\n");
      for (y = 0; y < y_dim; y++)
	for (x = ((y == 0) ? -start : 0); x < x_dim +
	     ((y == y_dim - 1) ? end : 0); x++)
	{
	  t = 0;
	  if (((x == 0) && left) || ((y == 0) && top) ||
	      ((x == x_dim - 1) && right) || ((y == y_dim - 1) && bottom))
	    t += TYPE_GHOST;
	  v = (y + y_off) * gridsize[X_DIR] + (x + x_off);

	  /* check if current vertex is a source */
	  if ((x + x_off == 0) || (x + x_off == gridsize[X_DIR] - 1) ||
	      (y + y_off == 0) || (y + y_off == gridsize[Y_DIR] - 1))
	  {
	    t |= TYPE_SOURCE;
	    s_val = 0;
	  }
	  for (i = 0; (i < N_sources) && (source[i] != v); i++) ;
	  if (i < N_sources)
	  {
	    t |= TYPE_SOURCE;
	    s_val = source_val[i];
	  }

          if (do_adapt)
	    fprintf(f, "%i %20.15e %20.15e %i ", y * x_dim + x + start,
		    grid[x+x_off+(y+y_off)*gridsize[X_DIR]+1].xpos,
		    grid[x+x_off+(y+y_off)*gridsize[X_DIR]+1].ypos, t);
          else
	    fprintf(f, "%i %20.15e %20.15e %i ", y * x_dim + x + start,
		    ((float) x + x_off) / (gridsize[X_DIR] - 1),
		    ((float) y + y_off) / (gridsize[Y_DIR] - 1), t);

	  if (t & TYPE_SOURCE)
	    fprintf(f, "%20.15e\n", s_val);
	  else
	    fprintf(f, "%20.15e\n", 0.0);
	}

      /* print elements */
      N_elm = 2 * (x_dim - 1) * (y_dim - 1) + start + end;
      fprintf(f, "N_elm: %i\n", N_elm);
      fprintf(f, "id v1 v2 v3\n");
      for (y = 0; y < y_dim - 1; y++)
	for (x = 0; x < x_dim - 1; x++)
	{
	  if ((y != 0) || (x != 0) || (start == 0))
	    fprintf(f, "%i %i %i %i\n", (y * (x_dim - 1) + x) * 2 + start,
		    y * x_dim + x + start, y * x_dim + x + 1 + start,
		    (y + 1) * x_dim + x + start);
	  if ((y != y_dim - 2) || (x != x_dim - 2) || (end == 0))
	    fprintf(f, "%i %i %i %i\n", (y * (x_dim - 1) + x) * 2 + 1 + start,
		    y * x_dim + x + 1 + start, (y + 1) * x_dim + x + start,
		    (y + 1) * x_dim + x + 1 + start);
	}

      /* print neighbour connectivity */
      fprintf(f, "Neighbours: %i\n", top + bottom + left + right +
	      ((top && right) ? 1 : 0) + ((bottom && left) ? 1 : 0));
      if (top)
      {
	fprintf(f, "from %i :", (py - 1) * P_grid[X_DIR] + px);
	for (x = left; x < x_dim - right; x++)
	  fprintf(f, " %i", x + start);
	fprintf(f, "\n");
	fprintf(f, "to %i :", (py - 1) * P_grid[X_DIR] + px);
	for (x = left; x < x_dim - right; x++)
	  fprintf(f, " %i", x_dim + x + start);
	fprintf(f, "\n");
      }
      if (bottom)
      {
	fprintf(f, "from %i :", (py + 1) * P_grid[X_DIR] + px);
	for (x = left; x < x_dim - right; x++)
	  fprintf(f, " %i", (y_dim - 1) * x_dim + x + start);
	fprintf(f, "\n");
	fprintf(f, "to %i :", (py + 1) * P_grid[X_DIR] + px);
	for (x = left; x < x_dim - right; x++)
	  fprintf(f, " %i", (y_dim - 2) * x_dim + x + start);
	fprintf(f, "\n");
      }
      if (left)
      {
	fprintf(f, "from %i :", py * P_grid[X_DIR] + px - 1);
	for (y = top; y < y_dim - bottom; y++)
	  fprintf(f, " %i", y * x_dim + start);
	fprintf(f, "\n");
	fprintf(f, "to %i :", py * P_grid[X_DIR] + px - 1);
	for (y = top; y < y_dim - bottom; y++)
	  fprintf(f, " %i", y * x_dim + 1 + start);
	fprintf(f, "\n");
      }
      if (right)
      {
	fprintf(f, "from %i :", py * P_grid[X_DIR] + px + 1);
	for (y = top; y < y_dim - bottom; y++)
	  fprintf(f, " %i", (y + 1) * x_dim - 1 + start);
	fprintf(f, "\n");
	fprintf(f, "to %i :", py * P_grid[X_DIR] + px + 1);
	for (y = top; y < y_dim - bottom; y++)
	  fprintf(f, " %i", (y + 1) * x_dim - 2 + start);
	fprintf(f, "\n");
      }
      if (top && right)
      {
	fprintf(f, "from %i : %i\n", (py - 1) * P_grid[X_DIR] + px + 1,
		x_dim - 1 + start);
	fprintf(f, "to %i : %i\n", (py - 1) * P_grid[X_DIR] + px + 1,
		2 * x_dim - 2 + start);
      }
      if (bottom && left)
      {
	fprintf(f, "from %i : %i\n", (py + 1) * P_grid[X_DIR] + px - 1,
		(y_dim - 1) * x_dim);
	fprintf(f, "to %i : %i\n", (py + 1) * P_grid[X_DIR] + px - 1,
		(y_dim - 2) * x_dim + 1);
      }
      fclose(f);
    }
}

void Write_GraphMap()
{
  int px, py;
  int n = 0;
  FILE *f;
  char filename[25];

  Debug("Write_GraphMap", 0);

  sprintf(filename, "mapping%i.dat", P_grid[X_DIR] * P_grid[Y_DIR]);
  if ((f = fopen(filename, "w")) == NULL)
    Debug("Write_GraphMap: Could not open mapping outputfile", 1);

  printf(" map\n");

  fprintf(f, "N_proc : %i\n", P_grid[X_DIR] * P_grid[Y_DIR]);
  fprintf(f, "number of neighbours :\n");
  for (py = 0; py < P_grid[Y_DIR]; py++)
    for (px = 0; px < P_grid[X_DIR]; px++)
    {
      if (py > 0)
	n++;
      if (py < P_grid[Y_DIR] - 1)
	n++;
      if (px > 0)
	n++;
      if (px < P_grid[X_DIR] - 1)
	n++;
      if ((py > 0) && (px < P_grid[X_DIR] - 1))
	n++;
      if ((px > 0) && (py < P_grid[Y_DIR] - 1))
	n++;
      fprintf(f, "%i\n", n);
    }

  fprintf(f, "neighbours :\n");
  for (py = 0; py < P_grid[Y_DIR]; py++)
    for (px = 0; px < P_grid[X_DIR]; px++)
    {
      if (py > 0)
	fprintf(f, "%i\n", (py - 1) * P_grid[X_DIR] + px);
      if (py < P_grid[Y_DIR] - 1)
	fprintf(f, "%i\n", (py + 1) * P_grid[X_DIR] + px);
      if (px > 0)
	fprintf(f, "%i\n", py * P_grid[X_DIR] + px - 1);
      if (px < P_grid[X_DIR] - 1)
	fprintf(f, "%i\n", py * P_grid[X_DIR] + px + 1);
      if ((py > 0) && (px < P_grid[X_DIR] - 1))
	fprintf(f, "%i\n", (py - 1) * P_grid[X_DIR] + px + 1);
      if ((px > 0) && (py < P_grid[Y_DIR] - 1))
	fprintf(f, "%i\n", (py + 1) * P_grid[X_DIR] + px - 1);
    }
  fclose(f);
}

int main(int argc, char **argv)
{
  Setup_Grid(argc, argv);
  if (do_adapt)
  {
    adaptgrid();
    Write_Grid();
  }
  Write_Datafiles();
  Write_GraphMap();

  free(source);
  free(source_val);

  return 0;
}

