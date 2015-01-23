/***
 * This program divides a 2D domain [0..1] by [0..1] up
 * into nx by ny triangular connected gridpoints.
 ***/

int nx;
int ny;
int nm;

typedef struct
{
  int id;
  int fixed;
  float xpos, ypos;
  float fx, fy;
  int neighbour[6];
} gridpoint;

typedef struct
{
  float xpos, ypos;
  float value;
} sourcepoint;

gridpoint *grid;
sourcepoint *source2;
int ngrid;
int nsource;

float sqr(float a)
{
  return(a*a);
}

void gridgen()
{
  float dx, dy;
  int i, j, num, id;
  gridpoint *g;

  ngrid = nx*ny;
  num = nx*ny + 1;
  grid = malloc(num * sizeof(gridpoint));

  dx = 1.0 / (nx-1);
  dy = 1.0 / (ny-1);
  for(i=0; i<nx; i++)
  {
    for(j=0; j<ny; j++)
    {
      id = i + j*nx + 1;
      g = &grid[id];

      g->id = id;
      g->xpos = i*dx;
      g->ypos = j*dy;

      g->neighbour[0] = (i>0)              ? id-1 : 0;
      g->neighbour[1] = (i>0 && (j<ny-1))  ? id-1+nx : 0;
      g->neighbour[2] = (j<ny-1)           ? id+nx : 0;
      g->neighbour[3] = (i<nx-1)           ? id+1 : 0;
      g->neighbour[4] = ((i<nx-1) && j>0)  ? id+1-nx : 0;
      g->neighbour[5] = (j>0)              ? id-nx : 0;

      if( (i>0) && (i<nx-1) && (j>0) && (j<ny-1) )
	g->fixed = 0;
      else
	g->fixed = 1;
    }
  }
}

void sources_on_gridpoints()
{
  int i, j, minidx;
  float r, minr;

  for(i=0; i<nsource; i++)
  {
    minr = 10.0;
    minidx = 0;
    for(j=1; j<=ngrid; j++)
    {
      r = sqrt(
        sqr(source2[i].xpos-grid[j].xpos) +
        sqr(source2[i].ypos-grid[j].ypos));
      if((grid[j].fixed==0) && (r<minr))
      {
        minr = r;
        minidx = j;
      }
    }
    grid[minidx].xpos = source2[i].xpos;
    grid[minidx].ypos = source2[i].ypos;
    grid[minidx].fixed = 1;
  }
}

float calc_springc(float x, float y)
{
  int i, minidx;
  float dx, dy, r;
  float minr = 10.0;
  float rate;

  for(i=0; i<nsource; i++)
  {
    dx = source2[i].xpos - x;
    dy = source2[i].ypos - y;
    r = sqrt(sqr(dx)+sqr(dy));
    if(r<minr)
    {
      minr = r;
      minidx = i;
    }
  }
  rate = (minr<0.25) ? minr : 0.25;
  rate = 1.0 - (1.0/0.25) * rate;
  return( 0.1 * sqr(sqr(rate)) + 0.02 );
}

float gridmove()
{
  gridpoint *g, *g2;
  float dx, dy, r;
  float f, fmin;
  float centerx, centery;
  float springc;
  float maxdiff = 0;
  int i, j;

  for(i=1; i<ngrid; i++)
  {
    g = &grid[i];
    if(g->fixed == 0)
    {
      g->fx = 0;
      g->fy = 0;
      for(j=0; j<6; j++)
      {
        g2 = &grid[g->neighbour[j]];
        dx = g2->xpos - g->xpos;
        dy = g2->ypos - g->ypos;
        r = sqrt( dx*dx + dy*dy );
        centerx = g->xpos + 0.5 * dx;
        centery = g->ypos + 0.5 * dy;
        springc = calc_springc(centerx, centery);
        g->fx += springc * dx;
	g->fy += springc * dy;
      }
    }
  }

  fmin = 0.1 / (nx + ny);

  for(i=1; i<ngrid; i++)
  {
    g = &grid[i];
    if(g->fixed == 0)
    {
      f = sqrt(sqr(g->fx)+sqr(g->fy));
      if( f > fmin )
      {
        g->fx *= fmin/f;
        g->fy *= fmin/f;
        f = fmin;
      }
      g->xpos += g->fx;
      g->ypos += g->fy;
      if( f > maxdiff ) maxdiff = f;
    }
  }

  return(maxdiff);
}

void grid_deform(int count)
{
  int i;
  float maxdiff;

  for(i=0; i<count; i++)
  {
    maxdiff = gridmove();
    fprintf(stderr, "iter %3i: %10.2e\n", i, maxdiff);
    if( maxdiff < 1e-7 )
      break;
  }
}

void adaptgrid()
{
  gridgen();
  sources_on_gridpoints();
  grid_deform(nm);
}


