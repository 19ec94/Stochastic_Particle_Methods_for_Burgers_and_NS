#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <time.h>
# include <stdint.h>
# include "ziggurat.h"

void find_cell_x (double position_x, int *xl, int *xr, int *xm, double *x_domain);
int main (int argc, char *argv[])
{

  FILE *fp, *fp1;
  fp = fopen ("serial.txt", "w+");
  //fp1 = fopen ("C3.txt", "w+");
  srand (time (NULL));
  int nb_cells = 10, nb_particles = 1;
  int total_nb_particles = nb_particles * nb_cells;
  double dt = 0.01, final_time = 0.1;
  uint32_t nt = final_time / dt;
  double nu = 0.0;
  double xStart = 0, xEnd = 1;
  double dx = (xEnd - xStart) / (double) nb_cells;
  float fn[128],wn[128];
  uint32_t kn[128],seed=time(NULL);
  r4_nor_setup(kn,fn,wn);
  
  //double x[nb_cells + 1];
  double *x = (double *) malloc ((nb_cells + 1) * sizeof (double));
  //double cell_centre[nb_cells];
  double *cell_centre = (double *) malloc ((nb_cells) * sizeof (double));
  //double cell_vel[nb_cells], new_cell_vel[nb_cells], particles[nb_cells];
  double *cell_vel = (double *) malloc ((nb_cells) * sizeof (double));
  double *new_cell_vel = (double *) malloc ((nb_cells) * sizeof (double));
  double *particles = (double *) malloc ((nb_cells) * sizeof (double));
  //double U[nt+1][nb_cells];
  double **U = (double **) malloc ((nt + 1) * sizeof (double *));
  for (int i = 0; i < (nt + 1); i++)
    U[i] = (double *) malloc (nb_cells * sizeof (double));
  //double par_old[3][total_nb_particles];
  double **par_old = (double **) malloc (3 * sizeof (double *));
  for (int i = 0; i < 3; i++)
    par_old[i] = (double *) malloc (total_nb_particles * sizeof (double));
  //double par_new[3][total_nb_particles];
  double **par_new = (double **) malloc (3 * sizeof (double *));
  for (int i = 0; i < 3; i++)
    par_new[i] = (double *) malloc (total_nb_particles * sizeof (double));
  clock_t start=clock(),diff;


  for (int i = 0; i < nb_cells + 1; i++)
    {
      x[i] = i * dx;
    }

  for (int i = 0; i < nb_cells; i++)
    {
      cell_centre[i] = x[i] + (dx / 2);
    }

  for (int i = 0; i < nb_cells; i++)
    {
      cell_vel[i] = sin (2 * M_PI * cell_centre[i]);
    }



  for (int i = 0; i < nt + 1; i++)
    {
      for (int j = 0; j < nb_cells; j++)
	{
	  U[i][j] = 0;
	  //printf("U[%d][%d]=%f\t",i,j,U[i][j]);
	}
      //printf("\n");
    }


  for (int i = 0; i < nb_cells; i++)
    {
      U[0][i] = cell_vel[i];
      //fprintf(fp,"%f ",U[0][i]);
      for (int j = (i * nb_particles); j < ((i + 1) * nb_particles); j++)
	{
	  par_old[0][j] = cell_centre[i];
	  par_old[1][j] = cell_vel[i];
	  par_old[2][j] = i;
	}
    }
  //fprintf (fp,"\n");

  for (int i = 0; i < total_nb_particles; i++)
    {
      par_new[0][i] = par_old[0][i];
      par_new[1][i] = par_old[1][i];
      par_new[2][i] = par_old[2][i];
    }


  for (int time = 1; time <= nt; time++)
    {				//time loop- start
      //printf ("Time = %d \n", time);
      for (int i = 0; i < total_nb_particles; i++)
	{			//pos-update-start
	  par_new[0][i] =
	    par_old[0][i] + (dt * par_old[1][i]) +
	    (sqrt (2 * nu) * sqrt (dt) * r4_nor ( &seed, kn, fn, wn ));
	  if (par_new[0][i] > xEnd)
	    {
	      par_new[0][i] = par_new[0][i] - xEnd;
	    }
	  else if (par_new[0][i] < xStart)
	    {
	      par_new[0][i] = par_new[0][i] + xEnd;
	    }
	  else
	    {
	    }
	}			//pos-update-end

      for (int i = 0; i < total_nb_particles; i++)
	{
	  double position_x = par_new[0][i];
	  int xl = 0, xr = nb_cells + 1, xm = floor ((xl + xr) / 2);
	  find_cell_x (position_x, &xl, &xr, &xm, x);
	  par_new[2][i] = xl;
	}
      for (int i = 0; i < nb_cells; i++)
	{
	  new_cell_vel[i] = 0;
	  particles[i] = 0;
	}

      for (int i = 0; i < total_nb_particles; i++)
	{
	  int c = par_new[2][i];
	  new_cell_vel[c] = new_cell_vel[c] + par_new[1][i];
	  particles[c] = particles[c] + 1;
	}
      for (int i = 0; i < nb_cells; i++)
	{
	  if (particles[i] == 0)
	    {
	      new_cell_vel[i] = 0;
	    }
	  else
	    {
	      new_cell_vel[i] = new_cell_vel[i] / particles[i];
	    }
	}
      for (int i = 0; i < total_nb_particles; i++)
	{
	  int c = par_new[2][i];
	  par_new[1][i] = new_cell_vel[c];
	}

      for (int i = 0; i < total_nb_particles; i++)
	{
	  par_old[0][i] = par_new[0][i];
	  par_old[1][i] = par_new[1][i];
	  par_old[2][i] = par_new[2][i];
	}
      for (int i = 0; i < nb_cells; i++)
	{
	  U[time][i] = new_cell_vel[i];
	   //fprintf(fp,"%f ",U[time][i]);
	}
      //fprintf(fp,"\n");


    }

  for (int i = 0; i <= nt; i++)
    {
      for (int j = 0; j < nb_cells; j++)
	{
	  fprintf (fp, "%f ", U[i][j]);
	}
      fprintf (fp, "\n");
    }

  diff =clock() -start;
  double msec =diff * 1000 /CLOCKS_PER_SEC;
  printf("time %f seconds  %d milliseconds\n",msec/1000,(int)msec%1000);
  //time loop-end
  fclose (fp);
  free (x);
  free (cell_centre);
  free (cell_vel);
  free (new_cell_vel);
  free (particles);

  return 0;
}
void
find_cell_x (double position_x, int *xl, int *xr, int *xm, double *x_domain)
{
  if (position_x <= x_domain[(*xm)])
    {
      *xl = *xl;
      *xr = *xm;
      *xm = floor ((*xl + *xr) / 2);
      if ((*xr - *xl) == 1)
	{
	  *xl = *xl;
	  *xr = *xr;
	  return;
	}
      else
	{
	  find_cell_x (position_x, xl, xr, xm, x_domain);
	}
    }
  else
    {
      *xl = *xm;
      *xr = *xr;
      *xm = floor ((*xl + *xr) / 2);
      if ((*xr - *xl) == 1)
	{
	  *xl = *xl;
	  *xr = *xr;
	  return;
	}
      else
	{
	  find_cell_x (position_x, xl, xr, xm, x_domain);
	}
    }
  return;
}

