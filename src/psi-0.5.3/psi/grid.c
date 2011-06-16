/* grid.c -- grid for Gibbs sampler
   Copyright (C) 2006 Sang Chul Choi
  
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA.
*/

/** @start 1 */
#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "common.h"
#include "error.h"
#include "ds.h"
#include "grid.h"

double grid_put_out (double d)
{
  return d;
}

/*
#include "system.h"
#include "defs.h"
#include "structs.h"
*/
static int want_debug_gibbs = 0;
static int want_verbose = 0;

static int gibbs_size;
/* static double gridpoint_s_begin;
static double gridpoint_s_end;
static int gridpoint_s_number;
static double gridpoint_p_begin;
static double gridpoint_p_end;
static int gridpoint_p_number;
static double gridpoint_n_begin;
static double gridpoint_n_end;
static int gridpoint_n_number; */

static Gibbs_data *gibbs_list = NULL;
static double *gmel_grid_gridpoints[5] = {NULL, NULL, NULL, NULL, NULL};
static GibbsPart ******gmel_grid_info = NULL;
/* static char grid_name[5] = { 'S', 'P', 'A', 'C', 'T' }; */
static int gmel_grid_n[5] = { 0, 0, 0, 0, 0 };
static double grid_min[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
static double grid_max[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
static int number_multi_gibbs = 1;
static int current_gibbs = 0;

/* must be called before the gmel_grid_new function */
void
psi_grid_create_multi_gibbs (int n)
{
  assert (n > 0);
  int i;

  number_multi_gibbs = n;
  if (n > 1)
    {
      for (i = 0; i < n; i++)
        {
          gibbs_set (&gibbs_list, i, NULL);
        }
    }
    
  /* current_gibbs = 0; already set to be zero as default */
}

void
psi_grid_change_multi_gibbs (int n)
{
  assert (number_multi_gibbs > 1);
  if (n == current_gibbs)
    {
      psi_warning ("gibbs grid change: No change");
      return;
    }
  /* old one */
  Gibbs_data *g = gibbs_find (gibbs_list, current_gibbs);
  if (g == NULL)
    psi_fatal ("no such gibbs: %d", current_gibbs);
  g->info = gmel_grid_info;

  /* new one */
  g = gibbs_find (gibbs_list, n);
  if (g == NULL)
    psi_fatal ("no such gibbs: %d", n);
  gmel_grid_info = g->info;

  current_gibbs = n;
}

static void 
free_all_gmel_grid_gridpoints (void);

/*!\brief Locate a gridpoint for a given one parameter

   \param index the return value of gridpoint index
   \param param_val parameter value near indexed gridpoint's parameter value
   \param which_grid one of five GRIDs
   \return 0 for EXIT_SUCCESS
 */
static int gmel_grid_locate_a_parameter (int *index, double param_val, 
                                         int which_grid);
static int gmel_grid_check_sparse_a_parameter (double top, double bot, 
                                               int which_grid);
static int gmel_grid_check_sparse_locate (double val, int which_grid);
static int gmel_grid_print_gridpoints (int which_grid);

int
gmel_grid_new_grid (double *grid, double min, double max, int n_g)
{
  /*
           0       1       2       3       4       5       6   --> n : 7
  GRID b---|-------|-------|-------|-------|-------|-------|---e
       m                                                       M
  w        |<----->| or        |<----->|
  
  Below is WRONG!
               0       1       2       3       4       5       6   --> n : 7
  GRID b---m---|-------|-------|-------|-------|-------|-------|---M---e
  w        |<----->| or        |<----->|
  */

  int i, m, n;
  double b, e, w;

  if (grid == NULL) 
    {
       psi_fatal ("Gridpoints must have been allocated before");
    }

  w = (max - min) / n_g;
  b = min - w/2; 
  e = max + w/2; 
  for (i = 0; i < n_g; i++) {
     m = i + 1;
     n = n_g - i;
     grid[i] = ((m * e) + (n * b)) / (n_g + 1);
  }
 
  return EXIT_SUCCESS;
}

int 
gmel_grid_new (double min_s, double max_s, int n_s,
               double min_p, double max_p, int n_p,
               double min_a, double max_a, int n_a,
               double min_c, double max_c, int n_c,
               double min_g, double max_g, int n_g,
               int gs)
{
  int r = EXIT_SUCCESS;
  int i, j, k, l, m, n;
  Gibbs_data *g = NULL;

  assert (min_s < max_s);
  assert (min_p < max_p);
  assert (min_a < max_a);
  assert (min_c < max_c);
  assert (min_g < max_g);
  assert (n_s > 0);
  assert (n_p > 0);
  assert (n_a > 0);
  assert (n_c > 0);
  assert (n_g > 0);

  gibbs_size = gs;
  gmel_grid_n[PSI_GRID_S] = n_s;
  gmel_grid_n[PSI_GRID_P] = n_p;
  gmel_grid_n[PSI_GRID_A] = n_a;
  gmel_grid_n[PSI_GRID_C] = n_c;
  gmel_grid_n[PSI_GRID_G] = n_g;
  grid_min[PSI_GRID_S] = min_s;
  grid_min[PSI_GRID_P] = min_p;
  grid_min[PSI_GRID_A] = min_a;
  grid_min[PSI_GRID_C] = min_c;
  grid_min[PSI_GRID_G] = min_g;
  grid_max[PSI_GRID_S] = max_s;
  grid_max[PSI_GRID_P] = max_p;
  grid_max[PSI_GRID_A] = max_a;
  grid_max[PSI_GRID_C] = max_c;
  grid_max[PSI_GRID_G] = max_g;

  /* Grid Point Memory Allocation */
  for (i = 0; i < PSI_GRID_G + 1; i++) 
    {
      gmel_grid_gridpoints[i] = XMALLOC (double, gmel_grid_n[i] + 1);
      /* Grid Point Assigment */
      r = gmel_grid_new_grid (gmel_grid_gridpoints[i], 
                              grid_min[i], grid_max[i], gmel_grid_n[i]);
      if (r != EXIT_SUCCESS) 
        {
          psi_fatal ("Error in new grid");
        }
    }

  if (number_multi_gibbs == 1)
    {
      /* Gibbs Informatoion Memory Allocation */
      gmel_grid_info = XMALLOC (GibbsPart *****, n_s + 1);
      for (i = 0; i < (n_s + 1); i++)
        {
          gmel_grid_info[i] = XMALLOC (GibbsPart ****, n_p + 1);
          for (j = 0; j < (n_p + 1); j++)
            {
              gmel_grid_info[i][j] = XMALLOC (GibbsPart ***, n_a + 1);
              for (k = 0; k < (n_a + 1); k++)
                {
                  gmel_grid_info[i][j][k] = XMALLOC(GibbsPart **, n_c + 1);
                  for (l = 0; l < (n_c + 1); l++)
                    {
                      gmel_grid_info[i][j][k][l] = XMALLOC (GibbsPart *, n_g + 1);
                      for (m = 0; m < (n_g + 1); m++) 
                        gmel_grid_info[i][j][k][l][m] = NULL;
                    }
                }
            }
        }
    }
  else
    {
      for (n = 0; n < number_multi_gibbs; n++)
        {
          gmel_grid_info = XMALLOC (GibbsPart *****, n_s + 1);
          for (i = 0; i < (n_s + 1); i++)
            {
              gmel_grid_info[i] = XMALLOC (GibbsPart ****, n_p + 1);
              for (j = 0; j < (n_p + 1); j++)
                {
                  gmel_grid_info[i][j] = XMALLOC (GibbsPart ***, n_a + 1);
                  for (k = 0; k < (n_a + 1); k++)
                    {
                      gmel_grid_info[i][j][k] = XMALLOC(GibbsPart **, n_c + 1);
                      for (l = 0; l < (n_c + 1); l++)
                        {
                          gmel_grid_info[i][j][k][l] = XMALLOC (GibbsPart *, n_g + 1);
                          for (m = 0; m < (n_g + 1); m++) 
                            gmel_grid_info[i][j][k][l][m] = NULL;
                        }
                    }
                }
            }
          g = gibbs_find (gibbs_list, n);
          if (g == NULL)
            psi_fatal ("no such gibbs: %d", n);
          g->info = gmel_grid_info;
          /* gmel_grid_info = NULL; */
        }
      g = gibbs_find (gibbs_list, 0);
      if (g == NULL)
        psi_fatal ("no such gibbs: %d", 0);
      gmel_grid_info = g->info;
      current_gibbs = 0;
    }

  return r;
}

int 
gmel_grid_load (const char* gname)
{
  FILE *gibbsfile;
  GibbsPart *sampled_seqs = NULL;
  GibbsPart *point_gibbs= NULL;
  int i, j, k, l, m, n;
  int r;
  double d_v;
  int i_v;
  
  gibbsfile = fopen (gname, "r");
  if (gibbsfile == NULL) {
     return EXIT_SUCCESS;
     /* psi_fatal ("grid load: could not open gibbs file"); */
  }
  /* check the gibbs setting is the same */
  fscanf (gibbsfile, "%d", &i_v);
  if (i_v != gibbs_size) {
     psi_fatal ("grid load: wrong gibbs size: %d %d\n", i_v, gibbs_size);
  }
  fscanf (gibbsfile, "%lf", &d_v);
  if (d_v != grid_min[PSI_GRID_S]) {
     psi_fatal ("grid load: wrong begin of grid point S");
  }
  fscanf (gibbsfile, "%lf", &d_v);
  if (d_v != grid_max[PSI_GRID_S]) {
     psi_fatal ("grid load: wrong end of grid point S");
  }
  fscanf (gibbsfile, "%d", &i_v);
  if (i_v != gmel_grid_n[PSI_GRID_S]) {
     psi_fatal ("grid load: wrong number of grid points S");
  }
  fscanf (gibbsfile, "%lf", &d_v);
  if (d_v != grid_min[PSI_GRID_P]) {
     psi_fatal ("grid load: wrong begin of grid point P");
  }
  fscanf (gibbsfile, "%lf", &d_v);
  if (d_v != grid_max[PSI_GRID_P]) {
     psi_fatal ("grid load: wrong end of grid point P");
  }
  fscanf (gibbsfile, "%d", &i_v);
  if (i_v != gmel_grid_n[PSI_GRID_P]) {
     psi_fatal ("grid load: wrong number of grid points P");
  }
  fscanf (gibbsfile, "%lf", &d_v);
  if (d_v != grid_min[PSI_GRID_A]) {
     psi_fatal ("grid load: wrong begin of grid point N");
  }
  fscanf (gibbsfile, "%lf", &d_v);
  if (d_v != grid_max[PSI_GRID_A]) {
     psi_fatal ("grid load: wrong end of grid point N");
  }
  r = fscanf (gibbsfile, "%d", &i_v);
  if (i_v != gmel_grid_n[PSI_GRID_A]) {
     psi_fatal ("grid load: wrong number of grid points N");
  }

  /* read one line */  
  while (r != EOF) {
     r = fscanf (gibbsfile, "%d%d%d%d%d", &i, &j, &k, &l, &m);
     if (r == EOF) {
        break;
     }
     /* allocate a memory */
     sampled_seqs = XMALLOC (GibbsPart, gibbs_size);
     point_gibbs = sampled_seqs; 
     /* copy the line into the memory */
     for (n = 0; n < gibbs_size; n++) {
        r = fscanf (gibbsfile, "%lf%lf%d%d%d", 
                           &(point_gibbs->S), &(point_gibbs->P),
                           &(point_gibbs->A), &(point_gibbs->C),
                           &(point_gibbs->G));
        if (want_debug_gibbs == 1) {
           fprintf (stderr, "Gibbs Read: %lf\t%lf\t%d\t%d\t%d\t", 
                    (point_gibbs->S), (point_gibbs->P),
                    (point_gibbs->A), (point_gibbs->C),
                    (point_gibbs->G));
        }
        point_gibbs++;
     }
     gmel_grid_info[i][j][k][l][m] = sampled_seqs;
     if (want_debug_gibbs == 1) {
        fprintf (stderr, "\n");
     }
  }
  
  fclose (gibbsfile);
  gibbsfile = NULL;

  return EXIT_SUCCESS;
}

int 
gmel_grid_save (const char* gname)
{
  FILE *gibbsfile;
  GibbsPart *sampled_seqs = NULL;
  int i, j, k, l, m, n;
  gibbsfile = fopen (gname, "w");
  if (gibbsfile == NULL) {
     psi_fatal ("cannot open gibbsfile");
  }
  /* write the gibbs setting */
  fprintf (gibbsfile, "%d\t%lf\t%lf\t%d\t%lf\t%lf\t%d\t%lf\t%lf\t%d\n", 
                      gibbs_size,
                      grid_min[PSI_GRID_S], grid_max[PSI_GRID_S], gmel_grid_n[PSI_GRID_S],
                      grid_min[PSI_GRID_P], grid_max[PSI_GRID_P], gmel_grid_n[PSI_GRID_P],
                      grid_min[PSI_GRID_A], grid_max[PSI_GRID_A], gmel_grid_n[PSI_GRID_A]);
  
  /* write all gibbs sample */
  for (i = 0; i < gmel_grid_n[PSI_GRID_S]; i++) { 
  for (j = 0; j < gmel_grid_n[PSI_GRID_P]; j++) {
  for (k = 0; k < gmel_grid_n[PSI_GRID_A]; k++) {
  for (l = 0; l < gmel_grid_n[PSI_GRID_C]; l++) {
  for (m = 0; m < gmel_grid_n[PSI_GRID_G]; m++) {
     if ( gmel_grid_info[i][j][k][l][m] != NULL ) {
        sampled_seqs = gmel_grid_info[i][j][k][l][m];
        fprintf (gibbsfile, "%d\t%d\t%d\t%d\t%d", i, j, k, l, m);
        for (n = 0; n < gibbs_size; n++)
        {
           fprintf (gibbsfile, "\t%lf\t%lf\t%d\t%d\t%d",
                               sampled_seqs[n].S,
                               sampled_seqs[n].P,
                               sampled_seqs[n].A,
                               sampled_seqs[n].C,
                               sampled_seqs[n].G);
        }
        fprintf (gibbsfile, "\n");
     }
  }
  }
  }
  }
  }
  
  fclose (gibbsfile);
  gibbsfile = NULL;
  return EXIT_SUCCESS;
}

int 
gmel_grid_del ()
{
  int i, j, k, l, m, n;

  if (number_multi_gibbs == 1)
    { 
      /* Gibbs Informatoion Memory Deallocation */
      for (i = 0; i < (gmel_grid_n[PSI_GRID_S] + 1); i++) { 
      for (j = 0; j < (gmel_grid_n[PSI_GRID_P] + 1); j++) {
      for (k = 0; k < (gmel_grid_n[PSI_GRID_A] + 1); k++) {
      for (l = 0; l < (gmel_grid_n[PSI_GRID_C] + 1); l++) {
      for (m = 0; m < (gmel_grid_n[PSI_GRID_G] + 1); m++) {
          XFREE (gmel_grid_info[i][j][k][l][m]);
      }
          XFREE (gmel_grid_info[i][j][k][l]);
      }
          XFREE (gmel_grid_info[i][j][k]);
      }
          XFREE (gmel_grid_info[i][j]);
      }
          XFREE (gmel_grid_info[i]);
      }
      XFREE(gmel_grid_info);
    }
  else
    {
      for (n = 0; n < number_multi_gibbs; n++)
        {
          Gibbs_data *g = gibbs_find (gibbs_list, n);
          if (g == NULL)
            psi_fatal ("no such gibbs: %d", n);
          gmel_grid_info = g->info;

          for (i = 0; i < (gmel_grid_n[PSI_GRID_S] + 1); i++) { 
          for (j = 0; j < (gmel_grid_n[PSI_GRID_P] + 1); j++) {
          for (k = 0; k < (gmel_grid_n[PSI_GRID_A] + 1); k++) {
          for (l = 0; l < (gmel_grid_n[PSI_GRID_C] + 1); l++) {
          for (m = 0; m < (gmel_grid_n[PSI_GRID_G] + 1); m++) {
              XFREE (gmel_grid_info[i][j][k][l][m]);
          }
              XFREE (gmel_grid_info[i][j][k][l]);
          }
              XFREE (gmel_grid_info[i][j][k]);
          }
              XFREE (gmel_grid_info[i][j]);
          }
              XFREE (gmel_grid_info[i]);
          }
          XFREE(gmel_grid_info);
        }
      gibbs_remove (&gibbs_list);
    }

  /* Grid Point Memory Deallocation */
  free_all_gmel_grid_gridpoints ();

  return EXIT_SUCCESS;
}

static void 
free_all_gmel_grid_gridpoints (void) {
  XFREE (gmel_grid_gridpoints[PSI_GRID_S]);
  XFREE (gmel_grid_gridpoints[PSI_GRID_P]);
  XFREE (gmel_grid_gridpoints[PSI_GRID_A]);
  XFREE (gmel_grid_gridpoints[PSI_GRID_C]);
  XFREE (gmel_grid_gridpoints[PSI_GRID_G]);
}

static int 
gmel_grid_locate_a_parameter (int *index, double param_val, int which_grid)
{
  /*

   --- !!! R I G H T !!! ---

        min                                                             max
            0       1       2       3       4       5       6       7 --> n : 8
   GRID b---|-------|-------|-------|-------|-------|-------|-------|---e
   VALUE   ---(------v------O------v------O------v------)-------
   
   w        |<----->|
   W    |<---------->|

   --- !!! W R O N G !!! ---

                0       1       2       3       4       5       6    --> n : 7
   GRID b-------|-------|-------|-------|-------|-------|-------|-------e
   VALUE   ---(------v------O------v------O------v------)-------
   
   w            |<----->|
   W    |<---------->|

   The index of gridpoint near the first middle point (v) is
   W/w's quotient.
   */


  double w, W, b, v;
  int i, n;

  assert (which_grid >= 0);
  assert (which_grid < PSI_GRID_G + 1);
  assert (gmel_grid_n[which_grid] > 1);
  n = gmel_grid_n[which_grid];
  w = (gmel_grid_gridpoints[which_grid][n - 1] - gmel_grid_gridpoints[which_grid][0]) / (n - 1);
  b = gmel_grid_gridpoints[which_grid][0] - w/2;
  v = param_val;
  W = v - b; 
  i = (int) floor (W/w);

  if (i < 0) {
     if (want_verbose == 1) {
        fprintf (stderr, "%d's grid value of %lf is out of range: %d\n",
                 which_grid, param_val, i);
     }
     i = 0;
  } else if (i >= n) {
     if (want_verbose == 1) {
        fprintf (stderr, "%d's grid value of %lf is out of range: %d\n",
                 which_grid, param_val, i);
     }
     i = n - 1;
  }

  assert (i >= 0 && i < n);
  if (i >= n) {
     gmel_grid_print_grid ();
     psi_fatal ("grid point is out of bound");
  }
  *index = i;

  return 0;
}

int 
gmel_grid_locate (gridpoint *gp, parameter param)
{
  gmel_grid_locate_a_parameter (&(gp->s), param.s, PSI_GRID_S); 
  gmel_grid_locate_a_parameter (&(gp->p), param.p, PSI_GRID_P); 
  gmel_grid_locate_a_parameter (&(gp->a), param.a, PSI_GRID_A); 
  gmel_grid_locate_a_parameter (&(gp->c), param.c, PSI_GRID_C); 
  gmel_grid_locate_a_parameter (&(gp->g), param.g, PSI_GRID_G); 

  /* A, C, G Grid Point Boundary Check */
  while ( gmel_grid_gridpoints[PSI_GRID_A][gp->a] 
          + gmel_grid_gridpoints[PSI_GRID_C][gp->c]
          + gmel_grid_gridpoints[PSI_GRID_G][gp->g] >= 1.0 ) {                               
     double diff_A, diff_C, diff_G;                                         
     if (gp->a > 0) {
        diff_A = param.a - gmel_grid_gridpoints[PSI_GRID_A][gp->a - 1];
     } else {
        diff_A = 1;
     }
     if (gp->c > 0) {
           diff_C = param.c - gmel_grid_gridpoints[PSI_GRID_C][gp->c - 1]; 
     } else {
        diff_C = 1;
     }
     if (gp->g > 0) {
        diff_G = param.g - gmel_grid_gridpoints[PSI_GRID_G][gp->g - 1];
     } else {
        diff_G = 1;
     }
     /* which one is the smallest? */                                       
     if (diff_A < diff_C) {
        if (diff_A < diff_G) {                                              
           (gp->a)--;
        } else {
           (gp->g)--;
        }
     } else {
        if (diff_C < diff_G) {                                              
           (gp->c)--;
        } else { 
           (gp->g)--;
        }
     }
     assert (gp->a >= 0 && gp->a < gmel_grid_n[PSI_GRID_A]);
     assert (gp->c >= 0 && gp->c < gmel_grid_n[PSI_GRID_C]);
     assert (gp->g >= 0 && gp->g < gmel_grid_n[PSI_GRID_G]);
  }

  {
     double a, c, g;
     a = gmel_grid_gridpoints[PSI_GRID_A][gp->a];
     c = gmel_grid_gridpoints[PSI_GRID_C][gp->c];
     g = gmel_grid_gridpoints[PSI_GRID_G][gp->g];
     assert (a + c + g < 1.0);
  }
 
  return 0;
}

int 
gmel_grid_param_gridpoint (parameter *param, gridpoint gp)
{
  param->s = gmel_grid_gridpoints[PSI_GRID_S][gp.s];
  param->p = gmel_grid_gridpoints[PSI_GRID_P][gp.p];
  param->a = gmel_grid_gridpoints[PSI_GRID_A][gp.a];
  param->c = gmel_grid_gridpoints[PSI_GRID_C][gp.c];
  param->g = gmel_grid_gridpoints[PSI_GRID_G][gp.g];
  param->t = 1.0 - (param->a + param->c + param->g);
  assert (param->a + param->c + param->g < 1.0);
  return 0;
}

int
gmel_grid_fix (int grid_type, double value)
{
  gmel_grid_gridpoints[grid_type][0] = value;
  return EXIT_SUCCESS;
}

int 
gmel_grid_put_info (GibbsPart *gibbs, gridpoint gp)
{
  gmel_grid_info[gp.s][gp.p][gp.a][gp.c][gp.g] = gibbs;
  return EXIT_SUCCESS;
}

int 
gmel_grid_get_info (GibbsPart **gibbs, gridpoint gp)
{
  *gibbs = gmel_grid_info[gp.s][gp.p][gp.a][gp.c][gp.g];
  return EXIT_SUCCESS;
}

int
gmel_grid_check_sparse (parameter top, parameter bot)
{
  int r;
  r = gmel_grid_check_sparse_a_parameter(top.s, bot.s, PSI_GRID_S);
  if (r == 1) {
     return 1;
  }
  r = gmel_grid_check_sparse_a_parameter(top.p, bot.p, PSI_GRID_P);
  if (r == 1) {
     return 1;
  }
  r = gmel_grid_check_sparse_a_parameter(top.a, bot.a, PSI_GRID_A);
  if (r == 1) {
     return 1;
  }
  r = gmel_grid_check_sparse_a_parameter(top.c, bot.c, PSI_GRID_C);
  if (r == 1) {
     return 1;
  }
  r = gmel_grid_check_sparse_a_parameter(top.g, bot.g, PSI_GRID_G);
  if (r == 1) {
     return 1;
  }

  return 0;
}

static int
gmel_grid_check_sparse_a_parameter (double top, double bot, int which_grid)
{
  int top_i, bot_i;

  top_i = gmel_grid_check_sparse_locate (top, which_grid); 
  bot_i = gmel_grid_check_sparse_locate (bot, which_grid); 

  if (top_i == bot_i) {
     return 1;
  } else {
     return 0;
  } 
} 

static int
gmel_grid_check_sparse_locate (double val, int which_grid)
{
  double w, W, b, v;
  int i, n;

  n = gmel_grid_n[which_grid];
  w = (gmel_grid_gridpoints[which_grid][n - 1] - gmel_grid_gridpoints[which_grid][0]) / (n - 1);
  b = gmel_grid_gridpoints[which_grid][0];
  v = val;
  W = v - b; 
  i = floor (W/w);

  return i;
}

int 
gmel_grid_print_grid ()
{
  int i;
  for (i = 0; i < 5; i++) {
     gmel_grid_print_gridpoints (i);
  }
  return 0;
}

static int
gmel_grid_print_gridpoints (int which_grid)
{
  int i, n;
  n = gmel_grid_n[which_grid];
  switch (which_grid) {
  case PSI_GRID_S:
     fprintf(stderr, "S: ");
     break; 
  case PSI_GRID_P:
     fprintf(stderr, "P: ");
     break; 
  case PSI_GRID_A:
     fprintf(stderr, "A: ");
     break; 
  case PSI_GRID_C:
     fprintf(stderr, "C: ");
     break; 
  case PSI_GRID_G:
     fprintf(stderr, "G: ");
     break; 
  default:
     assert(0);
  }
  for (i = 0; i < n; i++) {
     fprintf(stderr, "%.3lf ", gmel_grid_gridpoints[which_grid][i]);
  } 
  fprintf(stderr, "\n");
  return 0;
}



