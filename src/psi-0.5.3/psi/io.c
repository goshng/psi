/* io.c -- I/O
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
#include "rng.h"
#include "io.h"

static void next_line (int upto, FILE *fp);

void
write_parameter_header (FILE *fp)
{
  fprintf (fp, "Gen\tA\tC\tG\tT\tS\tP\n");
}

void 
write_parameter (FILE *fp, parameter param, int gen) 
{
  fprintf (fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
           gen, param.a, param.c, param.g, param.t, param.s, param.p);
}

void 
write_parameter_alldraws (FILE *fp, parameter param, int gen) 
{
  static int c = 0; 
  if (gen == GEN_INIT) {
     c = 0;
     fprintf (fp, "Gen\tA\tC\tG\tT\tS\tP\n");
  }
  c++;
  fprintf (fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
           c - 1, param.a, param.c, param.g, param.t, param.s, param.p);
}


void 
debug_write_parameter (FILE *fp, parameter param) {
  static int c = 0; 
  c++;
  fprintf (fp, "D:%9d %+.6lf %+.6lf %+.6lf %+.6lf %+.6lf %+.6lf %+.6lf %+.6lf ", 
           c, param.a, param.c, param.g, param.t, param.s, param.p, 
           param.max_energy, param.gibbs_sum);
  debug_write_gp (fp, param.gp);
}

void 
debug_write_gp (FILE *fp, gridpoint gp) {
  fprintf (fp, "%3d %3d %3d %3d %3d\n", gp.a, gp.c, gp.g, gp.s, gp.p);
}

int fgetline (char *line[], int max, FILE *f) {
  if (fgets((char *)line, max, f) == NULL)
     return 0;
  else
     return strlen((char *)line);
}

static void next_line (int upto, FILE *fp) {
  int i;
  char *line[PSI_IO_GMEL_LINE_MAX];
  for (i = 0; i < upto; i++) {
     fgetline(line, PSI_IO_GMEL_LINE_MAX, fp);
 /* fprintf(stderr, "NEXT LINE: %s", line); */
  }
}

void locate_m1_block1 ( FILE *fp, int sample_size ) {
  int upto = (sample_size + 1) * 3 + 1;
  next_line (upto, fp);
}

void locate_m4_block1_e1 ( FILE *fp, int sample_size ) {
  int upto = (sample_size + 1) * 0 + 1;
  next_line (upto, fp);
}

void locate_m4_block1_e2 ( FILE *fp, int sample_size ) {
  int upto = (sample_size + 1) * 1 + 1;
  next_line (upto, fp);
}

void locate_m4_block2_e1 ( FILE *fp, int sample_size ) {
  int upto = (sample_size + 1) * 1 + 1;
  next_line (upto, fp);
}

void locate_m4_block2_e2 ( FILE *fp, int sample_size ) {
  int upto = (sample_size + 1) * 2 + 1;
  next_line (upto, fp);
}

void locate_m4_block3_e1 ( FILE *fp, int sample_size ) {
  int upto = (sample_size + 1) * 2 + 1;
  next_line (upto, fp);
}

void pass_lines (int n, FILE *fp) {
  next_line (n, fp);
}

int 
read_param_line (FILE *bf, int *gen, 
                 double *a, double *c, double *g, double *t, 
                 double *s, double *p) 
{
  int r = EXIT_SUCCESS;
  r = fscanf (bf, 
              "%d%lf%lf%lf%lf%lf%lf",  /* BUGGY '\n'? */
              gen, a, c, g, t, s, p);
  if (r == EOF || r != 7) {
     fprintf (stderr, "%s:%d:%s; %s %d %lf %lf %lf %lf %lf %lf - r: %d\n",
              __FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "error of reading posterior draws",
              *gen, *a, *c, *g, *t, *s, *p, r);
     /* r = PSI_ERR_READ; */
     r = 1008;
  } else {
     r = EXIT_SUCCESS;
  }
  return r;
}

int 
read_parameter (FILE *fp, parameter *p_param) 
{
  int r = EXIT_SUCCESS;
  int gen;
  double a, c, g, t, s, p;

  r = read_param_line (fp, &gen, &a, &c, &g, &t, &s, &p);
  if (r != EXIT_SUCCESS) {
     fprintf (stderr, "%s:%d:%s; %s %d %lf %lf %lf %lf %lf %lf - r: %d\n",
              __FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "error of reading posterior draws of block 1 in model 1",
              gen, a, c, g, t, s, p, r);
     psi_fatal ("error of reading posterior draws of block 1 in model 1");
  }
  p_param->a = a;
  p_param->c = c;
  p_param->g = g;
  p_param->t = t;
  p_param->s = s;
  p_param->p = p;

  return EXIT_SUCCESS;
}

void
verbose_likelihood_gibbs (FILE *fp, double likelihood, parameter top, parameter bot)
{
  fprintf (fp, "Ratio %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf vs. %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf : %.3lf\n", 
           top.s, top.p, top.a, top.c, top.g, top.t, 
           bot.s, bot.p, bot.a, bot.c, bot.g, bot.t, 
           likelihood);
  return;
}

void
verbose_gibbs (FILE *fp, double ratio, 
               int n_a, int n_c, int n_g, int n_t, 
               double LnDiff_A, double LnDiff_C, 
               double LnDiff_G, double LnDiff_T, 
               double solv, double pair,
               double gibbs_top, double gibbs_bot,
               double top_max_energy, double bot_max_energy,
               parameter top, parameter bot)
{
     fprintf (fp, "gibbs: %.3lf : %lf - %lf ; %lf ; \
                      %d %.3lf = %.3lf, %d %.3lf = %.3lf, %d %.3lf = %.3lf, \
                      %d %.3lf = %.3lf -> %.3lf ; \
                      S %.3lf * %.3lf - %.3lf = %.3lf, \
                      P %.3lf * %.3lf - %.3lf = %.3lf, \
                      Max %.3lf %.3lf = %.3lf\n",
              ratio,
              log(gibbs_top), log(gibbs_bot), log(gibbs_top) - log(gibbs_bot), 
              n_a, LnDiff_A, n_a * LnDiff_A, 
              n_c, LnDiff_C, n_c * LnDiff_C, 
              n_g, LnDiff_G, n_g * LnDiff_G, 
              n_t, LnDiff_T, n_t * LnDiff_T, 
              n_a * LnDiff_A + n_c * LnDiff_C + n_g * LnDiff_G + n_t * LnDiff_T, 
              solv, top.s, bot.s, solv * (top.s - bot.s),
              pair, top.p, bot.p, pair * (top.p - bot.p),
              top_max_energy, bot_max_energy, top_max_energy - bot_max_energy);
}

int
gmel_read_draws (char *chain_name, double *draws, int w, int sample_size)
{
  int r = EXIT_SUCCESS;
  int i;
  double draw[6];
  int gen;
  FILE *pfile = NULL;

  pfile = fopen (chain_name, "r");
  if (pfile == NULL) {
     fprintf (stderr, "%s:%d:%s; %s\n",
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
     psi_fatal ("could open file for read: gmel_read_draw");
  }
  locate_m4_block1_e1 (pfile, sample_size);
  for (i = 0; i < sample_size; i++) {
     r = fscanf (pfile, "%d%lf%lf%lf%lf%lf%lf", 
                 &gen, 
                 &draw[PSI_SAMPLE_A],
                 &draw[PSI_SAMPLE_C],
                 &draw[PSI_SAMPLE_G],
                 &draw[PSI_SAMPLE_T],
                 &draw[PSI_SAMPLE_S],
                 &draw[PSI_SAMPLE_P]);
     if (r == EOF || r != 7) {
        fprintf (stderr, "%s:%d:%s; %s %d %lf %lf %lf %lf %lf %lf - r: %d\n",
                 __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                 "error of reading posterior draws of block 1 in model 4",
                 gen, draw[PSI_SAMPLE_A], draw[PSI_SAMPLE_C], 
                 draw[PSI_SAMPLE_G], draw[PSI_SAMPLE_T], 
                 draw[PSI_SAMPLE_S], draw[PSI_SAMPLE_P], r);
        psi_fatal ("error in reading posterior draws of block 1 in model 4");
     }
     draws[i] = draw[w]; 
  }
  
  return EXIT_SUCCESS;
}

int
psi_shuffle_file (const char *old_name, const char *new_name, int sample_size)
{
  int r = EXIT_SUCCESS;
  FILE *old_file = NULL;
  FILE *new_file = NULL;
  parameter *list_params = NULL;
  int i;
  int j;
  double u;

  list_params = XMALLOC (parameter, sample_size);

  old_file = fopen (old_name, "r");
  if (old_file == NULL) {
     fprintf (stderr, "%s:%s:%d:%s; %s\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
     psi_fatal ("old_file open error in psi_shuffle_file");
  }

  new_file = fopen (new_name, "w");
  if (new_file == NULL) {
     fprintf (stderr, "%s:%s:%d:%s; %s\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
     psi_fatal ("new_file open error in psi_shuffle_file");
  }

  /* Model 4 - Block 1 */
  pass_lines (1, old_file);
  for (i = 0; i < sample_size; i++) {
     read_parameter (old_file, &list_params[i]);
  }

  i = 0;
  u = psi_rng ();
  j = (int) (sample_size * u);
  write_parameter (new_file, list_params[j], GEN_INIT);
  for (i = 1; i < sample_size; i++) {
     u = psi_rng ();
     j = (int) (sample_size * u);
     write_parameter (new_file, list_params[j], GEN_NEXT);
  }
  
  /* Model 4 - Block 2 */
  pass_lines (1 + 1, old_file);
  for (i = 0; i < sample_size; i++) {
     read_parameter (old_file, &list_params[i]);
  }

  i = 0;
  u = psi_rng ();
  j = (int) (sample_size * u);
  write_parameter (new_file, list_params[j], GEN_INIT);
  for (i = 1; i < sample_size; i++) {
     u = psi_rng ();
     j = (int) (sample_size * u);
     write_parameter (new_file, list_params[j], GEN_NEXT);
  }
  
  /* Model 4 - Block 3 */
  pass_lines (1 + 1, old_file);
  for (i = 0; i < sample_size; i++) {
     read_parameter (old_file, &list_params[i]);
  }

  i = 0;
  u = psi_rng ();
  j = (int) (sample_size * u);
  write_parameter (new_file, list_params[j], GEN_INIT);
  for (i = 1; i < sample_size; i++) {
     u = psi_rng ();
     j = (int) (sample_size * u);
     write_parameter (new_file, list_params[j], GEN_NEXT);
  }
 
  /* Model 1 - Block 1 */
  pass_lines (1 + 1, old_file);
  for (i = 0; i < sample_size; i++) {
     read_parameter (old_file, &list_params[i]);
  }

  i = 0;
  u = psi_rng ();
  j = (int) (sample_size * u);
  write_parameter (new_file, list_params[i], GEN_INIT);
  for (i = 1; i < sample_size; i++) {
     u = psi_rng ();
     j = (int) (sample_size * u);
     write_parameter (new_file, list_params[j], GEN_NEXT);
  }

  fclose (old_file);
  old_file = NULL;
  fclose (new_file);
  new_file = NULL;
  free (list_params);
  list_params = NULL;
  return r; 
}

void
psi_mcmc_write_posterior_sample (char *pn, const double ***post, int n, int s)
{
  int i, j;
  char *filename = NULL;
  int len_filename = 0;
  FILE *fp = NULL;
  parameter theta;
  len_filename = strlen (pn) + 1 + 1;
  filename = XMALLOC (char, len_filename);

  for (i = 0; i < n; i++)
    {
      sprintf (filename, "%s%d", pn, i + 1);
      fp = fopen (filename, "w");
      write_parameter_header (fp);
      for (j = 0; j < s; j++)
        {
          theta.a = post[i][PSI_SAMPLE_A][j];
          theta.c = post[i][PSI_SAMPLE_C][j];
          theta.g = post[i][PSI_SAMPLE_G][j];
          theta.t = post[i][PSI_SAMPLE_T][j];
          theta.s = post[i][PSI_SAMPLE_S][j];
          theta.p = post[i][PSI_SAMPLE_P][j];
          write_parameter (fp, theta, j);
        }
      fclose (fp);
    } 
  XFREE (filename);
}

int
psi_io_check_est_p (char *pname, int sample_size)
{
  int r = EXIT_SUCCESS;
  int i;
  double draw[6];
  int gen;
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    return EXIT_FAILURE; /* no output file, so it's okay to start */

  next_line (1, fp);
  for (i = 0; i < sample_size; i++) 
    {
      r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                  &gen, 
                  &draw[PSI_SAMPLE_A],
                  &draw[PSI_SAMPLE_C],
                  &draw[PSI_SAMPLE_G],
                  &draw[PSI_SAMPLE_T],
                  &draw[PSI_SAMPLE_S],
                  &draw[PSI_SAMPLE_P]);
      if (r == EOF || r != 7) 
        {
          fclose (fp);
          return EXIT_FAILURE; 
        }
    }
  fclose (fp);
  return EXIT_SUCCESS;
}

int
psi_io_check_bf_p (char *pname, int sample_size)
{
  int r = EXIT_SUCCESS;
  int i, j;
  double draw[6];
  int gen;
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    return EXIT_FAILURE; /* no output file, so it's okay to start */

  for (i = 0; i < 4; i++) 
    {
      next_line (1, fp);
      for (j = 0; j < sample_size; j++) 
        {
          r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                      &gen, 
                      &draw[PSI_SAMPLE_A],
                      &draw[PSI_SAMPLE_C],
                      &draw[PSI_SAMPLE_G],
                      &draw[PSI_SAMPLE_T],
                      &draw[PSI_SAMPLE_S],
                      &draw[PSI_SAMPLE_P]);
          if (r == EOF || r != 7) 
            {
              fclose (fp);
              return EXIT_FAILURE; 
            }
        }
    }
  fclose (fp);
  return EXIT_SUCCESS;
}

int
psi_io_check_bf_state (char *pname, int sample_size)
{
  int r = EXIT_SUCCESS;
  int i, j;
  double draw[6];
  int gen;
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    return 0;

  for (i = 0; i < 4; i++) 
    {
      next_line (1, fp);
      for (j = 0; j < sample_size; j++) 
        {
          r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                      &gen, 
                      &draw[PSI_SAMPLE_A],
                      &draw[PSI_SAMPLE_C],
                      &draw[PSI_SAMPLE_G],
                      &draw[PSI_SAMPLE_T],
                      &draw[PSI_SAMPLE_S],
                      &draw[PSI_SAMPLE_P]);
          if (r == EOF || r != 7) 
            {
              fclose (fp);
              return i;
            }
        }
    }
  fclose (fp);
  return i;
}

int
psi_bf_stage_remained_sample_size (char *pname, int sample_size)
{
  int r = EXIT_SUCCESS;
  int i, j;
  double draw[6];
  int gen;
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    return sample_size;

  for (i = 0; i < 4; i++) 
    {
      next_line (1, fp);
      for (j = 0; j < sample_size; j++) 
        {
          r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                      &gen, 
                      &draw[PSI_SAMPLE_A],
                      &draw[PSI_SAMPLE_C],
                      &draw[PSI_SAMPLE_G],
                      &draw[PSI_SAMPLE_T],
                      &draw[PSI_SAMPLE_S],
                      &draw[PSI_SAMPLE_P]);
          if (r == EOF || r != 7) 
            {
              fclose (fp);
              return (sample_size - j);
            }
        }
    }
  fclose (fp);
  return 0;
}

parameter 
psi_bf_stage_theta_star (char *pname, int sample_size)
{
  int r = EXIT_SUCCESS;
  int j;
  double draw[6];
  int gen;
  parameter theta;
  theta.s = 0.0;
  theta.p = 0.0;
  theta.a = 0.0;
  theta.c = 0.0;
  theta.g = 0.0;
  theta.t = 0.0;
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    psi_fatal ("there should be parameter file to get theta star, this function should be called only when there is a complete set of estimation");

  next_line (1, fp);
  for (j = 0; j < sample_size; j++) 
    {
      r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                  &gen, 
                  &draw[PSI_SAMPLE_A],
                  &draw[PSI_SAMPLE_C],
                  &draw[PSI_SAMPLE_G],
                  &draw[PSI_SAMPLE_T],
                  &draw[PSI_SAMPLE_S],
                  &draw[PSI_SAMPLE_P]);
      if (r == EOF || r != 7) 
        {
          if (j > 0) 
            {
              theta.a /= j;
              theta.c /= j;
              theta.g /= j;
              theta.t = 1.0 - theta.a - theta.c - theta.g;
              theta.s /= j;
              theta.p /= j;
            }
          fclose (fp);
          return theta;
        }
      theta.a += draw[PSI_SAMPLE_A];
      theta.c += draw[PSI_SAMPLE_C];
      theta.g += draw[PSI_SAMPLE_G];
      theta.t += draw[PSI_SAMPLE_T];
      theta.s += draw[PSI_SAMPLE_S];
      theta.p += draw[PSI_SAMPLE_P];
    }
  fclose (fp);
  theta.a /= sample_size;
  theta.c /= sample_size;
  theta.g /= sample_size;
  theta.t = 1.0 - theta.a - theta.c - theta.g;
  theta.s /= sample_size;
  theta.p /= sample_size;
  return theta;
}

parameter 
psi_bf_stage_theta_last (char *pname, int sample_size)
{
  int r = EXIT_SUCCESS;
  int i, j;
  double draw[6];
  int gen;
  parameter theta;
  theta.s = 0.0;
  theta.p = 0.0;
  theta.a = 0.25;
  theta.c = 0.25;
  theta.g = 0.25;
  theta.t = 0.25;
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    return theta;

  for (i = 0; i < 4; i++) 
    {
      next_line (1, fp);
      for (j = 0; j < sample_size; j++) 
        {
          r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                      &gen, 
                      &draw[PSI_SAMPLE_A],
                      &draw[PSI_SAMPLE_C],
                      &draw[PSI_SAMPLE_G],
                      &draw[PSI_SAMPLE_T],
                      &draw[PSI_SAMPLE_S],
                      &draw[PSI_SAMPLE_P]);
          if (r == EOF || r != 7) 
            {
              fclose (fp);
              return theta;
            }
          theta.a = draw[PSI_SAMPLE_A];
          theta.c = draw[PSI_SAMPLE_C];
          theta.g = draw[PSI_SAMPLE_G];
          theta.t = draw[PSI_SAMPLE_T];
          theta.s = draw[PSI_SAMPLE_S];
          theta.p = draw[PSI_SAMPLE_P];
        }
    }
  fclose (fp);
  return theta;
}

int 
psi_io_index_last_theta (char *pname, int sample_size)
{
  int r = EXIT_SUCCESS;
  int i, j;
  double draw[6];
  int gen;
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    return 0;

  for (i = 0; i < 4; i++) 
    {
      next_line (1, fp);
      for (j = 0; j < sample_size; j++) 
        {
          r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                      &gen, 
                      &draw[PSI_SAMPLE_A],
                      &draw[PSI_SAMPLE_C],
                      &draw[PSI_SAMPLE_G],
                      &draw[PSI_SAMPLE_T],
                      &draw[PSI_SAMPLE_S],
                      &draw[PSI_SAMPLE_P]);
          if (r == EOF || r != 7) 
            {
              fclose (fp);
              return j;
            }
        }
    }
  fclose (fp);
  psi_fatal ("should not be called");
}



int
psi_io_check_bf_b (char *bname)
{
  int r = EXIT_SUCCESS;
  double draw[14];
  FILE *fp = NULL;

  fp = fopen (bname, "r");
  
  if (fp == NULL)
    return EXIT_FAILURE; /* no output file, so it's okay to start */

  r = fscanf (fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n", 
              &draw[0], &draw[1], &draw[2], &draw[3], &draw[4],
              &draw[5], &draw[6], &draw[7], &draw[8], &draw[9],
              &draw[10], &draw[11], &draw[12], &draw[13]);
  if (r == EOF || r != 14) 
    {
      fclose (fp);
      return EXIT_FAILURE; 
    }
  fclose (fp);
  return EXIT_SUCCESS;
}

int
psi_io_write_number_burn (FILE *fn, int i)
{
  if (fn == NULL)
    psi_fatal ("no file pointer for the parameter");
  fprintf (fn, "### B U R N ###\n");
  return EXIT_SUCCESS;
}

int 
psi_io_read_gen_last_burn (char *pname, int *gen)
{
  int r;
  double draw[6];
  char line[PSI_IO_GMEL_LINE_MAX];
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    return EXIT_FAILURE;

  while (fgets(line, PSI_IO_GMEL_LINE_MAX, fp) != NULL)
    {
      /* fprintf (stderr, "%s", line); */
      if (!strcmp (line, "### B U R N ###\n"))
        {
           r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                      gen, 
                      &draw[PSI_SAMPLE_A],
                      &draw[PSI_SAMPLE_C],
                      &draw[PSI_SAMPLE_G],
                      &draw[PSI_SAMPLE_T],
                      &draw[PSI_SAMPLE_S],
                      &draw[PSI_SAMPLE_P]);
          if (r == EOF || r != 7) 
            {
              fclose (fp);
              return EXIT_FAILURE;
            }
          fclose (fp);
          return EXIT_SUCCESS;  
        }
    }

  fclose (fp);
  return EXIT_FAILURE;  
}

int
psi_io_read_theta_last_burn (char *pname, parameter *theta)
{
  int r;
  double draw[6];
  int gen;
  char line[PSI_IO_GMEL_LINE_MAX];
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    return EXIT_FAILURE;

  while (fgets(line, PSI_IO_GMEL_LINE_MAX, fp) != NULL)
    {
      /* fprintf (stderr, "%s", line); */
      if (!strcmp (line, "### B U R N ###\n"))
        {
           r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                      &gen, 
                      &draw[PSI_SAMPLE_A],
                      &draw[PSI_SAMPLE_C],
                      &draw[PSI_SAMPLE_G],
                      &draw[PSI_SAMPLE_T],
                      &draw[PSI_SAMPLE_S],
                      &draw[PSI_SAMPLE_P]);
          if (r == EOF || r != 7) 
            {
              fclose (fp);
              psi_fatal ("there should be burn state");
            }
          theta->a = draw[PSI_SAMPLE_A];
          theta->c = draw[PSI_SAMPLE_C];
          theta->g = draw[PSI_SAMPLE_G];
          theta->t = draw[PSI_SAMPLE_T];
          theta->s = draw[PSI_SAMPLE_S];
          theta->p = draw[PSI_SAMPLE_P];
          fclose (fp);
          return EXIT_SUCCESS;  
        }
    }

  fclose (fp);
  /* not find B U R N  */
  return EXIT_FAILURE;  
}

int
psi_io_remove_number_burn (char *pname)
{
  char c;
  FILE *fp = NULL;
  FILE *fp_t = NULL;
  char *pname_t = NULL; 
  pname_t = XMALLOC (char, strlen (pname) + 1);
  strcpy (pname_t, pname); 
  pname_t[strlen (pname) - 1] = 'x';
  fprintf (stderr, "%s\n", pname_t);
  fp_t = fopen (pname_t, "w");
  fp = fopen (pname, "r");
  while ((c = fgetc (fp)) != EOF)
    {
      fputc (c, fp_t);
    }
  fclose (fp);
  fclose (fp_t);

  fp_t = fopen (pname_t, "r");
  fp = fopen (pname, "w");
  while ((c = fgetc (fp_t)) != '#')
    {
      fputc (c, fp);
    }
  fclose (fp);
  fclose (fp_t);

  unlink (pname_t);
  XFREE (pname_t);
  return EXIT_SUCCESS;
}

int
psi_io_write_number_sample (FILE *fn, int i)
{
  if (fn == NULL)
    psi_fatal ("no file pointer for the parameter");
  fprintf (fn, "$$$ S A M P L E $$$\n");
  return EXIT_SUCCESS;
}


int 
psi_io_read_gen_last_sample (char *pname, int *gen)
{
  int r;
  double draw[6];
  char line[PSI_IO_GMEL_LINE_MAX];
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    return EXIT_FAILURE;

  while (fgets(line, PSI_IO_GMEL_LINE_MAX, fp) != NULL)
    {
      /* fprintf (stderr, "%s", line); */
      if (!strcmp (line, "$$$ S A M P L E $$$\n"))
        {
           r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                      gen, 
                      &draw[PSI_SAMPLE_A],
                      &draw[PSI_SAMPLE_C],
                      &draw[PSI_SAMPLE_G],
                      &draw[PSI_SAMPLE_T],
                      &draw[PSI_SAMPLE_S],
                      &draw[PSI_SAMPLE_P]);
          if (r == EOF || r != 7) 
            {
              fclose (fp);
              return EXIT_FAILURE;
            }
          fclose (fp);
          return EXIT_SUCCESS;  
        }
    }

  fclose (fp);
  return EXIT_FAILURE;  
}

int
psi_io_read_theta_last_sample (char *pname, parameter *theta)
{
  int r;
  double draw[6];
  int gen;
  char line[PSI_IO_GMEL_LINE_MAX];
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    return EXIT_FAILURE;

  while (fgets(line, PSI_IO_GMEL_LINE_MAX, fp) != NULL)
    {
      /* fprintf (stderr, "%s", line); */
      if (!strcmp (line, "$$$ S A M P L E $$$\n"))
        {
           r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                      &gen, 
                      &draw[PSI_SAMPLE_A],
                      &draw[PSI_SAMPLE_C],
                      &draw[PSI_SAMPLE_G],
                      &draw[PSI_SAMPLE_T],
                      &draw[PSI_SAMPLE_S],
                      &draw[PSI_SAMPLE_P]);
          if (r == EOF || r != 7) 
            {
              fclose (fp);
              psi_fatal ("there should be burn state");
            }
          theta->a = draw[PSI_SAMPLE_A];
          theta->c = draw[PSI_SAMPLE_C];
          theta->g = draw[PSI_SAMPLE_G];
          theta->t = draw[PSI_SAMPLE_T];
          theta->s = draw[PSI_SAMPLE_S];
          theta->p = draw[PSI_SAMPLE_P];
          fclose (fp);
          return EXIT_SUCCESS;  
        }
    }

  fclose (fp);
  /* not find B U R N  */
  return EXIT_FAILURE;  
}

int
psi_io_remove_number_sample (char *pname)
{
  char c;
  FILE *fp = NULL;
  FILE *fp_t = NULL;
  char *pname_t = NULL; 
  pname_t = XMALLOC (char, strlen (pname) + 1);
  strcpy (pname_t, pname); 
  pname_t[strlen (pname) - 1] = 'x';
  fprintf (stderr, "%s\n", pname_t);
  fp_t = fopen (pname_t, "w");
  fp = fopen (pname, "r");
  while ((c = fgetc (fp)) != EOF)
    {
      fputc (c, fp_t);
    }
  fclose (fp);
  fclose (fp_t);

  fp_t = fopen (pname_t, "r");
  fp = fopen (pname, "w");
  while ((c = fgetc (fp_t)) != '$')
    {
      fputc (c, fp);
    }
  fclose (fp);
  fclose (fp_t);

  unlink (pname_t);
  XFREE (pname_t);
  return EXIT_SUCCESS;
}



int
psi_io_mean_theta (char *pname, parameter *theta, int sample_size, int w)
{
  int r = EXIT_SUCCESS;
  int i, j;
  double draw[6];
  int gen;
  theta->s = 0.0;
  theta->p = 0.0;
  theta->a = 0.0;
  theta->c = 0.0;
  theta->g = 0.0;
  theta->t = 0.0;
  FILE *fp = NULL;

  fp = fopen (pname, "r");
  
  if (fp == NULL)
    psi_fatal ("there should be parameter file to get theta star, this function should be called only when there is a complete set of estimation");

  for (i = 0; i < 4; i++) 
    {
      next_line (1, fp);
      for (j = 0; j < sample_size; j++) 
        {
          r = fscanf (fp, "%d%lf%lf%lf%lf%lf%lf\n", 
                      &gen, 
                      &draw[PSI_SAMPLE_A],
                      &draw[PSI_SAMPLE_C],
                      &draw[PSI_SAMPLE_G],
                      &draw[PSI_SAMPLE_T],
                      &draw[PSI_SAMPLE_S],
                      &draw[PSI_SAMPLE_P]);
          if (r == EOF || r != 7) 
            {
              return EXIT_FAILURE;
            }
          if (i == w)
            {
          theta->a += draw[PSI_SAMPLE_A];
          theta->c += draw[PSI_SAMPLE_C];
          theta->g += draw[PSI_SAMPLE_G];
          theta->t += draw[PSI_SAMPLE_T];
          theta->s += draw[PSI_SAMPLE_S];
          theta->p += draw[PSI_SAMPLE_P];
            }
        }
      if (i == w)
        {
          fclose (fp);
          theta->a /= sample_size;
          theta->c /= sample_size;
          theta->g /= sample_size;
          theta->t = 1.0 - theta->a - theta->c - theta->g;
          theta->s /= sample_size;
          theta->p /= sample_size;
          return EXIT_SUCCESS;
        }
    }
  return EXIT_SUCCESS;
}

void 
psi_io_write_mcmc_acceptance (FILE *fn,
                              double delta_s_accepted,
                              double delta_p_accepted,
                              double delta_s_rejected,
                              double delta_p_rejected,
                              int gmel_news_accept_s,
                              int gmel_news_total_s,
                              int gmel_news_accept_p,
                              int gmel_news_total_p,
                              int gmel_news_accept_n,
                              int gmel_news_total_n)
{
  fprintf (fn, "*** A C C E P T A N C E ***\n");
  fprintf (fn, "%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\t%d\n",
           delta_s_accepted,
           delta_p_accepted,
           delta_s_rejected,
           delta_p_rejected,
           gmel_news_accept_s,
           gmel_news_total_s,
           gmel_news_accept_p,
           gmel_news_total_p,
           gmel_news_accept_n,
           gmel_news_total_n);
}

void 
psi_io_read_mcmc_acceptance (char *pname,
                              double *delta_s_accepted,
                              double *delta_p_accepted,
                              double *delta_s_rejected,
                              double *delta_p_rejected,
                              int *gmel_news_accept_s,
                              int *gmel_news_total_s,
                              int *gmel_news_accept_p,
                              int *gmel_news_total_p,
                              int *gmel_news_accept_n,
                              int *gmel_news_total_n)
{
  /* locate the A C ... */
  /* read them */
  int r = EXIT_SUCCESS;
  char line[PSI_IO_GMEL_LINE_MAX];
  FILE *fp = NULL;
  char c;
  FILE *fp_t = NULL;
  char *pname_t = NULL; 

  fp = fopen (pname, "r");
  
  if (fp == NULL) 
    {
      *delta_s_accepted = 0.0;
      *delta_p_accepted = 0.0;
      *delta_s_rejected = 0.0;
      *delta_p_rejected = 0.0;
      *gmel_news_accept_s = 0;
      *gmel_news_total_s = 0;
      *gmel_news_accept_p = 0;
      *gmel_news_total_p = 0;
      *gmel_news_accept_n = 0;
      *gmel_news_total_n = 0;
      return;
    }

  while (fgets(line, PSI_IO_GMEL_LINE_MAX, fp) != NULL)
    {
      if (!strcmp (line, "*** A C C E P T A N C E ***\n"))
        {
           r = fscanf (fp, "%lf%lf%lf%lf%d%d%d%d%d%d\n", 
                       delta_s_accepted,
                       delta_p_accepted,
                       delta_s_rejected,
                       delta_p_rejected,
                       gmel_news_accept_s,
                       gmel_news_total_s,
                       gmel_news_accept_p,
                       gmel_news_total_p,
                       gmel_news_accept_n,
                       gmel_news_total_n);
          fclose (fp);
          if (r == EOF || r != 10) 
            {
              psi_fatal ("there should be burn state");
            }
          break;  
        }
    }
  if (r != 10)
    psi_fatal ("there must be ACCEPTANCE");

  /* remove them */
  pname_t = XMALLOC (char, strlen (pname) + 1);
  strcpy (pname_t, pname); 
  pname_t[strlen (pname) - 1] = 'x';
  fprintf (stderr, "%s\n", pname_t);
  fp_t = fopen (pname_t, "w");
  fp = fopen (pname, "r");
  while ((c = fgetc (fp)) != EOF)
    {
      fputc (c, fp_t);
    }
  fclose (fp);
  fclose (fp_t);

  fp_t = fopen (pname_t, "r");
  fp = fopen (pname, "w");
  while ((c = fgetc (fp_t)) != '*')
    {
      fputc (c, fp);
    }
  fclose (fp);
  fclose (fp_t);

  unlink (pname_t);
  XFREE (pname_t);

}
