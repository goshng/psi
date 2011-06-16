/* rng.c -- random number generator
   Copyright (C) 2000, 2006 Gary V. Vaughan
  
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
#include "gslwrap.h"
#include "rng.h"

static const int debug_detail = 0;

double rng_put_out (double d)
{
  return d;
}

static gsl_rng *gsl_r;

int 
init_rng (unsigned long int s)
{
  const gsl_rng_type *T;
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  gsl_rng_default_seed = s;
  gsl_r = gsl_rng_alloc (T);
  /* gsl_rng_set (gsl_r, s); */

  return EXIT_SUCCESS;
}

int 
save_rng (const char *rngname)
{
  FILE *stream = NULL;
  stream = fopen (rngname, "w");
  gsl_rng_fwrite (stream, gsl_r); 
  fclose (stream);

  return EXIT_SUCCESS;
}

int 
load_rng (const char *rngname)
{
  FILE *stream = NULL;
  stream = fopen (rngname, "r");
  if (stream == NULL) 
    return EXIT_FAILURE;

  const gsl_rng_type *T;
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  gsl_r = gsl_rng_alloc (T);

  gsl_rng_fread (stream, gsl_r); 
  fclose (stream);

  return EXIT_SUCCESS;
}

void
psi_rng_info (FILE *fp)
{
  assert (fp != NULL);  
  fprintf (fp, "generator type: %s\n", gsl_rng_name (gsl_r));
  fprintf (fp, "seed = %lu\n", gsl_rng_default_seed);
}

double 
psi_rng ()
{
  return gsl_rng_uniform (gsl_r);
}

int 
fin_rng ()
{
  gsl_rng_free (gsl_r);
  return EXIT_SUCCESS;
}

int 
choose_one_nucleotide (double* here) 
{
  double uniform_variable = gsl_rng_uniform (gsl_r);

  return choose_one_nucleotide_with_fixed (here, uniform_variable);
}

int
choose_one_nucleotide_with_fixed (double *here, double uniform_variable)
{
  int i;
  double r;
  double sum[4];

  sum[0] = *here; here++;
  for (i=1; i < 4; i++) 
    {
      sum[i] = logsum (*here, sum[i-1]);
      here++;
    }
  here -= i;

  r = log (uniform_variable) + sum[3];

  for (i = 0; i < 4; i++) 
    {
      if (r < sum[i]) 
        {
          if (debug_detail == 1) 
            {
              fprintf (stderr, "choose %d: u = %lf; r = %lf; %lf %lf %lf %lf; %lf, %lf, %lf, %lf\n\n", 
                      i, uniform_variable, exp(r), exp(here[0]), exp(here[1]), exp(here[2]), exp(here[3]), exp(sum[0]), exp(sum[1]), exp(sum[2]), exp(sum[3]));
            }  
          return i;
        }
    }

  psi_fatal ("we must have chosen one of four nucleotides");
  return EXIT_FAILURE;
}

int 
choose_one_sequence_by_rate (double *here, int len_dna)
{
  double uniform_variable = gsl_rng_uniform (gsl_r);
  
  return choose_one_element_with_fixed (here, len_dna * 3, uniform_variable);
}

int
choose_one_element_with_fixed (double *here, int n, double uniform_variable)
{
  int i;
  double r;
  double *sum = NULL;
  sum = XMALLOC (double, n);

  sum[0] = *here; here++;
  for (i=1; i < n; i++) 
    {
      sum[i] = logsum (*here, sum[i - 1]);
      here++;
    }
  here -= i;

  r = log (uniform_variable) + sum[n - 1];

  for (i = 0; i < n; i++) 
    {
      if (r < sum[i]) 
        {
          XFREE (sum);
          return i;
        }
    }

  psi_fatal ("we must have chosen one of 3*N sequences");
  return EXIT_FAILURE;
}


int
psi_ran_pi (const double pi[], double r_pi[])
{
  gsl_ran_dirichlet (gsl_r, 4, pi, r_pi);
  return EXIT_SUCCESS; 
}

double 
psi_ran_pi_lnpdf (const double alpha[], const double theta[])
{
  return gsl_ran_dirichlet_lnpdf (4, alpha, theta);
}

/* BUGGY */
double
psi_ran_flat (double a, double b)
{
  return gsl_ran_flat (gsl_r, a, b);
}

double 
psi_sample_ran_flat (double v, double delta, double min, double max)
{
  double ub, lb;
  lb = v - delta > min ? v - delta : min;
  ub = v + delta < max ? v + delta : max;
  return gsl_ran_flat (gsl_r, lb, ub);
}

double
psi_ran_flat_pdf (double x, double mu, double delta, double min, double max)
{
  double lb, ub;
  lb = mu - delta > min ? mu - delta : min;
  ub = mu + delta < max ? mu + delta : max;
  if (lb < x && x < ub)
    return 1 / (ub - lb);
  else
    return 0.0;
}

int
choose_a_site (int len_dna)
{
  return gsl_rng_uniform_int (gsl_r, (unsigned long int) len_dna);
}

int
psi_rng_dna_seq (int *dna, int n)
{
  int i;
  for (i = 0; i < n; i++)
    {
      dna[i] = gsl_rng_uniform_int (gsl_r, 4);
    }
  return EXIT_SUCCESS;
}
