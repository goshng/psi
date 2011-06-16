/* rng.c -- test program of the random number generator module
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

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <signal.h>

#include <psi/psi.h>

/* static int psi_init  (Sic *psi); */

double 
put_out (double a)
{
  return a;
}

/** @start 1 */
int
/* main (int argc, char * const argv[]) */
main (void)
{
  int i;
  int result = EXIT_SUCCESS;
  int r1, r2;
  double u;
  double log_array[3] = { 1, 2, 3 };

  init_rng (0);
  u = psi_rng ();
  if (gsl_fcmp (u, 0.999742, 1e-6) != 0)
    {
      return EXIT_FAILURE;  
    } 

  double here[4] = { log(1.0), log(1.0), log(2.0), log(3.0) };
  
  double here2[4] = { log(1.0/2.0), log(1.0/2.0), log(2.0/2.0), log(3.0/2.0) };
  r1 = choose_one_nucleotide_with_fixed (here, 0.1);
  r2 = choose_one_nucleotide_with_fixed (here2, 0.1);
  if (r1 != r2)
    {
      return EXIT_FAILURE;
    }
  r1 = choose_one_nucleotide_with_fixed (here, 0.2);
  r2 = choose_one_nucleotide_with_fixed (here2, 0.2);
  if (r1 != r2)
    {
      return EXIT_FAILURE;
    }
  r1 = choose_one_nucleotide_with_fixed (here, 0.5);
  r2 = choose_one_nucleotide_with_fixed (here2, 0.5);
  if (r1 != r2)
    {
      return EXIT_FAILURE;
    }

  fin_rng ();

  init_rng (0);
  for (i = 0; i < 10; i++) 
    {
      u = psi_rng ();
    }
  save_rng ("rng.txt");
  fin_rng ();
  
  load_rng ("rng.txt");
  u = psi_rng ();
  if (gsl_fcmp (u, 0.759944, 1e-6) != 0)
    {
      return EXIT_FAILURE;  
    } 
  fin_rng ();

  init_rng (10);
  u = psi_rng ();
  if (gsl_fcmp (u, 0.771321, 1e-6) != 0)
    {
      return EXIT_FAILURE;  
    } 
  fin_rng ();

  u = psi_gaussian_lnpdf (0.0, 1.0, 0.0) * 1000; 
  if (trunc(u) != -918) 
    {
      return EXIT_FAILURE;  
    } 

  u = psi_sum_exp (log_array, 3) * 1000;
  if (trunc(u) != 3407) 
    {
      return EXIT_FAILURE;  
    } 
  
  u = logsum (1, 1) * 1000;
  if (trunc(u) != 1693) 
    {
      return EXIT_FAILURE;  
    } 
  u = logsum (10, 1) * 1000;
  if (trunc(u) != 10000) 
    {
      return EXIT_FAILURE;  
    } 
  u = logsum (-10, -1) * 1000;
  if (trunc(u) != -999) 
    {
      return EXIT_FAILURE;  
    } 
  u = logsum (-100, -1) * 1000;
  if (trunc(u) != -1000) 
    {
      return EXIT_FAILURE;  
    } 
  u = logsum (-1000, -1) * 1000;
  if (trunc(u) != -1000) 
    {
      return EXIT_FAILURE;  
    } 
  
  double alpha[4] = {21.0, 19.0, 26.00, 34.00};
  double theta[4] = {0.21, 0.19, 0.26, 0.34};
  u = psi_pi_lnpdf (alpha, theta) * 1000;
  if (trunc(u) != 6961)
    {
      return EXIT_FAILURE;  
    } 

  
  double **sample = NULL;
  sample = XMALLOC (double *, 1);
  sample[0] = XMALLOC (double, 30);
  for (i = 0; i < 30; i++)
    {
      sample[0][i] = (double) i;
    }
  take_subserial_sample (sample, 1, 30, 7);
  if (sample[0][0] != 0.0 
      || sample[0][1] != 4.0 
      || sample[0][2] != 8.0
      || sample[0][3] != 12.0
      || sample[0][4] != 17.0
      || sample[0][5] != 21.0
      || sample[0][6] != 25.0)
    {
      return EXIT_FAILURE;  
    }
  XFREE (sample[0]);
  XFREE (sample);
  
/*
  sample = XMALLOC (double *, 2);
  sample[0] = XMALLOC (double, 10);
  sample[1] = XMALLOC (double, 10);
  sample[0][0] = -0.708370061;
  sample[0][1] =  0.337756654; 
  sample[0][2] =  1.580066975;
  sample[0][3] = -0.174330851;  
  sample[0][4] =  0.914981058;
  sample[0][5] = -0.009618598;
  sample[0][6] =  1.363124301;
  sample[0][7] = -0.509507790;
  sample[0][8] = -0.885103555;
  sample[0][9] =  0.602615608;
  sample[1][0] = 1.0;
  sample[1][1] = 2.0;
  sample[1][2] = 3.0;
  sample[1][3] = 4.0;
  sample[1][4] = 5.0;
  sample[1][5] = 6.0;
  sample[1][6] = 7.0;
  sample[1][7] = 8.0;
  sample[1][8] = 9.0;
  sample[1][9] = 10.0;
  u = get_longest_autocorr_time_sample (sample, 2, 10);
  printf ("max: %lf\n", u); 

  XFREE (sample[0]);
  XFREE (sample[1]);
  XFREE (sample);
*/

  /* read xaa and print out the autocorrelation time */
  double *draws = XMALLOC (double, 100); 
  double tau;
  gmel_read_draws ("../data/xaa", draws, PSI_SAMPLE_A, 100);
  tau = get_autocorrelation_time (draws, 100);
  if (gsl_fcmp (tau, 1.338333, 1e-6) != 0)
    {
      return EXIT_FAILURE;  
    } 

  gmel_read_draws ("../data/xaa", draws, PSI_SAMPLE_C, 100);
  tau = get_autocorrelation_time (draws, 100);
  if (gsl_fcmp (tau, 1.0, 1e-6) != 0)
    {
      return EXIT_FAILURE;  
    } 

  gmel_read_draws ("../data/xaa", draws, PSI_SAMPLE_G, 100);
  tau = get_autocorrelation_time (draws, 100);
  if (gsl_fcmp (tau, 1.0, 1e-6) != 0)
    {
      return EXIT_FAILURE;  
    } 

  gmel_read_draws ("../data/xaa", draws, PSI_SAMPLE_T, 100);
  tau = get_autocorrelation_time (draws, 100);
  if (gsl_fcmp (tau, 1.354983, 1e-6) != 0)
    {
      return EXIT_FAILURE;  
    } 

  gmel_read_draws ("../data/xaa", draws, PSI_SAMPLE_S, 100);
  tau = get_autocorrelation_time (draws, 100);
  if (gsl_fcmp (tau, 1.0, 1e-6) != 0)
    {
      return EXIT_FAILURE;  
    } 

  gmel_read_draws ("../data/xaa", draws, PSI_SAMPLE_P, 100);
  tau = get_autocorrelation_time (draws, 100);
  if (gsl_fcmp (tau, 1.151005, 1e-6) != 0)
    {
      return EXIT_FAILURE;  
    } 

  
/*
  u = psi_gaussian_lnpdf (0.0, 0.9, 0.0); 
  printf ("ln_normal: %lf\n", u);
  u = psi_gaussian_lnpdf (0.0, 1.0, 0.0); 
  printf ("ln_normal: %lf\n", u);
*/

  exit (result);
}
/** @end 1 */
