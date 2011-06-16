/* mc.c -- Multiple Chains
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
#include "seq.h"
#include "rng.h"
#include "energy.h"
#include "gslwrap.h"
#include "grid.h"
#include "gibbs.h"
#include "io.h"
#include "mcmc.h"
#include "bf.h"
#include "mc.h"

static int is_fixed_run = 0;
static int number_chains = 0;
static double **init_param;

int
unsetup_mc ()
{
  int i;
  assert (number_chains > 1); 
  for (i = 0; i < number_chains; i++)
    XFREE (init_param[i]);
  XFREE (init_param);
  
  return EXIT_SUCCESS;
}

int 
setup_mc_chains (int n)
{
  int i;
  assert (n > 1);
  number_chains = n;
  init_param = XMALLOC (double *, n);
  for (i = 0; i < n; i++)
    init_param[i] = XMALLOC (double, 6);
  return EXIT_SUCCESS;
}

int
setup_mc_option (int i)
{
  is_fixed_run = i;
  return EXIT_SUCCESS;
}

int 
setup_mc_add_ip (double a, double c, double g, double s, double p)
{
  assert (number_chains > 1);
  assert (a > 0.0);
  assert (c > 0.0);
  assert (g > 0.0);
  assert (a + c + g < 1.0);
  static int i = 0;
  if (i == number_chains)
    {
      i = 0;
    }
  if (i < -1)
    psi_fatal ("setup_mc_add_ip: error - i is negative");

  init_param[i][PSI_SAMPLE_A] = a;
  init_param[i][PSI_SAMPLE_C] = c;
  init_param[i][PSI_SAMPLE_G] = g;
  init_param[i][PSI_SAMPLE_T] = 1.0 - a - c - g;
  init_param[i][PSI_SAMPLE_S] = s;
  init_param[i][PSI_SAMPLE_P] = p;

  i++;
  return EXIT_SUCCESS;
}

/*
  requirements are the following:

  1. init_rng
  2. setup_energy_jones (datname, energyname);
  3. initialize_drevol (type_energy);

  4. setup_psi ();
  initialize the starting values of parameters
  setup_mcmc_init (init_s, init_p, init_a, init_c, init_g);
  setup_mcmc_sample (sample_size, sample_freq, burnin_ratio, log_freq);

  setup_mcmc_file (ofile, pname, aname);
  setup_mcmc_want (want_disable_s, want_disable_p, want_disable_n,
                   want_save_alldraws);

  These functions are not called by this but one that call this function
  setup_io (ofile, sample_size, want_verbose);
  setup_gibbs (gibbs_size, gibbs_burn, gibbs_freq);
  setup_mcmc_data (n_a, n_c, n_g, n_t, solv, pair);
  setup_mcmc_grid (gridpoint_s_begin, gridpoint_s_end, gridpoint_s_number,
  setup_mcmc_prior (flatprior_s_begin, flatprior_s_end,
  setup_mcmc_delta (delta_s, delta_p, delta_n, delta_n_mu);
     ---> done in function that call this function

*/
int 
execute_multiple_chains ()
{
/*
   --- sample frequency check ---
   1. run a single chain until it has a sample of size \f$N\f$,
   2. calculate the autocorrelation time (\f$\tau\f$) of the sample,
   3. set \f$N_e\f$ to be \f$N/\tau\f$
   4. take \f$N_e\f$ out of the sample
   5. set \f$N_d\f$ to be \f$N-N_e\f$
   6. if \f$N_d == 0\f$ then stop the algorithm,
      else set \f$N\f$ to be \f$N_d\f$ and go to the step 1

   --- Different Burnin and Sample Frequency are determined ---

   --- convergence check --- *** NOTE: different Gibbs sampler each chain ***
   Once we find the appropriate sample frequency, then we run more
   chains to check the convergence. If the psr, the estimate of 
   convergence statistic, is not close to 1.0001, then we rerun all 
   of the chains. We repeat this convergence check procedure until
   the psr estimate is close enought to 1.0001.

   --- output the sample --- 

   --- run three more repetitions above for Bayes factor ---
 */

  if (is_fixed_run == 0)
    find_mcmc_option ();

  /*
   NOTE: different Gibbs sampler each chain ***
   */
  if (number_chains > 1)
    {
      run_mcmc_multichain (number_chains, init_param, is_fixed_run);
    }
  else
    {
      psi_warning ("there should be more than one chain and");
      psi_warning ("we skip the multichain procedure");
      setup_mcmc_init (0.0, 0.0, 0.25, 0.25, 0.25);
    }

  /* more chains should be run in different mcmc chains */
  /* setup_mcmc_init (0.0, 0.0, 0.25, 0.25, 0.25); */
  /* setup_mcmc_init in run_mcmc_multichain */
  execute_bayesfactor ();

  /* running on the .p output file */
  /* HERE WE ARE */

  return EXIT_SUCCESS;
}




