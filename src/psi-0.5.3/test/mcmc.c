/* mcmc.c -- test program of the MCMC module
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
static char *pname = GMEL_STR_PARAM;
static char *aname = GMEL_STR_ALLDRAWS;
static char *oname = NULL;
static char *bname = GMEL_STR_BFACTOR;
static char *qname = GMEL_STR_Q;
static FILE *ofile;

static int unsetup_psi ();
static int setup_psi ();

/** @start 1 */
int
/* main (int argc, char * const argv[]) */
main (void)
{
  int result = EXIT_SUCCESS;
/*
  parameter theta_dummy;
*/
 
  if (oname == NULL) 
    ofile = stdout;
  else 
    ofile = fopen (oname, "w");

  setup_psi ();

/*
  int want_disable_s = 0;
  int want_disable_p = 0;
  result = execute_estimation (want_disable_s, want_disable_p, 0, 0, 
                               &theta_dummy);
*/
  if (result != EXIT_SUCCESS) 
    return EXIT_FAILURE;

  unsetup_psi ();

  if (ofile != stdout)
    fclose (ofile);

  exit (result);
}
/** @end 1 */

int 
unsetup_psi ()
{
  delete_gridpoints ();

  /* unsetup_bf (); we don't call any function of bf.c */
  unsetup_mcmc ();

  /* random number generator and energy data */
  finalize_drevol (PSI_ENERGY_JONES);
  fin_rng ();
  return EXIT_SUCCESS;
}

int 
setup_psi ()
{
  /* random number generator and energy data */
  init_rng (0);
  /* setup_energy_jones ("../data/pdb1g90.A.prd", "../data/energy.int"); */
  setup_energy_jones ("../data/pdbtest.2.dat", "../data/energy.int"); 
  initialize_drevol (PSI_ENERGY_JONES);

  /* gibbs */
  int gibbs_size = 2;
  int gibbs_burn = 1;
  int gibbs_freq = 1;
  setup_gibbs (gibbs_size, gibbs_freq, gibbs_burn);

  /* sample */
  int sample_size = 5;
  int sample_freq = 1;
  int burnin_ratio = 0;
  setup_mcmc_sample (sample_size, sample_freq, burnin_ratio);

  /* grid */
  double gridpoint_s_begin = GMEL_GRIDPOINT_S_BEGIN; /* -5.4  */
  double gridpoint_s_end = GMEL_GRIDPOINT_S_END;     /*  5.4  */
  int gridpoint_s_number = GMEL_GRIDPOINT_S_NUMBER;  /* 51    */
  double gridpoint_p_begin = GMEL_GRIDPOINT_P_BEGIN; /* -1.0  */
  double gridpoint_p_end = GMEL_GRIDPOINT_P_END;     /*  1.0  */
  int gridpoint_p_number = GMEL_GRIDPOINT_P_NUMBER;  /* 41    */
  double gridpoint_n_begin = GMEL_GRIDPOINT_N_BEGIN; /*  0.0  */
  double gridpoint_n_end = GMEL_GRIDPOINT_N_END;     /*  1.0  */
  int gridpoint_n_number = GMEL_GRIDPOINT_N_NUMBER;  /* 11    */
  setup_mcmc_grid (gridpoint_s_begin, gridpoint_s_end, gridpoint_s_number,
                   gridpoint_p_begin, gridpoint_p_end, gridpoint_p_number,
                   gridpoint_n_begin, gridpoint_n_end, gridpoint_n_number);

  /* prior */
  double flatprior_s_begin = GMEL_FLATPRIOR_S_BEGIN; /* -5.0  */
  double flatprior_s_end = GMEL_FLATPRIOR_S_END;     /*  5.0  */
  double flatprior_p_begin = GMEL_FLATPRIOR_P_BEGIN; /* -1.0  */
  double flatprior_p_end = GMEL_FLATPRIOR_P_END;     /*  1.0  */
  double flatprior_n_begin = GMEL_FLATPRIOR_N_BEGIN; /*  0.0  */
  double flatprior_n_end =  GMEL_FLATPRIOR_N_END;    /*  1.0  */
  setup_mcmc_prior (flatprior_s_begin, flatprior_s_end,
                    flatprior_p_begin, flatprior_p_end,
                    flatprior_n_begin, flatprior_n_end);

  /* delta */
  double delta_s = DELTA_W_S;
  double delta_p = DELTA_W_P;
  double delta_n = DELTA_PI;
  setup_mcmc_delta (delta_s, delta_p, delta_n);

  /* init */
  double init_s = GMEL_INIT_S;
  double init_p = GMEL_INIT_S;
  double init_a = GMEL_INIT_N;
  double init_c = GMEL_INIT_N;
  double init_g = GMEL_INIT_N;
  setup_mcmc_init (init_s, init_p, init_a, init_c, init_g);

  /* data */ 
  setup_mcmc_data ();

  /* file */
  setup_mcmc_file (ofile, pname, aname);


  /* want */
  int want_no_data = 0;
  int want_disable_s = 0;
  int want_disable_p = 0;
  int want_disable_n = 0;
  int want_save_alldraws = 0;
  int log_freq = GMEL_LOG_FREQ; /* print the number of iterations 
                                   every 10th iteration */
  int want_debug_est = 0;
  setup_mcmc_want (want_no_data, 
                   want_disable_s, want_disable_p, want_disable_n,
                   want_save_alldraws, want_debug_est, log_freq);

  /* these two functions should be gotten from mcmc.c module */

  /* want and file is okay but,
     sample_size and 0.9 what is this. 
     ofile should be replaced by oname and I will switch the positions
     of arguments: qname and aname. sample_size is a part of mcmc.c
     module. access_mcmc_sample_size () should be called in bf. c module.
   */
  int want_load_param = 0;
  setup_bf_want (want_load_param, want_save_alldraws);
  setup_bf_file (ofile, pname, aname, qname, bname);
  setup_bf_tol (0.9);

  /* creation of grid points */
  create_gridpoints(5);

  return EXIT_SUCCESS;
}


