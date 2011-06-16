/* mcmc.c -- MCMC procedure of the single sequence analysis
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

/*
   s and p
   prior      proposal
 *  Flat       Flat      : okay
   Normal     Normal    : okay
   Normal     Flat      : not implemented
   Flat       Normal    : not implemented

   pi
   prior      proposal
   Flat       Flat      : okay
 *  Flat       Dirichlet : okay 
   Dirichlet  Flat      : not implemented
   Dirichlet  Dirichlet : not implemented

   We will fix the prior and proposal like this:
   for s and p, flat prior and flat proposal
   for pi, flat prior and Dirichlet proposal
*/

/*!\file mcmc.c
   \author Sang Chul Choi
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
#include "grid.h"
#include "gslwrap.h"
#include "gibbs.h"
#include "io.h"
#include "mcmc.h"
#include "nsv.h"

/* state */
static int want_load_state = 0;
static int limit_time = 0;
static time_t start_time;

static const int want_verbose = 1;
static const int want_log = 0;
static const int want_progress_log = 1;

static const int want_use_Dirichlet_proposal_n = 1;
/* debug */
static const int want_debug_gibbs = 0;
static const int want_debug_input = 0;
static const int want_debug_multi = 1;
static const int ERR_TOOMANYZEROS = 1;

/* data */
static int data_num_A;
static int data_num_C;
static int data_num_G;
static int data_num_T;
static double data_minus2E_solv;
static double data_minus2E_pair;

/* prior */
static double flatprior_s_begin;
static double flatprior_s_end;
static double log_flatprior_s;
static double prior_normal_sd_s;
static double flatprior_p_begin;
static double flatprior_p_end;
static double log_flatprior_p;
static double prior_normal_sd_p;
static double flatprior_n_begin;
static double flatprior_n_end;
  

/* grid */
static double gridpoint_s_begin;
static double gridpoint_s_end;
static int gridpoint_s_number;
static double gridpoint_p_begin;
static double gridpoint_p_end;
static int gridpoint_p_number;
static double gridpoint_n_begin;
static double gridpoint_n_end;
static int gridpoint_n_number;
static int gibbs_size = 0;

/* delta */
static double delta_s;
static double delta_s_sd;
static double delta_p;
static double delta_p_sd;
static double delta_n;
static double delta_n_mu;

/* init */
static double init_s;
static double init_p;
static double init_a;
static double init_c;
static double init_g;

/* file */
static char *pname = NULL;
static char *aname = NULL;
static FILE *ofile = NULL;

/* disable */
static int want_disable_s = 0;
static int want_disable_p = 0;
static int want_disable_n = 0;
static int want_save_alldraws = 1;
static int want_debug_est = 0;
static int log_freq = 10;
static int want_no_data = 0;

/* sample */
static int sample_size = 100;
static int sample_freq = 1;
static int sample_burn = 0;

static double limit_tau = 1.1;
static double limit_psr = 1.01;

/* static functions */
static void psi_mcmc_find_means (double ***posteriors, int number_chains);
static void copy_draws (double **draws, double ***posteriors, int k, int n_c);
static int check_time ();

#define ENUMERATION_GRIDPOINTS_M2_BEGIN         \
  for (i = 0; i < gridpoint_s_number; i++) { \
  for (j = 0; j < gridpoint_p_number; j++) { \
  for (k = 0; k < gridpoint_n_number; k++) { \
  for (l = 0; l < gridpoint_n_number; l++) { \
  for (m = 0; m < gridpoint_n_number; m++) {    
#define ENUMERATION_GRIDPOINTS_M2_END           \
  } \
  } \
  } \
  } \
  }

static double delta_s_accepted = 0;
static double delta_p_accepted = 0;
static double delta_s_rejected = 0;
static double delta_p_rejected = 0;
static int gmel_news_accept_s = 0;
static int gmel_news_total_s = 0;
static int gmel_news_accept_p = 0;
static int gmel_news_total_p = 0;
static int gmel_news_accept_n = 0;
static int gmel_news_total_n = 0;

void 
write_mcmc_info ()
{
  double t = 1.0 - init_a - init_c - init_g;
  fprintf (ofile, "sample:\t%d\t%d\t%d\n", 
           sample_size, sample_burn, sample_freq);
  fprintf (ofile, "delta:\t%lf\t%lf\t%lf\n", 
           delta_s, delta_p, delta_n);
  fprintf (ofile, "init:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
           init_s, init_p, init_a, init_c, init_g, t);
  fprintf (ofile, "grid:\t%lf\t%lf\t%d\t%lf\t%lf\t%d\t%lf\t%lf\t%d\n", 
           gridpoint_s_begin, gridpoint_s_end, gridpoint_s_number,
           gridpoint_p_begin, gridpoint_p_end, gridpoint_p_number,
           gridpoint_n_begin, gridpoint_n_end, gridpoint_n_number);
}

void
initialize_mcmc_acceptance ()
{
  delta_s_accepted = 0;
  delta_p_accepted = 0;
  delta_s_rejected = 0;
  delta_p_rejected = 0;
  gmel_news_accept_s = 0;
  gmel_news_total_s = 0;
  gmel_news_accept_p = 0;
  gmel_news_total_p = 0;
  gmel_news_accept_n = 0;
  gmel_news_total_n = 0;
}



void
read_mcmc_acceptance (char *pname)
{
  psi_io_read_mcmc_acceptance (pname,
                               &delta_s_accepted,
                               &delta_p_accepted,
                               &delta_s_rejected,
                               &delta_p_rejected,
                               &gmel_news_accept_s,
                               &gmel_news_total_s,
                               &gmel_news_accept_p,
                               &gmel_news_total_p,
                               &gmel_news_accept_n,
                               &gmel_news_total_n);
}

int
access_mcmc_sample_size ()
{
  if (sample_size == 0)
    psi_fatal ("sample size should be set before calling this access");
  return sample_size;
}

int
access_mcmc_sample_freq ()
{
  return sample_freq;
}

int
access_mcmc_sample_burn ()
{
  return sample_burn;
}

int 
unsetup_mcmc ()
{
  XFREE (pname);
  XFREE (aname);
  return EXIT_SUCCESS;
}

int 
setup_mcmc_data ()
{
  data_members_jones_measure (&data_num_A,
                              &data_num_C,
                              &data_num_G,
                              &data_num_T,
                              &data_minus2E_solv,
                              &data_minus2E_pair);
  return EXIT_SUCCESS;
}

int 
setup_mcmc_prior (double s_b, double s_e, double p_b, double p_e, 
                  double n_b, double n_e)
{
  assert (s_b < s_e);
  assert (p_b < p_e);
  assert (n_b < n_e);
  flatprior_s_begin = s_b;
  flatprior_s_end = s_e;
  flatprior_p_begin = p_b;
  flatprior_p_end = p_e;
  flatprior_n_begin = n_b;
  flatprior_n_end = n_e;

  log_flatprior_s = log (1 / (s_e - s_b));
  log_flatprior_p = log (1 / (p_e - p_b));
  prior_normal_sd_s = (s_e - s_b) * 3.4641; /* just the sd of uniform */
  prior_normal_sd_p = (p_e - p_b) * 3.4641; /* just the sd of uniform */
  return EXIT_SUCCESS;
}

int 
setup_mcmc_grid (double s_b, double s_e, int s_n, 
                 double p_b, double p_e, int p_n,
                 double n_b, double n_e, int n_n)
{
  gridpoint_s_begin = s_b;
  gridpoint_s_end = s_e;
  gridpoint_s_number = s_n;
  gridpoint_p_begin = p_b;
  gridpoint_p_end = p_e;
  gridpoint_p_number = p_n;
  gridpoint_n_begin = n_b;
  gridpoint_n_end = n_e;
  gridpoint_n_number = n_n;
  gibbs_size = access_gibbs_size();
  return EXIT_SUCCESS;
}

int 
setup_mcmc_delta (double s, double p, double n)
{
  delta_s = s;
  delta_p = p;
  delta_n = n;
  delta_s_sd = s / 1.732051;
  delta_p_sd = p / 1.732051;
  delta_n_mu = 1 / (n * n);
  return EXIT_SUCCESS;
}

int
setup_mcmc_init (double s, double p, double a, double c, double g)
{
  init_s = s;
  init_p = p;
  init_a = a;
  init_c = c;
  init_g = g;
  return EXIT_SUCCESS;
}

int
setup_mcmc_file (FILE *fp, char *pn, char *an)
{
  assert (pn != NULL);
  assert (an != NULL);
  pname = xstrdup (pn);
  aname = xstrdup (an);
  ofile = fp;

  return EXIT_SUCCESS;
}

int 
setup_mcmc_want (int no, int s, int p, int n, int d, int debug, int lf)
{
  assert (no == 0 || no == 1);
  assert (s == 0 || s == 1);
  assert (p == 0 || p == 1);
  assert (n == 0 || n == 1);
  assert (d == 0 || d == 1);
  assert (debug == 0 || debug == 1);
  assert (lf > 0);

  want_no_data = no;
  want_disable_s = s;
  want_disable_p = p;
  want_disable_n = n;
  want_save_alldraws = d;
  want_debug_est = debug;
  log_freq = lf;
  return EXIT_SUCCESS;
}

int
setup_mcmc_sample (int s, int f, int b)
{
  assert (s > 0);
  assert (f > 0);
  assert (b >= 0); 

  sample_size = s;
  sample_freq = f;
  sample_burn = b;
  return EXIT_SUCCESS;
}

int
setup_mcmc_limit (double t, double p)
{
  assert (t > 1.0);
  assert (p > 1.0);

  limit_tau = t;
  limit_psr = p;
  return EXIT_SUCCESS;
}

int 
setup_mcmc_state (int s, int l)
{
  want_load_state = s;
  limit_time = l;
  start_time = time (NULL);
  return EXIT_SUCCESS;
}


int
new_pi (parameter *p_theta)
{
  int r = EXIT_SUCCESS;

  double log_random;
  double lnR, random, jumping_term, prior_term, likelihood_term;
  parameter theta_prime;

  /* Dirichlet proposal */
  r = get_theta_prime (p_theta, &theta_prime, PSI_PI_D_MODE);  

  /* Likelihood Term */
  if (want_no_data == 0)
    get_likelihood (&likelihood_term, p_theta, &theta_prime);
  else
    likelihood_term = 0.0;  

  /* Prior Term */
  r = get_prior (&prior_term, *p_theta, theta_prime, PSI_PI_D_MODE);
 
  /* Jumping Term */
  /* Dirichlet proposal */
  r = get_jumping (&jumping_term, *p_theta, theta_prime, PSI_PI_D_MODE);

  theta_prime.w = PSI_PI_D_MODE;

  if (!finite(likelihood_term) && !finite(prior_term) 
      && !finite(jumping_term)) 
    {
      /* debug_write_parameter (*p_theta);
      debug_write_parameter (theta_prime); */
      psi_fatal ("not finite likelihood, prior, or jumping");
    }

  lnR = likelihood_term + prior_term + jumping_term;
  if (want_debug_est == 1) 
      fprintf (stderr, "New PI: L. %.3lf, P. %.3lf, J. %.3lf\t", 
               likelihood_term, prior_term, jumping_term);
 
  random = psi_rng ();
  log_random = log (random);

  gmel_news_total_n++;
  if (log_random < lnR) 
    {
      if (want_debug_est == 1) 
        {
          fprintf (stderr, "Accept N: %.3lf < %.3lf\n", log_random, lnR);
        }
      *p_theta = theta_prime;
      gmel_news_accept_n++;
    } 
  else 
    {
      if (want_debug_est == 1) 
        {
          fprintf(stderr, "Reject N: %.3lf < %.3lf\n", log_random, lnR);
        }
    }
 
  return EXIT_SUCCESS;
}

/*!\brief Sample \f$s\f$
 
   This function samples \f$s\f$ from posterior distribution. 

   \param curr current parameter values
   \param energy energy information
 */
int
new_s (parameter *p_theta)
{
  int r = EXIT_SUCCESS;
  double log_random;
  double lnR, random, jumping_term, prior_term, likelihood_term;
  parameter theta_prime;

  /* flat proposal */
  r = get_theta_prime( p_theta, &theta_prime, PSI_S_MODE );

  /* Likelihood Term */
  if (want_no_data == 0)
    get_likelihood (&likelihood_term, p_theta, &theta_prime);
  else
    likelihood_term = 0.0;

  /* Prior Term: flat prior */
  r = get_prior (&prior_term, *p_theta, theta_prime, PSI_S_MODE);

  /* Jumping Term: flat proposal */
  r = get_jumping (&jumping_term, *p_theta, theta_prime, PSI_S_MODE);

  theta_prime.w = p_theta->w;

  if (!finite(likelihood_term) && !finite(prior_term) 
      && !finite(jumping_term)) {
     /* debug_write_parameter (*p_theta);
     debug_write_parameter (theta_prime); */
     psi_fatal ("new_s: not finite likelihood, prior, or jumping");
  }

  lnR = likelihood_term + prior_term + jumping_term;
  if (want_debug_est == 1)
     fprintf (stderr, "New S: %.3lf -> %.3lf, L. %.3lf, P. %.3lf, J. %.3lf\t", 
              p_theta->s, theta_prime.s, likelihood_term, prior_term, jumping_term);

  random = psi_rng ();
  errno = 0;
  log_random = log (random);

  gmel_news_total_s++;
  if (log_random < lnR) 
    {
      if (want_debug_est == 1) 
        {
          fprintf (stderr, "Accept S: %.3lf < %.3lf\n", log_random, lnR);
        }
      delta_s_accepted += fabs (p_theta->s - theta_prime.s);
      *p_theta = theta_prime;
      gmel_news_accept_s++;
    } 
  else 
    {
      delta_s_rejected += fabs (p_theta->s - theta_prime.s);
      if (want_debug_est == 1) 
        {
          fprintf(stderr, "Reject S: %.3lf < %.3lf\n", log_random, lnR);
        }
    }

  return EXIT_SUCCESS;
}

/*!\brief Sample \f$p\f$
 
   This function samples \f$p\f$ from posterior distribution. 

   \param curr current parameter values
   \param energy energy information
 */
int
new_p (parameter *p_theta)
{
  int r = EXIT_SUCCESS;
  double log_random;
  double lnR, random, jumping_term, prior_term, likelihood_term;
  parameter theta_prime;

  /* proposal from flat */
  r = get_theta_prime( p_theta, &theta_prime, PSI_P_MODE );

  /* Likelihood Term */
  if (want_no_data == 0)
    get_likelihood (&likelihood_term, p_theta, &theta_prime);
  else
    likelihood_term = 0.0;
                      
  /* Prior Term: flat */
  r = get_prior( &prior_term, *p_theta, theta_prime, PSI_P_MODE );

  /* Jumping Term: flat */
  r = get_jumping( &jumping_term, *p_theta, theta_prime, PSI_P_MODE );

  theta_prime.w = p_theta->w;
  
  if (!finite(likelihood_term) && !finite(prior_term) 
      && !finite(jumping_term)) 
    {
      /* debug_write_parameter (*p_theta);
      debug_write_parameter (theta_prime); */
      psi_fatal ("new_p: not finite likelihood, prior, or jumping");
    }

  lnR = likelihood_term + prior_term + jumping_term;
  if (want_debug_est == 1)
     fprintf (stderr, "New P: %.3lf -> %.3lf, L. %.3lf, P. %.3lf, J. %.3lf\t", 
              p_theta->p, theta_prime.p, likelihood_term, prior_term, jumping_term);
 
  random = psi_rng ();
  errno = 0;
  log_random = log (random); /* BUGGY: what if random is close to zero */

  gmel_news_total_p++;
  if (log_random < lnR) 
    {
      if (want_debug_est == 1) 
        {
          fprintf(stderr, "Accept P: %.3lf < %.3lf\n", log_random, lnR);
        }
      delta_p_accepted += fabs (p_theta->p - theta_prime.p);
      *p_theta = theta_prime;
      gmel_news_accept_p++;
    } 
  else 
    {
      delta_p_rejected += fabs (p_theta->p - theta_prime.p);
      if (want_debug_est == 1) 
        {
          fprintf(stderr, "Reject P: %.3lf < %.3lf\n", log_random, lnR);
        }
    }

  return EXIT_SUCCESS;
}

int
gmel_bf_s_sample(double s1, double *s2) {
  *s2 = psi_sample_ran_flat (s1, delta_s, flatprior_s_begin, flatprior_s_end);
  assert (isfinite(*s2));
  return EXIT_SUCCESS;
}

int
gmel_bf_p_sample(double p1, double *p2) {
  *p2 = psi_sample_ran_flat (p1, delta_p, flatprior_p_begin, flatprior_p_end);
  assert (isfinite(*p2));
  return EXIT_SUCCESS;
}

int
gmel_bf_n_sample (double a1, double c1, double g1, 
                  double *a2, double *c2, double *g2, int *w) 
{
  *w = 0; /* just for removing warning */
  int r = EXIT_SUCCESS;
  double t1 = 1.0 - a1 - c1 - g1;
  double t2_value;
  double *t2 = &t2_value;
  static int which_nuc = -1; which_nuc++;

  r = psi_err_base (a1, c1, g1, t1);
  if (r != EXIT_SUCCESS) {
     psi_fatal ("%s:%s:%d:%s; %s\n",
                program_name,
                __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                "invalid base frequency");
  }
  if (want_use_Dirichlet_proposal_n == 1) {
     double pi[4] = {a1 * delta_n_mu, 
                     c1 * delta_n_mu, 
                     g1 * delta_n_mu,
                     t1 * delta_n_mu};
     double r_pi[4] = {0.25, 0.25, 0.25, 0.25};
     psi_ran_pi (pi, r_pi);
     /* gsl_ran_dirichlet(gmel_gsl_r, 4, pi, r_pi); */
     while ( r_pi[PSI_DNA_A] < GMEL_N_MIN || r_pi[PSI_DNA_A] > GMEL_N_MAX
             || r_pi[PSI_DNA_C] < GMEL_N_MIN || r_pi[PSI_DNA_C] > GMEL_N_MAX
             || r_pi[PSI_DNA_G] < GMEL_N_MIN || r_pi[PSI_DNA_G] > GMEL_N_MAX
             || r_pi[PSI_DNA_T] < GMEL_N_MIN || r_pi[PSI_DNA_T] > GMEL_N_MAX ) {
        psi_ran_pi (pi, r_pi);
        /* gsl_ran_dirichlet(gmel_gsl_r, 4, pi, r_pi); */
     }
     r = psi_err_base (r_pi[PSI_DNA_A], r_pi[PSI_DNA_C], r_pi[PSI_DNA_G], r_pi[PSI_DNA_T]);
     if (r != EXIT_SUCCESS) {
        psi_fatal ("%s:%s:%d:%s; %s\n",
                   program_name,
                   __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                   "invalid base frequency");
     }
     *a2 = r_pi[PSI_DNA_A]; 
     *c2 = r_pi[PSI_DNA_C]; 
     *g2 = r_pi[PSI_DNA_G]; 
     *t2 = r_pi[PSI_DNA_T]; 
  } else {
     assert (0);
     psi_fatal ("no implemenation: 'simple' flat prior for nucleotide");
  }
  assert (isfinite(*a2));
  assert (isfinite(*c2));
  assert (isfinite(*g2));
  assert (isfinite(*t2));
  return r;
}

/* PSI_D_MODE, PSI_S_MODE, PSI_P_MODE only */
int 
get_theta_prime (parameter *p_theta, 
                 parameter *p_theta_prime, 
                 int mode ) 
{
  double minValue = 0.0;
  double maxValue = 0.0; 
  double delta = 0.0;
  double param = 0.0;
  double param_prime = 0.0;
  double L_theta, H_theta;

  *p_theta_prime = *p_theta;


  if ( mode == PSI_PI_A_MODE ) {
     param = p_theta->a;
     minValue = flatprior_n_begin;
     maxValue = flatprior_n_end;
     delta = delta_n;
  } else if ( mode == PSI_PI_C_MODE ) {
     param = p_theta->c;
     minValue = flatprior_n_begin;
     maxValue = flatprior_n_end;
     delta = delta_n;
  } else if ( mode == PSI_PI_G_MODE ) {
     param = p_theta->g;
     minValue = flatprior_n_begin;
     maxValue = flatprior_n_end;
     delta = delta_n;
  } else if ( mode == PSI_PI_T_MODE ) {
     param = p_theta->t;
     minValue = flatprior_n_begin;
     maxValue = flatprior_n_end;
     delta = delta_n;
  } else if ( mode == PSI_PI_D_MODE ) { /* Dirichlet Proposal */
     minValue = flatprior_n_begin;
     maxValue = flatprior_n_end;
     delta = delta_n;
  } else if ( mode == PSI_S_MODE ) {
     param = p_theta->s;
     minValue = flatprior_s_begin;
     maxValue = flatprior_s_end;
     delta = delta_s;
  } else if ( mode == PSI_P_MODE ) {
     param = p_theta->p;
     minValue = flatprior_p_begin;
     maxValue = flatprior_p_end;
     delta = delta_p;
  } else if ( mode == PSI_S_MODE_NORMAL ) {/* Normal Prior doesn't care boundary */
     param = p_theta->s;
     delta = delta_s;
  } else if ( mode == PSI_P_MODE_NORMAL ) {/* Normal Prior doesn't care boundary */
     param = p_theta->p;
     delta = delta_p;
  } else {
     psi_fatal ("get_theta_prime: invalid mode");
  }

  if (mode == PSI_PI_D_MODE) {
     double pi[4] = {p_theta->a * delta_n_mu, 
                     p_theta->c * delta_n_mu, 
                     p_theta->g * delta_n_mu,
                     p_theta->t * delta_n_mu};
     double r_pi[4] = {0.25, 0.25, 0.25, 0.25};
     /* get_sensible_pi (pi); */
     /* gsl_ran_dirichlet(gmel_gsl_r, 4, pi, r_pi); */
     psi_ran_pi (pi, r_pi);

     while (r_pi[PSI_DNA_A] < GMEL_N_MIN || r_pi[PSI_DNA_A] > GMEL_N_MAX
            || r_pi[PSI_DNA_C] < GMEL_N_MIN || r_pi[PSI_DNA_C] > GMEL_N_MAX
            || r_pi[PSI_DNA_G] < GMEL_N_MIN || r_pi[PSI_DNA_G] > GMEL_N_MAX
            || r_pi[PSI_DNA_T] < GMEL_N_MIN || r_pi[PSI_DNA_T] > GMEL_N_MAX) {
        /* gsl_ran_dirichlet (gmel_gsl_r, 4, pi, r_pi); */
        psi_ran_pi (pi, r_pi);
     }
     p_theta_prime->a = r_pi[PSI_DNA_A];
     p_theta_prime->c = r_pi[PSI_DNA_C];
     p_theta_prime->g = r_pi[PSI_DNA_G];
     p_theta_prime->t = r_pi[PSI_DNA_T];
  } else if (mode == PSI_PI_A_MODE || mode == PSI_PI_C_MODE 
             || mode == PSI_PI_G_MODE || mode == PSI_PI_T_MODE) { 
     assert (0);
  } else if (mode == PSI_S_MODE || mode == PSI_P_MODE) {
/*
     H_theta = (delta < (maxValue - param)) ? delta : (maxValue - param);
     L_theta = (delta < (param - minValue)) ? delta : (param - minValue);
     param_prime = param - L_theta + psi_rng () * (H_theta + L_theta);
*/
     L_theta = (param - delta > minValue) ? param - delta : minValue;
     H_theta = (param + delta < maxValue) ? param + delta : maxValue;
     param_prime = psi_ran_flat (L_theta, H_theta);
     if (param_prime >= maxValue)
       psi_fatal ("param: %lf, param_prime: %lf, delta: %lf, L_theta: %lf, H_theta: %lf, min: %lf, max: %lf\n",
                  param, param_prime, delta, L_theta, H_theta, minValue, maxValue);
     assert (param_prime > minValue);
     assert (param_prime < maxValue);
  } else if (mode == PSI_S_MODE_NORMAL || mode == PSI_P_MODE_NORMAL) {
     assert (0);
     psi_fatal ("not implemented: normal proposal");
/*
     - sqrt((2*delta)*(2*delta)/12), sqrt(3) is 1.732051 -
     double normal_jump_sd = delta / 1.732051; 
     param_prime = param + gsl_ran_gaussian (gmel_gsl_r, normal_jump_sd);
*/
  } else {
     assert(0);
     psi_fatal ("invalid mode");
  } 

  switch (mode) {
  case PSI_PI_A_MODE:
     p_theta_prime->a = param_prime;
     p_theta_prime->c = p_theta->c * (1 - param_prime) / (1 - param);
     p_theta_prime->g = p_theta->g * (1 - param_prime) / (1 - param);
     p_theta_prime->t = 1 - p_theta_prime->a - p_theta_prime->c - p_theta_prime->g;
     break;
  case PSI_PI_C_MODE:
     p_theta_prime->c = param_prime;
     p_theta_prime->g = p_theta->g * (1 - param_prime) / (1 - param);
     p_theta_prime->t = p_theta->t * (1 - param_prime) / (1 - param);
     p_theta_prime->a = 1 - p_theta_prime->c - p_theta_prime->g - p_theta_prime->t;
     break;
  case PSI_PI_G_MODE:
     p_theta_prime->g = param_prime;
     p_theta_prime->t = p_theta->t * (1 - param_prime) / (1 - param);
     p_theta_prime->a = p_theta->a * (1 - param_prime) / (1 - param);
     p_theta_prime->c = 1 - p_theta_prime->g - p_theta_prime->t - p_theta_prime->a;
     break;
  case PSI_PI_T_MODE:
     p_theta_prime->t = param_prime;
     p_theta_prime->a = p_theta->a * (1 - param_prime) / (1 - param);
     p_theta_prime->c = p_theta->c * (1 - param_prime) / (1 - param);
     p_theta_prime->g = 1 - p_theta_prime->t - p_theta_prime->a - p_theta_prime->c;
     break;
  case PSI_PI_D_MODE:
     /* Done when sampling, no need any Jacobian adjustment */
     break;
  case PSI_S_MODE:
  case PSI_S_MODE_NORMAL:
     p_theta_prime->s = param_prime;
     break;
  case PSI_P_MODE:
  case PSI_P_MODE_NORMAL:
     p_theta_prime->p = param_prime;
     break;
  default:
     assert( 0 );
  }

  return EXIT_SUCCESS;
}

int 
get_jumping (double *p_term, 
             parameter theta, parameter theta_prime,
             int mode ) 
{
  double minValue = 0;
  double maxValue = 0;
  double delta = 0;
  double param = 0;
  double param_prime = 0;
  double L_theta, H_theta;

  switch (mode) {
  case PSI_PI_A_MODE:
     param = theta.a;
     param_prime = theta_prime.a;
     minValue = flatprior_n_begin;
     maxValue = flatprior_n_end;
     delta = delta_n;
     break;
  case PSI_PI_C_MODE:
     param = theta.c;
     param_prime = theta_prime.c;
     minValue = flatprior_n_begin;
     maxValue = flatprior_n_end;
     delta = delta_n;
     break;
  case PSI_PI_G_MODE:
     param = theta.g;
     param_prime = theta_prime.g;
     minValue = flatprior_n_begin;
     maxValue = flatprior_n_end;
     delta = delta_n;
     break;
  case PSI_PI_T_MODE:
     param = theta.t;
     param_prime = theta_prime.t;
     minValue = flatprior_n_begin;
     maxValue = flatprior_n_end;
     delta = delta_n;
     break;
  case PSI_PI_D_MODE:
     break;
  case PSI_S_MODE:
     param = theta.s;
     param_prime = theta_prime.s;
     minValue = flatprior_s_begin;
     maxValue = flatprior_s_end;
     delta = delta_s;
     break;
  case PSI_P_MODE:
     param = theta.p;
     param_prime = theta_prime.p;
     minValue = flatprior_p_begin;
     maxValue = flatprior_p_end;
     delta = delta_p;
     break;
  case PSI_S_MODE_NORMAL:
  case PSI_P_MODE_NORMAL:
     assert (0);
     /* Normal Jump Does Not Affect Jumping Ratio */
     /* THAT IS WRONG! */
     /* This is true only if the proposal and prior are both
        normal density */
     break;
  default:
     assert (0);
  }

  if (mode == PSI_PI_D_MODE) {
     double p_theta[4] = {theta.a, 
                          theta.c,
                          theta.g,
                          theta.t};
     double p_alpha[4] = {theta.a * delta_n_mu, 
                          theta.c * delta_n_mu,
                          theta.g * delta_n_mu,
                          theta.t * delta_n_mu};
     double p_star_theta[4] = {theta_prime.a, 
                               theta_prime.c,
                               theta_prime.g,
                               theta_prime.t};
     double p_star_alpha[4] = {theta_prime.a * delta_n_mu, 
                               theta_prime.c * delta_n_mu,
                               theta_prime.g * delta_n_mu,
                               theta_prime.t * delta_n_mu};

     /* if (theta.a != 0 || theta.c != 0 || theta.g != 0 || theta.t != 0) { */
     assert(theta.a != 0);
     assert(theta.c != 0);
     assert(theta.g != 0);
     assert(theta.t != 0);
     assert(theta_prime.a != 0);
     assert(theta_prime.c != 0);
     assert(theta_prime.g != 0);
     assert(theta_prime.t != 0);

     /* get_sensible_pi (p_alpha); */
     /* get_sensible_pi (p_star_alpha); */
     *p_term = psi_ran_pi_lnpdf (p_star_alpha, p_theta)
               - psi_ran_pi_lnpdf (p_alpha, p_star_theta);
  } else if (mode == PSI_PI_A_MODE || mode == PSI_PI_C_MODE 
             || mode == PSI_PI_G_MODE || mode == PSI_PI_T_MODE) {
     assert (0);
  } else if (mode == PSI_S_MODE || mode == PSI_P_MODE) {
/*
     H_theta = (delta < (maxValue - param)) ? delta : (maxValue - param);
     L_theta = (delta < (param - minValue)) ? delta : (param - minValue);
*/
     H_theta = (param + delta < maxValue) ? param + delta : maxValue;
     L_theta = (param - delta > minValue) ? param - delta : minValue;
     *p_term = log(H_theta - L_theta); /* metro_bot */

/*
     if ( mode == PSI_PI_A_MODE || mode == PSI_PI_C_MODE 
          || mode == PSI_PI_G_MODE || mode == PSI_PI_T_MODE ) { 
        *p_term -= log(1.0 - param);
        *p_term -= log(1.0 - param);
     }
     H_theta = delta < (maxValue - param_prime) ? delta : (maxValue - param_prime);
     L_theta = delta < (param_prime - minValue) ? delta : (param_prime - minValue);
*/
     H_theta = (param_prime + delta < maxValue) 
               ? param_prime + delta : maxValue;
     L_theta = (param_prime - delta > minValue) 
               ? param_prime - delta : minValue;
     *p_term -= log(H_theta - L_theta); /* metro_top */
/*
     if ( mode == PSI_PI_A_MODE || mode == PSI_PI_C_MODE 
          || mode == PSI_PI_G_MODE || mode == PSI_PI_T_MODE ) { 
        *p_term += log(1.0 - param_prime);
        *p_term += log(1.0 - param_prime);
     }
*/
  } else if (mode == PSI_S_MODE_NORMAL || mode == PSI_P_MODE_NORMAL) {
     assert (0);
     /* Normal Jump Does Not Affect Jumping Ratio */
     /* This is true only if the proposal and prior are both
        normal density */
     *p_term = 0;
  } else {
     assert (0);
  }
  return EXIT_SUCCESS;
}

int 
get_prior (double *p_term, 
           parameter theta, parameter theta_prime,
           int mode ) 
{
  /* theta and theta_prime warn but ignore them! */

  switch (mode) {
  case PSI_PI_A_MODE:
  case PSI_PI_C_MODE:
  case PSI_PI_G_MODE:
  case PSI_PI_T_MODE:
  case PSI_PI_D_MODE: /* flat prior and Dirichlet prior */
  case PSI_S_MODE:
  case PSI_P_MODE:
     *p_term = 0;
     break;
  case PSI_S_MODE_NORMAL:
     assert (0);
     psi_fatal ("not implemented: normal prior");
     break;
/*
     n = gmel_simple_gaussian_lnpdf (0, prior_normal_sd_s, theta_prime.s); 
     d = gmel_simple_gaussian_lnpdf (0, prior_normal_sd_s, theta.s); 
     *p_term = n - d;
*/
  case PSI_P_MODE_NORMAL:
     assert (0);
     psi_fatal ("not implemented: normal prior");
     break;
/*
     n = gmel_simple_gaussian_lnpdf (0, prior_normal_sd_p, theta_prime.p); 
     d = gmel_simple_gaussian_lnpdf (0, prior_normal_sd_p, theta.p); 
     *p_term = n - d;
*/
  default:
     assert (0);
     psi_fatal ("invalid mode in get_prior");
  }

  return EXIT_SUCCESS;
}

int 
get_likelihood (double *p_term,
                parameter *p_theta, parameter *p_theta_prime)
                          /* top */           /* bot */       /* -> WRONG! */
                          /* bot */           /* top */
{
  /*!\note
     We will calcualte 
     \f$\sum_{h=1}^M e^{-2(s-s^*)E_s(\eta^{(h)}) - 2(p-p^*)E_p(\eta^{(h)}) }
        \prod_{n=1}^N \frac{\pi_{\eta_{n}^{(h)}}}{\pi_{\eta_{n}^{(h)}}^*}\f$
   */
  int r = EXIT_SUCCESS;

  double gibbs_top, gibbs_bot;  // same as const_top and const_bot
  double max_energy_top, max_energy_bot;
  parameter theta_star;

  /* New gridpoint is assigned to p_theta_prime's gp member */
  /* BUGGY */
  r = locate_gridpoint (&theta_star, p_theta, p_theta_prime);
  
  /* The messy gibbs sampling and summing */
  /* top gibbs */
  if ( p_theta->gp.a != p_theta_prime->gp.a 
       || p_theta->gp.c != p_theta_prime->gp.c
       || p_theta->gp.g != p_theta_prime->gp.g
       || p_theta->gp.s != p_theta_prime->gp.s
       || p_theta->gp.p != p_theta_prime->gp.p) {
     r = ERR_TOOMANYZEROS;
     while (r == ERR_TOOMANYZEROS) { 
        r = gmel_bf_likelihood_sumup_gibbs (&gibbs_top, 
                                            &max_energy_top, 
                                            *p_theta, 
                                            theta_star);
        if (want_debug_gibbs == 1 && r == ERR_TOOMANYZEROS) { 
           fprintf (stderr, "loop again because of ERR_TOOMANYZEROS\n");
        }
     }
  } else {
/* HERE */
     /* BUGGY: the first run after init_theta or the call
               of gmel_bf_likelihood_sumup_gibbs in the init_theta
               should not enter here */
     gibbs_top = p_theta->gibbs_sum;
     max_energy_top = p_theta->max_energy; /* will be removed */
  }

  /* bot gibbs */
  r = ERR_TOOMANYZEROS;
  while (r == ERR_TOOMANYZEROS) { 
     r = gmel_bf_likelihood_sumup_gibbs (&gibbs_bot, 
                                         &max_energy_bot, 
                                         *p_theta_prime, 
                                         theta_star);
     if (want_debug_gibbs == 1 && r == ERR_TOOMANYZEROS) { 
        fprintf (stderr, "loop again because of ERR_TOOMANYZEROS\n");
     }
  }

  /* for later update if the prime param is accepted */
  p_theta_prime->gibbs_sum = gibbs_bot; 
  p_theta_prime->max_energy = max_energy_bot;;

  double LnDiff_A, LnDiff_C, LnDiff_G, LnDiff_T, LnDiff_S, LnDiff_P;

  assert(p_theta->a != 0);
  assert(p_theta->c != 0);
  assert(p_theta->g != 0);
  assert(p_theta->t != 0);
  assert(p_theta_prime->a != 0);
  assert(p_theta_prime->c != 0);
  assert(p_theta_prime->g != 0);
  assert(p_theta_prime->t != 0);

  r = psi_err_base (p_theta->a, p_theta->c, p_theta->g, p_theta->t);
  r = psi_err_base (p_theta_prime->a, p_theta_prime->c,
                    p_theta_prime->g, p_theta_prime->t);

  LnDiff_A = log (p_theta_prime->a) - log (p_theta->a);
  LnDiff_C = log (p_theta_prime->c) - log (p_theta->c);
  LnDiff_G = log (p_theta_prime->g) - log (p_theta->g);
  LnDiff_T = log (p_theta_prime->t) - log (p_theta->t);
  LnDiff_S = p_theta_prime->s - p_theta->s;
  LnDiff_P = p_theta_prime->p - p_theta->p;

  if (!finite(LnDiff_S) || !finite(LnDiff_P)
      || !finite(LnDiff_A) || !finite(LnDiff_C)
      || !finite(LnDiff_G) || !finite(LnDiff_T) ) {
     if (want_debug_gibbs == 1) {
/*
        debug_write_parameter (*p_theta);
        debug_write_parameter (*p_theta_prime);
*/
     }
     psi_fatal ("diff is not finite");
  }

  /* *p_term = log(gibbs_top/gibbs_bot) */
  *p_term = gibbs_top - gibbs_bot
            + data_num_A * LnDiff_A
            + data_num_C * LnDiff_C
            + data_num_G * LnDiff_G
            + data_num_T * LnDiff_T
            + data_minus2E_solv * LnDiff_S
            + data_minus2E_pair * LnDiff_P;
            /* + max_energy_top - max_energy_bot; */
  if (want_debug_est == 1)
    {
      fprintf (stderr, "\tL: %.3lf\t%.3lf - %.3lf, (%d) %.3lf, (%d) %.3lf, (%d) %.3lf, (%d) %.3lf, (%.3lf) %.3lf, (%.3lf) %.3lf\n",
            *p_term, gibbs_top, gibbs_bot
            , data_num_A, LnDiff_A
            , data_num_C, LnDiff_C
            , data_num_G, LnDiff_G
            , data_num_T, LnDiff_T
            , data_minus2E_solv, LnDiff_S
            , data_minus2E_pair, LnDiff_P);
    }

  if (!finite (*p_term)) 
    {
      psi_fatal ("likelihood result is not finite");
    }
  return EXIT_SUCCESS;
}

int 
locate_gridpoint (parameter *p_theta_star, 
                  parameter *p_theta, parameter *p_theta_prime)
{
  int r = EXIT_SUCCESS;
  /* Find a middle value and locate its gridpoint */
  p_theta_star->s = (p_theta->s + p_theta_prime->s)/2;
  p_theta_star->p = (p_theta->p + p_theta_prime->p)/2;
  p_theta_star->a = (p_theta->a + p_theta_prime->a)/2;
  p_theta_star->c = (p_theta->c + p_theta_prime->c)/2;
  p_theta_star->g = (p_theta->g + p_theta_prime->g)/2;
  p_theta_star->t = 1.0 - p_theta_star->a - p_theta_star->c - p_theta_star->g;

  gmel_grid_locate (&(p_theta_star->gp), *p_theta_star);
  gmel_grid_param_gridpoint (p_theta_star, p_theta_star->gp);
  p_theta_prime->gp = p_theta_star->gp; /* for update of current gp */

  assert( p_theta_star->gp.a >= 0 && p_theta_star->gp.a < gridpoint_n_number );
  assert( p_theta_star->gp.c >= 0 && p_theta_star->gp.c < gridpoint_n_number );
  assert( p_theta_star->gp.g >= 0 && p_theta_star->gp.g < gridpoint_n_number );
  assert( p_theta_star->gp.s >= 0 && p_theta_star->gp.s < gridpoint_s_number );
  assert( p_theta_star->gp.p >= 0 && p_theta_star->gp.p < gridpoint_p_number );

  psi_err_base (p_theta_star->a, p_theta_star->c,
                p_theta_star->g, p_theta_star->t);

  return r;
}

int 
create_gridpoints (int n)
{
  int r = EXIT_SUCCESS;

  /* Grid-Points Setup */
  psi_grid_create_multi_gibbs (n);
  r = gmel_grid_new (gridpoint_s_begin, gridpoint_s_end, gridpoint_s_number,
                     gridpoint_p_begin, gridpoint_p_end, gridpoint_p_number,
                     gridpoint_n_begin, gridpoint_n_end, gridpoint_n_number,
                     gridpoint_n_begin, gridpoint_n_end, gridpoint_n_number,
                     gridpoint_n_begin, gridpoint_n_end, gridpoint_n_number,
                     gibbs_size);

  return r;
}

int 
load_gridpoints (const char *fn)
{
  int r = EXIT_SUCCESS;
  r = gmel_grid_load (fn);
  return r;
}

int 
save_gridpoints (const char *fn)
{
  int r = EXIT_SUCCESS;
  r = gmel_grid_save (fn);
  return r;
}

int delete_gridpoints() 
{
  gmel_grid_del ();
  return EXIT_SUCCESS;
}

double Maxof3(double a, double b, double c)
{
  if (a < b)
  {
      if (b < c)
      {
          return c;
      }
      else
      {
          return b;
      }
  }
  else
  {
      if (a < c)
      {
          return c;
      }
      else
      {
          return a;
      }
  }
}

int 
init_theta (parameter *p_theta)
{
  int r = EXIT_SUCCESS;
  /* INIT ENERGY Calculation: done initialize_drevol */

  /* INIT VALUES */
  p_theta->s = init_s;
  p_theta->p = init_p;
  p_theta->a = init_a;
  p_theta->c = init_c;
  p_theta->g = init_g;
  p_theta->t = 1.0 - (init_a + init_c + init_g);

  assert( p_theta->a >= 0 && p_theta->a <= 1 );
  assert( p_theta->c >= 0 && p_theta->c <= 1 );
  assert( p_theta->g >= 0 && p_theta->g <= 1 );
  assert( p_theta->t >= 0 && p_theta->t <= 1 );

  return r;
}

/* set the gp members and gibbs sum */
int 
psi_mcmc_init_theta_gp_gibbs_sum (parameter *p_theta)
{
  int r = EXIT_SUCCESS;

  assert( p_theta->a >= 0 && p_theta->a <= 1 );
  assert( p_theta->c >= 0 && p_theta->c <= 1 );
  assert( p_theta->g >= 0 && p_theta->g <= 1 );
  assert( p_theta->t >= 0 && p_theta->t <= 1 );

  /* find the grid point based on the current values of theta */
  r = gmel_grid_locate (&(p_theta->gp), *p_theta);

  /* check if the gibbs sum is equal to log (gibbs_size) */
  parameter theta_star; /*  This theta_star is a dummy because D_MODE of 
                            sumup_gibbs deals with only current theta.
                            This sumup_gibbs fill theta with gibbs_top and
                            max_energy so that the first calling of 
                            sumup_gibbs can perform correctly and smoothly. */
  theta_star = *p_theta;
  r = ERR_TOOMANYZEROS;
/*
  fprintf (stderr, "p_theta: ");
  debug_write_parameter (stderr, *p_theta);
  fprintf (stderr, "theta_s: ");
  debug_write_parameter (stderr, theta_star);
*/
  while (r == ERR_TOOMANYZEROS) { 
     r = gmel_bf_likelihood_sumup_gibbs (&(p_theta->gibbs_sum), 
                                         &(p_theta->max_energy), 
                                         *p_theta, theta_star);
     if (want_debug_gibbs == 1 && r == ERR_TOOMANYZEROS) { 
        fprintf (stderr, "loop again because of ERR_TOOMANYZEROS\n");
     }
  }
  if (p_theta->gibbs_sum < log (gibbs_size) - 0.001)
    {
      psi_fatal ("psi_mcmc_init_theta_gp_gibbs_sum: %lf < %lf", p_theta->gibbs_sum, log (gibbs_size));
    }
  if (p_theta->gibbs_sum > log (gibbs_size) + 0.001)
    {
      psi_fatal ("psi_mcmc_init_theta_gp_gibbs_sum: %lf > %lf", p_theta->gibbs_sum, log (gibbs_size));
    }


  /* New gridpoint is assigned to p_theta_prime's gp member */
  /* BUGGY */
  /* find the middle point of p_theta and theta_prime */
/* 
  double gibbs_top;  
  double max_energy_top;
  parameter theta_prime;
  theta_prime = *p_theta;

  r = locate_gridpoint (&theta_star, p_theta, &theta_prime);
  
  fprintf (stderr, "####################\n");
  fprintf (stderr, "p_theta: ");
  debug_write_parameter (stderr, *p_theta);
  fprintf (stderr, "theta_s: ");
  debug_write_parameter (stderr, theta_star);
  fprintf (stderr, "theta_p: ");
  debug_write_parameter (stderr, theta_prime);

  r = ERR_TOOMANYZEROS;
  while (r == ERR_TOOMANYZEROS) { 
     r = gmel_bf_likelihood_sumup_gibbs (&gibbs_top, 
                                         &max_energy_top, 
                                         *p_theta, 
                                         theta_star);
     if (want_debug_gibbs == 1 && r == ERR_TOOMANYZEROS) { 
        fprintf (stderr, "loop again because of ERR_TOOMANYZEROS\n");
     }
  }
  if (gibbs_top < log (gibbs_size) - 0.001)
    {
      psi_fatal ("psi_mcmc_init_theta_gp_gibbs_sum 2: %lf < %lf", gibbs_top, log (gibbs_size));
    }
  if (gibbs_top > log (gibbs_size) + 0.001)
    {
      psi_fatal ("psi_mcmc_init_theta_gp_gibbs_sum 2: %lf > %lf", gibbs_top, log (gibbs_size));
    }
  p_theta->gibbs_sum = gibbs_top;
*/
  p_theta->gp.s = -1;
  p_theta->gp.p = -1;
  p_theta->gp.a = -1;
  p_theta->gp.c = -1;
  p_theta->gp.g = -1;
  return r;
}

int correct_boundary( int *gp, int max ) {
  if (*gp == max) { 
     *gp = max - 1;
  } else if (*gp < 0) { 
     assert( 0 );
     *gp = 0;
  }
  return EXIT_SUCCESS;
}

/*
int
gmel_gsl_init() {
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  gmel_gsl_r = gsl_rng_alloc (T);
  gsl_rng_set (gmel_gsl_r, time(NULL));
  return EXIT_SUCCESS;
}

int
gmel_gsl_fin() {
  gsl_rng_free (gmel_gsl_r);
  return EXIT_SUCCESS;
}
*/

int
write_acceptance_ratio ()
{
  fprintf (ofile, "Acceptance:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
          (double) gmel_news_accept_n*100/gmel_news_total_n,
          (double) gmel_news_accept_s*100/gmel_news_total_s,
          (double) gmel_news_accept_p*100/gmel_news_total_p,
          (double) delta_s_accepted/gmel_news_accept_s,
          (double) delta_p_accepted/gmel_news_accept_p,
          (double) delta_s_rejected/(gmel_news_total_s - gmel_news_accept_s),
          (double) delta_p_rejected/(gmel_news_total_p - gmel_news_accept_p));
    
  if (want_verbose == 1)
    {
      fprintf(stderr, "Acceptance Ratio: N - %.1lf %%, S - %.1lf %%, P - %.1lf %%\n", 
              (double) gmel_news_accept_n*100/gmel_news_total_n,
              (double) gmel_news_accept_s*100/gmel_news_total_s,
              (double) gmel_news_accept_p*100/gmel_news_total_p);
      fprintf(stderr, "Average Accepted Deltas: DS - %.3lf, DP - %.3lf\n", 
              (double) delta_s_accepted/gmel_news_accept_s,
              (double) delta_p_accepted/gmel_news_accept_p);
      fprintf(stderr, "Average Rejected Deltas: DS - %.3lf, DP - %.3lf\n", 
              (double) delta_s_rejected/(gmel_news_total_s - gmel_news_accept_s),
              (double) delta_p_rejected/(gmel_news_total_p - gmel_news_accept_p));
    }
  return EXIT_SUCCESS;
}

/*
int
get_sensible_pi (double *pi)
{
  int i;
  for (i = 0; i < NUM_NUCLEOTIDE; i++) 
    {
      if (pi[i] < 0.1 * delta_n_mu)
          pi[i] = 0.1 * delta_n_mu;
    }
  return EXIT_SUCCESS;
}
*/


/* need to be exported from MCMC module */
int
execute_estimation (int is_fixed_s, int is_fixed_p, double s, double p,
                    parameter *theta_star)
{
  int r = EXIT_SUCCESS;
  int i, j;
  int total_iteration;
  parameter theta;

/*
  theta_star->s = 0;
  theta_star->p = 0;
  theta_star->a = 0;
  theta_star->c = 0;
  theta_star->g = 0;
  theta_star->t = 0;
*/

  FILE *pfile = NULL;
  FILE *afile = NULL;

  int last_burn_i;
  int last_i;
  int r_burn;
  int last_freq_i = 0;
  int r_sample;

  if (is_fixed_s == 1) {
    want_disable_s = 1;
    init_s = s;
  } else {
    want_disable_s = 0;
  }

  if (is_fixed_p == 1) {
    want_disable_p = 1;
    init_p = p;
  } else {
    want_disable_p = 0;
  }

  r_burn = psi_io_read_gen_last_burn (pname, &last_burn_i);
  r_burn = psi_io_read_theta_last_burn (pname, &theta);
  if (r_burn == EXIT_SUCCESS) 
    {
      psi_io_remove_number_burn (pname);
      /* last_i = psi_io_index_last_theta (pname, sample_size); */
      last_i = 0;
      last_freq_i = 0;
/* fprintf (stderr, "burn found: %d\n", last_i); */
    }
  else if (r_burn == EXIT_FAILURE) 
    {
      r_sample = psi_io_read_gen_last_sample (pname, &last_freq_i);
      last_freq_i++;
      r_sample = psi_io_read_theta_last_sample (pname, &theta);
      if (r_sample == EXIT_SUCCESS)
        {
          psi_io_remove_number_sample (pname);
          last_i = psi_io_index_last_theta (pname, sample_size);
          last_burn_i = sample_burn;
          if (want_verbose == 1)
            {
              fprintf (stderr, "in reading sample state: last_i = %d\n", 
                       last_i);
            }
        }
      else /* this will be the first of firsts iteration for hpc */
        {
          last_burn_i = 0;
          last_freq_i = 0;
          last_i = 0;
          r = init_theta (&theta);
          if (r != EXIT_SUCCESS)
            psi_fatal ("init theta failure"); 
          r = psi_mcmc_init_theta_gp_gibbs_sum (&theta);
          if (r != EXIT_SUCCESS)
            psi_fatal ("init theta failure"); 
        }
/* fprintf (stderr, "burn not found: %d\n", last_i);
      last_i = psi_io_index_last_theta (pname, sample_size);
      if (last_i == 0)
        {
          last_burn_i = 0;
        }
      else
        {
          last_burn_i = sample_burn;
        }
*/
    }

  pfile = fopen (pname, "a");
  if (pfile == NULL) {
     psi_fatal ("could not open the parameter file: pfile");
  }
  if (want_save_alldraws == 1) {
     afile = fopen (aname, "a");
     if (afile == NULL) {
        psi_fatal ("could not open the parameter file: afile");
     }
  }

  if (want_verbose == 1)
    { 
      fprintf (stderr, "param: ");
      write_parameter (stderr, theta, 0);
    }

  /*!\note
     Burn-in
   */ 
  /* gibbs_burn = gibbs_burn * Dat->len_dna; */
  if (want_debug_est == 1) {
     /* debug_write_parameter (theta); */
  }
  
  if (want_verbose == 1)
    {
      fprintf (stderr, "sample burn: %d\n", sample_burn);
      fprintf (stderr, "sample freq: %d\n", sample_freq);
      fprintf (stderr, "sample size: %d\n", sample_size);
    }
  fprintf (ofile, "EXE_ESTIMATION:\t%d\t%d\t%d\n", sample_burn, sample_freq, sample_size);

  total_iteration = sample_burn + sample_freq * sample_size;
  if (want_debug_input == 1) {
     fprintf (ofile, "Burnin: %d\n", sample_burn);
  }

  int gen_alldraws = 0;
  if (want_save_alldraws == 1)
    write_parameter_header (afile);

  for ( i = last_burn_i; i < sample_burn; i++ ) {
     if (want_log == 1 && i % log_freq == 0) {
        time_t t = time(NULL);
        fprintf (stderr, "%s:\tburnin %d\n", asctime(localtime(&t)), i);
     }
     if (want_save_alldraws == 1) {
        write_parameter (afile, theta, gen_alldraws);
        gen_alldraws++;
     }
     if (want_disable_n == 0) { 
        r = new_pi (&theta); /* 1: first */
        if (want_debug_est == 1) {
           /* debug_write_parameter (theta); */
        }
     }
     if (want_disable_s == 0) { 
        r = new_s (&theta);
        if (want_debug_est == 1) {
           /* debug_write_parameter (theta); */
        }
     }
     if (want_disable_p == 0) { 
        r = new_p (&theta);
        if (want_debug_est == 1) {
           /* debug_write_parameter (theta); */
        }
     } 
     r = check_time ();
     if (r != EXIT_SUCCESS)
       {
         psi_io_write_number_burn (pfile, i);
         write_parameter (pfile, theta, i);
         psi_io_write_mcmc_acceptance (pfile,
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
         i = sample_burn;    
       }
  }

  if (r != EXIT_SUCCESS)
    {
      fclose (pfile);
      pfile = NULL;
      if (want_save_alldraws == 1) 
        {
          fclose (afile);
          afile = NULL;
        }
      return r;
    }

  /*!\note
     MCMC
   */
  for (i = last_i; i < sample_size; i++) {
     if (want_verbose == 1)
       fprintf (stderr, "%d: ", i); 
     for ( j = last_freq_i; j < sample_freq; j++ ) { 
        if (want_save_alldraws == 1) {
           write_parameter (afile, theta, gen_alldraws);
           gen_alldraws++;
        }
        if (want_disable_n == 0) { 
           r = new_pi (&theta); /* first: 1, what does this mean? */
           if (want_debug_est == 1) {
              debug_write_parameter (ofile, theta);
           }
        }
        if (want_disable_s == 0) { 
           r = new_s (&theta);
           if (want_debug_est == 1) {
              debug_write_parameter (ofile, theta);
           }
        }
        if (want_disable_p == 0) { 
           r = new_p (&theta);
           if (want_debug_est == 1) {
              debug_write_parameter (ofile, theta);
           }
        }
        /* fprintf (stderr, "."); */
        r = check_time ();
        if (r != EXIT_SUCCESS)
          {
            psi_io_write_number_sample (pfile, j);
            write_parameter (pfile, theta, j);
            psi_io_write_mcmc_acceptance (pfile,
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
            j = sample_freq;
            i = sample_size;
          }
     }
     last_freq_i = 0;
     /* fprintf (stderr, "\n"); */
     if (want_log == 1 && i % log_freq == 0) {
        time_t t = time(NULL);
        fprintf (stderr, "%s:\testimation %d\n", asctime(localtime(&t)), i);
     }
     if (i == 0) {
        write_parameter_header (pfile);
     } 
     if (i < sample_size) {
        write_parameter (pfile, theta, i);
     } 

/*
     theta_star->s += theta.s;
     theta_star->p += theta.p;
     theta_star->a += theta.a;
     theta_star->c += theta.c;
     theta_star->g += theta.g;
*/

  }

/*
  theta_star->s /= sample_size;
  theta_star->p /= sample_size;
  theta_star->a /= sample_size;
  theta_star->c /= sample_size;
  theta_star->g /= sample_size;
  theta_star->t = 1.0 - theta_star->a - theta_star->c - theta_star->g;
*/

  fclose (pfile);
  pfile = NULL;

  if (want_save_alldraws == 1) {
     fclose (afile);
     afile = NULL;
  }

  return r;
}

int 
run_mcmc_chain (double **posterior, int b, int f, int s, int o)
{
  /* posterior is assumed to be 6-by-(s + o) matrix */
  int r;
  int i, j;
  int total_iteration;
  parameter theta;

  fprintf (ofile, "RUN_MCMC_CHAIN:\t%d\t%d\t%d\n", b, f, s);

  /* check want_disable_s and want_disable_p */
  /* init values: BUGGY */
  r = init_theta (&theta);
  if (r != EXIT_SUCCESS)
    return r;
  r = psi_mcmc_init_theta_gp_gibbs_sum (&theta);
  if (r != EXIT_SUCCESS)
    return r;
  if (want_debug_est == 1)
    {
      fprintf (stderr, "param: ");
      write_parameter (stderr, theta, 0);
    }

  /*!\note
     Burn-in
   */ 
  if (want_debug_est == 1)
    {
      debug_write_parameter (ofile, theta);
      fprintf (stderr, "sample burn: %d\n", b);
      fprintf (stderr, "sample freq: %d\n", f);
      fprintf (stderr, "sample size: %d\n", s);
    }

  for ( i = 0; i < b; i++ ) 
    {
      if (want_log == 1 && i % log_freq == 0)
        {
          time_t t = time(NULL);
          fprintf (stderr, "%s:\tburnin %d\n", asctime(localtime(&t)), i);
        }
      if (want_disable_n == 0) 
        { 
          new_pi (&theta); /* 1: first */
          if (want_debug_est == 1)
            debug_write_parameter (ofile, theta);
        }
      if (want_disable_s == 0)  
        { 
          new_s (&theta);
          if (want_debug_est == 1) 
            debug_write_parameter (ofile, theta);
        }
      if (want_disable_p == 0) 
        { 
          r = new_p (&theta);
          if (want_debug_est == 1) 
            debug_write_parameter (ofile, theta);
        } 
    }
  total_iteration = i;

  /*!\note
     MCMC
   */
  for (i = 0; i < s; i++) 
    {
      if (want_log == 1 && i % log_freq == 0)
        {
          time_t t = time(NULL);
          fprintf (stderr, "%s:\titerations %d\n", asctime(localtime(&t)), i);
        }

      for ( j = 0; j < f; j++ ) 
        { 
          if (want_disable_n == 0) 
            { 
              r = new_pi (&theta); /* first: 1, what does this mean? */
              if (want_debug_est == 1)
                debug_write_parameter (ofile, theta);
            }
          if (want_disable_s == 0) 
            { 
              r = new_s (&theta);
              if (want_debug_est == 1)
                debug_write_parameter (ofile, theta);
            }
          if (want_disable_p == 0) 
            { 
              r = new_p (&theta);
              if (want_debug_est == 1)
                debug_write_parameter (ofile, theta);
            }
        }
      total_iteration += j;
      posterior[PSI_SAMPLE_A][o + i] = theta.a;
      posterior[PSI_SAMPLE_C][o + i] = theta.c;
      posterior[PSI_SAMPLE_G][o + i] = theta.g;
      posterior[PSI_SAMPLE_T][o + i] = theta.t;
      posterior[PSI_SAMPLE_S][o + i] = theta.s;
      posterior[PSI_SAMPLE_P][o + i] = theta.p;
    }

  return total_iteration;
}

int
find_mcmc_option ()
{
  int i, N, N_e, B, F, T;
  int accumulated_burnin;
  double tau = 1.5;
  double **posterior = NULL;

  posterior = XMALLOC (double *, 6);
  for (i = 0; i < 6; i++)
    posterior[i] = XMALLOC (double, sample_size);
  
  /* find the burnin period and sampling frequency */
  N = sample_size;
  N_e = 0;
  B = sample_burn;
  F = sample_freq;
  accumulated_burnin = B; 
  setup_mcmc_init (0.0, 0.0, 0.25, 0.25, 0.25);
  while (N > 0)
    {
      /* mcmc.c */
      if (want_verbose == 1)
        fprintf (stderr, "finding sampling freq. until we have a good tau: %lf\n", tau);
      T = run_mcmc_chain (posterior, B, F, N, N_e); /* N + N_e = sample_size */
      /* gslwrap.c */
      tau = get_longest_autocorr_time_sample (posterior, 6, N); 
      /* if (tau < 1.1)  */
      if (tau < limit_tau) 
        break;
      /* assert (tau > 1.0); */
      N_e = (int) sample_size/tau;
/* fprintf (stderr, "N_e: %d, tau: %lf\n", N_e, tau); */
      /* N_e = (int) N/tau; */
      F = F + tau;
      accumulated_burnin += T - B - N;
      /* gslwrap.c */
      take_subserial_sample (posterior, 6, sample_size, N_e);
      N = sample_size - N_e;
      B = 0;
      setup_mcmc_init (posterior[PSI_SAMPLE_S][sample_size - 1],
                       posterior[PSI_SAMPLE_P][sample_size - 1],
                       posterior[PSI_SAMPLE_A][sample_size - 1],
                       posterior[PSI_SAMPLE_C][sample_size - 1],
                       posterior[PSI_SAMPLE_G][sample_size - 1]);
    }

  sample_burn = accumulated_burnin;
  sample_freq = F;
  fprintf (ofile, "FIND_MCMC_OPTION:\t%d\t%d\n", sample_burn, sample_freq);
  
  if (want_verbose == 1)
    fprintf (stderr, "B: %d, F: %d\n", sample_burn, sample_freq);

  for (i = 0; i < 6; i++)
    XFREE (posterior[i]);
  XFREE (posterior);

  return EXIT_SUCCESS;
}

int
run_mcmc_multichain (int number_chains, double **ip, int is_fixed_run)
{
  int i, j;

  double psr[6]; 
  double psr_max = 100; /* just for passing the while statement */
  double ***posteriors = NULL;
  posteriors = XMALLOC (double **, number_chains);
  for (i = 0; i < number_chains; i++)
    { 
      posteriors[i] = XMALLOC (double *, 6);
      for (j = 0; j < 6; j++)
        { 
          posteriors[i][j] = XMALLOC (double, sample_size);
        } 
    }

  double **draws = NULL;
  draws = XMALLOC (double *, number_chains);
  for (i = 0; i < number_chains; i++) {
     draws[i] = XMALLOC (double, sample_size);
  }

  /* multi gibbs grid for multi-chains */
  int iter_conv = -1;
  while (psr_max > limit_psr) 
    {
      iter_conv++;
      for (i = 0; i < number_chains; i++)
        { 
          /* change the gibbs grid */
          psi_grid_change_multi_gibbs (i);
          setup_mcmc_init (ip[i][PSI_SAMPLE_S],
                           ip[i][PSI_SAMPLE_P],
                           ip[i][PSI_SAMPLE_A],
                           ip[i][PSI_SAMPLE_C],
                           ip[i][PSI_SAMPLE_G]);
          if (want_verbose == 1)
            fprintf (stderr, "multi chains: %d - %d\n", iter_conv, i);
          run_mcmc_chain (posteriors[i], sample_burn, sample_freq, sample_size, 0);
          ip[i][PSI_SAMPLE_S] = posteriors[i][PSI_SAMPLE_S][sample_size - 1];
          ip[i][PSI_SAMPLE_P] = posteriors[i][PSI_SAMPLE_P][sample_size - 1];
          ip[i][PSI_SAMPLE_A] = posteriors[i][PSI_SAMPLE_A][sample_size - 1];
          ip[i][PSI_SAMPLE_C] = posteriors[i][PSI_SAMPLE_C][sample_size - 1];
          ip[i][PSI_SAMPLE_G] = posteriors[i][PSI_SAMPLE_G][sample_size - 1];
        }
      copy_draws (draws, posteriors, PSI_SAMPLE_A, number_chains);
      gmel_conv_psr (&psr[PSI_SAMPLE_A], draws, number_chains, sample_size);
      copy_draws (draws, posteriors, PSI_SAMPLE_C, number_chains);
      gmel_conv_psr (&psr[PSI_SAMPLE_C], draws, number_chains, sample_size);
      copy_draws (draws, posteriors, PSI_SAMPLE_G, number_chains);
      gmel_conv_psr (&psr[PSI_SAMPLE_G], draws, number_chains, sample_size);
      copy_draws (draws, posteriors, PSI_SAMPLE_T, number_chains);
      gmel_conv_psr (&psr[PSI_SAMPLE_T], draws, number_chains, sample_size);
      copy_draws (draws, posteriors, PSI_SAMPLE_S, number_chains);
      gmel_conv_psr (&psr[PSI_SAMPLE_S], draws, number_chains, sample_size);
      copy_draws (draws, posteriors, PSI_SAMPLE_P, number_chains);
      gmel_conv_psr (&psr[PSI_SAMPLE_P], draws, number_chains, sample_size);
      psr_max = psi_max (psr, 6);
      if (is_fixed_run == 1)
        {
          fprintf (ofile, "PSR:\t%lf\n", psr_max);
          break;
        }
    }
  assert (iter_conv > -1);

  if (is_fixed_run == 0)
    sample_burn += sample_freq * sample_size * iter_conv; 

  for (i = 0; i < number_chains; i++) {
     XFREE (draws[i]);
  }
  XFREE (draws);

  /* find the means of the multi chains and set the init values to them */
  psi_mcmc_find_means (posteriors, number_chains);
  /* save posterior samples */
  psi_mcmc_write_posterior_sample (pname, posteriors, number_chains, sample_size);

  for (i = 0; i < number_chains; i++)
    { 
      for (j = 0; j < 6; j++)
        { 
          XFREE (posteriors[i][j]);
        } 
      XFREE (posteriors[i]);
    }
  XFREE (posteriors);
  return EXIT_SUCCESS;
}

static void
psi_mcmc_find_means (double ***posteriors, int number_chains)
{
  int i, j;
  double ip[PSI_NUM_SAMPLE];

  for (i = 0; i < PSI_NUM_SAMPLE; i++)
    {
      ip[i] = 0.0;
      /* find the means */
      for (j = 0; j < number_chains; j++)
        ip[i] += psi_stats_mean (posteriors[j][i], sample_size);
      ip[i] /= number_chains;
    }
  setup_mcmc_init (ip[PSI_SAMPLE_S],
                   ip[PSI_SAMPLE_P],
                   ip[PSI_SAMPLE_A],
                   ip[PSI_SAMPLE_C],
                   ip[PSI_SAMPLE_G]);
}

static void 
copy_draws (double **draws, double ***posteriors, int k, int n_c)
{
  int i, j;
  for (i = 0; i < n_c; i++) 
    {
      for (j = 0; j < sample_size; j++) 
        {
          draws[i][j] = posteriors[i][k][j];
        }
    }
}

/* EXIT_FAILURE: overtime */
static int 
check_time ()
{
  time_t t;
  double d;
  if (limit_time != 0)
    {
      t = time(NULL);
      d = difftime (t, start_time);
      d /= 60.0;
      /* fprintf (stderr, "%lf > %d\n", d, limit_time); */
      if (d > limit_time)
        {
          /* fprintf (stderr, "overtime!!!\n"); */
          return EXIT_OTHERS;
        }
    }
  return EXIT_SUCCESS;
}

/* check the header file for their existence */ 
double 
psi_mcmc_s_flat_proposal_pdf (double x, double mu_s)
{
  return psi_ran_flat_pdf (x, mu_s, delta_s, flatprior_s_begin, flatprior_s_end);
}

double 
psi_mcmc_p_flat_proposal_pdf (double x, double mu_p)
{
  return psi_ran_flat_pdf (x, mu_p, delta_p, flatprior_p_begin, flatprior_p_end);
}

double 
psi_mcmc_s_normal_proposal_lnpdf (double x, double mu_s) 
{
  return psi_gaussian_lnpdf (mu_s, delta_s_sd, x); 
}

double 
psi_mcmc_p_normal_proposal_lnpdf (double x, double mu_p) 
{
  return psi_gaussian_lnpdf (mu_p, delta_p_sd, x); 
}

double
psi_mcmc_s_flat_prior_lnpdf (double s)
{
  assert (flatprior_s_begin < s && s < flatprior_s_end);
  return log_flatprior_s;
} 

double
psi_mcmc_p_flat_prior_lnpdf (double p)
{
  assert (flatprior_p_begin < p && p < flatprior_p_end);
  return log_flatprior_p;
}

double
psi_mcmc_s_normal_prior_lnpdf (double s)
{
  return psi_gaussian_lnpdf (0.0, prior_normal_sd_s, s);
}

double
psi_mcmc_p_normal_prior_lnpdf (double p)
{
  return psi_gaussian_lnpdf (0.0, prior_normal_sd_p, p);
}

double
psi_mcmc_n_dirichlet_proposal_lnpdf (double alpha[], double theta[])
{
  int i;
  for (i = 0; i < 4; i++)
    {
      alpha[i] *= delta_n_mu;
    }
  return psi_pi_lnpdf (alpha, theta);
}

double
psi_mcmc_n_likelihood (double a, double c, double g, double t)
{
  return data_num_A * log (a) + data_num_C * log (c) 
         + data_num_G * log (g) + data_num_T * log (t); 
}
