/* bf.h -- Bayes factor estimation
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
#include "nsv.h"
#include "bf.h"

/* setup function is needed */
/* file */
static FILE *ofile;
static char *pname;
static char *qname;
static char *aname;
static char *bname;
static int sample_size;
static int limit_sample_size = 30;
static int log_sample_size;
/* want */
static int want_load_param = 0;
static int want_save_alldraws = 0;
/* state */
static int want_load_state = 0;

static int want_verbose = 1;

static int load_estimation (parameter *theta_star_model4,
                            parameter *theta_star_model1);
int
unsetup_bf ()
{
  XFREE (pname);
  XFREE (aname);
  XFREE (qname);
  XFREE (bname);
  return EXIT_SUCCESS;
}

int
setup_bf_tol (double t)
{
  assert (t <= 1.0 && t > 0.0);
  /* accessed through mcmc */
  sample_size = access_mcmc_sample_size ();
  limit_sample_size = sample_size * t;
  log_sample_size = log (sample_size);
  return EXIT_SUCCESS;
}

int 
setup_bf_file (FILE *fp, const char *pn, const char *an, const char *qn, 
               const char *bn)
{
  assert (pn != NULL);
  assert (an != NULL);
  assert (qn != NULL);
  assert (bn != NULL);
  ofile = fp;
  pname = xstrdup (pn);
  aname = xstrdup (an);
  qname = xstrdup (qn);
  bname = xstrdup (bn);
  return EXIT_SUCCESS;
}

/* into header yet */
int
/* setup_bf_want (int ls, int lp, int d) */
setup_bf_want (int lp, int d)
{
  /* want_load_state = ls; */
  want_load_param = lp;
  want_save_alldraws = d;
  return EXIT_SUCCESS;
}

int 
setup_bf_state (int s)
{
  want_load_state = s;
  return EXIT_SUCCESS;
}

static int want_debug_bf = 0;
static int want_use_more_than_one_gridpoint = 0;
static int want_use_Dirichlet_prior_n = 1;
static int want_use_normal_prior_s = 0;
static int want_use_normal_prior_p = 0;
static int want_use_Dirichlet_proposal_n = 1;
static int want_use_normal_proposal_s = 0;
static int want_use_normal_proposal_p = 0;
static int want_save_q = 1;
static int want_load_q = 0;

static int psi_bf_stage (char *fn);

static int
psi_st_parameter (parameter *param,
                   double s,
                   double p,
                   double a,
                   double c,
                   double g);

int
psi_st_parameter (parameter *param,
                   double s,
                   double p,
                   double a,
                   double c,
                   double g)
{
  param->s = s;
  param->p = p;
  param->a = a;
  param->c = c;
  param->g = g;
  param->t = 1.0 - a - c - g;
  return 0;
}

int 
psi_bf_stage (char *fn)
{
  return psi_io_check_bf_state (fn, sample_size);
}


/* BUGGY: check this and get_likelihood version give the same result */
/* these two functions are going to be replaced by those in mcmc.c */
/*
static int 
gmel_bf_likelihood_gibbs_ratio (double *ratio, parameter top, parameter bot);
static int
gmel_bf_likelihood_locate_middle_point (parameter *p_theta_star, 
                                        parameter top, parameter bot);
*/

static int
gmel_bf_load_q (FILE *fp, parameter *p_param, int n, int w_b);

/******************************************************************************
 Static Functions
******************************************************************************/

/*
static int
gmel_bf_likelihood_gibbs_ratio (double *ratio, parameter top, parameter bot)
{
  int r = EXIT_SUCCESS;
  double gibbs_top, gibbs_bot;
  double top_max_energy, bot_max_energy;
  parameter theta_star;
  double LnDiff_A, LnDiff_C, LnDiff_G, LnDiff_T;

  gmel_bf_likelihood_locate_middle_point(&theta_star, top, bot); 
  r = ERR_BF_TOOMANYZEROS;
  while (r == ERR_BF_TOOMANYZEROS) { 
     r = gmel_bf_likelihood_sumup_gibbs (&gibbs_top, &top_max_energy, 
                                         bot, theta_star);
  }
  if (r != EXIT_SUCCESS) {
     return r;     
  }
  r = ERR_BF_TOOMANYZEROS;
  while (r == ERR_BF_TOOMANYZEROS) { 
     r = gmel_bf_likelihood_sumup_gibbs (&gibbs_bot, &bot_max_energy, 
                                      top, theta_star);
  }
  if (r != EXIT_SUCCESS) {
     return r;     
  }

  r = psi_err_base (top.a, top.c, top.g, top.t);
  if (r != EXIT_SUCCESS) {
     return r;
  }
  r = psi_err_base (bot.a, bot.c, bot.g, bot.t);
  if (r != EXIT_SUCCESS) {
     return r;
  }
  LnDiff_A = log(top.a / bot.a);
  LnDiff_C = log(top.c / bot.c);
  LnDiff_G = log(top.g / bot.g);
  LnDiff_T = log(top.t / bot.t);

  *ratio = gibbs_top - gibbs_bot;
  assert (isfinite(*ratio));
  r = psi_err_finite (*ratio);
  if (r != EXIT_SUCCESS) {
     return r;
  }
  *ratio += data_num_A * LnDiff_A
            + data_num_C * LnDiff_C
            + data_num_G * LnDiff_G
            + data_num_T * LnDiff_T
            + data_minus2E_solv * (top.s - bot.s)
            + data_minus2E_pair * (top.p - bot.p)
            + (top_max_energy - bot_max_energy); 
  assert (isfinite(*ratio));
  r = psi_err_finite (*ratio);
  if (r != EXIT_SUCCESS) {
     return r;
  }

  if (want_verbose == 1)
    {
      verbose_gibbs (ofile, *ratio,   
                     data_num_A, data_num_C, data_num_G, data_num_T,
                     LnDiff_A, LnDiff_C, LnDiff_G, LnDiff_T, 
                     data_minus2E_solv, data_minus2E_pair,
                     gibbs_top, gibbs_bot,
                     top_max_energy, bot_max_energy,
                     top, bot);
    }
  return r;
}

static int
gmel_bf_likelihood_locate_middle_point (parameter *p_theta_star, 
                                        parameter top, parameter bot) 
{
  parameter middle_param;
  gridpoint gp;

  middle_param.s = (top.s + bot.s)/2;
  middle_param.p = (top.p + bot.p)/2;
  middle_param.a = (top.a + bot.a)/2;
  middle_param.c = (top.c + bot.c)/2;
  middle_param.g = (top.g + bot.g)/2;
  middle_param.t = 1.0 - middle_param.a - middle_param.c - middle_param.g;
  gmel_grid_locate (&gp, middle_param);
  p_theta_star->gp = gp;
  gmel_grid_param_gridpoint (p_theta_star, gp);
  
  return 0;
}
*/

/*
static int 
gmel_bf_likelihood_create_middle_points (parameter **middles,
                                         int *number_middles,
                                         parameter top, parameter bot)
{
  int r = EXIT_SUCCESS;

  return r;
}
*/

static int
gmel_bf_load_q (FILE *fp, parameter *p_param, int n, int w_b)
{
  int r = EXIT_SUCCESS;
  int gen;
  double a, c, g, t, s, p;

  /* move to the right position */
  if (n == 0) {
     if (w_b == MODEL1_BLOCK1) {
        locate_m1_block1 (fp, sample_size);
     } else if (w_b == MODEL4_BLOCK1) {
        locate_m4_block1_e1 (fp, sample_size);
     } else if (w_b == MODEL4_BLOCK2) {
        locate_m4_block2_e1 (fp, sample_size);
     } else if (w_b == MODEL4_BLOCK3) {
        locate_m4_block3_e1 (fp, sample_size);
     }
  } 

  /* read one line and store it into p_param */
  r = read_param_line (fp, &gen, &a, &c, &g, &t, &s, &p);
  if (r != EXIT_SUCCESS) {
     psi_fatal ("%s:%s:%d:%s; %s %d %lf %lf %lf %lf %lf %lf - r: %d\n",
                program_name,
                __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                "error of reading posterior draws of block 1 in model 1",
                gen, a, c, g, t, s, p, r);
  }
  p_param->a = a;
  p_param->c = c;
  p_param->g = g;
  p_param->t = t;
  p_param->s = s;
  p_param->p = p;

  return r;
}

/******************************************************************************
 BAYES FACTOR USER LEVEL Functions
******************************************************************************/

int 
gmel_dbf (double *bf, parameter param_top, parameter param_bot, FILE *bfp)
/* BUGGY remove bfp, the last argument and made it as file static */
{
  int r = EXIT_SUCCESS;
  double likelihood = 0;
  double prior = 0;
  double posterior = 0;

  /* make theta_start paramter valaues based on mcmc sample */
  /* get the theta_star that achieve the high posterior */
  
  /* calculate each term in logarithmic scale */
  r = gmel_bf_likelihood (&likelihood, 4, 1, param_top, param_bot);
  r = gmel_bf_prior (&prior, 4, 1, param_top, param_bot);
  r = gmel_bf_posterior (&posterior, 4, 1, param_top, param_bot);


  fprintf(bfp, "%lf\t%lf\t%lf\t", likelihood, prior, posterior);
  if (want_debug_bf == 1) {
     fprintf (ofile, "%s:%s:%d:%s; %s: %lf, %s: %lf, %s: %lf\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__,
              "likelihood", likelihood,
              "prior", prior,
              "posterior", posterior);
  }

  *bf = likelihood + prior - posterior;
  return r;
}

int 
gmel_bf_likelihood (double *likelihood, int model_top, int model_bot,
                    parameter p_top, parameter p_bot)
{
  int r = EXIT_SUCCESS;
  parameter ps[2];
  ps[0] = p_top;
  ps[1] = p_bot;
  int models[2] = { model_top, model_bot };
  double likelihoods[2] = { 0.0, 0.0 };
  int i;

  if ((models[0] == 1 || models[0] == 3) 
      && (models[1] == 1 || models[1] == 3)) {
     for (i = 0; i < 2; i++) { 
        switch (models[i]) {
        case 1:
           r = gmel_bf_likelihood_m1 (&likelihoods[i], 
                                      ps[i].a, ps[i].c, ps[i].g); 
           break;
        case 3:
           psi_fatal ("no implementation: model 3 comparison");
           break;
        default:
           psi_fatal ("no implementation: other model except 1 and 3");
        }
     }
     *likelihood = likelihoods[0] - likelihoods[1];
  } else {
     if (p_top.s == p_bot.s && p_top.p == p_bot.p 
         && p_top.a == p_bot.a && p_top.c == p_bot.c 
         && p_top.g == p_bot.g) {
        *likelihood = 0;
     } else {
        r = gmel_bf_likelihood_gibbs (likelihood, p_top, p_bot);
     }
  }
  return r;
}

int 
gmel_bf_likelihood_gibbs (double *log_likelihood, parameter top, parameter bot)
{
  int r = EXIT_SUCCESS;
  double log_m_likelihood;
  parameter m_top; 
  parameter m_bot; 
  *log_likelihood = 0;

/*
  parameter *m_middles = NULL;
  int n_middles = 0;
*/

  m_top = top;
  if (want_use_more_than_one_gridpoint == 1) {
     psi_fatal ("multipoint grid should be considered");
/*
     r = gmel_bf_likelihood_create_middle_points (&m_middles, &n_middles,
                                                  top, bot);
     for (i = 0; i < n_middles; i++) {
        m_bot = m_middles[i];
        r = get_likelihood (&log_m_likelihood, &m_bot, &m_top);
        *log_likelihood += log_m_likelihood; 
        if (want_debug_bf == 1) {
           verbose_likelihood_gibbs (ofile, log_m_likelihood, m_top, m_bot);
        }
     m_top = m_bot;
     }
*/
  }
  m_bot = bot; 
/*r = gmel_bf_likelihood_gibbs_ratio(&log_m_likelihood, m_top, m_bot); */
  r = get_likelihood (&log_m_likelihood, &m_bot, &m_top);
  *log_likelihood += log_m_likelihood; 
  if (want_debug_bf == 1) {
     verbose_likelihood_gibbs (ofile, log_m_likelihood, m_top, m_bot);
  }
  assert (isfinite(*log_likelihood));

  return r;
}

int 
gmel_bf_prior (double *prior, int model_top, int model_bot,
                         parameter p_top, parameter p_bot)
{
  int r = EXIT_SUCCESS;
  parameter ps[2];
  ps[0] = p_top;
  ps[1] = p_bot;
  int models[2] = { model_top, model_bot };
  double priors[2] = { 0, 0 };
  int i;

  for (i = 0; i < 2; i++) { 
     switch (models[i]) {
     case 1:
        r = gmel_bf_prior_m1 (&priors[i], ps[i].a, ps[i].c, ps[i].g); 
        break;
     case 2:
        psi_fatal ("no implementation: %s, %d, %s\n", 
                    __FILE__, __LINE__, __PRETTY_FUNCTION__);
        /* r = gmel_bf_prior_m2 (&priors[i], ps[i].p, ps[i].a, ps[i].c, ps[i].g); */
        break;
     case 3:
        psi_fatal ("no implementation: %s, %d, %s\n", 
                    __FILE__, __LINE__, __PRETTY_FUNCTION__);
        /* r = gmel_bf_prior_m3 (&priors[i], ps[i].s, ps[i].a, ps[i].c, ps[i].g); */
        break;
     case 4:
        r = gmel_bf_prior_m4 (&priors[i], 
                              ps[i].s, ps[i].p, ps[i].a, ps[i].c, ps[i].g); 
        break;
     default:
        psi_fatal ("no implementation: %s, %d, %s\n", 
                    __FILE__, __LINE__, __PRETTY_FUNCTION__);
        assert(0);
     }
     assert (isfinite(priors[i]));
  }

  *prior = priors[0] - priors[1];
  return r;
}

int 
gmel_bf_posterior (double *posterior, int model_top, int model_bot,
                   parameter p_top, parameter p_bot)
{
  parameter ps[2];
  ps[0] = p_top;
  ps[1] = p_bot;
  int models[2] = { model_top, model_bot };
  double posteriors[2] = { 0, 0 };
  int i;
  int r = EXIT_SUCCESS;

  for (i = 0; i < 2; i++) { 
     switch (models[i]) {
     case 1:
        r = gmel_bf_posterior_m1 (&posteriors[i], 
                                  ps[i].a, ps[i].c, ps[i].g); 
        break;
     case 2:
        psi_fatal ("no implementation: %s, %d, %s\n", 
                    __FILE__, __LINE__, __PRETTY_FUNCTION__);
/*
        r = gmel_bf_posterior_m2 (&posteriors[i], 
                                  ps[i].p, ps[i].a, ps[i].c, ps[i].g); 
*/
        break;
     case 3:
        psi_fatal ("no implementation: %s, %d, %s\n", 
                    __FILE__, __LINE__, __PRETTY_FUNCTION__);
/*
        r = gmel_bf_posterior_m3 (&posteriors[i], 
                                  ps[i].s, ps[i].a, ps[i].c, ps[i].g); 
*/
        break;
     case 4:
        r = gmel_bf_posterior_m4 (&posteriors[i], 
                                  ps[i].s, ps[i].p, ps[i].a, ps[i].c, ps[i].g); 
        break;
     default:
        assert(0);
        psi_fatal ("no implementation: %s, %d, %s\n", 
                    __FILE__, __LINE__, __PRETTY_FUNCTION__);
     }
     assert (isfinite(posteriors[i]));
  }

  *posterior = posteriors[0] - posteriors[1];
  assert (isfinite(*posterior));

  return r;
}

/******************************************************************************
 SAMPLING Functions
******************************************************************************/


/******************************************************************************
 MODEL 1 Functions
******************************************************************************/

/* use this function only when comparing model 1 and model 3 */
/* WRONG! This can be used as it is */
/* if s=p=0, then we no longer need Gibbs Sampler */
int 
gmel_bf_likelihood_m1 (double *likelihood, 
                       double a, double c, double g)
{
  int r = EXIT_SUCCESS;
  double t = 1.0 - a - c - g;
  psi_err_base (a, c, g, t);
  *likelihood = psi_mcmc_n_likelihood (a, c, g, t);
  assert (isfinite(*likelihood));
  if (isfinite(*likelihood) == 0) {
     psi_fatal ("model1: likelihood is not finite"); 
  }

  return r;
}

int 
gmel_bf_prior_m1 (double *prior, 
                  double a, double c, double g)
{
  int r = EXIT_SUCCESS;
  *prior = 0;
  if (want_use_Dirichlet_prior_n == 1) {
     /* this should not be zero because Dirichlet
        density would be a constant. 
        This is cancelled out with the other prior */
  } else {
     assert (0);
     psi_fatal ("no implementation: %s:%d:%s"
                 __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }
  
  a = c = g = 0.0; /* for removing warning */
  return r;
}

int 
gmel_bf_posterior_m1 (double *log_val, 
                      double a, double c, double g)
{
  int r = EXIT_SUCCESS;
  double val_block1 = 0.0;
  r = gmel_bf_posterior_m1_block1 (&val_block1, a, c, g);
  assert (isfinite(val_block1));
  *log_val = val_block1;

  return r;
}

int 
gmel_bf_posterior_m1_block1 (double *log_val, 
                             double a, double c, double g)
{
  int r = EXIT_SUCCESS;
  double val_e1 = 0.0;
  double val_e2 = 0.0;

  r = gmel_bf_posterior_m1_block1_e1 (&val_e1, a, c, g);
  r = gmel_bf_posterior_m1_block1_e2 (&val_e2, a, c, g);
  *log_val = val_e1 - val_e2;
  assert (isfinite(*log_val));

  return r;
}

int 
gmel_bf_posterior_m1_block1_e1
   (double *log_val, double a, double c, double g)
{
  int i;
  double ag, cg, gg, tg, sg, pg;
  int gen;
  int w = 0; 
  int r = EXIT_SUCCESS;
  double log_e1 = 0;
  double log_alpha, q;
  double *A_i = NULL;
  FILE *bf = NULL;

  bf = fopen (pname, "r");
  assert (bf != NULL);
  if (bf == NULL) {
     psi_fatal ("%s:%d:%s; %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
  }

  locate_m1_block1 (bf, sample_size);
  A_i = XMALLOC (double, sample_size);
  for (i = 0; i < sample_size; i++) {
     r = fscanf (bf, "%d%lf%lf%lf%lf%lf%lf\n", /* BUGGY: "\n" at the end */
                 &gen, &ag, &cg, &gg, &tg, &sg, &pg);
     assert (r != EOF && r == 7);
     if (r == EOF || r != 7) {
        psi_fatal ("%s:%s:%d:%s; %s %d %lf %lf %lf %lf %lf %lf - r: %d\n",
                 program_name,
                 __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                 "error of reading posterior draws of block 1 in model 1",
                 gen, ag, cg, gg, tg, sg, pg, r);
     }
     /* calculate each term of the sum */
     r = gmel_bf_posterior_m1_block1_alpha (&log_alpha, ag, cg, gg, a, c, g, w);
     r = gmel_bf_posterior_m1_block1_q (&q, ag, cg, gg, a, c, g, w);
     A_i[i] = log_alpha + q;
  }
  log_e1 = psi_sum_exp (A_i, sample_size);

  if (isfinite(log_e1) == 0) {
     psi_print_double_array (A_i, sample_size);
  }
  assert (isfinite(log_e1));
  psi_err_finite (log_e1);
  XFREE (A_i); 

  fclose (bf);
  bf = NULL;

  *log_val = log_e1 - log_sample_size;
  assert (isfinite(*log_val));
  r = gmel_err_same (program_name, 
                     __FILE__, __LINE__, __PRETTY_FUNCTION__,
                     *log_val, log_e1);
  return r;
}

int 
gmel_bf_posterior_m1_block1_e2 (double *log_val, double a, double c, double g)
{
  int i;
  int r = EXIT_SUCCESS;
  double aj, cj, gj;
  double log_e2 = 0;
  double log_alpha;
  double *A_i = NULL;
  int w;
  parameter q_param;
  FILE *qfile = NULL;

  if (want_save_q == 1) {
     qfile = fopen (qname, "a");
  } else if (want_load_q == 1) {
     qfile = fopen (qname, "r");
  }
  if ((want_save_q == 1 || want_load_q == 1) && qfile == NULL) {
     psi_fatal ("%s:%s:%d:%s; %s\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
  } 
 
  A_i = XMALLOC (double, sample_size);
  /* calculate each term of the sum */
  for (i = 0; i < sample_size; i++) { /* FIX */
     if (want_load_q == 1) {
        gmel_bf_load_q (qfile, &q_param, i, MODEL1_BLOCK1);
        aj = q_param.a; cj = q_param.c; gj = q_param.g;
     } else {
        gmel_bf_n_sample(a, c, g, &aj, &cj, &gj, &w);
        if (want_save_q == 1) {
           psi_st_parameter (&q_param, 0, 0, aj, cj, gj); 
           if (i == 0) {
              write_parameter (qfile, q_param, GEN_INIT);
           } else {
              write_parameter (qfile, q_param, GEN_NEXT);
           } 
        } 
     }
     assert (a != aj && c != cj && g != gj);
     gmel_bf_posterior_m1_block1_alpha (&log_alpha, a, c, g, 
                                        aj, cj, gj, w);
     A_i[i] = log_alpha;
  }

  log_e2 = psi_sum_exp (A_i, sample_size);
  XFREE (A_i); 
  fclose (qfile);
  qfile = NULL;
  assert (isfinite(log_e2));
  psi_err_finite (log_e2);
  
  *log_val = log_e2 - log_sample_size;
  assert (isfinite(*log_val));
  r = gmel_err_same (program_name, 
                     __FILE__, __LINE__, __PRETTY_FUNCTION__,
                     *log_val, log_e2);

  return r;
}

int 
gmel_bf_posterior_m1_block1_alpha (double *log_alpha, 
                                   double ag, double cg, double gg,
                                   double a, double c, double g, int w)
{
  int r = EXIT_SUCCESS;
  double likelihood_top = 0;
  double likelihood_bot = 0;
  double prior_top = 0;
  double prior_bot = 0;
  double jumping_top = 0;
  double jumping_bot = 0;

  r = gmel_bf_likelihood_m1 (&likelihood_top, a, c, g);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(likelihood_top));
  r = gmel_bf_likelihood_m1 (&likelihood_bot, ag, cg, gg);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(likelihood_bot));
  r = gmel_bf_prior_m1 (&prior_top, a, c, g);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(prior_top));
  r = gmel_bf_prior_m1 (&prior_bot, ag, cg, gg);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(prior_bot));
  r = gmel_bf_posterior_m1_block1_q (&jumping_top, a, c, g, ag, cg, gg, w);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(jumping_top));
  r = gmel_bf_posterior_m1_block1_q (&jumping_bot, ag, cg, gg, a, c, g, w);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(jumping_bot));

  *log_alpha = likelihood_top - likelihood_bot
               + prior_top - prior_bot
               + jumping_top - jumping_bot; 
  assert (isfinite(*log_alpha));
  if (isfinite(*log_alpha) == 0) {
     psi_fatal ("log_alpha is not finite: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  if (*log_alpha > 0) { /* min (something, 1) */
     *log_alpha = 0; 
  }
  return r;
}

int 
gmel_bf_posterior_m1_block1_q
   (double *q, double a1, double c1, double g1, 
                   double a2, double c2, double g2, int w)
{
  /* r = gmel_bf_posterior_m3_block2_q(q, a1, c1, g1, a2, c2, g2, 0, w); */
  int r = EXIT_SUCCESS;
  double t1 = 1.0 - a1 - c1 - g1;
  double t2 = 1.0 - a2 - c2 - g2;
  assert(a1 != 0);
  assert(c1 != 0);
  assert(g1 != 0);
  assert(t1 != 0);
  assert(a2 != 0);
  assert(c2 != 0);
  assert(g2 != 0);
  assert(t2 != 0);
  double p_alpha[4] = {a1, c1, g1, t1};
  double p_theta_star[4] = {a2, c2, g2, t2};

  if (want_use_Dirichlet_proposal_n == 1) {
     *q = psi_mcmc_n_dirichlet_proposal_lnpdf (p_alpha, p_theta_star);
     if (want_debug_bf == 1) {
        fprintf (ofile, "Dirichlet Proposal Density: %lf, %lf, %lf, %lf -> %lf, %lf, %lf, %lf: %lf\n",
                 a1, c1, g1, t1, a2, c2, g2, t2, *q);
     }
  } else {
     psi_fatal ("no implementation of not Dirichlet: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  assert (r == EXIT_SUCCESS);
  assert (isfinite(*q));
  if (isfinite(*q) == 0) {
     psi_fatal ("not finite: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  w = 0; /* just for removing warning */
  return r;
}

/******************************************************************************
 MODEL 4 Functions
******************************************************************************/

int 
gmel_bf_likelihood_m4 (double *likelihood, 
                       double s, double p, double a, double c, double g)
{
  // goshng: Model 4's Likelihood
  // REPLACED BY GIBBS SAMPLER
  // THIS IS NOT A FUNCTION, DO NOT CALL IT
  assert (0);
  psi_fatal ("no implementation of likelihood model 4: %s %d %s",
             __FILE__, __LINE__, __PRETTY_FUNCTION__);
  *likelihood = s = p = a = c = g = 0.0; /* just for removing warning */
  return EXIT_FAILURE;
}

int 
gmel_bf_prior_m4 (double *prior, 
                  double s, double p, double a, double c, double g)
{
  int r = EXIT_SUCCESS;
  double v = 0; 
  *prior = 0;
  if (want_use_Dirichlet_prior_n == 1) {
     /* this should not be zero because Dirichlet
        density would be a constant. 
        This is cancelled out with the other prior */
  } else {
     psi_fatal ("no implementation of not Dirichlet: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }
  
  if (want_use_normal_prior_s == 1) {
     v = psi_mcmc_s_normal_prior_lnpdf (s); 
     *prior += v;
  } else {
     v = psi_mcmc_s_flat_prior_lnpdf (s); 
     /* v = - log_flatprior_s; */
     /* v = - log(flatprior_s_end - flatprior_s_begin); */
     *prior += v;
  }
 
  if (want_use_normal_prior_p == 1) {
     v = psi_mcmc_p_normal_prior_lnpdf (p); 
     /* v = psi_gaussian_lnpdf (0, prior_normal_sd_p, p); */
     *prior += v;
  } else {
     v = psi_mcmc_p_flat_prior_lnpdf (p); 
     /* v = - log_flatprior_p; */
     /* v = - log(flatprior_p_end - flatprior_p_begin); */
     *prior += v;
  }

  assert (isfinite(*prior));
  psi_err_finite (*prior);

  a = c = g = 0.0; /* just for removing warning */
  return r;
}

int 
gmel_bf_posterior_m4 (double *log_val, 
                                double s, double p, 
                                double a, double c, double g)
{
  int r;
  double val_block1 = 0; 
  double val_block2 = 0; 
  double val_block3 = 0; 
  
  r = gmel_bf_posterior_m4_block1(&val_block1, s);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(val_block1));

  r = gmel_bf_posterior_m4_block2(&val_block2, s, p);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(val_block2));

  r = gmel_bf_posterior_m4_block3(&val_block3, s, p, a, c, g);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(val_block3));

  if (want_debug_bf == 1) {
     fprintf (ofile, "%s:%s:%d:%s; %s block1=%lf, block2=%lf, block3=%lf\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "log_val is not isfinite", val_block1, val_block2, val_block3);
  } 

  *log_val = val_block1 + val_block2 + val_block3;
  assert (isfinite(*log_val));
  if (isfinite(*log_val) == 0) {
     psi_fatal ("%s:%s:%d:%s; %s block1=%lf, block2=%lf, block3=%lf\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "log_val is not isfinite", val_block1, val_block2, val_block3);
  } 

  return r;
}

int 
gmel_bf_posterior_m4_block1 (double *log_val, double s)
{
  int r = EXIT_SUCCESS;
  double val_e1 = 0;
  double val_e2 = 0;
  r = gmel_bf_posterior_m4_block1_e1(&val_e1, s);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(val_e1));

  if (want_debug_bf == 1) {
     fprintf (ofile, "%s:%s:%d:%s; %s: %lf\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__,
              "block1_e1", val_e1);
  }

  r = gmel_bf_posterior_m4_block1_e2(&val_e2, s);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(val_e2));

  if (want_debug_bf == 1) {
     fprintf (ofile, "%s:%s:%d:%s; s val_e1=%lf, val_e2=%lf\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, 
              val_e1, val_e2);
  }

  *log_val = val_e1 - val_e2;
  assert (isfinite(*log_val));
  if (isfinite(*log_val) == 0) {
     psi_fatal ("%s:%s:%d:%s; %s val_e1=%lf, val_e2=%lf\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "log_val is not isfinite", val_e1, val_e2);
  } 

  return r;
}

int 
gmel_bf_posterior_m4_block2
   (double *log_val, double s, double p)
{
  int r = EXIT_SUCCESS;
  double val_e1 = 0;
  double val_e2 = 0;
  r = gmel_bf_posterior_m4_block2_e1(&val_e1, s, p);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(val_e1));

  r = gmel_bf_posterior_m4_block2_e2(&val_e2, s, p);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(val_e2));

  if (want_debug_bf == 1) {
     fprintf (ofile, "%s:%s:%d:%s; s val_e1=%lf, val_e2=%lf\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, 
              val_e1, val_e2);
  }

  *log_val = val_e1 - val_e2;
  assert (isfinite(*log_val));
  if (isfinite(*log_val) == 0) {
     fprintf (stderr, "%s:%s:%d:%s; %s val_e1=%lf, val_e2=%lf\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "log_val is not isfinite", val_e1, val_e2);
     return ERR_BF_NOTFINITE;
  } 

  return r;
}

int 
gmel_bf_posterior_m4_block3
   (double *log_val, 
              double s, double p, double a, double c, double g)
{
  int r = EXIT_SUCCESS;
  double val_e1 = 0;
  double val_e2 = 0;
  r = gmel_bf_posterior_m4_block3_e1(&val_e1, s, p, a, c, g);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(val_e1));

  r = gmel_bf_posterior_m4_block3_e2(&val_e2, s, p, a, c, g);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(val_e2));

  if (want_debug_bf == 1) {
     fprintf (ofile, "%s:%s:%d:%s; val_e1=%lf, val_e2=%lf\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, 
              val_e1, val_e2);
  }

  *log_val = val_e1 - val_e2;
  assert (isfinite(*log_val));
  if (isfinite(*log_val) == 0) {
     fprintf (stderr, "%s:%s:%d:%s; %s val_e1=%lf, val_e2=%lf\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "log_val is not isfinite", val_e1, val_e2);
     return ERR_BF_NOTFINITE;
  } 

  return r;
}

int 
gmel_bf_posterior_m4_block1_e1 (double *log_val, double s)
{
  int i;
  int gen;
  int r = EXIT_SUCCESS;
  double ag, cg, gg, tg, sg, pg;
  double log_e1 = 0;
  double log_alpha, q;
  double *A_i = NULL;
  FILE *bf = NULL;

  bf = fopen (pname, "r");
  if (bf == NULL) {
     psi_fatal ("%s:%s:%d:%s; %s",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
  }
  assert( bf != NULL );
  
  locate_m4_block1_e1 (bf, sample_size);
  A_i = XMALLOC (double, sample_size);
  int size = 0; 
  for (i = 0; i < sample_size; i++) {
     r = fscanf (bf, "%d%lf%lf%lf%lf%lf%lf",  /* BUGGY: "\n" */
                       &gen, &ag, &cg, &gg, &tg, &sg, &pg);
     assert (r != EOF && r == 7);
     if (r == EOF || r != 7) {
        psi_fatal ("%s:%s:%d:%s; %s %d %lf %lf %lf %lf %lf %lf - r: %d\n",
                 program_name,
                 __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                 "error of reading posterior draws of block 1 in model 4",
                 gen, ag, cg, gg, tg, sg, pg, r);
     }
 
     /* calculate each term of the sum */
     r = gmel_bf_posterior_m4_block1_alpha (&log_alpha, sg, s, pg, ag, cg, gg);
     if (r == ERR_BF_UNDEFINED_Q) { 
        continue;
     }
     r = gmel_bf_posterior_m4_block1_q(&q, sg, s, pg, ag, cg, gg);
     if (r == ERR_BF_UNDEFINED_Q) { 
        continue;
     }
     A_i[size] = log_alpha + q;
     size++;
  }
  fclose (bf);
  bf = NULL;

  fprintf (ofile, "BF B1:\t%d\t%d\t%d\n",
           size, limit_sample_size, sample_size);
  if (size < limit_sample_size) {
     fprintf (ofile, "BF Fatal B1: small sample size\t%d\t%d\t%d\n", 
              size, limit_sample_size, sample_size);
     psi_fatal ("sample size is too small for calculating E1 of block 1: \
                 size - %d, limit - %d, sample_size - %d", 
                size, limit_sample_size, sample_size);
  }

  if (size < sample_size) {
     fprintf (ofile, "BF info: size is smaller than sample size for calculating E1 of block 1: \
                 size - %d, limit - %d, sample_size - %d\n", 
                size, limit_sample_size, sample_size);
  }
 
  if (want_debug_bf == 1) {
     psi_print_double_array (A_i, size);
  }
  log_e1 = psi_sum_exp (A_i, size);

  XFREE (A_i); 
  assert (isfinite(log_e1));
  r = psi_err_finite (log_e1);
  if (r != EXIT_SUCCESS) {
     psi_fatal ("not finite log_e1: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  *log_val = log_e1 - log_sample_size;
  assert (isfinite(*log_val));
  r = gmel_err_same (program_name, 
                     __FILE__, __LINE__, __PRETTY_FUNCTION__,
                     *log_val, log_e1);
  return r;
}

int 
gmel_bf_posterior_m4_block1_e2 (double *log_val, double s)
{
  int i;
  int r = EXIT_SUCCESS;
  int gen;
  double aj, cj, gj, tj, sj, pj;
  double log_e2 = 0;
  double log_alpha;
  double *A_i = NULL;
  FILE *bf = NULL;
  parameter q_param;
  FILE *qfile = NULL;

  bf = fopen (pname, "r");
  assert( bf != NULL );
  if (bf == NULL) {
     psi_fatal ("%s:%s:%d:%s; %s\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
  }

  if (want_save_q == 1) {
     qfile = fopen (qname, "a");
  } else if (want_load_q == 1) {
     qfile = fopen (qname, "r");
  }
  if ((want_save_q == 1 || want_load_q == 1) && qfile == NULL) {
     psi_fatal ("%s:%s:%d:%s; %s\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
     fclose (bf);
     bf = NULL;
  }

  locate_m4_block1_e2 (bf, sample_size);
  A_i = XMALLOC (double, sample_size);
  for (i = 0; i < sample_size; i++) {
     r = fscanf(bf, "%d%lf%lf%lf%lf%lf%lf", 
                       &gen, &aj, &cj, &gj, &tj, &sj, &pj);
     assert (r != EOF && r == 7);
     if (r == EOF || r != 7) {
        psi_fatal ("%s:%s:%d:%s; %s\n",
                 program_name,
                 __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                 "error of reading posterior draws of block 2 in model 4");
     }

     /* calculate each term of the sum */
     if (want_load_q == 1) {
        r = gmel_bf_load_q (qfile, &q_param, i, MODEL4_BLOCK1);
        if (r != EXIT_SUCCESS) {
           psi_fatal ("%s:%s:%d:%s; %s\n",
                    program_name,
                    __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                    "error in bootstrap from loaded q");
        } else {
           sj = q_param.s; pj = q_param.p;
           aj = q_param.a; cj = q_param.c; gj = q_param.g;
        }
     } else {
        r = gmel_bf_s_sample(s, &sj);
        if (r != EXIT_SUCCESS) {
           psi_fatal ("%s:%s:%d:%s; %s\n",
                    program_name,
                    __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                    "error of s sampling of posterior draws of block 2 in model 4");
        }
        if (want_save_q == 1) {
           psi_st_parameter (&q_param, sj, pj, aj, cj, gj); 
           if (i == 0) {
              write_parameter (qfile, q_param, GEN_INIT);
           } else {
              write_parameter (qfile, q_param, GEN_NEXT);
           } 
        } 
     } 
     r = gmel_bf_posterior_m4_block1_alpha (&log_alpha, s, sj, pj, aj, cj, gj);
     A_i[i] = log_alpha;
  }
  fclose (qfile);
  qfile = NULL;
  fclose (bf);
  bf = NULL;
  if (want_debug_bf == 1) {
     psi_print_double_array (A_i, sample_size);
  }
  log_e2 = psi_sum_exp (A_i, sample_size);
  XFREE (A_i); 
  assert (isfinite(log_e2));
  r = psi_err_finite (log_e2);
  if (r != EXIT_SUCCESS) {
     psi_fatal ("not finite log_e1: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  *log_val = log_e2 - log_sample_size;
  assert (isfinite(*log_val));
  r = gmel_err_same (program_name, 
                     __FILE__, __LINE__, __PRETTY_FUNCTION__,
                     *log_val, log_e2);

  return r;
}

int 
gmel_bf_posterior_m4_block1_alpha (double *log_alpha, 
                                   double sg, double s, 
                                   double pg, double ag, double cg, double gg)
{
  int r = EXIT_SUCCESS;
  double prior_top = 0;
  double prior_bot = 0;
  double jumping_top = 0;
  double jumping_bot = 0;
  double likelihood = 0;

  /* Changed by Gibbs Sampler */
  parameter top, bot;
  psi_st_parameter(&top, s, pg, ag, cg, gg);
  psi_st_parameter(&bot, sg, pg, ag, cg, gg);

  /* for efficiency */
  r = gmel_bf_posterior_m4_block1_q(&jumping_top, s, sg, pg, ag, cg, gg);
  if (r == ERR_BF_UNDEFINED_Q) { 
     return ERR_BF_UNDEFINED_Q;
  }
  r = gmel_bf_posterior_m4_block1_q(&jumping_bot, sg, s, pg, ag, cg, gg);
  if (r == ERR_BF_UNDEFINED_Q) { 
     return ERR_BF_UNDEFINED_Q;
  }

  if (top.s == bot.s && top.p == bot.p 
      && top.a == bot.a && top.c == bot.c && top.g == bot.g) {
     likelihood = 0;
  } else {
     r = gmel_bf_likelihood_gibbs(&likelihood, top, bot);
     if (r != EXIT_SUCCESS) {
        fprintf (stderr, "%s:%s:%d:%s; %s\n",
                 program_name,
                 __FILE__, __LINE__, __PRETTY_FUNCTION__,
                 "error in likelihood of alpha of block 1 in model 4");
        return r;
     }
     assert (r == EXIT_SUCCESS);
     assert (isfinite(likelihood));
  }

  r = gmel_bf_prior_m4(&prior_top, s, pg, ag, cg, gg);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(prior_top));

  r = gmel_bf_prior_m4(&prior_bot, sg, pg, ag, cg, gg);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(prior_bot));


  *log_alpha = likelihood
               + prior_top - prior_bot
               + jumping_top - jumping_bot; 
  assert (isfinite(*log_alpha));
  if (isfinite(*log_alpha) == 0) {
     psi_fatal ("log_alpah is not finite: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }
  if (*log_alpha > 0) {
     *log_alpha = 0; 
  }

  return r;
}

int 
gmel_bf_posterior_m4_block1_q (double *q, double sg, double s, 
                               double pg, double ag, double cg, double gg)
{
  int r = EXIT_SUCCESS;
  double v;
  

  if (want_use_normal_proposal_s == 1) 
    {
      *q = psi_mcmc_s_normal_proposal_lnpdf (s, sg); /* first is x and second is mu */
    } 
  else 
    { 
      v = psi_mcmc_s_flat_proposal_pdf (s, sg);
      if (v > 0.0)
        {
          *q = log (v); /* minus */
        }
      else
        {
          *q = 0.0;
          return ERR_BF_UNDEFINED_Q;
        }
    }
  assert (isfinite(*q));
  r = psi_err_finite (*q);

  pg = ag = cg = gg = 0.0; /* just for removing warning */
  return r;
}

int 
gmel_bf_posterior_m4_block2_e1 (double *log_val, double s, double p)
{
  int i;
  int gen;
  int r = EXIT_SUCCESS;
  double aj, cj, gj, tj, sj, pj;
  double log_e1 = 0;
  double log_alpha, q;
  double *A_i = NULL;
  FILE *bf = NULL;
  bf = fopen (pname, "r");
  assert( bf != NULL );
  if (bf == NULL) {
     psi_fatal ("%s:%s:%d:%s; %s\n",
                program_name,
                __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
  }
  
  locate_m4_block2_e1 (bf, sample_size);
  A_i = XMALLOC (double, sample_size);
  int size = 0;
  for (i = 0; i < sample_size; i++) {
     r = fscanf(bf, "%d%lf%lf%lf%lf%lf%lf",  /* BUGGY: "\n" */
                       &gen, &aj, &cj, &gj, &tj, &sj, &pj);
     assert (r != EOF && r == 7);
     if (r == EOF || r != 7) {
        psi_fatal ("%s:%s:%d:%s; %s\n",
                 program_name,
                 __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                 "error of reading posterior draws of block 1 in model 4");
     }
     /* calculate each term of the sum */
     r = gmel_bf_posterior_m4_block2_alpha(&log_alpha, pj, p, s, aj, cj, gj);
     if (r == ERR_BF_UNDEFINED_Q) { 
        continue;
     }
     r = gmel_bf_posterior_m4_block2_q(&q, pj, p, s, aj, cj, gj);
     if (r == ERR_BF_UNDEFINED_Q) { 
        continue;
     }
     A_i[size] = log_alpha + q;
     size++;
  }
  fclose (bf);
  bf = NULL;

  fprintf (ofile, "BF B2:\t%d\t%d\t%d\n",
           size, limit_sample_size, sample_size);
  if (size < limit_sample_size) {
     fprintf (ofile, "BF Fatal B2: small sample size\t%d\t%d\t%d\n", 
              size, limit_sample_size, sample_size);
     psi_fatal ("sample size is too small for calculating E1 of block 2: \
                 size - %d, limit - %d, sample_size - %d", 
                size, limit_sample_size, sample_size);
  }

  if (size < sample_size) {
     fprintf (ofile, "BF info: size is smaller than sample size for calculating E1 of block 2: \
                 size - %d, limit - %d, sample_size - %d\n", 
                size, limit_sample_size, sample_size);
  }
 
  if (want_debug_bf == 1) {
     psi_print_double_array (A_i, size);
  }
  log_e1 = psi_sum_exp (A_i, size);

  XFREE (A_i); 
  assert (isfinite(log_e1));
  r = psi_err_finite (log_e1);
  if (r != EXIT_SUCCESS) {
     psi_fatal ("log_alpah is not finite: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  *log_val = log_e1 - log_sample_size;
  assert (isfinite(*log_val));
  r = gmel_err_same (program_name, 
                     __FILE__, __LINE__, __PRETTY_FUNCTION__,
                     *log_val, log_e1);
  return r;
}

int 
gmel_bf_posterior_m4_block2_e2
   (double *log_val, double s, double p)
{
  int i;
  int gen;
  int r = EXIT_SUCCESS;
  double ak, ck, gk, tk, sk, pk;
  double log_e2 = 0;
  double log_alpha;
  double *A_i = NULL;
  FILE *bf = NULL;
  FILE *qfile = NULL;
  parameter q_param;

  bf = fopen (pname, "r");
  assert( bf != NULL );
  if (bf == NULL) {
     psi_fatal ("%s:%s:%d:%s; %s",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
  }

  if (want_save_q == 1) {
     qfile = fopen (qname, "a");
  } else if (want_load_q == 1) {
     qfile = fopen (qname, "r");
  }
  if ((want_save_q == 1 || want_load_q == 1) && qfile == NULL) {
     psi_fatal ("%s:%s:%d:%s; %s\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
  }

  locate_m4_block2_e2 (bf, sample_size);
  A_i = XMALLOC (double, sample_size);
  for (i = 0; i < sample_size; i++) {
     r = fscanf (bf, "%d%lf%lf%lf%lf%lf%lf", 
                       &gen, &ak, &ck, &gk, &tk, &sk, &pk);
     assert (r != EOF && r == 7);
     if (r == EOF || r != 7) {
        psi_fatal ("%s:%s:%d:%s; %s",
                 program_name,
                 __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                 "error of reading posterior draws of block 3 in model 4");
     }

     /* calculate each term of the sum */
     if (want_load_q == 1) {
        r = gmel_bf_load_q (qfile, &q_param, i, MODEL4_BLOCK2);
        if (r != EXIT_SUCCESS) {
           psi_fatal ("%s:%s:%d:%s; %s",
                    program_name,
                    __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                    "error in bootstrap from loaded q");
        } else {
           sk = q_param.s; pk = q_param.p;
           ak = q_param.a; ck = q_param.c; gk = q_param.g;
        }
     } else {
        r = gmel_bf_p_sample(p, &pk);
        if (want_save_q == 1) {
           psi_st_parameter (&q_param, sk, pk, ak, ck, gk); 
           if (i == 0) {
              write_parameter (qfile, q_param, GEN_INIT);
           } else {
              write_parameter (qfile, q_param, GEN_NEXT);
           } 
        } 
     } 
     r = gmel_bf_posterior_m4_block2_alpha(&log_alpha, p, pk, s, ak, ck, gk);
     A_i[i] = log_alpha;
  }
  fclose (qfile);
  qfile = NULL;
  fclose (bf);
  bf = NULL;
  if (want_debug_bf == 1) {
     psi_print_double_array (A_i, sample_size);
  }
  log_e2 = psi_sum_exp (A_i, sample_size);
  XFREE (A_i); 
  assert (isfinite(log_e2));
  r = psi_err_finite (log_e2);
  if (r != EXIT_SUCCESS) {
     psi_fatal ("log_alpah is not finite: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  *log_val = log_e2 - log_sample_size;
  assert (isfinite(*log_val));
  r = gmel_err_same (program_name, 
                     __FILE__, __LINE__, __PRETTY_FUNCTION__,
                     *log_val, log_e2);

  return r;
}

int 
gmel_bf_posterior_m4_block2_alpha (double *log_alpha, 
                                   double pj, double p, double s, 
                                   double aj, double cj, double gj)
{
  int r = EXIT_SUCCESS;
  double prior_top = 0;
  double prior_bot = 0;
  double jumping_top = 0;
  double jumping_bot = 0;
  double likelihood = 0;

  // Changed by Gibbs Sampler
  parameter top, bot;
  psi_st_parameter(&top, s, p, aj, cj, gj);
  psi_st_parameter(&bot, s, pj, aj, cj, gj);

  /* for efficiency */
  r = gmel_bf_posterior_m4_block2_q(&jumping_top, p, pj, s, aj, cj, gj);
  if (r == ERR_BF_UNDEFINED_Q) {
     return ERR_BF_UNDEFINED_Q;
  }
  r = gmel_bf_posterior_m4_block2_q(&jumping_top, pj, p, s, aj, cj, gj);
  if (r == ERR_BF_UNDEFINED_Q) {
     return ERR_BF_UNDEFINED_Q;
  }

  if (top.s == bot.s && top.p == bot.p 
      && top.a == bot.a && top.c == bot.c && top.g == bot.g) {
     likelihood = 0;
  } else {
     r = gmel_bf_likelihood_gibbs (&likelihood, top, bot);
     assert (isfinite(likelihood));
  }

  r = gmel_bf_prior_m4 (&prior_top, s, p, aj, cj, gj);
  assert (isfinite(prior_top));

  r = gmel_bf_prior_m4(&prior_bot, s, pj, aj, cj, gj);
  assert (isfinite(prior_bot));

  *log_alpha = likelihood
               + prior_top - prior_bot
               + jumping_top - jumping_bot; 
  assert (isfinite(*log_alpha));
  if (isfinite(*log_alpha) == 0) {
     psi_fatal ("log_alpah is not finite: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }
  if (*log_alpha > 0) {
     *log_alpha = 0; 
  }

  return r;
}

int 
gmel_bf_posterior_m4_block2_q (double *q, double pj, double p,
                               double s, double aj, double cj, double gj)
{
  int r = EXIT_SUCCESS;
  double v;

  if (want_use_normal_proposal_p == 1) 
    {
      *q = psi_mcmc_p_normal_proposal_lnpdf (p, pj); /* first is x and second is mu */
    } 
  else 
    { 
      v = psi_mcmc_p_flat_proposal_pdf (p, pj);
      if (v > 0.0)
        {
          *q = log (v); /* minus */
        }
      else
        {
          *q = 0.0;
          return ERR_BF_UNDEFINED_Q;
        }
    }
  assert (isfinite(*q));
  r = psi_err_finite (*q);

  s = aj = cj = gj = 0.0; /* just for removing warning */
  return r;
}

int 
gmel_bf_posterior_m4_block3_e1 (double *log_val, 
                                double s, double p, 
                                double a, double c, double g)
{
  int i;
  int r = EXIT_SUCCESS;
  int gen;
  int w = 0;
  double ak, ck, gk, tk, sk, pk;
  double log_e1 = 0;
  double q = 0;
  double log_alpha = 0;
  double *A_i = NULL;
  FILE *bf;
  bf = fopen (pname, "r");
  assert( bf != NULL );
  if (bf == NULL) {
     psi_fatal ("%s:%s:%d:%s; %s\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
  }
  
  locate_m4_block3_e1 (bf, sample_size);
  A_i = XMALLOC (double, sample_size);
  for (i = 0; i < sample_size; i++) {
     r = fscanf(bf, "%d%lf%lf%lf%lf%lf%lf", 
                       &gen, &ak, &ck, &gk, &tk, &sk, &pk);
     assert (r != EOF && r == 7);
     if (r == EOF || r != 7) {
        psi_fatal ("%s:%s:%d:%s; %s\n",
                 program_name,
                 __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                 "error of reading posterior draws of block 3 in model 4");
     }
     /* calculate each term of the sum */
     r = gmel_bf_posterior_m4_block3_alpha (&log_alpha, ak, ck, gk,
                                            a, c, g, s, p, w);
     r = gmel_bf_posterior_m4_block3_q (&q, ak, ck, gk, 
                                           a, c, g, s, p, w);
     A_i[i] = log_alpha + q;
  }
  fclose (bf);
  bf = NULL;
  if (want_debug_bf == 1) {
     psi_print_double_array (A_i, sample_size);
  }
  log_e1 = psi_sum_exp (A_i, sample_size);
  XFREE (A_i); 
  assert (isfinite(log_e1));
  r = psi_err_finite (log_e1);
  if (r != EXIT_SUCCESS) {
     psi_fatal ("log_e1 is not finite: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  *log_val = log_e1 - log_sample_size;
  assert (isfinite(*log_val));
  r = gmel_err_same (program_name, 
                     __FILE__, __LINE__, __PRETTY_FUNCTION__,
                     *log_val, log_e1);
 
  return r;
}

int 
gmel_bf_posterior_m4_block3_e2 (double *log_val, 
                                double s, double p, 
                                double a, double c, double g)
{
  int i;
  int r = EXIT_SUCCESS;
  double al, cl, gl;
  double log_e2 = 0;
  double log_alpha = 0;
  double *A_i = NULL;
  int w;
  parameter q_param;
  FILE *qfile = NULL;

  if (want_save_q == 1) {
     qfile = fopen (qname, "a");
  } else if (want_load_q == 1) {
     qfile = fopen (qname, "r");
  }
  if ((want_save_q == 1 || want_load_q == 1) && qfile == NULL) {
     psi_fatal ("%s:%s:%d:%s; %s\n",
              program_name,
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno));
  }

  A_i = XMALLOC (double, sample_size);
  for (i = 0; i < sample_size; i++) {
     if (want_load_q == 1) {
        r = gmel_bf_load_q (qfile, &q_param, i, MODEL4_BLOCK3);
        if (r != EXIT_SUCCESS) {
           psi_fatal ("%s:%s:%d:%s; %s\n",
                    program_name,
                    __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                    "error in bootstrap from loaded q");
        } else {
           al = q_param.a; cl = q_param.c; gl = q_param.g;
        }
     } else {
        r = gmel_bf_n_sample(a, c, g, &al, &cl, &gl, &w);
        if (want_save_q == 1) {
           psi_st_parameter (&q_param, s, p, al, cl, gl); 
           if (i == 0) {
              write_parameter (qfile, q_param, GEN_INIT);
           } else {
              write_parameter (qfile, q_param, GEN_NEXT);
           } 
        } 
     } 
 
     r = gmel_bf_posterior_m4_block3_alpha (&log_alpha, a, c, g, al, cl, gl, 
                                            s, p, w);
     A_i[i] = log_alpha;
  }
  fclose (qfile);
  qfile = NULL;
  if (want_debug_bf == 1) {
     psi_print_double_array (A_i, sample_size);
  }
  log_e2 = psi_sum_exp (A_i, sample_size);
  XFREE (A_i); 
  assert (isfinite(log_e2));
  r = psi_err_finite (log_e2);
  if (r != EXIT_SUCCESS) {
     psi_fatal ("log_e2 is not finite: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  *log_val = log_e2 - log_sample_size;
  assert (isfinite(*log_val));
  r = gmel_err_same (program_name, 
                     __FILE__, __LINE__, __PRETTY_FUNCTION__,
                     *log_val, log_e2);

  return r;
}

int 
gmel_bf_posterior_m4_block3_alpha (double *log_alpha, 
                                   double ak, double ck, double gk, 
                                   double a, double c, double g,
                                   double s, double p, int w)
{
  int r = EXIT_SUCCESS;
  double prior_top = 0;
  double prior_bot = 0;
  double jumping_top = 0;
  double jumping_bot = 0;
  double likelihood = 0;

  // Changed by Gibbs Sampler
  parameter top, bot;
  psi_st_parameter(&top, s, p, a, c, g);
  psi_st_parameter(&bot, s, p, ak, ck, gk);

  if (top.s == bot.s && top.p == bot.p 
      && top.a == bot.a && top.c == bot.c && top.g == bot.g) {
     likelihood = 0;
  } else {
     r = gmel_bf_likelihood_gibbs(&likelihood, top, bot);
     assert (r == EXIT_SUCCESS);
     assert (isfinite(likelihood));
  }

  r = gmel_bf_prior_m4(&prior_top, s, p, a, c, g);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(prior_top));

  r = gmel_bf_prior_m4(&prior_bot, s, p, ak, ck, gk);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(prior_bot));

  r = gmel_bf_posterior_m4_block3_q (&jumping_top, a, c, g, 
                                     ak, ck, gk, s, p, w);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(jumping_top));

  r = gmel_bf_posterior_m4_block3_q (&jumping_top, ak, ck, gk, 
                                     a, c, g, s, p, w);
  assert (r == EXIT_SUCCESS);
  assert (isfinite(jumping_bot));

  *log_alpha = likelihood
               + prior_top - prior_bot
               + jumping_top - jumping_bot; 
  assert (isfinite(*log_alpha));
  if (isfinite(*log_alpha) == 0) {
     psi_fatal ("log_alpah is not finite: %s %d %s",
                __FILE__, __LINE__, __PRETTY_FUNCTION__);
  }
  if (*log_alpha > 0) {
     *log_alpha = 0; 
  }

  return r;
}

int 
gmel_bf_posterior_m4_block3_q (double *q, 
                               double ak, double ck, double gk, 
                               double a, double c, double g,
                               double s, double p, int w)
{
  int r = EXIT_SUCCESS;
  r = gmel_bf_posterior_m1_block1_q(q, ak, ck, gk, a, c, g, w);
  assert (isfinite(*q));
  psi_err_finite (*q);

  s = p = 0.0; /* just for removing warning */
  return r;
}

/* */


int 
execute_bayesfactor ()
{
  FILE *bfp;
  double bf = 0;
  int r = EXIT_SUCCESS;
  parameter theta_star_model4;
  parameter theta_star_model1;
  parameter theta_dummy;
  int which_stage = 0;
  int remained_sample_size = 0;
  int s = access_mcmc_sample_size ();
  int f = access_mcmc_sample_freq ();
  int b = access_mcmc_sample_burn ();

  /* run mcmc for two models and save them into a file */
  if (want_load_param == 1) 
    {
/* here */
      r = load_estimation (&theta_star_model4, &theta_star_model1);
    } 
  else 
    {
     /* check the acceptance */
     write_mcmc_info ();
     /* initialize_mcmc_acceptance (); BUGGY */
     read_mcmc_acceptance (pname);

      if (want_load_state == 1) 
        {
          /* figure out in which state you have to start */
          which_stage = psi_io_check_bf_state (pname, sample_size);

          if (want_verbose == 1)
            fprintf (stderr, "BF: state %d\n", which_stage);

          /* remained_sample_size = psi_bf_stage_remained_sample_size (pname, sample_size); */
          if (which_stage > 0) 
            theta_star_model4 = psi_bf_stage_theta_star (pname, sample_size);
          /* theta_dummy = psi_bf_stage_theta_last (pname, sample_size); --> update */
        } 
      else 
        {
          remove (pname);
          if (want_save_alldraws == 1) 
            remove (aname);
          remove (qname);
        }

     switch (which_stage) { 
        case 0: 
          /* last point of burn is dealt with in the following estimation function */
          r = execute_estimation (0, 0, 0, 0, &theta_star_model4);
          if (r != EXIT_SUCCESS)
            return r;
          r = psi_io_mean_theta (pname, &theta_star_model4, sample_size, 0);
          if (r != EXIT_SUCCESS)
            return r;
          psi_nsv_report (ofile, pname, sample_size);
          write_acceptance_ratio ();
        case 1:
          r = execute_estimation (1, 0, theta_star_model4.s, 0, &theta_dummy);
          if (r != EXIT_SUCCESS)
            return r;
        case 2:
          r = execute_estimation (1, 1, theta_star_model4.s, theta_star_model4.p, 
                                  &theta_dummy);
          if (r != EXIT_SUCCESS)
            return r;
        case 3:
          r = execute_estimation (1, 1, 0, 0, &theta_star_model1);
          r = psi_io_mean_theta (pname, &theta_star_model1, sample_size, 3);
          if (r != EXIT_SUCCESS)
            return r;
        case 4:
          break;
      }
    }
  if (want_verbose == 1) {
     fprintf (stderr, "BF: after estimation\n");
  }

  /* estimate bayes factor */
  bfp = fopen (bname, "a");  
  if (bfp == NULL) {
     fprintf (stderr, "%s:%d:%s; %s of %s\n",
              __FILE__, __LINE__, __PRETTY_FUNCTION__,
              "could not open file", bname);
     psi_fatal ("could not open the file: bfp");
  }
  fprintf (bfp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
                theta_star_model4.a,
                theta_star_model4.c,
                theta_star_model4.g,
                theta_star_model4.t,
                theta_star_model4.s,
                theta_star_model4.p,
                theta_star_model1.a,
                theta_star_model1.c,
                theta_star_model1.g,
                theta_star_model1.t);
  r = gmel_dbf (&bf, theta_star_model4, theta_star_model1, bfp);
  if (r != EXIT_SUCCESS) {
     fprintf (stderr, "%s:%d:%s; %s\n",
              __FILE__, __LINE__, __PRETTY_FUNCTION__,
              "error in gmel_dbf");
     psi_fatal ("error in gmel_dbf");
  }
  fprintf (bfp, "%lf\n", bf);
  fclose (bfp);
  bfp = NULL;

  return r;
}

static int
load_estimation (parameter *theta_star_model4,
                 parameter *theta_star_model1)
{
  int i;
  FILE *bf;
  int r;
  int gen;
  double a, c, g, t, s, p;
  bf = fopen(pname, "r");
  if (bf == NULL) {
     fprintf (stderr, "%s:%d:%s; %s %d\n",
              __FILE__, __LINE__, __PRETTY_FUNCTION__, strerror (errno), errno);
     psi_fatal ("could not open the file: load_estimation");
  }
  theta_star_model4->a = 0; 
  theta_star_model4->c = 0; 
  theta_star_model4->g = 0; 
  theta_star_model4->t = 0; 
  theta_star_model4->s = 0; 
  theta_star_model4->p = 0; 
  theta_star_model1->a = 0; 
  theta_star_model1->c = 0; 
  theta_star_model1->g = 0; 
  theta_star_model1->t = 0; 

  pass_lines (1, bf);
  for (i = 0; i < sample_size; i++) {
     r = read_param_line (bf, &gen, &a, &c, &g, &t, &s, &p);
     if (r != EXIT_SUCCESS) {
        fprintf (stderr, "%s:%d:%s; %s\n",
                 __FILE__, __LINE__, __PRETTY_FUNCTION__,
                 "the sample of model 4 could not be read in");
        psi_fatal ("the sample of model 4 could not be read in: load_estimation");
     }
     if (want_debug_bf == 1) {
        fprintf (ofile, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", r, a, c, g, t, s, p);
     }
     theta_star_model4->a += a; 
     theta_star_model4->c += c; 
     theta_star_model4->g += g; 
     theta_star_model4->t += t; 
     theta_star_model4->s += s; 
     theta_star_model4->p += p; 
  }
  pass_lines (1, bf);
  pass_lines ((sample_size + 1) * 2 + 1, bf);
  for (i = 0; i < sample_size; i++) {
     r = read_param_line (bf, &gen, &a, &c, &g, &t, &s, &p);
     if (r != EXIT_SUCCESS) {
        fprintf (stderr, "%s:%d:%s; %s\n",
                 __FILE__, __LINE__, __PRETTY_FUNCTION__,
                 "the sample of model 1 could not be read in");
        psi_fatal ("the sample of model 1 could not be read in");
     }
     if (want_debug_bf == 1) {
        fprintf (stderr, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", r, a, c, g, t, s, p);
     }
     theta_star_model1->a += a; 
     theta_star_model1->c += c; 
     theta_star_model1->g += g; 
     theta_star_model1->t += t; 
  }
 
  theta_star_model4->a /= sample_size; 
  theta_star_model4->c /= sample_size; 
  theta_star_model4->g /= sample_size; 
  theta_star_model4->t = 1 - theta_star_model4->a - theta_star_model4->c - theta_star_model4->g;
  theta_star_model4->s /= sample_size; 
  theta_star_model4->p /= sample_size; 
  theta_star_model1->a /= sample_size; 
  theta_star_model1->c /= sample_size; 
  theta_star_model1->g /= sample_size; 
  theta_star_model1->t = 1 - theta_star_model1->a - theta_star_model1->c - theta_star_model1->g;
  theta_star_model1->s = 0; 
  theta_star_model1->p = 0; 
  fclose (bf);
  bf = NULL;

  return EXIT_SUCCESS;
}

