/* gibbs.c -- Gibbs sampler
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
#include "grid.h"
#include "gslwrap.h"
#include "gibbs.h"

#define PSI_NO_CHECK 0
static int psi_no_check = PSI_NO_CHECK; /* in config */
static int gibbs_size = 0;
static int gibbs_freq;
static int gibbs_burn;
static int want_debug_est = 0;
static int want_debug_gibbs = 0;
static int gmel_gibbs_usage = 0;
static int gmel_gibbs_total = 0;
static int gmel_gibbs_instance = 0;
static double *terms_sum = NULL;
static int *sampled_dna = NULL;
static int *sampled_pro = NULL;
extern int Nuc2AATable[4][4][4];
extern int Nuc2AATableStart[4][4][4];

int 
access_gibbs_size ()
{
  if (gibbs_size == 0)
    psi_fatal ("Gibbs sample size is zero. You have to determine its size before calling this function");
  return gibbs_size;
}

int
setup_gibbs (int s, int f, int b)
{
  assert (s > 0);
  assert (f > 0);
  assert (b > 0); /* let's avoid zero-burn-in */
  gibbs_size = s;
  gibbs_freq = f;
  gibbs_burn = b;
  terms_sum = XMALLOC (double, s);
  sampled_dna = XMALLOC (int, len_dna_jones_measure ());
  sampled_pro = XMALLOC (int, len_pro_jones_measure ());
  return EXIT_SUCCESS;
}

int
unsetup_gibbs ()
{
  XFREE (terms_sum);
  XFREE (sampled_dna);
  XFREE (sampled_pro);
  return EXIT_SUCCESS;
}

int
gmel_bf_likelihood_gibbs_sampler (GibbsPart *p_sampled_seqs, 
                                  parameter theta_star)
{
  /* NOTE: We only send in the root node sequence */
  int i, j, k;
  double E_solv, E_pair;
  double log_theta_star[PSI_NUM_SAMPLE];

  /* Argument Assignment */
  int *protein = NULL;
  int *dna = dna_jones_measure ();
  int len_dna = len_dna_jones_measure ();
  int len_pro = len_dna / 3;

  protein = XMALLOC (int, len_pro);
  assert (len_dna % 3 == 0);
  
  GibbsPart *gibbsSeq = p_sampled_seqs;

  /* START OF WORK */
  /* Copying DNA and protein sequence for sequence sampling */
  /* psi_rng_dna_seq (sampled_dna, len_dna); */
  memcpy (sampled_dna, dna, sizeof (int) * len_dna);
  dna2protein (sampled_dna, sampled_pro, len_pro);

  log_theta_star[PSI_SAMPLE_A] = log (theta_star.a);
  log_theta_star[PSI_SAMPLE_C] = log (theta_star.c);
  log_theta_star[PSI_SAMPLE_G] = log (theta_star.g);
  log_theta_star[PSI_SAMPLE_T] = log (theta_star.t);
  log_theta_star[PSI_SAMPLE_S] = theta_star.s;
  log_theta_star[PSI_SAMPLE_P] = theta_star.p;

  char *str_dna = NULL;
  char *str_pro = NULL;
  if (want_debug_gibbs == 1)
    {
      str_dna = XMALLOC (char, len_dna + 1); 
      str_pro = XMALLOC (char, len_pro + 1);
    }

  int site;
  for (i = 0; i < gibbs_burn; i++)
    {
      for (j = 0; j < len_dna; j++)
        {
          site = choose_a_site (len_dna);
          /* we change the DNA and protein sequence */
          choose_a_nuc_stationary_fast (PSI_ENERGY_JONES, 
                                        sampled_dna, 
                                        sampled_pro, 
                                        site, 
                                        log_theta_star);
        }
    }

  for (i = 0; i < gibbs_size; i++)
    {
      for (j = 0; j < gibbs_freq; j++)
        {
          for (k = 0; k < len_dna; k++)
            {
              site = choose_a_site (len_dna);
              /* we change the DNA and protein sequence */
              choose_a_nuc_stationary_fast (PSI_ENERGY_JONES, 
                                            sampled_dna, 
                                            sampled_pro, 
                                            site, 
                                            log_theta_star);
            }
        }
      /* store the sampled sequence */
      score_drevol (PSI_ENERGY_JONES, sampled_pro, &E_solv, &E_pair);

      /* CHECK for debug: match the translation and protein sequence */
      dna2protein (sampled_dna, protein, len_pro);
      assert (psi_seq_diff_seqs (sampled_pro, protein, len_pro) == 0);

      gibbsSeq->S = E_solv;
      gibbsSeq->P = E_pair;
      gibbsSeq->A = num_nucleotide (sampled_dna, len_dna, PSI_DNA_A);
      gibbsSeq->C = num_nucleotide (sampled_dna, len_dna, PSI_DNA_C);
      gibbsSeq->G = num_nucleotide (sampled_dna, len_dna, PSI_DNA_G);
      gibbsSeq++;
    }

  XFREE (protein);

  if (want_debug_gibbs == 1)
    {
      XFREE (str_dna);
      XFREE (str_pro);
    }

  return EXIT_SUCCESS;
}

int
gmel_bf_likelihood_sumup_gibbs (double *sum, double *max_energy, 
                                parameter theta, parameter theta_star)
{
  int r = EXIT_SUCCESS;
  int i;
  GibbsPart *sampled_seqs = NULL;
  double diffS, diffP, diffA, diffC, diffG, diffT;
  int len_dna = len_dna_jones_measure ();

  diffS = theta.s - theta_star.s;
  diffP = theta.p - theta_star.p;

  if (psi_no_check == 0)
    {
      psi_err_base (theta.a, theta.c, theta.g, theta.t);
      psi_err_base (theta_star.a, theta_star.c, theta_star.g, theta_star.t);
    }

  diffA = log (theta.a) - log (theta_star.a);
  diffC = log (theta.c) - log (theta_star.c);
  diffG = log (theta.g) - log (theta_star.g);
  diffT = log (theta.t) - log (theta_star.t);
 
  if (psi_no_check == 0)
    {
      if (!isfinite(diffS) || !isfinite(diffP)
          || !isfinite(diffA) || !isfinite(diffC) 
          || !isfinite(diffG) || !isfinite(diffT) ) 
        {
          if (want_debug_gibbs == 1) 
            {
/*
        debug_write_parameter (theta);
        debug_write_parameter (theta_star);
*/
            }
          psi_fatal ("diff is not finite");
        } 
    } 

  /* theta_star is a vector of parameter values of a given gridpoint,
     or the gridpoint's parameter values are the same as those of parameter
     itself
   */
  gmel_grid_get_info (&sampled_seqs, theta_star.gp); 
  gmel_gibbs_total++; 
  if (sampled_seqs == NULL) 
    {
      gmel_gibbs_instance++;
      sampled_seqs = XMALLOC (GibbsPart,  gibbs_size);
      gmel_bf_likelihood_gibbs_sampler (sampled_seqs, theta_star);
      r = gmel_grid_put_info (sampled_seqs, theta_star.gp);
      assert (r == EXIT_SUCCESS);
    } 
  else 
    {
      gmel_gibbs_usage++;
    }

  for (i = 0; i < gibbs_size; i++) 
    {
      terms_sum[i]= sampled_seqs[i].S * diffS          /* S <- -2*E_s(h) */
                    + sampled_seqs[i].P * diffP        /* P <- -2*E_p(h) */
                    + sampled_seqs[i].A * diffA        /* A <- Number of A */
                    + sampled_seqs[i].C * diffC        /* C <- Number of C */
                    + sampled_seqs[i].G * diffG        /* G <- Number of G */
                    + (len_dna - (sampled_seqs[i].A + sampled_seqs[i].C + sampled_seqs[i].G)) * diffT;
      if (isfinite(terms_sum[i]) == 0) 
        {
          fprintf (stderr,  
                   "%s:%s:%d:%s\nV: %lf\nS: %lf * %lf, %lf %lf\nP: %lf * %lf, %lf %lf\nA: %d * %lf, %lf %lf\nC: %d * %lf, %lf %lf\nG: %d * %lf, %lf %lf\nT: %d * %lf, %lf %lf\n",
                   program_name,
                   __FILE__, __LINE__, __PRETTY_FUNCTION__, 
                   terms_sum[i], 
                   sampled_seqs[i].S, diffS, theta.s, theta_star.s,
                   sampled_seqs[i].P, diffP, theta.p, theta_star.p,
                   sampled_seqs[i].A, diffA, theta.a, theta_star.a,
                   sampled_seqs[i].C, diffC, theta.c, theta_star.c,
                   sampled_seqs[i].G, diffG, theta.g, theta_star.g,
                   len_dna - (sampled_seqs[i].A + sampled_seqs[i].C + sampled_seqs[i].G), diffT, theta.t, theta_star.t);
          psi_fatal ("gibbs sum is not finite");
        }
   
      if (want_debug_est == 1) 
        {
          fprintf (stderr, "%4d: %+.3lf - %+.3lf : %+.3lf * %+.3lf = %+.3lf, %+.3lf * %+.3lf = %+.3lf, %3d * %+.3lf = %+.3lf, %3d * %+.3lf = %+.3lf, %3d * %+.3lf = %+.3lf, %3d * %+.3lf = %+.3lf\n", 
                   i, *max_energy, terms_sum[i],
                   sampled_seqs[i].S, diffS, sampled_seqs[i].S * diffS,
                   sampled_seqs[i].P, diffP, sampled_seqs[i].P * diffP,
                   sampled_seqs[i].A, diffA, sampled_seqs[i].A * diffA,
                   sampled_seqs[i].C, diffC, sampled_seqs[i].C * diffC,
                   sampled_seqs[i].G, diffG, sampled_seqs[i].G * diffG,
                   (len_dna - (sampled_seqs[i].A + sampled_seqs[i].C + sampled_seqs[i].G)), diffT,
                   (len_dna - (sampled_seqs[i].A + sampled_seqs[i].C + sampled_seqs[i].G)) * (diffT));
        }
    }

  *sum = psi_sum_exp (terms_sum, gibbs_size);

  return EXIT_SUCCESS;
}

int
gmel_gibbs_write_usage (FILE *fp)
{
  fprintf(fp, "Gibbs Usage:\t%.1lf%%\t%d\n",
          (double) gmel_gibbs_usage*100/gmel_gibbs_total, gmel_gibbs_instance);
  return 0;
}
