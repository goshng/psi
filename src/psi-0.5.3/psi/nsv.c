/* nsv.c -- Delta Measure
   It calculates nonsynonymous rate due to structure, which will allow 
   the analysis of Nonsynonymous Rate Variation

   Copyright (C) 2005-2006 Sang Chul Choi
  
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/*
   Test file is the prd and p files of pdb1afp._
*/
#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "common.h"
#include "error.h"
#include "gslwrap.h"
#include "energy.h"

/* #include <termios.h> */
/* #include <gsl/gsl_fit.h> will be removed */
/* #include <getopt.h> */

/* #define N_SAMPLE 1000 */
extern int Nuc2AATable[4][4][4];
/* static int sample_size = N_SAMPLE; */
/* static double aa_codon_frequency[NUM_AMINOACID]; */

static double fitness_s = 0.0;
static double fitness_p = 0.0;

static double** psi_mem_double_2 (int n_1, int n_2);
static void psi_mem_double_free_2 (double ***m, int n_2);
static void psi_nsv_statistics (const double data[], size_t n, 
                                double *m, double *sd, double *median, 
                                double *q05, double *q25, 
                                double *q75, double *q95);
static double** psi_mem_double_2 (int n_1, int n_2)
{
  int i;
  double **m = XMALLOC (double *, n_1);
  for (i = 0; i < n_1; i++)
    {
      m[i] = XMALLOC (double, n_2);
    }
  return m; 
}
static void psi_mem_double_free_2 (double ***m, int n_2)
{
  int i;
  for (i = 0; i < n_2; i++)
    {
      XFREE ((*m)[i]);
    }
  XFREE (*m);
}

static void psi_nsv_statistics (const double data[], size_t n, 
                                double *m, double *sd, double *median, 
                                double *q05, double *q25, 
                                double *q75, double *q95)
{
  *m = psi_stats_mean (data, n);
  *sd = psi_stats_sd (data, n);
  *median = psi_stats_median (data, n);
  *q05 = psi_stats_quantile (data, n, 0.05);
  *q25 = psi_stats_quantile (data, n, 0.25);
  *q75 = psi_stats_quantile (data, n, 0.75);
  *q95 = psi_stats_quantile (data, n, 0.95);
}

static int
psi_nsv_delta (FILE *fp);

int
psi_nsv_report (FILE *fp, const char *pn, int sample_size)
{
  int i;
  int r;
  int gen;
  char c;
  double **post_sample = psi_mem_double_2 (PSI_NUM_SAMPLE, sample_size);

  FILE *pf = fopen (pn, "r");
  /* read one line */
  while ((c = fgetc (pf)) != '\n') {}
  for (i = 0; i < sample_size; i++)
    {
      r = fscanf (pf, "%d%lf%lf%lf%lf%lf%lf\n", 
                  &gen, 
                  &post_sample[0][i], &post_sample[1][i], &post_sample[2][i],
                  &post_sample[3][i], &post_sample[4][i], &post_sample[5][i]);
      assert (r == 7);
    }

/*
  for (i = 0; i < sample_size; i++)
    {
      for (j = 0; j < PSI_NUM_SAMPLE; j++)
        {
          fprintf (stderr, "%lf ", post_sample[j][i]);
        }
      fprintf (stderr, "\n");
    }
*/
 
  double m, sd, median, q05, q25, q75, q95;
  psi_nsv_statistics (post_sample[PSI_SAMPLE_A], sample_size,
                      &m, &sd, &median, &q05, &q25, &q75, &q95);
  fprintf (fp, "REPORT A:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
           m, sd, median, q05, q25, q75, q95);
  psi_nsv_statistics (post_sample[PSI_SAMPLE_C], sample_size,
                      &m, &sd, &median, &q05, &q25, &q75, &q95);
  fprintf (fp, "REPORT C:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
           m, sd, median, q05, q25, q75, q95);
  psi_nsv_statistics (post_sample[PSI_SAMPLE_G], sample_size,
                      &m, &sd, &median, &q05, &q25, &q75, &q95);
  fprintf (fp, "REPORT G:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
           m, sd, median, q05, q25, q75, q95);
  psi_nsv_statistics (post_sample[PSI_SAMPLE_T], sample_size,
                      &m, &sd, &median, &q05, &q25, &q75, &q95);
  fprintf (fp, "REPORT T:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
           m, sd, median, q05, q25, q75, q95);
  psi_nsv_statistics (post_sample[PSI_SAMPLE_S], sample_size,
                      &m, &sd, &median, &q05, &q25, &q75, &q95);
  fprintf (fp, "REPORT S:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
           m, sd, median, q05, q25, q75, q95);
  fitness_s = m;
  psi_nsv_statistics (post_sample[PSI_SAMPLE_P], sample_size,
                      &m, &sd, &median, &q05, &q25, &q75, &q95);
  fprintf (fp, "REPORT P:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
           m, sd, median, q05, q25, q75, q95);
  fitness_p = m;
  
  psi_mem_double_free_2 (&post_sample, PSI_NUM_SAMPLE);
  
  psi_nsv_delta (fp);

  return EXIT_SUCCESS;
}

static int
psi_nsv_delta (FILE *fp)
{
  int i, j;
  double solv, pair;
  double solv_i, pair_i;
  int codon[3];
  int nuc;
  int aa;
  int aa_replaced;
  int len_pro, len_dna;
  int *dna = NULL;
  int *protein = NULL;
  int *dat_dna = NULL;
  int *dat_pro = NULL;
  double *s_delta = NULL;
  double *p_delta = NULL;
  double *s_p_delta = NULL;
  int n = 0;
  double s_delta_mean, s_delta_var, p_delta_mean, p_delta_var;
  double s_delta_plus_p_delta_mean, s_delta_plus_p_delta_var;
  double s_p_covariance, s_p_correlation;

  len_dna = len_dna_jones_measure ();
  len_pro = len_pro_jones_measure ();
  dat_dna = dna_jones_measure ();
  dat_pro = pro_jones_measure ();
 
  dna = XMALLOC (int, len_dna);
  protein = XMALLOC (int, len_pro);
  s_delta = XMALLOC (double, len_dna * 3);
  p_delta = XMALLOC (double, len_dna * 3);
  s_p_delta = XMALLOC (double, len_dna * 3);

  for (i = 0; i < len_dna; i++) {
    dna[i] = dat_dna[i];
  }
  for (i = 0; i < len_pro; i++) {
    protein[i] = dat_pro[i];
  }

  score_drevol (PSI_ENERGY_JONES, protein, &solv_i, &pair_i);
  solv_i /= -2.0; /* This is dirty but we always did that */
  pair_i /= -2.0; /* This is dirty but we always did that */
 
/*
  if (dosnp == 1) {
    calculate_snp ();
  } 
  fprintf (ofile, "#%s\n", datname);
  fprintf (ofile, "splusp\t\ts\t\tp\tn\n");
*/
  
  /* fprintf(ofile, "SOLV: %lf, PAIR: %lf\n", solv, pair); */

  for (i = 0; i < len_dna; i++) 
    {
      nuc = dna[i]; 
      codon[0] = dna[i - i%3];
      codon[1] = dna[i + 1 - i%3];
      codon[2] = dna[i + 2 - i%3];
      for (j = 1; j < 4; j++) 
        {
          nuc++;
          nuc %= 4;
          codon[i%3] = nuc;
          aa = Nuc2AATable[codon[0]][codon[1]][codon[2]];  
          if (aa != protein[i/3] && aa != 20) 
            {
              aa_replaced = protein[i/3];
              protein[i/3] = aa;

/*
        fprintf (ofile, "%3d - %3d - %d : 01234567890123456789012345678901234567890123456789\n", i, i/3, j);
        fprintf (ofile, "%3d - %3d - %d : ", i, i/3, j);
        for (k = 0; k < Dat->len_pro; k++) {
          fprintf (ofile, "%c", AMINOACID[protein[k]]);
        }
        fprintf (ofile, "\n");
*/

              score_drevol (PSI_ENERGY_JONES, protein, &solv, &pair);
              solv /= -2.0;
              pair /= -2.0;
/*
        fprintf(ofile, "%lf\t%lf\t%lf\n", 
                fitness_s*(solv_i - solv) + fitness_p*(pair_i - pair),
                fitness_s*(solv_i - solv), 
                fitness_p*(pair_i - pair));
*/
              s_delta[n] = fitness_s*(solv_i - solv);
              p_delta[n] = fitness_p*(pair_i - pair);
              s_p_delta[n] = s_delta[n] + p_delta[n];
              n++;
              protein[i/3] = aa_replaced;
            }
        }
    }
/*
  fprintf (ofile, "\n");
*/

  XFREE (dna);
  XFREE (protein);

  /* find the best fitting line */
  /* a b c d r1 r2 a b c d r1 r2  c0 c1 cov00 cov01 cov11 chisq c0 c1 cov00 cov01 cov11 chisq c1 cov11 chisq c1 cov11 chisq s_delta_mean s_delta_sd p_delta_mean p_delta_sd s_delta_plus_p_delta_mean s_delta_plus_p_delta_sd s_p_correlation s_p_covariance s p*/
  /* s_delta_mean s_delta_sd p_delta_mean p_delta_sd s_delta_plus_p_delta_mean s_delta_plus_p_delta_sd s_p_correlation s_p_covariance s p*/

  s_delta_mean = psi_stats_mean (s_delta, n);
  s_delta_var = psi_stats_var (s_delta, n);
  p_delta_mean = psi_stats_mean (p_delta, n);
  p_delta_var = psi_stats_var (p_delta, n);
  s_p_covariance = psi_stats_covariance (s_delta, p_delta, n);
  s_p_correlation = s_p_covariance / (sqrt(s_delta_var) * sqrt(p_delta_var));
  s_delta_plus_p_delta_mean = s_delta_mean + p_delta_mean;
  s_delta_plus_p_delta_var = psi_stats_var (s_p_delta, n);

  fprintf (fp, "DELTA:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
           s_delta_mean,
           s_delta_var,
           p_delta_mean,
           p_delta_var,
           s_p_covariance,
           s_p_correlation,
           s_delta_plus_p_delta_mean,
           s_delta_plus_p_delta_var);
 
  XFREE (s_delta);
  XFREE (p_delta);
  XFREE (s_p_delta);

  return EXIT_SUCCESS;
}



