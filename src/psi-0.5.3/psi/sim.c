/* sim.c -- Generating simulated data with the parameters fixed
  memcpy (sampled_dna, t2NUC, sizeof(int)*len_dna);
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "common.h"
#include "error.h"
#include "seq.h"
#include "rng.h"
#include "energy.h"
#include "sim.h"

void 
psi_sim_gibbs (int iter)
{
  /* 0. start with a initial sequence
     1. choose a site `k' uniformly between 1 and N
     2. choose A or C or G or T at the site `k' 
        according to relative stationary probabilities 
   */
  int r;
  int i;
  int site;
  int nucleotide;
  int len_dna = len_dna_jones_measure ();
  int len_pro = len_pro_jones_measure ();
  int *dna = NULL;
  int *dna_prev = NULL;
  int *pro = NULL;
  int *pro_prev = NULL;
  char *dna_str = NULL;
  double p[6];
  p[PSI_SAMPLE_A] = 0.25;
  p[PSI_SAMPLE_C] = 0.25;
  p[PSI_SAMPLE_G] = 0.25;
  p[PSI_SAMPLE_T] = 0.25;
  p[PSI_SAMPLE_S] = 0.7;
  p[PSI_SAMPLE_P] = 0.005;
  
  int j = 0;
  int l = 0;
  int renewed = 0;
  int delta = 0;
  pro = XMALLOC (int, len_pro);
  pro_prev = XMALLOC (int, len_pro);
  dna = XMALLOC (int, len_dna);
  dna_prev = XMALLOC (int, len_dna);
  for (i = 0; i < len_dna; i++)
    {
      /* dna[i] = dat_dna[i]; */
      dna[i] = 0;
      dna_prev[i] = 3;
    }

  r = EXIT_FAILURE;
  while (r == EXIT_FAILURE)
    {
      psi_rng_dna_seq (dna, len_dna);
      r = dna2protein (dna, pro, len_pro);
    }
  
  double solv, pair;
  score_drevol (PSI_ENERGY_JONES, pro, &solv, &pair); 

  /* HEAD */
  fprintf (stderr, "gen\tEs\tEp\tA\tC\tG\tT\tchange\tdiffD\tdiffP\n");
          psi_energy_dna_report (PSI_ENERGY_JONES, dna, p);
          dna2protein (dna, pro, len_pro);
          dna2protein (dna_prev, pro_prev, len_pro);
          fprintf (stderr, "%.2lf\t%d\t%d\n", 
                   (double) 0,
                   psi_seq_diff_seqs (dna, dna_prev, len_dna),
                   psi_seq_diff_seqs (pro, pro_prev, len_pro));

  int burnin_sim = 10000;
  for (j = 0; j < burnin_sim; j++)
    {
      site = choose_a_site (len_dna);
      assert (site < len_dna);
      assert (site > -1); 
      nucleotide = choose_a_nuc_stationary (PSI_ENERGY_JONES, dna, site, p);
      dna[site] = nucleotide;
    }

  double **state_value = NULL;
  state_value = XMALLOC (double *, 6);
  for (i = 0; i < 6; i++)
    {
      state_value[i] = XMALLOC (double, iter);
    }

  delta = 2048;
  for (i = 0; i < iter; i++)
    {
      renewed = 0;
      for (j = 0; j < delta; j++)
        {
          site = choose_a_site (len_dna);
          assert (site < len_dna);
          assert (site > -1); 
          nucleotide = choose_a_nuc_stationary (PSI_ENERGY_JONES, dna, site, p);
          if (dna[site] != nucleotide)
            renewed++;
          dna[site] = nucleotide;
        }
      
      psi_energy_dna_report (PSI_ENERGY_JONES, dna, p);
      psi_energy_dna_state (PSI_ENERGY_JONES, dna, p, state_value, i);
      dna2protein (dna, pro, len_pro);
      dna2protein (dna_prev, pro_prev, len_pro);
      fprintf (stderr, "%.2lf\t%d\t%d\n", 
               (double) renewed / delta,
               psi_seq_diff_seqs (dna, dna_prev, len_dna),
               psi_seq_diff_seqs (pro, pro_prev, len_pro));

      for (l = 0; l < len_dna; l++)
        {
          dna_prev[l] = dna[l]; 
        }
    } 
  
     
  r = EXIT_FAILURE;
  while (r == EXIT_FAILURE)
    {
      site = choose_a_site (len_dna);
      assert (site < len_dna);
      assert (site > -1); 
      nucleotide = choose_a_nuc_stationary (PSI_ENERGY_JONES, dna, site, p);
      dna[site] = nucleotide;

      r = psi_energy_valid_state (PSI_ENERGY_JONES, dna, p, state_value, iter);
    }

  dna_str = XMALLOC (char, len_dna + 1);
  dna_int2char (dna, dna_str, len_dna);
  printf ("%s\n", dna_str);

  XFREE (pro);
  XFREE (pro_prev);
  XFREE (dna);
  XFREE (dna_prev);
  XFREE (dna_str);
}

void
psi_sim_rate (int iter)
{
  /* 0. start with a initial sequence
     1. sequence `i' is the current sequence
     2. calculate Rij for all neighbors j (3N of these neighbors)
     3. choose one of 3N neighbor sequences with a probability
        Rik/\sum_j Rij 
     4. go back to 1
   */
  int i;
  int site;
  int nucleotide;
  int len_dna = len_dna_jones_measure ();
  int *dna = NULL;
  char *dna_str = NULL;
  double p[8];
  p[PSI_SAMPLE_A] = 0.25;
  p[PSI_SAMPLE_C] = 0.25;
  p[PSI_SAMPLE_G] = 0.25;
  p[PSI_SAMPLE_T] = 0.25;
  p[PSI_SAMPLE_S] = 0.50;
  p[PSI_SAMPLE_P] = 0.05;
  p[PSI_SAMPLE_K] = 1.0;
  p[PSI_SAMPLE_W] = 1.0;

  dna = XMALLOC (int, len_dna);
  for (i = 0; i < len_dna; i++)
    dna[i] = 0;

  for (i = 0; i < iter; i++)
    {
      choose_a_nuc_rate (PSI_ENERGY_JONES, &site, &nucleotide, dna, p);
      assert (site < len_dna);
      assert (site > -1); 
      dna[site] = nucleotide;
    } 
  
  dna_str = XMALLOC (char, len_dna + 1);
  dna_int2char (dna, dna_str, len_dna);
  printf ("%s\n", dna_str);

  XFREE (dna);
  XFREE (dna_str);
}


