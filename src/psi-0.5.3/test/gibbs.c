/* grid.c -- test program of the grid module
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

#ifndef LEN_DNA
#  define LEN_DNA 192
#endif

#ifndef LEN_PROTEIN
#  define LEN_PROTEIN 64
#endif

/* static int psi_init  (Sic *psi); */

/** @start 1 */
int
/* main (int argc, char * const argv[]) */
main (void)
{
  int result = EXIT_SUCCESS;

  /* setup something like filenames before init*/
  /* setup_energy_jones ("../data/pdb1r48.A.dat", "../data/energy.int"); */
  setup_energy_jones ("../data/pdbtest.2.dat", "../data/energy.int");

  /* initialize the energy */
  initialize_drevol (PSI_ENERGY_JONES);

/*
  write_dat (stderr);
  write_energy (stderr);
*/

/* call the energy module 
  score_drevol (PSI_ENERGY_JONES, protein, &solv, &pair);
  printf ("solv: %lf, pair: %lf\n", solv, pair);
  energy_is (protein, len_protein, &solv, &pair);
  printf ("solv: %lf, pair: %lf\n", solv, pair);
*/

  init_rng (0);

  int gibbs_size = 2;
  int gibbs_burn = 1;
  int gibbs_freq = 1;
  setup_gibbs (gibbs_size, gibbs_burn, gibbs_freq);
  parameter theta;
  theta.s = 0.5; theta.p = 0.03; 
  theta.a = 0.2; theta.c = 0.3; theta.g = 0.1; theta.t = 0.4;

  GibbsPart *sampled_seqs = XMALLOC (GibbsPart,  gibbs_size);
  gmel_bf_likelihood_gibbs_sampler (sampled_seqs, theta);

  if (sampled_seqs[0].A != 1 || sampled_seqs[0].C != 3
      || sampled_seqs[0].G != 0 
      || sampled_seqs[1].A != 2 || sampled_seqs[1].C != 2 
      || sampled_seqs[1].G != 1
      || trunc (sampled_seqs[0].S*1000) != -1241
      || trunc (sampled_seqs[0].P*1000) != -106
      || trunc (sampled_seqs[1].S*1000) != 649
      || trunc (sampled_seqs[1].P*1000) != 630)
    {
      return EXIT_FAILURE; 
    }
  XFREE (sampled_seqs);
  fin_rng ();

  /* finalize the energy */
  finalize_drevol (PSI_ENERGY_JONES);

  exit (result);
}
/** @end 1 */
