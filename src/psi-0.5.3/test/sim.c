/* sim.c -- Simulation Test module
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

/* static const int total_iter = 1000000; */
/* static const int total_iter = 200000000; */
static const int total_iter = 100; 
static char *oname = NULL;
static FILE *ofile;

static int unsetup_psi ();
static int setup_psi ();

/** @start 1 */
int
/* main (int argc, char * const argv[]) */
main (void)
{
  int result = EXIT_SUCCESS;
 
  if (oname == NULL) 
    ofile = stdout;
  else 
    ofile = fopen (oname, "w");

  setup_psi ();

  psi_sim_gibbs (total_iter);
/*
  psi_sim_rate (total_iter);
*/

  unsetup_psi ();

  if (ofile != stdout)
    fclose (ofile);

  exit (result);
}
/** @end 1 */

int 
unsetup_psi ()
{
  /* delete_gridpoints (); 
     unsetup_mcmc ();
     we don't need grid or mcmc for Nonsynonymous Structure Variation module */

  finalize_drevol (PSI_ENERGY_JONES);
  fin_rng();
  return EXIT_SUCCESS;
}

int 
setup_psi ()
{
  /* We don't need a random number generator either */
  init_rng (0);

  /* We only need the energy module */
  setup_energy_jones ("../data/pdb1afp._.prd", "../data/energy.int"); 
  /* setup_energy_jones ("../data/pdb1afp.T.prd", "../data/energy.int"); */
  /* setup_energy_jones ("../data/pdb1afp.U.prd", "../data/energy.int"); */
  initialize_drevol (PSI_ENERGY_JONES);

  return EXIT_SUCCESS;
}


