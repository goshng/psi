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
  double solv, pair;
/* for pdb1hqb.A.prd
  int len_protein = 80;
  int protein[80] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
*/

  /* for pdb1r48.A.dat
     DNIEQKIDDIDHEIADLQAKRTRLVQQHPR 
   */

  /* setup something like filenames before init*/
  setup_energy_jones ("../data/pdb1r48.A.dat", "../data/energy.int");
  int protein[30] = { 3,  2,  9,  6,  5, 11,  9,  3,  3,  9,  
                      3,  8,  6,  9,  0,  3, 10,  5,  0, 11,  
                      1, 16,  1, 10, 19,  5,  5,  8, 14,  1 };

/* 
  setup_energy_jones ("../data/pdbtest.2.dat", "../data/energy.int"); 
  int len_protein = 2;
  int protein[2] = { 3,  2 };
*/

  /* initialize the energy */
  initialize_drevol (PSI_ENERGY_JONES);

/*
  write_dat (stderr);
  write_energy (stderr);
*/

  /* call the energy module */
  score_drevol (PSI_ENERGY_JONES, protein, &solv, &pair);
  solv *= 1000;
  pair *= 1000;
  if (trunc (solv) != 9896 || trunc (pair) != 31089)
    {
      return EXIT_FAILURE;
    }


/*
  write_dat (stderr);
  write_energy (stderr);
*/

  /* finalize the energy */
  finalize_drevol (PSI_ENERGY_JONES);

  exit (result);
}
/** @end 1 */

/*
  solv *= 1000;
  pair *= 1000;
  if (trunc (solv) != 9896 || trunc (pair) != 31089)
    {
      return EXIT_FAILURE;
    }

  energy_is (protein, len_protein, &solv, &pair);
  solv *= 1000;
  pair *= 1000;
  if (trunc (solv) != 9896 || trunc (pair) != 31089)
    {
      return EXIT_FAILURE;
    }
*/


