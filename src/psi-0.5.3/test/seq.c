/* seq.c -- test program of the seq module
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

#include "psi.h"

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
  int i, j, k;
  int r;
  int gc_prt[4] = { 3, 1, 0, 2 };

  /* read_translaton_table ("gc.prt"); only for the start */ 

  int codon[3];
  int dna[LEN_DNA];
  int protein[LEN_PROTEIN];
  char dna_str[LEN_DNA + 1];
  char protein_str[LEN_PROTEIN + 1];
  for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 4; j++)
        {
          for (k = 0; k < 4; k++)
            {
              dna[4*4*3*i + 4*3*j + 3*k + 0] = gc_prt[i];
              dna[4*4*3*i + 4*3*j + 3*k + 1] = gc_prt[j];
              dna[4*4*3*i + 4*3*j + 3*k + 2] = gc_prt[k];
            }
        }
    }
  dna_int2char (dna, dna_str, LEN_DNA);
  i = num_nucleotide (dna, LEN_DNA, PSI_DNA_A); 
  if (i != LEN_DNA/4)
    {
      return EXIT_FAILURE;
    }
  i = num_nucleotide (dna, LEN_DNA, PSI_DNA_C); 
  if (i != LEN_DNA/4)
    {
      return EXIT_FAILURE;
    }
  i = num_nucleotide (dna, LEN_DNA, PSI_DNA_G); 
  if (i != LEN_DNA/4)
    {
      return EXIT_FAILURE;
    }
  i = num_nucleotide (dna, LEN_DNA, PSI_DNA_T); 
  if (i != LEN_DNA/4)
    {
      return EXIT_FAILURE;
    }

  find_codon_in_dna (codon, dna, 0);
  if (codon[0] != PSI_DNA_T || codon[1] != PSI_DNA_T || codon[2] != PSI_DNA_T)
    {
      return EXIT_FAILURE;
    }
  find_codon_in_dna (codon, dna, 19);
  if (codon[0] != PSI_DNA_T || codon[1] != PSI_DNA_C || codon[2] != PSI_DNA_A)
    {
      return EXIT_FAILURE;
    }
  find_codon_in_dna (codon, dna, LEN_DNA - 1);
  if (codon[0] != PSI_DNA_G || codon[1] != PSI_DNA_G || codon[2] != PSI_DNA_G)
    {
      return EXIT_FAILURE;
    }

  choose_transl_table (1, 0);
  r = dna2protein (dna, protein, LEN_PROTEIN);
  if (r == EXIT_SUCCESS)
    exit (EXIT_FAILURE);
  pro_int2char (protein, protein_str, LEN_PROTEIN);
  if (strcmp (protein_str, 
              "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
    {
      printf ("Wrong translation for id 1!\n%s\n", protein_str);
      exit (EXIT_FAILURE);
    }

  choose_transl_table (4, 0);
  r = dna2protein (dna, protein, LEN_PROTEIN);
  if (r == EXIT_SUCCESS)
    exit (EXIT_FAILURE);
  pro_int2char (protein, protein_str, LEN_PROTEIN);
  if (strcmp (protein_str, 
              "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
    {
      printf ("Wrong translation for id 4!\n");
      exit (EXIT_FAILURE);
    }

  choose_transl_table (11, 0);
  r = dna2protein (dna, protein, LEN_PROTEIN);
  if (r == EXIT_SUCCESS)
    exit (EXIT_FAILURE);
  pro_int2char (protein, protein_str, LEN_PROTEIN);
  if (strcmp (protein_str, 
              "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
    {
      printf ("Wrong translation for id 11!\n");
      exit (EXIT_FAILURE);
    }

  exit (result);
}

/** @end 1 */
