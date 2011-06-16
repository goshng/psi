/* rng.h -- random number generator
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

/*!\file rng.h
   \author Sang Chul Choi
   \brief A wrapper module of GSL random number generator module
 */

/** @start 1 **/
#ifndef PSI_RNG_H
#define PSI_RNG_H 1

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_math.h>
#include <psi/common.h>

BEGIN_C_DECLS

extern double rng_put_out (double d);
extern int init_rng (unsigned long int s);
extern int save_rng (const char *rngname);
extern int load_rng (const char *rngname);
extern void psi_rng_info (FILE *fp);
extern double psi_rng ();
extern int fin_rng ();
extern int choose_one_nucleotide (double* here); 
extern int choose_one_nucleotide_with_fixed (double *here, double v);
extern int psi_ran_pi (const double pi[], double r_pi[]);
extern double psi_ran_pi_lnpdf (const double alpha[], const double theta[]);
extern double psi_ran_flat (double a, double b);
extern double psi_sample_ran_flat (double v, double delta, 
                                   double min, double max);
extern double psi_ran_flat_pdf (double x, double mu, double delta, 
                                double min, double max);
extern int choose_a_site (int len_dna);
extern int choose_one_sequence_by_rate (double *here, int len_dna);
extern int choose_one_element_with_fixed (double *here, int n, double uniform_variable);
extern int psi_rng_dna_seq (int *dna, int n);

END_C_DECLS

#endif /* !PSI_RNG_H */
/** @end 1 **/
