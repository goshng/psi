/* gslwrap.h -- GNU Scientific Library wrapper functions
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

/*!File gslwrap.h
   Author Sang Chul Choi
   Brief Energy Calculation Module

 */

/** @start 1 **/
#ifndef PSI_GSLWRAP_H
#define PSI_GSLWRAP_H 1

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_math.h>
#include <psi/common.h>

BEGIN_C_DECLS

extern int psi_int_function (int a);
extern double psi_function (double a);

extern double psi_gaussian_lnpdf (double mu, double sigma, double x);
extern double psi_sum_exp (double *log_i, int n);
extern void psi_print_double_array (double *A_i, int n);
extern double logsum_nocheck (double x, double y);
extern double logsum (double x, double y);
extern void loginc (double *x, double y);
extern double logdiff (double x, double y);
extern double psi_pi_lnpdf (const double alpha[], double theta[]);
extern int gmel_conv_psr (double *psr, double **draws, int m, int n);
extern int get_sample_variance (int n, double *v, double *x, double *y);
extern double get_autocorrelation_time (double *x, int n);
extern double get_longest_autocorr_time_sample (double **s, int w, int n); 
extern int take_subserial_sample (double **s, int w, int n, int n_e);
extern double psi_max (const double data[], size_t n);
extern double psi_stats_mean (const double data[], size_t n);
extern double psi_stats_median (const double data[], size_t n);
extern double psi_stats_quantile (const double data[], size_t n, double f);
extern double psi_stats_sd (const double data[], size_t n);
extern double psi_stats_var (const double data[], size_t n);
extern double psi_stats_covariance (const double data1[], const double data2[], size_t n);

END_C_DECLS

#endif /* !PSI_GSLWRAP_H */
/** @end 1 **/


