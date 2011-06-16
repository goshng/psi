/* gslwrap.c -- GNU Scientific Library wrapper functions
   Copyright (C) 2000, 2006 Gary V. Vaughan
  
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
#include "gslwrap.h"

static const int debug_detail = 0;
/*
extern int sample_size;
extern char *program_invocation_short_name;
extern double tolerance_zeros;
extern int want_debug_bf;
extern int want_debug_gibbs;
extern gsl_rng *gmel_gsl_r;
extern FILE *ofile;
*/

static const double max_float = 3.40282347e+38F;
static const double log_0 = -3.40282347e+38F;
static const double log_limit = -3.40282347e+36F;
/* const double log_limit = -max_float/100; */
/* const double log_0 = -max_float; */

// For Pentium I (exp & log are base 2)
// fadd = 3,1
// fmul = 3,1
// f2xm1 = 13-57
// fyl2x = 22-111
// fdiv  = 39
// [ and fscale... ?]

// NATS is 52*log(2) for 52 bits of precision
// HMM... the long doubles have 64 bits of precision...
static const double NATS = 40;

int psi_int_function (int a)
{
  return a * a;
}

double psi_function (double a)
{
  printf ("in func: %lf\n", a);
  return a;
}

/* BUGGY: x might be better to be moved to the first argument */
double 
psi_gaussian_lnpdf (double mu, double sigma, double x)
{
  double lnpdf = (-0.5) * log (2.0 * M_PI)                      
                 - log (sigma)                    
                 - (x - mu)*(x - mu) / (2.0 * sigma * sigma);
  return lnpdf;
}

/* Note: the output is in logarithm scale */
double
psi_sum_exp (double *log_i, int n)
{
  int i;
  double log_sum = 0;

  assert (log_i != NULL);

  log_sum = *log_i; 
  log_i++;
  for(i = 1; i < n; i++, log_i++) {
    log_sum = logsum (*log_i, log_sum);
  }

  return log_sum;; 
}

void
psi_print_double_array (double *A_i, int n)
{
  int i;

  for (i = 0; i < n; i++) {
     fprintf (stderr, "%lf ", A_i[i]);
  }
  fprintf (stderr, "\n");
}

double 
logsum_nocheck(double x, double y) {
  if (abs(x-y) > NATS)
    return ((x > y) ? x : y);
  else
    return (x + log1p(exp(y - x)));
}

double 
logsum(double x, double y)
{
  double temp = y-x;
  if (temp > NATS || x < log_limit)
    return y;
  else if (temp < -NATS || y < log_limit)
    return x;
  else
    return (x + log1p(exp(temp)));
}

void 
loginc(double *x, double y)
{
  double temp = y-*x;
  if (temp > NATS || *x < log_limit)
    *x=y;
  else if (temp < -NATS || y < log_limit)
    ;
  else
    *x += log1p(exp(temp));
}

double 
logdiff (double x, double y) 
{
  assert(x > y);
  double temp = y-x;
  if (temp < -NATS || y < log_limit)
    return x;
  else
    return (x + log1p(-exp(temp)));
}

double 
psi_pi_lnpdf (const double alpha[], double theta[])
{
  theta[3] = 1.0 - theta[0] - theta[1] - theta[2]; 
  double u = gsl_ran_dirichlet_lnpdf (4, alpha, theta);
  return u;
}


int
gmel_conv_psr (double *psr, double **draws, int m, int n)
{
  int r = EXIT_SUCCESS;

  double B = 0;
  double W = 0;
  double var_plus = 0;
  double m_d_d = 0;
  double *m_d_j = NULL;
  double *m_2_d_j = NULL;
  double *s_j = NULL;
  double d_f = 0; 
  double V_hat = 0;
  double var_hat_V_hat = 0;
  double var_hat_s = 0;
  double cov_hat_s_x_2 = 0;
  double cov_hat_s_x = 0;
  double dt_m = (double) m;
  double dt_n = (double) n;
  int i, j;
  m_d_j = (double *) malloc (sizeof (double) * m);
  m_2_d_j = (double *) malloc (sizeof (double) * m);
  s_j = (double *) malloc (sizeof (double) * m);

  m_d_d = 0;
  for (j = 0; j < m; j++) {
     m_d_j[j] = 0;
     for (i = 0; i < n; i++) {
        m_d_j[j] += draws[j][i];
     }
     m_d_j[j] /= dt_n;
     m_2_d_j[j] = (m_d_j[j] * m_d_j[j]);
     m_d_d += m_d_j[j];
  }
  m_d_d /= dt_m;

  B = 0;
  for (j = 0; j < m; j++) {
     B += ((m_d_j[j] - m_d_d) * (m_d_j[j] - m_d_d));
  }
  B *= dt_n;
  B /= (dt_m - 1.0); /* ---------> B has been calculated */
   
  W = 0;
  for (j = 0; j < m; j++) {
     s_j[j] = 0;
     for (i = 0; i < n; i++) {
        s_j[j] += (draws[j][i] - m_d_j[j]) * (draws[j][i] - m_d_j[j]);
     }
     s_j[j] /= (dt_n - 1.0); 
     W += s_j[j];
  }
  W /= dt_m;       /* ---------> W has bben calculated */

  var_plus = ((dt_n - 1) * W / dt_n) + (B / dt_n);
  V_hat = var_plus + (B / (dt_m * dt_n));

  get_sample_variance (m, &var_hat_s, s_j, s_j);
  get_sample_variance (m, &cov_hat_s_x, s_j, m_d_j);
  get_sample_variance (m, &cov_hat_s_x_2, s_j, m_2_d_j);
  var_hat_V_hat = ((((dt_n - 1)*(dt_n - 1))/(dt_n * dt_n))/dt_m)*var_hat_s
                  + ((((dt_m + 1)*(dt_m + 1))/(dt_m * dt_m * dt_n * dt_n))*2/(dt_m - 1.0)) * B * B /* BUGGY: need to check */
                  + (2 * (dt_m + 1.0) * (dt_n - 1.0) / (dt_m * dt_n * dt_n))
                  * (dt_n / dt_m)*(cov_hat_s_x_2 - 2 * m_d_d * cov_hat_s_x);
  /* see the ref: Gelman and Rubin 1992 */

  d_f = 2 * V_hat / var_hat_V_hat;             /* ref: Gelman & Rubin 1992 */
  *psr = sqrt (((d_f + 3) * V_hat) / ((d_f + 1) * W));
  /* OLD ONE: *psr = sqrt(var_plus / W); */

/*
fprintf (stderr, "m_d_d: %e\n", m_d_d);
fprintf (stderr, "var_plus: %e\n", var_plus);
fprintf (stderr, "W: %e\n", W);
fprintf (stderr, "B: %e\n", B);
fprintf (stderr, "var_hat_s: %e\ncov_hat_s_x: %e\ncov_hat_s_x_2: %e\n",
         var_hat_s, cov_hat_s_x, cov_hat_s_x_2);
fprintf (stderr, "var_hat_V_hat: %e\n", var_hat_V_hat);
fprintf (stderr, "V_hat: %e\nd_f: %e\n", V_hat, d_f);
fprintf (stderr, "psr: %lf\n\n", *psr);
*/

  free (s_j);
  s_j = NULL;
  free (m_d_j);
  m_d_j = NULL;
  free (m_2_d_j);
  m_2_d_j = NULL;

  return r;
}

int 
get_sample_variance (int n, double *v, double *x, double *y)
{
  int i;
  double m_x = 0;
  double m_y = 0;

  m_x = 0;
  for (i = 0; i < n; i++) {
     m_x += x[i];
  }
  m_x /= n;

  m_y = 0;
  for (i = 0; i < n; i++) {
     m_y += y[i];
  }
  m_y /= n;
  
  *v = 0;
  for (i = 0; i < n; i++) {
     *v += ((x[i] - m_x) * (y[i] - m_y));
  }
  *v /= (n - 1);

  return EXIT_SUCCESS;
}

  /* Ask Ben before coding this part */  
  /*
     N_e = N * (1 - r) / (1 + r)
     N / N_e = tau
     tau = (1 + r) / (1 - r)
  double tau, r;
  r = gsl_stats_lag1_autocorrelation (x, 1, n);
  tau = (1.0 + r) / (1.0 - r);
  */
double 
get_autocorrelation_time (double *x, int n)
{
  double V, sum, tau;
  int i;
  double *rho = NULL;
  int n_e;
  int n_c = 2 + n/4;
  rho = XMALLOC (double, n_c);
/* for (i = 0; i < n; i++) fprintf (stderr, "%lf ", x[i]); */
/*
  fprintf (stderr, "mean: %lf\n", gsl_stats_mean (x, 1, n));
*/

  n_e = n_c;
  for (i = 0; i < n_c; i++)
    {
      rho[i] = gsl_stats_covariance (&x[0], 1, &x[i], 1, n - i);
      if (rho[i] <= 0.0) 
        {
          n_e = i;
          break;
        } 
    } 
/*
  if (i == n_c)
    n_e = i;
*/
  V = rho[0];
  sum = 0.0;
  for (i = 1; i < n_e; i++)
    {
      sum += rho[i];
    }
  tau = 1.0 + 2.0 * sum / (V);
/*
fprintf (stderr, "tau: %lf, n_e: %d\n", tau, n_e);
*/

  XFREE (rho);
  return tau;
}

double 
get_longest_autocorr_time_sample (double **s, int w, int n)
{
  double tau;
  double *taus = NULL;
  int i;

  assert (w > 0);
  assert (n > 0);
  
  taus = XMALLOC (double, w);
  for (i = 0; i < w; i++)
    {
      taus[i] = get_autocorrelation_time (s[i], n);
    }
  tau = gsl_stats_max (taus, 1, w);
  XFREE (taus);

  return tau;
}

int 
take_subserial_sample (double **s, int w, int n, int n_e)
{
  int i, j, k;
  double tau;

  assert (w > 0);
  assert (n > 0);
  assert (n_e > 0);
  assert (n_e <= n);

  tau = ((double) n) / ((double) n_e);
  for (i = 0; i < w; i++)
    {
      for (j = 0; j < n_e; j++)
        {
          k = (int) (j * tau);
          assert (k < n);
          s[i][j] = s[i][k];
        }
    }
 
  return EXIT_SUCCESS;
}

double
psi_max (const double data[], size_t n)
{
  return gsl_stats_max (data, 1, n);
}

double
psi_stats_mean (const double data[], size_t n)
{
  return gsl_stats_mean (data, 1, n);
}

double 
psi_stats_median (const double data[], size_t n)
{
  size_t i;
  double v;
  double *sorted = XMALLOC (double, n);
  for (i = 0; i < n; i++)
    sorted[i] = data[i];
  gsl_sort (sorted, 1, n);
  v = gsl_stats_median_from_sorted_data (sorted, 1, n);
  XFREE (sorted);
  return v;
}

double 
psi_stats_quantile (const double data[], size_t n, double f)
{
  size_t i;
  double v;
  double *sorted = XMALLOC (double, n);
  for (i = 0; i < n; i++)
    sorted[i] = data[i];
  gsl_sort (sorted, 1, n);
  v = gsl_stats_quantile_from_sorted_data (sorted, 1, n, f);
  XFREE (sorted);
  return v;
}

double 
psi_stats_sd (const double data[], size_t n)
{
  double v;
  v = gsl_stats_variance (data, 1, n);
  return sqrt (v);
}

double 
psi_stats_var (const double data[], size_t n)
{
  return gsl_stats_variance (data, 1, n);
}

double 
psi_stats_covariance (const double data1[], const double data2[], size_t n)
{
  return gsl_stats_covariance (data1, 1, data2, 1, n);
}
