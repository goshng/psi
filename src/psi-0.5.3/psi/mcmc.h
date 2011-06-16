/* mcmc.h -- MCMC procedure of the single sequence analysis
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

/*!\file mcmc.h
   \author Sang Chul Choi
   \brief MCMC - propose a new value

 */

/** @start 1 **/
#ifndef PSI_MCMC_H
#define PSI_MCMC_H 1

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <psi/common.h>

BEGIN_C_DECLS

enum { 
  PSI_PI_A_MODE,
  PSI_PI_C_MODE,
  PSI_PI_G_MODE,
  PSI_PI_T_MODE,
  PSI_PI_D_MODE, /* Using Dirichlet Density */
  PSI_PI_MODE,
  PSI_S_MODE,
  PSI_S_MODE_NORMAL, /* Normal Prior */
  PSI_P_MODE,
  PSI_P_MODE_NORMAL, /* Normal Prior */
  PSI_D_MODE
};

/*
extern GibbsPart *GibbsInfo[SIZEOFGRID][SIZEOFGRID][FREQGRIDSIZE][FREQGRIDSIZE][FREQGRIDSIZE];
extern double S_POINTS[SIZEOFGRID];
extern double P_POINTS[SIZEOFGRID];
extern double A_POINTS[FREQGRIDSIZE];
extern double G_POINTS[FREQGRIDSIZE];
extern double C_POINTS[FREQGRIDSIZE];
extern double T_POINTS[FREQGRIDSIZE];
*/

/*!\brief Store number of nucleotides, solv and pair

    This will remove the uses of Data, possibly Interaction structures
 */
extern int access_mcmc_sample_size ();
extern int access_mcmc_sample_freq ();
extern int access_mcmc_sample_burn ();
extern int unsetup_mcmc ();
extern int setup_mcmc_data ();
extern int setup_mcmc_prior (double s_b, double s_e, double p_b, double p_e, 
                             double n_b, double n_e);
extern int setup_mcmc_grid (double s_b, double s_e, int s_n, 
                            double p_b, double p_e, int p_n,
                            double n_b, double n_e, int n_n);
extern int setup_mcmc_delta (double s, double p, double n);
extern int setup_mcmc_init (double s, double p, double a, double c, double g);
extern int setup_mcmc_file (FILE *fp, char *pn, char *an);
extern int setup_mcmc_want (int no, int s, int p, int n, int d, int debug, int lf);
extern int setup_mcmc_sample (int s, int f, int b);
extern int setup_mcmc_limit (double t, double p);
extern int setup_mcmc_state (int s, int l);

/*!\brief Sample \f$\pi\f$
 
   This function samples \f$pi\f$ from posterior distribution. 

   \param curr current parameter values
   \param energy energy information
   \param dat 
   \param whichnuc
 */
extern int new_pi (parameter *p_theta);

/*!\brief Sample \f$s\f$
 */
extern int new_s (parameter *p_theta);

/*!\brief Sample \f$p\f$
 */
extern int new_p (parameter *p_theta);

/*!\brief Propose \f$theta^*\f$ given \f$\theta\f$

   Like most of the codes, this functions is almost the same as Doug's
   original code. This function is another name of  NewParameter of DrEvol.
 */
extern int get_theta_prime (parameter *theta, 
                            parameter *p_theta_prime,
                            int mode); 

/*!\brief Calculate Jumping Term of r

   It calcuates jumping density term: 
   \f$\frac{J(\theta',\theta}{J(\theta,\theta')}\f$. It returns its 
   logarithm.
 */
extern int get_jumping (double *p_term, 
                        parameter theta, parameter theta_prime,
                        int mode); 

/*!\brief Calculate Prior Term of r

   It calculates prior density term: 
   \f$\frac{p(\theta'}{p(\theta)}\f$. It returns its logarithm.
 */
extern int get_prior (double *p_term, 
                      parameter theta, parameter theta_prime,
                      int mode); 

/*!\brief Calculate Likelihood Term of r

   It calculates likelihood density term: 
   \f$\frac{p(i|\theta'}{p(i|\theta)}\f$. It returns its logarithm.
 */
extern int get_likelihood (double *p_term, 
                           parameter *p_theta, parameter *p_theta_prime);

/*!\brief Locate the nearest grid-point of two thetas

   For given \f$\theta\f$ and \f$\theta^*\f$, it locates the nearest 
   grid-point.
 */ 
extern int locate_gridpoint (parameter *p_theta_star,
                             parameter *p_theta, 
                             parameter *p_theta_prime);

/*!\brief 

   We will calcualte 
   \f$\sum_{h=1}^M e^{-2(s-s^*)E_s(\eta^{(h)}) - 2(p-p^*)E_p(\eta^{(h)}) }
   \prod_{n=1}^N \frac{\pi_{\eta_{n}^{(h)}}}{\pi_{\eta_{n}^{(h)}}^*}\f$
extern int sumup_gibbs (double *sum, double *maximum,
                           parameter theta, parameter theta_prime, 
                           SAMPLE_MODE mode, 
                           Interaction *FO, dat_t *dat);
 */

/*!\brief 

   Now, let's sample DNA sequences using Gibbs sampler!
extern int gibbs_sampler (Interaction *FO, dat_t *dat, 
                             int solv, int pair, int A, int C, int G);
 */

/*!\brief 

   It's an energy calculator, which will be replaced by Sang Chul's version.
extern void SimSeqNRG (Interaction *FO, int *tmpAAseq, GibbsPart *gibbsSeq, 
                          int AAlen);
 */

/*!\brief Simple Function Group
 */
extern double Maxof3 (double a, double b, double c);

/*!\brief Simple Function Group
 */
extern int create_gridpoints (int n); 

/*!\brief Load Gibbs Sample
 */
extern int load_gridpoints (const char *fn); 

/*!\brief Save Gibbs Sample
 */
extern int save_gridpoints (const char *fn); 

/*!\brief Simple Function Group
 */
extern int delete_gridpoints ();

/*!\brief Simple Function Group
 */
extern int init_theta (parameter *p_theta);
extern int psi_mcmc_init_theta_gp_gibbs_sum (parameter *p_theta);

/*!\brief Check over- or under- boundary 
 */
extern int correct_boundary ( int *gp, int max );

/*!\brief GSL random number generator initialization
extern int gmel_gsl_init ();
 */

/*!\brief GSL random number finalization
extern int gmel_gsl_fin ();
 */

/*!\brief Print acceptance ratio
 */
extern int write_acceptance_ratio ();


/*!\brief For guarding the crazy Dirichlet proposal
   
   set the value greater than or equal to 1.0

   \param a double-array of size 4
extern int get_sensible_pi (double *pi);
 */

extern int execute_estimation (int is_fixed_s, int is_fixed_p, 
                               double s, double p,
                               parameter *theta_star);
extern int run_mcmc_chain (double **posterior, int b, int f, int s, int o);
extern int find_mcmc_option ();
extern int run_mcmc_multichain (int number_chains, double **ip, 
                                int is_fixed_run);
extern int gmel_bf_s_sample (double s1, double *s2);
extern int gmel_bf_p_sample (double p1, double *p2);
extern int gmel_bf_n_sample (double a1, double c1, double g1, 
                             double *a2, double *c2, double *g2, int *w);

extern double psi_mcmc_s_flat_proposal_pdf (double x, double mu_s);
extern double psi_mcmc_p_flat_proposal_pdf (double x, double mu_p);
extern double psi_mcmc_s_normal_proposal_lnpdf (double x, double mu_s);
extern double psi_mcmc_p_normal_proposal_lnpdf (double x, double mu_p);
extern double psi_mcmc_s_flat_prior_lnpdf (double s);
extern double psi_mcmc_p_flat_prior_lnpdf (double p);
extern double psi_mcmc_s_normal_prior_lnpdf (double s);
extern double psi_mcmc_p_normal_prior_lnpdf (double p);
extern double psi_mcmc_n_dirichlet_proposal_lnpdf (double alpha[], double theta[]);
extern double psi_mcmc_n_likelihood (double a, double c, double g, double t);

extern void initialize_mcmc_acceptance ();
extern void write_mcmc_info ();



END_C_DECLS

#endif /* !PSI_MCMC_H */
/** @end 1 **/
