/* bf.h -- Bayes factor estimation
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
/*!\file
   \author Sang Chul Choi
   \brief Estimating Bayes Factor

   \example 
   gmel_bf_prior_set(3,1); FIX : this can be called inside the file
   gmel_bf(3,1);
 */

/** @start 1 **/
#ifndef PSI_BF_H
#define PSI_BF_H 1

#include <psi/common.h>

#define ERR_BF_MEMORY         1000
#define ERR_BF_INVDNACHR      1001
#define ERR_BF_INVPROCHR      1002
#define ERR_BF_NINF           1003
#define ERR_BF_NAN            1004
#define ERR_BF_NOTFINITE      1005
#define ERR_BF_NOTIMPLEMENTED 1006
#define ERR_BF_FILE           1007
#define ERR_BF_READ           1008
#define ERR_BF_GRID           1009
#define ERR_BF_ASSERT         1010
#define ERR_BF_MATH           1011
#define ERR_BF_TOOMANYZEROS   1012
#define ERR_BF_LOGZERO        1013
#define ERR_BF_SAME           1014
#define ERR_BF_INVBASE        1015
#define ERR_BF_BLOCK          1016
#define ERR_BF_GIBBS          1017
#define ERR_BF_INVARG         1018
#define ERR_BF_TRANSLATION    1019
#define ERR_BF_UNDEFINED_Q    1020

BEGIN_C_DECLS

extern const char *program_name;

enum { 
  MODEL1_BLOCK1, 
  MODEL4_BLOCK1, 
  MODEL4_BLOCK2, 
  MODEL4_BLOCK3 
};

enum {
  model_top,
  model_bot
};

extern int execute_bayesfactor ();
extern int unsetup_bf ();
extern int setup_bf_file (FILE *fp, const char *pn, const char *an, 
                          const char *qn, const char *bn);
/* extern int setup_bf_want (int ls, int lp, int d); */
extern int setup_bf_want (int lp, int d);
extern int setup_bf_tol (double tolerate_s_ratio);

/* might be removed */
extern int setup_bf_prior (double s_b, double s_e, double p_b, double p_e, 
                           double n_b, double n_e);
extern int setup_bf_delta (double s, double p, double n);
extern int setup_bf_state (int s);

/******************************************************************************
 BAYES FACTOR USER LEVEL Functions
******************************************************************************/

/*!\defgroup bf_user Bayes Factor - User Level Functions
 */

/*@{*/

/*!\brief Calculates Bayes Factor
   
   Bayes Factor is a ratio of two marginal likelihoods;
   \f$BF_{41} = \frac{P(Data|Model 4)}{P(Data|Model 1}\f$. For my
   conveniency, model 4 is called model top, and model 1 is called model 
   bottom.  We have two models

   - Model 1: s=0 & p=0, \f$\pi\f$ is a free parameter
   - Model 4: s, p and \f$\pi\f$ are free parameters

   \param bf the Bayes Factor in logarithm 
   \param theta_star_top the parameter values that achieve the high posterior
   \param theta_star_bot the parameter values that achieve the high posterior
   \param bfp Bayes factor output
   \return 0 for success, error code otherwise 
 */ 
extern int 
gmel_dbf (double *bf, parameter theta_star_top, 
                                parameter theta_star_bot, FILE *bfp);

/*!\brief Calculates Likelihood Term of Bayes Factor
  
   Bayes Factor is a product of three terms: likelihood ratio, prior ratio, 
   and posterior ratio.  This function calculates the first ratio, or 
   likelihood ratio. For example, Bayes Factor of model 3 over model 1 is
   \f$BF_{31} = P(i|M_3)/P(i|M_1)
              = \frac{P(i|M_3,s^{m^3},\pi^{m^3})}{P(i|M_1,\pi^{m^1})} \times
                \frac{P(s^{m^3},\pi^{m^3})}{P(\pi^{m^1})} \times
                \frac{P(\pi^{m^1}|i)}{P(s^{m^3},\pi^{m^3}|i)}\f$
   The first ratio is the target of this function.

   \param likelihood the returning value in logarithm scale
   \param model_top model number of numerator of Bayes Factor
   \param model_bot model number of numerator of Bayes Factor
   \param p_top parameter of top model
   \param p_bot parameter of bottom model
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_likelihood (double *likelihood, int model_top, int model_bot,
                              parameter p_top, parameter p_bot);

/*!\brief Calculates Likelihood Term of Bayes Factor Using Gibbs Sampler for Model 2 or 4 Cases
  
   Bayes Factor is a product of three terms: likelihood ratio, prior ratio, 
   and posterior ratio.  This function calculates the first ratio, or 
   likelihood ratio. For example, Bayes Factor of model 3 over model 1 is
   \f$BF_{31} = P(i|M_3)/P(i|M_1)
              = \frac{P(i|M_3,s^{m^3},\pi^{m^3})}{P(i|M_1,\pi^{m^1})} \times
                \frac{P(s^{m^3},\pi^{m^3})}{P(\pi^{m^1})} \times
                \frac{P(\pi^{m^1}|i)}{P(s^{m^3},\pi^{m^3}|i)}\f$
   The first ratio is the target of this function.

   \param likelihood the returning value in logarithm scale
   \param p_top parameter of top model
   \param p_bot parameter of bottom model
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_likelihood_gibbs (double *likelihood, 
                                    parameter p_top, parameter p_bot);



/*!\brief Calculates Prior Term of Bayes Factor
  
   Bayes Factor is a product of three terms: likelihood ratio, prior ratio, 
   and posterior ratio.  This function calculates the second ratio, or 
   prior ratio. For example, Bayes Factor of model 3 over model 1 is
   \f$BF_{31} = P(i|M_3)/P(i|M_1)
              = \frac{P(i|M_3,s^{m^3},\pi^{m^3})}{P(i|M_1,\pi^{m^1})} \times
                \frac{P(s^{m^3},\pi^{m^3})}{P(\pi^{m^1})} \times
                \frac{P(\pi^{m^1}|i)}{P(s^{m^3},\pi^{m^3}|i)}\f$
   The second ratio is the target of this function.

   \param prior the returning value in logarithm scale
   \param model_top model number of numerator of Bayes Factor
   \param model_bot model number of numerator of Bayes Factor
   \param p_top parameter of top model
   \param p_bot parameter of bottom model
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_prior (double *prior, int model_top, int model_bot,
                         parameter p_top, parameter p_bot);

/*!\brief Calculates Posterior of Bayes Factor
  
   Bayes Factor is a product of three terms: likelihood ratio, prior ratio, 
   and posterior ratio.  This function calculates the second ratio, or 
   posterior ratio. For example, Bayes Factor of model 3 over model 1 is
   \f$BF_{31} = P(i|M_3)/P(i|M_1)
              = \frac{P(i|M_3,s^{m^3},\pi^{m^3})}{P(i|M_1,\pi^{m^1})} \times
                \frac{P(s^{m^3},\pi^{m^3})}{P(\pi^{m^1})} \times
                \frac{P(\pi^{m^1}|i)}{P(s^{m^3},\pi^{m^3}|i)}\f$
   The third ratio is the target of this function.

   \param posterior the returning value in logarithm scale
   \param model_top model number of numerator of Bayes Factor
   \param model_bot model number of numerator of Bayes Factor
   \param p_top parameter of top model
   \param p_bot parameter of bottom model
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior (double *posterior, int model_top, int model_bot,
                             parameter p_top, parameter p_bot);

/*!\brief Make a new grid points
 
   This is a just a wrapper of gmel_grid module's gmel_grid_new function.
   The following explanation is from that gmel_grid_new function:

   The module, or this file, has static variables called gmel_grid_gridpoints
   and gmel_grid_info. The fist variable spcifies which floating-point number
   corresponds to which grid-point. The second one stores a list of 
   Gibbs-sampled sequence information for each grid-point. That is whay we
   have one more dimension: five dimensions for five parameters
   (S, P, A, C, G) and the last dimension for an array of information of 
   DNA sequence sampled by Gibbs Sampler.

   \param min_s minimum of s parameter
   \param max_s maximum of s parameter 
   \param n_s number of desired gridpoints
   \param min_p minimum of p parameter
   \param max_p maximum of p parameter 
   \param n_p number of desired gridpoints
   \param min_a minimum of a parameter
   \param max_a maximum of a parameter 
   \param n_a number of desired gridpoints
   \param min_c minimum of c parameter
   \param max_c maximum of c parameter 
   \param n_c number of desired gridpoints
   \param min_g minimum of g parameter
   \param max_g maximum of g parameter 
   \param n_g number of desired gridpoints
   \return 0 for SUCCESS 
 */
extern int 
gmel_bf_set_grid (double min_s, double max_s, int n_s,
                            double min_p, double max_p, int n_p,
                            double min_a, double max_a, int n_a,
                            double min_c, double max_c, int n_c,
                            double min_g, double max_g, int n_g);

/*!\brief Deallocate memory used for two static variables

   \return 0 for SUCCESS 
 */
extern int 
gmel_bf_del_grid ();

/*@}*/

/******************************************************************************
 SAMPLING Functions
******************************************************************************/

/*!\defgroup bf_sampling Bayes Factor - Sampling Functions
 */

/******************************************************************************
 MODEL 1 Functions
******************************************************************************/

/*!\defgroup bf_model1 Bayes Factor - Model 1
 */

/*@{*/

/*!\brief Calculate likelihood, or stationary density values of model 1 

   The model 1 has the following stationary distribution density fuction:
   \f$P(i|M_1,\pi) = {\displaystyle \prod_{m=1}^N\pi_{i_m} 
                      \over
                      \displaystyle 
                      \sum_{k\in\mathcal{K}}\prod_{n=1}^N\pi_{k_m}}\f$

   \param likelihood the return value of this function
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_likelihood_m1 (double *likelihood, 
                                 double a, double c, double g);

/*!\brief Calculate prior density values of model 1

   We are mostly using uniform density for a prior. This function returns
   the value of prior. 

   \param log_val the return value of this function in logarithmic scal
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_prior_m1 (double *prior, 
                            double a, double c, double g);

/*!\brief Approximate posterior density values of model 1 using Chib's idea

   The key of this whole program is going to start with this function. 
   We are going to take advantage of Chib's idea to approximate this last
   term, posterior density approximation.

   \param log_val the return value of this function in logarithmic scal
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m1 (double *log_val, 
                                double a, double c, double g);

/*!\brief Approximate posterior density's first block  values of model 1 
          using Chib's idea

   The last term of marginal distribution is the posterior term, or
   \f$P(\pi^{m^1}|i)\f$. This probability cannot be splitted into any more,
   and we have only one block.
   \f$P(\pi^{m^1}|i) = \frac{ E_1\{\alpha(\pi,\pi^{m^1}|i)q(\pi,\pi^{m^1}|i)\} }
	                  { E_2\{\alpha(\pi^{m^1},\pi|i)\} }\f$

   \param log_val the return value of this function in logarithmic scal
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m1_block1 (double *log_val, 
                                       double a, double c, double g);

/*!\brief Approximate numerator of posterior density's first block values 
          of model 1 

   Chib's paper said that the numerator is approximated with
   \f$E_1 \approx M^{-1}\sum_{m=1}^M \alpha(\pi^{(g)},\pi^{m^1}|i)
	                               q(\pi^{(g)},\pi^{m^1}|i)\f$

   \param log_val the return value of this function in logarithmic scal
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m1_block1_e1
   (double *log_val, double a, double c, double g);

/*!\brief Approximate denominator of posterior density's first block values 
          of model 1 

   Chib's paper said that the denominator is approximated with
   \f$E_2 \approx J^{-1}\sum_{j=1}^J \alpha(\pi^{m^1}, \pi^{(j)}|i)\f$

   \param log_val the return value of this function in logarithmic scal
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m1_block1_e2
   (double *log_val, double a, double c, double g);

/*!\brief Calcualte probability move function of numerator of posterior 
          density's first block of model 1

   This function calculates the probability move function:
   \f$\alpha(\pi^{(g)},\pi^{m^1}|i)\f$

   \param log_val the return value of this function in logarithmic scal
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m1_block1_alpha
   (double *alpha, double ag, double cg, double gg,
                             double a, double c, double g, int w);

/*!\brief Calcualte jumping density ordinate of numerator of posterior 
          density's first block of model 1 

   This function calculates the probability move function:
   \f$q(\pi^{(g)},\pi^{m^1}|i)\f$

   \param log_val the return value of this function in logarithmic scal
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m1_block1_q
   (double *q, double ag, double cg, double gg,
                         double a, double c, double g, int w);

/*@}*/

/******************************************************************************
 MODEL 3 Functions
******************************************************************************/

/*!\defgroup bf_model3 Bayes Factor - Model 3
 */

/*@{*/

/*!\brief Calculate likelihood, or stationary density values of model 3

   The stationary distribution of Doug's MBE paper is almost impossible to
   calculate without Gibbs Sampling approximation when pairwise interaction
   parameter, \f$p\f$, is involved in that stationary distribution. 
   Model 3, however, assumes that \f$p=0\f$ which is great beneficial to 
   easy calculation of the stationary distribution's function value which
   would have been impossible if \f$p\f$ was not equal to zero.

   \param likelihood the return value of this function
   \param s a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_likelihood_m3 (double *likelihood, 
                                 double s, double a, double c, double g);

/*!\brief Calculate prior density values of model 3

   We are mostly using uniform density for a prior. This function returns
   the value of prior. 

   \param log_val the return value of this function in logarithmic scal
   \param p flat prior information
   \param s a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_prior_m3 (double *prior, 
                            double s, double a, double c, double g);

/*!\brief Approximate posterior density values of model 3 using Chib's idea

   The key of this whole program is going to start with this function. 
   We are going to take advantage of Chib's idea to approximate this last
   term, posterior density approximation.

   \param log_val the return value of this function in logarithmic scale
   \param s a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3 (double *log_val, 
                                double s, double a, double c, double g);

/*!\brief Approximate posterior density's first block values of model 3

   The last term of marginal distribution is the posterior term, or
   \f$P(s^{m^3},\pi^{m^3}|i) = P(s^{m^3}|i)P(\pi^{m^3}|i,s^{m^3})\f$.
   The first term of right hand side of the equation is the target of this
   function.

   \param log_val the return value of this function in logarithmic scale
   \param s a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3_block1 (double *log_val, double s);

/*!\brief Approximate posterior density's second block values of model 3

   The last term of marginal distribution is the posterior term, or
   \f$P(s^{m^3},\pi^{m^3}|i) = P(s^{m^3}|i)P(\pi^{m^3}|i,s^{m^3})\f$.
   The second term of right hand side of the equation is the target of this
   function. This term is 
   \f$P(\pi^{m^3}|i,s^{m^3}) 
      = \frac{E_1\{ \alpha(\pi,\pi^{m^3}|i,s^{m^3})q(\pi,\pi^{m^3}|i,s^{m^3}) \}}
             {E_2\{ \alpha(\pi^{m^3},\pi|i,s^{m^3}) \}}\f$.

   \param log_val the return value of this function in logarithmic scale
   \param s a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3_block2
   (double *log_val, double s, double a, double c, double g);

/*!\brief Approximate numerator of posterior density's first block values 
          of model 3

   We are going to approximate the numerator (\f$E_1\f$) of full marginal:
   \f$P(s^{m^3}|i) = \frac{E_1\{ \alpha(s,s^{m^3}|i,\pi)q(s,s^{m^3}|i,\pi) \}}
	                  {E_2\{ \alpha(s^{m^3},s|i,\pi) \}}\f$.
   Chib's paper said that the numerator is approximated with
   \f$E_1 \approx M^{-1}\sum_{g=1}^M \alpha(s^{(g)},s^{m^3}|i,\pi^{(g)})
	                             q(s^{(g)},s^{m^3}|i,\pi^{(g)})\f$

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3_block1_e1
   (double *log_val, double s);

/*!\brief Approximate denominator of posterior density's first block values 
          of model 3

   We are going to approximate the denominator (\f$E_2\f$) of full marginal:
   \f$P(s^{m^3}|i) = \frac{E_1\{ \alpha(s,s^{m^3}|i,\pi)q(s,s^{m^3}|i,\pi) \}}
	                  {E_2\{ \alpha(s^{m^3},s|i,\pi) \}}\f$.
   Chib's paper said that the denominator is approximated with
   \f$E_2 \approx J^{-1}\sum_{j=1}^J \alpha(s^{m^3}, s^{(j)}|i,\pi^{(j)})\f$

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3_block1_e2
   (double *log_val, double s);

/*!\brief Calcualte probability move function of numerator of posterior 
          density's first block of model 3

   This function calculates the probability move function:
   \f$\alpha(s^{(g)},s^{m^3}|i,\pi^{(g)})\f$

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \param sg a parameter value of solven accessibility
   \param ag a parameter value of adenine frequency 
   \param cg a parameter value of cytosine frequency 
   \param gg a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3_block1_alpha
   (double *alpha, 
              double s, double sg, double ag, double cg, double gg);

/*!\brief Calcualte jumping density ordinate of numerator of posterior 
          density's first block of model 3

   This function calculates the probability move function:
   \f$q(s^{(g)},s^{m^3}|i,\pi^{(g)})\f$

   \param log_val the return value of this function in logarithmic scal
   \param s1 a parameter value of solven accessibility
   \param s2 a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3_block1_q
   (double *q, 
              double s, double sg, double ag, double cg, double gg);

/*!\brief Approximate numerator of posterior density's second block values 
          of model 3

   We are going to approximate the numerator (\f$E_1\f$) of full conditional:
   \f$P(\pi^{m^3}|i,s^{m^3}) 
	= \frac{E_1\{ \alpha(\pi,\pi^{m^3}|i,s^{m^3})q(\pi,\pi^{m^3}|i,s^{m^3}) \}}
	       {E_2\{ \alpha(\pi^{m^3},\pi|i,s^{m^3}) \}}\f$
   Chib's paper said that the numerator is approximated with
   \f$E_1 \approx J^{-1}\sum_{j=1}^J \alpha(\pi^{(j)},\pi^{m^3}|i,s^{m^3})
	                               q(\pi^{(j)},\pi^{m^3}|i,s^{m^3})\f$

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3_block2_e1
   (double *log_val, double s, double a, double c, double g);

/*!\brief Approximate denominator of posterior density's second block values 
          of model 3

   We are going to approximate the numerator (\f$E_1\f$) of full conditional:
   \f$P(\pi^{m^3}|i,s^{m^3}) 
	= \frac{E_1\{ \alpha(\pi,\pi^{m^3}|i,s^{m^3})q(\pi,\pi^{m^3}|i,s^{m^3}) \}}
	       {E_2\{ \alpha(\pi^{m^3},\pi|i,s^{m^3}) \}}\f$
   Chib's paper said that the numerator is approximated with
   \f$E_2 \approx K^{-1}\sum_{k=1}^K \alpha(\pi^{m^3}, \pi^{(k)}|i,s^{m^3})\f$

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3_block2_e2
   (double *log_val, double s, double a, double c, double g);

/*!\brief Calcualte probability move function of numerator of posterior 
          density's second block of model 3

   This function calculates the probability move function:
   \f$\alpha(\pi^{(j)},\pi^{m^3}|i,s^{m^3})\f$

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3_block2_alpha
   (double *alpha, 
              double aj, double cj, double gj, double a, double c, double g, 
              double s, int w);

/*!\brief Calcualte jumping density ordinate of numerator of posterior 
          density's second block of model 3

   This function calculates the probability move function:
   \f$q(\pi^{(j)},\pi^{m^3}|i,s^{m^3})\f$

   \param log_val the return value of this function in logarithmic scal
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \param aj a parameter value of adenine frequency 
   \param cj a parameter value of cytosine frequency 
   \param gj a parameter value of guanine frequency 
   \param w which nuc was used to proposed a new value
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m3_block2_q
   (double *q, 
              double aj, double cj, double gj, double a, double c, double g, 
              double s, int w);

/*@}*/

/******************************************************************************
 MODEL 2 Functions
******************************************************************************/

/*!\defgroup bf_model2 Bayes Factor - Model 2
 */

/*@{*/

/*!\brief Calculate likelihood, or stationary density values of model 2 

   The model 2 has the following stationary distribution density fuction:
   \f$P(i|M_2,p,\pi) = {\displaystyle e^{-2pE_p(i)} \prod_{m=1}^N\pi_{i_m} 
                        \over
                        \displaystyle 
                        \sum_{k\in\mathcal{K}}
                        e^{-2pE_p(k)}\prod_{n=1}^N\pi_{k_m}}\f$

   \param likelihood the return value of this function
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_likelihood_m2 (double *likelihood, 
                                 double p, double a, double c, double g);

/*!\brief Calculate prior density values of model 2 

   We are mostly using uniform density for a prior. This function returns
   the value of prior. 

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_prior_m2 (double *prior, 
                            double p, double a, double c, double g);

/*!\brief Approximate posterior density values of model 2 using Chib's idea

   The key of this whole program is going to start with this function. 
   We are going to take advantage of Chib's idea to approximate this last
   term, posterior density approximation.

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2 (double *log_val, 
                                double p, double a, double c, double g);

/*!\brief Approximate posterior density's first block values of model 2 

   The last term of marginal distribution is the posterior term, or
   \f$P(p^{m^2},\pi^{m^2}|i) = P(p^{m^2}|i)P(\pi^{m^2}|i,p^{m^2})\f$.
   The first term of right hand side of the equation is the target of this
   function.

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2_block1 
   (double *log_val, double p);

/*!\brief Approximate posterior density's second block values of model 2 

   The last term of marginal distribution is the posterior term, or
   \f$P(p^{m^2},\pi^{m^2}|i) = P(p^{m^2}|i)P(\pi^{m^2}|i,p^{m^2})\f$.
   The second term of right hand side of the equation is the target of this
   function. This term is 
   \f$P(\pi^{m^2}|i,p^{m^2}) 
      = \frac{E_1\{ \alpha(\pi,\pi^{m^2}|i,p^{m^2})q(\pi,\pi^{m^2}|i,p^{m^2}) \}}
             {E_2\{ \alpha(\pi^{m^2},\pi|i,p^{m^2}) \}}\f$.

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2_block2
   (double *log_val, double p, double a, double c, double g);

/*!\brief Approximate numerator of posterior density's first block values 
          of model 2

   We are going to approximate the numerator (\f$E_1\f$) of full marginal:
   \f$P(p^{m^2}|i) = \frac{E_1\{ \alpha(p,p^{m^2}|i,\pi)q(p,p^{m^2}|i,\pi) \}}
	                  {E_2\{ \alpha(p^{m^2},p|i,\pi) \}}\f$.
   Chib's paper said that the numerator is approximated with
   \f$E_1 \approx M^{-1}\sum_{g=1}^M \alpha(p^{(g)},p^{m^3}|i,\pi^{(g)})
	                             q(p^{(g)},p^{m^2}|i,\pi^{(g)})\f$

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2_block1_e1
   (double *log_val, double p);

/*!\brief Approximate denominator of posterior density's first block values 
          of model 2

   We are going to approximate the denominator (\f$E_2\f$) of full marginal:
   \f$P(p^{m^2}|i) = \frac{E_1\{ \alpha(p,p^{m^2}|i,\pi)q(p,p^{m^2}|i,\pi) \}}
	                  {E_2\{ \alpha(p^{m^2},p|i,\pi) \}}\f$.
   Chib's paper said that the denominator is approximated with
   \f$E_2 \approx J^{-1}\sum_{j=1}^J \alpha(p^{m^2}, p^{(j)}|i,\pi^{(j)})\f$

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2_block1_e2
   (double *log_val, double p);

/*!\brief Calcualte probability move function of numerator of posterior 
          density's first block of model 2

   This function calculates the probability move function:
   \f$\alpha(p^{(g)},p^{m^2}|i,\pi^{(g)})\f$

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param pg a parameter value of pairwise interaction
   \param ag a parameter value of adenine frequency 
   \param cg a parameter value of cytosine frequency 
   \param gg a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2_block1_alpha
   (double *prior, 
              double pg, double p, 
              double ag, double cg, double gg);

/*!\brief Calcualte jumping density ordinate of numerator of posterior 
          density's first block of model 2

   This function calculates the probability move function:
   \f$q(s^{(g)},s^{m^3}|i,\pi^{(g)})\f$

   \param log_val the return value of this function in logarithmic scal
   \param p1 a parameter value of pairwise interaction
   \param p2 a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2_block1_q
   (double *q, 
              double pg, double p, 
              double ag, double cg, double gg);

/*!\brief Approximate numerator of posterior density's second block values 
          of model 2

   We are going to approximate the numerator (\f$E_1\f$) of full conditional:
   \f$P(\pi^{m^2}|i,p^{m^2}) 
	= \frac{E_1\{ \alpha(\pi,\pi^{m^2}|i,p^{m^2})q(\pi,\pi^{m^2}|i,p^{m^2}) \}}
	       {E_2\{ \alpha(\pi^{m^2},\pi|i,p^{m^2}) \}}\f$
   Chib's paper said that the numerator is approximated with
   \f$E_1 \approx J^{-1}\sum_{j=1}^J \alpha(\pi^{(j)},\pi^{m^2}|i,p^{m^2})
	                               q(\pi^{(j)},\pi^{m^2}|i,p^{m^2})\f$

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2_block2_e1
   (double *log_val, double p, double a, double c, double g);

/*!\brief Approximate denominator of posterior density's second block values 
          of model 2

   We are going to approximate the numerator (\f$E_1\f$) of full conditional:
   \f$P(\pi^{m^2}|i,p^{m^2}) 
	= \frac{E_1\{ \alpha(\pi,\pi^{m^2}|i,p^{m^2})q(\pi,\pi^{m^2}|i,p^{m^2}) \}}
	       {E_2\{ \alpha(\pi^{m^2},\pi|i,p^{m^2}) \}}\f$
   Chib's paper said that the numerator is approximated with
   \f$E_2 \approx K^{-1}\sum_{k=1}^K \alpha(\pi^{m^2}, \pi^{(k)}|i,p^{m^2})\f$

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2_block2_e2
   (double *log_val, double p, double a, double c, double g);

/*!\brief Calcualte probability move function of numerator of posterior 
          density's second block of model 2

   This function calculates the probability move function:
   \f$\alpha(\pi^{(j)},\pi^{m^2}|i,p^{m^2})\f$

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2_block2_alpha
   (double *alpha, 
              double aj, double cj, double gj, double a, double c, double g, 
              double p, int w);

/*!\brief Calcualte jumping density ordinate of numerator of posterior 
          density's second block of model 2

   This function calculates the probability move function:
   \f$q(\pi^{(j)},\pi^{m^2}|i,p^{m^2})\f$

   \param log_val the return value of this function in logarithmic scal
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \param aj a parameter value of adenine frequency 
   \param cj a parameter value of cytosine frequency 
   \param gj a parameter value of guanine frequency 
   \param w which nuc was used to proposed a new value
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m2_block2_q
   (double *q, 
              double aj, double cj, double gj, double a, double c, double g, 
              double p, int w);

/*@}*/

/******************************************************************************
 MODEL 4 Functions
******************************************************************************/

/*!\defgroup bf_model4 Bayes Factor - Model 4
 */

/*@{*/

/*!\brief Calculate likelihood, or stationary density values of model 4 

   The model 4 has the following stationary distribution density fuction:
   \f$P(i|M_2,p,\pi) = {\displaystyle e^{-2pE_p(i)} \prod_{m=1}^N\pi_{i_m} 
                        \over
                        \displaystyle 
                        \sum_{k\in\mathcal{K}}
                        e^{-2pE_p(k)}\prod_{n=1}^N\pi_{k_m}}\f$

   \param likelihood the return value of this function
   \param s a parameter value of solven accessibility
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_likelihood_m4 (double *likelihood, 
                                 double s, double p, 
                                 double a, double c, double g);

/*!\brief Calculate prior density values of model 4 

   We are mostly using uniform density for a prior. This function returns
   the value of prior. 

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_prior_m4 (double *prior, 
                            double s, double p, double a, double c, double g);

/*!\brief Approximate posterior density values of model 4 using Chib's idea

   The key of this whole program is going to start with this function. 
   We are going to take advantage of Chib's idea to approximate this last
   term, posterior density approximation.

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4 (double *posterior, 
                                double s, double p, 
                                double a, double c, double g);

/*!\brief Approximate posterior density's first block values of model 4 

   The last term of marginal distribution is the posterior term, or
   \f$P(p^{m^2},\pi^{m^2}|i) = P(p^{m^2}|i)P(\pi^{m^2}|i,p^{m^2})\f$.
   The first term of right hand side of the equation is the target of this
   function.

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block1 
   (double *log_val, double s);

/*!\brief Approximate posterior density's second block values of model 4 

   The last term of marginal distribution is the posterior term, or
   \f$P(p^{m^2},\pi^{m^2}|i) = P(p^{m^2}|i)P(\pi^{m^2}|i,p^{m^2})\f$.
   The second term of right hand side of the equation is the target of this
   function. This term is 
   \f$P(\pi^{m^2}|i,p^{m^2}) 
      = \frac{E_1\{ \alpha(\pi,\pi^{m^2}|i,p^{m^2})q(\pi,\pi^{m^2}|i,p^{m^2}) \}}
             {E_2\{ \alpha(\pi^{m^2},\pi|i,p^{m^2}) \}}\f$.

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \param p a parameter value of pairwise interaction
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block2
   (double *log_val, double s, double p);

/*!\brief Approximate posterior density's third block values of model 4 

   The last term of marginal distribution is the posterior term, or
   \f$P(p^{m^2},\pi^{m^2}|i) = P(p^{m^2}|i)P(\pi^{m^2}|i,p^{m^2})\f$.
   The second term of right hand side of the equation is the target of this
   function. This term is 
   \f$P(\pi^{m^2}|i,p^{m^2}) 
      = \frac{E_1\{ \alpha(\pi,\pi^{m^2}|i,p^{m^2})q(\pi,\pi^{m^2}|i,p^{m^2}) \}}
             {E_2\{ \alpha(\pi^{m^2},\pi|i,p^{m^2}) \}}\f$.

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \param p a parameter value of pairwise interaction
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block3
   (double *log_val, 
              double s, double p, double a, double c, double g);

/*!\brief Approximate numerator of posterior density's first block values 
          of model 4

   We are going to approximate the numerator (\f$E_1\f$) of full marginal:
   \f$P(p^{m^2}|i) = \frac{E_1\{ \alpha(p,p^{m^2}|i,\pi)q(p,p^{m^2}|i,\pi) \}}
	                  {E_2\{ \alpha(p^{m^2},p|i,\pi) \}}\f$.
   Chib's paper said that the numerator is approximated with
   \f$E_1 \approx M^{-1}\sum_{g=1}^M \alpha(p^{(g)},p^{m^3}|i,\pi^{(g)})
	                             q(p^{(g)},p^{m^2}|i,\pi^{(g)})\f$

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block1_e1
   (double *log_val, double s);

/*!\brief Approximate denominator of posterior density's first block values 
          of model 4

   We are going to approximate the denominator (\f$E_2\f$) of full marginal:
   \f$P(p^{m^2}|i) = \frac{E_1\{ \alpha(p,p^{m^2}|i,\pi)q(p,p^{m^2}|i,\pi) \}}
	                  {E_2\{ \alpha(p^{m^2},p|i,\pi) \}}\f$.
   Chib's paper said that the denominator is approximated with
   \f$E_2 \approx J^{-1}\sum_{j=1}^J \alpha(p^{m^2}, p^{(j)}|i,\pi^{(j)})\f$

   \param log_val the return value of this function in logarithmic scal
   \param s a parameter value of solven accessibility
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block1_e2
   (double *log_val, double s);

/*!\brief Calcualte probability move function of numerator of posterior 
          density's first block of model 4

   This function calculates the probability move function:
   \f$\alpha(p^{(g)},p^{m^2}|i,\pi^{(g)})\f$

   \param alpha the return value of this function
   \param s a parameter value of solven accessibility
   \param sg a value of MCMC sample
   \param pg a value of MCMC sample
   \param ag a value of MCMC sample
   \param cg a value of MCMC sample
   \param gg a value of MCMC sample
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block1_alpha
   (double *alpha, double sg, double s, double pg,
                             double ag, double cg, double gg);

/*!\brief Calcualte jumping density ordinate of numerator of posterior 
          density's first block of model 4

   This function calculates the probability move function:
   \f$q(s^{(g)},s^{m^3}|i,\pi^{(g)})\f$

   \param q the return value of this function
   \param s a parameter value of solven accessibility
   \param sg a value of MCMC sample
   \param pg a value of MCMC sample
   \param ag a value of MCMC sample
   \param cg a value of MCMC sample
   \param gg a value of MCMC sample
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block1_q
   (double *q, double sg, double s, double pg,
                         double ag, double cg, double gg);

/*!\brief Approximate numerator of posterior density's second block values 
          of model 4

   We are going to approximate the numerator (\f$E_1\f$) of full conditional:
   \f$P(\pi^{m^2}|i,p^{m^2}) 
	= \frac{E_1\{ \alpha(\pi,\pi^{m^2}|i,p^{m^2})q(\pi,\pi^{m^2}|i,p^{m^2}) \}}
	       {E_2\{ \alpha(\pi^{m^2},\pi|i,p^{m^2}) \}}\f$
   Chib's paper said that the numerator is approximated with
   \f$E_1 \approx J^{-1}\sum_{j=1}^J \alpha(\pi^{(j)},\pi^{m^2}|i,p^{m^2})
	                               q(\pi^{(j)},\pi^{m^2}|i,p^{m^2})\f$

   \param log_val the return value of this function
   \param s a parameter value of solven accessibility
   \param p a parameter value of pairwise interaction
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block2_e1
   (double *log_val, double s, double p);

/*!\brief Approximate denominator of posterior density's second block values 
          of model 4

   We are going to approximate the numerator (\f$E_1\f$) of full conditional:
   \f$P(\pi^{m^2}|i,p^{m^2}) 
	= \frac{E_1\{ \alpha(\pi,\pi^{m^2}|i,p^{m^2})q(\pi,\pi^{m^2}|i,p^{m^2}) \}}
	       {E_2\{ \alpha(\pi^{m^2},\pi|i,p^{m^2}) \}}\f$
   Chib's paper said that the numerator is approximated with
   \f$E_2 \approx K^{-1}\sum_{k=1}^K \alpha(\pi^{m^2}, \pi^{(k)}|i,p^{m^2})\f$

   \param log_val the return value of this function
   \param s a parameter value of solven accessibility
   \param p a parameter value of pairwise interaction
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block2_e2
   (double *log_val, double s, double p);

/*!\brief Calcualte probability move function of numerator of posterior 
          density's second block of model 4

   This function calculates the probability move function:
   \f$\alpha(\pi^{(j)},\pi^{m^2}|i,p^{m^2})\f$

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block2_alpha
   (double *alpha, double pj, double p, double s, 
                             double aj, double cj, double gj);

/*!\brief Calcualte jumping density ordinate of numerator of posterior 
          density's second block of model 4

   This function calculates the probability move function:
   \f$q(\pi^{(j)},\pi^{m^2}|i,p^{m^2})\f$

   \param log_val the return value of this function in logarithmic scal
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \param aj a parameter value of adenine frequency 
   \param cj a parameter value of cytosine frequency 
   \param gj a parameter value of guanine frequency 
   \param w which nuc was used to proposed a new value
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block2_q
   (double *q, double pj, double p, double s, 
                         double aj, double cj, double gj);

/*!\brief Approximate numerator of posterior density's second block values 
          of model 4

   We are going to approximate the numerator (\f$E_1\f$) of full conditional:
   \f$P(\pi^{m^2}|i,p^{m^2}) 
	= \frac{E_1\{ \alpha(\pi,\pi^{m^2}|i,p^{m^2})q(\pi,\pi^{m^2}|i,p^{m^2}) \}}
	       {E_2\{ \alpha(\pi^{m^2},\pi|i,p^{m^2}) \}}\f$
   Chib's paper said that the numerator is approximated with
   \f$E_1 \approx J^{-1}\sum_{j=1}^J \alpha(\pi^{(j)},\pi^{m^2}|i,p^{m^2})
	                               q(\pi^{(j)},\pi^{m^2}|i,p^{m^2})\f$

   \param log_val the return value of this function
   \param s a parameter value of solven accessibility
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block3_e1
   (double *log_val, 
              double s, double p, double a, double c, double g);

/*!\brief Approximate denominator of posterior density's third block values 
          of model 4

   We are going to approximate the numerator (\f$E_1\f$) of full conditional:
   \f$P(\pi^{m^2}|i,p^{m^2}) 
	= \frac{E_1\{ \alpha(\pi,\pi^{m^2}|i,p^{m^2})q(\pi,\pi^{m^2}|i,p^{m^2}) \}}
	       {E_2\{ \alpha(\pi^{m^2},\pi|i,p^{m^2}) \}}\f$
   Chib's paper said that the numerator is approximated with
   \f$E_2 \approx K^{-1}\sum_{k=1}^K \alpha(\pi^{m^2}, \pi^{(k)}|i,p^{m^2})\f$

   \param log_val the return value of this function
   \param s a parameter value of solven accessibility
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block3_e2
   (double *log_val, 
              double s, double p, double a, double c, double g);

/*!\brief Calcualte probability move function of numerator of posterior 
          density's second block of model 4

   This function calculates the probability move function:
   \f$\alpha(\pi^{(j)},\pi^{m^2}|i,p^{m^2})\f$

   \param log_val the return value of this function in logarithmic scal
   \param p a parameter value of pairwise interaction
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block3_alpha
   (double *alpha, 
              double ak, double ck, double gk, double a, double c, double g,
              double s, double p, int w);

/*!\brief Calcualte jumping density ordinate of numerator of posterior 
          density's second block of model 4

   This function calculates the probability move function:
   \f$q(\pi^{(j)},\pi^{m^2}|i,p^{m^2})\f$

   \param log_val the return value of this function in logarithmic scal
   \param ak a parameter value of adenine frequency 
   \param ck a parameter value of cytosine frequency 
   \param gk a parameter value of guanine frequency 
   \param a a parameter value of adenine frequency 
   \param c a parameter value of cytosine frequency 
   \param g a parameter value of guanine frequency 
   \param s a parameter value of solven accessibility
   \param p a parameter value of pairwise interaction
   \param w which nuc was used to proposed a new value
   \return 0 for success, error code otherwise 
 */
extern int 
gmel_bf_posterior_m4_block3_q
   (double *alpha, 
              double ak, double ck, double gk, double a, double c, double g,
              double s, double p, int w);

/*@}*/


END_C_DECLS

#endif /* !PSI_BF_H */
/** @end 1 **/


