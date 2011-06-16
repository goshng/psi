/* mc.h -- Bayes factor estimation
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
#ifndef PSI_MC_H
#define PSI_MC_H 1

#include <psi/common.h>

BEGIN_C_DECLS

/*!\brief Execute the mulitiple chains

   It uses more than one chain to check the convergence. It also
   take advantage of autocorrelation technique to have an "independent" sample
   this way.

    --- sample frequency check ---
   1. run a single chain until it has a sample of size \f$N\f$,
   2. calculate the autocorrelation time (\f$\tau\f$) of the sample,
   3. set \f$N_e\f$ to be \f$N/\tau\f$
   4. take \f$N_e\f$ out of the sample
   5. set \f$N_d\f$ to be \f$N-N_e\f$
   6. if \f$N_d == 0\f$ then stop the algorithm,
      else set \f$N\f$ to be \f$N_d\f$ and go to the step 1

   --- Different Burnin and Sample Frequency are determined ---

   --- convergence check --- *** NOTE: different Gibbs sampler each chain ***
   Once we find the appropriate sample frequency, then we run more
   chains to check the convergence. If the psr, the estimate of 
   convergence statistic, is not close to 1.0001, then we rerun all 
   of the chains. We repeat this convergence check procedure until
   the psr estimate is close enought to 1.0001.

   --- output the sample --- 

   --- run three more repetitions above for Bayes factor ---

   \return 0 for success, error code otherwise 
 */
extern int execute_multiple_chains ();

extern int unsetup_mc ();
extern int setup_mc_chains (int n);
extern int setup_mc_add_ip (double a, double c, double g, double s, double p);
extern int setup_mc_option (int i);

END_C_DECLS

#endif /* !PSI_MC_H */
/** @end 1 **/


