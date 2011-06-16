/* gibbs.h -- Gibbs sampler
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
/** @start 1 **/
#ifndef PSI_GIBBS_H
#define PSI_GIBBS_H 1

#include <psi/common.h>

BEGIN_C_DECLS

/*!\file
   \author Sang Chul Choi
   \brief Gibbs Sampler

 */

extern int setup_gibbs (int s, int f, int b);
extern int unsetup_gibbs ();
extern int access_gibbs_size ();
extern int gmel_bf_likelihood_gibbs_sampler (GibbsPart *p_sampled_seqs, 
                                             parameter theta_star);
extern int gmel_bf_likelihood_sumup_gibbs (double *sum, 
                                          double *max_energy, 
                                          parameter theta, 
                                          parameter theta_star);
extern int gmel_gibbs_write_usage (FILE *fp);

END_C_DECLS

#endif /* !PSI_GIBBS_H */
/** @end 1 **/


