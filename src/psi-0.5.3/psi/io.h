/* io.h -- I/O
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
   \brief all I/O functions

 */

/** @start 1 **/
#ifndef PSI_IO_H
#define PSI_IO_H 1

#include <psi/common.h>

BEGIN_C_DECLS

#define PSI_IO_GMEL_LINE_MAX 200

extern const char *program_name;
extern int psi_shuffle_file (const char *old_name, const char *new_name,
                             int sample_size);
extern void write_parameter_header (FILE *fp);
extern void write_parameter (FILE *fp, parameter p, int gen);
extern void write_parameter_alldraws (FILE *fp, parameter p, int gen);
extern void debug_write_parameter (FILE *fp, parameter p);
extern void debug_write_gp (FILE *fp, gridpoint gp);
extern int fgetline (char *line[], int max, FILE *f);
extern void locate_m1_block1 ( FILE *fp, int sample_size );
extern void locate_m4_block1_e1 ( FILE *fp, int sample_size );
extern void locate_m4_block1_e2 ( FILE *fp, int sample_size );
extern void locate_m4_block2_e1 ( FILE *fp, int sample_size );
extern void locate_m4_block2_e2 ( FILE *fp, int sample_size );
extern void locate_m4_block3_e1 ( FILE *fp, int sample_size );
extern void pass_lines ( int n, FILE *fp );
extern int read_param_line (FILE *bf, int *gen,
                               double *a, double *c, double *g, double *t, 
                               double *s, double *p);
extern int read_parameter (FILE *fp, parameter *p_param);

extern void verbose_likelihood_gibbs (FILE *fp, double likelihood, 
                                          parameter top, 
                                          parameter bot);
extern void verbose_gibbs (FILE *fp, double ratio, 
                               int n_a, int n_c, int n_g, int n_t, 
                               double LnDiff_A, double LnDiff_C, 
                               double LnDiff_G, double LnDiff_T, 
                               double solv, double pair,
                               double gibbs_top, double gibbs_bot,
                               double top_max_energy, double bot_max_energy,
                               parameter top, parameter bot);

extern int gmel_read_draws (char *chain_name, 
                                double *draws, 
                                int w, int sample_size);
extern void psi_mcmc_write_posterior_sample (char *pn, const double ***post, 
                                             int n, int s);
extern int psi_io_check_est_p (char *pname, int s);
extern int psi_io_check_bf_p (char *pname, int s);
extern int psi_io_check_bf_b (char *bname);
extern int psi_io_check_bf_state (char *pname, int sample_size);
extern int psi_bf_stage_remained_sample_size (char *pname, int sample_size);
extern parameter psi_bf_stage_theta_star (char *pname, int sample_size);
extern parameter psi_bf_stage_theta_last (char *pname, int sample_size);
extern int psi_io_index_last_theta (char *pname, int sample_size);
extern int psi_io_write_number_burn (FILE *fn, int i);
extern int psi_io_read_gen_last_burn (char *pname, int *gen);
extern int psi_io_read_theta_last_burn (char *pname, parameter *theta);
extern int psi_io_remove_number_burn (char *pname);
extern int psi_io_write_number_sample (FILE *fn, int i);
extern int psi_io_read_gen_last_sample (char *pname, int *gen);
extern int psi_io_read_theta_last_sample (char *pname, parameter *theta);
extern int psi_io_remove_number_sample (char *pname);
extern int psi_io_mean_theta (char *pname, parameter *theta, int sample_size, int w);
extern void psi_io_write_mcmc_acceptance (FILE *fn,
                              double delta_s_accepted,
                              double delta_p_accepted,
                              double delta_s_rejected,
                              double delta_p_rejected,
                              int gmel_news_accept_s,
                              int gmel_news_total_s,
                              int gmel_news_accept_p,
                              int gmel_news_total_p,
                              int gmel_news_accept_n,
                              int gmel_news_total_n);
extern void psi_io_read_mcmc_acceptance (char *pname,
                              double *delta_s_accepted,
                              double *delta_p_accepted,
                              double *delta_s_rejected,
                              double *delta_p_rejected,
                              int *gmel_news_accept_s,
                              int *gmel_news_total_s,
                              int *gmel_news_accept_p,
                              int *gmel_news_total_p,
                              int *gmel_news_accept_n,
                              int *gmel_news_total_n);

END_C_DECLS

#endif /* !PSI_IO_H */
/** @end 1 **/

