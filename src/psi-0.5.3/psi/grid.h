/* grid.h -- grid for Gibbs sampler
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
   \brief GMEL Grid Point Module

 */

/** @start 1 **/
#ifndef PSI_GRID_H
#define PSI_GRID_H 1

#include <psi/common.h>

/* #include "structs.h" */

BEGIN_C_DECLS

extern double grid_put_out (double d);

enum { 
  PSI_GRID_S, 
  PSI_GRID_P, 
  PSI_GRID_A, 
  PSI_GRID_C, 
  PSI_GRID_G
};

/*!\brief Set 

   This function is used in news.c file.
   \param index the return value of gridpoint index
   \param value value whose nearest gridpoint index is sought
   \param grid_begin lower boundary value of the grid-points
   \param grid_end upper boundary value of the grid-points
   \param grid_number number of grid-points
   \return 0 for SUCCESS
 */



/*!\brief Locate a gridpoint for a given one value

   This function is used in news.c file.
   \param index the return value of gridpoint index
   \param value value whose nearest gridpoint index is sought
   \param grid_begin lower boundary value of the grid-points
   \param grid_end upper boundary value of the grid-points
   \param grid_number number of grid-points
   \return 0 for SUCCESS
 */
extern int gmel_grid_locate_parameter (int *index, double value,
                                       double grid_begin, double grid_end, 
                                       double grid_number);

/*!\brief Make a new grid points for only one parameter
 
   \param grid gridpoint array already memory allocated
   \param min lower bound of the parameter
   \param max upper bound of the parameter
   \param n_g number of gridpoints
   \return 0 for SUCCESS 
 */
extern int gmel_grid_new_grid (double *grid, double min, double max, int n_g);

/*!\brief Make a new grid points
 
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
   \param n_s number of desired gridpoints
   \return 0 for SUCCESS 
 */
extern int gmel_grid_new (double min_s, double max_s, int n_s,
                          double min_p, double max_p, int n_p,
                          double min_a, double max_a, int n_a,
                          double min_c, double max_c, int n_c,
                          double min_g, double max_g, int n_g, int gs);

/*!\brief Load Gibbs Sample

   \return 0 for SUCCESS 
 */
extern int gmel_grid_load (const char* gname);

/*!\brief Save Gibbs Sample

   \return 0 for SUCCESS 
 */
extern int gmel_grid_save (const char* gname);

/*!\brief Deallocate memory used for two static variables

   \return 0 for SUCCESS 
 */
extern int gmel_grid_del ();

/*!\brief Locate a  gridpoint nearby a given parameter
 
   \param gp the return value of gridpoint
   \param param the input values of parameter
   \return 0 for SUCCESS
 */
extern int gmel_grid_locate (gridpoint *gp, parameter param);

/*!\brief Return parameter values of a given gridpoint

   \param param
   \param gp
   \return 0 for SUCCESS
 */
extern int gmel_grid_param_gridpoint (parameter *param, gridpoint gp);

/*!\brief Set the first value of gridpoint at a fixed one

   \param grid_type one of five types: A_GRID, C_GRID, ...
   \param value value for the gridpoint
   \return 0 for SUCCESS
 */
extern int gmel_grid_fix (int grid_type, double value);


/*!\brief Store Gibbs information for a gridpoint

   \param gibbs the information to be stored for a given gridpoint
   \param gp gridpoint
   \return 0 for SUCCESS
 */
extern int gmel_grid_put_info (GibbsPart *gibbs, gridpoint gp);

/*!\brief Retrieve Gibbs information for a gridpoint
 
   \param gibbs the return value of Gibbs information stored at the gridpoint
   \param gp gridpoint
   \return 0 for SUCCESS
 */
extern int gmel_grid_get_info (GibbsPart **gibbs, gridpoint gp);

/*!\brief Check if the gridpoints are too sparse

   \param top
   \param bot
   \return 0 if no, 1 otherwise
 */
extern int gmel_grid_check_sparse (parameter top, parameter bot);

extern int gmel_grid_print_grid ();

extern void psi_grid_create_multi_gibbs (int n);
extern void psi_grid_change_multi_gibbs (int n);

END_C_DECLS
 
#endif /* !PSI_GRID_H */
/** @end 1 **/

