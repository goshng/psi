/* error.h -- display formatted error diagnostics of varying severity
   Copyright (C) 2000 Gary V. Vaughan
  
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef PSI_SIM_H
#define PSI_SIM_H 1

#include <psi/common.h>

BEGIN_C_DECLS

extern void psi_sim_gibbs (int iter);
extern void psi_sim_rate (int iter);

END_C_DECLS

#endif /* !PSI_SIM_H */
/** @end 1 **/
