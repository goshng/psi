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

#ifndef PSI_ERROR_H
#define PSI_ERROR_H 1

#include <psi/common.h>

BEGIN_C_DECLS

extern const char *program_name;

extern void set_program_name (const char *argv0);
extern void psi_warning	     (const char *message, ...);
extern void psi_error	     (const char *message, ...);
extern void psi_fatal	     (const char *message, ...);
extern int psi_err_base (double a, double c, double g, double t);
extern int psi_err_finite (double v);
extern int gmel_err_math (const char *prog_name, 
               const char *file_name,
               const int line_number,
               const char *func_name,
               double v, double d, int e);
extern int gmel_err_same (const char *prog_name, 
               const char *file_name,
               const int line_number,
               const char *func_name,
               double one, double prev_one);
extern int gmel_err_base (const char *prog_name, 
               const char *file_name,
               const int line_number,
               const char *func_name,
               double a, double c, double g, double t);
extern int gmel_err_block (const char *prog_name, 
                const char *file_name,
                const int line_number,
                const char *func_name,
                const char *content,
                int r);
extern int gmel_err_finite (const char *prog_name, 
                 const char *file_name,
                 const int line_number,
                 const char *func_name,
                 double value);


END_C_DECLS

#endif /* !PSI_ERROR_H */
/** @end 1 **/
