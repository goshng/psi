/* error.c -- display formatted error diagnostics of varying severity
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

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "common.h"
#include "error.h"

static void error (int exit_status, const char *mode,
		   const char *message, va_list ap);

const char *program_name = NULL;

void
set_program_name (const char *path)
{
  assert (!program_name);
  program_name = xstrdup ((char *) basename (path));
}

static void
error (int exit_status, const char *mode, const char *message, va_list ap)
{
  fprintf (stderr, "%s: %s: ", program_name, mode);
  vfprintf (stderr, message, ap);
  fprintf (stderr, ".\n");

  if (exit_status >= 0)
    exit (exit_status);
}

void
psi_warning (const char *message, ...)
{
  va_list ap;
  va_start (ap, message);
  error (-1, "warning", message, ap);
  va_end (ap);
}

void
psi_error (const char *message, ...)
{
  va_list ap;
  va_start (ap, message);
  error (-1, "ERROR", message, ap);
  va_end (ap);
}

void
psi_fatal (const char *message, ...)
{
  va_list ap;
  va_start (ap, message);
  error (EXIT_FAILURE, "FATAL", message, ap);
  va_end (ap);
}

int
psi_err_base (double a, double c, double g, double t)
{
  if (a >= 1.0 || c >= 1.0 || g >= 1.0 || t >= 1.0 ||
      a <= 0.0 || c <= 0.0 || g <= 0.0 || t <= 0.0 ||
      a + c + g + t > 1.00001)
    {
      psi_fatal ("invalid base frequency, %lf, %lf, %lf, %lf",
                 a, c, g, t);
    }
  return EXIT_SUCCESS;
}

int 
psi_err_finite (double v)
{
  if (!finite(v)) 
    {
      psi_fatal ("not finite");
    }
  return EXIT_SUCCESS;
}

int
gmel_err_math (const char *prog_name, 
               const char *file_name,
               const int line_number,
               const char *func_name,
               double v, double d, int e) 
{
  int r = EXIT_SUCCESS;
  if (e == ERANGE) {
     if (v == HUGE_VAL) {
        fprintf (stderr, "%s:%s:%d:%s; %s %lf = exp(%lf)\n",
                 prog_name, file_name, line_number, func_name,
                 "Huge value", v, d);
        psi_fatal ("Huge value");
     } else if (v == -HUGE_VAL) {
        fprintf (stderr, "%s:%s:%d:%s; %s %lf = exp(%lf)\n",
                 prog_name, file_name, line_number, func_name,
                 "Negative huge value", v, d);
        psi_fatal ("Negative huge value");
     } else if (v == 0) {
        fprintf (stderr, "%s:%s:%d:%s; %s %lf = exp(%lf)\n",
                 prog_name, file_name, line_number, func_name,
                 "Too small value", v, d);
        psi_fatal ("Too small value");
     } else {
        fprintf (stderr, "%s:%s:%d:%s; %s %lf = exp(%lf)\n",
                 prog_name, file_name, line_number, func_name,
                 "Something else", v, d);
        psi_fatal ("Something else");
        assert (0);
     }
     /* r = ERR_MATH; */
     r = EXIT_FAILURE;
  } else if (e == EDOM) {
     fprintf (stderr, "%s:%s:%d:%s; %s %lf = exp(%lf)\n",
              prog_name, file_name, line_number, func_name,
              "Domain Error", v, d);
     psi_fatal ("Domain Error");
  } 

  return r;
}

int
gmel_err_same (const char *prog_name, 
               const char *file_name,
               const int line_number,
               const char *func_name,
               double one, double prev_one) 
{
  int r = EXIT_SUCCESS;
  if (one == prev_one) {
    fprintf (stderr, "%s:%s:%d:%s; %s, %lf == %lf\n",
             prog_name, file_name, line_number, func_name,
             "Two numbers are the same", one, prev_one);
    psi_fatal ("Two numbers are the same");
  }
 
  return r;
}

int 
gmel_err_base (const char *prog_name, 
               const char *file_name,
               const int line_number,
               const char *func_name,
               double a, double c, double g, double t) 
{
  if (a >= 1.0 || c >= 1.0 || g >= 1.0 || t >= 1.0 ||
      a <= 0.0 || c <= 0.0 || g <= 0.0 || t <= 0.0 ||
      a + c + g + t > 1.0000001) {
     fprintf (stderr, "%s:%s:%d:%s; %s, a=%lf, c=%lf, g=%lf, t=%lf\n",
              prog_name, file_name, line_number, func_name,
              "Invalid base frequency", a, c, g, t);
     psi_fatal ("Invalid base frequency");
  }
  return EXIT_SUCCESS;
}

int
gmel_err_block (const char *prog_name, 
                const char *file_name,
                const int line_number,
                const char *func_name,
                const char *content,
                int r)
{ 
  if (r != EXIT_SUCCESS) {
     fprintf (stderr, "%s:%s:%d:%s; %s error code - %d\n",
              prog_name, file_name, line_number, func_name,
              content, r);
     psi_fatal ("ERR_BLOCK");
  }
  return EXIT_SUCCESS;
}

int
gmel_err_finite (const char *prog_name, 
                 const char *file_name,
                 const int line_number,
                 const char *func_name,
                 double value) 
{
  int r = EXIT_SUCCESS;
  if (finite (value) == 0) {
     fprintf (stderr, "%s:%s:%d:%s; %s, %lf\n",
              prog_name, file_name, line_number, func_name,
              "value is not finite", value);
     psi_fatal ("value is not finite");
  }
  return r;
}






