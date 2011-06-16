/* psi.c -- read commands, evaluate them and print the results
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

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <signal.h>

#include "psi.h"

/* static int psi_init  (Sic *psi); */

/** @start 1 */
int
/* main (int argc, char * const argv[]) */
main (void)
{
  int result = EXIT_SUCCESS;


  FILE *fp = NULL;

/*
  fp = fopen ("output/pdb1g90.A.p", "a");  
  psi_io_write_number_burn (fp, 3);

  fclose (fp);
  int r; 
  parameter theta;
  int last_burn_i; 
  r = psi_io_read_gen_last_burn ("output/pdb1g90.A.p", &last_burn_i);
  if (r == EXIT_SUCCESS)
    fprintf (stderr, "last burn: %d\n", last_burn_i);
  r = psi_io_read_theta_last_burn ("output/pdb1g90.A.p", &theta);
  if (r == EXIT_SUCCESS)
    fprintf (stderr, "theta: %lf %lf %lf %lf %lf %lf\n", 
             theta.a, theta.c, theta.g, theta.t, theta.s, theta.p);
*/

  psi_io_remove_number_burn ("output/pdb1g90.A.p");
/*
  Sic *psi = psi_new ();
  
  if (psi_init (psi) != PSI_OKAY)
      psi_fatal ("psi initialisation failed");
  signal (SIGINT, SIG_IGN);
  setbuf (stdout, NULL);

  psistate_set (psi, "PS1", "] ", NULL);
  psistate_set (psi, "PS2", "- ", NULL);
  
  evalstream (psi, stdin);
*/

  exit (result);
}
/** @end 1 */
