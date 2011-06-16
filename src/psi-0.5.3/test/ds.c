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

/*
  double m1[2][2] = { { 0.0, 1.0 }, { 2.0, 3.0 } };
  double m2[2][2] = { { 10.0, 11.0 }, { 12.0, 13.0 } };
  double m3[2][2] = { { 20.0, 21.0 }, { 22.0, 23.0 } };

  Psi p;
  p.inter = NULL;
  p.pro = XMALLOC (int, 3);

  inter_set (&p, 0, 1, m1); 
  inter_set (&p, 1, 2, m2); 
  inter_set (&p, 1, 3, m3); 
  inter_set (&p, 1, 3, m2); 
  inter_set (&p, 2, 3, m2); 

  inter_print (&p);

  fprintf (stderr, "----- SEARCH -----\n");

  Interaction_data *m = inter_find (&p, 0, 1);
  matrix_print (m);
  m = inter_find (&p, 2, 1);
  matrix_print (m);
  delete_psi (&p);
*/
  /* inter_remove (&p); */

/*
  Gibbs_data *p = NULL;
  Gibbs_data *found_p = NULL;
  GibbsPart ******gp = NULL;
  gibbs_set (&p, 0, gp);
  gibbs_set (&p, 1, gp);
  gibbs_set (&p, 2, gp);
  gibbs_set (&p, 3, gp);
  gibbs_set (&p, 4, gp);
  gibbs_set (&p, 0, gp);
  gibbs_print (p);
  found_p = gibbs_find (p, 2);
  fprintf (stderr, "found p: %d\n", found_p->chain);
  gibbs_remove (&p);
*/

  exit (result);
}
/** @end 1 */
