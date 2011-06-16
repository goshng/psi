/* grid.c -- test program of the grid module
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

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <signal.h>

#include <psi/psi.h>

#ifndef LEN_DNA
#  define LEN_DNA 192
#endif

#ifndef LEN_PROTEIN
#  define LEN_PROTEIN 64
#endif

/* static int psi_init  (Sic *psi); */

/** @start 1 */
int
/* main (int argc, char * const argv[]) */
main (void)
{
  int result = EXIT_SUCCESS;

  gmel_grid_new (-5, 5, 19,
                 -1, 1, 19,
                  0, 1, 10,
                  0, 1, 10,
                  0, 1, 10, 100);

  /* gmel_grid_print_grid (); */

  gmel_grid_del ();

  exit (result);
}
/** @end 1 */
