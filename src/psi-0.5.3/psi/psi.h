/* psi.h -- create and maintain Psi ADT
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

/** @start 1 */
#ifndef PSI_PSI_H
#define PSI_PSI_H 1

#include <psi/common.h>
#include <psi/defs.h>
#include <psi/error.h>
#include <psi/seq.h>
#include <psi/grid.h>
#include <psi/energy.h>
#include <psi/rng.h>
#include <psi/gslwrap.h>
#include <psi/gibbs.h>
#include <psi/mcmc.h>
#include <psi/io.h>
#include <psi/bf.h>
#include <psi/ds.h>
#include <psi/mc.h>
#include <psi/nsv.h>
#include <psi/sim.h>

/** @end 1 */
/** @start 2 */
BEGIN_C_DECLS

enum {
  PSI_ERROR = -1,
  PSI_OKAY = 0,
  PSI_INCOMPLETE,
  PSI_BREAK,
  PSI_CONTINUE,
  PSI_EXIT
};

struct builtintab;

typedef struct statedata {
  struct statedata *next;       /* so they can be chained */
  char *key;			/* used as a key to find the right data */
  void *data;                   /* associated state data */
  void (*delete) (void *data);
} PsiState;


END_C_DECLS
#endif /* !PSI_PSI_H */

