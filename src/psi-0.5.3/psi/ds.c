/* ds.c -- data structures
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

/** @start 1 */
#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "common.h"
#include "ds.h"

int
inter_set (Psi *psi, int res1, int res2, double matrix[][NUM_ELEMENT])
{
  int i, j;
  Interaction_data *inter= inter_find (psi, res1, res2); 
 
  if (inter)
    {
      assert (res1 > res2);
      return EXIT_FAILURE;
    }
  else 
    {
      assert (res1 < res2);
      inter = XMALLOC (Interaction_data, 1);
      inter->next = psi->inter;
      psi->inter = inter;
      inter->res1 = res1;
      inter->res2 = res2;
      for (i = 0; i < NUM_ELEMENT; i++)
        {
          for (j = 0; j < NUM_ELEMENT; j++)
            {
              inter->matrix[i][j] = matrix[i][j];
/* fprintf (stderr, "%lf ", matrix[i][j]); */
            }
        }
/* fprintf (stderr, "\n"); */
/*
      memcpy (inter->matrix, matrix, 
              sizeof (double) * NUM_ELEMENT * NUM_ELEMENT );
*/
    }

  return EXIT_SUCCESS;
}

Interaction_data *
inter_find (Psi *psi, const int res1, const int res2) 
{
  Interaction_data *inter;

  for (inter = psi->inter; inter; inter = inter->next) 
    {
      if (inter->res1 == res1 && inter->res2 == res2)
        break;
      if (inter->res1 == res2 && inter->res2 == res1)
        break;
    }

  return inter;
}

void 
inter_print (Psi *psi)
{
  Interaction_data *inter;
  for (inter = psi->inter; inter; inter = inter->next) 
    {
      matrix_print (inter); 
    }
}

void 
inter_remove (Interaction_data **head)
{
  Interaction_data *g = NULL;
  while ((*head)->next != NULL)
    {
      g = (*head);
      (*head) = (*head)->next;
      XFREE (g);
    }
  XFREE (*head);
}
/*
void 
inter_remove (Psi **psi)
{
  Interaction_data *g = NULL;

  Interaction_data **head = &((*psi)->inter);
  while ((*head)->next != NULL)
    {
      g = (*head);
      head = (*head)->next;
      XFREE (g);
    }
  XFREE (*head);
}
*/

/*
void
gibbs_remove (Gibbs_data **head)
{
  Gibbs_data *g = NULL;
  while ((*head)->next != NULL)
    {
      g = *head;
      *head = (*head)->next;
      XFREE (g);
    }
  XFREE (*head);
}
*/


void 
matrix_print (Interaction_data *m)
{
  int i, j;
  fprintf (stderr, "%d - %d :\n",  m->res1, m->res2);
  for (i = 0; i < NUM_ELEMENT; i++)
    {
      for (j = 0; j < NUM_ELEMENT; j++)
        {
          fprintf (stderr, "%lf ", m->matrix[i][j]);
        }
      fprintf (stderr, "\n");
    }
  fprintf (stderr, "\n");
}

double
inter_score (Psi *psi, int *protein)
{
  double v = 0.0;
  int aa1, aa2;

  Interaction_data *inter;
  for (inter = psi->inter; inter; inter = inter->next) 
    {
      aa1 = protein[inter->res1];
      aa2 = protein[inter->res2];
      v += inter->matrix[aa1][aa2];
    }

  return v;
}

double
inter_score_res (Psi *psi, int *protein, int site, int aa)
{
  double v = 0.0;
  int aa1, aa2;

  Interaction_data *inter;
  for (inter = psi->inter; inter; inter = inter->next) 
    {
      if (inter->res1 == site)
        {
          aa1 = aa;
          aa2 = protein[inter->res2];
          v += inter->matrix[aa1][aa2];
        }
      else if (inter->res2 == site) 
        {
          aa1 = protein[inter->res1];
          aa2 = aa;
          v += inter->matrix[aa1][aa2];
        }
    }

  return v;
}

int 
gibbs_set (Gibbs_data **g, int chain, GibbsPart ******info)
{
  assert (chain >= 0);
  Gibbs_data *gn = gibbs_find (*g, chain); 
 
  if (gn)
    {
      return EXIT_FAILURE;
    }
  else 
    {
      gn = XMALLOC (Gibbs_data, 1);
      gn->info = info;
      gn->chain = chain;
      if (*g == NULL)
        {
          *g = gn;
          (*g)->next = NULL;
        }
      else
        {
          gn->next = *g;
          *g = gn;
        }
    }

  return EXIT_SUCCESS;
}

Gibbs_data* 
gibbs_find (Gibbs_data *head, int chain)
{
  Gibbs_data *g;

  for (g = head; g; g = g->next) 
    {
      if (g->chain == chain)
        break;
    }

  return g;
}

void 
gibbs_print (Gibbs_data *head)
{
  Gibbs_data *g;

  for (g = head; g; g = g->next) 
    {
      fprintf (stderr, "chain: %d\n", g->chain);
    }
}

/* call memory deallocation function before calling this */
void
gibbs_remove (Gibbs_data **head)
{
  Gibbs_data *g = NULL;
  while ((*head)->next != NULL)
    {
      g = *head;
      *head = (*head)->next;
      XFREE (g);
    }
  XFREE (*head);
}



