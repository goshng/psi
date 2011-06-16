/* ds.h -- data structures
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
/** @start 1 **/
#ifndef PSI_DS_H
#define PSI_DS_H 1

#include <psi/common.h>

#ifndef NUM_ELEMENT
#  define NUM_ELEMENT 20
#endif

BEGIN_C_DECLS

typedef struct gibbs_data {
  struct gibbs_data *next;
  int chain;
  GibbsPart ******info;
} Gibbs_data;

typedef struct interaction_data {
  struct interaction_data *next;
  int res1;
  int res2;
  double matrix[NUM_ELEMENT][NUM_ELEMENT];
} Interaction_data;

typedef struct psi {
   int init_codon;     /*!< init codon is starting one */
   int transl_table_id;/*!< translation table id */
   int n_seq;          /*!< number of sequences */
   int len_dna;        /*!< length of DNA sequence */
   char *str_dna;      /*!< DNA sequence itself */
   int *dna;           /*!< DNA sequence itself in number representation */
   int n_a;            /*!< number of nucleotide A */
   int n_c;            /*!< number of nucleotide C */
   int n_g;            /*!< number of nucleotide G */
   int n_t;            /*!< number of nucleotide T */
   int len_pro;        /*!< length of protein sequence */
   char *str_pro;      /*!< protein sequence itself */
   int *pro;           /*!< protein sequence itself in number representation */
   double solv;        /*!< solvent accessibility energy of the sequence */
   double pair;        /*!< pairwise interaction energy of the sequence */
   Interaction_data *inter;
   int *access;
   int *res_num;
   double solvent_data[5][20]; /* five categories and twenty amino acids */
} Psi;

extern int inter_set (Psi *psi, int res1, int res2, 
                      double matrix[][NUM_ELEMENT]);
extern Interaction_data* inter_find (Psi *psi, const int res1, const int res2);
extern void inter_print (Psi *psi);
extern void inter_remove (Interaction_data **inter);
/* extern void inter_remove (Psi **psi); */
extern void matrix_print (Interaction_data *m);
extern double inter_score (Psi *psi, int *protein);
extern double inter_score_res (Psi *psi, int *protein, int site, int aa);

extern int gibbs_set (Gibbs_data **g, int chain, GibbsPart ******info);
extern Gibbs_data* gibbs_find (Gibbs_data *head, int chain);
extern void gibbs_print (Gibbs_data *head);
extern void gibbs_remove (Gibbs_data **head);

END_C_DECLS

#endif /* !PSI_DS_H */
/** @end 1 **/
