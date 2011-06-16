/* energy.h -- energy calculation
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
   \brief Energy Calculation Module
 */

/** @start 1 **/
#ifndef PSI_ENERGY_H
#define PSI_ENERGY_H 1

#include <gsl/gsl_math.h>
#include <psi/common.h>
#include <psi/ds.h>

BEGIN_C_DECLS

extern double energy_put_out (double d);

enum { 
  PSI_ENERGY_JONES,
  PSI_ENERGY_VIENNA,
  PSI_ENERGY_VLMM, 
  PSI_ENERGY_IEDB,
  PSI_ENERGY_BASTOLLA
};

#define DECLARE_GLOBAL_DREVOL                      \
  char *datname= GMEL_STR_DAT; /* --dat */      \
  FILE *datfile;                                   \
  char *energyname= GMEL_STR_ENERGY; /* --energy */   \
  FILE *energyfile;                                \
  Interaction *FirstOrder = NULL;                  \
  dat_t *Dat = NULL;                               \
  long g_seed[6];

#define DECLARE_LOCAL_DREVOL                       \
  1;

/*!\brief DrEvol's DAT structure

   DrEvol's DAT file has information of DNA sequence and protein structure. 
 */ 
typedef struct data_t {
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
   /* AAsiteInfo *AAinfo; !< energy information of a given DrEvol's DAT file */
} Data;

/*!\brief statiatical potential information of a given amino acid site

   Solvent accessiblity and pairwise interaction energy are assigned
   to each amino acid site.
 */ 
typedef struct AAsiteInfo_t {
   double *solvent;    /*!< solvent accessibility */
   double **pairwise;  /*!< pairwise interaction */
   int numNeighbors;   /*!< number of neighboring amino acid */
   int *DDDneighbors;  /*!< neighboring site within a distance */
} AAsiteInfo;
        
/*!\brief huge statistical potential strucure

   DrEvol's energy file has David Jone's statistical potential information:
   solvent accessibility and pairwise interaction energy from huge prtein
   3-D tertiay structure. Now, this structure also has energy information
   of a given DrEvol's DAT. I'm not sure if this, AAsiteInfo, member should
   be taken out of it.
 */ 
typedef struct interaction_t {
   AAsiteInfo *AAinfo; /*!< energy information of a given DrEvol's DAT file */
   int ***energyAcc;   /*!< stores the pairwise distance between all CB atoms */
   int **nAccess;      /*!< int list that obtains correct index for neighbor matrix */
   int **neighList;    /*!< stores the AA neighbors, correct Nuc, and K param */
   double ****energy;  /*!< stores all matrices for pairpotential NRG */
} Interaction;

#ifndef LEN_FILENAME
#  define LEN_FILENAME 256
#endif

typedef struct psi_energy_jones {
  char *dat_filename;
  FILE *datfile;
  char *energy_filename;
  FILE *energyfile;
  Data data;
  Interaction info;
  Psi psi;
} PsiEnergyJones;

/*!\brief Accessors of jones_measure
 */
extern int* dna_jones_measure ();
extern int* pro_jones_measure ();
extern int len_dna_jones_measure ();
extern int len_pro_jones_measure ();
extern int init_codon_jones_measure ();

extern void data_members_jones_measure (int *n_a, int *n_c, int *n_g, int *n_t,
                                        double *solv, double *pair);

/*!\brief Prepare data files for the energy
 */
extern int
setup_energy_jones (const char *data_fn, const char *energy_fn);

/*!\brief Initialize dat and energy files
   
   It open the dat and energy files. It also read the dat file
   and energy file.

   \param e use GENE_THREADER
   \return 0 for success, error code otherwise 
 */ 
extern int initialize_drevol (int e);

/*!\brief Finalize dat and energy files
   
   It open the dat and energy files.

   \param e use GENE_THREADER
   \return 0 for success, error code otherwise 
 */
extern int finalize_drevol (int e);

/*!\brief Calculate energy
  
   It uses Interaction information to calculate a solvent accessibility
   score and pairwise interaction score of a given protein sequence. 
 */
extern int score_drevol (int e, ...);

/*!\brief Calculate energy
  
   It uses Interaction information to calculate a solvent accessibility
   score and pairwise interaction score of a given protein sequence. 
 */
extern void calc_energy (int *protein, 
                         double *solvNRG, double *pairNRG);

/*!\brief Calculate energy of a neighboring sequence
  
   It is a cousin of calc_energy function which is still going to calculate
   solvent accessibility score and pairwise interaction score of a given
   protein sequence. It uses neighboring sequence energy score to 
   calculate score of a given sequence. This won't dramatically decrease
   computing time, but it's computing time effect grows as more and more
   energy scores are calculated. So, hopefully, it will reduce a huge 
   amount of time so that we can have ice-cream sitting on table outside
   rather than on computer desk like me.
 */
extern void calc_energy_neigh (int *protein, 
                               double *solv, double *pair,
                               int site, int aa);

/*!\brief Calculate jones' measure energy score of a protein sequence
 */
extern int energy_is (int *protein, double *solv, double *pair );
/*!\brief Calculate solvent accessibility energy score of a protein sequence
 */
extern int energy_solv_is (int *protein, double *solv);

/*!\brief Calculate energy score of a site
 */
extern int energy_res_solv_is (int pos, int aa, double *solv);
extern int energy_res_pair_is (int *protein, int pos, int aa, double *solv);
extern int energy_res_is (int *protein, int pos, int aa, double *solv, double *pair);

/*!\brief Calculate pairwise interaction energy score of a protein sequence
 */
extern int energy_pair_is (int *protein, double *pair);

extern int energy_near_is (int *protein, int site, int aa,
                              double *solv, double *pair);
extern int energy_near_solv_is (int *protein, int site, int aa,
                                   double *solv);
extern int energy_near_pair_is (int *protein, int site, int aa,
                                   double *pair);

/* from create */


/* from delete */
extern void delete_psi (Psi *psi);
extern void delete_dat (Data **dat);
extern void delete_energy (Interaction **energy, int AAlen);

extern void write_dat (FILE *ofile);
extern void write_energy (FILE *ofile);
extern int choose_a_nuc_stationary (int e, int *dna, const int site, double parameter[]);
extern int choose_a_nuc_stationary_fast (int e, int *dna, int *pro, const int site, double log_theta_star[]);
extern int choose_a_nuc_rate (int e, int *site, int *nucleotide, int *dna, double parameter[]);
extern int psi_energy_dna_report (int e, int *dna, double parameter[]);
extern int psi_energy_dna_state (int e, int *dna, double parameter[], double **states, int n);
extern int psi_energy_valid_state (int e, int *dna, double parameter[], double **state_value, int n);

END_C_DECLS

#endif /* !PSI_ENERGY_H */
/** @end 1 **/
