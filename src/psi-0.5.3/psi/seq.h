/* seq.h -- biological sequence
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

/*!\file seq.h
   \author Sang Chul Choi
   \brief Biological sequence
 */

/* We got the original codon table from 
   ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
*/

/** @start 1 **/
#ifndef PSI_SEQ_H
#define PSI_SEQ_H 1

#include <psi/common.h>

BEGIN_C_DECLS

#ifndef LEN_BUFFER
#  define LEN_BUFFER 80
#endif

#ifndef NUM_CODON
#  define NUM_CODON 64
#endif

enum {
  PSI_DNA_A, 
  PSI_DNA_C, 
  PSI_DNA_G, 
  PSI_DNA_T 
};

enum { 
  PSI_CODON_FIRST, 
  PSI_CODON_SECOND, 
  PSI_CODON_THIRD
};

enum {
  PSI_AA_ALA, /* A */
  PSI_AA_ARG, /* R */
  PSI_AA_ASN, /* N */
  PSI_AA_ASP, /* D */
  PSI_AA_CYS, /* C */
  PSI_AA_GLN, /* Q */
  PSI_AA_GLU, /* E */
  PSI_AA_GLY, /* G */
  PSI_AA_HIS, /* H */
  PSI_AA_ILE, /* I */
  PSI_AA_LEU, /* L */
  PSI_AA_LYS, /* K */
  PSI_AA_MET, /* M */
  PSI_AA_PHE, /* F */
  PSI_AA_PRO, /* P */
  PSI_AA_SER, /* S */
  PSI_AA_THR, /* T */
  PSI_AA_TRP, /* W */
  PSI_AA_TYR, /* Y */
  PSI_AA_VAL, /* V */
  PSI_AA_STP  /* STOP */
};

enum {
  PSI_AA_A, /* A */
  PSI_AA_R, /* R */
  PSI_AA_N, /* N */
  PSI_AA_D, /* D */
  PSI_AA_C, /* C */
  PSI_AA_Q, /* Q */
  PSI_AA_E, /* E */
  PSI_AA_G, /* G */
  PSI_AA_H, /* H */
  PSI_AA_I, /* I */
  PSI_AA_L, /* L */
  PSI_AA_K, /* K */
  PSI_AA_M, /* M */
  PSI_AA_F, /* F */
  PSI_AA_P, /* P */
  PSI_AA_S, /* S */
  PSI_AA_T, /* T */
  PSI_AA_W, /* W */
  PSI_AA_Y, /* Y */
  PSI_AA_V, /* V */
  PSI_AA_X  /* STOP */
};

/* the numbers are from the "gc.prt" file */
enum { 
  PSI_CODON_TABLE_1, 
  PSI_CODON_TABLE_2, 
  PSI_CODON_TABLE_3, 
  PSI_CODON_TABLE_4, 
  PSI_CODON_TABLE_5, 
  PSI_CODON_TABLE_6, 
  PSI_CODON_TABLE_9, 
  PSI_CODON_TABLE_10, 
  PSI_CODON_TABLE_11, 
  PSI_CODON_TABLE_12, 
  PSI_CODON_TABLE_13, 
  PSI_CODON_TABLE_14, 
  PSI_CODON_TABLE_15, 
  PSI_CODON_TABLE_16,
  PSI_CODON_TABLE_21, 
  PSI_CODON_TABLE_22, 
  PSI_CODON_TABLE_23
};

/*!\brief Number of a nucleotide in a DNA sequence
 */
extern int num_nucleotide (int *dna, int len, int w);

/*!\brief Search a DNA sequence for a codon
 */
extern int find_codon_in_dna (int *codon, int *dna, int pos);

/*
extern const char *program_name;
extern void set_program_name (const char *argv0);

extern void psi_warning      (const char *message);
extern void psi_error        (const char *message);
extern void psi_fatal        (const char *message);
*/

/*!\brief Choose one of 17 codon tables
 */
extern int choose_transl_table (int id, int init);


/*!\brief Read the translation table file, "gc.prt"

    The table file was downloaded from the NCBI ftp server,
    ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
 */
extern int read_translaton_table (const char *fn);

/*!\brief Translator from DNA to protein

extern void dna_to_pro (int *dna, int *pro, int nuclen); not implemented
   From DNA, it spits out protein sequence using standard codon translation
   table.
 */

/*!\brief Convert char to int type of DNA
 */
extern int dna_char2int (char *str_dna, int *dna, int len);

/*!\brief Convert int to char type of DNA
 */
extern int dna_int2char (int *dna, char *str_dna, int len);

/*!\brief Convert char to int type of protein
 */
extern int pro_char2int (char *str_pro, int *pro, int len);

/*!\brief Convert int to char type of protein
 */
extern int pro_int2char (int *pro, char *str_pro, int len);

/*!\brief 
   
   \param dna [in] DNA squence in integer
   \param protein [out] Protein sequence in integer
   \param len [in] Length of protein sequnce
 */
extern int dna2protein (int *dna, int *protein, int len);

/*!\brief Calculate codon frequency using nucleotide one

   \param pi [in] nucleotide frequency already allocated
   \param c [out] amino acid frequency already allocated
 */
extern int codon_frequency (double *pi, double *c);

/*!\brief Return amino acid of a codon

   \param codon [in] codon
   \return amino acid 
 */
extern int aa_codon_is (int codon[]);

/*!\brief Change translation table

   \param i [in] codon table number
   \return SUCCESS
 */
extern int change_translation_table (int i);

extern int psi_seq_diff_seqs (int *seq_a, int *seq_b, int len);

extern double seq_put_out (double d);

END_C_DECLS

#endif /* !PSI_SEQ_H */
/** @end 1 **/
