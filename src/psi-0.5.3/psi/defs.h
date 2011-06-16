/* defs.h -- definitions
   Copyright (C) 2004 
                 Bioinformatics Research Center,
                 North Carolina State University

   I don't have any idea about distribution issue of this program because
   I am not a original author. Doug and Jeff would be able to talk about 
   that issue. I think that I should change the copyright statement above. */

/*!/file 
   /author Sang Chul Choi
   /brief program wide definition

   These definitions are for whole program configuration. I might find other
   way to configure the program, like using config file.
 */

#ifndef PSI_DEFS_H_
#define PSI_DEFS_H_ 1

#define GMEL_RNG_INIT 1
#define GMEL_NUMBER_CHAINS 5
#define GMEL_POS_INIT_S 2.5
#define GMEL_NEG_INIT_S -2.5
#define GMEL_POS_INIT_P 0.25
#define GMEL_NEG_INIT_P -0.25

#define GMEL_GRIDPOINT_S_BEGIN    -5.0
#define GMEL_GRIDPOINT_S_END       5.0
#define GMEL_GRIDPOINT_S_NUMBER   51
#define GMEL_GRIDPOINT_P_BEGIN    -0.3
#define GMEL_GRIDPOINT_P_END       0.3
#define GMEL_GRIDPOINT_P_NUMBER   51
#define GMEL_GRIDPOINT_N_BEGIN     0.0
#define GMEL_GRIDPOINT_N_END       1.0
#define GMEL_GRIDPOINT_N_NUMBER   11 
#define GMEL_FLATPRIOR_S_BEGIN    -5.0
#define GMEL_FLATPRIOR_S_END       5.0
#define GMEL_FLATPRIOR_P_BEGIN    -0.3
#define GMEL_FLATPRIOR_P_END       0.3
#define GMEL_FLATPRIOR_N_BEGIN     0.0
#define GMEL_FLATPRIOR_N_END       1.0

#define GMEL_INIT_S 0.0
#define GMEL_INIT_P 0.0
#define GMEL_INIT_N 0.25

#define GMEL_SAMPLE_SIZE 1000
#define GMEL_SAMPLE_FREQ 1
#define GMEL_SAMPLE_BURN 100
#define GMEL_LOG_FREQ 100000

#define DELTA_PI      0.1
#define DELTA_PI_MU   100    /* 1 / (DELTA_PI * DELTA_PI) */
#define DELTA_W_S     1.0    
#define DELTA_W_P     0.1 

#define GMEL_LINE_MAX 200    /* one line length of output file */

#define GMEL_STR_DAT               "../data/pdb1r48.A.dat"
#define GMEL_STR_ENERGY            "../data/energy.int"
#define GMEL_STR_STDOUT            "stdout"
#define GMEL_STR_STATE             "rname"
#define GMEL_STR_PARAM             "param"
#define GMEL_STR_ALLDRAWS          "alldraws"
#define GMEL_STR_GIBBS             "gibbs"
#define GMEL_STR_SETUP             "setup"
#define GMEL_STR_BFACTOR           "bfactor"
#define GMEL_STR_Q                 "q"
#define GMEL_STR_PARAM_PREFIX      "pdb1r48.A.b."
#define GMEL_STR_PARAM_SUFFIX      ".out.p"

#define GMEL_GIBBS_BURN      1
#define GMEL_GIBBS_FREQ      1
#define GMEL_GIBBS_SIZE     30

#define INTERACT 10.0       /* Distance at which AA's can biologically 
                               interact                                        */

#define ACCESS_CATEGORY 5   /* There are 5 categories of accessibility         */
#define NUM_INTPAIR     5   /* There are 5: 0:CBCB, 1:CBN, 2:CB0, 3:NCB, 4:OCB */
#define NUM_AMINOACID   20  /* There is assumed to be 20 amino acid            */
#define NUM_EAA_MATRIX 400  /* There is assumed to be 20 * 20 elements      */
#define NUM_NUCLEOTIDE  4   /* There is assumed to be 4 nucleotides            */
/* #define NUM_CODON       61   There is assumed to be 61 codons                */
#define NUM_COL_NEIGH   27  /* There is assumed to be 9, 9, 9 columns for each
                               codon                                           */
#define NUM_CORD_NEIGH  12  /* There is assumed to be 3, 3, 3, 3 cordinates for
                               each amino acid                                 */                               
#define NUM_ALL_CORD_AA 15  /* There is assumed to be 3, 3, 3, 3, 3 cordinates 
                               for each amino acid                             */
#define NUM_POISSON     16  /* There is assumed to be 16 Poisson Values        */                                                              
#define NUM_ROSETTA     93  /* There is assumed to be 93 rows of rosetta stone */
#define NUM_COL_SEQ_60  60  /* There is assumed to be 60 columns for each seq  */
#define NUM_COL_SEQ_30  30  /* There is assumed to be 60 columns for each seq  */
#define NUM_COL_SEQ_20  20  /* There is assumed to be 60 columns for each seq  */
#define NUM_COL_SEQ_15  15  /* There is assumed to be 60 columns for each seq  */

#define ERROR_STOPCODON 1   /* There is a stop codon in the given sequence     */
#define FUNCTION_BEGIN  0
#define FUNCTION_END    1

/*
   SANG CHUL CHOI: Some Definitions */
#define LEN_FILENAME 256
#define STRPDB "/pdb"
#define POTENTIAL "/drpro.test.int"
#define F84FREQSIZE    12   /* SEE the function F84_freq for the               */
                            /*                                  implementation */
#define CF_STR_SIZE 256

#define Dist(x1,y1,z1,x2,y2,z2) sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));

#define ERR_MEMORY         1000
#define ERR_INVDNACHR      1001
#define ERR_INVPROCHR      1002
#define ERR_NINF           1003
#define ERR_NAN            1004
#define ERR_NOTFINITE      1005
#define ERR_NOTIMPLEMENTED 1006
#define ERR_FILE           1007
#define ERR_READ           1008
#define ERR_GRID           1009
#define ERR_ASSERT         1010
#define ERR_MATH           1011
#define ERR_TOOMANYZEROS   1012
#define ERR_LOGZERO        1013
#define ERR_SAME           1014
#define ERR_INVBASE        1015
#define ERR_BLOCK          1016
#define ERR_GIBBS          1017
#define ERR_INVARG         1018
#define ERR_TRANSLATION    1019


enum { 
  PSI_PARAM_PI_A, 
  PSI_PARAM_PI_C,
  PSI_PARAM_PI_G,
  PSI_PARAM_PI_T,
  PSI_PARAM_S,
  PSI_PARAM_P 
};


#endif /* !PSI_DEFS_H_ */
