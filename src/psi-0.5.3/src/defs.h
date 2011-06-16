/* defs.h -- some important definitions
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

#ifndef _DEFS_H_
#define _DEFS_H_ 1

//#define BURNIN       10   /* This is how long we wait until we take a sample */
//#define SAMPLE_PATH  11   /* Number of sampled paths                         */
//#define SAMPFREQ     10   /*Number of MCMC cycles between samples            */

#define GIBBS_BURN     10
#define GIBBS_FREQ     10
#define GIBBS_SIZE    100  /* 100 */
#define BURNIN         50  /* This is how long we wait until we take a sample */
#define SAMPLE_FREQ   100  /* 100 */
#define SAMPLE_SIZE  1000

#define SAMPLE_PATH  1000   /* Number of sampled paths                         */
#define SAMPFREQ     500    /*Number of MCMC cycles between samples            */

#define MAXITERATION  (BURNIN +((SAMPLE_PATH-1)*SAMPFREQ)) 
                            /* This is the # of iterations of the MCMC routine */
#define NODE_UPDATE    1    /* how many cycles between node updates            */
#define DELTA_PATH     5    /* For fixed node sequences, how many times per    */
                            /* MCMC cycle should each branch have a randomly   */
                            /* selected site path proposed?                    */
#define ITERATE        1    /* The number of times that seq will be checked    */
#define SEQSET        100   /* How many sequences we simulate for S and P 
                               update                                          */
#define SIZEOFGRID     9    /* How big of a grid we will use for the 
                               simulation                                      */
#define LOW_A         0.14  /* Set these all to be the same now, but later     */
#define HIGH_A        0.36
#define LOW_C         0.14  /* Set these all to be the same now, but later     */
#define HIGH_C        0.36
#define LOW_G         0.14  /* Set these all to be the same now, but later     */
#define HIGH_G        0.36
#define LOW_T         0.14  /* Set these all to be the same now, but later     */
#define HIGH_T        0.36  /* we can have separate range values for each nuc  */

#define PI_MAX           1 
#define PI_MIN           0 
#define K1           1.857  
#define L1            1.0  
#define SIM_PI          1
/***************************/
#define FREQGRIDSIZE   12   /* SEE the function SeqFreq for the implementation */
//#define ADENINE    0.281118 /* If the user defines, then this is Pi_A          */
//#define CYTOSINE   0.204641 /* If the user defines, then this is Pi_C          */
//#define GUANINE    0.267932 /* If the user defines, then this is Pi_G          */
#define ADENINE      0.25   /* If the user defines, then this is Pi_A          */
#define CYTOSINE     0.25   /* If the user defines, then this is Pi_C          */
#define GUANINE      0.25   /* If the user defines, then this is Pi_G          */
/***************************/
#define SIM_T          0    /* [0,1] tells whether Time will be updated        */
#define INIT_T       1.0    /* The evolutionary time separating the sequences  */
#define DELTA_T      0.01   /* This is the step size for the time proposal     */
#define T_MIN        0.0    /* Uniform prior MIN value for time                */
#define T_MAX        1.0    /* Uniform prior MAX value for time                */
/***************************/
#define SIM_U          1    /* [0,1] tells whether Time will be updated        */
#define INIT_U       1.0    /* The evolutionary time searating the sequences   */
#define DELTA_U      0.2    /* This is the step size for the time proposal     */
#define U_MIN        0.0    /* Uniform prior MIN value for time                */
#define U_MAX        2.0    /* Uniform prior MAX value for time                */
/***************************/
#define SIM_K          1    /* [0,1] tells whether kappa will be updated       */
#define MULTI_K       -1    /* 0: One K on entire topo, 1:K_i per branch of 
                               topo                                            */
#define INIT_K       2.0    /* The initial seed value of kappa                 */
#define DELTA_K      1.0    /* This is the step size for the kappa proposal    */
#define K_MIN        0.0    /* Uniform prior MIN value for kappa               */
#define K_MAX       10.0    /* Uniform prior MAX value for kappa               */
/***************************/
#define SIM_W          1    /* [0,1] tells whether omega will be updated       */
#define MULTI_W       -1    /* 0: One W on entire topo, 1:W_i per branch of 
                               topo                                            */
#define INIT_W       1.0    /* Behaves like synonymous / nonsynonymous 
                               parameter                                       */
#define DELTA_W      0.5    /* This is the step size for the omega proposal    */
#define W_MIN        0.0    /* Uniform prior MIN value for omega               */
#define W_MAX        10.0   /* Uniform prior MAX value for omega               */
/***************************/
#define SIM_W_S        1    /* [0,1]tells whether omega(solvent) will be 
                               updated                                         */
#define INIT_W_S     0.0    /* Makes NRG's behave like a quantitative measure  */
#define W_S_MIN     -2.0    /* Uniform prior MIN value for omega_S             */
#define W_S_MAX      2.0    /* Uniform prior MAX value for omega_S             */
/***************************/
#define SIM_W_P        1    /* [0,1] tells whether omega(pairs) will be 
                               updated                                         */
#define INIT_W_P     0.0    /* Makes NRG's behave like a quantitative measure  */
#define W_P_MIN     -0.15   /* Uniform prior MIN value for omega_P             */
#define W_P_MAX      0.15   /* Uniform prior MAX value for omega_P             */
/***************************/

#define ANCESTRAL_NUC 1     /* 0: print out AA seqs, 1: print out nucleotide 
                               seqs                                            */
#define INTERACT 10.0       /* Distance at which AA's can biologically 
                               interact                                        */
#define BIGNUM 10000000.0   /* Used as an upper bound on the logarithm         */
#define INVBIGNUM 0.0000001 /* Used as a lowerbound on the logarithm           */
#define LOGBIGNUM 16.118095651 
                            /* This is the LN of the BIGNUM                    */
#define LOUD 0              /* How much STUFF do you want? 0 = little, 1 = 
                               TONS!                                           */
#define WHERE_ARE_WE 1
/***************************/

/*
   SANG CHUL CHOI: Some Definitions */
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
#define SUCCESS   0
#define FAIL     -1
#define CF_STR_SIZE 256
/* #define PACKAGE "drpro" */
/* #define VERSION "1.0"   */
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

#define GMEL_GRIDPOINT_S_BEGIN    -5.4
#define GMEL_GRIDPOINT_S_END       5.4
#define GMEL_GRIDPOINT_S_NUMBER   35
#define GMEL_GRIDPOINT_P_BEGIN    -0.33
#define GMEL_GRIDPOINT_P_END       0.33
#define GMEL_GRIDPOINT_P_NUMBER   35
#define GMEL_GRIDPOINT_N_BEGIN     0
#define GMEL_GRIDPOINT_N_END       1.0
#define GMEL_GRIDPOINT_N_NUMBER   11 /* 20 */
#define GMEL_FLATPRIOR_S_BEGIN    -5.0
#define GMEL_FLATPRIOR_S_END       5.0
#define GMEL_FLATPRIOR_P_BEGIN    -0.3
#define GMEL_FLATPRIOR_P_END       0.3
#define GMEL_FLATPRIOR_N_BEGIN       0
#define GMEL_FLATPRIOR_N_END         1

#define DELTA_PI      0.1
#define DELTA_PI_MU   100
#define DELTA_W_S     1.0    /* This is the step size for the omega_S proposal  */
#define DELTA_W_P     0.05   /* This is the step size for the omega_P proposal  */

/*
#define GMEL_N_MIN 0.000001
#define GMEL_N_MAX 0.999999
*/
#define GMEL_LINE_MAX 200    /* one line length of output file */


#define GMEL_STR_DAT               "../data/pdb1r48.A.dat"
#define GMEL_STR_ENERGY            "../data/energy.int"
#define GMEL_STR_STDOUT            "stdout"
#define GMEL_STR_PARAM             "param"
#define GMEL_STR_ALLDRAWS          "alldraws"
#define GMEL_STR_GIBBS             "gibbs"
#define GMEL_STR_SETUP             "setup"
#define GMEL_STR_BFACTOR           "bfactor"
#define GMEL_STR_Q                 "q"
#define GMEL_STR_PARAM_PREFIX      "pdb1r48.A.b."
#define GMEL_STR_PARAM_SUFFIX      ".out.p"


#endif /* !_DEFS_H_ */
