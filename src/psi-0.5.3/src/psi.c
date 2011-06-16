/* psi.c -- A Computer Program of Protein Structure Impact
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

extern char AMINOACID[21];

/* static functions */
static void usage (int status);
static int decode_switches (int argc, char **argv);

static void collect_garbage (void);
static int execute_bootstrap (void);
static int execute_convergence (void);
static int execute_scaled_selection_coefficient (void);
static int get_psr (double *psr, double **draws, int w);
static int setup_psi ();
static int unsetup_psi ();
static int prepare_filenames ();  

/*!\file
   \author Sang Chul Choi
   \date
 */

/*
#include "../nrgdrpro/gmel.h" 
*/

/* global variables */
/* The name the program was run with, stripped of any leading path. */
/*
char *program_name;
char *program_invocation_short_name;
*/

/* Option flags and variables */
static char *oname = NULL;                /* --output */
static FILE *ofile = NULL;
static char *datname = GMEL_STR_DAT;
static char *energyname = GMEL_STR_ENERGY;
static char *rname = GMEL_STR_STATE;
static char *pname = GMEL_STR_PARAM;                 /* --save-param */
static char *aname = GMEL_STR_ALLDRAWS;              /* --save-alldraws */
static char *gname = GMEL_STR_GIBBS;                 /* --save-gibbs */
static char *sname = GMEL_STR_SETUP;                 /* --save-setup */
static char *bname = GMEL_STR_BFACTOR;               /* --output-bf */
static char *qname = GMEL_STR_Q;                     /* --save-q */
static char *data_directory = NULL;                  /* --data-directory */
static char *output_directory = NULL;                /* --directory */
static char *auto_data_name = NULL;                  /* --auto-data */
static int want_interactive = 0;                     /* --interactive */
static int want_quiet = 0;                           /* --quiet, --silent */
static int want_brief = 0;                           /* --brief */
static int want_verbose = 0;                         /* --verbose */
static int want_dry_run = 0;                         /* --dry-run */
static int want_no_warn = 0;                         /* --no-warn */
static int want_progress = 0;                        /* --progress */

static int sample_burn = GMEL_SAMPLE_BURN;           /* --sample-burn */
static int sample_freq = GMEL_SAMPLE_FREQ;           /* --sample-freq */
static int sample_size = GMEL_SAMPLE_SIZE;           /* --sample-size */
static int want_is_fixed_run = 1;                    /* --is-fixed-run */
static int want_no_data = 0;                         /* --no-data */
static unsigned long int psi_seed = GMEL_RNG_INIT;   /* --fix-rng-seed */

static int gibbs_burn = GMEL_GIBBS_BURN;                  /* --gibbs-burn */
static int gibbs_freq = GMEL_GIBBS_FREQ;                  /* --gibbs-freq */
static int gibbs_size = GMEL_GIBBS_SIZE;                  /* --gibbs-size */
static double gridpoint_s_begin = GMEL_GRIDPOINT_S_BEGIN; /* -5.0  */
static double gridpoint_s_end = GMEL_GRIDPOINT_S_END;     /* 5.0   */
static int gridpoint_s_number = GMEL_GRIDPOINT_S_NUMBER;  /* 51    */
static double gridpoint_p_begin = GMEL_GRIDPOINT_P_BEGIN; /* -1.0; */
static double gridpoint_p_end = GMEL_GRIDPOINT_P_END;     /* 1.0;  */
static int gridpoint_p_number = GMEL_GRIDPOINT_P_NUMBER;  /* 41;   */
static double gridpoint_n_begin = GMEL_GRIDPOINT_N_BEGIN; /* 0;    */
static double gridpoint_n_end = GMEL_GRIDPOINT_N_END;     /* 1.0;  */
static int gridpoint_n_number = GMEL_GRIDPOINT_N_NUMBER;  /* 20;   */
static double flatprior_s_begin = GMEL_FLATPRIOR_S_BEGIN; /* -5.0; */
static double flatprior_s_end = GMEL_FLATPRIOR_S_END;     /* 5.0;  */
static double prior_normal_sd_s = (GMEL_FLATPRIOR_S_END - GMEL_FLATPRIOR_S_BEGIN) * 1 / 3.4641;
static double flatprior_p_begin = GMEL_FLATPRIOR_P_BEGIN; /* -1.0; */
static double flatprior_p_end = GMEL_FLATPRIOR_P_END;     /* 1.0;  */
static double prior_normal_sd_p = (GMEL_FLATPRIOR_P_END - GMEL_FLATPRIOR_P_BEGIN) * 1 / 3.4641;
static double flatprior_n_begin = GMEL_FLATPRIOR_N_BEGIN; /* 0;    */
static double flatprior_n_end =  GMEL_FLATPRIOR_N_END;    /* 1;    */
/* prior_dirichlet_alpha_n (1,1,1,1) */
static double delta_n = DELTA_PI;
static double delta_n_mu = DELTA_PI_MU;
static double delta_s = DELTA_W_S;
static double delta_s_sd = DELTA_W_S / 1.732051; 
static double delta_p = DELTA_W_P;
static double delta_p_sd = DELTA_W_P / 1.732051; 
static double init_s = GMEL_INIT_S;
static double init_p = GMEL_INIT_P;
static double init_a = GMEL_INIT_N;
static double init_c = GMEL_INIT_N;
static double init_g = GMEL_INIT_N;

static int want_disable_s = 0;                    /* --disable-s */
static int want_disable_p = 0;                    /* --disable-p */
static int want_disable_n = 0;                    /* --disable-n */
static int want_use_normal_prior_s = 0;           /* --use-normal-prior-s */
static int want_use_normal_prior_p = 0;           /* --use-normal-prior-p */
static int want_use_Dirichlet_prior_n = 1;        /* --use-normal-prior-p */
static int want_use_normal_proposal_s = 0;        /* --use-normal-proposal-s */
static int want_use_normal_proposal_p = 0;        /* --use-normal-proposal-p */
static int want_use_Dirichlet_proposal_n = 1;     /* --use-normal-proposal-p */
static int want_use_more_than_one_gridpoint = 0;  /* --use-more-than-one-gridpoint */
static int want_load_state = 0;                   /* --load-state */
static int limit_time_to_run = 0;                 /* --limit-time-to-run in minutes */
static int want_load_param = 0;                   /* --load-param */
static int want_save_alldraws = 0;                /* --save-alldraws */
static int want_load_gibbs = 0;                   /* --load-gibbs */
static int want_save_gibbs = 0;                   /* --save-gibbs */
static int want_load_setup = 0;                   /* --load-setup */
static int want_save_setup = 0;                   /* --save-setup */
static int want_load_q = 0;                       /* --load-q */
static int want_save_q = 0;                       /* --save-q */
static int want_debug_gibbs = 0;                  /* --debug-gibbs */
static int want_debug_bf = 0;                     /* --debug-bf */
static int want_debug_input = 0;                  /* --debug-input */
static int want_debug_est = 0;                    /* --debug-est */
static int want_debug_boot = 0;                   /* --debug-boot */
static int want_no_overwrite = 0;                 /* --no-overwrite */
static int bootstrapped_number = 100;             /* --bootstrapped_number */
static double tolerance_zeros = 0.5;              /* --tolerance-zeros */
static double limit_tau = 1.1;                    /* --limit-tau */
static double limit_psr = 1.01;                   /* --limit-psr */

static char *param_prefix = GMEL_STR_PARAM_PREFIX;/* --param-prefix */
static char *param_suffix = GMEL_STR_PARAM_SUFFIX;/* --param-suffix */
static int number_chains = GMEL_NUMBER_CHAINS;    /* --number-chains */
static int type_energy = PSI_ENERGY_JONES;        /* --type-energy */
static double nucleotide_frequency[NUM_NUCLEOTIDE] 
   = {0.295, 0.205, 0.205, 0.295};
   /* Asger's = {0.290, 0.193, 0.204, 0.313}; */
                                           /* --nuc-a-freq */
                                           /* --nuc-c-freq */
                                           /* --nuc-g-freq */
static double aa_codon_frequency[NUM_AMINOACID];
static int log_freq = GMEL_LOG_FREQ;              /* print the number of iterations */

static int translation_table = 1;                 /* --translation-table */

static int want_exe_all = 1;                      /* --exe-all */
static int want_exe_estimation = 0;               /* --exe-estimation */
static int want_exe_bayesfactor = 0;              /* --exe-bayesfactor */
static int want_exe_bootstrap = 0;                /* --exe-bootstrap */
static int want_exe_convergence = 0;              /* --exe-convergence */
static int want_exe_scaled_selection_coefficient = 0; 
                                           /* --exe-scaled-slection-coefficient */

/* getopt_long return codes */
enum { DUMMY_CODE=129
      ,BRIEF_CODE
      ,DRYRUN_CODE
      ,NOWARN_CODE
      ,DATA_DIRECTORY_CODE
      ,DIRECTORY_CODE
      ,GRIDPOINT_S_BEGIN_CODE
      ,GRIDPOINT_S_END_CODE
      ,GRIDPOINT_S_NUMBER_CODE
      ,GRIDPOINT_P_BEGIN_CODE
      ,GRIDPOINT_P_END_CODE
      ,GRIDPOINT_P_NUMBER_CODE
      ,GRIDPOINT_N_BEGIN_CODE
      ,GRIDPOINT_N_END_CODE
      ,GRIDPOINT_N_NUMBER_CODE
      ,FLATPRIOR_S_BEGIN_CODE
      ,FLATPRIOR_S_END_CODE
      ,FLATPRIOR_P_BEGIN_CODE
      ,FLATPRIOR_P_END_CODE
      ,FLATPRIOR_N_BEGIN_CODE
      ,FLATPRIOR_N_END_CODE
      ,DELTA_S_CODE
      ,DELTA_P_CODE
      ,DELTA_N_CODE
      ,INIT_S_CODE
      ,INIT_P_CODE
      ,INIT_A_CODE
      ,INIT_C_CODE
      ,INIT_G_CODE
      ,GIBBS_BURN_CODE
      ,GIBBS_FREQ_CODE
      ,GIBBS_SIZE_CODE
      ,PROGRESS_CODE
      ,FIXED_S_CODE
      ,FIXED_P_CODE
      ,FIXED_N_CODE
      ,DISABLE_S_CODE
      ,DISABLE_P_CODE
      ,USE_NORMAL_PRIOR_S_CODE
      ,USE_NORMAL_PRIOR_P_CODE
      ,USE_DIRICHLET_PRIOR_N_CODE
      ,USE_NORMAL_PROPOSAL_S_CODE
      ,USE_NORMAL_PROPOSAL_P_CODE
      ,USE_DIRICHLET_PROPOSAL_N_CODE
      ,EXE_ESTIMATION_CODE
      ,EXE_BAYESFACTOR_CODE
      ,GRIDPOINT_NUMBER_CODE
      ,LOAD_STATE_CODE
      ,LOAD_PARAM_CODE
      ,SAVE_PARAM_CODE
      ,SAVE_ALLDRAWS_CODE
      ,LOAD_GIBBS_CODE
      ,SAVE_GIBBS_CODE
      ,LOAD_SETUP_CODE
      ,SAVE_SETUP_CODE
      ,LOAD_Q_CODE
      ,SAVE_Q_CODE
      ,OUTPUT_BF_CODE
      ,DEBUG_GIBBS_CODE
      ,DEBUG_BF_CODE
      ,DEBUG_INPUT_CODE
      ,DEBUG_EST_CODE
      ,DEBUG_BOOT_CODE
      ,USE_MORE_THAN_ONE_GRIDPOINT_CODE
      ,TOLERANCE_ZEROS_CODE
      ,LIMIT_TAU_CODE
      ,LIMIT_PSR_CODE
      ,EXE_BOOTSTRAP_CODE
      ,BOOTSTRAPPED_NUMBER_CODE
      ,EXE_CONVERGENCE_CODE
      ,PARAM_PREFIX_CODE
      ,PARAM_SUFFIX_CODE
      ,NUMBER_CHAINS_CODE
      ,TYPE_ENERGY_CODE
      ,NUC_A_FREQ_CODE
      ,NUC_C_FREQ_CODE
      ,NUC_G_FREQ_CODE
      ,EXE_SCALED_SELECTION_COEFFICIENT_CODE
      ,TRANSLATION_TABLE_CODE
      ,EXE_ALL_CODE
      ,AUTO_DATA_CODE
      ,NO_DATA_CODE
      ,IS_FIXED_RUN_CODE
      ,FIX_RNG_SEED_CODE
      ,NO_OVERWRITE_CODE
      ,LIMIT_TIME_TO_RUN_CODE
};

static struct option const long_options[] =
{
  {"limit-time-to-run", required_argument, 0, LIMIT_TIME_TO_RUN_CODE},
  {"no-overwrite", no_argument, 0, NO_OVERWRITE_CODE},
  {"fix-rng-seed", required_argument, 0, FIX_RNG_SEED_CODE},
  {"is-fixed-run", no_argument, 0, IS_FIXED_RUN_CODE},
  {"no-data", no_argument, 0, NO_DATA_CODE},
  {"auto-data", required_argument, 0, AUTO_DATA_CODE},
  {"exe-all", no_argument, 0, EXE_ALL_CODE},
  {"translation-table", required_argument, 0, TRANSLATION_TABLE_CODE},
  {"exe-scaled-selection-coefficient", no_argument, 0, EXE_SCALED_SELECTION_COEFFICIENT_CODE},
  {"nuc-a-freq", required_argument, 0, NUC_A_FREQ_CODE},
  {"nuc-c-freq", required_argument, 0, NUC_C_FREQ_CODE},
  {"nuc-g-freq", required_argument, 0, NUC_G_FREQ_CODE},
  {"type-energy", required_argument, 0, TYPE_ENERGY_CODE},
  {"exe-convergence", no_argument, 0, EXE_CONVERGENCE_CODE},
  {"param-prefix", required_argument, 0, PARAM_PREFIX_CODE},
  {"param-suffix", required_argument, 0, PARAM_SUFFIX_CODE},
  {"number-chains", required_argument, 0, NUMBER_CHAINS_CODE},
  {"bootstrapped-number", required_argument, 0, BOOTSTRAPPED_NUMBER_CODE},
  {"exe-bootstrap", no_argument, 0, EXE_BOOTSTRAP_CODE},
  {"tolerance-zeros", required_argument, 0, TOLERANCE_ZEROS_CODE},
  {"limit-tau", required_argument, 0, LIMIT_TAU_CODE},
  {"limit-psr", required_argument, 0, LIMIT_PSR_CODE},
  {"use-more-than-one-gridpoint", no_argument, 0, USE_MORE_THAN_ONE_GRIDPOINT_CODE},
  {"debug-gibbs", no_argument, 0, DEBUG_GIBBS_CODE},
  {"debug-bf", no_argument, 0, DEBUG_BF_CODE},
  {"debug-input", no_argument, 0, DEBUG_INPUT_CODE},
  {"debug-est", no_argument, 0, DEBUG_EST_CODE},
  {"debug-boot", no_argument, 0, DEBUG_BOOT_CODE},
  {"output-bf", required_argument, 0, OUTPUT_BF_CODE},
  {"load-setup", required_argument, 0, LOAD_SETUP_CODE},
  {"save-setup", required_argument, 0, SAVE_SETUP_CODE},
  {"load-q", required_argument, 0, LOAD_Q_CODE},
  {"save-q", required_argument, 0, SAVE_Q_CODE},
  {"load-state", required_argument, 0, LOAD_STATE_CODE},
  {"load-param", required_argument, 0, LOAD_PARAM_CODE},
  {"save-param", required_argument, 0, SAVE_PARAM_CODE},
  {"save-alldraws", required_argument, 0, SAVE_ALLDRAWS_CODE},
  {"load-gibbs", required_argument, 0, LOAD_GIBBS_CODE},
  {"save-gibbs", required_argument, 0, SAVE_GIBBS_CODE},
  {"gridpoint-number", required_argument, 0, GRIDPOINT_NUMBER_CODE},
  {"exe-estimation", no_argument, 0, EXE_ESTIMATION_CODE},
  {"exe-bayesfactor", no_argument, 0, EXE_BAYESFACTOR_CODE},
  {"use-normal-prior-s", no_argument, 0, USE_NORMAL_PRIOR_S_CODE},
  {"use-normal-prior-p", no_argument, 0, USE_NORMAL_PRIOR_P_CODE},
  {"use-Dirichlet-prior-n", no_argument, 0, USE_DIRICHLET_PRIOR_N_CODE},
  {"use-normal-proposal-s", no_argument, 0, USE_NORMAL_PROPOSAL_S_CODE},
  {"use-normal-proposal-p", no_argument, 0, USE_NORMAL_PROPOSAL_P_CODE},
  {"use-Dirichlet-proposal-n", no_argument, 0, USE_DIRICHLET_PROPOSAL_N_CODE},
  {"disable-s", no_argument, 0, DISABLE_S_CODE},
  {"disable-p", no_argument, 0, DISABLE_P_CODE},
  {"gridpoint-s-begin", required_argument, 0, GRIDPOINT_S_BEGIN_CODE},
  {"gridpoint-s-end", required_argument, 0, GRIDPOINT_S_END_CODE},
  {"gridpoint-s-number", required_argument, 0, GRIDPOINT_S_NUMBER_CODE},
  {"gridpoint-p-begin", required_argument, 0, GRIDPOINT_P_BEGIN_CODE},
  {"gridpoint-p-end", required_argument, 0, GRIDPOINT_P_END_CODE},
  {"gridpoint-p-number", required_argument, 0, GRIDPOINT_P_NUMBER_CODE},
  {"gridpoint-n-begin", required_argument, 0, GRIDPOINT_N_BEGIN_CODE},
  {"gridpoint-n-end", required_argument, 0, GRIDPOINT_N_END_CODE},
  {"gridpoint-n-number", required_argument, 0, GRIDPOINT_N_NUMBER_CODE},
  {"flatprior-s-begin", required_argument, 0, FLATPRIOR_S_BEGIN_CODE},
  {"flatprior-s-end", required_argument, 0, FLATPRIOR_S_END_CODE},
  {"flatprior-p-begin", required_argument, 0, FLATPRIOR_P_BEGIN_CODE},
  {"flatprior-p-end", required_argument, 0, FLATPRIOR_P_END_CODE},
  {"flatprior-n-begin", required_argument, 0, FLATPRIOR_N_BEGIN_CODE},
  {"flatprior-n-end", required_argument, 0, FLATPRIOR_N_END_CODE},
  {"delta-s", required_argument, 0, DELTA_S_CODE},
  {"delta-p", required_argument, 0, DELTA_P_CODE},
  {"delta-n", required_argument, 0, DELTA_N_CODE},
  {"init-s", required_argument, 0, INIT_S_CODE},
  {"init-p", required_argument, 0, INIT_P_CODE},
  {"init-a", required_argument, 0, INIT_A_CODE},
  {"init-c", required_argument, 0, INIT_C_CODE},
  {"init-g", required_argument, 0, INIT_G_CODE},
  {"gibbs-burn", required_argument, 0, GIBBS_BURN_CODE},
  {"gibbs-freq", required_argument, 0, GIBBS_FREQ_CODE},
  {"gibbs-size", required_argument, 0, GIBBS_SIZE_CODE},
  {"progress", no_argument, 0, PROGRESS_CODE},
  {"sample-freq", required_argument, 0, 'f'},
  {"sample-size", required_argument, 0, 's'},
  {"sample-burn", required_argument, 0, 'b'},
  {"dat", required_argument, 0, 'd'},
  {"energy", required_argument, 0, 'e'},
  {"interactive", no_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},
  {"quiet", no_argument, 0, 'q'},
  {"silent", no_argument, 0, 'q'},
  {"brief", no_argument, 0, BRIEF_CODE},
  {"verbose", no_argument, 0, 'v'},
  {"dry-run", no_argument, 0, DRYRUN_CODE},
  {"no-warn", no_argument, 0, NOWARN_CODE},
  {"data-directory", required_argument, 0, DATA_DIRECTORY_CODE},
  {"directory", required_argument, 0, DIRECTORY_CODE},
  {"help", no_argument, 0, 'h'},
  {"version", no_argument, 0, 'V'},
  {NULL, 0, NULL, 0}
};

int
main (int argc, char **argv)
{
  int i;
  int r = EXIT_SUCCESS;
  parameter theta_dummy;
  time_t timer;
  time_t timer_end;
  set_program_name (argv[0]);
  timer = time(NULL);

  /* should be called before decode_switches */
  psi_seed = (unsigned long int) timer;

  atexit (collect_garbage);
  i = decode_switches (argc, argv);
  r = prepare_filenames ();  
  if (r != EXIT_SUCCESS)
    {
      /* exit without the analysis */
      exit (EXIT_SUCCESS);
    }
/* HERE WE WERE */
/*  psi_io_check_output (); */
  
  if (want_load_setup == 1) {
     psi_fatal ("no implementation of load_setup");
  }
  if (want_save_setup == 1) {
     psi_fatal ("no implementation of save_setup");
  }

  if (want_load_state == 1) {
     /* load random number generator, gibbs sample, parameters
     load_rng (rname);
     load_parameters (pname);  
     load_gridpoints (gname);
     */
  }

  if (want_dry_run == 1) {
     fprintf (stderr, "The current time is %s.\n", 
              asctime (localtime (&timer) ));
     exit (EXIT_SUCCESS);   
  }

/* This depends on the data file
  change_translation_table (translation_table);
*/

  if (oname == NULL)
    ofile = stdout;
  else 
    ofile = fopen (oname, "a");

  fprintf(ofile, "START: %s\n",asctime(localtime(&timer)));
  fprintf(stderr, "START: %s\n",asctime(localtime(&timer)));

  /* do the work */
  setup_psi ();

  if (want_load_gibbs == 1) {
     load_gridpoints (gname);
  }

  if (want_exe_all == 1) {
     setup_mc_option (want_is_fixed_run);
     setup_mc_chains (number_chains);
     if (number_chains == 2)
       {
         setup_mc_add_ip (0.25, 0.25, 0.25, 
                          GMEL_POS_INIT_S, GMEL_POS_INIT_P);
         setup_mc_add_ip (0.25, 0.25, 0.25, 
                          GMEL_NEG_INIT_S, GMEL_NEG_INIT_P);
       }
     else if (number_chains == 3)
       {
         setup_mc_add_ip (0.25, 0.25, 0.25, 
                          GMEL_POS_INIT_S, GMEL_POS_INIT_P);
         setup_mc_add_ip (0.25, 0.25, 0.25, 
                          GMEL_NEG_INIT_S, GMEL_NEG_INIT_P);
         setup_mc_add_ip (0.25, 0.25, 0.25,  0.0,  0.0);
       }
     else if (number_chains == 5)
       {
         setup_mc_add_ip (0.25, 0.25, 0.25, 
                          GMEL_POS_INIT_S, GMEL_POS_INIT_P);
         setup_mc_add_ip (0.25, 0.25, 0.25, 
                          GMEL_NEG_INIT_S, GMEL_NEG_INIT_P);
         setup_mc_add_ip (0.25, 0.25, 0.25,  0.0,  0.0);
         setup_mc_add_ip (0.25, 0.25, 0.25, 
                          GMEL_POS_INIT_S, GMEL_NEG_INIT_P);
         setup_mc_add_ip (0.25, 0.25, 0.25, 
                          GMEL_NEG_INIT_S, GMEL_POS_INIT_P);
       }
     else
       {
         psi_fatal ("number of chains should be 2, 3, or 5");
       }
     execute_multiple_chains ();
     unsetup_mc ();
  } else if (want_exe_estimation == 1) {
     r = remove (pname);
     if (r != EXIT_SUCCESS) {
        psi_warning ("could not remove the parameter file");
     }
     if (want_save_alldraws == 1) {
        remove (aname);
     }
     r = execute_estimation (want_disable_s, want_disable_p, 0, 0, 
                             &theta_dummy);
     if (r != EXIT_SUCCESS) {
        psi_fatal ("Could not execute the estimation procedure");
     }
  } else if (want_exe_bayesfactor == 1) {
     r = execute_bayesfactor ();
     if (r != EXIT_OTHERS) {
        if (r != EXIT_SUCCESS) {
           psi_fatal ("Could not execute the bayes factor procedure");
        }
     }
  } else if (want_exe_bootstrap == 1) {
     r == execute_bootstrap ();
     if (r != EXIT_SUCCESS) {
        psi_fatal ("Could not execute the bootstrap procedure");
     }
  } else if (want_exe_convergence == 1) {
     r == execute_convergence ();
     if (r != EXIT_SUCCESS) {
        psi_fatal ("Could not execute the convergence procedure");
     }
  } else if (want_exe_scaled_selection_coefficient == 1) {
     r == execute_scaled_selection_coefficient ();
     if (r != EXIT_SUCCESS) {
        psi_fatal ("Could not execute the scaled_selection_coefficient procedure");
     }
  } else {
     assert (0);
  }

  /* Energy Objects */
  if (want_save_gibbs == 1) {
     save_gridpoints (gname);
  }


  unsetup_psi ();

  timer_end = time (NULL);
  fprintf(ofile, "END: %s\n",asctime(localtime(&timer_end)));
  fprintf(stderr, "END: %s\n",asctime(localtime(&timer_end)));
  fprintf (ofile, "TIME: %lf\n", difftime (timer_end, timer));
  fprintf (stderr, "TIME: %lf\n", difftime (timer_end, timer));

  if (ofile != stdout)
    fclose (ofile);

  exit (EXIT_SUCCESS);
}

/* Set all the option flags according to the switches specified.
   Return the index of the first non-option argument.  */
static int
decode_switches (int argc, char **argv)
{
  static int i = 0;
  i++;
/*
  int i;
  fprintf (stderr, "argc=%d\n", argc);
  for (i = 0; i < argc; i++) {
    fprintf (stderr, "argv[%d]=%s ", i, argv[i]);
  }
  fprintf (stderr, "\n");
*/
  int c;
  /* struct stat ofile_stat; */
  /* int i; */

/*  ofile = stdout;
  char *str_option;
  char *token;
  double epsilon = 1; 
*/
  optind = 0;
  while ((c = getopt_long (argc, argv,
                           "i"  /* interactive */
                           "q"  /* quiet or silent */
                           "l:" /* lenght */
                           "v"  /* verbose */
                           "o:" /* output */
			   "b:"	/* sample-burn */
			   "f:"	/* sample-freq */
			   "s:"	/* sample-size */
                           "d:" /* dat */
                           "e:" /* energy */
                           "h"  /* help */
                           "V", /* version */
                           long_options, (int *) 0)) != EOF)
    {
      switch (c)
	{
        case LIMIT_TIME_TO_RUN_CODE:   /* --limit-time-to-run */
          limit_time_to_run = atoi(optarg);
          want_load_state = 1;
          break; 
        case NO_OVERWRITE_CODE:
          want_no_overwrite = 1;        /* --no-overwrite */
	  break;
        case FIX_RNG_SEED_CODE:         /* --fix-rng-seed */
          psi_seed = (unsigned long int) atoi(optarg);
          break; 
        case IS_FIXED_RUN_CODE:         /* --is-fixed-run */
          want_is_fixed_run = 1;
          break; 
        case NO_DATA_CODE:              /* --no-data */
          want_no_data = 1;
          break; 
        case AUTO_DATA_CODE:            /* --auto-data */
          auto_data_name = xstrdup(optarg);
          break; 
        case TRANSLATION_TABLE_CODE:    /* --translation-table */
          translation_table = atoi(optarg);
          break; 
        case NUC_A_FREQ_CODE:            /* --nuc-a-freq */
          nucleotide_frequency[PSI_DNA_A] = atof(optarg);
          nucleotide_frequency[PSI_DNA_T] = 1.0 - nucleotide_frequency[PSI_DNA_A] 
                                            - nucleotide_frequency[PSI_DNA_C] 
                                            - nucleotide_frequency[PSI_DNA_G];
          break; 
        case NUC_C_FREQ_CODE:            /* --nuc-c-freq */
          nucleotide_frequency[PSI_DNA_C] = atof(optarg);
          nucleotide_frequency[PSI_DNA_T] = 1.0 - nucleotide_frequency[PSI_DNA_A] 
                                            - nucleotide_frequency[PSI_DNA_C] 
                                            - nucleotide_frequency[PSI_DNA_G];
          break; 
        case NUC_G_FREQ_CODE:            /* --nuc-g-freq */
          nucleotide_frequency[PSI_DNA_G] = atof(optarg);
          nucleotide_frequency[PSI_DNA_T] = 1.0 - nucleotide_frequency[PSI_DNA_A] 
                                            - nucleotide_frequency[PSI_DNA_C] 
                                            - nucleotide_frequency[PSI_DNA_G];
          break; 
        case TYPE_ENERGY_CODE:           /* --type-energy */
          type_energy = atoi(optarg);
          break; 
        case NUMBER_CHAINS_CODE:         /* --number-chains */
          number_chains = atoi(optarg);                    
          break; 
        case PARAM_SUFFIX_CODE:          /* --param-suffix */
          param_suffix = xstrdup(optarg);
          break; 
        case PARAM_PREFIX_CODE:          /* --param-prefix */
          param_prefix = xstrdup(optarg);
          break; 
        case BOOTSTRAPPED_NUMBER_CODE:   /* --bootstrapped-number */
          bootstrapped_number = atoi(optarg);
          break; 
        case TOLERANCE_ZEROS_CODE:       /* --tolerance-zeros */
          tolerance_zeros = atof(optarg);
          break; 
        case LIMIT_TAU_CODE:             /* --limit-tau */
          limit_tau = atof(optarg);
          break; 
        case LIMIT_PSR_CODE:             /* --limit-psr */
          limit_psr = atof(optarg);
          break; 
        case USE_MORE_THAN_ONE_GRIDPOINT_CODE:/* --use-more-than-one-gridpoint */
          want_use_more_than_one_gridpoint = 1;
	  break;
        case DEBUG_GIBBS_CODE:
          want_debug_gibbs = 1;          /* --debug-gibbs */
	  break;
        case DEBUG_BF_CODE:
          want_debug_bf = 1;             /* --debug-bf */
	  break;
        case DEBUG_INPUT_CODE:
          want_debug_input = 1;          /* --debug-input */
	  break;
        case DEBUG_EST_CODE:
          want_debug_est = 1;            /* --debug-est */
	  break;
        case DEBUG_BOOT_CODE:
          want_debug_boot = 1;           /* --debug-boot */
	  break;
        case EXE_ALL_CODE:               /* --exe-all */
          want_exe_all = 1;
	  want_exe_estimation = 0;
          want_exe_bayesfactor = 0;
          want_exe_bootstrap = 0;
          want_exe_convergence = 0;          
          want_exe_scaled_selection_coefficient = 0;
	  break;
        case EXE_ESTIMATION_CODE:        /* --exe-estimation */
          want_exe_all = 0;
	  want_exe_estimation = 1;
          want_exe_bayesfactor = 0;
          want_exe_bootstrap = 0;
          want_exe_convergence = 0;          
          want_exe_scaled_selection_coefficient = 0;
          number_chains = 1;
	  break;
        case EXE_BAYESFACTOR_CODE:       /* --exe-bayesfactor */
          want_exe_all = 0;
	  want_exe_estimation = 0;
          want_exe_bayesfactor = 1;
          want_exe_bootstrap = 0;
          want_exe_convergence = 0;          
          want_exe_scaled_selection_coefficient = 0;
          number_chains = 1;
	  break;
/*
          want_use_normal_prior_s = 1; 
          want_use_normal_prior_p = 1;
          want_use_Dirichlet_prior_n = 1;
          want_use_normal_proposal_s = 1;
          want_use_normal_proposal_p = 1;
          want_use_Dirichlet_proposal_n = 1;
*/
        case EXE_BOOTSTRAP_CODE:         /* --exe-bootstrap */
          want_exe_all = 0;
	  want_exe_estimation = 0;
          want_exe_bayesfactor = 0;
          want_exe_bootstrap = 1;
          want_exe_convergence = 0;          
          want_exe_scaled_selection_coefficient = 0;
          want_load_param = 1;
          want_load_q = 1;
          want_save_q = 0;
	  break;
/*
          want_use_normal_prior_s = 1; 
          want_use_normal_prior_p = 1;
          want_use_Dirichlet_prior_n = 1;
          want_use_normal_proposal_s = 1;
          want_use_normal_proposal_p = 1;
          want_use_Dirichlet_proposal_n = 1;
*/
        case EXE_CONVERGENCE_CODE:       /* --exe-convergence */
          want_exe_all = 0;
	  want_exe_estimation = 0;
          want_exe_bayesfactor = 0;
          want_exe_bootstrap = 0;
          want_exe_convergence = 1;          
          want_exe_scaled_selection_coefficient = 0;
          break; 
        case EXE_SCALED_SELECTION_COEFFICIENT_CODE:       
                                         /* --exe-scaled-selection-coefficient */
          want_exe_all = 0;
	  want_exe_estimation = 0;
          want_exe_bayesfactor = 0;
          want_exe_bootstrap = 0;
          want_exe_convergence = 1;          
          want_exe_scaled_selection_coefficient = 1;          
          break; 
        case USE_NORMAL_PRIOR_S_CODE:    /* --use-normal-prior-s */
	  want_use_normal_prior_s = 1;
	  break;
        case USE_NORMAL_PRIOR_P_CODE:    /* --use-normal-prior-p */
	  want_use_normal_prior_p = 1;
	  break;
        case USE_DIRICHLET_PRIOR_N_CODE: /* --use-Dirichlet-prior-n */
	  want_use_Dirichlet_prior_n = 1;
	  break;
        case USE_NORMAL_PROPOSAL_S_CODE: /* --use-normal-proposal-s */
	  want_use_normal_proposal_s = 1;
	  break;
        case USE_NORMAL_PROPOSAL_P_CODE: /* --use-normal-proposal-p */
	  want_use_normal_proposal_p = 1;
	  break;
        case USE_DIRICHLET_PROPOSAL_N_CODE:   /* --use-Dirichlet-proposal-n */
	  want_use_Dirichlet_proposal_n = 1;
	  break;
	case 'i':                        /* --interactive */
	  want_interactive = 1;
	  break;
        case DISABLE_S_CODE:             /* --disable-s */
          want_disable_s = 1;
          break; 
        case DISABLE_P_CODE:             /* --disable-p */
          want_disable_p = 1;
          break; 
        case GRIDPOINT_S_BEGIN_CODE:     /* --gridpoint-s-begin */
          gridpoint_s_begin = atof(optarg);
          break; 
        case GRIDPOINT_S_END_CODE:       /* --gridpoint-s-end */
          gridpoint_s_end = atof(optarg);
          break; 
        case GRIDPOINT_S_NUMBER_CODE:    /* --gridpoint-s-number */
          gridpoint_s_number = atoi(optarg);
          break; 
        case GRIDPOINT_P_BEGIN_CODE:     /* --gridpoint-p-begin */
          gridpoint_p_begin = atof(optarg);
          break; 
        case GRIDPOINT_P_END_CODE:       /* --gridpoint-p-end */
          gridpoint_p_end = atof(optarg);
          break;
        case GRIDPOINT_P_NUMBER_CODE:    /* --gridpoint-p-number */
          gridpoint_p_number = atoi(optarg);
          break; 
        case GRIDPOINT_N_BEGIN_CODE:     /* --gridpoint-n-begin */
          gridpoint_n_begin = atof(optarg);
          break; 
        case GRIDPOINT_N_END_CODE:       /* --gridpoint-n-end */
          gridpoint_n_end = atof(optarg);
          break; 
        case GRIDPOINT_N_NUMBER_CODE:    /* --gridpoint-n-number */
          gridpoint_n_number = atoi(optarg);
          break; 
        case FLATPRIOR_S_BEGIN_CODE:     /* --flatprior-s-begin */
          flatprior_s_begin = atof(optarg);
          prior_normal_sd_s = (flatprior_s_end - flatprior_s_begin) * 1 / 3.4641;  /* (b - a)^2 / 12 */
          break; 
        case FLATPRIOR_S_END_CODE:       /* --flatprior-s-end */
          flatprior_s_end = atof(optarg);
          prior_normal_sd_s = (flatprior_s_end - flatprior_s_begin) * 1 / 3.4641;  /* (b - a)^2 / 12 */
          break; 
        case FLATPRIOR_P_BEGIN_CODE:     /* --flatprior-p-begin */
          flatprior_p_begin = atof(optarg);
          prior_normal_sd_p = (flatprior_p_end - flatprior_p_begin) * 1 / 3.4641;  /* (b - a)^2 / 12 */
          break; 
        case FLATPRIOR_P_END_CODE:       /* --flatprior-p-end */
          flatprior_p_end = atof(optarg);
          prior_normal_sd_p = (flatprior_p_end - flatprior_p_begin) * 1 / 3.4641;  /* (b - a)^2 / 12 */
          break; 
        case FLATPRIOR_N_BEGIN_CODE:     /* --flatprior-n-begin */
          flatprior_n_begin = atof(optarg);
          break; 
        case FLATPRIOR_N_END_CODE:       /* --flatprior-n-end */
          flatprior_n_end = atof(optarg);
          break; 
        case DELTA_S_CODE:               /* --delta-s */
          delta_s = atof(optarg);
          delta_s_sd = delta_s / 1.732051; 
          break; 
        case DELTA_P_CODE:               /* --delta-p */
          delta_p = atof(optarg);
          delta_p_sd = delta_p / 1.732051; 
          break; 
        case DELTA_N_CODE:               /* --delta-n */
          delta_n = atof(optarg);
          delta_n_mu = 1/(delta_n * delta_n);
          break; 
        case INIT_S_CODE:                /* --init-s */
          init_s = atof(optarg);
          break; 
        case INIT_P_CODE:                /* --init-p */
          init_p = atof(optarg);
          break; 
        case INIT_A_CODE:                /* --init-a */
          init_a = atof(optarg);
          break; 
        case INIT_C_CODE:                /* --init-c */
          init_c = atof(optarg);
          break; 
        case INIT_G_CODE:                /* --init-g */
          init_g = atof(optarg);
          break; 
        case GIBBS_BURN_CODE:            /* --gibbs-burn */
          gibbs_burn = atoi(optarg);
          break; 
        case GIBBS_FREQ_CODE:            /* --gibbs-freq */
          gibbs_freq = atoi(optarg);
          break; 
        case GIBBS_SIZE_CODE:            /* --gibbs-size */
          gibbs_size = atoi(optarg);
          break; 
        case 'f':                        /* --sample-freq */
          sample_freq = atoi(optarg);
          break; 
        case 's':                        /* --sample-size */
          sample_size = atoi(optarg);
          break; 
        case 'b':                        /* --sample-burn */
          sample_burn = atoi(optarg);
          break; 
        case LOAD_STATE_CODE:            /* --load-state */
          rname = xstrdup(optarg);
          want_load_state = 1;
          break;
        case LOAD_PARAM_CODE:            /* --load-param */
          pname = xstrdup(optarg);
          want_load_param = 1;
          break;
        case SAVE_PARAM_CODE:            /* --save-param */
          pname = xstrdup(optarg);
          break;
        case SAVE_ALLDRAWS_CODE:         /* --save-alldraws */
          aname = xstrdup(optarg);
          want_save_alldraws = 1;
          break;
        case LOAD_GIBBS_CODE:            /* --load-gibbs */
          gname = xstrdup(optarg);
          want_load_gibbs = 1;
          break;
        case SAVE_GIBBS_CODE:            /* --save-gibbs */
          gname = xstrdup(optarg);
          want_save_gibbs = 1;
          break;
        case LOAD_SETUP_CODE:            /* --load-setup */
          sname = xstrdup(optarg);
          want_load_setup = 1;
          break;
        case SAVE_SETUP_CODE:            /* --save-setup */
          sname = xstrdup(optarg);
          want_save_setup = 1;
          break;
        case LOAD_Q_CODE:                /* --load-q */
          qname = xstrdup(optarg);
          want_load_q = 1;
          want_save_q = 0;
          break;
        case SAVE_Q_CODE:                /* --save-q */
          qname = xstrdup(optarg);
          want_save_q = 1;
          want_load_q = 0;
          break;
        case OUTPUT_BF_CODE:             /* --output-bf */
          bname = xstrdup(optarg);
          break;
        case 'd':                        /* --dat */
          datname = xstrdup(optarg);
          break;
        case 'e':                        /* --energy */
          energyname = xstrdup(optarg);
          break;
	case 'o':		         /* --output */
	  oname = xstrdup(optarg);
	  break;
	case PROGRESS_CODE:              /* --progress */
	  want_progress = 1;
	  break;
	case 'q':                        /* --quiet, --silent */
	  want_quiet = 1;
	  break;
	case BRIEF_CODE:                 /* --brief */
	  want_brief = 1;
	  break;
	case 'v':                        /* --verbose */
	  want_verbose = 1;
	  break;
	case DRYRUN_CODE:                /* --dry-run */
	  want_dry_run = 1;
	  break;
	case NOWARN_CODE:                /* --no-warn */
	  want_no_warn = 1;
	  break;
	case DATA_DIRECTORY_CODE:        /* --data-directory */
	  data_directory = xstrdup(optarg);
	  break;
	case DIRECTORY_CODE:             /* --directory */
	  output_directory = xstrdup(optarg);
	  break;
	case 'V':
	  printf ("drprot %s - Sun Apr  9 22:34:45 EDT 2006\n", VERSION);
	  exit (0);
	case 'h':
fprintf (stderr, "usage: %d", c);
	  usage (0);
	default:
fprintf (stderr, "usage: %d", c);
	  usage (EXIT_FAILURE);
	}
    }

  return optind;
}

static void
usage (int status)
{
  printf ("%s - \
C Version Of DrPro/DrEvol\n", program_name);
  printf ("Usage: %s [OPTION]... [FILE]...\n", program_name);
  printf ("\
Options:\n\
\n\
Input files:\n\
  --no-overwrite             default is overwrite the output\n\
  --auto-data NAME           [N/A] 1hi2A or 1r48_\n\
  --directory NAME           output direcotry\n\
  --data-directory NAME      data directory\n\
  -d, --dat NAME             [%s] read DrEvol's DAT file\n\
  -e, --energy NAME          [%s] read DrEvol's statistical potential energy file\n\
  --type-energy NUM          [%d] 0: David Jones' method\n\
                                 1: Berjerano's VLMM method\n\
  --nuc-a-freq NUM           [%lf] A freq for mutation rate (VLMM only)\n\
  --nuc-c-freq NUM           [%lf] C freq for mutation rate (VLMM only)\n\
  --nuc-g-freq NUM           [%lf] G freq for mutation rate (VLMM only)\n\
  --limit-time-to-run MINs   [%d] time to run in minutes\n\
  --load-state NAME          [%s] load and save the states of the analysis\n\
  --load-param NAME          [%s] load stored parameter for Bayes factor\n\
  --load-gibbs NAME          [%s] load stored Gibbs sampled sequences\n\
  --load-q NAME              [%s] load q file and off saving q\n\
  --load-setup NAME          [%s] load setup file and ignore all options\n\
\n\
Output files:\n\
  -o, --output NAME          [%s] send output to NAME instead of standard output\n\
  --save-param NAME          [%s] load stored parameter for Bayes factor\n\
  --save-gibbs NAME          [%s] save Gibbs sampled sequences\n\
  --save-q NAME              [%s] save q in file and off loading q\n\
  --save-setup NAME          [%s] save setup file\n\
  --output-bf NAME           [%s] save Bayes factor\n\
\n\
Executes:\n\
  --exe-all                  [TRUE] find tau, check conv., and est. of Bayes factor\n\
  --exe-estimation           [FALSE] execute the estimation of S, P, and N\n\
  --exe-bayesfactor          [FALSE] execute the bayes factor analysis\n\
                                     will use noraml and Dirichlet prior/proposal\n\
  --exe-bootstrap            [FALSE] bootstrap BF, use by loading param and q,\n\
                                     optionally loading gibbs\n\
                                     will use noraml and Dirichlet prior/proposal\n\
  --bootstrapped-number NUM  [%d] number of replications of bootstrap\n\
  --exe-convergence          [FALSE] check the convergence with multiple chains\n\
  --number-chains NUM        [%d] the number of chains\n\
  --param-prefix NAME        [%s] the prefix of param of a chain\n\
  --param-suffix NAME        [%s] the suffix of param of a chain\n\
  --exe-scaled-selection-coefficient [FALSE] get the scaled selection coefficient\n\
                                             using VLMM\n\
  --translation-table NUM    [%d] specify if non-standard translation table(1,4,6)\n\
\n\
Chain length setup:\n\
  -b, --sample-burn NUMBER   [%d] burn-in period\n\
  -f, --sample-freq NUMBER   [%d] posterior sampling frequency\n\
  -s, --sample-size NUMBER   [%d] posterior sample size\n\
  --is-fixed-run                  no check of mcmc options and psi\n\
  --fix-rng-seed NUMBER           set to the current time in seconds if not specified\n\
\n\
Estimation setup:\n\
  --disable-s                [FALSE] disable estimation of S\n\
  --disable-p                [FALSE] disable estimation of P\n\
  --init-s NUM               [%lf] init value of S\n\
  --init-p NUM               [%lf] init value of P\n\
  --init-a NUM               [%lf] init value of A\n\
  --init-c NUM               [%lf] init value of C\n\
  --init-g NUM               [%lf] init value of G\n\
\n\
Prior:\n\
  --flatprior-s-begin NUM    [%lf] beginning of flatprior of S\n\
  --flatprior-s-end NUM      [%lf] ending of flatprior of S\n\
  --flatprior-p-begin NUM    [%lf] beginning of flatprior of P\n\
  --flatprior-p-end NUM      [%lf] ending of flatprior of P\n\
  --flatprior-n-begin NUM    [%lf] beginning of flatprior of A, C, G or T\n\
  --flatprior-n-end NUM      [%lf] ending of flatprior of A, C, G or T\n\
  --use-normal-prior-p       [FALSE] enable normal prior for P\n\
  --use-normal-prior-s       [FALSE] enable normal prior for S\n\
  --use-Dirichlet-prior-n    [TRUE] NOT IMPLELMENTED\n\
\n\
Proposal:\n\
  --use-Dirichlet-proposal-n [TRUE] enable Dirichlet proposal density for N\n\
  --use-normal-proposal-p    [FALSE] enable normal proposal density for P\n\
                                     Use only with normal prior\n\
  --use-normal-proposal-s    [FALSE] enable normal proposal density for S\n\
                                     Use only with normal prior\n\
  --delta-s NUM              [%lf] delta of S\n\
  --delta-p NUM              [%lf] delta of P\n\
  --delta-n NUM              [%lf] delta of A, C, G or T\n\
\n\
Gibbs and grid:\n\
  --gibbs-burn NUMBER        [%d] gibbs burn-in period = this number x DNA length\n\
  --gibbs-freq NUMBER        [%d] gibbs sampling frequency = this number x DNA length\n\
  --gibbs-size NUMBER        [%d] gibbs sample size\n\
  --gridpoint-n-begin NUM    [%lf] grid-point of beginning of A, C, G or T\n\
  --gridpoint-n-end NUM      [%lf] grid-point of ending of A, C, G or T\n\
  --gridpoint-n-number NUM   [%d] number of grid-points of A, C, G or T\n\
  --gridpoint-p-begin NUM    [%lf] grid-point of beginning of P\n\
  --gridpoint-p-end NUM      [%lf] grid-point of ending of P\n\
  --gridpoint-p-number NUM   [%d] number of grid-points of P\n\
  --gridpoint-s-begin NUM    [%lf] grid-point of beginning of S\n\
  --gridpoint-s-end NUM      [%lf] grid-point of ending of S\n\
  --gridpoint-s-number NUM   [%d] number of grid-points of S\n\
  --use-more-than-one-gridpoint [FALSE] use middle points between two values in BF\n\
  --tolerance-zeros NUM      [%lf] how many zeros are allowed in sum of exp in BF\n\
  --limit-tau NUM            [%lf] bigger than 1.0 how much strict to choose autocorrelation time\n\
  --limit-psr NUM            [%lf] bigger than and near 1.0, less than 100 how much strict to choose convergence\n\
\n\
  -i, --interactive          prompt for confirmation\n\
  --progress                 show progress bar\n\
  --dry-run                  take no real actions\n\
  --no-warn                  disable warnings\n\
  -q, --quiet, --silent      inhibit usual output\n\
  --brief                    shorten output\n\
  --verbose                  print more information\n\
  -h, --help                 display this help and exit\n\
  -V, --version              output version information and exit\n\
", datname, energyname, type_energy, 
    nucleotide_frequency[PSI_DNA_A],
    nucleotide_frequency[PSI_DNA_C],
    nucleotide_frequency[PSI_DNA_G],
    limit_time_to_run, 
    rname, pname, gname, qname, sname, 
    oname, pname, gname, qname, sname, bname,
    bootstrapped_number,
    number_chains, param_prefix, param_suffix, translation_table,
    sample_burn, sample_freq, sample_size,
    init_s, init_p, init_a, init_c, init_g,
    flatprior_s_begin, flatprior_s_end,
    flatprior_p_begin, flatprior_p_end,
    flatprior_n_begin, flatprior_n_end,
    delta_s, delta_p, delta_n, 
    gibbs_burn, gibbs_freq, gibbs_size, 
    gridpoint_n_begin, gridpoint_n_end, gridpoint_n_number,
    gridpoint_p_begin, gridpoint_p_end, gridpoint_p_number,
    gridpoint_s_begin, gridpoint_s_end, gridpoint_s_number,
    tolerance_zeros, limit_tau, limit_psr
);
  exit (status);
}

static void 
collect_garbage (void)
{
  if (strcmp(datname, GMEL_STR_DAT) != 0) {
    free (datname);
    datname = NULL;
  }
  if (strcmp(energyname, GMEL_STR_ENERGY) != 0) {
    free (energyname);
    energyname = NULL;
  }
  if (oname != NULL) {
    free (oname);
    oname = NULL;
  }
  if (strcmp(rname, GMEL_STR_STATE) != 0) {
    free (rname);
    rname = NULL;
  }
  if (strcmp(pname, GMEL_STR_PARAM) != 0) {
    free (pname);
    pname = NULL;
  }
  if (strcmp(aname, GMEL_STR_ALLDRAWS) != 0) {
    free (aname);
    aname = NULL;
  }
  if (strcmp(gname, GMEL_STR_GIBBS) != 0) {
    free (gname);
    gname = NULL;
  }
  if (strcmp(sname, GMEL_STR_SETUP) != 0) {
    free (sname);
    sname = NULL;
  }
  if (strcmp(bname, GMEL_STR_BFACTOR) != 0) {
    free (bname);
    bname = NULL;
  }
  if (strcmp(qname, GMEL_STR_Q) != 0) {
    free (qname);
    qname = NULL;
  }
  if (strcmp(param_prefix, GMEL_STR_PARAM_PREFIX) != 0) {
    free (param_prefix);
    param_prefix = NULL;
  }
  if (strcmp(param_suffix, GMEL_STR_PARAM_SUFFIX) != 0) {
    free (param_suffix);
    param_suffix = NULL;
  }
}

static int
execute_bootstrap (void)
{
  int r = EXIT_SUCCESS;
  int i;
  char new_pname[GMEL_LINE_MAX];
  char new_qname[GMEL_LINE_MAX];
  char subt_name[GMEL_LINE_MAX];

  strcpy (new_pname, pname);
  strcat (new_pname, ".shuffle"); 
  strcpy (new_qname, qname);
  strcat (new_qname, ".shuffle"); 
  
  for (i = 0; i < bootstrapped_number; i++) {
     /* suffle the two files: param and q */
     if (want_debug_boot == 1) {     
        fprintf (stderr, "shuffling: pname - %s, new_pname - %s\n", pname, new_pname);
        fprintf (stderr, "shuffling: qname - %s, new_qname - %s\n", qname, new_qname);
     }
     r = psi_shuffle_file (pname, new_pname, sample_size);
     if (r != EXIT_SUCCESS) {
        fprintf (stderr, "%s:%d:%s; %s %s -> %s\n",
                 __FILE__, __LINE__, __PRETTY_FUNCTION__,
                 "Could not shuffle the file", pname, new_pname);
        psi_fatal ("Could not shuffle the file: pname");
     }
     r = psi_shuffle_file (qname, new_qname, sample_size);
     if (r != EXIT_SUCCESS) {
        fprintf (stderr, "%s:%d:%s; %s %s -> %s\n",
                 __FILE__, __LINE__, __PRETTY_FUNCTION__,
                 "Could not shuffle the file", qname, new_qname);
        psi_fatal ("Could not shuffle the file: qname");
        return r;
     }

     /* substitutue the new pname and qname for the current ones */
     strcpy (subt_name, pname); 
     free (pname);
     pname = NULL;
     pname = xstrdup (new_pname);
     strcpy (new_pname, subt_name);

     strcpy (subt_name, qname); 
     free (qname);
     qname = NULL;
     qname = xstrdup (new_qname);
     strcpy (new_qname, subt_name);
     if (want_debug_boot == 1) {     
        fprintf (stderr, "shuffling: pname - %s, new_pname - %s\n", pname, new_pname);
        fprintf (stderr, "shuffling: qname - %s, new_qname - %s\n", qname, new_qname);
     }
     /* execute the bayes_factor analysis */
     r = execute_bayesfactor ();
     if (r != EXIT_SUCCESS) {
        fprintf (stderr, "%s:%d:%s; %s %d-th\n",
                 __FILE__, __LINE__, __PRETTY_FUNCTION__,
                 "Could not execute the bayes factor procedure", i);
        psi_fatal ("Could not execute the bayes factor procedure");
        return r;
     }
     strcpy (subt_name, pname); 
     free (pname);
     pname = NULL;
     pname = xstrdup (new_pname);
     strcpy (new_pname, subt_name);

     strcpy (subt_name, qname); 
     free (qname);
     qname = NULL;
     qname = xstrdup (new_qname);
     strcpy (new_qname, subt_name);
  }

  return r;
}

static int
execute_convergence (void)
{
  int r = EXIT_SUCCESS;
  int i;
  double **draws = NULL;
  double psr;
  
  draws = XMALLOC (double *, number_chains);
  for (i = 0; i < number_chains; i++) {
     draws[i] = XMALLOC (double, sample_size);
  }

  get_psr (&psr, draws, PSI_PARAM_PI_A);
  fprintf (ofile, "%lf\t", psr);
  get_psr (&psr, draws, PSI_PARAM_PI_C);
  fprintf (ofile, "%lf\t", psr);
  get_psr (&psr, draws, PSI_PARAM_PI_G);
  fprintf (ofile, "%lf\t", psr);
  get_psr (&psr, draws, PSI_PARAM_PI_T);
  fprintf (ofile, "%lf\t", psr);
  get_psr (&psr, draws, PSI_PARAM_S);
  fprintf (ofile, "%lf\t", psr);
  get_psr (&psr, draws, PSI_PARAM_P);
  fprintf (ofile, "%lf\n", psr);
 
  for (i = 0; i < number_chains; i++) {
     XFREE (draws[i]);
  }
  XFREE (draws);
  return r;
}

/* goshng */
static int 
execute_scaled_selection_coefficient (void)
{

  int r = EXIT_SUCCESS;
  double two_times_Ns = 0;
  double log_prob_i = 0;
  /* double sum_log_codon_i = 0; */
  double log_prob_j = 0;
  /* double sum_log_codon_j = 0; */
  int which_site_aa;
  int which_pos_codon;
  int codon[3];
  int i, j;
  int wildtype_j;
  int wildtype_aa_j;
  FILE *pfile = NULL;
  int *dna = dna_jones_measure ();
  int *pro = pro_jones_measure ();
  int len_dna = len_dna_jones_measure ();
  int len_pro = len_pro_jones_measure ();

  pfile = fopen (pname, "w");
  fprintf (pfile, "Ns\n");  

  /* int array: mutant_type j */
  int *mutant_j = XMALLOC (int, len_dna); 
  int *mutant_aa_j = XMALLOC (int, len_pro); 
  char *AAseq_str = XMALLOC (char, len_pro + 1);
  for (i = 0; i < len_dna; i++) {
     mutant_j[i] = dna[i];
  }
  for (i = 0; i < len_pro; i++) {
     mutant_aa_j[i] = pro[i];
  }
  r = pro_int2char (mutant_aa_j, AAseq_str, len_pro);
  if (r != EXIT_SUCCESS) {
    fprintf (stderr, "pro_int error: %s, %d\n", AAseq_str, len_pro);
  }
/*
  calc_energy_vlmm (AAseq_str, mutant_aa_j, &log_prob_i, 
                    &dummy_value, Dat->len_pro);
*/

  if (want_verbose == 1) {
     fprintf (ofile, "DNA Len: %d\n", len_dna);
     fprintf (ofile, "Pro Len: %d\n", len_pro);
     fprintf (ofile, "Protein: %s\n", AAseq_str);
     for (i = 0; i < NUM_AMINOACID; i++) {
        fprintf (ofile, "%c %lf\n", AMINOACID[i], aa_codon_frequency[i]);
     }
  }

  for (i = 0; i < len_dna; i++) {
              /* translate mutant_j to mutant_aa_j */
              which_site_aa = i / 3;
              which_pos_codon  = i % 3;
              switch (which_pos_codon) {
              // We need the correct full codon for this site
              case PSI_CODON_FIRST:
                 codon[1] = mutant_j[i + 1]; 
                 codon[2] = mutant_j[i + 2];
                 break;
              case PSI_CODON_SECOND:
                 codon[0] = mutant_j[i - 1]; 
                 codon[2] = mutant_j[i + 1];
                 break;
              case PSI_CODON_THIRD:
                 codon[0] = mutant_j[i - 2]; 
                 codon[1] = mutant_j[i - 1];
                 break;
              }
     for (j = 0; j < NUM_NUCLEOTIDE; j++) {
        if (mutant_j[i] != j) {
           wildtype_j = mutant_j[i];
           mutant_j[i] = j;
           /* calculate two_times_Ns */

           /* convert mutant_aa_j to AAseq_str */ 
           codon[which_pos_codon] = j;
           wildtype_aa_j = mutant_aa_j[which_site_aa];
           mutant_aa_j[which_site_aa] = aa_codon_is(codon);
           if (mutant_aa_j[which_site_aa] == wildtype_aa_j) {
              if (want_verbose == 1) {
                 fprintf (ofile, "synonymous\n");
              }
              mutant_j[i] = wildtype_j;
              continue;
           }
           if (mutant_aa_j[which_site_aa] == 20) {
              if (want_verbose == 1) {
                 fprintf (ofile, "stop codon\n");
              }
              mutant_j[i] = wildtype_j;
              mutant_aa_j[which_site_aa] = wildtype_aa_j;
              continue;
           }
           AAseq_str[which_site_aa] = AMINOACID[mutant_aa_j[which_site_aa]];

/*
           calc_energy_vlmm (AAseq_str, mutant_aa_j, &log_prob_j, 
                             &dummy_value, Dat->len_pro);
*/
           two_times_Ns = (log_prob_i/2.0 - log_prob_j/2.0) * (-0.5);
           fprintf (pfile, "%lf\n", two_times_Ns);  
           if (want_verbose == 1) {

              fprintf (ofile, "Mutant %3d %3d: %s\n", i, which_site_aa, AAseq_str);
           }
           mutant_j[i] = wildtype_j;
           mutant_aa_j[which_site_aa] = wildtype_aa_j;
           AAseq_str[which_site_aa] = AMINOACID[wildtype_aa_j];
        }
     }
  }

  XFREE (mutant_j);
  XFREE (mutant_aa_j);
  XFREE (AAseq_str);
  fclose (pfile);
  pfile = NULL;
  return r;
}

static int
get_psr (double *psr, double **draws, int w)
{
  int r = EXIT_SUCCESS;
  int i;
  char chain_name[GMEL_LINE_MAX];
  char str_number[GMEL_LINE_MAX];;
  
  for (i = 0; i < number_chains; i++) {
     strcpy (chain_name, param_prefix);
     sprintf (str_number, "%c", i + 48);
     strcat (chain_name, str_number);
     strcat (chain_name, param_suffix);
     r = gmel_read_draws (chain_name, draws[i], w, sample_size); 
     if (r != EXIT_SUCCESS) {
        fprintf (stderr, "%s:%d:%s; %s\n",
                 __FILE__, __LINE__, __PRETTY_FUNCTION__,
                 "Could not read the draws");
        psi_fatal ("Could not read the draws");
     }
  }

/* PRINT THE READ CHAINS for one PARAMETER
  int j;
  for (i = 0; i < number_chains; i++) {
     fprintf (stderr, "chains %d:\n", i);
     for (j = 0; j < sample_size; j++) {
        fprintf (stderr, "%lf ", draws[i][j]);
     }
     fprintf (stderr, "\n");
  }
*/

  r = gmel_conv_psr (psr, draws, number_chains, sample_size);
  if (r != EXIT_SUCCESS) {
     fprintf (stderr, "%s:%d:%s; %s\n",
              __FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Could not estimate the potential reduction scale");
     psi_fatal ("Could not estimate the potential reduction scale");
  }

  return r;
}

static int 
unsetup_psi ()
{
  if (want_load_state == 1) {
     /* save random number generator, gibbs sample, parameters
     save_parameters (); 
     */
     save_rng (rname);
     save_gridpoints (gname);
  }

  delete_gridpoints ();

  /* gibbs */
  gmel_gibbs_write_usage (ofile);
  unsetup_gibbs();

  /* unsetup_bf (); we don't call any function of bf.c */
  unsetup_mcmc ();
  unsetup_bf ();

  /* random number generator and energy data */
  finalize_drevol (PSI_ENERGY_JONES);
  fin_rng ();
  return EXIT_SUCCESS;
}

static int 
setup_psi ()
{
  int r;
  /* random number generator and energy data */
  if (want_load_state == 1) 
    {
      r = load_rng (rname);
      if (r == EXIT_FAILURE)
        init_rng (psi_seed);
    }
  else
    init_rng (psi_seed);

  psi_rng_info (ofile);
  setup_energy_jones (datname, energyname);
  initialize_drevol (PSI_ENERGY_JONES);

  /* gibbs */
  setup_gibbs (gibbs_size, gibbs_freq, gibbs_burn);
  fprintf (ofile, "GIBBS:\t%d\t%d\t%d\n", gibbs_burn, gibbs_freq, gibbs_size);

  /* sample */
  setup_mcmc_sample (sample_size, sample_freq, sample_burn);

  /* grid */
  setup_mcmc_grid (gridpoint_s_begin, gridpoint_s_end, gridpoint_s_number,
                   gridpoint_p_begin, gridpoint_p_end, gridpoint_p_number,
                   gridpoint_n_begin, gridpoint_n_end, gridpoint_n_number);

  /* prior */
  setup_mcmc_prior (flatprior_s_begin, flatprior_s_end,
                    flatprior_p_begin, flatprior_p_end,
                    flatprior_n_begin, flatprior_n_end);

  /* delta */
  setup_mcmc_delta (delta_s, delta_p, delta_n);

  /* init */
  setup_mcmc_init (init_s, init_p, init_a, init_c, init_g);

  /* data */ 
  setup_mcmc_data ();

  /* file */
  setup_mcmc_file (ofile, pname, aname);

  /* want */
  setup_mcmc_want (want_no_data, 
                   want_disable_s, want_disable_p, want_disable_n,
                   want_save_alldraws, want_debug_est, log_freq);
 
  setup_mcmc_limit (limit_tau, limit_psr);
  setup_mcmc_state (want_load_state, limit_time_to_run);
  setup_bf_state (want_load_state);

  /* these two functions should be gotten from mcmc.c module */

  /* want and file is okay but,
     sample_size and 0.9 what is this. 
     ofile should be replaced by oname and I will switch the positions
     of arguments: qname and aname. sample_size is a part of mcmc.c
     module. access_mcmc_sample_size () should be called in bf. c module.
   */
  setup_bf_want (want_load_param, want_save_alldraws);

  setup_bf_file (ofile, pname, aname, qname, bname);
  setup_bf_tol (tolerance_zeros);

  if (want_verbose == 1)
    {
      fprintf (stderr, "%s\n%s\n%s\n%s\n%s\n%s\n%s\n", datname, energyname, pname, aname, gname, qname, bname);
    }

  /* creation of grid points */
  create_gridpoints(number_chains);

  /* random number generator and energy data */
  if (want_load_state == 1) 
    load_gridpoints (gname);

  gmel_grid_print_grid ();

  return EXIT_SUCCESS;
}

static int
prepare_filenames ()
{
  int r = EXIT_SUCCESS;
  int r1 = EXIT_SUCCESS;
  int r2 = EXIT_SUCCESS;
  int len_data_directory;
  int len_output_directory;
  int len_data;
  char pdbid[5];
  char chainid;
  
  /* 1hi2A or 1hi2_ */
  if (auto_data_name != NULL && output_directory != NULL && data_directory != NULL)
    {
      len_data_directory = strlen (data_directory);
      len_output_directory = strlen (output_directory);
      len_data = strlen (auto_data_name);
      assert (len_data == 5);
      strncpy (pdbid, auto_data_name, 4);
      pdbid[4] = '\0'; /* just for making sure if the last is null terminated */
      chainid = auto_data_name[4];
      datname = XMALLOC (char, len_data_directory + 4 + 4 + 6 + 1);
      energyname = XMALLOC (char, len_data_directory + 11 + 1);
      oname = XMALLOC (char, len_output_directory + 4 + 4 + 4 + 1);
      rname = XMALLOC (char, len_output_directory + 4 + 4 + 4 + 1);
      pname = XMALLOC (char, len_output_directory + 4 + 4 + 4 + 1);
      aname = XMALLOC (char, len_output_directory + 4 + 4 + 4 + 1);
      gname = XMALLOC (char, len_output_directory + 4 + 4 + 4 + 1);
      qname = XMALLOC (char, len_output_directory + 4 + 4 + 4 + 1);
      bname = XMALLOC (char, len_output_directory + 4 + 4 + 4 + 1);
      sprintf (datname, "%s/pdb%s.%c.prd", data_directory, pdbid, chainid);
      sprintf (energyname, "%s/energy.int", data_directory); 
      sprintf (oname, "%s/pdb%s.%c.o", output_directory, pdbid, chainid);
      sprintf (rname, "%s/pdb%s.%c.r", output_directory, pdbid, chainid);
      sprintf (pname, "%s/pdb%s.%c.p", output_directory, pdbid, chainid);
      sprintf (aname, "%s/pdb%s.%c.a", output_directory, pdbid, chainid);
      sprintf (gname, "%s/pdb%s.%c.g", output_directory, pdbid, chainid);
      sprintf (qname, "%s/pdb%s.%c.q", output_directory, pdbid, chainid);
      sprintf (bname, "%s/pdb%s.%c.b", output_directory, pdbid, chainid);

/*
      fprintf (stderr, "%s\n%s\n%s\n%s\n%s\n%s\n%s\n", datname, energyname, pname, aname, gname, qname, bname);
*/
    }

  if (want_no_overwrite == 1)                 /* --no-overwrite */
    {
      if (want_exe_estimation == 1)
        {
          r = psi_io_check_est_p (pname, sample_size);
          if (r == EXIT_SUCCESS)
            return EXIT_FAILURE;
        }
      else if (want_exe_bayesfactor == 1)
        {
          r1 = psi_io_check_bf_b (bname);
          r2 = psi_io_check_bf_p (pname, sample_size);
          if (r1 == EXIT_SUCCESS && r2 == EXIT_SUCCESS)
            return EXIT_FAILURE;
        }
      else
        {
          psi_fatal ("no implementation");
        }
    }

  return EXIT_SUCCESS;
}
