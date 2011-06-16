/* energy.c -- energy caluculation
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

/* BUGGY:
    
   We do not need this neighbor information, which is based on
   a specific translation table. In the multiple sequence case,
   we need to have this function.
*/ 

/** @start 1 */
#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "common.h"
#include "error.h"
#include "rng.h"
#include "seq.h"
#include "gslwrap.h"
#include "ds.h"
#include "energy.h"

double energy_put_out (double d)
{
  return d;
}

#ifndef ACCESS_CATEGORY
#  define ACCESS_CATEGORY 5
#endif

#ifndef NUM_INTPAIR
#  define NUM_INTPAIR 5
#endif

#ifndef NUM_ALL_CORD_AA
#  define NUM_ALL_CORD_AA 15
#endif

#ifndef Dist
#  define Dist(x1,y1,z1,x2,y2,z2) \
          sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
#endif

/* goshng: NUM_CODON_EXCEPT_STOP should be other than 61 in other translation table */
static const int want_verbose = 1;
static const int want_debug = 0;
static const int NUM_AMINOACID = 20;
static const int NUM_NUCLEOTIDE = 4;
static const int NUM_CODON_EXCEPT_STOP = 61;
static const int NUM_COL_NEIGH = 27;
static const int NUM_CORD_NEIGH = 12;
static const double WITHIN_INTERACT_RANGE = 10.0;
static const int NUM_ROSETTA = 93; 
static const int NUM_EAA_MATRIX = 400; 
static PsiEnergyJones jones_measure;

extern const char NUCLEOTIDE[4];
extern const int INT_NUCLEOTIDE[128];
extern const char AMINOACID[21];
extern const int INT_AMINOACID[128];
extern int Nuc2AATable[4][4][4];
extern int Nuc2AATableStart[4][4][4];

/* simple accessors of jones_measure */
int*
dna_jones_measure ()
{
  /* return jones_measure.data.dna; */
  return jones_measure.psi.dna;
}

int*
pro_jones_measure ()
{
  /* return jones_measure.data.pro; */
  return jones_measure.psi.pro;
}

int
len_dna_jones_measure ()
{
  /* return jones_measure.data.len_dna; */
  return jones_measure.psi.len_dna;
}

int
len_pro_jones_measure ()
{
  /* return jones_measure.data.len_pro; */
  return jones_measure.psi.len_pro;
}

int 
init_codon_jones_measure ()
{
  return jones_measure.psi.init_codon;
}

void
data_members_jones_measure (int *n_a, int *n_c, int *n_g, int *n_t, 
                            double *solv, double *pair)
{
  *n_a = jones_measure.data.n_a;
  *n_c = jones_measure.data.n_c;
  *n_g = jones_measure.data.n_g;
  *n_t = jones_measure.data.n_t;
  *solv = jones_measure.data.solv;
  *pair = jones_measure.data.pair;
}

int
setup_energy_jones (const char *data_fn, const char *energy_fn)
{
  assert (data_fn != NULL);
  assert (energy_fn != NULL);
/*
  int len_data_fn = strlen (data_fn);
  int len_energy_fn = strlen (energy_fn);

  jones_measure.dat_filename = XMALLOC (char, len_data_fn + 1);
  strcpy (jones_measure.dat_filename, data_fn);
  jones_measure.energy_filename = XMALLOC (char, len_energy_fn + 1);
  strcpy (jones_measure.energy_filename, energy_fn);
*/

  jones_measure.dat_filename = xstrdup (data_fn);
  jones_measure.energy_filename = xstrdup (energy_fn);

  jones_measure.datfile = fopen (data_fn, "r");
  if (jones_measure.datfile == NULL)
    {
      psi_fatal ("no dat file for Jones' protein structure measure: %s", data_fn);
    }
  jones_measure.energyfile = fopen (energy_fn, "r");
  if (jones_measure.energyfile == NULL)
    {
      psi_fatal ("no energy file for Jones' protein structure measure: %s", energy_fn);
    }

  return EXIT_SUCCESS;
}

static int init_jones ();
static int init_vienna ();
static int init_vlmm ();
static int init_iedb ();
static int init_bastolla ();
static int fin_jones ();
static int fin_vienna ();
static int fin_vlmm ();
static int fin_iedb ();
static int fin_bastolla ();
static int score_jones_solv (int* const protein, double *solvNRG);
static int score_jones_pair (int* const protein, double *pairNRG);
static int score_jones (int* const protein, double *solvNRG, double *pairNRG);
/* This does not change the amino acid of protein at the site */
static void score_jones_solv_onlyone_aa (int *protein, double *solvNRG, 
                                         int site, int aa);
/* This does not change the amino acid of protein at the site */
static void score_jones_pair_onlyone_aa (int *protein, double *pairNRG, 
                                         int site, int aa);
static void score_jones_onlyone_aa (int *protein, 
                                    double *solvNRG, double *pairNRG, 
                                    int site, int aa);
static int score_vienna ();
static int score_vlmm();
static int score_iedb ();
static int score_bastolla ();

static void read_dat ();
static void create_interaction ();
static void SetUpInteraction ();


/*!\brief Setup solvent accessibility

   It reads in solvent accessibility information from DAT and energy files,
   and then assign Interaction->AAinfo[site_idx]->solvent[NUM_AMINOACID].
   For each amino acid site, there are twenty possible solvent accessibility
   scores that are used to calculate solvent accessbiltiy score, or 
   \f$ E_s(i) \f$.
 */
static void SetSolvent ();

/*!\brief Find a matrix for pairwise interaction score

   Taking into account separtion of 1-D along the sequence, find a right
   matrix among about 450 matrices.
 */
static int MatrixHunter (int seqsep, int pair, 
                         double **rosetta, double distance);

/*!\brief Setup codon's relation information

   It reads in 61*27 neighboring information which represents relations
   between a codon and its nine possible neighboring codons.
 */
static void SetNeighbor ();

/*!\brief Read pairwise info
   
    It reads in pairwise interaction information from DAT file.
 */
static void SetProt (double **prot);

/*!\brief Set distance 
 
 */
static void SetDDDneigh (double **prot, int *newOrder);
static void SetDDD (double **DDD, int *newOrder);

/*!\brief Create I's energyAcc member

   \param FirstOrder The only one energy strucutre
   \param AAlen length of protein

   It creates 3-D integer array of size AAlen*numNeigh*5. The function name
   is misleading and it should be something else, like CreateEnergyAcc 
   because this function does not set any values at all.
 */
static void SetEnergyAcc ();

/*!\brief Set magic stone

  It reads in a magic stone made by Doug which will be used to locate a proper
  matrix.
 */
static void SetRosetta (double **rosetta);

/*!\brief 

   It calculates five pairwise distances between atoms. For each pair,
   it assigns a proper(?) matrix index, which is one of 445 matrices.
 */
static void SetMatrixIndex (double **prot, double **rosetta, int *newOrder);

/*!\brief Reading in a huge matrices from energy file

   It reads in 445 matrices from energy file.
 */ 
static void SetNRGmats ();

/*!\brief Five Solvent Accessibility Categories

   For a given solvent accessibility in a DAT file, it spits out
   one of five numbers (0, 1, 2, 3, 4).
 */
static int SolventAccess (float acc);

static void 
read_dat () 
{
  /* copy of psi */
  Psi *psi = &jones_measure.psi;
  Data *dat = &(jones_measure.data);
  FILE *datfile = jones_measure.datfile;

  int i;
  char nucleotide;
  int int_nucleotide;
  int n_nucleotide[4] = { 0, 0, 0, 0 }; 
  char *line[80];
  int n_seq, len_dna, transl_table_id, init_codon;

  fscanf (datfile, "%d%d%d%d\n", 
          &n_seq, &len_dna, &transl_table_id, &init_codon);
  /* Translation table and init codon */
  choose_transl_table (transl_table_id, init_codon);
  psi->init_codon = init_codon;

  dat->n_seq = n_seq;
  dat->len_dna = len_dna;
  dat->transl_table_id = transl_table_id;
  dat->init_codon = init_codon;
  /* copy of psi */
  psi->n_seq = n_seq;
  psi->len_dna = len_dna;
  psi->transl_table_id = transl_table_id;
  psi->init_codon = init_codon;

  assert (dat->len_dna % 3 == 0);
  dat->len_pro = dat->len_dna / 3;
  dat->str_dna = XMALLOC (char, dat->len_dna + 1);
  dat->dna = XMALLOC (int, dat->len_dna);
  /* copy of psi */
  assert (psi->len_dna % 3 == 0);
  psi->len_pro = psi->len_dna / 3;
  psi->str_dna = XMALLOC (char, psi->len_dna + 1);
  psi->dna = XMALLOC (int, psi->len_dna);

  if (want_debug == 1) 
    {
      fprintf (stderr, "DNA: ");
    }
  for (i = 0; i < dat->len_dna; i++) {
      fscanf (datfile, "%c", &nucleotide);
      if (nucleotide == '\n') 
        {
          i--;
          continue; 
        }
      if (want_debug == 1) 
        {
          fprintf (stderr, "%c", nucleotide);
        }
      dat->str_dna[i] = nucleotide;
  /* copy of psi */
      psi->str_dna[i] = nucleotide;
      int_nucleotide = (int) nucleotide;
      dat->dna[i] = INT_NUCLEOTIDE[int_nucleotide];
  /* copy of psi */
      psi->dna[i] = INT_NUCLEOTIDE[int_nucleotide];
      n_nucleotide[INT_NUCLEOTIDE[int_nucleotide]]++;
  }
  if (want_debug == 1) 
    {
      fprintf (stderr, "\n");
    }
  dat->str_dna[i] = '\0';
  dat->n_a = n_nucleotide[PSI_DNA_A];
  dat->n_c = n_nucleotide[PSI_DNA_C];
  dat->n_g = n_nucleotide[PSI_DNA_G];
  dat->n_t = n_nucleotide[PSI_DNA_T];

  /* copy of psi */
  psi->str_dna[i] = '\0';
  psi->n_a = n_nucleotide[PSI_DNA_A];
  psi->n_c = n_nucleotide[PSI_DNA_C];
  psi->n_g = n_nucleotide[PSI_DNA_G];
  psi->n_t = n_nucleotide[PSI_DNA_T];
  /* (1); */
  /* fgets((char *)line, 80, datfile); */
  fscanf (datfile, "%s", &line); /* WARNING: check it out! */
}

static int 
init_jones ()
{
  int r; 
  read_dat ();
  create_interaction ();
  /* copy of psi: replace data to psi */
  int *pro = NULL;
  pro = XMALLOC (int, jones_measure.data.len_pro);
  r = dna2protein (jones_measure.data.dna, pro, jones_measure.data.len_pro);
  assert (r == EXIT_SUCCESS);

  score_jones (pro, &jones_measure.data.solv, &jones_measure.data.pair);
  if (want_debug == 1) 
    {
      fprintf (stderr, "%s\n", jones_measure.data.str_pro);
      fprintf (stderr, "solv: %lf, pair: %lf\n", jones_measure.data.solv,
                                                 jones_measure.data.pair);
    }
  
  XFREE (pro);
  return EXIT_SUCCESS;
}

static int 
init_vienna ()
{
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}

static int 
init_vlmm ()
{
  /* partially implemented in drprot, nrgdrpro and pst */
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}

static int 
init_iedb ()
{
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}

static int 
init_bastolla ()
{
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}

static int 
fin_jones ()
{
  Interaction *FirstOrder = &jones_measure.info;
  Data *dat = &jones_measure.data;
  int len = jones_measure.data.len_pro;

  /* goshng: no argument */
  delete_psi (&jones_measure.psi);
  delete_energy (&FirstOrder, len);
  delete_dat (&dat);
  XFREE (jones_measure.dat_filename);
  XFREE (jones_measure.energy_filename);
  fclose (jones_measure.datfile);
  fclose (jones_measure.energyfile);
  return EXIT_SUCCESS;
}

static int 
fin_vienna ()
{
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}

static int 
fin_vlmm ()
{
  /* partially implemented in drprot, nrgdrpro and pst */
  psi_fatal ("no implementaton");
/*  delete_dat( &Dat );
  free_AB (&s_pst_AB);
  free_pst (energyname, &s_pst_AB, &s_pst_absize, &s_pst_max_L, &s_pst_T); */
  return EXIT_SUCCESS;
}

static int 
fin_iedb ()
{
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}

static int 
fin_bastolla ()
{
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}


static int 
score_jones_solv (int* const protein, double *solvNRG)
{
  Psi *psi = &jones_measure.psi;
  Interaction *nrg= &jones_measure.info;
  int len = jones_measure.data.len_pro;

  int i_psi;
  double s_psi = 0.0;
  double s = 0.0;
  int aa1;
  AAsiteInfo *AAinfoPtr = nrg->AAinfo;
  int *initAA = protein;

  /* goshng: from here */
  for (aa1 = 0; aa1 < len; aa1++)
    {
      /* for a given amino acid site */
      s += AAinfoPtr->solvent[*initAA];
      
      i_psi = psi->access[aa1];
      s_psi += psi->solvent_data[i_psi][*initAA];

      if (s != s_psi)
        {
          fprintf (stderr, "s: %lf, s_psi: %lf, - %d, %lf %lf - %d\n", s, s_psi, aa1, AAinfoPtr->solvent[*initAA], psi->solvent_data[i_psi][*initAA], *initAA);
        }
      AAinfoPtr++;
      initAA++;
    }

  if (s != s_psi)
    {
      fprintf (stderr, "s: %lf, s_psi: %lf\n", s, s_psi);
    }
  assert (s == s_psi);
 
  assert (finite(s));

  /* NOTE: Only the energy associated with the ancestral sequence is kept */
  *solvNRG = -2.0 * s; /* Need this for updating w_solv */
 
  return EXIT_SUCCESS;

}

static int 
score_jones_pair (int* const protein, double *pairNRG)
{
  Psi *psi = &jones_measure.psi;
  Interaction *nrg= &jones_measure.info;
  int len = jones_measure.data.len_pro;

  double p_psi = 0.0;
  double p = 0.0;
  double **T2NRG = NULL;
  double *T1NRG = NULL; 
  int aa1, aa2;
  AAsiteInfo *AAinfoPtr = nrg->AAinfo;
  int *DDD = NULL;;
  int numNeigh;
  int *initAA = protein;

/* here */
  p_psi = inter_score (psi, protein);

  for (aa1 = 0; aa1 < len; aa1++)
    {
      /* for neighboring sites */ 
      T2NRG = AAinfoPtr->pairwise;
      numNeigh = AAinfoPtr->numNeighbors;
      DDD = AAinfoPtr->DDDneighbors;
      for (aa2 = 0; aa2 < numNeigh; aa2++)
        {
          /* 20-by-20 matrix T1NRG */
          T1NRG = *T2NRG;
          /* for all neighboring sites to the right of a current site */
          /* this if condition prevents duplicated summing of pairwise energy */
          if (*DDD > aa1)
            {
              p += T1NRG[((*initAA) * 20) + protein[*DDD]];
            }
          /* next neighbor */
          DDD++; 
          T2NRG++;
        }
      initAA++; 
      AAinfoPtr++;
    }

  if (gsl_fcmp(p, p_psi, 1e-6) != 0)
    {
      fprintf (stderr, "p: %e, p_psi: %e\n", p, p_psi);
    }
/*  assert (trunc(p*1000) == trunc(p_psi*1000)); */
  assert (finite(p));

  /* NOTE: Only the energy associated with the ancestral sequence is kept */
  *pairNRG = -2.0 * p; /* Need this for updating w_pair */
 
  return EXIT_SUCCESS;
}

static int 
score_jones (int* const protein, double *solvNRG, double *pairNRG)
{
  score_jones_solv (protein, solvNRG);
  score_jones_pair (protein, pairNRG);
  return EXIT_SUCCESS;

  /* NO CODE RUN */
  Interaction *nrg= &jones_measure.info;
  int len = jones_measure.data.len_pro;

  double s = 0.0;
  double p = 0.0;
  double **T2NRG = NULL;
  double *T1NRG = NULL; 
  int aa1, aa2;
  AAsiteInfo *AAinfoPtr = nrg->AAinfo;
  int *DDD = NULL;;
  int numNeigh;
  int *initAA = protein;

  if (want_debug == 1) 
    {
      fprintf (stderr, "EACH SOLV: ");
    }
  for (aa1 = 0; aa1 < len; aa1++)
    {
      /* for neighboring sites */ 
      T2NRG = AAinfoPtr->pairwise;
      numNeigh = AAinfoPtr->numNeighbors;
      DDD = AAinfoPtr->DDDneighbors;
      s += AAinfoPtr->solvent[*initAA];
      if (want_debug == 1) 
        {
          fprintf (stderr, "%lf ", s);
        }
      for (aa2 = 0; aa2 < numNeigh; aa2++)
        {
          /* 20-by-20 matrix T1NRG */
          T1NRG = *T2NRG;
          /* for all neighboring sites to the right of a current site */
          /* this if condition prevents duplicated summing of pairwise energy */
          if (*DDD > aa1)
            {
              p += T1NRG[((*initAA) * 20) + protein[*DDD]];
            }
          /* next neighbor */
          DDD++; T2NRG++;
        }
      initAA++; 
      AAinfoPtr++;
    }
  if (want_debug == 1) 
    {
      fprintf (stderr, "\n");
    }

  assert (finite(s));
  assert (finite(p));

  /* NOTE: Only the energy associated with the ancestral sequence is kept */
  *solvNRG = -2.0 * s; /* Need this for updating w_solv */
  *pairNRG = -2.0 * p; /* Need this for updating w_pair */
 
  return EXIT_SUCCESS;
}

/* This does not change the amino acid of protein at the site */
static void 
score_jones_solv_onlyone_aa (int *protein,
                              double *solvNRG, 
                              int site, int aa)
{
  Interaction *nrg= &jones_measure.info;

  double s = 0.0;
  int prev_aa = protein[site];

  AAsiteInfo *AAinfoPtr = &(nrg->AAinfo[site]);
  s = AAinfoPtr->solvent[aa] - AAinfoPtr->solvent[prev_aa];
      
  *solvNRG = (*solvNRG) - 2.0 * s; /* Need this for updating w_solv */
  assert (finite(*solvNRG));

  return;
}

/* This does not change the amino acid of protein at the site */
static void 
score_jones_pair_onlyone_aa (int *protein,
                             double *pairNRG, 
                             int site, int aa)
{
  Interaction *nrg= &jones_measure.info;

  double p = 0.0;
  double **T2NRG = NULL;
  double *T1NRG = NULL; 
  int aa2;
  int *DDD = NULL;;
  int numNeigh;
  int prev_aa = protein[site];
  

  AAsiteInfo *AAinfoPtr = &(nrg->AAinfo[site]);
      
  T2NRG = AAinfoPtr->pairwise;
  numNeigh = AAinfoPtr->numNeighbors;
  DDD = AAinfoPtr->DDDneighbors;
  for (aa2 = 0; aa2 < numNeigh; aa2++)
    {
      /* 20-by-20 matrix T1NRG */
      T1NRG = *T2NRG;
      /* for all neighboring sites to the right of a current site */
      /* this if condition prevents duplicated summing of pairwise energy */
      if (*DDD > site)
        {
          p += T1NRG[(aa * 20) + protein[*DDD]] 
               - T1NRG[(prev_aa * 20) + protein[*DDD]];
        } 
      else 
        {
          p += T1NRG[aa + protein[*DDD] * 20] 
               - T1NRG[prev_aa + protein[*DDD] * 20];
        }
        DDD++; T2NRG++;
    }
  *pairNRG = (*pairNRG) - 2.0 * p; /* Need this for updating w_pair */
  assert (finite(*pairNRG));

  return;
}

static void 
score_jones_onlyone_aa (int *protein,
                        double *solvNRG, double *pairNRG, 
                        int site, int aa)
{
  if (protein[site] == aa) 
    {
      return;
    }

  score_jones_solv_onlyone_aa (protein, solvNRG, site, aa);
  score_jones_pair_onlyone_aa (protein, pairNRG, site, aa);

  protein[site] = aa;
  return;
}

static int 
score_vienna ()
{
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}

static int 
score_vlmm ()
{
  /* partially implemented in drprot, nrgdrpro and pst */
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}

static int 
score_iedb ()
{
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}

static int 
score_bastolla ()
{
  psi_fatal ("no implementaton");
  return EXIT_SUCCESS;
}

int 
initialize_drevol (int e)
{
  switch (e)   
    {
      case PSI_ENERGY_JONES: 
        init_jones (); break;
      case PSI_ENERGY_VIENNA: 
        init_vienna (); break;
      case PSI_ENERGY_VLMM:
        init_vlmm (); break;
      case PSI_ENERGY_IEDB:
        init_iedb (); break;
      case PSI_ENERGY_BASTOLLA:
        init_bastolla (); break;
      default:
        assert (0);
    }
  return EXIT_SUCCESS;
}

int 
finalize_drevol (int e)
{
  switch (e)   
    {
      case PSI_ENERGY_JONES: 
        fin_jones (); break;
      case PSI_ENERGY_VIENNA: 
        fin_vienna (); break;
      case PSI_ENERGY_VLMM:
        fin_vlmm (); break;
      case PSI_ENERGY_IEDB:
        fin_iedb (); break;
      case PSI_ENERGY_BASTOLLA:
        fin_bastolla (); break;
      default:
        assert (0);
    }
  return EXIT_SUCCESS;
}

int 
score_drevol (int e, ...)
{
  va_list ap;
  
/*
  int protein[3] = {3, 4, 5};
  double solvNRG, pairNRG;
*/
  int *protein;
  double *solvNRG;
  double *pairNRG;

  va_start (ap, e);
  switch (e)
    {
      case PSI_ENERGY_JONES: 
        protein = va_arg (ap, int *);
        solvNRG = va_arg (ap, double *);
        pairNRG = va_arg (ap, double *);
        score_jones (protein, solvNRG, pairNRG);
        break;
      case PSI_ENERGY_VIENNA: 
        score_vienna (); break;
      case PSI_ENERGY_VLMM:
        score_vlmm (); break;
      case PSI_ENERGY_IEDB:
        score_iedb (); break;
      case PSI_ENERGY_BASTOLLA:
        score_bastolla (); break;
      default:
        assert (0);
    }

  va_end (ap);
  return EXIT_SUCCESS;
}

void 
calc_energy (int *protein, double *solvNRG, double *pairNRG)
{
  score_jones (protein, solvNRG, pairNRG);
}

void 
calc_energy_neigh (int *protein, 
                   double *solv, double *pair,
                   int site, int aa)
{
  score_jones_onlyone_aa (protein, solv, pair, site, aa);
}

int 
energy_is (int *protein, double *solv, double *pair )
{
  score_jones (protein, solv, pair);
  return EXIT_SUCCESS;
}

int 
energy_solv_is (int *protein, double *solv)
{
  score_jones_solv (protein, solv);
  return EXIT_SUCCESS;
}

int 
energy_res_solv_is (int site, int aa, double *solv)
{
  Psi *psi = &jones_measure.psi;
  Interaction *nrg = &jones_measure.info;
  double s;
  int i_psi;
  double s_psi = 0.0;

  AAsiteInfo *tempInfo = nrg->AAinfo;
  s = tempInfo[site].solvent[aa];

  i_psi = psi->access[site];
  s_psi = psi->solvent_data[i_psi][aa];

  if (s != s_psi)
    {
      fprintf (stderr, "s: %lf, s_psi: %lf\n", s, s_psi);
    }
  assert (s == s_psi);
/*
  AAsiteInfo *AAinfoPtr = &(nrg->AAinfo[site]);
  s = AAinfoPtr->solvent[aa];
*/
  *solv = -2.0 * s;

  return EXIT_SUCCESS;
}

int 
energy_res_pair_is (int *protein, int site, int aa, double *pair)
{
  Psi *psi = &jones_measure.psi;
  Interaction *nrg= &jones_measure.info;

  double p = 0.0;
  double **T2NRG = NULL;
  double *T1NRG = NULL; 
  int aa2;
  int *DDD = NULL;;
  int numNeigh;
  int initAA;

  double p_psi = 0.0;
  p_psi = inter_score_res (psi, protein, site, aa);

  AAsiteInfo *AAinfoPtr = &(nrg->AAinfo[site]);
  T2NRG = AAinfoPtr->pairwise;
  numNeigh = AAinfoPtr->numNeighbors;
  DDD = AAinfoPtr->DDDneighbors;
  for (aa2 = 0; aa2 < numNeigh; aa2++)
    {
      /* 20-by-20 matrix T1NRG */
      T1NRG = *T2NRG;
      /* for all neighboring sites to the right of a current site */
      /* this if condition prevents duplicated summing of pairwise energy */
      initAA = protein[*DDD];
      if (*DDD < site) 
        {
          p += T1NRG[initAA * 20 + aa];
/* fprintf (stderr, "T1: %lf\n", T1NRG[initAA * 20 + aa]); */
        } 
      else 
        {
          p += T1NRG[aa * 20 + initAA];
/* fprintf (stderr, "T1: %lf\n", T1NRG[aa * 20 + initAA]); */
        }
      DDD++; T2NRG++;
    }

/*
  if (trunc(p*1000) != trunc(p_psi*1000))
    {
      psi_warning ("p: %lf, p_psi: %lf - %s %d %s\n", p, p_psi,  __FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
*/
/*  assert (trunc(p*1000) == trunc(p_psi*1000)); */
  assert (finite(p));

  /* NOTE: Only the energy associated with the ancestral sequence is kept */
  *pair = -2.0 * p; /* Need this for updating w_pair */
 
  return EXIT_SUCCESS;

/*
                       E_pair = 0.0;
                       T2NRG = AAinfoPtr->pairwise;
                       FO_numNeighAAcol = AAinfoPtr->numNeighbors;
                       DDD = AAinfoPtr->DDDneighbors;
                       for (aa2 = 0; aa2 < FO_numNeighAAcol; aa2++)
                       {
                          T1NRG = *T2NRG;
                          initAA = tmpAAseq[*DDD];
                          if (*DDD < AAcol) {
                             E_pair += T1NRG[initAA * 20 + tempAA];
                          } else {
                             E_pair += T1NRG[tempAA * 20 + initAA];
                          }
                          DDD++; T2NRG++;
                       }
*/

}

int energy_res_is (int *protein, int site, int aa, double *solv, double *pair)
{
  energy_res_solv_is (site, aa, solv);
  energy_res_pair_is (protein, site, aa, pair);
  return EXIT_SUCCESS;
}

int 
energy_pair_is (int *protein, double *pair)
{
  score_jones_pair (protein, pair);
  return EXIT_SUCCESS;
}

int 
energy_near_is (int *protein, int site, int aa,
                double *solv, double *pair)
{
  score_jones_onlyone_aa (protein, solv, pair, site, aa);
  return EXIT_SUCCESS;
}

int 
energy_near_solv_is (int *protein, int site, int aa,
                     double *solv)
{
  if (protein[site] == aa) 
    {
      return EXIT_SUCCESS;
    }
  score_jones_solv_onlyone_aa (protein, solv, site, aa);
  protein[site] = aa;
  return EXIT_SUCCESS;
}

int 
energy_near_pair_is (int *protein, int site, int aa,
                     double *pair)
{
  if (protein[site] == aa) 
    {
      return EXIT_SUCCESS;
    }
  score_jones_pair_onlyone_aa (protein, pair, site, aa);
  protein[site] = aa;
  return EXIT_SUCCESS;
}

int 
check_translation ()
{
  Data *dat = &jones_measure.data;

  int i;
  int r = EXIT_SUCCESS;
  int codon[3];
  for (i = 0; i < dat->len_pro; i++) 
    {
      codon[PSI_CODON_FIRST] = dat->dna[i*3];
      codon[PSI_CODON_SECOND] = dat->dna[i*3 + 1];
      codon[PSI_CODON_THIRD] = dat->dna[i*3 + 2];
      if (aa_codon_is (codon) != dat->pro[i])
        {
          psi_fatal ("Translatoin does not match!");
          break;
        }
     }

  return r;
}

static void 
create_interaction ()
{
  /* from arguments of previous version */
  Interaction *I = &jones_measure.info;
  
  I->AAinfo = NULL;
  I->energyAcc = NULL;
  I->nAccess = NULL;
  I->neighList = NULL;
  I->energy = NULL;
  
  SetUpInteraction ();

  return;
} 

static void 
SetUpInteraction ()
{
  /* from arguments of previous version */
  Psi *psi = &jones_measure.psi;
  Interaction *FirstOrder = &jones_measure.info;
  int AAlen = jones_measure.data.len_pro;
  FILE *datfile = jones_measure.datfile;
  FILE *energyfile = jones_measure.energyfile;

  int length;
  int i, j, k, r, aa_i, aa_j, AAcol, aa2;

  assert (energyfile != NULL); 
  /* fscanf(datfile, "%d%d", &n_seq, &len_dna); */

  FirstOrder->AAinfo = XMALLOC (AAsiteInfo, AAlen); 
  SetSolvent ();
  SetNeighbor (); 

  fscanf (datfile, "%d\n", &length); /* reading a line */
  if (want_debug == 1)
    {
      fprintf(stderr, "%d %d\n", length, AAlen);
    }
  if (length != AAlen)
    {
      psi_fatal ("THIS FILE DOES NOT MATCH THE PROTEIN! lenth: %d, AAlen: %d", length, AAlen);
    }

  /* goshng: missing residue */
  /* int *newOrder = XMALLOC (int, AAlen); */
  int *newOrder = NULL;
  /* int *newOrder = new int[AAlen]; */
  double **prot = XMALLOC (double *, AAlen);
  /* double **prot = new double *[AAlen]; */
  for (i = 0; i < AAlen; i++)
    {
      prot[i] = XMALLOC(double, NUM_CORD_NEIGH);
      /* prot[i] = new double[NUM_CORD_NEIGH]; */
      /* newOrder[i] = 0; */
    }
  SetProt (prot);
  SetDDDneigh (prot, newOrder);
  /* We are done with SetDDDneigh */
  SetEnergyAcc ();
  /* We are done with SetEnergyAcc: only memory allocations */

  /* copy of psi: we need this rosetta  */
  double **rosetta = XMALLOC (double *, NUM_ROSETTA);
  /* double **rosetta = new double *[NUM_ROSETTA]; */
  for (r = 0; r < NUM_ROSETTA; r++)
  {
      rosetta[r] = XMALLOC (double, ACCESS_CATEGORY);
      /* rosetta[r] = new double[ACCESS_CATEGORY]; */
  }

  SetRosetta (rosetta);
  /* We are done with reading 93-5 matrix */

  /* copy of psi: energyAcc[AAlen][numNeigh][5] is set */
  SetMatrixIndex (prot, rosetta, newOrder);
  /* We are done with SetMatrixIndex */

  // We now delete the memory that stored rosetta    //
  for (r = 0; r < NUM_ROSETTA; r++)
  {
      XFREE (rosetta[r]);
      /* delete [] rosetta[r]; */
  }  //
  for (i = 0; i < AAlen; i++)
  {
      XFREE (prot[i]);
      /* delete [] prot[i]; */
  }  //
  XFREE (prot);
  /* delete[] prot; */
  XFREE (rosetta);
  /* delete[] rosetta; */
  /* XFREE (newOrder); */
  /* delete[] newOrder; */

  /* copy of psi: read all 445 matrices from energy.int file */
  SetNRGmats ();
  /*  We are done with SetNRGmats */

  int ***FO_energyAcc = FirstOrder->energyAcc;
  double **mPtr0, **mPtr1, **mPtr2, **mPtr3, **mPtr4, *tmpPtr;
  /* CHANGE by JEFF   double **mPtr0, **mPtr1, **mPtr2, **mPtr3, **mPtr4, **mPtr5, *tmpPtr; */
  int *access;
  /* CHANGE by JEFF   int colPtr, *access; */

  double ***NRG0 = FirstOrder->energy[0]; /* 20-by-20 matrices */
  double ***NRG1 = FirstOrder->energy[1]; /* 20-by-20 matrices */
  double ***NRG2 = FirstOrder->energy[2]; /* 20-by-20 matrices */
  double ***NRG3 = FirstOrder->energy[3]; /* 20-by-20 matrices */
  double ***NRG4 = FirstOrder->energy[4]; /* 20-by-20 matrices */

  /* copy of psi */
  int res1, res2;
  double matrix[NUM_ELEMENT][NUM_ELEMENT];

  for (AAcol = 0; AAcol < AAlen; AAcol++)
  {
      res1 = AAcol;
      int numNeigh = FirstOrder->AAinfo[AAcol].numNeighbors;
      FirstOrder->AAinfo[AAcol].pairwise = XMALLOC(double *, numNeigh);
      for (aa2 = 0; aa2 < numNeigh; aa2++)
      {
          /* copy of psi */
          /* WRONG: res2 = aa2; */
          res2 = FirstOrder->AAinfo[res1].DDDneighbors[aa2];

          FirstOrder->AAinfo[AAcol].pairwise[aa2] = XMALLOC(double, NUM_EAA_MATRIX);
          access = FO_energyAcc[AAcol][aa2];
/* fprintf(stderr, "%d %d: %d %d %d %d %d\n", AAcol, aa2, access[0],access[1],access[2],access[3],access[4]); */
          assert( access[0] != 150 );
          assert( access[1] != 150 );
          assert( access[2] != 150 );
          assert( access[3] != 150 );
          assert( access[4] != 150 );
          mPtr0 = NRG0[access[0]]; /* one of 93 matrices */
          mPtr1 = NRG1[access[1]]; /* one of 89 matrices */
          mPtr2 = NRG2[access[2]]; /* one of 93 matrices */
          mPtr3 = NRG3[access[3]]; /* one of 88 matrices */
          mPtr4 = NRG4[access[4]]; /* one of 82 matrices */
          tmpPtr = FirstOrder->AAinfo[AAcol].pairwise[aa2];
          int ptrCount = 0;
          for (aa_i = 0; aa_i < 20; aa_i++)
          {
              for (aa_j = 0; aa_j < 20; aa_j++)
              {
                  tmpPtr[ptrCount] = mPtr0[aa_i][aa_j]
                                   + mPtr1[aa_i][aa_j]
                                   + mPtr2[aa_i][aa_j]
                                   + mPtr3[aa_i][aa_j]
                                   + mPtr4[aa_i][aa_j];
                  /* copy of psi */
                  matrix[aa_i][aa_j] = mPtr0[aa_i][aa_j]
                                       + mPtr1[aa_i][aa_j]
                                       + mPtr2[aa_i][aa_j]
                                       + mPtr3[aa_i][aa_j]
                                       + mPtr4[aa_i][aa_j];
                  ptrCount++;
              }
          }
          /* copy of psi: That's it! This is all that about 
             pairwise interaction score setup 
          */
          inter_set (psi, res1, res2, matrix);
      }
  }
  int zero[NUM_INTPAIR];
  zero[0] = 93; zero[1] = 89; zero[2] = 93; zero[3] = 88; zero[4] = 82;
  for (i = 0; i < NUM_INTPAIR; i++)
  {
      for (j = 0; j < zero[i]; j++)
      {
          for (k = 0; k < NUM_AMINOACID; k++)
          {
              XFREE (FirstOrder->energy[i][j][k]);
              /* delete [] FirstOrder->energy[i][j][k]; */
          }
          XFREE (FirstOrder->energy[i][j]);
          /* delete [] FirstOrder->energy[i][j]; */
      }
      XFREE (FirstOrder->energy[i]);
      /* delete [] FirstOrder->energy[i]; */ 
  }
  XFREE (FirstOrder->energy);
  /* delete [] FirstOrder->energy; */
  FirstOrder->energy = NULL;
  /* We are done with Interaction */

  // Check # 2 (The Big One!)
#if LOUD
  // Here we do a MASSIVE Check on the system to see if it worked
  //   void CheckInteract(FirstOrder, AAlen, out);
  // END MASSIVE Check on the System!!!!!!
#endif

}  // End SetUpInteraction Subroutine

static void 
SetSolvent ()
{
  Psi *psi= &jones_measure.psi;
  Interaction *FO = &jones_measure.info;
  int AAlen = jones_measure.data.len_pro;
  FILE *datfile = jones_measure.datfile;
  FILE *energyfile = jones_measure.energyfile;

  ///////////////////////////////////////////////////////////////
  // Now we are going to set up the Solvent vector, which will //
  // tell us exatly how accessible each amino acid is to the   //
  // surrounding solvent. We also read in DJ's Solvent         //
  // Accessibilities for each AA as well.                      //
  ///////////////////////////////////////////////////////////////
  int i, j, k;
 
  int *solvAcc = XMALLOC (int, AAlen); 
  /* copy of psi */ 
  psi->access = XMALLOC (int, AAlen); 

  /* int *solvAcc = new int[AAlen]; // Initiallize ALWAYS! */
  // for(i = 0; i < AAlen; i++){FirstOrder->solvAcc[i] = 0;}
  // Like in the contact accessibility, we call this only once, but never
  // again.  Then the column number should be fixed from here on in.
  double *T1, *T2;
  float sol;
//  out << endl << "This is the solvent accessibility info " << endl;
  if (want_debug == 1) 
    {
      fprintf (stderr, "SOVLENT: ");
    }
  for (i = 0; i < AAlen; i++)
  {
      fscanf (datfile, "%f\n", &sol);
      if (want_debug == 1)
        {
          fprintf (stderr, "%lf ", sol);
        }
      solvAcc[i] = SolventAccess (sol);
  /* copy of psi */ 
  /* We choose one of five categories based on DSSP information.
     So this SolventAccess function is the key for the solvent 
     accessibility score.
   */
      psi->access[i] = SolventAccess (sol); 
  }
  if (want_debug == 1) 
    {
      fprintf (stderr, "\n");
    }
  double **solv = XMALLOC (double *, ACCESS_CATEGORY);
  /* double **solv = new double *[ACCESS_CATEGORY]; */
  double solvent[ACCESS_CATEGORY];
  for (i = 0; i < ACCESS_CATEGORY; i++)
  {
      solv[i] = XMALLOC (double, NUM_AMINOACID);
      /* solv[i] = new double[20]; */
  }
  for (i = 0; i < NUM_AMINOACID; i++)
  {
      fscanf (energyfile, "%lf%lf%lf%lf%lf", 
              &solvent[0], &solvent[1], &solvent[2], &solvent[3], &solvent[4]);
     /*   &solv[0][i], &solv[1][i], &solv[2][i], &solv[3][i], &solv[4][i] ); */
      for (j = 0; j < ACCESS_CATEGORY; j++) {
         solv[j][i] = solvent[j];
  /* copy of psi */ 
  /* That's it. We do not duplicate the solvent accessibility information */
         psi->solvent_data[j][i] = solvent[j];
      }
  }
  for (i = 0; i < AAlen; i++)
  {
      FO->AAinfo[i].solvent = XMALLOC (double, NUM_AMINOACID);
      /* FO->AAinfo[i].solvent = new double[NUM_AMINOACID]; */
      T1 = FO->AAinfo[i].solvent;
      T2 = solv[solvAcc[i]];
      for (j = 0; j < NUM_AMINOACID; j++)
      {
          *T1 = *T2; T1++; T2++;
      }
  }

  /* solvAcc: */
  XFREE (solvAcc);
  /* delete[] solvAcc; */
  for (k = 0; k < ACCESS_CATEGORY; k++)
  {
      XFREE (solv[k]);
      /* delete[] solv[k]; */
  }
  XFREE (solv);
  /* delete[] solv; */
}  // End SetSolvent Subroutine

static int 
SolventAccess(float acc)
{
  if (acc < 12.0)
  {
      return 0;
  }
  else if (acc < 36.0)
  {
      return 1;
  }
  else if (acc < 44.0)
  {
      return 2;
  }
  else if (acc < 87.0)
  {
      return 3;
  }
  else
  {
      return 4;
  }
} // End SolventAccess Subroutine

/* We do not need this neighbor information, which is based on
   a specific translation table. In the multiple sequence case,
   we need to have this function.
*/ 
static void 
SetNeighbor ()
{
  Interaction *FirstOrder = &jones_measure.info;
  FILE *energyfile = jones_measure.energyfile;

  /*******************************************************
     Now we are going to set up the Nearest Neighbor    
     information.  This includes setting up an index    
     matrix, as well as entering the nearest neighbors.
     The inNeigh file contains the AA neighbors, what 
     nucleotide caused it, and if it was a TS or TV. 
     These are constant data files that do not change
     even when the protein changes                  
  *******************************************************/
  int i;

  FirstOrder->nAccess = XMALLOC (int *, NUM_AMINOACID);
  /* FirstOrder->nAccess = new int * [NUM_AMINOACID]; */
  for (i = 0; i < NUM_AMINOACID; i++)
  {
      FirstOrder->nAccess[i] = XMALLOC(int, NUM_NUCLEOTIDE);
      /* FirstOrder->nAccess[i] = new int[NUM_NUCLEOTIDE]; */
      fscanf (energyfile, "%d%d%d%d", 
              &(FirstOrder->nAccess[i][0]), 
              &(FirstOrder->nAccess[i][1]), 
              &(FirstOrder->nAccess[i][2]), 
              &(FirstOrder->nAccess[i][3]) );
      ////////////////////////////////////////////
      //for (j = 0; j < NUM_NUCLEOTIDE; j++)//
      //{                                       //
      //   g_inter >> FirstOrder->nAccess[i][j];//
      //}                                       //
      ////////////////////////////////////////////
  }

  /* goshng: NUM_CODON_EXCEPT_STOP should be other than 61 in other translation table */
  FirstOrder->neighList = XMALLOC(int *, NUM_CODON_EXCEPT_STOP);
  /* FirstOrder->neighList = new int * [NUM_CODON_EXCEPT_STOP]; */
  for (i = 0; i < NUM_CODON_EXCEPT_STOP; i++)
  {
      FirstOrder->neighList[i] = XMALLOC(int, NUM_COL_NEIGH);
      /* FirstOrder->neighList[i] = new int[NUM_COL_NEIGH]; */
      fscanf (energyfile, 
              "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d", 
              &(FirstOrder->neighList[i][0]),  /* neighboring AA */
              &(FirstOrder->neighList[i][1]),  /* neighboring AA */
              &(FirstOrder->neighList[i][2]),  /* neighboring AA */
              &(FirstOrder->neighList[i][3]),  /* neighboring AA */
              &(FirstOrder->neighList[i][4]),  /* neighboring AA */
              &(FirstOrder->neighList[i][5]),  /* neighboring AA */
              &(FirstOrder->neighList[i][6]),  /* neighboring AA */
              &(FirstOrder->neighList[i][7]),  /* neighboring AA */
              &(FirstOrder->neighList[i][8]),  /* neighboring AA */
              &(FirstOrder->neighList[i][9]),  /* */
              &(FirstOrder->neighList[i][10]), /* */
              &(FirstOrder->neighList[i][11]), /* */
              &(FirstOrder->neighList[i][12]), /* */
              &(FirstOrder->neighList[i][13]), /* */
              &(FirstOrder->neighList[i][14]), /* */
              &(FirstOrder->neighList[i][15]), /* */
              &(FirstOrder->neighList[i][16]), /* */
              &(FirstOrder->neighList[i][17]), /* */
              &(FirstOrder->neighList[i][18]), /* TV or TS */
              &(FirstOrder->neighList[i][19]), /* TV or TS */
              &(FirstOrder->neighList[i][20]), /* TV or TS */
              &(FirstOrder->neighList[i][21]), /* TV or TS */
              &(FirstOrder->neighList[i][22]), /* TV or TS */
              &(FirstOrder->neighList[i][23]), /* TV or TS */
              &(FirstOrder->neighList[i][24]), /* TV or TS */
              &(FirstOrder->neighList[i][25]), /* TV or TS */
              &(FirstOrder->neighList[i][26]));/* TV or TS */

      ///////////////////////////////////////////////
      //for (j = 0; j < NUM_COL_NEIGH; j++)    //
      //{                                          //
      //    g_inter >> FirstOrder->neighList[i][j];//
      //}                                          //
      ///////////////////////////////////////////////
  }
}  // End SetNeighbor Subroutine

static void 
SetProt (double **prot)
{
  Psi *psi = &jones_measure.psi;
  Data *dat = &(jones_measure.data);
  FILE *datfile = jones_measure.datfile;

  int residue, i;
  /********************************************************
   Here we read in the coordinates from the data file.
   Note: we only read in N, CA, O, CB atom coordinates 
   initialize the prot structure                       
  ********************************************************/
  char aaRes; /* This holds the initial AA tag in the first column of file */
  int res_num;
  int count;
  int n_read;
  double pairwise[NUM_ALL_CORD_AA];
  assert (dat->len_pro > 0);
  dat->str_pro = XMALLOC (char, dat->len_pro + 1);
  dat->pro = XMALLOC (int, dat->len_pro);
  /* copy of psi */
  psi->str_pro = XMALLOC (char, psi->len_pro + 1);
  psi->pro = XMALLOC (int, psi->len_pro);
  psi->res_num = XMALLOC (int, psi->len_pro);
  for (residue = 0; residue < dat->len_pro; residue++)
    {
              /* "%c%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n", &aaRes, INPUT */
/*
      n_read = fscanf (datfile, 
              "%d%c%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n", &res_num, &aaRes,
              &pairwise[0], &pairwise[1], &pairwise[2],
              &pairwise[3], &pairwise[4], &pairwise[5],
              &pairwise[6], &pairwise[7], &pairwise[8],
              &pairwise[9], &pairwise[10], &pairwise[11],
              &pairwise[12], &pairwise[13], &pairwise[14]);
      assert (n_read == 17);
*/
       n_read = fscanf (datfile, 
              "%d ", &res_num);
      assert (n_read == 1);

       n_read = fscanf (datfile, 
              "%c%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n", &aaRes,
              &pairwise[0], &pairwise[1], &pairwise[2],
              &pairwise[3], &pairwise[4], &pairwise[5],
              &pairwise[6], &pairwise[7], &pairwise[8],
              &pairwise[9], &pairwise[10], &pairwise[11],
              &pairwise[12], &pairwise[13], &pairwise[14]);
      assert (n_read == 16);
 
      
  /*    while ((aaRes = fgetc (datfile)) != '\n') {}; */
      psi->res_num[residue] = res_num;
      count = 0;
      dat->pro[residue] = INT_AMINOACID[(int) aaRes];
      dat->str_pro[residue] = aaRes;
  /* copy of psi */
      psi->pro[residue] = INT_AMINOACID[(int) aaRes];
      psi->str_pro[residue] = aaRes;
      for (i = 0; i < NUM_ALL_CORD_AA; i++) 
        {
          /* prot[residue][i] = (double) pairwise[i]; */
          if ((i < 6) || (i > 8)) 
            {
              prot[residue][count] = pairwise[i];
              count++;
            }
        }
    }
  dat->str_pro[residue] = '\0';
  /* copy of psi */
  psi->str_pro[residue] = '\0';
}  // End SetProt Subroutine

void 
SetDDD (double **DDD, int *newOrder)
{
  Interaction *FO = &jones_measure.info;
  int AAlen = jones_measure.data.len_pro;
  if (newOrder != NULL)
    psi_fatal ("newOrder is replace by residue number");
  
  int aa1, aa2, k;
  /*******************************************************************
   Only accessible from SetDDDneigh using the matrix of CA atom  
   distances in DDD.  If the distance is less than 10.0, then it  
   is considered a neighbor in 3 dimensions.  Only these residues 
   will be considered in the subsequent energy calculations.      
  *******************************************************************/
  int *tmpNeighbor = XMALLOC (int, AAlen);
  /* int *tmpNeighbor = new int[AAlen]; */
  int count;
  int avg_count = 0;
  double distance;
  for (aa1 = 0; aa1 < AAlen; aa1++)
  {
      count = 0;
      for (aa2 = 0; aa2 < AAlen; aa2++)
      {
          if (aa1 != aa2)
          {
              if (aa1 < aa2)
              {
                  distance = DDD[aa1][aa2];
              }
              else
              {
                  distance = DDD[aa2][aa1];
              }
              if (distance < WITHIN_INTERACT_RANGE)
              {
                  tmpNeighbor[count] = aa2;
                  count++;
              }
          }
      }
      FO->AAinfo[aa1].numNeighbors = count;
      avg_count += count;
      FO->AAinfo[aa1].DDDneighbors = XMALLOC(int, count);
      for (k = 0; k < count; k++)
      {
          FO->AAinfo[aa1].DDDneighbors[k] = tmpNeighbor[k];
      }
  }
  avg_count /= AAlen;
  /* fprintf(stderr, "%d %d\n", avg_count, AAlen); */
  XFREE (tmpNeighbor);
}  // End SetDDD Subroutine

static void 
SetDDDneigh (double **prot, int *newOrder)
{
  int AAlen = jones_measure.data.len_pro;
  ////////////////////////////////////////////////
  // Now we are going to set up the Distance    //
  // matrix, in order to see how close the AA's //
  // are in 3-D space.                          //
  ////////////////////////////////////////////////
  int i, res_i, res_j;

  double **DDD = XMALLOC (double *, AAlen);
  /* double **DDD = new double *[AAlen]; */
  for (i = 0; i < AAlen; i++)
  {
      DDD[i] = XMALLOC(double, AAlen);
      /* DDD[i] = new double[AAlen]; */
  }

  /* int *missingSites; */
  int numMiss; 
  
  ////////////////////
  //g_in >> numMiss;//
  ////////////////////
  numMiss = 0;
  /* I'm not sure if I need to deal with this missing residue */
/*****************************************************************************
//  o << "This is the number of missing sites " << numMiss << endl;
  if (numMiss > 0) // If there are any missing sites in the protein, then this
  {
      // handles them by putting in a space there!!!
      missingSites = new int[numMiss];
      for (i = 0; i < numMiss; i++)
      {
          g_in >> missingSites[i];
//          o << "Those sites in order are (" << missingSites[i] << ")" << endl;
      }
      int count = 0;
      for (i = 0; i < AAlen; i++)
      {
          if (count == numMiss)
          {
              newOrder[i] = count;
          }
          else
          {
              if (i < missingSites[count])
              {
                  newOrder[i] = count;
              }
              else
              {
                  if (i == missingSites[count])
                  {
                      count++;
                      if (count < numMiss)
                      {
                          missingSites[count] -= count;
                      }
                  }
                  i--;
              }
          }
      }
      XFREE (missingSites);
      delete[] missingSites;
  }
*****************************************************************************/

  /***************************************************************************
   Initialize the distance matrix

   The critical part of this program is that we find the distance between 
   pairs of CB atoms only.  This set will be the only atoms that affect  
   the rate away from a particular atom.  The problem is that other atoms 
   could be within 10 A, but those will not be considered in future calcs.
  ***************************************************************************/
  int minusOne, resIplusOne; 
  minusOne = AAlen - 1;
  for (res_i = 0; res_i < minusOne; res_i++)
  {
      resIplusOne = res_i + 1;
      for (res_j = resIplusOne; res_j < AAlen; res_j++)
      {
          if (res_i < res_j)
          {
              /**************************************************************
                For this case we use the C_A atoms the center of the 10 A ball
                         DDD[res_i][res_j]=Dist(prot[res_i][3],prot[res_i][4],prot[res_i][5],
                              prot[res_j][3],prot[res_j][4],prot[res_j][5]);
                For this case we use the C_B atoms the center of the 10 A ball
              **************************************************************/
              DDD[res_i][res_j] = Dist (prot[res_i][9], prot[res_i][10],
                                        prot[res_i][11], prot[res_j][9],
                                        prot[res_j][10], prot[res_j][11]);
          }
      }
  }
  SetDDD (DDD, newOrder);
  // We now give back all of the memory that we took up until now ////////////
  for (i = 0; i < AAlen; i++)
  {
      XFREE (DDD[i]);
      /* delete[] DDD[i]; */
  }
  XFREE (DDD);
  /* delete[] DDD; */
}  // End SetDDDneigh Subroutine

static void 
SetEnergyAcc ()
{
  Interaction *FirstOrder = &jones_measure.info;
  int AAlen = jones_measure.data.len_pro;

  int w, x;
  // This structure will hold NRG matrix indices that are being used
  FirstOrder->energyAcc = XMALLOC (int **, AAlen);
  /* FirstOrder->energyAcc = new int **[AAlen]; */
  int n_neigh;
  for (w = 0; w < AAlen; w++)
  {
      n_neigh = FirstOrder->AAinfo[w].numNeighbors;
      FirstOrder->energyAcc[w] = XMALLOC (int *, n_neigh);
      /* FirstOrder->energyAcc[w] = new int * [t]; */
      for (x = 0; x < n_neigh; x++)
      {
          FirstOrder->energyAcc[w][x] = XMALLOC (int, ACCESS_CATEGORY);
          /* FirstOrder->energyAcc[w][x] = new int[ACCESS_CATEGORY]; */
      }
  }
}  // End SetEnergyAcc Subroutine

static void 
SetMatrixIndex (double **prot, double **rosetta, int *newOrder)
{
  Interaction *FirstOrder = &jones_measure.info;
  Psi *psi = &jones_measure.psi;
  int AAlen = jones_measure.data.len_pro;
  if (newOrder != NULL)
    psi_fatal ("newOrder is replace by residue number");

  int res, neigh, i;
  int targetRes, tempMat[NUM_INTPAIR], seqsep, first, next;
  double dddpair[NUM_INTPAIR]; // 0:CBCB, 1:CBN, 2:CB0, 3:NCB, 4:OCB
  /* These are the placements of all the zero matrices (Out of Range)
     zero[0] = 92; zero[1] = 88; zero[2] = 92; zero[3] = 87; zero[4] = 81; 
   */
  int ***FO_energyAcc, **FO_energyAccRes;
  /* CHANGE by JEFF   int *FO_numNeigh, *DDDneighbor, ***FO_energyAcc, **FO_energyAccRes; */
  FO_energyAcc = FirstOrder->energyAcc; /* int [AAlen][numNeigh][5] */
  int FO_numNeigh_res;
  for (res = 0; res < AAlen; res++)
  {
      FO_numNeigh_res = FirstOrder->AAinfo[res].numNeighbors;
      FO_energyAccRes = FO_energyAcc[res];
      //  DDDneighbor = FirstOrder->AAinfo[res].DDDneighbors;
      for (neigh = 0; neigh < FO_numNeigh_res; neigh++)
      {
          /* select one of neighboring residues */
          targetRes = FirstOrder->AAinfo[res].DDDneighbors[neigh];
          /*   targetRes = *DDDneighbor; */
          /* one resiue is preceding its neighboring one */
          if (res < targetRes)
          {
              /* how much are they separated in 1-D sequence */
              /* 
              seqsep = (targetRes - res)
                     + (newOrder[targetRes] - newOrder[res]);
              */
              seqsep = psi->res_num[targetRes] - psi->res_num[res];
              first = res; next = targetRes; 
          }
          else
          {
              /*
              seqsep = (res - targetRes)
                     + (newOrder[res] - newOrder[targetRes]);
              */
              seqsep = psi->res_num[res] - psi->res_num[targetRes];
              first = targetRes; next = res;
          }

          /* The following five pairwise distances are of interest to us */
          // CB x CB
          dddpair[0] = Dist(prot[first][9], prot[first][10], prot[first][11],
                            prot[next][9], prot[next][10], prot[next][11]);
          // CB x N
          dddpair[1] = Dist(prot[first][9], prot[first][10], prot[first][11],
                            prot[next][0], prot[next][1], prot[next][2]);
          // CB x O
          dddpair[2] = Dist(prot[first][9], prot[first][10], prot[first][11],
                            prot[next][6], prot[next][7], prot[next][8]);
          // N x CB
          dddpair[3] = Dist(prot[first][0], prot[first][1], prot[first][2],
                            prot[next][9], prot[next][10], prot[next][11]);
          // O x CB
          dddpair[4] = Dist(prot[first][6], prot[first][7], prot[first][8],
                            prot[next][9], prot[next][10], prot[next][11]);

          for (i = 0; i < NUM_INTPAIR; i++)
          {
              if (dddpair[i] < WITHIN_INTERACT_RANGE)
              {
                  tempMat[i] = MatrixHunter (seqsep, i, rosetta, dddpair[i]);
              }
              else
              {
                  tempMat[i] = MatrixHunter (seqsep, i, rosetta, 9.9999);
              }
/* goshng
   fprintf(ofile, "%d %d %e %d\n", seqsep, i, dddpair[i], tempMat[i]);
              else{tempMat[i] = zero[i];}
 */
              FO_energyAccRes[neigh][i] = tempMat[i];
              /* cout << tempMat[i] << " "; */
          }
          //    DDDneighbor++;
          //         cout << endl;
      } // End running through all of the AA's close in 3 dimensions
  } // End running through the residues in the protein
}// End Set MatrixIndex Subroutine

static void 
SetRosetta (double **rosetta)
{
  FILE *energyfile = jones_measure.energyfile;
  int r;
  // Set up the Rosetta matrix of break points
  for (r = 0; r < NUM_ROSETTA; r++)
  {
      fscanf (energyfile, "%lf%lf%lf%lf%lf", 
              &rosetta[r][0], &rosetta[r][1], &rosetta[r][2], 
              &rosetta[r][3], &rosetta[r][4]);
      /////////////////////////////////////////////
      //for (s = 0; s < ACCESS_CATEGORY; s++)//
      //{                                        //
          //g_inter >> rosetta[r][s];            //
      //}                                        //
      /////////////////////////////////////////////
  }
}  // End SetRosetta Subroutine

static void 
SetNRGmats ()
{
  Interaction *FirstOrder = &jones_measure.info;
  FILE *energyfile = jones_measure.energyfile;

  int i, j, k, l;
  /*******************************************************
   * Now we are going to set up the HUGE ENERGY Matrix. 
   * This MASSIVE structure includes 450 20x20 amino   
   * acid matrices of doubles.  This is for multiple  
   * distances for close, middle and distant         
   * interactions.                                  
  *******************************************************/
  int zero[NUM_INTPAIR];
  zero[0] = 93; zero[1] = 89; zero[2] = 93; zero[3] = 88; zero[4] = 82;
  // THIS IS IT!!!  We now set up the HUGE energy matrix

  FirstOrder->energy = XMALLOC(double ***, NUM_INTPAIR);
  /* FirstOrder->energy = new double ***[NUM_INTPAIR]; */
  double pairwise[NUM_AMINOACID]; 
  for (i = 0; i < NUM_INTPAIR; i++)
  {
      FirstOrder->energy[i] = XMALLOC(double **, zero[i]);
      /* FirstOrder->energy[i] = new double **[zero[i]]; */
      for (j = 0; j < zero[i]; j++)
      {
          FirstOrder->energy[i][j] = XMALLOC(double *, NUM_AMINOACID);
          /* FirstOrder->energy[i][j] = new double *[NUM_AMINOACID]; */
          for (k = 0; k < NUM_AMINOACID; k++)
          {
              FirstOrder->energy[i][j][k] = XMALLOC(double, NUM_AMINOACID);
              fscanf (energyfile, 
                      "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
                      &pairwise[0],  &pairwise[1],  &pairwise[2],  &pairwise[3], 
                      &pairwise[4],  &pairwise[5],  &pairwise[6],  &pairwise[7], 
                      &pairwise[8],  &pairwise[9],  &pairwise[10], &pairwise[11], 
                      &pairwise[12], &pairwise[13], &pairwise[14], &pairwise[15], 
                      &pairwise[16], &pairwise[17], &pairwise[18], &pairwise[19]);

              for (l = 0; l < NUM_AMINOACID; l++) {
                 FirstOrder->energy[i][j][k][l] = pairwise[l];
              }
          }
      }
  }
} // End SetNRGmats Subroutine

/*
void CheckInteract(Interaction *FirstOrder, int AAlen, ofstream &out)
{
  int i;
  // This checks to see if we have copied many of the input
  // files properly.  If they are wrong then this is NOT GOOD!
  out << endl << "Here is the Solvent accessibility information " << endl;
  // Check on the system

  out << endl << "Here is the Neighbor Access Matrix " << endl;
  for (i = 0; i < 20; i++)
  {
      for (j = 0; j < 4; j++)
      {
          out << FirstOrder->nAccess[i][j] << " ";
      }
      out << endl;
  }

  out << endl << "Here is the 9 nearest neighbor list" << endl;
  for (i = 0; i < 61; i++)
  {
      for (j = 0; j < 27; j++)
      {
          out << FirstOrder->neighList[i][j] << " ";
      }
      out << endl;
  }
}  // End CheckInteract Subroutine
*/

static int 
MatrixHunter (int seqsep, int pair, double **rosetta, double distance)
{
  int zero, preMats, mat, i;
  if (seqsep == 1)
  {
      if (pair == 0)
      {
          mat = 9;  zero = 92;
      }
      else if (pair == 1)
      {
          return 88;
      }
      else if (pair == 2)
      {
          mat = 12; zero = 92;
      }
      else if (pair == 3)
      {
          mat = 7;  zero = 87;
      }
      else
      {
          return 81;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[i][pair])
          {
              return i;
          }
      }
      return zero; // If it is not one of the above, then it = 0.0 matrix
  }
  else if (seqsep == 2)
  {
      if (pair == 0)
      {
          mat = 17; preMats = 9;  zero = 92;
      }
      else if (pair == 1)
      {
          mat = 13; preMats = 0;  zero = 88;
      }
      else if (pair == 2)
      {
          mat = 14; preMats = 12; zero = 92;
      }
      else if (pair == 3)
      {
          mat = 16; preMats = 7;  zero = 87;
      }
      else
      {
          mat = 12; preMats = 0;  zero = 81;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
      return zero; // If it is not one of the above, then it = 0.0 matrix
  }
  else if (seqsep == 3)
  {
      if (pair == 0)
      {
          mat = 11; preMats = 26;
      }
      else if (pair == 1)
      {
          mat = 14; preMats = 13;
      }
      else if (pair == 2)
      {
          mat = 10; preMats = 26;
      }
      else if (pair == 3)
      {
          mat = 10; preMats = 23;
      }
      else
      {
          mat = 11; preMats = 12;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  else if (seqsep == 4)
  {
      if (pair == 0)
      {
          mat = 8;  preMats = 37;
      }
      else if (pair == 1)
      {
          mat = 11; preMats = 27;
      }
      else if (pair == 2)
      {
          mat = 8;  preMats = 36;
      }
      else if (pair == 3)
      {
          mat = 8;  preMats = 33;
      }
      else
      {
          mat = 9;  preMats = 23;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  else if (seqsep == 5)
  {
      if (pair == 0)
      {
          mat = 7; preMats = 45;
      }
      else if (pair == 1)
      {
          mat = 8; preMats = 38;
      }
      else if (pair == 2)
      {
          mat = 6; preMats = 44;
      }
      else if (pair == 3)
      {
          mat = 7; preMats = 41;
      }
      else
      {
          mat = 7; preMats = 32;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  else if (seqsep == 6)
  {
      if (pair == 0)
      {
          mat = 6; preMats = 52;
      }
      else if (pair == 1)
      {
          mat = 7; preMats = 46;
      }
      else if (pair == 2)
      {
          mat = 6; preMats = 50;
      }
      else if (pair == 3)
      {
          mat = 6; preMats = 48;
      }
      else
      {
          mat = 6; preMats = 39;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  else if (seqsep == 7)
  {
      if (pair == 0)
      {
          mat = 6; preMats = 58;
      }
      else if (pair == 1)
      {
          mat = 6; preMats = 53;
      }
      else if (pair == 2)
      {
          mat = 6; preMats = 56;
      }
      else if (pair == 3)
      {
          mat = 5; preMats = 54;
      }
      else
      {
          mat = 6; preMats = 45;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  else if (seqsep == 8)
  {
      if (pair == 0)
      {
          mat = 5; preMats = 64;
      }
      else if (pair == 1)
      {
          mat = 5; preMats = 59;
      }
      else if (pair == 2)
      {
          mat = 5; preMats = 62;
      }
      else if (pair == 3)
      {
          mat = 4; preMats = 59;
      }
      else
      {
          mat = 5; preMats = 51;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  else if (seqsep == 9)
  {
      if (pair == 0)
      {
          mat = 4; preMats = 69;
      }
      else if (pair == 1)
      {
          mat = 4; preMats = 64;
      }
      else if (pair == 2)
      {
          mat = 5; preMats = 67;
      }
      else if (pair == 3)
      {
          mat = 4; preMats = 63;
      }
      else
      {
          mat = 5; preMats = 56;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  else if (seqsep == 10)
  {
      if (pair == 0)
      {
          mat = 4; preMats = 73;
      }
      else if (pair == 1)
      {
          mat = 4; preMats = 68;
      }
      else if (pair == 2)
      {
          mat = 4; preMats = 72;
      }
      else if (pair == 3)
      {
          mat = 4; preMats = 67;
      }
      else
      {
          mat = 4; preMats = 61;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  else if (seqsep == 11)
  {
      if (pair == 0)
      {
          mat = 3; preMats = 77;
      }
      else if (pair == 1)
      {
          mat = 4; preMats = 72;
      }
      else if (pair == 2)
      {
          mat = 4; preMats = 76;
      }
      else if (pair == 3)
      {
          mat = 4; preMats = 71;
      }
      else
      {
          mat = 4; preMats = 65;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  else if (seqsep < 23)
  {
      if (pair == 0)
      {
          mat = 6; preMats = 80;
      }
      else if (pair == 1)
      {
          mat = 6; preMats = 76;
      }
      else if (pair == 2)
      {
          mat = 6; preMats = 80;
      }
      else if (pair == 3)
      {
          mat = 6; preMats = 75;
      }
      else
      {
          mat = 6; preMats = 69;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  else
  {
      if (pair == 0)
      {
          mat = 6; preMats = 86;
      }
      else if (pair == 1)
      {
          mat = 6; preMats = 82;
      }
      else if (pair == 2)
      {
          mat = 6; preMats = 86;
      }
      else if (pair == 3)
      {
          mat = 6; preMats = 81;
      }
      else
      {
          mat = 6; preMats = 75;
      }

      for (i = 0; i < mat; i++)
      {
          if (distance < rosetta[preMats + i][pair])
          {
              return (preMats + i);
          }
      }
  }
  return 150;
}  // End MatrixHunter Subroutine

void 
delete_psi (Psi *psi) 
{
  XFREE (psi->str_dna);
  XFREE (psi->dna);
  XFREE (psi->str_pro);
  XFREE (psi->pro);
  XFREE (psi->access);
  XFREE (psi->res_num);
  inter_remove (&(psi->inter));
/*
  XFREE (*d);
  (*d) = NULL;
*/
}


void 
delete_dat (Data **d) 
{
  XFREE ((*d)->str_dna);
  XFREE ((*d)->dna);
  XFREE ((*d)->str_pro);
  XFREE ((*d)->pro);
  (*d)->str_dna = NULL;
  (*d)->dna = NULL;
  (*d)->str_pro = NULL;
  (*d)->pro = NULL;
/*
  XFREE (*d);
  (*d) = NULL;
*/
}

void 
delete_energy (Interaction **i, int AAlen)
{
  int j, w, x, AAcol, aa2;
  for (j = 0; j < NUM_AMINOACID; j++)
  {
      XFREE ((*i)->nAccess[j]);
      /* delete [] (*i)->nAccess[j]; */
      (*i)->nAccess[j] = NULL;
  }
  XFREE ((*i)->nAccess);
  /* delete [] (*i)->nAccess; */
  (*i)->nAccess = NULL;

  for (j = 0; j < NUM_CODON_EXCEPT_STOP; j++)
  {
      XFREE ((*i)->neighList[j]);
      /* delete [] (*i)->neighList[j]; */
      (*i)->neighList[j] = NULL;
  }
  XFREE ((*i)->neighList);
  /* delete [] (*i)->neighList; */
  (*i)->neighList = NULL;

  for (w = 0; w < AAlen; w++)
  {
      int t = (*i)->AAinfo[w].numNeighbors;
      for (x = 0; x < t; x++)
      {
          XFREE ((*i)->energyAcc[w][x]);
          /* delete [] (*i)->energyAcc[w][x]; */
          (*i)->energyAcc[w][x] = NULL;
      }
      XFREE ((*i)->energyAcc[w]);
      /* delete [] (*i)->energyAcc[w]; */
      (*i)->energyAcc[w] = NULL;
  }
  XFREE ((*i)->energyAcc);
  /* delete [] (*i)->energyAcc; */
  (*i)->energyAcc = NULL;

  for (AAcol = 0; AAcol < AAlen; AAcol++)
  {
      int numNeigh = (*i)->AAinfo[AAcol].numNeighbors;
      for (aa2 = 0; aa2 < numNeigh; aa2++)
      {
          XFREE ((*i)->AAinfo[AAcol].pairwise[aa2]);
          /* delete [] (*i)->AAinfo[AAcol].pairwise[aa2]; */
          (*i)->AAinfo[AAcol].pairwise[aa2] = NULL;
      }
      XFREE ((*i)->AAinfo[AAcol].pairwise);
      /* delete [] (*i)->AAinfo[AAcol].pairwise; */
      (*i)->AAinfo[AAcol].pairwise = NULL;
  }

  for (j = 0; j < AAlen; j++)
  {
      XFREE ((*i)->AAinfo[j].DDDneighbors);
      /* delete [] (*i)->AAinfo[j].DDDneighbors; */
      (*i)->AAinfo[j].DDDneighbors = NULL;
      XFREE ((*i)->AAinfo[j].solvent);
      /* delete [] (*i)->AAinfo[j].solvent; */
      (*i)->AAinfo[j].solvent = NULL;
  }
  XFREE ((*i)->AAinfo);
  /* delete [] (*i)->AAinfo; */
  (*i)->AAinfo = NULL;

/*
  XFREE ((*i));
  (*i) = NULL;
*/
  
  return;
}

void write_dat (FILE *ofile) 
{
  Data *dat = &jones_measure.data;
  int i;
  fprintf(ofile, "============================================================\n");
  fprintf(ofile, "DAT File: %d\n", dat->n_seq);
  fprintf(ofile, "          %d\n", dat->len_dna);
  fprintf(ofile, "          %d\n", dat->len_pro);
  fprintf(ofile, "          %s\n", dat->str_dna);
  fprintf(ofile, "          ");
  for ( i = 0; i < dat->len_dna; i++ ) fprintf(ofile, "%d", dat->dna[i]);
  fprintf(ofile, "\n");
  fprintf(ofile, "          %3d %3d %3d %3d\n", dat->n_a, dat->n_c, dat->n_g, dat->n_t);
  fprintf(ofile, "          %s\n", dat->str_pro);
  fprintf(ofile, "          ");
  for ( i = 0; i < dat->len_pro; i++ ) fprintf(ofile, "%2d ", dat->pro[i]);
  fprintf(ofile, "\n");
  fprintf(ofile, "          %+.3lf %+.3lf\n", dat->solv, dat->pair);
  fprintf(ofile, "============================================================\n");
}

void write_energy (FILE *ofile) 
{
  Interaction *energy = &jones_measure.info;
  int AAlen = jones_measure.data.len_pro;

  int i, j, k, l;
  fprintf(ofile, "============================================================\n");
  fprintf(ofile, "INT File:\n");
  
  /* Solvent Accessibilty */
  fprintf(ofile, "\nSolvent Accessibilty:\n");
  fprintf(ofile, "------------------------------------------------------------\n");
  for (i = 0; i < AAlen; i++) {
     fprintf(ofile, "%d: ", i);
     for (j = 0; j < NUM_AMINOACID; j++) {
        fprintf(ofile, "%+.3lf ", energy->AAinfo[i].solvent[j]);
     }
     fprintf(ofile, "\n");
  } 

  /* Number of Neighbors */
  fprintf(ofile, "\nNumber of Neighbors:\n");
  fprintf(ofile, "------------------------------------------------------------\n");
  for (i = 0; i < AAlen; i++) {
     fprintf(ofile, "%d: ", i);
     fprintf(ofile, "%d ", energy->AAinfo[i].numNeighbors);
     fprintf(ofile, "\n");
  } 

  /* DDD of Neighbors */
  fprintf(ofile, "\nDDD of Neighbors:\n");
  fprintf(ofile, "------------------------------------------------------------\n");
  for (i = 0; i < AAlen; i++) {
     fprintf(ofile, "%d: ", i);
     for (j = 0; j < energy->AAinfo[i].numNeighbors; j++) {
        fprintf(ofile, "%d ", energy->AAinfo[i].DDDneighbors[j]);
     }
     fprintf(ofile, "\n");
  } 

  /* Pairwise Interaction */
  fprintf(ofile, "\nPairwise Interaction:\n");
  fprintf(ofile, "------------------------------------------------------------\n");
  for (i = 0; i < AAlen; i++) {
     fprintf(ofile, "\n*****************************\n");
     fprintf(ofile, "%d - \n", i);
     for (j = 0; j < energy->AAinfo[i].numNeighbors; j++) {
        fprintf(ofile, "%d - %d:\n", i, energy->AAinfo[i].DDDneighbors[j]);
        for (k = 0; k < 20; k++) {
           for (l = 0; l < 20; l++) {
              fprintf(ofile, "%+.3lf ", energy->AAinfo[i].pairwise[j][k*20 + l]);
           }
           fprintf(ofile, "\n");
        }
        fprintf(ofile, "\n");
     }
     fprintf(ofile, "\n");
  } 

  /* energyAcc */
  fprintf(ofile, "\nenergyAcc:\n");
  fprintf(ofile, "------------------------------------------------------------\n");
  for (i = 0; i < AAlen; i++) {
     fprintf(ofile, "*****************************\n");
     fprintf(ofile, "%d - \n", i);
     for (j = 0; j < energy->AAinfo[i].numNeighbors; j++) {
        fprintf(ofile, "%d - %d:\n", i, energy->AAinfo[i].DDDneighbors[j]);
        for (k = 0; k < ACCESS_CATEGORY; k++) {
           fprintf(ofile, "%2d ", energy->energyAcc[i][j][k]);
        }
        fprintf(ofile, "\n");
     }
     fprintf(ofile, "\n");
  } 

  /* nAccess */
  fprintf(ofile, "\nnAccess:\n");
  fprintf(ofile, "------------------------------------------------------------\n");
  for (i = 0; i < NUM_AMINOACID; i++)
  {
      fprintf(ofile, "%2d : ", i);
      for (j = 0; j < NUM_NUCLEOTIDE; j++)
      {
          fprintf(ofile, "%2d ", energy->nAccess[i][j]);
      }
      fprintf(ofile, "\n");
  }

  /* neighList */
  fprintf(ofile, "\nneighList:\n");
  fprintf(ofile, "------------------------------------------------------------\n");
  for (i = 0; i < NUM_CODON_EXCEPT_STOP; i++)
  {
      for (j = 0; j < NUM_COL_NEIGH; j++)
      {
          fprintf(ofile, "%2d ", energy->neighList[i][j]);
      }
      fprintf(ofile, "\n");
  }

  fprintf(ofile, "============================================================\n");
}

int
choose_a_nuc_stationary (int e, int *dna, const int site, double parameter[])
{
  /* 1. translate the DNA to protein sequence 
     2. calculate energy for `different' proteins

     oh, we need to fix the parameters, s and p and a, c, g, t
   */

  int r;
  int i;
  int len_pro = jones_measure.psi.len_pro;
  int *protein = NULL;
  double solv, pair;
  double exp_part[4];
  double solv_part[4];
  double pair_part[4];
  int chosen_nuc;
  
  int site_nuc = dna[site]; 

  protein = XMALLOC (int, len_pro);

  do 
    {
  for (i = 0; i < 4; i++)
    {
      /* find_codon_in_dna (codon, dna, site); */
      dna[site] = i; 
      r = dna2protein (dna, protein, len_pro);
      if (r == EXIT_SUCCESS)
        {
          score_drevol (e, protein, &solv, &pair); 
          solv_part[i] = solv;
          pair_part[i] = pair;
          /* stop codon? */
          /* save them */
          /* PSI_SAMPLE_A: 0 */
          exp_part[i] = log (parameter[i]) 
                        + parameter[PSI_SAMPLE_S] * solv
                        + parameter[PSI_SAMPLE_P] * pair;
        }
      else
        {
          exp_part[i] = -DBL_MAX;
        }
    }

  /* choose one nucleotide */
  chosen_nuc = choose_one_nucleotide (exp_part);
    } while (exp_part[chosen_nuc] == -DBL_MAX);
  
  dna[site] = site_nuc;
  XFREE (protein);

  
  return chosen_nuc; 
}

int
choose_a_nuc_stationary_fast (int e, int *dna, int *pro, const int site, 
                              double log_theta_star[])
{
  /* 1. translate the DNA to protein sequence 
     2. calculate energy for `different' proteins

     oh, we need to fix the parameters, s and p and a, c, g, t
   */

  int i;
  int c[3];
  int chosen_nuc;
  double Energy[4];
  int toAA;
  int AAcol;
  int w;
  double E_solv, E_pair;
  int init_codon;

  w = site % 3;
  AAcol = site / 3;
  init_codon = init_codon_jones_measure ();
  find_codon_in_dna (c, dna, site);
  do 
    {
      for (i = 0; i < 4; i++)
        {
          c[w] = i; /* This completes the codon */

          /* init codon */
          if (AAcol != 0 || init_codon == 0)
            {
              /* whatAmino[nuc] = Nuc2AATable[c[0]][c[1]][c[2]];  */
              toAA = Nuc2AATable[c[0]][c[1]][c[2]]; 
            }
          else 
            {
              /* whatAmino[nuc] = Nuc2AATableStart[c[0]][c[1]][c[2]];  */
              toAA = Nuc2AATableStart[c[0]][c[1]][c[2]]; 
            }
    
          if (toAA == PSI_AA_STP)
            {
              /* whatNRG[nuc] = -DBL_MAX; Fix it */
              Energy[i] = -DBL_MAX;  /* Fix it */
              continue;
            }
          if (e == PSI_ENERGY_JONES)
            {
              /* avoid double calculation in the case of the same AA */
              energy_res_solv_is (AAcol, toAA, &E_solv);
              /* pro's seq should not change */
              energy_res_pair_is (pro, AAcol, toAA, &E_pair); 
            }
          else 
            {
              psi_fatal ("no implementation");
            }
 
          Energy[i] = log_theta_star[i] 
                      + E_solv * log_theta_star[PSI_SAMPLE_S] 
                      + E_pair * log_theta_star[PSI_SAMPLE_P];
        }

      /* choose one nucleotide */
      chosen_nuc = choose_one_nucleotide (Energy);
    } while (Energy[chosen_nuc] == -DBL_MAX);
  
  c[w] = chosen_nuc;
  
  if (AAcol != 0 || init_codon == 0)
    pro[AAcol] = Nuc2AATable[c[0]][c[1]][c[2]]; 
  else
    pro[AAcol] = Nuc2AATableStart[c[0]][c[1]][c[2]]; 
  dna[site] = chosen_nuc;

  return chosen_nuc; 
}

int
psi_energy_dna_report (int e, int *dna, double parameter[])
{
  /* 1. translate the DNA to protein sequence 
     2. calculate energy for `different' proteins

     oh, we need to fix the parameters, s and p and a, c, g, t
   */

  int r;
  int len_pro = jones_measure.psi.len_pro;
  int len_dna = jones_measure.psi.len_dna;
  int *protein = NULL;
  double solv, pair;
  

  protein = XMALLOC (int, len_pro);
  r = dna2protein (dna, protein, len_pro);
  assert (r == EXIT_SUCCESS);
  score_drevol (e, protein, &solv, &pair); 
  XFREE (protein);

  static int gen = 0;
  if (want_verbose == 1)
    {
      int na, nc, ng, nt;
      na = num_nucleotide (dna, len_dna, PSI_DNA_A);
      nc = num_nucleotide (dna, len_dna, PSI_DNA_C);
      ng = num_nucleotide (dna, len_dna, PSI_DNA_G);
      nt = len_dna - na - nc - ng;
      fprintf (stderr, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t", 
               gen,
               parameter[PSI_SAMPLE_S] * solv, 
               parameter[PSI_SAMPLE_P] * pair,
               (double) na / len_dna, 
               (double) nc / len_dna, 
               (double) ng / len_dna, 
               (double) nt / len_dna);
      gen++;
    }
  return EXIT_SUCCESS;
}

int 
psi_energy_dna_state (int e, int *dna, double parameter[], double **states, int n)
{
  int r;
  int len_pro = jones_measure.psi.len_pro;
  int len_dna = jones_measure.psi.len_dna;
  int *protein = NULL;
  double solv, pair;
  

  protein = XMALLOC (int, len_pro);
  r = dna2protein (dna, protein, len_pro);
  assert (r == EXIT_SUCCESS);
  score_drevol (e, protein, &solv, &pair); 
  XFREE (protein);

  int na, nc, ng, nt;
  na = num_nucleotide (dna, len_dna, PSI_DNA_A);
  nc = num_nucleotide (dna, len_dna, PSI_DNA_C);
  ng = num_nucleotide (dna, len_dna, PSI_DNA_G);
  nt = len_dna - na - nc - ng;
 
  states[PSI_SAMPLE_A][n] = (double) na / len_dna;
  states[PSI_SAMPLE_C][n] = (double) nc / len_dna; 
  states[PSI_SAMPLE_G][n] = (double) ng / len_dna;
  states[PSI_SAMPLE_T][n] = (double) nt / len_dna;
  states[PSI_SAMPLE_S][n] = parameter[PSI_SAMPLE_S] * solv;
  states[PSI_SAMPLE_P][n] = parameter[PSI_SAMPLE_P] * pair;

  return EXIT_SUCCESS;
}

int
psi_energy_valid_state (int e, int *dna, double parameter[], double **state_value, int n)
{
  int r;
  int len_pro = jones_measure.psi.len_pro;
  int len_dna = jones_measure.psi.len_dna;
  int *protein = NULL;
  double solv, pair;

  protein = XMALLOC (int, len_pro);
  r = dna2protein (dna, protein, len_pro);
  assert (r == EXIT_SUCCESS);
  score_drevol (e, protein, &solv, &pair); 
  XFREE (protein);

  int na, nc, ng, nt;
  na = num_nucleotide (dna, len_dna, PSI_DNA_A);
  nc = num_nucleotide (dna, len_dna, PSI_DNA_C);
  ng = num_nucleotide (dna, len_dna, PSI_DNA_G);
  nt = len_dna - na - nc - ng;
 
  double v_a = (double) na / len_dna;
  double v_c = (double) nc / len_dna; 
  double v_g = (double) ng / len_dna;
  double v_t = (double) nt / len_dna;
  double v_s = parameter[PSI_SAMPLE_S] * solv;
  double v_p = parameter[PSI_SAMPLE_P] * pair;

  double lbq = 0.4;
  double ubq = 0.6;
  double lb = psi_stats_quantile (state_value[PSI_SAMPLE_A], n, lbq);
  double ub = psi_stats_quantile (state_value[PSI_SAMPLE_A], n, ubq);
  if (v_a < lb || v_a > ub)
    {
      fprintf (stderr, "A: %lf - %lf: %lf\n", lb, ub, v_a);
      return EXIT_FAILURE;
    }

  lb = psi_stats_quantile (state_value[PSI_SAMPLE_C], n, lbq);
  ub = psi_stats_quantile (state_value[PSI_SAMPLE_C], n, ubq);
  if (v_c < lb || v_c > ub)
    {
      fprintf (stderr, "C: %lf - %lf: %lf\n", lb, ub, v_c);
      return EXIT_FAILURE;
    }

  lb = psi_stats_quantile (state_value[PSI_SAMPLE_G], n, lbq);
  ub = psi_stats_quantile (state_value[PSI_SAMPLE_G], n, ubq);
  if (v_g < lb || v_g > ub)
    {
      fprintf (stderr, "G: %lf - %lf: %lf\n", lb, ub, v_g);
      return EXIT_FAILURE;
    }

  lb = psi_stats_quantile (state_value[PSI_SAMPLE_T], n, lbq);
  ub = psi_stats_quantile (state_value[PSI_SAMPLE_T], n, ubq);
  if (v_t < lb || v_t > ub)
    {
      fprintf (stderr, "T: %lf - %lf: %lf\n", lb, ub, v_t);
      return EXIT_FAILURE;
    }

  lb = psi_stats_quantile (state_value[PSI_SAMPLE_S], n, lbq);
  ub = psi_stats_quantile (state_value[PSI_SAMPLE_S], n, ubq);
  if (v_s < lb || v_s > ub)
    {
      fprintf (stderr, "S: %lf - %lf: %lf\n", lb, ub, v_s);
      return EXIT_FAILURE;
    }

  lb = psi_stats_quantile (state_value[PSI_SAMPLE_P], n, lbq);
  ub = psi_stats_quantile (state_value[PSI_SAMPLE_P], n, ubq);
  if (v_p < lb || v_p > ub)
    {
      fprintf (stderr, "P: %lf - %lf: %lf\n", lb, ub, v_p);
      return EXIT_FAILURE;
    }

  fprintf (stderr, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
           parameter[PSI_SAMPLE_A],
           parameter[PSI_SAMPLE_C],
           parameter[PSI_SAMPLE_G],
           parameter[PSI_SAMPLE_T],
           parameter[PSI_SAMPLE_S],
           parameter[PSI_SAMPLE_P]);
  fprintf (stderr, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
           psi_stats_median (state_value[PSI_SAMPLE_A], n),
           psi_stats_median (state_value[PSI_SAMPLE_C], n),
           psi_stats_median (state_value[PSI_SAMPLE_G], n),
           psi_stats_median (state_value[PSI_SAMPLE_T], n),
           psi_stats_median (state_value[PSI_SAMPLE_S], n),
           psi_stats_median (state_value[PSI_SAMPLE_P], n));
  fprintf (stderr, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
           v_a, v_c, v_g, v_t, v_s, v_p);

  return EXIT_SUCCESS;
}

int
choose_a_nuc_rate (int e, int *site, int *nucleotide, int *dna, double parameter[])
{
  /* 1. translate the DNA to protein sequence 
     2. calculate energy for `different' proteins

     oh, we need to fix the parameters, s and p and a, c, g, t
   */

  int r;
  int i, j, k;
  int nuc_dna_i;
  int len_dna = jones_measure.psi.len_dna;
  int len_pro = jones_measure.psi.len_pro;
  int *protein = NULL;
  double solv, pair;
  double solv_i, pair_i;
  double *exp_part = NULL;
  
  exp_part = XMALLOC (double, len_dna * 3);
  protein = XMALLOC (int, len_pro);

  r = dna2protein (dna, protein, len_pro);
  score_drevol (e, protein, &solv_i, &pair_i); 

  for (i = 0; i < len_dna; i++)
    {
      nuc_dna_i = dna[i];
      k = 0;
      for (j = 0; j < 4; j++)
        {
          if (j == nuc_dna_i)
            continue;
          dna[i] = j;
          r = dna2protein (dna, protein, len_pro);
          if (r == EXIT_SUCCESS)
            {
              score_drevol (e, protein, &solv, &pair); 
              /* stop codon? */
              /* save them */
              /* PSI_SAMPLE_A: 0 */
              exp_part[i*3 + k] = log (parameter[j]) 
                            + parameter[PSI_SAMPLE_S] * (solv - solv_i) / 2.0
                            + parameter[PSI_SAMPLE_P] * (pair - pair_i) / 2.0;
            }
          else
            {
              exp_part[i*3 + k] = -DBL_MAX;
            }
          k++;
        }
      dna[i] = nuc_dna_i;
    }

  XFREE (protein);

  /* choose one nucleotide */
  int w = choose_one_sequence_by_rate (exp_part, len_dna);
  /* use w to determin site and nucleotide */
/* HERE WE WERE */
  /*
     ACGT ACGT ACGT ACGT
     ^012 0^12 01^2 012^ -> dna[w/3],w%3(=0),w%3(=1),w%3(=2)
     |+++  |++   |+    |
     |123 0|23 01|3 012|
    
     ^: dna[w/3]
   */
  *site = w / 3;
  if (dna[w / 3] <= w % 3)
    *nucleotide = (w % 3) + 1;
  else
    *nucleotide = (w % 3);
  XFREE (exp_part);
  return EXIT_SUCCESS; 
}
