#include <iostream.h> // This program has several parameters that can be updated
#include <stdlib.h>   // or fixed.  One parameter in particular is time.  The
#include <fstream.h>  // time should be fixed.  In this model, the overall rate
#include <iomanip.h>  // parameter U os added. This allows the model to estimate
                      // the true branch length between the sequences.This model
#include "tools.h"    // is time reversible, because the exponential term is now
#include "tools.c"    // squared.  I hope this will work.
#include "eigen.c"    //////////////////////////////////////////////////////////

// Options: Doug's paper
#define BURNIN            0 // This is how long we wait until we take a sample
#define SAMPLE_PATH    100000 // Number of sampled paths
#define SAMPFREQ          1 // Number of MCMC cycles between samples

//#define BURNIN       100000   // This is how long we wait until we take a sample
//#define SAMPLE_PATH  10000 //Number of sampled paths
//#define SAMPFREQ     500 //Number of MCMC cycles between samples
#define MAXITERATION  (BURNIN +((SAMPLE_PATH-1)*SAMPFREQ)) // This is the # of
                                            // iterations of the MCMC routine
#define NODE_UPDATE    1    // how many cycles between node updates
#define DELTA_PATH     5    // For fixed node sequences, how many times per
                            // MCMC cycle should each branch have a randomly
                            // selected site path proposed?
#define ITERATE        1    // The number of times that seq will be checked
#define SEQSET        1000  // How many sequences we simulate for S and P update
#define SIZEOFGRID     9   // How big of a grid we will use for the simulation
#define LOW_A         0.14  // Set these all to be the same now, but later
#define HIGH_A        0.36
#define LOW_C         0.14  // Set these all to be the same now, but later
#define HIGH_C        0.36
#define LOW_G         0.14  // Set these all to be the same now, but later
#define HIGH_G        0.36
#define LOW_T         0.14 // Set these all to be the same now, but later
#define HIGH_T        0.36  // we can have separate range values for each nuc
#define DELTA_PI      0.01
#define K1           1.857  //
#define L1            0.1   //
#define SIM_PI          1
//////////////////////////////
#define FREQGRIDSIZE   12   // SEE the function SeqFreq for the implementation
#define ADENINE    0.281118 // If the user defines, then this is Pi_A
#define CYTOSINE   0.204641 // If the user defines, then this is Pi_C
#define GUANINE    0.267932 // If the user defines, then this is Pi_G
//////////////////////////////
#define SIM_T          0    // [0,1] tells whether Time will be updated
#define INIT_T       1.0    // The evolutionary time separating the sequences
#define DELTA_T      0.01   // This is the step size for the time proposal
#define T_MIN        0.0    // Uniform prior MIN value for time
#define T_MAX        1.0    // Uniform prior MAX value for time
//////////////////////////////
#define SIM_U          1    // [0,1] will U (branch rate) be updated?
#define INIT_U       1.0    // The evolutionary time searating the sequences
#define DELTA_U      0.2    // This is the step size for the time proposal (it was 0.05)
#define U_MIN        0.0    // Uniform prior MIN value for time
#define U_MAX        4.0    // Uniform prior MAX value for time
//////////////////////////////
#define SIM_K          1    // [0,1] tells whether kappa will be updated
#define MULTI_K        0    // 0: One K on entire topo, 1:K_i per branch of topo
#define INIT_K       2.0    // The initial seed value of kappa
#define DELTA_K      1.0    // This is the step size for the kappa proposal
#define K_MIN        0.0    // Uniform prior MIN value for kappa
#define K_MAX       10.0    // Uniform prior MAX value for kappa
//////////////////////////////
#define SIM_W          0    // [0,1] tells whether omega will be updated
#define MULTI_W        0    // 0: One W on entire topo, 1:W_i per branch of topo
#define INIT_W       1.0    // Behaves like synonymous / nonsynonymous parameter
#define DELTA_W      0.5    // This is the step size for the omega proposal
#define W_MIN        0.0    // Uniform prior MIN value for omega
#define W_MAX       10.0    // Uniform prior MAX value for omega
//////////////////////////////
#define SIM_W_S        1    // [0,1]tells whether omega(solvent) will be updated
#define INIT_W_S    -1.0    // Makes NRG's behave like a quantitative measure
#define DELTA_W_S    0.4    // This is the step size for the omega_S proposal
#define W_S_MIN     -2.0    // Uniform prior MIN value for omega_S
#define W_S_MAX      2.0    // Uniform prior MAX value for omega_S
//////////////////////////////
#define SIM_W_P        1    // [0,1] tells whether omega(pairs) will be updated
#define INIT_W_P    -0.1    // Makes NRG's behave like a quantitative measure
#define DELTA_W_P    0.03   // This is the step size for the omega_P proposal
#define W_P_MIN     -0.15   // Uniform prior MIN value for omega_P
#define W_P_MAX      0.15   // Uniform prior MAX value for omega_P
//////////////////////////////

#define ANCESTRAL_NUC 1     // 0: print out AA seqs, 1: print out nucleotide seqs
#define INTERACT     10.0   // Distance at which AA's can biologically interact
#define BIGNUM     10000000.0  // Used as an upper bound on the logarithm
#define INVBIGNUM  0.0000001   // Used as a lowerbound on the logarithm
#define LOGBIGNUM 16.118095651 // This is the LN of the BIGNUM
#define LOUD         0      // How much STUFF do you want? 0 = little, 1 = TONS!
#define WHERE_ARE_WE 1
//////////////////////////////

//* needed by rnd() routine */
typedef long longer[6];
static longer seed;
long newseed0, newseed1, newseed2, newseed3, newseed4, newseed5, newseed6;

struct AAsiteInfo{
        double *solvent;
        double **pairwise;
        int numNeighbors;
        int *DDDneighbors;
};

struct Interaction{
        AAsiteInfo *AAinfo;
        int ***energyAcc; // stores the pairwise distance between all CB atoms
        int **nAccess; //int list that obtains correct index for neighbor matrix
        int **neighList; // stores the AA neighbors, correct Nuc, and K param
        double ****energy; // stores all matrices for pairpotential NRG
};

struct Substitution{ //stores all the substitution information
	double time; //exact times of the substitutions
        double genHood;
        double subHood;
        int column;  // identifies the location of the substitution
        int priorNuc; // What nucleotide the substitution was from
        int pathNuc; // What nucleotide was chosen for the sub
        int amino; // Holds the amino acid that this sub codes for
        int priorAA; // Holds the amino acid that this sub came from
        double nrgDiff; // Holds the specific rate away for this change
        double solvNRGdiff; // Holds the solvent energy difference
        double expNRGdiff; // Holds the exponential of the specific rate away
        Substitution *UpathPtr;//points to previous node in path(ordered in time)
        Substitution *DpathPtr;//points to next node in path (ordered in time)
        Substitution *UcolPtr;   // points to the previous node in the column
		// if it is the initial sequence node, then it points to NULL
		// we can use this to find previous nucleotide for testing
        Substitution *DcolPtr;   // points to the next node in the column
		// used for quickly traversing the column to delete
};

struct Path{
        int nuc; // This is the nucleotide in the initial sequence
        int numsub; // This is the number of subs in the column
        Substitution *firstColSub; // This points to the first substitution
};

struct Bayes{
        double *kappa;         //Will hold the kappa parameter for the simulation
        double *omega;         //Will hold the omega parameter for the simulation
        double solvent;       //Will hold the solv parameter for the simulation
        double pairwise;      //Will hold the pair parameter for the simulation
        double tempParameter; //Will be the proposed value for all parameters
        double *rate;         //Will hold the rates for all branches
        double *tempRate;     //Will be a place holder for the rate parameters
        double *yangBranch;   //Will hold the branch lengths in the Yang form
        double *branchLen;    //Will hold the branch lengths for all branches
        double *MCMCnf;       //Will hold the pi parameters for the simulation
        double *tempMCMCnf;   //Will hold the proposed pi parameters
        //These hold the energies of the ancestral sequence and the
        //      values that have to do with the Gibbs sampler
        //  double maxValue;  Holds the max prob of pi_i's of simulated seqs
        //  double nrgConst;  This is the illustrious C constant
        double     solvNRG,     pairNRG,     nrgConst,     maxValue;
        double tempsolvNRG, temppairNRG, tempNRGconst, tempMAXvalue;
        //These are the number of each nucleotide in the ancestral sequence
        int number_A, number_C, number_G;
        // Holds the solvent coordinate of the grid
        int    A_gridPt, tempA_gridPt; // Holds A nuc    coordinate of the grid
        int    C_gridPt, tempC_gridPt; // Holds C nuc    coordinate of the grid
        int    G_gridPt, tempG_gridPt; // Holds G nuc    coordinate of the grid
        int    S_gridPt, tempS_gridPt; // Holds solvent  coordinate of the grid
        int    P_gridPt, tempP_gridPt; // Holds pairwise coordinate of the grid
        int numSeq;           //Will hold the number sequences in the data set
        int numBranch;        //Will be a flag to say how many banches in tree
};

struct Information{ // Holds everything!
        int numSeq;     //Will hold the number sequences in the data set
        int numBranch;
        int numNodes;   // Holds the number of nodes in the tree
        int numChild;   // Holds how many children this node has
        int parent;     // Holds the location of the parent node
        int *child;     // Holds the list of children of this node
        int *parentSeq; // Points to the parent sequence of this node
        int *seq;        // nucleotide sequences at this particular node
        int len;          // nucleotide sequence length
        double invlen;    // 1.0 / len
        int AAlen;        // Amino Acid sequence length
        double invAAlen;   // 1.0 / AAlen
        int *AA_seq;      // initial constant Amino Acid sequence
        int *parentAAseq; // Points to the parent AA sequence (at top of branch)
        int *neighAcc;    // Holds the row of neighList for all sites
        double AAseqNRG;  // The energy of the initial sequence
        int totNumSub;    // Total number of subs in this path
        Substitution *startPath; // This is the start of the entire path
        double pathHood;  // holds the likelihood of the entire path
        double avgRate;   // Holds what will be the branch length
        double probLastEvent; // Holds prob of no event in last time interval
        Path *seqPath; // Holds the sequence path for this particular branch
};

struct GibbsPart{
   double S, P;
   int A, C, G;
};

   int Nuc2AATable[4][4][4]; //Global array used for hash-table
   GibbsPart *GibbsInfo[SIZEOFGRID][SIZEOFGRID][FREQGRIDSIZE][FREQGRIDSIZE][FREQGRIDSIZE];
   double S_POINTS[SIZEOFGRID];
   double P_POINTS[SIZEOFGRID];
   double A_POINTS[FREQGRIDSIZE];
   double G_POINTS[FREQGRIDSIZE];
   double C_POINTS[FREQGRIDSIZE];
   double T_POINTS[FREQGRIDSIZE];
   double F84_freq[12];
   double F84_RATE[4][4];
   double F84_PROB[4][4];
   double POISSON_G[16];
   double POISSON_W[16];
   int    COMPARE[4][4];

//   double S_NRGholder[SIZEOFGRID][SIZEOFGRID][SEQSET];
//   double P_NRGholder[SIZEOFGRID][SIZEOFGRID][SEQSET];
#define Dist(x1,y1,z1,x2,y2,z2) sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));

int ParseCommandLine(char *preferred_defaultfilename,
                     char **infile_filename,  char **inseed_filename,
                     char **interaction_filename,
                     char **infoout_filename, char **posterior_filename,
                     char **pathout_filename, char **roundout_filename,
                     int argc, char *argv[]);
//******************************************************************************
int main(int argc, char *argv[])
{
   char *infile_filename, *inseed_filename, *interaction_filename,
        *infoout_filename;
   char *posterior_filename, *pathout_filename, *roundout_filename;

   if (!ParseCommandLine("SeqData.dat",
      &infile_filename,  &inseed_filename,&interaction_filename,
      &infoout_filename, &posterior_filename,&pathout_filename,
      &roundout_filename,argc, argv))
   {
      exit(1);
   }
//ifstream in("1AQBe2h.dat");           // This is the input nucleotide sequences
//ifstream seedling("inseed.dat");      // Gives the random number seed
//ifstream inter("Interaction.dat");    // Holds all constant info + NRG's
//ofstream out("1AQBe2h.info.out");     // Random information output file
//ofstream out1("1AQBe2h.out");         // kappa + omega : check for convergence
//ofstream out2("1AQBe2h.NRGpic.out");
   ofstream outP1("X5noSnoP_1.out");
   ofstream outP2("X5noSnoP_2.out");
//   ofstream outP3("WWfast_3.out");
//   ofstream outP4("WWfast_1.NRG.out");
//   ofstream outP5("WWfast_2.NRG.out");
//   ofstream outP6("WWfast_3.NRG.out");

//   cerr << "in\t:" << infile_filename << endl;
//   cerr << "seed\t:" << inseed_filename << endl;
//   cerr << "inter\t:" << interaction_filename << endl;
//   cerr << "out\t:" << infoout_filename << endl;
//   cerr << "out1\t:" << posterior_filename << endl;
//   cerr << "out2\t:" << pathout_filename << endl;
//   cerr << "roundout\t:" << roundout_filename << endl;
//   exit (0);


   ifstream in(infile_filename);        // This is the input nucleotide sequences
   ifstream seedling(inseed_filename);  // Gives the random number seed
   ifstream inter(interaction_filename);// Holds all constant info + NRG's
   ofstream out(infoout_filename);      // Random information output file
   ofstream out1(posterior_filename);   // kappa + omega : check for convergence
   ofstream out2(pathout_filename);     // path energy terms

   if (!in.is_open() )      {cerr << "Error: Could not open infile '" << infile_filename << "' !" << endl; exit(1);}
   if (!seedling.is_open() ){cerr << "Error: Could not open inseedfile '" << inseed_filename << "' !" << endl; exit(1);}
   if (!inter.is_open() )   {cerr << "Error: Could not open interactionfile '" << interaction_filename << "' !" << endl; exit(1);}
   if (!out.is_open() )     {cerr << "Error: Could not open outfile '" << infoout_filename << "' !" << endl; exit(1);}
   if (!out1.is_open() )    {cerr << "Error: Could not open posteriorfile '" << posterior_filename << "' !" << endl; exit(1);}
   if (!out2.is_open() )    {cerr << "Error: Could not open pathenergytermsfile '" << pathout_filename << "' !" << endl; exit(1);}

   double rnd(long* seed);
   void nucCompare();
   int Nuc2AA(int *c);
   char Int2AA(int c);
   char IntToNuc(int m);
   void initBayes(Bayes *pE);
   void initInteraction(Interaction *I);
   void SetUpInfo(Bayes *pE, Information *C, Information *P, long* seed,
                ifstream &in, ofstream &out);
   void SetUpInteraction(Interaction *FirstOrder, int AAlen, ifstream& in,
                ifstream &inter, ofstream& out);
   void CalcInitSeqNRG(Bayes *pE, Information *C, Information *P,
                Interaction *FO, ofstream &o);
   void SetNRGholder(Bayes *pE);
   double CalcConst(Bayes *pE, Information &P, Interaction *FO, double &maximum,
                long *seed, int round, ofstream &o);
   void CalcNumMuts(Information &Current, Path *currentPath, long* seed,
                ofstream& out );
   void CopyInitColumn(Path *emptyPath, Path *fullPath, int column, int AAlen);
   void InsertFirstColumn(Information &Info, Path *columnPath, int col);
   void InsertColumn(Information &Info, Path *columnPath, int col);
   void FirstParam(Bayes *pE, Information *C, Information *P, Interaction *FO,
                ofstream &o);
   void NewKappa(Bayes *pE, Information *C, Information *P, Interaction *FO,
         double ****RateNRG, int &accept_K, double *cf, long *seed, ofstream &o);
   void NewOmega(Bayes *pE, Information *C,Information *P,Interaction *FO,
        double ****RateNRG, int &accept_W, double *cf, long *seed, ofstream &o);
   void NewOmega_S(Bayes *pE, Information *C,Information *P,Interaction *FO,
                double ****RateNRG, int &accept_WS, long *seed, ofstream &o);
   void NewOmega_P(Bayes *pE, Information *C,Information *P,Interaction *FO,
                double ****RateNRG, int &accept_WP, long *seed, ofstream &o);
   void NewRate(Bayes *pE, Information *C,Information *P,Interaction *FO,
                double ****RateNRG, int &accept_U,  long *seed, ofstream &o);
   void NewPath(Bayes *pE, Information *Ctot, Information *Ptot,Interaction *FO,
                long *seed, int &totP, int &yesP, int &same_P,ofstream& o);
   void NewPi(Bayes *pE, Information *C,Information *P,Interaction *FO,
                double ****RateNRG, int &accept_pi, int whichNuc, long *seed,
                ofstream &o);
   void NewNodeSeq(Bayes *pE,Information *Ctot,Information *Ptot,
        Interaction *FO,long *seed,int &totP,int &yesP,int &same_P,ofstream& o);
   void PathPosterior(Bayes *pE, Substitution *tmp, Interaction *FO,int count,
        int q, double *nE, double *mT, int *newOrder, int *picSeq, ofstream &o);
   void DeleteSome(Bayes *pE, Information *C, Information *P,
                 long* seed, ofstream& out);
   void DeleteAll(Bayes *pE, Interaction *FO, Information *C, Information *P,
                ofstream& o);
   void PrintPath(Substitution *start, ofstream& o);
   void PathProfile(Information *C,Interaction *FO,double **profile);
   void PathOutput(double ***myPath, int *PathSubsCount, ofstream &o);
   void NRGposterior(Bayes *PE,Information &C,Interaction *FO, int q, int count,
        double ***RateNRG, ofstream &o);
   void ErrorHandle(char *statement, Bayes *pE, Interaction *FirstOrder,
              Information *C, Information *P,ofstream& o);
   void PathPostFull(Bayes *pE, Substitution *tmp, Interaction *FO,int count,
         int q, double *nE, double *mT, int *newOrder, double **pathPtr);
   void codonFrequency(double *cf, double *nfreq);
   void YangCodonRate(double *cfreq, double *c, double k, double w);
   double CalcYangRate(double bl,double *cf, double *cr,int AAlen, int *neigh);
        //Build the hash-table for the fast version of function Nuc2AA
        int c[3];
        for (c[0]=0;c[0]<4;c[0]++){
            for (c[1]=0;c[1]<4;c[1]++){
               for (c[2]=0;c[2]<4;c[2]++){
                   Nuc2AATable[c[0]] [c[1]] [c[2]] = Nuc2AA(c);
               }
            }
         }

         double a,aaa,ccc,ggg,ttt,b;
         a = (W_S_MAX - W_S_MIN)/SIZEOFGRID;
         b = (W_P_MAX - W_P_MIN)/SIZEOFGRID;
         S_POINTS[0] = W_S_MIN + (a/2.0);
         P_POINTS[0] = W_P_MIN + (b/2.0);
         out << endl << "These are the grid points that will be used " << endl;
         out << "S = " << S_POINTS[0] << " and P = " << P_POINTS[0] << endl;
         for(int i = 1; i < SIZEOFGRID; i++)
         {
            S_POINTS[i] = S_POINTS[i-1] + a;
            P_POINTS[i] = P_POINTS[i-1] + b;
            out << "S = " << S_POINTS[i] << " and P = " << P_POINTS[i] << endl;
         }
         out << endl;
         // NULL out each of the GibbsInfo pointers in this massive for loop
         for(int i=0; i < SIZEOFGRID; i++){
           for(int j=0; j < SIZEOFGRID; j++){
             for(int k=0; k < FREQGRIDSIZE; k++){
               for(int l=0; l < FREQGRIDSIZE; l++){
                 for(int m=0; m < FREQGRIDSIZE; m++){
                    GibbsInfo[i][j][k][l][m] = NULL;
                 }
               }
             }
           }
         }
         aaa = (HIGH_A - LOW_A)/(FREQGRIDSIZE-1);
         ccc = (HIGH_C - LOW_C)/(FREQGRIDSIZE-1);
         ggg = (HIGH_G - LOW_G)/(FREQGRIDSIZE-1);
         ttt = (HIGH_T - LOW_T)/(FREQGRIDSIZE-1);
         A_POINTS[0]= LOW_A;
         C_POINTS[0]= LOW_C;
         G_POINTS[0]= LOW_G;
         T_POINTS[0]= LOW_T;
         out << "A C G T = " << A_POINTS[0] << " " << C_POINTS[0] << " " << G_POINTS[0] << " " << T_POINTS[0] << endl;
         for(int i = 1; i < FREQGRIDSIZE; i++)
         {
            A_POINTS[i] = A_POINTS[i-1] + aaa;
            C_POINTS[i] = C_POINTS[i-1] + ccc;
            G_POINTS[i] = G_POINTS[i-1] + ggg;
            T_POINTS[i] = T_POINTS[i-1] + ttt;
            out << "A C G T = " << A_POINTS[i] << " " << C_POINTS[i] << " " << G_POINTS[i] << " " << T_POINTS[i] << endl;
         }
         out << "(NOTE: T grid points are not used !!)" << endl;
         out << endl;

// Changes by JEFF around here
//         A_POINTS[0]=C_POINTS[0]=G_POINTS[0]=T_POINTS[0]= LOW_A;
//         out << endl << "These are the Frequency grid points that will be used " << endl;
//         out << "A C G T = " << A_POINTS[0] << endl;
//         for(int i = 1; i < FREQGRIDSIZE; i++){
//           C_POINTS[i]=G_POINTS[i]=T_POINTS[i]= A_POINTS[i] = A_POINTS[i-1] + a;
//            out << "A C G T = " << A_POINTS[i] << endl;
//         }
/****************************************************
>From Jeff Thorne and or Joe Felsenstein (PHYLIP)
program used in Thorne divtime program.  This is
a great random number generator that works extremely
well.  One must insert a seed first to get the
calculation working well.
****************************************************/
        long initseed;
        double lastrnd;
        initseed = 0;

        // This reads in the random seed from the file
        seedling >> initseed;
        out << "This is the initial seed : " << initseed << endl;
        //We now place it on top of our result file for future use
//        out1 << initseed << endl << endl;
        for (int i = 0; i <= 5; i++) {seed[i] = 0;}
        int i = 0;
        do {
        seed[i] = initseed & 63;
        initseed /= 64;
        i++;
        } while (initseed != 0);
        for(int i=1; i<= 100; i++){lastrnd = rnd(seed);}
        seedling.close();
/*********************************************************
End Random number generator Progrm: Now begin main thrust
of the program.  It is all set up, now we must use it! We
throw out the first 500 iterations as a burn in time, and
hopefully we will be at a higher portion of the Surface.
*********************************************************/
        Bayes *PostEstimate = new Bayes;
        initBayes(PostEstimate);
        in >> PostEstimate->numSeq;
        out << "This is numSeq " << PostEstimate->numSeq << endl;
        PostEstimate->numBranch = 2*PostEstimate->numSeq - 3; ;
        out << "This is numBranch " << PostEstimate->numBranch << endl;
        // The Current structure will hold the current step of the iteration
        Information *Current  = new Information [PostEstimate->numBranch+1];
//        Current[0]->seq = new int *[2];
        // The Proposed structure will hold the proposed step of the iteration
	Information *Proposed = new Information [PostEstimate->numBranch+1];
        SetUpInfo(PostEstimate, Current, Proposed, seed, in, out);
        out << "Done with setting up the info " << endl;
        nucCompare();

///////////////////////////////////////////////////////////////
// Now we are going to set up the Interaction Structure.     //
// This structure holds all of the information calcualted    //
// from David Jones' program including contact accessibility,//
// solvent accessibility, and pair potential energies.       //
///////////////////////////////////////////////////////////////
        Interaction *FirstOrder = new Interaction;
        initInteraction(FirstOrder);
        SetUpInteraction(FirstOrder,Current[0].AAlen,in,inter,out);
        inter.close();
        int nucSeqLen, AAseqLen, numSeq, numBranch, numNodes, numInternal;
        nucSeqLen = Current[0].len; AAseqLen = Current[0].AAlen;
        numSeq = PostEstimate->numSeq;
        numBranch = PostEstimate->numBranch;
        numNodes = numBranch+1;
        numInternal = numNodes - numSeq;
        out << "This is the sequence len = " << nucSeqLen;
        out << " and the aa len = " << AAseqLen << endl;

        // Now we check the amino acid sequences
        int codon[3], e;
        for(int h = 0; h < numNodes; h++)
        {
           for(int j = 0; j < nucSeqLen; j++){
              codon[j%3] = Current[h].seq[j];
              if(j%3 == 2){
                 e = Nuc2AATable[codon[0]] [codon[1]] [codon[2]];
                 if(e == 20){
                   ErrorHandle("There is a stop codon in your sequences!!!",
                   PostEstimate,FirstOrder,Current,Proposed,out);
                 }
              }
           }
           out << "No stop codon in sequence " << h << endl;
        }

//        out << "This is nucSeq Len " << nucSeqLen << endl;
//        for(int j = 0; j < nucSeqLen; j++)
//        {
//           out << Current[numBranch].seq[j];
//        }
//        out << endl;

        // This part of the code will hold the parameter posterior distributions
        // It will print out the entire path posterior only at the end of the
        //      entire simulation.
        double **paramOut = new double *[SAMPLE_PATH];
        int paramnum = 6;
        paramnum += (5*numBranch);
        if(MULTI_K){paramnum += numBranch;}
        else{paramnum++;}
        if(MULTI_W){paramnum += numBranch;}
        else{paramnum++;}
        for(int h = 0; h < SAMPLE_PATH; h++){paramOut[h] = new double[paramnum];}

        double **paramOut1, *parameterPtr, cf[61], c_rate[61];
        paramOut1 = paramOut;

       // Similar to above, this part will store all of the paths for
        // each branch of the tree.  It will print them out only after the
        // entire simulation is complete.

        //WE WILL HAVE TO SEE HOW THIS WORKS!!!

        double ****PathStorage = new double ***[numBranch];
        int **PathSubsCount = new int *[numBranch];
        for(int h =0; h < numBranch;h++)
        {
           PathStorage[h] = new double **[SAMPLE_PATH];
           PathSubsCount[h] = new int[SAMPLE_PATH];
        }



        // This function fills in the sequence energies that will remain fixed

        CalcInitSeqNRG(PostEstimate,Current,Proposed,FirstOrder, out);
        SetNRGholder(PostEstimate);

        // We now set the energy constant for the ancestral sequence as
        // well as the maximum value for

        PostEstimate-> tempNRGconst = PostEstimate->nrgConst =
        CalcConst(PostEstimate, Current[numBranch], FirstOrder,
                PostEstimate->maxValue, seed, 9, out);
        PostEstimate->tempMAXvalue = PostEstimate->maxValue;
        Substitution *tmp;

        double ****RateNRG = new double ***[numBranch];
        int doubleNucLen = 2*nucSeqLen;
        for (int q = 0; q < numBranch; q++)
        {
           RateNRG[q] = new double **[doubleNucLen];
           for(int r = 0; r < doubleNucLen; r++)
           {
              RateNRG[q][r] = new double *[AAseqLen];
              for(int s = 0; s < AAseqLen; s++)
              {
                 RateNRG[q][r][s] = new double[27];
              }
           }
        }
        for(int zz = 0; zz < 61; zz++){c_rate[zz] = cf[zz]=0.0;}
        codonFrequency(cf, PostEstimate->MCMCnf);
        YangCodonRate(c_rate,cf,PostEstimate->kappa[0],PostEstimate->omega[0]);
        for (int q = 0; q < numBranch; q++)
        {
           out << "This is init branch " << q << " "<< PostEstimate->yangBranch[q] << endl;
           PostEstimate->rate[q] = CalcYangRate(PostEstimate->yangBranch[q],cf,
                        c_rate,AAseqLen,Current[q].neighAcc);
           out << "This is the rate of sequence " << q << " "
                << PostEstimate->rate[q] << endl;
        }

////////////////////////////////////////////////////////////////////////////////
// WE ARE NOW READY TO BEGIN THE CALCULATION!!!!!!!/////////////////////////////
////////////////////////////////////////////////////////////////////////////////
   for(int trial = 0; trial < ITERATE; trial++)
   {
      //cout << "This is trial = " << (trial+1) << endl << endl;
      // This part of the code sets up the Current path completely
      // and then copies it completely into the Proposed path
     for(int i = 0; i < numBranch; i++)
     {
         CalcNumMuts(Current[i],Current[i].seqPath,seed,out);
         out << "This is how many subs in path " << i ;
         out << " := " << Current[i].totNumSub << endl;
         Proposed[i].totNumSub = Current[i].totNumSub;
         //We now have to copy everything over to the proposed structure
         //This function orders the columns in time
         Path *cPath, *pPath;
         cPath = Current[i].seqPath;
         pPath = Proposed[i].seqPath;

         for(int col = 0; col < nucSeqLen; col++)
         {
            if(cPath[col].numsub > 0)
            {
               CopyInitColumn(pPath,cPath,col,AAseqLen);
            }
         }
         //After the first, we can add all the rest in automatically
         // Now we must connect everything in time, one column at  time
         // Doug Changes August 2005
         int q,rrrr; q = rrrr = 0;
         while(rrrr == 0)
         {
            if(q < nucSeqLen)
            {
               if(pPath[q].numsub == 0){q++;}
               else
               {
                  rrrr = 1;
                  InsertFirstColumn(Proposed[i], pPath, q);
                  for(int col = q+1; col < nucSeqLen; col++)
                  {
                     if(pPath[col].numsub > 0)
                     {
                        InsertColumn(Proposed[i], pPath, col);
                     }
                  }
               }
            }
            else{rrrr = 1;}
         }
/*
         while(pPath[q].numsub == 0){q++;}
         InsertFirstColumn(Proposed[i], pPath, q);
         for(int col = q+1; col < nucSeqLen; col++)
         {
            if(pPath[col].numsub > 0)
            {
               InsertColumn(Proposed[i], pPath, col);
            }
         }
*/
         out << "This is the time of the first sub " << endl;
   //      out << "This is the time of the first sub " << Current[i].startPath->time << endl;
      }
//      UpdateProbMatrix(Current,Proposed);
#if LOUD//////////////////////////////////////////////////////////////////////
     for(int i = 0; i < numBranch; i++)
     {
        out << "Here is the Current Path now " << endl << endl;                 //
        tmp = Current[i].startPath;                                               //
        PrintPath(tmp, out);                                                    //
        out << "Here is the Proposed Path now " << endl << endl;                //
        tmp = Proposed[i].startPath;                                              //
        PrintPath(tmp, out);
      }                                                    //
#endif  ///////////////////////////////////////////////////////////////////////

      int accept_K, accept_W, accept_WS, accept_WP, accept_P, same_P, accept_pi;
      int  total_K, total_W, total_WS, total_WP, total_P, total_U, accept_U;
      int totNode, yesNode, samePath;
      same_P = accept_K = accept_W = accept_WS = accept_WP = accept_P =
      accept_pi = total_K = total_W = total_WS = total_WP =total_P = total_U =
      accept_U = totNode = yesNode = samePath = 0;
      out << "Ready to begin the calculation " << endl;
      static int count = 0;
      static double **numEvents, **meanTime;
      static int ***InternalSeq, **InternalSeqCount, **WhichInternalSeq;
      static int *NumInternalSeq = new int[numInternal];
      for(int qq = 0; qq < numInternal; qq++){NumInternalSeq[qq] = 0;}
      static int numMiss, *newOrder, *missingSites, *tmpskip;
      int *AAptr1, *AAptr2;
// CHANGES by JEFF around here
      for(int M = 0; M <= MAXITERATION; M++)
      {
        if((M>=BURNIN)&&(!((M-BURNIN)%SAMPFREQ)))
        {
         if(!count)
         {
            InternalSeq = new int **[numInternal];
            WhichInternalSeq = new int *[numInternal];
            InternalSeqCount = new int *[numInternal];
            for(int qq = 0; qq < numInternal; qq++)
            {
               InternalSeq[qq] = new int *[SAMPLE_PATH];
               WhichInternalSeq[qq] = new int [SAMPLE_PATH];
               InternalSeqCount[qq] = new int [SAMPLE_PATH];
               for(int q = 0; q < SAMPLE_PATH; q++){InternalSeqCount[qq][q] = 0;}
            }
            //cout << "Done initializing the internal sequence count " << endl;
            numEvents = new double *[numBranch];
            meanTime  = new double *[numBranch];
            for(int q = 0; q < numBranch; q++)
            {
               numEvents[q] = new double[nucSeqLen];
               meanTime[q] = new double[nucSeqLen];
               for(int r = 0; r < nucSeqLen; r++)
               {
                  meanTime[q][r] = numEvents[q][r] = 0.0;
               }
            }
            newOrder = new int[AAseqLen];
            for(int q = 0; q < AAseqLen; q++){newOrder[q] = 1;}
            in >> numMiss;    // cout << "This is numMiss " << numMiss << endl;
            if(numMiss > 0)
            {
               missingSites = new int[numMiss];
               tmpskip = new int[numMiss];
               for(int q = 0; q < numMiss; q++)
               {
                  in >> missingSites[q];   tmpskip[q] = missingSites[q];
               }
               int counting = 0;
               for(int q = 0; q < AAseqLen; q++)
               {
                 if(counting == numMiss){newOrder[q] = counting+1;}
                 else
                 {
                   if(q < missingSites[counting]){newOrder[q] = counting+1;}
                   else
                   {
                     if(q == missingSites[counting])
                     {
                       counting++;
                       if(counting < numMiss){missingSites[counting] -= counting;}
                     }
                     q--;
                   }
                 }
               }
            }
            in.close();
         } // END IF COUNT == 0
        }

#if WHERE_ARE_WE
         if(!(M % 10000)){out << "This is iteration " << M << endl;}
         if(!(M % 500)){
           //ofstream outCount("1ALA.round.out");
           ofstream outCount(roundout_filename);
	       outCount << "We are somewhere between iteration " << M;
	       outCount << " and iteration " << (M+500) << endl;
         }
#endif
         if(!M){out << "This is First Param" << endl; FirstParam(PostEstimate,Current,Proposed,FirstOrder,out);out << "This is after First Param "<< endl;}
#if SIM_K///////////////////////////////////////////////////////////////////////
         total_K++;
         for(int zz = 0; zz < 61; zz++){cf[zz]=0.0;}
         codonFrequency(cf, PostEstimate->MCMCnf);
         if(!M)
         {
           if(MULTI_K)
           {
              for(int x = 0; x < numBranch; x++)
              {
                out << "This is Init Kappa parameter: " << PostEstimate->kappa[x];
                out << "  This is deltaMult for kappa : " << DELTA_K <<
                " prior [" << K_MIN << "," << K_MAX << "]"<< endl;
              }
           }
           else
           {
                out << "This is Init Kappa parameter: " << PostEstimate->kappa[0];
                out << "  This is deltaMult for kappa : " << DELTA_K <<
                " prior [" << K_MIN << "," << K_MAX << "]"<< endl;
           }
         }
         NewKappa(PostEstimate,Current,Proposed,FirstOrder,RateNRG,accept_K,cf,seed,out);
//        out << "After kappa " << M << " " << PostEstimate->kappa[0] << endl;
#endif//////////////////////////////////////////////////////////////////////////
#if SIM_W///////////////////////////////////////////////////////////////////////
         total_W++;
         if(!M)
         {
           if(MULTI_W)
           {
              for(int x = 0; x < numBranch; x++)
              {
              out << "This is Init Omega parameter: " << PostEstimate->omega[x];
              out << "  This is deltaMult for Omega : " << DELTA_W <<
              " prior [" << W_MIN << "," << W_MAX << "]"<< endl;
              }
           }
           else
           {
              out << "This is Init Omega parameter: " << PostEstimate->omega[0];
              out << "  This is deltaMult for Omega : " << DELTA_W <<
              " prior [" << W_MIN << "," << W_MAX << "]"<< endl;
           }
         }
         NewOmega(PostEstimate,Current,Proposed,FirstOrder,RateNRG,accept_W,cf,seed,out);
//         out << "After omega " << PostEstimate->omega[0] << endl;
#endif//////////////////////////////////////////////////////////////////////////
#if SIM_U///////////////////////////////////////////////////////////////////////
         total_U +=numBranch;
         if(!M){
           for(int x = 0; x < numBranch; x++)
           {
              out << "This is Init Rate parameter: " << PostEstimate->rate[x];
              out << "  This is deltaMult for rate : " << DELTA_U <<
              " prior [" << U_MIN << "," << U_MAX << "]"<< endl;
            }
         }
     //    out << "Before rate " << endl;
         NewRate(PostEstimate,Current, Proposed, FirstOrder,RateNRG,accept_U,seed,out);
#endif//////////////////////////////////////////////////////////////////////////
#if SIM_PI//////////////////////////////////////////////////////////////////////
         if(!M){out << "We are about to update the pi's right now " << endl;}
         NewPi(PostEstimate,Current,Proposed,FirstOrder, RateNRG, accept_pi,(M%4),seed,out);
#endif//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Put the energy stuff in here????





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#if SIM_W_S/////////////////////////////////////////////////////////////////////
         total_WS++;
         if(!M)
         {
            out << "This is Init Omega solv: " << PostEstimate->solvent << " ";
            out << "This is deltaMult for w_solv : " << DELTA_W_S <<
            " prior [" << W_S_MIN << "," << W_S_MAX << "]"<< endl;
         }
       //  out << "Before solvent " << endl;
         NewOmega_S(PostEstimate,Current,Proposed,FirstOrder,RateNRG,accept_WS,seed,out);
 //        out << "After solvent " << PostEstimate->solvent << endl;
#endif//////////////////////////////////////////////////////////////////////////
#if SIM_W_P/////////////////////////////////////////////////////////////////////
         total_WP++;
         if(!M){
            out << "This is Init Omega pair: " << PostEstimate->pairwise << " ";
            out << "This is deltaMult for w_pair  : " << DELTA_W_P <<
            " prior [" << W_P_MIN << "," << W_P_MAX << "]"<< endl;
         }
     //    out << "Before pairwise " << endl;
         NewOmega_P(PostEstimate,Current,Proposed,FirstOrder,RateNRG,accept_WP,seed,out);
//         out << "After pairwise " << PostEstimate->pairwise << endl;
#endif//////////////////////////////////////////////////////////////////////////
         if((M >= BURNIN)&&(!((M-BURNIN)%SAMPFREQ)))
         {
            // Here we deal with the parameter posterior distribution
            parameterPtr = *paramOut1;
            int nonsyn[3];
            for(int v = 0; v < numBranch; v++)
            {
               nonsyn[v] = 0;
               tmp = Current[v].startPath;
               while(tmp != NULL)
               {
                  if(tmp->amino != tmp->priorAA){nonsyn[v]++;}
                  tmp = tmp->DpathPtr;
               }
            }

            if(MULTI_K)
            {
               for(int v = 0; v < numBranch; v++)
               {
                 *parameterPtr = PostEstimate->kappa[v];
                 parameterPtr++;
               }
            }
            else{*parameterPtr=PostEstimate->kappa[0];parameterPtr++;}
            if(MULTI_W)
            {
               for(int v = 0; v < numBranch; v++)
               {
                 *parameterPtr = PostEstimate->omega[v];
                 parameterPtr++;
               }
            }
            else{*parameterPtr=PostEstimate->omega[0];parameterPtr++;}

            *parameterPtr=PostEstimate->solvent;parameterPtr++;
            *parameterPtr=PostEstimate->pairwise;parameterPtr++;
            *parameterPtr=PostEstimate->MCMCnf[0];parameterPtr++;
            *parameterPtr=PostEstimate->MCMCnf[1];parameterPtr++;
            *parameterPtr=PostEstimate->MCMCnf[2];parameterPtr++;
            *parameterPtr=PostEstimate->MCMCnf[3]; parameterPtr++;
           // cout << "We are done with the parameters " <<*parameterPtr << endl;

            for(int v = 0; v < numBranch; v++)
            {
               *parameterPtr=PostEstimate->rate[v];parameterPtr++;
               *parameterPtr=nonsyn[v];parameterPtr++;
               *parameterPtr=Current[v].totNumSub;parameterPtr++;
               *parameterPtr=PostEstimate->branchLen[v]*Current[0].invlen;parameterPtr++;
               *parameterPtr=PostEstimate->yangBranch[v];parameterPtr++;
            }
            paramOut1++;

            // Now we are going to deal with the Path Posterior!!
            for(int q = 0; q < numBranch; q++)
            {
              tmp = Current[q].startPath;
              int tSub;
              PathSubsCount[q][count] = tSub = Current[q].totNumSub;
              PathStorage[q][count] = new double *[tSub];
              for(int x = 0; x < tSub; x++)
              {
                 PathStorage[q][count][x] = new double[18];
              }

              double **pathPtr = PathStorage[q][count];
              PathPostFull(PostEstimate,tmp,FirstOrder,count,q,numEvents[q],
                     meanTime[q], newOrder, pathPtr);
              if(!(count%10))
              {
                 NRGposterior(PostEstimate,Current[q], FirstOrder,q,count,
                        RateNRG[q],outP2);
              }
            }
//            out << "After priniting everything out to the files " << endl;
            #if ANCESTRAL_NUC
            for(int internal = 0; internal < numInternal; internal++)
            {
               if(!NumInternalSeq[internal])
               {
                  InternalSeqCount[internal][0]++;
                  WhichInternalSeq[internal][count] = count;
                  InternalSeq[internal][0] = new int[nucSeqLen];
                  AAptr1 = InternalSeq[internal][0];
                  AAptr2 = Current[internal+numSeq].seq;
                  for(int u = 0; u < nucSeqLen; u++)
                  {
                     *AAptr1 = *AAptr2; AAptr1++; AAptr2++;
                  }
                  NumInternalSeq[internal]++;
               }
               else
               {
                  int same = 1;
                  for(int u = 0; u < NumInternalSeq[internal]; u++)
                  {
                     same = 1;
                     AAptr1 = InternalSeq[internal][u];
                     AAptr2 = Current[internal+numSeq].seq;
                     for(int v = 0; v < nucSeqLen; v++)
                     {
                        if(*AAptr1 != *AAptr2){same = 0; v = nucSeqLen+1;}
                        AAptr1++; AAptr2++;
                     }
                     if(same)
                     {
                        InternalSeqCount[internal][u]++;
                        WhichInternalSeq[internal][count] = u;
                        u = NumInternalSeq[internal]+1;
                     }
                  }
                  if(!same)
                  {
                     InternalSeqCount[internal][NumInternalSeq[internal]]++;
                     WhichInternalSeq[internal][count]=NumInternalSeq[internal];
                     InternalSeq[internal][NumInternalSeq[internal]] = new int[nucSeqLen];
                     AAptr1 = InternalSeq[internal][NumInternalSeq[internal]];
                     AAptr2 = Current[internal+numSeq].seq;
                     for(int u = 0; u < nucSeqLen; u++)
                     {
                        *AAptr1 = *AAptr2; AAptr1++; AAptr2++;
                     }
                     NumInternalSeq[internal]++;
                  }
               } // END else statement
            } // END looping through the internal nodes
            #endif
            #if (ANCESTRAL_NUC == 0)
            for(int internal = 0; internal < numInternal; internal ++)
            {
              if(!NumInternalSeq[internal])
              {
                InternalSeqCount[internal][0]++;
                WhichInternalSeq[internal][count] = count;
                InternalSeq[internal][0] = new int[AAseqLen];
                AAptr1 = InternalSeq[internal][0];
                AAptr2 = Current[internal+numSeq].AA_seq;
                for(int u = 0; u < AAseqLen; u++)
                {
                  *AAptr1 = *AAptr2; AAptr1++; AAptr2++;
                }
                NumInternalSeq[internal]++;
              }
              else
              {
                 int same = 1;
                 for(int u = 0; u < NumInternalSeq[internal]; u++)
                 {
                   same = 1;
                   AAptr1 = InternalSeq[internal][u];
                   AAptr2 = Current[internal+numSeq].AA_seq;
                   for(int v = 0; v < AAseqLen; v++)
                   {
                     if(*AAptr1 != *AAptr2){same = 0; v = AAseqLen+1;}
                     AAptr1++; AAptr2++;
                   }
                   if(same)
                   {
                     InternalSeqCount[internal][u]++;
                     WhichInternalSeq[internal][count] = u;
                     u = NumInternalSeq[internal]+1;
                   }
                 }
                 if(!same)
                 {
                   InternalSeqCount[internal][NumInternalSeq[internal]]++;
                   WhichInternalSeq[internal][count] = NumInternalSeq[internal];
                   InternalSeq[internal][NumInternalSeq[internal]] = new int[AAseqLen];
                   AAptr1 = InternalSeq[internal][NumInternalSeq[internal]];
                   AAptr2 = Current[internal+numSeq].AA_seq;
                   for(int u = 0; u < AAseqLen; u++)
                   {
                      *AAptr1 = *AAptr2; AAptr1++; AAptr2++;
                   }
                   NumInternalSeq[internal]++;
                 }
              } // END else statement
            } // END looping through the internal nodes
            #endif
            count++;
            if(count == SAMPLE_PATH)
            {
              // This is the output of the parameter posterior distribution
              paramOut1 = paramOut;
//              cout << "We are about to start the Node Update " << endl;

              for(int ww = 0; ww < SAMPLE_PATH; ww++)
              {
//                 cout << "This is the iteration number " << ww  << endl;
                 parameterPtr = *paramOut1;
                 for(int qqq = 0; qqq < paramnum; qqq++)
                 {
                    out1 << *parameterPtr << " ";    parameterPtr++;
                 }
                 out1 << endl; paramOut1++;
              }
              for(int ww = 0; ww < SAMPLE_PATH; ww++){delete[] paramOut[ww];}
              delete[] paramOut;

              out << "We are done with the parameter out" << endl;
              double ***branchPathPtr;
              for(int q = 0; q < numBranch; q++)
              {
                 branchPathPtr = PathStorage[q];
                 PathOutput(branchPathPtr,PathSubsCount[q],outP1);
                 out << "Done with branch paths " << q << endl;
              }

              int *pathCountHolder;
              for(int ww = 0; ww < numBranch; ww++)
              {
                 pathCountHolder = PathSubsCount[ww];
                 for(int qqq = 0; qqq < SAMPLE_PATH; qqq++)
                 {
                    for(int qwqw = 0; qwqw < *pathCountHolder; qwqw++)
                    {
                       delete[] PathStorage[ww][qqq][qwqw];
                    }
                    delete[] PathStorage[ww][qqq];
                    pathCountHolder++;
                 }
                 delete[] PathStorage[ww];
              }
              delete[] PathStorage;
              for(int qwqw =0; qwqw < numBranch;qwqw++)
              {
                 delete[] PathSubsCount[qwqw];
              }
              delete[] PathSubsCount;


              //Complete all functions and erase all of the memory we borrowed!
              /////////////////////////////////////////////////////////////////
              // This is the code to print out the distribution of the internal
              // node sequence ////////////////////////////////////////////////
              /////////////////////////////////////////////////////////////////
              if(numMiss > 0)
              {
                 for(int q = 0; q < numMiss; q++){missingSites[q] = tmpskip[q];}
              }


              #if ANCESTRAL_NUC
              for(int internal = 0; internal < numInternal; internal ++)
              {
                  out2 << "Here is the distribution of the internal sequence "<<
                        internal+numSeq << endl << endl<<endl<<endl<<endl<<endl<<endl<<endl;

                for(int u = 0; u < NumInternalSeq[internal]; u++)
                {
                  out2 << "---------- Internal node " << internal+numSeq <<
                  " Sequence number " << u << " ---------- " <<
                  " with frequency " << InternalSeqCount[internal][u] <<
                  "  ----------" << endl;
                  AAptr1 = InternalSeq[internal][u];
                  if(numMiss == 0)
                  {
                    for(int v = 0; v < nucSeqLen; v++)
                    {
                      out2 << IntToNuc(*AAptr1); AAptr1++;
                    }
                  }
                  else
                  {
                    int skipsite = 0;
                    for(int v = 0; v < numMiss; v++){tmpskip[v] = 3*missingSites[v];}
                    for(int q = 0; q < nucSeqLen; q++)
                    {
                      if(skipsite == numMiss){out2 << IntToNuc(*AAptr1); AAptr1++;}
                      else
                      {
                        if(q < tmpskip[skipsite]){out2 << IntToNuc(*AAptr1); AAptr1++;}
                        else
                        {
                          if(q == tmpskip[skipsite])
                          {
                            out2 << "---";
                            skipsite++;
                            if(skipsite < numMiss){tmpskip[skipsite] -= 3*skipsite;}
                          }
                          q--;
                        }
                      }
                    } // End printing out an internal node sequence
                  } // End else statement
                  out2 << endl << endl;
                } // End printing our all internal node sequences
              }// END looping through the internal nodes
              #endif
              #if (ANCESTRAL_NUC == 0)
              out2 << "Here is the distribution of the internal sequence "<<
              internal+numSeq<<endl<<endl<<endl<<endl<<endl<<endl<<endl<<endl;

              for(int internal = 0; internal < numInternal; internal ++)
              {
                for(int u = 0; u < NumInternalSeq[internal]; u++)
                {
                  out2 << "---------- Internal node " << internal+numSeq <<
                  " Sequence number " << u << " ---------- " <<
                  " with frequency " << InternalSeqCount[internal][u] <<
                  "  ----------" << endl;
                  AAptr1 = InternalSeq[internal][u];
                  if(numMiss == 0)
                  {
                    for(int v = 0; v < AAseqLen; v++)
                    {
                       out2 << Int2AA(*AAptr1); AAptr1++;
                    }
                  }
                  else
                  {
                    int skipsite = 0;
                    for(int v = 0; v < numMiss; v++){tmpskip[v] = missingSites[v];}
                    for(int q = 0; q < AAseqLen; q++)
                    {
                      if(skipsite == numMiss){out2 << Int2AA(*AAptr1); AAptr1++;}
                      else
                      {
                        if(q < tmpskip[skipsite]){out2 << Int2AA(*AAptr1); AAptr1++;}
                        else
                        {
                          if(q == tmpskip[skipsite])
                          {
                            out2 << "-";
                            skipsite++;
                            if(skipsite < numMiss){tmpskip[skipsite] -= skipsite;}
                          }
                          q--;
                        }
                      }
                    } // End printing out an internal node sequence
                  } // End else statement
                  out2 << endl << endl;
                } // End printing our all internal node sequences
              }// END looping through the internal nodes
              #endif
              for(int internal = 0; internal < numInternal; internal ++)
              {
              out2 << "******************************************************* " << endl;
              out2 << "This is the same info as above, but in a readable form " << endl;
              out2 << "To be read as Sampled Sequence # has frequence ## " << endl;
              out2 << "******************************************************* " << endl;
                for(int u = 0; u < NumInternalSeq[internal]; u++)
                {
                  out2<<internal<<"\t"<<u<<"\t"<<InternalSeqCount[internal][u] <<endl;
                }
                out2 << endl << endl << endl <<

              "******************************************************* " << endl <<
              "Here we show which paths have the same ancestral sequence" << endl <<
              "To be read as Ancestral Sequence # is identical to  ## " << endl <<
              "******************************************************* " << endl << endl;
                for(int u = 0; u < SAMPLE_PATH; u++)
                {
                  out2<<internal<<"\t"<<u<<"\t"<<WhichInternalSeq[internal][u]<<endl;
                }
                for(int u = 0; u < NumInternalSeq[internal]; u++)
                {
                   delete[] InternalSeq[internal][u];
                }
                delete[] WhichInternalSeq[internal];
                delete[] InternalSeq[internal];
                delete[] InternalSeqCount[internal];
              } // END looping through the internal nodes
              delete[] WhichInternalSeq;  delete[] InternalSeq;
              delete[] InternalSeqCount;  delete[] missingSites;
              delete[] newOrder;          delete[] tmpskip;
              delete[] NumInternalSeq;
              /////////////////////////////////////////////////////////////////
              // This is where we plot the mean time of events in each col
              /////////////////////////////////////////////////////////////////
               out2 << endl << endl <<
               "Here is the mean time of the sub events in each column" << endl <<
               "col#    mean b1    mean b2    mean b3" << endl << endl;
           for(int internal = 0; internal < numInternal; internal ++)
           {

               for(int r = 0; r < nucSeqLen; r++)
               {
                  out2<<internal<<"\t"<<r<<"\t";
                  for(int q = 0; q < numBranch; q++)
                  {
                     if(numEvents[q][r] > 0.0)
                     {
                         out2 << (meanTime[q][r] / numEvents[q][r]) << "\t";
                     }
                     else{out2 << 0.0 << "\t";}
                  }
                  out2 << endl;
               }
               out2 << endl << endl << endl << endl << endl;
               /////////////////////////////////////////////////////////////////
               // This is where we plot the average number of events in each col
               /////////////////////////////////////////////////////////////////
               out2 << "Here is the absolute and average number of sub events per column" << endl;
               out2 << "col#    avg b1           avg b2           avg b3" << endl << endl;
               for(int r = 0; r < nucSeqLen; r++)
               {
                  out2<<internal<<"\t"<<r<<"\t";
                  for(int q = 0; q < numBranch; q++)
                  {
                      out2<<(numEvents[q][r])<<"\t"<<(numEvents[q][r]/count)<<"\t";

                  }
                  out2 << endl;
               }
            }
               for(int q = 0; q < numBranch; q++)
               {
                  delete[] numEvents[q];
                  delete[] meanTime[q];
               }
               delete[] numEvents;
               delete[] meanTime;
            }
         }


         if(!M){out << "We are about to update the path right now " << endl;}
//         out << "Before path " << endl;
         NewPath(PostEstimate, Current, Proposed, FirstOrder, seed, total_P,
                        accept_P, same_P, out);
//         out << "After path ";
//        for(int BMS = 0; BMS < numBranch; BMS++){out << Current[BMS].pathHood << " ";}
 //        out << endl;
////////////////////////////////////////////////////////////////////////////////
// Here is where we put the output! ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//         cout << "We are about to start the Node Update " << M << endl;
//         out << "Before node " << M;
         if(!(M%NODE_UPDATE) && (numSeq > 2))
         {
            //ofstream outJunk("junk.out");
            //outJunk << "Before NewNodeSeq " << endl;
            if(!M){out << "We are about to update the node right now " << endl;}
            NewNodeSeq(PostEstimate, Current, Proposed,FirstOrder, seed,
                        totNode, yesNode, samePath,out);
         }
//        out << "After node " ;
//        for(int BMS = 0; BMS < numBranch; BMS++){out << Current[BMS].pathHood << " ";}
//        out << endl;

////////////////////////////////////////////////////////////////////////////////
#if LOUD////////////////////////////////////////////////////////////////////////
         out1 << "This is the current path after round " << M << endl << endl;//
         tmp = Current[0].startPath;                                          //
         PrintPath(tmp, out1);                                                //
         out1 << M << endl;                                                   //
         out1<< "This is the proposed path after round " << M << endl << endl;//
         tmp = Proposed[0].startPath;                                         //
         PrintPath(tmp, out1);                                                //
#endif//////////////////////////////////////////////////////////////////////////
      }  // This ends the HUGE loop through the max number of interations///////
////////////////////////////////////////////////////////////////////////////////
// This is where we output all of the acceptance / rejection data for MCMC sim//
      double JT,LT;
      JT = total_K; LT = accept_K;
      if(MULTI_K){total_K *= numBranch;}
      out << "Total number of kappa = " << total_K;
      if(JT < 1){out << ", number accept = " << accept_K << "  =  0"<< endl;}     //
      else{out << ", number accept = " << accept_K << "  =  "<< LT/JT << endl;}//
      JT = total_P; LT = accept_P;
      out << "Total number of Path = " << total_P;
      if(JT < 1){out << ", number accept = " << accept_P << "  =  0"<< endl;}     //
      else{out << ", number accept = " << accept_P << "  =  "<< LT/JT << endl;}//
      out << "Total number of same paths = " << same_P << endl;               //
      JT = total_W; LT = accept_W;
      if(MULTI_W){total_W *= numBranch;}
      out << "Total number of omega = " << total_W;
      if(JT < 1){out << ", number accept = " << accept_W << "  =  0"<< endl;}     //
      else{out << ", number accept = " << accept_W<< "  =  "<< LT/JT << endl;}//
      JT = total_WS; LT = accept_WS;
      out << "Total number of w_solv = " << total_WS;
      if(JT < 1){out << ", number accept = " << accept_WS << "  =  0"<< endl;}     //
      else{out << ", number accept = " << accept_WS << "  =  "<< LT/JT << endl;}//                        //
      JT = total_WP; LT = accept_WP;
      out << "Total number of w_pair = " << total_WP;
      if(JT < 1){out << ", number accept = " << accept_WP << "  =  0"<< endl;}     //
      else{out << ", number accept = " << accept_WP << "  =  "<< LT/JT << endl;}//                     //
      JT = total_U; LT = accept_U;
      out << "Total number of U = " << total_U;
      if(JT < 1){out << ", number accept = " << accept_U << "  =  0"<< endl;}     //
      else{out << ", number accept = " << accept_U << "  =  "<< LT/JT << endl;}//                            //
      JT = totNode; LT = yesNode;
      out << "Total number of Node Updates = " << totNode;
      if(JT < 1){out << ", number accept = " << yesNode << "  =  0"<< endl;}     //
      else{out << ", number accept = " << yesNode << "  =  "<< LT/JT << endl;}//
                           //
////////////////////////////////////////////////////////////////////////////////
      if(trial < (ITERATE - 1))
      {
         DeleteSome(PostEstimate,Current,Proposed,seed,out);
      }
      out1 << endl << endl << endl;
   } // This ends the MASSIVE loop through all of the Sequences/////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

   //for (int q = 0; q < numSeq; q++)
   for (int q = 0; q < numBranch; q++)
   {
      for(int r = 0; r < doubleNucLen; r++)
      {
         for(int s = 0; s < AAseqLen; s++){delete[] RateNRG[q][r][s];}
         delete[] RateNRG[q][r];
      }
      delete[] RateNRG[q];
   }
   delete[] RateNRG;

   outP1.close();
   outP2.close();
//   outP3.close();
//   outP4.close();
//   outP5.close();
//   outP6.close();
   out1.close();
   out2.close();
   DeleteAll(PostEstimate, FirstOrder,Current,Proposed,out);
   out << "This is DrEVOL over and out! " << endl;
   out.close();
}
//******************************************************************************
/* From Joe Felsenstein's PHYLIP package ...
   check for "needed by rnd" for other lines of code this routine needs.

   Performance-optimized by R.Hopperger (office@hope-soft.co.at)
*/
double rnd(long* seed)
{
  newseed0 = 13 * seed[0];
  //for (j=0; j <= 4; j++) {newseed[j+1] += newseed[j] / 64;   newseed[j] &= 63;}
  newseed1 = newseed0>>6;   newseed0 &= 63;
  newseed2 = newseed1>>6;   newseed1 &= 63;
  newseed3 = newseed2>>6;   newseed2 &= 63;
  newseed4 = newseed3>>6;   newseed3 &= 63;
  newseed5 = newseed4>>6;   newseed4 &= 63;


  newseed1 += (13 * seed[1]) + (24 * seed[0]);
  //for (j=1; j <= 4; j++) {newseed[j+1] += newseed[j] / 64;   newseed[j] &= 63;}
  newseed2 += newseed1>>6;   newseed1 &= 63;
  newseed3 += newseed2>>6;   newseed2 &= 63;
  newseed4 += newseed3>>6;   newseed3 &= 63;
  newseed5 += newseed4>>6;   newseed4 &= 63;


  newseed2 += (13 * seed[2]) + (24 * seed[1]) + (22 * seed[0]);
  //for (j=2; j <= 4; j++) {newseed[j+1] += newseed[j] / 64;   newseed[j] &= 63;}
  newseed3 += newseed2>>6;   newseed2 &= 63;
  newseed4 += newseed3>>6;   newseed3 &= 63;
  newseed5 += newseed4>>6;   newseed4 &= 63;


  newseed3 += (13 * seed[3]) + (24 * seed[2]) + (22 * seed[1]) + (6  * seed[0]);
  //for (j=3; j <= 4; j++) {newseed[j+1] += newseed[j] / 64;   newseed[j] &= 63;}
  newseed4 += newseed3>>6;   newseed3 &= 63;
  newseed5 += newseed4>>6;   newseed4 &= 63;


  newseed4 += (13 * seed[4]) + (24 * seed[3]) + (22 * seed[2]) + (6  * seed[1]);
  //for (j=4; j <= 4; j++) {newseed[j+1] += newseed[j] / 64;   newseed[j] &= 63;}
  newseed5 += newseed4>>6;   newseed4 &= 63;


  newseed5 += (13 * seed[5]) + (24 * seed[4]) + (22 * seed[3]) + (6  * seed[2]);
  //for (j=5; j <= 4; j++) {newseed[j+1] += newseed[j] / 64;   newseed[j] &= 63;}

  return ((((( (seed[0]=newseed0) / 64.0 + (seed[1]=newseed1) ) / 64.0 +
  (seed[2]=newseed2) ) / 64.0 + (seed[3]=newseed3) ) / 64.0 +
  (seed[4]=newseed4) ) / 64.0 + (seed[5]=newseed5 & 3) ) / 4.0;
}
//******************************************************************************
/* From Joe Felsenstein's PHYLIP package ...
   check for "needed by rnd" for other lines of code this routine needs
double rnd(long* seed)
{
  /* random number generator -- slow but machine independent
  long i, j, k, sum;
  longer mult, newseed;
  double x;
  mult[0] = 13;
  mult[1] = 24;
  mult[2] = 22;
  mult[3] = 6;
  for (i = 0; i <= 5; i++){newseed[i] = 0;} // initialize the components
  for (i = 0; i <= 5; i++)
  {   sum = newseed[i];
      k = i;
      if (i > 3)
      k = 3;
      for (j = 0; j <= k; j++){sum += mult[j] * seed[i - j];}
      newseed[i] = sum;
      for (j = i; j <= 4; j++) {
          newseed[j + 1] += newseed[j] / 64;
          newseed[j] &= 63;
      }
  }
  memcpy(seed, newseed, sizeof(longer));
  seed[5] &= 3;
  x = 0.0;
  for (i = 0; i <= 5; i++)
    x = x / 64.0 + seed[i];
  x /= 4.0;
  return x;
}    rnd */
//******************************************************************************
double Maxof3(double a, double b, double c)
{
   if(a < b){if(b < c){return c;}else{return b;}}
   else{if(a < c){return c;}else{return a;}}
}
//******************************************************************************
void initPath(Path *c, Path *p, int seqlen, int *seq)
// Initializes the Path Struct
{
        int *firstSeq = seq;
        for(int i = 0; i < seqlen; i++)
        {
           c[i].nuc    = p[i].nuc    = *firstSeq;
           c[i].numsub = p[i].numsub = 0;
           c[i].firstColSub = NULL;
           p[i].firstColSub = NULL;
           firstSeq++;
        }
} // End initPath Subroutine
//******************************************************************************
void initInteraction(Interaction *I) // Initializes the Interaction structure
{
        I->energyAcc   = NULL;
        I->nAccess     = NULL;
        I->neighList   = NULL;
        I->energy      = NULL;
}  // End initInteraction Subroutine
//******************************************************************************
void initBayes(Bayes *pE)
{
        pE->solvent = pE->pairwise = pE->solvNRG =
        pE->pairNRG = pE->nrgConst = pE->maxValue = pE->tempsolvNRG =
        pE->temppairNRG=pE->tempNRGconst=pE->tempMAXvalue=pE->tempParameter = 0.0;
        pE->A_gridPt = pE->tempA_gridPt = pE->C_gridPt = pE->tempC_gridPt =
        pE->G_gridPt = pE->tempG_gridPt = pE->S_gridPt = pE->tempS_gridPt =
        pE->P_gridPt = pE->tempP_gridPt = pE->numBranch = 0;
        pE->kappa = NULL;
        pE->omega = NULL;
        pE->rate = NULL;
        pE->tempRate = NULL;
        pE->yangBranch = NULL;
        pE->branchLen = NULL;
        pE->MCMCnf = NULL;
        pE->tempMCMCnf = NULL;
}
//******************************************************************************
void initI(Information &I) // Initializes the information matrix
{
        I.seq         = NULL;
        I.AA_seq      = NULL;
        I.neighAcc    = NULL;
        I.startPath   = NULL;
        I.totNumSub = 0;
        I.AAseqNRG = I.pathHood = I.avgRate = 0.0;
}  // End initI Subroutine
//******************************************************************************
double minof2(double x, double y)  // Returns the min of two doubles
{
        if(x < y) {return x;}  else {return y;}
}  // End minof2 Subroutine
//******************************************************************************
void nucCompare()
{
   /* here 	0 = same               This matrix gives a quick reference
   //		1 = Transition         to see how related nucleotides
   //		2 = Transversion       are in an alignment
   */
   for(int i = 0; i < 4; i++)
   {
      for(int j = 0; j < 4; j++)
      {
         if(i == j)                {COMPARE[i][j] = 0;}
	 else if((i + j) % 2 == 0) {COMPARE[i][j] = 1;}
	 else                      {COMPARE[i][j] = 2;}
      }
   }
}  // End nucCompare Subroutine
//******************************************************************************
char IntToNuc(int m)
{
   char nuc;
   switch(m){
      case 0: nuc = 'A'; break;
      case 1: nuc = 'C'; break;
      case 2: nuc = 'G'; break;
      case 3: nuc = 'T'; break;
   }
   return nuc;
}
//******************************************************************************
int NucToInt(char m)  //Converts nucs to integer representation
{                                 // Also calculates freqs of nucleotides
	int num;
	switch ( m ){
        case 'A':  num = 0; F84_freq[0]+=1.0; break;
        case 'a':  num = 0; F84_freq[0]+=1.0; break;
        case 'C':  num = 1; F84_freq[1]+=1.0; break;
        case 'c':  num = 1; F84_freq[1]+=1.0; break;
        case 'G':  num = 2; F84_freq[2]+=1.0; break;
        case 'g':  num = 2; F84_freq[2]+=1.0; break;
        case 'T':  num = 3; F84_freq[3]+=1.0; break;
        case 't':  num = 3; F84_freq[3]+=1.0; break;
        case 'U':  num = 3; F84_freq[3]+=1.0; break;
        case 'u':  num = 3; F84_freq[3]+=1.0; break;
	   default: cerr << m << "is not a valid character " << endl;
        }
        return num;
}  // End NucToInt Subroutine
//******************************************************************************
int Nuc2AA(int *c)
{// 0 = ALA     5 = GLN     10 = LEU     15 = SER             //
 // 1 = ARG     6 = GLU     11 = LYS     16 = THR             //
 // 2 = ASN     7 = GLY     12 = MET     17 = TRP   20 = STOP //
 // 3 = ASP     8 = HIS     13 = PHE     18 = TYR             //
 // 4 = CYS     9 = ILE     14 = PRO     19 = VAL             //
 ///////////////////////////////////////////////////////////////
   if(c[0] == 1)
   {
      if(c[1] == 1) {return 14;}           // PRO
      else if(c[1] == 2) {return 1;}       // ARG
      else if(c[1] == 3) {return 10;}      // LEU
      else if((c[2]%2) == 0) {return 5;}   // GLN
      else {return 8;}                     // HIS
   }
   else if(c[0] == 2)
   {
      if(c[1] == 1) {return 0;}            // ALA
      else if(c[1] == 2) {return 7;}       // GLY
      else if(c[1] == 3) {return 19;}      // VAL
      else if((c[2]%2) == 0) {return 6;}   // GLU
      else {return 3;}                     // ASP
   }
   else if(c[0] == 3)
   {
      if(c[1] == 1) {return 15;}           // SER
      else if(c[1] == 3)
      {
         if((c[2]%2) == 0) {return 10;}    // LEU
         else {return 13;}                 // PHE
      }
      else if(c[1] == 0)
      {
         if((c[2]%2) == 1) {return 18;}    // TYR
         else {return 20;}                 // STOP
      }
      else if((c[2]%2) == 1) {return 4;}   // CYS
      else if(c[2] == 0) {return 20;}      // STOP
      else {return 17;}                    // TRP
   }
   else
   {
      if(c[1] == 1) {return 16;}           // THR
      else if(c[1] == 0)
      {
         if((c[2]%2) == 0) {return 11;}    // LYS
         else {return 2;}                  // ASN
      }
      else if(c[1] == 2)
      {
         if((c[2]%2) == 0) {return 1;}     // ARG
         else {return 15;}                 // SER
      }
      else if((c[2]%2) == 1) {return 9;}   // ILE
      else if(c[2] == 0) {return 9;}       // ILE
      else {return 12;}                    // MET
   }
}  // End Nuc2AA Subroutine
//******************************************************************************
char Int2AA(int c)
{// 0 = ALA     5 = GLN     10 = LEU     15 = SER             //
 // 1 = ARG     6 = GLU     11 = LYS     16 = THR             //
 // 2 = ASN     7 = GLY     12 = MET     17 = TRP   20 = STOP //
 // 3 = ASP     8 = HIS     13 = PHE     18 = TYR             //
 // 4 = CYS     9 = ILE     14 = PRO     19 = VAL             //
 ///////////////////////////////////////////////////////////////
   char AA;
   switch(c)
   {
      case 0: AA = 'A'; break;
      case 1: AA = 'R'; break;
      case 2: AA = 'N'; break;
      case 3: AA = 'D'; break;
      case 4: AA = 'C'; break;
      case 5: AA = 'Q'; break;
      case 6: AA = 'E'; break;
      case 7: AA = 'G'; break;
      case 8: AA = 'H'; break;
      case 9: AA = 'I'; break;
      case 10: AA = 'L'; break;
      case 11: AA = 'K'; break;
      case 12: AA = 'M'; break;
      case 13: AA = 'F'; break;
      case 14: AA = 'P'; break;
      case 15: AA = 'S'; break;
      case 16: AA = 'T'; break;
      case 17: AA = 'W'; break;
      case 18: AA = 'Y'; break;
      case 19: AA = 'V'; break;
   }
   return AA;
}  // End Int2AA Subroutine
//******************************************************************************
void PupkoNucChoice(Bayes *pE, Information *C, Information *P, long *seed)
{
   double rnd(long* seed);
   int numNodes, numSeq, numInternal, rootChoice, seqlen, tempNuc;
   int nodesMinus1, nodesMinus2, codon[3], Ci_numChild;
// CHANGE made by JEFF
   double tempProb[4], totalProb, targetProb;
//+ Sang Chul : Jiaye's Revised Pupko
   double sum;
//- Sang Chul : Jiaye's Revised Pupko

   seqlen   = C[0].len;
   numNodes = C[0].numNodes;
   numSeq   = C[0].numSeq;
   numInternal = numNodes - numSeq;
   nodesMinus1 = numNodes - 1;
   nodesMinus2 = numNodes - 2;

   int **tempSeq = new int *[numInternal];
   for(int i = 0; i < numInternal; i++){tempSeq[i] = new int[seqlen];}
   int *nucChoice = new int[numNodes];
   double **probHolder = new double *[numNodes];
   int **tempChoice = new int *[numNodes];
   for(int i = 0; i < numNodes; i++)
   {
      tempChoice[i] = new int[4];
      probHolder[i] = new double[4];
   }

   for(int nuc = 0; nuc < seqlen; nuc++)
   {
//      cout << "This is nuc " << nuc << endl;
      for(int i = 0; i < numNodes; i++)
      {
         Ci_numChild = C[i].numChild;
         if(i < numSeq)// it is a tip node
         {
            tempNuc = C[i].seq[nuc];
            for(int j = 0; j < 4; j++){probHolder[i][j] = F84_PROB[j][tempNuc];}
         }
         else
         {
            if(i < nodesMinus1) //if not root node
            {
               for(int j = 0; j < 4; j++)
               {
                  totalProb = 0.0;
                  for(int k = 0; k < 4; k++)
                  {
                     tempProb[k] = F84_PROB[j][k];
                     //+ Sang Chul : Jiaye's Revised Pupko
                     for(int l = 0; l < Ci_numChild; l++)
                     {
                        sum = 0;
                        for(int m = 0; m < 4; m++) {
                           sum += probHolder[C[i].child[l]][m];
                        }
                        tempProb[k] *= sum;
                     }
                     //- Sang Chul : Jiaye's Revised Pupko
                     //+ Sang Chul : Old version of Pupko
                     //for(int l = 0; l < Ci_numChild; l++)
                     //{
                     //   tempProb[k] *= probHolder[C[i].child[l]][k];
                     //}
                     //- Sang Chul : Old version of Pupko
                     totalProb +=tempProb[k];
                  }
                  targetProb = rnd(seed)*totalProb;
                  if(targetProb < tempProb[0])
                  {
                     probHolder[i][j] = tempProb[0];
                     tempChoice[i][j] = 0;
                  }
                  else if(targetProb < (tempProb[0]+tempProb[1]))
                  {
                     probHolder[i][j] = tempProb[1];
                     tempChoice[i][j] = 1;
                  }
                  else if(targetProb < (tempProb[0]+tempProb[1]+tempProb[2]))
                  {
                     probHolder[i][j] = tempProb[2];
                     tempChoice[i][j] = 2;
                  }
                  else
                  {
                     probHolder[i][j] = tempProb[3];
                     tempChoice[i][j] = 3;
                  }
               }
            }
            else // We are now working with the root node
            {
               totalProb = 0.0;
               for(int j = 0; j < 4; j++)
               {
                  tempProb[j] = F84_freq[j];
                  //+ Sang Chul : Jiaye's Revised Pupko
                  for(int k = 0; k < Ci_numChild; k++)
                  {
                     sum = 0;
                     for(int m = 0; m < 4; m++) {
                        sum += probHolder[C[i].child[k]][m];
                     }
                     tempProb[j] *= sum;
                  }
                  //- Sang Chul : Jiaye's Revised Pupko
                  //for(int k = 0; k < Ci_numChild; k++)
                  //{
                  //   tempProb[j] *= probHolder[C[i].child[k]][j];
                  //}
                  //- Sang Chul : Old version of Pupko
                  totalProb +=tempProb[j];
               }
               targetProb = rnd(seed)*totalProb;
               if(targetProb < tempProb[0]){rootChoice = 0;}
               else if(targetProb < (tempProb[0]+tempProb[1])){rootChoice = 1;}
               else if(targetProb < (tempProb[0]+tempProb[1]+tempProb[2]))
               {
                  rootChoice = 2;
               }
               else{rootChoice = 3;}
            }
         }
      }
      // We now look at the second part of the Pupko Algorithm
//      cout << "At the second part " << endl;
      tempSeq[nodesMinus1-numSeq][nuc] = nucChoice[numNodes-1] = rootChoice;
      for(int j = nodesMinus2; j >= numSeq; j--)
      {
         int j_seq = j - numSeq;
         tempSeq[j_seq][nuc]=nucChoice[j]=tempChoice[j][nucChoice[C[j].parent]];
      }

      int good, seqCount;
// CHANGE by JEFF      int good, seqCount, stop;
      seqCount = 0; good = 1;
      if((nuc%3) == 2)
      {
//         cout << "After the second part " << numInternal << endl;
         while(seqCount < numInternal)
         {
            if(good)
            {
               codon[0] = tempSeq[seqCount][nuc-2];
               codon[1] = tempSeq[seqCount][nuc-1];
               codon[2] = tempSeq[seqCount][nuc];
               if(Nuc2AATable[codon[0]][codon[1]][codon[2]]==20)
               {
                  nuc -= 2; good = 0; seqCount = numInternal + 10;
               }
               seqCount++;
            }
//            cout << "This is good " << good << " and seqCount " << seqCount << endl;
         }
      }
   }
   /* Now we have to fill in the internal node sequences properly
        to all the nodes in the Information structure.
   */
   int *Pseq, *Cseq, *Tseq;
   for(int i = 0; i < numInternal; i++)
   {
      int true_i = i + numSeq;
      Pseq = P[true_i].seq;  Cseq = C[true_i].seq; Tseq = tempSeq[i];
      int *PaaSeq = P[true_i].AA_seq;
      int *CaaSeq = C[true_i].AA_seq;
      for(int nuc = 0; nuc < seqlen; nuc++)
      {
         codon[nuc%3] = *Pseq = *Cseq = *Tseq; Pseq++; Cseq++; Tseq++;
         if((nuc%3) == 2)
         {  // This will be the ancestral amino acid sequence as well!
            *PaaSeq=*CaaSeq=Nuc2AATable[codon[0]][codon[1]][codon[2]];
            PaaSeq++; CaaSeq++;
         }
      }
   }
   for(int i = 0; i < numInternal; i++){delete[] tempSeq[i];}
   for(int i = 0; i < numNodes; i++)
   {
      delete[] tempChoice[i];
      delete[] probHolder[i];
   }
   delete[] probHolder;
   delete[] tempSeq;
   delete[] nucChoice;
   delete[] tempChoice;

}  // End PupkoNucChoice Algorithm
//******************************************************************************
// We now begin the tree toplology part of the program
//******************************************************************************
int TaxaToInt(char z, char y)
{
	int t;

	switch (z){
		case '1': t = 1; break;
		case '2': t = 2; break;
		case '3': t = 3; break;
		case '4': t = 4; break;
		case '5': t = 5; break;
		case '6': t = 6; break;
		case '7': t = 7; break;
		case '8': t = 8; break;
		case '9': t = 9; break;
        }
        switch (y){
                case ',': break;
                case ')': break;
                case ':': break;
                case '0': t = t*10; break;
                case '1': t = t*10+1; break;
                case '2': t = t*10+2; break;
                case '3': t = t*10+3; break;
                case '4': t = t*10+4; break;
                case '5': t = t*10+5; break;
                case '6': t = t*10+6; break;
                case '7': t = t*10+7; break;
                case '8': t = t*10+8; break;
                case '9': t = t*10+9; break;
        }
//   		cerr << " This is t " << t << endl;
        return t;
}
//******************************************************************************
void Topology(ifstream& in, int **Topo, int n)
{
       // Topo is the structure that holds the tree topology
       // n is the number of sequences
       int TaxaToInt(char z, char y);
       char y,z;
       int *child = new int[2*n-1];
       int taxa, dadnum, childnum, level;
//CHANGE by JEFF       double length;
       // Initialize the Topology

       for(int i = 0; i < 2*n-1; i++){
           child[i] = 0;  // init the child Vector
           for(int j = 0; j < n+2; j++){
               Topo[i][j] = 0; // init the child Vector
           }
       }
       // Begin reading the topology
       dadnum = n; // will keep track of what ) you are on
       level = (-1);
       in >> z;
       while(z != ';'){
      // outf << child << endl;
           switch(z){
           case '(':
              level++;
              in >> z ;
              for (int i=(2*n-3); i >=0; i--) {
                  child[i+1] = child[i];
              }
              child[0] = 0;
              break;
/*
   This will not be needed in this algorithm.
   The branch lengths are always going to be estimated
           case ':':
               in >> length;
               branch[taxa] = length;
               in >> z; break;
*/
           case ',':
               in >> z; break;

           case ')':
               dadnum++;
               taxa = dadnum;
               childnum=1;
               while (childnum <= Topo[level][n]) {
                   Topo[dadnum-1][childnum-1] = child[0];
                   Topo[child[0]-1][n+1] = dadnum;
                   for (int i=0; i < 2*n-2; i++) {
                       child[i] = child[i+1]; //shifts the entire child nodes up
                   }
                   child[2*n-2]=0; //re-init the last child slot
                   childnum++;
               }
               if (child[0] == 0) {
                   child[0] = dadnum;
               }
               else{
                   for (int i=2*n-3; i>=0; i--) {
                       child[i+1] = child[i];
                   } // this shifts everything down
                   child[0] = dadnum;
               }
               if (level > 0) {
                   Topo[level][n]=0;
                   Topo[level-1][n]++;
               }
               level--;
               in >> z; break;
           default:
               Topo[level][n]++;
               in >> y;
               taxa = TaxaToInt(z,y);
               if (y == ',' || y == ')' || y == ':') {
                  z = y;
               }
	           else{
	       	       in >> z;
	           }
               if (child[0]== 0) {
                   child[0] = taxa;
               }
               else{
                   for (int i=2*n-3; i>=0; i--) {
                       child[i+1] = child[i];
                   }
                   child[0] = taxa;
               }
               break;
           } // end the switch statement
       } // end the while statemnt
       Topo[0][n] = dadnum-1;
}
//******************************************************************************
void FillTopology(ifstream& in, Information *C, Information *P)
{
   void Topology(ifstream& in, int **Topo, int n);

//We now have to read in the topology
   int numSeq, numNodes;
//CHANGE by JEFF   int numSeq, numNodes, childHolder;
   numSeq = C[0].numSeq;
   numNodes = C[0].numNodes;

   int **Topo = new int *[2*numSeq - 1];
   for(int i = 0; i < (2*numSeq - 1); i++){Topo[i] = new int[numSeq + 2];}
   Topology(in, Topo, numSeq);
   int *tempChild = new int[numSeq];

   for(int node = 0; node < numNodes; node++)
   {
      if(node < numSeq)
      {
         if((Topo[node][numSeq+1]) >= numNodes)
         {
            C[node].parent = P[node].parent = numNodes - 1;
         }
         else
         {
            C[node].parent = P[node].parent = Topo[node][numSeq+1] - 1;
         }
         C[node].numChild  = P[node].numChild  = 0;
      }
      else
      {
         if(node < (numNodes-1))
         {
            int childCount = 0;
            int childHolder = Topo[node][childCount];
            if((Topo[node][numSeq+1]) >= numNodes)
            {
               C[node].parent = P[node].parent = numNodes - 1;
            }
            else
            {
               C[node].parent = P[node].parent = Topo[node][numSeq+1] - 1;
            }
            while(childHolder != 0)
            {
               childCount++;
               childHolder = Topo[node][childCount];
            }
            C[node].numChild  = P[node].numChild  = childCount;
            C[node].child = new int[childCount];
            P[node].child = new int[childCount];
            for(int ch = 0; ch < childCount; ch++)
            {
               C[node].child[ch] = P[node].child[ch] = (Topo[node][ch] - 1);
            }
         }
         else
         {
            int childCount = 0;
            int childHolder = Topo[node][childCount];
            while(childHolder != 0)
            {
               tempChild[childCount] = childHolder;
               childCount++;
               childHolder = Topo[node][childCount];
            }
            int cCount1 = childCount;
            childCount = 0;
            childHolder = Topo[node+1][childCount];
            while(childHolder != 0)
            {
               if(childHolder <= node)
               {
                  tempChild[cCount1]=childHolder;
                  cCount1++;
               }
               childCount++;
               childHolder = Topo[node+1][childCount];
            }
            C[node].numChild  = P[node].numChild  = cCount1;
            C[node].child = new int[cCount1];
            P[node].child = new int[cCount1];
            for(int ch = 0; ch < cCount1; ch++)
            {
               C[node].child[ch] = P[node].child[ch] = (tempChild[ch]-1);
            }
         }
      }
   }
   delete[] tempChild;
   for(int i = 0; i < (2*numSeq - 1); i++){delete[] Topo[i];}
   delete[] Topo;
}
//******************************************************************************
void SeqFreq(Bayes *pE,Information *C,Information *P,long *seed,ifstream& in,
                        ofstream& o)
{
// This function reads in the sequences, and calculates the relative
// frequecies of the nucleotides, and also fills purPyr correctly
   int NucToInt(char m);
   void PupkoNucChoice(Bayes *pE, Information *C, Information *P,
               long *seed);
   void RateMatrix(double G, double W);
   void EigenCalc(double time);
   double Pois(double rate, int numMuts);
   void FillTopology(ifstream& in, Information *C, Information *P);
   void initPath(Path *C, Path *P, int seqlen, int *seq);
   double total;
   char tmpnuc;
   int c[3], nucSeqLen, numSeq;
   nucSeqLen = C[0].len;
   numSeq = pE->numSeq;
/********************************************************
We need something here to say that we have the internal
node sequence, and we can simulate from it to the tip
sequences.  So at this point we should have all of the
sequences at the internal nodes.
   For multiple sequences, we should have a clear
understanding of how we are going to traverse any
tree that we are given.  These internal node sequences are
crucial to how the calculation will work and proceed.
It is in this function that the initial "sample" of
internal node sequences will be set and used as a
starting point for the whole program.
********************************************************/
   for(int f = 0; f < 12; f++){F84_freq[f] = 0.0;}
   int *N1, *N2;
   if(numSeq == 2)  // If there are only 2 seqs, thus one branch between them
   {
      for(int i = 1; i >= 0; i--)
      {
         N1 = P[i].seq;
         N2 = C[i].seq;
         for(int j = 0; j < nucSeqLen; j++)
         {
            in >> tmpnuc;   // Read in from the file
            // Converts char to int
            c[j%3] = *N1 = *N2 = NucToInt(tmpnuc);  N1++; N2++;
            if((j%3) == 2)
            {  // This will be the ancestral amino acid sequence as well!
               P[i].AA_seq[j/3]=C[i].AA_seq[j/3]=Nuc2AATable[c[0]][c[1]][c[2]];
            }
         }
         if(!i)
         {
            P[0].numChild   =   C[0].numChild = 0;
            P[0].parent = 1;    C[0].parent = 1;
            P[0].child = NULL;  C[0].child = NULL;
            P[0].parentSeq   = P[1].seq;
            C[0].parentSeq   = C[1].seq;
            P[0].parentAAseq = P[1].AA_seq;
            C[0].parentAAseq = C[1].AA_seq;
         }
         else
         {
            P[1].numChild = C[1].numChild = 1;
// CHANGE below by Jeff
            P[1].parent = 0;     C[1].parent = 0;
            P[1].child = new int[1];
            C[1].child = new int[1];
            P[1].child[0] = C[1].child[0] = 0;
            P[1].parentSeq = NULL;  C[1].parentSeq = NULL;
            P[1].parentAAseq = NULL;
            C[1].parentAAseq = NULL;

         }
      }
   }  // This does it for the two sequence case!!!
   else // If there are greater than two sequences, then we have internal nodes
   {    // This one is going to be vastly more complicated
      for(int i = 0; i < numSeq; i++)
      {
         N1 = P[i].seq;
         N2 = C[i].seq;
         for(int j = 0; j < nucSeqLen; j++)
         {
            in >> tmpnuc;    // Read in from the file
            // Converts char to int
            c[j%3] = *N1 = *N2 = NucToInt(tmpnuc);
            N1++; N2++;
            if((j%3) == 2)
            {  // This will be the ancestral amino acid sequence as well!
               P[i].AA_seq[j/3]=C[i].AA_seq[j/3]=Nuc2AATable[c[0]][c[1]][c[2]];
            }
         }
      }
   }

   total = F84_freq[0] + F84_freq[1] + F84_freq[2] + F84_freq[3];
   F84_freq[0] /= total;  F84_freq[1] /= total;
   F84_freq[2] /= total;  F84_freq[3] /= total;
   F84_freq[4]  = F84_freq[0] + F84_freq[1];    // Pi_A + Pi_C
   F84_freq[5]  = F84_freq[4] + F84_freq[2];    // Pi_A + Pi_C + Pi_G
   F84_freq[6]  = F84_freq[0] + F84_freq[2];    // Pi_A + Pi_G  PURINES
   F84_freq[7]  = F84_freq[1] + F84_freq[3];    // Pi_C + Pi_T  PYRIMIDINES
   F84_freq[8]  = F84_freq[0] / F84_freq[6];    // Pi_A / (Pi_A + Pi_G)
   F84_freq[9]  = F84_freq[1] / F84_freq[7];    // Pi_C / (Pi_C + Pi_T)
   F84_freq[10] = 1.0 - F84_freq[8];         // Pi_G / (Pi_A + Pi_G)
   F84_freq[11] = 1.0 - F84_freq[9];         // Pi_T / (Pi_C + Pi_T)
   o << "Here are the F84 nucfreqs: A:" << F84_freq[0] << " C:" <<
   F84_freq[1] << " G:" << F84_freq[2] << " T:" << F84_freq[3] << endl << endl;

   double t1, h1, G, W; t1 = 0.0;
   for(int q = 0; q < 4; q++){
      if(q%2 == 0){h1 = F84_freq[6];}
      else{h1 = F84_freq[7];}
      t1 += F84_freq[q]*((1.0 - F84_freq[q]) +
            K1*(1.0 - (F84_freq[q] / h1)));
   }
   G = 1.0 / t1;
   W = G * K1;
   G *=L1;
   W *=L1;
   o << "For the F84 parameters, ";
   o << "This is G = " << G << "  and this is W " << W << endl;
//***********************
//***Begin Rate Matrix***     This routine calculates the rate matrix
//***********************
   RateMatrix(G,W);// Because this updates both g and w

   EigenCalc(1.0);
//   Change by JEFF EigenCalc(L1);
   for(int i = 0; i < 16; i++)
   {
      POISSON_G[i]  = Pois(G, i);
      POISSON_W[i]  = Pois(W, i);
      if(i == 0){o << POISSON_G[i] << "  " << POISSON_W[i] << endl;}
   }
   //We now fill the topology and our information structures
   FillTopology(in, C, P);
   o << "We are done filling the topology" << endl;
   for(int i = 0; i <= pE->numBranch; i++)
   {
      o<< "Node " << i << " has parent " << C[i].parent << endl;
      if(C[i].numChild > 0)
      {
         for(int s = 0; s < C[i].numChild; s++)
         {
            o << "Node " << i << " has " << C[i].numChild << " children " << C[i].child[s] << endl;
         }
      }
   }
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//  We must figure out what to do here for > 2 sequences
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
   if(numSeq > 2)
   {
     PupkoNucChoice(pE,C,P,seed); // This is how we choose all internal nodes
     o << "We are done with Pupko NucChoice" << endl;
     // We now fill in all of the parent nucleotide and amino acid sequences
     // We point each node to its correct parent.  We should not need to
     //   touch this again for the rest of the program.  The parent sequences
     //   start each of the branches.
     for(int i = 0; i < pE->numBranch; i++)
     {
        P[i].parentSeq   = P[P[i].parent].seq;
        C[i].parentSeq   = C[C[i].parent].seq;
        P[i].parentAAseq = P[P[i].parent].AA_seq;
        C[i].parentAAseq = C[C[i].parent].AA_seq;
     }
   }
   /////////////////////////////////////////////////////////////////
   // This is where we create and fill all of the sequence paths!!!/
   /////////////////////////////////////////////////////////////////
   for(int i = 0; i < pE->numBranch; i++)
   {
      C[i].seqPath = new Path [nucSeqLen];
      P[i].seqPath = new Path [nucSeqLen];
      initPath(C[i].seqPath, P[i].seqPath, nucSeqLen, C[i].parentSeq);
   }


   // We now calculate the information we need about the
   //   internal node sequence for the pi updates.
   // Because these values can change, we must treat them carefully.

   int tempA, tempC, tempG;   tempA = tempC = tempG = 0;
   int *C_seq = C[C[0].numBranch].seq;// There are numNodes Nodes, but they are
                                      //  numbered 0 - numBranch
   for(int i = 0; i < nucSeqLen; i++)
   {
      switch(*C_seq){
         case 0: tempA++; break;
         case 1: tempC++; break;
         case 2: tempG++; break;
         case 3:          break;
      }
      C_seq++;
   }

   pE->number_A = tempA;
   pE->number_C = tempC;
   pE->number_G = tempG;

   o << "Number A: " << tempA << " Number C: " << tempC;
   o << " Number G: " << tempG << " Number T: ";
   o << (nucSeqLen - tempA - tempC - tempG) << endl;

} // End SeqFreq Subroutine
//******************************************************************************
void SetParams(Bayes *pE, long *seed, ofstream &out)
{
        double rnd(long* seed);


        pE->kappa[0]    = INIT_K;
        pE->omega[0]    = INIT_W;
        pE->solvent  = INIT_W_S;
        pE->pairwise = INIT_W_P;
        pE->rate[0]  = INIT_U;// Should fix this so we can input branch lengths
/* #if SIM_U                    // and fix them at these values for the calculation
        for(int i = 0; i < pE->numBranch; i++)
        {
           pE->rate[i]  = U_MIN + rnd(seed)*(U_MAX - U_MIN);
        }
#endif
*/
#if SIM_U                    // and fix them at these values for the calculation
        for(int i = 0; i < pE->numBranch; i++)
        {
           pE->yangBranch[i]  = U_MIN + rnd(seed)*(U_MAX - U_MIN);
        }
#endif
#if SIM_K
        if(MULTI_K)
        {
           for(int i = 0; i < pE->numBranch; i++)
           {
              pE->kappa[i]  = K_MIN + rnd(seed)*(K_MAX - K_MIN);
           }
        }
        else
        {
              pE->kappa[0]  = K_MIN + rnd(seed)*(K_MAX - K_MIN);
        }
#endif
#if SIM_W
        if(MULTI_W)
        {
           for(int i = 0; i < pE->numBranch; i++)
           {
              pE->omega[i]  = W_MIN + rnd(seed)*(W_MAX - W_MIN);
           }
        }
        else
        {
              pE->omega[0]  = W_MIN + rnd(seed)*(W_MAX - W_MIN);
        }
#endif
#if SIM_W_S
        pE->solvent  = W_S_MIN + rnd(seed)*(W_S_MAX - W_S_MIN);
#endif
#if SIM_W_P
        pE->pairwise  = W_P_MIN + rnd(seed)*(W_P_MAX - W_P_MIN);
#endif
        pE->MCMCnf[0] = pE->tempMCMCnf[0] = ADENINE;
        pE->MCMCnf[1] = pE->tempMCMCnf[1] = CYTOSINE;
        pE->MCMCnf[2] = pE->tempMCMCnf[2] = GUANINE;
        pE->MCMCnf[3] = pE->tempMCMCnf[3] = 1.0-(ADENINE+CYTOSINE+GUANINE);
        out << "Done initializing everything " << endl;
        if(MULTI_K)
        {
           for(int q = 0; q < pE->numBranch; q++)
           {
              out << "Init K(" << q << ") = " << pE->kappa[q] << " ";
           }
        }
        else{out << "Init K(0) = " << pE->kappa[0] << " ";}
        if(MULTI_W)
        {
           for(int q = 0; q < pE->numBranch; q++)
           {
              out << "Init W(" << q << ") = " << pE->omega[q] << " ";
           }
        }
        else{out << "Init W(0) = " << pE->omega[0] << " ";}
        out << endl << ", w_S = " << pE->solvent << " Init w_p = " <<
        pE->pairwise;
        for(int q = 0; q < pE->numBranch; q++)
        {
           out << " Init U(" << q << ") = " << pE->rate[q] << " ";
        }
        out << endl << "and freqs A: " << pE->MCMCnf[0] <<
        " C: " << pE->MCMCnf[1] << " G: " << pE->MCMCnf[2] <<" T: "<<
        pE->MCMCnf[3] << endl;
}
//******************************************************************************
/* void UpdateProbMatrix(Information *C, Information *P)
{
        void EigenCalc(double** RateMat, double** ProbMat, double time);
        double Pois(double rate, int numMuts);
        double t1, h1, G, W; t1 = 0.0;
        for(int q = 0; q < 4; q++){
           if(q%2 == 0){h1 = F84_freq[6];}
           else{h1 = F84_freq[7];}
           t1 += F84_freq[q]*((1.0 -F84_freq[q]) +
                 K1*(1.0 - (F84_freq[q] / h1)));
        }
        P->g = C->g = L1 / t1;
        P->w = C->w = C->g * K1;
        EigenCalc(C->F84Rate, C->F84Prob, L1);
        G = C->g; W = C->w;
        for(int i = 0; i < 16; i++)
        {
           C->sumG[i]  = Pois(G, i);
           C->sumW[i]  = Pois(W, i);
           P->sumG[i] = C->sumG[i];
           P->sumW[i] = C->sumW[i];
        }
        for(int i = 0; i < 4; i++){
           for(int j = 0; j < 4; j++){
              P->F84Prob[i][j] = C->F84Prob[i][j];
           }
        }
} // End UpdateProbMatrix
*/
//******************************************************************************
void RateMatrix(double G, double W)
{  // This creates the F84 instantaneous rate matrix

   for(int i = 0; i < 4; i++){
     for(int j = 0; j < 4; j++){
       if(i != j){
         if(((i+j)%2) == 1){F84_RATE[i][j] = G*F84_freq[j];}
         else{
           F84_RATE[i][j] = G*F84_freq[j] + W*F84_freq[j + 8];
         }
       } // Initialize before the calculation
       else{F84_RATE[i][j] = 0.0;}
     }
   }
   for(int i = 0; i < 4; i++){
      for(int j = 0; j < 4; j++){
         if(i!=j){F84_RATE[i][i] -= F84_RATE[i][j];}
      }
   }
   double temp; temp = 0.0;
//   void Normalize(double* relfreq, double& temp);
//   Normalize(F84_freq, temp);
}  // End RateMatrix Subroutine
//******************************************************************************
void Normalize(double* relfreq, double& temp)
{  //This functions to normalize any matrix that is given
   // To normalize we must SUM: -1*pi_i*R_{ii}
   for(int i = 0; i < 4; i++){temp += -1.0*relfreq[i]*F84_RATE[i][i];}
   for(int i = 0; i < 4; i++){
      for(int j = 0; j < 4; j++){F84_RATE[i][j] = F84_RATE[i][j] / temp;}
   }
}  // End Normalize Subroutine
//******************************************************************************
void EigenCalc(double time)
{   // This function uses Yang's routine to calculate the eigen-values
   void ProbMatrix(double **ev,double **evi, double **eval, double time);

   double ri[4], A[16], rr[4], vr[16], vi[16], q[8];
   for(int j=0; j < 4; j++){
      for(int k=0; k < 4; k++){
         A[4*j+k] = F84_RATE[j][k];
      }
   } // We now find the eigenvalues and eigenvectors using Yang routine
   eigen(1,A,4,rr,ri,vr,vi,q);

   double **ev, **evi, **eval;
   ev = new double *[4];  evi = new double *[4];  eval = new double *[4];

   for(int j=0; j < 4; j++){
      ev[j] = new double[4]; evi[j] = new double[4]; eval[j] = new double[4];
      for(int k=0; k < 4; k++){
         eval[j][k] = 0.0;
         ev[j][k] = vr[4*j+k];       //Puts eigenvectors into the matrix
         if(j==k){eval[j][j] = rr[j];} //Puts eigenvalues into my matrix
      }
   }

   //We now calculate the inverse of a 4 x 4 matrix directly
   double ev00,ev01,ev02,ev03,ev10,ev11,ev12,ev13,ev20,ev21,ev22,ev23,ev30,ev31,
          ev32,ev33, dnom;

   ev00 = ev[0][0]; ev01 = ev[0][1]; ev02 = ev[0][2]; ev03 = ev[0][3];
   ev10 = ev[1][0]; ev11 = ev[1][1]; ev12 = ev[1][2]; ev13 = ev[1][3];
   ev20 = ev[2][0]; ev21 = ev[2][1]; ev22 = ev[2][2]; ev23 = ev[2][3];
   ev30 = ev[3][0]; ev31 = ev[3][1]; ev32 = ev[3][2]; ev33 = ev[3][3];
   dnom = ev00*(-ev11*ev22*ev33 + ev11*ev23*ev32 + ev21*ev12*ev33 -
                 ev21*ev13*ev32 - ev31*ev12*ev23 + ev31*ev13*ev22) +
          ev10*(ev01*ev22*ev33 - ev01*ev23*ev32 - ev21*ev02*ev33 +
                ev21*ev03*ev32 + ev31*ev02*ev23 - ev31*ev03*ev22) +
          ev20*(-ev01*ev12*ev33 + ev01*ev13*ev32 + ev11*ev02*ev33 -
                 ev11*ev03*ev32 - ev31*ev02*ev13 + ev31*ev03*ev12) +
          ev30*(ev01*ev12*ev23 - ev01*ev13*ev22 - ev11*ev02*ev23 +
                ev11*ev03*ev22 + ev21*ev02*ev13 - ev21*ev03*ev12);
   evi[0][0] = -(-ev11*ev22*ev33 + ev11*ev23*ev32 - ev21*ev13*ev32 -
                  ev31*ev12*ev23 + ev31*ev13*ev22 + ev21*ev12*ev33)*dnom;
   evi[0][1] =  (-ev01*ev22*ev33 + ev01*ev23*ev32 + ev21*ev02*ev33 -
                  ev21*ev03*ev32 - ev31*ev02*ev23 + ev31*ev03*ev22)*dnom;
   evi[0][2] = -(-ev01*ev12*ev33 + ev01*ev13*ev32 + ev11*ev02*ev33 -
                  ev11*ev03*ev32 - ev31*ev02*ev13 + ev31*ev03*ev12)*dnom;
   evi[0][3] =  (-ev01*ev12*ev23 + ev01*ev13*ev22 + ev11*ev02*ev23 -
                  ev11*ev03*ev22 - ev21*ev02*ev13 + ev21*ev03*ev12)*dnom;
   evi[1][0] =  (-ev10*ev22*ev33 + ev10*ev23*ev32 + ev20*ev12*ev33 -
                  ev20*ev13*ev32 - ev30*ev12*ev23 + ev30*ev13*ev22)*dnom;
   evi[1][1] = -(-ev00*ev22*ev33 + ev00*ev23*ev32 + ev20*ev02*ev33 -
                  ev20*ev03*ev32 - ev30*ev02*ev23 + ev30*ev03*ev22)*dnom;
   evi[1][2] =  -(ev00*ev12*ev33 - ev00*ev13*ev32 - ev10*ev02*ev33 +
                  ev10*ev03*ev32 + ev30*ev02*ev13 - ev30*ev03*ev12)*dnom;
   evi[1][3] =   (ev00*ev12*ev23 - ev00*ev13*ev22 - ev10*ev02*ev23 +
                  ev10*ev03*ev22 + ev20*ev02*ev13 - ev20*ev03*ev12)*dnom;
   evi[2][0] = -(-ev10*ev21*ev33 + ev10*ev23*ev31 + ev20*ev11*ev33 -
                  ev20*ev13*ev31 - ev30*ev11*ev23 + ev30*ev13*ev21)*dnom;
   evi[2][1] =  -(ev00*ev21*ev33 - ev00*ev23*ev31 - ev20*ev01*ev33 +
                  ev20*ev03*ev31 + ev30*ev01*ev23 - ev30*ev03*ev21)*dnom;
   evi[2][2] =   (ev00*ev11*ev33 - ev00*ev13*ev31 - ev10*ev01*ev33 +
                  ev10*ev03*ev31 + ev30*ev01*ev13 - ev30*ev03*ev11)*dnom;
   evi[2][3] =  -(ev00*ev11*ev23 - ev00*ev13*ev21 - ev10*ev01*ev23 +
                  ev10*ev03*ev21 + ev20*ev01*ev13 - ev20*ev03*ev11)*dnom;
   evi[3][0] =  (-ev10*ev21*ev32 + ev10*ev22*ev31 + ev20*ev11*ev32 -
                  ev20*ev12*ev31 - ev30*ev11*ev22 + ev30*ev12*ev21)*dnom;
   evi[3][1] =   (ev00*ev21*ev32 - ev00*ev22*ev31 - ev20*ev01*ev32 +
                  ev20*ev02*ev31 + ev30*ev01*ev22 - ev30*ev02*ev21)*dnom;
   evi[3][2] =  -(ev00*ev11*ev32 - ev00*ev12*ev31 - ev10*ev01*ev32 +
                  ev10*ev02*ev31 + ev30*ev01*ev12 - ev30*ev02*ev11)*dnom;
   evi[3][3] =   (ev00*ev11*ev22 - ev00*ev12*ev21 - ev10*ev01*ev22 +
                  ev10*ev02*ev21 + ev20*ev01*ev12 - ev20*ev02*ev11)*dnom;

   ProbMatrix(ev,evi,eval,time);
   for(int j=0; j < 4; j++){delete[] ev[j];delete[] evi[j];delete[] eval[j];}
   delete[] ev;     delete[] evi;        delete[] eval;
}  // End EigenCalc Subroutine
//******************************************************************************
void MatrixMult(double **A, double **B,double **D)
{
   D[0][0] = A[0][0]*B[0][0]+A[0][1]*B[1][0]+A[0][2]*B[2][0]+A[0][3]*B[3][0];
   D[0][1] = A[0][0]*B[0][1]+A[0][1]*B[1][1]+A[0][2]*B[2][1]+A[0][3]*B[3][1];
   D[0][2] = A[0][0]*B[0][2]+A[0][1]*B[1][2]+A[0][2]*B[2][2]+A[0][3]*B[3][2];
   D[0][3] = A[0][0]*B[0][3]+A[0][1]*B[1][3]+A[0][2]*B[2][3]+A[0][3]*B[3][3];

   D[1][0] = A[1][0]*B[0][0]+A[1][1]*B[1][0]+A[1][2]*B[2][0]+A[1][3]*B[3][0];
   D[1][1] = A[1][0]*B[0][1]+A[1][1]*B[1][1]+A[1][2]*B[2][1]+A[1][3]*B[3][1];
   D[1][2] = A[1][0]*B[0][2]+A[1][1]*B[1][2]+A[1][2]*B[2][2]+A[1][3]*B[3][2];
   D[1][3] = A[1][0]*B[0][3]+A[1][1]*B[1][3]+A[1][2]*B[2][3]+A[1][3]*B[3][3];

   D[2][0] = A[2][0]*B[0][0]+A[2][1]*B[1][0]+A[2][2]*B[2][0]+A[2][3]*B[3][0];
   D[2][1] = A[2][0]*B[0][1]+A[2][1]*B[1][1]+A[2][2]*B[2][1]+A[2][3]*B[3][1];
   D[2][2] = A[2][0]*B[0][2]+A[2][1]*B[1][2]+A[2][2]*B[2][2]+A[2][3]*B[3][2];
   D[2][3] = A[2][0]*B[0][3]+A[2][1]*B[1][3]+A[2][2]*B[2][3]+A[2][3]*B[3][3];

   D[3][0] = A[3][0]*B[0][0]+A[3][1]*B[1][0]+A[3][2]*B[2][0]+A[3][3]*B[3][0];
   D[3][1] = A[3][0]*B[0][1]+A[3][1]*B[1][1]+A[3][2]*B[2][1]+A[3][3]*B[3][1];
   D[3][2] = A[3][0]*B[0][2]+A[3][1]*B[1][2]+A[3][2]*B[2][2]+A[3][3]*B[3][2];
   D[3][3] = A[3][0]*B[0][3]+A[3][1]*B[1][3]+A[3][2]*B[2][3]+A[3][3]*B[3][3];
}
//******************************************************************************
void ProbMatrix(double **ev,double **evi, double **eval, double time)
{ // This takes the eigenvalues and eigenvectors and calculates the prob
   void MatrixMult(double **A, double **B,double **D);
   double **tempProb = new double *[4];
   double **temp2  = new double *[4];

   for(int i = 0; i < 4; i++)
   {
      tempProb[i] = new double[4];
      temp2[i]    = new double[4];
      eval[i][i]  = exp(time * eval[i][i]);
   }

   MatrixMult(ev,eval,tempProb);
   MatrixMult(tempProb,evi,temp2);
   if(temp2[0][0] < 0.0)
   {
      for(int i = 0; i < 4; i++){for(int j = 0; j < 4; j++){
         F84_PROB[i][j] = -temp2[i][j];}
      }
   }
   else
   {
      for(int i = 0; i < 4; i++){for(int j = 0; j < 4; j++){
         F84_PROB[i][j] = temp2[i][j];}
      }
   }

   for(int i = 0; i < 4; i++){delete []tempProb[i];delete[] temp2[i];}
   delete[] tempProb;   delete[] temp2;
}  // End ProbMatrix Subroutine
//******************************************************************************
void SetUpBayes(Bayes *pE)
{
   int numBranch = pE->numBranch;
   if(MULTI_K){pE->kappa = new double[numBranch];}
   else{pE->kappa = new double[1];}
   if(MULTI_W){pE->omega = new double[numBranch];}
   else{pE->omega = new double[1];}
   pE->rate       = new double[numBranch];
   pE->tempRate   = new double[numBranch];
   pE->yangBranch = new double[numBranch];
   pE->branchLen  = new double[numBranch];
   pE->MCMCnf     = new double[4];
   pE->tempMCMCnf = new double[4];
}
//******************************************************************************
void SetUpInfo(Bayes *pE, Information *C, Information *P, long* seed,
                ifstream &in, ofstream &out)
{
        void initI(Information &I);
        void SeqFreq(Bayes *pE, Information *C, Information *P, long *seed,
                        ifstream& in, ofstream& o);
        void SetUpBayes(Bayes *pE);
        void SetParams(Bayes *pE, long *seed, ofstream &out);

        int numBranch, nucSeqLen, AAseqLen;
        double invlen, invAAlen;
        numBranch = pE->numBranch;
        in >> nucSeqLen; // get length of nucleotide sequence
        invlen = nucSeqLen;
        invlen   = 1.0 / invlen;         // inverse of the nucleotide length
        invAAlen = AAseqLen = nucSeqLen / 3;   // Amino acid sequence length
        invAAlen = 1.0 / invAAlen;  // The inverse of AA length
        out << "This is the nuc length = " << nucSeqLen <<
                 ", and AA " << AAseqLen << endl;
        out << "This is the seed " << *seed << endl;
        int numNodes = numBranch+1;
        for(int doug = 0; doug < numNodes; doug++)
        {  //We now initialize all of the information Structures
           initI(C[doug]);
           initI(P[doug]);
           P[doug].numSeq   = C[doug].numSeq   = pE->numSeq;
           P[doug].numBranch= C[doug].numBranch= numBranch;
           P[doug].numNodes = C[doug].numNodes = numNodes;
           P[doug].len      = C[doug].len      = nucSeqLen;
           P[doug].AAlen    = C[doug].AAlen    = AAseqLen;
           P[doug].invlen   = C[doug].invlen   = invlen;
           P[doug].invAAlen = C[doug].invAAlen = invAAlen;
           C[doug].seq  = new int[nucSeqLen]; // This is the Current nuc sequence
           P[doug].seq  = new int[nucSeqLen]; // This is the Proposed nuc sequence
           C[doug].AA_seq   = new int[AAseqLen]; // This is the AA sequence
           P[doug].AA_seq   = new int[AAseqLen]; // This is the AA sequence
           C[doug].neighAcc = new int[AAseqLen];
           P[doug].neighAcc = new int[AAseqLen];
           out << "We are done init a lot of C[" << doug << "]" << endl;
        }  // Done initializing most of the information matrices,  HUGE JOB!!
        out << "We are about to set up bayes " << endl;
        SetUpBayes(pE);
// This function reads in the sequences, and calculates the relative
// frequencies of the nucleotides, and also fills purPyr correctly
        out << "We are about to set SeqFreq " << endl;
         out << "This is the seed " << *seed << endl;
        SeqFreq(pE, C, P, seed, in, out);
         out << "This is the seed " << *seed << endl;
        out << "We are about to SetParams " << endl;
        SetParams(pE,seed,out);
//**********************
//***End Rate Matrix ***
//***Begin Prob Matrix**  Given the rate matrix, this calcs the prob matrix
//**********************
        //To calculate the eigenvalues and eigenvectors we need to
        // use PAML and convert our matrix into a usable form
        //For this the matrix class library is the best choice around
        //This routine also calculates the F84 probability matrix
        // Now we have the probability of change in our F84 format
        out << endl << "Here is the Rate matrix: " << endl << endl ;
        for(int i = 0; i < 4; i++){
           for(int j = 0; j < 4; j++){
              out << F84_RATE[i][j] << " ";
           }
           out << endl;
        }

        out << endl << "Here is the Probability matrix: " << endl << endl ;
        for(int i = 0; i < 4; i++){
           for(int j = 0; j < 4; j++){
              out << F84_PROB[i][j] << " ";
           }
           out << endl;
        }
}  // End SetUpInfo Subroutine
//******************************************************************************
void SetSolvent(Interaction *FO, int AAlen, ifstream& in,
        ifstream &inter, ofstream& out)
{
   int SolventAccess(double acc);
///////////////////////////////////////////////////////////////
// Now we are going to set up the Solvent vector, which will //
// tell us exatly how accessible each amino acid is to the   //
// surrounding solvent. We also read in DJ's Solvent         //
// Accessibilities for each AA as well.                      //
///////////////////////////////////////////////////////////////

   int *solvAcc = new int[AAlen]; // Initiallize ALWAYS!
   // for(int i = 0; i < AAlen; i++){FirstOrder->solvAcc[i] = 0;}
   // Like in the contact accessibility, we call this only once, but never
   // again.  Then the column number should be fixed from here on in.
   double s, *T1, *T2;
   out << endl << "This is the solvent accessibility info "<< endl;
   for(int i = 0; i < AAlen; i++)
   {
      in >> s; //      out << s << " ";
      solvAcc[i] = SolventAccess(s);
   }
   double **solv = new double *[5];
   for(int i = 0; i < 5; i++){solv[i] = new double[20];}
   for(int i = 0; i < 20; i++)
   {
      for(int j = 0; j < 5; j++){inter >> solv[j][i];}
   }
/*   for(int i = 0; i < 5; i++)
   {
      for(int j = 0; j < 20; j++){out << solv[i][j] << " ";}
      out << endl;
   }
*/
   for(int i = 0; i < AAlen; i++)
   {
      FO->AAinfo[i].solvent = new double[20];
      T1 = FO->AAinfo[i].solvent;
      T2 = solv[solvAcc[i]];
      for(int j = 0; j < 20; j++){*T1 = *T2; T1++; T2++;}
   }


   delete[] solvAcc;
   for(int k = 0; k < 5; k++){delete[] solv[k];}
   delete[] solv;
}  // End SetSolvent Subroutine
//******************************************************************************
void SetNeighbor(Interaction *FirstOrder,ifstream& inter)
{
////////////////////////////////////////////////////////
// Now we are going to set up the Nearest Neighbor    //
// information.  This includes setting up an index    //
// matrix, as well as entering the nearest neighbors. //
// The inNeigh file contains the AA neighbors, what   //
// nucleotide caused it, and if it was a TS or TV.    //
// These are constant data files that do not change   //
// even when the protein changes                      //
////////////////////////////////////////////////////////

   FirstOrder->nAccess = new int *[20];
   for(int i = 0; i < 20; i++)
   {
      FirstOrder->nAccess[i] = new int[4];
      for(int j = 0; j < 4; j++){inter >> FirstOrder->nAccess[i][j];}
   }

   FirstOrder->neighList = new int *[61];
   for(int i = 0; i < 61; i++)
   {
      FirstOrder->neighList[i] = new int[27];
      for(int j = 0; j < 27; j++){inter >> FirstOrder->neighList[i][j];}
   }
}  // End SetNeighbor Subroutine
//******************************************************************************
void SetProt(int AAlen, double **prot, ifstream& in, ofstream& out)
{
/////////////////////////////////////////////////////////
// Here we read in the coordinates from the data file. //
// Note: we only read in N, CA, O, CB atom coordinates //
// initialize the prot structure                       //
/////////////////////////////////////////////////////////
   char aaRes; // This holds the initial AA tag in the first column of file
   int count;
   double z;
   for(int residue = 0; residue < AAlen; residue++)
   {
      in >> aaRes;
      count = 0;
      for(int i = 0; i < 15; i++)
      {
         in >> z;
         if((i < 6) || (i > 8))
         {
            prot[residue][count] = z;
            count++;
         }
      }
   }

//   Check # 1 Must keep this here or else it will be deleted
#if LOUD
   out << "Here is the protein 3D coordinates " << endl << endl;
   for(int i = 0; i < AAlen; i++){
      for(int j = 0; j < 12; j++){
         out << prot[i][j] << "  ";
      }
      out << endl;
   }
   out << endl << endl;
#endif
}  // End SetProt Subroutine
//******************************************************************************
void SetDDD(Interaction *FO,int AAlen,double **DDD,int *newOrder,ofstream& out)
{
////////////////////////////////////////////////////////////////////
// Only accessible from SetDDDneigh using the matrix of CA atom   //
// distances in DDD.  If the distance is less than 10.0, then it  //
// is considered a neighbor in 3 dimensions.  Only these residues //
// will be considered in the subsequent energy calculations.      //
////////////////////////////////////////////////////////////////////
   int *tmpNeighbor = new int[AAlen];
   int count;
   double distance;
   for(int aa1 = 0; aa1 < AAlen; aa1++)
   {
      out << "col = " << aa1+1+newOrder[aa1] << " : ";
      count = 0;
      for(int aa2 = 0; aa2 < AAlen; aa2++)
      {
         if(aa1 != aa2)
         {
            if(aa1 < aa2){distance = DDD[aa1][aa2];}
            else {distance = DDD[aa2][aa1];}
            if(distance < INTERACT)
            {
               tmpNeighbor[count] = aa2;
               count++;
               out << aa2+1+newOrder[aa2] <<  " ";
            }
         }
      }
      FO->AAinfo[aa1].numNeighbors = count;
      FO->AAinfo[aa1].DDDneighbors = new int[count];
      out << " total = " << count << endl;
      for(int k = 0; k < count; k++){FO->AAinfo[aa1].DDDneighbors[k] = tmpNeighbor[k];}
   }
   delete[] tmpNeighbor;
}  // End SetDDD Subroutine
//******************************************************************************
void SetDDDneigh(Interaction *FirstOrder,int AAlen,double **prot,
                int *newOrder, ifstream &in, ofstream& o)
{
   void SetDDD(Interaction *FirstOrder,int AAlen,double **DDD,int *newOrder,
				ofstream& o);
////////////////////////////////////////////////
// Now we are going to set up the Distance    //
// matrix, in order to see how close the AA's //
// are in 3-D space.                          //
////////////////////////////////////////////////
   double **DDD = new double *[AAlen];
   for(int i =0;i < AAlen;i++){DDD[i] = new double[AAlen];}

   int *missingSites;
   int numMiss;  in >> numMiss;
   o << "This is the number of missing sites " << numMiss << endl;
   if(numMiss > 0) // If there are any missing sites in the protein, then this
   {               // handles them by putting in a space there!!!
      missingSites = new int[numMiss];
      for(int i = 0; i < numMiss; i++)
      {
         in >> missingSites[i];
         o << "Those sites in order are (" << missingSites[i] << ")" << endl;
      }
      int count = 0;
      for(int i = 0; i < AAlen; i++)
      {
        if(count == numMiss){newOrder[i] = count;}
        else
        {
          if(i < missingSites[count]){newOrder[i] = count;}
          else
          {
            if(i == missingSites[count])
            {
              count++;
              if(count < numMiss){missingSites[count] -= count;}
            }
            i--;
          }
        }
      }
      delete[] missingSites;
   }
   // Initialize the distance matrix
////////////////////////////////////////////////////////////////////////////
// The critical part of this program is that we find the distance between //
// pairs of CB atoms only.  This set will be the only atoms that affect   //
// the rate away from a particular atom.  The problem is that other atoms //
// could be within 10 A, but those will not be considered in future calcs.//
////////////////////////////////////////////////////////////////////////////
   int minusOne, resIplusOne; minusOne = AAlen-1;
   for(int res_i = 0; res_i < minusOne; res_i++){
     resIplusOne = res_i+1;
     for(int res_j = resIplusOne; res_j < AAlen; res_j++){
       if(res_i < res_j){
// For this case we use the C_A atoms the center of the 10 A ball
//         DDD[res_i][res_j]=Dist(prot[res_i][3],prot[res_i][4],prot[res_i][5],
//              prot[res_j][3],prot[res_j][4],prot[res_j][5]);
// For this case we use the C_B atoms the center of the 10 A ball
         DDD[res_i][res_j]=Dist(prot[res_i][9],prot[res_i][10],prot[res_i][11],
              prot[res_j][9],prot[res_j][10],prot[res_j][11]);
       }
     }
   }
   SetDDD(FirstOrder,AAlen,DDD,newOrder,o);
   // We now give back all of the memory that we took up until now ////////////
   for(int i = 0; i < AAlen; i++){delete[] DDD[i];}
   delete[] DDD;
}  // End SetDDDneigh Subroutine
//******************************************************************************
void SetEnergyAcc(Interaction *FirstOrder, int AAlen)
{
   // This structure will hold NRG matrix indices that are being used
   FirstOrder->energyAcc = new int **[AAlen];
   int t;
   for(int w = 0; w < AAlen; w++)
   {
      t = FirstOrder->AAinfo[w].numNeighbors;
      FirstOrder->energyAcc[w] = new int *[t];
      for(int x = 0; x < t; x++)
      {
         FirstOrder->energyAcc[w][x] = new int[5];
      }
   }
}  // End SetEnergyAcc Subroutine
//******************************************************************************
void SetRosetta(double **rosetta, ifstream& inter)
{
   // Set up the Rosetta matrix of break points
   for(int r = 0; r < 93; r++){
      for(int s = 0; s < 5; s++){inter >> rosetta[r][s];}
   }
}  // End SetRosetta Subroutine
//******************************************************************************
void SetMatrixIndex(Interaction *FirstOrder, int AAlen, double **prot,
		double **rosetta, int* newOrder, ofstream &o)
{
   int MatrixHunter(int seqsep, int pair, double **rosetta, double distance);

   int targetRes, tempMat[5], seqsep, first, next;
   double dddpair[5]; // 0:CBCB, 1:CBN, 2:CB0, 3:NCB, 4:OCB
   //These are the placements of all the zero matrices (Out of Range)
   //zero[0] = 92; zero[1] = 88; zero[2] = 92; zero[3] = 87; zero[4] = 81;
   int ***FO_energyAcc, **FO_energyAccRes;
// CHANGE by JEFF   int *FO_numNeigh, *DDDneighbor, ***FO_energyAcc, **FO_energyAccRes;
   FO_energyAcc   = FirstOrder->energyAcc;
   int FO_numNeigh_res;
   for(int res = 0; res < AAlen; res++){
      FO_numNeigh_res = FirstOrder->AAinfo[res].numNeighbors;
      FO_energyAccRes = FO_energyAcc[res];
    //  DDDneighbor = FirstOrder->AAinfo[res].DDDneighbors;
      for(int neigh = 0; neigh < FO_numNeigh_res; neigh++){
         targetRes = FirstOrder->AAinfo[res].DDDneighbors[neigh];
      //   targetRes = *DDDneighbor;
         if(res < targetRes){
            seqsep = (targetRes - res)+(newOrder[targetRes] - newOrder[res]);
            first = res; next = targetRes;
         }
         else
         {
            seqsep = (res - targetRes)+(newOrder[res] - newOrder[targetRes]);
            first = targetRes; next = res;
         }

         // CB x CB
         dddpair[0] = Dist(prot[first][9],prot[first][10],prot[first][11],
               prot[next][9],prot[next][10],prot[next][11]);
         // CB x N
         dddpair[1] = Dist(prot[first][9],prot[first][10],prot[first][11],
               prot[next][0],prot[next][1],prot[next][2]);
         // CB x O
         dddpair[2] = Dist(prot[first][9],prot[first][10],prot[first][11],
               prot[next][6],prot[next][7],prot[next][8]);
         // N x CB
         dddpair[3] = Dist(prot[first][0],prot[first][1],prot[first][2],
               prot[next][9],prot[next][10],prot[next][11]);
         // O x CB
         dddpair[4] = Dist(prot[first][6],prot[first][7],prot[first][8],
               prot[next][9],prot[next][10],prot[next][11]);

         for(int i = 0; i < 5; i++){
            if(dddpair[i] < INTERACT){
               tempMat[i] = MatrixHunter(seqsep, i, rosetta, dddpair[i]);
            }
            else{tempMat[i] = MatrixHunter(seqsep, i, rosetta, 9.9999);}
            //else{tempMat[i] = zero[i];}
            FO_energyAccRes[neigh][i] = tempMat[i];
//            cout << tempMat[i] << " ";
         }
     //    DDDneighbor++;
//         cout << endl;
      } // End running through all of the AA's close in 3 dimensions
   } // End running through the residues in the protein
}// End Set MatrixIndex Subroutine
//******************************************************************************
void SetNRGmats(Interaction *FirstOrder, ifstream& inter)
{
////////////////////////////////////////////////////////
// Now we are going to set up the HUGE ENERGY Matrix. //
// This MASSIVE structure includes 450 20x20 amino    //
// acid matrices of doubles.  This is for multiple    //
// distances for close, middle and distant            //
// interactions.                                      //
////////////////////////////////////////////////////////
   int zero[5];
   zero[0] = 93; zero[1] = 89; zero[2] = 93; zero[3] = 88; zero[4] = 82;
   // THIS IS IT!!!  We now set up the HUGE energy matrix

   FirstOrder->energy = new double ***[5];
   for(int i = 0; i < 5; i++){
     FirstOrder->energy[i] = new double **[zero[i]];
     for(int j = 0; j < zero[i]; j++){
       FirstOrder->energy[i][j] = new double *[20];
       for(int k = 0; k < 20; k++){
         FirstOrder->energy[i][j][k] = new double[20];
         for(int l = 0; l < 20; l++){inter >> FirstOrder->energy[i][j][k][l];}
       }
     }
   }
} // End SetNRGmats Subroutine
//******************************************************************************
void CheckInteract(Interaction *FirstOrder, int AAlen, ofstream& out)
{
      // This checks to see if we have copied many of the input
      // files properly.  If they are wrong then this is NOT GOOD!
      out << endl << "Here is the Solvent accessibility information " << endl;
      // Check on the system

      out << endl << "Here is the Neighbor Access Matrix " << endl;
      for(int i = 0; i < 20; i++){
         for(int j = 0; j < 4; j++){
            out << FirstOrder->nAccess[i][j] << " ";
         }
         out << endl;
      }

      out << endl << "Here is the 9 nearest neighbor list" << endl;
      for(int i = 0; i < 61; i++){
         for(int j = 0; j < 27; j++){
            out << FirstOrder->neighList[i][j] << " ";
         }
         out << endl;
      }
}  // End CheckInteract Subroutine
//******************************************************************************
void SetUpInteraction(Interaction *FirstOrder, int AAlen, ifstream& in,
        ifstream &inter, ofstream& out)
{
   void SetSolvent(Interaction *FirstOrder, int AAlen, ifstream& in,
                ifstream &inter, ofstream& out);
   void SetNeighbor(Interaction *FirstOrder,ifstream& inter);
   void SetProt(int AAlen, double **prot, ifstream& in, ofstream& out);
   void SetDDDneigh(Interaction *FirstOrder,int AAlen,double **prot,
                int *NewOrder, ifstream &in, ofstream& o);
   void SetEnergyAcc(Interaction *FirstOrder, int AAlen);
   void SetRosetta(double **rosetta,ifstream& inter);
   void SetMatrixIndex(Interaction *FirstOrder, int AAlen, double **prot,
		double **rosetta, int* newOrder, ofstream& o);
   void SetNRGmats(Interaction *FirstOrder, ifstream& inter);
   void CheckInteract(Interaction *FirstOrder, int AAlen, ofstream& out);

   FirstOrder->AAinfo = new AAsiteInfo[AAlen];
      out << "We are setting up interaction " << endl;
   SetSolvent(FirstOrder, AAlen, in, inter, out);
   SetNeighbor(FirstOrder,inter);
   int length;
   in >> length;
   if(length != AAlen){cerr << "THIS FILE DOES NOT MATCH THE PROTEIN!" << endl;}


   int *newOrder = new int[AAlen];
   double **prot = new double *[AAlen];
   for(int i = 0; i < AAlen; i++)
   {
      prot[i] = new double[12];
      newOrder[i] = 0;
   }
   SetProt(AAlen,prot,in,out);
   SetDDDneigh(FirstOrder,AAlen,prot,newOrder,in,out);
   out << "We are done with SetDDDneigh " << endl;

   SetEnergyAcc(FirstOrder,AAlen);
   out << "We are done with SetEnergyAcc " << endl;

   double **rosetta = new double *[93];
   for(int r = 0; r < 93; r++){rosetta[r] = new double[5];}

   SetRosetta(rosetta,inter);
   SetMatrixIndex(FirstOrder,AAlen,prot,rosetta,newOrder,out);
   out << "We are done with SetMatrixIndex " << endl;
   // We now delete the memory that stored rosetta    //
   for(int i = 0; i < AAlen; i++){delete[] prot[i];}  //
   for(int r = 0; r < 93; r++){delete[] rosetta[r];}  //
   delete[] prot; delete[] rosetta;delete[] newOrder; //
   /////////////////////////////////////////////////////
   SetNRGmats(FirstOrder,inter);
   out << "We are done with SetNRGmats " << endl;

   int ***FO_energyAcc = FirstOrder->energyAcc;
   double **mPtr0, **mPtr1, **mPtr2, **mPtr3, **mPtr4, *tmpPtr;
//CHANGE by JEFF   double **mPtr0, **mPtr1, **mPtr2, **mPtr3, **mPtr4, **mPtr5, *tmpPtr;
   int *access;
//CHANGE by JEFF   int colPtr, *access;

   double ***NRG0 = FirstOrder->energy[0];
   double ***NRG1 = FirstOrder->energy[1];
   double ***NRG2 = FirstOrder->energy[2];
   double ***NRG3 = FirstOrder->energy[3];
   double ***NRG4 = FirstOrder->energy[4];


   for(int AAcol = 0; AAcol < AAlen; AAcol++)
   {
      //out << "This is AAcol " << AAcol << endl;
      int numNeigh = FirstOrder->AAinfo[AAcol].numNeighbors;
    //  out << "This is numNeigh " << numNeigh << endl;
      FirstOrder->AAinfo[AAcol].pairwise = new double *[numNeigh];
    //  out << "This is AAcol 1: " << AAcol << endl;
      for(int aa2 = 0; aa2 < numNeigh; aa2++)
      {
         FirstOrder->AAinfo[AAcol].pairwise[aa2] = new double [400];
         access = FO_energyAcc[AAcol][aa2];
      //   out << "This is Access0 " << access[0] << endl;
         mPtr0 = NRG0[access[0]];  mPtr1 = NRG1[access[1]];
         mPtr2 = NRG2[access[2]];  mPtr3 = NRG3[access[3]];
         mPtr4 = NRG4[access[4]];
         tmpPtr = FirstOrder->AAinfo[AAcol].pairwise[aa2];
         int ptrCount = 0;
         for(int aa_i = 0; aa_i < 20; aa_i++)
         {
            for(int aa_j = 0; aa_j < 20; aa_j++)
            {
               tmpPtr[ptrCount] = mPtr0[aa_i][aa_j] + mPtr1[aa_i][aa_j] +
                                  mPtr2[aa_i][aa_j] + mPtr3[aa_i][aa_j] +
                                  mPtr4[aa_i][aa_j];
               ptrCount++;
            }
         }
      }
   }
   int zero[5];
   zero[0] = 93; zero[1] = 89; zero[2] = 93; zero[3] = 88; zero[4] = 82;
   for(int i = 0; i < 5; i++){
      for(int j = 0; j < zero[i]; j++){
         for(int k = 0; k < 20; k++){
            delete[] FirstOrder->energy[i][j][k];
         }
         delete[] FirstOrder->energy[i][j];
      }
      delete[] FirstOrder->energy[i];
   }
   delete[] FirstOrder->energy;
   out << "We are done with Interaction " << endl;

   // Check # 2 (The Big One!)
#if LOUD
// Here we do a MASSIVE Check on the system to see if it worked
//   void CheckInteract(FirstOrder, AAlen, out);
// END MASSIVE Check on the System!!!!!!
#endif
}  // End SetUpInteraction Subroutine
//******************************************************************************
int NeighborAccess(int AA, int *c, int **nAccess)
{
   int index;

   if(AA == 1){
      if(c[0] == 0){
         if(c[2] == 0){index = 4;}
         else {index = 5;}
      }
      else{index = nAccess[AA][c[2]];}
   }
   else if(AA == 10){
      if(c[0] == 3){
         if(c[2] == 0){index = 33;}
         else {index = 34;}
      }
      else{index = nAccess[AA][c[2]];}
   }
   else if(AA == 15){
      if(c[0] == 0){
         if(c[2] == 1){index = 44;}
         else {index = 45;}
      }
      else{index = nAccess[AA][c[2]];}
   }
   else{index = nAccess[AA][c[2]];}
   if(index == 65)
   {
      cerr << "WARNING!!!  Neighbor Access out of range!!! " << endl;
      cerr << "For amino acid = " << AA << ", does not match codon = ";
      cerr << c[0] << c[1] << c[2] << endl << "Correct amino acid is = ";
      cerr << Nuc2AATable[c[0]] [c[1]] [c[2]] <<
      " , return incorrect of 65 as safguard!" << endl;
   }
   return index;
} // End NeighborAccess Subroutine
//******************************************************************************
int SolventAccess(double acc)
{
   if(acc < 12.0) {return 0;}
   else if(acc < 36.0) {return 1;}
   else if(acc < 44.0) {return 2;}
   else if(acc < 87.0) {return 3;}
   else {return 4;}
} // End SolventAccess Subroutine
//******************************************************************************
int MatrixHunter(int seqsep, int pair, double **rosetta, double distance)
{
   int zero, preMats, mat;
   if(seqsep == 1)
   {
      if(pair == 0)      {mat = 9;  zero = 92;}
      else if(pair == 1) {return 88;}
      else if(pair == 2) {mat = 12; zero = 92;}
      else if(pair == 3) {mat = 7;  zero = 87;}
      else               {return 81;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[i][pair]){return i;}
      }
      return zero; // If it is not one of the above, then it = 0.0 matrix
   }
   else if(seqsep == 2)
   {
      if(pair == 0)      {mat = 17; preMats = 9;  zero = 92;}
      else if(pair == 1) {mat = 13; preMats = 0;  zero = 88;}
      else if(pair == 2) {mat = 14; preMats = 12; zero = 92;}
      else if(pair == 3) {mat = 16; preMats = 7;  zero = 87;}
      else               {mat = 12; preMats = 0;  zero = 81;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
      return zero; // If it is not one of the above, then it = 0.0 matrix
   }
   else if(seqsep == 3)
   {
      if(pair == 0)      {mat = 11; preMats = 26;}
      else if(pair == 1) {mat = 14; preMats = 13;}
      else if(pair == 2) {mat = 10; preMats = 26;}
      else if(pair == 3) {mat = 10; preMats = 23;}
      else               {mat = 11; preMats = 12;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   else if(seqsep == 4)
   {
      if(pair == 0)      {mat = 8;  preMats = 37;}
      else if(pair == 1) {mat = 11; preMats = 27;}
      else if(pair == 2) {mat = 8;  preMats = 36;}
      else if(pair == 3) {mat = 8;  preMats = 33;}
      else               {mat = 9;  preMats = 23;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   else if(seqsep == 5)
   {
      if(pair == 0)      {mat = 7; preMats = 45;}
      else if(pair == 1) {mat = 8; preMats = 38;}
      else if(pair == 2) {mat = 6; preMats = 44;}
      else if(pair == 3) {mat = 7; preMats = 41;}
      else               {mat = 7; preMats = 32;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   else if(seqsep == 6)
   {
      if(pair == 0)      {mat = 6; preMats = 52;}
      else if(pair == 1) {mat = 7; preMats = 46;}
      else if(pair == 2) {mat = 6; preMats = 50;}
      else if(pair == 3) {mat = 6; preMats = 48;}
      else               {mat = 6; preMats = 39;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   else if(seqsep == 7)
   {
      if(pair == 0)      {mat = 6; preMats = 58;}
      else if(pair == 1) {mat = 6; preMats = 53;}
      else if(pair == 2) {mat = 6; preMats = 56;}
      else if(pair == 3) {mat = 5; preMats = 54;}
      else               {mat = 6; preMats = 45;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   else if(seqsep == 8)
   {
      if(pair == 0)      {mat = 5; preMats = 64;}
      else if(pair == 1) {mat = 5; preMats = 59;}
      else if(pair == 2) {mat = 5; preMats = 62;}
      else if(pair == 3) {mat = 4; preMats = 59;}
      else               {mat = 5; preMats = 51;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   else if(seqsep == 9)
   {
      if(pair == 0)      {mat = 4; preMats = 69;}
      else if(pair == 1) {mat = 4; preMats = 64;}
      else if(pair == 2) {mat = 5; preMats = 67;}
      else if(pair == 3) {mat = 4; preMats = 63;}
      else               {mat = 5; preMats = 56;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   else if(seqsep == 10)
   {
      if(pair == 0)      {mat = 4; preMats = 73;}
      else if(pair == 1) {mat = 4; preMats = 68;}
      else if(pair == 2) {mat = 4; preMats = 72;}
      else if(pair == 3) {mat = 4; preMats = 67;}
      else               {mat = 4; preMats = 61;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   else if(seqsep == 11)
   {
      if(pair == 0)      {mat = 3; preMats = 77;}
      else if(pair == 1) {mat = 4; preMats = 72;}
      else if(pair == 2) {mat = 4; preMats = 76;}
      else if(pair == 3) {mat = 4; preMats = 71;}
      else               {mat = 4; preMats = 65;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   else if(seqsep < 23)
   {
      if(pair == 0)      {mat = 6; preMats = 80;}
      else if(pair == 1) {mat = 6; preMats = 76;}
      else if(pair == 2) {mat = 6; preMats = 80;}
      else if(pair == 3) {mat = 6; preMats = 75;}
      else               {mat = 6; preMats = 69;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   else
   {
      if(pair == 0)      {mat = 6; preMats = 86;}
      else if(pair == 1) {mat = 6; preMats = 82;}
      else if(pair == 2) {mat = 6; preMats = 86;}
      else if(pair == 3) {mat = 6; preMats = 81;}
      else               {mat = 6; preMats = 75;}

      for(int i = 0; i < mat; i++)
      {
         if(distance < rosetta[preMats + i][pair]){return (preMats + i);}
      }
   }
   return 150;

}  // End MatrixHunter Subroutine
//******************************************************************************
void SimSeqNRG(Interaction *FO, int *tmpAAseq, GibbsPart *gibbsSeq, int AAlen)
{
////////////////////////////////////////////////////
// Here gridPoint is the solvent gridpoint, and   //
// seqnum is the  ((pairwise gridpoint * SEQSET)  //
// + what ever seq we are on).  The 0 th point is //
// the solvent term, and the first point is the   //
// pairwise term. To be used later in the program.//
////////////////////////////////////////////////////
   AAsiteInfo *tempInfo; tempInfo = FO->AAinfo;

   int  *DDD, *initAA, numNeigh;
   double s,p, **T2NRG, *T1NRG; s = p = 0.0;
   initAA = tmpAAseq;

   for(int aa1 = 0; aa1 < AAlen; aa1++)
   {
      AAsiteInfo &AAinfoPtr = *tempInfo;
      T2NRG = AAinfoPtr.pairwise;
      numNeigh = AAinfoPtr.numNeighbors;
      DDD = AAinfoPtr.DDDneighbors;
      s += AAinfoPtr.solvent[*initAA];

      for(int aa2 = 0; aa2 < numNeigh; aa2++)
      {
         T1NRG = *T2NRG;
         if(*DDD > aa1){p += T1NRG[(*initAA)*20+tmpAAseq[*DDD]];}
         DDD++; T2NRG++;
      }
      initAA++; tempInfo++;
   }
   gibbsSeq->S = -2.0*s;
   gibbsSeq->P = -2.0*p;
}  // End SimSeqNRG Subroutine
//******************************************************************************
void GibbsSampler(Information &P, Interaction *FO, long *seed, int solv,
                int pair, int A, int C, int G)
{  //NOTE: We only send in the root node sequence
   void SimSeqNRG(Interaction *FO,int *tmpAAseq,GibbsPart *gibbsSeq, int AAlen);

   int *DDD;
// CHANGE by JEFF   int *DDD, *access;
   int killCol, c[3], whichNuc, tempAA, AAcol, initAA;
   double Energy[4], E_total, E_solv, E_pair, whatNRG[4];
   double solv_pt,pair_pt,nuc_pt[4];
   int burn, upperBound, nucLen, count; count = 0;
   int AAseqLen, separation, FO_numNeighAAcol, whatAmino[4], flag;
// CHANGE by JEFF   int AAseqLen, separation, FO_numNeighAAcol, whatAmino[4], flag, nucCount[4];
   AAsiteInfo *tempInfo; tempInfo = FO->AAinfo;

   nucLen = P.len; AAseqLen = P.AAlen; burn = 10;//burn = 1000;
   separation = 10; // This tells us how many seqs we generate between samples

   GibbsInfo[solv][pair][A][C][G] = new GibbsPart[SEQSET];
   GibbsPart *gibbsSeq = GibbsInfo[solv][pair][A][C][G];

   static int *tmpNUCseq;
   static int *tmpAAseq;
   static int z = 0;
   if(z==0)
   {
      tmpNUCseq = new int[nucLen];
      tmpAAseq  = new int[AAseqLen];
      z = 1;
   }

/*
   double ***FOenergy0      = FO->energy[0];
   double ***FOenergy1      = FO->energy[1];
   double ***FOenergy2      = FO->energy[2];
   double ***FOenergy3      = FO->energy[3];
   double ***FOenergy4      = FO->energy[4];
   int    ***FO_energyAcc   = FO->energyAcc;
   int      *FO_numNeigh    = FO->numNeigh;
   int     **FO_DDDneighbor = FO->DDDneighbor;
   double  **FO_solvent     = FO->solvent;
   int      *FO_solvAcc     = FO->solvAcc;
//   int      *PAAseq         = P.AA_seq;
//   int     **Pseq           = P.seq;
*/

   // Here we have a 2000, burn in, and then we take every 10th sample point
   // in the chain.  I think the chains are too related, and this
   upperBound = burn + (separation * SEQSET); // is one way of dealing with it!!!
   int numA,numC,numG,numT,NA,NC,NG,NT;     NA = NC = NG = NT = 0;


   int *t1AA, *t2AA, *t1NUC, *t2NUC;
   t1AA = tmpAAseq; t2AA = P.AA_seq;
   t1NUC = tmpNUCseq; t2NUC = P.seq;
   for(int nuc = 0; nuc < nucLen; nuc++){
      if(nuc < AAseqLen){*t1AA = *t2AA; t1AA++; t2AA++;}
      *t1NUC = *t2NUC;
      switch(*t2NUC){
          case 0: NA++; break;
          case 1: NC++; break;
          case 2: NG++; break;
          case 3: NT++; break;
      }
      t1NUC++; t2NUC++;
   }
   numA = NA; numC = NC; numG = NG; numT = NT;

   solv_pt = S_POINTS[solv];
   pair_pt = P_POINTS[pair];
   nuc_pt[0] = A_POINTS[A];   nuc_pt[1] = C_POINTS[C];
   nuc_pt[2] = G_POINTS[G];   nuc_pt[3] =1.0-(nuc_pt[0]+nuc_pt[1] + nuc_pt[2]);

   for(int seq = 0; seq < upperBound; seq++){
     for(int i = 0; i < nucLen; i++){
       killCol = (int)(nucLen*rnd(seed));//Choose the column we wish to kill
       whichNuc = killCol%3; // calculate the associated AA column
       AAcol = killCol/3;
       AAsiteInfo &AAinfoPtr = tempInfo[AAcol];
       switch(whichNuc){ // We need the correct full codon for this site
         case 0:c[1]=tmpNUCseq[killCol+1];c[2]=tmpNUCseq[killCol+2];break;
         case 1:c[0]=tmpNUCseq[killCol-1];c[2]=tmpNUCseq[killCol+1];break;
         case 2:c[0]=tmpNUCseq[killCol-2];c[1]=tmpNUCseq[killCol-1];break;
       }
       for(int nuc = 0; nuc < 4; nuc++){
         c[whichNuc] = nuc; // This completes the codon
         whatAmino[nuc] = tempAA = Nuc2AATable[c[0]][c[1]][c[2]]; // AA of that nuc
         flag = 0;
         switch(nuc){
         case 1: if(tempAA == whatAmino[0])
                 {
                   flag=1;
                   Energy[1] = nuc_pt[1] * whatNRG[0];
                 } break;
         case 2: if(tempAA == whatAmino[0])
                 {
                   flag=1;
                   Energy[2] = nuc_pt[2] * whatNRG[0];
                 }
                 else{
                   if((!flag)&&(tempAA == whatAmino[1]))
                   {
                     flag = 1;
                     Energy[2] = nuc_pt[2] * whatNRG[1];
                   }
                 } break;
         case 3: if(tempAA == whatAmino[0])
                 {
                   flag=1;
                   Energy[3] = nuc_pt[3] * whatNRG[0];
                 }
                 else if((!flag)&&(tempAA == whatAmino[1]))
                 {
                   flag = 1;
                   Energy[3] = nuc_pt[3] * whatNRG[1];
                 }
                 else{
                   if((!flag)&&(tempAA == whatAmino[2]))
                   {
                     flag = 1;
                     Energy[3] = nuc_pt[3] * whatNRG[2];
                   }
                } break;
         } // end switch statement
         double **T2NRG, *T1NRG;

         if(!flag){
           if(tempAA == 20){whatNRG[nuc] = Energy[nuc] = 0.0;} // If it is a stop codon
           else{
             // Here we calculate the Solvent energy
 	     E_solv = -2.0 * AAinfoPtr.solvent[tempAA];
             // Here we calculate the Pairwise energy
             E_pair = 0.0;
             T2NRG = AAinfoPtr.pairwise;
             FO_numNeighAAcol = AAinfoPtr.numNeighbors;
             DDD = AAinfoPtr.DDDneighbors;
             for(int aa2 = 0; aa2 < FO_numNeighAAcol; aa2++){
                T1NRG = *T2NRG;
                initAA = tmpAAseq[*DDD];
                if(*DDD < AAcol){
                   E_pair += T1NRG[initAA*20+tempAA];
                }
                else{
                   E_pair += T1NRG[tempAA*20+initAA];
                }
                DDD++; T2NRG++;
             }
             E_pair *= -2.0;
             whatNRG[nuc] = exp((E_solv*solv_pt) + (E_pair*pair_pt));
             Energy[nuc] = nuc_pt[nuc] * whatNRG[nuc];
           } // end if not a stop codon
         } // End not flag
       } // end looping through the 4 nucleotide possibilities
       E_total = (Energy[0] + Energy[1] + Energy[2] + Energy[3])*rnd(seed);

       switch(tmpNUCseq[killCol])
       {
          case 0: numA--; break;
          case 1: numC--; break;
          case 2: numG--; break;
          case 3: numT--; break;
       }

       if(E_total < Energy[0]){
         c[whichNuc] = tmpNUCseq[killCol] = 0;   numA++;
         tmpAAseq[AAcol] = whatAmino[0];
       }
       else if(E_total < (Energy[0]+Energy[1])){
         c[whichNuc] = tmpNUCseq[killCol] = 1;   numC++;
         tmpAAseq[AAcol] = whatAmino[1];
       }
         else if(E_total < (Energy[0]+Energy[1]+Energy[2])){
         c[whichNuc] = tmpNUCseq[killCol] = 2;   numG++;
         tmpAAseq[AAcol] = whatAmino[2];
       }
       else{
         c[whichNuc] = tmpNUCseq[killCol] = 3;   numT++;
         tmpAAseq[AAcol] = whatAmino[3];
       }
     } // End making the C->len updates to the sequence
     if((seq >= burn) && (!(seq%separation))){
       //nucCount[0] = nucCount[1] = nucCount[2] = nucCount[3] = 0;
       SimSeqNRG(FO,tmpAAseq,gibbsSeq,AAseqLen);
       //for(int q=0; q < nucLen; q++){nucCount[tmpNUCseq[q]]++;}
       gibbsSeq->A = numA;
       gibbsSeq->C = numC;
       gibbsSeq->G = numG;
       count++; //if(count%1000 == 0){cout << "Another thousand! " << endl;}
       gibbsSeq++;
     }
   } // End sampling all of the sequences in the SEQSET!!!
#if WHERE_ARE_WE
   static int heather = 0;
   if(!(heather%1000))
   {
     cout << "We simulated " << count << " sequences " << heather << " times " << endl;
   }
   heather++;
#endif
   //delete[] tmpNUCseq;   delete[] tmpAAseq;
}
//******************************************************************************
double CalcConst(Bayes *pE, Information &P, Interaction *FO, double &maximum,
                long *seed, int round, ofstream &o)
{   // WARNING: THIS NEEDS UPDATING IN THE CODE TO ACCOUNT FOR BAYES
   void GibbsSampler(Information &P, Interaction *FO, long *seed, int solv,
        int pair, int A, int C, int G);
   static double *tmpSeqProb;
   static int z = 0;
   if(z==0){tmpSeqProb = new double[SEQSET]; z = 1;}

   double sum, solvent, pairwise;
   double sim_s, sim_p, sim_A, sim_C, sim_G, diffA, diffC, diffG, diffT;
   int tempA, tempC, tempG, tempS, tempP, nucLen; nucLen = P.len;
   GibbsPart *tempSeq;

   maximum = -99999999999999.0;
   sum = 0.0;

   switch(round){
   case 3: tempA = pE->A_gridPt;        tempC = pE->C_gridPt;
           tempG = pE->G_gridPt;        tempS = pE->tempS_gridPt;
           tempP = pE->P_gridPt;        solvent = pE->tempParameter;
           pairwise = pE->pairwise;     sim_s = S_POINTS[tempS];
           sim_p = P_POINTS[tempP];     sim_A = A_POINTS[tempA];
           sim_G = G_POINTS[tempG];     sim_C = C_POINTS[tempC];
           diffA = log(pE->MCMCnf[0]/sim_A);
           diffC = log(pE->MCMCnf[1]/sim_C);
           diffG = log(pE->MCMCnf[2]/sim_G);
           diffT = log(pE->MCMCnf[3]/(1.0-(sim_A + sim_C + sim_G))); break;

   case 4: tempA = pE->A_gridPt;        tempC = pE->C_gridPt;
           tempG = pE->G_gridPt;        tempS = pE->S_gridPt;
           tempP = pE->tempP_gridPt;    solvent = pE->solvent;
           pairwise = pE->tempParameter;    sim_s = S_POINTS[tempS];
           sim_p = P_POINTS[tempP];         sim_A = A_POINTS[tempA];
           sim_G = G_POINTS[tempG];         sim_C = C_POINTS[tempC];

           diffA = log(pE->MCMCnf[0]/sim_A);
           diffC = log(pE->MCMCnf[1]/sim_C);
           diffG = log(pE->MCMCnf[2]/sim_G);
           diffT = log(pE->MCMCnf[3]/(1.0-(sim_A + sim_C + sim_G))); break;

   case 5: tempA = pE->tempA_gridPt;    tempC = pE->tempC_gridPt;
           tempG = pE->tempG_gridPt;    tempS = pE->S_gridPt;
           tempP = pE->P_gridPt;        solvent = pE->solvent;
           pairwise = pE->pairwise;     sim_s = S_POINTS[tempS];
           sim_p = P_POINTS[tempP];     sim_A = A_POINTS[tempA];
           sim_G = G_POINTS[tempG];     sim_C = C_POINTS[tempC];

           diffA = log(pE->tempMCMCnf[0]/sim_A);
           diffC = log(pE->tempMCMCnf[1]/sim_C);
           diffG = log(pE->tempMCMCnf[2]/sim_G);
           diffT = log(pE->tempMCMCnf[3]/(1.0-(sim_A + sim_C + sim_G))); break;

   case 9: tempA = pE->A_gridPt;        tempC = pE->C_gridPt;
           tempG = pE->G_gridPt;        tempS = pE->S_gridPt;
           tempP = pE->P_gridPt;        solvent = pE->solvent;
           pairwise = pE->pairwise;     sim_s = S_POINTS[tempS];
           sim_p = P_POINTS[tempP];     sim_A = A_POINTS[tempA];
           sim_G = G_POINTS[tempG];     sim_C = C_POINTS[tempC];

           diffA = log(pE->MCMCnf[0]/sim_A);
           diffC = log(pE->MCMCnf[1]/sim_C);
           diffG = log(pE->MCMCnf[2]/sim_G);
           diffT = log(pE->MCMCnf[3]/(1.0-(sim_A + sim_C + sim_G))); break;
   }

// CHANGE by JEFF   static int doug = 0;

   if(GibbsInfo[tempS][tempP][tempA][tempC][tempG] == NULL)
   {
     GibbsSampler(P, FO, seed, tempS,tempP,tempA,tempC,tempG);
     #if WHERE_ARE_WE
     static int doug = 0;
     if(!(doug%1000))
     {
       o << "We have entered the Gibbs Sampler "  << doug << " times " << endl;
     }
     doug++;
     #endif
   }

   tempSeq = GibbsInfo[tempS][tempP][tempA][tempC][tempG];
   double *tempProbPtr = tmpSeqProb;
   for(int seqnum = 0; seqnum < SEQSET; seqnum++){
      *tempProbPtr =
                tempSeq->S*(solvent - sim_s) + tempSeq->P*(pairwise - sim_p)
                +tempSeq->A*(diffA)+tempSeq->C*(diffC) + tempSeq->G*(diffG)
                +(nucLen-(tempSeq->A+tempSeq->C+tempSeq->G))*(diffT);
      if(*tempProbPtr > maximum){maximum = *tempProbPtr;}
      tempSeq++; tempProbPtr++;
   }
//   cout << "We just calcualted the constant of " << count << " seqs " << endl;
   tempProbPtr = tmpSeqProb;
   for(int i = 0; i < SEQSET; i++)
   {
      sum += exp(*tempProbPtr - maximum);
      tempProbPtr++;
   }

//   P->maxValue = maximum;
   return sum;
}  // End CalcConst Subroutine
//******************************************************************************
void SetNRGholder(Bayes *pE)
{
   double A_MULT, HALF_A, C_MULT, HALF_C, G_MULT, HALF_G;

   A_MULT = (FREQGRIDSIZE-1) / (HIGH_A - LOW_A);
   HALF_A = 0.5/A_MULT;
   C_MULT = (FREQGRIDSIZE-1) / (HIGH_C - LOW_C);
   HALF_C = 0.5/C_MULT;
   G_MULT = (FREQGRIDSIZE-1) / (HIGH_G - LOW_G);
   HALF_G = 0.5/G_MULT;

   pE->S_gridPt = pE->tempS_gridPt =
      ((int) floor( ((pE->solvent - W_S_MIN)*SIZEOFGRID) / (W_S_MAX - W_S_MIN)) );
   pE->P_gridPt = pE->tempP_gridPt =
      ((int) floor( ((pE->pairwise - W_P_MIN)*SIZEOFGRID) / (W_P_MAX - W_P_MIN)) );
   pE->A_gridPt = pE->tempA_gridPt =
      (int) (1.0 +(floor((ADENINE - (LOW_A + HALF_A))*A_MULT)));
   pE->C_gridPt = pE->tempC_gridPt =
      (int) (1.0 +(floor((CYTOSINE - (LOW_C + HALF_C))*C_MULT)));
   pE->G_gridPt = pE->tempG_gridPt =
      (int) (1.0 +(floor((GUANINE - (LOW_G + HALF_G))*G_MULT)));
}  // End SetNRGholder Subroutine
//******************************************************************************
void CalcNeighNRG(Bayes *pE, Interaction *FO,Information &C,Information &P,
                int *tmpAAseq,int*tmpNUCseq)
{
   int NeighborAccess(int AA, int *c, int **nAccess);

   int codon[3], spacer, AAseqLen;
// CHANGE by JEFF   int codon[3], nIndex, spacer, AAseqLen;
   AAseqLen = C.AAlen;
 //  cout << "Inside CalcNeighNRG  " << C.AAlen << endl;
   int *N1, *N2, *AAptr;
   N1 = P.neighAcc;
   N2 = C.neighAcc;
   AAptr = tmpAAseq;


   for(int AAcol = 0; AAcol < AAseqLen; AAcol++)
   {
      spacer = 3*AAcol;
      for(int nuc = 0; nuc < 3; nuc++){codon[nuc] = tmpNUCseq[spacer + nuc];}
      // Obtain the correct row for the Neighbor list in FirstOrder
     // P.neighAcc[AAcol] = C.neighAcc[AAcol] =
       *N1 = *N2 = NeighborAccess(*AAptr, codon, FO->nAccess);
       N1++; N2++; AAptr++;
   }
}  // End CalcNeighNRG Subroutine
//******************************************************************************
void SeqEnergy(Interaction *FO, int *tmpAAseq, double &solvNRG, double &pairNRG,
                        int AAseqLen, ofstream &out)
{
   AAsiteInfo *tempInfo; tempInfo = FO->AAinfo;
   int  *DDD, *initAA, numNeigh;
   double s,p, **T2NRG, *T1NRG; s = p = 0.0;
   initAA = tmpAAseq;

   for(int aa1 = 0; aa1 < AAseqLen; aa1++)
   {
      AAsiteInfo &AAinfoPtr = *tempInfo;
      T2NRG = AAinfoPtr.pairwise;
      numNeigh = AAinfoPtr.numNeighbors;
      DDD = AAinfoPtr.DDDneighbors;
      s += AAinfoPtr.solvent[*initAA];
      for(int aa2 = 0; aa2 < numNeigh; aa2++){
         T1NRG = *T2NRG;
         if(*DDD > aa1){p += T1NRG[((*initAA)*20)+tmpAAseq[*DDD]];}
         DDD++; T2NRG++;
      }
      initAA++; tempInfo++;
   }
   // NOTE: Only the energy associated with the ancestral sequence is kept
   solvNRG = -2.0*s; // Need this for updating w_solv
   pairNRG = -2.0*p; // Need this for updating w_pair
}
//******************************************************************************
void CalcInitSeqNRG(Bayes *pE, Information *C, Information *P, Interaction *FO,
                ofstream &o)
{   // These are constant throughout the calculation
   void SeqEnergy(Interaction *FO, int *tmpAAseq, double &solv, double &pair,
                        int AAseqLen, ofstream &out);
   void CalcNeighNRG(Bayes *pE, Interaction *FO,Information &C,Information &P,
                   int *tmpAAseq, int*tmpNUCseq);


   int nucLen, AAseqLen; nucLen = C[0].len; AAseqLen = C[0].AAlen;
// CHANGE by JEFF   int nucLen, AAseqLen,q[3]; nucLen = C[0].len; AAseqLen = C[0].AAlen;
   int *tmpNUCseq = new int[nucLen];
   int *tmpAAseq  = new int[AAseqLen];
   int numBranch = pE->numBranch;

   double solv, pair;

   int *I1, *I2, *I3, *I4;
   I1 = tmpAAseq;       I2 = C[numBranch].AA_seq;
   I3 = tmpNUCseq;      I4 = C[numBranch].seq;

   for(int col = 0; col < nucLen; col++)
   {
      if(col < AAseqLen){*I1 = *I2; I1++; I2++;}
      *I3 = *I4; I3++; I4++;
   }

   o << "This is the Root Node energy of the sequence " << endl;
   SeqEnergy(FO,tmpAAseq,solv,pair,AAseqLen,o);
   pE->solvNRG = solv;
   pE->pairNRG = pair;
   solv /=(-2.0);
   pair /=(-2.0);
   o << "Total Energy Solvent = " << solv << ", and pairs = ";
   o << pair << " and total = " << (solv + pair) << endl;

   o << "This is the number of nodes that we have " << C[0].numNodes << endl;
   for(int i = 0; i < C[0].numNodes; i++)
   {
      o << endl <<"This is the energy associated with sequence " << i << endl;

      for(int s = 0; s < AAseqLen; s++){tmpAAseq[s] = C[i].AA_seq[s];}
      SeqEnergy(FO,tmpAAseq,solv,pair,AAseqLen,o);
      solv /= (-2.0);
      pair /= (-2.0);
      o << "Total Energy Solvent = " << solv << ", and pairs = ";
      o << pair << " and total = " << (solv + pair) << endl;
   }

   // run through all sites in protein
   for(int i = 0; i < C[0].numBranch; i++)
   {
      I1 = tmpAAseq;       I2 = C[i].parentAAseq;
      I3 = tmpNUCseq;      I4 = C[i].parentSeq;

      for(int col = 0; col < nucLen; col++)
      {
         if(col < AAseqLen){*I1 = *I2; I1++; I2++;}
         *I3 = *I4; I3++; I4++;
      }
      CalcNeighNRG(pE,FO,C[i],P[i],tmpAAseq,tmpNUCseq);
   }
   o << "We are done with CalcInitSeqEnergy " << endl;
   delete[] tmpNUCseq;     delete[] tmpAAseq;
}  // End CalcInitSeqNRG Subroutine
//******************************************************************************
double EasyProb(Information *I, double** ProbMat)
{    // Simple F84 maximum likelihood calculation
        double prob;         prob = 0.0;
        for(int i = 0; i < I->len; i++)
        {
           prob += log(F84_freq[(I->parentSeq[i])]) +
               log(ProbMat[(I->parentSeq[i])][(I->seq[i])]);
        }
        return prob;
}  // End EasyProb Subroutine
//******************************************************************************
double fact(int n) // To be used for the Poisson distribution   FACTORIAL
{
   double f;
   if(n == 0){return 1.0;}
   else
   {
       f = 1.0;
       for(int i = 1; i <= n; i++){f *= i;}
       return f;
   }
}   // End fact Subroutine
//******************************************************************************
double Pois(double rate, int numMuts) // Poisson dist given rates
{
     double fact(int n);
     double temp;
     temp = (exp(-rate)* pow((rate),numMuts))/(fact(numMuts));
     return temp;
}  // End Pois Subroutine
//******************************************************************************
void Case2Num(double rand1, int &nw)
{             // The case where  ng = 0, nw > = 1
     double right;
     int stop, count;  stop = 0; count = 1;

     right = POISSON_W[count]; // For nw = 1

     if(rand1 < right){nw = count;}
     else{
        while(stop == 0)    // Go thru and create intervals for each case nw > 0
        {
           // A check on the system to see when we are getting close to overflow
           if(count > 13)
           {
              cerr << "WARNING IN CASE 2!! NW = "<< count << " rand ";
              cerr << rand1 << " right = " << right << endl;
           } // This should never occur!
           count++;                    // Increase nw by 1
           right += POISSON_W[count];     // Add the resultant to the sum
           if(rand1 < right){
              stop = 1;   // This is a flag that stops the loop from executing
              nw = count; // This copies the result into the returned element
           }
        }
     }
}   // End Case2Num Subroutine
//******************************************************************************
void Case3Num(double rand1, int& ng, int& nw)
{                // case for ng > 0, nw > 0
     double right; right = 0.0;
     int numg, numw, numberg, stop;
     numg = 0; stop = 0;  // min value ng can take : this resets the value
     while (stop == 0)
     {
        numw = 0; // This always starts at zero
        numg++;       // This updates the column that we are on
        numberg = numg;// this also ensures that we do not repeat a diagonal
        while(numw < numberg)  // We calculate the prob of each situation
        {                   // Then create intervals with that probability
           if((numberg > 13) || (numw > 13))
           {
             cerr << "WARNING IN CASE 3. NG =  " << numberg;
             cerr << " NW = " << numw << endl;
             cerr << "This is rand , right " << rand1 << "," << right << endl;
           }
           right += POISSON_G[numberg]*POISSON_W[numw]; // updates right
           if(rand1 < right){
              stop = 1;   // This signals end of while loop
              ng = numberg;
              nw = numw;
              numberg = -1;  //This forces the end of the loop!
           }
           numw ++;
        }
        while(numberg > 0)
        {
           right += POISSON_G[numberg]*POISSON_W[numw]; // updates right
           if(rand1 < right){
              stop = 1;   // This signals end of while loop
              ng = numberg;
              nw = numw;
              numberg = -1;  //This forces the end of the loop!
           }
           numberg --;
        }
     }
}  // End Case3Num Subroutine
//******************************************************************************
void CalcMuts(int col, int& ng, int& nw, Information &I, long* seed,
                ofstream& out)
{
/************************************************
This subroutine calculates the specific number
of substitutions conditional on the starting and
ending nucleotide in the given column.
************************************************/
     double rnd(long*);
     void Case2Num(double rand1, int &nw);
     void Case3Num(double rand1, int& ng, int& nw);
     double rand1, first, second, secondInt;
     rand1 = rnd(seed);
     int *I_seq, *I_parentSeq;
     I_seq = I.seq;
     I_parentSeq = I.parentSeq;

     if(COMPARE[I_parentSeq[col]][I_seq[col]] == 0)  // If SAME
     {
        rand1 *= F84_PROB[I_parentSeq[col]][I_seq[col]];
        first = POISSON_G[0] * POISSON_W[0];
        if(rand1 < first){ng = 0; nw = 0;}
        else{
           rand1 -= first; // rand1 = rand1 - (expG*expW / P_{ij}(T)
           second = POISSON_G[0] * F84_freq[I_seq[col] + 8];
           secondInt = second * (1.0 - POISSON_W[0]);

           if(rand1 < secondInt){
              rand1/=second;  Case2Num(rand1,nw);  ng = 0;
           }
           else{
              rand1 -= secondInt;  rand1 /= F84_freq[I_seq[col]];
              Case3Num(rand1,ng,nw);
           }
        }
#if LOUD
   out << "SAME: for column " << col << ": This is ng, and nw ";
   out << ng << " " << nw << endl;
#endif
     }
     else if(COMPARE[I_parentSeq[col]][I_seq[col]] == 1)  // if TS
     {
        rand1 *= F84_PROB[I_parentSeq[col]][I_seq[col]];
        second = POISSON_G[0] * F84_freq[I_seq[col] + 8];
        secondInt = second * (1.0 - POISSON_W[0]);
        if(rand1 < secondInt){
           rand1/=second;  Case2Num(rand1,nw);  ng = 0;
        }
        else{
           rand1 -= secondInt;  rand1 /= F84_freq[I_seq[col]];
           Case3Num(rand1,ng,nw);
        }
#if LOUD
   out << "TS: for column " << col << ": This is ng, and nw ";
   out << ng << " " << nw << endl;
#endif
     }
     else{  // Else it is a TV   ng > 0 & nw > = 0
        rand1 *= F84_PROB[I_parentSeq[col]][I_seq[col]];
        rand1 /= F84_freq[I_seq[col]];
        Case3Num(rand1,ng,nw);
#if LOUD
   out << "TV: for column " << col << ": This is ng, and nw ";
   out << ng << " " << nw << endl;
#endif
     }
}   // End CalcMuts Subroutine
//******************************************************************************
void StopCheck(Path *P, int col, int &stop)
{// This program checks the three column codon path for stop codons. //
 // If no stop codons, then we can use it in our path calculation.   //
 // We must start from the third position of the codon in order for  //
 // this function to work properly.                                  //
 //////////////////////////////////////////////////////////////////////
   int codon[3], new_AA, zero, first;
   zero = col - 2; first = (col - 1);
   Substitution *subCol0, *subCol1, *subCol2;

   double t[3]; // To hold the times for each column
   subCol0 = P[zero].firstColSub;
   subCol1 = P[first].firstColSub;
   subCol2 = P[col].firstColSub;

   codon[0] = P[zero].nuc;
   codon[1] = P[first].nuc;
   codon[2] = P[col].nuc;
   // Find out what the top AA is in the codon-Path

   while((subCol0 != NULL) || (subCol1 != NULL) || (subCol2 != NULL))
   {
      if(subCol0 != NULL){t[0] = subCol0->time;}
      else {t[0] = 2.0;}
      if(subCol1 != NULL){t[1] = subCol1->time;}
      else {t[1] = 2.0;}
      if(subCol2 != NULL){t[2] = subCol2->time;}
      else {t[2] = 2.0;}

      if(t[0] <= t[1])
      {
         if(t[0] <= t[2])
         {
//            cerr << "SubCol 0 is the min " << endl;
            codon[0] = subCol0->pathNuc;  // Copy nucleotide into the codon
            subCol0 = subCol0->DcolPtr;
         }
         else
         {
//            cerr << "SubCol 2 is the min - 1 " << endl;
            codon[2] = subCol2->pathNuc;  // Copy nucleotide into the codon
            subCol2 = subCol2->DcolPtr;
         }
      }
      else // This implies that  if(t[1] <= t[0])
      {
         if(t[1] <= t[2])
         {
//            cerr << "SubCol 1 is the min " << endl;
            codon[1] = subCol1->pathNuc;  // Copy nucleotide into the codon
            subCol1 = subCol1->DcolPtr;
         }
         else
         {
//            cerr << "SubCol 2 is the min - 2" << endl;
            codon[2] = subCol2->pathNuc;  // Copy nucleotide into the codon
            subCol2 = subCol2->DcolPtr;
         }
      }
      new_AA = Nuc2AATable[codon[0]] [codon[1]] [codon[2]];
      if(new_AA == 20) {stop = 1;}
   } // end while loop
}  // end StopCheck Subroutine
//******************************************************************************
void SetCodonPath(Path *P, int col)
{// This program sets what amino acids are created in the three column   //
 // codon path.  Because we are certain there are no stop codons, we can //
 // simply assign the correct AA's to their substitution.                //
 //////////////////////////////////////////////////////////////////////////
   int current_AA, zero, first; zero = col - 2; first = col - 1;
   double t[3]; // To hold the times for each column
   Substitution *nowSub, *subCol0, *subCol1, *subCol2;
   subCol0 = P[zero].firstColSub;
   subCol1 = P[first].firstColSub;
   subCol2 = P[col].firstColSub;

   int codon[3];
   codon[0] = P[zero].nuc;
   codon[1] = P[first].nuc;
   codon[2] = P[col].nuc;
   //  cerr << "This is codon for col (" << col << ") = ";
   //  cerr << codon[0] << codon[1] << codon[2] << endl;
   // Find out what the top AA is in the codon-Path
   current_AA = Nuc2AATable[codon[0]] [codon[1]] [codon[2]];

   while((subCol0 != NULL) || (subCol1 != NULL) || (subCol2 != NULL)){
      if(subCol0 != NULL){t[0] = subCol0->time;}
      else {t[0] = 2.0;}
      if(subCol1 != NULL){t[1] = subCol1->time;}
      else {t[1] = 2.0;}
      if(subCol2 != NULL){t[2] = subCol2->time;}
      else {t[2] = 2.0;}

      if(t[0] < t[1]){
         if(t[0] <= t[2]){
            codon[0] = subCol0->pathNuc;  // Copy nucleotide into the codon
	    nowSub = subCol0;
            subCol0 = subCol0->DcolPtr;
         }
         else{
            codon[2] = subCol2->pathNuc;  // Copy nucleotide into the codon
	    nowSub = subCol2;
            subCol2 = subCol2->DcolPtr;
         }
      }
      else{ //  This implies that if(t[1] <= t[0])
         if(t[1] <= t[2]){
            codon[1] = subCol1->pathNuc;  // Copy nucleotide into the codon
	    nowSub = subCol1;
            subCol1 = subCol1->DcolPtr;
         }
         else{
            codon[2] = subCol2->pathNuc;  // Copy nucleotide into the codon
	    nowSub = subCol2;
            subCol2 = subCol2->DcolPtr;
         }
      }
      nowSub->priorAA = current_AA; // set the priorAA term
      //Check what the new AA formed is
      current_AA = nowSub->amino = Nuc2AATable[codon[0]] [codon[1]] [codon[2]];
       // set the amino term and Update the current AA in the loop to the new_AA
   } // end while loop
}  // end SetCodonPath Subroutine
//******************************************************************************
void CalcNumMuts(Information &C,Path *cP, long* seed,ofstream& out)
{
/***************************************************
This is the big dog subroutine part of the program.
It calculates the path information initially.  With
the MCMC, this will be the only time we use this,
and all other times we will just update what we have.
****************************************************/
     void CalcMuts(int col, int& ng, int& nw, Information &I, long* seed,
                ofstream& out);
     void CalcTimes(Path *columnPath, long* seed, int colnum, int AAlen,
                      int ng, int *type, int *force);
     void ChoosingNuc(Information &Info, Path *columnPath, int col,
                long* seed, int *type, int *force);
     void InsertFirstColumn(Information &Info, Path *columnPath, int col);
     void InsertColumn(Information &Info,Path *columnPath, int col);
     void PruneColumn(Information &Info, Path *columnPath, int col);
     void DeleteStopColumn(Information &Info, Path *columnPath, int col);
     void StopCheck(Path *P, int col, int &stop);
     void SetCodonPath(Path *P, int col);
/**********************************************************
Here we calculate the specific number of mutations for each
column of the alignment.
Here we calculate the times uniformly over the interval
for each column of the alignment.  With each partially
ordered column, we also calculate the specific nucleotides
in the path.  This is all done in the one for loop below.
**********************************************************/
    int *type = new int[25];
    int *force = new int[25];

    int ng,nw,stop,nucLen,AAseqLen;  stop = 0;
    nucLen = C.len; AAseqLen = C.AAlen;
    for(int col = 0; col < nucLen; col++)
    {
       CalcMuts(col,ng,nw,C,seed,out);
       if((ng+nw) > 0){
          C.totNumSub += cP[col].numsub = ng + nw;
          CalcTimes(cP,seed,col,AAseqLen,ng,type,force);
          ChoosingNuc(C,cP,col,seed,type,force);
          PruneColumn(C,cP,col);
       }
//     out << "This is column =" << col << ", and col % 3 = " << col%3 << endl;

       if((col%3) == 2)
       {
          StopCheck(cP, col, stop);
        //    out << "We are now done with the first set of changes.";
        //    out <<  Stop = " << stop << endl;
          while(stop == 1) //if it is a stop codon, then we must re-choose path
          {
             for(int tmpCol = col - 2; tmpCol <= col; tmpCol++)
             {
                if(cP[tmpCol].numsub > 0)
                { // In order to replace a column, we must delete it first!
                   C.totNumSub -= cP[tmpCol].numsub;
                   DeleteStopColumn(C,cP,tmpCol);
                }
                CalcMuts(tmpCol,ng,nw,C,seed,out);
                if((ng+nw) > 0){
                   C.totNumSub += cP[tmpCol].numsub = ng + nw;
                   CalcTimes(cP,seed,tmpCol,AAseqLen,ng,type,force);
                   ChoosingNuc(C,cP,tmpCol,seed,type,force);
                   PruneColumn(C,cP,tmpCol);
                } // end if statement
             } // end tmpCol loop
             stop = 0;
             StopCheck(cP, col, stop);
          }// This sets up the priorAA and amino in each substituion
          SetCodonPath(cP, col);
       } // end working with that particular codon in the path
    } // end looping through the entire sequence
////////////////////////////////////////
// Now we order the sub events inside //
// of the path structure!             //
////////////////////////////////////////
    delete[] type;   delete[] force;
    //After the first, we can add all the rest in automatically
     // Doug Changes August 2005 in order to counteract initial paths with
     // zero substitutions.
    int q,r; q = 0; r = 0;
    while(r == 0)
    {
       if(q < nucLen)
       {
          if(cP[q].numsub == 0){q++;}
          else
          {
             r = 1;
             InsertFirstColumn(C, cP, q);
             for(int col = q+1; col < nucLen; col++)
             {
                if(cP[col].numsub > 0){InsertColumn(C, cP, col);}
             }
          }
       }
       else{r = 1;}
    }

    /*
    while(cP[q].numsub == 0){q++;}
    InsertFirstColumn(C, cP, q);
    for(int col = q+1; col < nucLen; col++){
       if(cP[col].numsub > 0){InsertColumn(C, cP, col);}
    }
*/
} // End CalcNumMuts Subroutine
//******************************************************************************
void SetSubTimes(Path *columnPath, int colnum, int AAlen, double *orderTime)
{
     // This sets up the first node at the top
     int numsub; numsub = columnPath[colnum].numsub; // Number of subs in colnum
     Substitution *temp;
     Substitution *cSub = new Substitution;
     cSub->UpathPtr    = NULL;     cSub->DpathPtr    = NULL;
     cSub->UcolPtr     = NULL;     cSub->DcolPtr     = NULL;

     // Initialize everything BABY!
     cSub->column  = colnum;
     cSub->time    = orderTime[0];
     // This connects everything to the path
     columnPath[colnum].firstColSub = cSub;
     // This sets up the remaining substitutions in the column
     for(int i = 1; i < numsub; i++)
     {
        temp = cSub;
        cSub = new Substitution;
        cSub->UpathPtr    = NULL;     cSub->DpathPtr    = NULL;
        cSub->DcolPtr     = NULL;

        cSub->column  = colnum;
        cSub->time    = orderTime[i];
        temp->DcolPtr = cSub;
        cSub->UcolPtr = temp;
     }
}// End SetSubTimes Subroutine
//******************************************************************************
void CalcTimes(Path *columnPath, long* seed, int colnum, int AAlen,
                        int ng, int *type, int *force)
{
/**************************************
Given the number of substitutions, this
routine calculates the times uniformly
over the time interval.
**************************************/
      void SetSubTimes(Path *columnPath,int colnum,int AAlen,double *orderTime);
      double rnd(long*);
      int numsub, flag;
      numsub = columnPath[colnum].numsub;    // This is how many column subs

      static double *tempTime;   // Holds the initial times
      static double *orderTime;  // Holds the ordered times
      static int z = 0;
      if(z==0){tempTime = new double[25]; orderTime= new double[25]; z = 1;}

      for(int i = 0; i < numsub; i++){
         tempTime[i] = rnd(seed); // Place the times in the unordered holder
         force[i] = 0;
      }

      double min;
      for(int i = 0; i < numsub; i++) // Bubble sort method (slow)
      {
         min = 100.0;
         int k;
         for(int j = 0; j < numsub; j++)
         {
            if(tempTime[j] < min)
            {
               k = j;
               min = tempTime[j];
            }
         }
         orderTime[i] = tempTime[k]; // Choose min and order it in time
         if(k < ng) {type[i] = 0;}   // It is a general change
	 else {type[i] = 1;}         // It is a within group change
         tempTime[k] = 100.0;  // update it so it can not be chosen again!
      }
      for(int i = (numsub-1); i >= 0; i--)
      {  // if it is the last substitution, then it is forced
         if(i == (numsub-1)){force[i] = 1;}
         else
         {  // if it is a within group, then it is unforced
            if(type[i] == 1){force[i] = 0;}
            else  // it is a general change
            {
               flag = 0;
               // Are there any other general changes after this one?
               for(int j = i+1; j < numsub; j++)
               {
                  if(!type[j]){flag++;}
               }
               // if no, then it is forced to the group of the last nuc
               if(!flag){force[i] = 2;}
               //if yes, then the general change is unforced
               else {force[i] = 0;}
            }
         }
      }
      SetSubTimes(columnPath, colnum, AAlen, orderTime);
//      delete[] tempTime;     // return memory that we borrowed
//      delete[] orderTime;    // return memory that we borrowed
} // End CalcTimes Subroutine
//******************************************************************************
int Intervals(long* seed, int targetNuc, int type, int force)
{  // This calculates the intervals for the prob of the chosen nucleotide
   // The different cases are for the the different situations that occur
       double rnd(long*);
       double rand1;    rand1 = rnd(seed);
       //force = 2 -> general substitutions that are forced to group of TARGET
       if(force == 2){
          if(!(targetNuc%2)){ //Purines
             if(rand1 <= F84_freq[8]) {return 0;}
             else                    {return 2;}
          }
          else
          {  //Pyrimidines
             if(rand1 <= F84_freq[9]) {return 1;}
             else                    {return 3;}
          }
       }
       else
       {  //force = 0 for general: Prob of choice = relative freq of nucleotide
          if(!type)
          {
             if(rand1 <= F84_freq[0])      {return 0;}  // # <= Pi_A
             else if(rand1 <= F84_freq[4]) {return 1;}
             else if(rand1 <= F84_freq[5]) {return 2;}
             else                         {return 3;}    // Pi_G < # <= 1.0
          }
          else
          { //force = 0 for within: prob of choice remains in group only
             if(!(targetNuc%2)){ //Purines
                if(rand1 <= F84_freq[8]) {return 0;}
                else                    {return 2;}
             }
             else
             {  //Pyrimidines
                if(rand1 <= F84_freq[9]) {return 1;}
                else                    {return 3;}
             }
          }
       }
} // End Intervals Subroutine
//******************************************************************************
void ChoosingNuc(Information &Info, Path *columnPath, int col, long* seed,
                int *type, int *force)
{    // This program chooses the specific nucleotide given the circumstances
       int Intervals(long* seed,int targetNuc, int type, int force);

       int tempNuc, colsub;
       int *Info_seq = Info.seq;
       tempNuc = columnPath[col].nuc;
       colsub = columnPath[col].numsub;
       Substitution *temp;
       temp = columnPath[col].firstColSub;

       for(int sub = 0; sub < colsub; sub++)
       {
          if(force[sub] == 1)
          {
             temp->priorNuc = tempNuc;
             temp->pathNuc  = Info_seq[col];
          }
          //if it is a general forced to the group of the target nucleotide
          else if(force[sub] == 2)
          {
             //Here we choose our nucleotide based on the calculated interval
             temp->priorNuc = tempNuc;
             tempNuc=Intervals(seed,Info_seq[col],type[sub],force[sub]);
             // if it is accepted, we put it into our structure
             temp->pathNuc = tempNuc;
          } //end FORCE = 2
          //if the substitution is unforced, then it can be of two types
          else
          {  //The first type is a general unforced change
             if(type[sub] == 0)
             {//In this case the prob of sub depends on the nuc freqs themselves
                temp->priorNuc = tempNuc;
                tempNuc = Intervals(seed,tempNuc,type[sub],force[sub]);
                temp->pathNuc = tempNuc;
             } //end FORCE = 0, SUBTYPE = 0

             else
             {//The last type is the within group unforced substitution
                temp->priorNuc = tempNuc;
                tempNuc = Intervals(seed,tempNuc,type[sub],force[sub]);
                temp->pathNuc = tempNuc;
             } //end FORCE = 0, SUBTYPE = 1
          } //End FORCE = 0
          temp = temp->DcolPtr;
        }// End total substitution loop
} // End ChoosingNuc Subroutine
//******************************************************************************
void PruneColumn(Information &Info, Path *columnPath, int col)
{    // This function removes all non-event changes from a column
   Substitution *tmp, *deltmp, *previous;
   int count, colsub, AAseqLen; count = 0; colsub = columnPath[col].numsub;
   AAseqLen = Info.AAlen;
   tmp      = columnPath[col].firstColSub;
   for(int i = 0; i < colsub; i++)
   {
      if(tmp->priorNuc == tmp->pathNuc)
      {
         previous = tmp->UcolPtr;
         deltmp = tmp;
         tmp = tmp->DcolPtr;
         if(previous == NULL)
         {
            if(tmp == NULL) // if it is the only node left in the column
            {
               columnPath[col].firstColSub = NULL;
               delete[] deltmp;
            }
            else // if it is the first node with subsequent nodes
            {
               columnPath[col].firstColSub = tmp;
               tmp->UcolPtr = NULL;
               delete[] deltmp;
            }
         }
         else
         {
            if(tmp == NULL) // if it's the last node in column, but not the only
            {
               previous->DcolPtr = NULL;
               delete[] deltmp;
            }
            else // if the node has nodes both before and after
            {
               previous->DcolPtr = tmp;
               tmp->UcolPtr = previous;
               delete[] deltmp;
            }
         }
         count++;
      }
      else {tmp = tmp->DcolPtr;}
   }
   Info.totNumSub -= count;
   columnPath[col].numsub -= count;
} // End PruneColumn Subroutine
//******************************************************************************
void InsertFirstColumn(Information &Info, Path *columnPath, int col)
{  // Used only twice: once when path is first created
   //  then used when after we copy everything into the Proposed structure

   Substitution *temp, *next;

   temp = columnPath[col].firstColSub; // This sets temp to first sub in the col
   Info.startPath = temp;            // Temp is also the first sub in path

   // We now want to add the rest of the column subs to the path
   int colsub = columnPath[col].numsub;
   for(int i = 1; i < colsub; i++)
   {
         next = temp->DcolPtr;
         temp->DpathPtr = next;
         next->UpathPtr = temp;
         temp = next;
   }
} // End InsertFirstColumn Subroutine
//******************************************************************************
void InsertColumn(Information &Info, Path *columnPath, int col)
{  // used whenerver we create a new column for the path
  int colsub = columnPath[col].numsub;
  if(colsub > 0)
  {
   Substitution *column, *temp, *next;
   column = columnPath[col].firstColSub; // This is the column to be inserted
   temp = Info.startPath; // Temp is initially the first sub in the path
   next = Info.startPath; // Next will be the substitution following temp

   // We only want to loop through as many times as we have substitutions
   for(int i = 0; i < colsub; i++)
   {
      while(column->time > next->time) // advances until temp < col < next
      {
         if(next->DpathPtr == NULL) // if we add a node to the end of the path
         {
            next->DpathPtr = column; // connects next to column
            column->UpathPtr = next; // connects column to next
            temp = next;             // advances temp to the next
            next = column;   // next = column makes it bounces out of while loop
         }
         else
         {  // advances the temp and next node down the path one by one
            temp = next;
            next = temp->DpathPtr;
         }
      }
      // For case when we add a node that will be the first node of the path
      if(next->UpathPtr == NULL)
      {
         Info.startPath = column;   // connects column as start of the path
         column->DpathPtr = next;    // Makes the old 1st = 2nd sub of the path
         next->UpathPtr = column;    // Connects 2nd node to 1st node of path
         temp = column;              // sets column as the new temp node
         column = column->DcolPtr;   // advances to the next node of the column
      }
      else
      {  // For case where node is in the middle of the path somewhere
         if((next->DpathPtr != NULL) || (column->time < next->time))
         {  // it may be before last, or somewhere in the middle
            column->DpathPtr = next;  // places column before next
            column->UpathPtr = temp;  // places column after temp
            next->UpathPtr = column;  // connects next to column
            temp->DpathPtr = column;  // connects temp to column
            temp = column;            // advances temp to the column node
            column = column->DcolPtr; // advances to the next node of the column
         }
         else
         {  // This advances column down its own column until the next node
            column = column->DcolPtr;
         }
      }  // This continues until the last column node = NULL
   }     // This will then "kick" out of the program
  }
} // End InsertColumn Subroutine
//******************************************************************************
void CopyInitColumn(Path *emptyPath, Path *fullPath, int column, int AAlen)
{    // used only once to copy the entire Current column into the Proposed
	Substitution *next, *temp, *full;
        full = fullPath[column].firstColSub;  // starts with first substitution
	temp = new Substitution;              // creates a new substitution
        temp->UpathPtr    = NULL;     temp->DpathPtr    = NULL;
        temp->UcolPtr     = NULL;     temp->DcolPtr     = NULL;
	emptyPath[column].numsub = fullPath[column].numsub;
        temp->column   = full->column;
	temp->time     = full->time;     // copies all info from full to temp
	temp->priorNuc = full->priorNuc;
	temp->pathNuc  = full->pathNuc;
        temp->amino    = full->amino;
        temp->priorAA  = full->priorAA;

        // We now connect temp to the Path structure in the correct location
	emptyPath[column].firstColSub = temp;
        // We now go through and add the rest of the subs in the column
	for(int i = 1; i < fullPath[column].numsub; i++)
	{
	   next = new Substitution;
           next->UpathPtr    = NULL;     next->DpathPtr    = NULL;
           next->DcolPtr     = NULL;
	   temp->DcolPtr = next;   // connects next to temp
           next->UcolPtr = temp;   // connects temp to next
	   temp = next;            // advances temp
	   full = full->DcolPtr;   // advances full
           temp->column   = full->column;    // copies all of the information
	   temp->time     = full->time;      //  from full to temp
	   temp->priorNuc = full->priorNuc;
	   temp->pathNuc  = full->pathNuc;
           temp->amino    = full->amino;
           temp->priorAA  = full->priorAA;
	}
} // End CopyInitColumn Subroutine
//******************************************************************************
void SetPathNodes(Path* emptyPath,Substitution *prePath, Substitution *postPath,
        Substitution *emptySub, Substitution* &startPath)
{   // This fn sets up the path pointers for the empty path
     Substitution *temp, *downPath;
     // This first part connects the empty node's Upath Pointer
     if(prePath == NULL)
     {   // if the first sub is the first sub in the path as well
        startPath = emptySub;  // connects the empty to the starting sub
        if(postPath != NULL)
        {
           downPath = emptyPath[postPath->column].firstColSub;
           // Use this line only when 2nd sub is in the same col as startPath
           if(downPath->time != postPath->time){downPath = downPath->DcolPtr;}
        }
        else{downPath = NULL;}
     }   // finds the next sub in the path
     else
     {  // temp is the prior path node, but at the top of the col
        temp = emptyPath[prePath->column].firstColSub; // finds correct node
        while(temp->time != prePath->time){temp = temp->DcolPtr;}
        if(temp->DpathPtr != NULL)
        {
          if(temp->DpathPtr->time != emptySub->time){downPath = temp->DpathPtr;}
          else //This is for the unique case when the first two subs are in
          {    // the same column.  Used very rarely!!
             //cout << "This is the rare case " << endl;
             if(postPath == NULL){downPath = NULL;}
             else
             {
                downPath = emptyPath[postPath->column].firstColSub;
                while(downPath->time != postPath->time)
                   {downPath = downPath->DcolPtr;}
             }
          }
        }
        else {downPath = temp->DpathPtr;}
        temp->DpathPtr = emptySub; // connects the path to the empty node
        emptySub->UpathPtr = temp;  // connects the empty to the path
     }
     // This part connects the empty node's Dpath Pointer
     if(downPath != NULL) // If it is the last node, then it points to NULL
     {  // If not, we did a neat trick above to get the node quickly!!!!
        emptySub->DpathPtr = downPath; // Connects the node to the Down path
        downPath->UpathPtr = emptySub; // Connects the path to the node
     }
} // End SetPathNodes Subroutine
//******************************************************************************
void CopyUpdateColumn(Path *emptyPath,Path *fullPath,int col, Information &Info)
{   // all inclusive function that not only copies the intended column, but
    //    uses the path information to quickly order them in the path
    void SetPathNodes(Path* emptyPath, Substitution *prePath,
            Substitution *postPath,Substitution *emptySub,Substitution* &start);

    Substitution *prePath;   // Used as the path sub prior to the current
    Substitution *postPath;  // Used as the path sub after the current
    Substitution *nowPath;   // is the current sub that is being copied
    Substitution *emptySub;  // is the sub that is being filled
    Substitution *next;      // this is the next one in the empty column

    nowPath  = fullPath[col].firstColSub;  // places the first into NOW
    prePath  = nowPath->UpathPtr;         // finds the pre path
    postPath = nowPath->DpathPtr;         // finds the post path
    emptySub = new Substitution;          // creates a new empty
    emptySub->UpathPtr    = NULL;     emptySub->DpathPtr    = NULL;
    emptySub->UcolPtr     = NULL;     emptySub->DcolPtr     = NULL;
    emptyPath[col].firstColSub = emptySub;    // connects empty to path

    // copies the number of subs
    emptyPath[col].numsub = fullPath[col].numsub;
    emptySub->column      = nowPath->column;
    emptySub->time        = nowPath->time;      // This copies all of the info
    emptySub->priorNuc    = nowPath->priorNuc;
    emptySub->pathNuc     = nowPath->pathNuc;

    // now we step through the rest of them and add them one by one
    int colsub = emptyPath[col].numsub;
    if(colsub == 1)
    {
       SetPathNodes(emptyPath,prePath,postPath,emptySub,Info.startPath);
    }
    else
    {
       for(int i = 1; i < colsub; i++)
       {
          next = new Substitution;
          next->UpathPtr    = NULL;     next->DpathPtr    = NULL;
          next->DcolPtr     = NULL;
          emptySub->DcolPtr = next;
          next->UcolPtr     = emptySub;
          nowPath           = nowPath->DcolPtr;
          next->column      = nowPath->column;
          next->time        = nowPath->time;      // This copies all of the info
          next->priorNuc    = nowPath->priorNuc;
          next->pathNuc     = nowPath->pathNuc;
          SetPathNodes(emptyPath,prePath,postPath,emptySub,Info.startPath);
          emptySub = next;
          prePath  = nowPath->UpathPtr;
          postPath = nowPath->DpathPtr;
          if(i == (colsub - 1))
          {
             SetPathNodes(emptyPath,prePath,postPath,emptySub,Info.startPath);
          }
       }
    }
}  // End CopyUpdateColumn Subroutine
//******************************************************************************
void DeleteStopColumn(Information &Info, Path *columnPath, int col)
{  // deletes the chosen Codon Column from the path,
   // We must treat this differently than the normal delete column

  Substitution *column, *tmpColumn;   // to be deleted in most all cases
   // This is to be used to step along the column
   column = columnPath[col].firstColSub;
   int colsub = columnPath[col].numsub;
   for(int i = 0; i < colsub; i++)
   {
      tmpColumn = column;
      column = column->DcolPtr;
      delete[] tmpColumn;
   }
   columnPath[col].firstColSub = NULL; // We must Null out the path as well
   columnPath[col].numsub = 0;//after deleting, it resets the number of subs = 0
}  // End DeleteStopColumn Subroutine
//******************************************************************************
void UpdateStopColumn(Path *loserPath,Path *winnerPath, int killCol,
                        Information &Info)
{
    Substitution *loser;
    Substitution *winner;
    loser  =  loserPath[killCol].firstColSub;
    winner = winnerPath[killCol].firstColSub;
    int colsub = winnerPath[killCol].numsub;
    for(int i = 0; i < colsub; i++)
    {
       loser->priorAA = winner->priorAA;
       loser->amino   = winner->amino;
       loser->genHood = winner->genHood;
       loser->subHood = winner->subHood;
       loser->nrgDiff     = winner->nrgDiff;
       loser->expNRGdiff  = winner->expNRGdiff;
       loser->solvNRGdiff = winner->solvNRGdiff;
       loser  = loser->DcolPtr;     winner = winner->DcolPtr;
    }
} // End UpdateStopColumn Subroutine
//******************************************************************************
void DeleteColumn(Information &Info, Path *columnPath, int col)
{  // deletes the chosen column from the path, and connects the path nodes
        //  around the now missing node

   Substitution *column;
   Substitution *tmpColumn;   // to be deleted in most all cases
   Substitution *temp;
   Substitution *next;

   // This is to be used to step along the column
   column = columnPath[col].firstColSub;
// CHANGE by JEFF   int AAseqLen = Info.AAlen;
   int colsub = columnPath[col].numsub;
   for(int i = 0; i < colsub; i++)
   {
      temp = column->UpathPtr; // Temp is substitution before column in the path
      next = column->DpathPtr; // next is substitution after column in the path

      if(temp == NULL) // If it is the first node in the entire path
      {
        if(next == NULL) // This is for rare case of path with one substitution
        {
         Info.startPath = NULL;
         delete[] column;
        }
        else
        {
         Info.startPath = next; // then the next becomes the first substitution
         next->UpathPtr = NULL;  // Then we must point the UP to NULL
         tmpColumn = column;     // Setting up so we can DELETE
         column = column->DcolPtr; // move the column down the column pointer
         delete[] tmpColumn;       // DELETE the temporary node
        }
      }
      else if(next == NULL) // If next is the last node in the path
      {                     // This means that temp is the last node in the path
         temp->DpathPtr = NULL;
         delete[] column;     // Now we can delete the column sub safely
      }
      else   // if it is somewhere in the middle of the path
      {
         temp->DpathPtr = next;  // Then we connect the temp to next in the path
         next->UpathPtr = temp;  // Then we connect the next to temp in the path
         tmpColumn = column;     // Like above, we put in the temporary
         column = column->DcolPtr; // Advance the column to the next node
         delete[] tmpColumn;       // DELETE the temporary node
      }
   }
   columnPath[col].firstColSub = NULL; // We must Null out the path as well
   columnPath[col].numsub = 0;//after deleting, it resets the number of subs = 0
   if(Info.totNumSub == 0){Info.startPath = NULL;}
}  // End DeleteColumn Subroutine
//******************************************************************************
double MetroParameter(double k1, double k2, double interval)
{
/****************************************************
This subroutine calculates the metropolis hastings
portion of the Kappa update.  Because the prior on
kappa is uniform, we just take the difference in the
kappa values, and multiply by the height of the
rectangle.  These will be the numerator and
denominator in the calculation
*****************************************************/
    double tmp;
    if((k1 - k2) > 0.0) {tmp = (k1 - k2)/interval;}
    else {tmp = (k2 - k1)/interval;}
    return tmp;
}  // End MetroParameter Subroutine
//******************************************************************************
void NewParameter(long *seed, double& metro_top, double& metro_bot,
          double parameter,   double& tempParam, double maxValue,
          double minValue,    double delta)
{
     double rnd(long*);

     double L_theta, H_theta;

     // We use the following to calculate metro_bot in MCMC routine
     if(delta < (maxValue - parameter)) {H_theta = delta;}
     else {H_theta=(maxValue - parameter);}

     if(delta < (parameter - minValue)){L_theta = delta;}
     else {L_theta = (parameter - minValue);}

     tempParam = parameter - L_theta + rnd(seed)*(H_theta + L_theta);
     // After choosing a new Kappa, we find the prob of the New given the old
//     metro_bot = MetroParameter(parameter, oldParameter, (H_theta + L_theta));
     metro_bot = 1.0 /(H_theta + L_theta);

     // We use the following to calculate metro_top in MCMC routine
     if(delta < (maxValue - tempParam)){H_theta = delta;}
     else {H_theta=(maxValue - tempParam);}

     if(delta < (tempParam - minValue)) {L_theta = delta;}
     else {L_theta = (tempParam - minValue);}

     // After finding the new, we can find the prob of the old given the New
//     metro_top = MetroParameter(parameter, oldParameter, (H_theta + L_theta));
     metro_top = 1.0 /(H_theta + L_theta);
}  // End NewParameterSubroutine
//******************************************************************************
double Metro_Path(Path *columnPath, int column)
{  ////////////////////////////////////////////////////////
   // This subroutine calculates the metropolis hastings //
   // portion of the Path update.  The columns are       //
   // independent, thus we test the original column path //
   // and then the updated column path.  These will be   //
   // the numerator and denominator in the calculation   //
   ////////////////////////////////////////////////////////
     double tmp, W, time;  tmp = 0.0;
     Substitution *temp;  temp = columnPath[column].firstColSub;
     Substitution *next;
     int colsub = columnPath[column].numsub;
     for(int sub = 0; sub < colsub; sub++)
     {
        W = F84_RATE[temp->priorNuc][temp->priorNuc];
        if(!sub) {time = temp->time;}  // T_1 = t_1
        else
        {
           next = temp->UcolPtr;
           time = temp->time - next->time;
        }           // t_i = t_i - t_{i-1}
        if(F84_RATE[temp->priorNuc][temp->pathNuc] < 0.0)
        {
           cerr << "IN Metro Path rate matrix has negative entry " << endl;
        }
        tmp +=  W*time + log(F84_RATE[temp->priorNuc][temp->pathNuc]);
        if(sub < (colsub - 1)){temp = temp->DcolPtr;}
     }
     W = F84_RATE[temp->pathNuc][temp->pathNuc];

     time = 1.0 - temp->time; // T - T_n
     tmp += W * time; // ln ( e^(-WT))
     return tmp;
}  // End Metro_Path Subroutine
//******************************************************************************
int NewHistory(Information &P, Path *cP, Path *pP, long *seed, double &metro_top,
	double &metro_bot, int& stop, ofstream &out)
{
     double rnd(long*);
     void CalcMuts(int col, int& ng, int& nw, Information &I, long* seed,
                ofstream& out);
     void CalcTimes(Path *columnPath, long* seed, int colnum, int AAlen,
                      int ng, int *type, int *force);
     void ChoosingNuc(Information &Info, Path *columnPath, int col,
                long* seed, int *type, int *force);
     void PruneColumn(Information &Info, Path *columnPath, int col);
     void InsertFirstColumn(Information &Info, Path *columnPath, int col);
     void InsertColumn(Information &Info, Path *columnPath, int col);
     void DeleteColumn(Information &Info, Path *columnPath, int col);
     double Metro_Path(Path *columnPath, int col);
     void StopCheck(Path *P, int col, int &stop);
     void SetCodonPath(Path *P, int col);

     int *P_seq, *P_parentSeq;
     P_seq = P.seq;
     P_parentSeq = P.parentSeq;
     int killCol, p_numsub;
     // Fast method for choosing the correct column to delete
     killCol = (int) (P.len * rnd(seed));
     static int z = 0;
     static int *type;
     static int *force;
     if(z == 0)
     {
        type  = new int[25];
        force = new int[25];
        z = 1;
     }

     p_numsub = pP[killCol].numsub;
     if(p_numsub > 0)
     {
        // If this column has substitutions, subtract that number from the total
        P.totNumSub -= p_numsub;
        // delete this column in P path, making it clear for a new column
        DeleteColumn(P,pP,killCol);
     }
     int ng, nw;       ng = nw = 0;    // samples ng and nw from dist.
     CalcMuts(killCol,ng,nw,P,seed,out);
     // if ng + nw = 0, then we would be wasting our times here
     int zeroflag = 0;
     if(P.totNumSub == 0){zeroflag = 1;}
     if((ng+nw) > 0)
     {
        int subs = ng + nw;
        pP[killCol].numsub = subs;  // set numsub at top of col path
        P.totNumSub += subs; // add to total
        CalcTimes(pP,seed,killCol,P.AAlen,ng,type,force);
        ChoosingNuc(P,pP,killCol,seed,type,force);
        PruneColumn(P,pP,killCol);
        if(zeroflag == 1){InsertFirstColumn(P,pP,killCol);}
        else{InsertColumn(P,pP, killCol);}
     }

     // Now we must check for the creation of a stop codon in the mix
     int tmpCol; // if so, stop will equal 1, and we have to deal with it.
     stop = 0;  // We initially think that it will be zero
     int killMod3 = killCol%3;
     if(!killMod3)
     {
        tmpCol = killCol+2; // We must start from position 3 of the codon
        StopCheck(pP,tmpCol, stop);
        if(!stop){SetCodonPath(pP, tmpCol);}
     }
     else if(killMod3 == 1)
     {
        tmpCol = killCol+1; // We must start from position 3 of the codon
        StopCheck(pP,tmpCol, stop);
        if(!stop){SetCodonPath(pP, tmpCol);}
     }
     else
     {
        StopCheck(pP,killCol, stop);
        if(!stop){SetCodonPath(pP, killCol);}
     }

     // This is the probability of the P column given the original.
     // This takes care of no subs in the original col path
     if(stop == 0)
     {
       if(cP[killCol].numsub == 0)
       {  // ln ( e^(-WT))     or substitutions to oneself
          metro_top = F84_RATE[P_parentSeq[killCol]][P_seq[killCol]];
       } // This takes care of no subs in the new col path
       else{metro_top = Metro_Path(cP, killCol);}

       // If proposed path has a diff # of subs, then we must recalc Rate Matrix
       if(pP[killCol].numsub == 0)
       {  // ln ( e^(-WT))      or substitutions to oneself
          metro_bot = F84_RATE[P_parentSeq[killCol]][P_seq[killCol]];
       }
       else{metro_bot = Metro_Path(pP, killCol);}
     } // End if (stop = 0) statement
     else{metro_bot = metro_top = 0.0;}// This is case for existence of stop
     return killCol;
}  // End NewHistory Subroutine
/*******************************************************************************
double SpecificRateAway(Information *Info, Substitution *temp)
      // calculates the specific rate away for substitution in question
{
        return Info->F84Rate[temp->priorNuc][temp->pathNuc];
}
//******************************************************************************
double GeneralRateAway(int *sequence, Information *Info)
       // General rate away from the current sequence in question
{
        double tmp; tmp = 0.0;

        for(int i = 0; i < Info->len; i++) // run through each col and sum over
        {                                // all possibilities
           tmp += -1*Info->F84Rate[sequence[i]][sequence[i]];
        }
        return tmp;
}
*******************************************************************************/
//******************************************************************************
double FastSpecificRateAway(Bayes *pE,Substitution *kill,int whichBranch,int round)
{ // This function is a streamlined version of the Specific Rate away. //
  // I now save the energies associated with the substitution, and     //
  // Kappa, I use the exponential of the energy away, but for          //
  // Omega, I use just the sum of the energies and use omega to get    //
  //              the exponential (according to the formula)           //
  ///////////////////////////////////////////////////////////////////////
   void PrintPath(Substitution *start, ofstream& o);
   double rate, NRG_diff,K, U, W, S, P;
   double *I_nf = pE->MCMCnf;
   if(MULTI_K){K = pE->kappa[whichBranch];}
   else{K = pE->kappa[0];}
   if(MULTI_W){W = pE->omega[whichBranch];}
   else{W = pE->omega[0];}
   U = pE->rate[whichBranch];
   S = pE->solvent;  P = pE->pairwise;

   switch(round)
   {
      case 0: K = pE->tempParameter; break;
      case 1: U = pE->tempParameter; break;
      case 2: W = pE->tempParameter; break;
      case 3: S = pE->tempParameter; break;
      case 4: P = pE->tempParameter; break;
      case 5: I_nf = pE->tempMCMCnf; break;
   }

   if(kill->amino == kill->priorAA) // Substitution is synonymous
   {
      if(((kill->priorNuc + kill->pathNuc) % 2) == 0)
      {  // The substitution is caused by a transition
         rate = (I_nf[kill->pathNuc]) * K;
      } // The substitution is caused by a transversion
      else{rate = I_nf[kill->pathNuc];}
      kill->solvNRGdiff = kill->nrgDiff = 0.0;
      kill->expNRGdiff = 1.0;
   } // End the synonymous if statment
   else // Substitution is non - synonymous
   {  // the if statement is for iterations involving kappa
      // the else statement is for iterations involving omega

      if((round == 3)||(round == 4)){  // if updating w_solv or w_pair
        kill->expNRGdiff = NRG_diff =
        exp(kill->solvNRGdiff * S + kill->nrgDiff * P);
      }
      else{NRG_diff = kill->expNRGdiff;}
      // Now we calculate the specific rate away from this site!
      if(((kill->priorNuc + kill->pathNuc) % 2) == 0)
      {  // The substitution is caused by a transition
         rate = (I_nf[kill->pathNuc]) * K * W * (NRG_diff);
      }
      else // The substitution is caused by a transversion
      {
         rate = (I_nf[kill->pathNuc]) * W * (NRG_diff);
      }
   } // End nonsynonymous else statement
   if(rate <= 0.0)
   {
      ofstream o("junk.out");
      o << "Rate = " << rate;
      o << ", Prior nuc = " << kill->priorNuc;
      o << ", chosen = " << kill->pathNuc;
      o << " ,Prior AA  = " << kill->priorAA;
      o  << ", chosen = " << kill->amino << endl;
      o << " N freq = " << I_nf[kill->pathNuc];
      o << ", Kappa = " << K << " , omega " << W << endl;
      o << ", NRG_diff = " << NRG_diff << " and " << kill->expNRGdiff;
      o << " nrg = " << kill->nrgDiff << endl;
      o << " e^nrg = " << kill->expNRGdiff;
      o << " col = " << kill->column;
      o << " time = " << kill->time;
      o << "  round " << round << endl;
      cerr << "Prior nuc = " << kill->priorNuc;
      cerr << ", chosen = " << kill->pathNuc;
      cerr << " ,Prior AA  = " << kill->priorAA;
      cerr  << ", chosen = " << kill->amino << endl;
      cerr << " N freq = " << I_nf[kill->pathNuc];
      cerr << ", Kappa = " << K << " , omega " << W << endl;
      cerr << ", NRG_diff = " << NRG_diff << " and " << kill->expNRGdiff;
      cerr << " nrg = " << kill->nrgDiff << endl;
      cerr << " e^nrg = " << kill->expNRGdiff;
      cerr << " col = " << kill->column;
      cerr << "  round " << round << endl;
      PrintPath(kill, o);
   }

   return U * rate;
}  // End FastSpecificRateAway Subroutine
//******************************************************************************
double SlowSpecificRateAway(Bayes *pE, Interaction * FO, Substitution *kill,
                int whichBranch, int* tmpAAseq)
{ // This function not only calculates the energy of the site,     //
  // but it also updates every site along the protein accordingly. //
  // In a sense, it initializes the next rate away calculation.    //
  // It also updates the saved energy differences and exponential  //
  // of the energy difference in order to speed up the parameter   //
  // update procedure.                                             //
  ///////////////////////////////////////////////////////////////////

   double *I_nf = pE->MCMCnf;
   double rate, K, W;
   int killAACol = kill->column / 3;
   if(MULTI_K){K = pE->kappa[whichBranch];}
   else{K = pE->kappa[0];}
   if(MULTI_W){W = pE->omega[whichBranch];}
   else{W = pE->omega[0];}

   AAsiteInfo &AAinfoPtr = FO->AAinfo[killAACol];

   if(kill->amino == kill->priorAA) // Substitution is synonymous
   {
      if(((kill->priorNuc + kill->pathNuc) % 2) == 0)
      {  // The substitution is caused by a transition
         rate = (I_nf[kill->pathNuc]) * K;
      } // The substitution is caused by a transversion
      else{rate = I_nf[kill->pathNuc];}
      kill->solvNRGdiff = kill->nrgDiff = 0.0;
      kill->expNRGdiff = 1.0;
   } // End the synonymous if statment
   else // Substitution is non - synonymous
   {
// CHANGE by JEFF      int sAcc, targetAA, prior, nowAA;
      int targetAA, prior, nowAA;
      double NRG_diff;
      prior = kill->priorAA;
      nowAA = kill->amino;
      NRG_diff = kill->nrgDiff = kill->solvNRGdiff = 0.0;
      #if SIM_W_S
         kill->solvNRGdiff= AAinfoPtr.solvent[prior] - AAinfoPtr.solvent[nowAA];
      #endif
      // The below calculates the pairwise potential, now we add the
      // contact and solvent energies to the total energy of this site.
      // Calculate the NRG of the affected site and update all other sites
#if SIM_W_P
      double **T2NRG, *T1NRG;
      int  *DDD = AAinfoPtr.DDDneighbors;
      int FO_numneigh = AAinfoPtr.numNeighbors;
      T2NRG = AAinfoPtr.pairwise;
      for(int AAcol = 0; AAcol < FO_numneigh; AAcol++){
          T1NRG = *T2NRG;
          // Calculate the energy of the actual change that occurred
          // Find & update the affect of change on all sites in the protein
          // Notice that this is arranged as BEFORE - AFTER
          targetAA = tmpAAseq[*DDD];
          if(*DDD < killAACol){
             NRG_diff += (T1NRG[targetAA*20+prior] -  T1NRG[targetAA*20+nowAA]);
          }
          else{
             NRG_diff += (T1NRG[prior*20+targetAA] - T1NRG[nowAA*20+targetAA]);
          }
          DDD++; T2NRG++;
      } // end looping through all of the 3D nearest neighbors
      kill->nrgDiff = NRG_diff;
#endif

      kill->expNRGdiff = NRG_diff = 1.0;
#if (SIM_W_S || SIM_W_P)
      kill->expNRGdiff = NRG_diff =
      exp(kill->solvNRGdiff * pE->solvent + kill->nrgDiff*pE->pairwise);
#endif
      // Now we calculate the specific rate away from this site!
      if(((kill->priorNuc + kill->pathNuc) % 2) == 0)
      {  // The substitution is caused by a transition
         rate = (I_nf[kill->pathNuc])* K * W *(NRG_diff);
      }
      else // The substitution is caused by a transversion
      {
         rate = (I_nf[kill->pathNuc])* W *(NRG_diff);
      }
   } // End nonsynonymous else statement
   return (pE->rate[whichBranch])*(rate);
}  // End SlowSpecificRateAway Subroutine
//******************************************************************************
double SetRateAway(Bayes *pE, Interaction *FO, double *rateAway,int *tmpAAseq,
		int whichBranch, int *siteNeigh, int AAseqLen)
{  // The general rate away is NOT from kill, but from the substitution prior //
   //  to kill.  We name this substitution UpSub. The specific rate away is   //
   //  determined by kill.  This also updates the energies in the NRGaway. We //
   // use this to seriously speed our calculations in the FastGeneralRate!    //
   /////////////////////////////////////////////////////////////////////////////
     double SiteNRGCalc(Interaction *FO, double &sNRG, int AAcol,
                int testAA, int *tmpAAseq);
     //ofstream o("zzzz.out");
     int *nIndex, siteAA, siteNuc, siteKappa, *currentAA;
     double siteNRG, sNRG, rate, totalRate, siteRate;  totalRate = 0.0;

     int *aaRowPtr, *nucRowPtr, *kappaRowPtr;


     double K,W,S,P,U;    U = pE->rate[whichBranch];
     if(MULTI_K){K = pE->kappa[whichBranch];}
     else{K = pE->kappa[0];}
     if(MULTI_W){W = pE->omega[whichBranch];}
     else{W = pE->omega[0];}

     S = pE->solvent; P = pE->pairwise;
     double *I_nf = pE->MCMCnf;
     int **FO_neighList = FO->neighList;
     double *RA = rateAway;
     currentAA = tmpAAseq;
     nIndex = siteNeigh;
     // run through all sites in protein
    for(int AAcol = 0; AAcol < AAseqLen; AAcol++){
       siteRate = 0.0;
       // Obtain the correct row from Neighbor list
       // run through the 9 nearest neighbors
       aaRowPtr = FO_neighList[*nIndex];
       nucRowPtr = aaRowPtr + 9;
       kappaRowPtr = nucRowPtr + 9;
       for(int neighbor = 0; neighbor < 9; neighbor++)
       {
          rate = 0.0; //re-initialize the rate variable
          siteAA    = *aaRowPtr;
          siteNuc   = *nucRowPtr;
          siteKappa = *kappaRowPtr;
          if(siteAA == *currentAA) // If the change is a synonymous one
          {
             if(siteKappa == 1){rate = I_nf[siteNuc];}
             else{rate = I_nf[siteNuc] * K;}
          }
          else  // it is a nonsynonymous change
          {
             if(siteAA != 20) // if not a stop
             {
                siteNRG = 1.0;
                #if (SIM_W_S || SIM_W_P)
                siteNRG = SiteNRGCalc(FO,sNRG,AAcol,siteAA,tmpAAseq);
                siteNRG = exp(sNRG*S + siteNRG*P);
                #endif
                if(siteKappa == 1){rate = I_nf[siteNuc] * W * siteNRG;}
                else{rate = I_nf[siteNuc] * K * W * siteNRG;}
             }
          }

          siteRate += rate;
          //update the pointers to the structure to advance to the next entry
          aaRowPtr++; nucRowPtr++; kappaRowPtr++;
       } // end running through all of the nearest neighbors
       *RA = siteRate; RA++;
       //o << "This is the rate for column " << AAcol << " now  " << rateAway[AAcol] << endl;
       totalRate += siteRate;
       nIndex++; currentAA++;
    } // End running through the entire sequence

    return U * totalRate;
}  // End SetRateAway Subroutine
//******************************************************************************
double FirstSetRateAway(Bayes *pE, Interaction *FO, double *rateAway,int *tmpAAseq,
              int whichBranch, int *siteNeigh, int AAseqLen, double **AllRates)
{
     // This was written to handle the Kappa parameter
     double SiteNRGCalc(Interaction *FO, double &sNRG, int AAcol,
                int testAA, int *tmpAAseq);

     int *nIndex, siteAA, siteNuc, siteKappa, *currentAA;
     double *RA, siteNRG, sNRG, rate, totalRate, siteRate;  totalRate = 0.0;

     int *aaRowPtr, *nucRowPtr, *kappaRowPtr;
     RA = rateAway;

     double K,W,S,P,U;    U = pE->rate[whichBranch];
     K = pE->tempParameter;
//     if(MULTI_K){K = pE->kappa[whichBranch];}
//     else{K = pE->kappa[0];}
     if(MULTI_W){W = pE->omega[whichBranch];}
     else{W = pE->omega[0];}

     S = pE->solvent; P = pE->pairwise;
     double *I_nf = pE->MCMCnf;
     int **FO_neighList = FO->neighList;

     double *solvPtr, *pairwisePtr, *expPtr;
     currentAA = tmpAAseq;
     nIndex = siteNeigh;

     // run through all sites in protein
    for(int AAcol = 0; AAcol < AAseqLen; AAcol++)
    {
       solvPtr = *AllRates;
       pairwisePtr = solvPtr + 9;
       expPtr = pairwisePtr + 9;
       siteRate = 0.0;
       // Obtain the correct row from Neighbor list
       // run through the 9 nearest neighbors
       aaRowPtr = FO_neighList[*nIndex];
       nucRowPtr = aaRowPtr + 9;
       kappaRowPtr = nucRowPtr + 9;
       for(int neighbor = 0; neighbor < 9; neighbor++)
       {
          rate = 0.0; //re-initialize the rate variable
          siteAA    = *aaRowPtr;
          siteNuc   = *nucRowPtr;
          siteKappa = *kappaRowPtr;
          if(siteAA == *currentAA) // If the change is a synonymous one
          {
             if(siteKappa == 1){rate = I_nf[siteNuc];}
             else{rate = I_nf[siteNuc] * K;}
          }
          else  // it is a nonsynonymous change
          {
             if(siteAA != 20) // if not a stop
             {
                siteNRG = 1.0;
                #if (SIM_W_S || SIM_W_P)
                *pairwisePtr = siteNRG = SiteNRGCalc(FO,sNRG,AAcol,siteAA,tmpAAseq);
                *solvPtr = sNRG;
                *expPtr = siteNRG = exp(sNRG*S + siteNRG*P);
                #endif
                if(siteKappa == 1){rate = I_nf[siteNuc] * W * siteNRG;}
                else{rate = I_nf[siteNuc] * K * W * siteNRG;}
               //o << "This is rate " << rate << " and exp Ptr = " << *expPtr << endl;
             }
          }
//          o << "This is rate " << rate << endl;
          siteRate += rate;
          //update the pointers to the structure to advance to the next entry
          solvPtr++; pairwisePtr++; expPtr++; //siteNeigh++;
          aaRowPtr++;  nucRowPtr++;   kappaRowPtr++; //tmpAAseq++;
       } // end running through all of the nearest neighbors
       *RA = siteRate;  RA++;
       //o << "This is the rate for column " << AAcol << " now  " << rateAway[AAcol] << endl;
       totalRate += siteRate;
       AllRates++; nIndex++;  currentAA++;

    } // End running through the entire sequence
   // o << endl << endl;
    return U * totalRate;
}  // End FirstSetRateAway Subroutine
//******************************************************************************
double FastSetRateAway(Bayes *pE, Interaction *FO, double *rateAway,int *tmpAAseq,
    int whichBranch, int *siteNeigh, int AAseqLen, double **AllRates, int round)
{
     // This was written to handle the W, U, S, P parameters
     double SiteNRGCalc(Interaction *FO, double &sNRG, int AAcol,
                int testAA, int *tmpAAseq);

     int *nIndex, siteAA, siteNuc, siteKappa, *currentAA;
     double *RA, siteNRG, rate, totalRate, siteRate;  totalRate = 0.0;
// CHANGE by JEFF     double *RA, siteNRG, sNRG, rate, totalRate, siteRate;  totalRate = 0.0;
     int *aaRowPtr, *nucRowPtr, *kappaRowPtr;
     RA = rateAway;
     currentAA = tmpAAseq;
     nIndex = siteNeigh;

     double K,W,S,P,U;    U = pE->rate[whichBranch];
     if(MULTI_K){K = pE->kappa[whichBranch];}
     else{K = pE->kappa[0];}
     if(MULTI_W){W = pE->omega[whichBranch];}
     else{W = pE->omega[0];}

     S = pE->solvent; P = pE->pairwise;
     double *I_nf = pE->MCMCnf;
     int **FO_neighList = FO->neighList;

     double *solvPtr, *pairwisePtr, *expPtr;

     switch(round)
     {
        case 0: K = pE->tempParameter; break;
        case 1: U = pE->tempParameter; break;
        case 2: W = pE->tempParameter; break;
        case 3: S = pE->tempParameter; break;
        case 4: P = pE->tempParameter; break;
        case 5: I_nf = pE->tempMCMCnf; break;
     }
     // run through all sites in protein
    for(int AAcol = 0; AAcol < AAseqLen; AAcol++){
       solvPtr = *AllRates;
       pairwisePtr = solvPtr + 9;
       expPtr  = pairwisePtr + 9;
       siteRate = 0.0;
       // Obtain the correct row from Neighbor list
       // run through the 9 nearest neighbors
       aaRowPtr = FO_neighList[*nIndex];
       nucRowPtr = aaRowPtr + 9;
       kappaRowPtr = nucRowPtr + 9;
       for(int neighbor = 0; neighbor < 9; neighbor++)
       {
          rate = 0.0; //re-initialize the rate variable
          siteAA    = *aaRowPtr;
          siteNuc   = *nucRowPtr;
          siteKappa = *kappaRowPtr;
          if(siteAA == *currentAA) // If the change is a synonymous one
          {
             if(siteKappa == 1){rate = I_nf[siteNuc];}
             else{rate = I_nf[siteNuc] * K;}
          }
          else  // it is a nonsynonymous change
          {
             if(siteAA != 20) // if not a stop
             {
                siteNRG = 1.0;
                #if (SIM_W_S || SIM_W_P)
                if((round == 3) || (round == 4))
                {
                   *expPtr = siteNRG = exp((*solvPtr)*S + (*pairwisePtr)*P);
                }
                else{siteNRG = *expPtr;}
                #endif
                if(siteKappa == 1){rate = I_nf[siteNuc] * W * siteNRG;}
                else{rate = I_nf[siteNuc] * K * W * siteNRG;}
             }
          }

          siteRate += rate;
          //update the pointers to the structure to advance to the next entry
         solvPtr++; pairwisePtr++; expPtr++; //siteNeigh++;
         aaRowPtr++;  nucRowPtr++;   kappaRowPtr++; //tmpAAseq++;
       } // end running through all of the nearest neighbors
       *RA = siteRate; RA++;
       //o << "This is the rate for column " << AAcol << " now  " << rateAway[AAcol] << endl;
       totalRate += siteRate;
       AllRates++; nIndex++; currentAA++;

    } // End running through the entire sequence
    //o << endl << endl;
    return U * totalRate;
}  // End FastSetRateAway Subroutine
//******************************************************************************
void FastCalcSingleSub(Bayes *pE, Interaction *FO, double *rateAway,int *tmpAAseq,
     int *siteNeigh,int whichBranch, int AAcol, double &totalRate, double *solvPtr,
     int round)
{
     // This is for calculating the rates away for the Kappa parameter

     double SiteNRGCalc(Interaction *FO, double &sNRG, int AAcol,
                int testAA, int *tmpAAseq);

// CHANGE by JEFF     int nIndex, currentAA, siteAA, siteNuc, siteKappa, killAAcol;
     int nIndex, currentAA, siteAA, siteNuc, siteKappa;
     double siteNRG, rate, siteRate, K, W, S, P, U;
// cHANGE by JEFF     double siteNRG, sNRG, rate, siteRate, K, W, S, P, U;
     int *aaRowPtr, *nucRowPtr, *kappaRowPtr;

     if(MULTI_K){K = pE->kappa[whichBranch];}
     else{K = pE->kappa[0];}
     if(MULTI_W){W = pE->omega[whichBranch];}
     else{W = pE->omega[0];}

     U = pE->rate[whichBranch];
     S = pE->solvent;            P = pE->pairwise;

     double *I_nf = pE->MCMCnf;
     int **FO_neighList = FO->neighList;

     switch(round)
     {
        case 0: K = pE->tempParameter; break;
        case 1: U = pE->tempParameter; break;
        case 2: W = pE->tempParameter; break;
        case 3: S = pE->tempParameter; break;
        case 4: P = pE->tempParameter; break;
        case 5: I_nf = pE->tempMCMCnf; break;
     }
     double *pairwisePtr, *expPtr;
     pairwisePtr = solvPtr + 9;
     expPtr  = pairwisePtr + 9;

     siteRate = 0.0;
     nIndex = siteNeigh[AAcol];// Obtain the correct row from Neighbor list
     // run through the 9 nearest neighbors
     aaRowPtr = FO_neighList[nIndex];
     nucRowPtr = aaRowPtr + 9;
     kappaRowPtr = nucRowPtr + 9;
     currentAA = tmpAAseq[AAcol];

     for(int neighbor = 0; neighbor < 9; neighbor++)
     {
        rate = 0.0; //re-initialize the rate variable
        siteAA    = *aaRowPtr;
        siteNuc   = *nucRowPtr;
        siteKappa = *kappaRowPtr;
        if(siteAA == currentAA) // If the change is a synonymous one
        {
           if(siteKappa == 1){rate = I_nf[siteNuc];}
           else{rate = I_nf[siteNuc] * K;}
        }
        else  // it is a nonsynonymous change
        {
           if(siteAA != 20) // if not a stop
           {
              siteNRG = 1.0;
              #if (SIM_W_S || SIM_W_P)
              if((round == 3) || (round == 4))
              {
                 *expPtr = siteNRG = exp((*solvPtr)*S + (*pairwisePtr)*P);
              }
              else{siteNRG = (*expPtr);}
              #endif
              if(siteKappa == 1){rate = I_nf[siteNuc] * W * siteNRG;}
              else{rate = I_nf[siteNuc] * K * W * siteNRG;}
           }
        }
        siteRate += rate;
        //update the pointers to the structure to advance to the next entry
        solvPtr++;  pairwisePtr++; expPtr++;
        aaRowPtr++; nucRowPtr++;   kappaRowPtr++;
     } // end running through all of the nearest neighbors
     totalRate += (U * (siteRate - rateAway[AAcol]));
     rateAway[AAcol] = siteRate;
} // End FastCalcSingleSub Subroutine
//******************************************************************************
void FirstCalcSingleSub(Bayes *pE, Interaction *FO, double *rateAway,int *tmpAAseq,
     int *siteNeigh,int whichBranch, int AAcol, double &totalRate, double *solvPtr)
{
     // This is for calculating the rates away for the Kappa parameter

     double SiteNRGCalc(Interaction *FO, double &sNRG, int AAcol,
                int testAA, int *tmpAAseq);

     int nIndex, currentAA, siteAA, siteNuc, siteKappa;
//CHANGE by JEFF     int nIndex, currentAA, siteAA, siteNuc, siteKappa, killAAcol;
     double siteNRG, sNRG, rate, siteRate, K, W, S, P, U;
     int *aaRowPtr, *nucRowPtr, *kappaRowPtr;

//     if(MULTI_K){K = pE->kappa[whichBranch];}
//     else{K = pE->kappa[0];}
     if(MULTI_W){W = pE->omega[whichBranch];}
     else{W = pE->omega[0];}
     K = pE->tempParameter;
     U = pE->rate[whichBranch];
     S = pE->solvent;            P = pE->pairwise;

     double *I_nf = pE->MCMCnf;
     int **FO_neighList = FO->neighList;

     double *pairwisePtr, *expPtr;
     pairwisePtr = solvPtr + 9;
     expPtr  = pairwisePtr + 9;

     siteRate = 0.0;
     nIndex = siteNeigh[AAcol];// Obtain the correct row from Neighbor list
     // run through the 9 nearest neighbors
     aaRowPtr = FO_neighList[nIndex];
     nucRowPtr = aaRowPtr + 9;
     kappaRowPtr = nucRowPtr + 9;
     currentAA = tmpAAseq[AAcol];

     for(int neighbor = 0; neighbor < 9; neighbor++)
     {
        rate = 0.0; //re-initialize the rate variable
        siteAA    = *aaRowPtr;
        siteNuc   = *nucRowPtr;
        siteKappa = *kappaRowPtr;
        if(siteAA == currentAA) // If the change is a synonymous one
        {
           if(siteKappa == 1){rate = I_nf[siteNuc];}
           else{rate = I_nf[siteNuc] * K;}
        }
        else  // it is a nonsynonymous change
        {
           if(siteAA != 20) // if not a stop
           {
              siteNRG = 1.0;
              #if (SIM_W_S || SIM_W_P)
              *pairwisePtr = siteNRG = SiteNRGCalc(FO,sNRG,AAcol,siteAA,tmpAAseq);
              *solvPtr = sNRG;
              *expPtr  = siteNRG = exp(sNRG*S + siteNRG*P);
              #endif
              if(siteKappa == 1){rate = I_nf[siteNuc] * W * siteNRG;}
              else{rate = I_nf[siteNuc] * K * W * siteNRG;}
           }
        }
        siteRate += rate;
        //update the pointers to the structure to advance to the next entry
        solvPtr++;  pairwisePtr++; expPtr++;
        aaRowPtr++; nucRowPtr++;   kappaRowPtr++;
     } // end running through all of the nearest neighbors
     totalRate += (U * (siteRate - rateAway[AAcol]));
     rateAway[AAcol] = siteRate;
} // End FirstCalcSingleSub Subroutine
//******************************************************************************
void CalcSingleSub(Bayes *pE, Interaction *FO, double *rateAway, int *tmpAAseq,
       int *siteNeigh, int whichBranch, int AAcol, double &totalRate)
{
     double SiteNRGCalc(Interaction *FO, double &sNRG, int AAcol,
                int testAA, int *tmpAAseq);

// CHANGE by JEFF     int nIndex, currentAA, siteAA, siteNuc, siteKappa, killAAcol;
     int nIndex, currentAA, siteAA, siteNuc, siteKappa;
     double siteNRG, sNRG, rate, siteRate, K, W, S, P, U;
     int *aaRowPtr, *nucRowPtr, *kappaRowPtr;

     if(MULTI_K){K = pE->kappa[whichBranch];}
     else{K = pE->kappa[0];}
     if(MULTI_W){W = pE->omega[whichBranch];}
     else{W = pE->omega[0];}

     U = pE->rate[whichBranch];
     S = pE->solvent;            P = pE->pairwise;

     double *I_nf = pE->MCMCnf;
     int **FO_neighList = FO->neighList;

     siteRate = 0.0;
     nIndex = siteNeigh[AAcol];// Obtain the correct row from Neighbor list
     // run through the 9 nearest neighbors
     aaRowPtr = FO_neighList[nIndex];
     nucRowPtr = aaRowPtr + 9;
     kappaRowPtr = nucRowPtr + 9;
     currentAA = tmpAAseq[AAcol];

     for(int neighbor = 0; neighbor < 9; neighbor++)
     {
        rate = 0.0; //re-initialize the rate variable
        siteAA    = *aaRowPtr;
        siteNuc   = *nucRowPtr;
        siteKappa = *kappaRowPtr;
        if(siteAA == currentAA) // If the change is a synonymous one
        {
           if(siteKappa == 1){rate = I_nf[siteNuc];}
           else{rate = I_nf[siteNuc] * K;}
        }
        else  // it is a nonsynonymous change
        {
           if(siteAA != 20) // if not a stop
           {
              siteNRG = 1.0;
              #if (SIM_W_S || SIM_W_P)
              siteNRG = SiteNRGCalc(FO,sNRG,AAcol,siteAA,tmpAAseq);
              siteNRG = exp(sNRG*S + siteNRG*P);
              #endif
              if(siteKappa == 1){rate = (I_nf[siteNuc] * W * siteNRG);}
              else{rate = (I_nf[siteNuc] * K * W * siteNRG);}
           }
        }
        siteRate += rate;
        aaRowPtr++; nucRowPtr++; kappaRowPtr++;
     } // end running through all of the nearest neighbors
     totalRate += (U * (siteRate - rateAway[AAcol]));
     rateAway[AAcol] = siteRate;
} // End CalcSingleSub Subroutine
//******************************************************************************
double SiteNRGCalc(Interaction *FO, double &sNRG, int AAcol, int testAA,
                int *tmpAAseq)
{
   AAsiteInfo &AAinfoPtr = FO->AAinfo[AAcol];
   int siteAA, originalAA;
// CHANGE by JEFF   int sAcc, siteAA, originalAA;
   double NRG_diff; NRG_diff = sNRG = 0.0;
   originalAA = tmpAAseq[AAcol];

//#if SIM_W_S
      sNRG = AAinfoPtr.solvent[originalAA] - AAinfoPtr.solvent[testAA];
//#endif

//#if SIM_W_P  // Calculate the NRG of the affected site
      double **T2NRG, *T1NRG;
      int  *DDD = AAinfoPtr.DDDneighbors;
      int FO_numneigh = AAinfoPtr.numNeighbors;
      T2NRG = AAinfoPtr.pairwise;

      for(int column = 0; column < FO_numneigh; column++){
         T1NRG = *T2NRG;
         siteAA = tmpAAseq[*DDD];
         // Calculate the energy of the actual change that occurred
         // Calculate the difference this change makes on all other sites
         // This is calculated by taking (Before - After(proposed))
         if(*DDD < AAcol){
             NRG_diff +=(T1NRG[siteAA*20+originalAA] - T1NRG[siteAA*20+testAA]);
         }
         else{
             NRG_diff +=(T1NRG[originalAA*20+siteAA] - T1NRG[testAA*20+siteAA]);
         }
         DDD++; T2NRG++;
      } // End looping through the 3D nearest neighbors
//#endif // end calculating the NRG of the affected site
   return NRG_diff;
/*******************************************************************
// Now something must be thought of to handle the energies.  Now  //
// that I have them, the thought is to take the negative          //
// exponential of their difference.  Instead, we will reverse the //
// order and take the positive exponential instead.  The other    //
// is some variant of this thought, but I will try that later.    //
// exp (w *G(S1 - S2)) IDEA 1: will implement this one now        //
// or    Let ES* = median NRG of initial sequence one and two     //
// exp (w |ES* - ES1| - |ES* - ES2|)  IDEA 2: will try this later //
*******************************************************************/
}  // End SiteNRGCalc Subroutine
//******************************************************************************
void GeneralRateAway(Bayes *pE,Interaction *FO,double *rateAway,int *tmpAAseq,
     int whichBranch, int *siteNeigh, double &totalRate, int AAcol, int synFlag)
{
    void CalcSingleSub(Bayes *pE,Interaction *FO,double *rateAway,int *tmpAAseq,
      int *siteNeigh, int whichBranch, int AAcol, double &totalRate);

    // NOTE: At this point in the process, the change has already
    //  been introduced into the nucleotide sequence, the
    //   amino acid sequence, and the neighAcc holders.

    // synFlag = 0: nonsynonymous
    // synFlag = 1: synonymouns

    if(synFlag) // SYNONYMOUS SUBSTITUTION EVENT
    {
       CalcSingleSub(pE, FO, rateAway, tmpAAseq, siteNeigh, whichBranch, AAcol,
                        totalRate);
    }
    else // NONSYNONYMOUS SUBSTITUTION EVENT!!!
    {
       CalcSingleSub(pE, FO, rateAway, tmpAAseq, siteNeigh, whichBranch, AAcol,
                        totalRate);

       AAsiteInfo &AAinfoPtr = FO->AAinfo[AAcol];

       int *DDD = AAinfoPtr.DDDneighbors;
       int FO_numberNeigh = AAinfoPtr.numNeighbors;

       for(int AAcolumn = 0; AAcolumn < FO_numberNeigh; AAcolumn++)
       {
         CalcSingleSub(pE, FO, rateAway, tmpAAseq, siteNeigh, whichBranch,
                                *DDD, totalRate);
         DDD++;
       }
    }
}  // End GeneralRateAway Subroutine
//******************************************************************************
void FastGeneralRateAway(Bayes *pE,Interaction *FO,double *rateAway,int *tmpAAseq,
     int whichBranch, int *siteNeigh, double &totalRate, int AAcol,
     double **SpeedRate,int synFlag, int round)
{
     void FirstCalcSingleSub(Bayes *pE, Interaction *FO, double *rateAway,
        int *tmpAAseq,int *siteNeigh,int whichBranch, int AAcol,
        double &totalRate, double *AllRates);
     void FastCalcSingleSub(Bayes *pE, Interaction *FO, double *rateAway,
        int *tmpAAseq, int *siteNeigh,int whichBranch, int AAcol,
        double &totalRate, double *AllRates, int round);

    // NOTE: At this point in the process, the change has already
    //  been introduced into the nucleotide sequence, the
    //   amino acid sequence, and the neighAcc holders.

    // synFlag = 0: nonsynonymous
    // synFlag = 1: synonymouns

    if(synFlag) // SYNONYMOUS SUBSTITUTION EVENT
    {
       if(!round)
       {
          FirstCalcSingleSub(pE,FO,rateAway,tmpAAseq,siteNeigh,whichBranch,AAcol,
                                totalRate, SpeedRate[AAcol]);
       }
       else
       {
	  FastCalcSingleSub(pE,FO,rateAway,tmpAAseq, siteNeigh,whichBranch,AAcol,
				totalRate, SpeedRate[AAcol], round);
       }
    }
    else // NONSYNONYMOUS SUBSTITUTION EVENT!!!
    {
       if(!round)
       {
          FirstCalcSingleSub(pE,FO,rateAway,tmpAAseq,siteNeigh,whichBranch,AAcol,
                                totalRate, SpeedRate[AAcol]);

          AAsiteInfo &AAinfoPtr = FO->AAinfo[AAcol];

          int *DDD = AAinfoPtr.DDDneighbors;
          int FO_numberNeigh = AAinfoPtr.numNeighbors;

          for(int AAcolumn = 0; AAcolumn < FO_numberNeigh; AAcolumn++)
          {
            FirstCalcSingleSub(pE,FO,rateAway,tmpAAseq,siteNeigh,whichBranch,
                            *DDD,totalRate, SpeedRate[*DDD]);
            DDD++;
          }
       }
       else
       {
	  FastCalcSingleSub(pE,FO,rateAway,tmpAAseq, siteNeigh,whichBranch,AAcol,
				totalRate, SpeedRate[AAcol], round);

          AAsiteInfo &AAinfoPtr = FO->AAinfo[AAcol];

          int *DDD = AAinfoPtr.DDDneighbors;
          int FO_numberNeigh = AAinfoPtr.numNeighbors;

          for(int AAcolumn = 0; AAcolumn < FO_numberNeigh; AAcolumn++)
          {
            FastCalcSingleSub(pE,FO,rateAway,tmpAAseq,siteNeigh,whichBranch,
                            *DDD,totalRate, SpeedRate[*DDD], round);
            DDD++;
          }
       }
    }
}  // End FastGeneralRateAway Subroutine
//******************************************************************************
double ProbPathWTimes(Bayes *pE, Interaction *FO,Information &Info,
     Substitution *start, int *tmpAAseq,int *tmpNUCseq,int *siteNeigh,double gH,
     double sH, int count, int which, int LASTFLAG, double maxTime)
{
   double SetRateAway(Bayes *pE, Interaction *FO, double *rateAway,int *tmpAAseq,
		int whichBranch, int *siteNeigh, int AAseqLen);
   void GeneralRateAway(Bayes *pE,Interaction *FO,double *rateAway,int *tmpAAseq,
     int whichBranch, int *siteNeigh, double &totalRate, int AAcol, int synFlag);

   double SlowSpecificRateAway(Bayes *pE, Interaction * FO, Substitution *kill,
                int whichBranch, int* tmpAAseq);
   int NeighborAccess(int AA, int *c, int **nAccess);

   double W, specificRate, time, W_time;
   int c[3], synFlag, AAseqLen; AAseqLen = Info.AAlen;
   Substitution *next; next = NULL;
   Substitution *kill; kill = Info.startPath;
   if(kill != NULL)
   {
      if(start->time != kill->time){
         kill=start->DpathPtr;
         if(kill == NULL){kill = start;}
      }
   }
   static int z = 0;
   static double *rateAway;
   if(z == 0){rateAway = new double[AAseqLen]; z = 1;}

   W = SetRateAway(pE,FO,rateAway,tmpAAseq,which,siteNeigh,AAseqLen);

   if(LASTFLAG)
   {
      while(kill != NULL){
         Substitution *previous = kill->UpathPtr;
         if(previous == NULL){time = kill->time;}  // T_1 = t_1
         else{time = kill->time - previous->time;} // t_i = t_i - t_{i-1}
         sH *= kill->subHood = specificRate = SlowSpecificRateAway(pE,FO,kill,which,tmpAAseq);
         gH += W_time = -W*time;
         kill->genHood = W_time;
         if(sH > BIGNUM){sH *= INVBIGNUM; count++;}
         else if(sH < INVBIGNUM){sH *= BIGNUM; count--;}
         int NUCcol = kill->column;
         int AAcol = NUCcol / 3;
         tmpNUCseq[NUCcol] = kill->pathNuc; // updates the nuc sequence
         tmpAAseq[AAcol] = kill->amino; // updates the AA sequence
         switch(NUCcol % 3){ // We need the correct full codon for this site
            case 0: c[0] = tmpNUCseq[NUCcol]; c[1] = tmpNUCseq[NUCcol+1];
                    c[2] = tmpNUCseq[NUCcol+2]; break;
            case 1: c[0] = tmpNUCseq[NUCcol-1]; c[1] = tmpNUCseq[NUCcol];
                    c[2] = tmpNUCseq[NUCcol+1]; break;
            case 2: c[0] = tmpNUCseq[NUCcol-2]; c[1] = tmpNUCseq[NUCcol-1];
                    c[2] = tmpNUCseq[NUCcol];   break;
         } // siteNeigh is truly an intense speed up in the code!
         siteNeigh[AAcol]=NeighborAccess(kill->amino,c,FO->nAccess);
         if(kill->priorAA == kill->amino){synFlag = 1;}//SYNONYMOUS
         else{synFlag = 0;}//NONSYNONYMOUS
         GeneralRateAway(pE,FO,rateAway,tmpAAseq,which,siteNeigh,W,AAcol,synFlag);
         if(kill->DpathPtr == NULL){next = kill; kill = kill->DpathPtr;}
         else{kill = kill->DpathPtr;}
      }
      if(next != NULL){kill = next;}
      if(Info.totNumSub != 0){time = 1.0 - kill->time;}  // T - T_n
      else{time = 1.0;}
      gH += Info.probLastEvent = -W*time; // ln ( e^(-WT))
      if(sH < 0.0){cerr << "WARNING: sH < 0.0 " << endl;}
      Info.pathHood = gH+log(sH)+count*LOGBIGNUM;
   }
   else
   {
       while(kill->time <= maxTime){
         Substitution *previous = kill->UpathPtr;
         if(previous == NULL){time = kill->time;}  // T_1 = t_1
         else{time = kill->time - previous->time;} // t_i = t_i - t_{i-1}
         sH *= kill->subHood = specificRate = SlowSpecificRateAway(pE,FO,kill,which,tmpAAseq);
         gH += W_time = -W*time;
         kill->genHood = W_time;
         if(sH > BIGNUM){sH *= INVBIGNUM; count++;}
         else if(sH < INVBIGNUM){sH *= BIGNUM; count--;}
         int NUCcol = kill->column;
         int AAcol = NUCcol / 3;
         tmpNUCseq[NUCcol] = kill->pathNuc; // updates the nuc sequence
         tmpAAseq[AAcol] = kill->amino; // updates the AA sequence
         switch(NUCcol % 3){ // We need the correct full codon for this site
            case 0: c[0] = tmpNUCseq[NUCcol]; c[1] = tmpNUCseq[NUCcol+1];
                    c[2] = tmpNUCseq[NUCcol+2]; break;
            case 1: c[0] = tmpNUCseq[NUCcol-1]; c[1] = tmpNUCseq[NUCcol];
                    c[2] = tmpNUCseq[NUCcol+1]; break;
            case 2: c[0] = tmpNUCseq[NUCcol-2]; c[1] = tmpNUCseq[NUCcol-1];
                    c[2] = tmpNUCseq[NUCcol];   break;
         } // siteNeigh is truly an intense speed up in the code!
         siteNeigh[AAcol]=NeighborAccess(kill->amino,c,FO->nAccess);
         if(kill->priorAA == kill->amino){synFlag = 1;}//SYNONYMOUS
         else{synFlag = 0;}//NONSYNONYMOUS
         GeneralRateAway(pE,FO,rateAway,tmpAAseq,which,siteNeigh,W,AAcol,synFlag);
         if(kill->DpathPtr == NULL){next = kill; kill = kill->DpathPtr;}
         else{kill = kill->DpathPtr;}
       }
      while(kill != NULL)
      {
         gH += kill->genHood;
         sH *= kill->subHood;
         if(sH > BIGNUM){sH *= INVBIGNUM; count++;}
         else if(sH < INVBIGNUM){sH *= BIGNUM; count--;}
         if(sH < 0.0){cerr << "WARNING: sH < 0.0  in probPathWTimes update" << endl;}
         if(kill->DpathPtr == NULL){next = kill; kill = kill->DpathPtr;}
         else{kill = kill->DpathPtr;}
      }
      gH += Info.probLastEvent;
      Info.pathHood = gH+log(sH)+count*LOGBIGNUM;
   }
   return Info.pathHood;
}  // End ProbPathWTimes Subroutine
//******************************************************************************
double FastProbParam(Bayes *pE, Interaction *FO, Information &Info, int round,
                int which, ofstream &o, double *genHood, double *subHood,
                double ***AllRates)
{    // Calculates the probability of the independent site model without
        // integrating over all times
   double FirstSetRateAway(Bayes *pE, Interaction *FO, double *rateAway,
         int *tmpAAseq,int whichBranch, int *siteNeigh, int AAseqLen,
         double **AllRates);

   double FastSetRateAway(Bayes *pE, Interaction *FO, double *rateAway,
         int *tmpAAseq, int whichBranch, int *siteNeigh, int AAseqLen,
         double **AllRates, int round);

   void FastGeneralRateAway(Bayes *pE,Interaction *FO,double *rateAway,
         int *tmpAAseq,int whichBranch, int *siteNeigh, double &totalRate,
         int AAcol, double **SpeedRate,int synFlag, int round);

   double FastSpecificRateAway(Bayes *pE,Substitution *kill,int which,int round);

   int NeighborAccess(int AA, int *c, int **nAccess);

   double W, tmpGenHood, tmpSiteHood, time, specificRate, *gH, *sH, W_time;
   Info.avgRate = tmpGenHood = 0.0; tmpSiteHood = 1.0;
   gH = genHood;   sH = subHood;
   Substitution *next; // next = NULL;
   Substitution *kill; kill = Info.startPath;
   int c[3], InflateCount, InfoAAlen, Infolen, synFlag; InflateCount = 0;
   //Speed up access to the Info-object...
   InfoAAlen = Info.AAlen;  Infolen = Info.len;

   static int z = 0;
   static int *tmpNUCseq, *tmpAAseq, *siteNeigh;
   static double *rateAway;

   if(z == 0){
      tmpNUCseq = new int[Infolen];
      tmpAAseq  = new int[InfoAAlen];
      siteNeigh = new int[InfoAAlen];
      rateAway  = new double[InfoAAlen];
      z = 1;
   }
   int *t1AA, *t2AA, *t1NUC, *t2NUC, *t1NEIGH, *t2NEIGH;
   t1AA = tmpAAseq; t2AA = Info.parentAAseq;
   t1NUC = tmpNUCseq; t2NUC = Info.parentSeq;
   t1NEIGH = siteNeigh; t2NEIGH = Info.neighAcc;
   for(int col = 0; col < Infolen; col++)
   {// To protect the sequence, we copy it into temporary and work on that
      if(col < InfoAAlen)
      {
         *t1AA = *t2AA; *t1NEIGH = *t2NEIGH;
         t1AA++; t2AA++; t1NEIGH++; t2NEIGH++;
      }
      *t1NUC = *t2NUC; t1NUC++; t2NUC++;
   }

   if(!round)
   {
//      o << "about to do FirstSet RateAway " << endl;
      W = FirstSetRateAway(pE,FO,rateAway,tmpAAseq,which,siteNeigh,InfoAAlen,
		AllRates[0]);
   }
   else
   {
      W = FastSetRateAway(pE,FO,rateAway,tmpAAseq,which, siteNeigh, InfoAAlen,
		AllRates[0], round);
   }
//   W = SetRateAway(pE,FO,rateAway,tmpAAseq,which,siteNeigh,InfoAAlen,round);
   int I_totsub = Info.totNumSub;

   for(int sub = 0; sub < I_totsub; sub++)
   {
      if(sub == 0) {time = kill->time;}  // T_1 = t_1
      else{
         next = kill->UpathPtr;
         time = kill->time - next->time; // t_i = t_i - t_{i-1}
      }
      tmpSiteHood *= *sH = specificRate = FastSpecificRateAway(pE,kill,which,round);
      if(specificRate <= 0)
      {
         cerr << "In FastProbParam, the specific rate < 0 " << endl;
         o << "In FastProbParam, the specific rate < 0 " << endl;
         o << "param  " << round << endl;
      }
      Info.avgRate += W_time = W*time;
      tmpGenHood += *gH = -W_time;
      if(tmpSiteHood > BIGNUM){
         tmpSiteHood *= INVBIGNUM;
         InflateCount++;
      }
      else if(tmpSiteHood < INVBIGNUM){
         tmpSiteHood *= BIGNUM;
         InflateCount--;
      }

      int NUCcol = kill->column;
      int AAcol = NUCcol / 3;
      tmpNUCseq[NUCcol] = kill->pathNuc; // updates the nuc sequence
      tmpAAseq[AAcol] = kill->amino; // updates the AA sequence
      switch(NUCcol % 3){
      case 0: c[0] = tmpNUCseq[NUCcol]; c[1] = tmpNUCseq[NUCcol+1];
              c[2] = tmpNUCseq[NUCcol+2]; break;
      case 1: c[0] = tmpNUCseq[NUCcol-1]; c[1] = tmpNUCseq[NUCcol];
              c[2] = tmpNUCseq[NUCcol+1]; break;
      case 2: c[0] = tmpNUCseq[NUCcol-2]; c[1] = tmpNUCseq[NUCcol-1];
              c[2] = tmpNUCseq[NUCcol];   break;
      }
      siteNeigh[AAcol]=NeighborAccess(kill->amino,c,FO->nAccess);
      if(kill->priorAA == kill->amino){synFlag = 1;}//SYNONYMOUS
      else{synFlag = 0;}//NONSYNONYMOUS
      FastGeneralRateAway(pE,FO,rateAway,tmpAAseq,which, siteNeigh, W, AAcol,
                 AllRates[sub+1], synFlag, round);
      if(sub < (I_totsub - 1)){kill = kill->DpathPtr;}
      gH++; sH++;
   }

   if(I_totsub != 0){time = 1.0 - kill->time;}  // T - T_n
   else{time = 1.0;}
   tmpGenHood += Info.probLastEvent = W_time = -W*time;// ln ( e^(-WT))
   Info.avgRate -= W_time;
   Info.pathHood = tmpGenHood+log(tmpSiteHood)+InflateCount*LOGBIGNUM;
   return Info.pathHood;
}  // End FastProbParam Subroutine
//******************************************************************************
double InitProbParam(Bayes *pE, Interaction *FO, Information &Info, int which,
                ofstream &o)
{    // Calculates the probability of the independent site model without
        // integrating over all times
   double SetRateAway(Bayes *pE, Interaction *FO, double *rateAway,int *tmpAAseq,
		int whichBranch, int *siteNeigh, int AAseqLen);
   double SlowSpecificRateAway(Bayes *pE, Interaction * FO,
             Substitution *kill, int whichBranch, int* tmpAAseq);
   void GeneralRateAway(Bayes *pE,Interaction *FO,double *rateAway,int *tmpAAseq,
     int whichBranch, int *siteNeigh, double &totalRate, int AAcol, int synFlag);
   int NeighborAccess(int AA, int *c, int **nAccess);

   double W, tmpGenHood, tmpSiteHood, time, specificRate, W_time;
   W = tmpGenHood = Info.avgRate = 0.0; tmpSiteHood = 1.0;
   int c[3];
   int InflateCount = 0;
   int nucLen, AAseqLen, totsub, synFlag;
   nucLen = Info.len;  AAseqLen = Info.AAlen;  totsub = Info.totNumSub;


   int *tmpNUCseq     = new int[nucLen];
   int *tmpAAseq      = new int[AAseqLen];
   int *siteNeigh     = new int[AAseqLen];
   double *rateAway   = new double[AAseqLen];
   int **FO_nAccess = FO->nAccess;
   int *t1AA, *t2AA, *t1NUC, *t2NUC, *t1NEIGH, *t2NEIGH;
   t1AA = tmpAAseq; t2AA = Info.parentAAseq;
   t1NUC = tmpNUCseq; t2NUC = Info.parentSeq;
   t1NEIGH = siteNeigh; t2NEIGH = Info.neighAcc;
   for(int col = 0; col < nucLen; col++)
   {// To protect the sequence, we copy it into temporary and work on that
      if(col < AAseqLen)
      {
         *t1AA = *t2AA; *t1NEIGH = *t2NEIGH;
         t1AA++; t2AA++; t1NEIGH++; t2NEIGH++;
      }
      *t1NUC = *t2NUC; t1NUC++; t2NUC++;
   }
   Substitution *next;
   Substitution *kill = Info.startPath;

   W = SetRateAway(pE,FO,rateAway,tmpAAseq,which,siteNeigh,AAseqLen);
   for(int sub = 0; sub < totsub; sub++)
   {
      if(sub == 0) {time = kill->time;}  // T_1 = t_1
      else{
         next = kill->UpathPtr;
         time = kill->time - next->time; // t_i = t_i - t_{i-1}
      }
      tmpSiteHood *= kill->subHood=specificRate=SlowSpecificRateAway(pE,FO,kill,which,tmpAAseq);
      tmpGenHood += kill->genHood = W_time = -W*time;
      Info.avgRate += -W_time;
      if(tmpSiteHood > BIGNUM){
         tmpSiteHood *= INVBIGNUM;
         InflateCount++;
      }
      else if(tmpSiteHood < INVBIGNUM){
         tmpSiteHood *= BIGNUM;
         InflateCount--;
      }
      int NUCcol = kill->column;
      int AAcol = NUCcol / 3;
      tmpNUCseq[NUCcol] = kill->pathNuc; // updates the nuc sequence
      tmpAAseq[AAcol] = kill->amino; // updates the AA sequence
      switch(NUCcol % 3){
      case 0: c[0] = tmpNUCseq[NUCcol]; c[1] = tmpNUCseq[NUCcol+1];
              c[2] = tmpNUCseq[NUCcol+2]; break;
      case 1: c[0] = tmpNUCseq[NUCcol-1]; c[1] = tmpNUCseq[NUCcol];
              c[2] = tmpNUCseq[NUCcol+1]; break;
      case 2: c[0] = tmpNUCseq[NUCcol-2]; c[1] = tmpNUCseq[NUCcol-1];
              c[2] = tmpNUCseq[NUCcol];   break;
      }
      siteNeigh[AAcol]=NeighborAccess(kill->amino,c,FO_nAccess);
      if(kill->priorAA == kill->amino){synFlag = 1;}//SYNONYMOUS
      else{synFlag = 0;}//NONSYNONYMOUS
      GeneralRateAway(pE,FO,rateAway,tmpAAseq,which,siteNeigh,W,AAcol,synFlag);
      if(sub < (totsub - 1)){kill = kill->DpathPtr;}
   }

   //time = 1.0 - kill->time;  Had to update to allow 0 substitutions
   // Doug Changes August 2005

   if(totsub != 0){time = 1.0 - kill->time;}  // T - T_n
   else{time = 1.0;}

   tmpGenHood += Info.probLastEvent = W_time = -W*time;// ln ( e^(-WT))
   Info.avgRate -= W_time;
   Info.pathHood = tmpGenHood+log(tmpSiteHood)+InflateCount*LOGBIGNUM;
//   o << "This is the likelihood " << Info.pathHood << endl;
//   o << "tmpgen = " << tmpGenHood << "  log site = " << log(tmpSiteHood);
//   o << "  count  = " <<  InflateCount*LOGBIGNUM << " TOTAL = " <<    Info.pathHood << endl;
//    o << endl << endl;
   delete[] tmpNUCseq;   delete[] tmpAAseq; // return the memory we borrowed
   delete[] siteNeigh; delete[] rateAway;
   return Info.pathHood;
}  // End InitProbParam Subroutine
//******************************************************************************
/*double ProfileNRG(Interaction *FO, double &sNRG, int AAcol, int testAA,
                int *tmpAAseq)
{
   int *access;
   double ***NRG0, ***NRG1, ***NRG2, ***NRG3, ***NRG4;
   double **NRGmat0, **NRGmat1, **NRGmat2, **NRGmat3, **NRGmat4;
   int ***FO_energyAcc   = FO->energyAcc;
   int *DDD = FO->DDDneighbors[AAcol];
   double **FO_solvent = FO->solvent;

   int sAcc, siteAA, originalAA; originalAA = tmpAAseq[AAcol];
   double NRG_diff = 0.0;

   sAcc = FO->solvAcc[AAcol];
   sNRG = FO_solvent[originalAA][sAcc] - FO_solvent[testAA][sAcc];


   NRG0 = FO->energy[0];  NRG1 = FO->energy[1];  NRG2 = FO->energy[2];
   NRG3 = FO->energy[3];  NRG4 = FO->energy[4];
   int FO_numneigh = FO->numNeigh[AAcol];

   // Calculate the NRG of the affected site
   for(int column = 0; column < FO_numneigh; column++)
   {
      access = FO_energyAcc[AAcol][column]; // obtain correct matrix
      siteAA = tmpAAseq[*DDD];
      // Calculate the energy of the actual change that occurred
      // Calculate the difference this change makes on all other sites
      // This is calculated by taking (Before - After(proposed))
      NRGmat0 = NRG0[access[0]];
      NRGmat1 = NRG1[access[1]];
      NRGmat2 = NRG2[access[2]];
      NRGmat3 = NRG3[access[3]];
      NRGmat4 = NRG4[access[4]];

      if(*DDD < AAcol){
         NRG_diff +=
                NRGmat0[siteAA][originalAA] + NRGmat1[siteAA][originalAA] +
                NRGmat2[siteAA][originalAA] + NRGmat3[siteAA][originalAA] +
                NRGmat4[siteAA][originalAA] - NRGmat0[siteAA][testAA] -
                NRGmat1[siteAA][testAA] - NRGmat2[siteAA][testAA] -
                NRGmat3[siteAA][testAA] - NRGmat4[siteAA][testAA];
       }
      else{
         NRG_diff +=
                NRGmat0[originalAA][siteAA] + NRGmat1[originalAA][siteAA] +
                NRGmat2[originalAA][siteAA] + NRGmat3[originalAA][siteAA] +
                NRGmat4[originalAA][siteAA] - NRGmat0[testAA][siteAA] -
                NRGmat1[testAA][siteAA] - NRGmat2[testAA][siteAA] -
                NRGmat3[testAA][siteAA] - NRGmat4[testAA][siteAA];
      }
      DDD++;
   }

   NRG_diff = (sNRG + NRG_diff);
   return -1.0*NRG_diff;
}  // End ProfileNRG Subroutine
//******************************************************************************
void PathProfile(Information *C,Interaction *FirstOrder,double **profile)
{
     double ProfileNRG(Interaction *FirstOrder, double &sNRG,int AAcol,
                int testAA, int *tmpAAseq);
     static int z = 0;
     int col, count; count = 0;
     int *tmpAAseq = new int[C->AAlen];
     for(int i = 0; i < C[0].AAlen; i++){tmpAAseq[i] = C[0].AA_seq[i];}
     double seqNRG,sNRG,siteNRG, timeStep, tempTime; seqNRG = C[0].AAseqNRG;
     timeStep = 1.0 / 20.0;
     tempTime = timeStep;
     Substitution *temp; temp = C[0].startPath;
     profile[z][count] = C[0].AAseqNRG;

     while(temp != NULL)
     {
        if(temp->time > tempTime)
        {
           count++;
           profile[z][count] = seqNRG;
           tempTime+= timeStep;
        }
        if(temp->amino != temp->priorAA)
        {
           col = temp->column/3;
           siteNRG = ProfileNRG(FirstOrder,sNRG,col,temp->amino,tmpAAseq);
           seqNRG += siteNRG;
           tmpAAseq[col] = temp->amino;
        }
        temp = temp->DpathPtr;
     } // end while loop
     count++;
     profile[z][count] = seqNRG;
//     o << seqNRG << endl;
     delete[] tmpAAseq;
     z++;
}  // End PathProfile Subroutine
*/
//******************************************************************************
void PathPosterior(Bayes *pE, Substitution *tmp, Interaction *FO,int count,
         int q, double *nE, double *mT, int *newOrder, int *picSeq, ofstream &o)
{
   double SiteNRGCalc(Interaction *FO, double &sNRG, int AAcol, int testAA,
                int *tmpAAseq);
   double sNRG, pNRG, eNRG;

   if(MULTI_W)
   {
     while(tmp!=NULL)
     {
       nE[tmp->column]+=1.0;
       mT[tmp->column]+=tmp->time;
//       if(tmp->amino != tmp->priorAA)
//       {

         int tmpCol = tmp->column/3;
         pNRG = SiteNRGCalc(FO, sNRG, tmpCol, tmp->amino,picSeq);
         eNRG = exp(sNRG*pE->solvent + pNRG*pE->pairwise);
         o << count << " " << (tmp->column-2)+(3*newOrder[tmpCol]) << " " <<
         tmpCol + newOrder[tmpCol] << " " << tmp->time << " " << tmp->priorNuc <<
         " " << tmp->pathNuc << " " << tmp->priorAA << " " << tmp->amino <<" "<<
         pE->omega[q] << " " << pE->solvent << " " << pE->pairwise << " " <<
         sNRG << " " << pNRG << " " << pE->solvent*sNRG << " " <<
         pNRG*pE->pairwise << " " << eNRG << " " << eNRG*pE->omega[q] << endl;
         picSeq[tmpCol] = tmp->amino;
//       }
       tmp = tmp->DpathPtr;
     }
   }
   else
   {
     while(tmp!=NULL)
     {
       nE[tmp->column]+=1.0;
       mT[tmp->column]+=tmp->time;
//       if(tmp->amino != tmp->priorAA)
//       {
         int tmpCol = tmp->column/3;
         pNRG = SiteNRGCalc(FO, sNRG, tmpCol, tmp->amino,picSeq);
         eNRG = exp(sNRG*pE->solvent + pNRG*pE->pairwise);
         o << count << " " << (tmp->column-2)+(3*newOrder[tmpCol]) << " " <<
         tmpCol + newOrder[tmpCol] << " " << tmp->time << " " << tmp->priorNuc <<
         " " << tmp->pathNuc << " " << tmp->priorAA << " " << tmp->amino <<" "<<
         pE->omega[0] << " " << pE->solvent << " " << pE->pairwise << " " <<
         sNRG << " " << pNRG << " " << pE->solvent*sNRG << " " <<
         pNRG*pE->pairwise << " " << eNRG << " " << eNRG*pE->omega[0] << endl;
         picSeq[tmpCol] = tmp->amino;
 //      }
       tmp = tmp->DpathPtr;
     }
   }
}
//******************************************************************************
//******************************************************************************
void PathPostFull(Bayes *pE, Substitution *tmp, Interaction *FO,int count,
         int q, double *nE, double *mT, int *newOrder, double **pathPtr)
{
   double **PathPtr1, *pathRowPtr;
   PathPtr1 = pathPtr;

   if(MULTI_W)
   {
     while(tmp!=NULL)
     {
         pathRowPtr = *PathPtr1;
         nE[tmp->column]+=1.0;
         mT[tmp->column]+=tmp->time;
//       if(tmp->amino != tmp->priorAA)
//       {

         int tmpCol = tmp->column/3;
         pathRowPtr[0] = q;
         pathRowPtr[1] = count;
         pathRowPtr[2] = (tmp->column-2)+(3*newOrder[tmpCol]);
         pathRowPtr[3] = tmpCol + newOrder[tmpCol];
         pathRowPtr[4] = tmp->time;
         pathRowPtr[5] = tmp->priorNuc;
         pathRowPtr[6] = tmp->pathNuc;
         pathRowPtr[7] = tmp->priorAA;
         pathRowPtr[8] = tmp->amino;
         pathRowPtr[9] = pE->omega[q];
         pathRowPtr[10] = pE->solvent;
         pathRowPtr[11] = pE->pairwise;
         pathRowPtr[12] = tmp->solvNRGdiff;
         pathRowPtr[13] = tmp->nrgDiff;
         pathRowPtr[14] = pE->solvent*tmp->solvNRGdiff;
         pathRowPtr[15] = tmp->nrgDiff*pE->pairwise;
         pathRowPtr[16] = tmp->expNRGdiff;
         pathRowPtr[17] = tmp->expNRGdiff*pE->omega[q];


//       }
       tmp = tmp->DpathPtr;  PathPtr1++;
     }
   }
   else
   {
     while(tmp!=NULL)
     {
         pathRowPtr = *PathPtr1;
         nE[tmp->column]+=1.0;
         mT[tmp->column]+=tmp->time;
//       if(tmp->amino != tmp->priorAA)
//       {
         int tmpCol = tmp->column/3;
         pathRowPtr[0] = q;
         pathRowPtr[1] = count;
         pathRowPtr[2] = (tmp->column-2)+(3*newOrder[tmpCol]);
         pathRowPtr[3] = tmpCol + newOrder[tmpCol];
         pathRowPtr[4] = tmp->time;
         pathRowPtr[5] = tmp->priorNuc;
         pathRowPtr[6] = tmp->pathNuc;
         pathRowPtr[7] = tmp->priorAA;
         pathRowPtr[8] = tmp->amino;
         pathRowPtr[9] = pE->omega[0];
         pathRowPtr[10] = pE->solvent;
         pathRowPtr[11] = pE->pairwise;
         pathRowPtr[12] = tmp->solvNRGdiff;
         pathRowPtr[13] = tmp->nrgDiff;
         pathRowPtr[14] = pE->solvent*tmp->solvNRGdiff;
         pathRowPtr[15] = tmp->nrgDiff*pE->pairwise;
         pathRowPtr[16] = tmp->expNRGdiff;
         pathRowPtr[17] = tmp->expNRGdiff*pE->omega[0];

//       }
       tmp = tmp->DpathPtr;  PathPtr1++;
     }
   }
}
//******************************************************************************
void NRGposterior(Bayes *pE, Information &C, Interaction *FO, int q, int count,
        double ***RateNRG, ofstream &o)
{
   void PrintRateAway(Bayes *pE, Interaction *FO, int *tmpAAseq,int whichBranch,
        int *siteNeigh, int AAseqLen,ofstream& o1);
   int NeighborAccess(int AA, int *c, int **nAccess);
   int AAseqLen = C.AAlen;
   int nucSeqLen = C.len;
   static int *nucRate, *AArate, *neighRate;
   if(!count)
   {
     nucRate   = new int[nucSeqLen];
     AArate    = new int[AAseqLen];
     neighRate = new int[AAseqLen];
   }
   ////////////////////////////////////////////////////////////////////
   // Here we talk about the NRG associated with the posterior paths //
   ////////////////////////////////////////////////////////////////////

     //cout << " we made it through " << count%10 << endl;
     int **FO_nAccess = FO->nAccess;
     int *S1, *S2, *S3, *S4, *S5, *S6;
     S1 = nucRate; S2 = C.parentSeq;
     S3 = AArate; S4 = C.parentAAseq;
     S5 = neighRate; S6 = C.neighAcc;
     for(int d = 0; d < nucSeqLen; d++)
     {
       *S1 = *S2; S1++; S2++;
       if(d < AAseqLen)
       {
         *S3 = *S4; S3++; S4++;
         *S5 = *S6; S5++; S6++;
       }
     }
     double rateTime = 0.100;
     Substitution *rateSub = C.startPath;
     o << q << " " << count/10 << " " << 0.0 << " ";

     PrintRateAway(pE, FO, AArate, q, neighRate, AAseqLen, o);
     while(rateSub!=NULL)
     {
       if(rateSub->time <= rateTime)
       {
         int NUCcol = rateSub->column;
         int AAcol = NUCcol / 3;
         int c8[3];
         nucRate[NUCcol] = rateSub->pathNuc;// updates the nuc sequence
         AArate[AAcol] = rateSub->amino; // updates the AA sequence
         switch(NUCcol % 3)
         {
           case 0: c8[0] = nucRate[NUCcol]; c8[1] = nucRate[NUCcol+1];
                   c8[2] = nucRate[NUCcol+2]; break;
           case 1: c8[0] = nucRate[NUCcol-1]; c8[1] = nucRate[NUCcol];
                   c8[2] = nucRate[NUCcol+1]; break;
           case 2: c8[0] = nucRate[NUCcol-2]; c8[1] = nucRate[NUCcol-1];
                   c8[2] = nucRate[NUCcol];   break;
         }
         neighRate[AAcol]=NeighborAccess(rateSub->amino,c8,FO_nAccess);
         rateSub = rateSub->DpathPtr;
       }
       else
       {
         o << q << " " << count/10 << " " << rateTime << " ";

         PrintRateAway(pE, FO, AArate, q, neighRate, AAseqLen, o);
         rateTime += 0.100;
       }
     }
     while(rateTime < 1.05)
     {
         o << count/10 << " " << rateTime << " ";
         PrintRateAway(pE, FO, AArate, q, neighRate, AAseqLen, o);
         rateTime += 0.100;
     }

   ////////////////////////////////////////////////////////////////////
   // Here we STOP talking about the energies of the sequence  paths //
   ////////////////////////////////////////////////////////////////////
   if(count == SAMPLE_PATH)
   { //Complete all functions and erase all of the memory we borrowed!
     delete[] nucRate;   delete[] AArate;   delete[] neighRate;
   }
}
//******************************************************************************
void PrintRateAway(Bayes *pE, Interaction *FO, int *tmpAAseq, int whichBranch,
		int *siteNeigh, int AAseqLen, ofstream& o)
{
     double SiteNRGCalc(Interaction *FO, double &sNRG, int AAcol,
                int testAA, int *tmpAAseq);
     int *nIndex, siteAA, siteNuc, siteKappa, *currentAA;
     double siteNRG, sNRG, rate, siteRate, K, W, S, P, U;
     int *aaRowPtr, *nucRowPtr, *kappaRowPtr;

     U = pE->rate[whichBranch];
     S = pE->solvent; P = pE->pairwise;
     if(MULTI_K){K = pE->kappa[whichBranch];}
     else{K = pE->kappa[0];}
     if(MULTI_W){W = pE->omega[whichBranch];}
     else{W = pE->omega[0];}

     double *I_nf = pE->MCMCnf;
     int **FO_neighList = FO->neighList;
     currentAA = tmpAAseq;
     nIndex = siteNeigh;
    // run through all sites in protein
    for(int AAcol = 0; AAcol < AAseqLen; AAcol++)
    {
       siteRate = 0.0;
       // Obtain the correct row from Neighbor list
       // run through the 9 nearest neighbors
       aaRowPtr = FO_neighList[*nIndex];
       nucRowPtr = aaRowPtr + 9;
       kappaRowPtr = nucRowPtr + 9;
       for(int neighbor = 0; neighbor < 9; neighbor++)
       {
          rate = 0.0; //re-initialize the rate variable
          siteAA    = *aaRowPtr;
          siteNuc   = *nucRowPtr;
          siteKappa = *kappaRowPtr;
          if(siteAA == *currentAA) // If the change is a synonymous one
          {
             if(siteKappa == 1){rate = I_nf[siteNuc];}
             else{rate = I_nf[siteNuc] * K;}
          }
          else  // it is a nonsynonymous change
          {
             if(siteAA != 20) // if not a stop
             {
                siteNRG = 1.0;
                #if (SIM_W_S || SIM_W_P)
                  siteNRG = SiteNRGCalc(FO,sNRG,AAcol,siteAA,tmpAAseq);
                  siteNRG = exp(sNRG*S + siteNRG*P);
                //  siteNRG = *expPtr;
                #endif
                if(siteKappa == 1){rate = I_nf[siteNuc] * W * siteNRG;}
                else{rate = I_nf[siteNuc] * K * W * siteNRG;}
             }
          }
          siteRate += rate;
          //update the pointers to the structure to advance to the next entry
          aaRowPtr++;  nucRowPtr++;   kappaRowPtr++;
       } // end running through all of the nearest neighbors
       currentAA++;  nIndex++;
       o << U*siteRate << " ";
    } // End running through the entire sequence
    o << endl;
}  // End PrintRateAway Subroutine
//******************************************************************************
void FirstParam(Bayes *pE, Information *C, Information *P, Interaction *FO,
                        ofstream &o)
{
   double InitProbParam(Bayes *pE, Interaction *FO, Information &C, int which,
                ofstream &o);

   Substitution *tmp, *tmp1;
   int numBranch = pE->numBranch;
// CHANGE by JEFF   int AAseqLen = C[0].AAlen;

   for(int i = 0; i < numBranch; i++)
   {
      o << "Inside FirstParam, this is iteration " << i <<  endl;
      P[i].pathHood = C[i].pathHood = InitProbParam(pE, FO, C[i], i, o);
      o << "We are done with InitProbParam for path " << i << endl;
      P[i].probLastEvent = C[i].probLastEvent;
      pE->branchLen[i] = P[i].avgRate = C[i].avgRate;
      tmp = C[i].startPath;  tmp1 = P[i].startPath;
      // We now copy all of the sequence information into Proposed subs
      while(tmp != NULL){
         tmp1->nrgDiff     = tmp->nrgDiff;
         tmp1->expNRGdiff  = tmp->expNRGdiff;
         tmp1->solvNRGdiff = tmp->solvNRGdiff;
         tmp1->genHood = tmp->genHood;
         tmp1->subHood = tmp->subHood;
         tmp = tmp->DpathPtr; tmp1 = tmp1->DpathPtr;
      }
      o << "This is the likelihood of the path " << P[i].pathHood << endl;
   }
}  // End FirstParam Subroutine
//******************************************************************************
// codonFrequency fills the cf vector with the correct codon frequnencies
//      accounting for the stop codons.
void codonFrequency(double *cf, double *nfreq)
{
double stopfreq;

stopfreq = 1.0/(1.0 - (nfreq[3]*nfreq[0]*nfreq[0]+nfreq[3]*nfreq[0]*nfreq[2]+nfreq[3]*nfreq[2]*nfreq[0]));


   cf[0]=nfreq[2]*nfreq[1]*nfreq[0]*stopfreq;
   cf[1]=nfreq[2]*nfreq[1]*nfreq[1]*stopfreq;
   cf[2]=nfreq[2]*nfreq[1]*nfreq[2]*stopfreq;
   cf[3]=nfreq[2]*nfreq[1]*nfreq[3]*stopfreq;
   cf[4]=nfreq[0]*nfreq[2]*nfreq[0]*stopfreq;
   cf[5]=nfreq[0]*nfreq[2]*nfreq[2]*stopfreq;
   cf[6]=nfreq[1]*nfreq[2]*nfreq[0]*stopfreq;
   cf[7]=nfreq[1]*nfreq[2]*nfreq[1]*stopfreq;
   cf[8]=nfreq[1]*nfreq[2]*nfreq[2]*stopfreq;
   cf[9]=nfreq[1]*nfreq[2]*nfreq[3]*stopfreq;
   cf[10]=nfreq[0]*nfreq[0]*nfreq[1]*stopfreq;
   cf[11]=nfreq[0]*nfreq[0]*nfreq[3]*stopfreq;
   cf[12]=nfreq[2]*nfreq[0]*nfreq[1]*stopfreq;
   cf[13]=nfreq[2]*nfreq[0]*nfreq[3]*stopfreq;
   cf[14]=nfreq[3]*nfreq[2]*nfreq[1]*stopfreq;
   cf[15]=nfreq[3]*nfreq[2]*nfreq[3]*stopfreq;
   cf[16]=nfreq[1]*nfreq[0]*nfreq[0]*stopfreq;
   cf[17]=nfreq[1]*nfreq[0]*nfreq[2]*stopfreq;
   cf[18]=nfreq[2]*nfreq[0]*nfreq[0]*stopfreq;
   cf[19]=nfreq[2]*nfreq[0]*nfreq[2]*stopfreq;
   cf[20]=nfreq[2]*nfreq[2]*nfreq[0]*stopfreq;
   cf[21]=nfreq[2]*nfreq[2]*nfreq[1]*stopfreq;
   cf[22]=nfreq[2]*nfreq[2]*nfreq[2]*stopfreq;
   cf[23]=nfreq[2]*nfreq[2]*nfreq[3]*stopfreq;
   cf[24]=nfreq[1]*nfreq[0]*nfreq[1]*stopfreq;
   cf[25]=nfreq[1]*nfreq[0]*nfreq[3]*stopfreq;
   cf[26]=nfreq[0]*nfreq[3]*nfreq[0]*stopfreq;
   cf[27]=nfreq[0]*nfreq[3]*nfreq[1]*stopfreq;
   cf[28]=nfreq[0]*nfreq[3]*nfreq[3]*stopfreq;
   cf[29]=nfreq[1]*nfreq[3]*nfreq[0]*stopfreq;
   cf[30]=nfreq[1]*nfreq[3]*nfreq[1]*stopfreq;
   cf[31]=nfreq[1]*nfreq[3]*nfreq[2]*stopfreq;
   cf[32]=nfreq[1]*nfreq[3]*nfreq[3]*stopfreq;
   cf[33]=nfreq[3]*nfreq[3]*nfreq[0]*stopfreq;
   cf[34]=nfreq[3]*nfreq[3]*nfreq[2]*stopfreq;
   cf[35]=nfreq[0]*nfreq[0]*nfreq[0]*stopfreq;
   cf[36]=nfreq[0]*nfreq[0]*nfreq[2]*stopfreq;
   cf[37]=nfreq[0]*nfreq[3]*nfreq[2]*stopfreq;
   cf[38]=nfreq[3]*nfreq[3]*nfreq[1]*stopfreq;
   cf[39]=nfreq[3]*nfreq[3]*nfreq[3]*stopfreq;
   cf[40]=nfreq[1]*nfreq[1]*nfreq[0]*stopfreq;
   cf[41]=nfreq[1]*nfreq[1]*nfreq[1]*stopfreq;
   cf[42]=nfreq[1]*nfreq[1]*nfreq[2]*stopfreq;
   cf[43]=nfreq[1]*nfreq[1]*nfreq[3]*stopfreq;
   cf[44]=nfreq[0]*nfreq[2]*nfreq[1]*stopfreq;
   cf[45]=nfreq[0]*nfreq[2]*nfreq[3]*stopfreq;
   cf[46]=nfreq[3]*nfreq[1]*nfreq[0]*stopfreq;
   cf[47]=nfreq[3]*nfreq[1]*nfreq[1]*stopfreq;
   cf[48]=nfreq[3]*nfreq[1]*nfreq[2]*stopfreq;
   cf[49]=nfreq[3]*nfreq[1]*nfreq[3]*stopfreq;
   cf[50]=nfreq[0]*nfreq[1]*nfreq[0]*stopfreq;
   cf[51]=nfreq[0]*nfreq[1]*nfreq[1]*stopfreq;
   cf[52]=nfreq[0]*nfreq[1]*nfreq[2]*stopfreq;
   cf[53]=nfreq[0]*nfreq[1]*nfreq[3]*stopfreq;
   cf[54]=nfreq[3]*nfreq[2]*nfreq[2]*stopfreq;
   cf[55]=nfreq[3]*nfreq[0]*nfreq[1]*stopfreq;
   cf[56]=nfreq[3]*nfreq[0]*nfreq[3]*stopfreq;
   cf[57]=nfreq[2]*nfreq[3]*nfreq[0]*stopfreq;
   cf[58]=nfreq[2]*nfreq[3]*nfreq[1]*stopfreq;
   cf[59]=nfreq[2]*nfreq[3]*nfreq[2]*stopfreq;
   cf[60]=nfreq[2]*nfreq[3]*nfreq[3]*stopfreq;
} // End of Codon Frequency function
//******************************************************************************
void YangCodonRate(double *cfreq, double *c, double k, double w)
{   // This functions computes the rate away for each of the 61 codons according
//      to the Goldman and Yang rate matrix.  It uses the kappa, omega and the
//      codon frequencies computed in the previous function.
   cfreq[0]=  c[1]+c[3]+k*(c[2])+w*(c[18]+c[20]+c[40]+c[46]+k*(c[50]+c[57]));
   cfreq[1]=  c[0]+c[2]+k*(c[3])+w*(c[12]+c[21]+c[41]+c[47]+k*(c[51]+c[58]));
   cfreq[2]=  c[1]+c[3]+k*(c[0])+w*(c[19]+c[22]+c[42]+c[48]+k*(c[52]+c[59]));
   cfreq[3]=  c[0]+c[2]+k*(c[1])+w*(c[13]+c[23]+c[43]+c[49]+k*(c[53]+c[60]));
   cfreq[4]=  c[6]+k*(c[5])+w*(c[26]+c[44]+c[45]+c[50]+k*(c[20]+c[35]));
   cfreq[5]=  c[8]+k*(c[4])+w*(c[37]+c[44]+c[45]+c[52]+c[54]+k*(c[22]+c[36]));
   cfreq[6]=  c[4]+c[7]+c[9]+k*(c[8])+w*(c[20]+c[29]+c[40]+k*(c[16]));
   cfreq[7]=  c[6]+c[8]+k*(c[9])+w*(c[21]+c[30]+c[41]+c[44]+k*(c[14]+c[24]));
   cfreq[8]=  c[5]+c[7]+c[9]+k*(c[6])+w*(c[22]+c[31]+c[42]+k*(c[17]+c[54]));
   cfreq[9]=  c[6]+c[8]+k*(c[7])+w*(c[23]+c[32]+c[43]+c[45]+k*(c[15]+c[25]));
   cfreq[10]=  k*(c[11])+w*(c[24]+c[27]+c[35]+c[36]+c[51]+c[55]+k*(c[12]+c[44]));
   cfreq[11]=  k*(c[10])+w*(c[25]+c[28]+c[35]+c[36]+c[53]+c[56]+k*(c[13]+c[45]));
   cfreq[12]=  k*(c[13])+w*(c[1]+c[18]+c[19]+c[24]+c[55]+c[58]+k*(c[10]+c[21]));
   cfreq[13]=  k*(c[12])+w*(c[3]+c[18]+c[19]+c[25]+c[56]+c[60]+k*(c[11]+c[23]));
   cfreq[14]=  k*(c[15])+w*(c[21]+c[38]+c[44]+c[47]+c[54]+k*(c[7]+c[55]));
   cfreq[15]=  k*(c[14])+w*(c[23]+c[39]+c[45]+c[49]+c[54]+k*(c[9]+c[56]));
   cfreq[16]=  k*(c[17])+w*(c[18]+c[24]+c[25]+c[29]+c[35]+c[40]+k*(c[6]));
   cfreq[17]=  k*(c[16])+w*(c[19]+c[24]+c[25]+c[31]+c[36]+c[42]+k*(c[8]));
   cfreq[18]=  k*(c[19])+w*(c[0]+c[12]+c[13]+c[16]+c[57]+k*(c[20]+c[35]));
   cfreq[19]=  k*(c[18])+w*(c[2]+c[12]+c[13]+c[17]+c[59]+k*(c[22]+c[36]));
   cfreq[20]=  c[21]+c[23]+k*(c[22])+w*(c[0]+c[6]+c[57]+k*(c[4]+c[18]));
   cfreq[21]=  c[20]+c[22]+k*(c[23])+w*(c[1]+c[7]+c[14]+c[58]+k*(c[12]+c[44]));
   cfreq[22]=  c[21]+c[23]+k*(c[20])+w*(c[2]+c[8]+c[54]+c[59]+k*(c[5]+c[19]));
   cfreq[23]=  c[20]+c[22]+k*(c[21])+w*(c[3]+c[9]+c[15]+c[60]+k*(c[13]+c[45]));
   cfreq[24]=  k*(c[25])+w*(c[10]+c[12]+c[16]+c[17]+c[30]+c[41]+k*(c[7]+c[55]));
   cfreq[25]=  k*(c[24])+w*(c[11]+c[13]+c[16]+c[17]+c[32]+c[43]+k*(c[9]+c[56]));
   cfreq[26]=  c[27]+c[28]+w*(c[4]+c[29]+c[33]+c[35]+k*(c[37]+c[50]+c[57]));
   cfreq[27]=  c[26]+k*(c[28])+w*(c[10]+c[30]+c[37]+c[38]+c[44]+k*(c[51]+c[58]));
   cfreq[28]=  c[26]+k*(c[27])+w*(c[11]+c[32]+c[37]+c[39]+c[45]+k*(c[53]+c[60]));
   cfreq[29]=  c[30]+c[32]+k*(c[31]+c[33])+w*(c[6]+c[16]+c[26]+c[57]+k*(c[40]));
   cfreq[30]=  c[29]+c[31]+k*(c[32])+w*(c[7]+c[24]+c[27]+c[58]+k*(c[38]+c[41]));
   cfreq[31]=  c[30]+c[32]+k*(c[29]+c[34])+w*(c[8]+c[17]+c[37]+c[59]+k*(c[42]));
   cfreq[32]=  c[29]+c[31]+k*(c[30])+w*(c[9]+c[25]+c[28]+c[60]+k*(c[39]+c[43]));
   cfreq[33]=  k*(c[29]+c[34])+w*(c[26]+c[38]+c[39]+c[57]+k*(c[46]));
   cfreq[34]=  k*(c[33]+c[31])+w*(c[37]+c[38]+c[39]+c[54]+c[59]+k*(c[48]));
   cfreq[35]=  k*(c[36])+w*(c[10]+c[11]+c[16]+c[26]+c[50]+k*(c[4]+c[18]));
   cfreq[36]=  k*(c[35])+w*(c[10]+c[11]+c[17]+c[37]+c[52]+k*(c[5]+c[19]));
   cfreq[37]=  w*(c[5]+c[27]+c[28]+c[31]+c[34]+c[36]+k*(c[26]+c[52]+c[59]));
   cfreq[38]=  k*(c[39])+w*(c[14]+c[27]+c[33]+c[34]+c[55]+c[58]+k*(c[30]+c[47]));
   cfreq[39]=  k*(c[38])+w*(c[15]+c[28]+c[33]+c[34]+c[56]+c[60]+k*(c[32]+c[49]));
   cfreq[40]=  c[41]+c[43]+k*(c[42])+w*(c[0]+c[6]+c[16]+c[50]+k*(c[29]+c[46]));
   cfreq[41]=  c[40]+c[42]+k*(c[43])+w*(c[1]+c[7]+c[24]+c[51]+k*(c[30]+c[47]));
   cfreq[42]=  c[41]+c[43]+k*(c[40])+w*(c[2]+c[8]+c[17]+c[52]+k*(c[31]+c[48]));
   cfreq[43]=  c[40]+c[42]+k*(c[41])+w*(c[3]+c[9]+c[25]+c[53]+k*(c[32]+c[49]));
   cfreq[44]=  k*(c[45])+w*(c[4]+c[7]+c[5]+c[14]+c[27]+c[51]+k*(c[10]+c[21]));
   cfreq[45]=  k*(c[44])+w*(c[4]+c[9]+c[5]+c[15]+c[28]+c[53]+k*(c[11]+c[23]));
   cfreq[46]=  c[47]+c[49]+k*(c[48])+w*(c[0]+c[50]+k*(c[33]+c[40]));
   cfreq[47]=  c[46]+c[48]+k*(c[49])+w*(c[1]+c[14]+c[51]+c[55]+k*(c[38]+c[41]));
   cfreq[48]=  c[47]+c[49]+k*(c[46])+w*(c[2]+c[52]+c[54]+k*(c[34]+c[42]));
   cfreq[49]=  c[46]+c[48]+k*(c[47])+w*(c[3]+c[15]+c[53]+c[56]+k*(c[39]+c[43]));
   cfreq[50]=  c[51]+c[53]+k*(c[52])+w*(c[4]+c[35]+c[40]+c[46]+k*(c[0]+c[26]));
   cfreq[51]=  c[50]+c[52]+k*(c[53])+w*(c[10]+c[41]+c[44]+c[47]+k*(c[1]+c[27]));
   cfreq[52]=  c[51]+c[53]+k*(c[50])+w*(c[5]+c[36]+c[42]+c[48]+k*(c[2]+c[37]));
   cfreq[53]=  c[50]+c[52]+k*(c[51])+w*(c[11]+c[43]+c[45]+c[49]+k*(c[3]+c[28]));
   cfreq[54]=  w*(c[5]+c[14]+c[15]+c[22]+c[34]+c[48]+k*(c[8]));
   cfreq[55]=  k*(c[56])+w*(c[10]+c[12]+c[38]+c[47]+k*(c[14]+c[24]));
   cfreq[56]=  k*(c[55])+w*(c[11]+c[13]+c[39]+c[49]+k*(c[15]+c[25]));
   cfreq[57]=  c[58]+c[60]+k*(c[59])+w*(c[18]+c[20]+c[29]+c[33]+k*(c[0]+c[26]));
   cfreq[58]=  c[57]+c[59]+k*(c[60])+w*(c[12]+c[21]+c[30]+c[38]+k*(c[1]+c[27]));
   cfreq[59]=  c[58]+c[60]+k*(c[57])+w*(c[19]+c[22]+c[31]+c[34]+k*(c[2]+c[37]));
   cfreq[60]=  c[57]+c[59]+k*(c[58])+w*(c[13]+c[23]+c[32]+c[39]+k*(c[3]+c[28]));
} // This ends Yang Codon Rate function
//******************************************************************************
double CalcYangRate(double bl,double *cf, double *cr,int AAlen, int *neigh)
{
   static int z = 0;
   static double *counter;
   if(z == 0)
   {
      counter=new double[61];
      z = 1;
   }

   double temp,rateAway;
   temp = rateAway = 0.0;

   for(int i = 0; i < 61; i++){counter[i]=0.0;}
   for(int i = 0; i < AAlen; i++){counter[*neigh]+=1.0; neigh++;}
   for(int i = 0; i < 61; i++){rateAway += (counter[i]*cf[i]*cr[i]);}
   temp = bl/rateAway;

   return temp;
}
//******************************************************************************
void Update_KWU(Information &P, Information &C, double *gP, double *sP)
{
   Substitution *L, *W;

   C.pathHood   = P.pathHood;
   L = C.startPath; W = P.startPath;

   while(L != NULL) // nothing with energies updates in this process
   {
      L->genHood = W->genHood = *gP;
      L->subHood = W->subHood = *sP;
      L = L->DpathPtr; W = W->DpathPtr; gP++; sP++;
   } // end while statement
}// End updating KAPPA PARAMETER
//*****************************************************************************
void NewKappa(Bayes *pE, Information *C, Information *P, Interaction *FO,
         double ****RateNRG, int &accept_K, double *cf, long *seed, ofstream &o)
{
   void NewParameter(long *seed, double& m_top, double& m_bot, double param,
               double& tempParam, double maxValue,double minValue,double delta);
   double FastProbParam(Bayes *pE, Interaction *FO, Information &Info, int round,
                int which, ofstream &o, double *genHood, double *subHood,
                double ***AllRates);
   void YangCodonRate(double *cfreq, double *c, double k, double w);
   double CalcYangRate(double bl,double *cf, double *cr,int AAlen, int *neigh);
// CHANGE above by JEFF
   void Update_KWU(Information &P, Information &C, double *gP, double *sP);
   double rnd(long* seed);

   double lnProb_bot, lnProb_top, metro_top, metro_bot, R, random;
   int numBranch = pE->numBranch;

   static int z = 0;
   static double **GEN_PROB;
   static double **SPEC_PROB;

   if(z == 0)
   {
      GEN_PROB = new double *[numBranch];
      SPEC_PROB= new double *[numBranch];
      z = 1;
   }

   if(MULTI_K)
   {
     for(int i = 0; i < numBranch; i++)
     {
       // The proposed kappa path has to be completely recalculated
       NewParameter(seed, metro_top, metro_bot, pE->kappa[i], pE->tempParameter,
                            K_MAX, K_MIN, DELTA_K);
       lnProb_bot = C[i].pathHood;
       GEN_PROB[i]  = new double[P[i].totNumSub];
       SPEC_PROB[i] = new double[P[i].totNumSub];
       // lnProb_top only depends on the new value of Kappa and not Rate Matrix
       lnProb_top = FastProbParam(pE,FO,P[i],0,i,o, GEN_PROB[i], SPEC_PROB[i],RateNRG[i]);
       // R is the probability that we accept such a proposal
       if((lnProb_top - lnProb_bot) >= 100.0){R = exp(100.0);}
       else{R = exp(lnProb_top - lnProb_bot);}
       R *= metro_top/metro_bot;
       random = rnd(seed);
       if(random < R){   // The proposed value is chosen as the winner
          accept_K++;    // Tells us that Proposed Kappa is better
          pE->kappa[i] = pE->tempParameter;
          pE->branchLen[i] = C[i].avgRate  = P[i].avgRate;
          C[i].probLastEvent = P[i].probLastEvent;
          Update_KWU(P[i], C[i], GEN_PROB[i], SPEC_PROB[i]);
       }
       else // The current value is chosen as the winner in this case
       {   // All is clear on the Western Front
          P[i].pathHood = C[i].pathHood;
          P[i].avgRate  = C[i].avgRate;
          P[i].probLastEvent = C[i].probLastEvent;
       }
     }
   }
   else // Then there is only one kappa for all branches of the topology
   {
     // The proposed kappa path has to be completely recalculated
     NewParameter(seed, metro_top, metro_bot, pE->kappa[0], pE->tempParameter,
                         K_MAX, K_MIN, DELTA_K);
     double c_rate[61]; for(int zz = 0; zz < 61; zz++){c_rate[zz]=0.0;}
     YangCodonRate(c_rate,cf,pE->tempParameter,pE->omega[0]);

     double lnR = 0.0;
     for(int i = 0; i < numBranch; i++)
     {
       lnProb_bot = C[i].pathHood;
       GEN_PROB[i]  = new double[P[i].totNumSub];
       SPEC_PROB[i] = new double[P[i].totNumSub];
       // lnProb_top only depends on the new value of Kappa and not Rate Matrix
       pE->tempRate[i] = pE->rate[i];
       pE->rate[i] = CalcYangRate(pE->yangBranch[i],cf,c_rate,P[i].AAlen,P[i].neighAcc);
       lnProb_top = FastProbParam(pE,FO,P[i],0,i,o, GEN_PROB[i], SPEC_PROB[i],RateNRG[i]);
       lnR += lnProb_top - lnProb_bot;
//       if((lnProb_top - lnProb_bot) >= 100.0){R *= exp(100.0);}
//       else{R *= exp(lnProb_top - lnProb_bot);}
     }
     // R is the probability that we accept such a proposal in Kappa
//     R *= (metro_top/metro_bot);
// Change by DOUG and JEFF
     lnR += log(metro_top/metro_bot);
     if(lnR > 0.0)
     {
        R = 1.0;
     }
     else
     {
        R = exp(lnR);
     }
     random = rnd(seed);
     if(random < R){   // The proposed value is chosen as the winner
        accept_K++;    // Tells us that Proposed Kappa is better
        pE->kappa[0] = pE->tempParameter;
        for(int i = 0; i < numBranch; i++)
        {
           pE->branchLen[i] = C[i].avgRate  = P[i].avgRate;
           C[i].probLastEvent = P[i].probLastEvent;
           Update_KWU(P[i], C[i], GEN_PROB[i], SPEC_PROB[i]);
        }
     }
     else{ // The current value is chosen as the winner in this case
        for(int i = 0; i < numBranch; i++) // All is clear on the Western Front
        {
           P[i].pathHood = C[i].pathHood;
           P[i].avgRate  = C[i].avgRate;
           P[i].probLastEvent = C[i].probLastEvent;
           pE->rate[i] = pE->tempRate[i];
        }
     }
   }

   for(int i = 0; i < numBranch; i++)
   {
      delete[] GEN_PROB[i];   delete[] SPEC_PROB[i];
   }
}  // End NewKappa Subroutine
//******************************************************************************
void NewOmega(Bayes *pE, Information *C,Information *P,Interaction *FO,
         double ****RateNRG, int &accept_W, double *cf, long *seed, ofstream &o)
{
   void NewParameter(long *seed, double& m_top, double& m_bot, double param,
               double& tempParam, double maxValue,double minValue,double delta);
   double FastProbParam(Bayes *pE, Interaction *FO, Information &Info, int round,
                int which, ofstream &o, double *genHood, double *subHood,
                double ***AllRates);
   void Update_KWU(Information &P, Information &C, double *gP, double *sP);
   double rnd(long* seed);
   void YangCodonRate(double *cfreq, double *c, double k, double w);
   double CalcYangRate(double bl,double *cf, double *cr,int AAlen, int *neigh);

   double lnProb_bot, lnProb_top, metro_top, metro_bot, R, random;
   int numBranch = pE->numBranch;

   static int z = 0;
   static double **GEN_PROB;
   static double **SPEC_PROB;

   if(z == 0)
   {
      GEN_PROB = new double *[numBranch];
      SPEC_PROB= new double *[numBranch];
      z = 1;
   }
   if(MULTI_W)
   {
     for(int i = 0; i < numBranch; i++)
     {
       // The proposed kappa path has to be completely recalculated
       NewParameter(seed, metro_top, metro_bot, pE->omega[i], pE->tempParameter,
                            W_MAX, W_MIN, DELTA_W);
       lnProb_bot = C[i].pathHood;
       GEN_PROB[i]  = new double[P[i].totNumSub];
       SPEC_PROB[i] = new double[P[i].totNumSub];
       // lnProb_top only depends on the new value of Omega and not Rate Matrix
       lnProb_top = FastProbParam(pE,FO,P[i],2,i,o,GEN_PROB[i],SPEC_PROB[i],RateNRG[i]);
       // R is the probability that we accept such a proposal
       if((lnProb_top - lnProb_bot) >= 100.0){R = exp(100.0);}
       else{R = exp(lnProb_top - lnProb_bot);}
       R *= metro_top/metro_bot;
       random = rnd(seed);
       if(random < R){   // The proposed value is chosen as the winner
          accept_W++;    // Tells us that Proposed Omega is better
          pE->omega[i] = pE->tempParameter;
          pE->branchLen[i] = C[i].avgRate  = P[i].avgRate;
          C[i].probLastEvent = P[i].probLastEvent;
          Update_KWU(P[i], C[i], GEN_PROB[i], SPEC_PROB[i]);
       }
       else // The current value is chosen as the winner in this case
       {   // All is clear on the Western Front
          P[i].pathHood = C[i].pathHood;
          P[i].avgRate  = C[i].avgRate;
          P[i].probLastEvent = C[i].probLastEvent;
       }
     }
   }
   else // Then there is only one Omega for all branches of the topology
   {
     // The proposed Omega path has to be completely recalculated
     NewParameter(seed, metro_top, metro_bot, pE->omega[0], pE->tempParameter,
                         W_MAX, W_MIN, DELTA_W);
     double c_rate[61]; for(int zz = 0; zz < 61; zz++){c_rate[zz]=0.0;}
     YangCodonRate(c_rate,cf,pE->kappa[0],pE->tempParameter);

//     R = 1.0;
     double lnR = 0.0;
     for(int i = 0; i < numBranch; i++)
     {
       lnProb_bot = C[i].pathHood;
       GEN_PROB[i]  = new double[P[i].totNumSub];
       SPEC_PROB[i] = new double[P[i].totNumSub];
       // lnProb_top only depends on the new value of Omega and not Rate Matrix
       pE->tempRate[i] = pE->rate[i];
       pE->rate[i] = CalcYangRate(pE->yangBranch[i],cf,c_rate,P[i].AAlen,P[i].neighAcc);
       lnProb_top = FastProbParam(pE,FO,P[i],2,i,o, GEN_PROB[i], SPEC_PROB[i],RateNRG[i]);
       lnR += lnProb_top - lnProb_bot;
//       if((lnProb_top - lnProb_bot) >= 100.0){R *= exp(100.0);}
//       else{R *= exp(lnProb_top - lnProb_bot);}
     }
     // R is the probability that we accept such a proposal in Omega
//     R *= (metro_top/metro_bot);
// Change by DOUG and JEFF
     lnR += log(metro_top/metro_bot);
     if(lnR > 0.0)
     {
        R = 1.0;
     }
     else
     {
        R = exp(lnR);
     }
     random = rnd(seed);
     if(random < R){   // The proposed value is chosen as the winner
        accept_W++;    // Tells us that Proposed Omega is better
        pE->omega[0] = pE->tempParameter;
        for(int i = 0; i < numBranch; i++)
        {
           pE->branchLen[i] = C[i].avgRate  = P[i].avgRate;
           C[i].probLastEvent = P[i].probLastEvent;
           Update_KWU(P[i], C[i], GEN_PROB[i], SPEC_PROB[i]);
        }
     }
     else{ // The current value is chosen as the winner in this case
        for(int i = 0; i < numBranch; i++) // All is clear on the Western Front
        {
           P[i].pathHood = C[i].pathHood;
           P[i].avgRate  = C[i].avgRate;
           P[i].probLastEvent = C[i].probLastEvent;
           pE->rate[i] = pE->tempRate[i];
        }
     }
   }

   for(int i = 0; i < numBranch; i++)
   {
      delete[] GEN_PROB[i];   delete[] SPEC_PROB[i];
   }
}  // End NewOmega Subroutine
//*****************************************************************************/
void UpdateOmega_SP(Bayes *pE, Information *Winner, Information *Loser,
        int f, double **gP, double **sP)
{
   double *general, *specific;
//CHANGE by JEFF   double *L_EXP, *W_EXP, *general, *specific;
   Substitution *L, *W;
   int numBranch, AAseqLen; AAseqLen = Winner[0].AAlen;
   numBranch = pE->numBranch;


   for(int i = 0; i < numBranch; i++)
   {
      pE->branchLen[i] = Loser[i].avgRate = Winner[i].avgRate;
      Loser[i].pathHood = Winner[i].pathHood;
      Loser[i].probLastEvent = Winner[i].probLastEvent;

      L = Loser[i].startPath; W = Winner[i].startPath;
// CHANGE by JEFF      int count = 0;
      general = gP[i];
      specific = sP[i];

      if(!f)
      {
         while(L != NULL)
         {
            L->expNRGdiff  = W->expNRGdiff;
            L = L->DpathPtr; W = W->DpathPtr;
         }
      }
      else
      {
         while(L != NULL)  // For the case when f == 1
         { // nrgDiff does NOT change, do NOT update
           // expNRGdiff only changes with Omega because it is in the formulation
            L->expNRGdiff  = W->expNRGdiff;
            L->genHood = W->genHood = *general;
            L->subHood = W->subHood = *specific;
            general++; specific++;
            L = L->DpathPtr; W = W->DpathPtr;
         }
      }
   }
}// End updating OMEGA SOLVENT OR PAIRWISE PARAMETER
//******************************************************************************
void NewOmega_S(Bayes *pE, Information *C, Information *P, Interaction *FO,
                double ****RateNRG, int &accept_WS, long *seed, ofstream &o)
{
   void NewParameter(long *seed, double& m_top, double& m_bot, double param,
               double& tempParam, double maxValue,double minValue,double delta);
   double FastProbParam(Bayes *pE, Interaction *FO, Information &Info, int round,
                int which, ofstream &o, double *genHood, double *subHood,
                double ***AllRates);
   void UpdateOmega_SP(Bayes *pE, Information *Winner, Information *Loser,
        int f, double **gP, double **sP);
   double CalcConst(Bayes *pE, Information &P, Interaction *FO, double &maximum,
                long *seed, int round, ofstream &o);
   double rnd(long* seed);

   static double S_MULT = SIZEOFGRID / (W_S_MAX - W_S_MIN);
   double lnProb_bot,lnProb_top,metro_top,metro_bot,R,random,const_bot;
   double const_top, mean_S, originalConst, originalMax;
   int originalPoint, tmpGridPt;
   int numBranch = pE->numBranch;

   static int z = 0;
   static double **GEN_PROB;
   static double **SPEC_PROB;

   if(z == 0)
   {
      GEN_PROB = new double *[numBranch];
      SPEC_PROB= new double *[numBranch];
      z = 1;
   }

   NewParameter(seed, metro_top, metro_bot, pE->solvent, pE->tempParameter,
                W_S_MAX, W_S_MIN, DELTA_W_S);

   originalPoint = pE->S_gridPt;
   originalConst = pE->nrgConst;
   originalMax   = pE->maxValue;
   mean_S = 0.5 * (pE->tempParameter + pE->solvent);
   tmpGridPt = (int)(floor((mean_S - W_S_MIN) * S_MULT));

 //  R = 1.0;
   double lnR = 0.0;
   for(int i = 0; i < numBranch; i++)
   {
      lnProb_bot = C[i].pathHood;
      GEN_PROB[i]  = new double[P[i].totNumSub];
      SPEC_PROB[i] = new double[P[i].totNumSub];
   // The lnProb_top only depends on the new value of Kappa and not Rate Matrix
      lnProb_top =  FastProbParam(pE,FO,P[i],3,i, o, GEN_PROB[i],
                        SPEC_PROB[i], RateNRG[i]);
   // R is the probability that we accept such a proposal in solvent
//      if((lnProb_top - lnProb_bot) >= 100.0){R *= exp(100.0);}
//      else{R *= exp(lnProb_top - lnProb_bot);}
      lnR += lnProb_top - lnProb_bot;
   }
   if(tmpGridPt != pE->S_gridPt){
      pE->S_gridPt = tmpGridPt;// Now the grid point is possibly different
      // CalcConst changes the nrgConst as well as the maxvalue
      const_top = pE->nrgConst = CalcConst(pE,C[C[0].numBranch],FO,pE->maxValue,seed,9,o);
   }
   else{const_top = pE->nrgConst;}

   pE->tempS_gridPt = tmpGridPt;
   const_bot=pE->tempNRGconst = CalcConst(pE,P[C[0].numBranch],FO,pE->tempMAXvalue,seed,3,o);

//   R *=((metro_top*const_top)/(metro_bot*const_bot))*
//         exp(pE->solvNRG*(pE->tempParameter - pE->solvent) +
//         (pE->maxValue - pE->tempMAXvalue));

// Change by DOUG and JEFF
    lnR += log((metro_top*const_top)/(metro_bot*const_bot));
    lnR += pE->solvNRG*(pE->tempParameter - pE->solvent) +
         (pE->maxValue - pE->tempMAXvalue);

    if(lnR > 0.0)
    {
       R = 1.0;
    }
    else
    {
       R = exp(lnR);
    }

   random = rnd(seed);
   int f; // The proposed value is chosen as the winner
   if(random < R)
   {
      accept_WS++; f=1;
      pE->solvent  = pE->tempParameter;
      pE->nrgConst = pE->tempNRGconst;
      pE->maxValue = pE->tempMAXvalue;
      UpdateOmega_SP(pE ,P,C,f,GEN_PROB,SPEC_PROB);

   }
   else{// The current value is chosen as the winner in this case
      pE->nrgConst = pE->tempNRGconst = originalConst;
      pE->maxValue = pE->tempMAXvalue = originalMax;
      pE->S_gridPt = pE->tempS_gridPt = originalPoint;
      f = 0;
      UpdateOmega_SP(pE,C, P, f,GEN_PROB,SPEC_PROB);

   }
   for(int i = 0; i < numBranch; i++)
   {
      delete[] GEN_PROB[i];   delete[] SPEC_PROB[i];
   }
}  // End NewOmega_S Subroutine
//******************************************************************************
void NewOmega_P(Bayes *pE, Information *C, Information *P, Interaction *FO,
                double ****RateNRG, int &accept_WP, long *seed, ofstream &o)
{
   void NewParameter(long *seed, double& m_top, double& m_bot, double param,
               double& tempParam, double maxValue,double minValue,double delta);
   double FastProbParam(Bayes *pE, Interaction *FO, Information &Info, int round,
                int which, ofstream &o, double *genHood, double *subHood,
                double ***AllRates);
   void UpdateOmega_SP(Bayes *pE, Information *Winner, Information *Loser,
        int f, double **gP, double **sP);
   double CalcConst(Bayes *pE, Information &P, Interaction *FO, double &maximum,
                long *seed, int round, ofstream &o);
   double rnd(long* seed);

   static double P_MULT = SIZEOFGRID / (W_P_MAX - W_P_MIN);
   double lnProb_bot,lnProb_top,metro_top,metro_bot,R,random,const_bot;
   double const_top, mean_P, originalConst, originalMax;
   int originalPoint, tmpGridPt;
   int numBranch = pE->numBranch;

   static int z = 0;
   static double **GEN_PROB;
   static double **SPEC_PROB;

   if(z == 0)
   {
      GEN_PROB = new double *[numBranch];
      SPEC_PROB= new double *[numBranch];
      z = 1;
   }

   NewParameter(seed,metro_top,metro_bot,pE->pairwise, pE->tempParameter,
                W_P_MAX, W_P_MIN, DELTA_W_P);

   originalPoint = pE->P_gridPt;
   originalConst = pE->nrgConst;
   originalMax   = pE->maxValue;
   mean_P = 0.5 * (pE->tempParameter + pE->pairwise);
   tmpGridPt = (int)(floor((mean_P - W_P_MIN) * P_MULT));

//   R = 1.0;
   double lnR= 0.0;
   for(int i = 0; i < numBranch; i++)
   {
      lnProb_bot = C[i].pathHood;
      GEN_PROB[i]  = new double[P[i].totNumSub];
      SPEC_PROB[i] = new double[P[i].totNumSub];
   // The lnProb_top only depends on the new value of Kappa and not Rate Matrix
      lnProb_top =  FastProbParam(pE,FO,P[i],4,i, o, GEN_PROB[i],
                        SPEC_PROB[i], RateNRG[i]);
   // R is the probability that we accept such a proposal in solvent
//      if((lnProb_top - lnProb_bot) >= 100.0){R *= exp(100.0);}
//      else{R *= exp(lnProb_top - lnProb_bot);}
      lnR += lnProb_top - lnProb_bot;

   }
   if(tmpGridPt != pE->P_gridPt){
      pE->P_gridPt = tmpGridPt;// Now the grid point is possibly different
      // CalcConst changes the nrgConst as well as the maxvalue
      const_top = pE->nrgConst = CalcConst(pE,C[C[0].numBranch],FO,pE->maxValue, seed,9,o);
   }
   else{const_top = pE->nrgConst;}

   pE->tempP_gridPt = tmpGridPt;
   const_bot = pE->tempNRGconst = CalcConst(pE,P[C[0].numBranch],FO,pE->tempMAXvalue,seed,4,o);

// Change by DOUG and JEFF
     lnR += log((metro_top*const_top)/(metro_bot*const_bot));
     lnR += pE->pairNRG*(pE->tempParameter - pE->pairwise) +
             (pE->maxValue - pE->tempMAXvalue);
     if(lnR > 0.0)
     {
        R = 1.0;
     }
     else
     {
        R = exp(lnR);
     }

//   R *=((metro_top*const_top)/(metro_bot*const_bot))*
//         exp(pE->pairNRG*(pE->tempParameter - pE->pairwise) +
//         (pE->maxValue - pE->tempMAXvalue));

   random = rnd(seed);
   int f; // The proposed value is chosen as the winner
   if(random < R)
   {
      accept_WP++; f=1;
      pE->pairwise  = pE->tempParameter;
      pE->nrgConst  = pE->tempNRGconst;
      pE->maxValue  = pE->tempMAXvalue;
      UpdateOmega_SP(pE,P,C,f,GEN_PROB,SPEC_PROB);

   }
   else{// The current value is chosen as the winner in this case
      pE->nrgConst = pE->tempNRGconst = originalConst;
      pE->maxValue = pE->tempMAXvalue = originalMax;
      pE->P_gridPt = pE->tempP_gridPt = originalPoint;
      f = 0;
      UpdateOmega_SP(pE,C, P, f,GEN_PROB,SPEC_PROB);
   }
   for(int i = 0; i < numBranch; i++)
   {
      delete[] GEN_PROB[i];   delete[] SPEC_PROB[i];
   }
}  // End NewOmega_P Subroutine
//******************************************************************************
void UpdatePi(Bayes *pE,Information *P,Information *C,double **gP,double **sP)
{
   pE->MCMCnf[0] = pE->tempMCMCnf[0];
   pE->MCMCnf[1] = pE->tempMCMCnf[1];
   pE->MCMCnf[2] = pE->tempMCMCnf[2];
   pE->MCMCnf[3] = pE->tempMCMCnf[3];
   pE->nrgConst = pE->tempNRGconst;
   pE->maxValue = pE->tempMAXvalue;

   int count;
   Substitution *L;   Substitution *W;
   double *general, *specific;

   for(int i = 0; i < pE->numBranch; i++)
   {
      count = 0;
      C[i].pathHood = P[i].pathHood;
      pE->branchLen[i] = C[i].avgRate  = P[i].avgRate;
      C[i].probLastEvent = P[i].probLastEvent;

      L = C[i].startPath; W = P[i].startPath;
      general  = gP[i];
      specific = sP[i];

      while(L != NULL) // nothing with energies updates in this process
      {
         L->genHood = W->genHood = *general;
         L->subHood = W->subHood = *specific;
         L = L->DpathPtr; W = W->DpathPtr; general++; specific++;
      } // end while statement
   } // end looping through the branches
}// End updating PI PARAMETERS
//******************************************************************************
void NewPi(Bayes *pE, Information *C,Information *P,Interaction *FO,
      double ****RateNRG, int &accept_pi, int whichNuc, long *seed, ofstream &o)
{
   void NewParameter(long *seed, double& m_top, double& m_bot, double param,
               double& tempParam, double maxValue,double minValue,double delta);
   double FastProbParam(Bayes *pE, Interaction *FO, Information &Info, int round,
                int which, ofstream &o, double *genHood, double *subHood,
                double ***AllRates);
   void UpdatePi(Bayes *pE,Information *P,Information *C,double **gP,double **sP);
   double CalcConst(Bayes *pE, Information &P, Interaction *FO, double &maximum,
                long *seed, int round, ofstream &o);
   double Maxof3(double a, double b, double c);
   double rnd(long* seed);
   void codonFrequency(double *cf, double *nfreq);
   void YangCodonRate(double *cfreq, double *c, double k, double w);
   double CalcYangRate(double bl,double *cf, double *cr,int AAlen, int *neigh);

   static double A_MULT = (FREQGRIDSIZE-1) / (HIGH_A - LOW_A);
   static double HALF_A = 0.5/A_MULT;
   static double C_MULT = (FREQGRIDSIZE-1) / (HIGH_C - LOW_C);
   static double HALF_C = 0.5/C_MULT;
   static double G_MULT = (FREQGRIDSIZE-1) / (HIGH_G - LOW_G);
   static double HALF_G = 0.5/G_MULT;
   int number_T = C[0].len - pE->number_A - pE->number_C - pE->number_G;

   double lnProb_bot,lnProb_top,metro_top,metro_bot,R,random,const_bot,max_freq,
      const_top, mean_A,mean_C,mean_G, originalConst, originalMax, nfdivision,
      cf[61], c_rate[61];
   int originalPoint[3], chosen_A, chosen_C, chosen_G;
   int numBranch = pE->numBranch;

   static int z = 0;
   static double **GEN_PROB;
   static double **SPEC_PROB;

   if(z == 0)
   {
      GEN_PROB = new double *[numBranch];
      SPEC_PROB= new double *[numBranch];
      z = 1;
   }

  // We iterate this 5 times to "help" it converge faster

   NewParameter(seed,metro_top,metro_bot,pE->MCMCnf[whichNuc],
                pE->tempMCMCnf[whichNuc],1.0, 0.0,DELTA_PI);
   // UPDATE ALL OTHER FREQUENCIES
//   metro_top /= ((1.0-pE->MCMCnf[whichNuc])*(1.0-pE->MCMCnf[whichNuc]));
//   metro_bot /= ((1.0-pE->tempMCMCnf[whichNuc])*(1.0-pE->tempMCMCnf[whichNuc]));
////////////////////////////////////////////////////////
//   Above is the correct implementation using        //
//   division.  Taking out the division and replacing //
//   it with multiplication (faster) is given below.  //
//              Both should be equivalent             //
////////////////////////////////////////////////////////
   metro_bot *= ((1.0-pE->MCMCnf[whichNuc])*(1.0-pE->MCMCnf[whichNuc]));
   metro_top *= ((1.0-pE->tempMCMCnf[whichNuc])*(1.0-pE->tempMCMCnf[whichNuc]));

   double pi1,pi2,pi3;
   pi1 = pE->tempMCMCnf[whichNuc];
   nfdivision =(1.0-pE->tempMCMCnf[whichNuc]) / (1.0-pE->MCMCnf[whichNuc]);
   int dummy = (whichNuc+1)%4;
   pi2 = pE->tempMCMCnf[dummy] = pE->MCMCnf[dummy] * nfdivision;
   dummy = (whichNuc+2)%4;
   pi3 = pE->tempMCMCnf[dummy] = pE->MCMCnf[dummy] * nfdivision;
   dummy = (whichNuc+3)%4;
   pE->tempMCMCnf[dummy] = 1.0 - (pi1 + pi2 + pi3);

   originalPoint[0] = pE->A_gridPt;
   originalPoint[1] = pE->C_gridPt;
   originalPoint[2] = pE->G_gridPt;
   originalConst = pE->nrgConst;
   originalMax   = pE->maxValue;

   mean_A = 0.5 * (pE->tempMCMCnf[0] + pE->MCMCnf[0]);
   mean_C = 0.5 * (pE->tempMCMCnf[1] + pE->MCMCnf[1]);
   mean_G = 0.5 * (pE->tempMCMCnf[2] + pE->MCMCnf[2]);

   if(mean_A < (LOW_A+HALF_A)){chosen_A = 0;}
   else if(mean_A > (HIGH_A-HALF_A)){chosen_A = FREQGRIDSIZE-1;}
   else{chosen_A = (int)(1.0+(floor((mean_A - (LOW_A+HALF_A)) * A_MULT)));}

   if(mean_C < (LOW_C+HALF_C)){chosen_C = 0;}
   else if(mean_C > (HIGH_C-HALF_C)){chosen_C = FREQGRIDSIZE-1;}
   else{chosen_C = (int)(1.0+(floor((mean_C - (LOW_C+HALF_C)) * C_MULT)));}

   if(mean_G < (LOW_G+HALF_G)){chosen_G = 0;}
   else if(mean_G > (HIGH_G-HALF_G)){chosen_G = FREQGRIDSIZE-1;}
   else{chosen_G = (int)(1.0+(floor((mean_G - (LOW_G+HALF_G)) * G_MULT)));}

   double pi_A, pi_C, pi_G;
   pi_A = A_POINTS[chosen_A];
   pi_C = C_POINTS[chosen_C];
   pi_G = G_POINTS[chosen_G];

   while((pi_A+pi_C+pi_G) > 0.99)
   {
      max_freq = Maxof3(pi_A,pi_C,pi_G);
      if(max_freq == pi_A){chosen_A--; pi_A = A_POINTS[chosen_A];}
      else{
         if(max_freq == pi_C){chosen_C--; pi_C = C_POINTS[chosen_C];}
         else{chosen_G--; pi_G = G_POINTS[chosen_G];}
      }
   }

   pE->A_gridPt = chosen_A; // Now the grid point is possibly different
   pE->C_gridPt = chosen_C;
   pE->G_gridPt = chosen_G;
   if((pE->A_gridPt != originalPoint[0]) || (pE->C_gridPt != originalPoint[1])
              || (pE->G_gridPt != originalPoint[2]))
   {
      // CalcConst changes the nrgConst as well as the maxvalue
      const_top = pE->nrgConst = CalcConst(pE,C[C[0].numBranch],FO,pE->maxValue,seed,9,o);
   }
   else{const_top = pE->nrgConst;}

   pE->tempA_gridPt = chosen_A;
   pE->tempC_gridPt = chosen_C;
   pE->tempG_gridPt = chosen_G;
   const_bot = pE->tempNRGconst = CalcConst(pE,P[C[0].numBranch],FO,pE->tempMAXvalue,seed,5,o);

   double LnDiff_A, LnDiff_C, LnDiff_G, LnDiff_T;

   LnDiff_A = log(pE->tempMCMCnf[0] / pE->MCMCnf[0]);
   LnDiff_C = log(pE->tempMCMCnf[1] / pE->MCMCnf[1]);
   LnDiff_G = log(pE->tempMCMCnf[2] / pE->MCMCnf[2]);
   LnDiff_T = log(pE->tempMCMCnf[3] / pE->MCMCnf[3]);

   for(int zz = 0; zz < 61; zz++){c_rate[zz] = cf[zz]=0.0;}
   codonFrequency(cf, pE->tempMCMCnf);
   YangCodonRate(c_rate,cf,pE->kappa[0],pE->omega[0]);


//   R = 1.0;
   double lnR = 0.0;
   for(int i = 0; i < numBranch; i++)
   {
      lnProb_bot = C[i].pathHood;
      GEN_PROB[i]  = new double[P[i].totNumSub];
      SPEC_PROB[i] = new double[P[i].totNumSub];
   // The lnProb_top only depends on the new value of pi and not Rate Matrix
      pE->tempRate[i] = pE->rate[i];
      pE->rate[i] = CalcYangRate(pE->yangBranch[i],cf,c_rate,P[i].AAlen,P[i].neighAcc);
      lnProb_top =  FastProbParam(pE,FO,P[i],5,i, o, GEN_PROB[i],
                SPEC_PROB[i], RateNRG[i]);
   // R is the probability that we accept such a proposal in pi
//      if((lnProb_top - lnProb_bot) >= 100.0){R *= exp(100.0);}
//      else{R *= exp(lnProb_top - lnProb_bot);}
      lnR += lnProb_top - lnProb_bot;
   }

//   R *=((metro_top*const_top)/(metro_bot*const_bot))*
//      exp(pE->number_A*LnDiff_A + pE->number_C*LnDiff_C +
//      pE->number_G*LnDiff_G + number_T*LnDiff_T +
//      (pE->maxValue - pE->tempMAXvalue));

   lnR += log((metro_top*const_top)/(metro_bot*const_bot));
   lnR += pE->number_A*LnDiff_A + pE->number_C*LnDiff_C +
     pE->number_G*LnDiff_G + number_T*LnDiff_T +
     (pE->maxValue - pE->tempMAXvalue);
    if(lnR > 0.0)
    {
       R = 1.0;
    }
    else
    {
       R = exp(lnR);
    }
   random = rnd(seed);
    // The proposed value is chosen as the winner
   if(random < R){accept_pi++;UpdatePi(pE, P, C, GEN_PROB, SPEC_PROB);}
   else{// The current value is chosen as the winner in this case
      pE->tempMCMCnf[0] = pE->MCMCnf[0];
      pE->tempMCMCnf[1] = pE->MCMCnf[1];
      pE->tempMCMCnf[2] = pE->MCMCnf[2];
      pE->tempMCMCnf[3] = pE->MCMCnf[3];
      pE->A_gridPt = pE->tempA_gridPt = originalPoint[0];
      pE->C_gridPt = pE->tempC_gridPt = originalPoint[1];
      pE->G_gridPt = pE->tempG_gridPt = originalPoint[2];
      pE->nrgConst = pE->tempNRGconst = originalConst;
      pE->maxValue = pE->tempMAXvalue = originalMax;

      for(int i = 0; i < numBranch; i++)
      {
         P[i].pathHood = C[i].pathHood;
         P[i].avgRate  = C[i].avgRate;
         P[i].probLastEvent = C[i].probLastEvent;
         pE->rate[i] = pE->tempRate[i];
      }
   }
   for(int i = 0; i < numBranch; i++)
   {
      delete[] GEN_PROB[i];   delete[] SPEC_PROB[i];
   }
}  // End NewOmega_S Subroutine
/******************************************************************************
void UpdateRate(Information &P, Information &C, double *gP, double *sP)
{
   int count = 0;
   Substitution *L, *W;

   C.pathHood   = P.pathHood;
   L = C.startPath; W = P.startPath;

   while(L != NULL) // nothing with energies updates in this process
   {
      L->genHood = W->genHood = gP[count];
      L->subHood = W->subHood = sP[count];
      L = L->DpathPtr; W = W->DpathPtr; count++;
   } // end while statement
}// End updating RATE PARAMETER
//*****************************************************************************/
void NewRate(Bayes *pE, Information *C, Information *P, Interaction *FO,
                double ****RateNRG, int &accept_U, long *seed, ofstream &o)
{
   void NewParameter(long *seed, double& m_top, double& m_bot, double param,
               double& tempParam, double maxValue,double minValue,double delta);
   double FastProbParam(Bayes *pE, Interaction *FO, Information &Info, int round,
                int which, ofstream &o, double *genHood, double *subHood,
                double ***AllRates);
   void Update_KWU(Information &P, Information &C, double *gP, double *sP);
   double rnd(long* seed);
   void codonFrequency(double *cf, double *nfreq);
   void YangCodonRate(double *cfreq, double *c, double k, double w);
   double CalcYangRate(double bl,double *cf, double *cr,int AAlen, int *neigh);

   double lnProb_bot,lnProb_top,metro_top,metro_bot,R,random,cf[61],c_rate[61],u;
   int numBranch = pE->numBranch;

   static int z = 0;
   static double **GEN_PROB;
   static double **SPEC_PROB;

   if(z == 0)
   {
      GEN_PROB = new double *[numBranch];
      SPEC_PROB= new double *[numBranch];
      z = 1;
   }
   for(int zz = 0; zz < 61; zz++){c_rate[zz] = cf[zz]=0.0;}
   codonFrequency(cf, pE->MCMCnf);
   YangCodonRate(c_rate,cf,pE->kappa[0],pE->omega[0]);

   for(int i = 0; i < numBranch; i++)
   {
      // The proposed kappa path has to be completely recalculated
      NewParameter(seed, metro_top, metro_bot, pE->yangBranch[i],
                u, U_MAX, U_MIN, DELTA_U);
      lnProb_bot = C[i].pathHood;
      GEN_PROB[i]  = new double[P[i].totNumSub];
      SPEC_PROB[i] = new double[P[i].totNumSub];
   // The lnProb_top only depends on the new value of Kappa and not Rate Matrix
      pE->tempRate[i] = pE->rate[i];
      pE->tempParameter = CalcYangRate(u,cf,c_rate,P[i].AAlen,P[i].neighAcc);
      lnProb_top =  FastProbParam(pE,FO,P[i],1,i, o, GEN_PROB[i],
                        SPEC_PROB[i], RateNRG[i]);

      // R is the probability that we accept such a proposal
      if((lnProb_top - lnProb_bot) >= 100.0){R = exp(100.0);}
      else{R = exp(lnProb_top - lnProb_bot);}

      R *= metro_top/metro_bot;

      random = rnd(seed);

      if(random < R){   // The proposed value is chosen as the winner
         accept_U++;    // Tells us that Proposed Kappa is better
         pE->yangBranch[i] = u;
         pE->rate[i] = pE->tempParameter;
         pE->branchLen[i] = C[i].avgRate  = P[i].avgRate;
         C[i].probLastEvent = P[i].probLastEvent;
         Update_KWU(P[i], C[i], GEN_PROB[i], SPEC_PROB[i]);
      }
      else{ // The current value is chosen as the winner in this case
         // All is clear on the Western Front
         P[i].pathHood = C[i].pathHood;
         P[i].avgRate  = C[i].avgRate;
         P[i].probLastEvent = C[i].probLastEvent;
      }
   }
   for(int i = 0; i < numBranch; i++)
   {
      delete[] GEN_PROB[i];   delete[] SPEC_PROB[i];
   }
}  // End NewRate Subroutine
//******************************************************************************
void UpdatePath(Information &Win, Path *wP, Information &Lose, Path *lP,
      Substitution *start, int killCol, int flag, int LASTFLAG, double maxTime,
      ofstream &o)
{  // After a winner is chosen, we must replace the loser column with the winner
   void DeleteColumn(Information &Info, Path *columnPath, int col);
   void CopyUpdateColumn(Path *emptyP,Path *fullP,int col,Information &Empty);
   void UpdateStopColumn(Path *emptyP,Path *fullP,int col,Information &Empty);
   void SetCodonPath(Path *P, int col);

   Substitution *L, *W;

   int lNumsub, wNumsub, f, AAseqLen; f = 0;
   AAseqLen = Win.AAlen;
   lNumsub = lP[killCol].numsub;
   wNumsub = wP[killCol].numsub;
   // if we start from beginning
   if(start != NULL){if(start->UpathPtr == NULL){f = 1;}}
   // We only have to delete it if there are substitutions there
   if(lNumsub > 0){DeleteColumn(Lose,lP,killCol);}
   //Also we only have to add it if there are substitutions there
   if(wNumsub > 0){CopyUpdateColumn(lP,wP,killCol,Lose);}
   // If we have a stop codon, then if we update carefully we save time
   // we do not have to update everything, just things in this column
   if(flag == 3)
   {
      UpdateStopColumn(lP,wP,killCol,Lose);
      Lose.totNumSub = Win.totNumSub;
   }
   else if((lNumsub > 0)||(wNumsub > 0))
   {  ////////////////////////////////////////////////////////////////////
      // Because there is a difference in loser vs winner columns       //
      // This means the amino acid choice is also different.  Thus we   //
      // must update the amino acid "priorAA & amino" for these columns //
      ////////////////////////////////////////////////////////////////////
      Lose.totNumSub = Win.totNumSub;
      int tmpCol;
      // We must start from position 3 of the codon
      if((killCol%3) == 0){tmpCol = killCol+2;       SetCodonPath(lP, tmpCol);}
      else if((killCol%3) == 1){tmpCol = killCol+1;  SetCodonPath(lP, tmpCol);}
      else{SetCodonPath(lP, killCol);}
      // Now we have to copy the winner liklihoods to the loser's
      // Also the sequence energies must be copied exactly as well.
      if(start != NULL)
      {
         if(f == 1){// Here we start from the beginning of the path
            L = Lose.startPath;
            W = Win.startPath;
         }
         else  // Here we start from a node != start of path
         {
           L = lP[start->column].firstColSub;
           W = wP[start->column].firstColSub;
           while(L->time != start->time){
              L = L->DcolPtr; W = W->DcolPtr;
              if(L == NULL){
                 cerr << "Warning in Update Path!!! ";
                 cerr << " Could not find start" << endl;
              }
           }
         }
         Lose.pathHood = Win.pathHood;
         Lose.probLastEvent = Win.probLastEvent;

         if(LASTFLAG)
         {
            while(L != NULL){
               L->genHood = W->genHood;
               L->subHood = W->subHood;
               L->nrgDiff     = W->nrgDiff;
               L->expNRGdiff  = W->expNRGdiff;
               L->solvNRGdiff = W->solvNRGdiff;
               L = L->DpathPtr; W = W->DpathPtr;
            }
         }
         else
         {
            while(L != NULL)
            {
               L->genHood = W->genHood;
               L->subHood = W->subHood;
               L->nrgDiff     = W->nrgDiff;
               L->expNRGdiff  = W->expNRGdiff;
               L->solvNRGdiff = W->solvNRGdiff;
               L = L->DpathPtr; W = W->DpathPtr;
               if(L->time > maxTime){while(L!=NULL){L = L->DcolPtr;}}
            }
         }
      }
   }
} // End updating the Path structure
//******************************************************************************
void NewPath(Bayes *pE, Information *C, Information *P, Interaction *FO,
        long *seed, int &totP, int &yesP,int &same_P, ofstream& o)
{
   int NewHistory(Information &P, Path *cP, Path *pP,long *seed,
     double &metro_top,double &metro_bot,int& stop,ofstream &out);
   double ProbPathWTimes(Bayes *pE, Interaction *FO,Information &Info,
     Substitution *start, int *tmpAAseq,int *tmpNUCseq,int *siteNeigh,double gH,
     double sH, int count, int which, int LAST, double maxtime);
   void UpdatePath(Information &Win, Path *wP, Information &Lose, Path *lP,
      Substitution *start, int killCol, int flag, int LASTFLAG, double maxTime,
      ofstream &o);
   int NeighborAccess(int AA, int *c, int **nAccess);
   double minof2(double x, double y);
   double rnd(long *seed);

   Substitution *c_tmp, *p_tmp;
   int stop, killCol, same,InflateCount, numBranch;
   double R, metro_top, metro_bot, random, lnProb_top, lnProb_bot;
   double minTargetTime, maxTargetTime, tmpGenHood, tmpSiteHood;
   Substitution *start;  // We must NULL the start pointer!

   numBranch = pE->numBranch;
   int NucLen, AAseqLen, LASTFLAG;
   NucLen = C[0].len; AAseqLen = C[0].AAlen;

   static int z = 0;
   static int *tmpNUCseq;
   static int *tmpAAseq;
   static int *siteNeigh;
   if(z == 0){
      tmpNUCseq = new int[NucLen];    // Temporary nuc sequence
      tmpAAseq  = new int[AAseqLen];  // Temporary AA sequence
      siteNeigh = new int[AAseqLen];
      z = 1;
   }
   int *t1AA, *t2AA, *t1NUC, *t2NUC, *t1NEIGH, *t2NEIGH;
   Path *cP1, *pP1;
for(int i = 0; i < numBranch; i++)
{
   cP1 = C[i].seqPath;  pP1 = P[i].seqPath;
   for(int pathUpdate = 0; pathUpdate < DELTA_PATH; pathUpdate++)
   {
      LASTFLAG=stop=same=InflateCount=0; tmpGenHood = 0.0; tmpSiteHood = 1.0;
      maxTargetTime = 2.0;
      totP++; start = NULL;
      killCol = NewHistory(P[i],cP1,pP1,seed,metro_top,metro_bot,stop,o);
      if(stop == 1){ // We have got a stop codon, and we must deal with it!
         R = 0.0; // We must go back to the current path because of stop
      }

      else if((cP1[killCol].numsub == 0) && (pP1[killCol].numsub == 0))
      { // We do this because the current and proposed are identical
         R = 0.0; // Thus if we don't, then we will be recalculating
         same = 1; same_P++;
      }  // the same thing twice, so to save time we avoid it all

      else{ // if the paths are different, we have to test them both
          // Paths are similar up to a point, thus we do not recalculate
          // Since we just calculated this, we do not need to do it again
         lnProb_bot = C[i].pathHood; // This allows us to not recalc
         //Find out when the first sub in current and proposed occurs
         c_tmp = cP1[killCol].firstColSub;
         p_tmp = pP1[killCol].firstColSub;
         //If there are no substitutions in the Current col
         if(c_tmp == NULL)
         {
            minTargetTime = p_tmp->time;
            int p_count = pP1[killCol].numsub-1;
            for(int qq = 0; qq < p_count; qq++){p_tmp = p_tmp->DcolPtr;}
            Substitution *maxSub = p_tmp->DpathPtr;
            if(maxSub == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
            else
            {
               if(maxSub->DpathPtr == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
               else{maxTargetTime = maxSub->time;}
            }
         }
         //If there are no substitutions in the Proposed col
         else
         {
            if(p_tmp == NULL)
            {
               minTargetTime = c_tmp->time;
               int c_count = cP1[killCol].numsub-1;
               for(int qq = 0; qq < c_count; qq++){c_tmp = c_tmp->DcolPtr;}
               Substitution *maxSub = c_tmp->DpathPtr;
               if(maxSub == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
               else
               {
                  if(maxSub->DpathPtr == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
                  else{maxTargetTime = maxSub->time;}
               }
            }
         //Both will never be true simultaneously, else we get the first case
         // minTargetTime is the first place where the paths differ
            else
            {
               Substitution *maxSub;
               minTargetTime = minof2(c_tmp->time, p_tmp->time);
               int p_count = pP1[killCol].numsub-1;
               for(int qq = 0; qq < p_count; qq++){p_tmp = p_tmp->DcolPtr;}
               int c_count = cP1[killCol].numsub-1;
               for(int qq = 0; qq < c_count; qq++){c_tmp = c_tmp->DcolPtr;}
               if(c_tmp->time > p_tmp->time)
               {
                  maxSub = c_tmp->DpathPtr;
                  if(maxSub == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
                  else
                  {
                     if(maxSub->DpathPtr == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
                     else{maxTargetTime = maxSub->time;}
                  }
               }
               else
               {
                  maxSub = p_tmp->DpathPtr;
                  if(maxSub == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
                  else
                  {
                     if(maxSub->DpathPtr == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
                     else{maxTargetTime = maxSub->time;}
                  }
               }
            }
         }
///////////////////////////////////////////////////////////////////
// Up to Min target time, the two paths are identical, thus we   //
// do not want to waste valuable time recalculating information. //
// This is the basis of our efficiency algorithm here.           //
///////////////////////////////////////////////////////////////////
         t1AA = tmpAAseq; t2AA = C[i].parentAAseq;
         t1NUC = tmpNUCseq; t2NUC = C[i].parentSeq;
         t1NEIGH = siteNeigh; t2NEIGH = C[i].neighAcc;

         for(int col = 0; col < NucLen; col++)
         {// To protect the sequence, we copy it into temporary and work on that
            if(col < AAseqLen){
               *t1AA = *t2AA; *t1NEIGH = *t2NEIGH;
               t1AA++; t2AA++; t1NEIGH++; t2NEIGH++;
            }
            *t1NUC = *t2NUC; t1NUC++; t2NUC++;
         }
         p_tmp  = P[i].startPath;// We now set tmp to the start of the path
         c_tmp  = C[i].startPath;
         int switchFlag = 0;
         if(p_tmp == NULL){p_tmp = C[i].startPath; switchFlag = 1;}
         int minFlag = 0;
         while(p_tmp->time < minTargetTime){
           minFlag = 1;
           int c[3];
           int NUCcol = p_tmp->column;
           int AAcol = NUCcol / 3;
           tmpNUCseq[NUCcol] = p_tmp->pathNuc; // updates the nuc sequence
           tmpAAseq[AAcol] = p_tmp->amino; // updates the AA sequence
           tmpGenHood += p_tmp->genHood;
           tmpSiteHood *= p_tmp->subHood;
           if(tmpSiteHood > BIGNUM){
              tmpSiteHood *= INVBIGNUM;
              InflateCount++;
           }
           else if(tmpSiteHood < INVBIGNUM){
              tmpSiteHood *= BIGNUM;
              InflateCount--;
           }

           switch(NUCcol % 3){
           case 0: c[0] = tmpNUCseq[NUCcol]; c[1] = tmpNUCseq[NUCcol+1];
                   c[2] = tmpNUCseq[NUCcol+2]; break;
           case 1: c[0] = tmpNUCseq[NUCcol-1]; c[1] = tmpNUCseq[NUCcol];
                   c[2] = tmpNUCseq[NUCcol+1]; break;
           case 2: c[0] = tmpNUCseq[NUCcol-2]; c[1] = tmpNUCseq[NUCcol-1];
                   c[2] = tmpNUCseq[NUCcol];   break;
           }
           siteNeigh[AAcol]=NeighborAccess(p_tmp->amino,c,FO->nAccess);
	   p_tmp = p_tmp->DpathPtr;
           if(p_tmp == NULL){
              minTargetTime = tmpGenHood = 0.0;
              p_tmp = P[i].startPath;
              tmpSiteHood = 1.0;  InflateCount = 0;
           }
         }
         p_tmp = p_tmp->UpathPtr;
         if(p_tmp == NULL)
         {
           start = P[i].startPath;
           if(minFlag)
           {
              t1AA = tmpAAseq; t2AA = C[i].parentAAseq;
              t1NUC = tmpNUCseq; t2NUC = C[i].parentSeq;
              t1NEIGH = siteNeigh; t2NEIGH = C[i].neighAcc;

              for(int col = 0; col < NucLen; col++)
              {
                 if(col < AAseqLen){
                   *t1AA = *t2AA; *t1NEIGH = *t2NEIGH;
                    t1AA++; t2AA++; t1NEIGH++; t2NEIGH++;
                 }
                 *t1NUC = *t2NUC; t1NUC++; t2NUC++;
              }
           }
         }
         else if(p_tmp->UpathPtr == NULL)
         {
            t1AA = tmpAAseq; t2AA = C[i].parentAAseq;
            t1NUC = tmpNUCseq; t2NUC = C[i].parentSeq;
            t1NEIGH = siteNeigh; t2NEIGH = C[i].neighAcc;
            for(int col = 0; col < NucLen; col++)
            {
               if(col < AAseqLen)
               {
                 *t1AA = *t2AA; *t1NEIGH = *t2NEIGH;
                 t1AA++; t2AA++; t1NEIGH++; t2NEIGH++;
               }
               *t1NUC = *t2NUC; t1NUC++; t2NUC++;
             }
             minTargetTime = tmpGenHood = 0.0;
             tmpSiteHood = 1.0;  InflateCount = 0;
             start = P[i].startPath;
         }
         else
         {
           start = p_tmp;
           if(switchFlag){start = P[i].startPath;}
         }
         lnProb_top = ProbPathWTimes(pE,FO,P[i],start,tmpAAseq,tmpNUCseq,siteNeigh,
		tmpGenHood,tmpSiteHood,InflateCount,i,LASTFLAG,maxTargetTime);
         // R is the acceptance probability of the path
         R = exp(lnProb_top - lnProb_bot + metro_top - metro_bot);
      } // End path update if there is no stop codon or if killCol is full
      if(!same)
      {
// CHANGE by JEFF         R = minof2(1.0,R);
         random = rnd(seed);
         if(random < R)   // The proposed value is chosen as the winner
         {
            yesP++;
            //int flag = 1;
            UpdatePath(P[i], P[i].seqPath, C[i], C[i].seqPath, start, killCol, 1, LASTFLAG,
                maxTargetTime, o);
         }
         else
         {  // The current value is chosen as the winner in this case
            int flag;
            if(stop){flag = 3;} // cerr << "STOP CODON " << endl;}
            // Houston, we have got problems!
            else {flag = 0;}  // All is clear on the Western Front
            // Above we set start to Proposed->startPath, but because the Current
            // path wins here, we must put start equal to the Current start now.
            if(p_tmp == NULL){start = C[i].startPath;}
            UpdatePath(C[i], C[i].seqPath, P[i], P[i].seqPath, start, killCol, flag, LASTFLAG,
                maxTargetTime, o);
         }
      } // end same loop
   } // End looping through all of the path updates
}
//   delete[] tmpNUCseq; delete[] tmpAAseq; delete[] siteNeigh;
} // End NewPath Subroutine
//******************************************************************************
void ChooseNewNuc(Information *C, int killCol, int *AA, int *nucChoice,
                long *seed, int &stop)
{
   double rnd(long* seed);
   int numNodes, numSeq, numInternal, rootChoice, seqlen, tempNuc;
   int nodesMinus1, nodesMinus2, codon[3];


   double tempProb[4], totalProb, targetProb;
//+ Sang Chul : Jiaye's Revised Pupko
   double sum;
//- Sang Chul : Jiaye's Revised Pupko

   seqlen   = C[0].len;
   numNodes = C[0].numNodes;
   numSeq   = C[0].numSeq;
   numInternal = numNodes - numSeq;
   nodesMinus1 = numNodes - 1;
   nodesMinus2 = numNodes - 2;

   static int z = 0;
   static double **probHolder;
   static int **tempChoice;
   if(z == 0)
   {
     probHolder = new double *[numNodes];
     tempChoice = new int*[numNodes];
     for(int i = 0; i < numNodes; i++)
     {
        tempChoice[i] = new int[4];
        probHolder[i] = new double[4];
     }
     z = 1;
   }

//   for(int i = numSeq; i < numNodes; i++){currentChoice[i] = C[i].seq[killCol];}
   for(int i = 0; i < numNodes; i++)
   {
      int Ci_numChild = C[i].numChild;
      if(i < numSeq)// it is a tip node
      {
         tempNuc = C[i].seq[killCol];
         for(int j = 0; j < 4; j++){probHolder[i][j] = F84_PROB[j][tempNuc];}
      }
      else
      {
         if(i < nodesMinus1) //if not root node
         {
            for(int j = 0; j < 4; j++)
            {
               totalProb = 0.0;
               for(int k = 0; k < 4; k++)
               {
                  tempProb[k] = F84_PROB[j][k];
                  //+ Sang Chul : Jiaye's Revised Pupko
                  for(int l = 0; l < Ci_numChild; l++)
                  {
                     sum = 0;
                     for(int m = 0; m < 4; m++) {
                        sum += probHolder[C[i].child[l]][m];
                     }
                     tempProb[k] *= sum;
                  }
                  //- Sang Chul : Jiaye's Revised Pupko
                  //+ Sang Chul : Old version of Pupko
                  //for(int l = 0; l < Ci_numChild; l++)
                  //{
                  //   tempProb[k] *= probHolder[C[i].child[l]][k];
                  //}
                  //- Sang Chul : Old version of Pupko
                  totalProb +=tempProb[k];
               }
               targetProb = rnd(seed)*totalProb;
               if(targetProb < tempProb[0])
               {
                  probHolder[i][j] = tempProb[0];
                  tempChoice[i][j] = 0;
               }
               else if(targetProb < (tempProb[0]+tempProb[1]))
               {
                  probHolder[i][j] = tempProb[1];
                  tempChoice[i][j] = 1;
               }
               else if(targetProb < (tempProb[0]+tempProb[1]+tempProb[2]))
               {
                  probHolder[i][j] = tempProb[2];
                  tempChoice[i][j] = 2;
               }
               else
               {
                  probHolder[i][j] = tempProb[3];
                  tempChoice[i][j] = 3;
               }
            }
         }
         else // We are now working with the root node
         {
            totalProb = 0.0;
            for(int j = 0; j < 4; j++)
            {
               tempProb[j] = F84_freq[j];

               //+ Sang Chul : Jiaye's Revised Pupko
               for(int k = 0; k < Ci_numChild; k++)
               {
                  sum = 0;
                  for(int m = 0; m < 4; m++) {
                     sum += probHolder[C[i].child[k]][m];
                  }
                  tempProb[j] *= sum;
               }
               //- Sang Chul : Jiaye's Revised Pupko
               //+ Sang Chul : Old version of Pupko
               //for(int k = 0; k < Ci_numChild; k++)
               //{
               //   tempProb[j] *= probHolder[C[i].child[k]][j];
               //}
               //- Sang Chul : Old version of Pupko
               totalProb +=tempProb[j];
            }
            targetProb = rnd(seed)*totalProb;
            if(targetProb < tempProb[0]){rootChoice = 0;}
            else if(targetProb < (tempProb[0]+tempProb[1])){rootChoice = 1;}
            else if(targetProb < (tempProb[0]+tempProb[1]+tempProb[2]))
            {
               rootChoice = 2;
            }
            else{rootChoice = 3;}
         }
      }
   }
   // We now look at the second part of the Pupko Algorithm
   nucChoice[numNodes-1] = rootChoice;
   for(int j = nodesMinus2; j >= numSeq; j--)
   {
      nucChoice[j] = tempChoice[j][nucChoice[C[j].parent]];
   }

   int killMod3, seqCount; killMod3 = killCol%3; seqCount = 0;
   stop = 0;
//   cout << "This is numInternal " << numInternal << " " << killMod3 << endl;
   while(seqCount < numInternal)
   {
  //    cout << "This is seqCount " << seqCount << endl;
      int s_seq = seqCount+numSeq;
      if(!stop)
      {
         int *C_seq = C[s_seq].seq;
         if(!killMod3)
         {
            codon[0] = nucChoice[s_seq];
      //      cout << codon[0] << " " << rootChoice << " Stop codon? " << endl;
            codon[1] = C_seq[killCol+1];
            codon[2] = C_seq[killCol+2];
            AA[s_seq] = Nuc2AATable[codon[0]][codon[1]][codon[2]];
          //  cout << "This is the AA " <<  AA[s_seq] << endl;
            if(AA[s_seq] == 20)
            {
               stop = 1;seqCount = numInternal;
             //  cout <<  " Hi 0 " << codon[0] << " " << codon[1] << " " << codon[2] << endl;
            }
         }
         else if(killMod3 == 1)
         {
            codon[0] = C_seq[killCol-1];
            codon[1] = nucChoice[s_seq];
            codon[2] = C_seq[killCol+1];
            AA[s_seq] = Nuc2AATable[codon[0]][codon[1]][codon[2]];
            if(AA[s_seq] == 20)
            {
               stop = 1;seqCount = numInternal;
            //   cout << " Hi 1 " <<codon[0] << " " << codon[1] << " " << codon[2] << endl;
            }
         }
         else
         {
            codon[0] = C_seq[killCol-2];
            codon[1] = C_seq[killCol-1];
            codon[2] = nucChoice[s_seq];
            AA[s_seq] = Nuc2AATable[codon[0]][codon[1]][codon[2]];
            if(AA[s_seq] == 20)
            {
               stop = 1;seqCount = numInternal;
             //  cout << " Hi 2 " <<codon[0] << " " << codon[1] << " " << codon[2] << endl;
            }
         }
      } // End !stop
      seqCount++;
   } //END while loop
}    // End ChooseNewNuc function
//******************************************************************************
/*
int ChooseNewNuc(Information *C, int killCol, int &AA, long *seed, int numBranch,
                        int &stop)
{  // With a randomly chosen column, this function calculates a new nucleotide
   //  based on the tip sequences, also the measure of how well the new nuc
   //  compares to the old, as well as a check to see if the new nucleotide is
   //  forms a premature stop codon
   double rnd(long* seed);
   double P_A, P_C, P_G, P_T, total, probs;
   int codon[3], tempNuc;

   P_A = F84_freq[0];   P_C = F84_freq[1];
   P_G = F84_freq[2];   P_T = F84_freq[3];

   for(int j = 0; j < numBranch; j++)
   {
      int tempHolder = C[j].seq[1][killCol];
      P_A *= F84_PROB[0][tempHolder];
      P_C *= F84_PROB[1][tempHolder];
      P_G *= F84_PROB[2][tempHolder];
      P_T *= F84_PROB[3][tempHolder];
   }

   total = P_A + P_C + P_G + P_T;
   P_A /= total;
   P_C /= total;
   P_G /= total;
   P_T /= total;
   probs = rnd(seed);

   if(probs < P_A){tempNuc = 0;}
   else if(probs < (P_A + P_C)){tempNuc= 1;}
   else if(probs < (P_A + P_C + P_G)){tempNuc = 2;}
   else{tempNuc = 3;}

   stop = 0;  // We initially think that it will be zero
   int killMod3 = killCol%3;
   int *tempSEQ = C[0].seq[0];
   if(!killMod3)
   {
      codon[0] = tempNuc;
      codon[1] = tempSEQ [killCol+1];
      codon[2] = tempSEQ [killCol+2];
      AA = Nuc2AATable[codon[0]] [codon[1]] [codon[2]];
      if(AA == 20){stop = 1;}
   }
   else if(killMod3 == 1)
   {
      codon[0] = tempSEQ [killCol-1];
      codon[1] = tempNuc;
      codon[2] = tempSEQ [killCol+1];
      AA = Nuc2AATable[codon[0]] [codon[1]] [codon[2]];
      if(AA == 20){stop = 1;}
   }
   else
   {
      codon[0] = tempSEQ [killCol-2];
      codon[1] = tempSEQ [killCol-1];
      codon[2] = tempNuc;
      AA = Nuc2AATable[codon[0]] [codon[1]] [codon[2]];
      if(AA == 20){stop = 1;}
   }
   return tempNuc;
} // End ChooseNewNuc function     */
//******************************************************************************
void NewNodeSitePath(Information &C, Information &P, Path *cP, Path *pP,
        long *seed, double &metro_top, double &metro_bot, int& stop,int killCol,
        ofstream &out)
{    // This function proposes a new site path for each path in the topology
     double rnd(long*);
     void CalcMuts(int col, int& ng, int& nw, Information &I, long* seed,
                ofstream& out);
     void CalcTimes(Path *columnPath, long* seed, int colnum, int AAlen,
                      int ng, int *type, int *force);
     void ChoosingNuc(Information &Info, Path *columnPath, int col,
                long* seed, int *type, int *force);
     void PruneColumn(Information &Info, Path *columnPath, int col);
     void InsertFirstColumn(Information &Info, Path *columnPath, int col);
     void InsertColumn(Information &Info, Path *columnPath, int col);
     void DeleteColumn(Information &Info, Path *columnPath, int col);
     void CopyUpdateColumn(Path *emptyP,Path *fullP,int col,Information &Empty);
     void UpdateStopColumn(Path *emptyP,Path *fullP,int col,Information &Empty);
     double Metro_Path(Path *columnPath, int col);
     void StopCheck(Path *P, int col, int &stop);
     void SetCodonPath(Path *P, int col);

     // Fast method for choosing the correct column to delete

     static int z = 0;
     static int *type;
     static int *force;
     if(z == 0)
     {
        type  = new int[25];
        force = new int[25];
        z = 1;
     }

     if(pP[killCol].numsub > 0)
     {
        // If this column has substitutions, subtract that number from the total
        P.totNumSub -= pP[killCol].numsub;
        // delete this column in P path, making it clear for a new column
        DeleteColumn(P,pP,killCol);
     }
     int ng, nw;       ng = nw = 0;    // samples ng and nw from dist.
     CalcMuts(killCol,ng,nw,P,seed,out);
     // if ng + nw = 0, then we would be wasting our times here
     int zeroflag = 0;
     if(P.totNumSub == 0){zeroflag = 1;}
     if((ng+nw) > 0)
     {
        int subs = ng + nw;
        pP[killCol].numsub = subs;  // set numsub at top of col path
        P.totNumSub += subs; // add to total
        CalcTimes(pP,seed,killCol,P.AAlen,ng,type,force);
        ChoosingNuc(P,pP,killCol,seed,type,force);
        PruneColumn(P,pP,killCol);
        if(zeroflag == 1){InsertFirstColumn(P,pP,killCol);}
        else{InsertColumn(P,pP, killCol);}
     }

     // Now we must check for the creation of a stop codon in the mix
     int tmpCol; // if so, stop will equal 1, and we have to deal with it.
     // We initially set stop = 0 before this function call
     int killMod3 = killCol%3;
     if(!killMod3)
     {
        tmpCol = killCol+2; // We must start from position 3 of the codon
        StopCheck(pP, tmpCol, stop);
        if(!stop){SetCodonPath(pP, tmpCol);}
     }
     else if(killMod3 == 1)
     {
        tmpCol = killCol+1; // We must start from position 3 of the codon
        StopCheck(pP, tmpCol, stop);
        if(!stop){SetCodonPath(pP, tmpCol);}
     }
     else
     {
        StopCheck(pP, killCol, stop);
        if(!stop){SetCodonPath(pP, killCol);}
     }

     // This is the probability of the P column given the original.
     // This takes care of no subs in the original col path
     if(!stop)
     {
       if(cP[killCol].numsub == 0)
       {  // ln ( e^(-WT))     or substitutions to oneself
          metro_top = F84_RATE[C.parentSeq[killCol]][C.seq[killCol]];
       } // This takes care of no subs in the new col path
       else{metro_top = Metro_Path(cP, killCol);}

       // If proposed path has a diff # of subs, then we must recalc Rate Matrix
       if(pP[killCol].numsub == 0)
       {  // ln ( e^(-WT))      or substitutions to oneself
          metro_bot = F84_RATE[P.parentSeq[killCol]][P.seq[killCol]];
       }
       else{metro_bot = Metro_Path(pP, killCol);}
     } // End if (stop = 0) statement
     else{metro_bot = metro_top = 0.0;}// This is case for existence of stop
}  // End NewNodeSitePath Subroutine
//******************************************************************************
void UpdateNodeSeq(Bayes *pE, Information *Win, Information *Lose,
        Interaction *FO, int killCol, int winFlag)
{  // After a winner is chosen, we must replace the loser column with the winner
   void DeleteColumn(Information &Info, Path *columnPath, int col);
   void CopyUpdateColumn(Path *emptyP,Path *fullP,int col,Information &Empty);
   void UpdateStopColumn(Path *emptyP,Path *fullP,int col,Information &Empty);
   void SetCodonPath(Path *P, int col);

//***********FROM HERE***************
   int numBranch = pE->numBranch;
   int initNuc, newNuc, killAAcol, AAseqLen,numSeq;
   numSeq = pE->numSeq;
   killAAcol = killCol / 3;
   AAseqLen = Win[0].AAlen;
   if(winFlag)
   {
      initNuc = Lose[numBranch].seq[killCol];
      newNuc  = Win[numBranch].seq[killCol];
   }
   else
   {
      initNuc = Win[numBranch].seq[killCol];
      newNuc  = Lose[numBranch].seq[killCol];
   }

   // we have to update all of the energies
   // as well in order for the program to work !!!
   for(int i = 0; i < numBranch; i++)
   {
      Lose[i].seqPath[killCol].nuc = Win[Lose[i].parent].seq[killCol];
      Lose[i].neighAcc[killAAcol] = Win[i].neighAcc[killAAcol];
      if(i >= numSeq)
      {
         Lose[i].seq[killCol]      =  Win[i].seq[killCol];
         Lose[i].AA_seq[killAAcol] =  Win[i].AA_seq[killAAcol];
      }
   }

   Lose[numBranch].seq[killCol]      =  Win[numBranch].seq[killCol];
   Lose[numBranch].AA_seq[killAAcol] =  Win[numBranch].AA_seq[killAAcol];


      // If Propsed wins: This is what we need to do with the Bayes
   if(winFlag )  // winFlag == 1 means Proposed won!
   {
      pE->solvNRG = pE->tempsolvNRG;
      pE->pairNRG = pE->temppairNRG;
      switch(initNuc){
      case 0: pE->number_A--; break;
      case 1: pE->number_C--; break;
      case 2: pE->number_G--; break;
      case 3: break;
      }

      switch(newNuc){
      case 0: pE->number_A++; break;
      case 1: pE->number_C++; break;
      case 2: pE->number_G++; break;
      case 3: break;
      }
   }
//***********TO HERE****************************
// only needs to  be done if newNuc != initNuc *
//**********************************************
//******************************************************************************
// This is what we need to do if Either wins regardless
//******************************************************************************

   Substitution *L, *W;
   Path *WinPath, *LosePath;
   for(int x = 0; x < numBranch; x++)
   {
      WinPath = Win[x].seqPath; LosePath = Lose[x].seqPath;
      // We only have to delete it if there are substitutions there
      if(LosePath[killCol].numsub > 0){DeleteColumn(Lose[x],LosePath,killCol);}

      //Also we only have to add it if there are substitutions there
      if(WinPath[killCol].numsub > 0)
      {
         CopyUpdateColumn(LosePath,WinPath,killCol,Lose[x]);
      }

      Lose[x].totNumSub = Win[x].totNumSub;
      Lose[x].pathHood  = Win[x].pathHood;
      Lose[x].probLastEvent = Win[x].probLastEvent;

      int tmpCol;

      // We must start from position 3 of the codon
      if((killCol%3) == 0){tmpCol = killCol+2;     SetCodonPath(LosePath,tmpCol);}
      else if((killCol%3) == 1){tmpCol = killCol+1;SetCodonPath(LosePath,tmpCol);}
      else{SetCodonPath(LosePath, killCol);}

      L = Lose[x].startPath;
      W = Win[x].startPath;

      while(W != NULL){
         L->genHood = W->genHood;
         L->subHood = W->subHood;
         L->nrgDiff     = W->nrgDiff;
         L->expNRGdiff  = W->expNRGdiff;
         L->solvNRGdiff = W->solvNRGdiff;
         L = L->DpathPtr; W = W->DpathPtr;
      }
   }
}
//******************************************************************************
// This ends what we need to do if either wins regardless
//******************************************************************************
//******************************************************************************
void NewNodeSeq(Bayes *pE, Information *C, Information *P, Interaction *FO,
        long *seed, int &totP, int &yesP,int &same_P, ofstream& o)
{
   void ChooseNewNuc(Information *C, int killCol, int *AA, int *nucChoice,
                long *seed, int &stop);
   double ProbPathWTimes(Bayes *pE, Interaction *FO,Information &Info,
     Substitution *start, int *tmpAAseq,int *tmpNUCseq,int *siteNeigh,double gH,
     double sH, int count, int which, int LAST, double maxtime);
   void SeqEnergy(Interaction *FO, int *tmpAAseq, double &solv, double &pair,
                        int AAlen, ofstream &out);
   void NewNodeSitePath(Information &C, Information &P, Path *cP, Path *pP,
        long *seed, double &metro_top, double &metro_bot, int &stop,int killCol,
        ofstream &o);
   void UpdateNodeSeq(Bayes *pE, Information *C, Information *P,
        Interaction *FO, int killCol, int winFlag);
   void UpdatePath(Information &Win, Path *wP, Information &Lose, Path *lP,
      Substitution *start, int killCol, int flag, int LASTFLAG, double maxTime,
      ofstream &o);
   double CalcConst(Bayes *pE, Information &P, Interaction *FO, double &maximum,
                long *seed, int round, ofstream &o);
   void DeleteColumn(Information &Info, Path *columnPath, int col);
   void CopyUpdateColumn(Path *emptyP,Path *fullP,int col,Information &Empty);
   void UpdateStopColumn(Path *emptyP,Path *fullP,int col,Information &Empty);
   int NeighborAccess(int AA, int *c, int **nAccess);
   double minof2(double x, double y);
   double rnd(long *seed);

   Path *cP, *pP;
   Substitution *c_tmp, *p_tmp;
   int stop, killCol, killAAcol, same, InflateCount, numBranch, sameNucFlag;
//CHANGE by JEFF   int newAA, NucLen, AAseqLen, stopFlag, Heather, initNuc, newNeigh, codon[3];
   int NucLen, AAseqLen, stopFlag, Heather, initNuc, codon[3];
   int whichBranchStop, LASTFLAG,numSeq,numNodes;

   double R, metro_top, metro_bot, random, lnProb_top, lnProb_bot;
   double minTargetTime, maxTargetTime, tmpGenHood, tmpSiteHood;
   Substitution *start;  // We must NULL the start pointer!
   numSeq = pE->numSeq;
   numBranch = pE->numBranch;
   numNodes = C[0].numNodes;
   metro_top = metro_bot = 1.0;
   NucLen = C[0].len; AAseqLen = C[0].AAlen;
   totP++;
   static int z = 0;
   static int *tmpNUCseq,*tmpAAseq,*siteNeigh;
   static int *nucChoice,*oldNucs,*AAchoice,*oldAAchoice;

   if(z == 0){
      tmpNUCseq = new int[NucLen];    // Temporary nuc sequence
      tmpAAseq  = new int[AAseqLen];  // Temporary AA sequence
      siteNeigh = new int[AAseqLen];
      nucChoice = new int[numNodes];
      oldNucs   = new int[numNodes];
      AAchoice  = new int[numNodes];
      oldAAchoice  = new int[numNodes];
      z = 1;
   }

   whichBranchStop = Heather = stop = stopFlag = sameNucFlag = 0;
   // Here we choose the column in which we are going to kill
   killCol = (int) (NucLen * rnd(seed));
//   o << "This is the killCol " << killCol << endl;
   nucChoice[numBranch]=20;
   killAAcol = killCol / 3;
   initNuc = C[numBranch].seq[killCol];
   AAchoice[numBranch] = 20;
   ChooseNewNuc(C, killCol, AAchoice, nucChoice, seed, stop);
   int newNuc = nucChoice[numBranch];
//   o << " newNuc= " << newNuc << " at killcol = " << killCol << " AA choice " << AAchoice[numBranch] << endl;
//   o << "This is aachoice " << AAchoice[numBranch] << endl;
   pE->tempsolvNRG = pE->solvNRG;
   pE->temppairNRG = pE->pairNRG;

   if(stop){//Then we do not have to update anything here because
        // everything is still the same.
        stopFlag = Heather = 1; //o << "Heather = 1: WOW!" << endl;
   }
   else
   {
      for(int a = 0; a < numBranch; a++)
      {
         P[a].seqPath[killCol].nuc = nucChoice[P[a].parent];
         if(a >= numSeq)
         {
            oldNucs[a] = C[a].seq[killCol];
            P[a].seq[killCol] = nucChoice[a];
            P[a].AA_seq[killAAcol] = AAchoice[a];
         }
      }
      oldNucs[numBranch] = C[numBranch].seq[killCol];
      P[numBranch].seq[killCol] = nucChoice[numBranch];
      P[numBranch].AA_seq[killAAcol] = AAchoice[numBranch];

      int *t1AA, *t2AA, *t1NUC, *t2NUC, *t1NEIGH, *t2NEIGH;
// CHANGE by JEFF      R = 1.0;
      double lnR = 0.0;
// CHANGE by JEFF      int energyFlag = 0;
      LASTFLAG = 1;
      for(int i = 0; i < numBranch; i++)//must create a new path for each branch
      {
        //o << "This is start of update for branch " << i << endl;
        cP = C[i].seqPath;     pP = P[i].seqPath;
        same = InflateCount = 0; tmpGenHood = 0.0; tmpSiteHood = 1.0;
        start = NULL; maxTargetTime = 2.0;
        t1AA = tmpAAseq; t2AA = P[i].parentAAseq;
        t1NUC = tmpNUCseq; t2NUC = P[i].parentSeq;
        t1NEIGH = siteNeigh; t2NEIGH = P[i].neighAcc;

        for(int col = 0; col < NucLen; col++)
        {// To protect the sequence, we copy it into temporary and work on that
           if(col < AAseqLen)
           {
              *t1AA = *t2AA; *t1NEIGH = *t2NEIGH;
              t1AA++; t2AA++; t1NEIGH++; t2NEIGH++;
           }
           *t1NUC = *t2NUC; t1NUC++; t2NUC++;
        }
        if(nucChoice[P[i].parent] != oldNucs[P[i].parent]){sameNucFlag = 0;}
        else{sameNucFlag = 1;}

        if(!sameNucFlag)
        {
           LASTFLAG = 1;
           if((i == (numBranch-1))&&(AAchoice[numBranch]!=C[numBranch].AA_seq[killAAcol]))
           {
              SeqEnergy(FO,tmpAAseq,pE->tempsolvNRG,pE->temppairNRG,AAseqLen,o);
           }

           switch(killCol % 3){
           case 0: codon[0] = tmpNUCseq[killCol];
                   codon[1] = tmpNUCseq[killCol+1];
                   codon[2] = tmpNUCseq[killCol+2]; break;
           case 1: codon[0] = tmpNUCseq[killCol-1];
                   codon[1] = tmpNUCseq[killCol];
                   codon[2] = tmpNUCseq[killCol+1]; break;
           case 2: codon[0] = tmpNUCseq[killCol-2];
                   codon[1] = tmpNUCseq[killCol-1];
                   codon[2] = tmpNUCseq[killCol];   break;
           }
           P[i].neighAcc[killAAcol] = siteNeigh[killAAcol] =
                      NeighborAccess(AAchoice[P[i].parent],codon,FO->nAccess);
           // HERE IS A CHECK ON THE PROGRAM!!!!
//           if(siteNeigh[killAAcol] != P[i].neighAcc[killAAcol])
//           {
//              cout << "Warning!! in new node, site neigh not equal!!!" << endl;
//           }
           //o << "This is the new sightNeigh " << siteNeigh[killAAcol] << endl;
           //o << "First sightNeigh " << siteNeigh[0] << endl;
           lnProb_bot = C[i].pathHood;
          //o << "Nucs not the same : current pathhood = " << lnProb_bot << endl;
        }
        if(!stop)
        {
           //o << "before new node site path " << endl;
           NewNodeSitePath(C[i],P[i],cP,pP,seed,metro_top,metro_bot,stop,killCol,o);
           start = P[i].startPath;
          //o << "metro_top = " << metro_top << " metro_bot = " << metro_bot << endl;
        }

        if(stop){ // We have got a stop codon, and we must deal with it!
           if(!stopFlag){whichBranchStop = i+1;}
           stopFlag = 1; // We must go back to the current path because of stop
           i = numBranch + 1;  // This stops the iteration through the branches
           //o << "stop codon in sequence  " << i-1 << endl;
        }
        else
        {
          if(sameNucFlag)//This can only be used if nucleotides are the same
          {
            LASTFLAG = 0;
            if((cP[killCol].numsub == 0) && (pP[killCol].numsub == 0))
            { // We do this because the current and proposed are identical
              // Thus if we don't, then we will be recalculating
               //o << "branch " << i << " both have 0 subs " << endl;
               same = 1; same_P++;
            }  // the same thing twice, so to save time we avoid it all

            else // if the paths are different, we have to test them both
            {    // Paths are similar up to a point, thus we do not recalculate
                // Since we just calculated this, we do not need to do it again
               //o << "The two branches are not the same " << endl;
               lnProb_bot = C[i].pathHood; // This allows us to not recalc
               //o << "This is prob bot " << lnProb_bot << endl;
               //Find out when the first sub in current and proposed occurs
               c_tmp = cP[killCol].firstColSub;
               p_tmp = pP[killCol].firstColSub;
               //If there are no substitutions in the Current col
               if(i < numSeq)
               {
                 if(c_tmp == NULL)
                 {
                    minTargetTime = p_tmp->time;
                    int p_count = pP[killCol].numsub-1;
                    for(int qq = 0; qq < p_count; qq++){p_tmp = p_tmp->DcolPtr;}
                    Substitution *maxSub = p_tmp->DpathPtr;
                    if(maxSub == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
                    else
                    {
                       if(maxSub->DpathPtr == NULL)
                       {
                          LASTFLAG = 1; maxTargetTime = 2.0;
                       }
                       else{maxTargetTime = maxSub->time;}
                    }
                 }
                 //If there are no substitutions in the Proposed col
                 else
                 {
                    if(p_tmp == NULL)
                    {
                       minTargetTime = c_tmp->time;
                       int c_count = cP[killCol].numsub-1;
                       for(int qq = 0; qq < c_count;qq++){c_tmp=c_tmp->DcolPtr;}
                       Substitution *maxSub = c_tmp->DpathPtr;
                       if(maxSub == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
                       else
                       {
                          if(maxSub->DpathPtr == NULL)
                          {
                             LASTFLAG = 1; maxTargetTime = 2.0;
                          }
                          else{maxTargetTime = maxSub->time;}
                       }
                    }
               //Both will never be true simultaneously, else we get the first case
               // minTargetTime is the first place where the paths differ
                    else
                    {
                       Substitution *maxSub;
                       minTargetTime = minof2(c_tmp->time, p_tmp->time);
                       int p_count = pP[killCol].numsub-1;
                       for(int qq = 0; qq < p_count;qq++){p_tmp=p_tmp->DcolPtr;}
                       int c_count = cP[killCol].numsub-1;
                       for(int qq = 0; qq < c_count;qq++){c_tmp=c_tmp->DcolPtr;}
                       if(c_tmp->time > p_tmp->time)
                       {
                          maxSub = c_tmp->DpathPtr;
                          if(maxSub == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
                          else
                          {
                             if(maxSub->DpathPtr == NULL)
                             {
                                LASTFLAG = 1; maxTargetTime = 2.0;
                             }
                             else{maxTargetTime = maxSub->time;}
                          }
                       }
                       else
                       {
                          maxSub = p_tmp->DpathPtr;
                          if(maxSub == NULL){LASTFLAG = 1; maxTargetTime = 2.0;}
                          else
                          {
                             if(maxSub->DpathPtr == NULL)
                             {
                                LASTFLAG = 1; maxTargetTime = 2.0;
                             }
                             else{maxTargetTime = maxSub->time;}
                          }
                       }
                    }
                 } // End else statement
               } // End if(i < numSeq)
               else  // i >= numSeq
               {
                  LASTFLAG = 1;
                  if(c_tmp == NULL){minTargetTime = p_tmp->time;}
                  else if(p_tmp == NULL){minTargetTime = c_tmp->time;}
                  else{minTargetTime = minof2(c_tmp->time, p_tmp->time);}
               }  // end  i >= numSeq
///////////////////////////////////////////////////////////////////
// Up to Min target time, the two paths are identical, thus we   //
// do not want to waste valuable time recalculating information. //
// This is the basis of our efficiency algorithm here.           //
///////////////////////////////////////////////////////////////////
               // We now set tmp to the start of the path
               p_tmp  = P[i].startPath;
               c_tmp  = C[i].startPath;
               int switchFlag = 0;
               if(p_tmp == NULL){p_tmp = C[i].startPath; switchFlag = 1;}
               int minFlag = 0;
               while(p_tmp->time < minTargetTime)
               {
                   minFlag = 1;
                   int c[3];
                   int NUCcol = p_tmp->column;
                   int AAcol = NUCcol / 3;
                   tmpNUCseq[NUCcol] = p_tmp->pathNuc;//updates the nuc sequence
                   tmpAAseq[AAcol] = p_tmp->amino; // updates the AA sequence
                   tmpGenHood += p_tmp->genHood;
                   tmpSiteHood *= p_tmp->subHood;
                   if(tmpSiteHood > BIGNUM){
                     tmpSiteHood *= INVBIGNUM;
                     InflateCount++;
                   }
                   else if(tmpSiteHood < INVBIGNUM){
                     tmpSiteHood *= BIGNUM;
                     InflateCount--;
                   }

                   switch(NUCcol % 3){
                   case 0: c[0] = tmpNUCseq[NUCcol]; c[1] = tmpNUCseq[NUCcol+1];
                           c[2] = tmpNUCseq[NUCcol+2]; break;
                   case 1: c[0] = tmpNUCseq[NUCcol-1]; c[1] = tmpNUCseq[NUCcol];
                           c[2] = tmpNUCseq[NUCcol+1]; break;
                   case 2: c[0] = tmpNUCseq[NUCcol-2]; c[1]=tmpNUCseq[NUCcol-1];
                           c[2] = tmpNUCseq[NUCcol];   break;
                   }
                   siteNeigh[AAcol]=NeighborAccess(p_tmp->amino,c,FO->nAccess);
        	   p_tmp = p_tmp->DpathPtr;
                   if(p_tmp == NULL){
                      minTargetTime = tmpGenHood = 0.0;
                      p_tmp = P[i].startPath;
                      tmpSiteHood = 1.0;  InflateCount = 0;
                   }
               } // End while loop
               p_tmp = p_tmp->UpathPtr;
               if(p_tmp == NULL)  // If we are at the beginning of the path
               {

                 start = P[i].startPath;  // Start at the beginning
                 if(minFlag)  // if minFlag == 1, then we recopy everything
                 {
                    t1AA = tmpAAseq; t2AA = P[i].parentAAseq;
                    t1NUC = tmpNUCseq; t2NUC = P[i].parentSeq;
                    t1NEIGH = siteNeigh; t2NEIGH = P[i].neighAcc;

                    for(int col = 0; col < NucLen; col++)
                    {
                       if(col < AAseqLen){
                         *t1AA = *t2AA; *t1NEIGH = *t2NEIGH;
                          t1AA++; t2AA++; t1NEIGH++; t2NEIGH++;
                       }
                       *t1NUC = *t2NUC; t1NUC++; t2NUC++;
                    }
                 }
               }
               else if(p_tmp->UpathPtr == NULL)
               {
                  t1AA = tmpAAseq; t2AA = P[i].parentAAseq;
                  t1NUC = tmpNUCseq; t2NUC = P[i].parentSeq;
                  t1NEIGH = siteNeigh; t2NEIGH = P[i].neighAcc;

                  for(int col = 0; col < NucLen; col++)
                  {
                     if(col < AAseqLen){
                       *t1AA = *t2AA; *t1NEIGH = *t2NEIGH;
                        t1AA++; t2AA++; t1NEIGH++; t2NEIGH++;
                     }
                     *t1NUC = *t2NUC; t1NUC++; t2NUC++;
                  }
                  minTargetTime = tmpGenHood = 0.0;
                  tmpSiteHood = 1.0;  InflateCount = 0;
                  start = P[i].startPath;
               }

               else
               {
                 start = p_tmp;   // Start from correct place in path
                 if(switchFlag){start = P[i].startPath;}// start = beginning
               }
            } // end else statement for paths not being identical
          }  // end if(sameNucFlag) statement for one path
        }   // end if(stop) else statement

        if((!stopFlag)&&(!same))
        {
           //o << "We now calc lnProb_top " << endl;
           lnProb_top = ProbPathWTimes(pE,FO,P[i],start,tmpAAseq,tmpNUCseq,
                        siteNeigh,tmpGenHood,tmpSiteHood,InflateCount,i,
                        LASTFLAG, maxTargetTime);
          //o << "LnProb_bot= " <<lnProb_bot<<" lnProb_top "<<lnProb_top<< endl;
           //o << "We now calc lnProb_top  = " << lnProb_top << endl;
           // R is the acceptance probability of the path
           lnR += lnProb_top - lnProb_bot + metro_top - metro_bot;
           //o << "This is R so far " << R << endl;
        } // End path update if there is no stop codon or if killCol is full

        //We have to find an easy way to add the energy of the seqs
        // this gets at the probability of each sequence.
      } // end running through all of the branches in the tree topology
      //o << "this is before R " << endl;
      lnR += (pE->solvent*(pE->tempsolvNRG - pE->solvNRG)) +
              (pE->pairwise*(pE->temppairNRG - pE->pairNRG)) ;
      R = exp(lnR)*
        ((pE->MCMCnf[newNuc]*F84_freq[initNuc])/(pE->MCMCnf[initNuc]*F84_freq[newNuc]));
        //o << "  This is after R " << R << endl;
   } // End huge else statement if(stop) ELSE statement
   if(!stopFlag)
   {
// CHANGE by JEFF      R = minof2(1.0,R);
      random = rnd(seed);
      if(random < R)   // The proposed value is chosen as the winner
      {
  //      o << "Proposed new path won and rnd = " << random << endl;
         yesP++;
         // int winFlag = 1;
         UpdateNodeSeq(pE, P, C, FO, killCol, 1);
      }
      else // int winFlag = 0;
      {  // The current value is chosen as the winner in this case
  //       o << "Current new path won and rnd = " << random << endl;
         //Here instead of newAA, we send in the original AA so it can be copied
         UpdateNodeSeq(pE,C,P,FO,killCol, 0);
      }
   } // end same loop
   else
   {  //if Heather == 1: Then we do nothing because nothing was changed
 //    o << "There was a stop codon" << endl;
      if(!Heather) //if Heather == 0
      { // We have to redo everything because all paths were already made
        // So this means we have a stop codon in one of the paths, and
        // we have to update all of the columns
         //o << "One of the sequences had a stop codon " << whichBranchStop << endl;
         //If there is a stop, then we have got to update a whole lot!!
   //      o << "We are in Not Heather ! " << endl;
         LASTFLAG = 1;
         maxTargetTime = 2.0;
         for(int a = 0; a < numBranch; a++)
         {
            P[a].seqPath[killCol].nuc = C[P[a].parent].seq[killCol];
            P[a].neighAcc[killAAcol] = C[a].neighAcc[killAAcol];
            if(a >= numSeq)
            {
               P[a].seq[killCol]      = C[a].seq[killCol];
               P[a].AA_seq[killAAcol] = C[a].AA_seq[killAAcol];
            }

         }
         P[numBranch].seq[killCol]      = C[numBranch].seq[killCol];
         P[numBranch].AA_seq[killAAcol] = C[numBranch].AA_seq[killAAcol];

         for(int q = 0; q < whichBranchStop; q++)
         {
            cP = C[q].seqPath;  pP = P[q].seqPath;
            if(q == (whichBranchStop - 1))
            {
               UpdatePath(C[q], cP, P[q], pP, C[q].startPath, killCol, 3,
                        LASTFLAG, maxTargetTime, o);
            }
            else
            {
               UpdatePath(C[q], cP, P[q], pP, C[q].startPath, killCol, 0,
                        LASTFLAG, maxTargetTime, o);
            }
         }
      } // end if not heather statement
   } // end else statement
//   delete[] tmpNUCseq; delete[] tmpAAseq; delete[] siteNeigh;
} // End NewNodeSeq Subroutine
//******************************************************************************
void PathOutput(double ***myPath, int *PathSubsCount, ofstream &o)
{

   double **pathPtr, *pathRowPtr;
   int *numCount = PathSubsCount;

   for(int qq = 0; qq < SAMPLE_PATH; qq++)
   {
      pathPtr = *myPath;
      for(int j = 0; j < *numCount; j++)
      {
         pathRowPtr = *pathPtr;
         for(int qqq = 0; qqq < 17; qqq++){o<< *pathRowPtr << " ";pathRowPtr++;}
         o << endl; pathPtr++;
      }
      myPath++; numCount++;
   }
}
//******************************************************************************

void DeleteSome(Bayes *pE, Information *C,Information *P, long* seed,
        ofstream& out)
{
     void DeleteColumn(Information &Info, Path *columnPath, int col);
     double rnd(long *seed);
     int Nuclen = C[0].len;

     int numBranch = pE->numBranch;


   for(int i = 0; i < numBranch; i++)
   {
     Path *cP = C[i].seqPath;
     Path *pP = P[i].seqPath;
     for(int col = 0; col < Nuclen; col++)
     {
        if(cP[col].numsub > 0){DeleteColumn(C[i],cP,col);}
        if(pP[col].numsub > 0){DeleteColumn(P[i],pP,col);}
     }
     C[i].startPath = NULL; P[i].startPath = NULL;
     //The next line is temporary, and will be
     //updated in the maximization routine
//        Current->kappa  = INIT_K;   Proposed->kappa  = INIT_K;
//        Current->omega  = INIT_W;   Proposed->omega  = INIT_W;
//        Current->w_solv = INIT_W_S; Proposed->w_solv = INIT_W_S;
//        Current->w_pair = INIT_W_P; Proposed->w_pair = INIT_W_P;
        C[i].totNumSub = 0;        P[i].totNumSub = 0;
   }

#if SIM_U
   for(int i = 0; i < numBranch; i++)
   {
      pE->rate[i]  = U_MIN + rnd(seed)*(U_MAX - U_MIN);
   }
#endif
#if SIM_K
        if(MULTI_K)
        {
           for(int i = 0; i < pE->numBranch; i++)
           {
              pE->kappa[i]  = K_MIN + rnd(seed)*(K_MAX - K_MIN);
           }
        }
        else
        {
              pE->kappa[0]  = K_MIN + rnd(seed)*(K_MAX - K_MIN);
        }
#endif
#if SIM_W
        if(MULTI_W)
        {
           for(int i = 0; i < pE->numBranch; i++)
           {
              pE->omega[i]  = W_MIN + rnd(seed)*(W_MAX - W_MIN);
           }
        }
        else
        {
              pE->omega[0]  = W_MIN + rnd(seed)*(W_MAX - W_MIN);
        }
#endif
#if SIM_W_S
   pE->solvent  = W_S_MIN + rnd(seed)*(W_S_MAX - W_S_MIN);
#endif
#if SIM_W_P
   pE->pairwise  = W_P_MIN + rnd(seed)*(W_P_MAX - W_P_MIN);
#endif

        out << "Done initializing everything " << endl;
        out << "Init K = " << pE->kappa[0] << ", and W ";
        out << pE->omega[0] << ", w_S = " << pE->solvent;
        out << ", w_p = " << pE->pairwise;
        out << ", and T = " << 1.0 << endl;
        out << ", and U = " << pE->rate[0] << " " << pE->rate[1];
        out << " " << pE->rate[2] ;
        // Put the matrices into their respective structures
//**********************
//***End Rate Matrix ***
//***Begin Prob Matrix**  Given the rate matrix, this calcs the prob matrix
//**********************
        //To calculate the eigenvalues and eigenvectors we need to
        // use PAML and convert our matrix into a usable form
        //For this the matrix class library is the best choice around
        //This routine also calculates the F84 probability matrix
      //  EigenCalc(Current->F84Rate, Current->F84Prob, t);
     // This copies the matrices into the proposed structure

}  // End DeleteSome Subroutine
//******************************************************************************
void DelInteract(Interaction *FirstOrder, int AAlen)
{

   for(int i = 0; i < 20 ; i++){delete[] FirstOrder->nAccess[i];}
   delete[] FirstOrder->nAccess;

   for(int i = 0; i < 61 ; i++){delete[] FirstOrder->neighList[i];}
   delete[] FirstOrder->neighList;

   for(int AAcol = 0; AAcol < AAlen; AAcol++)
   {
      delete[] FirstOrder->AAinfo[AAcol].DDDneighbors;
      delete[] FirstOrder->AAinfo[AAcol].solvent;
      int numNeigh = FirstOrder->AAinfo[AAcol].numNeighbors;
      for(int aa2 = 0; aa2 < numNeigh; aa2++)
      {
         delete[] FirstOrder->AAinfo[AAcol].pairwise[aa2];
      }
      delete[] FirstOrder->AAinfo[AAcol].pairwise;
    }
    delete[] FirstOrder->AAinfo;
}
//******************************************************************************
void DeleteInfo(Information &I)
{
   Substitution *temp;
   Substitution *next;
   temp = I.startPath;
   int totsub, AAseqLen; totsub = I.totNumSub; AAseqLen = I.AAlen;
   for(int i = 0; i < totsub; i++)
   {
      next = temp->DpathPtr;
      delete[] temp;
      if(i < (I.totNumSub - 1)){temp = next;}
   }
   delete[] I.seq;
   delete[] I.neighAcc;
   delete[] I.AA_seq;
}
//******************************************************************************
void DeleteBayes(Bayes *pE)
{
   delete[] pE->kappa;
   delete[] pE->omega;
   delete[] pE->rate;
   delete[] pE->branchLen;
   delete[] pE->MCMCnf;
   delete[] pE->tempMCMCnf;
   delete[] pE;
}
//******************************************************************************
void DeleteAll(Bayes *pE, Interaction *FirstOrder, Information *Current,
               Information *Proposed, ofstream& o)
{
///////////////////////////////////////////////////////////
// At last we are finished.  We borrowed a great deal of //
// memory, and now we must return all of it to the       //
// computer.  I now finish by DELETING IT ALL!!!!        //
// DELETE ALL I COMMAND!!                                //
///////////////////////////////////////////////////////////
   void DelInteract(Interaction *FO, int AAlen);
   void DeleteInfo(Information &I);
   void DeleteBayes(Bayes *pE);
   int count = 0;
   for(int i=0; i < SIZEOFGRID; i++){
      for(int j=0; j < SIZEOFGRID; j++){
         for(int k=0; k < FREQGRIDSIZE; k++){
            for(int l=0; l < FREQGRIDSIZE; l++){
               for(int m=0; m < FREQGRIDSIZE; m++){
                  if(GibbsInfo[i][j][k][l][m] != NULL)
                     {delete[] GibbsInfo[i][j][k][l][m]; count++;}
               }
            }
         }
      }
   }

   o << "We have entered the Gibbs sampler a total of " << count << " times. ";
   DelInteract(FirstOrder,Current[0].AAlen); delete[] FirstOrder;
   for(int i = 0; i < pE->numBranch; i++)
   {
      delete[] Current[i].seqPath;
      delete[] Proposed[i].seqPath;
      DeleteInfo(Current[i]);
      DeleteInfo(Proposed[i]);
   }
   delete[] Current;
   delete[] Proposed;
   DeleteBayes(pE);
}  // End DeleteAll Subroutine
//******************************************************************************
void PrintPath(Substitution *start, ofstream& o)
{  // Simple function that prints out everything about the path to locate errors
        int i; i = 0;
        while(start != NULL)
        {
           o << "Sub (" << i << "), col = " << start->column <<
           ", time = " << start->time <<  ", prior = " << start->priorNuc <<
           ", chosen = " << start->pathNuc << ", up AA = " << start->priorAA <<
           ", d AA = " << start->amino << ", NRG = " << start->nrgDiff <<
           " , exp NRG = " << start->expNRGdiff << endl;
           start = start->DpathPtr;
           i++;
        }
}   // End PrintPath Subroutine
//******************************************************************************
void ErrorHandle(char *statement, Bayes *pE, Interaction *FirstOrder,
          Information *Current, Information *Proposed, ofstream& o)
{
   void DeleteAll(Bayes *pE, Interaction *FirstOrder, Information *Current,
               Information *Proposed, ofstream& o);

   cerr << statement << endl;
   DeleteAll(pE, FirstOrder,Current, Proposed, o);
   void abort(void);
}  // End ErrorHandle Subroutine
//******************************************************************************


/******************************************************************************
Function:        ParseCommandLine
Description:     This function evaluates the command-line-parameters
                 2 styles are supported.
                 old-style (which has a fixed order of parameters)
                 new-style (which receives only a few arguments and defaults the rest)
                 If any new-style-parameters are recognized (infile=... ) then new-style is
                 assumed for all and evaluation is done accordingly.
//******************************************************************************/
int ParseCommandLine(char *preferred_defaultfilename,
                     char **infile_filename,  char **inseed_filename,
                     char **interaction_filename,
                     char **infoout_filename, char **posterior_filename,
                     char **pathout_filename, char **roundout_filename,
                     int argc, char *argv[])
{
  static char staticfilename0[260];
  static char staticfilename1[260];
  static char staticfilename2[260];
  static char staticfilename3[260];
  static char staticfilename4[260];

  int CommandLineType=0, validcommand;
  char *cp, errbuf[200];
  int DebugShowParameters;

  //If no parameters given show usage of the program
  if (argc<2)
  {
     cerr << "\nUsage1:     DrEvol [parameters]";
     cerr << "\nUsage2:     DrEvol +     uses the defaults and shows them on the screen.";
     cerr << "\n";
     cerr << "\nParameters:";
     cerr << "\ninfile=...               filename for sequence data + information";
     cerr << "\ninseed=...               filename for integer seed for rnd";
     cerr << "\ninteraction=...          filename for file that holds all constant info + NRG's";
     cerr << "\noutfile=...              filename (info.out)";
     cerr << "\npostout=...              filename for posterior sample points";
     cerr << "\npathout=...              filename for path energy terms";
     cerr << "\nroundout=...             filename round.out";
     cerr << "\n";
     cerr << "\nshow_commandline=...     yes OR no  displays commandline-parameters";
     cerr << "\n+                        a single + as any parameter shows the parameters too";
     cerr << "\n";
     return 0;
  }


  //PRF_PUSH("calc---ParseCommandLine");

  //----------Protocoll the used parameters to logfile
  /*
  #if DBG_LOG_COMMANDLINE==DBG_ON
      char tmp_buffer[300];
      strcpy(tmp_buffer, "Calculate_likelihood ");
      for (int i=1; i<argc; i++)
      {
          sprintf(&tmp_buffer[strlen(tmp_buffer)], " %s", argv[i]);
      }
      sprintf(&tmp_buffer[strlen(tmp_buffer)], "%c", 10);
      ofstream commandline_logfile("commands.log",ios::app);
      commandline_logfile.write(tmp_buffer,strlen(tmp_buffer));
      commandline_logfile.close();
  #endif
  */
  //----------Protocoll the used parameters to logfile



  //----------Initialize all parameters to null-values
  *infile_filename       =NULL;
  *inseed_filename       =NULL;
  *interaction_filename  =NULL;
  *infoout_filename      =NULL;
  *posterior_filename    =NULL;
  *pathout_filename      =NULL;
  *roundout_filename     =NULL;
  DebugShowParameters    =0;
  //----------Initialize all parameters to null-values


  //----------loop over all arguments
  //Do we have any new-style-arguments ("default" or "infile=..."
  //if so then we assume that all our arguments are new-style
  strcpy(errbuf,"");
  for (int i=1; i<argc; i++)
  {
    if (i==1)
    {
         if ((strstr(argv[i], "usedef")!=NULL) || (strstr(argv[i], "default")!=NULL) || (!strcmp(argv[i], "+")))
            CommandLineType=1;
    }
    cp=strstr(argv[i], "=");

    //This is either a parameter-error or we are dealing with old-style-params, we don't know yet... 2Apr24
    if (cp==NULL)
    {
       if ((argv[i][0]=='?') || (argv[i][0]=='+'))                                  //2Apr25
       {
          DebugShowParameters=1; continue;                                          //2Apr25
       }
    }
    else
    {
        //was there an old-style-parameter before this new-style-param ?
        if (errbuf[0]!=0)
        {
            cerr << "\n" << errbuf; getchar(); return(0);
        }
    }

    if (cp!=NULL)
    {
        validcommand=0;
        cp=strstr(argv[i], "infile=");      if (cp!=NULL) {*infile_filename=&cp[7];         CommandLineType=1; validcommand=1;}
        cp=strstr(argv[i], "inseed=");      if (cp!=NULL) {*inseed_filename=&cp[7];         CommandLineType=1; validcommand=1;}
        cp=strstr(argv[i], "interaction="); if (cp!=NULL) {*interaction_filename=&cp[12];   CommandLineType=1; validcommand=1;}
        cp=strstr(argv[i], "outfile=");     if (cp!=NULL) {*infoout_filename=&cp[8];        CommandLineType=1; validcommand=1;}
        cp=strstr(argv[i], "postout=");     if (cp!=NULL) {*posterior_filename=&cp[8];      CommandLineType=1; validcommand=1;}
        cp=strstr(argv[i], "pathout=");     if (cp!=NULL) {*pathout_filename=&cp[8];        CommandLineType=1; validcommand=1;}
        cp=strstr(argv[i], "roundout=");    if (cp!=NULL) {*roundout_filename=&cp[9];       CommandLineType=1; validcommand=1;}
        cp=strstr(argv[i], "show_commandline=");
        if (cp!=NULL)
        {
            validcommand=1;
            if ((cp[17]=='y') || (cp[17]=='Y') || (cp[17]=='1'))
                {DebugShowParameters=1; CommandLineType=1;}
            else
                {DebugShowParameters=0; CommandLineType=1;}
        }

        //unknown command, issue a warning !
        if (!validcommand)
        {
            sprintf(errbuf,"Parameter error !  '%s' is not a valid parameter !",argv[i]);
            cerr << "\n" << errbuf; getchar();
        }
     }
  }
  //----------loop over all arguments



  //----------new style parameters
  //DEFAULTS:

  //We have found new-style-parameters, so we fill all omitted parameters with defaults
  if (CommandLineType==1)
  {
        if (*infile_filename==NULL)
        {
            strcpy(staticfilename0,preferred_defaultfilename);
            *infile_filename=staticfilename0;
        }
        if (*interaction_filename==NULL)  *interaction_filename="Interaction.dat";

        //Generate the basename from the infile
        char basename[260];
        strcpy (basename, *infile_filename);
        cp=strstr(basename, ".dat");
        if (cp!=NULL)
            cp[0]=0;
        else
        {
            cp=strstr(basename, ".Dat");
            if (cp!=NULL)
                cp[0]=0;
            else
            {
               cp=strstr(basename, ".DAT");
               if (cp!=NULL)
                   cp[0]=0;
               else
               {
                   for (int i=strlen(basename);i>0;i--)
                   {
                      if (basename[i]=='.') {basename[i]=0;break;}
                   }
               }
            }
        }


        if (*inseed_filename==NULL)             *inseed_filename="inseed.dat";

        if (*infoout_filename==NULL)
        {
            strcpy(staticfilename1,basename); strcat(staticfilename1,".info.out");
            *infoout_filename=staticfilename1;
        }
        if (*posterior_filename==NULL)
        {
            strcpy(staticfilename2,basename); strcat(staticfilename2,".out");
            *posterior_filename=staticfilename2;
        }
        if (*pathout_filename==NULL)
        {
            strcpy(staticfilename3,basename); strcat(staticfilename3,".NRGpic.out");
            *pathout_filename=staticfilename3;
        }
        if (*roundout_filename==NULL)
        {
            strcpy(staticfilename4,basename); strcat(staticfilename4,".round.out");
            *roundout_filename=staticfilename4;
        }
  }
  //----------new style parameters

  else

  //----------old style parameters
  {
     //We have found no new-style-parameters so we just fill the pointers
     //in old-style (job-to-do=argv[1], ...)
     if (argc >= 7)
     {
        *infile_filename     =argv[1];
        *inseed_filename     =argv[2];
        *interaction_filename=argv[3];
        *infoout_filename    =argv[4];
        *posterior_filename  =argv[5];
        *pathout_filename    =argv[6];
        *roundout_filename   =argv[7];
     }
     else
     {
        sprintf(errbuf,"Parameter error !  Old-style parameters, - not enough params given !");
        cerr << "\n" << errbuf; getchar();
        return(0);
     }
  }
  //----------old style parameters


  //----------validate plausible values
  //No matter if we used new-style or old-style, now is the time to check the
  //validity of the data.  We return 0 if we find a problem.
  if (
     (!argc) ||
     (*infile_filename==NULL)    || (*inseed_filename==NULL)    || (interaction_filename==NULL) ||
     (*infoout_filename==NULL)   || (*posterior_filename==NULL) || (*pathout_filename==NULL)    ||
     (*roundout_filename==NULL)
     )
  {
    //PRF_POP(0L);
    return 0;
  }
  //----------validate plausible values


  //----------debug-print all parameters
  //Just for debug-purposes, show all command-line-parameters we recognized
  char dbuf[200];
  if (DebugShowParameters==1)             //2Mar19
  {
     cout << endl << endl << endl;
     cout << "------------------------------------------------------------------" << endl;
     cout << "DrEvol.cpp  called with the following command-line-parameters:" << endl;
     cout << endl;
     for (int i=0; i<argc; i++) {cout << argv[i] << " ";} cout << endl << endl;
     sprintf(dbuf,"p1:   infile=      %s",*infile_filename);      cout << dbuf << endl;
     sprintf(dbuf,"p2:   inseed=      %s",*inseed_filename);      cout << dbuf << endl;
     sprintf(dbuf,"p3:   interaction= %s",*interaction_filename); cout << dbuf << endl;
     sprintf(dbuf,"p4:   outfile=     %s",*infoout_filename);     cout << dbuf << endl;
     sprintf(dbuf,"p5:   postout=     %s",*posterior_filename);   cout << dbuf << endl;
     sprintf(dbuf,"p6:   pathout=     %s",*pathout_filename);     cout << dbuf << endl;
     sprintf(dbuf,"p7:   roundout=    %s",*roundout_filename);    cout << dbuf << endl;

     cout << "Please press any key...";
     getchar();
     cout << endl << endl << endl;
  }
  //----------debug-print all parameters

  //PRF_POP(0L);
  return 1;
}

