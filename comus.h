#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


#define NS 5000
#define NBRANCH (NS*2-2)
#define MAXNSONS 20
#define LSPNAME 50
#define NCODE 64
#define NCATG 40
#define ADDSIZE 50
#define SITESINC 10 
#define MAXTREESTRING 10000
#define FILESIZE 256
#define MAXTRIALS 500


enum {PrBranch=1, PrNodeNum=2, PrLabel=4, PrAge=8, PrOmega=16} OutTreeOptions;

struct CommonInfo {
  char *z[2*NS-1], spname[NS][LSPNAME+1], daafile[96], cleandata, readpattern;
  int ns, ls, npatt, np, ntime, ncode, clock, rooted, model, icode;
  int seqtype, *pose, ncatG, NSsites;
  double *fpatt, kappa, omega, alpha, pi[64], *conP, daa[20*20];
  double freqK[NCATG], rK[NCATG];
  char *siteID;    /* used if ncatG>1 */
  double *siterates;   /* rates for gamma or omega for site or branch-site models */
  double *omegaBS, *QfactorBS;     /* omega IDs for branch-site models */
}  com;

struct devent {
  double time;
  int popi;
  int popj;
  int isSpeciation;
  double paramv;
  double **mat ;
  char detype ;
  struct devent *nextde;
} ;

struct deventArray {
  int type;
  int nevents;
  double time;
  int sp1;
  int sp2;
  int popi;
  int popj;
  int isSpeciation;
  double paramv;
  double **mat;
  char detype;
  int id;
} *pastEvents;

struct c_params {
  int npop;
  
  int nsam;

  int *config;
  
  double *samplingTime;

  double **mig_mat;

  int **ignoreSpeciationJoints;

  double r;
  int nsites;
  int msites;
  double f;
  double track_len;
  double *size;
  double *alphag;
  struct devent *deventlist ;
  
} ;


struct m_params {
  double theta;
  int segsitesin;
  int treeflag;
  int timeflag;
  int mfreq;
} ;


struct params { 
  struct c_params cp;
  struct m_params mp;
  int commandlineseedflag ;
  int output_precision;
};


struct initial_devent {
  double time;
  int popi;
  int popj;
  double paramv;
  double **mat ;
  char detype ;
  struct devent *nextde;
} ;

struct initial_c_params {
  int npop;
  int nsam;
  int *config;
  double **mig_mat;
  int **ignoreSpeciationJoints;
  double r;
  int nsites;
  double f;
  double track_len;
  double *size;
  double *alphag;
  struct devent *deventlist ;
} ;

struct initial_m_params {
  double theta;
  int segsitesin;
  int treeflag;
  int timeflag;
  int mfreq;
} ;

struct initial_params { 
  struct c_params cp;
  struct m_params mp;
  int commandlineseedflag ;
  int output_precision;
};
	

struct TREEB {
  int nbranch, nnode, root, branches[NBRANCH][2], leaves, internalNodes;
}  tree;


struct TREEN {
  int father, nson, sons[MAXNSONS], ibranch, nodeID, spID, originalIndex,
    sampleSize, *samplePops, npop, offset, *pops, maxIndPop, sumofnodespop, isLeaf, nseqs;
  
  double branch, age, omega, label, *conP, coalAge, isolationTime, 
    theta, rho, Ne, partialDuration, partialMaxMigration;
  
  char *nodeStr, fossil;
  
}  *nodes, *initialNodes;


/* struct SPECIES { */
/*   int sampleSize, /\* total sample size for the species *\/ */
/*     *samplePops, /\* sample size of each pop of the species *\/ */
/*     npop, /\* the number of populations that the species contains *\/ */
/*     offset; /\* from which population we should start to count the species *\/ */
  
/*   double theta,  /\* theta = 4Nmu per species *\/ */
/*     rho,  /\* rho = 4Nr per species *\/ */
/*     Ne; /\* relative population sizes for species *\/ */
      
/*   //struct devent *deventlist ; */
/* } *spConf; */



/* struct MIG{ */
/*   /\* stores the migration events set by the command line *\/ */
/*   int events,  */
/*     *fromSpecies,  */
/*     *toSpecies, */
/*     *fromPop,  */
/*     *toPop; */
/* } cmdMigEvents; */
  
struct GLOBAL{
  
  int multipleSpecies;

  int partialIsolationModel;
  
  int readMigration;

  int readSampling;

  int readN;

  int readG;

  int readeG;

  int readej;
  
  int smt;

  int readeN;
  
  int totalPops;

  int currentEvents;

  int maxEventsSize;

  int finiteSiteModel;

  int printSeparator;

  int printInSeparateFiles;
  
} globalVar;



struct segl {
  int beg;
  int length;
  int recLength;
  struct node *ptree;
  int next;
}  ;


/* source code based on Fletcher and Yang's INDELible */
void randomTreeGeneration(int ntaxa, double birth, double death, double sample, double mut, int option, char *outputTree, int mode, double torigin, FILE *phyloTreeOut, double oldestOrigin, double *samplingTime );

void speciationToCoalescentTimes(double theta, double scale); 

int OutTreeN (FILE *fout, int spnames, int printopt);

struct params pars ;	

struct params initial_pars;

void seedit( const char * ) ;

void getpars( int argc, char *argv[], int *howmany, FILE **CoalescentFile, char CoalescentFileName[FILESIZE]  )  ;

void getparsCMD_multiSpeciesSimulate(int argc, char *argv[], int *phowmany, FILE **phyloOut, char phyloFileName[FILESIZE], FILE **CoalescentFile, char CoalescentFileName[FILESIZE], double *partIsoPeriod, double *partMaxMigration, double *samplingRate, double *mu, double *lambda, double *mut, char phyloInputFileName[FILESIZE], FILE **phyloInputFile, int *phyloMode, double *torigin, double *oldestOrigin , unsigned int *timeSeed);

void addtoelist( struct devent *pt, struct devent *elist ); 

int gensam( char **list, double *probss, double *ptmrca, double *pttot, struct segl **seglst, FILE *ResultFile, int *nsegs, FILE* InfoFile, int simindex) ;

void constructUpdatedCommandLine();

void addtoelist( struct devent *pt, struct devent *elist ); 

void argcheck( int arg, int argc, char ** ) ;


int commandlineseed( char ** ) ;

void free_eventlist( struct devent *pt, int npop );

void phyTreeToCoalEvents(int *current_node, int *smallerSon);

void phyTreeToCoalEvents2();

void nodesToCoalEvents();

int pick2(int n, int *i, int *j);

int isseg(int start, int c, int *psg);

int links(int c);

int xover(int nsam,int ic, int is);

int ca(int nsam, int nsites, int c1, int c2);

int pick2_chrom(int pop, int config[], int *pc1, int *pc2);

int ranvec(int n, double pbuf[]);

void free_streec_variables();

void free_rand1_file();

int  readTreeFromFile(FILE *infile, int nspecies, int *positionsArray);

int readNumberOfReplications(FILE *infile);

int checkTreeSamplingConsistency(double *samplingTime, int n);
