/************* comus ********************



CoMuS: Simulating coalescent histories and polymorphic data from multiple species

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

CoMuS is heavily based on Hudson's ms and Rambaut's Seq-Gen. Thus, we have also put
the preambles of those programs. 


----------------------- ms.c preamble ----------------------------

*
*       Generates samples of gametes ( theta given or fixed number 
*						of segregating sites.)
	Usage is shown by typing ms without arguments.   
 usage: ms nsam howmany  -t  theta  [options]
 or
 ms nsam howmany -s segsites  [options] 

 nsam is the number of gametes per sample.
 howmany is the number of samples to produce.
 With -t the numbers of segregating sites will randomly vary 
 from one sample to the next.
 with -s segsites,  the number of segregating sites will be
 segsites in each sample.

 Other options: See msdoc.pdf or after downloading and compiling, type ms<CR>.


 *Arguments of the options are explained here:

 npop:  Number of subpopulations which make up the total population
 ni:  the sample size from the i th subpopulation (all must be 
 specified.) The output will have the gametes in order such that
 the first n1 gametes are from the first island, the next n2 are
 from the second island, etc.
 nsites: number of sites between which recombination can occur.
 theta: 4No times the neutral mutation rate 
 rho: recombination rate between ends of segment times 4No
 f: ratio of conversion rate to recombination rate. (Wiuf and Hein model.)
 track_len:  mean length of conversion track in units of sites.  The 
 total number of sites is nsites, specified with the -r option.
 mig_rate: migration rate: the fraction of each subpop made up of
 migrants times 4No. 
 howmany: howmany samples to generate.

 Note:  In the above definition, No is the total diploid population if
 npop is one, otherwise, No is the diploid population size of each
 subpopulation. 
 A seed file called "seedms" will be created  if it doesn't exist. The
 seed(s) in this file will be modified by the program. 
 So subsequent runs
 will produce new output.  The initial contents of seedms will be
 printed on the second line of the output.
 Output consists of one line with the command line arguments and one
 line with the seed(s).
 The samples appear sequentially following that line.
 Each sample begins with "//", then the number of segregating sites, the positions
 of the segregating sites (on a scale of 0.0 - 1.0). On the following
 lines are the sampled gametes, with mutants alleles represented as
 ones and ancestral alleles as zeros.
 To compile:  cc -o ms  ms.c  streec.c  rand1.c -lm
 or:  cc -o ms ms.c streec.c rand2.c -lm
 (Of course, gcc would be used instead of cc on some machines.  And -O3 or 
 some other optimization switches might be usefully employed with some 
 compilers.) ( rand1.c uses drand48(), whereas rand2.c uses rand() ).

 *
 *   Modifications made to combine ms and mss on 25 Feb 2001
 *	Modifications to command line options to use switches  25 Feb 2001
 *	Modifications to add // before each sample  25 Feb 2001
 Modifications to add gene conversion 5 Mar 2001
 Added demographic options -d  13 Mar 2001
 Changed ran1() to use rand(). Changed seed i/o to accomodate this change. 20 April.
 Changed cleftr() to check for zero rand() .13 June 2001
 Move seed stuff to subroutine seedit()  11 July 2001
 Modified streec.c to handle zero length demographic intervals 9 Aug 2001
 Corrected problem with negative growth rates (Thanks to D. Posada and C. Wiuf) 13 May 2002
 Changed sample_stats.c to output thetah - pi rather than pi - thetah.  March 8 2003.
 Changed many command line options, allowing arbitrary migration matrix, and subpopulation
 sizes.  Also allows parameters to come from a file. Option to output trees.  Option to
 split and join subpopulations.   March 8, 2003. (Old versions saved in msold.tar ).
 !!! Fixed bug in -en option.  Earlier versions may have produced garbage when -en ... used. 9 Dec 2003
 Fixed bug which resulted in incorrect results for the case where
 rho = 0.0 and gene conversion rate > 0.0. This case was not handled
 correctly in early versions of the program. 5 Apr 2004.  (Thanks to
 Vincent Plagnol for pointing out this problem.) 
 Fixed bug in prtree().  Earlier versions may have produced garbage when the -T option was used.
 1 Jul 2004.
 Fixed bug in -e. options that caused problems with -f option  13 Aug 2004.
 Fixed bug in -es option, which was a problem when used with -eG. (Thanks again to V. Plagnol.) 6 Nov. 2004
 Added -F option:  -F minfreq  produces output with sites with minor allele freq < minfreq filtered out.  11 Nov. 2004.
 Fixed bug in streec.c (isseg() ).  Bug caused segmentation fault, crash on some machines. (Thanks
 to Melissa Jane Todd-Hubisz for finding and reporting this bug.)
 Added -seeds option 4 Nov 2006
 Added "tbs" arguments feature 4 Nov 2006
 Added -L option.  10 May 2007
 Changed -ej option to set Mki = 0 pastward of the ej event.  See msdoc.pdf.  May 19 2007.
 fixed bug with memory allocation when using -I option. This caused problems expecially on Windows
 machines.  Thanks to several people, including Vitor Sousa and Stephane De Mita for help on this one.
 Oct. 17, 2007.
 Modified pickb() and pickbmf() to eliminate rare occurrence of fixed mutations Thanks to J. Liechty and K. Thornton. 4 Feb 2010.
 Added -p n  switch to allow position output to be higher precision.  10 Nov 2012.

------------------------ seq-gen preamble -----------------------------------

Sequence Generator - seq-gen, version 1.3.3
Copyright (c)1996-2011, Andrew Rambaut & Nick Grassly
Institute of Evolutionary Biology, University of Edinburgh			
All rights reserved.                          

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. The names of its contributors may not be used to endorse or promote 
products derived from this software without specific prior written 
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


For Seq-Gen:

   Any feedback is very welcome.
   http://tree.bio.ed.ac.uk/software/seqgen/
   email: a.rambaut@ed.ac.uk


For CoMuS:

email: pavlidisp@gmail.com
web:   http://pop-gen.eu/wordpress/software/comus-coalescent-of-multiple-species

***************************************************************************/

#include "comus.h"
#include "tree.h"
#include "model.h"
#include "comus_evolve.h"
#include "global.h"
#include "nucmodels.h"
#include "aamodels.h"
#include "treefile.h"
#include "twister.h"

#define MAXSIZEMODELNAME  100

#define MINTIME 0.00000000001

/* prototypes */

TNode *NewNode(TTree *tree);

TTree *NewTree();

static unsigned maxsites = SITESINC ;

struct node{
  int abv;
  int ndes;
  float time;
};


static double *posit ;

static double segfac ;

static int count, ntbs, nseeds ;

/* seq-gen related variables */
int scaleTrees, scaleBranches, ancestorSeq, writeAncestors, writeRates, numDatasets;

int *partitionLengths;

double *partitionRates;

double treeScale, branchScale;

double *mutationTimes;

FILE* mutationTimesFile;

/******************************/


double gettime(void)
{
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

void speciationIsolationTimes()
{
  int i;
  
  for( i = 0; i < tree.nnode; ++i)
    nodes[i].isolationTime = (nodes[i].coalAge - nodes -> partialDuration ) > 0 ? nodes[i].coalAge - nodes -> partialDuration : MINTIME;


}


/* the function scales the node times to coalescent time units */
void speciationToCoalescentTimes(double theta, double scale)
{
  int i;

  for(i = 0; i < tree.nnode; ++i)
    {
  
      nodes[i].coalAge = (nodes[i].age * scale)/theta;
      
      /* printf("age of node %d is %.10f, scales to %.10f, %f, %f\n", i, nodes[i].age, nodes[i].coalAge, scale, theta); */
      
    }
}

void setInternalFlag()
{
  int i;
  for( i = 0; i < 2 * tree.leaves - 1; ++i)
    {
      if( nodes[i].nson > 0)
	nodes[i].isLeaf = 0;
      else
	nodes[i].isLeaf = 1;
    }
}

int compare_struct_type( const void *a_, const void *b_)
{
  struct deventArray *a = (struct deventArray*)a_;
  
  struct deventArray *b = (struct deventArray*)b_;

  if( a->type > b->type ) return 1;
  
  if( b->type > a->type ) return -1;

  return 0;

}

int compare_struct( const void *a_, const void *b_)
{
  struct deventArray *a = (struct deventArray*)a_;
  struct deventArray *b = (struct deventArray*)b_;

  if( a->time > b->time ) return 1;
  if( b->time > a->time ) return -1;

  return 0;

}


void allocateParValues()
{
  int i;
  
  pars.cp.size = calloc( globalVar.totalPops, sizeof( double ) );
  
  assert( pars.cp.size != NULL);

  pars.cp.mig_mat = calloc( globalVar.totalPops, sizeof(double * ) );
  
  assert( pars.cp.mig_mat != NULL);

  pars.cp.ignoreSpeciationJoints = calloc( globalVar.totalPops, sizeof( int*) );

  assert( pars.cp.ignoreSpeciationJoints != NULL);

  for( i = 0; i < globalVar.totalPops; ++i)
    {
      pars.cp.mig_mat[i] = calloc( globalVar.totalPops, sizeof(double) );
      assert( pars.cp.mig_mat[i] != NULL);
      pars.cp.ignoreSpeciationJoints[i] = calloc( globalVar.totalPops, sizeof(int) );
      assert( pars.cp.ignoreSpeciationJoints[i] != NULL );
      
    }

  pars.cp.config = calloc( globalVar.totalPops, sizeof( int ) );

  pars.cp.samplingTime = calloc( globalVar.totalPops, sizeof(double) );

  pars.cp.alphag = calloc( globalVar.totalPops, sizeof( double ) );

}


void setParsFromInitialValues()
{
  int i, j, k;
  pars.cp.nsam = initial_pars.cp.nsam;

  
  for( i = 0; i < globalVar.totalPops; ++i)
    pars.cp.alphag[i] = initial_pars.cp.alphag[i];

  
  for( i = 0; i < globalVar.totalPops; ++i)
    pars.cp.samplingTime[i] = initial_pars.cp.samplingTime[i];

  
  for( i = 0; i < globalVar.totalPops; ++i)
    pars.cp.size[i] = initial_pars.cp.size[i];

  for( i = 0; i < globalVar.totalPops; ++i)
    for( j = 0; j < globalVar.totalPops; ++j)
      {
	pars.cp.mig_mat[i][j] = initial_pars.cp.mig_mat[i][j];
	pars.cp.ignoreSpeciationJoints[i][j] = initial_pars.cp.ignoreSpeciationJoints[i][j];
      }

  pars.commandlineseedflag = initial_pars.commandlineseedflag;

  pars.output_precision = initial_pars.output_precision;

  pars.cp.r = initial_pars.cp.r;

  pars.mp.theta = initial_pars.mp.theta;

  pars.cp.f = initial_pars.cp.f;

  pars.cp.npop = initial_pars.cp.npop;

  k = 0;
  
  for( i = 0; i < tree.leaves; ++i)
    for( j = 0; j < nodes[i].npop; ++j)
      {
	pars.cp.config[k] = nodes[i].samplePops[j];
	k++;
      }
  
  pars.cp.track_len = initial_pars.cp.track_len;

  pars.cp.npop = initial_pars.cp.npop;

  pars.mp.segsitesin =  initial_pars.mp.segsitesin;

  pars.mp.treeflag = initial_pars.mp.treeflag;

  pars.mp.mfreq = initial_pars.mp.mfreq;

  pars.cp.nsites = initial_pars.cp.nsites;

  pars.cp.msites = initial_pars.cp.msites;

}

void initialize_global_values()
{
  globalVar.partialIsolationModel = 0;

  globalVar.readMigration = 0;

  globalVar.readSampling = 0;
  
  globalVar.readN = 0;
  
  globalVar.currentEvents = 0;
  
  globalVar.maxEventsSize = ADDSIZE;

  globalVar.readeN = 0;

  globalVar.smt = 0;

  globalVar.finiteSiteModel = 1;

  globalVar.printSeparator = 1;

  globalVar.printInSeparateFiles = 0;

  initial_pars.cp.msites = pars.cp.msites = 100;

}

void createEventListFromEventArray(struct deventArray *pastEvents, struct devent *list)
{

  int i;
  
  struct devent *pt, *ptemp;
  
  /* put events in the list */
  for( i = 0; i < globalVar.currentEvents; ++i)
    {
      pt = (struct devent *) malloc( sizeof( struct devent) );
      
      pt->detype = pastEvents[i].detype;

      pt->isSpeciation = pastEvents[i].isSpeciation;
      
      pt->time = pastEvents[i].time;

      pt->popi = pastEvents[i].popi;

      pt->popj = pastEvents[i].popj;

      pt->paramv = pastEvents[i].paramv;

      pt->nextde = NULL;

      /* set the matrix */
      if( pastEvents[i].detype == 'a')
	{
	  
	  int pop, pop2;

	  pt->mat = calloc( globalVar.totalPops, sizeof( double* ) );
	  for( pop = 0; pop < globalVar.totalPops; ++pop)
	    {
	      pt->mat[pop] = calloc( globalVar.totalPops, sizeof(double));

	      for(pop2 = 0; pop2 < globalVar.totalPops; ++pop2)
		pt->mat[pop][pop2] = pastEvents[i].mat[pop][pop2];
	    }
	}
      
      if( pars.cp.deventlist == NULL)
	pars.cp.deventlist = pt;

      else if( pt->time < pars.cp.deventlist -> time)
	{
	  ptemp = pars.cp.deventlist;
	  pars.cp.deventlist = pt;
	  pt->nextde = ptemp;
	}
      else
	{
	  addtoelist ( pt, pars.cp.deventlist );
	}
      
    }

}

int getMaxInt(int *a, int n)
{
  if( n == 0)
    return 0;
  
  int i, max = a[0];

  assert(n > 0);
  
  for( i = 1; i < n; ++i)
    max = (a[i] > max) ? a[i] : max;

  return max;

}

void setInitialNodes()
{
  
  int i,j;

  initialNodes = malloc((tree.leaves * 2 -1)*sizeof(struct TREEN) );

  for(i = 0; i < tree.nnode; ++i)
    {
      initialNodes[i].offset = nodes[i].offset;
      initialNodes[i].sampleSize = nodes[i].sampleSize;
      initialNodes[i].npop = nodes[i].npop;
      initialNodes[i].nson = 0;
      initialNodes[i].maxIndPop = nodes[i].maxIndPop;
      initialNodes[i].pops = NULL;
      initialNodes[i].samplePops = NULL;
      initialNodes[i].nseqs = 0;

    }

  for( i =0; i < tree.leaves; ++i)
    {
      
      initialNodes[i].samplePops = calloc( nodes[i].npop, sizeof(int) );
      
      for(j = 0; j < nodes[i].npop; ++j)
	initialNodes[i].samplePops[j] = nodes[i].samplePops[j];
      
    }
      
}


void updateTreePops(struct TREEN *nodes)
{
 
  int i,j, k, s, son;

  nodes->sumofnodespop = 0;

  for( i = tree.leaves; i < tree.nnode; ++i)
    nodes[i].npop = 0;
  
  for( i =0; i < 2*tree.leaves - 1; ++i)
    {
      /* here we do not enter in the leaves */
      for(j = 0; j < nodes[i].nson; ++j)
	{
	  assert( i > tree.leaves - 1);
	  nodes[i].npop += nodes[ nodes[i].sons[j] ].npop;
	}
    }
  

  /* I need to update the pops because 
     it matters when the tree has been simulated
     or it has been read from the file
  */
  j = 0;
  for( i = 0; i < tree.leaves; ++i)
    {
      
      /* printf("*i: %d, npop: %d\n", i, nodes[i].npop); */
      if(nodes[i].pops == NULL )
  	nodes[i].pops = calloc( nodes[i].npop, sizeof(int) );
      else
  	{
  	  free( nodes[i].pops );
  	  nodes[i].pops = calloc( nodes[i].npop, sizeof(int) );
  	}
      
      for( k = 0; k < nodes[i].npop; ++k)
  	nodes[i].pops[k] = j++;
      
    }

  
  for( i = tree.leaves; i < 2*tree.leaves - 1; ++i)
    {
      nodes->sumofnodespop += nodes[i].npop - 1;
      
      if(nodes[i].pops == NULL )
	nodes[i].pops = calloc( nodes[i].npop, sizeof(int) );
      else
	{
	  free(nodes[i].pops);
	  
	  nodes[i].pops = calloc( nodes[i].npop , sizeof(int) );
	}
          
      j = 0;
      for( s = 0; s < nodes[i].nson; ++s)
	{
	  son = nodes[i].sons[s];
	  
	  for( k = 0; k < nodes[ son ].npop; ++k)
	    {
	      assert( j < nodes[i].npop);

	      /* printf("pops1: %d\n", nodes[i].pops[j]); */

	      /* printf("pops2*: son: %d, k: %d\n", son, k); */

	      /* printf("pops2: %d, son: %d, k: %d\n", nodes[son].pops[k], son, k); */
	      
	      nodes[i].pops[j] = nodes[ son ].pops[k];

	      j += 1;
	    }
	}
    }

  
  for( i = 0; i < 2 * tree.leaves - 1; ++i)
    {
      nodes[i].maxIndPop = getMaxInt( nodes[i].pops, nodes[i].npop );
    }

}



int order(int n, double pbuf[])
{
  int gap, i, j;
  double temp;

  for( gap= n/2; gap>0; gap /= 2)
    for( i=gap; i<n; i++)
      for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
	temp = pbuf[j];
	pbuf[j] = pbuf[j+gap];
	pbuf[j+gap] = temp;
      }

  return 0;
}

/****  tdesn : returns 1 if tip is a descendant of node in *ptree, otherwise 0. **/

int tdesn(struct node *ptree, int tip, int node )
{
  int k;

  for( k= tip ; k < node ; k = (ptree+k)->abv ) ;
  if( k==node ) return(1);
  else return(0);
}


void free_global()
{
  
  int i;

  free(pars.cp.config);
  free( pars.cp.samplingTime);
  free(pars.cp.alphag);
  
  free( initial_pars.cp.size);
  free( initial_pars.cp.alphag);
  free( initial_pars.cp.samplingTime);

  for( i = 0; i < globalVar.totalPops; ++i)
    {
      free( initial_pars.cp.mig_mat[i]);
      free( initial_pars.cp.ignoreSpeciationJoints[i]);
      free( pars.cp.ignoreSpeciationJoints[i]);
      free( pars.cp.mig_mat[i]);
    }
  
  free( initial_pars.cp.mig_mat);
  free( pars.cp.mig_mat);
  
  
  free( initial_pars.cp.ignoreSpeciationJoints);
  free( pars.cp.ignoreSpeciationJoints);
  
  for(i = 0; i < globalVar.totalPops; ++i)
    {

    }

  for( i = 0; i <2*tree.leaves - 1; ++i)
    {
      if( nodes[i].samplePops != NULL)
	free( nodes[i].samplePops );
      
      if( nodes[i].pops != NULL )
	free(nodes[i].pops);
    }

  /* free(pars.cp.deventlist); */

  free( pars.cp.size);

  free_streec_variables();

  free(nodes);

  nodes = NULL;
    
  free( posit );

  posit = NULL;

  free( pastEvents );

  pastEvents = NULL;

  free_rand1_file();
}

	
int usage()
{

  fprintf(stderr, "\n##################################################\n");
  fprintf(stderr, "#       comus v 2.0                              #\n");
  fprintf(stderr, "#      ------------                              #\n");
  fprintf(stderr, "# January 2016                                   #\n");
  fprintf(stderr, "# Developed by P. Pavlidis and S. Papadantonakis #\n");
  fprintf(stderr, "##################################################\n");
  fprintf(stderr, "\nProperties:\n");
  fprintf(stderr, "\t--infinite or finite site model.\n");
  fprintf(stderr, "\t--single or multiple species polymorphic data\n");
  fprintf(stderr, "\t--various birth-death models for speciation\n");
  fprintf(stderr, "\t--ancestral sampling + present day sampling (e.g. co-analysis of fossils with modern samples)\n");
  fprintf(stderr, "\t--substitution rate heterogeneity (gamma categories) between sites\n");
  fprintf(stderr, "\t--gradual speciation through gradual isolation\n\n\n");
  
  fprintf(stderr, "Infinite site model and single species simulations are similar to Hudson's ms\n");
  //fprintf(stderr, "The command line is the same:\n\n");
  fprintf(stderr,"------------------Infinite Site Model, Single Species---------------\nusage: comus nsam howmany \n");
  fprintf(stderr,"  Options: \n"); 
  fprintf(stderr,"\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
  fprintf(stderr,"\t -s segsites   ( fixed number of segregating sites)\n");
  fprintf(stderr,"\t -T          (Output gene tree.)\n");
  fprintf(stderr,"\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
  fprintf(stderr,"\t -r rho nsites     (rho here is 4Nc)\n");
  fprintf(stderr,"\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
  fprintf(stderr,"\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n"); 
  fprintf(stderr,"\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");      
  fprintf(stderr,"\t -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");    
  fprintf(stderr,"\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n"); 
  fprintf(stderr,"\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n"); 
  fprintf(stderr,"\t\t -n i size_i   (popi has size set to size_i*N0 \n");
  fprintf(stderr,"\t\t -g i alpha_i  (If used must appear after -M option.)\n"); 
  fprintf(stderr,"\t   The following options modify parameters at the time 't' specified as the first argument:\n");
  fprintf(stderr,"\t -eG t alpha  (Modify growth rate of all pop's.)\n");     
  fprintf(stderr,"\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");    
  fprintf(stderr,"\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n"); 
  fprintf(stderr,"\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n"); 
  fprintf(stderr,"\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");  
  fprintf(stderr,"\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");    
  fprintf(stderr,"\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
  fprintf(stderr,"\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");    
  fprintf(stderr,"\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
  fprintf(stderr,"\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
  fprintf(stderr,"\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
  fprintf(stderr,"\t\t  size, alpha and M are unchanged.\n");  
  fprintf(stderr,"\t  -f filename     ( Read command line arguments from file filename.)\n"); 
  fprintf(stderr,"\t  -p n ( Specifies the precision of the position output.  n is the number of digits after the decimal.)\n");
  fprintf(stderr," See msdoc.pdf for explanation of these parameters.\n");
  

  /************* multiple species **************/

  fprintf(stderr, "---------------------- Multiple species or Finite site model -------------------------\n\n");

  fprintf(stderr, "Options: \n\n");

  fprintf(stderr, "For INFINITE site model do NOT provide -mm <MUTATION MODEL> e.g.\ncomus <nspecies> <n_sp1> <n_sp2> ... <replications> -t <float> .... \n\n");

  fprintf(stderr, "For FINITE site model DO provide -mm <MUTATION MODEL> and the -msites <int> e.g.\ncomus <nspecies> <n_sp1> <n_sp2> ... <replications> -t <float> -mm hky -msites 10000 .... \n\n");

  fprintf(stderr, "-seed <int>:\tThe seed for random numbers. If absent the program is seeded by the time clock of the PC\n");

  fprintf(stderr, "-iphylo <file>:\tFile name with the newick tree. IMPORTANT: use ultrametric rooted tree... Otherwise results will be meaningless\n");

  fprintf(stderr, "-oldestOrigin <float>:\tdelimit the TMRCA of the phylogenetic tree to NOT be older than <float>\n");

  fprintf(stderr, "-phyloMode <int>:\n");
  fprintf(stderr, "\t0:ala Yang and Rannala, i.e. fix the TMRCA to 1.0, condition on the number of species and specify the birth rate (default:100), death rate(default: 0), and sampling proportion(default 1.0)\n");
  fprintf(stderr, "\t1:Condition on the number of taxa and fix the TMRCA to a user defined value. It's similar to mode 0, but now the user can define the TMRCA. For birth rate (default:100), death rate(default: 0), and sampling proportion(default 1.0) and TMRCA (default: 1.0)\n");
  fprintf(stderr, "\t2:Same as mode 1, but now condition on the time of the origin of the process, not on the time of the most recent common ancestor\n");
  fprintf(stderr, "\t3:Condition only on the number of taxa\n");
  fprintf(stderr, "\t4:Condition on the number of taxa and on the fact that the process cannot be older than a certain time\n");
  fprintf(stderr, "-torigin <float>:\t the time of origin or the TMRCA. This will depend whether you use phylomode 1 or 2\nIMPORTANT: USE THE -phyloMode FIRST. IN THIS CASE IT CAN BE EITHER 1 OR 2\n");
  fprintf(stderr, "-samplingrate (case insensitive) <float>: sets the sampling rate to <float> (default 1.0)\n");
  fprintf(stderr, "-partisolation\tUse the gradual isolation model. This flag make the partisolation = TRUE\n");
  fprintf(stderr, "-partmaxmigration <float>:\tSpecify the maximum migration rate for the gradual isolation model (default: 1.0). This flag also makes the partisolation = TRUE\n");
  fprintf(stderr, "-partisoperiod <float>:\tSpecify the time for which sister species remain in partial contact after a speciation event. Units are in coalescent terms (default: 0.0)\n");
  fprintf(stderr, "-oPhylo or -oPhylogeny:\tOutput the phylogenetic guide tree\n");
  fprintf(stderr, "-oCoalescent:\tOutput the Coalescent trees\n");
  fprintf(stderr, "-noSeparator:\tDo NOT print the '\\' between the fasta alignments. This is useful when only one alignment is simulated and you wish to view your alignment in other programs that do not like the '\\'\n");
  fprintf(stderr, "-oSeparateFiles:\tSeparate fasta alignments to different files\n");
  fprintf(stderr, "-t <float>:\tthe mutation parameter value\n");
  fprintf(stderr, "-msites <int>:\tNumber of sites. Should be greater than rsites. This parameters swithches on the finite site model\n");
  fprintf(stderr, "-r <float>: recombination rate\n");
  fprintf(stderr, "-rsites <int>: number of segmentes that can recombine\n");
  fprintf(stderr, "-G:\n");
  fprintf(stderr, "\t<float>:\tSets the growth rate to <float> for all populations\n");
  fprintf(stderr, "\t<int> <float>:\tSets the growth rate for species <int> to <float>\n");
  fprintf(stderr, "\t<int1> <int2> <float>:\tSets the growth rate for the population int2 of species int1 to float\n");
  fprintf(stderr, "-N:\n");
  fprintf(stderr, "\t<float>:\tSets the popsize to <float> for all populations(relative to 1)\n");
  fprintf(stderr, "\t<int> <float>:\tSets the popsize for species <int> to <float> (relative to 1)\n");
  fprintf(stderr, "\t<int1> <int2> <float>:\tSets the posize for the population int2 of species int1 to float (relative to 1)\n");
  fprintf(stderr, "-eN:\n");
  fprintf(stderr, "\t<float1> <float2>:\tSets the population size of all populations at time <float1> to <float2>\n");
  fprintf(stderr, "\t<float1> <int> <float2>:\tSets the population size of all populations of species <int> at time <float1> to <float2>\n");
  fprintf(stderr, "\t<float1> <int1> <int2> <float2>:\tSets at time <float1>, the population size of population <int2> of species <int1> to <float2>\n");
  fprintf(stderr, "-eN:\n");
  fprintf(stderr, "\t<float1> <float2>:\tSets the growth rate of all populations at time <float1> to <float2>\n");
  fprintf(stderr, "\t<float1> <int> <float2>:\tSets the growth rate of all populations of species <int> at time <float1> to <float2>\n");
  fprintf(stderr, "\t<float1> <int1> <int2> <float2>:\tSets at time <float1>, the growth rate of population <int2> of species <int1> to <float2>\n");
  fprintf(stderr, "-I <species int1>  <npop int2>  <ss int3> <ss int4> ...:\tDenotes that species int1 has int2 populations with sample sizes int3, int4, etc. IMPORTANT: do not specify migration rates here\n");
  fprintf(stderr,"-migration:\n");
  fprintf(stderr, "\t<float>:\tSet the migration rate between all pops to <float>\n");
  fprintf(stderr, "\t<int> <float>:\tSet the migration rate between all populations of species <int> to <float>\n");
  fprintf(stderr, "\t<int1> <int2> <float>:\tSets the migration rate between all populations of species int1 and species int2 to float\n");
  fprintf(stderr, "\t<int1> <int2> <int3> <float>:\tSets the migration rate between population int2 and int3 of species int1 to float\n");
  fprintf(stderr, "\t<int1> <int2> <int3> <int4> <float>:\tSets the migration rate between population int2 of species int1 and population int4 of species int3 to float\n");
  fprintf(stderr, "-ej:\n");
  fprintf(stderr, "\t<float> <int1>  <int2> <int3>:\tJoins at time <float> population int2 and int3 from species int1\n");
  fprintf(stderr, "\t<float> <int1> <int2> <int3> <int4>:\tJoins at time <float> population int2 from species int1 with population int4 from species int3\n");
  fprintf(stderr, "-birth:\tbirth rate for the speciation process\n");
  fprintf(stderr, "-death:\tdeath rate for the speciation process\n");
  

  fprintf(stderr, "\n\n---- the following flags refer to the finite site model\n\n");
  fprintf(stderr, "-name <string>:\tSpecify the name of the run\n");
  fprintf(stderr, "-oFormat <string>:\tSpecify the output format. Either FASTA or PHYLIP\n");
  fprintf(stderr, "-mm <STRING>:\tSpecify the mutation model (JC, HKY, GTR, F84)\n");
  fprintf(stderr, "-msites <int>:\tThe number of sites for the finite site model\n");
  fprintf(stderr, "-frequencies <float> <float> ...:\tThe frequencies of the nucleotides\n");
  fprintf(stderr, "-rates <float> <float> ...:\tThe relative rates of nucleotide substitution\n");
  fprintf(stderr, "-titv <float>:\ttransition to transversion rate\n");
  fprintf(stderr, "-gammaCategories <int>:\tSets the number of Gamma categories for taking into heterogeneity in nucleotide substitutions\nb");
  fprintf(stderr, "-alpha <float>:\tParameter of the gamma distribution for drawing heterogeinity rates. It is the shape parameter of the gamma distribution\b");
  fprintf(stderr, "-invariable <float>:\tProportion of invariable sites\n");
  
  // sampling time

  fprintf(stderr,"-samplingtime (case insensitive):\n");
  fprintf(stderr, "\t<float>:\tSet the sampling time for ALL pops to <float>\n");
  fprintf(stderr, "\t<int> <float>:\tSet the sampling time for ALL populations OF SPECIES  <int> to <float>\n");
  fprintf(stderr, "\t<int1> <int2> <float>:\tSet the sampling time for population <int2> of species <int1> to <float>\n");  
  
  exit(1);
}



void stringToUpper(char *text, char *nText){
  
  int i;
  
  for(i=0; i<=strlen(text); i++){
    if( (text[i] > 96 ) && (text[i] < 123) ) // is the char lower case
      nText[i] = text[i] - 'a' + 'A';   //make upper
    else
      nText[i] = text[i]; //do nothing
  }   
  
}


static int checkAllCommandLine(int argc, char **argv, int *unsupported)
{
  int success = 1,  i = 0;

  
  char stringTemp[MAX_NAME_LEN];

  for( i = 1; i < argc; ++i)
    {
      stringToUpper( argv[i], stringTemp );

      

      if(argv[i][0] != '-')
	continue;

      if( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 )
	{
	  usage();
	  exit(1);
	}

      if( strcmp(stringTemp, "-SAMPLINGTIME") == 0 )
	{
	  continue;
	}

      if( strcmp(argv[i], "-t") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-s") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-T") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-F") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-r") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-G") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-eG") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-I") == 0)
	{
	  continue;
	}
      
      if( strcmp(argv[i], "-ma") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-n") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-g") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-m") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-eG") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-eg") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-eM") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-ema") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-en") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-es") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-ej") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-p") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-SEED") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-IPHYLO") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-OLDESTORIGIN") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-PHYLOMODE") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-PARTISOPERIOD") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-PARTMAXMIGRATION") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-PARTISOLATION") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-TORIGIN") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-OPHYLO") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-OCOALESCENT") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-NOSEPARATOR") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-OFORMAT") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-OSEPARATEFILES") == 0)
	{
	  
	  continue;
	}

       if( strcmp("-MIGRATION", stringTemp) == 0 )
	 {
	   continue;
	 }

       
       if( strcmp("-ANCESTRALMIGRATION", stringTemp) == 0 )
	 {
	   continue;
	 }


      if( strcmp(argv[i], "-msites") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-rsites") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-N") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-eN") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-en") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-n") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-birth") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-mm") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-death") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-frequencies") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-name") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-frequencies") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-rates") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-titv") == 0)
	{
	  continue;
	}
      if( strcmp(stringTemp, "-GAMMACATEGORIES") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-alpha") == 0)
	{
	  continue;
	}
      if( strcmp(argv[i], "-invariable") == 0)
	{
	  continue;
	}

      success = 0;

      *unsupported = i;
    }
  return success;
  
  
}


static void resizePastEvents(int addsize)
{
  assert(addsize > 0);
  
  globalVar.maxEventsSize += addsize;
  
  pastEvents = realloc( pastEvents, globalVar.maxEventsSize * sizeof( struct deventArray ) );
  
  assert( pastEvents != NULL );
  
}


static void checkRequiredArguments(char name[MAX_NAME_LEN])
{
  if(name[0] == 0)
    {
      fprintf(stderr, "\n\nPlease provide a name, -name, for this run\n\n");
      assert(name[0] != 0);
    }
  
}

static int getparsMutationModel( int argc, char *argv[], char name[MAX_NAME_LEN] );

void printHeader(FILE *fout)
{
  fprintf(fout, "\n\n==================== comus v2.0================================================\n");
  fprintf(fout,     "An open source (GPLv3) software to simulate polymorphic samples from multiple species\n");
  fprintf(fout,     "under the infinite/finite site model and various demographic scenarios\n\n");
  fprintf(fout,     "Contacts:\n");
  fprintf(fout,     "Stefanos Papadantonakis:\tspapadadonakis@gmail.com\n");
  fprintf(fout,     "Pavlos Pavlidis:\t\tpavlidisp@gmail.com\n\n");
  fprintf(fout,     "Please cite: comus manuscript XXXX\n\n");

  fprintf(fout,     "comus is based on the following publications. Please cite them as well:\n\n");
  fprintf(fout,     "\t(1) Seq-Gen: an application for the Monte Carlo simulation of DNA sequence\n\tevolution along phylogenetic trees, Comput Appl Biosci. 1997 Jun;13(3):235-8.\n\n");
  fprintf(fout,     "\t(2) Generating samples under a Wright-Fisher neutral model of genetic\n\tvariation, Bioinformatics. 2002 Feb;18(2):337-8.\n");
  fprintf(fout,     "================================================================================\n\n");  
}

void printIncrementer(FILE *fout, int repl, int total, double currentTimeSeconds)
{
  int barsize = 50;
  
  static char bar1[1000] = "=";

  static char incrementer[256];
  
  int howmany = total/barsize;
  
  if( howmany > 0 && (repl % howmany ) == 0 )
    {
      sprintf(bar1, "%s", bar1);
  
      sprintf(incrementer, "ET: %3.2f, CT: %3.2f ... %s %2.1f%%", currentTimeSeconds * total/repl, currentTimeSeconds, bar1, (double)repl/total * 100.);

      fprintf(fout, "%s\r", incrementer);

      fflush(fout);
    }
  else if(howmany == 0)
    {

      sprintf(bar1, "%s", bar1);
  
      sprintf(incrementer, "ET: %3.2f, CT: %3.2f ... %s %2.1f%%", currentTimeSeconds * total/repl, currentTimeSeconds, bar1, (double)repl/total * 100.);

      fprintf(fout, "%s\r", incrementer);

      fflush(fout);

    }
  
  if(repl == total)
    fprintf(fout, "\n\n");
  


}


TNode *getNodeFromMS(struct node *ptree, TTree *tree, TNode *parent, int* descl, int *descr, int noden)
{
  TNode *node, *node2;
  
  //double len; //, param = 0.;

  char name[256];

  if( (node = NewNode(tree) ) == NULL )
    return NULL;
  
  if(descl[noden] == -1 )
    {
      /* node->length0 = (ptree + (ptree + noden)->abv ) -> time */
      
      node->length0 = (ptree + (ptree + noden)->abv ) -> time - (ptree+noden)->time ;

      /* printf("coal age of node %d is %f\n", noden, node->length0); */
      
      /* node->length2 = node2->length0 = (ptree + (ptree + descr[noden])->abv)->time - (ptree + descr[noden])->time  ; */

      node->branch1 = node->branch2 = NULL;
      node->length1 = node->length2 = 0.;
      node->branch0 = parent;
      node->tipNo = noden;

      
     
      assert( sprintf(name, "seq%d", noden) > 0);
 
      if (tree->names[node->tipNo]==NULL) {
	if ( (tree->names[node->tipNo]=(char *)malloc(MAX_NAME_LEN+1))==NULL )
	  {
	    fprintf(stderr, "Out of memory creating name.");
	    return NULL;
	  }
      }
      strcpy(tree->names[node->tipNo], name);

      tree->tips[node->tipNo] = node;
      
      tree->numTips++;
      
    }
  else
    {

      node2 = getNodeFromMS( ptree, tree, node, descl, descr, descl[noden]);
      
      node->branch1 = node2;
      
      node2->branch0 = node;
      node->length1 = node2->length0 = (ptree + (ptree + descl[noden])->abv)->time - (ptree + descl[noden])->time  ;
      
      
      node2 = getNodeFromMS( ptree, tree, node, descl, descr, descr[noden]);
      node->branch2 = node2;
      node2->branch0 = node;
      node->length2 = node2->length0 = (ptree + (ptree + descr[noden])->abv)->time - (ptree + descr[noden])->time  ;
      
      
      if( (ptree + noden)->abv == 0 )
	{
	  tree->root = node;
	  
	}
      
    }

  node->sampledNode = 1;

  return node;
}



void copyTree(struct node *ptree, TTree *tree, int nsam)
{

  
  int i, *descl, *descr ;

  descl = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
  descr = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );

  for( i=0; i<2*nsam-1; i++) 
    descl[i] = descr[i] = -1 ;

  for( i = 0; i< 2*nsam-2; i++){
    if( descl[ (ptree+i)->abv ] == -1 ) 
      descl[(ptree+i)->abv] = i ;
    else 
      descr[ (ptree+i)->abv] = i ;
  }
  
  getNodeFromMS(ptree, tree, NULL, descl, descr, 2*nsam - 2);
  
  free( descl ) ;
  free( descr ) ;
  
}


int getSequencesFiniteModel(struct segl *seglst, int nsegs, int nsam, FILE *ResultFile, TTree **treeSet, int *partitionLengths, double *partitionRates)
{
  int i, j, seg, k, start, end, len = 0, numTrees = nsegs, sumLength = 0, currentLength = 0, tempSum = 0, pop=0; //, nt;
            
  char *ancestor = NULL;
      
  double ratioMSitesNSites = 1., remainingRatio = 0., timeScaleChange = 1., scale = 1.;
  
  numTaxa = nsam;
      
  numPartitions = nsegs;

  CreateRates();
    
  ratioMSitesNSites = 1.;

  remainingRatio = 0.;
      
  assert(pars.cp.msites >= pars.cp.nsites);

  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {

     partitionRates[k] = 1.;
	
    /* if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) ){ */
	  
      end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : pars.cp.nsites-1 );
	  
      start = seglst[seg].beg ;
	  
      len = end - start + 1 ;

      /* printf("start: %d, end: %d, len: %d\n", start, end, len); */

      partitionLengths[k] = (len * pars.cp.msites)/pars.cp.nsites;
	  
      ratioMSitesNSites = (1. * len * pars.cp.msites)/(double)pars.cp.nsites;
	  
      remainingRatio += ratioMSitesNSites - partitionLengths[k];

      if(remainingRatio >= 0.99999999999999)
	{
	  remainingRatio -= 0.99999999999999;
	  partitionLengths[k]++;
	}

      tempSum += partitionLengths[k];

      sumLength += len;
	
      seglst[seg].length = partitionLengths[k];

      seglst[seg].recLength = len;
      
      int tmpi = 0;
      for( pop = 0; pop < pars.cp.npop; pop++)
	{
	  if( pars.cp.samplingTime[ pop ] > 0 )
	    {
	      for( i = 0; i < pars.cp.config[ pop ]; ++i)
		{
		  (seglst[seg].ptree + tmpi)->time = pars.cp.samplingTime[ pop ];
		  /* printf("sampling time of pop: %d is: %f\n", pop, pars.cp.samplingTime[ pop ]); */
		  tmpi++;
		}
	    }
	else
	  tmpi += pars.cp.config[pop];
      }

      
    copyTree(seglst[seg].ptree, treeSet[k], nsam);

    /* /\* decide whether tip nodes are sampled or not...  */
    /*    if not, then terminal branches will not be mutated. */
    /* *\/ */
    /* int tmpi = 0; */
    /* for( pop = 0; pop < pars.cp.npop; pop++) */
    /*   { */
    /* 	if( pars.cp.samplingTime[ pop ] > 0 ) */
    /* 	  { */
    /* 	    for( i = 0; i < pars.cp.config[ pop ]; ++i) */
    /* 	      { */
    /* 		treeSet[k]->tips[tmpi]->sampledNode = 0; */
    /* 		tmpi++; */
    /* 	      } */
    /* 	  } */
    /* 	else */
    /* 	  tmpi += pars.cp.config[pop]; */
    /*   } */

    /* for( i = 0; i < nsam; ++i) */
    /*   { */
    /* 	//treeSet[k]->tips[i]->sampledNode=0; */
    /* 	//printf("TIME of NODE %d is %f\n", i, (seglst[seg].ptree + i )->time); */
    /*   } */

    
    
    treeSet[k]->root->sequence = calloc(partitionLengths[k], sizeof(char) );

    treeSet[k]->rooted = 1;
      
    partitionRates[k] *= pars.mp.theta/pars.cp.msites;

    /* printf("rate for seg %d is %f\n", k, partitionRates[k]); */
      
    currentLength += len;
  }

  timeScaleChange = pars.mp.theta/pars.cp.msites;
      
  for (i = 0; i < nsegs; i++)
    CreateSequences(treeSet[i], partitionLengths[i]);
      
  if (nsegs > 1) 
    {
      double sumRates = 0.0;
      for (i = 0; i < nsegs; i++)
	sumRates += partitionRates[i] * partitionLengths[i];
      
      for (i = 0; i < nsegs; i++)
	{

	  partitionRates[i] *= (double)numSites / sumRates;

	  partitionRates[i] *= timeScaleChange;

	  /* printf("rate for seg %d is %f\n", i, partitionRates[i]); */

	}

    }
	
  /* per simulated dataset per tree, i.e. you can
     use the same tree (phylo+coal) for multiple Simulations
  */
  for (i=0; i<numDatasets; i++) {
    
    SetCategories();
	
    k = 0;
	
    for (j = 0; j < nsegs; j++) {
	 
      scale = partitionRates[j];
	  
      if (scaleTrees) 
	scale *= treeScale/treeSet[j]->totalLength;
      else if (scaleBranches)
	scale *= branchScale;
      
      EvolveSequences(treeSet[j], k, partitionLengths[j], scale, ancestor); 
	  	  
      k += partitionLengths[j]; 
	 	     
    }	
    
    WriteSequences(ResultFile, (numTrees > 1 ? nsegs+1 : -1), (numDatasets > 1 ? numDatasets+1 : -1), treeSet, partitionLengths, count);
		
  }

  freeRates();

  for ( i = 0; i < nsegs; ++i)
    freeSequences(treeSet[i]);

  for (i = 0; i < nsegs; i++) 
    DisposeTree(treeSet[i]);
  
  return 1;

}


void prtree( FILE *fout, struct node *ptree, int nsam, double scale)
{
  
  int i, *descl, *descr ;

  void parens(FILE *fout,  struct node *ptree, int *descl, int *descr, int noden, double scale );

  descl = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
  descr = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
  
  for( i=0; i<2*nsam-1; i++) 
    descl[i] = descr[i] = -1 ;
  
  for( i = 0; i< 2*nsam-2; i++){

    if( descl[ (ptree+i)->abv ] == -1 ) 
      descl[(ptree+i)->abv] = i ;
    
    else descr[ (ptree+i)->abv] = i ;
    
  }
  
  parens(fout,  ptree, descl, descr, 2*nsam-2, scale);
  
  free( descl ) ;
  free( descr ) ;
}


void printCoalescentTrees(FILE *fout, struct segl *seglst, int nsegs, int nsam, double scale)
{
  int seg, k;
  static int counter = 0;

  if(counter++ > 0 && globalVar.printSeparator  == 1) 
    {
      fprintf(fout, "\n//\n");
    }
  else if( globalVar.printSeparator  == 1)
    fprintf(fout, "//\n");

  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) 
    {
      if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) )
	{	  
	  fprintf(fout,"[%d]", seglst[seg].length);
	}
      prtree(fout,  seglst[seg].ptree, nsam, scale ) ;
      
    }

}


int setPartialIsolationEvents(struct deventArray *pastEvents, double duration)
{
  
  int i, counter = 0;
  for(i = 0; i < pastEvents -> nevents; ++i)
    {
      if( pastEvents[i].isSpeciation == 1)
	counter++;
    }

  
  if(globalVar.maxEventsSize < globalVar.currentEvents + counter )
    resizePastEvents( counter );
  
  for( i = 0; i < pastEvents -> nevents; ++i)
    {

      if( pastEvents[i].isSpeciation == 1)
	
	{
	  
	  ++globalVar.currentEvents;
	  
	  pastEvents[ globalVar.currentEvents - 1].detype = 'i';
	  
	  pastEvents[ globalVar.currentEvents - 1].paramv = 0.;
	  
	  pastEvents[ globalVar.currentEvents - 1].isSpeciation = 0;

	  pastEvents[ globalVar.currentEvents - 1].time = (pastEvents[ i ].time - duration ) > 0 ? pastEvents[ i ].time - duration : MINTIME;
	  	    
	}
    }
  
  pastEvents ->nevents = globalVar.currentEvents;

  return 1;
}


void freeList(struct devent *pt)
{
  struct devent *tmp;

  while (pt != NULL)
    {
      tmp = pt;
      pt = pt->nextde;
      free(tmp);
    }

}

void updateTreeFromFile(int *positionsArray)
{
  int i = 0, j;
  

  for( i = 0; i < tree.nnode; ++i)
    {
      
      nodes[i].offset = initialNodes[i].offset; 
      
    }


  
  for( i = 0; i < tree.nnode; ++i)
    {
      
      nodes[i].npop = initialNodes[ i ].npop;

    }
  
  for( i = 0; i < tree.nnode; ++i)
    {
      
      nodes[i].Ne = initialNodes[i].Ne; //tempDouble[ ind ];
    }

  
  for( i = 0; i < tree.nnode; ++i)
    {
      nodes[i].sampleSize = initialNodes[i].sampleSize; //tempInt[ ind ];
    }

  
  for( i = 0; i < tree.nnode; ++i)
    {
      nodes[i].nseqs = initialNodes[i].nseqs;
  
      if( i < tree.leaves)
	nodes[i].isLeaf = 1;
      else
	nodes[i].isLeaf =0;
    }
  
  for( i = 0; i < tree.nnode; ++i)
    {

      nodes[i].samplePops = realloc( nodes[i].samplePops, ( nodes[i].npop) * sizeof(int) );
      
      for( j = 0; j < nodes[i].npop; ++j)
	nodes[i].samplePops[j] = initialNodes[i].samplePops[j]; 
    }
}


void getPhyloTreeANDJointEvents( FILE *phyloInputFile, int ntaxa, double birth, double death, double sample, double mut, int option, char *outputTree, int mode, double torigin, FILE *phyloTreeOut, double theta, double msites, double duration, double maxMigration, double oldestOrigin, int *replicationsPerTree, double *speciesSamplingTime)
{
  
  /* Either generate a random phylogenetic tree or read a tree from a file */
  /* These trees describe the speciation history */
  assert(nodes != NULL);

  static int counter = 0;
  
  if(counter == 0 && phyloInputFile != NULL)
    setInitialNodes();

  int *positionsArray;

  *replicationsPerTree = 0;

  counter++;

     
  
  if(phyloInputFile != NULL)
    {

      if(counter == 1)
	fprintf(stderr, "Tree file will be read...\n");
      
      positionsArray = calloc( tree.nnode, sizeof(int) );
      
      *replicationsPerTree = readNumberOfReplications(phyloInputFile); 

      
      if( readTreeFromFile(phyloInputFile, tree.leaves, positionsArray) == 1 )
	{
	  updateTreeFromFile(positionsArray);
	}
      else
	{
	  fprintf(stderr, "\n\nWarning... Unexpected EOF was found.\n\nPerhaps more trees are needed\n\n");

	  assert(0);
	  
	}
      
      

      if(phyloTreeOut != NULL)
	{
	  
	  if(phyloTreeOut != NULL && com.ns<20) 
	    {
	      OutTreeN(phyloTreeOut,0,1);
	    }
	  else if(phyloTreeOut != NULL)
	    {
	      OutTreeN(phyloTreeOut,1,1);
	    }
	}

      free(positionsArray);
      
    }
  
  else
    {
      
      fprintf(stderr, "Tree will be simulated....\n");

      randomTreeGeneration(tree.leaves, birth, death, sample, mut, 2, outputTree, mode, torigin, phyloTreeOut, oldestOrigin, speciesSamplingTime);

    }


    
  setInternalFlag();
  
  /* set Partial Migration values and rescale time */
  nodes -> partialDuration = duration * (msites / theta);

  /* fprintf(stderr, "isoPeriod: %f\n", nodes->partialDuration); */
  
  nodes -> partialMaxMigration = maxMigration;
    
  /* based on the phylogenetic tree and the command line,
     update the populations that are present on the tree
  */
  updateTreePops(nodes);
  
  if(counter == 1 && phyloInputFile != NULL)
    updateTreePops(initialNodes);
  
  /* int i; */
  /* for(i = 0; i < tree.nnode; ++i) */
  /*   printf("-repl: %d, node %d, time: %f, npop: %d, samplesize: %d\n", count, i, nodes[i].age, nodes[i].npop, nodes[i].sampleSize); */

  
  /* transform phylogenetic time (expected number of substitutions)
     to coalescent time units 
  */
  speciationToCoalescentTimes(theta, (double)msites); 
  
  /* turn nodes -- speciation events -- to coalescent events */
  nodesToCoalEvents();
  
  /* if partial Isolation, then set the partial IsolationEvents */
  
  pastEvents -> nevents = globalVar.currentEvents;  
    
  assert( pastEvents -> nevents > 0);
  
  if( globalVar.partialIsolationModel == 1 && pastEvents -> nevents > 0)
    {
      setPartialIsolationEvents( pastEvents, nodes -> partialDuration );
      
      speciationIsolationTimes();
      
    }
    
  /* sort the past events from the most recent to the oldest */
  qsort( pastEvents, globalVar.currentEvents, sizeof( struct deventArray ), compare_struct ); 
  
  pastEvents -> nevents = globalVar.currentEvents;

  /* printf("pastEvents: %d\n", pastEvents->nevents); */
  
  /* from the sorted array of events create a list of events
     this list can now be used from the core of ms 
  */
  createEventListFromEventArray(pastEvents, pars.cp.deventlist);          
  
}

int copyCMDEventsToPastEvents(int nevents, int *cmdEvents)
{
  int i;

  qsort( pastEvents, globalVar.currentEvents, sizeof( struct deventArray ), compare_struct_type); 

  pastEvents->nevents = nevents;

  for( i =0; i < nevents; ++i)
    {
      if(pastEvents[i].type == 1)
	break;
    }
  
  *cmdEvents = i;
      
  return 1;
}

void getSpeciesSamplingTime( int leaves, int totalPops, double *popsSamplingTime, double *speciesSamplingTime )
{
  int i,j, pop;

  double maxSamplingTime = 0.;

  for( i = 0; i < leaves; ++i )
    {
      
      maxSamplingTime = 0.;
      
      for( j = 0; j < nodes[i].npop; ++j)
	{
	  pop = nodes[i].pops[j];

	  assert( pop < totalPops);

	  if( maxSamplingTime < popsSamplingTime[ pop ] )
	    maxSamplingTime = popsSamplingTime[ pop ];
	  
	}

      speciesSamplingTime[ i ] = maxSamplingTime;

      /* printf("Sampling time for species %d is %f\n", i, maxSamplingTime); */
    }
}

void setLeavesPops()
{
  
  int i = 0, j = 0, k = 0;
  
  for( i = 0; i < tree.leaves; ++i)
    {
      
      /* printf("i: %d, NPOP: %d\n", i, nodes[i].npop); */
      
      if(nodes[i].pops == NULL )
	nodes[i].pops = calloc( nodes[i].npop, sizeof(int) );
      else
	{
	  free(nodes[i].pops);
	  nodes[i].pops = calloc( nodes[i].npop , sizeof(int) );
	}

      assert( nodes[i].pops != NULL);

      /* printf("node %d has : ", i); */
      
      for( k = 0; k < nodes[i].npop; ++k)
	{
	  nodes[i].pops[k] = j;
	  j += 1;
	  
	  /* printf(" %d ", nodes[i].pops[k] ); */
	 
	}
      
      /* printf("\n"); */
      
    }

}

int main(int argc, char **argv)
{
    
  int i, k, howmany, segsites, nsegs, seg, maxSegs = 0, *partitionLengths = NULL, treeGenerationMode=3, totalPastEvents = 0, cmdEvents = 0, replicationsPerInputTree = 0, coalPerPhylo = 0, CumulativecoalPerPhylo = 0, unsupportedArgument = -1; 
  
  unsigned int timeSeed; 
  
  char **list = NULL, 
    **cmatrix(), 
    **tbsparamstrs, 
    name[MAX_NAME_LEN], 
    InfoFileName[FILESIZE], 
    ResultFileName[FILESIZE], 
    phyloFileName[FILESIZE], 
    CoalescentFileName[FILESIZE],
    phyloInputFileName[FILESIZE];

  name[0] = '\0';
  
  FILE *pf, 
    *fopen(), 
    *phyloFile = NULL, *nfo = stderr, 
    *InfoFile = NULL, 
    *ResultFile=NULL, 
    *CoalescentFile = NULL, 
    *phyloInputFile = NULL ;

   /* the main structure for the total tree: coalescent + phylogenetics */
  TTree **treeSet = NULL;
    
  double *partitionRates = NULL, 
    partMaxMigration =0,
    partIsoPeriod = 0., 
    probss, 
    tmrca, 
    ttot,  
    time0 = gettime(), 
    mut = 1.0, /* mutation rate on phylogenetic tree*/
    samplingRate = 1.0, /* sampling rate */
    lambda = 100.0, /* birth rate */
    mu = 0.0, /* death rate */
    torigin = 1.0, /* time of origin of process */
    oldestOrigin = -1., 
    scale = 1.,
    *speciesSamplingTime = NULL;

  struct timeval tv;

  static struct segl *seglst = NULL;
  
  printHeader(stderr);
  
  initialize_global_values();
  
  /* seeding using the microseconds */
  gettimeofday( &tv, NULL);

  timeSeed = tv.tv_usec * tv.tv_sec;

  SetSeed(timeSeed);

  
  ntbs = 0 ;   /* these next few lines are for reading in parameters from a file (for each sample) */
  tbsparamstrs = (char **)malloc( argc*sizeof(char *) ) ;

  for( i=0; i<argc; i++) fprintf(stderr, "%s ",argv[i]);
  
  fprintf(stderr, "\n");
  
  for( i =0; i<argc; i++) tbsparamstrs[i] = (char *)malloc(30*sizeof(char) ) ;
  
  for( i = 1; i<argc ; i++)
    if( strcmp( argv[i],"tbs") == 0 )  argv[i] = tbsparamstrs[ ntbs++] ;

  
  tree.leaves = 5;

  /* get the phylogenetic tree */
  char outputTree[MAXTREESTRING];
  
  pastEvents = calloc( globalVar.maxEventsSize, sizeof( struct deventArray ) );
  
  count=0;
  
  if( ntbs > 0 )  
    for( k=0; k<ntbs; k++)  
      if(scanf(" %s", tbsparamstrs[k] ) != 1)
	{
	  fprintf(stderr, "ERROR in reading the tbsparamstrs\n");
	  assert(0);
	}
  

  /* read mutation model parameters */
  if ( getparsMutationModel(argc, argv, name) == 0 )
    {
      globalVar.finiteSiteModel = 0;
    }
  else
    {
      globalVar.finiteSiteModel = 1;
      
      if (rateHetero == CodonRates && invariableSites) {
	fprintf(stderr, "Invariable sites model cannot be used with codon rate heterogeneity.\n");
	exit(0);
      }
      
    }

  if(globalVar.finiteSiteModel == 1 )
    
    {

      checkRequiredArguments(name);
      
      sprintf(InfoFileName, "comus_Info.%s", name);
      
      InfoFile = fopen(InfoFileName, "w");

      nfo = InfoFile;
      
      
      
            
      sprintf(phyloFileName, "comus_Phylo.%s", name);
      
      sprintf(CoalescentFileName, "comus_Coalescent.%s", name);
      
      printHeader(InfoFile);

    }

    
  if(globalVar.finiteSiteModel == 0)
    {
      fprintf(nfo, "Infinite site model -- ");
    }
  else
    {      
      fprintf(stderr, "Finite site model -- ");
    }
  
  /*********************************/

  

  if( argc > 1 && checkAllCommandLine(argc, argv, &unsupportedArgument) == 0)
    {
      
      fprintf(stderr, "\n\nError!\n\n");
      fprintf(stderr, "Argument %s is unsupported. Check usage by typing ./comus -h or comus -h \n", argv[unsupportedArgument]);
      exit(1);
    }
  
  if( argc > 3 && argv[3][0] != '-')
    {

      globalVar.multipleSpecies = 1;
      
      fprintf(nfo, "Multiple Species model\n\n");
      
      fprintf(stderr, "Multiple Species model\n\n");

      getparsCMD_multiSpeciesSimulate( argc, argv, &howmany, &phyloFile, phyloFileName, &CoalescentFile, CoalescentFileName,  &partIsoPeriod, &partMaxMigration, &samplingRate, &mu, &lambda, &mut, phyloInputFileName, &phyloInputFile, &treeGenerationMode, &torigin, &oldestOrigin, &timeSeed);

        
      SetSeed( timeSeed );
      
      assert ( treeGenerationMode < 5 && treeGenerationMode > -1);

      assert ( mu >= 0. );

      assert ( lambda >= 0.);

      assert ( mut >= 0. );

      assert ( oldestOrigin == -1. || oldestOrigin > 0 );

      assert( howmany > 0 );

      /* allocate space for the pars.cp parameters 
	 These parameters are used in the coalescent itself 
      */
      allocateParValues();
      
      /* 
	 set pars.cp values from initial_pars.cp 
      */
      setParsFromInitialValues();

      setLeavesPops();

      speciesSamplingTime = calloc( tree.leaves, sizeof(double) );

      if( globalVar.readSampling == 1)
	{
	  getSpeciesSamplingTime( tree.leaves, globalVar.totalPops,  pars.cp.samplingTime, speciesSamplingTime);

	  /* rescale sampling times */
	  double scale = pars.cp.msites / pars.mp.theta;
	  
	  for( i = 0; i < pars.cp.npop; ++i)
	    pars.cp.samplingTime[i] *= scale;
	}

      /* for( i = 0; i < globalVar.totalPops; ++i) */
      /* 	{ */
      /* 	  printf("samplng time for population %d is %f\n", i, pars.cp.samplingTime[i]); */
      /* 	} */


      if(oldestOrigin > 0 && treeGenerationMode != 4 )
	{
	  fprintf(stderr, "WARNING! If you specify the oldestOrigin parameter you should choose mode == 4. Tree Generation mode becomes 4 automatically.\n\n");

	  treeGenerationMode = 4;
	  
	}

      
      if( globalVar.printInSeparateFiles == 0 )
	{
	  sprintf(ResultFileName, "comus_Results.%s", name);
	  ResultFile = fopen(ResultFileName, "w");
	}
      else 
	ResultFile = NULL;

      
      if(tree.leaves > 1 )
	{
	  
	  getPhyloTreeANDJointEvents( phyloInputFile, tree.leaves, lambda, mu, samplingRate, mut, 2, outputTree, treeGenerationMode, torigin, phyloFile, pars.mp.theta,  (double)pars.cp.msites, partIsoPeriod, partMaxMigration, oldestOrigin, &replicationsPerInputTree, speciesSamplingTime );
            
	  totalPastEvents = pastEvents->nevents;

	  
	  if( checkTreeSamplingConsistency(speciesSamplingTime, tree.leaves) == 0)
	    {
	      assert(0);
	    }
  
	}

      
    }
  else
    {
      
      if(globalVar.finiteSiteModel == 1)
	{
	  fprintf(nfo, "Error please specify the parameters to make make species, even if you have just a single species\n");
	  assert( globalVar.finiteSiteModel = 0);
	}

      fprintf(nfo, "Single Species model\n\n");
      
      globalVar.multipleSpecies = 0;
      
      getpars( argc, argv, &howmany) ;   /* results are stored in global variable, pars */
    }

  
  for( i=0; i<argc; i++) 
    fprintf(nfo, "%s ",argv[i]);
  
  fprintf(nfo, "\n");
  
  
  if( !pars.commandlineseedflag ) 
    seedit( "s");
  
  
  if( ResultFile != NULL && globalVar.printInSeparateFiles == 0 )
    pf = ResultFile;
  else if( globalVar.printInSeparateFiles == 0 )
    pf = stdout;

  if( pars.mp.segsitesin ==  0 ) {
    list = cmatrix(pars.cp.nsam, maxsites+1);
    posit = (double *)malloc( (unsigned)( maxsites*sizeof( double)) ) ;
  }
  else {
    list = cmatrix(pars.cp.nsam, pars.mp.segsitesin+1 ) ;
    posit = (double *)malloc( (unsigned)( pars.mp.segsitesin*sizeof( double)) ) ;
    if( pars.mp.theta > 0.0 ){
      segfac = 1.0 ;
      for(  i= pars.mp.segsitesin; i > 1; i--) segfac *= i ;
    }
  }
  
  
  /**** finite mutation model *****/
  
  if(globalVar.finiteSiteModel == 1)
    {
      SetModel(model);
      numSites = pars.cp.msites = initial_pars.cp.msites;
    }

  /*******************************/

  if(phyloInputFile != NULL)
    coalPerPhylo = 0;

  if( globalVar.multipleSpecies == 1)
    {
      CumulativecoalPerPhylo = replicationsPerInputTree;
      
      treeSet = malloc( (pars.cp.nsites + 10  ) * sizeof(TTree*) ); /* TODO debugging 20140113 */
    }


  fprintf(nfo, "\ntimeSeed:%u\n", timeSeed);


  mutationTimesFile = fopen("mutationTimes.txt", "w");
  
  while( howmany-count++ ) {

    fprintf(mutationTimesFile, "//\n");
    
    /* for(i = 0; i < pars.cp.npop; ++i) */
    /*   printf("COUNT: %d, pop %d has %d samples\n", count, i, pars.cp.config[i]); */
    
    
    if( tree.leaves > 1 && (coalPerPhylo > 0) && ((count-1) % CumulativecoalPerPhylo) == 0)
      {

	copyCMDEventsToPastEvents(totalPastEvents, &cmdEvents);

	globalVar.currentEvents = cmdEvents;
	
	/*
	  set pars.cp values from initial_pars.cp
	*/
	setParsFromInitialValues();
	
	/* for(i = 0; i < pars.cp.npop; ++i) */
	/*   printf("1.COUNT: %d, pop %d has %d samples\n", count, i, pars.cp.config[i]); */
    	
	if(oldestOrigin > 0 && treeGenerationMode != 4 )
	  {
	    fprintf(stderr, "WARNING! If you specify the oldestOrigin parameter you should choose mode == 4. Tree Generation mode becomes 4 automatically.\n\n");
	    
	    treeGenerationMode = 4;
	    
	  }
	
	getPhyloTreeANDJointEvents( phyloInputFile, tree.leaves, lambda, mu, samplingRate, mut, 2, outputTree, treeGenerationMode, torigin, phyloFile, pars.mp.theta,  (double)pars.cp.msites, partIsoPeriod, partMaxMigration, oldestOrigin, &replicationsPerInputTree, speciesSamplingTime);

	CumulativecoalPerPhylo += replicationsPerInputTree;

      }
    
    /* for(i = 0; i < pars.cp.npop; ++i) */
    /*   printf("COUNT: %d, pop %d has %d samples\n", count, i, pars.cp.config[i]); */
    

    if(replicationsPerInputTree > 0)
      {
	coalPerPhylo = replicationsPerInputTree;	
      }
    else
      coalPerPhylo = 0;

    
    /* for(i = 0; i < tree.nnode; ++i) */
    /*   printf("repl: %d, node %d, time: %f, npop: %d, samplesize: %d\n", count, i, nodes[i].age, nodes[i].npop, nodes[i].sampleSize); */
    
    
    /* for( i =0; i < totalPastEvents; ++i) */
    /*   printf("events: %c, type: %d, time: %f\n", pastEvents[i].detype, pastEvents[i].type, pastEvents[i].time); */

    
    if( (ntbs > 0) && (count >1 ) ){
      for( k=0; k<ntbs; k++){ 
	if( scanf(" %s", tbsparamstrs[k]) == EOF ){
	  if( !pars.commandlineseedflag ) seedit( "end" );
	  exit(0);
	}
      }
      
      getpars( argc, argv, &howmany) ;
      
    }

    if(globalVar.printInSeparateFiles == 1)
      {
	
	if(ResultFile != NULL)
	  {
	    fclose(ResultFile);
	    ResultFile =  NULL;
	  }
	
	sprintf(ResultFileName, "comus_Results_alignment_%d.%s", count, name);
      
	ResultFile = fopen(ResultFileName, "w");

	pf = ResultFile;
	
      }
    
    if( globalVar.printSeparator )
      fprintf(pf,"\n//\n");
    else if(count > 1)
      fprintf(pf, "\n");
    
    if( ntbs >0 ){
      for(k=0; k< ntbs; k++) fprintf(ResultFile, "\t%s", tbsparamstrs[k] ) ;
    }
        
    segsites = gensam( list, &probss, &tmrca, &ttot, &seglst, ResultFile, &nsegs, InfoFile, count) ; 

    if(globalVar.finiteSiteModel == 1)
      {
	
	if(nsegs > maxSegs )
	  {
	    
	    if(maxSegs == 0)
	      {
		
		partitionLengths = malloc(sizeof(int) * nsegs);  
		
		partitionRates = (double *)malloc(sizeof(double) * nsegs);

		
		for (i = 0; i < nsegs; i++) {
		  
		  if ( (treeSet[i]=NewTree())==NULL ) {
		    fprintf(stderr, "Out of memory\n");
		    exit(0);
		  }
		  
		}

	      }
	    else
	      {
		treeSet = realloc( treeSet, nsegs * sizeof(TTree*) );

		partitionLengths = realloc( partitionLengths, nsegs * sizeof(int) );
		partitionRates = realloc(partitionRates, sizeof(double) * nsegs);
				
		for (i = maxSegs; i < nsegs; i++) {
		  
		  if ( (treeSet[i]=NewTree())==NULL ) {
		    fprintf(stderr, "Out of memory\n");
		    exit(0);
		  }
		  
		}

	      }

	    
	    if (partitionLengths==NULL) {
	      fprintf(stderr, "Out of memory\n");
	      exit(0);
	    }
	    
	    if (partitionRates==NULL) {
	      fprintf(stderr, "Out of memory\n");
	      exit(0);
	    }
	          
	    if (treeSet==NULL) {
	      fprintf(stderr, "Out of memory\n");
	      exit(0);
	    }

	    maxSegs = nsegs;
	  }
	
	getSequencesFiniteModel( seglst, nsegs, pars.cp.nsam, ResultFile, treeSet, partitionLengths, partitionRates);
	

      }
    else if(globalVar.finiteSiteModel == 0)
      {
    
	if( pars.mp.timeflag ) fprintf(pf,"time:\t%lf\t%lf\n",tmrca, ttot ) ;
	if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
	  if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 )) 
	    fprintf(pf,"prob: %g\n", probss ) ;
	  fprintf(pf,"segsites: %d\n",segsites);
	  if( segsites > 0 )	fprintf(pf,"positions: ");
	  for( i=0; i<segsites; i++)
	    fprintf(pf,"%6.*lf ", pars.output_precision,posit[i] );
	  fprintf(pf,"\n");
	  if( segsites > 0 ) 
	    for(i=0;i<pars.cp.nsam; i++) { fprintf(pf,"%s\n", list[i] ); }
	}
      }
    
    if( CoalescentFile != NULL)
      {
        if(globalVar.finiteSiteModel == 1 && globalVar.multipleSpecies == 1)
          scale = pars.mp.theta/pars.cp.msites;
	
	printCoalescentTrees(CoalescentFile, seglst, nsegs, pars.cp.nsam, scale);
      }
   
    for(seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, ++k)
      {
	free(seglst[seg].ptree);
      }

    
    double time1 = gettime();

    printIncrementer(stderr, count, howmany, time1 - time0);
    
    printIncrementer(nfo, count, howmany, time1 - time0);
 
  }

  if(mutationTimesFile != NULL)
    fclose(mutationTimesFile);
  
  
  if( !pars.commandlineseedflag ) seedit( "end" );
  
  if(phyloInputFile != NULL)
    fclose(phyloInputFile);

  /* free tbsparamstrs */  
  for( i =0; i<argc; i++) 
    free(tbsparamstrs[i]);
  free(tbsparamstrs);
  
  /* free the list of SNPs */
  if(list != NULL)
    {
      for( i = 0; i < pars.cp.nsam; ++i)
	free(list[i]);
      free(list);
    }
  
  //freeGlobalEvolve();
  if(globalVar.finiteSiteModel == 1)
    free_global();

  
  if(pastEvents != NULL)
    free(pastEvents);


  if(phyloFile != NULL)
    fclose(phyloFile);

  if(CoalescentFile != NULL)
    fclose(CoalescentFile);

  if(InfoFile != NULL)
    fclose(InfoFile);

  freeList(pars.cp.deventlist);

  if(globalVar.finiteSiteModel == 1)
    {
      
      if(partitionLengths != NULL)
	free( partitionLengths );
      
      if(partitionRates != NULL)
	free( partitionRates );
    }
  
  printf("Total time: %lf\n", gettime() - time0);
    
  return 1;
  
}



int mnmial(n,nclass,p,rv)
     int n, nclass, rv[];
     double p[];
{
  double x, s;
  int i, j = 0;

  for(i=0; i<nclass; i++) 
    rv[i]=0;
  
  for(i=0; i<n ; i++) {
    x = rndu();
    
    j=0;
    
    s = p[0];
    
    while( (x > s) && ( j<(nclass-1) ) )  
      s += p[++j];
    
    rv[j]++;
    
  }
  
  return(j);
}




/**** ordran.c  ***/
int ordran(int n, double pbuf[])
{
  ranvec(n,pbuf);
  order(n,pbuf);
  return 0;
}



int locate(int n, double beg, double len, double *ptr)
{
  int i;
  ordran(n,ptr);
  for(i=0; i<n; i++)
    ptr[i] = beg + ptr[i]*len ;

  return 0;

}


int biggerlist(int nsam,  char **list )
{
  int i;

  /*  fprintf(stderr,"maxsites: %d\n",maxsites);  */	
  for( i=0; i<nsam; i++){
    list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
    if( list[i] == NULL ) perror( "realloc error. bigger");
  }

  return 0;
}


int poisso(u)
     double u;
{
  double  cump, ru,  p, gasdev(double, double) ;
  int i=1;

  if( u > 30. ){
    i =  (int)(0.5 + gasdev(u,u)) ;
    if( i < 0 ) return( 0 ) ;
    else return( i ) ;
  }
	 
  ru = rndu();
  p = exp(-u);
  if( ru < p) return(0);
  cump = p;
	
  while( ru > ( cump += (p *= u/i ) ) )
    i++;
  return(i);
}



int  gensam( char **list, double *pprobss, double *ptmrca, double *pttot, struct segl **seglst, FILE *ResultFile, int *numberofSegments, FILE* infofile, int simindex) 
{
  int nsegs, i, k, seg, ns = 0, start, end, len, segsit ;
  
  struct segl *segtre_mig(struct c_params *p, int *nsegs, FILE *infofile, int simindex ), 
    *segtre_mig_multSpecies(struct c_params *p, int *nsegs, struct TREEN *nodes, struct TREEB *tree); /* used to be: [MAXSEG];  */

  double nsinv,  tseg, tt, ttime(struct node *, int nsam), ttimemf(struct node *, int nsam, int mfreq) ;
  double *pk;
  int *ss;
  int segsitesin,nsites;
  double theta, es ;
  int nsam, mfreq ;
  void prtree(FILE *fout,  struct node *ptree, int nsam, double scale);
  int make_gametes(int nsam, int mfreq,  struct node *ptree, double tt, int newsites, int ns, char **list );
  void ndes_setup( struct node *, int nsam );

  nsites = pars.cp.nsites ;
  
  nsinv = 1./nsites;
  
  if(globalVar.partialIsolationModel == 1)
    (*seglst) = segtre_mig_multSpecies( &(pars.cp), &nsegs, nodes, &tree);
  else
    (*seglst) = segtre_mig(&(pars.cp),  &nsegs, infofile, simindex ) ;

  *numberofSegments = nsegs;
	
  nsam = pars.cp.nsam;

  segsitesin = pars.mp.segsitesin ;

  theta = pars.mp.theta ;

  mfreq = pars.mp.mfreq ;


  if( pars.mp.timeflag ) {
    tt = 0.0 ;
    for( seg=0, k=0; k<nsegs; seg=(*seglst)[seg].next, k++) { 
      if( mfreq > 1 ) ndes_setup( (*seglst)[seg].ptree, nsam );
      end = ( k<nsegs-1 ? (*seglst)[(*seglst)[seg].next].beg -1 : nsites-1 );
      start = (*seglst)[seg].beg ;
      if( (nsegs==1) || ( ( start <= nsites/2) && ( end >= nsites/2 ) ) )
	*ptmrca = ( (*seglst)[seg].ptree + 2*nsam-2) -> time ;
      len = end - start + 1 ;
      tseg = len/(double)nsites ;
      if( mfreq == 1 ) tt += ttime( (*seglst)[seg].ptree,nsam)*tseg ;
      else tt += ttimemf( (*seglst)[seg].ptree,nsam, mfreq)*tseg ;
      
      /* if( (segsitesin == 0) && ( theta == 0.0 )  )  */
      /* 	free( (*seglst)[seg].ptree) ; */
    }
    *pttot = tt ;
  }	

  /* printf("theta: %f\n", theta); */
	
  if( globalVar.finiteSiteModel == 0 && (segsitesin == 0) && ( theta > 0.0)   ) {
    ns = 0 ;


    
    for( seg=0, k=0; k<nsegs; seg= (*seglst)[seg].next, k++) { 
      if( mfreq > 1 ) ndes_setup( (*seglst)[seg].ptree, nsam );
      end = ( k<nsegs-1 ? (*seglst)[(*seglst)[seg].next].beg -1 : nsites-1 );
      start = (*seglst)[seg].beg ;
      len = end - start + 1 ;
      tseg = len*(theta/nsites) ;
      if( mfreq == 1) tt = ttime((*seglst)[seg].ptree, nsam);
      else tt = ttimemf( (*seglst)[seg].ptree, nsam, mfreq );
      segsit = poisso( tseg*tt );

      /* printf("tt: %f, segsit: %d, tseg: %f, theta: %f, nsites: %d\n", tt, segsit, tseg, theta, nsites); */
      /* exit(-1); */

      if( (segsit + ns) >= maxsites ) {
	maxsites = segsit + ns + SITESINC ;
	posit = (double *)realloc(posit, maxsites*sizeof(double) ) ;
	biggerlist(nsam, list) ; 
      }

      make_gametes(nsam,mfreq, (*seglst)[seg].ptree,tt, segsit, ns, list );
      
      //free( seglst[seg].ptree) ;
      
      locate(segsit,start*nsinv, len*nsinv,posit+ns);   
      
      ns += segsit;
    }
    fprintf(mutationTimesFile, "\n");
  }
  else if( globalVar.finiteSiteModel == 0 && segsitesin > 0 ) {

    pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
    ss = (int *)malloc((unsigned)(nsegs*sizeof(int)));
    if( (pk==NULL) || (ss==NULL) ) perror("malloc error. gensam.2");


    tt = 0.0 ;
    for( seg=0, k=0; k<nsegs; seg=(*seglst)[seg].next, k++) { 
      if( mfreq > 1 ) ndes_setup( (*seglst)[seg].ptree, nsam );
      end = ( k<nsegs-1 ? (*seglst)[(*seglst)[seg].next].beg -1 : nsites-1 );
      start = (*seglst)[seg].beg ;
      len = end - start + 1 ;
      tseg = len/(double)nsites ;
      if( mfreq == 1 ) pk[k] = ttime( (*seglst)[seg].ptree,nsam)*tseg ;
      else pk[k] = ttimemf( (*seglst)[seg].ptree,nsam, mfreq)*tseg ;
      tt += pk[k] ;
    }
    if( theta > 0.0 ) { 
      es = theta * tt ;
      *pprobss = exp( -es )*pow( es, (double) segsitesin) / segfac ;
    }
    if( tt > 0.0 ) {
      for (k=0;k<nsegs;k++) pk[k] /= tt ;
      mnmial(segsitesin,nsegs,pk,ss);
    }
    else
      for( k=0; k<nsegs; k++) ss[k] = 0 ;
    ns = 0 ;
    for( seg=0, k=0; k<nsegs; seg=(*seglst)[seg].next, k++) { 
      end = ( k<nsegs-1 ? (*seglst)[(*seglst)[seg].next].beg -1 : nsites-1 );
      start = (*seglst)[seg].beg ;
      len = end - start + 1 ;
      tseg = len/(double)nsites;
      make_gametes(nsam,mfreq, (*seglst)[seg].ptree,tt*pk[k]/tseg, ss[k], ns, list);

      //free(seglst[seg].ptree) ;
      locate(ss[k],start*nsinv, len*nsinv,posit+ns);   
      ns += ss[k] ;
    }
    free(pk);
    free(ss);

  }

  /* printf("nsam: %d, segsites: %d\n", nsam, ns); */

  /* exit(-1); */
  if(globalVar.finiteSiteModel == 0)
    for(i=0;i<nsam;i++) list[i][ns] = '\0' ;
  
  return( ns ) ;
}

void 
ndes_setup(struct node *ptree, int nsam )
{
  int i ;

  for( i=0; i<nsam; i++) (ptree+i)->ndes = 1 ;
  for( i = nsam; i< 2*nsam -1; i++) (ptree+i)->ndes = 0 ;
  for( i= 0; i< 2*nsam -2 ; i++)  (ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes ;

}

	   


/* allocates space for gametes (character strings) */
char **cmatrix(int nsam, int len)
{
  
  int i;
  char **m;
  
  if( ! ( m = (char **) malloc( (unsigned) nsam*sizeof( char* ) ) ) )
    perror("alloc error in cmatrix") ;
  for( i=0; i<nsam; i++) {
    if( ! ( m[i] = (char *) malloc( (unsigned) len*sizeof( char ) )))
      perror("alloc error in cmatric. 2");
  }
  return( m );
}



int NSEEDS = 3 ;


void
getpars(int argc, char *argv[], int *phowmany )
{
  int arg, i, j, sum , pop , argstart, npop , npop2, pop2 ;
  double migr, mij, psize, palpha ;
  void addtoelist( struct devent *pt, struct devent *elist ); 
  void argcheck( int arg, int argc, char ** ) ;
  int commandlineseed( char ** ) ;
  void free_eventlist( struct devent *pt, int npop );
  struct devent *ptemp , *pt ;
  FILE *pf ;
  char ch3 ;
	

  if( count == 0 ) {
    if( argc < 4 ){ fprintf(stderr,"Too few command line arguments.\n\nUse ./comus -h for a full list of options\n\n");
      exit(0);
    }
    pars.cp.nsam = atoi( argv[1] );
    if( pars.cp.nsam <= 0 ) { fprintf(stderr,"First argument error. nsam <= 0. \n"); assert( pars.cp.nsam > 0 ); }
    *phowmany = atoi( argv[2] );
    if( *phowmany  <= 0 ) { fprintf(stderr,"Second argument error. howmany <= 0. \n"); assert( *phowmany > 0 ); }
    pars.commandlineseedflag = 0 ;
    pars.output_precision = 4 ;
    pars.cp.r = pars.mp.theta =  pars.cp.f = 0.0 ;
    pars.cp.track_len = 0. ;
    pars.cp.npop = npop = 1 ;
    pars.cp.mig_mat = (double **)malloc( (unsigned) sizeof( double *) );
    pars.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double ));
    pars.cp.mig_mat[0][0] =  0.0 ;
    pars.mp.segsitesin = 0 ;
    pars.mp.treeflag = 0 ;
    pars.mp.timeflag = 0 ;
    pars.mp.mfreq = 1 ;
    
    pars.cp.config = (int *) malloc( (unsigned)(( pars.cp.npop +1 ) *sizeof( int)) );
    pars.cp.samplingTime = calloc( pars.cp.npop +1, sizeof( double) );
    
    
    (pars.cp.config)[0] = pars.cp.nsam ;
    
    pars.cp.size= (double *) malloc( (unsigned)( pars.cp.npop * sizeof( double )) );
    (pars.cp.size)[0] = 1.0  ;
    pars.cp.alphag = (double *) malloc( (unsigned)(( pars.cp.npop ) *sizeof( double )) );
    (pars.cp.alphag)[0] = 0.0  ;
    pars.cp.nsites = 2 ;
  }
  else{
    npop = pars.cp.npop ;
    free_eventlist( pars.cp.deventlist, npop );
  }
  pars.cp.deventlist = NULL ;

  arg = 3 ;

  while( arg < argc ){
    if( argv[arg][0] != '-' ) { fprintf(stderr," argument should be -%s ?\n", argv[arg]); assert(argv[arg][0] == '-'); }
    
    switch ( argv[arg][1] ){
   	     
    case 'f' :
      if( ntbs > 0 ) { fprintf(stderr," can't use tbs args and -f option.\n"); exit(1); }
      arg++;
      argcheck( arg, argc, argv);
      pf = fopen( argv[arg], "r" ) ;
      if( pf == NULL ) {fprintf(stderr," no parameter file %s\n", argv[arg] ); exit(0);}
      arg++;
      argc++ ;
      argv = (char **)malloc(  (unsigned)(argc+1)*sizeof( char *) ) ;
      argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
      argstart = arg ;
      while( fscanf(pf," %s", argv[arg]) != EOF ) {
	arg++;
	argc++;
	argv = (char **)realloc( argv, (unsigned)argc*sizeof( char*) ) ;
	argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
      }
      fclose(pf);
      argc--;
      arg = argstart ;
      break;
    case 'r' : 
      arg++;
      argcheck( arg, argc, argv);
      pars.cp.r = atof(  argv[arg++] );
      argcheck( arg, argc, argv);
      pars.cp.nsites = atoi( argv[arg++]);
      if( pars.cp.nsites <2 ){
	fprintf(stderr,"with -r option must specify both rec_rate and nsites>1\n");
	assert( pars.cp.nsites > 1); 
      }
      break;	
    case 'p' :
      arg++;
      argcheck(arg,argc,argv);
      pars.output_precision = atoi( argv[arg++] ) ;
      break;
    case 'c' : 
      arg++;
      argcheck( arg, argc, argv);
      pars.cp.f = atof(  argv[arg++] );
      argcheck( arg, argc, argv);
      pars.cp.track_len = atof( argv[arg++]);
      if( pars.cp.track_len <1. ){
	fprintf(stderr,"with -c option must specify both f and track_len>0\n");
	assert( pars.cp.track_len >= 1. );
      }
      break;		
    case 't' : 
      arg++;
      argcheck( arg, argc, argv);
      pars.mp.theta = atof(  argv[arg++] );
      break;
    case 's' : 
      arg++;
      argcheck( arg, argc, argv);
      if( argv[arg-1][2] == 'e' ){  /* command line seeds */
	pars.commandlineseedflag = 1 ;
	if( count == 0 ) nseeds = commandlineseed(argv+arg );
	arg += nseeds ;
      }
      else {
	pars.mp.segsitesin = atoi(  argv[arg++] );
      }
      break;
    case 'F' : 
      arg++;
      argcheck( arg, argc, argv);
      pars.mp.mfreq = atoi(  argv[arg++] );
      if( (pars.mp.mfreq < 2 ) || (pars.mp.mfreq > pars.cp.nsam/2 ) ){
	fprintf(stderr," mfreq must be >= 2 and <= nsam/2.\n");
	assert( pars.mp.mfreq >= 2);
      }
      break;
    case 'T' : 
      pars.mp.treeflag = 1 ;
      arg++;
      break;
    case 'L' : 
      pars.mp.timeflag = 1 ;
      arg++;
      break;
    case 'I' : 
      arg++;
      if( count == 0 ) {
	argcheck( arg, argc, argv);
	pars.cp.npop = atoi( argv[arg]);
	pars.cp.config = (int *) realloc( pars.cp.config, (unsigned)( pars.cp.npop*sizeof( int)));
	pars.cp.samplingTime = realloc( pars.cp.samplingTime, (unsigned)( pars.cp.npop * sizeof( double) ) );
	npop = pars.cp.npop ;
      }
      arg++;
      for( i=0; i< pars.cp.npop; i++) {
	argcheck( arg, argc, argv);
	pars.cp.config[i] = atoi( argv[arg++]);
      }
      if( count == 0 ){
	
	pars.cp.mig_mat = 
	  (double **)realloc(pars.cp.mig_mat, (unsigned)(pars.cp.npop*sizeof(double *) )) ;
	pars.cp.mig_mat[0] = 
	  (double *)realloc(pars.cp.mig_mat[0], (unsigned)( pars.cp.npop*sizeof(double)));
	for(i=1; i<pars.cp.npop; i++) pars.cp.mig_mat[i] = 
					(double *)malloc( (unsigned)( pars.cp.npop*sizeof(double)));
	pars.cp.size = (double *)realloc( pars.cp.size, (unsigned)( pars.cp.npop*sizeof( double )));
	pars.cp.alphag = 
	  (double *) realloc( pars.cp.alphag, (unsigned)( pars.cp.npop*sizeof( double )));
	for( i=1; i< pars.cp.npop ; i++) {
	  (pars.cp.size)[i] = (pars.cp.size)[0]  ;
	  (pars.cp.alphag)[i] = (pars.cp.alphag)[0] ;
	}
      }
      
      if( (arg <argc) && ( argv[arg][0] != '-' ) ) {
	argcheck( arg, argc, argv);
	migr = atof(  argv[arg++] );
      }
      else migr = 0.0 ;
      for( i=0; i<pars.cp.npop; i++) 
	for( j=0; j<pars.cp.npop; j++) pars.cp.mig_mat[i][j] = migr/(pars.cp.npop-1) ;
      for( i=0; i< pars.cp.npop; i++) pars.cp.mig_mat[i][i] = migr ;
      break;
    case 'm' :
      if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); assert( npop >= 2);}
      if( argv[arg][2] == 'a' ) {
	arg++;
	for( pop = 0; pop <npop; pop++)
	  for( pop2 = 0; pop2 <npop; pop2++){
	    argcheck( arg, argc, argv);
	    pars.cp.mig_mat[pop][pop2]= atof( argv[arg++] ) ;
	  }
	for( pop = 0; pop < npop; pop++) {
	  pars.cp.mig_mat[pop][pop] = 0.0 ;
	  for( pop2 = 0; pop2 < npop; pop2++){
	    if( pop2 != pop ) pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2] ;
	  }
	}	
      }
      else {
	arg++;
	argcheck( arg, argc, argv);
	i = atoi( argv[arg++] ) -1;
	argcheck( arg, argc, argv);
	j = atoi( argv[arg++] ) -1;
	argcheck( arg, argc, argv);
	mij = atof( argv[arg++] );
	pars.cp.mig_mat[i][i] += mij -  pars.cp.mig_mat[i][j]  ;
	pars.cp.mig_mat[i][j] = mij;
      }
      break;
    case 'n' :
      if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); assert( npop >= 2);}
      arg++;
      argcheck( arg, argc, argv);
      pop = atoi( argv[arg++] ) -1;
      argcheck( arg, argc, argv);
      psize = atof( argv[arg++] );
      pars.cp.size[pop] = psize ;
      break;
    case 'g' :
      if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); assert(npop >= 2); }
      arg++;
      argcheck( arg, argc, argv);
      pop = atoi( argv[arg++] ) -1;
      if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); assert( arg < argc);  }
      palpha = atof( argv[arg++] );
      pars.cp.alphag[pop] = palpha ;
      break;
    case 'G' :
      arg++;
      if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); assert( arg < argc);  }
      palpha = atof( argv[arg++] );
      for( i=0; i<pars.cp.npop; i++) 
	pars.cp.alphag[i] = palpha ;
      break;
    case 'e' :
      pt = (struct devent *)malloc( sizeof(struct devent) ) ;
      
      pt->detype = argv[arg][2] ;
      
      ch3 = argv[arg][3] ;
      
      arg++;
      
      argcheck( arg, argc, argv);
      
      pt->time = atof( argv[arg++] ) ;
      
      pt->nextde = NULL ;
      
      if( pars.cp.deventlist == NULL ) 
	pars.cp.deventlist = pt ;
      
      else if ( pt->time < pars.cp.deventlist->time ) { 
	ptemp = pars.cp.deventlist ;
	pars.cp.deventlist = pt ;
	pt->nextde = ptemp ;	
      }	
      else
	addtoelist( pt, pars.cp.deventlist ) ;
      switch( pt->detype ) {
      case 'N' :
	argcheck( arg, argc, argv);
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'G' :
	if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eG.\n"); assert( arg < argc);  }
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'M' :
	argcheck( arg, argc, argv);
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'n' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	argcheck( arg, argc, argv);
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'g' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eg.\n"); assert( arg < argc); }
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 's' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	argcheck( arg, argc, argv);
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'm' :
	if( ch3 == 'a' ) {
	  pt->detype = 'a' ;
	  argcheck( arg, argc, argv);
	  npop2 = atoi( argv[arg++] ) ;
	  pt->mat = (double **)malloc( (unsigned)npop2*sizeof( double *) ) ;
	  for( pop =0; pop <npop2; pop++){
	    (pt->mat)[pop] = (double *)malloc( (unsigned)npop2*sizeof( double) );
	    for( i=0; i<npop2; i++){
	      if( i == pop ) arg++;
	      else {
		argcheck( arg, argc, argv); 
		(pt->mat)[pop][i] = atof( argv[arg++] ) ;
	      }
	    }
	  }
	  
	  for( pop = 0; pop < npop2; pop++) {
	    (pt->mat)[pop][pop] = 0.0 ;
	    for( pop2 = 0; pop2 < npop2; pop2++){
	      if( pop2 != pop ) (pt->mat)[pop][pop] += (pt->mat)[pop][pop2] ;
	    }
	  }
	}
	else {
	  argcheck( arg, argc, argv);
	  pt->popi = atoi( argv[arg++] ) -1 ;
	  argcheck( arg, argc, argv);
	  pt->popj = atoi( argv[arg++] ) -1 ;
	  argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	}
	break;
      case 'j' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	argcheck( arg, argc, argv);
	pt->popj = atoi( argv[arg++] ) -1 ;
	break;
      default: 
	fprintf(stderr,"e event\n");  
	assert(0);
      }
      break;
    default: 
      fprintf(stderr," option default\n");
      assert(0);
    }
  }
  if( (pars.mp.theta == 0.0) && ( pars.mp.segsitesin == 0 ) && ( pars.mp.treeflag == 0 ) && (pars.mp.timeflag == 0) ) {
    fprintf(stderr," either -s or -t or -T option must be used. \n");
    exit(1);
  }
  sum = 0 ;
  for( i=0; i< pars.cp.npop; i++) sum += (pars.cp.config)[i] ;
  if( sum != pars.cp.nsam ) {
    fprintf(stderr," sum sample sizes != nsam\n");
    exit(1);
  }
}


void
argcheck( int arg, int argc, char *argv[] )
{
  if( (arg >= argc ) || ( argv[arg][0] == '-') ) {
    fprintf(stderr,"not enough arguments after %s\n", argv[arg-1] ) ;
    fprintf(stderr,"For usage type: comus -h <return>\n");
    exit(0);
  }
}




/* void constructUpdatedCommandLine() */
/* { */
/*   int i; */
/*   for(i = 0; i < tree.nnode; ++i) */
/*     printf("node: %i, name: %s, id: %d, nson: %d, age: %f, branch: %f\n", i, nodes[i].nodeStr, nodes[i].nodeID, nodes[i].nson, nodes[i].age, nodes[i].branch); */
   
/* } */


void printError(char error[1000])
{
  fprintf(stderr, "\n\n%s\n\n", error);
}

void checkArguments( int arg, int argc, char **argv )
{
  if( arg == argc - 1 )
    {
      fprintf(stderr,  "Argument %s should be followed by a value\n", 
	       argv[ arg ] );

      assert( arg < argc - 1 );
    }

  if( arg < argc -1 )
    {
      if( argv[ arg + 1 ][0] == '-' )
	{
	  fprintf(stderr,  "After argument %s a value is expected. However there is: %s\n", argv[ arg ], argv[ arg + 1 ]);
	  
	  assert( argv[ arg + 1 ][0] != '-' );
	}
    }
  
}


int getTheta(int *arg, int argc, char **argv)
{
  int i;
  
  if(argc == *arg + 2)
    {
      for( i = 0; i < tree.leaves; ++i)
	nodes[i].theta = atof(argv[*arg + 1] ); /* spConf[i].theta =*/
      ++(*arg);
    }
  else if( argc > *arg + 2 && argv[*arg + 2][0] == '-')
    {
      for( i = 0; i < tree.leaves; ++i)
	{
	  assert( argv[*arg + 1][0] != '-' );
	  nodes[i].theta = atof(argv[*arg + 1] ); //  spConf[i].theta =
	}
      ++(*arg);
    }
  else if( argc > *arg + tree.leaves )
    {
      for( i = 0; i < tree.leaves; ++i)
	{
	  assert(argv[*arg+1][0] != '-');
	  nodes[i].theta =  atof(argv[++(*arg)]);
		  
	}
    }
  else
    {
      return 0;
    }
	  
  return 1;
}



int getNe(int *arg, int argc, char **argv)
{
  int i;
  
  if(argc == *arg + 2)
    {
      for( i = 0; i < tree.leaves; ++i)
	nodes[i].Ne = atof(argv[*arg + 1] );
      ++(*arg);
    }
  else if( argc > *arg + 2 && argv[*arg + 2][0] == '-')
    {
      for( i = 0; i < tree.leaves; ++i)
	nodes[i].Ne =  atof(argv[*arg + 1] );
      ++(*arg);
    }
  else if( argc > *arg + tree.leaves )
    {
      for( i = 0; i < tree.leaves; ++i)
	{
	  assert(argv[*arg+1][0] != '-');
	  nodes[i].Ne  = atof(argv[++(*arg)]);
		  
	}
    }
  else
    {
      return 0;
    }
	  
  return 1;
}



int getRho(int *arg, int argc, char **argv)
{
  int i;
  

  /* if rho is the last argument */
  if(argc == *arg + 2)
    {
      
      for( i = 0; i < tree.leaves; ++i)
	{
	  assert( argv[*arg + 1][0] != '-' );
	  nodes[i].rho = atof(argv[*arg + 1] );
	}
            
      (*arg) += 1;
    }
  
  /* if rho is in the middle of the command line and one rho is given */
  else if( argc > *arg + 2 && argv[*arg + 2][0] == '-')
    {
      
      for( i = 0; i < tree.leaves; ++i)
	{
	  assert( argv[*arg + 1][0] != '-' );

	  initial_pars.cp.r = nodes[i].rho =  atof(argv[*arg + 1] );
	}
            
      (*arg) += 1;
    }
  /* if rho is in the middle of the command line and more rho are given */
  else if( argc > *arg + tree.leaves + 1 )
    {
      
      for( i = 0; i < tree.leaves; ++i)
	{
	  assert(argv[*arg+1][0] != '-');
	  
	  nodes[i].rho = atof(argv[++(*arg)]);
		  
	}
      
    }
  else
    {
      return 0;
    }
	  
  return 1;
}


unsigned long getFactorial(int n)
{
  unsigned long  i;
  unsigned long f = 1;
  
  for(i = 2; i <= n; ++i)
    {
      f = f * i;
    }
  
  return f;

}

unsigned long getCombination(int n, int x)
{
  unsigned long numerator = getFactorial(n);

  unsigned long den1 = getFactorial(x);
  
  unsigned long den2 = getFactorial( n - x );

  unsigned long res = numerator / den1;
  
  res /= den2;

  return res;
}

/* int updateMigrationStructure(int nvalues, int dim, int arg, int argc, char **argv) */
/* { */
/*   int i; */
/*   int sp1, sp2, pop1, pop2; */
/*   double migrationRate; */


  

/*   if( nvalues == 2 ) */
/*     { *//*       /\* 1 species, set value for all pop pairs *\/ */
/*       sp1 = atoi( argv[++arg] ); */
  
/*       assert( sp1 > 0 ); */
/*       assert( sp1 < dim + 1); // species number starts from 1 */
    
/*       migrationRate = atof( argv[ ++arg ] ); */
/*       assert( migrationRate > 0.); */

/*       cmdMigEvents.events += 2 * getCombination( spConf[sp1 - 1].npop, 2 ); */

/*       for(i = 0; i < spConf[sp1 - 1].npop; ++i) */
/* 	for( j = 0; j < spConf[sp1 -1].npop; ++j) */
/* 	  { */
	    
/* 	    if( i == j) */
/* 	      continue; */
	    
/* 	    cmdMigEvents.events++; */
	    
/* 	    cmdMigEvents */
	    
/* 	  } */
      
/*     } */
/*   else if( nvalues == 3 ) */
/*     { */
/*       /\* two species set values for all pairs of pops *\/ */
/*     } */


/* } */



int getValuesAfterFlag(int arg, int argc, char **argv)
{

  int valArgs = 0, i;
  
  for(i = arg; i < argc; ++i)
    {
      
      if( argv[i][0] != '-')
	++valArgs;
      else
	break;
      
    }
  
  return valArgs;
}


int getPopIndex(int species, int pop)
{
  return ( nodes[ species - 1].offset + pop - 1);
  //return ( spConf[ species - 1].offset + pop - 1);
}

int setEffectivePopulationSize( int nvalues, int arg, int argc, char **argv, double* popArray)
{
  int i,sp, pop, popIndex;
  double Ne; 
  
  if( nvalues  < 1 )
    {
      fprintf(stderr, " \n\n -- number is required after -N\n\n");
      return(0);
    }

  if( nvalues == 1 )
    {
      Ne = atof( argv[++arg] );
      
      for( i = 0; i < globalVar.totalPops; ++i)
	popArray[i] = Ne;
      
      return 1;
      
    }

  if( nvalues == 2 )
    {
      sp = atoi( argv[ ++arg ] );

      assert( sp > 0 && sp  < tree.leaves + 1);

      Ne = atof( argv[ ++arg ] );

      assert( Ne > 0);

            
      for( i = nodes[ sp - 1].offset; i < nodes[ sp- 1].offset + nodes[sp-1].npop; ++i)
	popArray[i] = Ne;


      return 1;
    }

  if( nvalues == 3 )
    {

      sp = atoi( argv[ ++arg ] );
      
      assert( sp > 0 && sp < tree.leaves + 1);

      pop = atoi( argv[++arg ] );
      
      assert( pop > 0 && pop < nodes[ sp - 1].npop+1 );

      Ne = atof( argv[ ++arg ]);
      
      popIndex = getPopIndex( sp, pop );

      popArray[ popIndex ] = Ne;

      return 1;
    }

  return 0; 
  
}


int setSamplingTimes( int nvalues, int arg, int argc, char **argv, double *samplingTime)
{

  int i, sp, pop, popIndex;

  double s = 0.;

  if(nvalues < 1)
    {
      fprintf( stderr, "\n\n --- number is required after -samplingTime\n\n");
      
      return 0;
    }

  if( nvalues == 1 )
    {
      s = atof( argv[++arg]);

      for( i = 0; i < globalVar.totalPops; ++i)
	samplingTime[i] = s;

      return 1;
      
    }

  if( nvalues == 2)
    {
      sp = atoi( argv[ ++arg ] );

      assert( sp > 0 && sp < tree.leaves + 1);

      s = atof( argv[ ++arg ] );

      /* printf("SAMPLINGTIME sampling time for %d is %f\n", sp, s); */

      for( i = nodes[ sp - 1].offset; i < nodes[ sp- 1].offset + nodes[sp-1].npop; ++i)
	{
	  assert( i >= 0);
	  samplingTime[i] = s;
	  /* printf("***SAMPLINGTIME sampling time for %d is %f\n", i, s); */
	}

      return 1;
    }

    if( nvalues == 3 )
    {

      sp = atoi( argv[ ++arg ] );
      
      assert( sp > 0 && sp < tree.leaves + 1);

      pop = atoi( argv[++arg ] );
      
      assert( pop > 0 && pop < nodes[ sp - 1].npop+1 );

      s = atof( argv[ ++arg ]);
      
      popIndex = getPopIndex( sp, pop );

      samplingTime[ popIndex ] = s;

      return 1;
    }

  return 0; 

}

int setGrowthRates( int nvalues, int arg, int argc, char **argv, double* popArray)
{
  int i,sp, pop, popIndex;
  double G; 
  
  if( nvalues  < 1 )
    {
      fprintf(stderr, " \n\n -- number is required after -G\n\n");
      return(0);
    }

  if( nvalues == 1 )
    {
      G = atof( argv[++arg] );
      
      for( i = 0; i < globalVar.totalPops; ++i)
	popArray[i] = G;
      
      return 1;
      
    }

  if( nvalues == 2 )
    {
      sp = atoi( argv[ ++arg ] );

      assert( sp > 0 && sp  < tree.leaves + 1);

      G = atof( argv[ ++arg ] );
            
      for( i = nodes[ sp - 1].offset; i < nodes[ sp- 1].offset + nodes[sp-1].npop; ++i)
	popArray[i] = G;

      return 1;
    }

  if( nvalues == 3 )
    {

      sp = atoi( argv[ ++arg ] );
      
      assert( sp > 0 && sp < tree.leaves + 1);

      pop = atoi( argv[++arg ] );
      
      assert( pop > 0 && pop < nodes[ sp - 1].npop+1 );

      G = atof( argv[ ++arg ]);
      
      popIndex = getPopIndex( sp, pop );

      popArray[ popIndex ] = G;

      return 1;
    }

  return 0; 
  
}


int setMigrationMatrix( int nvalues, int arg, int argc, char ** argv, double** migMat)
{
  
  int i,j, sp, sp2, minCoord, 
    maxCoord, minCoord2, maxCoord2, 
    c1, c2, pop1, pop2;

  double mig = 0.;

  if( nvalues == 1)
    {
      
      mig = atof( argv[++arg] );

      for( i = 0; i < globalVar.totalPops-1; ++i)
	{

	  for( j = i+1; j < globalVar.totalPops; ++j)
	    {
	      /* attention... this is different than ms. 
		 in ms the migration rate in the command 
		 line is divided by the numbe rof populations. 
		 TODO.... put that in the manual 
	      */
	      migMat[i][j] = mig; // /( globalVar.totalPops - 1);
	      migMat[j][i] = mig; // /( globalVar.totalPops - 1);
	    }
	}

    }

  if(nvalues == 2)
    {
      /* assumed scenario: 1 species 
	 total migration rate for the species, 
	 which is equally distributed among all populations
      */
      
      sp = atoi( argv[ ++arg ] );

      assert( sp > 0 && sp < tree.leaves + 1);

      mig = atof( argv[ ++arg ] );

      minCoord = nodes[ sp - 1].offset;

      maxCoord = nodes[ sp - 1].offset + nodes[ sp - 1].npop - 1;

      for( i = minCoord; i < maxCoord; ++i)
	{
	  for( j = i + 1; j < maxCoord + 1; ++j)
	    {
	      /* attention... this is different than ms. 
		 in ms the migration rate in the command 
		 line is divided by the numbe rof populations. 
		 TODO.... put that in the manual 
	      */
	      migMat[i][j] = mig; // /(nodes[sp - 1].npop - 1);
	      migMat[j][i] = mig; // /(nodes[sp - 1].npop - 1);
	    }
	}
	  
      return 1;
    }

  if( nvalues == 3 )
    {
      sp = atoi( argv[ ++arg ] );
      
      sp2 = atoi( argv[ ++arg ] );

      mig = atof( argv [++arg ] );
      
      assert( sp > 0 && sp < tree.leaves + 1);

      assert( sp2 > 0 && sp < tree.leaves + 1);

      minCoord = nodes[ sp - 1 ].offset;
      
      minCoord2 = nodes[ sp2 - 1].offset;

      maxCoord2 = nodes[ sp2 - 1 ].offset + 
	/* nodes[ sp2 - 1].offset  + */ nodes[ sp2 - 1].npop -1; 

      maxCoord = nodes[ sp - 1].offset +
	/*nodes[ sp - 1].offset + */nodes[ sp - 1].npop - 1;

      for( i = minCoord; i < maxCoord + 1; ++i)
	for( j = minCoord2 ; j < maxCoord2 + 1; ++j)
	  {
	    migMat[i][j] = mig;
	    //migMat[j][i] = mig;
	  }
      
      return 1;
    }

  if( nvalues == 4 )
    {
      /* 1 species, two specific populations */
      sp = atoi( argv[ ++arg ] );
      
      assert( sp > 0 && sp < tree.leaves + 1);

      pop1 = atoi( argv[ ++arg] ) - 1;

      pop2 = atoi( argv[ ++arg] ) - 1;
      
      assert(pop1 > -1 && pop1 < nodes[sp].npop);
      
      assert(pop2 > -1 && pop2 < nodes[sp].npop);

      c1 = nodes[sp - 1].offset + pop1;
      c2 = nodes[sp - 1].offset + pop2;

      migMat[c1][c2] = atof(argv[++arg]);

      return 1;      
    }

  if( nvalues == 5 )
    {
      sp = atoi( argv[ ++arg ] );

      assert( sp > 0 && sp < tree.leaves + 1);
      
      pop1 = atoi( argv[ ++arg ] ) - 1;

      sp2 = atoi( argv[ ++arg ] );
      
      assert( sp2 > 0 && sp < tree.leaves + 1);

      pop2 = atoi( argv[ ++arg ] ) - 1;

      mig = atof( argv[ ++arg ] );

      c1 = nodes[sp-1].offset + pop1;

      c2 = nodes[sp2-1].offset + pop2;

      assert( pop1 > -1 && pop1 < nodes[sp-1].npop);
      
      if( pop2 < 0 || pop2 > nodes[sp2-1].npop - 1)
	{
	  fprintf(stderr, "sp2: %d, pop2: %d, nodes[sp2].npop: %d\n", 
		  sp2, pop2, nodes[sp2-1].npop);
	  
	  assert( pop2 > -1 && pop2 < nodes[sp2-1].npop);
	  
	}

      assert( c1 > -1 && c1 < globalVar.totalPops);
      
      assert( c2 > -1 && c2 < globalVar.totalPops);
      
      migMat[c1][c2] = mig;

      return 1;
      
    }

  return 0;  
}


int setJointEvent( int nvalues, int arg, int argc, char **argv, struct deventArray *pastEvents)
{

  int sp1, sp2, pop1 = -1, pop2 = -1, success = 0;
  
  pastEvents[globalVar.currentEvents - 1].detype = 'j';

  pastEvents[globalVar.currentEvents - 1].isSpeciation = 0;
  
  pastEvents[globalVar.currentEvents - 1].time = atof( argv[++arg] );

  pastEvents[globalVar.currentEvents -1].type = 0;
  
  sp1 = pastEvents[globalVar.currentEvents - 1].sp1 = atoi( argv[++arg] ) - 1;

  if( nvalues == 4 )
    {
      pop1 = pastEvents[globalVar.currentEvents - 1].popi = atoi(argv[++arg] ) + nodes[sp1].offset - 1;
      
      pop2 = pastEvents[globalVar.currentEvents - 1].popj = atoi(argv[++arg] ) + nodes[sp1].offset - 1;

      success = 1;
    }

  if( nvalues == 5 )
    {

      pop1 = pastEvents[globalVar.currentEvents - 1].popi = atoi(argv[++arg] ) + nodes[sp1].offset - 1;
       
      sp2 = pastEvents[globalVar.currentEvents - 1].sp2 = atoi(argv[++arg] ) - 1;
      
      pop2 = pastEvents[globalVar.currentEvents - 1].popj = atoi(argv[++arg] ) + nodes[sp2].offset - 1;
      
      success = 1;
    }

  assert( pop1 >= 0 && pop2 >= 0);

  initial_pars.cp.ignoreSpeciationJoints[pop1][pop2] = 1;
    
  return(success);

}


int setPopSizeChange( int nvalues, int arg, int argc, char **argv, struct deventArray *pastEvents)
{
  int sp1,  success = 0, i, offset, popi;
  double time;
  
  
  /* time, sizechange */
  if ( nvalues == 2 )
    {
      
      pastEvents[globalVar.currentEvents - 1].time = atof( argv[++arg] );  
    
      pastEvents[globalVar.currentEvents - 1].detype = 'N';
      
      pastEvents[globalVar.currentEvents - 1].paramv = atof( argv[++arg] );

      pastEvents[globalVar.currentEvents - 1].type = 0;

      success = 1;
    }

  /* time, speceis, sizechange */
  else if( nvalues == 3 )
    {
      time = pastEvents[globalVar.currentEvents - 1].time = atof( argv[++arg] );  
            
      sp1 = atoi(argv[++arg]);

      assert(sp1 > 0);
      
      if(globalVar.maxEventsSize  < globalVar.currentEvents + nodes[sp1 -1].npop)
	{
	  resizePastEvents( nodes[sp1 - 1].npop );
	}  

      offset = nodes[sp1 - 1].offset;

      assert( sp1 < tree.leaves && sp1 > 0 );

   
      
      for( i = 0; i < nodes[sp1 - 1].npop; ++i)
	{
	  	  
	  pastEvents[globalVar.currentEvents - 1].time = time;
	  
	  
	  pastEvents[globalVar.currentEvents -1].detype = 'n';

	  pastEvents[globalVar.currentEvents - 1].paramv = atof( argv[arg + 1]);
	  
	  pastEvents[globalVar.currentEvents - 1].popi = i + offset;

	  ++globalVar.currentEvents;
	}
      --globalVar.currentEvents;

      success = 1;
    }

  else if( nvalues == 4 )
    {
      pastEvents[globalVar.currentEvents - 1].time = atof( argv[++arg] );  

      sp1 = atoi(argv[++arg]);
      
      assert(sp1 > -1 && sp1 < tree.leaves );
      
      offset = nodes[sp1 - 1].offset;

      pastEvents[globalVar.currentEvents - 1].detype = 'n';
      
      popi = pastEvents[globalVar.currentEvents - 1].popi = atoi(argv[++arg]) + offset - 1;

      assert( popi > -1 && popi < offset + nodes[sp1 - 1].npop);

      pastEvents[globalVar.currentEvents - 1].paramv = atof(argv[++arg]);
      
      success = 1;
      
    }

  return success;
}



int setPopGrowthRate( int nvalues, int arg, int argc, char **argv, struct deventArray *pastEvents)
{
  int sp1,  success = 0, i, offset, popi;
  double time;
    
  /* time, sizechange */
  if ( nvalues == 2 )
    {
      
      pastEvents[globalVar.currentEvents - 1].time = atof( argv[++arg] );  
    
      pastEvents[globalVar.currentEvents - 1].detype = 'G';
      
      pastEvents[globalVar.currentEvents - 1].paramv = atof( argv[++arg] );

      pastEvents[globalVar.currentEvents - 1].type = 0;

      success = 1;
    }

  /* time, speceis, sizechange */
  else if( nvalues == 3 )
    {
      time = pastEvents[globalVar.currentEvents - 1].time = atof( argv[++arg] );  
            
      sp1 = atoi(argv[++arg]);
      
      if(globalVar.maxEventsSize  < globalVar.currentEvents + nodes[sp1 -1].npop)
	{
	  resizePastEvents( nodes[sp1 - 1].npop );
	}  

      offset = nodes[sp1 - 1].offset;

      assert( sp1 < tree.leaves && sp1 >= 0 );
      
      for( i = 0; i < nodes[sp1 - 1].npop; ++i)
	{
	  	  
	  pastEvents[globalVar.currentEvents - 1].time = time;
	  
	  
	  pastEvents[globalVar.currentEvents -1].detype = 'g';

	  pastEvents[globalVar.currentEvents - 1].paramv = atof( argv[arg + 1]);
	  
	  pastEvents[globalVar.currentEvents - 1].popi = i + offset;

	  ++globalVar.currentEvents;
	}
      --globalVar.currentEvents;

      success = 1;
    }

  else if( nvalues == 4 )
    {
      pastEvents[globalVar.currentEvents - 1].time = atof( argv[++arg] );  

      sp1 = atoi(argv[++arg]);
      
      assert(sp1 > -1 && sp1 < tree.leaves );
      
      offset = nodes[sp1 - 1].offset;

      pastEvents[globalVar.currentEvents - 1].detype = 'g';
      
      popi = pastEvents[globalVar.currentEvents - 1].popi = atoi(argv[++arg]) + offset - 1;

      assert( popi > -1 && popi < offset + nodes[sp1 - 1].npop);

      pastEvents[globalVar.currentEvents - 1].paramv = atof(argv[++arg]);
      
      success = 1;
      
    }

  return success;
}

void checkAndRelabelEvents(int e, int source, int target, double time)
{
  int i;
  
  /* if there is an event at a later time point 
     that involves the source population,
     then the source population should be relabeled
  */

  for( i = 0; i < e; ++i)
    {
      /* joint event */
      if( pastEvents[i].detype == 'j' && 
	  pastEvents[i].time > time )
	{
	  if( pastEvents[i].popi == source )
	    pastEvents[i].popi = target;
	  else if( pastEvents[i].popj == source )
	    pastEvents[i].popj = target;
	}
    }
  
}


void nodesToCoalEvents()
{
  
  
  int i, numberOfEvents, j, targetPop, sourcePop, x;
  
  /* for(i = 0; i < globalVar.totalPops; ++i) */
  /*   { */
  /*     for( j = 0; j < globalVar.totalPops; ++j) */
  /* 	printf(" %d ", pars.cp.ignoreSpeciationJoints[i][j] ); */
  /*     printf("\n"); */
  /*   } */
  
  if( globalVar.currentEvents + nodes->sumofnodespop > globalVar.maxEventsSize)
    {
      globalVar.maxEventsSize += nodes->sumofnodespop;
      
      pastEvents = realloc( pastEvents, globalVar.maxEventsSize * sizeof( struct deventArray ) );
      
      assert( pastEvents != NULL);
    }

  
  /* internal nodes */
  for( i = tree.leaves; i < 2 * tree.leaves - 1; ++i)
    {
      
      targetPop = nodes[i].maxIndPop;
      
      /* printf("targetPop: %d\n", targetPop); */
      
      /* create a joint event for all pops of the node. Each pop goes to the population with the largest index */
      numberOfEvents = nodes[i].npop - 1;
      
      /* I must make one more comparison because I don't know which pop has the maxIndex */
      for( j = 0; j < numberOfEvents + 1; ++j)
	{
	  
	  /* for( x = 0; x < globalVar.totalPops; ++x) */
	  /*   { */
	  /*     for( y = 0; y < globalVar.totalPops; ++y) */
	  /* 	printf("%d,", pars.cp.ignoreSpeciationJoints[x][y]); */
	  /*     printf("\n"); */
	  /*   } */
	  
	  
	  sourcePop = nodes[i].pops[j];

	  /* printf("sourcePop: %d\n", sourcePop); */
	  
	  /* printf("set event: %d --> %d\n", nodes[i].pops[j], nodes[i].maxIndPop); */
	  
	  if( targetPop == sourcePop )
	    continue;
	  
	  assert( sourcePop < globalVar.totalPops );
	  
	  assert( targetPop < globalVar.totalPops );
	  
	  /* printf(" \nmeta %d %d %d\n", nodes[i].pops[j], nodes[i].maxIndPop, globalVar.totalPops); */
	  
	  if( pars.cp.ignoreSpeciationJoints[ sourcePop ] [ targetPop ] == 1 ||
	      pars.cp.ignoreSpeciationJoints[ targetPop ] [ sourcePop ] == 1 )
	    continue;
	  
	  for( x = 0; x < globalVar.totalPops; ++x)
	    {
	      pars.cp.ignoreSpeciationJoints[targetPop][x] = pars.cp.ignoreSpeciationJoints[targetPop][x] | pars.cp.ignoreSpeciationJoints[sourcePop][x];
	      pars.cp.ignoreSpeciationJoints[x][targetPop] = pars.cp.ignoreSpeciationJoints[x][targetPop] | pars.cp.ignoreSpeciationJoints[x][sourcePop];
	    }
	  
	  
	  /* for( x = 0; x < globalVar.totalPops; ++x) */
	  /*   { */
	  /*     for( y = 0; y < globalVar.totalPops; ++y) */
	  /* 	printf("%d,", pars.cp.ignoreSpeciationJoints[x][y]); */
	  /*     printf("\n"); */
	  /*   } */
	  /* printf("\n==========\n\n"); */
	  
	  	  
	  
	  /* if( pars.cp.ignoreSpeciationJoints[ sourcePop ] [ targetPop ] == 1 || */
	  /*     pars.cp.ignoreSpeciationJoints[ targetPop ] [ sourcePop ] == 1 ) */
	  /*   continue; */
	  
	  
	  checkAndRelabelEvents(globalVar.currentEvents, sourcePop, targetPop, nodes[i].coalAge);
	  
	  
	  pastEvents[globalVar.currentEvents].detype='j';

	  pastEvents[globalVar.currentEvents].type = 1;
	  
	  pastEvents[globalVar.currentEvents].isSpeciation = 1;
	 
	  
	  pastEvents[globalVar.currentEvents].time = nodes[i].coalAge;
	  
	  pastEvents[globalVar.currentEvents].popi = nodes[i].pops[j];
	  
	  pastEvents[globalVar.currentEvents].popj = nodes[i].maxIndPop;
	  
	  ++globalVar.currentEvents;
	  
	}
    }
}



int getparsMutationModel(int argc, char *argv[], char name[MAX_NAME_LEN])
{

  int finiteSiteModel = 0;

  int i,j;

  model = NONE;
    
  scaleTrees=0;
  treeScale=0.0;
  scaleBranches=0;
  branchScale=0.0;
  
  maxPartitions=1;
  numPartitions=1;
  
  userSeed = 0;
  
  numCats=1;
  rateHetero=NoRates;
  catRate[0]=1.0;
  gammaShape=1.0;
  
  invariableSites=0;
  proportionInvariable = 0.0;
  
  equalFreqs = 1;
  equalTstv = 1; /* equal transition transversion */
  tstv=0.50002; 
  
  for (i = 0; i < NUM_AA_REL_RATES; i++) {
    aaRelativeRate[i] = 1.0;
  }
  
  for (i = 0; i < NUM_AA; i++) {
    aaFreq[i] = 1.0;
  }
  
  aaFreqSet = 0;

  /* number of sites equals to the number of sites that can be mutated */
  numSites = initial_pars.cp.msites;

  /* number of datasets per tree */
  /*TODO put 1 now, but later this can be from a command line parameter */
  numDatasets = 1;

  /* is an ancestor Seq given */
  ancestorSeq = 0;

  
  writeAncestors=0;
  
  writeRates=0;
  
  verbose=1;
  
  fileFormat = PHYLIPFormat;
  
  quiet=0;

  char *stringModel = calloc( MAXSIZEMODELNAME, sizeof(char) );

  char stringTemp[MAX_NAME_LEN];

  for( i = 1; i < argc; ++i)
    {
      stringToUpper( argv[i], stringTemp );

      if(strcmp("-name", argv[i]) == 0)
	{
	  
	  checkArguments(i, argc, argv);
	  
	  strcpy(name, argv[++i]);

	  continue;
	}
      
      if (strcmp("-OFORMAT" , stringTemp) == 0 ) 
	{
	  
	  assert( i < argc-1 );
	  
	  stringToUpper(argv[++i], stringTemp);
	  
	  if (strcmp("FASTA", stringTemp) == 0)
	    fileFormat=FASTAFormat;
	  
	  else if(strcmp("PHYLIP" , stringTemp) == 0)
	    fileFormat = RelaxedFormat;
	  else 
	    {
	      fprintf(stderr, "File format %s does not exist... Exit!!\n", stringTemp);
	      assert(0);
	    }
	  
	  continue;
	}
    
      if(strcmp("-mm", argv[i]) == 0 )
	{
	  
	  finiteSiteModel = 1;

	  globalVar.finiteSiteModel = 1;

	  stringToUpper(argv[++i], stringModel);
	  
	  for( j = JC; j < numModels; ++j)
	    {
	  
	      if( strcmp( modelNames[j], stringModel ) == 0)
		{
	
		  model = j;

		  if( model <= GTR )
		    {
		      isNucModel = 1;
		      numStates = 4;
		    }
		  else
		    {
		      isNucModel = 0;
		      numStates = 20;
		    }

		  break;
		}
	      
	      else if (strncmp(stringModel, "REV", 3)==0) 
		{
		  model=GTR;
		  isNucModel = 1;
		  numStates = 4;
		  break;
		}
	    }
	  
	  if( j == numModels )
	    {
	      fprintf(stderr, "(j: %d, numModels: %d)...model %s is not yet impelemnted\n", j, numModels, argv[i]);
	      assert( j < numModels );
	    }

	  if( model == -1)
	    {
	      fprintf(stderr, "Unknown model: %s\n\n", argv[i]);
	      exit(-1);
	    }
	  
	  continue;
	}

      if(strcmp("-length", argv[i]) == 0)
	{
	  fprintf(stderr, "length should be provided with the -msites parameter for the finite site model\n");

	  assert(0);
	  
	}

      if(strcmp("-frequencies", argv[i]) == 0 )
	{
	  double freqsum = 0.;
	  int nvalues = getValuesAfterFlag(i+1, argc, argv );
	  if(nvalues != NUM_NUC)
	    {
	      fprintf(stderr, "After -frequencies you should provide 4 arguments that sum up to 1.0\n");
	      assert(nvalues == NUM_NUC);
	    }

	  for(j = 0; j < NUM_NUC; ++j)
	    freqsum += atof(argv[i + j + 1]);

	  
	  for(j = 0; j < NUM_NUC; ++j)
	    nucFreq[j] = atof(argv[i + j + 1])/freqsum;

	  i += NUM_NUC;

	  continue;
	}

      if(strcmp("-rates", argv[i]) == 0)
	{
	  int nvalues = getValuesAfterFlag(i+1, argc, argv);
	  
	  if(nvalues !=  NUM_NUC_REL_RATES )
	    {
	      fprintf(stderr, "Rates should be %d. You provide however %d\n", 
		      (NUM_NUC * (NUM_NUC - 1) )/2, nvalues);
	      assert(0);
	    }

	  for( j = 0; j < nvalues; ++j)
	    {
	      nucRelativeRates[j] = atof(argv[i + j + 1]);
	    }

	  /* rates are relative to the last one */
	  for( j = 0; j < nvalues; ++j)
	    nucRelativeRates[j] /= nucRelativeRates[NUM_NUC_REL_RATES - 1];

	  i += nvalues;
	  
	  continue;

	}

      if(strcmp("-titv", argv[i]) == 0)
	{
	  equalTstv = 0;
	  tstv = atof(argv[++i]);
	  continue;
	}

    

      if(strcmp( "-GAMMACATEGORIES", stringTemp) == 0 )
	 
	{
	  rateHetero = DiscreteGammaRates;
	  numCats = atoi(argv[++i]);
	  assert(numCats > 1);
	  assert(numCats < MAX_RATE_CATS + 1);
	  continue;
	}

	if(strcmp("-alpha", argv[i]) == 0)
	  {
	    assert(rateHetero != CodonRates );
	    if(rateHetero == NoRates)
	      rateHetero = GammaRates;

	    gammaShape = atof(argv[++i]);

	    assert(gammaShape > 0.);
	    continue;
	  }

      

	if(strcmp("-invariable", argv[i]) == 0)
	  {
	    proportionInvariable = atof(argv[++i]);
	    assert( proportionInvariable <= 1.);
	    continue;
	  }
	
    } /* end of arguments loop */
    
  free(stringModel);  
  return finiteSiteModel;

}


void getparsCMD_multiSpeciesSimulate(int argc, char *argv[], int *phowmany, FILE **phyloFile, char phyloFileName[FILESIZE], FILE **CoalescentFile, char CoalescentFileName[FILESIZE], double *partIsoPeriod, double *partMaxMigration, double *samplingRate, double *mu, double *lambda, double *mut, char phyloInputFileName[FILESIZE], FILE **phyloInputFile, int *phyloMode, double *torigin, double *oldestOrigin , unsigned int *timeSeed)
  
{
  
  int arg, i, j, npop, currentArg = 0, tempSpeciesTotalSamples = 0, nvalues = 0, *tempINT, **tempIntMat; 
  
  double thetaBasis, **tempDoubleMat = NULL, *tempDoubleArray = NULL; 
  
  char stringTemp[MAX_NAME_LEN];

  if( count == 0 ) {

    initial_pars.cp.nsam = 0;
    
    if( argc < 4 )
      { 
	fprintf(stderr,"Too few command line arguments\n"); 
	assert(argc >= 4);
      }
    
    tree.leaves = atoi( argv[1] );

    globalVar.totalPops = tree.leaves;

    /* allocate intial memory for the population sizes */
    initial_pars.cp.size = calloc( tree.leaves, sizeof( double ) );

    for( i = 0; i < tree.leaves; ++i)
      initial_pars.cp.size[i] = 1.;

    /* allocate initial memory for the samplingTimes */
    initial_pars.cp.samplingTime = calloc( tree.leaves, sizeof( double ) );
    
    for( i = 0; i < tree.leaves; ++i)
      initial_pars.cp.samplingTime[i] = 0.;
    
    /* allocate initial memory for the migration matrix */
    initial_pars.cp.mig_mat = calloc( tree.leaves, 
				      sizeof(double*) );

    assert( initial_pars.cp.mig_mat != NULL );

    for( i = 0; i < tree.leaves; ++i)
      {
	initial_pars.cp.mig_mat[i] = calloc( tree.leaves, 
					     sizeof(double) );
	
	assert( initial_pars.cp.mig_mat[i] != NULL );
	
      }

    /***************************************************/

    initial_pars.cp.ignoreSpeciationJoints = calloc( globalVar.totalPops, 
						     sizeof(int*) );
    
    assert( initial_pars.cp.ignoreSpeciationJoints != NULL);
    

    for( i = 0; i < globalVar.totalPops; ++i)
      {
	initial_pars.cp.ignoreSpeciationJoints[i] = calloc( globalVar.totalPops, sizeof(int) );

	assert(initial_pars.cp.ignoreSpeciationJoints[i] != NULL);
      }
          
    assert( tree.leaves > 0);
    
    /* speciesSamples = malloc( tree.leaves * sizeof(int) ); */
    
    
    /* allocate memory for the tree  */
    nodes = malloc( (tree.leaves * 2 -1)*sizeof(struct TREEN) );

    tree.nnode = tree.leaves * 2 -1;    

    tree.internalNodes = tree.leaves - 1;

    /* spConf = malloc( tree.leaves * sizeof( struct SPECIES ) ); */

    
    /** set the offset for each species
	this will be used to access the migration matrix correctly */
    for(i = 0; i < tree.leaves; ++i)
      {
		
	nodes[2 * tree.leaves - 2 - i].offset = -1; 
	nodes[i].offset = i;

	
	nodes[2 * tree.leaves - 2 -i].nseqs = 0;
	nodes[i].nseqs = 0;
      }
    
    assert( argc > 1 + tree.leaves);

    for( i = 0; i < tree.nnode; ++i)
      {
	nodes[i].sampleSize = 0;
      }
    
    for(i = 0; i < tree.leaves; ++i)
      {
	
	if( argv[2+i][0] == '-' )
	  {
	    
	    fprintf(stderr, "Argument %d should not start with -. The sample size of species %d is expected\n\n", 2+i, i+1);
	    
	    assert(argv[2+i][0] != '-');
	    
	  }
	
	nodes[i].sampleSize = atoi(argv[2+i]);
	
	initial_pars.cp.nsam += nodes[i].sampleSize;

	nodes[i].Ne = 1.;
	
	//nodes[2*tree.leaves - 2 - i ].npop = 1;
	
	nodes[2*tree.leaves - 2 -i].npop = 0;

	nodes[i].npop = 1;

	nodes[i].pops = NULL;

	nodes[2*tree.leaves - 2 - i].pops = NULL;

	if( i != tree.leaves - 1)
	  {
	    nodes[2*tree.leaves - 2 - i].samplePops = calloc( 1, sizeof(int) );
	
	  }
	
	nodes[i].samplePops = calloc( 1, sizeof(int) );

	nodes[2*tree.leaves - 2 - i].samplePops[0] = nodes[2*tree.leaves - 2 - i].sampleSize;

	nodes[i].samplePops[0] = nodes[i].sampleSize;

	nodes[2*tree.leaves -2-i].maxIndPop = 0;

	nodes[i].maxIndPop = i;
	
      }
        
    currentArg = 2 + tree.leaves;
        
    if( initial_pars.cp.nsam <= 0 ) 
      { 
      	fprintf(stderr,"Total sample size should be positive. Now it is: %d\n", initial_pars.cp.nsam); 
	assert( initial_pars.cp.nsam > 0 );
      }
    
    *phowmany = atoi( argv[ currentArg++ ] );
    
    if( *phowmany  <= 0 ) 
      { 
	fprintf(stderr,"Second argument error. howmany <= 0. \n"); 
	assert( *phowmany > 0 );
      }

    /* printf("howmany: %d\n", *phowmany); */

    initial_pars.commandlineseedflag = 0 ;
    
    initial_pars.output_precision = 7;
    
    initial_pars.cp.r = initial_pars.mp.theta =  initial_pars.cp.f = 0.0 ;
    
    initial_pars.cp.track_len = 0. ;
    
    initial_pars.cp.npop = npop = 1 ;
    
    initial_pars.mp.segsitesin = 0 ;
    
    pars.mp.treeflag = initial_pars.mp.treeflag = 0 ;
    
    initial_pars.mp.timeflag = 0 ;
    
    initial_pars.mp.mfreq = 1 ;
    
    
    /* initial growth rates */
    initial_pars.cp.alphag = calloc( tree.leaves , sizeof( double ) );
    
    initial_pars.cp.nsites = 2 ;
    
  } // if count == 0
   
  pars.cp.deventlist = NULL ;
    
  arg = currentArg ;

    
  /* printf("*Total number of populations is: %d\n", initial_pars.cp.npop); */
	
  
  for(arg = currentArg; arg < argc; ++arg)
    {

      stringToUpper( argv[arg], stringTemp );

      if(strcmp( stringTemp, "-SAMPLINGTIME") == 0)
	{
	  globalVar.readSampling = 1;

	  nvalues = getValuesAfterFlag( arg+1, argc, argv);

	  if( setSamplingTimes(nvalues, arg, argc, argv, initial_pars.cp.samplingTime) != 1)
	    {
	      fprintf(stderr, " --- setting the sampling times for populations ... ERROR\n\n\n");

	      assert(0);
	      
	      
	    }

	  
	  continue;
	}

      if(strcmp( stringTemp, "-SAMPLINGRATE") == 0)
	{
	  checkArguments(arg, argc, argv);

	  *samplingRate = atof(argv[++arg]);
	  
	  continue;
	}

      if(strcmp( stringTemp, "-SEED") == 0)
	{
	  checkArguments (arg, argc, argv);

	  *timeSeed = atoll(argv[++arg]);
	  
	  continue;
	}

      if(strcmp( stringTemp, "-IPHYLO") == 0 )
	{
	  checkArguments( arg, argc, argv);

	  strcpy(phyloInputFileName, argv[++arg]);

	  *phyloInputFile = fopen( phyloInputFileName, "r");

	  if( *phyloInputFile == NULL)
	    {
	      fprintf(stderr, "Cannot open properly the file %s: null pointer\n", phyloInputFileName);
	      assert(*phyloInputFile != NULL);
	    }

	  

	  continue;
	}

      if(strcmp( stringTemp, "-OLDESTORIGIN") == 0 )
	{
	  checkArguments ( arg, argc, argv );

	  *oldestOrigin = atof(argv[++arg] );

	  continue;
	}

      if(strcmp(stringTemp, "-PHYLOMODE") == 0 )
	{
	  checkArguments( arg, argc, argv );
	  
	  *phyloMode = atoi(argv[++arg]);

	  if( !(*phyloMode == 1|| *phyloMode ==2 || *phyloMode == 3 || *phyloMode ==4 || *phyloMode == 0) )
	    {
	      fprintf(stderr, "PhyloMode can be either 0,1,2,3,4\n");

	      assert(*phyloMode == 1|| *phyloMode ==2 || *phyloMode == 3 || *phyloMode ==4 || *phyloMode == 0);
	    }

	  continue;	  
	}

      if(strcmp("-TORIGIN", stringTemp ) == 0)
	{
	  checkArguments( arg, argc, argv );

	  *torigin = atof( argv[++arg] );

	  if( *phyloMode != 1 && *phyloMode != 2)
	    {
	      fprintf(stderr, "\n\nPlease specify -phyloMode first\n");
	      fprintf(stderr, "\nIt must be either 1 or 2\n\n");
	      exit(-1);
	    }
	  
	  continue;
	}

            
      if(strcmp(stringTemp, "-PARTISOLATION") == 0)
	{
	  globalVar.partialIsolationModel = 1;
	  continue;
	}
      
      if(strcmp( stringTemp,  "-PARTMAXMIGRATION" ) == 0 )
	{
	  checkArguments(arg, argc, argv);
	  *partMaxMigration = atof(argv[++arg]);
	  globalVar.partialIsolationModel = 1;
	  continue;
	}

      if(strcmp( stringTemp, "-PARTISOPERIOD") == 0 )
	{
	  checkArguments( arg, argc, argv);
	  *partIsoPeriod = atof(argv[++arg]);
	  globalVar.partialIsolationModel = 1;
	  continue;
	}

      if( strcmp( stringTemp, "-OPHYLO") == 0 || strcmp( stringTemp, "-OPHYLOGENY") == 0 )
	{
	  *phyloFile = fopen( phyloFileName, "w" );
	  continue;
	}

      if(strcmp( stringTemp, "-OCOALESCENT") == 0)
	{
	  *CoalescentFile = fopen(CoalescentFileName, "w");
	  pars.mp.treeflag = 1;
	  continue;
	}

      if(strcmp( stringTemp, "-NOSEPARATOR") == 0)
	{
	  if( globalVar.printInSeparateFiles ==  0 )
	    globalVar.printSeparator = 0;
	  continue;
	}

      if(strcmp( stringTemp, "-OSEPARATEFILES") == 0 )
	{
	  globalVar.printInSeparateFiles = 1;
	  globalVar.printSeparator = 0;
	  continue;
	}

      if(strcmp(argv[arg], "-T") == 0)
	{
	  pars.mp.treeflag = 1;
	  continue;
	}
      
      if(strcmp(argv[arg], "-t" ) == 0 )
	{
	  
	  if( getTheta( &arg, argc, argv) == 0)
	    {

	      fprintf(stderr, "\nYou have given wrong number of values for the parameter THETA.\nNumber of species: %d, argc: %d, current_arg: %d\n\n", tree.leaves, argc, arg);
	      	      
	      assert(0);
	    }

	  initial_pars.mp.theta = thetaBasis = nodes[0].theta;
	  	  
	  continue;
	  
	}
           
      if( strcmp(argv[arg], "-r" )  == 0 )
	{
	  
	  if( getRho( &arg, argc, argv) == 0 )
	    {
	      
	      fprintf(stderr, "\nYou have given wrong number of values for the parameter RHO.\nNumber of species: %d, argc: %d, current_arg: %d\n\n", tree.leaves, argc, arg);
	      	      
	      assert(0);

	    }
	  continue;
	}

      if( strcmp(argv[arg], "-rsites" ) == 0)
	{
	  checkArguments(arg, argc, argv);

	  initial_pars.cp.nsites = atoi( argv[ ++arg ] );

	  continue;
	}

      if( strcmp( argv[arg], "-msites" ) == 0)
	{
	  checkArguments(arg, argc, argv);

	  initial_pars.cp.msites = atoi( argv[ ++arg ] );
	  
	  continue;
	}

      
      if( strcmp( argv[arg], "-G" ) == 0 )
	{

	  globalVar.readG = 1;
	  
	  nvalues = getValuesAfterFlag( arg+1, argc, argv );
	  
	  if( setGrowthRates(nvalues, arg, argc, argv, initial_pars.cp.alphag) != 1)
	    {
	      fprintf( stderr, " --- setting the growth rate... error\n\n");
	      assert ( 0 );
	    }

	  /* for( i = 0; i < globalVar.totalPops; ++i) */
	  /*   printf("pop%d size: %f\n", i, initial_pars.cp.size[i]); */
	  /* printf("\n\n"); */
	  
	  continue;	  
	}

      if( strcmp( argv[arg], "-N" ) == 0 )
	{

	  globalVar.readN = 1;
	  
	  nvalues = getValuesAfterFlag( arg+1, argc, argv );
	  
	  if( setEffectivePopulationSize(nvalues, arg, argc, argv, initial_pars.cp.size) != 1)
	    {
	      fprintf( stderr, " --- setting the effective population size... error\n\n");
	      assert ( 0 );
	    }

	  for( i = 0; i < globalVar.totalPops; ++i)
	    printf("pop%d size: %f\n", i, initial_pars.cp.size[i]);
	  printf("\n\n");
	  
	  continue;
	  
	  
	}

      if( strcmp( argv[ arg ], "-eN" ) == 0 )
	{
	  /* make sure that -I has been read before */
	  globalVar.readeN = 1;
	  
	  /* the number of arguments after the -eN flag */
	  nvalues = getValuesAfterFlag( arg+1, argc, argv);

	  ++globalVar.currentEvents;

	  
	  if( globalVar.currentEvents > globalVar.maxEventsSize )
	    {
	      resizePastEvents(ADDSIZE);
	    }	  

	  setPopSizeChange( nvalues, arg, argc, argv, pastEvents);	  

	  arg += nvalues;
	  
	  continue;
	}

      if( strcmp( argv[ arg ], "-eG" ) == 0 )
	{
	  /* make sure that -I has been read before */
	  globalVar.readeG = 1;
	  
	  /* the number of arguments after the -eN flag */
	  nvalues = getValuesAfterFlag( arg+1, argc, argv);

	  ++globalVar.currentEvents;
	  
	  if( globalVar.currentEvents > globalVar.maxEventsSize )
	    {
	      resizePastEvents(ADDSIZE);
	    }	  

	  setPopGrowthRate( nvalues, arg, argc, argv, pastEvents);	  

	  arg += nvalues;
	  
	  continue;
	}


      
      if( strcmp( argv[ arg], "-I" ) == 0 )
	{
	  /* format expected 
	     -I SPEC NPOP SS1 SS2 ... 
	     Different than ms. I do not expect to find 
	     any migration rate here
	  */
	  
	  assert( globalVar.readMigration == 0);
	  assert( globalVar.readN == 0 );
	  assert( globalVar.readeN == 0 );
	  assert( globalVar.readG == 0);
	  assert( globalVar.readeG == 0);
	  assert( globalVar.readSampling == 0);
	  
	  /* assert( migrationHasBeenSet == 0 ); */

	  int i, species;
	  
	  checkArguments( arg, argc, argv );
	  
	  species = atoi( argv[ ++arg] );

	  assert ( species <= tree.leaves );

	  checkArguments( arg, argc, argv );
	  
	  nodes[species - 1].npop =  atoi( argv[++arg] );	  

	  nodes[species - 1].maxIndPop += nodes[species -1].npop - 1;

	  for( i = species; i < tree.leaves; ++i)
	    {
	      nodes[i].offset += nodes[species - 1].npop - 1;
	      
	      nodes[i].maxIndPop += (nodes[species - 1].npop - 1);
	    }
	  
	  globalVar.totalPops += nodes[species - 1].npop - 1; 
	  
	  tempDoubleMat = realloc( initial_pars.cp.mig_mat, 
				   globalVar.totalPops * 
				   sizeof( double * ) );

	  if( tempDoubleMat == NULL )
	    {
	      fprintf(stderr, "Error in migration matrix realloc\n" );
	      assert( tempDoubleMat != NULL );
	      
	    }

	  initial_pars.cp.mig_mat = tempDoubleMat;

	  tempIntMat = realloc( initial_pars.cp.ignoreSpeciationJoints, 
				globalVar.totalPops *
				sizeof( int* ) );

	  if( tempIntMat == NULL)
	    {
	      fprintf(stderr, "Error in ignoreSpeciation matrix realloc\n");

	      assert( tempIntMat != NULL);
	    }

	  initial_pars.cp.ignoreSpeciationJoints = tempIntMat;

	  for( i = 0; i < globalVar.totalPops - nodes[ species - 1].npop; ++i)
	    {
	      tempDoubleArray = realloc( initial_pars.cp.mig_mat[i], globalVar.totalPops * sizeof( double ) );
	      
	      assert( tempDoubleArray != NULL);

	      initial_pars.cp.mig_mat[i] = tempDoubleArray;

	      tempINT = realloc( initial_pars.cp.ignoreSpeciationJoints[i], globalVar.totalPops * sizeof( int) );

	      assert( tempINT != NULL);

	      initial_pars.cp.ignoreSpeciationJoints[i] = tempINT;

	      for( j = globalVar.totalPops - nodes[species - 1].npop; j < globalVar.totalPops; ++j)
		{
		  initial_pars.cp.mig_mat[i][j] = 0.;

		  initial_pars.cp.ignoreSpeciationJoints[i][j] = 0;
		}

	    }

	  for( i = globalVar.totalPops - nodes[ species - 1].npop; i < globalVar.totalPops; ++i)
	    {
	      initial_pars.cp.mig_mat[i] = calloc( globalVar.totalPops, sizeof( double ) );
	      
	      assert( initial_pars.cp.mig_mat[i] != NULL);


	      initial_pars.cp.ignoreSpeciationJoints[i] = calloc( globalVar.totalPops, sizeof(int) );

	      assert( initial_pars.cp.ignoreSpeciationJoints[i] != NULL);
	      
	    }
	  
	  	  
	  if( nodes[ species - 1 ].npop > 1 )
	    {
	  
	      /* get sure that tempINT points at different 
		 memory location every time
		 -- it seems to work fine
	      */

	      tempINT = NULL;
	      
	      tempINT = realloc( nodes[ species -1].samplePops, nodes[ species - 1].npop * sizeof(int) );
	      
	      printf("memory address of tempINT is: %p and stores the value %d\n", tempINT, tempINT[0]);

	      if( tempINT == NULL )
		{
		  assert(tempINT != NULL);
		}
	      nodes[ species-1].samplePops = tempINT;
	      
	      /* increase the size of the initial_pars.cp.nsize 
	       */

	      tempDoubleArray = NULL;

	      tempDoubleArray = realloc( initial_pars.cp.size, globalVar.totalPops * sizeof( double ) );

	      if( tempDoubleArray == NULL )
		{
		  fprintf(stderr,  "\n\n -- memory realloc error pars.cp.size\n");
		  assert( tempDoubleArray != NULL );
		  
		}
	      

	      initial_pars.cp.size = tempDoubleArray;

	      for( i = 0; i < globalVar.totalPops; ++i)
		{
		  initial_pars.cp.size[i] = 1.;
		}
	      
	      tempDoubleArray = NULL;
	      
	      tempDoubleArray = realloc( initial_pars.cp.alphag, globalVar.totalPops * sizeof( double ) );
	      
	      if( tempDoubleArray == NULL )
		{
		  fprintf(stderr,  "\n\n -- memory realloc error pars.cp.size\n");
		  assert( tempDoubleArray != NULL );
		  
		}
	      
	      initial_pars.cp.alphag = tempDoubleArray;

	      for( i = 0; i < globalVar.totalPops; ++i)
		{
		  initial_pars.cp.alphag[i] = 0.;
		}


	      tempDoubleArray = NULL;
	      
	      tempDoubleArray = realloc( initial_pars.cp.samplingTime, globalVar.totalPops * sizeof( double ) );
	      
	      if( tempDoubleArray == NULL )
		{
		  fprintf(stderr,  "\n\n -- memory realloc error pars.cp.size\n");
		  assert( tempDoubleArray != NULL );
		  
		}
	      
	      initial_pars.cp.samplingTime = tempDoubleArray;

	      for( i = 0; i < globalVar.totalPops; ++i)
		{
		  initial_pars.cp.samplingTime[i] = 0.;
		}

	    }
	  

	  tempSpeciesTotalSamples = 0;
	  
	  for( i = 0; i < nodes[ species - 1].npop; ++i )
	    {
	      
	      printf("currentLocation: %d\n", arg + 1);
	      
	      if( arg + 1 >= argc)
		{
		  fprintf(stderr, "\n\nPerhaps not enough sample sizes for the sub-pops of species %d are given\n\n\n", species);
		  
		  assert( arg + 1  < argc );
		}
	      
	      if( argv[arg + 1][0] == '-' )
		{
		  
		  fprintf(stderr, "argument %s (location: %d) is expected to be an integer\n", argv[arg + 1], arg + 1);
		  
		  for( j = arg + 1; j < argc; ++j)
		    fprintf(stderr, "%s ", argv[j]);
		  
		  fprintf(stderr, "\n");
		      
		  assert( argv[ arg + 1 ][0] != '-' );
		  
		}
	      

	      nodes[species - 1]. samplePops[i] = atoi(argv[++arg]);	      
	      
	      tempSpeciesTotalSamples += nodes[species - 1].samplePops[i]; //spConf[ species - 1].samplePops[i];
	      
	    }

	  if( tempSpeciesTotalSamples != nodes[species - 1].sampleSize )
	    {
	      fprintf(stderr, "Total sample size of species %d does not match with the sum of sample sizes for its subpopulations (total: %d, sum: %d)\n", species, nodes[species -1].sampleSize, tempSpeciesTotalSamples);
	      
	      assert( tempSpeciesTotalSamples == nodes[species - 1].sampleSize );
	      
	    }

	  continue;
	  
	}// up to here read -I

    
      if( strcmp("-MIGRATION", stringTemp) == 0 )
	{

	  globalVar.readMigration = 1;

	  nvalues = getValuesAfterFlag(arg+1, argc, argv);

	  if( nvalues < 1 || nvalues > 5)
	    {
	      printf("  --- nvalues: %d ", nvalues);
	      
	      for( i = arg; i < arg+nvalues; ++i)
		printf(" %s ", argv[i]);
	      
	      printf("\n");
	      
	      assert( nvalues > 0 );
	      
	      assert( nvalues < 6 );

	    }
	  
	  /* printf("nvalues: %d\n", nvalues); */

	  setMigrationMatrix( nvalues, arg, argc, argv, initial_pars.cp.mig_mat );
	  
	  continue;
	  
	}/* up to here -migration */


      
      if( strcmp("-ANCESTRALMIGRATION", stringTemp) == 0 )
	{

	  int pop2, pop;

	  globalVar.readMigration = 1;

	  ++globalVar.currentEvents;
	  
	  if( globalVar.currentEvents > globalVar.maxEventsSize )
	    {
	      resizePastEvents(ADDSIZE);
	    }	  
	  	  
	  pastEvents[globalVar.currentEvents - 1].time = atof( argv[++arg] );  
	  
	  pastEvents[globalVar.currentEvents - 1].detype = 'a';
	  
	  nvalues = getValuesAfterFlag(arg+1, argc, argv);
	  
	  if( nvalues < 1 || nvalues > 5)
	    {
	      printf("  --- nvalues: %d ", nvalues);
	      
	      for( i = arg; i < arg+nvalues; ++i)
		printf(" %s ", argv[i]);
	      
	      printf("\n");
	      
	      assert( nvalues > 0 );
	      
	      assert( nvalues < 6 );

	    }
	  
	  /* allocate memory for matrix */
	  pastEvents[globalVar.currentEvents - 1].mat = calloc( globalVar.totalPops, sizeof(double*) );

	  for(pop = 0; pop < globalVar.totalPops; ++pop)
	    pastEvents[globalVar.currentEvents - 1].mat[pop] = calloc( globalVar.totalPops, sizeof(double) );
	  
	  setMigrationMatrix( nvalues, arg, argc, argv, pastEvents[globalVar.currentEvents - 1].mat);
	  	  	  
	  for( pop = 0; pop < globalVar.totalPops; pop++) 
	    {
	    
	      pastEvents[globalVar.currentEvents - 1].mat[pop][pop] = 0.0 ;
	      
	      for( pop2 = 0; pop2 < globalVar.totalPops; pop2++)
		{
		  if( pop2 != pop ) 
		    pastEvents[globalVar.currentEvents - 1].mat[pop][pop] += pastEvents[globalVar.currentEvents - 1].mat[pop][pop2] ;
		}
	      
	    }
	  
	  continue;
	  
	}/* up to here -ancestralmigration */


      if( strcmp( argv[arg], "-ej" ) == 0 )
	{
	  globalVar.readej = 1;

	  nvalues = getValuesAfterFlag(arg+1, argc, argv);
	  
	  if( nvalues < 4 || nvalues > 5 )
	    {
	      fprintf(stderr, "-ej flags accepts 4 or 5 values\n");
	      
	      assert( nvalues > 3 && nvalues < 6 );
	    }

	  ++globalVar.currentEvents;

	  
	  if( globalVar.currentEvents > globalVar.maxEventsSize )
	    {
	      resizePastEvents(ADDSIZE);
	    }	  	  
	  	  

	  setJointEvent(nvalues, arg, argc, argv, pastEvents);

	  arg += nvalues;

	  continue;

	}/* up to here, -ej */

      if(strcmp("-PHYLOMUT", stringTemp) == 0 )
	{
	  checkArguments( arg, argc, argv );
	  *mut = atof(argv[++arg]);
	  continue;
	}

      if( strcmp("-BIRTH", stringTemp) == 0 )
	{
	  checkArguments( arg, argc, argv );
	  *lambda = atof(argv[++arg]);
	  continue;
	}

      if(strcmp("-DEATH", stringTemp) == 0 )
	{
	  checkArguments( arg, argc, argv );
	  *mu = atof(argv[++arg]);
	  continue;
	}
            
    } // read all command line arguments

  

  /* for( i = 0; i < globalVar.currentEvents; ++i) */
  /*   { */
  /*     printf("event: %d -- %c %f, pops: popi %d and popj %d\n", i, pastEvents[i].detype, pastEvents[i].time, pastEvents[i].popi, pastEvents[i].popj); */
      
  /*   } */

  
  /* set diagonal values to total input mig rate for each subpopulation for the migration matrix 
   */
  for( i = 0; i < globalVar.totalPops; ++i)
    {
      
      initial_pars.cp.mig_mat[i][i] = 0.0;
      
      for( j = 0; j < globalVar.totalPops; ++j)
	{
	  if( i == j )
	    continue;
	  
	  initial_pars.cp.mig_mat[i][i] += 
	    initial_pars.cp.mig_mat[i][j];
	}
    }
 
  /* get the total number of populations */
  initial_pars.cp.npop = globalVar.totalPops;
  
}


 void addtoelist( struct devent *pt, struct devent *elist ) 
{
  struct devent *plast, *pevent, *ptemp  ;

  plast = pevent = elist ;
  while(  (pevent != NULL ) && ( pevent->time <= pt->time ) )  {
    plast = pevent ;
    pevent = pevent->nextde ;
  }
  ptemp = plast->nextde ;
  plast->nextde = pt ;
  pt->nextde = ptemp ;
}

void 
free_eventlist( struct devent *pt, int npop )
{
  struct devent *next ;
  int pop ;
   
  while( pt != NULL){
    next = pt->nextde ;
    if( pt->detype == 'a' ) {
      for( pop = 0; pop < npop; pop++) free( (pt->mat)[pop] );
      free( pt->mat );
    }
    free(pt);
    pt = next ;
  }
}

	
/************ make_gametes.c  *******************************************
 *
 *
 *****************************************************************************/

#define STATE1 '1'
#define STATE2 '0'

int   make_gametes(int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list )
{
  int  tip, j,  node ;
  int pickb(int nsam, struct node *ptree, double tt), 
    pickbmf(int nsam, int mfreq, struct node *ptree, double tt) ;

  
  for(  j=ns; j< ns+newsites ;  j++ ) {
    if( mfreq == 1 )
      node = pickb(  nsam, ptree, tt);
    else
      node = pickbmf(  nsam, mfreq, ptree, tt);

    
    
  fprintf(mutationTimesFile, "%.8f ", rndu() * (ptree[ ptree[node].abv ].time - ptree[node].time) + ptree[node].time);
    
    for( tip=0; tip < nsam ; tip++) {
      if( tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
      else list[tip][j] = STATE2 ;
    }
  }

  return 0;
}


/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/

double
ttime( ptree, nsam)
     struct node *ptree;
     int nsam;
{
  double t;
  int i;

  t = (ptree + 2*nsam-2) -> time ;
  for( i=nsam; i< 2*nsam-1 ; i++)
    t += (ptree + i)-> time ;
  return(t);
}


double
ttimemf( ptree, nsam, mfreq)
     struct node *ptree;
     int nsam, mfreq;
{
  double t;
  int i;

  t = 0. ;
  for( i=0;  i< 2*nsam-2  ; i++)
    if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
      t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
  return(t);
}


void parens(FILE *fout,  struct node *ptree, int *descl, int *descr,  int noden, double scale)
{
  double time ;

  if( descl[noden] == -1 ) {
    fprintf(fout, "seq%d:%.10lf", noden+1, ( (ptree+ ((ptree+noden)->abv))->time - (ptree + noden)->time) * scale );
  }
  else{
    fprintf(fout, "(");
    parens(fout,  ptree, descl,descr, descl[noden], scale ) ;
    fprintf(fout, ",");
    parens(fout, ptree, descl, descr, descr[noden], scale ) ;
    if( (ptree+noden)->abv == 0 ) fprintf(fout, ");\n"); 
    else {
      time = (ptree + (ptree+noden)->abv )->time - (ptree+noden)->time ;
      fprintf(fout, "):%.10lf", time * scale );
    }
  }
}

/***  pickb : returns a random branch from the tree. The probability of picking
      a particular branch is proportional to its duration. tt is total
      time in tree.   ****/

int
pickb(nsam, ptree, tt)
     int nsam;
     struct node *ptree;
     double tt;
{
  double x, y, rndu();
  int i;

  x = rndu()*tt;
  for( i=0, y=0; i < 2*nsam-2 ; i++) {
    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
    if( y >= x ) return( i ) ;
  }
  return( 2*nsam - 3  );  /* changed 4 Feb 2010 */
}

int
pickbmf(nsam, mfreq, ptree, tt )
     int nsam, mfreq;
     struct node *ptree;
     double tt;
{
  double x, y;
  int i, lastbranch = 0 ;

  x = rndu()*tt;
  for( i=0, y=0; i < 2*nsam-2 ; i++) {
    if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ){
      y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
      lastbranch = i ;    /* changed 4 Feb 2010 */
    }
    if( y >= x ) return( i ) ;
  }
  return( lastbranch );   /*  changed 4 Feb 2010 */
}


/* pick2()  */
int pick2(int n, int *i, int *j)
{
    
  *i = n * rndu() ;
  while( ( *j = n * rndu() ) == *i )
    ;
  return(0) ;
}



int ranvec(int n, double pbuf[])
{
  int i;
  for(i=0; i<n; i++)
    pbuf[i] = rndu();

  return 0;
}


/* a slight modification of crecipes version */

double gasdev(m,v)
     double m, v;
{
  static int iset=0;
  static float gset;
  float fac,r,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*rndu()-1.0;
      v2=2.0*rndu()-1.0;
      r=v1*v1+v2*v2;
    } while (r >= 1.0);
    fac=sqrt(-2.0*log(r)/r);
    gset= v1*fac;
    iset=1;
    return( m + sqrt(v)*v2*fac);
  } else {
    iset=0;
    return( m + sqrt(v)*gset ) ;
  }
}

