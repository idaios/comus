/**********  segtre_mig.c **********************************
 *
 *	This subroutine uses a Monte Carlo algorithm described in
 *	Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201, to produce
 *	a history of a random sample of gametes under a neutral
 *	Wright-Fisher model with recombination and geographic structure. 
 *	Input parameters
 *	are the sample size (nsam), the number of sites between
 *	which recombination can occur (nsites), and the recombination
 *	rate between the ends of the gametes (r). The function returns
 *	nsegs, the number of segments the gametes were broken into
 *	in tracing back the history of the gametes.  The histories of
 *	these segments are passed back to the calling function in the
 *	array of structures seglst[]. An element of this array,  seglst[i],
 * 	consists of three parts: (1) beg, the starting point of
 *	of segment i, (2) ptree, which points to the first node of the
 *	tree representing the history of the segment, (3) next, which
 *	is the index number of the next segment.
 *	     A tree is a contiguous set of 2*nsam nodes. The first nsam
 *	nodes are the tips of the tree, the sampled gametes.  The other
 *	nodes are the nodes ancestral to the sampled gametes. Each node
 *	consists of an "abv" which is the number of the node ancestral to
 *	that node, an "ndes", which is not used or assigned to in this routine,
 *	and a "time", which is the time (in units of 4N generations) of the
 *	node. For the tips, time equals zero.
 *	Returns a pointer to an array of segments, seglst.

 **************************************************************************/


#include "comus.h"
#include "twister.h"
#include <float.h>

#define NL putchar('\n')
#define size_t unsigned

#define MIN(x, y) ( (x)<(y) ? (x) : (y) )

#define ERROR(message) fprintf(stderr,message),NL,exit(1)

#define SEGINC 80 

#define TOL 0.00000000001

extern int flag;

int nchrom, begs, nsegs;
long nlinks;

static int *nnodes = NULL ;  

double t, cleft , pc, lnpc ;

static unsigned seglimit = SEGINC ;

static unsigned maxchr ;

struct seg{
  int beg;
  int end;
  int desc;
};

struct chromo{
  int nseg;
  int pop;
  struct seg  *pseg;
  double samplingTime;
} ;

static struct chromo *chrom = NULL ;

struct node{
  int abv;
  int ndes;
  float time;
} *ptree1, *ptree2;


static struct segl *seglst = NULL ;


void free_streec_variables()
{
  free( nnodes );
  free( chrom );
  free( seglst );
}




int  randPoissonLinear(double tis, double tsp, double currentTime, double maxMig, double *value)
{

  int cnt = 0, maxCNT=1000, gotNumber = 1;
  
  double t , r1, r2;

  if(currentTime < tis)
    currentTime = tis;
  
  if( tsp < currentTime ) 
    return 0;
  
  do 
    {      
     
      if( cnt > maxCNT )
	{
	  gotNumber = 0;
	  break;
	}

      r1 = rndu();

      t = currentTime - 1./maxMig * log( 1. - r1 );

      r2 = rndu ();
            
      if( t > tsp )
	{
	  gotNumber = 0;
	  break;
	}

      cnt++;

    } while(r2 > t / (tsp - tis )  - tis / (tsp - tis ) );
    
    
    *value = t-currentTime;  
    
  assert( gotNumber == 0 || t - currentTime >= 0.);
  
  return gotNumber;
}

int checkPartialMigration(struct TREEN *phyloNodes, int *activePopsOnNodes, int nnodes, double currentTime, double duration, int *config)
{


  int i, isActive = 0, j = -1, activePops, pop;
  
  for( i = 0; i < nnodes; ++i)
    {

      activePops = 0;
      
      if(phyloNodes[i].isLeaf == 1)
	continue;

      ++j;

      for( pop = 0; pop < phyloNodes[i].npop; ++pop)
	if( config[phyloNodes[i].pops[pop]] > 0 )
	  activePops++;

      if( phyloNodes[i].nseqs > 0 && 
	  phyloNodes[i].isolationTime < currentTime + TOL && 
	  currentTime < phyloNodes[i].coalAge )
	{
	  
	  activePopsOnNodes[j] = activePops - 1;
	  
	  if( activePopsOnNodes[j] < 0 )
	    {
	      assert(activePopsOnNodes[j] == -1 );

	      activePopsOnNodes[j] = 0;
	    }
	  
	  if(activePops - 1 > 0)
	    isActive = 1;
	}
      else
	activePopsOnNodes[j] = 0;
      
    }

  return isActive;

}


int seqsPerNode( struct TREEN node, int *config)
{
  int popIndex, i, nseqs = 0;

  for( i = 0 ; i < node.npop; ++i)
    {
      popIndex = node.pops[i];

      nseqs += config[ popIndex ];

    }

  return nseqs;  
}


int sequenceBelongsToNode(int seq, struct TREEN node, struct chromo *chrom)
{
  int i;
  for( i = 0; i < node.npop; ++i)
    {
      if(node.pops[i] == chrom[seq].pop)
	return 1;
    }

  return 0;
}


int chooseRandomSequenceAtNode( struct TREEN node, int *pop, int *config, struct chromo *chrom, int nchrom)
{
  int nseqs = 0, i, popIndex, randomSeq = 0, j;

  for(i = 0; i < node.npop; ++i)
    {
      popIndex = node.pops[i];

      nseqs += config[ popIndex ];

    }

  assert( nseqs > 0);

  randomSeq = rndu() * nseqs;

  assert( randomSeq < nseqs );

  nseqs = 0;
  for( i = 0; i < node.npop; ++i)
    {
      popIndex = node.pops[i];

      nseqs += config[ popIndex ];

      if(randomSeq < nseqs )
	{
	  *pop = popIndex;
	  break;
	}	  
    }
  
  /* find the right index of the randomSeq */
  j = 0;
  for(i = 0; i < nchrom; ++i)
    {

      if( sequenceBelongsToNode(i, node, chrom )== 1)
	j++;
      
      if( j == randomSeq+1)
	break;
      
    }

  assert( i < nchrom);
  
  return i;
  
}

void printConfig(int *config, struct TREEB tree, struct TREEN *nodes)
{
  int i, j; //, k =0;

  for( i = 0; i < tree.nnode; ++i)
    {
      for( j = 0; j < nodes[i].npop; ++j)
	{
	  printf("population %d of node %d has %d sequences\n", nodes[i].pops[j], i, config[ nodes[i].pops[j] ] );
	}
    }
  
}

void printAllPopsAtNodes(struct TREEB tree, struct TREEN *nodes)
{
  int i,j;
  
  for(i = 0; i < tree.nnode; ++i)
    {
      printf("***node %d contains %d populations: \t", i, nodes[i].npop);
      for( j = 0; j < nodes[i].npop; ++j)
	{
	  printf("%d, ", nodes[i].pops[j]);
	}
      printf("\n");
    }
}

void printAllSequencesPops(int nchrom, struct chromo *chrom)
{
  int i;
  for(i = 0; i < nchrom; ++i)
    {
      printf("@@@sequence %d belongs to population %d\n", i, chrom[i].pop);
    }
  
}

int chooseRandomPopulation ( int *array, int n)
{
  int i;

  i = rndu() * n;

  assert( i < n );

  return array[i];
  
}




struct segl * segtre_mig(struct c_params *cp, 
			 int *pnsegs, 
			 FILE *infofile,
			 int simindex
			 ) 
{

  int i, j, k, /*seg,*/ dec, pop, pop2, c1, c2, ind, rchrom, ic, is; //, eventOccured = 0; /*, intn  ; */
  int migrant, source_pop, *config; /*, flagint ;  */
  double  sum, x, /*tcoal,*/ ttemp, rft, clefta,  tmin = DBL_MAX, p  ;
  double prec, cin,  prect, /*nnm1,  nnm0,*/ mig, ran, coal_prob, /*prob,*/ rdum , arg ;
  char /*c,*/ event = 'x' ;
  int re(), re2(int nsam, int ic, int is), getSpot(), xover(), cinr(), cleftr(), eflag, cpop = 0;
  void getRecPosition(int nlinks, int *ic, int *is);
  int nsam, npop, nsites, /*nintn,*/ *inconfig, migrationEventsNumber = 0 ;
  double r,  f, rf,  track_len, /**nrec, *npast, *tpast,*/ **migm ;
  double *size, *alphag, *tlast, *samplingTime, *inSamplingTime ;
  struct devent *nextevent ;

  nsam = cp->nsam;
  npop = cp->npop;
  nsites = cp->nsites;
  inconfig = cp->config;
  inSamplingTime = cp->samplingTime;
  r = cp->r ;
  f = cp->f ;
  track_len = cp->track_len ;
  migm = (double **)malloc( (unsigned)npop*sizeof(double *) ) ;
  for( i=0; i<npop; i++) {
    migm[i] = (double *)malloc( (unsigned)npop*sizeof( double) ) ;
    for( j=0; j<npop; j++) migm[i][j] = (cp->mig_mat)[i][j] ;
  }
  nextevent = cp->deventlist ;



  /* printf("\nMigration Matrix\n"); */
  /* for( i=0; i<npop; i++) { */
    
  /*   for( j=0; j<npop; j++) */
  /*     printf("%f\t", migm[i][j]); */

  /*   printf("\n"); */
  /* } */
	
  /* Initialization */

  if( chrom == NULL ) {
    maxchr = nsam + 20 ;
    chrom = calloc( (unsigned)( maxchr), sizeof( struct chromo) ) ;
    if( chrom == NULL ) perror( "malloc error. segtre");
  }
  if( nnodes == NULL ){
    nnodes = (int*) malloc((unsigned)(seglimit*sizeof(int)))  ;
    if( nnodes == NULL ) perror("malloc error. segtre_mig");
  }
  if( seglst == NULL ) {
    seglst = (struct segl *)malloc((unsigned)(seglimit*sizeof(struct segl)) ) ;
    if( seglst == NULL ) perror("malloc error. segtre_mig.c 2");
  }

  config = (int *)malloc( (unsigned) ((npop+1)*sizeof(int) )) ;
  if( config == NULL ) perror("malloc error. segtre.");

  samplingTime = calloc( npop+1, sizeof(double) );
  if( samplingTime == NULL) perror("Malloc error... samplingTime.");

  size = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
  if( size == NULL ) perror("malloc error. segtre.");

  alphag = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
  if( alphag == NULL ) perror("malloc error. segtre.");

  tlast = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
  if( alphag == NULL ) perror("malloc error. segtre.");


  for(pop=0;pop<npop;pop++) {
    config[pop] = inconfig[pop] ;
    samplingTime[pop] = inSamplingTime[pop];
    size[pop] = (cp->size)[pop] ;
    alphag[pop] = (cp->alphag)[pop] ;
    tlast[pop] = 0.0 ;
  }


  for(pop=ind=0;pop<npop;pop++)
    for(j=0; j<inconfig[pop];j++,ind++) {
			
      chrom[ind].nseg = 1;
      if( !(chrom[ind].pseg = (struct seg*)malloc((unsigned)sizeof(struct seg)) ))
	ERROR("calloc error. se1");

      (chrom[ind].pseg)->beg = 0;
      (chrom[ind].pseg)->end = nsites-1;
      (chrom[ind].pseg)->desc = ind ;
      chrom[ind].pop = pop ;
    }
  

  seglst[0].beg = 0;
  
  if( !(seglst[0].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node)) ))
    perror("calloc error. se2");


  nnodes[0] = nsam - 1 ;
  nchrom=nsam;
  nlinks = ((long)(nsam))*(nsites-1) ;
  nsegs=1;
  t = 0.;
  r /= (nsites-1);
  if( f > 0.0 ) 	pc = (track_len -1.0)/track_len ;
  else pc = 1.0 ;
  lnpc = log( pc ) ;
  cleft = nsam* ( 1.0 - pow( pc, (double)(nsites-1) ) ) ;
  if( r > 0.0 ) rf = r*f ;
  else rf = f /(nsites-1) ;
  rft = rf*track_len ;
  /*flagint = 0 ; */

  /* Main loop */

  while( nchrom > 1 ) 
    {
      
      prec = nlinks*r;
      cin = nlinks*rf ;
      clefta = cleft*rft ;
      prect = prec + cin + clefta ;
      mig = 0.0;
      
      for( i=0; i<npop; i++) 
	{ 
	  mig += config[i]*migm[i][i] ;
	}
      
      if( (npop > 1) && ( mig == 0.0) && ( nextevent == NULL)) 
	{
	  i = 0;
	  for( j=0; j<npop; j++) 
	    if( config[j] > 0 ) i++;
	  if( i > 1 ) {
	    fprintf(stderr," Infinite coalescent time. No migration.\n");
	    exit(1);
	  }
	}
      
          
      eflag = 0 ;
      
      if( prect > 0.0 ) {      /* cross-over or gene conversion */

	while( (rdum = rndu() )  == 0.0 ) ;
	ttemp = -log( rdum)/prect ;
	if( (eflag == 0) || (ttemp < tmin ) ){
	  tmin = ttemp;
	  event = 'r' ;
	  eflag = 1;
	}
      }

      
      if( mig > 0.0 ){         /* migration   */
	
	while( (rdum = rndu() ) == 0.0 ) ;
	ttemp = -log( rdum)/mig ;
	
	if( (eflag == 0) || (ttemp < tmin ) ){
	  tmin = ttemp;
	  event = 'm' ;
	  eflag = 1 ;
	}
      }
      
      for(pop=0; pop<npop ; pop++) {     /* coalescent */
	coal_prob = ((double)config[pop])*(config[pop]-1.) ;
	if( coal_prob > 0.0 ) {
	  while( ( rdum = rndu() )  == .0 ) ;
	  if( alphag[pop] == 0 ){
	    ttemp = -log( rdum )*size[pop] /coal_prob ;
	    if( (eflag == 0) || (ttemp < tmin ) ){
	      tmin = ttemp;
	      event = 'c' ;
	      eflag = 1 ;
	      cpop = pop;
	    }
	  }
	  else {
	    arg  = 1. - alphag[pop]*size[pop]*exp(-alphag[pop]*(t - tlast[pop] ) )* log(rdum) / coal_prob     ;
	    if( arg > 0.0 ) {                          /*if arg <= 0,  no coalescent within interval */ 
	      ttemp = log( arg ) / alphag[pop]  ;
	      if( (eflag == 0) || (ttemp < tmin ) ){
		tmin = ttemp;
		event = 'c' ;
		eflag = 1 ;
		cpop = pop ;
	      }
	    }
	  }
	 
	}
	
      }// here we know what event will take place
      
      

      
      
      if( (eflag == 0) && ( nextevent == NULL) ) {
	fprintf(stderr,
		" infinite time to next event. Negative growth rate in last time interval or non-communicating subpops, tmin: %f, eflag: %d, nextevent: %p\n", tmin, eflag, nextevent);
	assert( 0);
      }
      
      if( ( ( eflag == 0) && (nextevent != NULL))|| ( (nextevent != NULL) &&  ( (t+tmin) >=  nextevent->time)) )  {
	
	t = nextevent->time ;

	/* printf("next event: %c\n", nextevent->detype); */
	
	switch(  nextevent->detype ) {
	  
	case 'N' :

	  for(pop =0; pop <npop; pop++)
	    {
	      
	      size[pop]= nextevent->paramv ;
	      
	      alphag[pop] = 0.0 ;
	    }
	  
	  nextevent = nextevent->nextde ;
	  
	  break;
	  
	case 'n' :
	  size[nextevent->popi]= nextevent->paramv ;
	  
	  alphag[nextevent->popi] = 0.0 ;
	  
	  nextevent = nextevent->nextde ;
	  
	  break;
	case 'G' :
	  for(pop =0; pop <npop; pop++)
	    {
	      size[pop] = size[pop]*exp( -alphag[pop]*(t - tlast[pop]) ) ;
	      alphag[pop]= nextevent->paramv ;
	      tlast[pop] = t ;
	    }
	  nextevent = nextevent->nextde ;
	  
	  break;
	  
	case 'g' :
	  pop = nextevent->popi ;
	  size[pop] = size[pop]*exp( - alphag[pop]*(t-tlast[pop]) ) ;
	  alphag[pop]= nextevent->paramv ;
	  tlast[pop] = t ;
	  nextevent = nextevent->nextde ;
	  break;
	case 'M' :
	  for(pop =0; pop <npop; pop++)
	    for( pop2 = 0; pop2 <npop; pop2++) migm[pop][pop2] = (nextevent->paramv) /(npop-1.0) ;
	  for( pop = 0; pop <npop; pop++)
	    migm[pop][pop]= nextevent->paramv ;
	  nextevent = nextevent->nextde ;
	  break;
	case 'a' :
	  for(pop =0; pop <npop; pop++)
	    for( pop2 = 0; pop2 <npop; pop2++)
	      {
		migm[pop][pop2] = (nextevent->mat)[pop][pop2]  ;
		
	      }
	  nextevent = nextevent->nextde ;
	  break;
	case 'm' :
	  i = nextevent->popi ;
	  j = nextevent->popj ;
	  migm[i][i] += nextevent->paramv - migm[i][j];
	  migm[i][j]= nextevent->paramv ;
	  nextevent = nextevent->nextde ;
	  break;
	case 'j' :         /* merge pop i into pop j  (join) */
	
	  i = nextevent->popi ;
	  
	  j = nextevent->popj ;

		  
	  config[j] += config[i] ;
	  config[i] = 0 ;
	  for( ic = 0; ic<nchrom; ic++) if( chrom[ic].pop == i ) chrom[ic].pop = j ;
	  /*  the following was added 19 May 2007 */
	  for( k=0; k < npop; k++){
	    if( k != i) {
	      migm[k][k] -= migm[k][i] ;
	      migm[k][i] = 0. ;
	    }
	  }

	  /* end addition */
	  nextevent = nextevent->nextde ;
	  break;
	case 's' :         /*split  pop i into two;p is the proportion from pop i, and 1-p from pop n+1  */
	  i = nextevent->popi ;
	  p = nextevent->paramv ;
	  npop++;
	  config = (int *)realloc( config, (unsigned)(npop*sizeof( int) )); 
	  size = (double *)realloc(size, (unsigned)(npop*sizeof(double) ));
	  alphag = (double *)realloc(alphag, (unsigned)(npop*sizeof(double) ));
	  tlast = (double *)realloc(tlast,(unsigned)(npop*sizeof(double) ) ) ;
	  tlast[npop-1] = t ;
	  size[npop-1] = 1.0 ;
	  alphag[npop-1] = 0.0 ;
	  migm = (double **)realloc(migm, (unsigned)(npop*sizeof( double *)));
	  for( j=0; j< npop-1; j++)
	    migm[j] = (double *)realloc(migm[j],(unsigned)(npop*sizeof(double)));
	  migm[npop-1] = (double *)malloc( (unsigned)(npop*sizeof( double) ) ) ;
	  for( j=0; j<npop; j++) migm[npop-1][j] = migm[j][npop-1] = 0.0 ;
	  config[npop-1] = 0 ;
	  config[i] = 0 ;
	  for( ic = 0; ic<nchrom; ic++){
	    if( chrom[ic].pop == i ) {
	      if( rndu() < p ) config[i]++;
	      else {
		chrom[ic].pop = npop-1 ;
		config[npop-1]++;
	      }
	    }
	  }
	  nextevent = nextevent->nextde ;
	  break;
	}
      }
      else {
	
	t += tmin ;

	

	//eventOccured = 0;
	
	if( event == 'r' ) {   

	  if( (ran = rndu()) < ( prec / prect ) ){ /*recombination*/
	    
	    getRecPosition(nlinks, &ic, &is);

	    
	    if( samplingTime[ chrom[ ic ].pop ] <= t )
	      {
		//eventOccured = 1;
		rchrom = re2(nsam, ic, is);
		//rchrom = re2(nsam, ic, is);	    	    
		config[ chrom[rchrom].pop ] += 1 ;
	      }
	  }
	  else if( ran < (prec + clefta)/(prect) ){    /*  cleft event */
	    rchrom = cleftr(nsam);
	    config[ chrom[rchrom].pop ] += 1 ;
	  }
	  else  {         /* cin event */
	    rchrom = cinr(nsam,nsites);
	    if( rchrom >= 0 ) config[ chrom[rchrom].pop ] += 1 ;
	  }
	}
	else if ( event == 'm' ) {  /* migration event */
	  x = mig*rndu();
	  sum = 0.0 ;
	  for( i=0; i<nchrom; i++) {
	    sum += migm[chrom[i].pop][chrom[i].pop] ;
	    if( x <sum ) break;
	  }

	  migrant = i ;

	  x = rndu()*migm[chrom[i].pop][chrom[i].pop];

	  sum = 0.0;
	  for(i=0; i<npop; i++){
	    if( i != chrom[migrant].pop ){
	      sum += migm[chrom[migrant].pop][i];
	      if( x < sum ) break;
	    }
	  }
	  source_pop = i;
	  
	  
	  if( samplingTime[ chrom[ migrant ].pop ] <= t &&
	      samplingTime[ source_pop ] <= t )
	    {
	      //eventOccured = 1;
	      if( infofile != NULL)
		fprintf(infofile, "sim: %d: migration from source %d to pop %d at time %f\n", simindex, source_pop, chrom[migrant].pop, t);
	 	      
	      config[chrom[migrant].pop] -= 1;
	      
	      config[source_pop] += 1;

	      migrationEventsNumber++;
	      
	      //fprintf(stderr, "source: %d, dest: %d time: %f\n", source_pop, chrom[migrant].pop, t);
	      
	      chrom[migrant].pop = source_pop ;
	    }
	}
	else { 								 /* coalescent event */
	  /* pick the two, c1, c2  */

	  /* printf("samplingTime of %d is %f\n", cpop, samplingTime[cpop]); */
	  
	  if( samplingTime[ cpop ] <= t )
	    {
	      
	      
	      pick2_chrom( cpop, config, &c1,&c2);  /* c1 and c2 are chrom's to coalesce */
	      dec = ca(nsam,nsites,c1,c2 );
	      config[cpop] -= dec ;
	    }
	}
      }

    }

  //fprintf(stdout, "number of migration events: %d\n", migrationEventsNumber);

  *pnsegs = nsegs ;
  free(config); 
  free( size ) ;
  free( alphag );
  free( tlast );
  for( i=0; i<npop; i++) free ( migm[i] ) ;
  free( migm ) ;
  return( seglst );
}


void getSequencesOnNodes(struct TREEN *phyloNodes, struct TREEB *tree, int *config, int npops)
{
  
  int i = 0, j = 0; 

  for( i = 0; i < tree->nnode; ++i)
    {
      phyloNodes[i].nseqs = 0;
      for ( j = 0; j < phyloNodes[i].npop; ++j)
	phyloNodes[i].nseqs += config[ phyloNodes[i].pops[j] ];
    }
    
}


int nodeOfPop(int pop, struct TREEN *nodes, int nnodes)
{
  int i, j;

  for(i = 0; i < nnodes; ++i)
    for( j = 0; j < nodes[i].npop; ++j)
      if( nodes[i].pops[j] == pop )
	return i;

  fprintf(stderr, "pop %d was not found.... error\n", pop);
  assert( i != nnodes );

  return 0;
}



struct segl * segtre_mig_multSpecies(struct c_params *cp, 
				     int *pnsegs,
				     struct TREEN *phyloNodes,
				     struct TREEB *treeStructure
				     
			 ) 
{


  /* phylogenetic variables */
  int totalNodes = treeStructure->nnode, 
    internalPhyloNodes = treeStructure -> internalNodes,
    *activePopsOnNodes = calloc( internalPhyloNodes, sizeof(int) );

  double duration = phyloNodes -> partialDuration, tis = 0.;;

  int partialIsActive = 0, activeNode = -1;;
  
    
  /*************************/


  int i, j, k, /*seg,*/ dec, pop, pop2, c1, c2, ind, rchrom, ic, is; //, eventOccured = 0; /*, intn  ; */
  int migrant, source_pop, *config; /*, flagint ;  */
  double  ran1(), sum, x, /*tcoal,*/ ttemp,ttemp2,  rft, clefta,  tmin = DBL_MAX, p  ;
  double prec, cin,  prect, /*nnm1,  nnm0,*/ mig, ran, coal_prob, /*prob,*/ rdum , arg ;
  char /*c,*/ event='x' ;
  int re(), re2(int nsam, int ic, int is), xover(), cinr(), cleftr(), eflag, cpop = 0;
  void getRecPosition(int nlinks, int *ic, int *is);
  int nsam, npop, nsites, /*nintn,*/ *inconfig ;
  double r,  f, rf,  track_len, /**nrec, *npast, *tpast,*/ **migm ;
  double *size, *alphag, *tlast, *samplingTime, *inSamplingTime ;
  struct devent *nextevent ;

  /* int *popsToNodes; */

  nsam = cp->nsam;
  
  npop = cp->npop;
  
  nsites = cp->nsites;
  
  inconfig = cp->config;
  
  inSamplingTime = cp->samplingTime;
  r = cp->r ;
  f = cp->f ;
  track_len = cp->track_len ;
  migm = (double **)malloc( (unsigned)npop*sizeof(double *) ) ;
  for( i=0; i<npop; i++) {
    migm[i] = (double *)malloc( (unsigned)npop*sizeof( double) ) ;
    for( j=0; j<npop; j++) migm[i][j] = (cp->mig_mat)[i][j] ;
  }
  nextevent = cp->deventlist ;
	
  /* Initialization */

  

  if( chrom == NULL ) {
    maxchr = nsam + 20 ;
    chrom = (struct chromo *)malloc( (unsigned)( maxchr*sizeof( struct chromo) )) ;
    if( chrom == NULL ) perror( "malloc error. segtre");
  }
  if( nnodes == NULL ){
    nnodes = (int*) malloc((unsigned)(seglimit*sizeof(int)))  ;
    if( nnodes == NULL ) perror("malloc error. segtre_mig");
  }
  if( seglst == NULL ) {
    seglst = (struct segl *)malloc((unsigned)(seglimit*sizeof(struct segl)) ) ;
    if( seglst == NULL ) perror("malloc error. segtre_mig.c 2");
  }

  config = (int *)malloc( (unsigned) ((npop+1)*sizeof(int) )) ;
  if( config == NULL ) perror("malloc error. segtre.");

  samplingTime = calloc( npop+1, sizeof(double) );

  if( samplingTime == NULL)
    perror( "Malloc error...sampling time\n\n");

  /* popsToNodes = calloc( (npop + 1), sizeof(int) ); */

  assert( config != NULL );

  size = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
  if( size == NULL ) perror("malloc error. segtre.");
  alphag = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
  if( alphag == NULL ) perror("malloc error. segtre.");
  tlast = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
  if( alphag == NULL ) perror("malloc error. segtre.");

  


  for(pop=0;pop<npop;pop++) {
    config[pop] = inconfig[pop] ;
    samplingTime[pop] = inSamplingTime[pop];
    size[pop] = (cp->size)[pop] ;
    alphag[pop] = (cp->alphag)[pop] ;
    tlast[pop] = 0.0 ;
  }


  for(pop=ind=0;pop<npop;pop++)
    for(j=0; j<inconfig[pop];j++,ind++) {
			
      chrom[ind].nseg = 1;
      if( !(chrom[ind].pseg = (struct seg*)malloc((unsigned)sizeof(struct seg)) ))
	ERROR("calloc error. se1");

      (chrom[ind].pseg)->beg = 0;
      (chrom[ind].pseg)->end = nsites-1;
      (chrom[ind].pseg)->desc = ind ;
      chrom[ind].pop = pop ;
    }
  seglst[0].beg = 0;
  if( !(seglst[0].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node)) ))
    perror("calloc error. se2");



  nnodes[0] = nsam - 1 ;
  nchrom=nsam;
  nlinks = ((long)(nsam))*(nsites-1) ;
  nsegs=1;
  t = 0.;
  r /= (nsites-1);
  if( f > 0.0 ) 	pc = (track_len -1.0)/track_len ;
  else pc = 1.0 ;
  lnpc = log( pc ) ;
  cleft = nsam* ( 1.0 - pow( pc, (double)(nsites-1) ) ) ;
  if( r > 0.0 ) rf = r*f ;
  else rf = f /(nsites-1) ;
  rft = rf*track_len ;
  /*flagint = 0 ; */

  
  partialIsActive = checkPartialMigration( phyloNodes, activePopsOnNodes, totalNodes, 0., duration, config );

  /* for(pop = 0; pop < npop; ++pop) */
  /*   { */
  /*     popsToNodes[pop] = nodeOfPop(pop, phyloNodes, treeStructure->nnode); */
  /*   } */

  getSequencesOnNodes(phyloNodes, treeStructure, config, npop);
  
  for( i = 0; i < treeStructure->nnode; ++i)
    phyloNodes[i].nseqs = 0;


  /* Main loop */

  while( nchrom > 1 ) 
    {
      
      getSequencesOnNodes(phyloNodes, treeStructure, config, npop);
  
      /* printAllPopsAtNodes(tree, phyloNodes); */

      
      /* printAllSequencesPops(nchrom, chrom); */
      
      /* printConfig(config, tree, phyloNodes); */
      
      prec = nlinks*r;
      cin = nlinks*rf ;
      clefta = cleft*rft ;
      prect = prec + cin + clefta ;
      mig = 0.0;
      
      for( i=0; i<npop; i++) 
	{ 
	  mig += config[i]*migm[i][i] ;
	}
      
      if( (npop > 1) && ( mig == 0.0) && ( nextevent == NULL)) 
	{
	  i = 0;
	  for( j=0; j<npop; j++) 
	    if( config[j] > 0 ) i++;
	  if( i > 1 ) {
	    fprintf(stderr," Infinite coalescent time. No migration.\n");
	    exit(1);
	  }
	}
      
    
      
      eflag = 0 ;
      
      if( prect > 0.0 ) {      /* cross-over or gene conversion */
	while( (rdum = rndu() )  == 0.0 ) ;
	ttemp = -log( rdum)/prect ;
	if( (eflag == 0) || (ttemp < tmin ) ){
	  tmin = ttemp;
	  event = 'r' ;
	  eflag = 1;
	}
      }
      if( partialIsActive == 1 && /* partial migration */
	  phyloNodes -> partialMaxMigration > 0 && 
	  duration > 0)
	{
	  ttemp = DBL_MAX;
	  
	  for( i = 0; i < internalPhyloNodes; ++i)
	    {
	      
	      j = i + treeStructure -> leaves;

	      if( activePopsOnNodes[i] == 0)
		continue;
	      
	      tis = phyloNodes[j].isolationTime; 

	      	      
	      if( randPoissonLinear( tis, phyloNodes[j].coalAge, t, activePopsOnNodes[i] * phyloNodes -> partialMaxMigration, &ttemp2) == 1  &&
		  ttemp2 < ttemp)
		{
		  
		  ttemp = ttemp2;
		  activeNode = j;
		}
	      
	    }

	  if( eflag == 0 || ttemp < tmin )
	    {
	      tmin = ttemp;
	      event = 'i';
	      eflag = 1;
	    }
	}
      

      if( mig > 0.0 ){         /* migration   */
	
	while( (rdum = rndu() ) == 0.0 ) ;
	
	ttemp = -log( rdum)/mig ;
	
	if( (eflag == 0) || (ttemp < tmin ) )
	  {
	    
	    tmin = ttemp;
	    event = 'm' ;
	    eflag = 1 ;
	  }
      }
      
      for(pop=0; pop<npop ; pop++) {     /* coalescent */
	coal_prob = ((double)config[pop])*(config[pop]-1.) ;
	if( coal_prob > 0.0 ) {
	  while( ( rdum = rndu() )  == .0 ) ;
	  if( alphag[pop] == 0 ){
	    ttemp = -log( rdum )*size[pop] /coal_prob ;
	    if( (eflag == 0) || (ttemp < tmin ) ){
	      tmin = ttemp;
	      event = 'c' ;
	      eflag = 1 ;
	      cpop = pop;
	    }
	  }
	  else {
	    arg  = 1. - alphag[pop]*size[pop]*exp(-alphag[pop]*(t - tlast[pop] ) )* log(rdum) / coal_prob     ;
	    if( arg > 0.0 ) {                          /*if arg <= 0,  no coalescent within interval */ 
	      ttemp = log( arg ) / alphag[pop]  ;
	      if( (eflag == 0) || (ttemp < tmin ) ){
		tmin = ttemp;
		event = 'c' ;
		eflag = 1 ;
		cpop = pop ;
	      }
	  }
	  }
	}		
      }


      
      if( (eflag == 0) && ( nextevent == NULL) ) {
	fprintf(stderr,
		" infinite time to next event. Negative growth rate in last time interval or non-communicating subpops.\n");
	assert(0);
      }
      
      if( ( ( eflag == 0) && (nextevent != NULL))|| ( (nextevent != NULL) &&  ( (t+tmin) >=  nextevent->time)) )  {
	
	t = nextevent->time ;
	
	switch(  nextevent->detype ) {

	  printf("Next event is %c\n", nextevent->detype);
	  
	case 'i' : /* new event, partial isolation after (forward) a speciation event */

	  partialIsActive = checkPartialMigration( phyloNodes, activePopsOnNodes, totalNodes, t, duration, config);

	  
	  for( i = 0; i < internalPhyloNodes; ++i)
	    {
	      
	      j = i + treeStructure -> leaves;

	    }
	  
	  nextevent = nextevent -> nextde;
	  
	  break;
	  
	case 'N' :

	  for(pop =0; pop <npop; pop++)
	    {
	      
	      size[pop]= nextevent->paramv ;
	      
	      alphag[pop] = 0.0 ;
	    }
	  
	  nextevent = nextevent->nextde ;
	  
	  break;
	  
	case 'n' :
	  size[nextevent->popi]= nextevent->paramv ;
	  
	  alphag[nextevent->popi] = 0.0 ;
	  
	  nextevent = nextevent->nextde ;
	  
	  break;
	case 'G' :
	  for(pop =0; pop <npop; pop++)
	    {
	      size[pop] = size[pop]*exp( -alphag[pop]*(t - tlast[pop]) ) ;
	      alphag[pop]= nextevent->paramv ;
	      tlast[pop] = t ;
	    }
	  nextevent = nextevent->nextde ;
	  
	  break;
	  
	case 'g' :
	  pop = nextevent->popi ;
	  size[pop] = size[pop]*exp( - alphag[pop]*(t-tlast[pop]) ) ;
	  alphag[pop]= nextevent->paramv ;
	  tlast[pop] = t ;
	  nextevent = nextevent->nextde ;
	  break;
	case 'M' :
	  for(pop =0; pop <npop; pop++)
	    for( pop2 = 0; pop2 <npop; pop2++) 
	      {
		
		migm[pop][pop2] = (nextevent->paramv) /(npop-1.0) ;
		//printf("migration at %f is %f\n", t, migm[pop][pop2]);
	      }
	  for( pop = 0; pop <npop; pop++)
	    {
	      
	      migm[pop][pop]= nextevent->paramv ;
	      
	    }
	  nextevent = nextevent->nextde ;
	  break;
	case 'a' :
	  for(pop =0; pop <npop; pop++)
	    for( pop2 = 0; pop2 <npop; pop2++) migm[pop][pop2] = (nextevent->mat)[pop][pop2]  ;
	  nextevent = nextevent->nextde ;
	  break;
	case 'm' :
	  
	  if(partialIsActive == 1)
	    {
	      assert(0);
	    }
	  else if( mig > 0)
	    {

	      i = nextevent->popi ;
	      j = nextevent->popj ;
	      migm[i][i] += nextevent->paramv - migm[i][j];
	      migm[i][j]= nextevent->paramv ;
	      nextevent = nextevent->nextde ;
	    }
	  break;
	case 'j' :         /* merge pop i into pop j  (join) */
	  
	  if( nextevent -> isSpeciation == 1)
	    {
	      partialIsActive = checkPartialMigration( phyloNodes, activePopsOnNodes, totalNodes, t, duration, config);
	      
	      
	      for( i = 0; i < internalPhyloNodes; ++i)
		{
		  
		  j = i + treeStructure -> leaves;
		  
		}


	    }

	  i = nextevent->popi ;
	  
	  j = nextevent->popj ;
		  
	  config[j] += config[i] ;

	  /* phyloNodes[ popsToNodes[j] ].nseqs += phyloNodes[ popsToNodes[ i ] ].nseqs; */
	  
	  config[i] = 0 ;

	  /* phyloNodes[ popsToNodes[i] ].nseqs = 0; */
	  
	  for( ic = 0; ic<nchrom; ic++) 
	    if( chrom[ic].pop == i ) 
	      chrom[ic].pop = j ;
	  
	  /*  the following was added 19 May 2007 */
	  for( k=0; k < npop; k++){
	    if( k != i) {
	      migm[k][k] -= migm[k][i] ;
	      migm[k][i] = 0. ;
	    }
	  }

	  /* end addition */
	  nextevent = nextevent->nextde ;
	  
	  break;
	  
	case 's' : 
	  
	  /*split  pop i into two;p is the proportion from pop i, and 1-p from pop n+1  */
	  fprintf(stderr, "split pops not implemented yet... exit\n");
	  assert(0);
	  
	  i = nextevent->popi ;
	  p = nextevent->paramv ;
	  npop++;
	  config = (int *)realloc( config, (unsigned)(npop*sizeof( int) )); 
	  size = (double *)realloc(size, (unsigned)(npop*sizeof(double) ));
	  alphag = (double *)realloc(alphag, (unsigned)(npop*sizeof(double) ));
	  tlast = (double *)realloc(tlast,(unsigned)(npop*sizeof(double) ) ) ;
	  tlast[npop-1] = t ;
	  size[npop-1] = 1.0 ;
	  alphag[npop-1] = 0.0 ;
	  migm = (double **)realloc(migm, (unsigned)(npop*sizeof( double *)));
	  for( j=0; j< npop-1; j++)
	    migm[j] = (double *)realloc(migm[j],(unsigned)(npop*sizeof(double)));
	  migm[npop-1] = (double *)malloc( (unsigned)(npop*sizeof( double) ) ) ;
	  for( j=0; j<npop; j++) migm[npop-1][j] = migm[j][npop-1] = 0.0 ;
	  config[npop-1] = 0 ;
	  config[i] = 0 ;
	  for( ic = 0; ic<nchrom; ic++){
	    if( chrom[ic].pop == i ) {
	      if( rndu() < p ) config[i]++;
	      else {
		chrom[ic].pop = npop-1 ;
		config[npop-1]++;
	      }
	    }
	  }
	  nextevent = nextevent->nextde ;
	  break;
	}
      }
      
      else {
	
	t += tmin ;	
	
	//eventOccured = 0;

	if( event == 'r' ) {   
	  if( (ran = rndu()) < ( prec / prect ) ){ /*recombination*/
	    
	    getRecPosition(nlinks, &ic, &is);

	    if( samplingTime[ chrom[ ic ].pop ] <= t )
	      {
		//eventOccured = 1;
		rchrom = re2(nsam, ic, is);
		config[ chrom[rchrom].pop ] += 1 ;
	      }
	    
	    //rchrom = re(nsam);
	    
	    /* phyloNodes[ popsToNodes[ chrom[rchrom].pop ] ].nseqs += 1; */
	    
	  }
	  else if( ran < (prec + clefta)/(prect) ){    /*  cleft event */
	    rchrom = cleftr(nsam);
	    
	    config[ chrom[rchrom].pop ] += 1 ;
	    /* phyloNodes[ popsToNodes[ chrom[rchrom].pop ] ].nseqs += 1; */
	  }
	  else  {         /* cin event */
	    rchrom = cinr(nsam,nsites);
	    if( rchrom >= 0 ) 
	      {
		config[ chrom[rchrom].pop ] += 1 ;
		/* phyloNodes[ popsToNodes[ chrom[rchrom].pop ] ].nseqs += 1; */
	      }
	  }
	}
	else if ( event == 'm' ) {  /* migration event */
	  
	  x = mig*rndu();
	  
	  sum = 0.0 ;
	  for( i=0; i<nchrom; i++) {
	    
	    sum += migm[chrom[i].pop][chrom[i].pop] ;
	    
	    if( x <sum ) break;
	  }
	  
	  

	  migrant = i ;
	  x = rndu()*migm[chrom[i].pop][chrom[i].pop];
	  sum = 0.0;
	  for(i=0; i<npop; i++){
	    if( i != chrom[migrant].pop ){
	      sum += migm[chrom[migrant].pop][i];
	      if( x < sum ) break;
	    }
	  }
	  source_pop = i;

	  if( samplingTime[ chrom[ migrant ].pop ] <= t &&
	      samplingTime[ source_pop ] <= t )
	    {
	      //eventOccured = 1;
	  
	      config[chrom[migrant].pop] -= 1;

	      /* phyloNodes[ popsToNodes[ chrom[migrant].pop ] ].nseqs -= 1; */
	      
	      config[source_pop] += 1;
	      
	      /* phyloNodes[ popsToNodes[ source_pop ] ].nseqs += 1; */
	      
	      chrom[migrant].pop = source_pop ;
	      
	      partialIsActive = checkPartialMigration( phyloNodes, activePopsOnNodes, totalNodes, t, duration, config);
	      
	      
	      for( i = 0; i < internalPhyloNodes; ++i)
		{
		  
		  j = i + treeStructure -> leaves;
		  
		}
	      
	    }
	}
	
	else if( event == 'i' ) 
	  {
	    
	    int seqPop;

	    assert(activeNode > -1);
	    
	    assert(activeNode < treeStructure -> nnode);
	    	  
	    /* choose a sequence from this node */
	    migrant = chooseRandomSequenceAtNode( phyloNodes[ activeNode ], &seqPop, config, chrom, nchrom);
	    
	    source_pop = chooseRandomPopulation( phyloNodes[ activeNode ].pops, phyloNodes[ activeNode ].npop);

	    
	    /* check whether the event has taken place after sampling */
	    if( samplingTime[ source_pop ] <= t &&
		samplingTime[ chrom[ migrant ].pop ] <= t)
	      {
    
		//eventOccured = 1;
		
		config[ chrom[migrant].pop ] -= 1;

		
		/* phyloNodes[ popsToNodes[ chrom[migrant].pop ] ].nseqs -= 1; */
		
		config[ source_pop ] += 1;
	    
	    
		/* phyloNodes[ popsToNodes[ source_pop ] ].nseqs += 1; */
		
		chrom[ migrant].pop = source_pop;
		
		partialIsActive = checkPartialMigration( phyloNodes, activePopsOnNodes, totalNodes, t, duration, config);
		
		
		for( i = 0; i < internalPhyloNodes; ++i)
		  {
		    
		    j = i + treeStructure -> leaves;
		    
		  }
		
		assert( config[ chrom[migrant].pop ] >= 0 );
	    
	      }
	  }
	else { 	
	  /* coalescent event */
	  
	  /* pick the two, c1, c2  */
	  
	  /* printAllSequencesPops(nchrom, chrom); */

	  /* printConfig(config, tree, phyloNodes); */

	  /* printf("samplingTime of %d is %f\n", cpop, samplingTime[cpop]); */

	  if( samplingTime[ cpop ] <= t )
	    {
	      
	      /* c1 and c2 are chrom's to coalesce */
	      pick2_chrom( cpop, config, &c1,&c2);  
	      
	      assert( c1 < nchrom);
	      
	      assert( c2 < nchrom);
	      
	      dec = ca(nsam,nsites,c1,c2 );
	      
	      config[cpop] -= dec ;
	      
	      /* phyloNodes[ popsToNodes[ cpop ] ].nseqs -= dec; */
	    }
	  
	}

      }
      
      
      /* printAllSequencesPops(nchrom, chrom); */


    }
  
  *pnsegs = nsegs ;
  free(config); 
  free( size ) ;
  free( alphag );
  free( tlast );
  for( i=0; i<npop; i++) free ( migm[i] ) ;
  free( migm ) ;
  free(activePopsOnNodes);
  return( seglst );
}




/******  recombination subroutine ***************************
	 Picks a chromosome and splits it in two parts. If the x-over point
	 is in a new spot, a new segment is added to seglst and a tree set up
	 for it.   ****/


int
re(nsam)
     int nsam;
{
  struct seg *pseg = NULL;
  
  int  el, lsg, lsgm1,  ic,  is, spot;
    /* First generate a random x-over spot, then locate it as to chrom and seg. */

  spot = nlinks*rndu() + 1.;
  
  /* get chromosome # (ic)  */

  for( ic=0; ic<nchrom ; ic++) {
    lsg = chrom[ic].nseg ;
    lsgm1 = lsg - 1;
    pseg = chrom[ic].pseg;
    el = ( (pseg+lsgm1)->end ) - (pseg->beg);
    if( spot <= el ) break;
    spot -= el ;
  }
  
  assert(pseg != NULL);

  is = pseg->beg + spot -1;
  xover(nsam, ic, is);
  return(ic);	
}

void getRecPosition(int nlinks, int *ic, int *is)
{
  
  int  el, lsg, lsgm1,  spot;

  struct seg *pseg = NULL;

  spot = nlinks*rndu() + 1;

  /* get chromosome # (ic)  */

  for( *ic=0; *ic<nchrom ; (*ic)++) {
    lsg = chrom[*ic].nseg ;
    lsgm1 = lsg - 1;
    pseg = chrom[*ic].pseg;
    el = ( (pseg+lsgm1)->end ) - (pseg->beg);
    if( spot <= el ) break;
    spot -= el ;
  }
  
  assert(pseg != NULL);

  *is = pseg->beg + spot -1;

}


int re2(int nsam, int ic, int is)
{
  
  xover(nsam, ic, is);

  return(ic);	
}


int
cleftr( int nsam)
{
  struct seg *pseg ;
  int   /*lsg, lsgm1,*/  ic,  is; /*, in, spot; */
  double x, sum, len  ;
  int xover(int, int, int);

  while( (x = cleft*rndu() )== 0.0 )  ;
  sum = 0.0 ;
  ic = -1 ;
  while ( sum < x ) {
    sum +=  1.0 - pow( pc, links(++ic) )  ;
  }
  pseg = chrom[ic].pseg;
  len = links(ic) ;
  is = pseg->beg + floor( 1.0 + log( 1.0 - (1.0- pow( pc, len))*rndu() )/lnpc  ) -1  ;
  xover( nsam, ic, is);
  return( ic) ;
}

int
cinr( int nsam, int nsites)
{
  struct seg *pseg = NULL;
  int len,  el, lsg, lsgm1,  ic,  is, spot, endic ;
  double rndu();
  int xover(), ca() ;


  /* First generate a random x-over spot, then locate it as to chrom and seg. */

  spot = nlinks*rndu() + 1.;

  /* get chromosome # (ic)  */

  for( ic=0; ic<nchrom ; ic++) {
    lsg = chrom[ic].nseg ;
    lsgm1 = lsg - 1;
    pseg = chrom[ic].pseg;
    el = ( (pseg+lsgm1)->end ) - (pseg->beg);
    if( spot <= el ) break;
    spot -= el ;
  }

  assert(pseg != NULL);

  is = pseg->beg + spot -1;
  endic = (pseg+lsgm1)->end ;
  xover(nsam, ic, is);

  len = floor( 1.0 + log( rndu() )/lnpc ) ;
  if( is+len >= endic ) return(ic) ;  
  if( is+len < (chrom[nchrom-1].pseg)->beg ){
    ca( nsam, nsites, ic, nchrom-1);
    return(-1) ;
  }
  xover( nsam, nchrom-1, is+len ) ;
  ca( nsam,nsites, ic,  nchrom-1);
  return(ic);	

}

int xover(int nsam,int ic, int is)
{
  struct seg *pseg, *pseg2;
  int i,  lsg, lsgm1, newsg,  jseg, k,  in;
  double len ;


  pseg = chrom[ic].pseg ;
  lsg = chrom[ic].nseg ;
  len = (pseg + lsg -1)->end - pseg->beg ;
  cleft -= 1 - pow(pc,len) ;
  /* get seg # (jseg)  */

  for( jseg=0; is >= (pseg+jseg)->end ; jseg++) ;
  if( is >= (pseg+jseg)->beg ) in=1;
  else in=0;
  newsg = lsg - jseg ;

  /* copy last part of chrom to nchrom  */

  nchrom++;
  if( nchrom >= maxchr ) {
    maxchr += 20 ;
    chrom = (struct chromo *)realloc( chrom, (unsigned)(maxchr*sizeof(struct chromo))) ;
    if( chrom == NULL ) perror( "malloc error. segtre2");
  }
  if( !( pseg2 = chrom[nchrom-1].pseg = (struct seg *)calloc((unsigned)newsg,sizeof(struct seg)) ) )
    ERROR(" alloc error. re1");
  chrom[nchrom-1].nseg = newsg;
  chrom[nchrom-1].pop = chrom[ic].pop ;
  pseg2->end = (pseg+jseg)->end ;
  if( in ) {
    pseg2->beg = is + 1 ;
    (pseg+jseg)->end = is;
  }
  else pseg2->beg = (pseg+jseg)->beg ;
  pseg2->desc = (pseg+jseg)->desc ;
  for( k=1; k < newsg; k++ ) {
    (pseg2+k)->beg = (pseg+jseg+k)->beg;
    (pseg2+k)->end = (pseg+jseg+k)->end;
    (pseg2+k)->desc = (pseg+jseg+k)->desc;
  }

  lsg = chrom[ic].nseg = lsg-newsg + in ;
  lsgm1 = lsg - 1 ;
  nlinks -= pseg2->beg - (pseg+lsgm1)->end ;
  len = (pseg+lsgm1)->end - (pseg->beg) ;
  cleft += 1.0 - pow( pc, len) ;
  len = (pseg2 + newsg-1)->end - pseg2->beg ;
  cleft += 1.0 - pow(pc, len) ;
  if( !(chrom[ic].pseg = 
	(struct seg *)realloc(chrom[ic].pseg,(unsigned)(lsg*sizeof(struct seg)) )) )
    perror( " realloc error. re2");
  if( in ) {
    begs = pseg2->beg;
    for( i=0,k=0; (k<nsegs-1)&&(begs > seglst[seglst[i].next].beg-1);
	 i=seglst[i].next, k++) ;
    if( begs != seglst[i].beg ) {
      /* new tree  */

      if( nsegs >= seglimit ) {  
	seglimit += SEGINC ;
	nnodes = (int *)realloc( nnodes,(unsigned)(sizeof(int)*seglimit)) ; 
	if( nnodes == NULL) perror("realloc error. 1. segtre_mig.c");
	seglst =
	  (struct segl *)realloc( seglst,(unsigned)(sizeof(struct segl)*seglimit));
	if(seglst == NULL ) perror("realloc error. 2. segtre_mig.c");

      } 
      seglst[nsegs].next = seglst[i].next;
      seglst[i].next = nsegs;
      seglst[nsegs].beg = begs ;
      if( !(seglst[nsegs].ptree = (struct node *)calloc((unsigned)(2*nsam), sizeof(struct
										   node)) )) perror("calloc error. re3.");
      nnodes[nsegs] = nnodes[i];
      ptree1 = seglst[i].ptree;
      ptree2 = seglst[nsegs].ptree;
      nsegs++ ;
      for( k=0; k<=nnodes[i]; k++) {
	(ptree2+k)->abv = (ptree1+k)->abv ;
	(ptree2+k)->time = (ptree1+k)->time;
      }
    }
  }
  return(ic) ;
}

/***** common ancestor subroutine **********************
       Pick two chromosomes and merge them. Update trees if necessary. **/

int ca(int nsam, int nsites, int c1, int c2)
{
  int yes1, yes2, seg1, seg2, seg ;
  int tseg, start, end, desc, k;
  struct seg *pseg;
  struct node *ptree;

  seg1=0;
  seg2=0;

  if( !(pseg = (struct seg *)calloc((unsigned)nsegs,sizeof(struct seg) ))) 
    perror("alloc error.ca1");

  tseg = -1 ;

  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
    start = seglst[seg].beg;

    yes1 = isseg(start, c1, &seg1);
    yes2 = isseg(start, c2, &seg2);
    if( yes1 || yes2 ) {
      tseg++;
      (pseg+tseg)->beg=seglst[seg].beg;
      end = ( k< nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1 ) ;
      (pseg+tseg)->end = end ;

      if( yes1 && yes2 ) {
	nnodes[seg]++;
	if( nnodes[seg] >= (2*nsam-2) ) tseg--;
	else
	  (pseg+tseg)->desc = nnodes[seg];
	ptree=seglst[seg].ptree;
	desc = (chrom[c1].pseg + seg1) ->desc;
	(ptree+desc)->abv = nnodes[seg];
	desc = (chrom[c2].pseg + seg2) -> desc;
	(ptree+desc)->abv = nnodes[seg];
	(ptree+nnodes[seg])->time = t;

      }
      else {
	(pseg+tseg)->desc = ( yes1 ?
			      (chrom[c1].pseg + seg1)->desc :
			      (chrom[c2].pseg + seg2)->desc);
      }
    }
  }
  nlinks -= links(c1);
  
  cleft -= 1.0 - pow(pc, (double)links(c1));
    
  free(chrom[c1].pseg) ;
  
  if( tseg < 0 ) {
    free(pseg) ;
    chrom[c1].pseg = chrom[nchrom-1].pseg;
    chrom[c1].nseg = chrom[nchrom-1].nseg;
    chrom[c1].pop = chrom[nchrom-1].pop ;
    if( c2 == nchrom-1 ) c2 = c1;
    nchrom--;
  }
  else {
    if( !(pseg = realloc(pseg,(unsigned)((tseg+1)*sizeof(struct seg)))))
      perror(" realloc error. ca1");
    chrom[c1].pseg = pseg;
    chrom[c1].nseg = tseg + 1 ;
    nlinks += links(c1);
    cleft += 1.0 - pow(pc, (double)links(c1));
  }
  nlinks -= links(c2);
  
  cleft -= 1.0 - pow(pc, (double)links(c2));
  
  free(chrom[c2].pseg) ;
  
  chrom[c2].pseg = chrom[nchrom-1].pseg;
  chrom[c2].nseg = chrom[nchrom-1].nseg;
  chrom[c2].pop = chrom[nchrom-1].pop ;
  nchrom--;
  if(tseg<0) return( 2 );  /* decrease of nchrom is two */
  else return( 1 ) ;
}

/*** Isseg: Does chromosome c contain the segment on seglst which starts at
     start? *psg is the segment of chrom[c] at which one is to begin 
     looking.  **/

int isseg(int start, int c, int *psg)
{
  int ns;
  struct seg *pseg;

  ns = chrom[c].nseg;
  pseg = chrom[c].pseg;

  

  /*  changed order of test conditions in following line on 6 Dec 2004 */
  for(  ; ((*psg) < ns ) && 
	  ( (pseg+(*psg))->beg <= start ) ; 
	++(*psg) )
    {
      if( (pseg+(*psg))->end >= start ) return(1);
    }
  return(0);
}
	


int pick2_chrom(int pop, int config[], int *pc1, int *pc2)
{
  int c1, c2, cs,cb,i, count;
	
  pick2(config[pop],&c1,&c2);

  cs = (c1>c2) ? c2 : c1;
  cb = (c1>c2) ? c1 : c2 ;
  i=count=0;

  for(;;){
    while( chrom[i].pop != pop ) i++;
    if( count == cs ) break;
    count++;
    i++;
  }
  *pc1 = i;
  i++;
  count++;
  for(;;){
    while( chrom[i].pop != pop ) i++;
    if( count == cb ) break;
    count++;
    i++;
  }
  *pc2 = i ;

  return 1;
}
	
	

/****  links(c): returns the number of links between beginning and end of chrom **/

int links(int c)
{
  int ns;

  ns = chrom[c].nseg - 1 ;

  return( (chrom[c].pseg + ns)->end - (chrom[c].pseg)->beg);
}

