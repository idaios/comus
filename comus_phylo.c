#include <float.h>
#include "comus.h"
#include "twister.h"

#define cntMAX 10000


#define NAMESIZE 100

#define TOKENSIZE 10024

#define NAMESIZE 100



/*******************************/


char readChar(FILE *stream, int *depth)
{
    char chr;
    
    do {
        if (fread(&chr, sizeof(char), 1, stream) != 1) {
	  return '\0';
        }
    } while (chr == ' ' || chr == '\n' || chr == 9 || chr == 13);
    
    // keep track of paren depth
    if (chr == '(') (*depth)++;
    if (chr == ')') (*depth)--;

    return chr;
}



char readLastChar(FILE *stream, int *depth)
{
  char chr;
  long int lastPos = 0;
  do {
    lastPos = ftell(stream);    
    if (fread(&chr, sizeof(char), 1, stream) != 1) {
      return '\0';
    }

    if( chr == ':' || chr == ';' ) // go back and exit
      {
	fseek(stream, lastPos, SEEK_SET);
	break;
      }

    // just exit
    if( chr == ')' || chr == '('){
      break;
    }
    
  } while (chr == ' ' || chr == '\n' || chr == 9 || chr == 13 || (chr > 64 && chr < 91) || // capital letter
	   (chr > 96 && chr < 123) || // small letter
	   (chr > 47 && chr < 58) || //integer
	   chr == 95 || chr == 46  ); // period, underscore

   // keep track of paren depth
  if (chr == '(') (*depth)++;
  if (chr == ')') (*depth)--;
  
  return chr;
}



char readUntil(FILE *stream, char token[TOKENSIZE], const char *stops, int *depth)
{
    char chr;

    int i = 1;

    const char *j = stops;

    token[0] = token[1] = '\0';
        
    while (1) {
      
        chr = readChar(stream, depth);

	/* if end character */
        if (!chr)
            return chr;
        
        // compare char to stop characters
	while( *j )
	  {
	  
	    if (chr == *j)
	      return chr;
	    
	    ++i;
	    ++j;
        }

	assert( i + 1 < TOKENSIZE );
	
	token[i] = chr;
	
	token[++i] = '\0';
    }
}


int nodeTimeNameCmp(const void *_a, const void *_b)
{
    struct TREEN *a = (struct TREEN*) _a;
    struct TREEN *b = (struct TREEN*) _b;

    /* printf("a->age: %f, b->age: %f\n", a->age, b->age); */
    
    if(a->age < b->age)
      return -1;
    
    else if(a->age > b->age)
      return 1;
    
    else if(a->age != 0. && b->age != 0.)
      return 0;
    
    else if(a->age == 0 && b->age == 0.)
      {
	if(a->nodeID < b->nodeID)
	  return -1;
	else if(a->nodeID > b->nodeID)
	  return 1;
	else
	  return 0;
      }
    
    return 0;
}


int readNumberOfReplications(FILE *infile)
{
  char chr;
  int nrepl = 0;

  chr = fgetc(infile);
  
  while( chr == 9 || chr == 32 || chr == 13 || chr == 10 )
    chr = fgetc(infile);

  if(chr == EOF)
    {
      /* go one position back. The readTree function will then report the enfOfFile */
      ungetc(chr, infile);
      return 0;
    }

  if(chr != '[' && chr != '(')
    {
      fprintf(stderr, "Either [ or ( is expected; however %c was found.\n", chr);
      assert(chr == '[' || chr == '(' );
    }

  if( chr == '[' )
    {
      
      if( fscanf( infile, "%d", &nrepl ) != 1 )
	{
	  fprintf( stderr, "\nError while expecting to read the number of replications for the tree\n\n");
	  assert( 0 );
	}

      while( (chr = fgetc(infile) ) != '(')
	{
	  if(chr != ']' && chr != 32 && chr != 9 && chr != 10 && chr != 13)
	    {
	      fprintf( stderr, "\n] or SPACE or TAB is expected before the ( symbol. However %c was found\n\n", chr );

	      assert(0);
	    }
	}
    }

  else
    {
      assert(chr == '(');
      nrepl = 0;
    }

  if(chr != '(' )
    {
      fprintf(stderr, "character ( is expected but %c was read\n\n", chr);
      assert(chr == '(');
    }

  assert( chr == '(');

  /* go one position back so that next function will read ( for sure */
  ungetc(chr, infile);
  
  return nrepl;

}


/* code from Fletcher and Yang's INDELible */
void ClearNode (int inode)
{
  
  nodes[inode].father=nodes[inode].ibranch=-1;

  nodes[inode].nson=0;
  
  nodes[inode].branch=nodes[inode].age=0;
  
}


float readDist(FILE *infile)
{
    float dist = 0;
    
    fscanf(infile, "%f", &dist);
    
    return dist;
}

int assignChildrenToNodes(struct TREEN *nodes, int nnodes, int *positionsArray)
{
  int i, k;
  int *index = calloc( nnodes, sizeof(int) );

  assert(index != NULL);
  
  for(i = 0; i < nnodes; ++i)
    {
      if(nodes[i].father > -1)
	{
	  assert(nodes[i].father < nnodes);
	  k = nodes[i].father;
	  nodes[i].father = positionsArray[k];
	  //	  fprintf(stderr, "i: %d, nodes[i].father: %d\n", i, nodes[i].father);
	  nodes[ nodes[i].father ].sons[index[ nodes[i].father ] ] = i;
	  index[ nodes[i].father ]++;
	}
    }
  free( index );
  return 1;
}


int readNewickNode(FILE *infile, int parentIndex, int *depth, int nspecies, int *nodesInFile)
{

  static int node = -1;

  if(parentIndex == -1)
    node = -1;

  int currentNode; 
  
  int id2;
  
  char chr, chr1;

  int depth2;

  char token[TOKENSIZE];
  double dist = 0.;
  
  /* read the next character and keep track of the parenthesis depth */
  if (!(chr1  = readLastChar(infile, depth))) {
    fprintf(stderr, "unexpected end of file");
    assert(0);
  }
  /* chr1 can be either (, or ) */

  assert( *depth >= 0 );

  
  if(chr1 == '(')
    {

      ++(*nodesInFile);

      node++;

      currentNode = node;
      
      depth2 = *depth;

      nodes[currentNode].father = parentIndex;

      nodes[currentNode].originalIndex = node;

      /* nodes[currentNode].nodeStr = calloc(256, sizeof(char)); */
      /* strcpy(nodes[currentNode].nodeStr, "internal"); */
      
      nodes[currentNode].nodeID = -1;
      
      nodes[currentNode].nson = 0;
      
      if(parentIndex > -1)      
	++nodes[parentIndex].nson;
      
      nodes[currentNode].age = 0.;

      while( *depth == depth2 )
	{
	  
	  id2 = readNewickNode(infile, currentNode, depth, nspecies, nodesInFile);
	  
	  if(id2 < 0)
	    {
	      fprintf(stderr, "Warning ... comus.phylo %d\n", 1);
	      return -1;
	    }
	  
	}
 
      chr = readUntil(infile, token, "):,;", depth);

      if( chr == ':' )
	{

	  dist = readDist(infile);

	  if( nodes[parentIndex].age == 0 ||
	      nodes[parentIndex].age < nodes[currentNode].age + dist )
	    {
	      nodes[parentIndex].age = nodes[currentNode].age + dist;
	      //fprintf(stderr, "parent: %d, age: %f, child: %d, age: %f\n", parentIndex, nodes[parentIndex].age, currentNode, nodes[currentNode].age);
	      
	    }
	  assert(nodes[parentIndex].age > nodes[currentNode].age);
	  
	  nodes[currentNode].branch = dist;
	  	  
	  if(!(chr = readUntil(infile, token, "):,", depth)))
	    {
	      fprintf(stderr, "Warning ... comus.phylo %d\n", 2);
	      return -1;
	    }
	}
      
      return (currentNode);
      
    }
  
  else
    {
      ++(*nodesInFile);
      
      node++;

      currentNode = node;

      nodes[currentNode].originalIndex = node;
      
      nodes[currentNode].father = parentIndex;

      nodes[currentNode].nson = 0;

      nodes[currentNode].age = 0.;
      
      /* nodes[currentNode].nodeStr = calloc(256, sizeof(char)); */
      
      /* strcpy(nodes[currentNode].nodeStr, "external"); */
     
      ++nodes[parentIndex].nson;
      
      if(!(chr = readUntil(infile, token, ":),", depth)))
	{
	  fprintf(stderr, "Warning ... comus.phylo %d\n", 3);
	  return -1;
	}
      
      token[0] = chr1;
      
      nodes[currentNode].nodeID = node; //atoi(token);

      if( chr == ':' )
	{
	  
	  dist = readDist(infile);
	  
	  //nodes[parentIndex].age = nodes[currentNode].age + dist;
	  
	  if( nodes[parentIndex].age == 0 ||
	      nodes[parentIndex].age < nodes[currentNode].age + dist )
	    {
	      nodes[parentIndex].age = nodes[currentNode].age + dist;
	    }
	  
	  if( nodes[parentIndex].age <= nodes[currentNode].age)
	    {
	      fprintf(stderr, "2. parent: %d, age: %f, child: %d, age: %f\n", parentIndex, nodes[parentIndex].age, currentNode, nodes[currentNode].age);
	      assert(nodes[parentIndex].age > nodes[currentNode].age);
	    }
	  
	  if(! (chr = readUntil(infile, token, ":),", depth)))
	    {
	      fprintf(stderr, "%d\n", 4);
	      return -1;
	    }
	  
	}
      
      return (currentNode);
      
    }
}

char skipSpaces(FILE *infile)
{

  char chr = fgetc(infile);

  while ( chr == 9 || chr == 13 || chr == 10 || chr == 32 )
    chr = fgetc(infile);

  return chr;
  
}


int readTreeFromFile(FILE *infile, int nspecies, int *positionsArray)
{

  char chr;
  
  if( (chr = skipSpaces(infile) ) != EOF)
    ungetc(chr, infile);

  if( chr == EOF)
    return 0;

  int nnodes = 2 * nspecies - 1;
  
  if(nodes == NULL)
    nodes = calloc(nnodes, sizeof(struct TREEN) );

  int *originalIndex = calloc( nnodes, sizeof(int) );

  assert(nodes != NULL);

  int nodesInFile = 0;
  
  int i, rootNode; //, it = 0, *nodea = (int*)space;

  int depth = 0, parentIndex = -1;

  for( i = 0; i < nnodes; ++i)
    {
      nodes[ i ].nson = 0;
    }

  nodes[0].branch = 0.;

  tree.root = rootNode = readNewickNode(infile, parentIndex, &depth, nspecies, &nodesInFile);

  
  if(rootNode != 0)
    {
      printf("Error in reading the input tree. Root Node is: %d, 0 is expected\n", rootNode);
      assert(rootNode == 0);
    }

  if(nodesInFile != nnodes)
    {
      fprintf(stderr, "The number of leaves in the tree ( %d ) does not match with the number of species in the command line ( %d )\n", (nodesInFile+1)/2, nspecies);
      assert(nodesInFile == nnodes);
    }

  qsort(nodes, nodesInFile, sizeof(struct TREEN), nodeTimeNameCmp);

  for ( i = 0; i < nnodes; ++i)
    {
      originalIndex[i] = nodes[i].originalIndex;
      positionsArray[ originalIndex[i] ] = i;
      /* printf("originalIndex %i is %d\n",  i, originalIndex[i]); */
    }

  assignChildrenToNodes(nodes, nodesInFile, positionsArray);

  /*set the tree.root */
  tree.root = -1;
  for( i = 0; i < nnodes; ++i)
    if( nodes[i].father == -1 )
      {
	tree.root = i;
	break;
      }
  assert(tree.root != -1 );

  tree.leaves = nspecies; 

  tree.internalNodes = nnodes - nspecies;

  tree.nbranch = tree.nnode + 1;

  /* set the branch */
  for( i = 0; i < nnodes; ++i )
    if( i != tree.root )
      nodes[i].branch = nodes[ nodes[i].father ].age - nodes[i].age;
    
  
  free(originalIndex);

  return 1;
}



/*******************************/


int OutSubTreeN (FILE *fout, int inode, int spnames, int printopt, char *labelfmt)
{
  int i,j, dad = nodes[inode].father, nsib = (inode==tree.root ? 0 : nodes[dad].nson);

  if(inode != tree.root && inode == nodes[dad].sons[0])
    fputc ('(', fout);

  for(i=0; i<nodes[inode].nson; i++)
    OutSubTreeN(fout, nodes[inode].sons[i], spnames, printopt, labelfmt);
  
  if(nodes[inode].nson==0) { 
    
    fprintf(fout, "%d", inode+1);
  }
  
  if((printopt & PrNodeNum) && nodes[inode].nson) 
    fprintf(fout," %d ", inode+1);
  if((printopt & PrLabel) && nodes[inode].label>0)
    fprintf(fout, labelfmt, nodes[inode].label);
  if((printopt & PrAge) && nodes[inode].age) 
    fprintf(fout, " @%.3f", nodes[inode].age);

  /*  Add branch labels to be read by Rod Page's TreeView. */
#if (defined CODEML)
  if((printopt & PrOmega) && inode != tree.root)
    fprintf(fout," '#%.4f' ", nodes[inode].omega);
#elif (defined (EVOLVER) || defined (MCMCTREE))
  if((printopt & PrLabel) && nodes[inode].nodeStr && nodes[inode].nodeStr[0])
    fprintf(fout," '%s'", nodes[inode].nodeStr);
#endif
  
  if((printopt & PrBranch) &&  
     (inode!=tree.root || 
      nodes[inode].branch>0))
    {
      fprintf(fout,":%.16f:%d", nodes[inode].branch, nodes[inode].nson);
      if(nodes[inode].nson > 0)
	{
	  fprintf(fout, ":[");
	  for(j = 0; j < nodes[inode].nson; ++j)
	    {
	      fprintf(fout, "%d ", nodes[inode].sons[j]+1);
	    }
	  fprintf(fout, "]");
	}
    }

  if(nsib == 0)            /* root */
    fprintf(fout, ";\n");
  else if (inode == nodes[dad].sons[nsib-1])  /* last sib */
    fputc(')', fout);
  else                     /* not last sib */
    fprintf(fout, ", ");

  return (0);
}


int OutTreeN (FILE *fout, int spnames, int printopt)
{
  /* print the current tree.
     Can the block of print statements be moved inside the recursive function?
  */
  int i, intlabel=1;
  char* labelfmt[2]={"'#%.5f'", "'#%.0f'"};

  if(printopt & PrLabel) {
    for(i=0; i<tree.nnode; i++) 
      if(nodes[i].label-(int)nodes[i].label != 0) intlabel=0;
  }

  OutSubTreeN(fout, tree.root, spnames, printopt, labelfmt[intlabel]);

  return(0);
}


/* code from Fletcher and Yang's INDELible */
void NodeToBranchSub (int inode)
{
  int i, ison;
  
  for(i=0; i<nodes[inode].nson; i++) {
    
    tree.branches[tree.nbranch][0] = inode;
    tree.branches[tree.nbranch][1] = ison = nodes[inode].sons[i];
    nodes[ison].ibranch = tree.nbranch++;
    
    if(nodes[ison].nson>0)  
      NodeToBranchSub(ison);
  }
  
}

/* code from Fletcher and Yang's INDELible */
void NodeToBranch (void)
{
  tree.nbranch=0;
  
  NodeToBranchSub (tree.root);

  /* printf("\n\ntree.nnode: %d, tree.nbranch+1: %d\n",  */
  /* 	 tree.nnode, tree.nbranch + 1); */
  

  if(tree.nnode != tree.nbranch+1)
    {
      fprintf(stderr, "\n\nERROR!!!    tree.nnode: %d, tree.nbranch+1: %d\n", 
	      tree.nnode, tree.nbranch + 1);
      
      assert( tree.nnode == tree.nbranch+1);
    }
  
}





/* code from Fletcher and Yang's INDELible */
int RandomLHistory (int ns, int rooted, double space[])
{
    
  /* 
     random coalescence tree, with each labeled history having equal probability. interior nodes are numbered ns, ns+1, ..., 2*ns-1-!rooted
  */

  int i, j, it=0, 
    *nodea=(int*)space;
  
  for (i=0; i<2*ns-1-!rooted; i++) 
    ClearNode(i);

  for (i=0; i<ns; i++) 
    nodea[i]=i;
  
  for (i=ns; i>(1+!rooted); i--) {
    
    nodes[ it=2*ns-i ].nson=2;

    nodes[it].nodeID = -1;
    
    j=(int)(i*rndu());
    
    nodes[ nodea[j] ].father=it; 

    nodes[it].sons[0] = nodea[j];

    nodes[ nodea[j] ].nodeID = nodea[j];
    
    nodea[j]=nodea[i-1];
    
    j=(int)((i-1)*rndu());
    
    nodes[nodea[j]].father=it; 

    nodes[it].sons[1]=nodea[j];

    nodes[nodea[j]].nodeID = nodea[j];
    
    nodea[j]=it;
    
    if (!rooted && i==3) {
      nodes[it].nson++;
      nodes[nodea[1-j]].father=it; 
      nodes[it].sons[2]=nodea[1-j];
    }
    
  }
  tree.root=it;  tree.nnode=ns*2-1-!rooted;
  
  NodeToBranch();
  
  return (0);
}



void BranchLengthBD_RannalaYang(int ns, int rooted, double birth, double death, double sample, double mut)
{

  
  /* Generate random branch lengths (nodes[].branch) using the birth and
     death process with species sampling, or the Yule (coalescent?) process
     if sample=0, when only parameter mut is used.
     Note: older interior nodes have larger node numbers, so root is at
     node com.ns*2-2 with time t[ns-2], while the youngest node is at 
     node com.ns with time t[0].  When unrooted=0, the root is removed with
     branch lengths adjusted.
     This works with the tree generated from RandomLHistory().
  */

  int i,j, it, imin, fixt0=1;
  
  double la=birth, rho=sample, tmin, r, t[NS-1], dif = birth - death;
  double phi, eml, y, temp;
  

  //printf("birth: %f, death: %f, mutation: %e, sample: %f\n", birth, death, mut, rho);
  
  if (sample==0)  /* coalescent model.  Check this!!!  */
    for (i=ns,y=0; i>1; i--)
      nodes[ns*2-i].age = y += -log(rndu())/(i*(i-1.)/2.)*mut/2;
  else  
    {         /* BD with sampling */
      if (fixt0) 
	t[ns-2]=1;
      
      if (fabs(dif)>1e-6) {
	
	eml=exp(-dif);  
	
	phi=( rho * la * (eml - 1) - dif * eml ) / (eml-1);
	
	for (i=0; i<ns-1-(fixt0); i++) 
	  {
	    r=rndu(); 

	    temp = phi-r*rho*la;
	    
	    t[i] = ( log( temp ) - log(temp + r * dif) );
	    
	    t[i] = -t[i]/ dif;
	    
	  }
      }
      else
	for (i=0; i<ns-1-(fixt0); i++)
	  { 
	    r=rndu();  
	    
	    t[i]=r/(1+la*rho*(1-r)); 
	  }
      
      /* bubble sort */
      for (i=0; i<ns-1-1; i++) 
	{
	  for (j=i+1,tmin=t[i],imin=i; j<ns-1; j++)
	    if (tmin>t[j]) 
	      { 
		tmin=t[j]; 
		imin=j; 
	      }
	  t[imin]=t[i];  
	  t[i]=tmin;
	}
      
      for (i=ns; i>1; i--) 
	nodes[ns*2-i].age = t[ns-i] * mut;
    }
  
  for(i = 0; i < ns; ++i) 
    nodes[i].age = 0;
  
  for (i=0; i<tree.nnode; i++)
    {
      if (i!=tree.root)
	nodes[i].branch=nodes[nodes[i].father].age-nodes[i].age;
    }
  
  if (!rooted) 
    {
      it=nodes[tree.root].sons[2];
      nodes[it].branch =
	2*nodes[2*ns-2].age-nodes[tree.root].age-nodes[it].age;
    }

}



void BranchLengthBD_age_taxa_mrca(int ns, int rooted, double birth, double death, double sample, double mut, double tmrca)
{
  /* generate node ages by using the theory developed by Tanja Stadler
     It is the same as Ziheng's Yang method. It just extends a bit Ziheng's code 
     for arbitrary t-origin
  */
  
  int i,j, it, imin;
  
  double la=birth, mu=death, rho=sample, tmin, r, lamb1, mu1, *t = malloc( (ns - 1) * sizeof( double ) ), y;
  
  
  if (sample==0)  /* coalescent model.  Check this!!!  */
    
    for (i=ns,y=0; i>1; i--)
      
      nodes[ns*2-i].age = y += -log(rndu())/(i*(i-1.)/2.)*mut/2;
  
  else  
    {         /* BD with sampling */
      t[ns-2]=tmrca;
      
      
      lamb1 = rho * la;
      
      mu1 = mu - la * ( 1. - rho );
      
      if( la > mu )
	{
	  
	  for(i = 0; i < ns-2; ++i)
	    {
	      r = rndu();
	      
	      t[i] = 1/(lamb1-mu1)*log((lamb1-mu1* exp((-lamb1+mu1)*tmrca) -mu1*(1-exp((-lamb1+mu1)*tmrca)) *r )/(lamb1-mu1* exp((-lamb1+mu1)*tmrca) - lamb1*(1-exp((-lamb1+mu1)*tmrca)) *r )   )  ;
	      
            } 
	}
      else 
	{
	  for(i = 0; i < ns - 2; ++i)
	    {
	      r = rndu();
	      
	      t[i] =  -((tmrca* r)/(-1 - la * rho* tmrca + la* rho* tmrca*r ));
	    }
	}
	  
                
      /* bubble sort */
      for (i=0; i<ns-1-1; i++) 
	{
	  for (j=i+1,tmin=t[i],imin=i; j<ns-1; j++)
	    if (tmin>t[j]) 
	      { 
		tmin=t[j]; 
		imin=j; 
	      }
	  t[imin]=t[i];  
	  t[i]=tmin;
	}
      
      for (i=ns; i>1; i--) 
	nodes[ns*2-i].age = t[ns-i] * mut;
    }
  
  for(i = 0; i < ns; ++i) 
    nodes[i].age = 0;
  
  for (i=0; i<tree.nnode; i++)
    {
      if (i!=tree.root)
	nodes[i].branch=nodes[nodes[i].father].age-nodes[i].age;
    }
  
  if (!rooted) 
    {
      it=nodes[tree.root].sons[2];
      nodes[it].branch =
	2*nodes[2*ns-2].age-nodes[tree.root].age-nodes[it].age;
    }

  free(t);

}


void BranchLengthBD_Stadler_age_taxa_origin(int ns, int rooted, double birth, double death, double sample, double mut, double torigin)
{
  /* generate node ages by using the theory developed by Tanja Stadler
     It is the same as Ziheng's Yang method. It just extends a bit Ziheng's code 
     for arbitrary t-origin
  */
  
  int i,j, it, imin;
  
  double la=birth, mu=death, rho=sample, tmin, r, lamb1, mu1,  y;

  /* the size of this array is ns to include the origin */
  double *t = malloc( (ns ) * sizeof( double ) ); 
  
  if (sample==0)  /* coalescent model.  Check this!!!  */
    for (i=ns,y=0; i>1; i--)
      nodes[ns*2-i].age = y += -log(rndu())/(i*(i-1.)/2.)*mut/2;
  else  
    {         /* BD with sampling */
      t[ns-1]=torigin;

      lamb1 = rho * la;
      
      mu1 = mu - la * ( 1. - rho );
      
      if( la > mu )
	{
	  
	  /* compared to tmrca we sample one more value */
	  for(i = 0; i < ns-1; ++i)
	    {
	      
	      r = rndu();
	      
	      t[i] = 1/(lamb1-mu1)*log((lamb1-mu1* exp((-lamb1+mu1)*torigin) -mu1*(1-exp((-lamb1+mu1)*torigin)) *r )/(lamb1-mu1* exp((-lamb1+mu1)*torigin) - lamb1*(1-exp((-lamb1+mu1)*torigin)) *r ) );
	      
	      //fprintf(stderr, "t[%d]: %f, lamb1: %f, mu1: %f, dif: %f, torigin: %f\n", i, t[i], lamb1, mu1, lamb1-mu1, torigin);

            } 
	}
      else 
	{
	  /* compared to tmrca we sample one more value */
	  for(i = 0; i < ns - 1; ++i)
	    {
	      r = rndu();
	      
	      t[i] =  -((torigin * r)/(-1 - la * rho* torigin + la* rho* torigin * r ));
	    }
	}
	  
                
      /* bubble sort for all but the origin*/
      for (i=0; i<ns-1; i++) 
	{
	  for (j=i+1,tmin=t[i],imin=i; j<ns-1; j++)
	    
	    if (tmin>t[j]) 
	      { 
		tmin=t[j]; 
		imin=j; 
	      }
	  t[imin]=t[i];  
	  t[i]=tmin;
	}
      
      /* rescale based on the mutation rate */
      for (i=ns; i>1; i--) 
	nodes[ns*2-i].age = t[ns-i] * mut;
    }
  
  for(i = 0; i < ns; ++i) 
    nodes[i].age = 0;
  
  for (i=0; i<tree.nnode; i++)
    {
      if (i!=tree.root)
	nodes[i].branch=nodes[nodes[i].father].age-nodes[i].age;
    }
  
  if (!rooted) 
    {
      it=nodes[tree.root].sons[2];
      nodes[it].branch =
	2*nodes[2*ns-2].age-nodes[tree.root].age-nodes[it].age;
    }

  free(t);

}


/* 
   generate node ages by using the theory developed by Tanja Stadler 
   In detail, this function corresponds to 
   sim2.bd.fast.single from TreeSim
   or equivalently to sim.bd.taxa.loop when complete==FALSE and stochsampling ==TRUE

   As Tanja describes in the manual these settings correspond to 
   "If stochsampling=TRUE: Each tip is included into the final tree with
   probability frac."
*/
   
void BranchLengthBD_Stadler_ntaxa(int ns, int rooted, double birth, double death, double sample, double mut)
{
    
  int i,j, imin;
  
  double la=birth, mu=death, rho=sample, tmin, r, lamb1, mu1, y, enadiani, rtoenadiani, torigin;
  
  /* the size of this array is ns to include the origin */
  double *t = malloc( (ns ) * sizeof( double ) ); 
  
  if (sample==0)  /* coalescent model.  Check this!!!  */
    {
      for (i=ns,y=0; i>1; i--)
	
	nodes[ns*2-i].age = y += -log( rndu() )/(i*(i-1.)/2.)*mut/2;
      
    }
  else  
    {
      
      lamb1 = rho * la;
      
      mu1 = mu - la * ( 1. - rho );
      
      enadiani = 1./ns;
      
      /* sample torigin */
      r = rndu();

      rtoenadiani = pow(r, enadiani);
      
      if( la > mu )
	{
	  
	  torigin = log((-lamb1 - la * r * enadiani + mu * r * enadiani + lamb1 * r * enadiani) / (lamb1 * ( -1 + rtoenadiani) ) ) / (la - mu);

	  
	  
	}
      else
	torigin = -rtoenadiani / (la * ( -1 + rtoenadiani * rho ) );
      
      t[ns-1]=torigin;
      
      if( la > mu )
	{
	  
	  /* compared to tmrca we sample one more value */
	  for(i = 0; i < ns-1; ++i)
	    {
	      r = rndu();
	      
	      t[i] = 1/(lamb1-mu1)*log((lamb1-mu1* exp((-lamb1+mu1)*torigin) -mu1*(1-exp((-lamb1+mu1)*torigin)) *r )/(lamb1-mu1* exp((-lamb1+mu1)*torigin) - lamb1*(1-exp((-lamb1+mu1)*torigin)) *r )   )  ;

	      //fprintf(stderr, "t[%d]: %f, lamb1: %f, mu1: %f, dif: %f\n", i, t[i], lamb1, mu1, lamb1-mu1);
	      
	    } 
	}
      else 
	{
	  /* compared to tmrca we sample one more value */
	  for(i = 0; i < ns - 1; ++i)
	    {
	      r = rndu();
	      
	      t[i] =  -((torigin * r)/(-1 - la * rho* torigin + la* rho* torigin * r ));
	    }
	}
      
      /* bubble sort for all but the origin*/
      for (i=0; i<ns-1; i++) 
	{
	  for (j=i+1,tmin=t[i],imin=i; j<ns-1; j++)
	    
	    if (tmin>t[j]) 
	      { 
		tmin=t[j]; 
		imin=j; 
	      }
	  t[imin]=t[i];  
	  t[i]=tmin;
	}
      
      /* rescale based on the mutation rate */
      for (i=ns; i>1; i--) 
	nodes[ns*2-i].age = t[ns-i] * mut;
      
    }
  
  for(i = 0; i < ns; ++i) 
    nodes[i].age = 0;
  
  for (i=0; i<tree.nnode; i++)
    {
      if (i!=tree.root)
	nodes[i].branch=nodes[nodes[i].father].age-nodes[i].age;
    }
  
  free(t);

}
   
void BranchLengthBD_Stadler_ntaxa_ProcessOldestStart(int ns, int rooted, double birth, double death, double sample, double mut, double oldestOrigin)
{
  
    
  int i,j, imin, cnt =0; 
  
  double la=birth, mu=death, rho=sample, tmin, r, lamb1, mu1, enadiani, rtoenadiani, torigin;
  
  /* the size of this array is ns to include the origin */
  double *t = malloc( (ns ) * sizeof( double ) ); 
  
  assert( sample > 0 );
  
  lamb1 = rho * la;
  
  mu1 = mu - la * ( 1. - rho );
  
  enadiani = 1./ns;
  
  /* sample torigin */

  torigin = DBL_MAX;

  cnt = 0;
  
  if( la > mu )
    {

      while( torigin > oldestOrigin && cnt < cntMAX)
	{
	  cnt++;
	  
	  r = rndu();
	  
	  rtoenadiani = pow(r, enadiani);
	  
	  torigin = log((-lamb1 - la * r * enadiani + mu * r * enadiani + lamb1 * r * enadiani) / (lamb1 * ( -1 + rtoenadiani) ) ) / (la - mu);

	}
      
    }
  else
    {
      while( torigin > oldestOrigin && cnt < cntMAX)
      {
	cnt++;
		
	r = rndu();
	  
	rtoenadiani = pow(r, enadiani);
	
	torigin = -rtoenadiani / (la * ( -1 + rtoenadiani * rho ) );
	
      }
    }

  if( cnt >= cntMAX )
    {
      printf("\nCannot get a tree with taxa: %d conditioning that origin is younger than %f after %d trials\n\n", ns, oldestOrigin, cntMAX);

      assert( cnt < cntMAX );
    }
  
  t[ns-1]=torigin;
  
  if( la > mu )
    {
      
      /* compared to tmrca we sample one more value */
      for(i = 0; i < ns-1; ++i)
	{
	  r = rndu();
	  
	  t[i] = 1/(lamb1-mu1)*log((lamb1-mu1* exp((-lamb1+mu1)*torigin) -mu1*(1-exp((-lamb1+mu1)*torigin)) *r )/(lamb1-mu1* exp((-lamb1+mu1)*torigin) - lamb1*(1-exp((-lamb1+mu1)*torigin)) *r )   )  ;
	  
	} 
    }
  else 
    {
      /* compared to tmrca we sample one more value */
      for(i = 0; i < ns - 1; ++i)
	{
	  r = rndu();
	  
	  t[i] =  -((torigin * r)/(-1 - la * rho* torigin + la* rho* torigin * r ));
	}
    }
  
  /* bubble sort for all but the origin*/
  for (i=0; i<ns-1; i++) 
    {
      for (j=i+1,tmin=t[i],imin=i; j<ns-1; j++)
	
	if (tmin>t[j]) 
	  { 
	    tmin=t[j]; 
	    imin=j; 
	  }
      t[imin]=t[i];  
      t[i]=tmin;
    }
  
  /* rescale based on the mutation rate */
  for (i=ns; i>1; i--) 
    nodes[ns*2-i].age = t[ns-i] * mut;
  

  
  for(i = 0; i < ns; ++i) 
    nodes[i].age = 0;
  
  for (i=0; i<tree.nnode; i++)
    {
      if (i!=tree.root)
	nodes[i].branch=nodes[nodes[i].father].age-nodes[i].age;
    }
  
  free(t);

}



void BranchLengthBDProcess( int  ns, int rooted, double birth, double death, double sample, double mut, double torigin, int mode, double oldestOrigin)

{
  
  if(mode == 0)
    {
      if(torigin != 1.0)
	{
	  fprintf(stderr, "Warning: you have specified an origin for the process and this is ignored\n");
	}
      BranchLengthBD_RannalaYang( ns, rooted, birth, death, sample, mut);
    }
  else if( mode == 1)
    {
      BranchLengthBD_age_taxa_mrca( ns, rooted, birth, death, sample, mut, torigin);
    }
  else if( mode == 2)
    {
      
      assert(torigin > 0);

      BranchLengthBD_Stadler_age_taxa_origin( ns, rooted, birth, death, sample, mut, torigin);
    }
  
  else if( mode == 3)
    {
      
      BranchLengthBD_Stadler_ntaxa(ns, rooted, birth, death, sample, mut);
      
    }

  else if (mode == 4 )
    {
      BranchLengthBD_Stadler_ntaxa_ProcessOldestStart(ns, rooted, birth, death, sample, mut, oldestOrigin);
    }
  
}

int speciationNodesYoungerThanSampling( double *speciesSamplingTime, int n)
{
  
  int i, father;

  for( i = 0; i < n; ++i)
    {
      father = nodes[i].father;
      
      if( nodes[father].age < speciesSamplingTime[ i ] )
	return 1;

    }

  return 0;
}

int checkTreeSamplingConsistency(double *samplingTime, int n)
{
  
  int i, father;

  for( i = 0; i < n; ++i)
    {
      father = nodes[i].father;
      
      if( nodes[father].age < samplingTime[ i ] )
	{
	  fprintf(stderr, "Error... Sampling time (%f) for species %d is older than its speciation event time (%f)\n", samplingTime[i], i, nodes[father].age);
	  return 0;
	}
    }
  
  return 1;
}


void randomTreeGeneration(int ntaxa, double birth, double death, double sample, double mut, int option, char *outputTree, int mode, double torigin, FILE *phyloTreeOut, double oldestOrigin, double *samplingTime )
{
  
  // Function written by W.F. adapting code from the main body of evolver.c in PAML
  // produces random trees using the birth-death process
  // Yang and Rannala. (1997) Bayesian phylogenetic inference using DNA sequences: a Markov Chain Monte Carlo Method. Mol. Biol. Evol. 14:717-724.  
  
  double *space;
  
  int ntree=1, i, ns = com.ns = ntaxa;
  
  if(ns>NS) fprintf(stderr, "Too many species.  Raise NS.");
  
  if((space=(double*)malloc(10000*sizeof(double)))==NULL) 
    fprintf(stderr, "oom");

  int rooted=!(option%2);

  int BD=1; //use birth death process for branch lengths
  
  if( nodes == NULL )
    fprintf(stderr, "oom");
  
  for(i=0; i<ns; i++)          /* default spname */
    sprintf(com.spname[i],"S%d",i+1);

  for(i=0;i<ntree;i++) {
    
    /* TODO not very sure whehter this should be in or out the do-while
       loop
    */
    //RandomLHistory (com.ns, rooted, space);

    int trials = 0; 
  
    do
      {
	
	trials++;

	if(trials > MAXTRIALS && trials % 10000000 == 0)
	  {
	    fprintf(stderr, "\nWarning... with the current phylogenetics parameters (birth rate, sampling times)\nyour simulations will take some time (if they ever finish  --- (current trials %d)\nEither decrease the speciation rate, or decreease the sampling time\n", trials);
	    
	  }

	//if(trials % 1 == 0)
	RandomLHistory (com.ns, rooted, space);
	
	BranchLengthBDProcess(com.ns, 1, birth, death, sample, mut, torigin, mode, oldestOrigin);

      }while( speciationNodesYoungerThanSampling( samplingTime, ntaxa ) == 1 );

    fprintf(stderr, "Tree was successful after %d trials (due to ancestral sampling)\n", trials);
      if(phyloTreeOut != NULL && com.ns<20 && ntree<10) 
      {
	OutTreeN(phyloTreeOut,0,BD);
      }
    else if(phyloTreeOut != NULL)
      {
	OutTreeN(phyloTreeOut,1,BD);
      }
  }

  free(space);  
}



