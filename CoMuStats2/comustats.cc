#include <iostream>
#include <vector>
#include <Sequence/PolySNP.hpp>
#include <sstream>
#include <Sequence/Fasta.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/FST.hpp>
#include <Sequence/SimData.hpp>
#include <numeric>
#include <functional>
#include <stdexcept>
#include <ctime>
#include <limits>
#include <cmath>
#include <sys/time.h>
#include <Sequence/Kimura80.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include <cstring>

#include <cassert>

#include <string>



#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/SeqUtilities.hpp>
#include <Sequence/CountingOperators.hpp>
#include <algorithm>
#include <vector>
#include <iostream>
#include <functional>
#include <cctype>

//so that you know where things come from
using std::map;
using std::vector;
using std::cout;
using std::cerr;
using Sequence::Fasta;
using Sequence::PolySites;
using Sequence::makeCountList; // <Sequence/SeqUtilities.hpp>
using Sequence::operator+=;    //
using Sequence::Alignment::GetData;
using Sequence::Alignment::IsAlignment;

using Sequence::Alignment::GetData;
using Sequence::Alignment::IsAlignment;


using namespace std;
using namespace Sequence;



bool validStates(const std::map<char,unsigned> & counts)
/*
  counts is assumed to be a map of nucleotides and their number of occurences.
  The function returns false if any of the nucleotides are not in the set {A,G,C,T,N},
  and toupper is used so that matching is case-insensitive.
 */
{
  for( map<char,unsigned>::const_iterator i = counts.begin() ;
       i != counts.end() ;
       ++i )
    {
      char ch = toupper(i->first);
      if ( ch != 'A' && ch != 'G' && ch != 'C' && ch != 'T' && ch != 'N' )
	return false;
    }
  return true;
}


/*
  gets the current time. It is used to
  measure the time of different parts of the code
*/
double gettime(void)
{
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
}


/*
  split a string to a vector of strings 
*/
vector<string> split(const string &s, char delim)
{
  std::stringstream ss(s);
  
  string item; 
  
  vector<string> tokens;

  while( getline(ss, item, delim) ) {

    /* push back only non-empty items */
    if(item.length() > 0 )
      tokens.push_back(item);
  }
  return tokens;
    
}


/* 
   split a string to a vector of doubles
*/
vector<double> splitDouble(const string &s, char delim)
{
  stringstream ss(s);

  string item; 
  
  vector<double> tokens;
  
  while( getline(ss, item, delim) ) {
    
    /* push back only non-empty items */
    if(item.length() > 0 )
      tokens.push_back(atof(item.c_str()));
  }
  
  return tokens;

}



/* 
   get the data for one dataset from an ms file
*/
int getMSData( ifstream &inputFile, vector <double > &pos, vector <string> &data, bool haveOutgroup, unsigned outgroupSeq ){
  
  /*!
    The function starts reading every line in a ms file.
    It searches for the *s* of segsites. This triggers reading
    a new dataset. Then, positions and the actual data is read
  **/
  
  string str, seq, name;

  static int cnt = 0;
  
  unsigned int newDataset = 0, j = 0, seqCounter = 0;
  // numberSegsites = 0, 
    
  
  vector<string> tempStringVector; 

  /* clear both the pos and the data vector. 
     This will guarantee that no previous data will remain there
  */
  pos.clear();
  data.clear();
  
  
  /* read until the end of the file 
     There is, however, a break statement below if a '//' is found
   */
  while( !inputFile.eof() )
    {
      
      getline( inputFile, str);
      
      j = 0; 
      
      /**
	 skip the white spaces
      **/
      while(str[j] == 32 || str[j] == 9 || str[j] == 13 || str[j] == 10)
	++j;

      /* this means that the line contains only white spaces */
      if(j == str.length())
	continue;

      /* if we find a // string then we have found a new dataset */
      if( str[j] == '/' && ( str.length() > j+1 && str[j+1] == '/' ) )
	{
	  
	  cnt++;

	  newDataset = 1;

	  /* break when the next dataset starts
	     for the first dataset this does not apply
	  */
	  if( cnt > 1 )
	    break;
	  
	}
      
      
      if(str[j] == 's')
	{
	  
	  tempStringVector = split(str, ' ');
	  
	  //numberSegsites = atoi( tempStringVector[1].c_str() );
	  	  
	  newDataset = 1;

	  continue;
	  
	}

      /* this is meant to get the *p*ositions line */
      if( str[j] == 'p' && newDataset == 1)
	{
	  tempStringVector = split(str, ' ');

	  for( unsigned int i = 0; i < tempStringVector.size()-1; ++i)
	    {
	      pos.push_back( atof(tempStringVector[i+1].c_str()) );
	      
	    }

	  newDataset = 2;
	  continue;
	}

      
      /* 
	 the lines that start with 1 or 0 contain the actual data
      */
      if( (str[j] == '1' || str[j] == '0') && (newDataset ==2 || newDataset == 3) )
	{
	  if(haveOutgroup == false || seqCounter != outgroupSeq)
	    data.push_back(str);
	  
	  newDataset = 3;
	}
    
    }
  
  // for( int i = 0; i < pos.size(); ++i )
  //   cout << pos[i] << "\t";
  // cout << endl;

  // for( int i = 0; i < data.size(); ++i)
  //   cout << data[i] << endl;

  if( pos.size() > 0 && data.size() > 0)
    return 1;

  return 0;
  
}



int getAlignment( ifstream &inputFile, vector<Fasta> &fst)
{
  string str, seq, name;

  unsigned int j = 0, length = 0, i = 0, cnt = 0;

  while( ! inputFile.eof() )
    {
      
      getline( inputFile, str);

      // cerr << "*** " << str << endl;

      j = 0;

      while(str[j] == 32 || str[j] == 9 || str[j] == 13 || str[j] == 10)
	++j;

      if(j == str.length())
	continue;
      
      if( str[j] == '>' && cnt == 0 ) 
	{
	  cnt++;
	  name = str;
	  i = 1;
	  seq = "";
	}

      else if( str[j] == '>' && cnt > 0 )
	{

	  //cerr << seq.length() << "\t" << length <<endl;

	  if( length != 0 && seq.length() != length )
	    {
	      cerr << "Sequence " <<  cnt << " has length " << seq.length() << ". Length should be " << length << endl;

	      
	      assert( length == seq.length() );
	    }

	  
	  length = seq.length();

	  
	  Sequence::Fasta newFasta(name, seq);
	  
	  fst.push_back(newFasta);

	  cnt++;
	  
	  name = str;

	  seq = "";

	  i = 1;
	}
      
      else if(str[j] != '/' && (i == 1 || i == 2 ) )
	{
	  seq += str;
	  
	  i = 2;
	}

      else if(str[j] == '/' && name.length()  < 1 )
	continue;
      
      else if(str[j] == '/' && str.length() > j+1 && str[j+1] == '/')
	break;
      
      else if(str[j] == '/')
	continue;
      
      else
	{
	  cerr << "Something is wrong... unexpected line: " << str << endl;
	  exit(-1);
	}
      
    }


  if( i == 2 )
    {

      
      if( length != 0 && seq.length() != length )
	{
	  cerr << "Sequence " <<  cnt << " has length " << seq.length() << ". Length should be " << length << endl;
	  
	  
	  assert( length == seq.length() );
	  
	  length = seq.length();
	}
      
      
      Sequence::Fasta newFasta(name, seq);
      
      fst.push_back(newFasta);
      
    }
  
  assert(cnt == fst.size());

  if( cnt > 0 )
    {
      assert( length > 0 );
      cerr << "Alignment contains " << cnt << " sequences, with length " << length << "." << endl;
    }

  
  if(inputFile.eof() == 0 || fst.size() > 0)
    return 1;
  
  return 0;
  
}

void printVersions()
{
  fprintf(stderr, "\nCoMuStats version 3.1, March 2016\n");
  fprintf(stderr, "\n\n------------------ COMMENTS --------------------\n");
  fprintf(stderr, "Please cite the manuscript of CoMuS, as well as K. Thornton's Libsequence\n");
  fprintf(stderr, "The SFS code is based on ufs.cc of K.Thornton (see the examples of libsequence\n");	  
  fprintf(stderr, "\n\n===================================\n\n\n");
}


unsigned int getDif( const std::string &seq1, const std::string &seq2)
{
  
  unsigned int dif;
  
  dif = inner_product( seq1.begin(), seq1.end(), seq2.begin(), 0, std::plus<unsigned int>(), std::not2(std::equal_to<std::string::value_type>()) );
  
  return dif;
}

unsigned int getMaxDistance( vector <Fasta> &data, bool haveOutgroup, int outgroup )
{
  unsigned int maxDis = 0, dis;
  string seq1, seq2;

  if(haveOutgroup == false)
    {

      for(unsigned int i = 0; i < data.size(); ++i)
	{
	  seq1 = data[i].GetSeq();
	  for(unsigned int j = i + 1; j < data.size(); ++j)
	    {
	      seq2 = data[j].GetSeq();
	      
	      dis = getDif( seq1, seq2 );
	      
	      if(maxDis <  dis )
		maxDis = dis;
	      
	    }
	}
    }

  else 
    {
      assert( outgroup > -1 && (unsigned)outgroup < data.size() );

      for(unsigned int i = 0; i < data.size(); ++i)
	{
	  
	  if(i == (unsigned)outgroup)
	    continue;
	  
	  seq1 = data[i].GetSeq();
	  for(unsigned int j = i + 1; j < data.size(); ++j)
	    {
	      
	      if(j == (unsigned)outgroup)
		continue;
	  
	      seq2 = data[j].GetSeq();
	      
	      dis = getDif( seq1, seq2 );
	      
	      if(maxDis <  dis )
		maxDis = dis;
	    }
	}
            
    }
  return maxDis;
}



int getSubSet(unsigned int k, unsigned int n, bool haveOutgroup, int outgroup, const unsigned int *config, const vector<Fasta> &data, vector<Fasta> &out)
{

  
  assert(k >= 0 && k < n);
  
  unsigned int start = 0, i = 0, j, jj;

  for(i = 0; i < k; ++i)
    start += config[i];

  if( haveOutgroup && start >= (unsigned)outgroup )
    start++;

  jj = j = start;
  
  while( j  < start + config[k])
    {
      
      if( haveOutgroup == false || (unsigned)outgroup != jj)
	{
	  out.push_back( data[jj] );
	  ++j;
	  ++jj;
	}
      else if( haveOutgroup == true && (unsigned)outgroup == jj)
	{
	  ++jj;
	  continue;
	}
    }

  /* if there is outgroup, this is in the end of each subset */
  if( haveOutgroup == true )
    out.push_back( data[ outgroup ] );

  return 1;

}



int getSubSetMS(unsigned int k, unsigned int n, const unsigned int *config, const vector<string> &data, vector<string> &out, const bool haveOutgroup, const int outgroupSeq)
{
  assert(k >= 0 && k < n);

  out.clear();

  if(haveOutgroup == true )
    {
      assert( (unsigned)outgroupSeq >= 0 );
      assert( (unsigned)outgroupSeq < data.size() );
    }

  
  unsigned int start = 0, i = 0, j;

  for(i = 0; i < k; ++i)
    start += config[i];

  /* if there is an outgroup sequence before the "supposed start" then
     we should actually move on sequence below
  HERE HERE HERE
  */
  if( haveOutgroup == true && (unsigned)outgroupSeq < start)
    ++start;

  /* j is the index of the first sequence of the subset */
  j = start;
  
  while( j  < start + config[k])
    {
      out.push_back( data[j] );
      ++j;
    }

  if( haveOutgroup == true)
    out.push_back( data[outgroupSeq] );
    
  return 1;

}


int getSubSetGroups(unsigned *k, // which groups we'd like to get
		    unsigned int kn, // how many we 'd like to get
		    unsigned int n, // total sample size
		    bool haveOutgroup, // is there any outgroup
		    int outgroup, // the index of the outgroup
		    const unsigned int *config, // configuration of samples within the populations
		    const vector <Fasta> &data, // all fasta data
		    vector <Fasta> &out, // output fasta vector
		    unsigned int **outconf // configuration of the output vector
		    )
{

  /* empty the vector */
  
  if( out.empty() == false)
    out.erase( out.begin(), out.end());

  unsigned int i = 0, j = 0;

  *outconf = (unsigned int*)calloc( kn, sizeof(unsigned int) );

  int *start = (int*)calloc( n, sizeof(int) );

  vector <Fasta> tmp;

  start[0] = 0;
  for( i = 1; i < n; ++i)
    start[i] = start[ i  - 1 ] + config[i - 1];
    

  assert( *outconf != NULL);

  for( i = 0; i < kn; ++i )
    {
      assert( k[i] >= 0 && k[i] < n );

      (*outconf)[i] = config[ k[i] ];
      
    }

  for( i = 0; i < kn; ++i)
    {
      if( getSubSet( k[i], n, haveOutgroup, outgroup, config, data, tmp ) == 0 )
	{
	  fprintf(stderr, "\nERROR Cannot get subset %d\n\n", k[i]);

	  assert(0);
	}

      for(j = 0; j < tmp.size(); ++j)
	{
	  out.push_back( tmp[j] );
	}
      
      tmp.erase( tmp.begin(), tmp.end() );
	
    }

  return 1;

}




int getSubSetGroupsMS(unsigned int *k, // which groups we'd like to get
		      unsigned int kn, // how many we 'd like to get
		      unsigned int n, // total sample size
		      const unsigned int *config, // configuration of samples within the populations
		      const vector <string> &data, // all fasta data
		      vector <string> &out, // output fasta vector
		      unsigned int **outconf // configuration of the output vector
		    )
{

  /* empty the vector */
  
  if( out.empty() == false)
    out.erase( out.begin(), out.end());

  unsigned int i = 0, j = 0;

  *outconf = (unsigned int*)calloc( kn, sizeof(unsigned int) );

  int *start = (int*)calloc( n, sizeof(int) );

  vector <string> tmp;

  start[0] = 0;
  for( i = 1; i < n; ++i)
    start[i] = start[ i  - 1 ] + config[i - 1];
    

  assert( *outconf != NULL);

  for( i = 0; i < kn; ++i )
    {
      assert( k[i] >= 0 && k[i] < n );
      
      (*outconf)[i] = config[ k[i] ];
      
    }

  for( i = 0; i < kn; ++i)
    {
      if( getSubSetMS( k[i], n, config, data, tmp, false, 0 ) == 0 )
	{
	  fprintf(stderr, "\nERROR Cannot get subset %d\n\n", k[i]);

	  assert(0);
	}

      for(j = 0; j < tmp.size(); ++j)
	{
	  out.push_back( tmp[j] );
	}
      
      tmp.erase( tmp.begin(), tmp.end() );
	
    }

  return 1;

}


/* this part of the code uses many lines from the ufs.cc of K. Thornton */
vector< unsigned > ufs ( vector<Fasta> &f, unsigned outgroup){

  PolySites SNPtable(f, true, 1, false, false, 0);
  
  vector < unsigned > u (f.size(), 0) ;

  int polsites = 0;
  
  for( PolySites::const_site_iterator i = SNPtable.sbegin(); i != SNPtable.send(); ++i)
    {
      polsites++;
      const char ancstate = i->second[outgroup];
      
      std::map<char, unsigned> counts = makeCountList(i->second.begin(), i->second.begin() + outgroup);
      
      counts += makeCountList(i->second.begin()+outgroup+1,i->second.end());
      
      if( ! validStates(counts) )
	{
	  cerr << "site " << i->first << " contains characters other than A,G,C,T,N\n";
	}
      else
	{
	  int counter = 0;
	  
	  for( map<char,unsigned>::const_iterator i = counts.begin() ;
	       i != counts.end() ;
	       ++i )
	    {
	      counter++;
	      if ( toupper(i->first) != ancstate )
		{
		  u[i->second]++;
		}

	    }
	  
	}
    }

  u[0] = f[0].length() - polsites;
  
  return(u);
  
}


int main(int argc, char** argv)

{

  double time0 = gettime();

  char inputFileName[1024];
  
  strcpy(inputFileName, "XXXX");

  int header = 1, npop = 1;
  
  unsigned int outgroupSeq = 0;
  
  unsigned int *config = NULL, *sfs = NULL;
  
  vector<Fasta> data, subset, pairsData;

  vector<string> msdata, subsetms, pairsmsdata;
  vector<double> mspositions;

  bool haveOutgroup = false, sepPops = false, pairwiseFst = false, divergence = false, 
    lastIsOutgroup = false, windowsAnalysis = false,  slidingWindowMode = false, msflag = false, maxDisRel = false, sfsflag = false;

  int slidingWindowOffset = 0, slidingWindowLength = 0, totalSampleSize = 0;
  

  if(argc < 2)
    {
      printVersions();
      
      fprintf(stderr, "./CoMuStats -input <FileName> -npop <int> <int>...-outgroup <int>\n\n");
      fprintf(stderr, "It prints out summary statistics for fasta files\n");
      fprintf(stderr, "-v: prints out the version number\n\n");
      fprintf(stderr, "-input <filename>:\nreads in the file name. It can have multiple fasta alignments, separated by //. Note that the whole alignment should be separated by //.\n\n");
      fprintf(stderr, "-npop <int1> <int2> ...:\nint1 is the number of species, int2 ... the sample size per species, in the right5B5B5B order\n\n");
      fprintf(stderr, "-outgroup <int>: the outgroup sequence (1 to the number of sequences in the alignment)\n\n");
      fprintf(stderr, "-div : it calculates the divergence from the outgroup and normalizes several statistics against divergence\n\n");
      fprintf(stderr, "-maxDisRel : calculates the relative maximum distance\n\n");
      fprintf(stderr, "-sepPops : calculates summary statistics also in separate groups (defined by npops)\n");
      fprintf(stderr, "-pairwiseFst : calculates the Fst for all the population pairs\n");
      fprintf(stderr, "-slidwin <int: window length> <int: window offset> : calculate statistics in a sliding window fashion. You should provide the length and the offset for the window\n");
      fprintf(stderr, "-ms : data are in ms format (default is fasta)\n");
      fprintf(stderr, "-sfs : get the sfs (currently it works only with ms format -- January 2016\n");
      fprintf(stderr, "IMPORTANT: When outgroup is specified, and it is the ith sequence, the program will read all sequences and will set the ith as the outgroup. For FST calcuations, the npop must be the populations of the ingroup, and the sum of sample sizes should be the total sample size of the INGROUP, i.e. total seqs - 1\n\n");
      exit(-1);
    }

  for(int i = 1; i < argc; ++i)
    {

      /* turns on the sfs flag */
      if(strcmp(argv[i], "-sfs") == 0 )
	{
	  sfsflag = true;
	  continue;
	}
	  
      
      /* switches on the calculations for an ms-like file */
      if(strcmp(argv[i], "-ms") == 0 )
	{
	  msflag = true;
	  continue;
	}


      /* switches on the sliding window mode
	 additional flags are required: the window length and the window offset
      */
      if(strcmp(argv[i], "-slidwin") == 0 )
	{
	  windowsAnalysis = true;
	  slidingWindowMode = true;
	  slidingWindowLength = atoi(argv[++i]);
	  slidingWindowOffset = atoi(argv[++i]);
	  
	  continue;
	}

      /*
	prints the version of the program
      */
      if(strcmp(argv[i], "-v") == 0)
	{
	  printVersions();
	  exit(1);
	}


      /* 
	 the input file
      */
      if( strcmp(argv[ i ] , "-input" ) == 0 )
	{
	  strcpy(inputFileName, argv[++i] );
	  
	  continue;
	}

      
      /*
	don't print the header. This is useful when later I'd like to cat the various files
      */
      if( strcmp( argv[ i ], "-noheader") == 0 )
	{
	  header = 0;
	  
	  continue;
	}

      /*
	calculate the maximum relative distance. TODO: what exactly is this?
      */
      if( strcmp( argv[ i], "-maxdis") == 0)
	{
	  maxDisRel = true;
	  continue;
	}

      /* From how many populations do the data come? */
      if( strcmp( argv[ i ], "-npop") == 0)
	{
	  npop = atoi( argv[ ++i ]);

	  assert( npop > 0 );

	  config = (unsigned int* )calloc( npop, sizeof(unsigned int) );
	  
	  	  
	  if( i + npop >= argc )
	    {
	      
	      fprintf(stderr, "npop: %d, i: %d, argc: %d\n", npop, i, argc);
	      assert( i + npop  < argc );
	      
	    }	  
	  

	  for( int j = 0; j < npop; ++j)
	    {
	      
	      

	      assert( argv[ i+1 ][0] != '-' );
	      config[ j ] = atoi(argv[ ++i ]);
	    }

	  continue;
	}

      /*
	which sequence is the outgroup?
	Here, you must provide the index of the sequence, start counting from 1
      */
      if( strcmp( argv[ i ], "-outgroup") == 0 )
	{
	  haveOutgroup = true;

	  /* often the last sequence is the outgroup
	     then, just write "last"
	  */
	  if( strcmp( argv[i+1], "last") == 0 )
	    {
	      ++i;
	      lastIsOutgroup = true;
	    }
	  else
	    {
	      outgroupSeq = atoi( argv[ ++i]) - 1;
	      assert(outgroupSeq >= 0);
	    }
	  
	  continue;
	}

      /* should I calculate divergence? */
      if( strcmp( argv[ i ], "-div") == 0 )
	{
	  divergence = true; 
	  
	  if( haveOutgroup == false)
	    {
	      cerr << "Please provide an outgroup first using the -outgroup" << endl;
	      assert(haveOutgroup == true);
	    }
	  continue;
	  
	}

      /* calculate statistics separately in populations */
      if( strcmp( argv[ i ], "-sepPops") == 0)
	{
	  sepPops = true;
	  continue;
	}

      /* calcuate pairwise FST */
      if( strcmp( argv[ i ], "-pairwiseFst") == 0 )
	{
	  pairwiseFst = true;
	  continue;
	}
     
      else
	{
	  fprintf(stderr, "Argument %s is not supported\n", argv[i]);
	  exit(0);
	}
    }


  /* open the file. It may contain several fasta alignments separated by  // 
   */

  // ifstream inputFile("test.txt");
  
  ifstream inputFile(inputFileName, ifstream::in);

  if(!inputFile.is_open())
    {
      cerr << "Cannot open file " << inputFileName << endl;
      exit(0);
    }

  double time1 = 0, averageTime = 0;

  int cnt = 0;

  double div = 0., avDiv = 0.;

  SimData *sd = NULL;

  
  
  /* we either read FASTA (msflag == FALSE) or ms data (msflag == TRUE) */
  while( (msflag == false && (getAlignment(inputFile, data) == 1) ) || ( msflag == true && ( getMSData(inputFile, mspositions, msdata, haveOutgroup, outgroupSeq) == 1 ) ) )
    {

      time1 = gettime();

      ++cnt;

      if( msflag == true )
	{
	  
	  assert( mspositions.size() > 0 );

	  assert( msdata.size() > 0 );
	  
	  sd = new SimData(mspositions, msdata);
	  
	  assert(sd != NULL);
	  
	}
      else /* if msflag == false */
	{

	  if(lastIsOutgroup )
	    outgroupSeq = data.size() - 1;
	  
	  assert( outgroupSeq < data.size() );
	  
	  assert(Alignment::IsAlignment(data));
	  
	  if (Alignment::Gapped(data))
	    {
	      Alignment::RemoveTerminalGaps(data);
	    }
	  
	}
      
      
      
      if(sepPops && config != NULL)
	{

	  for( int pop = 0; pop < npop; ++pop)
	    {
	      
	      
	      if(header == 1)
		{

		  printf("ThetaPi%d\tThetaW%d\tVarPi%d\tStochasticVarPi%d\tSamplingVarPi%d\tVarThetaW%d\tNumPoly%d\tNumMutations%d\tNumSingletons%d\tTajimaD%d\tHprime%d\tDnominator%d\tFuLiDStar%d\tFuLiFStar%d\tDandVH%d\tDandVK%d\tWallsB%d\tWallsBprime%d\tWallsQ%d\tHudsonsC%d", pop, pop, pop, pop, pop, pop,pop, pop, pop,pop, pop, pop,pop, pop, pop,pop, pop, pop, pop, pop) ;

		  if( divergence == true )
		    printf( "\tDivergence%d\tAverageDivergence%d", pop, pop );

		  
		  if(haveOutgroup)
		    printf( "\tThetaH%d\tThetaL%d\tFuLiD%d\tFuLiF%d", pop, pop, pop, pop); 

		  
		  if(maxDisRel == true)
		    printf("\tmaxDisRel%d", pop);

		  if(sfsflag == true )
		    {
		      
		      for( int j = 0; j < config[pop] + 1; ++j)
			printf("\tsfs_%d_%d", pop, j);
		    }
		    

		  printf("\t");		  
		  
		}
	      
	    }
	}

      
      if(header == 1)
	{

	  if( config != NULL && pairwiseFst == true)
	    {

	      for(int pop1 = 0; pop1 < npop-1; ++pop1)
		{
		  for( int pop2 = pop1 + 1; pop2 < npop; ++pop2)
		    {
		      cout << "Fst_" << pop1 << "_" << pop2 << "\t";
		    }
		}
	    }
	}


      
      if(header == 1)
	{

	  if(windowsAnalysis == true)
	    cout << "window\t";
	  
	  cout << "ThetaPi\tThetaW\tVarPi\tStochasticVarPi\tSamplingVarPi\tVarThetaW\tNumPoly\tNumMutations\tNumSingletons\tTajimaD\tHprime\tDnominator\tFuLiDStar\tFuLiFStar\tDandVH\tDandVK\tWallsB\tWallsBprime\tWallsQ\tHudsonsC";

	  if(divergence == true )
	    {
	      cout << "\tDivergence\tAverageDivergence";
	    }

	  if(haveOutgroup)
	    cout << "\tThetaH\tThetaL\tFuLiD\tFuLiF";
	  
	  if(config != NULL)
	    cout << "\tFst";

	  if( maxDisRel == true)
	    cout << "\tmaxDisRel";

	  if(sfsflag == true)
	    {
	      
	      totalSampleSize = data.size();
	      
	      if( haveOutgroup )
		totalSampleSize--;	      
	      	      
	      sfs = (unsigned int *)calloc( totalSampleSize + 1, sizeof(unsigned int) );
	      
	      for(int i = 0; i <totalSampleSize + 1; ++i)
		printf("\tsfs_%d", i);
	      
	    }

	  cout <<endl;
	  
	  ++header;
	}
      
      
      
      if(sepPops && config != NULL)
	{
	  for( int pop = 0; pop < npop; ++pop)
	    {

	      PolySites *polytable = NULL;
	      
	      PolySNP *analyze = NULL;

	      SimData *simDataPop = NULL;


	      if( msflag == false )
		{
		  /* if there is outgroup, this is in the end of each subset */
		  getSubSet( pop, npop, haveOutgroup, outgroupSeq, config, data, subset);
		  
		  
		  if (Alignment::Gapped(subset))
		    Alignment::RemoveTerminalGaps(subset);
		  
		  if( sfsflag == true && haveOutgroup)
		    {

		      vector< unsigned > unfoldedSFS = ufs(subset, subset.size() - 1);
		      for( int i = 0; i < unfoldedSFS.size(); ++i)
			{
			  sfs[i] = unfoldedSFS[i];
			  //cerr << i << "\t" << sfs[i]<< endl;
			}
		      
		    }
		  
		  if( divergence == true)
		    {
		      div = 0.;
		      
		      Kimura80 *k80 = NULL;
		      
		      unsigned int randomSequence = rand() % (subset.size() - 1);
		      
		      k80 = new Kimura80( &subset[ randomSequence ], &subset[ subset.size() - 1] );
		      
		      double div = k80 -> K();
		      
		      
		      if(div == 999.0)
			{
			  div = NAN;
			}
		      
		      delete k80;
		      
		      
		      avDiv = 0.;
		      
		      for(unsigned int i = 0; i < subset.size() - 1; ++i )
			{
			  Kimura80 *k80 = new Kimura80( &subset[i], &subset[subset.size() - 1] );
			  avDiv += k80 -> K();
			}
		      
		      avDiv /= ( subset.size() - 1 );
		      
		      
		    }

		  polytable = new PolySites(subset);
		  
		  
		  if( haveOutgroup == true )
		    analyze = new PolySNP(polytable, true, config[pop]);
		  else
		    analyze = new PolySNP(polytable, false, 0);
		  
		  assert(polytable != NULL);

		}
	      else // if we have ms data
		{
		  getSubSetMS(pop, npop, config, msdata, subsetms, false, 0);

		  fprintf(stderr, "msdata size: %lu\n", subsetms.size());
		  
		  simDataPop = new SimData( mspositions, subsetms);
		  
		  analyze = new PolySNP( simDataPop );

		  // calculate the sfs now
		  
		  SimData::const_site_iterator itr = simDataPop -> sbegin();
		  
		  if( sfsflag == true)
		    {
		      // reset the sfs
		      for( int i = 0; i < totalSampleSize + 1; ++i)
		  	sfs[i] = 0;

		      while( itr < simDataPop -> send() )
		  	{
			  
		  	  unsigned ones = count(itr->second.begin(), itr->second.end(), '1');

			  assert(ones >= 0);

			  assert(ones < totalSampleSize + 1);
		  	  sfs[ones]++;

		  	  ++itr;
		  	}

		    }
		 
		}
	    
	      
	      //assert(polytable != NULL);
		  
      	      
	      
	      unsigned int maxDistance = 0;

	      if(maxDisRel == true && msflag == false)
		maxDistance = getMaxDistance( subset, false, -1);

	      
	      double thetaPi = analyze->ThetaPi();

	    	      
	      cout << thetaPi << "\t";
	      
	      cout << analyze->ThetaW() << "\t";
	      	      
	      cout << analyze->VarPi() << "\t";
	      
	      cout << analyze->StochasticVarPi() << "\t";
	      
	      cout << analyze->SamplingVarPi() << "\t";
	      
	      cout << analyze->VarThetaW() << "\t";
	      
	      cout << analyze->NumPoly() << "\t";
	      
	      cout << analyze->NumMutations() << "\t";
	      
	      cout << analyze->NumSingletons() << "\t";
	      
	      cout << analyze->TajimasD() << "\t";
	      
	      cout << analyze->Hprime() << "\t";
	      
	      cout << analyze->Dnominator() << "\t";
	      	      
	      cout << analyze->FuLiDStar() << "\t";
	      
	      cout << analyze->FuLiFStar() << "\t";
	      
	      cout << analyze->DandVH() << "\t";
	      
	      cout << analyze->DandVK() << "\t";
	      
	      cout << analyze->WallsB() << "\t";
	      
	      cout << analyze->WallsBprime() << "\t";
	      
	      cout << analyze->WallsQ() << "\t";
	      
	      cout << analyze->HudsonsC() << "\t";

	      if( divergence == true )
		cout << div << "\t" << avDiv << "\t";

	      
	      if(haveOutgroup)
		{
		  
		  cout  <<  analyze->ThetaH() << "\t" ; 
		  
		  cout  << analyze->ThetaL() << "\t" ;
	  
		  cout  << analyze->FuLiD() << "\t" ;
		  
		  cout  << analyze->FuLiF() << "\t" ;
		}
      
	      
	      if(maxDisRel == true)
		cout << maxDistance/thetaPi << "\t";

	      if(sfsflag == true)
		for( int i = 0; i < config[pop] + 1; ++i)
		  {
		    cout << sfs[i] << "\t";
		    sfs[i] = 0;
		  }

	      subset.erase(subset.begin(), subset.end());

	      if(simDataPop != NULL)
		delete simDataPop;
	      
	      if( polytable != NULL)
		delete polytable;
	      
	      if( analyze != NULL)
		delete analyze;
	      
	      
	    }
	  
	  
	}


      if(pairwiseFst == true && config != NULL)
	{
	  unsigned int *pairs = (unsigned int*)calloc(2, sizeof(int));
	  
	  unsigned int *pairsConf = NULL;

	  FST *pairFST = NULL;

	  
	  for( int pop1 = 0; pop1 < npop-1; ++pop1)
	    {
	      pairs[0] = pop1;
	      
	      for( int pop2 = pop1+1; pop2 < npop; ++pop2)
		{

		  		  
		  pairs[1] = pop2;

		  PolySites *polytable = NULL;

		  SimData *simDataPairs = NULL;

		  if( msflag == false )
		    {
		      getSubSetGroups( pairs, 2, npop, false, 0, config, data, pairsData, &pairsConf); 
		      
		      
		      if (Alignment::Gapped(pairsData))
			Alignment::RemoveTerminalGaps(pairsData);
		      
		      polytable = new PolySites(pairsData);
		    }
		  else
		    {
		      getSubSetGroupsMS( pairs, 2, npop, config, msdata, pairsmsdata, &pairsConf);
		      
		      simDataPairs = new SimData( mspositions, pairsmsdata);

		      
		      
		      
		    }
		  
		  if( msflag == false)
		    pairFST = new FST(polytable, 2, pairsConf, NULL, false, 0);
		  else
		    pairFST = new FST(simDataPairs, 2, pairsConf, NULL, false, 0);

		  assert( pairsConf[0] + pairsConf[1] == pairsData.size());

		  cout << pairFST ->HBK() << "\t";

		  delete pairFST;

		  pairsData.erase( pairsData.begin(), pairsData.end());

		  delete polytable;
		  
		  
		}
	    }
	  
	  
	}
      

      PolySNP *analyze = NULL;
      
      PolySites *polytable = NULL;
      
      
      if(msflag == false )
	{
	  polytable = new PolySites(data);
        
	  analyze = new PolySNP(polytable, haveOutgroup, outgroupSeq);

	}
      else
	{
	  polytable = new PolySites( mspositions, msdata);
	  
	  analyze = new PolySNP(sd);
	}
      
      unsigned int maxDistance = 0;

      
      PolySites::const_site_iterator itr = polytable -> sbegin();
      
      if( sfsflag == true && msflag == true)
	{
	  // reset the sfs
	  for( int i = 0; i < totalSampleSize + 1; ++i)
	    sfs[i] = 0;
	  
	  while( itr < polytable -> send() )
	    {
	      
	      unsigned ones = count(itr->second.begin(), itr->second.end(), '1');
	      sfs[ones]++;
	      ++itr;

	    }
	}
      else if( sfsflag == true && msflag == false)
	{
	  
	  vector< unsigned > unfoldedSFS = ufs(data, outgroupSeq);
	  for( int i = 0; i < unfoldedSFS.size(); ++i)
	    {
	      sfs[i] = unfoldedSFS[i];
	      //cerr << i << "\t" << sfs[i]<< endl;
	    }
	  
	  
	}


      if( msflag == false && divergence == true )
	{
	  
	  int cnt1 = 0;
	  
	  for( unsigned int i = 0; i < data.size(); ++i)
	    {
	      if( haveOutgroup && i == outgroupSeq )
		continue;

	      cnt1++;

	      Kimura80 *k80A = new Kimura80( &data[i], &data[ outgroupSeq ] );
	      avDiv += k80A->K();

	      delete k80A;
	    }

	  avDiv /= cnt1;

	  //cerr << "average div: " << avDiv << endl;


	  unsigned int randomSequence = outgroupSeq;

	  
	  while( randomSequence == outgroupSeq )
	    {
	      // get a random sequence 
	      randomSequence = rand() % (data.size());
	    }
	  
	  Kimura80 *k80A = new Kimura80( &data[randomSequence], &data[ outgroupSeq ] );
	  
	  div = k80A -> K();
	  
	  //int divsites = k80A->sites();

	  //cerr << "div: " << div << endl;

	  delete k80A;

	  
	}
      
      if(msflag == false && maxDisRel == true)
	maxDistance = getMaxDistance( data, haveOutgroup, outgroupSeq);
      
      FST *fst = NULL; 
      
      if( msflag == false && config != NULL)
	fst = new FST(polytable, npop, config, NULL, haveOutgroup, outgroupSeq);
      else if( msflag == true && config != NULL)
	fst = new FST(polytable, npop, config, NULL, false, 0);


      PolyTableSlice<PolySites> *windows = NULL; 
      if ( slidingWindowMode == true )
	{
	  windows = new PolyTableSlice<PolySites> ( polytable->sbegin(), 
						    polytable->send(), 
						    10000, 
						    1000, 
						    1000000., 
						    1000000.);
	  
	  
      
      cerr << "windows size: " << windows->size() << endl;
	}
      
      
      if(slidingWindowMode == true && windowsAnalysis == true)
	cout << 0 << "\t";

      
      double thetaPi = analyze->ThetaPi();

      cout << analyze->ThetaPi() << "\t";
      
      cout << analyze->ThetaW() << "\t";
            
      cout << analyze->VarPi() << "\t";
      
      cout << analyze->StochasticVarPi() << "\t";
      
      cout << analyze->SamplingVarPi() << "\t";
      
      cout << analyze->VarThetaW() << "\t";
      
      cout << analyze->NumPoly() << "\t";
      
      cout << analyze->NumMutations() << "\t";
      
      cout << analyze->NumSingletons() << "\t";
      
      cout << analyze->TajimasD() << "\t";
      
      cout << analyze->Hprime() << "\t";
      
      cout << analyze->Dnominator() << "\t";
      
      cout << analyze->FuLiDStar() << "\t";
      
      cout << analyze->FuLiFStar() << "\t";
      
      cout << analyze->DandVH() << "\t";
      
      cout << analyze->DandVK() << "\t";
      
      cout << analyze->WallsB() << "\t";
      
      cout << analyze->WallsBprime() << "\t";
      
      cout << analyze->WallsQ() << "\t";
      
      cout << analyze->HudsonsC();

      
      if( divergence == true )
	cout << div << "\t" << avDiv << "\t";

      
      if(haveOutgroup)
	{
            
	  cout <<"\t" <<  analyze->ThetaH() ;
	  
	  cout << "\t" << analyze->ThetaL() ;
	  
	  cout << "\t" << analyze->FuLiD() ;
	  
	  cout << "\t" << analyze->FuLiF() ;
	}
            
      if(config != NULL)
	cout << "\t" << fst -> HBK() << "\t";

      if(maxDisRel  == true )
	cout << "\t" << maxDistance/thetaPi;
      
      
      
      if(sfsflag == true)
	for( int i = 0; i < totalSampleSize + 1; ++i)
	  {
	    cout << "\t" << sfs[i];
	    sfs[i] = 0;
	  }

  
      
      cout << endl;

      


      for( unsigned i = 0; slidingWindowMode == true && i < windows -> size(); ++i)
	{
	  
	  
	  cout << i+1 << "\t";

	  PolySites w( (*windows)[i] );
	  PolySNP an(&w);
	  
	  cout << an.ThetaPi() << "\t";
	  
	  cout << an.ThetaW() << "\t";
	  
	  cout << an.VarPi() << "\t";
	  
	  cout << an.StochasticVarPi() << "\t";
	  
	  cout << an.SamplingVarPi() << "\t";
	  
	  cout << an.VarThetaW() << "\t";
	  
	  cout << an.NumPoly() << "\t";
	  
	  cout << an.NumMutations() << "\t";
	  
	  cout << an.NumSingletons() << "\t";
	  
	  cout << an.TajimasD() << "\t";
	  
	  cout << an.Hprime() << "\t";
	  
	  cout << an.Dnominator() << "\t";
	  
	  cout << an.FuLiDStar() << "\t";
	  
	  cout << an.FuLiFStar() << "\t";
	  
	  cout << an.DandVH() << "\t";
	  
	  cout << an.DandVK() << "\t";
	  
	  cout << an.WallsB() << "\t";
	  
	  cout << an.WallsBprime() << "\t";
	  
	  cout << an.WallsQ() << "\t";
	  
	  cout << an.HudsonsC() << endl;
	  

	}


      delete polytable;
      
      delete analyze;
      
      data.erase( data.begin(), data.end());

      averageTime += (gettime() - time1);

      if( cnt % 10 == 0)
	fprintf(stderr, "current average time (%d): %f\n", cnt,  averageTime/cnt);
    }

  fprintf(stderr, "Total time: %f, average time: %f\n", gettime() - time0, averageTime/cnt);

  exit(1);
}
