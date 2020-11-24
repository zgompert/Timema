// file: main.C

// estpEM - a program to estimate allele frequencies from
// pre-calculated genotype likelihoods using a soft EM algorithm

// Time-stamp: <Friday, 03 October 2014, 08:16 MDT -- zgompert>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <time.h>
#include <getopt.h>
#include "estpEM.H"

using namespace std;

gsl_rng * r;  /* global state variable for random number generator */

/* ----------------- */
/* beginning of main */
/* ----------------- */

int main(int argc, char *argv[]) {

  time_t start = time(NULL);
  time_t end;
  int rng_seed = 0;
  int ch, i;
  int nh = 0; // number of header lines, excluding the required, nind nlocus, header

  double tolerance = 0.001; // maximum allele frequency difference for convergence
  int max = 20; // maximum number of EM iterations before giving up

  string infile = "undefined";
  string outfile = "out.txt";

  datacont data;

  gsl_matrix * alleleFreqs;

  FILE * OUT;

  // get command line arguments
  if (argc < 2) {
    usage(argv[0]);
  }
  
  while ((ch = getopt(argc, argv, "i:o:e:m:h:")) != -1){
    switch(ch){
    case 'i':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'e':
      tolerance = atof(optarg);
      break;
    case 'm':
      max = atoi(optarg);
      break;
    case 'h':
      nh = atoi(optarg);
      break;
    case '?':
    default:
      usage(argv[0]);
    }
  }

  // open outfile
  OUT = fopen(outfile.c_str(),"w");

  // set up gsl random number generation 
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);
  srand(time(NULL));
  rng_seed = rand();
  gsl_rng_set(r, rng_seed); /* seed gsl_rng with output of rand, which
			       was seeded with result of time(NULL) */

  // read genotype likelihood formatted data
  cerr << "Reading data from " << infile << endl;
  readData(infile, &data, nh);

  // allocate memory for allele frequency estimate matrix
  alleleFreqs = gsl_matrix_calloc(data.nloci, 2); // initial and EM estimates, respectively

  // use EM algorithm to estimate allele frequencies
  cerr << "Using EM algorithm to estimate allele frequencies" << endl;
  for(i=0; i<data.nloci; i++){
    expMax(&data,alleleFreqs,tolerance,max,i);
  }

  // write results
  cerr << "Writing results to " << outfile << endl;
  write(&data,alleleFreqs,OUT);
  fclose(OUT);

  // prints run time
  end = time(NULL);
  cerr << "Runtime: " << (end-start)/3600 << " hr " << (end-start)%3600/60 << " min ";
  cerr << (end-start)%60 << " sec" << endl;
  return 0;
}
