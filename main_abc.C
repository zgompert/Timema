// file: main.C

// wgs_abc = approximate Bayesian inference of direct selection on SNPs from the
// Timema cristinae within-generation selection experiment
// infile contains non-integer genotype estimates for one treatment (A or C)

// Time-stamp: <Saturday, 24 October 2020, 20:38 MDT -- zgompert>

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
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <time.h>
#include <getopt.h>
#include "wgs_abc.H"

using namespace std;

gsl_rng * r;  /* global state variable for random number generator */

/* ----------------- */
/* beginning of main */
/* ----------------- */

int main(int argc, char *argv[]) {

  time_t start = time(NULL);
  time_t end;
  int rng_seed = 0, ch = 0, step = 0, nsteps = 1000;
  int j;
  
  string infile = "undefined";
  string designfile = "undefined";
  string prefix = "out_abc";
  paramcont params;
  datacont data;

  FILE * OUT; // outfile
  FILE * OUTNG; // out for qtl
  FILE * OUTS; // out for s
  FILE * OUTFIT; // optional outfile for expected fitness
  
  // hard-code # of blocks
  data.nblock = 5;

  // default dimension of data, need to enter
  data.nloci = 100;
  data.nind = 100;
  data.mdist = 0.5; // minimum distance to write

  // set default values for priors
  params.alpha = 10;
  params.beta = 90; // beta prios on S
  params.mu = 3; // prior mean for ngamma or prior max
  params.cv = 0; // default to no cross-validation
  params.wfit = 0; // default to no expected fitness outfile
  
  // get command line arguments
  if (argc < 2) {
    usage(argv[0]);
  }
  
  while ((ch = getopt(argc, argv, "n:i:j:g:d:o:a:b:m:l:c:f:")) != -1){
    switch(ch){
    case 'n':
      nsteps = atoi(optarg);
      break;
    case 'i':
      data.nloci = atoi(optarg);
      break;
    case 'j':
      data.nind = atoi(optarg);
      break;
    case 'g':
      infile = optarg;
      break;
    case 'd':
      designfile = optarg; // block and survival
      break;
    case 'o':
      prefix = optarg;
      break;
    case 'a':
      params.alpha = atof(optarg);
      break;
    case 'b':
      params.beta = atof(optarg);
      break;
    case 'm':
      params.mu = atof(optarg);
      break;
    case 'l':
      data.mdist = atof(optarg);
      break;
    case 'c':
      params.cv = atoi(optarg);
      break;
    case 'f':
      params.wfit = atoi(optarg);
      break;
    case '?':
    default:
      usage(argv[0]);
    }
  }

  // outfiles
  string outmain = prefix + ".hyp.txt";
  string outng = prefix + ".ngamma.txt";
  string outs = prefix + ".s.txt";
  string outfit = prefix + ".fit.txt";
  
  // open outfile
  OUT = fopen(outmain.c_str(),"w");
  OUTNG = fopen(outng.c_str(),"w");
  OUTS = fopen(outs.c_str(),"w");
  if(params.wfit == 1)
    OUTFIT = fopen(outfit.c_str(),"w");
  
  // set up gsl random number generation 
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);
  srand(time(NULL));
  rng_seed = rand();
  gsl_rng_set(r, rng_seed); /* seed gsl_rng with output of rand, which
			       was seeded with result of time(NULL) */

  // read data files
  cerr << "Reading in genotype file" << endl;
  readData(infile, &data);
  cerr << "Reading in design file" << endl;
  readSB(designfile, &data);

  // allocate memory
  data.lvec = gsl_vector_uint_calloc(data.nloci);
  
  // allocate memory for params 
  params.w = gsl_vector_calloc(data.nind);
  params.wsum = gsl_vector_calloc(data.nblock);
  params.ssurv = gsl_vector_uint_calloc(data.nind);
  if(params.cv == 1){
    params.perm = gsl_permutation_calloc(data.nind);
    gsl_permutation_init(params.perm);
    params.sets = gsl_vector_int_calloc(data.nind);
    for(j=0; j<floor(data.nind * 0.8); j++){
      gsl_vector_int_set(params.sets, j, 1);
    }
  }
    
  cerr << "Begining ABC simulations" << endl;
  
  // main loop for simulations
  for(step=0;step<nsteps;step++){ // mcmc loop
    oneSim(&data, &params, OUT, OUTNG, OUTS, OUTFIT);
      if((step % 1000) == 0)
	cerr << "Sim. " << step << endl;
  }
  fclose(OUT);
  fclose(OUTNG);
  fclose(OUTS);
  if(params.wfit == 1)
    fclose(OUTFIT);
  
  // prints run time
  end = time(NULL);
  cerr << "Runtime: " << (end-start)/3600 << " hr " << (end-start)%3600/60 << " min ";
  cerr << (end-start)%60 << " sec" << endl;
  return 0;
}
 
