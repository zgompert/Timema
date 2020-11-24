#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort_vector.h>

#define VERSION "0.1 -- Sept. 2020"

extern gsl_rng * r;
// function declarations
using namespace std;


// struct definitions

struct datacont{
  gsl_matrix * G; // genotype matrix, rows = nloci, cols =  nind
  gsl_vector_uint * surv; // boolean, survival
  gsl_vector_uint * block; // block number for each ind.
  int nloci;
  int nind;
  int nblock;
  gsl_vector_int * bsurv; // number surviving for each block, length nblock
  gsl_vector_uint * lvec; // vector of integers 0 ..  nloci-1, used for sampling
  double mdist; // minimum distance to write
  
};

struct paramcont{
  gsl_vector * w; // fitness values, ninds
  gsl_vector * wsum; // sum survival
  gsl_vector_uint * ssurv; // simulated survival

  int ngamma; // number of qtl
  double sbar; // mean selection coefficient

  unsigned int mu; // prior mean for ngamma (exponential) or max of uniform
  double alpha; // a and b for beta prior on sbar
  double beta;

  double dist; // distance (prop.) of differences in survival
  double distcv; // distrance (prop.) for test set, used only in CV mode

  int cv; // boolean, run cross-validation mode
  int wfit; // boolean, write exptected fitness
  gsl_permutation * perm; // permutation used to shuffle test vs. training set
  gsl_vector_int * sets; // vector of 0 and 1 specifiying training (1) and test (0) sets
  // allocation is 80% train and 20% test
};

// function declarations
void usage(char * name);
void readData(string infile, datacont * data);
void readSB(string infile, datacont * data);
void oneSim(datacont * data, paramcont * params, FILE * OUT, FILE * OUTNG, FILE * OUTS, FILE * OUTFIT);
void calcFitness(datacont * data, paramcont * params, gsl_vector * s, gsl_vector_uint * qtl);
void simSurv(datacont * data, paramcont * params);
double calcDist(datacont * data, paramcont * params);
double calcDistCv(datacont * data, paramcont * params);
void writeOut(datacont * data, paramcont * params, gsl_vector_uint * qtl, gsl_vector * s,
	      FILE * OUT, FILE * OUTNG, FILE * OUTS, FILE * OUTFIT);
