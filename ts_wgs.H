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

};

struct paramcont{
  gsl_vector * S; // selection coefficient, vector of length MCMC steps for a single locus
  gsl_vector * w; // fitness values, ninds
  gsl_vector * wsum; // sum survival
  
  double Scur; // current value of S
  
  double alpha; // a and b for beta prior on S
  double beta;
  double prop; // +- uniform proposal for S
  
};

// function declarations
void usage(char * name);
void readData(string infile, datacont * data);
void readSB(string infile, datacont * data);
void initScur(datacont * data, paramcont *params, int locus);
void updateS(datacont * data, paramcont * params, int locus);
void calcFitness(datacont * data, paramcont * params, int locus, double s);
double calcLProb(datacont * data, paramcont * params, double s);
void writeS(paramcont * params, FILE * OUT, int locus);
