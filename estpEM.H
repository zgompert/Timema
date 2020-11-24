#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

// macro definitions
#define VERSION "0.1 -- 2 October 2014"
#define NGEN 3 // number of genotypes, set to 3 for diploid inds. and bi-allelic loci

extern gsl_rng * r;
// function declarations
using namespace std;


// struct definitions

struct datacont{
  gsl_matrix ** gl; // genotype likelihoods; locus x ind x genotype; Hom Ref, Het, Hom Alt
  string * lIds;
  int nloci;
  int nind;
};


// function definitions
void usage(char *);
void readData(string infile, datacont * data, int nh);
void expMax(datacont * data, gsl_matrix * af, double tol, int max, int i);
void write(datacont * data, gsl_matrix * af, FILE * OUT);
