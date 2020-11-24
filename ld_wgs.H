#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_multifit.h>


#define VERSION "0.1 -- Oct. 2020"

// function declarations
using namespace std;


// struct definitions

struct datacont{
  gsl_matrix * G; // genotype matrix, rows = nloci, cols =  nind
  gsl_matrix * Gf; // genotype matrix for focal snps
  gsl_vector * g1;
  gsl_vector * g2; // data for one SNP
  gsl_matrix * X; // transpose of Gf, used for multiple regression
  gsl_vector_int * sub; // binary, use individual in LD calculation
  int nloci;
  int nfocal;
  int nind;
  int subset;
  

};

struct paramcont{
  gsl_vector * mld; // max ld
  gsl_vector_int * wmld; // SNP number for max LD
  gsl_vector_int * sign; // pos = 1, neg = -1 sign of LD
  gsl_vector * betas; // regression coefficient esimtaes
  gsl_matrix * covar; // regression covariance matrix, p x p, where p = nfocal
  double chi2; // regression chi2
  double tss; // total sum of squares
  double r2; // r2
  
};

// function declarations
void usage(char * name);
void readData(string infile, datacont * data);
void readFocal(string infile, datacont * data);
void readSub(string infile, datacont * data);
void subG(datacont * data);
void calcLd(datacont * data, paramcont * params, int locus);
void calcMLd(datacont * data, paramcont * params, int locus);
void calcMfitLd(datacont * data, paramcont * params, int locus, FILE * OUT);
