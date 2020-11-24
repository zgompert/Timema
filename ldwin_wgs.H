#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>


#define VERSION "0.1 -- Oct. 2020"

// function declarations
using namespace std;


// struct definitions

struct datacont{
  gsl_matrix * G; // genotype matrix, rows = nloci, cols =  nind
  gsl_matrix_view Gwin; // genotype submatrix for one window, as matrix view
  gsl_vector_int * pos; // vector of SNP positions

  gsl_vector * g1;
  gsl_vector * g2; // data for one SNP
  
  int nloci;
  int nind;
  int nsnps;
  int winsize;
  int max; 
  

};

struct paramcont{
  gsl_vector * ld; // ld for
  double mld; // mean LD for window
  
};

// function declarations
void usage(char * name);
void readData(string infile, datacont * data);
void readPos(string infile, datacont * data);
int subG(datacont * data, int lb);
void calcLd(datacont * data, paramcont * params);
