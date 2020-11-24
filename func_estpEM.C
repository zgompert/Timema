#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <float.h>

#include <math.h>

#include "estpEM.H"

using namespace std;

void usage(char * name){
  fprintf(stderr,"\n%s version %s\n\n", name, VERSION); 
  fprintf(stderr, "Usage: popmod -i infile.txt [options]\n");
  fprintf(stderr, "-i Infile with genetic data for the population\n");
  fprintf(stderr, "-o Outfile for allele frequency estimates [default = out.txt]\n");
  fprintf(stderr, "-e Tolerance for EM convergence [default = 0.001]\n");
  fprintf(stderr, "-m Maximum number of EM iterations [default = 20]\n");
  fprintf(stderr, "-h Number of header lines [default = 0]\n");
  exit(1);
}

// function to read the genotype likelihoods
void readData(string infile, datacont * data, int nh){
  int i, j, k, l;
  ifstream file;
  istringstream stream;
  string line, element;
  double genotypes[NGEN];
  double gensum = 0;

  // open genotype likelihood file
  file.open(infile.c_str());
  if (!file){
    cerr << "Cannot open file " << infile << endl;
    exit(1);
  }

  // read line with number of loci and individuals
  getline(file, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of individuals
  data->nind = atoi(element.c_str());  
  stream >> element; // number of loci
  data->nloci = atoi(element.c_str()); 

  cerr << "Number of loci: " << data->nloci << endl;
  cerr << "Number of individuals: " << data->nind << endl;

  // remove excess header lines
  for(l=0; l<nh; l++){
    getline(file,line);
  }

  // allocate memory for data
  data->gl = new gsl_matrix * [data->nloci];
  for(i=0; i<data->nloci; i++){
    data->gl[i] = gsl_matrix_calloc(data->nind, NGEN);
  }
  data->lIds = new string[data->nloci];

  // genotype likelihood lines, this are phred-scaled
  for(i=0; i<data->nloci; i++){
    getline(file, line); // data for one locus
    stream.str(line);
    stream.clear();
    stream >> element; // locus id, this is not retained
    data->lIds[i] = element;
    for(j=0; j<data->nind; j++){
      gensum = 0;
      for(k=0; k<NGEN; k++){ // store genotype likelihoods
	stream >> element; // need to convert from phred scale: phred
	// = -10 log10(prob), prob = 10^(phred/-10)
	genotypes[k] = pow(10, (atof(element.c_str())/-10.0));
	gensum += genotypes[k];
      }
      for(k=0; k<NGEN; k++){ // normalize genotype likelihoods
	genotypes[k] = genotypes[k]/gensum;
	if (genotypes[k] == 0) // set to DBL_MIN if 0
	  genotypes[k] = DBL_MIN;
	gsl_matrix_set(data->gl[i], j, k, genotypes[k]);
      }
    }
  } 
}

// soft EM algorithm for allele frequency estimation
void expMax(datacont * data, gsl_matrix * af, double tol, int max, int i){
  
  int j, k, n;
  double p = 0, pml = 0, gml;
  double dif = 1;
  
  // initilize using simple point estimate
  for(j=0; j<data->nind; j++){
    for(k=1; k<NGEN; k++){ // can ignore 0 b/c 0 * gl = 0
      p += k * gsl_matrix_get(data->gl[i], j, k);
     }
  }
  p /= (2.0 * data->nind);
  gsl_matrix_set(af, i, 0, p);
  //cerr << endl << p << " :: ";
  // iterate EM algorithm
  for(n=0; (n<max) && (dif > tol); n++){
    for(j=0; j<data->nind; j++){
      gml = 0 + 2 * p * (1-p) * gsl_matrix_get(data->gl[i], j, 1) +
	p * p * gsl_matrix_get(data->gl[i], j, 2) * 2;
      gml /= (1-p) * (1-p) * gsl_matrix_get(data->gl[i], j, 0) + 
	2 * p * (1-p) * gsl_matrix_get(data->gl[i], j, 1) +
	p * p * gsl_matrix_get(data->gl[i], j, 2);
      pml += (gml/2.0);
    }
    if (n == max)
      cerr << "failed to converge: locus " << i << endl;
    pml /= data->nind;
    dif = pml - p;
    if (dif < 0)
      dif = dif * -1;
    p = pml;
    pml = 0;
    //cerr << p << ", ";
  }   
  //cerr << endl << endl;
  gsl_matrix_set(af, i, 1, p); // EM estimate

}

// write results to a text file
void write(datacont * data, gsl_matrix * af, FILE * OUT){

  int i;

  for(i=0; i<data->nloci; i++){
    fprintf(OUT, "%s %.4f %.4f\n", data->lIds[i].c_str(), gsl_matrix_get(af, i, 0),
	    gsl_matrix_get(af, i, 1));
  }
  
}
