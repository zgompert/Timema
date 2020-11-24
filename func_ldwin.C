#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics_double.h>
#include <float.h>
#include <math.h>
#include "ldwin_wgs.H"

using namespace std;

void usage(char * name){
  fprintf(stderr,"\n%s version %s\n\n", name, VERSION); 
  fprintf(stderr, "Usage: wgtcr [options]\n");
  fprintf(stderr, "-g Infile, matrix of genotype estimates [default = undefined]\n");
  fprintf(stderr, "-s Infile, position for each SNP [default = undefined]\n");
  fprintf(stderr, "-i Total number of genetic loci [default = 100]\n");
  fprintf(stderr, "-j Total number of individuals [default = 100]\n");
  fprintf(stderr, "-w Window size in bp [default = 10000]\n");
  fprintf(stderr, "-o Text outfile for LD estimates [default = out_ld.txt]\n");
 

  exit(1);
}

//function to read the genotype data
void readData(string infile, datacont * data){
  int i, j;
  ifstream file;
  istringstream stream;
  string line, element;

  // open file with genotype matrix
  file.open(infile.c_str());
  if (!file){
    cerr << "Cannot open the genotype matrix infile " << infile << endl;
    exit(1);
  }

  // allocate data matrix
  data->G = gsl_matrix_calloc(data->nloci, data->nind);
    
  // read file, use user input number of loci and individuals
  // store in matrix G
  for(i=0; i<data->nloci; i++){
    getline(file, line); // data for one locus
    stream.str(line);
    stream.clear(); 
    for(j=0; j<data->nind; j++){
      stream >> element;
      gsl_matrix_set(data->G, i, j, atof(element.c_str()));
    }     
  }
}

//function to the SNP positions data
void readPos(string infile, datacont * data){
  int i;
  ifstream file;
  istringstream stream;
  string line, element;

  // open file with genotype matrix
  file.open(infile.c_str());
  if (!file){
    cerr << "Cannot open the genotype matrix infile " << infile << endl;
    exit(1);
  }

  // allocate data matrix
  data->pos = gsl_vector_int_calloc(data->nloci);
    
  // read file, store binary in sub
  for(i=0; i<data->nloci; i++){
    getline(file, line); // data for one locus
    stream.str(line);
    stream.clear(); 
    stream >> element;
    gsl_vector_int_set(data->pos, i, atoi(element.c_str()));
  }
  data->max = atoi(element.c_str());
}

//function to subset genotype matrixes and redefine nind
int subG(datacont * data, int lb){
  int ub;
  int a = -1, b = 0;
  int i;
  double val;
  int minSnps = 4;
  int inrange = 0;
  
  ub = lb + data->winsize - 1; // define upper bound of window
  for(i=0;i<data->nloci;i++){
    val = gsl_vector_int_get(data->pos, i);
    if(val < lb)
      a = i;
    if(val < ub)
      b = i;
    if(val >= lb && val <= ub)
      inrange = 1;
  }
  a++; // add one so that you are at first SNP of window
  b = b-a+1; // convert ub to size

  data->nsnps = b;
  
  if(b > minSnps && inrange ==1){ // must have at least minSnps SNP to create the view
    data->Gwin = gsl_matrix_submatrix(data->G, a, 0, b, data->nind);
    return(1);
  }
  else{
    return(0);
  }
}


// function to compute the max LD for a SNP
void calcLd(datacont * data, paramcont * params){
  int i, k, x;
  double ld;
  int d;

  d = (data->nsnps * (data->nsnps-1))/2.0; // number of unique pairwise combos
  
  params->ld = gsl_vector_calloc(d);
  
  x = 0;
  for(i=0;i<(data->nsnps-1);i++){
    gsl_matrix_get_row(data->g1, &data->Gwin.matrix, i);
    for(k=(i+1);k<data->nsnps;k++){
      gsl_matrix_get_row(data->g2, &data->Gwin.matrix, k);
      
      ld = gsl_stats_correlation(data->g1->data, data->g1->stride, data->g2->data,
				 data->g2->stride, data->g2->size);
      ld = gsl_pow_2(ld);
      gsl_vector_set(params->ld, x, ld);
      x++;
    }
  }
  params->mld = gsl_stats_mean(params->ld->data, params->ld->stride, params->ld->size);
    
  gsl_vector_free(params->ld);

}

