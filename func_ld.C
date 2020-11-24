#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sort_vector.h>
#include <float.h>
#include <math.h>
#include "ld_wgs.H"

using namespace std;

void usage(char * name){
  fprintf(stderr,"\n%s version %s\n\n", name, VERSION); 
  fprintf(stderr, "Usage: ld [options]\n");
  fprintf(stderr, "-g Infile, matrix of genotype estimates [default = undefined]\n");
  fprintf(stderr, "-f Infile, matrix of focal genotypes [default = undefined]\n");
  fprintf(stderr, "-i Total number of genetic loci [default = 100]\n");
  fprintf(stderr, "-n Number of focal genetic loci [default = 100]\n");
  fprintf(stderr, "-j Total number of individuals [default = 100]\n");
  fprintf(stderr, "-o Text outfile for LD estimates [default = out_ld.txt]\n");
  fprintf(stderr, "-s Binary, use only subset of inds. [default = 0]\n");
  fprintf(stderr, "-b Infile with binary var. specifying inds. to keep [default = undefined]\n");
  fprintf(stderr, "-m Binary, use multiple regression model [default = 0]\n");
  fprintf(stderr, "-e Binary, use expected fitness model [default = 0]\n");
  
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

//function to read the genotype data from the focal SNPs file
void readFocal(string infile, datacont * data){
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
  data->Gf = gsl_matrix_calloc(data->nfocal, data->nind);
    
  // read file, use user input number of loci and individuals
  // store in matrix Gf
  for(i=0; i<data->nfocal; i++){
    getline(file, line); // data for one locus
    stream.str(line);
    stream.clear(); 
    for(j=0; j<data->nind; j++){
      stream >> element;
      gsl_matrix_set(data->Gf, i, j, atof(element.c_str()));
    }     
  }
}

//function to read binary indicator variable, retain ind. in calc of LD
void readSub(string infile, datacont * data){
  int j;
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
  data->sub = gsl_vector_int_calloc(data->nind);
    
  // read file, store binary in sub
  for(j=0; j<data->nind; j++){
    getline(file, line); // data for one locus
    stream.str(line);
    stream.clear(); 
    stream >> element;
    gsl_vector_int_set(data->sub, j, atoi(element.c_str()));
  }
}

//function to subset genotype matrixes and redefine nind
void subG(datacont * data){
  gsl_matrix * tempG;
  gsl_matrix * tempGf;
  int Nkeep;
  int i, j, x;
  double val;
  
  Nkeep = 0;
  for(j=0;j<data->nind;j++){
    Nkeep += gsl_vector_int_get(data->sub, j);
  }
  
  // create and fill temp matrixes
  tempG = gsl_matrix_calloc(data->nloci, data->nind);
  tempGf = gsl_matrix_calloc(data->nfocal, data->nind);  
  gsl_matrix_memcpy(tempG, data->G);
  gsl_matrix_memcpy(tempGf, data->Gf);
  // redefine main matrixes
  gsl_matrix_free(data->G);
  gsl_matrix_free(data->Gf);
  data->G = gsl_matrix_calloc(data->nloci, Nkeep);
  data->Gf = gsl_matrix_calloc(data->nfocal, Nkeep);

  // loops to fill matrixes
  x = 0;
  for(j=0;j<data->nind;j++){
    if(gsl_vector_int_get(data->sub, j)==1){ // keep this one
      for(i=0;i<data->nloci;i++){
	val = gsl_matrix_get(tempG, i, j);
	gsl_matrix_set(data->G, i, x, val);
      }
      for(i=0;i<data->nfocal;i++){
	val = gsl_matrix_get(tempGf, i, j);
	gsl_matrix_set(data->Gf, i, x, val);
      }
      x++;
    }
  }
  data->nind = Nkeep;
  gsl_matrix_free(tempG);
  gsl_matrix_free(tempGf);
  
}


// function to compute the max LD for a SNP
void calcLd(datacont * data, paramcont * params, int locus){
  int k;
  double maxld = 0;
  double ld;
  int sig = 0;
  int msig = 0;
  int wld = 0;

  gsl_matrix_get_row(data->g1, data->G, locus);
  
  for(k=0;k<data->nfocal;k++){
    gsl_matrix_get_row(data->g2, data->Gf, k);
    ld = gsl_stats_correlation(data->g1->data, data->g1->stride, data->g2->data,
			       data->g2->stride, data->g2->size);
    // set sign
    if(ld > 0)
      sig = 1;
    else
      sig = -1;
    
    ld = gsl_pow_2(ld);
    
    if(ld > maxld){
      maxld = ld;
      wld = k;
      msig = sig;
    }
    
  }
  gsl_vector_set(params->mld, locus, maxld);
  gsl_vector_int_set(params->wmld, locus, wld);
  gsl_vector_int_set(params->sign, locus, msig);
}

// function to compute multiple regression r2 for each snps
void calcMLd(datacont * data, paramcont * params, int locus){

  gsl_multifit_linear_workspace * work;
  
  gsl_matrix_get_row(data->g1, data->G, locus);

  // multiple regression
  // allocate n * p workspace
  work = gsl_multifit_linear_alloc(data->nind, data->nfocal);
  // ols fit
  gsl_multifit_linear(data->X, data->g1, params->betas,
			   params->covar, &params->chi2, work);

  // compute the total sum of squares
  params->tss = gsl_stats_tss(data->g1->data, data->g1->stride, data->g1->size);
  
  // computer r2
  params->r2 = 1.0 - (params->chi2/params->tss);

  // free memory
  gsl_multifit_linear_free(work);
}

// function to compute multiple regression r2 for each snps
void calcMfitLd(datacont * data, paramcont * params, int locus, FILE * OUT){

  int x;
  double val;
  gsl_multifit_linear_workspace * work;
  gsl_vector * vec;
  gsl_vector * r2pp;

  vec = gsl_vector_calloc(data->nind);
  r2pp = gsl_vector_calloc(data->nfocal);

  gsl_matrix_get_row(data->g1, data->G, locus);

  // multiple regression
  // allocate n * p workspace
  work = gsl_multifit_linear_alloc(data->nind, 2);

  // loop over post
  for(x=0;x<data->nfocal;x++){
    gsl_matrix_get_row(vec, data->Gf, x); // grab xth sample from post
    gsl_matrix_set_col(data->X, 1, vec); // set to column 1 in data/regression vec
  
    // ols fit
    gsl_multifit_linear(data->X, data->g1, params->betas,
			params->covar, &params->chi2, work);
    // compute the total sum of squares
    params->tss = gsl_stats_tss(data->g1->data, data->g1->stride, data->g1->size);
  
    // computer r2
    params->r2 = 1.0 - (params->chi2/params->tss);
    gsl_vector_set(r2pp, x, params->r2);
  }

  // point estimate and quantiles
  gsl_sort_vector(r2pp);
  val = gsl_stats_median_from_sorted_data(r2pp->data, r2pp->stride, r2pp->size);
  fprintf(OUT, "%.4f", val);
  // 5, 95, 2.5, 97.5
  val = gsl_stats_quantile_from_sorted_data(r2pp->data, r2pp->stride, r2pp->size, 0.05);
  fprintf(OUT, " %.4f", val);
  val = gsl_stats_quantile_from_sorted_data(r2pp->data, r2pp->stride, r2pp->size, 0.1);
  fprintf(OUT, " %.4f", val);
  val = gsl_stats_quantile_from_sorted_data(r2pp->data, r2pp->stride, r2pp->size, 0.025);
  fprintf(OUT, " %.4f", val);
  val = gsl_stats_quantile_from_sorted_data(r2pp->data, r2pp->stride, r2pp->size, 0.975);
  fprintf(OUT, " %.4f\n", val);

  
  // free memory
  gsl_multifit_linear_free(work);
  gsl_vector_free(vec);
  gsl_vector_free(r2pp);
}
