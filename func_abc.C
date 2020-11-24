#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <float.h>
#include <math.h>
#include "wgs_abc.H"

using namespace std;

void usage(char * name){
  fprintf(stderr,"\n%s version %s\n\n", name, VERSION); 
  fprintf(stderr, "Usage: wgtcr [options]\n");
  fprintf(stderr, "-g Infile, matrix of genotype estimates [default = undefined]\n");
  fprintf(stderr, "-d Infile with survival and block data [default = undefined]\n");
  fprintf(stderr, "-n Number of simulations [default = 1000]\n");
  fprintf(stderr, "-i Number of genetic loci [default = 100]\n");
  fprintf(stderr, "-j Number of individuals [default = 100]\n");
  fprintf(stderr, "-o Prefix for outfiles [default = out_abc]\n");
  fprintf(stderr, "-a Parameter alpha for beta prior on sbar [default = 10]\n");
  fprintf(stderr, "-b Parameter beta for beta prior on sbar [default = 90]\n");
  fprintf(stderr, "-m Parameter max number of causal variants for uniform [default = 3]\n");
  fprintf(stderr, "-c Run cross-validation (train:test = 80:20) [default = 0]\n");
  fprintf(stderr, "-f Write outfile with expected fitness [default = 0]\n");
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

//function to read survival and block data
void readSB(string infile, datacont * data){
  int j;
  ifstream file;
  istringstream stream;
  string line, element;
  double val;

  // open file with genotype matrix
  file.open(infile.c_str());
  if (!file){
    cerr << "Cannot open the survival/block infile " << infile << endl;
    exit(1);
  }

  // allocate data vectors
  data->surv = gsl_vector_uint_calloc(data->nind);
  data->block = gsl_vector_uint_calloc(data->nind);
  data->bsurv = gsl_vector_int_calloc(data->nblock);
  // read file and store survival (column 0, bin) and block (column 1, 1-5)
  for(j=0; j<data->nind; j++){
    getline(file, line); // data for one individual
    stream.str(line);
    stream.clear(); 
    stream >> element;
    gsl_vector_uint_set(data->surv, j, atoi(element.c_str()));
    stream >> element;
    gsl_vector_uint_set(data->block, j, atoi(element.c_str()));

    val = gsl_vector_uint_get(data->surv,j) +
      gsl_vector_int_get(data->bsurv, atoi(element.c_str()));
    gsl_vector_int_set(data->bsurv, atoi(element.c_str()), val);
  }
}

//function for running a ABC simulation
void oneSim(datacont * data, paramcont * params, FILE * OUT, FILE * OUTNG, FILE * OUTS,
	    FILE * OUTFIT){
  gsl_vector_uint * qtl;
  gsl_vector * s; 
  double val;
  int * lvec;
  lvec = new int [data->nloci]; 
  int i;

  for(i=0; i<data->nloci; i++){
    lvec[i] = i;
  }
  
  // sample number of causal variants
  if(gsl_rng_uniform(r) < 0.5)
    params->ngamma = gsl_rng_uniform_int(r, params->mu);
  else
    params->ngamma = 0;
  
  // sample locations for causal variants
  if(params->ngamma > 0){
    qtl = gsl_vector_uint_calloc(params->ngamma);
    gsl_ran_shuffle(r, lvec, data->nloci, sizeof (int));
    for(i=0; i<params->ngamma; i++){
      gsl_vector_uint_set(qtl, i, lvec[i]);
    }
    // sample fitness effects
    // average for exp. dist of S
    params->sbar = gsl_ran_beta(r, params->alpha, params->beta);
    s = gsl_vector_calloc(params->ngamma);
    for(i=0;i<params->ngamma;i++){
      val = gsl_ran_exponential(r, params->sbar);
      if(val > 0.5)
	val = 0.5;
      if(gsl_rng_uniform(r) < 0.5)
	val = val * -1.0;
      gsl_vector_set(s, i, val);
    }
  }
  else{ // no causal variants, thus average effect = 0 by def.
    params->sbar = 0;
    qtl = gsl_vector_uint_calloc(1);
    s = gsl_vector_calloc(1);
  }
  // compute fitness, multiplicative model
  calcFitness(data, params, s, qtl);
  
  // simulate survival
  simSurv(data, params);
  
  // calculate summary statistic
  if(params->cv == 0)
    params->dist = calcDist(data, params);
  else
    params->dist = calcDistCv(data, params); // version with cross-validation
  
  // write results
  if(params->dist <= data->mdist)
    writeOut(data, params, qtl, s, OUT, OUTNG, OUTS, OUTFIT);
  
  // free dynamic memory
  gsl_vector_uint_free(qtl);
  gsl_vector_free(s);
  delete[] lvec;

}

// function to compute absolute fitness, mult. model
void calcFitness(datacont * data, paramcont * params, gsl_vector * s, gsl_vector_uint * qtl){
  int j, i, l, bl;
  double val, cur;

  // zero wsum
  gsl_vector_set_zero(params->wsum);
  
  // relative fitnesses
  for(j=0;j<data->nind;j++){
    if(params->ngamma > 0){
      val = 1;
      for(i=0;i<params->ngamma;i++){
	l = gsl_vector_uint_get(qtl, i);
	val *= 1 - (gsl_matrix_get(data->G, l, j) * gsl_vector_get(s,i));
      }
      gsl_vector_set(params->w, j, val);
    }
    else{
      gsl_vector_set(params->w, j, 1);
    }
    bl = gsl_vector_uint_get(data->block, j);
    cur = gsl_vector_get(params->wsum, bl);
    val = cur+gsl_vector_get(params->w, j);
    gsl_vector_set(params->wsum, bl, val);
  }
  // obtain absolute surv. probs
  for(j=0;j<data->nind;j++){
    val = gsl_vector_get(params->w, j);
    bl = gsl_vector_uint_get(data->block, j);
    val = (val/gsl_vector_get(params->wsum,bl)) * gsl_vector_int_get(data->bsurv, bl);
    gsl_vector_set(params->w, j, val);
  }
}

// function to randomly sample survivors based on abs. fitness
void simSurv(datacont * data, paramcont * params){
  int j;
  double val;

  for(j=0; j<data->nind; j++){
    val = gsl_ran_bernoulli(r, gsl_vector_get(params->w, j));
    gsl_vector_uint_set(params->ssurv, j, val);
  }

}

// function to compute the distance between the obs. and sim. survivor vectors
double calcDist(datacont * data, paramcont * params){
  int j;
  double d = 0.0;
  
  for(j=0; j<data->nind; j++){
    if(gsl_vector_uint_get(data->surv, j) != gsl_vector_uint_get(params->ssurv, j))
      d += 1.0;
  }
  d = d/data->nind;
  return d;
}

// function to compute the distance between the obs. and sim. survivor vectors CV version
double calcDistCv(datacont * data, paramcont * params){
  int j;
  double d = 0.0;
  double dcv = 0.0;
  double c = 0, ccv = 0;
  
  // permutation for train and test sets
  gsl_ran_shuffle(r, params->perm->data, data->nind, sizeof(size_t));
  gsl_permute_vector_int(params->perm, params->sets);
  
  for(j=0; j<data->nind; j++){
    if(gsl_vector_int_get(params->sets, j) == 1){ // training
      if(gsl_vector_uint_get(data->surv, j) != gsl_vector_uint_get(params->ssurv, j))
	d += 1.0;
      c++;
    }
    else{ // test
      if(gsl_vector_uint_get(data->surv, j) != gsl_vector_uint_get(params->ssurv, j))
	dcv += 1.0;
      ccv++;
    }
  }
  d = d/c;
  params->distcv = dcv/ccv;
  return d;
}

// function to write to outfile
void writeOut(datacont * data, paramcont * params, gsl_vector_uint * qtl, gsl_vector * s,
	      FILE * OUT, FILE * OUTNG, FILE * OUTS, FILE * OUTFIT){
  int i, j;
  double val, bl;

  // hyperparms
  if(params->cv == 0)
    fprintf(OUT, "%i %.4f %.4f\n", params->ngamma, params->sbar, params->dist);
  else
    fprintf(OUT, "%i %.4f %.4f %.4f\n", params->ngamma, params->sbar, params->dist, params->distcv);
  
  if(params->ngamma == 0){
    fprintf(OUTNG, "NA\n");
    fprintf(OUTS, "NA\n");
  }
  else{
    fprintf(OUTNG, "%i", params->ngamma);
    fprintf(OUTS, "%i", params->ngamma);
    for(i=0;i<params->ngamma;i++){
      fprintf(OUTNG, " %i", gsl_vector_uint_get(qtl, i));
      fprintf(OUTS, " %.4f", gsl_vector_get(s, i));
    }
    fprintf(OUTNG, "\n");
    fprintf(OUTS, "\n");
  }
  if(params->wfit == 1){ // write exptected fitness
    for(j=0;j<data->nind;j++){
      val = gsl_vector_get(params->w, j);
      bl = gsl_vector_uint_get(data->block, j);
      val = val / gsl_vector_int_get(data->bsurv, bl);
      if(j==0)
	fprintf(OUTFIT, "%.4f", val);
      else
	fprintf(OUTFIT, " %.4f", val);
    }
    fprintf(OUTFIT, "\n");
  }
  
}
