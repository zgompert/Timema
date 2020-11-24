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
#include <float.h>
#include <math.h>
#include "ts_wgs.H"

using namespace std;

void usage(char * name){
  fprintf(stderr,"\n%s version %s\n\n", name, VERSION); 
  fprintf(stderr, "Usage: wgtcr [options]\n");
  fprintf(stderr, "-g Infile, matrix of genotype estimates [default = undefined]\n");
  fprintf(stderr, "-d Infile with survival and block data [default = undefined]\n");
  fprintf(stderr, "-m Number of MCMC steps for the analysis [default = 1000]\n");
  fprintf(stderr, "-n Discard the first n MCMC samples as a burn-in [default = 100]\n");
  fprintf(stderr, "-t Thin MCMC samples by recording every nth value [default = 2]\n");
  fprintf(stderr, "-i Total number of genetic loci [default = 100]\n");
  fprintf(stderr, "-j Total number of individuals [default = 100]\n");
  fprintf(stderr, "-o Text outfile for parameter estimates [default = out_wgtcr.txt]\n");
  fprintf(stderr, "-a Parameter alpha for beta prior on S [default = 10]\n");
  fprintf(stderr, "-b Parameter beta for beta prior on S [default = 90]\n");
  fprintf(stderr, "-l Lower bound of range of loci to analyze [default = 0]\n");
  fprintf(stderr, "-u Upper bound of range of loci to analyze [default = 1000]\n");
  fprintf(stderr, "-p Proposal dist. +- for s [default = 0.2]\n");
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

//function to initialize Scur to a value near the point estimate, s = dp/(p*(1-p))
void initScur(datacont * data, paramcont *params, int locus){
  double p0 = 0, p1 = 0;
  double N2 = 0, N21 = 0;
  double dp, s;
  int j;

  for(j=0;j<data->nind;j++){
    if(gsl_vector_uint_get(data->surv, j) == 1){
      p1 += gsl_matrix_get(data->G, locus, j);
      N21 += 2;
    }
    p0 += gsl_matrix_get(data->G, locus, j);
    N2 +=  2;
  }
  p0 = p0/N2;
  p1 = p1/N21;
  dp = p1 - p0;
  s = -1 * dp/(p0 * (1-p0));
  //cerr << dp << ", " << s << ",";
  s += gsl_ran_flat(r,-0.02,0.02);
  
  if(s > 0.2)
    s = 0.2;
  if(s < -0.2)
    s = -0.2;
  params->Scur = s;
  
}

//function for MCMC update of S
void updateS(datacont * data, paramcont * params, int locus){
  double Sprop;
  double lpr1, lpr0; // new and old log prob
  
  // propose new value of S from beta prior
  Sprop = params->Scur + gsl_ran_flat(r, -1.0 * params->prop, params->prop);
  if(abs(Sprop) <= 0.5){ // otherwise Pr(S) = 0, don't need to do anything
   
    // calculate survival probs. under current S
    gsl_vector_set_zero(params->wsum);
    calcFitness(data, params, locus, params->Scur);
    lpr0 = calcLProb(data, params, params->Scur);

    // calculate survival probs. under proposed S
    gsl_vector_set_zero(params->wsum);
    calcFitness(data, params, locus, Sprop);
    lpr1 = calcLProb(data, params, Sprop);

    //cerr << Sprop << ":" << lpr1 << "    " << params->Scur << ":" << lpr0 << endl;
    // calc MH ratio and accept/reject
    if((lpr1-lpr0) > log(gsl_rng_uniform(r)))
      params->Scur = Sprop;
  }
}

// function to compute absolute fitness
void calcFitness(datacont * data, paramcont * params, int locus, double s){
  int j, bl;
  double val, cur;
  
  for(j=0;j<data->nind;j++){
    val = 1 - (gsl_matrix_get(data->G, locus, j) * s);
    gsl_vector_set(params->w, j, val);
    bl = gsl_vector_uint_get(data->block, j);
    cur = gsl_vector_get(params->wsum, bl);
    val = cur+val;
    //cerr << val << " " << bl << " " << gsl_vector_get(params->wsum, bl) << " xx " << endl;
    gsl_vector_set(params->wsum, bl, val);
  }
  // obtain absolute surv. probs
  for(j=0;j<data->nind;j++){
    val = gsl_vector_get(params->w, j);
    bl = gsl_vector_uint_get(data->block, j);
    // cerr << val <<  " " << bl << endl;
    // cerr << gsl_vector_uint_get(data->surv, j) << " " << gsl_vector_int_get(data->bsurv, bl) << endl;
    val = (val/gsl_vector_get(params->wsum,bl)) * gsl_vector_int_get(data->bsurv, bl);
    //cerr << "::: " << val << endl;
    gsl_vector_set(params->w, j, val);
  }
}

// function to compute log Pr(Y|S) Pr(S)
double calcLProb(datacont * data, paramcont * params, double s){
  int j;
  double pr = 0;
  double tpr, val;
  
  // lprob (Y|S)
  for(j=0;j<data->nind;j++){
    val = gsl_vector_get(params->w, j);
    tpr = gsl_ran_bernoulli_pdf(gsl_vector_uint_get(data->surv, j), val);
    if(tpr == 0)
      tpr = DBL_MIN;
    if(gsl_finite(tpr) == 0)
      tpr = DBL_MAX;
    pr += log(tpr);
  }


  // lprob (S)
  val = abs(s*2.0); // drop minus and double to get back to beta distribution
  tpr = gsl_ran_beta_pdf(val, params->alpha, params->beta);
  if(tpr == 0)
    tpr = DBL_MIN;
  if(gsl_finite(tpr) == 0)
    tpr = DBL_MAX;
  pr += log(tpr);
  return pr;
}

// function to write to outfile
void writeS(paramcont * params, FILE * OUT, int locus){
  double mn;
  double q50, q5, q95;
  double pp = 0;
  double n = 0;
  int x;
  
  mn = gsl_stats_mean(params->S->data, params->S->stride, params->S->size);
//  cerr << " " << mn << endl;
  gsl_sort_vector(params->S);
  q50 = gsl_stats_median_from_sorted_data(params->S->data, params->S->stride, params->S->size);
  q5 = gsl_stats_quantile_from_sorted_data(params->S->data, params->S->stride, params->S->size,0.05);
  q95 = gsl_stats_quantile_from_sorted_data(params->S->data, params->S->stride, params->S->size,0.95);
  for(x=0; x<int(params->S->size); x++){
    n++;
    if(gsl_vector_get(params->S, x) > 0)
      pp++;
  }
  pp = pp/n;
  
  fprintf(OUT, "%i %.4f %.4f %.4f %.4f %.3f\n", locus, mn, q50, q5, q95, pp);
  
}
