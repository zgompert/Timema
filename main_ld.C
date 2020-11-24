// file: main.C

// ld_wgs = LD (r2) between genome wide SNPs and a set of focal SNPs
// for use with Timema cristinae within-generation selection experiment
// infiles contains non-integer genotype estimates for genome-wide SNPs and
// focal SNPs

// oct 22. added an option to compute r2 from multiple regression of each SNP
// on all causal variants

// Time-stamp: <Wednesday, 14 October 2020, 13:35 MDT -- zgompert>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <time.h>
#include <getopt.h>
#include "ld_wgs.H"

using namespace std;

/* ----------------- */
/* beginning of main */
/* ----------------- */

int main(int argc, char *argv[]) {

  time_t start = time(NULL);
  time_t end;
  int i, ch;
  int mregmod, efitmod;
  
  string infile = "undefined";
  string focalfile = "undefined";
  string subfile = "undefined";
  string outfile = "out_ld.txt";
  paramcont params;
  datacont data;

  FILE * OUT; // outfile
  

  // default dimension of data, need to enter
  data.nloci = 100;
  data.nfocal = 100;
  data.nind = 100;
  data.subset = 0; // binary, only work with subset

  // set model, default is max ld, set to 1 to run multiple regression
  mregmod = 0;
  efitmod = 0; // set to 1 to run regression over posterior for expected fitness
  
  // get command line arguments
  if (argc < 2) {
    usage(argv[0]);
  }
  
  while ((ch = getopt(argc, argv, "i:j:n:g:o:f:s:b:m:e:")) != -1){
    switch(ch){
    case 'i':
      data.nloci = atoi(optarg);
      break;
    case 'j':
      data.nind = atoi(optarg);
      break;
    case 'n':
      data.nfocal = atoi(optarg);
      break;
    case 'g':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'f':
      focalfile = optarg;
      break;
    case 's':
      data.subset = atoi(optarg);
      break;
    case 'b':
      subfile = optarg;
      break;
    case 'm':
      mregmod = atoi(optarg);
      break;
    case 'e':
      efitmod = atoi(optarg);
      break;
    case '?':
    default:
      usage(argv[0]);
    }
  }

  // open outfile
  OUT = fopen(outfile.c_str(),"w");
  
  // read data files
  cerr << "Reading in genome-wide SNPs file" << endl;
  readData(infile, &data);
  cerr << "Reading in focal SNP set file" << endl;
  readFocal(focalfile, &data);
  if(data.subset == 1){
    cerr << "Reading subset to use for calculations" << endl;
    readSub(subfile, &data);
    // shrink matrixes and reset nind
    subG(&data);
  }
    
  // allocate memory
  data.g1 = gsl_vector_calloc(data.nind);

  if(mregmod == 0){
    data.g2 = gsl_vector_calloc(data.nind);
    params.mld = gsl_vector_calloc(data.nloci);
    params.wmld = gsl_vector_int_calloc(data.nloci);
    params.sign = gsl_vector_int_calloc(data.nloci);
  }
  else if(mregmod == 1 && efitmod == 0){
    data.X = gsl_matrix_calloc(data.nind, data.nfocal);
    params.betas = gsl_vector_calloc(data.nfocal);
    params.covar = gsl_matrix_calloc(data.nfocal, data.nfocal);
    gsl_matrix_transpose_memcpy(data.X, data.Gf);
  }
  else if(efitmod == 1){
    data.X = gsl_matrix_calloc(data.nind, 2);
    params.betas = gsl_vector_calloc(2);
    params.covar = gsl_matrix_calloc(2, 2);
    gsl_matrix_set_all(data.X, 1); // set to 1, leave first row at 1 as intercepts
  }
  
  // loop over genome-wide SNPs
  cerr << "Computing LD" << endl;
  for(i=0;i<data.nloci;i++){
    if((i % 10000) == 0)
      cerr << "Working on SNP: " << i << endl;
    if(mregmod == 0){// max LD model
      calcLd(&data, &params, i);
      fprintf(OUT, "%.4f %d %d\n", gsl_vector_get(params.mld, i), gsl_vector_int_get(params.wmld, i),
	    gsl_vector_int_get(params.sign, i));
    }
    else if(mregmod == 1 && efitmod == 0){// multiple regression r2
      calcMLd(&data, &params, i);
      fprintf(OUT, "%.4f\n", params.r2);
    }
    else if(efitmod == 1){// multiple regression r2 post on efit
      calcMfitLd(&data, &params, i, OUT);

    }
  }
  fclose(OUT);
  
  // prints run time
  end = time(NULL);
  cerr << "Runtime: " << (end-start)/3600 << " hr " << (end-start)%3600/60 << " min ";
  cerr << (end-start)%60 << " sec" << endl;
  return 0;
}
 
