// file: main.C

// ld_wgs = LD (r2) between genome wide SNPs and a set of focal SNPs
// for use with Timema cristinae within-generation selection experiment
// infiles contains non-integer genotype estimates for genome-wide SNPs and
// focal SNPs

// Time-stamp: <Thursday, 08 October 2020, 15:45 MDT -- zgompert>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <time.h>
#include <getopt.h>
#include "ldwin_wgs.H"

using namespace std;

/* ----------------- */
/* beginning of main */
/* ----------------- */

int main(int argc, char *argv[]) {

  time_t start = time(NULL);
  time_t end;
  int win, ch, o;
  
  string infile = "undefined";
  string posfile = "undefined";
  string outfile = "out_ld.txt";
  paramcont params;
  datacont data;

  FILE * OUT; // outfile
  

  // default dimension of data, need to enter
  data.nloci = 100;
  data.nind = 100;
  data.winsize = 10000;
  
  // get command line arguments
  if (argc < 2) {
    usage(argv[0]);
  }
  
  while ((ch = getopt(argc, argv, "i:j:g:o:s:w:")) != -1){
    switch(ch){
    case 'i':
      data.nloci = atoi(optarg);
      break;
    case 'j':
      data.nind = atoi(optarg);
      break;
    case 'g':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 's':
      posfile = optarg;
      break;
    case 'w':
      data.winsize = atoi(optarg);
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
  cerr << "Reading SNP position file" << endl;
  readPos(posfile, &data);

    
  // allocate memory
  data.g1 = gsl_vector_calloc(data.nind);
  data.g2 = gsl_vector_calloc(data.nind);
  
  // loop over genome-wide SNPs
  cerr << "Computing LD" << endl;
  win = 0;
  cerr << win << " " << data.winsize << " " << data.max << endl;
  while((win+data.winsize-1)<data.max){
    cerr << "Working on window: " << win << endl;
    o = subG(&data, win);
    if(o==1){
      calcLd(&data, &params);
      fprintf(OUT, "%i %.4f\n", win, params.mld); 
    }
    else{
      fprintf(OUT, "%i NA\n", win);
    }
    win += data.winsize;

  }

  fclose(OUT);
  
  // prints run time
  end = time(NULL);
  cerr << "Runtime: " << (end-start)/3600 << " hr " << (end-start)%3600/60 << " min ";
  cerr << (end-start)%60 << " sec" << endl;
  return 0;
}
 
