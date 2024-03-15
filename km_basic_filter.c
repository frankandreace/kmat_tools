#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdbool.h>
#include <string.h>

/*
THIS FUNCTION TAKES AS INPUT A MATRIX OF KMERS (KMERS ARE ROWS, SAMPLES ARE COLUMNS) AND 
RETURN A NEW MATRIX THAT SATISFIES THE DESIRED CONDITIONS:

1 - A MINIMUM ABUNDANCE OF A KMER TO BE CONSIDERED PRESENT IN A SAMPLE
2 - A MINIMUM NUMBER OF SAMPLES IN WHICH A KMER HAS TO BE ABSENT FROM
3 - A MINIMUM NUMBER OF SAMPLES IN WHICH A KMER HAS TO HAVE ABUNDANCE > THAN THE ONE SET IN THE 1ST CONDITION

CONDITION 2 AND 3 CAN BE SPECIFIED AS ABSOLUTE VALUE OR FRACTION OF THE SAMPLES.
INPUT AND OUTPUT CAN BE SPECIFIED OR BE STDIN/OUT
*/


int main(int argc, char **argv) {

  int min_zeros=10, min_nz=10, min_abund=10;
  double min_zero_frac=0.5, min_nz_frac=0.1;
  char *out_fname = NULL;
  bool verbose_opt=false, help_opt=false;
  
  bool min_zero_frac_opt=false, min_nz_frac_opt=false;

  int c;
  while ((c = getopt(argc, argv, "a:f:F:n:N:o:vh")) != -1) {
    switch (c) {
      case 'a':
        min_abund = strtol(optarg, NULL, 10);
        break;
      case 'n':
        min_zeros = strtol(optarg, NULL, 10);
        break;
      case 'N':
        min_nz = strtol(optarg, NULL, 10);
        break;
      case 'o':
        out_fname = optarg;
        break;
      case 'f':
        min_zero_frac_opt = true;
        min_zero_frac = atof(optarg);
        break;
      case 'F':
        min_nz_frac_opt = true;
        min_nz_frac = atof(optarg);
        break;
      case 'v':
        verbose_opt = true;
        break;
      case 'h':
        help_opt = true;
        break;
      case '?':
        return 1;
      default:
        abort();
    }
  }

  if(min_zero_frac_opt && (min_zero_frac < 0.01 || min_zero_frac >0.99)) {
    fprintf(stderr, "[error] -f must be in the [0.01,0.99] interval.\n");
    return 1;
  }
  if(min_nz_frac_opt && (min_nz_frac < 0.01 || min_nz_frac > 0.95)) {
    fprintf(stderr, "[error] -F must be in the [0.01,0.95] interval.\n");
    return 1;
  }

  if(argc-optind != 1 || help_opt) {
    fprintf(stdout, "Usage: km_basic_filter [options] <in.mat>\n\n");

    fprintf(stdout, "Filter a matrix by selecting k-mers that are potentially differential.\n\n");
    
    fprintf(stdout, "Options:\n");
    fprintf(stdout, "  -a INT    min abundance to define a k-mer as present in a sample [10]\n");
    fprintf(stdout, "  -n INT    min number of samples for which a k-mer should be absent [10]\n");
    fprintf(stdout, "  -f FLOAT  fraction of samples for which a k-mer should be absent (overrides -n)\n");
    fprintf(stdout, "  -N INT    min number of samples for which a k-mer should be present [10]\n");
    fprintf(stdout, "  -F FLOAT  fraction of samples for which a k-mer should be present (overrides -N)\n");
    fprintf(stdout, "  -o FILE   output filtered matrix to FILE [stdout]\n");
    fprintf(stdout, "  -v        verbose output\n");
    fprintf(stdout, "  -h        print this help message\n");
    return 0;
  }

  // set the matrix file pointer to the standard input or to the path provided in the arguments (non option)
  FILE *matfile = strcmp(argv[optind],"-") ? fopen(argv[optind],"r") : stdin;
  if(matfile == NULL) { 
    fprintf(stderr,"[error] cannot open file \"%s\"\n",argv[optind]);
    return 1;
  }

  // set the output file pointer to stdout or to the path provided in the options
  FILE *outfile = out_fname ? fopen(out_fname,"w") : stdout;
  if(outfile != stdout && outfile == NULL) {
    if(matfile != stdin){ fclose(matfile); }
    fprintf(stderr,"[error] cannot open output file \"%s\"\n",out_fname);
    return 1;
  }

  size_t n_samples = 0, n_kmers = 0, n_retrieved = 0;

  // use two char arrays to store the line to be parsed. the 1st to parse the line, the second to write it on the output if to be retained
  char *line = NULL, *line_cpy = NULL;
  size_t line_size = 0, line_cpy_size = 0;

  // read a line and place into the char buffer
  ssize_t ch_read = getline(&line, &line_size, matfile);
  while(ch_read >= 0) {  // while char in the line 

    // if the line copy array is smaller than the line one, resize it. Needed for the strcpy
    if(line_cpy_size < line_size) {
      line_cpy_size = line_size;
      line_cpy = (char*)realloc(line_cpy,line_cpy_size);
    }
    // copy line into line_cpy
    line_cpy = strcpy(line_cpy,line);

    // tockenize (split) the string based on ' '(space), '\t'(tab) or '\n'(newline) 
    char *elem = strtok(line," \t\n"); 
    if(elem == NULL){ continue; } // skip empty lines (if no token returned)
    
    // IF LINE NOT EMPTY 
    
    // record seeing a new kmer 
    ++n_kmers;

    size_t n_zeros = 0, n_present = 0;
    
    // While there are elements in the current line
    while((elem = strtok(NULL," \t\n")) != NULL) { 

      // when reading the first kmer, recod how many columns (samples) are in the matrix
      if(n_kmers == 1){ 
        ++n_samples;  
      }

      // convert element into a base_10 long integer
      long val = strtol(elem,NULL,10);
      
      // add zeros (if == 0 ) or present (if > minimum abundance)
      if(val == 0){ ++n_zeros; } else if(val >= min_abund){ ++n_present; }
    }

    // check conditions to record the kmer row
    
    // check if there were enough kmers with zero in the row (based on frequency or absolute count)
    bool enough_zeros = (min_zero_frac_opt && n_zeros >= min_zero_frac*n_samples) || (!min_zero_frac_opt && n_zeros >= min_zeros);

    // check if there were enough kmers with values grater than the minimum abundance in the row (based on frequency or absolute count)
    bool enough_nz = (min_nz_frac_opt && n_present >= min_nz_frac*n_samples) || (!min_nz_frac_opt && n_present >= min_nz);

    // in case both options are true, the row is retrieved and written to the output filestream
    if(enough_zeros && enough_nz) {
      ++n_retrieved;
      fputs(line_cpy, outfile);
    }

    if(verbose_opt && (n_kmers & ((1U<<20)-1)) == 0) {
      fprintf(stderr, "%lu k-mers processed, %lu retrieved\n", n_kmers, n_retrieved);
    }

    // get the next line to be processed
    ch_read = getline(&line, &line_size, matfile);
  }

  fprintf(stderr, "[info] %lu\tsamples\n", n_samples);
  fprintf(stderr, "[info] %lu\ttotal k-mers\n", n_kmers);
  fprintf(stderr, "[info] %lu\tretained k-mers\n", n_retrieved);

  // cleaning up - freeing the memory for the arrays and closing the filestreams
  free(line); free(line_cpy);
  if(matfile != stdin){ fclose(matfile); }
  if(outfile != stdout){ fclose(outfile); }

  return 0;
}
