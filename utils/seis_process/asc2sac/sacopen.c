
#include "formats.h"

FILE *
sacopen(char *infile, char * mode) {
  FILE *fp;
  if((fp = fopen(infile, mode)) == NULL) {
    fprintf(stderr, "SacOpen: cant open file: %s\n", infile);
    return(NULL);
  }
  return(fp);
}
