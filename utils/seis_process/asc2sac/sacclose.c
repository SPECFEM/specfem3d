
#include "formats.h"

int
sacclose(FILE *fp) {
  if(fclose(fp) != 0) {
  fprintf(stderr,"Sacclose: Could not open file\n");
    return(-1);
  }
  return(0);
}
