
#include "formats.h"

int
sacdatawrite(FILE *fpout, float *data_ptr, int npts) {
  if(fwrite(data_ptr, npts*sizeof(float), 1, fpout) != 1) {
    fprintf(stderr, "SacDataWrite: could not write sac data\n");
    return(-1);
  }
  return(1);
}
