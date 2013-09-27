

#include "formats.h"

int
sachdrwrite(FILE *fpout, struct sac *sac_hdr) {
  if(fwrite(sac_hdr, sizeof(struct sac), 1, fpout) != 1) {
    fprintf(stderr, "SacHeaderWrite: could not write sac header\n");
    return(-1);
  }
  return(1);
}
