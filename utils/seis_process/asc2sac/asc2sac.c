
#include <stdio.h>
#include "formats.h"

int main(int argc, char*argv[]) {
  struct sac *sac_hdr;
  int npts;
  float time, amp, *data, dt;
  float bval;
  FILE *fp;
  int i;
  char *filename;
  char *comp[8];

  if(argc < 4) {
    fprintf(stderr, "%s ascii-file npts sac-file\n", argv[0]);
    exit(0);
  }
  npts = atoi(argv[2]);
  filename = argv[3];
  if((fp = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Error opening ascii-file: %s\n", argv[1]);
    exit(-1);
  }
  data = (float *) malloc(sizeof(float) * npts);
  i = 0;
  while(fscanf(fp, "%f %f", &time, &amp) == 2) {
    if(i > npts) {
      fprintf(stderr, "Number of points read in (%d) is greater than specified (%d)\n",
        i, npts);
      exit(-1);
    }
    if(i == 0) {
      bval = time;
    }
    if(i == 1) {
       dt = time - bval;
    }
    data[i] = amp;
    i++;
  }
  if(i != npts) {
    fprintf(stderr, "Number of points read in (%d) is less than specified (%d)\n",
      i, npts);
    exit(-1);
  }
  fclose(fp);
  /* Must set parameters */
  sac_hdr = &sac_null;      /* Sets header to a null header */
  sac_hdr->npts = npts; /* number of points */
  sac_hdr->delta = dt; /* Time between points */
  sac_hdr->leven = 1;     /* Evenly spaced data */
  sac_hdr->iftype = 1;    /* Type of File */
  sac_hdr->internal4 = 6; /* Internal Variable */
  sac_hdr->b = bval;         /* Starting time */
  /* Open and write a newly created  Sac file */
  /*      fprintf(stderr, "%s\n", filename); */
  if((fp = sacopen(filename, "w")) == NULL) {
    fprintf(stderr, "Error opening sac-file: %s\n", filename);
    exit(-1);
  }
  sachdrwrite(fp, sac_hdr);
  sacdatawrite(fp, data, npts);
  sacclose(fp);
  return(0);
}


