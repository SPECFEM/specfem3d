#include <stdio.h>
#include <stdlib.h>
//#include "drw_ascfile.h"
// Qinya Liu, Caltech, May 2007


void dread_ascfile_c_(const char *ascfile,
		   double *t0, double *dt, int *n,
		   double *data) 

{
  FILE *fd;
  double junk,junk1;
  int i;

  if ((fd = fopen(ascfile,"r")) == NULL) {
    printf(" file %s cannot be opened\n",ascfile);
    exit(1);
  }
  i = 0;
  while ( fscanf(fd,"%lf  %lf\n",&junk, data+i) != EOF ) {
    if (i == 0) junk1 = junk;
    if (i == 1) *dt = junk - junk1;
    i++;
  }
  *t0 = junk1;
  *n = i;
  if (fclose(fd) != 0) {
    printf(" file %s cannot be closed\n",ascfile);
    exit(1);
  }

}

void dwrite_ascfile_c_(const char *ascfile,
		    double *t0, double *dt, int *n,
		    const double *data) 

{
  FILE *fd;
  int i;
  double tmp,tmp0;
  
  //printf("t0: %f ,dt: %f ,n: %i \n",*t0,*dt,*n);
  
  if ((fd = fopen(ascfile,"w")) == NULL) {
    printf(" file %s cannot be opened to write\n",ascfile);
    exit(1);
  }
  i = 0;
  tmp = *dt;
  tmp0 = *t0;
  for (i=0; i< *n; i++) {  
    fprintf(fd,"%14.7g %18.7g\n", tmp0+i*tmp, data[i]);
  }
  if (fclose(fd) != 0) {
    printf("file %s cannot be closed\n",ascfile);
    exit(1);
  }

}
