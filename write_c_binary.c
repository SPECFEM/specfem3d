/*
!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================
*/

// after Brian's function

#include "config.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static int fd;

void
FC_FUNC_(open_file,OPEN_FILE)(char *file) {
  /*    fprintf(stderr, "Opening file: %s\n", file); */
  fd = open(file, O_WRONLY | O_CREAT, 0644);
  if(fd == -1) {
    fprintf(stderr, "Error opening file: %s exiting\n", file);
    exit(-1);
  }
}

void
FC_FUNC_(close_file,CLOSE_FILE)() {
  /*    fprintf(stderr, "Closing file\n"); */
  close(fd);
}

void
FC_FUNC_(write_integer,WRITE_INTEGER)(int *z) {
  write(fd, z, sizeof(int));
}

void
FC_FUNC_(write_real,WRITE_REAL)(float *z) {
  write(fd, z, sizeof(float));
}

