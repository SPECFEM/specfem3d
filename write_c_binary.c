/*!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 2
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================*/

// after Brian's function

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static int fd;

void
open_file_(char *file) {
  /*    fprintf(stderr, "Opening file: %s\n", file); */
  fd = open(file, O_WRONLY | O_CREAT, 0644);
  if(fd == -1) {
    fprintf(stderr, "Error opening file: %s exiting\n", file);
    exit(-1);
  }
}

void
close_file_() {
  /*    fprintf(stderr, "Closing file\n"); */
  close(fd);
}

void
write_integer_(int *z) {
  write(fd, z, sizeof(int));
}

void
write_real_(float *z) {
  write(fd, z, sizeof(float));
}

