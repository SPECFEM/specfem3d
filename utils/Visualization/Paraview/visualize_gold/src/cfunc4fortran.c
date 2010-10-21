/* Functions below are intended to be called from fortran routines to write 
binary files. 

- First, any number of files (usually 1) are opened using:
  call open_file(filename, fd)
- Any number and any kind (presently string, int, and float) of data are then
  written to the file/s corresponding to the file descriptor fd using, e.g.,
  call write_string(string,fd) or
  call write_integer(integer,fd) or
  call write_float(float,fd) 
- Then the file/s are closed using:
  close_file(fd)
! COMPILE
!   >> gcc -c cfunc4fortran.c
! HISTORY:
!   Hom Nath Gharti, NORSAR
!   Mar 10,2010 (NORSAR)
! FEEDBACK:
!   homnath_AT_norsar_DOT_no
!------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

/* fname: file name */
/* fd   : file descriptor */
 
/* open file to read as a file descriptor (not FILE pointer!)*/
void open_file2read_(char *fname, int *fd) {
  *fd = open(fname, O_RDONLY, 0);
  if(*fd == -1) {
    fprintf(stderr, "ERROR: file %s cannot be opened!\n", fname);
    exit(-1);
  }  
}

/* open file to write as a file descriptor (not FILE pointer!)*/
void open_file2write_(char *fname, int *fd) {
  *fd = open(fname, O_WRONLY | O_APPEND | O_CREAT | O_TRUNC, 0644);
  if(*fd == -1) {
    fprintf(stderr, "ERROR: file %s cannot be opened!\n", fname);
    exit(-1);
  }  
}

/* open file to append as a file descriptor (not FILE pointer!)*/
void open_file2append_(char *fname, int *fd) {
  *fd = open(fname, O_WRONLY | O_APPEND, 0644);
  if(*fd == -1) {
    fprintf(stderr, "ERROR: file %s cannot be opened!\n", fname);
    exit(-1);
  }  
}

/* close file */
void close_file_(int *fd) {
  close(*fd);
}

/* delete file */
void delete_file_(char *fname) {
  unlink(fname);
}

/* close and delete file */
void close_delete_file_(char *fname, int *fd) {
  close(*fd);
  unlink(fname);
}

/* write a string */
void write_string_(char *string, int *fd) {  
  int slen;
  slen=strlen(string);
  write(*fd, string, slen);
}


/* write a integer */
void write_integer_(int *z, int *fd) {
  write(*fd, z, sizeof(int));
}

/* write a float */
void write_float_(float *z, int *fd) {
  write(*fd, z, sizeof(float));
}

/* write a real */
void write_real_(float *z, int *fd) {
  write(*fd, z, sizeof(float));
}

/* read a integer */
void read_integer_(int *z, int *fd) {
  read(*fd, z, sizeof(int));
}

/* read a float */
void read_float_(float *z, int *fd) {
  read(*fd, z, sizeof(float));
}

/* This function determines the byte order of the processor architecture
source: http://www.ibm.com/developerworks/aix/library/au-endianc/index.html?ca=drs-
Use a character pointer to the bytes of an int and then check its first byte to 
see if it is 0 or 1.
HISTORY
  Hom Nath Gharti
  Mar 18,2009 (Princeton University)
*/
#define LE 0 /* Little Endian */
#define BE 1 /* Big Endian */

void get_endian_(int *z) {
  int i = 1;
  char *p = (char *)&i;

  if (p[0] == 1)
    *z = LE;
  else
    *z = BE;
}

