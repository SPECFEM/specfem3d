/*
 * =====================================================================================
 *
 *       Filename:  enc.c
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  09/28/2010 03:58:30 PM JST
 *       Revision:  none
 *       Compiler:  mpicc
 *
 *         Author:  Leonardo BAUTISTA GOMEZ (leobago@matsulab.is.titech.com),
 *        Company:  Tokyo Institue of Technology
 *
 * =====================================================================================
 */

#include <python2.4/Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "mpi.h"

#include "galois.h"
#include "fti.h"

void enc_region_xor(char *r1, char *r2, char *r3, int nbytes)
{
  long *l1, *l2, *l3, *ltop;
  char *ctop;

  ctop = r1 + nbytes;
  ltop = (long *) ctop;
  l1 = (long *) r1;
  l2 = (long *) r2;
  l3 = (long *) r3;

  while (l1 < ltop) {
    *l3 = ((*l1)  ^ (*l2));
    l1++;
    l2++;
    l3++;
  }
}

void enc_w08_region_multiply(char *region, int multby, int nbytes, char *r2, int add)
{
  unsigned char *ur1, *ur2, *cp;
  unsigned char prod;
  int i, srow, j;
  unsigned long l, *lp2;
  unsigned char *lp;
  int sol;

  ur1 = (unsigned char *) region;
  ur2 = (r2 == NULL) ? ur1 : (unsigned char *) r2;

  if (galois_mult_tables[8] == NULL) {
    if (galois_create_mult_tables(8) < 0) {
      fprintf(stderr, "galois_08_region_multiply -- couldn't make multiplication tables\n");
      exit(1);
    }
  }
  srow = multby * nw[8];
  if (r2 == NULL || !add) {
    for (i = 0; i < nbytes; i++) {
      prod = galois_mult_tables[8][srow+ur1[i]];
      ur2[i] = prod;
    }
  } else {
    sol = sizeof(long);
    lp2 = &l;
    lp = (unsigned char *) lp2;
    for (i = 0; i < nbytes; i += sol) {
      cp = ur2+i;
      lp2 = (unsigned long *) cp;
      for (j = 0; j < sol; j++) {
        prod = galois_mult_tables[8][srow+ur1[i+j]];
        lp[j] = prod;
      }
      *lp2 = (*lp2) ^ l;
    }
  }
  return;
}

void enc_w16_region_multiply(char *region, int multby, int nbytes, char *r2, int add)
{
  unsigned short *ur1, *ur2, *cp;
  int prod;
  int i, log1, j, log2;
  unsigned long l, *lp2, *lptop;
  unsigned short *lp;
  int sol;

  ur1 = (unsigned short *) region;
  ur2 = (r2 == NULL) ? ur1 : (unsigned short *) r2;
  nbytes /= 2;

  if (multby == 0) {
    if (!add) {
      lp2 = (unsigned long *) ur2;
      ur2 += nbytes;
      lptop = (unsigned long *) ur2;
      while (lp2 < lptop) { *lp2 = 0; lp2++; }
    }
    return;
  }

  if (galois_log_tables[16] == NULL) {
    if (galois_create_log_tables(16) < 0) {
      fprintf(stderr, "galois_16_region_multiply -- couldn't make log tables\n");
      exit(1);
    }
  }
  log1 = galois_log_tables[16][multby];

  if (r2 == NULL || !add) {
    for (i = 0; i < nbytes; i++) {
      if (ur1[i] == 0) {
        ur2[i] = 0;
      } else {
        prod = galois_log_tables[16][ur1[i]] + log1;
        ur2[i] = galois_ilog_tables[16][prod];
      }
    }
  } else {
    sol = sizeof(long)/2;
    lp2 = &l;
    lp = (unsigned short *) lp2;
    for (i = 0; i < nbytes; i += sol) {
      cp = ur2+i;
      lp2 = (unsigned long *) cp;
      for (j = 0; j < sol; j++) {
        if (ur1[i+j] == 0) {
          lp[j] = 0;
        } else {
          log2 = galois_log_tables[16][ur1[i+j]];
          prod = log2 + log1;
          lp[j] = galois_ilog_tables[16][prod];
        }
      }
      *lp2 = (*lp2) ^ l;
    }
  }
  return;
}


void FTI_Encode(char *filename, int fs, int maxFs) {

    char *myData, *data, *coding, fn[25], efn[25], mfn[25];
    int cnt, pos, ps, l, i, j, init, src, offset, dest;
    MPI_Request reqSend, reqRecv;
    MPI_Status status;
    FILE *fd, *efd;
    double tim;
    tim = MPI_Wtime();

    // Calculating padding size
    ps = ((maxFs/FTI_Bksz))*FTI_Bksz;
    if (ps < maxFs)
        ps = ps + FTI_Bksz;

    // Initializing encoded variables
    myData = talloc(char, FTI_Bksz);
    data = talloc(char, 2*FTI_Bksz);
    coding = talloc(char, FTI_Bksz);

    // Openning files
    sprintf(fn,"%s/%s",FTI_Cdir, filename);
    sprintf(efn,"%s/%s.enc",FTI_Cdir, filename);
    if (truncate(fn,ps) == -1)
        FTI_Prerror("Error with truncate on checkpoint file");
    fd = fopen(fn, "rb");
    efd = fopen(efn, "wb");
    if (fd == NULL || efd == NULL)
        FTI_Prerror("Impossible to open the checkpoint file or encoded checkpoint file");

    pos = 0;
    while(pos < ps) {

        // Reading checkpoint files
        fread(myData, sizeof(char), FTI_Bksz, fd);

        // Encoding work
        dest = FTI_Idgr;
        i = FTI_Idgr;
        offset = 0;
        init = 0;
        cnt = 0;
        while(cnt < FTI_Grsz) {

            // Overlaping calculation and communication
            if (cnt == 0) {
                // The first loop we just copy the read data
                memcpy(&(data[offset*FTI_Bksz]), myData, sizeof(char)*FTI_Bksz);
            } else {
                // After first loop we need to receive the data
                MPI_Wait(&reqSend, &status);
                MPI_Wait(&reqRecv, &status);
            }

            if (cnt != FTI_Grsz-1) {
                // At every loop but the last one we send the data
                dest = (dest+FTI_Grsz-1)%FTI_Grsz;
                src = (i+1)%FTI_Grsz;
                //printf("proc %d i %d src %d dest %d\n",FTI_Idgr, i, src, dest);
                MPI_Isend(myData, FTI_Bksz, MPI_CHAR, dest, FTI_Mtag, FTI_Cenc, &reqSend);
                MPI_Irecv(&(data[(1-offset)*FTI_Bksz]), FTI_Bksz, MPI_CHAR, src, FTI_Mtag, FTI_Cenc, &reqRecv);
            }

            // First copy or xor any data that does not need to be multiplied by a factor
            if (FTI_Mtrx[FTI_Idgr*FTI_Grsz+i] == 1) {
                if (init == 0) {
                    memcpy(coding, &(data[offset*FTI_Bksz]), FTI_Bksz);
                    init = 1;
                } else {
                    enc_region_xor(&(data[offset*FTI_Bksz]), coding, coding, FTI_Bksz);
                }
            }
            // Now do the data that needs to be multiplied by a factor
            if (FTI_Mtrx[FTI_Idgr*FTI_Grsz+i] != 0 && FTI_Mtrx[FTI_Idgr*FTI_Grsz+i] != 1) {
                switch (FTI_Wdsz) {
                    case 8:  enc_w08_region_multiply(&(data[offset*FTI_Bksz]), FTI_Mtrx[FTI_Idgr*FTI_Grsz+i], FTI_Bksz, coding, init); break;
                    case 16: enc_w16_region_multiply(&(data[offset*FTI_Bksz]), FTI_Mtrx[FTI_Idgr*FTI_Grsz+i], FTI_Bksz, coding, init); break;
                }
                init = 1;
            }

            offset = 1 - offset;
            i = (i+1)%FTI_Grsz;
            cnt++;
        }

        // Writting encoded checkpoints
        fwrite(coding, sizeof(char), FTI_Bksz, efd);
        pos = pos + FTI_Bksz;
    }

    // Closing and resizing files
    fclose(fd);
    if (truncate(fn,fs) == -1)
        FTI_Prerror("Error with truncate on ckpt file");
    fclose(efd);

    // Freeing memory
    free(data);
    free(coding);
    free(myData);

    // Printing preformance and finalizing
    MPI_Barrier(FTI_COMM_WORLD);
    if (FTI_Rank == 0) {
        printf("[FTI] : Encoding time : %f\n", MPI_Wtime() - tim);
    }
}


