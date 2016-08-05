/*
 * =====================================================================================
 *
 *       Filename:  dec.c
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
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "mpi.h"
#include "fti.h"
#include "galois.h"

void region_xor(char *r1, char *r2, char *r3, int nbytes)
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

void w08_region_multiply(char *region, int multby, int nbytes, char *r2, int add)
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

void w16_region_multiply(char *region, int multby, int nbytes, char *r2, int add)
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


void matrix_dotprod(int k, int w, int *matrix_row, int *src_ids, int dest_id,
                          char **data_ptrs, char **coding_ptrs, int size)
{
  int init;
  char *dptr, *sptr;
  int i;

  if (w != 1 && w != 8 && w != 16) {
    FTI_Prerror("matrix_dotprod() called and w is not 1, 8, 16 or 32\n");
  }
  init = 0;
  dptr = (dest_id < k) ? data_ptrs[dest_id] : coding_ptrs[dest_id-k];

  /* First copy or xor any data that does not need to be multiplied by a factor */
  for (i = 0; i < k; i++) {
    if (matrix_row[i] == 1) {
      if (src_ids == NULL) {
        sptr = data_ptrs[i];
      } else if (src_ids[i] < k) {
        sptr = data_ptrs[src_ids[i]];
      } else {
        sptr = coding_ptrs[src_ids[i]-k];
      }
      if (init == 0) {
        memcpy(dptr, sptr, size);
        init = 1;
      } else {
        region_xor(sptr, dptr, dptr, size);
      }
    }
  }

  /* Now do the data that needs to be multiplied by a factor */
  for (i = 0; i < k; i++) {
    if (matrix_row[i] != 0 && matrix_row[i] != 1) {
      if (src_ids == NULL) {
        sptr = data_ptrs[i];
      } else if (src_ids[i] < k) {
        sptr = data_ptrs[src_ids[i]];
      } else {
        sptr = coding_ptrs[src_ids[i]-k];
      }
      switch (w) {
        case 8:  w08_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
        case 16: w16_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
      }
      init = 1;
    }
  }
}

void FTI_Decode(char *filename, int fs, int maxFs) {
    int i, k, m, *erased, ps, j, l;
    int *tmpids, *decMatrix, *dm_ids, *tmpmat, pos;
    char **coding, **data, *dataTmp, fn[25], efn[25], mfn[25];
    FILE *mfd, *fd, *efd;
    double tim;

    tim = MPI_Wtime();
    k = FTI_Grsz;
    m = k;

    // Calculating padding size
    ps = ((maxFs/FTI_Bksz))*FTI_Bksz;
    if (ps < maxFs)
        ps = ps + FTI_Bksz;

    sprintf(fn,"%s/%s",FTI_Cdir, filename);
    sprintf(efn,"%s/%s.enc",FTI_Cdir, filename);

    // Initializing variables
    data = talloc(char *, k);
    coding = talloc(char *, m);
    dm_ids = talloc(int, k);
    decMatrix = talloc(int, k*k);
    tmpmat = talloc(int, k*k);
    erased = talloc(int, (k+m));
    dataTmp = talloc(char, FTI_Bksz*k);
    for (i = 0; i < m; i++) {
        coding[i] = talloc(char, FTI_Bksz);
        data[i] = talloc(char, FTI_Bksz);
        erased[i] = 0;
        erased[i+k] = 0;
    }

    // Reading erasures locally
    fd = fopen(fn, "rb");
    if(fd == NULL)
        erased[FTI_Idgr] = 1;
    else
        fclose(fd);
    fd = fopen(efn, "rb");
    if(fd == NULL)
        erased[FTI_Idgr+k] = 1;
    else
        fclose(fd);

    // Distributing erasures globally
    MPI_Allgather(erased+FTI_Idgr, 1, MPI_INT, erased, 1, MPI_INT, FTI_Cenc);
    MPI_Allgather(erased+FTI_Idgr+k, 1, MPI_INT, erased+k, 1, MPI_INT, FTI_Cenc);

    // Counting the number of erasures
    l = 0;
    for(i = 0; i < k+m; i++)
        if(erased[i])
            l++;

    if (l > m)
        FTI_Prerror("Too much erasures");

    if (l > 0) {

        if (FTI_Rank == 0) {
            printf("[FTI] : There are %d (encoded)-checkpoint files missing\n",l);
            printf("[FTI] : Regenerating the lost data...\n");
        }

        // Getting ids for the decoding matrix
        j = 0;
        for (i = 0; j < k; i++) {
            if (erased[i] == 0) {
                dm_ids[j] = i;
                j++;
            }
        }

        // Building the matrix
        for (i = 0; i < k; i++) {
            if (dm_ids[i] < k) {
                for (j = 0; j < k; j++) tmpmat[i*k+j] = 0;
                tmpmat[i*k+dm_ids[i]] = 1;
            } else {
                for (j = 0; j < k; j++) {
                    tmpmat[i*k+j] = FTI_Mtrx[(dm_ids[i]-k)*k+j];
                }
            }
        }

        // Inversing the matrix
        if (FTI_InvertMatrix(tmpmat, decMatrix, k, FTI_Wdsz) < 0) {
            free(tmpmat);
            free(erased);
            free(dm_ids);
            free(decMatrix);
            FTI_Prerror("Error inversing matrix");
        }

        // Openning files
        if(erased[FTI_Idgr] == 0) {
            if (truncate(fn,ps) == -1)
                FTI_Prerror("Error with truncate on checkpoint file");
            fd = fopen(fn, "rb");
            efd = fopen(efn, "rb");
        } else {
            fd = fopen(fn, "wb");
            efd = fopen(efn, "wb");
        }
        if (fd == NULL || efd == NULL)
            FTI_Prerror("Impossible to open checkpoint file or encoded checkpoint file");

        // Main loop, block by block
        pos = 0;
        while(pos < ps) {

            // Reading the data
            if(erased[FTI_Idgr] == 0) {
                fread(data[FTI_Idgr]+0, sizeof(char), FTI_Bksz, fd);
                fread(coding[FTI_Idgr]+0, sizeof(char), FTI_Bksz, efd);
            } else {
                // Erasure found
                bzero(data[FTI_Idgr], FTI_Bksz);
                bzero(coding[FTI_Idgr], FTI_Bksz);
            }

            // Gathering the data in each process
            MPI_Allgather(data[FTI_Idgr]+0, FTI_Bksz, MPI_CHAR, dataTmp, FTI_Bksz, MPI_CHAR, FTI_Cenc);
            for (i = 0; i < k; i++) {
                memcpy(data[i]+0, &(dataTmp[i*FTI_Bksz]), sizeof(char)*FTI_Bksz);
            }

            // Gathering the encoded data in each process
            MPI_Allgather(coding[FTI_Idgr]+0, FTI_Bksz, MPI_CHAR, dataTmp, FTI_Bksz, MPI_CHAR, FTI_Cenc);
            for (i = 0; i < k; i++) {
                memcpy(coding[i]+0, &(dataTmp[i*FTI_Bksz]), sizeof(char)*FTI_Bksz);
            }

            // Decoding the lost data work
            if (erased[FTI_Idgr]) {
                matrix_dotprod(k, FTI_Wdsz, decMatrix+(FTI_Idgr*k), dm_ids, FTI_Idgr, data, coding, FTI_Bksz);
                l--;
            }

            // Gathering the recovered data in each process
            MPI_Allgather(data[FTI_Idgr]+0, FTI_Bksz, MPI_CHAR, dataTmp, FTI_Bksz, MPI_CHAR, FTI_Cenc);
            for (i = 0; i < k; i++) {
                memcpy(data[i]+0, &(dataTmp[i*FTI_Bksz]), sizeof(char)*FTI_Bksz);
            }

            // Finally, re-encode any erased encoded checkpoint file
            if (erased[FTI_Idgr + k]) {
                matrix_dotprod(k, FTI_Wdsz, FTI_Mtrx+(FTI_Idgr*k), NULL, FTI_Idgr+k, data, coding, FTI_Bksz);
            }

            // Writing restored checkpoints
            if (erased[FTI_Idgr])
                fwrite(data[FTI_Idgr]+0, sizeof(char), FTI_Bksz, fd);

            // Writing restored encoded checkpoints
            if (erased[FTI_Idgr + k])
                fwrite(coding[FTI_Idgr]+0, sizeof(char), FTI_Bksz, efd);

            pos = pos + FTI_Bksz;
        }

        // Closing files
        fclose(fd);
        fclose(efd);
        if (truncate(fn,fs) == -1)
            FTI_Prerror("Error with truncate on ckpt file");
    }

    // Freeing memory
    free(tmpmat);
    free(erased);
    free(dm_ids);
    free(decMatrix);
    free(data);
    free(dataTmp);
    free(coding);

    // Printing decoding time and finalizing
    MPI_Barrier(FTI_COMM_WORLD);
    if (FTI_Rank == 0) {
        printf("[FTI] : Decoding time : %f\n", MPI_Wtime() - tim);
    }
}

