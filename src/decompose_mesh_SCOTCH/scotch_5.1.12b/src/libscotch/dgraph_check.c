/* Copyright 2007-2010 ENSEIRB, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : dgraph_check.c                          **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Francois CHATENET (P0.0)                **/
/**                Sebastien FOUCAULT (P0.0)               **/
/**                Nicolas GICQUEL (P0.1)                  **/
/**                Jerome LACOSTE (P0.1)                   **/
/**                Cedric CHEVALIER                        **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel static mapper.       **/
/**                This module contains the distributed    **/
/**                graph consistency checking routine.     **/
/**                                                        **/
/**                # Version P0.0 : from : 01 apr 1997     **/
/**                                 to     20 jun 1997     **/
/**                # Version P0.1 : from : 14 apr 1998     **/
/**                                 to     20 jun 1998     **/
/**                # Version P0.2 : from : 11 may 1999     **/
/**                                 to     02 feb 2000     **/
/**                # Version 5.0  : from : 04 jul 2005     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 20 nov 2008     **/
/**                                 to   : 30 jul 2010     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH_CHECK

#include "module.h"
#include "common.h"
#include "dgraph.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This function checks the consistency
** of the given distributed graph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
dgraphCheck (
const Dgraph * restrict const grafptr)
{
  MPI_Comm            proccomm;                   /* Graph communicator                           */
  int                 procglbnbr;                 /* Number of processes sharing graph data       */
  int                 proclocnum;                 /* Number of this process                       */
  int                 procrcvnum;                 /* Number of process from which to receive      */
  int                 procsndnum;                 /* Number of process to which to send           */
  int                 procngbsel;                 /* Value of the currently used neighbor buffers */
  int                 procngbnum;                 /* Number of current neighbor process           */
  Gnum *              procngbtab;                 /* Array of neighbor vertex ranges              */
  Gnum                vertlocnum;
  Gnum                vertngbmin;                 /* Smallest vertex number of neighbor process   */
  Gnum                vertngbmax;                 /* Largest vertex number of neighbor process    */
  int                 vertngbnbr[2];              /* Size of the neighbor vertex arrays           */
  Gnum * restrict     vertngbtab[2];              /* Array of two neighbor vertex arrays          */
  Gnum * restrict     vendngbtab[2];              /* Array of two neighbor end vertex arrays      */
  Gnum                edgelocnbr;                 /* Local number of edges                        */
  int                 edgengbnbr[2];              /* Size of the neighbor vertex arrays           */
  Gnum * restrict     edgengbtab[2];              /* Array of two neighbor edge arrays            */
  Gnum * restrict     edlongbtab[2];              /* Array of two neighbor edge load arrays       */
  Gnum                edlolocsiz;                 /* Size of neighbor edge load array (if any)    */
  int                 cheklocval;                 /* Local consistency flag                       */
  int                 chekglbval;                 /* Global consistency flag                      */
  Gnum                reduloctab[20];             /* Arrays for reductions                        */
  Gnum                reduglbtab[20];
  MPI_Request         requloctab[8];              /* Arrays for pipelined communications          */
  MPI_Status          statloctab[8];

  proccomm = grafptr->proccomm;                   /* Simplify */

  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {    /* Synchronize */
    errorPrint ("dgraphCheck: communication error (1)");
    return     (1);
  }

  cheklocval =                                    /* Assume everything is all right */
  chekglbval = 0;
  MPI_Comm_size (proccomm, &procglbnbr);          /* Get communicator data */
  MPI_Comm_rank (proccomm, &proclocnum);

  if ((grafptr->procglbnbr != procglbnbr) ||
      (grafptr->proclocnum != proclocnum) ||
      ((grafptr->procdsptab == NULL) && (grafptr->vertlocnbr != 0))) {
    errorPrint ("dgraphCheck: inconsistent communication data (1)");
    cheklocval = 1;
  }

  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCheck: communication error (2)");
    return     (1);
  }
  if (chekglbval != 0)
    return (1);
  reduloctab[0] = (grafptr->procdsptab == NULL) ? 0 : 1; /* If private data not initialized */
  if (MPI_Allreduce (reduloctab, reduglbtab, 1, GNUM_MPI, MPI_SUM, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCheck: communication error (3)");
    return     (1);
  }
  if (reduglbtab[0] == 0)                         /* If distributed graph is empty         */
    return (0);                                   /* Do not go any further                 */
  if (reduglbtab[0] != procglbnbr) {              /* If private data not consistently here */
    errorPrint ("dgraphCheck: inconsistent communication data (2)");
    return     (1);
  }

  for (procrcvnum = 0; procrcvnum < grafptr->procglbnbr; procrcvnum ++) {
    if ((grafptr->proccnttab[procrcvnum] < 0)                                                                        ||
        (grafptr->proccnttab[procrcvnum] != (grafptr->procdsptab[procrcvnum + 1] - grafptr->procdsptab[procrcvnum])) ||
        (grafptr->proccnttab[procrcvnum] >  (grafptr->procvrttab[procrcvnum + 1] - grafptr->procvrttab[procrcvnum]))) {
      errorPrint ("dgraphCheck: inconsistent communication data (3)");
      cheklocval = 1;
    }
  }
  if (grafptr->proccnttab[proclocnum] != grafptr->vertlocnbr) {
    errorPrint ("dgraphCheck: inconsistent communication data (4)");
    cheklocval = 1;
  }

  procrcvnum = (proclocnum + 1) % procglbnbr;     /* Compute indices of neighbors */
  procsndnum = (proclocnum - 1 + procglbnbr) % procglbnbr;

  if ((procngbtab = (Gnum *) memAlloc ((procglbnbr + 1) * sizeof (Gnum))) == NULL) {
    errorPrint ("dgraphCheck: out of memory (1)");
    cheklocval = 1;
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCheck: communication error (4)");
    return     (1);
  }
  if (chekglbval != 0) {
    if (procngbtab != NULL)
      memFree (procngbtab);
    return (1);
  }

  MPI_Sendrecv (grafptr->procdsptab, procglbnbr + 1, GNUM_MPI, procsndnum, TAGPROCVRTTAB, /* Check vertex range array */
                procngbtab,          procglbnbr + 1, GNUM_MPI, procrcvnum, TAGPROCVRTTAB,
                proccomm, &statloctab[0]);
  for (procngbnum = 0; procngbnum <= procglbnbr; procngbnum ++) {
    if (grafptr->procdsptab[procngbnum] != procngbtab[procngbnum]) {
      errorPrint ("dgraphCheck: inconsistent communication data (5)");
      cheklocval = 1;
      break;
    }
  }

  MPI_Sendrecv (grafptr->procvrttab, procglbnbr + 1, GNUM_MPI, procsndnum, TAGPROCVRTTAB, /* Check vertex range array */
                procngbtab,          procglbnbr + 1, GNUM_MPI, procrcvnum, TAGPROCVRTTAB,
                proccomm, &statloctab[0]);
  for (procngbnum = 0; procngbnum <= procglbnbr; procngbnum ++) {
    if (grafptr->procvrttab[procngbnum] != procngbtab[procngbnum]) {
      errorPrint ("dgraphCheck: inconsistent communication data (6)");
      cheklocval = 1;
      break;
    }
  }
  memFree (procngbtab);

  if ((grafptr->baseval < 0) ||                   /* Elementary constraints on graph fields */
      (grafptr->baseval > 1) ||                   /* Strong limitation on base value        */
      (grafptr->vertlocnbr < 0) ||
      (grafptr->vertlocnnd != (grafptr->vertlocnbr + grafptr->baseval)) ||
      (((grafptr->flagval & DGRAPHHASEDGEGST) != 0) &&
       ((grafptr->vertgstnbr < grafptr->vertlocnbr) ||
        (grafptr->vertgstnnd != (grafptr->vertgstnbr + grafptr->baseval)))) ||
      (grafptr->edgelocnbr < 0) ||
      (grafptr->edgelocsiz < grafptr->edgelocnbr)) {
    errorPrint ("dgraphCheck: inconsistent local graph data");
    cheklocval = 1;
  }

  reduloctab[ 0] =   grafptr->flagval;
  reduloctab[ 1] = - grafptr->flagval;
  reduloctab[ 2] =   grafptr->baseval;
  reduloctab[ 3] = - grafptr->baseval;
  reduloctab[ 4] =   grafptr->vertglbnbr;
  reduloctab[ 5] = - grafptr->vertglbnbr;
  reduloctab[ 6] =   grafptr->vertglbmax;
  reduloctab[ 7] = - grafptr->vertglbmax;
  reduloctab[ 8] =   grafptr->vertlocnbr;
  reduloctab[ 9] =   grafptr->edgeglbnbr;
  reduloctab[10] = - grafptr->edgeglbnbr;
  reduloctab[11] =   grafptr->edgeglbmax;
  reduloctab[12] = - grafptr->edgeglbmax;
  reduloctab[13] =   grafptr->edgelocnbr;
  reduloctab[14] =   grafptr->edgelocsiz;
  reduloctab[15] =   grafptr->edgeglbsmx;
  reduloctab[16] = - grafptr->edgeglbsmx;
  reduloctab[17] =   grafptr->degrglbmax;
  reduloctab[18] = - grafptr->degrglbmax;
  reduloctab[19] = (Gnum) cheklocval;
  if (MPI_Allreduce (reduloctab, reduglbtab, 20, GNUM_MPI, MPI_MAX, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCheck: communication error (5)");
    return     (1);
  }
  if (reduglbtab[19] != 0)
    return (1);
  if ((reduglbtab[ 1] != - reduglbtab[ 0]) ||       /* Check if global graph data match */
      (reduglbtab[ 3] != - reduglbtab[ 2]) ||
      (reduglbtab[ 5] != - reduglbtab[ 4]) ||
      (reduglbtab[ 7] != - reduloctab[ 6]) ||
      (reduglbtab[ 8] !=   reduloctab[ 6]) ||
      (reduglbtab[10] != - reduglbtab[ 9]) ||
      (reduglbtab[12] != - reduglbtab[11]) ||
      (reduglbtab[13] !=   reduloctab[11]) ||     /* Recompute and test maximum number of local edges    */
      (reduglbtab[14] !=   reduloctab[15]) ||     /* Recompute and test maximum size of local edge array */
      (reduglbtab[16] != - reduloctab[15]) ||
      (reduglbtab[18] != - reduloctab[17])) {
    errorPrint ("dgraphCheck: inconsistent global graph data (1)");
    cheklocval = 1;
  }
  reduloctab[0] = (grafptr->veloloctax != NULL) ? 1 : 0; /* Check consistency */
  reduloctab[1] = (grafptr->edgegsttax != NULL) ? 1 : 0;
  reduloctab[2] = (grafptr->edloloctax != NULL) ? 1 : 0;
  reduloctab[3] = (grafptr->vnumloctax != NULL) ? 1 : 0;
  reduloctab[4] = grafptr->vertlocnbr;            /* Recompute local sizes */
  reduloctab[5] = grafptr->edgelocnbr;
  reduloctab[6] = (Gnum) cheklocval;
  if (MPI_Allreduce (reduloctab, reduglbtab, 7, GNUM_MPI, MPI_SUM, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCheck: communication error (6)");
    return     (1);
  }
  if (reduglbtab[6] != 0)
    return (1);
  if (((reduglbtab[0] != 0) && (reduglbtab[0] != procglbnbr)) ||
      ((reduglbtab[1] != 0) && (reduglbtab[1] != procglbnbr)) ||
      ((reduglbtab[2] != 0) && (reduglbtab[2] != procglbnbr)) ||
      ((reduglbtab[3] != 0) && (reduglbtab[3] != procglbnbr)) ||
      (reduglbtab[4] != grafptr->vertglbnbr)                  ||
      (reduglbtab[5] != grafptr->edgeglbnbr)) {
    errorPrint ("dgraphCheck: inconsistent global graph data (2)");
    cheklocval = 1;
  }

  for (vertlocnum = grafptr->baseval, edgelocnbr = 0; vertlocnum < grafptr->vertlocnnd; vertlocnum ++) {
    Gnum                edgelocnum;

    if ((grafptr->vendloctax[vertlocnum] < grafptr->vertloctax[vertlocnum]) ||
        (grafptr->vendloctax[vertlocnum] > (grafptr->edgelocsiz + grafptr->baseval))) {
      errorPrint ("dgraphCheck: inconsistent local vertex arrays");
      edgelocnbr = grafptr->edgelocnbr;           /* Avoid unwanted cascaded error messages */
      cheklocval = 1;
      break;
    }
    edgelocnbr += grafptr->vendloctax[vertlocnum] - grafptr->vertloctax[vertlocnum];

    if ((grafptr->flagval & DGRAPHHASEDGEGST) != 0) { /* If ghost edge array is valid */
      for (edgelocnum = grafptr->vertloctax[vertlocnum]; edgelocnum < grafptr->vendloctax[vertlocnum]; edgelocnum ++) {
        if ((grafptr->edgegsttax[edgelocnum] < grafptr->baseval) ||
            (grafptr->edgegsttax[edgelocnum] >= grafptr->vertgstnnd)) {
          errorPrint ("dgraphCheck: inconsistent ghost edge array");
          edgelocnbr = grafptr->edgelocnbr;       /* Avoid unwanted cascaded error messages */
          cheklocval = 1;
          break;
        }
      }
    }
  }
  if (edgelocnbr != grafptr->edgelocnbr) {
    errorPrint ("dgraphCheck: invalid local number of edges");
    cheklocval = 1;
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCheck: communication error (7)");
    return     (1);
  }
  if (chekglbval != 0)
    return (1);

  if (grafptr->veloloctax != NULL) { /* We must check that load data are consistent */
    Gnum                velolocsum;
    Gnum                veloglbsum;

    for (vertlocnum = grafptr->baseval, velolocsum = 0; vertlocnum < grafptr->vertlocnnd; vertlocnum ++)
      velolocsum += grafptr->veloloctax[vertlocnum];

    MPI_Allreduce (&velolocsum, &veloglbsum, 1, GNUM_MPI, MPI_SUM, proccomm);

    cheklocval = 0;
    if (velolocsum != grafptr->velolocsum) {
      errorPrint ("dgraphCheck: invalid local vertex load sum");
      cheklocval = 1;
    }
    if (veloglbsum != grafptr->veloglbsum) {
      errorPrint ("dgraphCheck: invalid global vertex load sum");
      cheklocval = 1;
    }
    MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, proccomm);
    if (chekglbval != 0)
      return (1);
  }

  edlolocsiz = (grafptr->edloloctax != NULL) ? grafptr->edgeglbsmx : 0;
  if (memAllocGroup ((void **) (void *)
                     &vertngbtab[0], (size_t) (grafptr->vertglbmax * sizeof (Gnum)), /* Send vertex and vertex end arrays, even when they are compact */
                     &vertngbtab[1], (size_t) (grafptr->vertglbmax * sizeof (Gnum)),
                     &vendngbtab[0], (size_t) (grafptr->vertglbmax * sizeof (Gnum)),
                     &vendngbtab[1], (size_t) (grafptr->vertglbmax * sizeof (Gnum)),
                     &edgengbtab[0], (size_t) (grafptr->edgeglbsmx * sizeof (Gnum)),
                     &edgengbtab[1], (size_t) (grafptr->edgeglbsmx * sizeof (Gnum)),
                     &edlongbtab[0], (size_t) (edlolocsiz          * sizeof (Gnum)),
                     &edlongbtab[1], (size_t) (edlolocsiz          * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphCheck: out of memory (2)");
    cheklocval = 1;
  }
  if (grafptr->edloloctax == NULL) {              /* If graph edges are not weighted */
    edlongbtab[0] =                               /* Edge load arrays are fake       */
    edlongbtab[1] = NULL;
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCheck: communication error (8)");
    return     (1);
  }
  if (chekglbval != 0) {
    if (cheklocval == 0)
      memFree (vertngbtab[0]);                    /* Free group leader */
    return (1);
  }

  MPI_Irecv (vertngbtab[0], grafptr->vertglbmax, GNUM_MPI, procrcvnum, TAGVERTLOCTAB, proccomm, &requloctab[0]);
  MPI_Irecv (vendngbtab[0], grafptr->vertglbmax, GNUM_MPI, procrcvnum, TAGVENDLOCTAB, proccomm, &requloctab[1]);
  MPI_Irecv (edgengbtab[0], grafptr->edgeglbsmx, GNUM_MPI, procrcvnum, TAGEDGELOCTAB, proccomm, &requloctab[2]);
  if (grafptr->edloloctax != NULL)
    MPI_Irecv (edlongbtab[0], grafptr->edgeglbsmx, GNUM_MPI, procrcvnum, TAGEDLOLOCTAB, proccomm, &requloctab[3]);

  MPI_Send (grafptr->vertloctax + grafptr->baseval, grafptr->vertlocnbr, GNUM_MPI, procsndnum, TAGVERTLOCTAB, proccomm);
  MPI_Send (grafptr->vendloctax + grafptr->baseval, grafptr->vertlocnbr, GNUM_MPI, procsndnum, TAGVENDLOCTAB, proccomm);
  MPI_Send (grafptr->edgeloctax + grafptr->baseval, grafptr->edgelocsiz, GNUM_MPI, procsndnum, TAGEDGELOCTAB, proccomm);
  if (grafptr->edloloctax != NULL) {              /* Send synchronously vertloctab and vendloctab to avoid MPI Isend problem if same array */
    MPI_Send (grafptr->edloloctax + grafptr->baseval, grafptr->edgelocsiz, GNUM_MPI, procsndnum, TAGEDLOLOCTAB, proccomm);
    MPI_Waitall (4, &requloctab[0], &statloctab[0]);
  }
  else
    MPI_Waitall (3, &requloctab[0], &statloctab[0]);

  MPI_Get_count (&statloctab[0], GNUM_MPI, &vertngbnbr[0]);
  MPI_Get_count (&statloctab[2], GNUM_MPI, &edgengbnbr[0]);

  for (procngbnum  = (proclocnum + 1) % procglbnbr, procngbsel = 0; /* Loop on all other processes */
       procngbnum != proclocnum;
       procngbnum  = (procngbnum + 1) % procglbnbr, procngbsel ^= 1) {
    Gnum                vertlocnum;
    Gnum                vertglbnum;
    Gnum * restrict     edgengbtax;
    Gnum * restrict     edlongbtax;

    vertngbmin = grafptr->procvrttab[procngbnum]; /* Get neighbor vertex number range */
    vertngbmax = grafptr->procvrttab[procngbnum + 1];

#ifdef SCOTCH_DEBUG_DGRAPH2
    if ((vertngbnbr[procngbsel] != grafptr->proccnttab[procngbnum]) ||
        (vertngbnbr[procngbsel]  > (vertngbmax - vertngbmin))) {
      errorPrint ("dgraphCheck: internal error (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

    if (((procngbnum + 1) % procglbnbr) != proclocnum) {
      MPI_Irecv (vertngbtab[1 - procngbsel], grafptr->vertglbmax, GNUM_MPI, procrcvnum, TAGVERTLOCTAB, proccomm, &requloctab[4]);
      MPI_Irecv (vendngbtab[1 - procngbsel], grafptr->vertglbmax, GNUM_MPI, procrcvnum, TAGVENDLOCTAB, proccomm, &requloctab[5]);
      MPI_Irecv (edgengbtab[1 - procngbsel], grafptr->edgeglbsmx, GNUM_MPI, procrcvnum, TAGEDGELOCTAB, proccomm, &requloctab[6]);
      MPI_Isend (vertngbtab[procngbsel], vertngbnbr[procngbsel], GNUM_MPI, procsndnum, TAGVERTLOCTAB, proccomm, &requloctab[1]);
      MPI_Isend (vendngbtab[procngbsel], vertngbnbr[procngbsel], GNUM_MPI, procsndnum, TAGVENDLOCTAB, proccomm, &requloctab[2]);
      MPI_Isend (edgengbtab[procngbsel], edgengbnbr[procngbsel], GNUM_MPI, procsndnum, TAGEDGELOCTAB, proccomm, &requloctab[3]);
      if (grafptr->edloloctax != NULL) {
        MPI_Irecv (edlongbtab[1 - procngbsel], grafptr->edgeglbsmx, GNUM_MPI, procrcvnum, TAGEDLOLOCTAB, proccomm, &requloctab[7]);
        MPI_Isend (edlongbtab[procngbsel], edgengbnbr[procngbsel], GNUM_MPI, procsndnum, TAGEDLOLOCTAB, proccomm, &requloctab[0]);
      }
    }
    edgengbtax = edgengbtab[procngbsel] - grafptr->baseval;
    edlongbtax = (grafptr->edloloctax != NULL) ? (edlongbtab[procngbsel] - grafptr->baseval) : NULL;

    if (((procngbnum + 1) % procglbnbr) != proclocnum) { /* Before MPI 2.2, complete communications before accessing arrays being sent */
      if (grafptr->edloloctax != NULL)
        MPI_Waitall (8, &requloctab[0], &statloctab[0]);
      else
        MPI_Waitall (6, &requloctab[1], &statloctab[1]);
      MPI_Get_count (&statloctab[4], GNUM_MPI, &vertngbnbr[1 - procngbsel]);
      MPI_Get_count (&statloctab[6], GNUM_MPI, &edgengbnbr[1 - procngbsel]);
    }

    for (vertlocnum = grafptr->baseval, vertglbnum = grafptr->procvrttab[proclocnum];
         vertlocnum < grafptr->vertlocnnd; vertlocnum ++, vertglbnum ++) {
      Gnum                edgelocnum;

      for (edgelocnum = grafptr->vertloctax[vertlocnum];
           edgelocnum < grafptr->vendloctax[vertlocnum]; edgelocnum ++) {
        Gnum                vertglbend;

        vertglbend = grafptr->edgeloctax[edgelocnum];
        if ((vertglbend >= vertngbmin) &&         /* If end vertex belongs to current neighbor process */
            (vertglbend <  vertngbmax)) {
          Gnum                edgengbnum;
          Gnum                edgengbnnd;
          Gnum                edgengbcnt;

          for (edgengbnum = vertngbtab[procngbsel][vertglbend - vertngbmin],
               edgengbnnd = vendngbtab[procngbsel][vertglbend - vertngbmin], edgengbcnt = 0;
               edgengbnum < edgengbnnd; edgengbnum ++) {
            if (edgengbtax[edgengbnum] == vertglbnum) { /* If matching edge found */
              edgengbcnt ++;                      /* Account for it               */
              if ((edlongbtax != NULL) &&         /* If edge weights do not match */
                  (edlongbtax[edgengbnum] != grafptr->edloloctax[edgelocnum]))
                cheklocval = 3;
            }
          }
          if (edgengbcnt < 1)                     /* If matching edge not found */
            cheklocval = 1;
          else if (edgengbcnt > 1)                /* If duplicate edge */
            cheklocval = 2;
        }
      }
    }
    MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, proccomm);

    if (chekglbval != 0) {                        /* Error number ranges from 1 to 3 */
      if (chekglbval == 1)
        errorPrint ("dgraphCheck: arc data do not match");
      else if (chekglbval == 2)
        errorPrint ("dgraphCheck: duplicate arc");
      else
        errorPrint ("dgraphCheck: arc load data do not match");
      memFree (vertngbtab[0]);                    /* Free group leader */
      return  (1);
    }
  }
  memFree (vertngbtab[0]);                        /* Free group leader */

  return (0);
}
