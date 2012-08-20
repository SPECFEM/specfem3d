/* Copyright 2007-2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_band.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a distributed band **/
/**                graph from the given frontier array.    **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 11 nov 2007     **/
/**                                 to   : 20 feb 2011     **/
/**                                                        **/
/**   NOTES      : # This code derives from the code of    **/
/**                  vdgraph_separate_bd.c in version 5.0. **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DGRAPH_BAND

#include "module.h"
#include "common.h"
#include "dgraph.h"

/**********************************/
/*                                */
/* Distance computation routines. */
/*                                */
/**********************************/

/* This routine computes a distributed index array
** of given width around the current separator,
** using collective communication.
** It returns:
** - 0   : if the index array could be computed.
** - !0  : on error.
*/

int
dgraphBandColl (
Dgraph * restrict const           grafptr,        /*+ Distributed graph                                        +*/
const Gnum                        queulocnbr,     /*+ Number of frontier vertices, start size for vertex queue +*/
Gnum * restrict const             queuloctab,     /*+ Array of frontier vertices, re-used as queue array       +*/
const Gnum                        distmax,        /*+ Maximum distance from separator vertices                 +*/
Gnum * restrict * restrict const  vnumgstptr,     /*+ Pointer to vnumgsttax                                    +*/
Gnum * restrict const             bandvertlvlptr, /*+ Pointer to based start index of last level               +*/
Gnum * restrict const             bandvertlocptr, /*+ Pointer to bandvertlocnnd                                +*/
Gnum * restrict const             bandedgelocptr) /*+ Pointer to bandedgelocnbr                                +*/
{
  Gnum                    queulocnum;
  Gnum                    vertlocnnd;
  Gnum * restrict         vnumgsttax;             /* Index array for vertices kept in band graph */
  Gnum                    vnumgstsiz;             /* Size of vnumgsttax; TRICK: re-use           */
  Gnum *                  vrcvdattab;             /* Data arrays [norestrict:async]              */
  Gnum * restrict         vsnddattab;
  Gnum * restrict         procvgbtab;             /* Array of neighbor bounds [+1]               */
  int                     procngbnbr;
  int                     procngbnum;
  int * restrict          vrcvcnttab;
  int * restrict          vsndcnttab;
  int * restrict          vrcvdsptab;
  int * restrict          vsnddsptab;
  int * restrict          nsndidxtab;
  int                     vrcvdspnum;
  int                     vsnddspnum;
  Gnum                    queuheadidx;            /* Index of head of queue                      */
  Gnum                    queutailidx;            /* Index of tail of queue                      */
  Gnum                    bandvertlvlnum;
  Gnum                    bandvertlocnnd;
  Gnum                    bandedgelocnbr;
  Gnum                    distval;
#ifdef SCOTCH_DEBUG_DGRAPH1
  Gnum                    reduloctab[3];
  Gnum                    reduglbtab[3];
#else /* SCOTCH_DEBUG_DGRAPH1 */
  Gnum                    reduglbtab[1];
#endif /* SCOTCH_DEBUG_DGRAPH1 */

  const Gnum * restrict const vertloctax = grafptr->vertloctax;
  const Gnum * restrict const vendloctax = grafptr->vendloctax;
  const Gnum * restrict const edgegsttax = grafptr->edgegsttax;
  const Gnum * restrict const edgeloctax = grafptr->edgeloctax;

  procngbnbr = grafptr->procngbnbr;

  reduglbtab[0] = 0;                             /* Assume everything is all right */
  vnumgstsiz    = MAX ((grafptr->vertgstnbr * sizeof (Gnum)), (grafptr->procglbnbr * sizeof (int))); /* TRICK: re-use array for further error collective communications */
  if (((vnumgsttax = memAlloc (vnumgstsiz)) == NULL) ||
      (memAllocGroup ((void **) (void *)
                      &procvgbtab, (size_t) ((procngbnbr + 1) * sizeof (Gnum)),
                      &nsndidxtab, (size_t) (procngbnbr * sizeof (int)),
                      &vrcvcnttab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                      &vsndcnttab, (size_t) (grafptr->procglbnbr * sizeof (int)), /* TRICK: vsndcnttab, vrcvdsptab, vrcvdattab, vrcvdattab joined */
                      &vrcvdsptab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                      &vsnddsptab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                      &vrcvdattab, (size_t) (grafptr->procsndnbr * sizeof (Gnum)), /* Senders and receivers inverted because we send local, not halo vertices */
                      &vsnddattab, (size_t) ((grafptr->vertgstnbr - grafptr->vertlocnbr) * sizeof (Gnum)), NULL) == NULL)) {
    errorPrint ("dgraphBandColl: out of memory (1)");
    reduglbtab[0] = 1;
  }
#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication not needed and not absorbable by the algorithm */
  reduloctab[0] = reduglbtab[0];
  reduloctab[1] =   distmax;
  reduloctab[2] = - distmax;
  if (MPI_Allreduce (reduloctab, reduglbtab, 3, GNUM_MPI, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphBandColl: communication error (1)");
    return     (1);
  }
  if (reduglbtab[1] != - reduglbtab[2]) {
    errorPrint ("dgraphBandColl: invalid parameters");
    reduglbtab[0] = 1;
  }
#endif /* SCOTCH_DEBUG_DGRAPH1 */
  if (reduglbtab[0] != 0) {
    if (vnumgsttax != NULL) {
      if (procvgbtab != NULL)
        memFree (procvgbtab);                     /* Free group leader */
      memFree (vnumgsttax);
    }
    return (1);
  }

  memSet (vsndcnttab, 0, (((int *) vrcvdattab) - vsndcnttab) * sizeof (int)); /* TRICK: vsndcnttab, vrcvdsptab, vrcvdattab, vrcvdattab joined */

  for (procngbnum = 0, vrcvdspnum = vsnddspnum = 0; /* Build communication index arrays */
       procngbnum < procngbnbr; procngbnum ++) {
    int                 procglbnum;

    procglbnum = grafptr->procngbtab[procngbnum];
    procvgbtab[procngbnum] = grafptr->procvrttab[procglbnum];
    vrcvdsptab[procglbnum] = vrcvdspnum;
    vsnddsptab[procglbnum] = vsnddspnum;
    vrcvdspnum += grafptr->procsndtab[procglbnum]; /* Senders and receivers are reversed */
    vsnddspnum += grafptr->procrcvtab[procglbnum];
  }
  procvgbtab[procngbnum] = grafptr->procvrttab[grafptr->procglbnbr];

  bandvertlvlnum =                                /* Start index of last level is start index */
  bandvertlocnnd = grafptr->baseval;              /* Reset number of band vertices, plus base */
  bandedgelocnbr = 0;
  memSet (vnumgsttax, ~0, grafptr->vertgstnbr * sizeof (Gnum)); /* Reset part array */
  vnumgsttax -= grafptr->baseval;

  vertlocnnd = grafptr->vertlocnnd;
  for (queulocnum = 0; queulocnum < queulocnbr; queulocnum ++) { /* All frontier vertices will be first vertices of band graph */
    Gnum              vertlocnum;

    vertlocnum = queuloctab[queulocnum];
    vnumgsttax[vertlocnum] = bandvertlocnnd ++;   /* Keep frontier vertex in band               */
    bandedgelocnbr += vendloctax[vertlocnum] - vertloctax[vertlocnum]; /* Account for its edges */
  }
  queuheadidx = 0;                                /* No queued vertex read yet                  */
  queutailidx = queulocnbr;                       /* All frontier vertices are already in queue */

  for (distval = 0; ++ distval <= distmax; ) {
    Gnum              queunextidx;                /* Tail index for enqueuing vertices of next band */
    int               procngbnum;

    bandvertlvlnum = bandvertlocnnd;              /* Save start index of current level, based */
    *bandvertlvlptr = bandvertlvlnum;

    for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) /* Build communication index arrays */
      nsndidxtab[procngbnum] = vsnddsptab[grafptr->procngbtab[procngbnum]];

    for (queunextidx = queutailidx; queuheadidx < queutailidx; ) { /* For all vertices in queue */
      Gnum              vertlocnum;
      Gnum              edgelocnum;

      vertlocnum = queuloctab[queuheadidx ++];    /* Dequeue vertex */
      for (edgelocnum = vertloctax[vertlocnum]; edgelocnum < vendloctax[vertlocnum]; edgelocnum ++) {
        Gnum              vertlocend;

        vertlocend = edgegsttax[edgelocnum];
        if (vnumgsttax[vertlocend] != ~0)         /* If end vertex has already been processed */
          continue;                               /* Skip to next vertex                      */

        if (vertlocend < vertlocnnd) {            /* If end vertex is local                           */
          vnumgsttax[vertlocend] = bandvertlocnnd ++; /* Label vertex as enqueued                     */
          queuloctab[queunextidx ++] = vertlocend; /* Enqueue vertex for next pass                    */
          bandedgelocnbr += vendloctax[vertlocend] - vertloctax[vertlocend]; /* Account for its edges */
        }
        else {                                    /* End vertex is a ghost */
          Gnum              vertglbend;
          int               procngbnum;
          int               procngbmax;

          vnumgsttax[vertlocend] = 0;             /* Label ghost vertex as enqueued  */
          vertglbend = edgeloctax[edgelocnum];    /* Get global number of end vertex */

          procngbnum = 0;
          procngbmax = procngbnbr;
          while ((procngbmax - procngbnum) > 1) { /* Find owner process by dichotomy on procvgbtab */
            int                 procngbmed;

            procngbmed = (procngbmax + procngbnum) / 2;
            if (procvgbtab[procngbmed] > vertglbend)
              procngbmax = procngbmed;
            else
              procngbnum = procngbmed;
          }
#ifdef SCOTCH_DEBUG_DGRAPH2
          if ((grafptr->procvrttab[grafptr->procngbtab[procngbnum]]     >  vertglbend) ||
              (grafptr->procvrttab[grafptr->procngbtab[procngbnum] + 1] <= vertglbend) ||
              (nsndidxtab[procngbnum] >= (vsnddsptab[grafptr->procngbtab[procngbnum]] + grafptr->procrcvtab[grafptr->procngbtab[procngbnum]]))) {
            errorPrint ("dgraphBandColl: internal error (1)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
          vsnddattab[nsndidxtab[procngbnum] ++] = vertglbend - procvgbtab[procngbnum] + grafptr->baseval; /* Buffer local value on neighbor processor */
        }
      }
    }

    for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) {
      int                 procglbnum;

      procglbnum = grafptr->procngbtab[procngbnum];
      vsndcnttab[procglbnum] = nsndidxtab[procngbnum] - vsnddsptab[procglbnum];
    }

    if (MPI_Alltoall (vsndcnttab, 1, MPI_INT, vrcvcnttab, 1, MPI_INT, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphBandColl: communication error (2)");
      return     (1);
    }
    if (MPI_Alltoallv (vsnddattab, vsndcnttab, vsnddsptab, GNUM_MPI,
                       vrcvdattab, vrcvcnttab, vrcvdsptab, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphBandColl: communication error (3)");
      return     (1);
    }

    for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) { /* For all receive buffers */
      Gnum * restrict   vrcvdatptr;
      int               vertrcvnum;
      int               procglbnum;
      int               statsiz;

      procglbnum = grafptr->procngbtab[procngbnum];
      vrcvdatptr = vrcvdattab + vrcvdsptab[procglbnum];
      statsiz    = vrcvcnttab[procglbnum];
      for (vertrcvnum = 0; vertrcvnum < statsiz; vertrcvnum ++) {
        Gnum              vertlocend;

        vertlocend = vrcvdatptr[vertrcvnum];      /* Get local vertex from message */
#ifdef SCOTCH_DEBUG_DGRAPH2
        if ((vertlocend < grafptr->baseval) || (vertlocend >= vertlocnnd)) {
          errorPrint ("dgraphBandColl: internal error (2)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
        if (vnumgsttax[vertlocend] != ~0)         /* If end vertex has already been processed */
          continue;                               /* Skip to next vertex                      */

        vnumgsttax[vertlocend] = bandvertlocnnd ++; /* Label vertex as enqueued                     */
        queuloctab[queunextidx ++] = vertlocend;  /* Enqueue vertex for next pass                   */
        bandedgelocnbr += vendloctax[vertlocend] - vertloctax[vertlocend]; /* Account for its edges */
      }
    }

    queutailidx = queunextidx;                    /* Prepare queue for next sweep */
  }

  memFree (procvgbtab);                           /* Free group leader */

  *vnumgstptr     = vnumgsttax;
  *bandvertlocptr = bandvertlocnnd - grafptr->baseval;
  *bandedgelocptr = bandedgelocnbr;

  return (0);
}

/* This routine computes a distributed index array
** of given width around the current separator,
** using point-to-point communication.
** It returns:
** - 0   : if the index array could be computed.
** - !0  : on error.
*/

int
dgraphBandPtop (
Dgraph * restrict const           grafptr,        /*+ Distributed graph                                        +*/
const Gnum                        queulocnbr,     /*+ Number of frontier vertices, start size for vertex queue +*/
Gnum * restrict const             queuloctab,     /*+ Array of frontier vertices, re-used as queue array       +*/
const Gnum                        distmax,        /*+ Maximum distance from separator vertices                 +*/
Gnum * restrict * restrict const  vnumgstptr,     /*+ Pointer to vnumgsttax                                    +*/
Gnum * restrict const             bandvertlvlptr, /*+ Pointer to based start index of last level               +*/
Gnum * restrict const             bandvertlocptr, /*+ Pointer to bandvertlocnnd                                +*/
Gnum * restrict const             bandedgelocptr) /*+ Pointer to bandedgelocnbr                                +*/
{
  Gnum                    queulocnum;
  Gnum                    vertlocnnd;
  Gnum * restrict         vnumgsttax;             /* Index array for vertices kept in band graph */
  Gnum                    vnumgstsiz;             /* Size of vnumgsttax; TRICK: re-use           */
  Gnum *                  vrcvdattab;             /* Data arrays [norestrict:async]              */
  Gnum *                  vsnddattab;
  Gnum * restrict         procvgbtab;             /* Array of neighbor bounds [+1]               */
  int                     procngbnbr;
  int                     procngbnum;
  int                     procngbnxt;
  Gnum * restrict         nrcvdsptab;
  Gnum * restrict         nsnddsptab;
  Gnum                    nrcvdspnum;
  Gnum                    nsnddspnum;
  Gnum * restrict         nsndidxtab;
  Gnum                    queuheadidx;            /* Index of head of queue                      */
  Gnum                    queutailidx;            /* Index of tail of queue                      */
  MPI_Request *           nrcvreqtab;             /* Array of receive requests                   */
  MPI_Request *           nsndreqtab;             /* Array of receive requests                   */
  Gnum                    bandvertlvlnum;
  Gnum                    bandvertlocnnd;
  Gnum                    bandedgelocnbr;
  Gnum                    distval;
#ifdef SCOTCH_DEBUG_DGRAPH1
  Gnum                    reduloctab[3];
  Gnum                    reduglbtab[3];
#else /* SCOTCH_DEBUG_DGRAPH1 */
  Gnum                    reduglbtab[1];
#endif /* SCOTCH_DEBUG_DGRAPH1 */

  const Gnum * restrict const vertloctax = grafptr->vertloctax;
  const Gnum * restrict const vendloctax = grafptr->vendloctax;
  const Gnum * restrict const edgegsttax = grafptr->edgegsttax;
  const Gnum * restrict const edgeloctax = grafptr->edgeloctax;

  procngbnbr = grafptr->procngbnbr;

  reduglbtab[0] = 0;                             /* Assume everything is all right */
  vnumgstsiz    = MAX ((grafptr->vertgstnbr * sizeof (Gnum)), (grafptr->procglbnbr * sizeof (int))); /* TRICK: re-use array for further error collective communications */
  if (((vnumgsttax = memAlloc (vnumgstsiz)) == NULL) ||
      (memAllocGroup ((void **) (void *)
                      &procvgbtab, (size_t) ((procngbnbr + 1) * sizeof (Gnum)),
                      &nrcvdsptab, (size_t) ((procngbnbr + 1) * sizeof (Gnum)), /* +1 to check against end of array */
                      &nsnddsptab, (size_t) ((procngbnbr + 1) * sizeof (Gnum)), /* +1 to check against end of array */
                      &nsndidxtab, (size_t) (procngbnbr * sizeof (Gnum)), /* Here Gnum's since point-to-point       */
                      &nrcvreqtab, (size_t) (procngbnbr * sizeof (MPI_Request)),
                      &nsndreqtab, (size_t) (procngbnbr * sizeof (MPI_Request)),
                      &vrcvdattab, (size_t) (grafptr->procsndnbr * sizeof (Gnum)), /* Senders and receivers inverted because we send local, not halo vertices */
                      &vsnddattab, (size_t) ((grafptr->vertgstnbr - grafptr->vertlocnbr) * sizeof (Gnum)), NULL) == NULL)) {
    errorPrint ("dgraphBandPtop: out of memory (1)");
    reduglbtab[0] = 1;
  }
#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication not needed and not absorbable by the algorithm */
  reduloctab[0] = reduglbtab[0];
  reduloctab[1] =   distmax;
  reduloctab[2] = - distmax;
  if (MPI_Allreduce (reduloctab, reduglbtab, 3, GNUM_MPI, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphBandPtop: communication error (1)");
    return     (1);
  }
  if (reduglbtab[1] != - reduglbtab[2]) {
    errorPrint ("dgraphBandPtop: invalid parameters");
    reduglbtab[0] = 1;
  }
#endif /* SCOTCH_DEBUG_DGRAPH1 */
  if (reduglbtab[0] != 0) {
    if (vnumgsttax != NULL) {
      if (procvgbtab != NULL)
        memFree (procvgbtab);                     /* Free group leader */
      memFree (vnumgsttax);
    }
    return (1);
  }

  for (procngbnum = 0, nrcvdspnum = nsnddspnum = procngbnxt = 0; /* Build communication index arrays */
       procngbnum < procngbnbr; procngbnum ++) {
    int                 procglbnum;

    procglbnum = grafptr->procngbtab[procngbnum];
    if ((procngbnxt == 0) && (procglbnum > grafptr->proclocnum)) /* Find index of first neighbor of higher rank */
      procngbnxt = procngbnum;
    procvgbtab[procngbnum] = grafptr->procvrttab[procglbnum];
    nrcvdsptab[procngbnum] = nrcvdspnum;          /* Arrays are indexed per neighbor since we are doing point-to-point communication */
    nsnddsptab[procngbnum] = nsnddspnum;
    nrcvdspnum += grafptr->procsndtab[procglbnum]; /* Senders and receivers are reversed */
    nsnddspnum += grafptr->procrcvtab[procglbnum];
  }
  procvgbtab[procngbnum] = grafptr->procvrttab[grafptr->procglbnbr];
  nrcvdsptab[procngbnum] = nrcvdspnum;            /* Mark end of communication index arrays */
  nsnddsptab[procngbnum] = nsnddspnum;

  procngbnum = procngbnxt;                        /* Create receive requests in descending order */
  do {
    procngbnum = (procngbnum + (procngbnbr - 1)) % procngbnbr; /* Pre-decrement neighbor rank */
    if (MPI_Recv_init (vrcvdattab + nrcvdsptab[procngbnum], (int) (nrcvdsptab[procngbnum + 1] - nrcvdsptab[procngbnum]), GNUM_MPI,
                       grafptr->procngbtab[procngbnum], TAGBAND,
                       grafptr->proccomm, nrcvreqtab + procngbnum) != MPI_SUCCESS) {
      errorPrint ("dgraphBandPtop: communication error (2)");
      return     (1);
    }
  } while (procngbnum != procngbnxt);

  bandvertlvlnum =                                /* Start index of last level is start index */
  bandvertlocnnd = grafptr->baseval;              /* Reset number of band vertices, plus base */
  bandedgelocnbr = 0;
  memSet (vnumgsttax, ~0, grafptr->vertgstnbr * sizeof (Gnum)); /* Reset part array */
  vnumgsttax -= grafptr->baseval;

  vertlocnnd = grafptr->vertlocnnd;
  for (queulocnum = 0; queulocnum < queulocnbr; queulocnum ++) { /* All frontier vertices will be first vertices of band graph */
    Gnum              vertlocnum;

    vertlocnum = queuloctab[queulocnum];
    vnumgsttax[vertlocnum] = bandvertlocnnd ++;   /* Keep frontier vertex in band               */
    bandedgelocnbr += vendloctax[vertlocnum] - vertloctax[vertlocnum]; /* Account for its edges */
  }
  queuheadidx = 0;                                /* No queued vertex read yet                  */
  queutailidx = queulocnbr;                       /* All frontier vertices are already in queue */

  for (distval = 0; ++ distval <= distmax; ) {
    Gnum              queunextidx;                /* Tail index for enqueuing vertices of next band */
    int               vrcvreqnbr;

    if (MPI_Startall (procngbnbr, nrcvreqtab) != MPI_SUCCESS) { /* Start all receive operations from neighbors */
      errorPrint ("dgraphBandPtop: communication error (3)");
      return     (1);
    }

    bandvertlvlnum = bandvertlocnnd;              /* Save start index of current level, based */
    *bandvertlvlptr = bandvertlvlnum;

    memCpy (nsndidxtab, nsnddsptab, procngbnbr * sizeof (Gnum)); /* Reset send buffer indices */

    for (queunextidx = queutailidx; queuheadidx < queutailidx; ) { /* For all vertices in queue */
      Gnum              vertlocnum;
      Gnum              edgelocnum;

      vertlocnum = queuloctab[queuheadidx ++];    /* Dequeue vertex */
      for (edgelocnum = vertloctax[vertlocnum]; edgelocnum < vendloctax[vertlocnum]; edgelocnum ++) {
        Gnum              vertlocend;

        vertlocend = edgegsttax[edgelocnum];
        if (vnumgsttax[vertlocend] != ~0)         /* If end vertex has already been processed */
          continue;                               /* Skip to next vertex                      */

        if (vertlocend < vertlocnnd) {            /* If end vertex is local                           */
          vnumgsttax[vertlocend] = bandvertlocnnd ++; /* Label vertex as enqueued                     */
          queuloctab[queunextidx ++] = vertlocend; /* Enqueue vertex for next pass                    */
          bandedgelocnbr += vendloctax[vertlocend] - vertloctax[vertlocend]; /* Account for its edges */
        }
        else {                                    /* End vertex is a ghost */
          Gnum              vertglbend;
          int               procngbnum;
          int               procngbmax;

          vnumgsttax[vertlocend] = 0;             /* Label ghost vertex as enqueued  */
          vertglbend = edgeloctax[edgelocnum];    /* Get global number of end vertex */

          procngbnum = 0;
          procngbmax = procngbnbr;
          while ((procngbmax - procngbnum) > 1) { /* Find owner process by dichotomy on procvgbtab */
            int                 procngbmed;

            procngbmed = (procngbmax + procngbnum) / 2;
            if (procvgbtab[procngbmed] > vertglbend)
              procngbmax = procngbmed;
            else
              procngbnum = procngbmed;
          }
#ifdef SCOTCH_DEBUG_DGRAPH2
          if ((grafptr->procvrttab[grafptr->procngbtab[procngbnum]]     >  vertglbend) ||
              (grafptr->procvrttab[grafptr->procngbtab[procngbnum] + 1] <= vertglbend) ||
              (nsndidxtab[procngbnum] >= nsnddsptab[procngbnum + 1])) {
            errorPrint ("dgraphBandPtop: internal error (1)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
          vsnddattab[nsndidxtab[procngbnum] ++] = vertglbend - procvgbtab[procngbnum] + grafptr->baseval; /* Buffer local value on neighbor processor */
        }
      }
    }

    procngbnum = procngbnxt;                      /* Send all buffers to neighbors */
    do {
      int               procglbnum;

      procglbnum = grafptr->procngbtab[procngbnum];

      if (MPI_Isend (vsnddattab + nsnddsptab[procngbnum], nsndidxtab[procngbnum] - nsnddsptab[procngbnum],
                     GNUM_MPI, grafptr->procngbtab[procngbnum], TAGBAND, grafptr->proccomm,
                     nsndreqtab + procngbnum) != MPI_SUCCESS) {
        errorPrint ("dgraphBandPtop: communication error (4)");
        return     (1);
      }
      procngbnum = (procngbnum + 1) % procngbnbr; /* Post-increment neighbor rank */
    } while (procngbnum != procngbnxt);

    for (vrcvreqnbr = procngbnbr; vrcvreqnbr > 0; vrcvreqnbr --) { /* For all pending receive requests */
      Gnum * restrict   vrcvdatptr;
      int               vertrcvnum;
      MPI_Status        statdat;
      int               statsiz;
      int               o;

#ifdef SCOTCH_DETERMINISTIC
      procngbnum = vrcvreqnbr - 1;
      o = MPI_Wait (&nrcvreqtab[procngbnum], &statdat);
#else /* SCOTCH_DETERMINISTIC */
      o = MPI_Waitany (procngbnbr, nrcvreqtab, &procngbnum, &statdat);
#endif /* SCOTCH_DETERMINISTIC */
      if ((o != MPI_SUCCESS) ||
          (MPI_Get_count (&statdat, GNUM_MPI, &statsiz) != MPI_SUCCESS)) {
        errorPrint ("dgraphBandPtop: communication error (5)");
        return     (1);
      }
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (statdat.MPI_SOURCE != grafptr->procngbtab[procngbnum]) {
        errorPrint ("dgraphBandPtop: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

      vrcvdatptr = vrcvdattab + nrcvdsptab[procngbnum];
      for (vertrcvnum = 0; vertrcvnum < statsiz; vertrcvnum ++) {
        Gnum              vertlocend;

        vertlocend = vrcvdatptr[vertrcvnum];      /* Get local vertex from message */
#ifdef SCOTCH_DEBUG_DGRAPH2
        if ((vertlocend < grafptr->baseval) || (vertlocend >= vertlocnnd)) {
          errorPrint ("dgraphBandPtop: internal error (3)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
        if (vnumgsttax[vertlocend] != ~0)         /* If end vertex has already been processed */
          continue;                               /* Skip to next vertex                      */

        vnumgsttax[vertlocend] = bandvertlocnnd ++; /* Label vertex as enqueued                     */
        queuloctab[queunextidx ++] = vertlocend;  /* Enqueue vertex for next pass                   */
        bandedgelocnbr += vendloctax[vertlocend] - vertloctax[vertlocend]; /* Account for its edges */
      }
    }

    queutailidx = queunextidx;                    /* Prepare queue for next sweep */

    if (MPI_Waitall (procngbnbr, nsndreqtab, MPI_STATUSES_IGNORE) != MPI_SUCCESS) { /* Wait until all send operations completed */
      errorPrint ("dgraphBandPtop: communication error (6)");
      return     (1);
    }
  }
  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) { /* Free persistent receive requests */
    if (MPI_Request_free (nrcvreqtab + procngbnum) != MPI_SUCCESS) {
      errorPrint ("dgraphBandPtop: communication error (7)");
      return     (1);
    }
  }

  memFree (procvgbtab);                           /* Free group leader */

  *vnumgstptr     = vnumgsttax;
  *bandvertlocptr = bandvertlocnnd - grafptr->baseval;
  *bandedgelocptr = bandedgelocnbr;

  return (0);
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine computes a distributed index array
** of given width around the current separator.
** It returns:
** - 0   : if the index array could be computed.
** - !0  : on error.
*/

int
dgraphBand (
Dgraph * restrict const             grafptr,      /*+ Distributed graph                                      +*/
const Gnum                          fronlocnbr,   /*+ Number of frontier vertices                            +*/
Gnum * restrict const               fronloctab,   /*+ Array of frontier vertices, re-used as queue array     +*/
const GraphPart * restrict const    partgsttax,   /*+ Part array for original graph ({0,1} or {0,1,2})       +*/
const Gnum                          complocload0, /*+ Load in part 0 or {0,2}                                +*/
const Gnum                          complocload1, /*+ Load in part 1                                         +*/
Gnum                                distmax,      /*+ Maximum distance from separator vertices               +*/
Dgraph * restrict const             bandgrafptr,  /*+ Pointer to band graph structure to fill                +*/
Gnum * restrict * const             bandfronlocptr, /*+ Pointer to bandfronloctab                            +*/
GraphPart * restrict * const        bandpartgstptr, /*+ Pointer to bandpartgsttax                            +*/
Gnum * const                        bandvertlvlptr, /*+ Pointer to based start index of last level           +*/
Gnum * const                        bandvertlocptr1, /*+ Pointer to number of band vertices in part 1        +*/
Gnum * const                        bandvertlocancptr) /*+ Pointer to flag set if anchor vertices overloaded +*/
{
  Gnum                    bandvertlocnnd;         /* End of local band vertex array, (without anchor vertices) */
  Gnum                    bandvertlocnbr;         /* Number of local band vertices (including anchor vertices) */
  Gnum                    bandvertlocnbr1;        /* Number of band graph vertices in part 1 except anchor 1   */
  Gnum                    bandvertlvlnum;         /* Index of first band vertex belonging to last level        */
  Gnum * restrict         bandvertloctax;
  Gnum                    bandvertlocadj;         /* Ajust value for local-to-global band vertex indices       */
  Gnum                    bandvertlocancadj;      /* Flag set when anchor(s) represent unexistent vertices     */
  Gnum                    bandvertlocnum;
  Gnum                    bandvelolocsum;
  Gnum                    bandvelolocsum1;
  Gnum * restrict         bandedgeloctax;
  Gnum                    bandedgelocnum;
  Gnum                    bandedgeloctmp;
  Gnum                    bandedgelocnbr;         /* Number of local edges in band graph                       */
  Gnum * restrict         bandedloloctax;
  Gnum                    bandedlolocnbr;         /* Size of local band edge load array                        */
  Gnum * restrict         bandfronloctab;
  GraphPart * restrict    bandpartgsttax;
  Gnum * restrict         bandvnumgsttax;         /* Indices of selected band vertices in band graph           */
  Gnum                    banddegrlocmax;
  Gnum                    degrval;
  Gnum                    veloval;
  const Gnum * restrict   edgegsttax;
  Gnum                    fronlocnum;
  int                     cheklocval;
  int                     procngbnum;

  if (dgraphGhst (grafptr) != 0) {                /* Compute ghost edge array if not already present */
    errorPrint ("dgraphBand: cannot compute ghost edge array");
    return     (1);
  }

  if ((((grafptr->flagval & DGRAPHCOMMPTOP) != 0) ? dgraphBandPtop : dgraphBandColl)
      (grafptr, fronlocnbr, fronloctab, distmax, &bandvnumgsttax, &bandvertlvlnum, &bandvertlocnbr, &bandedgelocnbr) != 0)
    return (1);

  if (bandvertlvlptr != NULL)
    *bandvertlvlptr = bandvertlvlnum;

  bandedgelocnbr += 2 * ((bandvertlocnbr + grafptr->baseval - bandvertlvlnum) + (grafptr->procglbnbr - 1)); /* Add edges to and from anchors */
  bandvertlocnbr += 2;                            /* Add anchor vertices */
  bandedlolocnbr  = (grafptr->edloloctax != NULL) ? bandedgelocnbr : 0;

  dgraphInit (bandgrafptr, grafptr->proccomm);
  bandgrafptr->flagval = (DGRAPHFREEALL ^ DGRAPHFREECOMM) | DGRAPHVERTGROUP | DGRAPHEDGEGROUP; /* Arrays created by the routine itself */
  bandgrafptr->baseval = grafptr->baseval;

  cheklocval = 0;
  if (memAllocGroup ((void **) (void *)           /* Allocate distributed graph private data */
                     &bandgrafptr->procdsptab, (size_t) ((grafptr->procglbnbr + 1) * sizeof (Gnum)),
                     &bandgrafptr->proccnttab, (size_t) (grafptr->procglbnbr       * sizeof (Gnum)),
                     &bandgrafptr->procngbtab, (size_t) (grafptr->procglbnbr       * sizeof (int)),
                     &bandgrafptr->procrcvtab, (size_t) (grafptr->procglbnbr       * sizeof (int)),
                     &bandgrafptr->procsndtab, (size_t) (grafptr->procglbnbr       * sizeof (int)), NULL) == NULL) {
    errorPrint ("dgraphBand: out of memory (1)");
    cheklocval = 1;
  }
  else if (memAllocGroup ((void **) (void *)      /* Allocate distributed graph public data */
                          &bandgrafptr->vertloctax, (size_t) ((bandvertlocnbr + 1) * sizeof (Gnum)), /* Compact vertex array */
                          &bandgrafptr->vnumloctax, (size_t) (bandvertlocnbr       * sizeof (Gnum)),
                          &bandgrafptr->veloloctax, (size_t) (bandvertlocnbr       * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphBand: out of memory (2)");
    cheklocval = 1;
  }
  else if (bandgrafptr->vertloctax -= bandgrafptr->baseval,
           bandgrafptr->veloloctax -= bandgrafptr->baseval,
           bandgrafptr->vnumloctax -= bandgrafptr->baseval,
           (memAllocGroup ((void **) (void *)
                           &bandedgeloctax, (size_t) (bandedgelocnbr * sizeof (Gnum)),
                           &bandedloloctax, (size_t) (bandedlolocnbr * sizeof (Gnum)), NULL) == NULL)) {
    errorPrint ("dgraphBand: out of memory (3)");
    cheklocval = 1;
  }
  else {
    bandedgeloctax -= bandgrafptr->baseval;
    bandedloloctax  = (grafptr->edloloctax != NULL) ? (bandedloloctax - bandgrafptr->baseval) : NULL;

    if ((bandfronloctab = memAlloc (bandvertlocnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("dgraphBand: out of memory (4)");
      cheklocval = 1;
    }
    else if ((bandpartgsttax = memAlloc ((bandvertlocnbr + bandedgelocnbr) * sizeof (GraphPart))) == NULL) { /* Upper bound on number of ghost vertices */
      errorPrint ("dgraphBand: out of memory (5)");
      cheklocval = 1;
    }
    else
      bandpartgsttax -= bandgrafptr->baseval;
  }

  if (cheklocval != 0) {                          /* In case of memory error */
    bandgrafptr->procdsptab[0] = -1;
    if (MPI_Allgather (&bandgrafptr->procdsptab[0], 1, GNUM_MPI, /* Send received data to dummy array */
                       bandvnumgsttax + bandgrafptr->baseval, 1, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphBand: communication error (1)");
      return     (1);
    }
    if (bandfronloctab != NULL)
      memFree (bandfronloctab);
    dgraphExit (bandgrafptr);
    memFree    (bandvnumgsttax + bandgrafptr->baseval);
    return     (1);
  }
  else {
    bandgrafptr->procdsptab[0] = bandvertlocnbr;
    if (MPI_Allgather (&bandgrafptr->procdsptab[0], 1, GNUM_MPI,
                       &bandgrafptr->procdsptab[1], 1, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphBand: communication error (2)");
      return     (1);
    }
  }
  bandgrafptr->procdsptab[0] = bandgrafptr->baseval; /* Build vertex-to-process array */
#ifdef SCOTCH_DEBUG_DGRAPH2
  memSet (bandgrafptr->vnumloctax + bandgrafptr->baseval, ~0, (bandvertlocnbr * sizeof (Gnum)));
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  for (procngbnum = 1; procngbnum <= grafptr->procglbnbr; procngbnum ++) { /* Process potential error flags from other processes */
    if (bandgrafptr->procdsptab[procngbnum] < 0) { /* If error notified by another process                                       */
      if (bandpartgsttax != NULL)
        memFree (bandpartgsttax + bandgrafptr->baseval);
      if (bandfronloctab != NULL)
        memFree (bandfronloctab);
      dgraphExit (bandgrafptr);
      memFree    (bandvnumgsttax + bandgrafptr->baseval);
      return     (1);
    }
    bandgrafptr->procdsptab[procngbnum]    += bandgrafptr->procdsptab[procngbnum - 1];
    bandgrafptr->proccnttab[procngbnum - 1] = bandgrafptr->procdsptab[procngbnum] - bandgrafptr->procdsptab[procngbnum - 1];
  }

  for (fronlocnum = 0, bandvertlocnum = bandgrafptr->baseval, bandvertlocadj = bandgrafptr->procdsptab[grafptr->proclocnum] - bandgrafptr->baseval;
       fronlocnum < fronlocnbr; fronlocnum ++, bandvertlocnum ++) { /* Turn all graph frontier vertices into band frontier vertices */
    Gnum              vertlocnum;

    bandfronloctab[fronlocnum] = bandvertlocnum;  /* All frontier vertices are first vertices of band graph */
    vertlocnum = fronloctab[fronlocnum];
    bandgrafptr->vnumloctax[bandvertlocnum] = vertlocnum;
    bandvnumgsttax[vertlocnum] += bandvertlocadj; /* Turn local indices in band graph into global indices */
  }
  for (bandvertlocnnd = bandvertlocnbr + bandgrafptr->baseval - 2; /* Pick selected band vertices from rest of frontier array without anchors */
       bandvertlocnum < bandvertlocnnd; fronlocnum ++, bandvertlocnum ++) {
    Gnum              vertlocnum;

    vertlocnum = fronloctab[fronlocnum];
    bandgrafptr->vnumloctax[bandvertlocnum] = vertlocnum;
    bandvnumgsttax[vertlocnum] += bandvertlocadj; /* Turn local indices in band graph into global indices */
  }
  bandgrafptr->vnumloctax[bandvertlocnnd]     =   /* Prevent Valgrind from yelling when centralizing band graphs */
  bandgrafptr->vnumloctax[bandvertlocnnd + 1] = -1;

  if (dgraphHaloSync (grafptr, (byte *) (bandvnumgsttax + bandgrafptr->baseval), GNUM_MPI) != 0) { /* Share global indexing of halo vertices */
    errorPrint ("dgraphBand: cannot perform halo exchange");
    return     (1);
  }

  edgegsttax = grafptr->edgegsttax;

  veloval = 1;
  bandvertloctax  = bandgrafptr->vertloctax;
  bandvertlocnbr1 = 0;
  bandvelolocsum  = 0;
  bandvelolocsum1 = 0;
  banddegrlocmax  = 0;
  for (bandvertlocnum = bandedgelocnum = bandgrafptr->baseval; /* Build global vertex array of band graph             */
       bandvertlocnum < bandvertlvlnum; bandvertlocnum ++) { /* For all vertices that do not belong to the last level */
    Gnum              vertlocnum;
    Gnum              edgelocnum;
    Gnum              degrval;
    GraphPart         partval;
    Gnum              partval1;

    vertlocnum = bandgrafptr->vnumloctax[bandvertlocnum];
    partval    = partgsttax[vertlocnum];
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (partval > 2) {
      errorPrint ("dgraphBand: internal error (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    partval1 = partval & 1;
    bandvertlocnbr1 += partval1;                  /* Count vertices in part 1 */
    bandpartgsttax[bandvertlocnum] = partval;
    bandvertloctax[bandvertlocnum] = bandedgelocnum;
    if (grafptr->veloloctax != NULL) {
      veloval = grafptr->veloloctax[vertlocnum];
      bandvelolocsum  += veloval;
      bandvelolocsum1 += veloval & (- partval1);  /* Sum vertex load if (partval == 1) */
    }
    bandgrafptr->veloloctax[bandvertlocnum] = veloval;

    degrval = grafptr->vendloctax[vertlocnum] - grafptr->vertloctax[vertlocnum];
    if (banddegrlocmax < degrval)
      banddegrlocmax = degrval;

    for (edgelocnum = grafptr->vertloctax[vertlocnum]; /* For all original edges */
         edgelocnum < grafptr->vendloctax[vertlocnum]; edgelocnum ++) {
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (bandvnumgsttax[edgegsttax[edgelocnum]] == ~0) { /* All ends should belong to the band graph too */
        errorPrint ("dgraphBand: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      bandedgeloctax[bandedgelocnum ++] = bandvnumgsttax[edgegsttax[edgelocnum]];
    }
  }
  for ( ; bandvertlocnum < bandvertlocnnd; bandvertlocnum ++) { /* For all vertices that belong to the last level except anchors */
    Gnum              vertlocnum;
    Gnum              edgelocnum;
    Gnum              degrval;
    GraphPart         partval;
    Gnum              partval1;

    vertlocnum = bandgrafptr->vnumloctax[bandvertlocnum];
    partval    = partgsttax[vertlocnum];
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (partval > 2) {
      errorPrint ("dgraphBand: internal error (3)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    partval1 = partval & 1;
    bandvertlocnbr1 += partval1;                  /* Count vertices in part 1 */
    bandpartgsttax[bandvertlocnum] = partval;
    bandvertloctax[bandvertlocnum] = bandedgelocnum;
    if (grafptr->veloloctax != NULL) {
      veloval = grafptr->veloloctax[vertlocnum];
      bandvelolocsum  += veloval;
      bandvelolocsum1 += veloval & (- partval1);  /* Sum vertex load if (partval == 1) */
    }
    bandgrafptr->veloloctax[bandvertlocnum] = veloval;

    for (edgelocnum = grafptr->vertloctax[vertlocnum]; /* For all original edges */
         edgelocnum < grafptr->vendloctax[vertlocnum]; edgelocnum ++) {
      Gnum              bandvertlocend;

      bandvertlocend = bandvnumgsttax[edgegsttax[edgelocnum]];
      if (bandvertlocend != ~0) {                 /* If end vertex belongs to band graph  */
        if (bandedloloctax != NULL)               /* If graph has edge weights, copy load */
          bandedloloctax[bandedgelocnum] = grafptr->edloloctax[edgelocnum];
        bandedgeloctax[bandedgelocnum ++] = bandvertlocend;
      }
    }
    if (bandedloloctax != NULL)                   /* If graph has edge weights */
      bandedloloctax[bandedgelocnum] = 1;         /* Edge to anchor has load 1 */
    bandedgeloctax[bandedgelocnum ++] = bandvertlocnnd + bandvertlocadj + partval1; /* Add edge to anchor of proper part */

    degrval = bandedgelocnum - bandvertloctax[bandvertlocnum];
    if (banddegrlocmax < degrval)
      banddegrlocmax = degrval;
  }

  memFree (bandvnumgsttax + bandgrafptr->baseval);  /* Free useless space */

  bandpartgsttax[bandvertlocnnd]     = 0;         /* Set parts of anchor vertices */
  bandpartgsttax[bandvertlocnnd + 1] = 1;
  bandvertloctax[bandvertlocnum] = bandedgelocnum; /* Process anchor vertex in part 0                 */
  for (procngbnum = 0; procngbnum < grafptr->proclocnum; procngbnum ++) /* Build clique with anchors of part 0 */
    bandedgeloctax[bandedgelocnum ++] = bandgrafptr->procdsptab[procngbnum + 1] - 2;
  for (procngbnum ++; procngbnum < grafptr->procglbnbr; procngbnum ++) /* Build clique with anchors of part 0 */
    bandedgeloctax[bandedgelocnum ++] = bandgrafptr->procdsptab[procngbnum + 1] - 2;
  bandedgeloctmp = bandedgelocnum + (bandvertlocnnd - bandvertlvlnum);
  for (procngbnum = 0; procngbnum < grafptr->proclocnum; procngbnum ++) /* Build clique with anchors of part 1 */
    bandedgeloctax[bandedgeloctmp ++] = bandgrafptr->procdsptab[procngbnum + 1] - 1;
  for (procngbnum ++; procngbnum < grafptr->procglbnbr; procngbnum ++) /* Build clique with anchors of part 1 */
    bandedgeloctax[bandedgeloctmp ++] = bandgrafptr->procdsptab[procngbnum + 1] - 1;
  bandvertloctax[bandvertlocnnd + 2] = bandedgeloctmp;
  bandedgelocnbr = bandedgeloctmp - bandgrafptr->baseval; /* Set real number of edges */
  for (bandvertlocnum = bandvertlvlnum, bandedgeloctmp = bandedgelocnum + (bandvertlocnnd - bandvertlvlnum); /* Link vertices of last level to anchors */
       bandvertlocnum < bandvertlocnnd; bandvertlocnum ++) {
    if (bandpartgsttax[bandvertlocnum] == 0)
      bandedgeloctax[bandedgelocnum ++] = bandvertlocnum + bandvertlocadj;
    else
      bandedgeloctax[-- bandedgeloctmp] = bandvertlocnum + bandvertlocadj;
  }
  bandvertloctax[bandvertlocnnd + 1] = bandedgeloctmp;
  degrval = bandvertloctax[bandvertlocnnd + 1] - bandvertloctax[bandvertlocnnd];
  if (banddegrlocmax < degrval)
    banddegrlocmax = degrval;
  degrval = bandvertloctax[bandvertlocnnd + 2] - bandvertloctax[bandvertlocnnd + 1];
  if (banddegrlocmax < degrval)
    banddegrlocmax = degrval;

  if (bandedloloctax != NULL) {                   /* If graph has edge weights */
    Gnum              edgelocnum;
    Gnum              edgelocnnd;

    for (bandvertlocnum = bandgrafptr->baseval;   /* For all vertices that do not belong to the last level */
         bandvertlocnum < bandvertlvlnum; bandvertlocnum ++) { 
      Gnum              vertlocnum;
      Gnum              bandedgelocnum;

      vertlocnum     = bandgrafptr->vnumloctax[bandvertlocnum];
      bandedgelocnum = bandvertloctax[bandvertlocnum];
      memCpy (bandedloloctax + bandedgelocnum,    /* Copy edge load array */
              &grafptr->edloloctax[grafptr->vertloctax[vertlocnum]],
              (bandvertloctax[bandvertlocnum + 1] - bandedgelocnum) * sizeof (Gnum));
    }                                             /* Vertices of last level have been processed before */
    for (edgelocnum = bandvertloctax[bandvertlocnnd], /* Loads of anchor edges are all 1's too         */
         edgelocnnd = bandvertloctax[bandvertlocnnd + 2];
         edgelocnum < edgelocnnd; edgelocnum ++)
      bandedloloctax[edgelocnum] = 1;
  }

  if (grafptr->veloloctax == NULL) {              /* If original graph is not weighted */
    bandgrafptr->veloloctax[bandvertlocnnd]     = complocload0 + bandvertlocnbr1 - bandvertlocnbr + 2; /* Plus 2 for anchors */
    bandgrafptr->veloloctax[bandvertlocnnd + 1] = complocload1 - bandvertlocnbr1;
  }
  else {
    bandgrafptr->veloloctax[bandvertlocnnd]     = complocload0 + bandvelolocsum1 - bandvelolocsum;
    bandgrafptr->veloloctax[bandvertlocnnd + 1] = complocload1 - bandvelolocsum1;
  }
  bandvertlocancadj = 0;
  if ((bandgrafptr->veloloctax[bandvertlocnnd]     == 0) || /* If at least one anchor is empty */
      (bandgrafptr->veloloctax[bandvertlocnnd + 1] == 0)) {
    bandvertlocancadj = 1;
    bandgrafptr->veloloctax[bandvertlocnnd] ++;   /* Increase weight of both anchors to keep balance */
    bandgrafptr->veloloctax[bandvertlocnnd + 1] ++;
  }

  bandgrafptr->procvrttab = bandgrafptr->procdsptab; /* Graph does not have holes */
  bandgrafptr->vertlocnbr = bandvertlocnbr;
  bandgrafptr->vertlocnnd = bandvertlocnbr + bandgrafptr->baseval;
  bandgrafptr->vendloctax = bandvertloctax + 1;   /* Band graph is compact */
  bandgrafptr->velolocsum = grafptr->velolocsum + 2 * bandvertlocancadj;
  bandgrafptr->edgeloctax = bandedgeloctax;
  bandgrafptr->edloloctax = bandedloloctax;
  bandgrafptr->edgelocnbr = bandedgelocnbr;
  bandgrafptr->edgelocsiz = bandedgelocnbr;
  bandgrafptr->degrglbmax = banddegrlocmax;       /* Local maximum degree will be turned into global maximum degree */
  if (dgraphBuild4 (bandgrafptr) != 0) {
    errorPrint ("dgraphBand: cannot build band graph");
    return     (1);
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphCheck (bandgrafptr) != 0) {
    errorPrint ("dgraphBand: internal error (4)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  *bandfronlocptr    = bandfronloctab;
  *bandpartgstptr    = bandpartgsttax;
  *bandvertlocptr1   = bandvertlocnbr1;
  *bandvertlocancptr = bandvertlocancadj;

  return (0);
}
