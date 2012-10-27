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
/**   NAME       : dgraph_induce.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Jun-Ho HER                              **/
/**                                                        **/
/**   FUNCTION   : This module handles the source graph    **/
/**                subgraph-making functions.              **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 08 apr 2006     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 31 mar 2008     **/
/**                                 to   : 30 jul 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DGRAPH
#define DGRAPH_INDUCE

#include "module.h"
#include "common.h"
#include "dgraph.h"

/****************************************/
/*                                      */
/* These routines handle source graphs. */
/*                                      */
/****************************************/

/* This routine builds the graph induced
** by the original graph and the list of
** selected vertices.
** The induced vnumtab array is the global
** translation of the list array if the
** original graph does not have a vnumtab,
** or the proper subset of the original
** vnumtab else.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphInduceList (
Dgraph * restrict const       orggrafptr,
const Gnum                    indlistnbr,
const Gnum * restrict const   indlisttab,         /* Local list of kept vertices */
Dgraph * restrict const       indgrafptr)
{
  Gnum * restrict       orgindxgsttax;            /* Based access to vertex translation array       */
  Gnum                  indvertlocnnd;            /* Based index of end of local vertex array       */
  Gnum                  indvertlocnum;            /* Number of current vertex in induced graph      */
  Gnum                  indvertglbnum;            /* Number of current vertex in global ordering    */
  Gnum                  indvelolocnbr;            /* Size of local vertex load array                */
  Gnum                  indvelolocsum;            /* Sum of vertex loads                            */
  Gnum                  indedgelocmax;            /* (Approximate) number of edges in induced graph */
  Gnum                  indedgelocnbr;            /* Real number of edges in induced graph          */
  Gnum                  indedgelocnum;
  Gnum * restrict       indedloloctax;
  Gnum                  inddegrlocmax;            /* Local maximum degree                           */
  const Gnum * restrict indlisttax;
  int                   cheklocval;
  int                   chekglbval;

  if (dgraphGhst (orggrafptr) != 0) {             /* Compute ghost edge array if not already present */
    errorPrint ("dgraphInduceList: cannot compute ghost edge array");
    return     (1);
  }

  dgraphInit (indgrafptr, orggrafptr->proccomm);
  indgrafptr->flagval = (DGRAPHFREEALL ^ DGRAPHFREECOMM) | DGRAPHVERTGROUP | DGRAPHEDGEGROUP;

  if (orggrafptr->veloloctax != NULL) {
    indvelolocnbr = indlistnbr;
    indvelolocsum = 0;
  }
  else {
    indvelolocnbr = 0;
    indvelolocsum = indlistnbr;
  }
  indedgelocmax = orggrafptr->edgelocnbr;         /* Choose best upper bound on number of edges (avoid multiply overflow) */
  if ((orggrafptr->degrglbmax > 0) && (indlistnbr < (indedgelocmax / orggrafptr->degrglbmax)))
    indedgelocmax = indlistnbr * orggrafptr->degrglbmax;
  if (orggrafptr->edloloctax != NULL)             /* If graph has edge weights */
    indedgelocmax *= 2;                           /* Account for edge weights  */

  cheklocval =
  chekglbval = 0;
  if (memAllocGroup ((void **) (void *)           /* Allocate distributed graph private data */
                     &indgrafptr->procdsptab, (size_t) ((orggrafptr->procglbnbr + 1) * sizeof (Gnum)),
                     &indgrafptr->proccnttab, (size_t) (orggrafptr->procglbnbr       * sizeof (Gnum)),
                     &indgrafptr->procngbtab, (size_t) (orggrafptr->procglbnbr       * sizeof (int)),
                     &indgrafptr->procrcvtab, (size_t) (orggrafptr->procglbnbr       * sizeof (int)),
                     &indgrafptr->procsndtab, (size_t) (orggrafptr->procglbnbr       * sizeof (int)), NULL) == NULL) {
    errorPrint ("dgraphInduceList: out of memory (1)");
    cheklocval = 1;
  }
  else if (memAllocGroup ((void **) (void *)      /* Allocate distributed graph public data */
                          &indgrafptr->vertloctax, (size_t) ((indlistnbr + 1) * sizeof (Gnum)), /* Compact vertex array */
                          &indgrafptr->vnumloctax, (size_t) (indlistnbr       * sizeof (Gnum)),
                          &indgrafptr->veloloctax, (size_t) (indvelolocnbr    * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphInduceList: out of memory (2)");
    cheklocval = 1;
  }
  else if (indgrafptr->vertloctax -= orggrafptr->baseval,
           indgrafptr->vnumloctax -= orggrafptr->baseval,
           indgrafptr->veloloctax  = (orggrafptr->veloloctax != NULL) ? indgrafptr->veloloctax - orggrafptr->baseval : NULL,
           memAllocGroup ((void **) (void *)
                          &indgrafptr->edgeloctax, (size_t) (indedgelocmax          * sizeof (Gnum)), /* Pre-allocate space for edgetab (and edlotab) */
                          &orgindxgsttax,          (size_t) (orggrafptr->vertgstnbr * sizeof (Gnum)), NULL) == NULL) { /* orgindxgsttab is at the end */
    errorPrint ("dgraphInduceList: out of memory (3)");
    cheklocval = 1;
  }
  else
    indgrafptr->edgeloctax -= orggrafptr->baseval;

  if (cheklocval != 0) {                          /* In case of memory error */
    Gnum                procngbnum;
    Gnum                dummyval;

    dummyval   = -1;
    chekglbval = 1;
    if (MPI_Allgather (&dummyval, 1, GNUM_MPI,    /* Use proccnttab of orggraf as dummy receive array (will be regenerated) */
                       orggrafptr->proccnttab, 1, GNUM_MPI, indgrafptr->proccomm) != MPI_SUCCESS)
      errorPrint ("dgraphInduceList: communication error (1)");

    for (procngbnum = 1; procngbnum <= orggrafptr->procglbnbr; procngbnum ++) /* Rebuild proccnttab of orggraf */
      orggrafptr->proccnttab[procngbnum - 1] = orggrafptr->procdsptab[procngbnum] - orggrafptr->procdsptab[procngbnum - 1];
  }
  else {
    indgrafptr->procdsptab[0] = indlistnbr;
    if (MPI_Allgather (&indgrafptr->procdsptab[0], 1, GNUM_MPI,
                       &indgrafptr->proccnttab[0], 1, GNUM_MPI, indgrafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphInduceList: communication error (2)");
      chekglbval = 1;
    }
    else {
      Gnum                procngbnum;

      indgrafptr->procdsptab[0] = orggrafptr->baseval; /* Build vertex-to-process array                                                */
      for (procngbnum = 0; procngbnum < indgrafptr->procglbnbr; procngbnum ++) { /* Process potential error flags from other processes */
        if (indgrafptr->procdsptab[procngbnum] < 0) { /* If error notified by another process                                          */
          chekglbval = 1;
          break;
        }
        indgrafptr->procdsptab[procngbnum + 1] = indgrafptr->procdsptab[procngbnum] + indgrafptr->proccnttab[procngbnum];
      }
    }
    indgrafptr->procvrttab = indgrafptr->procdsptab; /* Graph does not have holes */
  }
  if (chekglbval != 0) {                          /* If something went wrong in all of the above */
    dgraphExit (indgrafptr);
    return     (1);
  }

  memSet (orgindxgsttax, ~0, orggrafptr->vertlocnbr * sizeof (Gnum)); /* Preset index array */
  orgindxgsttax -= orggrafptr->baseval;

  indlisttax    = indlisttab - orggrafptr->baseval;
  indvertlocnnd = indlistnbr + orggrafptr->baseval;
  for (indvertlocnum = orggrafptr->baseval, indvertglbnum = indgrafptr->procdsptab[indgrafptr->proclocnum], indedgelocmax = 0; /* Fill index array while recomputing tighter upper bound on arcs */
       indvertlocnum < indvertlocnnd; indvertlocnum ++, indvertglbnum ++) {
    Gnum                orgvertlocnum;

    orgvertlocnum = indlisttax[indvertlocnum];
    orgindxgsttax[orgvertlocnum] = indvertglbnum; /* Mark selected vertices */
    indedgelocmax += orggrafptr->vendloctax[orgvertlocnum] - orggrafptr->vertloctax[orgvertlocnum];
  }

  if (dgraphHaloSync (orggrafptr, (byte *) (orgindxgsttax + orggrafptr->baseval), GNUM_MPI) != 0) { /* Share global indexing of subgraph vertices */
    errorPrint ("dgraphInduceList: cannot perform halo exchange");
    dgraphExit (indgrafptr);
    return     (1);
  }

  indedloloctax = (orggrafptr->edloloctax != NULL) ? indgrafptr->edgeloctax + indedgelocmax : NULL;
  inddegrlocmax = 0;
  for (indvertlocnum = indedgelocnum = orggrafptr->baseval;
       indvertlocnum < indvertlocnnd; indvertlocnum ++) {
    Gnum                orgvertlocnum;
    Gnum                orgedgelocnum;

    orgvertlocnum = indlisttax[indvertlocnum];
    indgrafptr->vertloctax[indvertlocnum] = indedgelocnum;
    if (indgrafptr->veloloctax != NULL) {         /* If graph has vertex weights */
      indvelolocsum +=                            /* Accumulate vertex loads     */
      indgrafptr->veloloctax[indvertlocnum] = orggrafptr->veloloctax[orgvertlocnum];
    }

    for (orgedgelocnum = orggrafptr->vertloctax[orgvertlocnum];
         orgedgelocnum < orggrafptr->vendloctax[orgvertlocnum]; orgedgelocnum ++) {
      Gnum                indvertgstend;

      indvertgstend = orgindxgsttax[orggrafptr->edgegsttax[orgedgelocnum]];
      if (indvertgstend != ~0) {                  /* If edge should be kept */
        indgrafptr->edgeloctax[indedgelocnum] = indvertgstend;
        if (indedloloctax != NULL)
          indedloloctax[indedgelocnum] = orggrafptr->edloloctax[orgedgelocnum];
        indedgelocnum ++;
      }
    }
    if (inddegrlocmax < (indedgelocnum - indgrafptr->vertloctax[indvertlocnum]))
      inddegrlocmax = (indedgelocnum - indgrafptr->vertloctax[indvertlocnum]);
  }
  indgrafptr->vertloctax[indvertlocnum] = indedgelocnum; /* Mark end of edge array */
  indedgelocnbr = indedgelocnum - orggrafptr->baseval;

  if (indedloloctax != NULL) {                    /* Re-allocate arrays and delete orgindxtab             */
    size_t              indedlooftval;            /* Offset of edge load array with respect to edge array */

    indedlooftval = indedloloctax - indgrafptr->edgeloctax;
    indgrafptr->edgeloctax  = memRealloc (indgrafptr->edgeloctax + orggrafptr->baseval, (indedlooftval + indedgelocnbr) * sizeof (Gnum));
    indgrafptr->edgeloctax -= orggrafptr->baseval;
    indedloloctax = indgrafptr->edgeloctax + indedlooftval; /* Use old index into old array as new index to avoid alignment problems */
  }
  else {
    indgrafptr->edgeloctax  = memRealloc (indgrafptr->edgeloctax + orggrafptr->baseval, indedgelocnbr * sizeof (Gnum));
    indgrafptr->edgeloctax -= orggrafptr->baseval;
  }

  if (orggrafptr->vnumloctax != NULL) {           /* Adjust vnumloctax */
    for (indvertlocnum = orggrafptr->baseval; indvertlocnum < indvertlocnnd; indvertlocnum ++)
      indgrafptr->vnumloctax[indvertlocnum] = orggrafptr->vnumloctax[indlisttax[indvertlocnum]];
  }
  else {
    Gnum                orgvertglbadj;

    orgvertglbadj = orggrafptr->procdsptab[orggrafptr->proclocnum] - orggrafptr->baseval; /* Set adjustement for global ordering */
    for (indvertlocnum = orggrafptr->baseval; indvertlocnum < indvertlocnnd; indvertlocnum ++)
      indgrafptr->vnumloctax[indvertlocnum] = indlisttax[indvertlocnum] + orgvertglbadj;
  }

  indgrafptr->baseval    = orggrafptr->baseval;
  indgrafptr->vertlocnbr = indlistnbr;
  indgrafptr->vertlocnnd = indlistnbr + orggrafptr->baseval;
  indgrafptr->vendloctax = indgrafptr->vertloctax + 1; /* Induced graph is compact */
  indgrafptr->velolocsum = indvelolocsum;
  indgrafptr->edgelocnbr = indedgelocnbr;
  indgrafptr->edgelocsiz = indedgelocnbr;
  indgrafptr->edloloctax = indedloloctax;
  indgrafptr->degrglbmax = inddegrlocmax;         /* Local maximum degree will be turned into global maximum degree */
  if (dgraphBuild4 (indgrafptr) != 0) {
    errorPrint ("dgraphInduceList: cannot build induced graph");
    return     (1);
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphCheck (indgrafptr) != 0) {            /* Check graph consistency */
    errorPrint ("dgraphInduceList: inconsistent graph data");
    dgraphExit (indgrafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (0);
}

/* This routine builds the graph induced
** by the original graph and the vector of
** selected vertices.
** The induced vnumtab array is the global
** translation of the list array if the
** original graph does not have a vnumtab,
** or the proper subset of the original
** vnumtab else.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphInducePart (
Dgraph * restrict const           orggrafptr,     /* Pointer to original distributed graph          */
const GraphPart * restrict const  orgpartgsttax,  /* Based array of local vertex partition flags    */  
const Gnum                        indvertnbr,     /* Number of local vertices in selected part      */
const GraphPart                   indpartval,
Dgraph * restrict const           indgrafptr)
{
  Gnum * restrict       orgindxgsttax;            /* Based access to vertex translation array       */
  Gnum                  orgvertlocnum;            /* Number of current vertex in original graph     */
  Gnum                  indvertlocnnd;            /* Based index of end of local vertex array       */
  Gnum                  indvertlocnum;            /* Number of current vertex in induced graph      */
  Gnum                  indvertglbnum;            /* Number of current vertex in global ordering    */
  Gnum                  indvelolocnbr;            /* Size of local vertex load array                */
  Gnum                  indvelolocsum;            /* Sum of vertex loads                            */
  Gnum                  indedgelocmax;            /* (Approximate) number of edges in induced graph */
  Gnum                  indedgelocnbr;            /* Real number of edges in induced graph          */
  Gnum                  indedgelocnum;
  Gnum * restrict       indedloloctax;
  Gnum                  inddegrlocmax;            /* Local maximum degree                           */
  int                   cheklocval;
  int                   chekglbval;

  if (dgraphGhst (orggrafptr) != 0) {             /* Compute ghost edge array if not already present */
    errorPrint ("dgraphInducePart: cannot compute ghost edge array");
    return     (1);
  }

  dgraphInit (indgrafptr, orggrafptr->proccomm);
  indgrafptr->flagval = (DGRAPHFREEALL ^ DGRAPHFREECOMM) | DGRAPHVERTGROUP | DGRAPHEDGEGROUP;
  indvertlocnnd = orggrafptr->baseval + indvertnbr;
  if (orggrafptr->veloloctax != NULL) {
    indvelolocnbr = indvertnbr;
    indvelolocsum = 0;
  }
  else {
    indvelolocnbr = 0;
    indvelolocsum = indvertnbr;
  }
  indedgelocmax = orggrafptr->edgelocnbr;         /* Choose best upper bound on number of edges (avoid multiply overflow) */
  if ((orggrafptr->degrglbmax > 0) && (indvertnbr < (indedgelocmax / orggrafptr->degrglbmax)))
    indedgelocmax = indvertnbr * orggrafptr->degrglbmax;
  if (orggrafptr->edloloctax != NULL)             /* If graph has edge weights */
    indedgelocmax *= 2;                           /* Account for edge weights  */

  cheklocval =
  chekglbval = 0;
  if (memAllocGroup ((void **) (void *)           /* Allocate distributed graph private data */
                     &indgrafptr->procdsptab, (size_t) ((orggrafptr->procglbnbr + 1) * sizeof (Gnum)),
                     &indgrafptr->proccnttab, (size_t) (orggrafptr->procglbnbr       * sizeof (Gnum)),
                     &indgrafptr->procngbtab, (size_t) (orggrafptr->procglbnbr       * sizeof (int)),
                     &indgrafptr->procrcvtab, (size_t) (orggrafptr->procglbnbr       * sizeof (int)),
                     &indgrafptr->procsndtab, (size_t) (orggrafptr->procglbnbr       * sizeof (int)), NULL) == NULL) {
    errorPrint ("dgraphInducePart: out of memory (1)");
    cheklocval = 1;
  }
  else if (memAllocGroup ((void **) (void *)      /* Allocate distributed graph public data */
                          &indgrafptr->vertloctax, (size_t) ((indvertnbr + 1) * sizeof (Gnum)), /* Compact vertex array */
                          &indgrafptr->vnumloctax, (size_t) (indvertnbr       * sizeof (Gnum)),
                          &indgrafptr->veloloctax, (size_t) (indvelolocnbr    * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphInducePart: out of memory (2)");
    cheklocval = 1;
  }
  else if (indgrafptr->vertloctax -= orggrafptr->baseval,
           indgrafptr->vnumloctax -= orggrafptr->baseval,
           indgrafptr->veloloctax  = (orggrafptr->veloloctax != NULL) ? indgrafptr->veloloctax - orggrafptr->baseval : NULL,
           memAllocGroup ((void **) (void *)
                          &indgrafptr->edgeloctax, (size_t) (indedgelocmax          * sizeof (Gnum)), /* Pre-allocate space for edgetab (and edlotab) */
                          &orgindxgsttax,          (size_t) (orggrafptr->vertgstnbr * sizeof (Gnum)), NULL) == NULL) { /* orgindxgsttab is at the end */
    errorPrint ("dgraphInducePart: out of memory (3)");
    cheklocval = 1;
  }
  else
    indgrafptr->edgeloctax -= orggrafptr->baseval;

  if (cheklocval != 0) {                          /* In case of memory error */
    Gnum                procngbnum;
    Gnum                dummyval;

    dummyval   = -1;
    chekglbval = 1;
    if (MPI_Allgather (&dummyval, 1, GNUM_MPI,    /* Use proccnttab of orggraf as dummy receive array (will be regenerated) */
                       orggrafptr->proccnttab, 1, GNUM_MPI, indgrafptr->proccomm) != MPI_SUCCESS)
      errorPrint ("dgraphInducePart: communication error (1)");

    for (procngbnum = 1; procngbnum <= orggrafptr->procglbnbr; procngbnum ++) /* Rebuild proccnttab of orggraf */
      orggrafptr->proccnttab[procngbnum - 1] = orggrafptr->procdsptab[procngbnum] - orggrafptr->procdsptab[procngbnum - 1];
  }
  else {
    indgrafptr->procdsptab[0] = indvertnbr;
    if (MPI_Allgather (&indgrafptr->procdsptab[0], 1, GNUM_MPI,
                       &indgrafptr->proccnttab[0], 1, GNUM_MPI, indgrafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphInducePart: communication error (2)");
      chekglbval = 1;
    }
    else {
      Gnum                procngbnum;

      indgrafptr->procdsptab[0] = orggrafptr->baseval; /* Build vertex-to-process array                                                */
      for (procngbnum = 0; procngbnum < indgrafptr->procglbnbr; procngbnum ++) { /* Process potential error flags from other processes */
        if (indgrafptr->procdsptab[procngbnum] < 0) { /* If error notified by another process                                          */
          chekglbval = 1;
          break;
        }
        indgrafptr->procdsptab[procngbnum + 1] = indgrafptr->procdsptab[procngbnum] + indgrafptr->proccnttab[procngbnum];
      }
    }
    indgrafptr->procvrttab = indgrafptr->procdsptab; /* Graph does not have holes */
  }
  if (chekglbval != 0) {                          /* If something went wrong in all of the above */
    dgraphExit (indgrafptr);
    return     (1);
  }

  memSet (orgindxgsttax, ~0, orggrafptr->vertlocnbr * sizeof (Gnum)); /* Preset index array */
  orgindxgsttax -= orggrafptr->baseval;

  for (indvertlocnum = orggrafptr->baseval, indvertglbnum = indgrafptr->procdsptab[indgrafptr->proclocnum], 
       indedgelocmax = 0, orgvertlocnum = orggrafptr->baseval; /* Fill index array while recomputing tighter upper bound on arcs */
       orgvertlocnum < orggrafptr->vertlocnnd; orgvertlocnum ++) {
    if (orgpartgsttax[orgvertlocnum] == indpartval) {
      orgindxgsttax[orgvertlocnum] = indvertglbnum; /* Mark selected vertices */
      indgrafptr->vnumloctax[indvertlocnum] = orgvertlocnum;
      indedgelocmax += orggrafptr->vendloctax[orgvertlocnum] - orggrafptr->vertloctax[orgvertlocnum];
      indvertlocnum ++, indvertglbnum ++;
    }
    else
      orgindxgsttax[orgvertlocnum] = ~0;
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if ((indvertlocnum - orggrafptr->baseval) != indvertnbr) {
    errorPrint ("dgraphInducePart: inconsistent data (1)");
    dgraphExit (indgrafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if (dgraphHaloSync (orggrafptr, (byte *) (orgindxgsttax + orggrafptr->baseval), GNUM_MPI) != 0) { /* Share global indexing of subgraph vertices */
    errorPrint ("dgraphInducePart: cannot perform halo exchange");
    dgraphExit (indgrafptr);
    return     (1);
  }

  indedloloctax = (orggrafptr->edloloctax != NULL) ? indgrafptr->edgeloctax + indedgelocmax : NULL;
  inddegrlocmax = 0;
  for (indvertlocnum = indedgelocnum = orggrafptr->baseval;
       indvertlocnum < indvertlocnnd; indvertlocnum ++) {
    Gnum                orgvertlocnum;
    Gnum                orgedgelocnum;

    orgvertlocnum = indgrafptr->vnumloctax[indvertlocnum];
    indgrafptr->vertloctax[indvertlocnum] = indedgelocnum;
    if (indgrafptr->veloloctax != NULL) {         /* If graph has vertex weights */
      indvelolocsum +=                            /* Accumulate vertex loads     */
      indgrafptr->veloloctax[indvertlocnum] = orggrafptr->veloloctax[orgvertlocnum];
    }

    for (orgedgelocnum = orggrafptr->vertloctax[orgvertlocnum];
         orgedgelocnum < orggrafptr->vendloctax[orgvertlocnum]; orgedgelocnum ++) {
      Gnum                indvertgstend;

      indvertgstend = orgindxgsttax[orggrafptr->edgegsttax[orgedgelocnum]];
      if (indvertgstend != ~0) {                  /* If edge should be kept */
        indgrafptr->edgeloctax[indedgelocnum] = indvertgstend;
        if (indedloloctax != NULL)
          indedloloctax[indedgelocnum] = orggrafptr->edloloctax[orgedgelocnum];
        indedgelocnum ++;
      }
    }
    if (inddegrlocmax < (indedgelocnum - indgrafptr->vertloctax[indvertlocnum]))
      inddegrlocmax = (indedgelocnum - indgrafptr->vertloctax[indvertlocnum]);
  }
  indedgelocnbr = indedgelocnum - orggrafptr->baseval;
  indgrafptr->vertloctax[indvertlocnum] = indedgelocnum; /* Mark end of edge array */
  indgrafptr->baseval    = orggrafptr->baseval;
  indgrafptr->vertlocnbr = indvertnbr;
  indgrafptr->vertlocnnd = indvertnbr + indgrafptr->baseval;
  indgrafptr->vendloctax = indgrafptr->vertloctax + 1; /* Induced graph is compact */
  indgrafptr->velolocsum = indvelolocsum;
  indgrafptr->edgelocnbr = indedgelocnbr;
  indgrafptr->edgelocsiz = indedgelocnbr;

  if (indedloloctax != NULL) {                    /* Re-allocate arrays and delete orgindxtab             */
    size_t              indedlooftval;            /* Offset of edge load array with respect to edge array */

    indedlooftval = indedloloctax - indgrafptr->edgeloctax;
    indgrafptr->edgeloctax  = memRealloc (indgrafptr->edgeloctax + orggrafptr->baseval, (indedlooftval + indedgelocnbr) * sizeof (Gnum));
    indgrafptr->edgeloctax -= orggrafptr->baseval;
    indedloloctax = indgrafptr->edgeloctax + indedlooftval; /* Use old index into old array as new index to avoid alignment problems */
  }
  else {
    indgrafptr->edgeloctax  = memRealloc (indgrafptr->edgeloctax + orggrafptr->baseval, indedgelocnbr * sizeof (Gnum));
    indgrafptr->edgeloctax -= orggrafptr->baseval;
  }

  indvertlocnum = indgrafptr->baseval;
  if (orggrafptr->vnumloctax != NULL) {           /* Adjust vnumloctax */
    for ( ; indvertlocnum < indvertlocnnd; indvertlocnum ++)
      indgrafptr->vnumloctax[indvertlocnum] = orggrafptr->vnumloctax[indgrafptr->vnumloctax[indvertlocnum]];
  }
  else {
    Gnum                orgvertglbadj;

    orgvertglbadj = orggrafptr->procvrttab[orggrafptr->proclocnum] - orggrafptr->baseval; /* Set adjustement for global indexing */
    for ( ; indvertlocnum < indvertlocnnd; indvertlocnum ++)
      indgrafptr->vnumloctax[indvertlocnum] += orgvertglbadj;
  }

  indgrafptr->edloloctax = indedloloctax;
  indgrafptr->degrglbmax = inddegrlocmax;         /* Local maximum degree will be turned into global maximum degree */
  if (dgraphBuild4 (indgrafptr) != 0) {
    errorPrint ("dgraphInducePart: cannot build induced graph");
    return     (1);
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphCheck (indgrafptr) != 0) {            /* Check graph consistency */
    errorPrint ("dgraphInducePart: inconsistent graph data (2)");
    dgraphExit (indgrafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (0);
}
