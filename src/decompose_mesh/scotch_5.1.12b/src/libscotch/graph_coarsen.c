/* Copyright 2004,2007,2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : graph_coarsen.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the source graph   **/
/**                coarsening functions.                   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to     18 may 1993     **/
/**                # Version 1.3  : from : 30 apr 1994     **/
/**                                 to     18 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     31 oct 1994     **/
/**                # Version 3.0  : from : 07 jul 1995     **/
/**                                 to     28 sep 1995     **/
/**                # Version 3.1  : from : 28 nov 1995     **/
/**                                 to     08 jun 1996     **/
/**                # Version 3.2  : from : 07 sep 1996     **/
/**                                 to     17 sep 1998     **/
/**                # Version 4.0  : from : 13 dec 2001     **/
/**                                 to     31 aug 2005     **/
/**                # Version 5.0  : from : 13 dec 2006     **/
/**                                 to     24 mar 2008     **/
/**                # Version 5.1  : from : 30 oct 2009     **/
/**                                 to     30 oct 2009     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH_COARSEN

#include "module.h"
#include "common.h"
#include "graph.h"
#include "graph_coarsen.h"

/*
**  The static variables.
*/

static Gnum              (* graphCoarsenFuncTab[GRAPHCOARNBR]) (const Graph * const, Gnum *, const Gnum, const Gnum) = { /* Tables of edge-matching routines */
                              graphCoarsenMatchHy,
                              graphCoarsenMatchSc,
                              graphCoarsenMatchCs,
                              graphCoarsenMatchCh };

/***************************/
/*                         */
/* The coarsening routine. */
/*                         */
/***************************/

/* This routine coarsens the given "finegraph" into
** "coargraph", as long as the coarsening ratio remains
** below some threshold value and the coarsened graph
** is not too small.
** It returns:
** - 0  : if the graph has been coarsened.
** - 1  : if the graph could not be coarsened.
** - 2  : on error.
*/

int
graphCoarsen (
const Graph * restrict const          finegrafptr, /*+ Graph to coarsen                    +*/
Graph * restrict const                coargrafptr, /*+ Coarse graph to build               +*/
GraphCoarsenMulti * restrict * const  coarmultptr, /*+ Pointer to multinode table to build +*/
const Gnum                            coarnbr,    /*+ Minimum number of coarse vertices    +*/
const double                          coarrat,    /*+ Maximum contraction ratio            +*/
const GraphCoarsenType                coartype)   /*+ Edge matching type                   +*/
{
  Gnum                          coarhashnbr;      /* Size of the hash table                   */
  Gnum                          coarhashmsk;      /* Mask for access to hash table            */
  GraphCoarsenHash * restrict   coarhashtab;      /* Table of edges to other multinodes       */
  Gnum                          coarvertnbr;      /* Number of coarse vertices                */
  Gnum                          coarvertnum;      /* Number of current multinode vertex       */
  Gnum                          coarvertmax;      /* Maximum number of multinode vertices     */
  Gnum                          coarvelomax;      /* Maximum vertex weight allowed            */
  GraphCoarsenMulti * restrict  coarmulttax;      /* Multinode array                          */
  Gnum * restrict               finecoartax;      /* Based access to finecoartab              */
  Gnum                          finevertnum;      /* Number of currently selected fine vertex */
  size_t                        coarmultoftval;
  size_t                        coarvelooftval;
  size_t                        coaredgeoftval;
  size_t                        coaredlooftval;

#ifdef SCOTCH_DEBUG_GRAPH2
  if (coartype >= GRAPHCOARNBR) {
    errorPrint ("graphCoarsen: invalid parameter");
    return     (2);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */
#ifdef SCOTCH_DEBUG_GRAPH1
  if (coarrat < 0.5L)                             /* If impossible coarsening ratio wanted */
    return (1);                                   /* We will never succeed                 */
#endif /* SCOTCH_DEBUG_GRAPH1 */

  coarvertmax = (Gnum) ((double) finegrafptr->vertnbr * coarrat); /* Maximum number of coarse vertices */
  if (coarvertmax < coarnbr)                      /* If there will be too few vertices in graph        */
    return (1);                                   /* It is useless to go any further                   */

  if ((finecoartax = (Gnum *) memAlloc (finegrafptr->vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("graphCoarsen: out of memory (1)"); /* Allocate coarse graph uncoarsening array */
    return     (2);
  }
  memSet (finecoartax, ~0, finegrafptr->vertnbr * sizeof (Gnum));
  finecoartax -= finegrafptr->baseval;            /* Set based access to finecoartab */

  coarvelomax = (3 * finegrafptr->velosum) / (2 * coarnbr) + 1;

  coarvertnbr = graphCoarsenFuncTab[coartype] (finegrafptr, finecoartax, coarvertmax, coarvelomax); /* Call proper matching function */

  if (coarvertnbr >= coarvertmax) {               /* If coarsened graph too large */
    memFree (finecoartax + finegrafptr->baseval); /* Do not proceed any further   */
    return  (1);
  }

#ifdef SCOTCH_DEBUG_GRAPH2
  for (finevertnum = finegrafptr->baseval; finevertnum < finegrafptr->vertnnd; finevertnum ++) {
    if (finecoartax[finevertnum] <= ~0) {         /* If coarsening not aborted, this should not happen */
      errorPrint ("graphCoarsen: internal error (1)");
      memFree    (finecoartax + finegrafptr->baseval);
      return     (2);
    }
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  memSet (coargrafptr, 0, sizeof (Graph));        /* Initialize coarse graph */
  coargrafptr->flagval = GRAPHFREEVERT | GRAPHVERTGROUP | GRAPHEDGEGROUP;
  coargrafptr->baseval = finegrafptr->baseval;
  coargrafptr->vertnbr = coarvertnbr;
  coargrafptr->vertnnd = coarvertnbr + coargrafptr->baseval;
  coargrafptr->velosum = finegrafptr->velosum;    /* Keep load of finer graph */

  for (coarhashmsk = 31; coarhashmsk < finegrafptr->degrmax; coarhashmsk = coarhashmsk * 2 + 1) ;
  coarhashmsk = coarhashmsk * 4 + 3;
  coarhashnbr = coarhashmsk + 1;

  if (memAllocGroup ((void **) (void *)
                     &coargrafptr->verttax, (size_t) ((coarvertnbr + 1)    * sizeof (Gnum)),
                     &coargrafptr->velotax, (size_t) (coarvertnbr          * sizeof (Gnum)),
                     &coarmulttax,          (size_t) (coarvertnbr          * sizeof (GraphCoarsenMulti)),
                     &coargrafptr->edgetax, (size_t) (finegrafptr->edgenbr * sizeof (Gnum)), /* Pre-allocate space for edge arrays */ 
                     &coargrafptr->edlotax, (size_t) (finegrafptr->edgenbr * sizeof (Gnum)),
                     &coarhashtab,          (size_t) (coarhashnbr          * sizeof (GraphCoarsenHash)), NULL) == NULL) {
    errorPrint ("graphCoarsen: out of memory (2)"); /* Allocate coarser graph structure */
    memFree    (finecoartax + finegrafptr->baseval);
    return     (2);
  }
  coargrafptr->verttax -= coargrafptr->baseval;   /* Base coarse graph arrays */
  coargrafptr->velotax -= coargrafptr->baseval;
  coargrafptr->edgetax -= coargrafptr->baseval;
  coargrafptr->edlotax -= coargrafptr->baseval;
  coarmulttax          -= coargrafptr->baseval;

  for (finevertnum = finegrafptr->baseval, coarvertnum = coargrafptr->baseval; /* Finalize finecoartab array */
       finevertnum < finegrafptr->vertnnd; finevertnum ++) {
    Gnum                finematenum;              /* Number of current mate vertex */

    finematenum = finecoartax[finevertnum];       /* Get mate number                               */
    if (finematenum >= finevertnum) {             /* If mate has larger number                     */
      coarmulttax[coarvertnum].vertnum[0] = finevertnum; /* Build new multinode                    */
      coarmulttax[coarvertnum].vertnum[1] = finematenum; /* Second index always biggest            */
      finecoartax[finematenum] =                  /* Point to coarse vertex                        */
      finecoartax[finevertnum] = coarvertnum;     /* Always valid since coarvertnum <= finevertnum */
      coarvertnum ++;                             /* One more multinode created                    */
    }
  }
#ifdef SCOTCH_DEBUG_GRAPH2
  if ((coarvertnum - coargrafptr->baseval) != coarvertnbr) {
    errorPrint ("graphCoarsen: internal error (2)");
    graphFree  (coargrafptr);
    memFree    (finecoartax + finegrafptr->baseval);
    return     (2);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  if (finegrafptr->velotax != NULL) {             /* If fine graph is weighted */
    for (coarvertnum = coargrafptr->baseval; coarvertnum < coargrafptr->vertnnd; coarvertnum ++) {
      Gnum                coarveloval;

      coarveloval = finegrafptr->velotax[coarmulttax[coarvertnum].vertnum[0]];
      if (coarmulttax[coarvertnum].vertnum[0] != coarmulttax[coarvertnum].vertnum[1])
        coarveloval += finegrafptr->velotax[coarmulttax[coarvertnum].vertnum[1]];
      coargrafptr->velotax[coarvertnum] = coarveloval;
    }
  }
  else {                                          /* Fine graph is not weighted */
    for (coarvertnum = coargrafptr->baseval; coarvertnum < coargrafptr->vertnnd; coarvertnum ++)
      coargrafptr->velotax[coarvertnum] = (coarmulttax[coarvertnum].vertnum[0] != coarmulttax[coarvertnum].vertnum[1]) ? 2 : 1;
  }

  memSet (coarhashtab, ~0, coarhashnbr * sizeof (GraphCoarsenHash));
  if (finegrafptr->edlotax != NULL)               /* If edge loads available */
    graphCoarsenEdgeLl (finegrafptr, finecoartax, coarmulttax, coargrafptr, coarhashtab, coarhashmsk);
  else                                            /* Fine edges not weighted */
    graphCoarsenEdgeLu (finegrafptr, finecoartax, coarmulttax, coargrafptr, coarhashtab, coarhashmsk);
  coargrafptr->edgenbr = coargrafptr->verttax[coargrafptr->vertnnd] - coargrafptr->baseval; /* Set exact number of edges */

  memFree (finecoartax + finegrafptr->baseval);

  coarvelooftval = coargrafptr->velotax - coargrafptr->verttax;
  coarmultoftval = (Gnum *) coarmulttax - coargrafptr->verttax;
  coaredgeoftval = coargrafptr->edgetax - coargrafptr->verttax;
  coaredlooftval = coargrafptr->edlotax - coargrafptr->verttax;
  memReallocGroup ((void *) (coargrafptr->verttax + coargrafptr->baseval), /* Re-allocate data, wiping temporary arrays */
                   &coargrafptr->verttax, (size_t) ((coarvertnbr + 1)    * sizeof (Gnum)),
                   &coargrafptr->velotax, (size_t) (coarvertnbr          * sizeof (Gnum)),
                   &coarmulttax,          (size_t) (coarvertnbr          * sizeof (GraphCoarsenMulti)),
                   &coargrafptr->edgetax, (size_t) (finegrafptr->edgenbr * sizeof (Gnum)),
                   &coargrafptr->edlotax, (size_t) (coargrafptr->edgenbr * sizeof (Gnum)), NULL);
  coargrafptr->verttax -= coargrafptr->baseval;
  coargrafptr->vendtax  = coargrafptr->verttax + 1; /* Use compact representation of arrays */
  coargrafptr->velotax  = coargrafptr->verttax + coarvelooftval;
  coargrafptr->edgetax  = coargrafptr->verttax + coaredgeoftval;
  coargrafptr->edlotax  = coargrafptr->verttax + coaredlooftval;
  coarmulttax           = (GraphCoarsenMulti *) (coargrafptr->verttax + coarmultoftval);
  *coarmultptr          = coarmulttax;            /* Return pointer to multinode array */

#ifdef SCOTCH_DEBUG_GRAPH2
  if (graphCheck (coargrafptr) != 0) {            /* Check graph consistency */
    errorPrint ("graphCoarsen: inconsistent graph data");
    graphFree  (coargrafptr);
    return     (2);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  return (0);
}

/****************************************/
/*                                      */
/* The edge array building subroutines. */
/*                                      */
/****************************************/

#define GRAPHCOARSENEDGENAME        graphCoarsenEdgeLl
#define GRAPHCOARSENEDGEINIT        const Gnum * restrict const fineedlotax = finegrafptr->edlotax
#define GRAPHCOARSENEDGEEDLOINIT    coaredlotax[coaredgenum] = fineedlotax[fineedgenum]
#define GRAPHCOARSENEDGEEDLOADD     coaredlotax[coarhashtab[h].edgenum] += fineedlotax[fineedgenum]
#define GRAPHCOARSENEDGEEDLOSUB     coaredlosum -= finegrafptr->edlotax[fineedgenum]
#include "graph_coarsen_edge.c"
#undef GRAPHCOARSENEDGENAME
#undef GRAPHCOARSENEDGEINIT
#undef GRAPHCOARSENEDGEEDLOINIT
#undef GRAPHCOARSENEDGEEDLOADD
#undef GRAPHCOARSENEDGEEDLOSUB

#define GRAPHCOARSENEDGENAME        graphCoarsenEdgeLu
#define GRAPHCOARSENEDGEINIT
#define GRAPHCOARSENEDGEEDLOINIT    coaredlotax[coaredgenum] = 1
#define GRAPHCOARSENEDGEEDLOADD     coaredlotax[coarhashtab[h].edgenum] ++
#define GRAPHCOARSENEDGEEDLOSUB     coaredlosum --
#include "graph_coarsen_edge.c"
#undef GRAPHCOARSENEDGENAME
#undef GRAPHCOARSENEDGEINIT
#undef GRAPHCOARSENEDGEEDLOINIT
#undef GRAPHCOARSENEDGEEDLOADD
#undef GRAPHCOARSENEDGEEDLOSUB

/*****************************/
/*                           */
/* The matching subroutines. */
/*                           */
/*****************************/

static
Gnum
graphCoarsenMatchHy (
const Graph * restrict const  finegrafptr,        /* Fine graph to perform matching on */
Gnum * restrict               finecoartax,        /* Fine to coarse vertex index array */
const Gnum                    coarvertmax,        /* Maximum number of vertices to get */
const Gnum                    coarvelomax)        /* Maximum vertex weight allowed     */
{
  Gnum                  coarvertnum;              /* Number of current multinode vertex       */
  Gnum                  finepertbas;              /* Index of base of perturbation area       */
  Gnum                  finepertnbr;              /* Size of perturbation area                */
  const Gnum * restrict fineverttax;              /* Based access to vertex array             */
  const Gnum * restrict finevendtax;              /* Based access to end vertex array         */
  const Gnum * restrict finevelotax;              /* Based access to end vertex array         */
  const Gnum * restrict fineedgetax;              /* Based access to end vertex array         */
  const Gnum * restrict fineedlotax;              /* Based access to end vertex array         */
  Gnum                  finevertnnd;              /* Current end of vertex array              */
  Gnum                  finevertnum;              /* Number of currently selected fine vertex */

  if (finegrafptr->edlotax == NULL)               /* If no edge loads, perform scan matching instead */
    return (graphCoarsenMatchSc (finegrafptr, finecoartax, coarvertmax, coarvelomax));

  fineverttax = finegrafptr->verttax;
  finevendtax = finegrafptr->vendtax;
  finevelotax = finegrafptr->velotax;
  fineedgetax = finegrafptr->edgetax;
  fineedlotax = finegrafptr->edlotax;
  finevertnnd = finegrafptr->vertnnd;
  coarvertnum = 0;

  if (finegrafptr->velotax != NULL) {
    Gnum                finevelodlt;              /* Minimum load of neighbor */

    finevelodlt = (3 * finegrafptr->velosum) / (5 * (finevertnnd - finegrafptr->baseval));

    for (finevertnum = finegrafptr->baseval;      /* Pre-selection loop for isolated and lightest vertices */
         finevertnum < finevertnnd; finevertnum ++) {
      if (fineverttax[finevertnum] == finevendtax[finevertnum]) { /* If isolated vertex      */
        while (finecoartax[-- finevertnnd] != ~0) ; /* Search for first matchable "neighbor" */

        finecoartax[finevertnum] = finevertnnd;   /* At worst we will stop at finevertnum */
        finecoartax[finevertnnd] = finevertnum;
        coarvertnum ++;                           /* One more coarse vertex created */
      }
      else {                                      /* Vertex has neighbors */
        if ((finevelotax[finevertnum] < finevelodlt) &&
            (finecoartax[finevertnum] == ~0)) {   /* If vertex is too light on average  */
          Gnum                finevertbst;        /* Number of current best neighbor    */
          Gnum                fineedlobst;        /* Edge load of current best neighbor */
          Gnum                fineedgenum;

          if (coarvertnum >= coarvertmax)         /* If coarse graph is too large       */
            return (coarvertmax);                 /* Return that we cannot coarsen more */

          finevertbst = finevertnum;              /* No matching neighbor found yet */
          fineedlobst = 0;
          for (fineedgenum = fineverttax[finevertnum]; /* For all adjacent vertices */
               fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
            if ((finecoartax[fineedgetax[fineedgenum]] == ~0) && /* If unmatched vertex */
                (fineedlotax[fineedgenum] > fineedlobst)) { /* And is better candidate  */
              fineedlobst = fineedlotax[fineedgenum];
              finevertbst = fineedgetax[fineedgenum];
            }
          }

          finecoartax[finevertnum] = finevertbst;
          finecoartax[finevertbst] = finevertnum;
          coarvertnum ++;                         /* One more coarse vertex created */
        }
      }
    }
  }

  finepertnbr = 2 + intRandVal (GRAPHCOARPERTPRIME - 2); /* Compute perturbation area size (avoid DIV0 in random) */
  for (finepertbas = finegrafptr->baseval; finepertbas < finevertnnd; /* Run cache-friendly perturbation          */
       finepertbas += finepertnbr) {
    Gnum                finepertval;              /* Current index in perturbation area */

    if (finepertbas + finepertnbr > finevertnnd)
      finepertnbr = finevertnnd - finepertbas;

    finepertval = 0;                              /* Start from first perturbation vertex */
    do {                                          /* Loop on perturbation vertices        */
      finevertnum = finepertbas + finepertval;    /* Compute corresponding vertex number  */

      if (finecoartax[finevertnum] == ~0) {       /* If vertex has not been picked already */
        Gnum                finevertbst;          /* Number of current best neighbor       */
        Gnum                fineedlobst;          /* Edge load of current best neighbor    */
        Gnum                finevelodlt;          /* Maximum load of neighbor              */
        Gnum                fineedgenum;

        if (coarvertnum >= coarvertmax)           /* If coarse graph too large       */
          return (coarvertmax);                   /* Return that cannot coarsen more */

        finevertbst = finevertnum;                /* No matching vertex found yet */
        fineedlobst = 0;
        finevelodlt = coarvelomax - ((finevelotax != NULL) ? finevelotax[finevertnum] : 1);

        for (fineedgenum = fineverttax[finevertnum]; /* For all adjacent vertices */
             fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
          if ((finecoartax[fineedgetax[fineedgenum]] == ~0) && /* If unmatched vertex */
              (fineedlotax[fineedgenum] > fineedlobst) && /* And better candidate     */
              ((finevelotax == NULL) ||  /* And does not create overloads             */
               (finevelodlt >= finevelotax[fineedgetax[fineedgenum]]))) {
            fineedlobst = fineedlotax[fineedgenum];
            finevertbst = fineedgetax[fineedgenum];
          }
        }

        finecoartax[finevertnum] = finevertbst;
        finecoartax[finevertbst] = finevertnum;
        coarvertnum ++;                           /* One more coarse vertex created */
      }

      finepertval = (finepertval + GRAPHCOARPERTPRIME) % finepertnbr; /* Compute next perturbation index */
    } while (finepertval != 0);
  }

  return (coarvertnum);                           /* Return number of coarse vertices */
}

static
Gnum
graphCoarsenMatchSc (
const Graph * restrict const  finegrafptr,        /* Fine graph to perform matching on */
Gnum * restrict               finecoartax,        /* Fine to coarse vertex index array */
const Gnum                    coarvertmax,        /* Maximum number of vertices to get */
const Gnum                    coarvelomax)        /* Maximum vertex weight allowed     */
{
  Gnum                  coarvertnum;              /* Number of current multinode vertex       */
  Gnum                  finepertbas;              /* Index of base of perturbation area       */
  const Gnum * restrict fineverttax;              /* Based access to vertex array             */
  Gnum                  finepertnbr;              /* Size of perturbation area                */
  Gnum                  finevertnnd;              /* Current end of vertex array              */
  Gnum                  finevertnum;              /* Number of currently selected fine vertex */

  fineverttax = finegrafptr->verttax;

  coarvertnum = 0;
  for (finepertbas = finegrafptr->baseval, finevertnnd = finegrafptr->vertnnd;
       finepertbas < finevertnnd; finepertbas += finepertnbr) { /* Run cache-friendly perturbation */
    Gnum                finepertval;              /* Current index in perturbation area            */

    finepertnbr = finegrafptr->degrmax * 2 + intRandVal (finegrafptr->degrmax + 1) + 1; /* Compute perturbation area size (avoid DIV0 in random) */
    if (finepertnbr >= GRAPHCOARPERTPRIME)
      finepertnbr = 32 + intRandVal (GRAPHCOARPERTPRIME - 34);

    if (finepertbas + finepertnbr > finevertnnd)
      finepertnbr = finevertnnd - finepertbas;

    finepertval = 0;                              /* Start from first perturbation vertex */
    do {                                          /* Loop on perturbation vertices        */
      finevertnum = finepertbas + finepertval;    /* Compute corresponding vertex number  */

      if (finecoartax[finevertnum] == ~0) {       /* If vertex has not been picked already  */
        Gnum                finevertbst;          /* Number of current best matching vertex */

        if (coarvertnum >= coarvertmax)           /* If coarse graph is too large          */
          return (coarvertmax);                   /* Return that we cannot coarsen more    */

        if (fineverttax[finevertnum] == finegrafptr->vendtax[finevertnum]) { /* If isolated vertex */
          while (finecoartax[-- finevertnnd] != ~0) ; /* Search for first matchable "neighbor"     */
          finevertbst = finevertnnd;              /* Unmatched vertex will act as neighbor         */
        }
        else {                                    /* Vertex has neighbors */
          Gnum                finevelodlt;        /* Overload limit       */
          Gnum                fineedgenum;        /* Current edge number  */

          finevertbst = finevertnum;              /* No matching vertex found yet */
          finevelodlt = coarvelomax - ((finegrafptr->velotax != NULL) ? finegrafptr->velotax[finevertnum] : 1);

          for (fineedgenum = fineverttax[finevertnum]; /* For all adjacent vertices */
               fineedgenum < finegrafptr->vendtax[finevertnum]; fineedgenum ++) {
            if ((finecoartax[finegrafptr->edgetax[fineedgenum]] == ~0) && /* If unmatched vertex */
                ((finegrafptr->velotax == NULL) || /* And does not create overloads              */
                 (finevelodlt >= finegrafptr->velotax[finegrafptr->edgetax[fineedgenum]]))) {
              finevertbst = finegrafptr->edgetax[fineedgenum];
              break;
            }
          }
        }

        finecoartax[finevertnum] = finevertbst;
        finecoartax[finevertbst] = finevertnum;
        coarvertnum ++;                           /* One more coarse vertex created */
      }

      finepertval = (finepertval + GRAPHCOARPERTPRIME) % finepertnbr; /* Compute next perturbation index */
    } while (finepertval != 0);
  }

  return (coarvertnum);                           /* Return number of coarse vertices */
}

static
Gnum
graphCoarsenMatchCs (                             /* Crystallographic scan             */
const Graph * restrict const  finegrafptr,        /* Fine graph to perform matching on */
Gnum * restrict               finecoartax,        /* Fine to coarse vertex index array */
const Gnum                    coarvertmax,        /* Maximum number of vertices to get */
const Gnum                    coarvelomax)        /* Maximum vertex weight allowed     */
{
  Gnum                  coarvertnum;              /* Number of current multinode vertex       */
  const Gnum * restrict fineverttax;              /* Based access to vertex array             */
  const Gnum * restrict finevendtax;              /* Based access to end vertex array         */
  const Gnum * restrict finevelotax;              /* Based access to vertex load array        */
  const Gnum * restrict fineedgetax;              /* Based access to edge array               */
  Gnum                  finevertnum;              /* Number of currently selected fine vertex */
  Gnum * restrict       finequeutab;
  Gnum                  finequeuheadval;
  Gnum                  finequeutailval;
  Gnum                  finepermnum;              /* Permutation number for finding connected components */

  fineverttax = finegrafptr->verttax;
  finevendtax = finegrafptr->vendtax;
  finevelotax = finegrafptr->velotax;
  fineedgetax = finegrafptr->edgetax;
  if ((finequeutab = memAlloc (finegrafptr->vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("graphCoarsenMatchCs: out of memory");
    return     (graphCoarsenMatchSc (finegrafptr, finecoartax, coarvertmax, coarvelomax)); /* Fallback strategy */
  }

  coarvertnum     = 0;
  finequeuheadval = 1;
  finequeutailval = 0;
  finequeutab[0] = finegrafptr->baseval + intRandVal (finegrafptr->vertnbr); /* Start from random vertex */
  finecoartax[finequeutab[0]] = -2;               /* Set vertex as enqueued */
  
  for (finepermnum = finegrafptr->baseval; finequeutailval < finegrafptr->vertnbr; ) {
    if (finequeutailval < finequeuheadval) {      /* If vertices remain in queue */
      Gnum                finevertbst;            /* Best vertex found till now  */
      Gnum                finevelodlt;            /* Overload limit              */
      Gnum                fineedgenum;            /* Current edge number         */

      finevertnum = finequeutab[finequeutailval ++]; /* Select a vertex from the queue */

      if (finecoartax[finevertnum] >= 0) {        /* If selected vertex already matched */
        for (fineedgenum = fineverttax[finevertnum]; /* For all adjacent vertices       */
             fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
          Gnum                finevertend;

          finevertend = fineedgetax[fineedgenum];
          if (finecoartax[finevertend] == ~0) {
            finequeutab[finequeuheadval ++] = finevertend;
            finecoartax[finevertend] = -2;
          }
        }
        continue;                                 /* Skip to next vertex */
      }

      if (coarvertnum >= coarvertmax)             /* If coarse graph is too large       */
        break;                                    /* Return that we cannot coarsen more */

      finevertbst = finevertnum;                  /* No matching vertex found yet */
      finevelodlt = coarvelomax - ((finegrafptr->velotax != NULL) ? finegrafptr->velotax[finevertnum] : 1);

      for (fineedgenum = fineverttax[finevertnum]; /* For all adjacent vertices */
           fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
        Gnum                finevertend;
        Gnum                finecoarval;

        finevertend = fineedgetax[fineedgenum];
        finecoarval = finecoartax[finevertend];

        if (finecoarval < 0) {                    /* If vertex not matched  */
          if (finecoartax[finevertend] == ~0) {   /* If vertex not enqueued */
            finequeutab[finequeuheadval ++] = finevertend; /* Enqueue it    */
            finecoartax[finevertend] = -2;
          }
          if ((finevelotax == NULL) ||            /* And does not create overloads */
              (finevelodlt >= finevelotax[finevertend])) {
            finevertbst = finevertend;            /* Get matching vertex */

            while (++ fineedgenum < finevendtax[finevertnum]) { /* Scan and enqueue remaining neighbors */
              finevertend = fineedgetax[fineedgenum];
              if (finecoartax[finevertend] == ~0) {
                finequeutab[finequeuheadval ++] = finevertend;
                finecoartax[finevertend] = -2;
              }
            }
          }
        }
      }

      finecoartax[finevertnum] = finevertbst;     /* Match both vertices */
      finecoartax[finevertbst] = finevertnum;
      coarvertnum ++;                             /* One more coarse vertex created */
    }
    else {                                        /* Search for other connected component */
      Gnum                finevertbst;

      for ( ; finecoartax[finepermnum] >= 0; finepermnum ++) { /* Scan vertices in ascending order */
#ifdef SCOTCH_DEBUG_GRAPH2
        if (finepermnum >= finegrafptr->vertnnd) {
          errorPrint ("graphCoarsenMatchCs: internal error (1)");
          memFree    (finequeutab);
          return     (finegrafptr->vertnbr);      /* Coarsening aborted */
        }
#endif /* SCOTCH_DEBUG_GRAPH2 */
      }
#ifdef SCOTCH_DEBUG_GRAPH2
      if (finecoartax[finepermnum] != ~0) {
        errorPrint ("graphCoarsenMatchCs: internal error (2)");
        memFree    (finequeutab);
        return     (finegrafptr->vertnbr);        /* Coarsening aborted */
      }
#endif /* SCOTCH_DEBUG_GRAPH2 */
      finevertnum = finepermnum ++;               /* Start from found vertex */

      if (fineverttax[finevertnum] != finevendtax[finevertnum]) { /* If vertex not isolated */
        finequeutab[finequeuheadval ++] = finevertnum; /* Enqueue it for normal processing  */
        continue;                                 /* Skip to main loop to process it        */
      }

      finequeuheadval = ++ finequeutailval;       /* One more vertex enqueued-edqueued */

      if (coarvertnum >= coarvertmax)             /* If coarse graph is too large       */
        break;                                    /* Return that we cannot coarsen more */

      if (finequeutailval >= finegrafptr->vertnbr) /* If isolated vertex is last available vertex */
        finevertbst = finevertnum;
      else {
        for ( ; finecoartax[finepermnum] >= 0; finepermnum ++) {
#ifdef SCOTCH_DEBUG_GRAPH2
          if (finepermnum >= finegrafptr->vertnnd) {
            errorPrint ("graphCoarsenMatchCs: internal error (3)");
            memFree    (finequeutab);
            return     (finegrafptr->vertnbr);    /* Coarsening aborted */
          }
#endif /* SCOTCH_DEBUG_GRAPH2 */
        }
#ifdef SCOTCH_DEBUG_GRAPH2
        if (finecoartax[finepermnum] != ~0) {
          errorPrint ("graphCoarsenMatchCs: internal error (4)");
          memFree    (finequeutab);
          return     (finegrafptr->vertnbr);      /* Coarsening aborted */
        }
#endif /* SCOTCH_DEBUG_GRAPH2 */
        finevertbst     = finepermnum ++;         /* Get found vertex                  */
        finequeuheadval = ++ finequeutailval;     /* One more vertex enqueued-edqueued */
      }

      finecoartax[finevertnum] = finevertbst;     /* Match both vertices */
      finecoartax[finevertbst] = finevertnum;
      coarvertnum ++;                             /* One more coarse vertex created */
    }
  }

  memFree (finequeutab);

  return (coarvertnum);                           /* Return number of coarse vertices */
}

static
Gnum
graphCoarsenMatchCh (                             /* Crystallographic heavy edge       */
const Graph * restrict const  finegrafptr,        /* Fine graph to perform matching on */
Gnum * restrict               finecoartax,        /* Fine to coarse vertex index array */
const Gnum                    coarvertmax,        /* Maximum number of vertices to get */
const Gnum                    coarvelomax)        /* Maximum vertex weight allowed     */
{
  Gnum                  coarvertnum;              /* Number of current multinode vertex       */
  const Gnum * restrict fineverttax;              /* Based access to vertex array             */
  const Gnum * restrict finevendtax;              /* Based access to end vertex array         */
  const Gnum * restrict finevelotax;              /* Based access to vertex load array        */
  const Gnum * restrict fineedgetax;              /* Based access to edge array               */
  const Gnum * restrict fineedlotax;              /* Based access to edge load array          */
  Gnum                  finevertnum;              /* Number of currently selected fine vertex */
  Gnum * restrict       finequeutab;
  Gnum                  finequeuheadval;
  Gnum                  finequeutailval;
  Gnum                  finepermnum;              /* Permutation number for finding connected components */

  if (finegrafptr->edlotax == NULL)               /* If no edge loads, perform scan matching instead */
    return (graphCoarsenMatchCs (finegrafptr, finecoartax, coarvertmax, coarvelomax));

  fineverttax = finegrafptr->verttax;
  finevendtax = finegrafptr->vendtax;
  finevelotax = finegrafptr->velotax;
  fineedgetax = finegrafptr->edgetax;
  fineedlotax = finegrafptr->edlotax;
  if ((finequeutab = memAlloc (finegrafptr->vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("graphCoarsenMatchCh: out of memory");
    return     (graphCoarsenMatchSc (finegrafptr, finecoartax, coarvertmax, coarvelomax)); /* Fallback strategy */
  }

  coarvertnum     = 0;
  finequeuheadval = 1;
  finequeutailval = 0;
  finequeutab[0] = finegrafptr->baseval + intRandVal (finegrafptr->vertnbr); /* Start from random vertex */
  finecoartax[finequeutab[0]] = -2;               /* Set vertex as enqueued */
  
  for (finepermnum = finegrafptr->baseval; finequeutailval < finegrafptr->vertnbr; ) {
    if (finequeutailval < finequeuheadval) {      /* If vertices remain in queue        */
      Gnum                finevertbst;            /* Best vertex found till now         */
      Gnum                fineedlobst;            /* Edge load of current best neighbor */
      Gnum                finevelodlt;            /* Overload limit                     */
      Gnum                fineedgenum;            /* Current edge number                */

      finevertnum = finequeutab[finequeutailval ++]; /* Select a vertex from the queue */

      if (finecoartax[finevertnum] >= 0) {        /* If selected vertex already matched */
        for (fineedgenum = fineverttax[finevertnum]; /* For all adjacent vertices       */
             fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
          Gnum                finevertend;

          finevertend = fineedgetax[fineedgenum];
          if (finecoartax[finevertend] == ~0) {
            finequeutab[finequeuheadval ++] = finevertend;
            finecoartax[finevertend] = -2;
          }
        }
        continue;                                 /* Skip to next vertex */
      }

      if (coarvertnum >= coarvertmax)             /* If coarse graph is too large       */
        break;                                    /* Return that we cannot coarsen more */

      finevertbst = finevertnum;                  /* No matching vertex found yet */
      fineedlobst = 0;
      finevelodlt = coarvelomax - ((finegrafptr->velotax != NULL) ? finegrafptr->velotax[finevertnum] : 1);

      for (fineedgenum = fineverttax[finevertnum]; /* For all adjacent vertices */
           fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
        Gnum                finevertend;
        Gnum                finecoarval;

        finevertend = fineedgetax[fineedgenum];
        finecoarval = finecoartax[finevertend];

        if (finecoarval < 0) {                    /* If vertex not matched */
          Gnum                fineedloval;

          fineedloval = fineedlotax[fineedgenum];
          if (finecoartax[finevertend] == ~0) {   /* If vertex not enqueued */
            finequeutab[finequeuheadval ++] = finevertend; /* Enqueue it    */
            finecoartax[finevertend] = -2;
          }
          if (((finevelotax == NULL) ||            /* And does not create overloads */
               (finevelodlt >= finevelotax[finevertend])) &&
              (fineedloval > fineedlobst)) {
            finevertbst = finevertend;            /* Get matching vertex */
            fineedlobst = fineedloval;
          }
        }
      }

      finecoartax[finevertnum] = finevertbst;     /* Match both vertices */
      finecoartax[finevertbst] = finevertnum;
      coarvertnum ++;                             /* One more coarse vertex created */
    }
    else {                                        /* Search for other connected component */
      Gnum                finevertbst;

      for ( ; finecoartax[finepermnum] >= 0; finepermnum ++) { /* Scan vertices in ascending order */
#ifdef SCOTCH_DEBUG_GRAPH2
        if (finepermnum >= finegrafptr->vertnnd) {
          errorPrint ("graphCoarsenMatchCh: internal error (1)");
          memFree    (finequeutab);
          return     (finegrafptr->vertnbr);      /* Coarsening aborted */
        }
#endif /* SCOTCH_DEBUG_GRAPH2 */
      }
#ifdef SCOTCH_DEBUG_GRAPH2
      if (finecoartax[finepermnum] != ~0) {
        errorPrint ("graphCoarsenMatchCh: internal error (2)");
        memFree    (finequeutab);
        return     (finegrafptr->vertnbr);        /* Coarsening aborted */
      }
#endif /* SCOTCH_DEBUG_GRAPH2 */
      finevertnum = finepermnum ++;               /* Start from found vertex */

      if (fineverttax[finevertnum] != finevendtax[finevertnum]) { /* If vertex not isolated */
        finequeutab[finequeuheadval ++] = finevertnum; /* Enqueue it for normal processing  */
        continue;                                 /* Skip to main loop to process it        */
      }

      finequeuheadval = ++ finequeutailval;       /* One more vertex enqueued-edqueued */

      if (coarvertnum >= coarvertmax)             /* If coarse graph is too large       */
        break;                                    /* Return that we cannot coarsen more */

      if (finequeutailval >= finegrafptr->vertnbr) /* If isolated vertex is last available vertex */
        finevertbst = finevertnum;
      else {
        for ( ; finecoartax[finepermnum] >= 0; finepermnum ++) {
#ifdef SCOTCH_DEBUG_GRAPH2
          if (finepermnum >= finegrafptr->vertnnd) {
            errorPrint ("graphCoarsenMatchCh: internal error (3)");
            memFree    (finequeutab);
            return     (finegrafptr->vertnbr);    /* Coarsening aborted */
          }
#endif /* SCOTCH_DEBUG_GRAPH2 */
        }
#ifdef SCOTCH_DEBUG_GRAPH2
        if (finecoartax[finepermnum] != ~0) {
          errorPrint ("graphCoarsenMatchCh: internal error (4)");
          memFree    (finequeutab);
          return     (finegrafptr->vertnbr);      /* Coarsening aborted */
        }
#endif /* SCOTCH_DEBUG_GRAPH2 */
        finevertbst     = finepermnum ++;         /* Get found vertex                  */
        finequeuheadval = ++ finequeutailval;     /* One more vertex enqueued-edqueued */
      }

      finecoartax[finevertnum] = finevertbst;     /* Match both vertices */
      finecoartax[finevertbst] = finevertnum;
      coarvertnum ++;                             /* One more coarse vertex created */
    }
  }

  memFree (finequeutab);

  return (coarvertnum);                           /* Return number of coarse vertices */
}
