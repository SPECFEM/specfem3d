/* Copyright 2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : kgraph_map_ml.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module performs the Dual Recursive **/
/**                Bipartitioning mapping algorithm.       **/
/**                It is now a branching routine.          **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 13 jul 2010     **/
/**                                 to     14 jul 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_ML

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "graph_coarsen.h"
#include "arch.h"
#include "mapping.h"
#include "kgraph.h"
#include "kgraph_map_ml.h"
#include "kgraph_map_st.h"

/*********************************************/
/*                                           */
/* The coarsening and uncoarsening routines. */
/*                                           */
/*********************************************/

/* This routine builds a coarser graph from the
** graph that is given on input. The coarser
** graphs differ at this stage from classical
** active graphs as their internal gains are not
** yet computed.
** It returns:
** - 0  : if the coarse graph has been built.
** - 1  : if threshold reached or on error.
*/

static
int
kgraphMapMlCoarsen (
const Kgraph * restrict const         finegrafptr, /*+ Finer graph                         +*/
Kgraph * restrict const               coargrafptr, /*+ Coarser graph to build              +*/
GraphCoarsenMulti * restrict * const  coarmultptr, /*+ Pointer to multinode table to build +*/
const KgraphMapMlParam * const        paraptr)    /*+ Method parameters                    +*/
{
  Gnum                vertmax;

  vertmax = (archVar (&finegrafptr->m.archdat) ? 1 : finegrafptr->m.domnnbr) * paraptr->coarnbr;

  if (graphCoarsen (&finegrafptr->s, &coargrafptr->s, coarmultptr, vertmax,
                    paraptr->coarrat, paraptr->coartype) != 0)
    return (1);                                   /* Return if coarsening failed */

  coargrafptr->m.baseval   = finegrafptr->m.baseval;
  coargrafptr->m.vertnbr   = coargrafptr->s.vertnbr;
  coargrafptr->m.parttax   = NULL;                /* Do not allocate part array yet */
  coargrafptr->m.domntab   = finegrafptr->m.domntab; /* Re-use domain array         */
  coargrafptr->m.domnnbr   = finegrafptr->m.domnnbr;
  coargrafptr->m.domnmax   = finegrafptr->m.domnmax;
  coargrafptr->m.archdat   = finegrafptr->m.archdat;
  coargrafptr->m.domnorg   = finegrafptr->m.domnorg;
  coargrafptr->frontab     = finegrafptr->frontab; /* Re-use frontier array if it exists */
  coargrafptr->comploadavg = finegrafptr->comploadavg; /* Re-use existing load arrays    */
  coargrafptr->comploaddlt = finegrafptr->comploaddlt;
  coargrafptr->levlnum     = finegrafptr->levlnum + 1; /* Graph level is coarsening level */

  return (0);
}

/* This routine propagates the separation of the
** coarser graph back to the finer graph, according
** to the multinode table of collapsed vertices.
** After the separation is propagated, it finishes
** to compute the parameters of the finer graph that
** were not computed at the coarsening stage.
** It returns:
** - 0   : if coarse graph data has been propagated to fine graph.
** - !0  : on error.
*/

static
int
kgraphMapMlUncoarsen (
Kgraph * restrict const                   finegrafptr, /*+ Finer graph     +*/
Kgraph * restrict const                   coargrafptr, /*+ Coarser graph   +*/
const GraphCoarsenMulti * restrict const  coarmulttax) /*+ Multinode array +*/
{
  Gnum                          coarvertnum;
  const Anum * restrict         coarparttax;
  Gnum *                        coarfrontab;      /* [norestrict] */
  Gnum                          coarfronnum;
  Anum * restrict               fineparttax;
  Gnum *                        finefrontab;      /* [norestrict] */
  Gnum                          finefronnbr;

  const Gnum * restrict const fineverttax = finegrafptr->s.verttax;
  const Gnum * restrict const finevendtax = finegrafptr->s.vendtax;
  const Gnum * restrict const fineedgetax = finegrafptr->s.edgetax;

  if (finegrafptr->m.parttax == NULL) {           /* If partition array not yet allocated */
    if ((finegrafptr->m.parttax = (Anum *) memAlloc (finegrafptr->s.vertnbr * sizeof (Anum))) == NULL) {
      errorPrint ("kgraphMapMlUncoarsen: out of memory (1)");
      return     (1);
    }
    finegrafptr->s.flagval |= KGRAPHFREEPART;     /* Allocated data will be freed along with graph structure */
    finegrafptr->m.parttax -= finegrafptr->s.baseval;
  }
  if (finegrafptr->frontab == NULL) {             /* If frontier array not yet allocated */
    if ((finegrafptr->frontab = (Gnum *) memAlloc (finegrafptr->s.vertnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("kgraphMapMlUncoarsen: out of memory (2)");
      return     (1);
    }
  }

  if (coargrafptr == NULL) {                      /* If no coarse graph provided            */
    kgraphFrst (finegrafptr);                     /* Assign all vertices to first subdomain */
    return     (0);
  }

  finegrafptr->m.domntab = coargrafptr->m.domntab; /* Get pointer to domain array again in case it was reallocated */
  finegrafptr->m.domnnbr = coargrafptr->m.domnnbr;
  finegrafptr->m.domnmax = coargrafptr->m.domnmax;
  coargrafptr->m.domntab = NULL;                  /* Do not doubly free array */

  finegrafptr->comploadavg = coargrafptr->comploadavg; /* Get pointer to load array again in case it was reallocated */
  finegrafptr->comploaddlt = coargrafptr->comploaddlt;
  coargrafptr->comploadavg = NULL;                /* Do not doubly free array */

  coarparttax = coargrafptr->m.parttax;
  fineparttax = finegrafptr->m.parttax;
  for (coarvertnum = coargrafptr->s.baseval; coarvertnum < coargrafptr->s.vertnnd; coarvertnum ++) {
    Gnum                finevertnum0;             /* First multinode vertex  */
    Gnum                finevertnum1;             /* Second multinode vertex */
    Anum                partval;

    finevertnum0 = coarmulttax[coarvertnum].vertnum[0];
    finevertnum1 = coarmulttax[coarvertnum].vertnum[1];
    partval      = coarparttax[coarvertnum];

    fineparttax[finevertnum0] = partval;
    if (finevertnum0 != finevertnum1)
      fineparttax[finevertnum1] = partval;
  }

  coarfrontab = coargrafptr->frontab;             /* TRICK: may also be equal to finefrontab */
  finefrontab = finegrafptr->frontab;
  for (coarfronnum = 0, finefronnbr = coargrafptr->fronnbr; /* Re-cycle frontier array from coarse to fine graph */
       coarfronnum < coargrafptr->fronnbr; coarfronnum ++) {
    Gnum                coarvertnum;
    Gnum                finevertnum0;             /* First multinode vertex  */
    Gnum                finevertnum1;             /* Second multinode vertex */

    coarvertnum  = coarfrontab[coarfronnum];
    finevertnum0 = coarmulttax[coarvertnum].vertnum[0];
    finevertnum1 = coarmulttax[coarvertnum].vertnum[1];
      
    if (finevertnum0 != finevertnum1) {           /* If multinode si made of two distinct vertices */
      Anum                coarpartval;
      Gnum                fineedgenum;

      coarpartval = coarparttax[coarvertnum];

#ifdef SCOTCH_DEBUG_KGRAPH2
      finefrontab[coarfronnum] = ~0;
#endif /* SCOTCH_DEBUG_KGRAPH2 */

      for (fineedgenum = fineverttax[finevertnum0];
           fineedgenum < finevendtax[finevertnum0]; fineedgenum ++) {
        if (fineparttax[fineedgetax[fineedgenum]] != coarpartval) { /* If first vertex belongs to frontier */
          finefrontab[coarfronnum] = finevertnum0; /* Record it in lieu of the coarse frontier vertex      */
          break;
        }
      }
      if (fineedgenum >= finegrafptr->s.vendtax[finevertnum0]) { /* If first vertex not in frontier */
        finefrontab[coarfronnum] = finevertnum1;  /* Then second vertex must be in frontier         */
        continue;                                 /* Skip to next multinode                         */
      }

      for (fineedgenum = fineverttax[finevertnum1]; /* Check if second vertex belong to frontier too */
           fineedgenum < finevendtax[finevertnum1]; fineedgenum ++) {
        if (fineparttax[fineedgetax[fineedgenum]] != coarpartval) { /* If second vertex belongs to frontier      */
          finefrontab[finefronnbr ++] = finevertnum1; /* Record it at the end of the (recycled ?) frontier array */
          break;
        }
      }

#ifdef SCOTCH_DEBUG_KGRAPH2
      if (finefrontab[coarfronnum] == ~0) {
        errorPrint ("kgraphMapMlUncoarsen: internal error");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    }
    else                                          /* If coarse vertex is single node */
      finefrontab[coarfronnum] = finevertnum0;    /* Then it belongs to the frontier */
  }
  finegrafptr->fronnbr = finefronnbr;

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (kgraphCheck (finegrafptr) != 0) {
    errorPrint ("kgraphMapMlUncoarsen: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  return (0);
}

/* This routine recursively performs the partitioning.
** It returns:
** - 0   : if mapping could be computed.
** - !0  : on error.
*/

static
int
kgraphMapMl2 (
Kgraph * restrict const         grafptr,
const KgraphMapMlParam * const  paraptr)
{
  Kgraph                        coargrafdat;
  GraphCoarsenMulti * restrict  coarmulttax;
  int                           o;

  if (kgraphMapMlCoarsen (grafptr, &coargrafdat, &coarmulttax, paraptr) == 0) {
    if (((o = kgraphMapMl2         (&coargrafdat, paraptr)) == 0)              &&
        ((o = kgraphMapMlUncoarsen (grafptr, &coargrafdat, coarmulttax)) == 0) &&
	((o = kgraphMapSt          (grafptr, paraptr->stratasc)) != 0)) /* Apply ascending strategy */
      errorPrint ("kgraphMapMl2: cannot apply ascending strategy");
    kgraphExit (&coargrafdat);
  }
  else {                                          /* Cannot coarsen due to lack of memory or error */
    if (((o = kgraphMapMlUncoarsen (grafptr, NULL, NULL)) == 0) && /* Finalize graph               */
        ((o = kgraphMapSt          (grafptr, paraptr->stratlow)) != 0)) /* Apply low strategy      */
      errorPrint ("kgraphMapMl2: cannot apply low strategy");
  }

  return (o);
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the muti-level mapping.
** It returns:
** - 0 : if separator could be computed.
** - 1 : on error.
*/

int
kgraphMapMl (
Kgraph * const                  grafptr,          /*+ Graph to map      +*/
const KgraphMapMlParam * const  paraptr)          /*+ Method parameters +*/
{
  Gnum                levlnum;                    /* Save value for graph level */
  int                 o;

  levlnum = grafptr->levlnum;                     /* Save graph level            */
  grafptr->levlnum = 0;                           /* Initialize coarsening level */
  o = kgraphMapMl2 (grafptr, paraptr);            /* Perform multi-level mapping */
  grafptr->levlnum = levlnum;                     /* Restore graph level         */

  return (o);
}
