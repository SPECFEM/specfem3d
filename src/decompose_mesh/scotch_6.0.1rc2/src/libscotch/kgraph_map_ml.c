/* Copyright 2010,2011,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module maps an active graph        **/
/**                to a specific architecture graph        **/
/**                using a multi-level scheme.             **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 13 jul 2010     **/
/**                                 to     14 jul 2010     **/
/**   DATES      : # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     09 oct 2012     **/
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
#include "arch.h"
#include "mapping.h"
#include "graph_coarsen.h" 
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
  GraphCoarsenMulti * restrict  coarmulttab;
  Gnum                          coarvertnum;      /* Number of current multinode vertex */

  const Anum * restrict const finepfixtax = finegrafptr->pfixtax;

  if (graphCoarsen (&finegrafptr->s, &coargrafptr->s, coarmultptr, paraptr->coarnbr, paraptr->coarval, finegrafptr->r.m.parttax, finegrafptr->pfixtax, finegrafptr->vfixnbr, &coargrafptr->vfixnbr) != 0)
    return (1);                                   /* Return if coarsening failed */

  coargrafptr->frontab    = finegrafptr->frontab; /* Use frontier array of finer graph as coarse frontier array */
  coargrafptr->a          = finegrafptr->a;
  coargrafptr->s.flagval &= ~KGRAPHFREEFRON;      /* Be sure to teep top frontab for all levels */
  coargrafptr->m.parttax  = NULL;                 /* Do not allocate partition data yet         */
  coargrafptr->m.domntab  = NULL;                 /* Do not allocate domain data yet            */
  coargrafptr->m.archptr  = &coargrafptr->a;
  coargrafptr->m.grafptr  = &coargrafptr->s;
  coargrafptr->m.flagval  = finegrafptr->m.flagval;
  coargrafptr->m.domnorg  = finegrafptr->m.domnorg;
  coargrafptr->m.domnnbr  = 0;                    /* Number of domains not yet known */
  coargrafptr->m.domnmax  = finegrafptr->m.domnmax; /* Propagate relevant estimation */

  coarmulttab = *coarmultptr;

  if (finegrafptr->r.m.parttax != NULL) {
    Gnum *                fineparotax;
    Gnum *                coarparotab;
    Gnum *                coarvmlotab;
    const Gnum * restrict finevmlotax; 
    Gnum                  coarvertnbr;
 
    coargrafptr->r.m = finegrafptr->r.m;          /* Clone old mapping; never free */

    coarvertnbr = coargrafptr->s.vertnbr;
    if ((coarparotab = (Anum *) memAlloc (coarvertnbr * sizeof (Anum))) == NULL) {
      errorPrint ("kgraphMapMlCoarsen: out of memory (1)");
      return     (1);
    }

    if ((coarvmlotab = (Gnum *) memAlloc (coarvertnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("kgraphMapMlCoarsen: out of memory (2)");
      return     (1);
    }
    coargrafptr->s.flagval |= KGRAPHFREEVMLO;

    fineparotax = finegrafptr->r.m.parttax;

    finevmlotax = finegrafptr->r.vmlotax;
    for (coarvertnum = 0; coarvertnum < coarvertnbr; coarvertnum ++) {
      Gnum                finevertnum0;
      Gnum                finevertnum1;

      finevertnum0 = coarmulttab[coarvertnum].vertnum[0];
      finevertnum1 = coarmulttab[coarvertnum].vertnum[1];
      coarvmlotab[coarvertnum] = finevmlotax[finevertnum0] + ((finevertnum0 == finevertnum1) ? 0 : finevmlotax[finevertnum1]);
      coarparotab[coarvertnum] = fineparotax[finevertnum0];
#ifdef SCOTCH_DEBUG_KGRAPH2
      if (((fineparotax[finevertnum1] != fineparotax[finevertnum0]) && /* Vertices were not in the same part */
          ((finegrafptr->pfixtax == NULL) ||  
          ((finegrafptr->pfixtax[finevertnum1] == -1) &&               /* And they are not fixed */
           (finegrafptr->pfixtax[finevertnum0] == -1)))) ||
          ((finegrafptr->pfixtax != NULL) &&
          ((finegrafptr->pfixtax[finevertnum1] !=                      /* Or they are fixed, but not in the same part */
            finegrafptr->pfixtax[finevertnum0])))) {
        errorPrint ("kgraphMapMlCoarsen: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    }

    coargrafptr->r.m.parttax = coarparotab - finegrafptr->s.baseval;
    coargrafptr->r.vmlotax   = coarvmlotab - finegrafptr->s.baseval;
  }
  else {
    coargrafptr->r.m.parttax = NULL;
    coargrafptr->r.vmlotax   = NULL;
  }

  if (finepfixtax != NULL) {                      /* If we have fixed vertices */  
    Gnum                coarvertnbr;
    Anum * restrict     coarpfixtab;
    Gnum                coarvfixnbr;

    coarvertnbr = coargrafptr->s.vertnbr;
    if ((coarpfixtab = (Anum *) memAlloc (coarvertnbr * sizeof (Anum))) == NULL) {
      errorPrint ("kgraphMapMlCoarsen: out of memory (3)");
      return     (1);
    }
    for (coarvertnum = 0, coarvfixnbr = 0;
         coarvertnum < coarvertnbr; coarvertnum ++) {
      coarpfixtab[coarvertnum] = finepfixtax[coarmulttab[coarvertnum].vertnum[0]];
      if (coarpfixtab[coarvertnum] >= 0)
        coarvfixnbr ++;
    }

    coargrafptr->pfixtax = coarpfixtab - coargrafptr->s.baseval;
    coargrafptr->vfixnbr = coarvfixnbr;
  }
  else {
    coargrafptr->pfixtax = NULL;
    coargrafptr->vfixnbr = 0;
  }

#ifdef SCOTCH_DEBUG_KGRAPH2
  if ((finegrafptr->comploadavg == NULL) || (finegrafptr->comploaddlt == NULL)) {
    errorPrint ("kgraphMapMlCoarsen: internal error (3)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  coargrafptr->comploadavg = finegrafptr->comploadavg; /* Use target average loads array of finer graph as coarse one */  
  coargrafptr->comploaddlt = finegrafptr->comploaddlt; /* Use target imbalances array of finer graph as coarse one    */
  coargrafptr->comploadrat = finegrafptr->comploadrat;
  coargrafptr->r.cmloval   = finegrafptr->r.cmloval;
  coargrafptr->r.crloval   = finegrafptr->r.crloval;
  coargrafptr->kbalval     = finegrafptr->kbalval;
  coargrafptr->levlnum     = finegrafptr->levlnum + 1;

  return (0);
}

/* This routine propagates the partitioning of the
** coarser graph back to the finer graph, according
** to the multinode table of collapsed vertices.
** After the partitioning is propagated, it finishes
** to compute the parameters of the finer graph that
** were not computed at the coarsening stage.
** It returns:
** - 0   : if coarse graph data has been propagated to fine graph.
** - !0  : on error.
*/

static
int
kgraphMapMlUncoarsen (
Kgraph * restrict const         finegrafptr,      /*+ Finer graph                +*/
const Kgraph * const            coargrafptr,      /*+ Coarser graph              +*/
const GraphCoarsenMulti * const coarmulttax)      /*+ Pointer to multinode array +*/
{
  Anum * restrict       fineparttax;
  Gnum                  finefronnum;
  const Anum * restrict coarparttax;
  Gnum                  coarvertnum;
  Gnum                  coarfronnum;
  Gnum * restrict       coarfrontab;

  const Gnum * restrict const fineverttax = finegrafptr->s.verttax; /* Fast accesses */
  const Gnum * restrict const finevendtax = finegrafptr->s.vendtax;
  const Gnum * restrict const fineedgetax = finegrafptr->s.edgetax;

  if (finegrafptr->m.parttax == NULL) {             /* If partition array not yet allocated */
    if ((finegrafptr->m.parttax = (Anum *) memAlloc (finegrafptr->s.vertnbr * sizeof (Anum))) == NULL) {
      errorPrint ("kgraphMapMlUncoarsen: out of memory");
      return     (1);                               /* Allocated data will be freed along with graph structure */
    }
    finegrafptr->m.parttax -= finegrafptr->s.baseval;
  }
  if (finegrafptr->m.domntab == NULL) {
    if ((finegrafptr->m.domntab = (ArchDom *) memAlloc (finegrafptr->m.domnmax * sizeof (ArchDom))) == NULL) {
      errorPrint ("kdgraphMapMlUncoarsen: out of memory (2)");
      return     (1);
    }
    finegrafptr->m.flagval |= MAPPINGFREEDOMN;
  }

  if (coargrafptr == NULL) {                      /* If no coarse graph provided                */
    kgraphFrst (finegrafptr);                     /* Assign all vertices to the first subdomain */
    return     (0);
  }
  else {
    finegrafptr->m.domnnbr = coargrafptr->m.domnnbr;
    memCpy (finegrafptr->m.domntab, coargrafptr->m.domntab, finegrafptr->m.domnnbr * sizeof (ArchDom)); 
  }

  coarparttax = coargrafptr->m.parttax;
  coarfrontab = coargrafptr->frontab;
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

  finegrafptr->commload = coargrafptr->commload;

  for (coarfronnum = 0, finefronnum = coargrafptr->fronnbr; /* Re-cycle frontier array from coarse to fine graph */
       coarfronnum < coargrafptr->fronnbr; coarfronnum ++) {
    Gnum                coarvertnum;
    Gnum                finevertnum0;             /* First multinode vertex  */
    Gnum                finevertnum1;             /* Second multinode vertex */

    coarvertnum  = coarfrontab[coarfronnum];
    finevertnum0 = coarmulttax[coarvertnum].vertnum[0];
    finevertnum1 = coarmulttax[coarvertnum].vertnum[1];
      
    if (finevertnum0 != finevertnum1) {           /* If multinode si made of two distinct vertices */
      Gnum                coarpartval;
      Gnum                fineedgenum;

      coarpartval = coarparttax[coarvertnum];

#ifdef SCOTCH_DEBUG_KGRAPH2
      coarfrontab[coarfronnum] = ~0;
#endif /* SCOTCH_DEBUG_KGRAPH2 */

      for (fineedgenum = fineverttax[finevertnum0];
           fineedgenum < finevendtax[finevertnum0]; fineedgenum ++) {
        if (fineparttax[fineedgetax[fineedgenum]] != coarpartval) { /* If first vertex belongs to frontier */
          coarfrontab[coarfronnum] = finevertnum0; /* Record it in lieu of the coarse frontier vertex      */
          break;
        }
      }
      if (fineedgenum >= finegrafptr->s.vendtax[finevertnum0]) { /* If first vertex not in frontier */
        coarfrontab[coarfronnum] = finevertnum1;  /* Then second vertex must be in frontier         */
        continue;                                 /* Skip to next multinode                         */
      }

      for (fineedgenum = fineverttax[finevertnum1]; /* Check if second vertex belong to frontier too */
           fineedgenum < finevendtax[finevertnum1]; fineedgenum ++) {
        if (fineparttax[fineedgetax[fineedgenum]] != coarpartval) { /* If second vertex belongs to frontier  */
          coarfrontab[finefronnum ++] = finevertnum1; /* Record it at the end of the recycled frontier array */
          break;
        }
      }

#ifdef SCOTCH_DEBUG_KGRAPH2
      if (coarfrontab[coarfronnum] == ~0) {
        errorPrint ("kgraphMapMlUncoarsen: internal error");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    }
    else                                          /* If coarse vertex is single node */
      coarfrontab[coarfronnum] = finevertnum0;    /* Then it belongs to the frontier */
  }
  finegrafptr->fronnbr = finefronnum;

  if (coargrafptr->r.m.parttax != NULL)
    memFree (coargrafptr->r.m.parttax + coargrafptr->s.baseval);

  if ((finegrafptr->r.m.parttax != NULL) && (coargrafptr->r.m.parttax == NULL)) /* If repartitioning data has been removed */
    finegrafptr->r.m.parttax = NULL;              /* Propagate this information */

  if (coargrafptr->pfixtax != NULL)
    memFree (coargrafptr->pfixtax + coargrafptr->s.baseval);
  memFree (coarparttax + coargrafptr->s.baseval);

  kgraphCost (finegrafptr);
#ifdef SCOTCH_DEBUG_KGRAPH2
  if (kgraphCheck (finegrafptr) != 0) {
    errorPrint ("kgraphMapMlUncoarsen: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  return (0);
}

/* This routine performs the
** partitioning recursion.
** It returns:
** - 0 : if partitioning could be computed.
** - 1 : on error.
*/

static
int
kgraphMapMl2 (
Kgraph * restrict const           grafptr,        /*+ Active graph      +*/
const KgraphMapMlParam * const    paraptr)        /*+ Method parameters +*/
{
  Kgraph              coargrafdat;
  GraphCoarsenMulti * coarmultptr;
  int                 o;

  if (kgraphMapMlCoarsen (grafptr, &coargrafdat, &coarmultptr, paraptr) == 0) {
    if (((o = kgraphMapMl2         (&coargrafdat, paraptr))              == 0) &&
        ((o = kgraphMapMlUncoarsen (grafptr, &coargrafdat, coarmultptr)) == 0) &&
        ((o = kgraphMapSt          (grafptr, paraptr->stratasc))         != 0)) /* Apply ascending strategy */
      errorPrint ("kgraphMapMl2: cannot apply ascending strategy");
    kgraphExit (&coargrafdat);
  }
  else {                                          /* Cannot coarsen due to lack of memory or error */
    if (((o = kgraphMapMlUncoarsen (grafptr, NULL, NULL))        == 0) && /* Finalize graph        */
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

/* This routine performs the multi-level mapping.
** It returns:
** - 0 : if mapping could be computed.
** - 1 : on error.
*/

int
kgraphMapMl (
Kgraph * const                  grafptr,          /*+ Active graph      +*/
const KgraphMapMlParam * const  paraptr)          /*+ Method parameters +*/
{
  Gnum                levlnum;                    /* Save value for graph level */
  int                 o;

  levlnum = grafptr->levlnum;                     /* Save graph level                */
  grafptr->levlnum = 0;                           /* Initialize coarsening level     */
  o = kgraphMapMl2 (grafptr, paraptr);            /* Perform multi-level mapping     */
  grafptr->levlnum = levlnum;                     /* Restore graph level             */

  return (o);
}

