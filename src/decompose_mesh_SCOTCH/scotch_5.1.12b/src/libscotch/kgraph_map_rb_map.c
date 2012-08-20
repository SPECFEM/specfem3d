/* Copyright 2004,2007-2009,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : kgraph_map_rb_map.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module performs the Dual Recursive **/
/**                Bipartitioning mapping algorithm.       **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 31 mar 1993     **/
/**                                 to     31 mar 1993     **/
/**                # Version 1.0  : from : 04 oct 1993     **/
/**                                 to     06 oct 1993     **/
/**                # Version 1.1  : from : 15 oct 1993     **/
/**                                 to     15 oct 1993     **/
/**                # Version 1.3  : from : 09 apr 1994     **/
/**                                 to     11 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     17 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     18 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     19 oct 1995     **/
/**                # Version 3.1  : from : 30 oct 1995     **/
/**                                 to     14 jun 1996     **/
/**                # Version 3.2  : from : 23 aug 1996     **/
/**                                 to     07 sep 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to     08 dec 1998     **/
/**                # Version 3.4  : from : 01 jun 2001     **/
/**                                 to     07 nov 2001     **/
/**                # Version 4.0  : from : 12 jan 2004     **/
/**                                 to     06 mar 2005     **/
/**                # Version 5.1  : from : 22 nov 2007     **/
/**                                 to     31 aug 2011     **/
/**                                                        **/
/**   NOTES      : # This code is a complete rewrite of    **/
/**                  the original code of kgraphMapRb(),   **/
/**                  hence the kept history.               **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_RB_MAP

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "bgraph.h"
#include "bgraph_bipart_st.h"
#include "kgraph.h"
#include "kgraph_map_rb.h"
#include "kgraph_map_rb_map.h"

/*
**  The static variables.
*/

static KgraphMapRbMapPoolLink kgraphmaprbmappooldummy; /* Dummy links for pool routines; TRICK */

/************************************/
/*                                  */
/* These routines handle job pools. */
/*                                  */
/************************************/

/* This routine initializes the job pool
** structures.
** It returns:
** - VOID  : in all cases.
*/

static
int
kgraphMapRbMapPoolInit (
KgraphMapRbMapPoolData * restrict const poolptr,
Kgraph * restrict const                 grafptr,
const KgraphMapRbParam * restrict const paraptr)
{
  int                 flagval;

  memSet (grafptr->m.parttax + grafptr->m.baseval, 0, grafptr->s.vertnbr * sizeof (ArchDomNum)); /* Initialize partition data */
  grafptr->m.domnnbr    = 1;                      /* Only one valid domain to date */
  grafptr->m.domntab[0] = grafptr->m.domnorg;     /* All vertices are mapped to it */

  flagval = 0;
  poolptr->grafptr = NULL;                        /* Assume we don't need top level graph data */
  if (archVar (&grafptr->m.archdat) != 0)
    flagval |= KGRAPHMAPRBMAPARCHVAR;
  if (archPart (&grafptr->m.archdat) != 0)
    flagval |= KGRAPHMAPRBMAPARCHCMPLT;
  else
    poolptr->grafptr = &grafptr->s;               /* We will need top-level graph data */

  poolptr->linktab[0].prev = 
  poolptr->linktab[0].next =
  poolptr->linktab[1].prev = 
  poolptr->linktab[1].next = &kgraphmaprbmappooldummy;
  poolptr->pooltab[0] = &poolptr->linktab[0];
  poolptr->pooltab[1] = (paraptr->flagjobtie != 0) ? &poolptr->linktab[0] : &poolptr->linktab[1];
  poolptr->polival = ((flagval & KGRAPHMAPRBMAPARCHCMPLT) != 0) ? KGRAPHMAPRBPOLILEVEL : paraptr->polival;
  poolptr->mappptr = &grafptr->m;                 /* Will be used at exiting time */
  if ((poolptr->jobtab = (KgraphMapRbMapJob *) memAlloc (grafptr->m.domnmax * sizeof (KgraphMapRbMapJob))) == NULL) {
    errorPrint ("kgraphMapRbMapPoolInit: out of memory (1)");
    return     (1);
  }

  poolptr->domntmp = grafptr->m.domntab;          /* Keep track of original domain array         */
  if (paraptr->flagmaptie != 0) {                 /* If mappings are tied, use same domain array */
    poolptr->domntab = grafptr->m.domntab;
    flagval |= KGRAPHMAPRBMAPPARTHALF;
  }
  else {
    if ((poolptr->domntab = (ArchDom *) memAlloc (grafptr->m.domnmax * sizeof (ArchDom))) == NULL) {
      errorPrint ("kgraphMapRbMapPoolInit: out of memory (2)");
      memFree    (poolptr->jobtab);
      return     (1);
    }
  }

  poolptr->flagval = flagval;

  return (0);
}

/* This routine frees all of the internal arrays
** involved in the DRB algorithms. Great care
** should be taken that this routine always
** succeeds, whatever part of the algorithm it
** is called from.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolExit (
KgraphMapRbMapPoolData * restrict const poolptr)
{
  int                 jobnum;

  for (jobnum = 0; jobnum < poolptr->mappptr->domnnbr; jobnum ++) { /* For all potential jobs in both pools */
    if (poolptr->jobtab[jobnum].poolflag != 0) /* If job slot is active                                     */
      graphFree (&poolptr->jobtab[jobnum].grafdat); /* Free job graph, if not clone of original graph       */
  }

  if (poolptr->mappptr->domntab != poolptr->domntmp) { /* If current mapping domain array is not original domain array      */
    if (poolptr->domntab == poolptr->domntmp) {   /* If original domain array area was preserved                            */
      memCpy (poolptr->domntmp, poolptr->mappptr->domntab, poolptr->mappptr->domnnbr * sizeof (ArchDom)); /* Just update it */
      poolptr->domntab          = poolptr->mappptr->domntab; /* Prepare to free the other domain array */
      poolptr->mappptr->domntab = poolptr->domntmp; /* Point to original domain array                  */
    }
  }
  if (poolptr->domntab != poolptr->mappptr->domntab) /* If mappings were not tied, free second mapping array */
    memFree (poolptr->domntab);

  memFree (poolptr->jobtab);
}

/* This routine adds a job to pool 1 of the
** given pool data structure.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolAdd (
KgraphMapRbMapPoolLink * restrict const linkptr,
KgraphMapRbMapJob * const               jobptr)
{
  jobptr->poollink.prev = linkptr;                /* Link job in pool: TRICK */
  jobptr->poollink.next = linkptr->next;
  jobptr->poolflag      = 1;                      /* Job is in pool    */
  jobptr->poolptr       = linkptr;                /* Point to the pool */
  linkptr->next->prev   = &jobptr->poollink;
  linkptr->next         = &jobptr->poollink;
}

/* This routine gets the best job available from
** the given pool, according to the given policy.
** It returns:
** - !NULL  : pointer to the job.
** - NULL   : if the pool is empty.
*/

static
KgraphMapRbMapJob *
kgraphMapRbMapPoolGet (
KgraphMapRbMapPoolData * const  poolptr)
{
  KgraphMapRbMapJob * jobbest;                    /* Best job found */
  KgraphMapRbMapJob * jobptr;

  jobbest = (KgraphMapRbMapJob *) poolptr->pooltab[0]->next;  /* Get first job in pool */
  for (jobptr  = jobbest;                         /* For all jobs in pool              */
       jobptr != (KgraphMapRbMapJob *) (void *) &kgraphmaprbmappooldummy;
       jobptr  = (KgraphMapRbMapJob *) jobptr->poollink.next) {
    if (jobptr->priolvl > jobbest->priolvl)       /* If the current job has stronger priority */
      jobbest = jobptr;                           /* Select it as the best job                */
  }

  if (jobbest != (KgraphMapRbMapJob *) (void *) &kgraphmaprbmappooldummy) { /* If job found */
    jobbest->poollink.next->prev = jobbest->poollink.prev; /* Remove it from pool           */
    jobbest->poollink.prev->next = jobbest->poollink.next; /* But do not mark it unused     */
  }
  else                                            /* Dummy job means no job found */
    jobbest = NULL;

  return (jobbest);
}

/* This routine adds a job to the given pool
** as the first bipartitioning job.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolFrst (
KgraphMapRbMapPoolData * const  poolptr,
KgraphMapRbMapJob * const       jobptr)           /* Job to be added */
{
  switch (poolptr->polival) {                     /* Set job priority value */
    case KGRAPHMAPRBPOLIRANDOM :
      jobptr->prioval =
      jobptr->priolvl = intRandVal (INTVALMAX);
      break;
    case KGRAPHMAPRBPOLILEVEL   :
    case KGRAPHMAPRBPOLINGLEVEL :
      jobptr->prioval = jobptr->grafdat.vertnbr;
      jobptr->priolvl = 0;
      break;
    case KGRAPHMAPRBPOLISIZE   :
    case KGRAPHMAPRBPOLINGSIZE :
      jobptr->prioval =
      jobptr->priolvl = jobptr->grafdat.vertnbr;
      break;
#ifdef SCOTCH_DEBUG_KGRAPH2
    default :
      errorPrint ("kgraphMapRbMapPoolFrst: unknown job selection policy");
      jobptr->prioval = 0;
      jobptr->priolvl = 0;
      return;
#endif /* SCOTCH_DEBUG_KGRAPH2 */
  }

  kgraphMapRbMapPoolAdd (poolptr->pooltab[0], jobptr); /* Add job to pool */
}

/* This routine updates the given job
** table with both of the given subjob
** data.
** This routine can be called only if
** the parent jobs of the vertices to
** be updated still exist.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolUpdt1 (
KgraphMapRbMapPoolData * const  poolptr,
const KgraphMapRbMapJob * const joboldptr,        /* Job to be removed      */
const GraphPart * const         parttax,
KgraphMapRbMapJob * const       jobnewptr,        /* Its only active subjob */
const GraphPart                 partval)
{
  KgraphMapRbMapJob * restrict  jobtab;
  const Anum * restrict         mapparttax;       /* Based pointer to mapping part array */
  const Gnum * restrict         topverttax;
  const Gnum * restrict         topvendtax;
  const Gnum * restrict         topedgetax;
  Gnum                          prioval;
  Gnum                          priolvl;

  priolvl = 0;                                    /* Prepare for neighbor updating methods */

  switch (poolptr->polival) {                     /* Set job priority value */
    case KGRAPHMAPRBPOLIRANDOM :
      prioval =
      priolvl = intRandVal (INTVALMAX);
      break;
    case KGRAPHMAPRBPOLILEVEL :
      priolvl = joboldptr->priolvl + 1;
    case KGRAPHMAPRBPOLINGLEVEL :
      prioval = joboldptr->prioval - 1;
      break;
    case KGRAPHMAPRBPOLISIZE :
      priolvl = jobnewptr->grafdat.vertnbr;
    case KGRAPHMAPRBPOLINGSIZE :
      prioval = jobnewptr->grafdat.vertnbr;
      break;
#ifdef SCOTCH_DEBUG_KGRAPH2
    default :
      errorPrint ("kgraphMapRbMapPoolUpdt1: unknown job selection policy");
      jobnewptr->prioval = 0;
      jobnewptr->priolvl = 0;
      return;
#endif /* SCOTCH_DEBUG_KGRAPH2 */
  }

  jobnewptr->prioval = prioval;

  if (poolptr->polival >= KGRAPHMAPRBPOLINEIGHBOR) { /* If neighbors have to be updated */
    Gnum                jobvertnum;
    Gnum                prioold;

    jobtab     = poolptr->jobtab;
    mapparttax = poolptr->mappptr->parttax;
    topverttax = poolptr->grafptr->verttax;       /* Point to top-level graph arrays */
    topvendtax = poolptr->grafptr->vendtax;
    topedgetax = poolptr->grafptr->edgetax;
  
    prioold = joboldptr->prioval;

    if (joboldptr->grafdat.vertnbr < poolptr->grafptr->vertnbr) { /* If subgraph is not top graph */
      const Gnum * restrict jobvnumtax;
      const Gnum * restrict jobverttax;
      const Gnum * restrict jobvendtax;

      jobvnumtax = joboldptr->grafdat.vnumtax;    /* Change priority of neighboring jobs of old job */
      jobverttax = joboldptr->grafdat.verttax;
      jobvendtax = joboldptr->grafdat.vendtax;

      jobnewptr->poolflag = 0;                    /* TRICK: avoid new job being considered for update */

      for (jobvertnum = joboldptr->grafdat.baseval; jobvertnum < joboldptr->grafdat.vertnnd; jobvertnum ++) {
        Gnum                topvertnum;
        Gnum                topedgenum;

        if (parttax[jobvertnum] == partval)       /* If vertex belongs to part which is still alive */
          continue;                               /* Do not consider update part as removed         */

        topvertnum = jobvnumtax[jobvertnum];      /* If graph is smaller than top graph, then vnumtax must exist */

        if ((topvendtax[topvertnum] - topverttax[topvertnum]) == /* If vertex is internal, skip it */
            (jobvendtax[jobvertnum] - jobverttax[jobvertnum]))
          continue;

        for (topedgenum = topverttax[topvertnum]; topedgenum < topvendtax[topvertnum]; topedgenum ++) {
          KgraphMapRbMapJob * restrict  jobnghbptr; /* (Old ?) job of neighbor vertex */

          jobnghbptr = &jobtab[mapparttax[topedgetax[topedgenum]]]; /* Get pointer to neighboring job */

          if ((jobnghbptr->poolflag != 0) &&      /* If neighbor is active                   */
              (jobnghbptr->prioval <= prioold))   /* And had not already a stronger priority */
            jobnghbptr->priolvl ++;               /* Update neighbor priority                */
        }
      }

      jobnewptr->poolflag = 1;                    /* TRICK: new job is active again */
    }

    if (jobnewptr->grafdat.vertnbr < poolptr->grafptr->vertnbr) { /* If subgraph is not top graph, update priority of neighbors of new job only */
      const Gnum * restrict jobvnumtax;
      const Gnum * restrict jobverttax;
      const Gnum * restrict jobvendtax;

      jobvnumtax = jobnewptr->grafdat.vnumtax;
      jobverttax = jobnewptr->grafdat.verttax;
      jobvendtax = jobnewptr->grafdat.vendtax;

      for (jobvertnum = jobnewptr->grafdat.baseval; jobvertnum < jobnewptr->grafdat.vertnnd; jobvertnum ++) {
        Gnum                          topvertnum;
        Gnum                          topedgenum;

        topvertnum = jobvnumtax[jobvertnum];      /* For subjobs jobvnumtax always exists */

        if ((topvendtax[topvertnum] - topverttax[topvertnum]) == /* If vertex is internal, skip it */
            (jobvendtax[jobvertnum] - jobverttax[jobvertnum]))
          continue;

        for (topedgenum = topverttax[topvertnum]; topedgenum < topvendtax[topvertnum]; topedgenum ++) {
          KgraphMapRbMapJob * restrict  jobnghbptr; /* (Old ?) job of neighbor vertex */

          jobnghbptr = &jobtab[mapparttax[topedgetax[topedgenum]]]; /* Get pointer to neighboring job   */
          if (jobnghbptr == jobnewptr)            /* If it is the current job, do not consider the edge */
            continue;

          if ((jobnghbptr->poolflag == 0) ||      /* If neighbor is not active            */
              (prioval > jobnghbptr->prioval))    /* Or if we have higher priority        */
            priolvl ++;                           /* Increase our priority                */
          else if ((prioval <  jobnghbptr->prioval) && /* Else if neighbor has higher one */
                   (prioold >= jobnghbptr->prioval)) /* Which it did not already have     */
            jobnghbptr->priolvl ++;               /* Update neighbor priority             */
        }
      }
    }
  }

  jobnewptr->priolvl = priolvl;

  kgraphMapRbMapPoolAdd (poolptr->pooltab[1], jobnewptr); /* Add job to pool */
}

static
void
kgraphMapRbMapPoolUpdt2 (
KgraphMapRbMapPoolData * const  poolptr,
const KgraphMapRbMapJob * const joboldptr,        /* Job to be removed */
const GraphPart * const         parttax,
KgraphMapRbMapJob * const       jobnewptr0,       /* Its two subjobs   */
KgraphMapRbMapJob * const       jobnewptr1)
{
  KgraphMapRbMapJob * restrict  jobtab;
  const Anum * restrict         mapparttax;       /* Based pointer to mapping part array */
  const Gnum * restrict         jobvnumtax;
  const Gnum * restrict         jobverttax;
  const Gnum * restrict         jobvendtax;
  const Gnum * restrict         topverttax;
  const Gnum * restrict         topvendtax;
  const Gnum * restrict         topedgetax;
  KgraphMapRbMapJob * restrict  jobnewtab[2];
  int                           i;

  jobnewtab[0] = jobnewptr0;
  jobnewtab[1] = jobnewptr1;

  for (i = 1; i >= 0; i --) {
    KgraphMapRbMapJob * jobnewptr;
    Gnum                prioval;
    Gnum                priolvl;

    jobnewptr = jobnewtab[i];                     /* Get concerned subjob */

    priolvl = 0;                                  /* Prepare for neighbor updating methods */

    switch (poolptr->polival) {                   /* Set job priority value */
      case KGRAPHMAPRBPOLIRANDOM :
        prioval =
        priolvl = intRandVal (INTVALMAX);
        break;
      case KGRAPHMAPRBPOLILEVEL :
        priolvl = joboldptr->priolvl + 1;
      case KGRAPHMAPRBPOLINGLEVEL :
        prioval = joboldptr->prioval - 1;
        break;
      case KGRAPHMAPRBPOLISIZE :
        priolvl = jobnewptr->grafdat.vertnbr;
      case KGRAPHMAPRBPOLINGSIZE :
        prioval = jobnewptr->grafdat.vertnbr;
        break;
#ifdef SCOTCH_DEBUG_KGRAPH2
      default :
        errorPrint ("kgraphMapRbMapPoolUpdt2: unknown job selection policy");
        jobnewptr->prioval = 0;
        jobnewptr->priolvl = 0;
        return;
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    }

    jobnewptr0->prioval = prioval + 1;            /* TRICK: when processing subdomain 1, subdomain 0 has higher priority value */
    jobnewptr->prioval  = prioval;                /* Then in its turn subdomain 0 will have its proper priority value          */

    if (poolptr->polival >= KGRAPHMAPRBPOLINEIGHBOR) { /* If neighbors have to be updated */
      Gnum                jobvertnum;
      Gnum                prioold;

      jobtab     = poolptr->jobtab;
      mapparttax = poolptr->mappptr->parttax;
      topverttax = poolptr->grafptr->verttax;     /* Point to top-level graph arrays */
      topvendtax = poolptr->grafptr->vendtax;
      topedgetax = poolptr->grafptr->edgetax;
      jobvnumtax = jobnewptr->grafdat.vnumtax;
      jobverttax = jobnewptr->grafdat.verttax;
      jobvendtax = jobnewptr->grafdat.vendtax;

      prioold = joboldptr->prioval;
      for (jobvertnum = jobnewptr->grafdat.baseval; jobvertnum < jobnewptr->grafdat.vertnnd; jobvertnum ++) {
        Gnum                          topvertnum;
        Gnum                          topedgenum;

        topvertnum = jobvnumtax[jobvertnum];      /* For subjobs jobvnumtax always exists */

        if ((topvendtax[topvertnum] - topverttax[topvertnum]) == /* If vertex is internal, skip it */
            (jobvendtax[jobvertnum] - jobverttax[jobvertnum]))
          continue;

        for (topedgenum = topverttax[topvertnum]; topedgenum < topvendtax[topvertnum]; topedgenum ++) {
          KgraphMapRbMapJob * jobnghbptr;         /* (Old ?) job of neighbor vertex */

          jobnghbptr = &jobtab[mapparttax[topedgetax[topedgenum]]]; /* Get pointer to neighboring job */

          if ((jobnghbptr->poolflag != 0)      && /* If neighbor is in active job  */
              (jobnghbptr->prioval >  prioval) && /* Which gained priority over us */
              (jobnghbptr->prioval <= prioold)) {
            jobnghbptr->priolvl ++;               /* Update neighbor priority */
          }
          if ((jobnghbptr->poolflag == 0) ||      /* If neighbor is fully known    */
              (jobnghbptr->prioval < prioval))    /* Or has smaller priority value */
            priolvl ++;                           /* Then we should be processed   */
        }
      }
    }

    jobnewptr->priolvl = priolvl;                 /* Set new priority          */
    kgraphMapRbMapPoolAdd (poolptr->pooltab[1], jobnewptr); /* Add job to pool */
  }
}

/*
** This routine removes the influence of the
** given job from its neighbor jobs.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolRemv (
KgraphMapRbMapPoolData * const  poolptr,
const KgraphMapRbMapJob * const joboldptr)        /* Job to be removed */
{
  KgraphMapRbMapJob * restrict  jobtab;
  const Anum * restrict         mapparttax;       /* Based pointer to mapping part array */
  const Gnum * restrict         jobvnumtax;
  const Gnum * restrict         jobverttax;
  const Gnum * restrict         jobvendtax;
  const Gnum * restrict         topverttax;
  const Gnum * restrict         topvendtax;
  const Gnum * restrict         topedgetax;

  if (poolptr->polival >= KGRAPHMAPRBPOLINEIGHBOR) { /* If neighbors have to be modified */
    Gnum                jobvertnum;
    Gnum                prioold;

    jobtab     = poolptr->jobtab;
    mapparttax = poolptr->mappptr->parttax;
    topverttax = poolptr->grafptr->verttax;       /* Point to top-level graph arrays */
    topvendtax = poolptr->grafptr->vendtax;
    topedgetax = poolptr->grafptr->edgetax;
    jobvnumtax = joboldptr->grafdat.vnumtax;
    jobverttax = joboldptr->grafdat.verttax;
    jobvendtax = joboldptr->grafdat.vendtax;

    prioold = joboldptr->prioval;

    for (jobvertnum = joboldptr->grafdat.baseval; jobvertnum < joboldptr->grafdat.vertnnd; jobvertnum ++) {
      Gnum                topvertnum;             /* Source graph vertex number */
      Gnum                topedgenum;             /* Source graph edge number   */

      topvertnum = (jobvnumtax == NULL) ? jobvertnum : jobvnumtax[jobvertnum];

      if ((topvendtax[topvertnum] - topverttax[topvertnum]) == /* If vertex is internal, skip it */
          (jobvendtax[jobvertnum] - jobverttax[jobvertnum]))
        continue;

      for (topedgenum = topverttax[topvertnum]; topedgenum < topvendtax[topvertnum]; topedgenum ++) {
        KgraphMapRbMapJob * jobnghbptr;           /* (Old ?) job of neighbor vertex */

        jobnghbptr = &jobtab[mapparttax[topedgetax[topedgenum]]]; /* Get pointer to neighboring job */

        if ((jobnghbptr->poolflag != 0) &&        /* If neighbor job is active                       */
            (jobnghbptr->prioval <= prioold))     /* And had not already a stronger priority         */
          jobnghbptr->priolvl ++;                 /* Increase its priority since we are now inactive */
      }
    }
  }
}

/**********************************************/
/*                                            */
/* These routines handle the pool part array. */
/*                                            */
/**********************************************/

static
void
kgraphMapRbMapPartBoth (
KgraphMapRbMapPoolData * restrict const poolptr,
const Bgraph * restrict const           bipgrafptr,
const Anum * restrict const             jobsubnum)
{
  Gnum                       bipvertnum;
  const GraphPart * restrict bipparttax;
  Anum * restrict            mapparttax;
  Anum                       mappartval1;
  Anum                       mappartdlt;

  bipparttax  = bipgrafptr->parttax;
  mapparttax  = poolptr->mappptr->parttax;
  mappartval1 = jobsubnum[1];
  mappartdlt  = jobsubnum[0] - jobsubnum[1];

  if (bipgrafptr->s.vnumtax != NULL) {
    const Gnum * restrict      bipvnumtax;

    bipvnumtax = bipgrafptr->s.vnumtax;
    for (bipvertnum = bipgrafptr->s.baseval; bipvertnum < bipgrafptr->s.vertnnd; bipvertnum ++)
      mapparttax[bipvnumtax[bipvertnum]] = mappartval1 + ((((Anum) bipparttax[bipvertnum]) - 1) & mappartdlt);
  }
  else {
    for (bipvertnum = bipgrafptr->s.baseval; bipvertnum < bipgrafptr->s.vertnnd; bipvertnum ++)
      mapparttax[bipvertnum] = mappartval1 + ((((Anum) bipparttax[bipvertnum]) - 1) & mappartdlt);
  }
}

static
void
kgraphMapRbMapPartOne (
KgraphMapRbMapPoolData * restrict const poolptr,
const Bgraph * restrict const           bipgrafptr,
const Anum                              jobsubnum1)
{
  Gnum                       bipvertnum;
  const GraphPart * restrict bipparttax;
  Anum * restrict            mapparttax;

  bipparttax  = bipgrafptr->parttax;
  mapparttax  = poolptr->mappptr->parttax;

  if (bipgrafptr->s.vnumtax != NULL) {
    const Gnum * restrict      bipvnumtax;

    bipvnumtax = bipgrafptr->s.vnumtax;
    for (bipvertnum = bipgrafptr->s.baseval; bipvertnum < bipgrafptr->s.vertnnd; bipvertnum ++)
      if (bipparttax[bipvertnum] == 1)
        mapparttax[bipvnumtax[bipvertnum]] = jobsubnum1;
  }
  else {
    for (bipvertnum = bipgrafptr->s.baseval; bipvertnum < bipgrafptr->s.vertnnd; bipvertnum ++)
      if (bipparttax[bipvertnum] == 1)
        mapparttax[bipvertnum] = jobsubnum1;
  }
}

/********************************************/
/*                                          */
/* This is the entry point for the Dual     */
/* Recursive Bipartitioning mapping method. */
/*                                          */
/********************************************/

/* This routine runs the Dual Recursive
** Bipartitioning algorithm.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kgraphMapRbMap (
Kgraph * restrict const                 grafptr,
const KgraphMapRbParam * restrict const paraptr)
{
  KgraphMapRbMapPoolData  pooldat;                /* Data for handling jobs and job pools */
  ArchDom                 domsubtab[2];           /* Subdomains of current job domain     */
  KgraphMapRbMapJob       joborgdat;              /* Aera to save original job data       */
  Anum                    jobsubnum[2];           /* Number of subjob slots in job array  */
  Gnum                    jobsubsiz[2];           /* Sizes of subjobs                     */
  Bgraph                  bipgrafdat;             /* Bipartition graph                    */
  double                  comploadmin;            /* Minimum vertex load per target load  */
  double                  comploadmax;            /* Maximum vertex load per target load  */
  int                     i;

#ifdef SCOTCH_DEBUG_KGRAPH2
  grafptr->m.domnmax = 1;                         /* Force resizing of job arrays, for debugging */
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  if (kgraphMapRbMapPoolInit (&pooldat, grafptr, paraptr) != 0) /* Initialize pool data; done first for kgraphMapRbMapPoolExit() to succeed afterwards */
    return (1);

  if ((((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) == 0) && (archDomSize (&grafptr->m.archdat, &grafptr->m.domnorg) <= 1)) || /* If single-vertex domain */
      (((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) != 0) && (grafptr->s.vertnbr <= 1))) { /* Or if variable-sized architecture with single vertex graph  */
    kgraphMapRbMapPoolExit (&pooldat);
    return (0);                                   /* Job already done */
  }

  pooldat.jobtab[0].domnorg = grafptr->m.domnorg; /* Build first job                         */
  pooldat.jobtab[0].grafdat = grafptr->s;         /* Clone original graph as first job graph */
  pooldat.jobtab[0].grafdat.flagval &= ~GRAPHFREETABS; /* Do not free its arrays on exit     */
  kgraphMapRbMapPoolFrst (&pooldat, &pooldat.jobtab[0]); /* Add initial job                  */

  comploadmin = (1.0 - paraptr->kbalval) * grafptr->comploadrat; /* Ratio can have been tilted when working on subgraph */
  comploadmax = (1.0 + paraptr->kbalval) * grafptr->comploadrat;

  while (! kgraphMapRbMapPoolEmpty (&pooldat)) {  /* For all non-empty pools */
    KgraphMapRbMapJob * joborgptr;                /* Pointer to current job  */

    while ((joborgptr = kgraphMapRbMapPoolGet (&pooldat)) != NULL) { /* For all jobs in pool */
      int                 partval;

      jobsubnum[0] = joborgptr - pooldat.jobtab;  /* Get current (and first son) job slot number before possible move of pointers */
      joborgdat = *joborgptr;                     /* Save current job data (clone graph)                                          */

      if (archDomBipart (&grafptr->m.archdat, &joborgdat.domnorg, &domsubtab[0], &domsubtab[1]) != 0) {
        errorPrint ("kgraphMapRbMap: cannot bipartition domain");
        kgraphMapRbMapPoolExit (&pooldat);        /* Copied graph will be freed as not yet removed */
        return (1);
      }

      if (bgraphInit (&bipgrafdat, &joborgdat.grafdat, /* Create bipartition graph */
                      pooldat.grafptr, pooldat.mappptr, domsubtab) != 0) {
        errorPrint ("kgraphMapRbMap: cannot create bipartition graph");
        kgraphMapRbMapPoolExit (&pooldat);        /* Copied graph will be freed as not yet removed */
        return (1);
      }
      bipgrafdat.s.flagval |= (joborgdat.grafdat.flagval & GRAPHFREETABS); /* Bipartition graph is responsible for freeing the cloned graph data fields */
      joborgptr->poolflag = 0;                    /* Original slot is now considered unused so that cloned graph data will not be freed twice           */

      if ((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) == 0) { /* If not variable-sized, impose constraints on bipartition */
        double              comploadavg;

        comploadavg = (double) bipgrafdat.s.velosum / (double) archDomWght (&grafptr->m.archdat, &joborgdat.domnorg);
        bipgrafdat.compload0min = bipgrafdat.compload0avg -
                                  (Gnum) MIN ((comploadmax - comploadavg) * (double) bipgrafdat.domwght[0],
                                              (comploadavg - comploadmin) * (double) bipgrafdat.domwght[1]);
        bipgrafdat.compload0max = bipgrafdat.compload0avg +
                                  (Gnum) MIN ((comploadavg - comploadmin) * (double) bipgrafdat.domwght[0],
                                              (comploadmax - comploadavg) * (double) bipgrafdat.domwght[1]);
      }

      if (bgraphBipartSt (&bipgrafdat, paraptr->strat) != 0) { /* Perform bipartitioning */
        errorPrint             ("kgraphMapRbMap: cannot bipartition job");
        bgraphExit             (&bipgrafdat);
        kgraphMapRbMapPoolExit (&pooldat);
        return                 (1);
      }
      memFree (bipgrafdat.frontab);               /* Frontier array of bipartitioning graph is no longer necessary */
      bipgrafdat.frontab = NULL;

      if ((partval = 1, bipgrafdat.compload0 == 0) || /* If no bipartition found */
          (partval = 0, bipgrafdat.compload0 == bipgrafdat.s.velosum)) {
        if (((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) != 0) || /* If architecture is variable-sized     */
            (archDomSize (&grafptr->m.archdat, &domsubtab[partval]) <= 1)) { /* Or domain is terminal    */
          pooldat.domntab[jobsubnum[0]] = joborgdat.domnorg; /* Update domain in next pool               */
          kgraphMapRbMapPoolRemv (&pooldat, &joborgdat); /* Remove job from pool as long as graph exists */
          bgraphExit (&bipgrafdat);               /* Free bipartitioning data as well as current graph   */
          continue;                               /* Process next job in current pool                    */
        }
        else {
          joborgptr->domnorg =                    /* New job takes same graph and non-empty subdomain                */
          pooldat.domntab[jobsubnum[0]] = domsubtab[partval]; /* Update domain in next pool                          */
          kgraphMapRbMapPoolUpdt1 (&pooldat, &joborgdat, bipgrafdat.parttax, joborgptr, partval); /* Add job to pool */
          continue;                               /* Process next job in current pool */
        }
      }

      if ((pooldat.mappptr->domnnbr == pooldat.mappptr->domnmax) && /* If all job slots busy and if cannot resize */
          (kgraphMapRbMapPoolResize (&pooldat) != 0)) {
        errorPrint             ("kgraphMapRbMap: cannot resize structures");
        kgraphMapRbMapPoolExit (&pooldat);
        return                 (1);
      }

      jobsubnum[1] = pooldat.mappptr->domnnbr ++; /* Get slot number of new subdomain */
      jobsubsiz[1] = bipgrafdat.s.vertnbr - bipgrafdat.compsize0;
      jobsubsiz[0] = bipgrafdat.compsize0;

      pooldat.jobtab[jobsubnum[1]].poolflag = 0;  /* Assume that new job is inactive in case of premature freeing                                 */
      pooldat.mappptr->domntab[jobsubnum[1]] = joborgdat.domnorg; /* Copy original domain to new subdomain as old mapping shares parttax with new */
      pooldat.domntab[jobsubnum[0]] = domsubtab[0]; /* Set subdomains of second mapping before relinking subjobs in pool                          */
      pooldat.domntab[jobsubnum[1]] = domsubtab[1];

      if ((pooldat.flagval & KGRAPHMAPRBMAPPARTHALF) != 0) /* If can only update second half */
        kgraphMapRbMapPartOne (&pooldat, &bipgrafdat, jobsubnum[1]);
      else
        kgraphMapRbMapPartBoth (&pooldat, &bipgrafdat, jobsubnum);

      for (i = 1; i >= 0; i --) {                 /* For both subdomains */
        KgraphMapRbMapJob * jobsubptr;

        jobsubptr = &pooldat.jobtab[jobsubnum[i]]; /* Point to subdomain job slot                                */
        jobsubptr->poollink.prev =                /* Prevent Valgrind from yelling in kgraphMapRbMapPoolResize() */
        jobsubptr->poollink.next = NULL;
        jobsubptr->prioval =                      /* Prevent Valgrind from yelling in kgraphMapRbMapPoolRemv()/Updt1()/Updt2() */
        jobsubptr->priolvl = 0;

        if ((((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) == 0) && (archDomSize (&grafptr->m.archdat, &domsubtab[i]) <= 1)) || /* If single-vertex domain */
            (((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) != 0) && (jobsubsiz[i] <= 1))) { /* Or if variable-sized architecture with single vertex graph  */
          jobsubsiz[i] = 0;                       /* Cancel subjob */
          continue;
        }

        partval = i;                              /* At least this subjob works   */

        if (graphInducePart (&bipgrafdat.s, bipgrafdat.parttax, jobsubsiz[i], (GraphPart) i, &jobsubptr->grafdat) != 0) {
          errorPrint             ("kgraphMapRbMap: cannot create induced subgraph");
          bgraphExit             (&bipgrafdat);
          kgraphMapRbMapPoolExit (&pooldat);
          return                 (1);
        }
        jobsubptr->poolflag = 1;                  /* So that graph is freed in case of error on other part */
        jobsubptr->domnorg  = domsubtab[i];
      }

      if ((jobsubsiz[0] | jobsubsiz[1]) == 0)     /* If both subjobs do not need further processing */
        kgraphMapRbMapPoolRemv (&pooldat, &joborgdat);
      else if (jobsubsiz[1 - partval] == 0)       /* If one of the subjobs only needs further processing */
        kgraphMapRbMapPoolUpdt1 (&pooldat, &joborgdat, bipgrafdat.parttax, &pooldat.jobtab[jobsubnum[partval]], (GraphPart) partval);
      else
        kgraphMapRbMapPoolUpdt2 (&pooldat, &joborgdat, bipgrafdat.parttax, &pooldat.jobtab[jobsubnum[0]], &pooldat.jobtab[jobsubnum[1]]);

      bgraphExit (&bipgrafdat);                   /* Free bipartition graph data */
    }

    kgraphMapRbMapPoolSwap (&pooldat);            /* Swap current and next levels */
  }

  kgraphMapRbMapPoolExit (&pooldat);              /* Free internal structures and propagate back new partition */

  return (0);
}

/**********************************************/
/*                                            */
/* These routines handle internal structures. */
/*                                            */
/**********************************************/

/* This routine doubles the size all of the arrays
** involved in handling the target architecture,
** to make room for new domains of variable-sized
** architectures.
** It returns:
** - 0   : if resize succeeded.
** - !0  : if out of memory.
*/

static
int
kgraphMapRbMapPoolResize (
KgraphMapRbMapPoolData * restrict const poolptr)
{
  KgraphMapRbMapJob * restrict  joboldtab;        /* Pointer to old job array          */
  ArchDom *                     domnoldtab;       /* Temporary pointer to domain array */
  Anum                          domnmax;

  domnmax = poolptr->mappptr->domnmax * 2;

  joboldtab = poolptr->jobtab;                    /* Save old job table address */
  if ((poolptr->jobtab = (KgraphMapRbMapJob *) memRealloc (joboldtab, domnmax * sizeof (KgraphMapRbMapJob))) == NULL) {
    errorPrint ("kgraphMapRbMapPoolResize: out of memory (1)");
    poolptr->jobtab = joboldtab;
    return (1);
  }

  if (poolptr->jobtab != joboldtab) {             /* If job array moved              */
    KgraphMapRbMapJob *       joboldtnd;          /* Pointer to end of old job array */
    Anum                      jobnum;             /* Temporary job index             */
    size_t                    jobdlt;             /* Address delta value             */

    joboldtnd = joboldtab + poolptr->mappptr->domnmax;
    jobdlt = (byte *) poolptr->jobtab - (byte *) joboldtab; /* Compute delta between addresses */

    for (jobnum = 0; jobnum < poolptr->mappptr->domnmax; jobnum ++) {
      if ((poolptr->jobtab[jobnum].poollink.prev >= (KgraphMapRbMapPoolLink *) joboldtab) && /* If old pointers within bounds of old array, adjust them */
          (poolptr->jobtab[jobnum].poollink.prev <  (KgraphMapRbMapPoolLink *) joboldtnd))
        poolptr->jobtab[jobnum].poollink.prev = (KgraphMapRbMapPoolLink *) ((byte *) poolptr->jobtab[jobnum].poollink.prev + jobdlt);
      if ((poolptr->jobtab[jobnum].poollink.next >= (KgraphMapRbMapPoolLink *) joboldtab) &&
          (poolptr->jobtab[jobnum].poollink.next <  (KgraphMapRbMapPoolLink *) joboldtnd))
        poolptr->jobtab[jobnum].poollink.next = (KgraphMapRbMapPoolLink *) ((byte *) poolptr->jobtab[jobnum].poollink.next + jobdlt);
    }
    if (poolptr->linktab[0].next != &kgraphmaprbmappooldummy) /* Update first pool pointer */
      poolptr->linktab[0].next = (KgraphMapRbMapPoolLink *) ((byte *) poolptr->linktab[0].next + jobdlt);
    if (poolptr->pooltab[0] != poolptr->pooltab[1]) { /* If job pools not tied                */
      if (poolptr->linktab[1].next != &kgraphmaprbmappooldummy) /* Update second pool pointer */
        poolptr->linktab[1].next = (KgraphMapRbMapPoolLink *) ((byte *) poolptr->linktab[1].next + jobdlt);
    }
  }

  domnoldtab = poolptr->mappptr->domntab;
  if ((poolptr->mappptr->domntab = (ArchDom *) memRealloc (domnoldtab, domnmax * sizeof (ArchDom))) == NULL) {
    errorPrint ("kgraphMapRbMapPoolResize: out of memory (2)");
    poolptr->mappptr->domntab = domnoldtab;
    return (1);
  }

  if (poolptr->domntab != domnoldtab) {           /* If mappings not tied */
    domnoldtab = poolptr->domntab;
    if ((poolptr->domntab = (ArchDom *) memRealloc (domnoldtab, domnmax * sizeof (ArchDom))) == NULL) {
      errorPrint ("kgraphMapRbMapPoolResize: out of memory (3)");
      poolptr->domntab = domnoldtab;
      return (1);
    }
  }
  else
    poolptr->domntab = poolptr->mappptr->domntab; /* Update second domain array in case first one moved */

  poolptr->mappptr->domnmax = domnmax;            /* Double job slot limit */

  return (0);
}
