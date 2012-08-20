/* Copyright 2008,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : kgraph_map_rb_part.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module performs the Dual Recursive **/
/**                Bipartitioning mapping algorithm for    **/
/**                (eventually weighted) complete graph    **/
/**                target architectures.                   **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 16 sep 2008     **/
/**                                 to     31 aug 2011     **/
/**                                                        **/
/**   NOTES      : # This is a rewrite of kgraphMapRb()    **/
/**                  for complete-graph target topologies. **/
/**                  Its advantage over kgraphMapRbMap()   **/
/**                  is that no job arrays are allocated,  **/
/**                  which can save space for instance for **/
/**                  using the variable-sized complete     **/
/**                  graph architecture.                   **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_RB_PART

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
#include "kgraph_map_rb_part.h"

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

static
int
kgraphMapRbPart3 (
const Graph * restrict const      orggrafptr,     /* Graph to induce and bipartition         */
const GraphPart * restrict const  orgparttax,     /* Part array of original graph            */
const GraphPart                   indpartval,     /* Part of graph to consider               */
const int                         domnnum,        /* Index of domain onto which map the part */
Mapping * restrict const          mappptr)        /* Final mapping                           */
{
  Gnum               vertnum;

  if (orgparttax == NULL) {                       /* If graph is full graph */
#ifdef SCOTCH_DEBUG_KGRAPH2
    if ((orggrafptr->vnumtax != NULL) || (domnnum != 0)) {
      errorPrint ("kgraphMapRbPart3: internal error");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    memSet (mappptr->parttax + mappptr->baseval, 0, orggrafptr->vertnbr * sizeof (ArchDomNum));
  }
  else {                                          /* Graph to consider is a subgraph of the original graph       */
    if (orggrafptr->vnumtax == NULL) {            /* If original graph is not itself a subgraph                  */
      for (vertnum = orggrafptr->baseval; vertnum < orggrafptr->vertnnd; vertnum ++) { /* For all graph vertices */
        if (orgparttax[vertnum] == indpartval)    /* If vertex belongs to the right part                         */
          mappptr->parttax[vertnum] = domnnum;
      }
    }
    else {
      for (vertnum = orggrafptr->baseval; vertnum < orggrafptr->vertnnd; vertnum ++) { /* For all graph vertices */
        if (orgparttax[vertnum] == indpartval)    /* If vertex belongs to the right part                         */
          mappptr->parttax[orggrafptr->vnumtax[vertnum]] = domnnum;
      }
    }
  }

  return (0);
}

static
int
kgraphMapRbPart2 (
KgraphMapRbPartData * restrict const  topdataptr, /* Top-level graph and partition data       */
Graph * restrict const                orggrafptr, /* Graph to induce and bipartition          */
const GraphPart * restrict const      orgparttax, /* Part array of original graph to consider */
const GraphPart                       indpartval, /* Part of graph to consider                */
const Gnum                            indvertnbr, /* Number of vertices in part or in graph   */
const Anum                            domnnum)    /* Index of domain onto which map the part  */
{
  Graph               indgrafdat;
  Graph *             indgrafptr;
  Bgraph              actgrafdat;
  Anum                domnsubidx;
  Anum                domnsubdlt;
  ArchDom             domnsubtab[2];              /* Target subdomains              */
  Anum                domnsubnum[2];              /* Index of subdomains in mapping */
  Gnum                grafsubsiz[2];
  Mapping * restrict  mappptr;
  int                 avarval;                    /* Flag set if variable-sized */
  int                 i;
  int                 o;

  mappptr = topdataptr->mappptr;
  avarval = archVar (&mappptr->archdat);
  o = (avarval &&                                 /* If architecture is variable-sized   */
       (indvertnbr <= 1))                         /* And source subgraph of minimal size */
      ? 1                                         /* Then do not bipartition target more */
      : archDomBipart (&mappptr->archdat, &mappptr->domntab[domnnum], &domnsubtab[0], &domnsubtab[1]);

  switch (o) {
    case 1 :                                      /* If target domain is terminal */
      return (kgraphMapRbPart3 (orggrafptr, orgparttax, indpartval, domnnum, mappptr)); /* Update mapping and return */
    case 2 :                                      /* On error */
      errorPrint ("kgraphMapRbPart2: cannot bipartition domain");
      return     (1);
  }

  indgrafptr = orggrafptr;                        /* Assume we will work on the original graph */
  if (orgparttax != NULL) {                       /* If not the case, build induced subgraph   */
    indgrafptr = &indgrafdat;
    if (graphInducePart (orggrafptr, orgparttax, indvertnbr, indpartval, &indgrafdat) != 0) {
      errorPrint ("kgraphMapRbPart2: cannot induce graph");
      return     (1);
    }
  }

  if (bgraphInit (&actgrafdat, indgrafptr, NULL, mappptr, domnsubtab) != 0) { /* Create active graph */
    errorPrint ("kgraphMapRbPart2: cannot create bipartition graph");
    return     (1);
  }

  if (! avarval) {                                /* If not variable-sized, impose constraints on bipartition */
    double              comploadavg;

    comploadavg = (double) actgrafdat.s.velosum / (double) archDomWght (&mappptr->archdat, &mappptr->domntab[domnnum]);
    actgrafdat.compload0min = actgrafdat.compload0avg -
                              (Gnum) MIN ((topdataptr->comploadmax - comploadavg) * (double) actgrafdat.domwght[0],
                                          (comploadavg - topdataptr->comploadmin) * (double) actgrafdat.domwght[1]);
    actgrafdat.compload0max = actgrafdat.compload0avg +
                              (Gnum) MIN ((comploadavg - topdataptr->comploadmin) * (double) actgrafdat.domwght[0],
                                          (topdataptr->comploadmax - comploadavg) * (double) actgrafdat.domwght[1]);
  }

  if (bgraphBipartSt (&actgrafdat, topdataptr->paraptr->strat) != 0) { /* Perform bipartitioning */
    errorPrint ("kgraphMapRbPart2: cannot bipartition graph");
    bgraphExit (&actgrafdat);
    return     (1);
  }
  memFree (actgrafdat.frontab);                   /* Frontier array of bipartitioning graph is no longer necessary */
  actgrafdat.s.flagval &= ~BGRAPHFREEFRON;

  if (archVar (&mappptr->archdat)) {              /* If architecture is variable-sized */
    if ((actgrafdat.compload0 == 0) ||            /* If bipartition failed             */
        (actgrafdat.compload0 == actgrafdat.s.velosum))
      return (kgraphMapRbPart3 (orggrafptr, orgparttax, indpartval, domnnum, mappptr)); /* Update mapping with original domain and return */
  }

  domnsubdlt = mappptr->domnnbr - domnnum;        /* Increment in domain number */
  domnsubidx = domnnum - domnsubdlt;              /* Place where to insert subdomain */
  mappptr->domnnbr --;                            /* One less subdomain as for now   */
  grafsubsiz[0] = actgrafdat.compsize0;
  grafsubsiz[1] = actgrafdat.s.vertnbr - actgrafdat.compsize0;

  o = 0;
  for (i = 1; i >= 0; i --) {                     /* For all subparts             */
    if (grafsubsiz[i] <= 0)                       /* If subpart is empty, skip it */
      continue;

    mappptr->domnnbr ++;                          /* One more subdomain to account for */
    if (mappptr->domnnbr > mappptr->domnmax) {
      Anum                      domnmax;
      ArchDom * restrict        domntmp;

      domnmax = mappptr->domnmax + (mappptr->domnmax >> 2) + 8; /* Increase size by 25% */
      if ((domntmp = memRealloc (mappptr->domntab, domnmax * sizeof (ArchDom))) == NULL) {
        errorPrint ("kgraphMapRbPart: cannot resize structures");
        o = 1;
        break;
      }
      mappptr->domnmax = domnmax;
      mappptr->domntab = domntmp;
    }
    domnsubidx   += domnsubdlt;                   /* Compute location of subdomain */
    domnsubnum[i] = domnsubidx;                   /* Record it before recursion    */
    mappptr->domntab[domnsubidx] = domnsubtab[i]; /* Write it at this place        */
  }

  if (o == 0) {
    for (i = 1; i >= 0; i --) {                   /* For all subparts             */
      if (grafsubsiz[i] <= 0)                     /* If subpart is empty, skip it */
        continue;

      if ((o = kgraphMapRbPart2 (topdataptr, indgrafptr, actgrafdat.parttax, (GraphPart) i, grafsubsiz[i], domnsubnum[i])) != 0)
        return (1);                               /* If problem in recursion, stop */
    }
  }

  bgraphExit (&actgrafdat);                       /* Free bipartition graph (that is, parttax) */
  if (indgrafptr == &indgrafdat)                  /* If an induced subgraph had been created   */
    graphExit (indgrafptr);                       /* Free it                                   */

  return (o);
}

int
kgraphMapRbPart (
Kgraph * restrict const                 grafptr,
const KgraphMapRbParam * restrict const paraptr)
{
  KgraphMapRbPartData   topdatadat;

  topdatadat.topgrafptr  = &grafptr->s;
  topdatadat.topfrontab  = NULL;
  topdatadat.topfronnbr  = 0;
  topdatadat.mappptr     = &grafptr->m;
  topdatadat.paraptr     = paraptr;
  topdatadat.comploadmin = (1.0 - paraptr->kbalval) * grafptr->comploadrat; /* Ratio can have been tilted when working on subgraph */
  topdatadat.comploadmax = (1.0 + paraptr->kbalval) * grafptr->comploadrat;

  grafptr->m.domnnbr    = 1;                      /* Reset mapping to original (sub)domain */
  grafptr->m.domntab[0] = grafptr->m.domnorg;

  return (kgraphMapRbPart2 (&topdatadat, &grafptr->s, NULL, 0, grafptr->s.vertnbr, 0));
}
