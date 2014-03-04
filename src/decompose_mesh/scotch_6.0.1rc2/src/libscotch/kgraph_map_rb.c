/* Copyright 2004,2007,2008,2011,2013,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_rb.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module performs the Dual Recursive **/
/**                Bipartitioning mapping algorithm.       **/
/**                It is now a branching routine.          **/
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
/**                                 to     07 oct 2008     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     07 feb 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_RB

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "arch_dist.h"
#include "mapping.h"
#include "bgraph.h"
#include "bgraph_bipart_st.h"
#include "kgraph.h"
#include "kgraph_map_rb.h"
#include "kgraph_map_rb_map.h"
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

int
kgraphMapRb (
Kgraph * const                          grafptr,
const KgraphMapRbParam * restrict const paraptr)
{
  Graph                         savgrafdat;       /* Current graph with fixed vertices             */
  int                           o;
  Gnum                          vertnum;
  Arch                          usrarch;          /* Architecture given by the user                */
  Arch                          migarch;          /* Pseudo-architecture to handle migration edges */
  ArchDist *                    migarchdataptr;   /* Pointer to pseudo-architecture data           */
  Gnum *                        rvnutax;          /* Ancestor to current vertex number array       */
  Gnum *                        vnumtax;          /* Current to ancestor vertex number array       */
  Gnum                          fronnum;

  grafptr->kbalval = paraptr->kbalval;            /* Store last k-way imbalance ratio */

  rvnutax = NULL;                                 /* We do not need rvnutax if we do not have fixed vertices */
  if (grafptr->pfixtax != NULL) {                 /* We have fixed vertices */
    VertList                    indlist;
    Gnum                        vnfinum;
    Gnum                        partnum;

    indlist.vnumnbr = grafptr->s.vertnbr - grafptr->vfixnbr;
    if (memAllocGroup ((void **) (void *) 
                       &rvnutax,         (Gnum) (grafptr->s.vertnbr * sizeof (Gnum)),
                       &indlist.vnumtab, (Gnum) (indlist.vnumnbr    * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("kgraphMapRb: out of memory (1)");
      return     (1);
    }
    rvnutax -= grafptr->s.baseval;

    vnfinum = 0;
    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
      partnum = grafptr->pfixtax[vertnum];
      if (partnum == -1)                          /* If vertex not fixed */
        indlist.vnumtab[vnfinum ++] = vertnum;
    }
    
    savgrafdat = grafptr->s;
    graphInit (&grafptr->s);
    graphInduceList (&savgrafdat, &indlist, &grafptr->s, rvnutax);
    vnumtax = grafptr->s.vnumtax;
    grafptr->s.vnumtax = NULL;                    /* Since a global allocation is done vnumtax will be correctly freed */
    grafptr->m.grafptr = &grafptr->s;
  }

  if (grafptr->r.m.parttax != NULL) {             /* We are doing a repartitioning                           */
    usrarch = grafptr->a;                         /* Save the user given arch                                */
    archDistArchBuild (&grafptr->a, &usrarch, (Anum) grafptr->r.crloval); /* Build the distance architecture */
  }

  o = (archPart (grafptr->m.archptr)
       ? kgraphMapRbPart (grafptr, paraptr, &savgrafdat, rvnutax)
       : kgraphMapRbMap  (grafptr, paraptr, &savgrafdat, rvnutax));

  if (grafptr->r.m.parttax != NULL)               /* If we are doing a repartitioning */
    grafptr->a = usrarch;                         /* Restore original architecture    */

  if (rvnutax != NULL) {                          /* If we have fixed vertices */
    Gnum                orgvertnum;
    Arch *              tgtarchptr;
    Anum *              trmdomntab;
    Anum                trmdomnnbr;
    ArchDom             fstdomdat;
    Anum                domnnum;

    tgtarchptr = grafptr->m.archptr;
    archDomFrst (tgtarchptr, &fstdomdat);         /* Get first domain                                */
    trmdomnnbr = archVar (tgtarchptr)             /* Get (upper bound on) number of terminal domains */
                 ? grafptr->s.vertnbr
                 : archDomSize (tgtarchptr, &fstdomdat);

    if ((trmdomntab = memAlloc (trmdomnnbr * sizeof (Anum))) == NULL) {
      errorPrint ("kgraphMapRb: out of memory (2)");
      return     (1);
    }
    memSet (trmdomntab, ~0, trmdomnnbr * sizeof (Anum));
    for (domnnum = 0; domnnum < grafptr->m.domnnbr; domnnum ++) {
      ArchDom *                 domnptr;

      domnptr = &grafptr->m.domntab[domnnum];
#ifdef SCOTCH_DEBUG_KGRAPH2
      if (archDomSize (tgtarchptr, domnptr) != 1) {
        errorPrint ("kgraphMapRb: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      trmdomntab[archDomNum (tgtarchptr, domnptr)] = domnnum;
    }

    for (orgvertnum = savgrafdat.vertnnd - 1, vertnum = grafptr->s.vertnnd - 1; /* Prolong back part array */
         orgvertnum >= grafptr->s.baseval; ) {
      Gnum                fixpartnum;

      fixpartnum = grafptr->pfixtax[orgvertnum];
      if (fixpartnum < 0) {                       /* If original vertex was not fixed */
#ifdef SCOTCH_DEBUG_KGRAPH2
	if (vnumtax[vertnum] != orgvertnum) {     /* If vertex does not match induced vertex */
          errorPrint ("kgraphMapRb: internal error (2)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
        grafptr->m.parttax[orgvertnum --] = grafptr->m.parttax[vertnum --];
      }
      else {                                      /* Original vertex is fixed */
        Anum                fixdomnnum;

#ifdef SCOTCH_DEBUG_KGRAPH2
	if ((vertnum >= grafptr->s.baseval) &&    /* If vertex does not match induced vertex */
            (vnumtax[vertnum] >= orgvertnum)) {
          errorPrint ("kgraphMapRb: internal error (3)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

        fixdomnnum = trmdomntab[fixpartnum];      /* Get fixed domain number                */
        if (fixdomnnum == -1) {                   /* If fixed domain not already in parttab */
#ifdef SCOTCH_DEBUG_KGRAPH1
          int                 o;

          o =
#endif /* SCOTCH_DEBUG_KGRAPH1 */
          archDomTerm (tgtarchptr, &grafptr->m.domntab[grafptr->m.domnnbr], fixpartnum);
#ifdef SCOTCH_DEBUG_KGRAPH1
          if (o != 0) {
            errorPrint ("kgraphMapRb: invalid part array");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH1 */

          trmdomntab[fixpartnum] =               /* This is the new domain number */
          fixdomnnum             = grafptr->m.domnnbr ++;
        }
        grafptr->m.parttax[orgvertnum --] = fixdomnnum; /* Fill empty vertex with fixed vertex part */

        if (orgvertnum == vertnum)                /* If all remaining vertices are not fixed   */
          break;                                  /* Then beginning of array need not be moved */
      }
    }
#ifdef SCOTCH_DEBUG_KGRAPH2
    if ((orgvertnum != vertnum) && (vertnum != (grafptr->s.baseval - 1))) {
      errorPrint ("kgraphMapRb: internal error (4)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

    graphExit (&grafptr->s);
    grafptr->s = savgrafdat;                      /* Restore original graph */
    grafptr->m.grafptr = &grafptr->s;

    memFree (trmdomntab);
    memFree (rvnutax + grafptr->s.baseval);       /* Free group leader */
  }
      
  for (vertnum = grafptr->s.baseval, fronnum = 0; vertnum < grafptr->s.vertnnd; vertnum ++) {
    Gnum                partval;
    Gnum                edgenum;

    partval = grafptr->m.parttax[vertnum];

    for (edgenum = grafptr->s.verttax[vertnum]; edgenum < grafptr->s.vendtax[vertnum]; edgenum ++) {
      if (grafptr->m.parttax[grafptr->s.edgetax[edgenum]] != partval) { /* If first vertex belongs to frontier */
        grafptr->frontab[fronnum] = vertnum;
        fronnum ++;
        break;
      }
    }
  }
  grafptr->fronnbr = fronnum;

  kgraphCost (grafptr);

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (kgraphCheck (grafptr) != 0) {
    errorPrint ("kgraphMapRb: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  return (o);
}
