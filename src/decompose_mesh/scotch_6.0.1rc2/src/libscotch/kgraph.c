/* Copyright 2004,2007,2008,2010-2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : Part of a bipartitioning mapper.        **/
/**                This module handles the k-way active    **/
/**                graph and save data structure handling  **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 03 oct 1997     **/
/**                                 to     26 may 1998     **/
/**                # Version 3.4  : from : 30 oct 2001     **/
/**                                 to     30 oct 2001     **/
/**                # Version 4.0  : from : 24 jun 2004     **/
/**                                 to     16 feb 2005     **/
/**                # Version 5.1  : from : 28 sep 2008     **/
/**                                 to     31 aug 2011     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     17 oct 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "kgraph.h"

/***********************************/
/*                                 */
/* Active graph handling routines. */
/*                                 */
/***********************************/

/* This routine builds the active graph
** corresponding to the given k-way
** partition parameters.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kgraphInit (
Kgraph * restrict const         actgrafptr,       /*+ Active graph                         +*/
const Graph * restrict const    srcgrafptr,       /*+ Source graph                         +*/
const Arch * restrict const     archptr,          /*+ Target architecture                  +*/
const ArchDom * restrict const  archdomptr,       /*+ Target architecture initial domain   +*/
Anum * restrict const           parttax,          /*+ Partition array                      +*/
Anum * restrict const           parotax,          /*+ Old partition array                  +*/
const Gnum                      crloval,          /*+ Coefficient load for regular edges   +*/
const Gnum                      cmloval,          /*+ Coefficient load for migration edges +*/
const Gnum *                    vmlotax,          /*+ Vertex migration cost array          +*/
Anum *                          pfixtax)          /*+ Fixed partition array                +*/
{
  ArchDom                   domfrst;              /* First, largest domain */
  const ArchDom * restrict  domfrstptr;
  Gnum                      domfrstload;          /* Load of first domain  */
  Gnum                      vfixnbr;
  Anum                      termnum;

  if (srcgrafptr != &actgrafptr->s) {             /* If graph not self-provided, clone source graph */
    actgrafptr->s          = *srcgrafptr;
    actgrafptr->s.flagval &= (GRAPHBITSUSED & ~GRAPHFREETABS); /* Remove extended graph class flags and do not allow freeing */
  }
  if (archptr != &actgrafptr->a)                  /* If architecture not self-provided, clone architecture */
    actgrafptr->a = *archptr;
  actgrafptr->pfixtax = pfixtax;

  archDomFrst (archptr, &domfrst);                /* Get first, largest domain          */
  domfrstptr = (archdomptr == NULL) ? &domfrst : archdomptr; /* Use it if none provided */

  mapInit (&actgrafptr->m, &actgrafptr->s, &actgrafptr->a, domfrstptr, parttax, 0);
  if (parotax != NULL) {
    mapInit (&actgrafptr->r.m, &actgrafptr->s, &actgrafptr->a, domfrstptr, parotax, 1);

    if (vmlotax != NULL)
      actgrafptr->r.vmlotax = vmlotax;
    else {
      Gnum *              vmlotab;
      Gnum                vertnum;

      if ((vmlotab = memAlloc (actgrafptr->s.vertnbr * sizeof (Gnum))) == NULL) {
        errorPrint ("kgraphInit: out of memory (1)");
        return     (1);
      }
      actgrafptr->s.flagval |= KGRAPHFREEVMLO;
      for (vertnum = 0; vertnum < actgrafptr->s.vertnbr; vertnum ++)
        vmlotab[vertnum] = 1;

      actgrafptr->r.vmlotax = vmlotab - actgrafptr->s.baseval;
    }
  }
  else {
    memSet (&actgrafptr->r.m, 0, sizeof (Mapping));
    actgrafptr->r.vmlotax = NULL;
  }
#ifdef SCOTCH_DEBUG_KGRAPH2
  if ((crloval < 1) || (cmloval < 1)) {
    errorPrint ("kgraphInit: invalid parameters");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
  actgrafptr->r.crloval = crloval;
  actgrafptr->r.cmloval = cmloval;

  vfixnbr = 0;                                    /* Assume there are no fixed vertices              */
  if (pfixtax != NULL) {                          /* If we have fixed vertices, compute their number */
    Gnum                vertnum;

    for (vertnum = srcgrafptr->baseval; vertnum < srcgrafptr->vertnnd; vertnum ++)
      if (actgrafptr->pfixtax[vertnum] != -1)
        vfixnbr ++;
  }
  actgrafptr->vfixnbr = vfixnbr;
  
  if (memAllocGroup ((void **) (void *)           /* Allocation and initialization of imbalance arrays */
                     &actgrafptr->frontab,     (size_t) (actgrafptr->s.vertnbr * sizeof (Gnum)),
                     &actgrafptr->comploadavg, (size_t) (actgrafptr->m.domnmax * sizeof (Gnum)), /* TRICK: can send both compload arrays in one piece */
                     &actgrafptr->comploaddlt, (size_t) (actgrafptr->m.domnmax * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("kgraphInit: out of memory (2)");
    return     (1);
  }

  actgrafptr->s.flagval |= KGRAPHFREECOMP | KGRAPHFREEFRON | KGRAPHFREERMAP | KGRAPHFRCOGROUP; /* comploadavg and comploaddlt are always grouped */

  domfrstload = archDomWght (archptr, &domfrst);  /* Get load of whole architecture */

  actgrafptr->comploadavg[0] = actgrafptr->s.velosum;
  actgrafptr->comploaddlt[0] = 0;

  actgrafptr->fronnbr     = 0;                    /* No frontier yet */
  actgrafptr->comploadrat = (double) srcgrafptr->velosum / (double) domfrstload; /* Always balance with respect to original graph */
  actgrafptr->commload    = 0;
  actgrafptr->levlnum     = 0;
  actgrafptr->kbalval     = 1;                    /* No information on imbalance yet */

  return (0);
}

/* This routine frees the contents
** of the given active graph and
** updates the mapping data accordingly.
** It returns:
** - VOID  : in all cases.
*/

void
kgraphExit (
Kgraph * restrict const     grafptr)
{
  mapExit (&grafptr->m);
  if ((grafptr->r.m.parttax != NULL) &&
      ((grafptr->s.flagval & KGRAPHFREERMAP) != 0))
    mapExit (&grafptr->r.m);
  kgraphFree (grafptr);                           /* Free kgraph data */
  graphExit (&grafptr->s);

#ifdef SCOTCH_DEBUG_KGRAPH2
  memSet (grafptr, ~0, sizeof (Kgraph));          /* Purge kgraph fields */
#endif /* SCOTCH_DEBUG_KGRAPH2 */
}

/* This routine frees the graph data and
** updates the mapping data accordingly.
** It returns:
** - VOID  : in all cases.
*/

void
kgraphFree (
Kgraph * restrict const     grafptr)
{
  if (((grafptr->s.flagval & KGRAPHFREEFRON) != 0) && /* If frontab must be freed */
      (grafptr->frontab != NULL))                 /* And if it exists             */
    memFree (grafptr->frontab);                   /* Free it                      */

  if (((grafptr->s.flagval & KGRAPHFREECOMP) != 0) && /* If comptabs must be freed */
      (grafptr->comploadavg != NULL) &&           /* And if it exists              */
      ((grafptr->s.flagval & KGRAPHFRCOGROUP) == 0)) /* And not grouped            */
    memFree (grafptr->comploadavg);               /* Free it                       */

  if (((grafptr->s.flagval & KGRAPHFREEVMLO) != 0) && /* If vmlotax must be freed */ 
      (grafptr->r.vmlotax != NULL))               /* And if it exists             */
    memFree (grafptr->r.vmlotax + grafptr->s.baseval); /* Free it                 */
}

/* This routine moves all of the graph
** vertices to the first subdomain, and
** computes the resulting gains.
** It returns:
** - VOID  : in all cases.
*/

void
kgraphFrst (
Kgraph * restrict const     grafptr)
{
  archDomFrst (grafptr->m.archptr, &grafptr->m.domntab[0]);
  grafptr->m.domnnbr = 1;

  memSet (grafptr->m.parttax + grafptr->s.baseval, 0, grafptr->s.vertnbr * sizeof (Anum)); /* Set all vertices to subdomain 0 */
  memSet (grafptr->comploadavg + 1, 0, (2 * grafptr->m.domnmax - 1) * sizeof (Gnum));

  grafptr->comploadavg[0] = grafptr->s.velosum;
  grafptr->commload       = 0;
  grafptr->m.domntab[0]   = grafptr->m.domnorg;   /* Point to first domain */
  grafptr->fronnbr        = 0;                    /* No frontier vertices  */
}

/* This routine computes the cost of the
** current partition.
** It returns:
** - VOID  : in all cases.
*/

void
kgraphCost (
Kgraph * restrict const     grafptr)
{
  Gnum                          vertnum;
  Gnum * restrict               compload;
  Gnum                          commload;
  double                        fdomwgt;
  Gnum                          fvelsum;
  Gnum                          velosum;
  Anum                          domnnum;
  ArchDom                       domndat;
  double                        domnrat;

  const Arch * restrict const     archptr = &grafptr->a;
  const ArchDom * restrict const  domntab = grafptr->m.domntab;
  Anum * restrict const           parttax = grafptr->m.parttax;
  const Gnum * restrict const     verttax = grafptr->s.verttax;
  const Gnum * restrict const     velotax = grafptr->s.velotax;
  const Gnum * restrict const     vendtax = grafptr->s.vendtax;
  const Gnum * restrict const     edlotax = grafptr->s.edlotax;
  const Gnum * restrict const     edgetax = grafptr->s.edgetax;
  const Anum                      domnnbr = grafptr->m.domnnbr;

  commload = 0;
  compload = grafptr->comploaddlt;                   /* Use delta array as temporary storage */
  memSet (compload, 0, domnnbr * sizeof (Gnum));
  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
    Gnum                edgenum;
    Gnum                edgennd;
    Anum                partval;                  /* Part of current vertex                                */
    Anum                partlst;                  /* Part of last vertex for which a distance was computed */
    Anum                distlst;                  /* Last distance computed                                */
    Gnum                veloval;

    partval = parttax[vertnum];
    partlst = -1;                                 /* Invalid part to recompute distance */
    distlst = -1;                                 /* To prevent compiler from yielding  */

#ifdef SCOTCH_DEBUG_KGRAPH2
    if ((partval < 0) || (partval >= domnnbr)) {
      errorPrint ("kgraphCost: invalid parameters (1)");
      return;
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    veloval = (velotax != NULL) ? velotax[vertnum] : 1;
    compload[partval] += veloval;

    for (edgenum = verttax[vertnum], edgennd = vendtax[vertnum];
         edgenum < edgennd; edgenum ++) {
      Gnum                vertend;
      Anum                partend;

      vertend = edgetax[edgenum];
      if (vertend > vertnum)                      /* Compute loads only once */
        continue;

      partend = parttax[vertend];
#ifdef SCOTCH_DEBUG_KGRAPH2
      if ((partend < 0) || (partend >= domnnbr)) {
        errorPrint ("kgraphCost: invalid parameters (2)");
        return;
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      if (partval != partend) {
        Anum                distval;

        distval = (partend != partlst) ? archDomDist (archptr, &domntab[partval], &domntab[partend]) : distlst;
        distlst = distval;
        partlst = partend;
          
        commload += (Gnum) distval * ((edlotax != NULL) ? edlotax[edgenum] : 1);
      }
    }
  }
  grafptr->commload = commload;

  fdomwgt = 0;
  fvelsum = 0;
  if ((grafptr->s.flagval & KGRAPHHASANCHORS) != 0) {
    const Gnum                  vertancnnd = grafptr->s.vertnnd - domnnbr;
    Gnum                        veloval;

    for (domnnum = 0; domnnum < domnnbr; domnnum ++)
      if ((grafptr->s.verttax[vertancnnd + domnnum + 1] - grafptr->s.verttax[vertancnnd + domnnum]) != 0)
        continue;

    if (domnnum != domnnbr) {
      for (domnnum = 0; domnnum < domnnbr; domnnum ++) {
        if ((grafptr->s.verttax[vertancnnd + domnnum + 1] - grafptr->s.verttax[vertancnnd + domnnum]) == 0) {
          veloval = grafptr->s.velotax[vertancnnd + domnnum];
  
          fdomwgt += (double) archDomWght (archptr, &grafptr->m.domntab[domnnum]);
          fvelsum += veloval;
          compload[domnnum] -= 
          grafptr->comploadavg[domnnum] = veloval;
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (compload[domnnum] != 0) {
            errorPrint ("kgraphCost: invalid comploaddlt");
            return;
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
        }
      }
    }
  }
  archDomFrst (archptr, &domndat);
  domnrat = (double) archDomWght (archptr, &domndat);
  domnrat -= fdomwgt;
  velosum = grafptr->s.velosum - fvelsum;
  for (domnnum = 0; domnnum < domnnbr; domnnum ++) {
    compload[domnnum]            -=
    grafptr->comploadavg[domnnum] = (Gnum) ((double) velosum * ((double) archDomWght (archptr, &grafptr->m.domntab[domnnum]) / domnrat));
  }
}

/* This routine computes the frontier
** array of the current partition.
** It returns:
** - VOID  : in all cases.
*/

void
kgraphFron (
Kgraph * restrict const     grafptr)
{
  Gnum                          vertnum;
  Gnum                          vertnnd;
  Gnum                          fronnbr;

  const Gnum * restrict const verttax = grafptr->s.verttax;
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  const Anum * restrict const parttax = grafptr->m.parttax;
  Gnum * restrict const       frontab = grafptr->frontab;

  for (vertnum = grafptr->s.baseval, vertnnd = grafptr->s.vertnnd, fronnbr = 0;
       vertnum < vertnnd; vertnum ++) {
    Anum                partval;                  /* Part of current vertex */
    Gnum                edgenum;                  /* Number of current edge */

    partval = parttax[vertnum];
    for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
      Anum                partend;

      partend = parttax[edgetax[edgenum]];
      if (partend != partval) {                   /* If neighbor belongs to another part */
        frontab[fronnbr ++] = vertnum;
        break;
      }
    }
  }

  grafptr->fronnbr = fronnbr;
}
