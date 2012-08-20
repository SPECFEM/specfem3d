/* Copyright 2004,2007,2008,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : bgraph.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the bipartition    **/
/**                graph data structure handling           **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to     12 may 1993     **/
/**                # Version 1.3  : from : 06 apr 1994     **/
/**                                 to     09 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     01 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     30 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     15 aug 1995     **/
/**                # Version 3.1  : from : 15 nov 1995     **/
/**                                 to     16 nov 1995     **/
/**                # Version 3.2  : from : 24 aug 1996     **/
/**                                 to   : 14 oct 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     19 oct 1998     **/
/**                # Version 4.0  : from : 18 dec 2001     **/
/**                                 to     31 aug 2004     **/
/**                # Version 5.0  : from : 17 dec 2006     **/
/**                                 to     10 sep 2007     **/
/**                # Version 5.1  : from : 08 oct 2008     **/
/**                                 to     18 mar 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BGRAPH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "bgraph.h"

/*************************/
/*                       */
/* These routines handle */
/* bipartition graphs.   */
/*                       */
/*************************/

/* This routine builds the active graph
** corresponding to the given bipartitioning
** job parameters.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
bgraphInit (
Bgraph * restrict const         actgrafptr,       /* Active graph                     */
const Graph * restrict const    indgrafptr,       /* Induced source subgraph          */
const Graph * restrict const    srcgrafptr,       /* Original source graph            */
const Mapping * restrict const  mapptr,           /* Current mapping of halo vertices */
const ArchDom                   domsubtab[])      /* Subdomains                       */
{
  Anum                domdist;                    /* Distance between both subdomains   */
  Anum                domwght0;                   /* Processor workforce in each domain */
  Anum                domwght1;

  domdist  = archDomDist (&mapptr->archdat, &domsubtab[0], &domsubtab[1]); /* Get distance between subdomains */
  domwght0 = archDomWght (&mapptr->archdat, &domsubtab[0]); /* Get weights of subdomains                      */
  domwght1 = archDomWght (&mapptr->archdat, &domsubtab[1]);

  actgrafptr->s         = *indgrafptr;            /* Get source graph data */
  actgrafptr->s.flagval = (indgrafptr->flagval & ~GRAPHFREETABS) | BGRAPHFREEFRON | BGRAPHFREEPART; /* Graph is a clone with own grouped bipartitioning arrays */
  actgrafptr->s.vlbltax = NULL;                   /* Remove vertex labels    */
  actgrafptr->veextax   = NULL;                   /* No external gains (yet) */

  if (((actgrafptr->parttax = memAlloc (actgrafptr->s.vertnbr * sizeof (GraphPart))) == NULL) ||
      ((actgrafptr->frontab = memAlloc (actgrafptr->s.vertnbr * sizeof (Gnum)))      == NULL)) {
    errorPrint ("bgraphInit: out of memory");
    if (actgrafptr->parttax != NULL)
      memFree (actgrafptr->parttax);
    return (1);
  }
  actgrafptr->parttax -= actgrafptr->s.baseval;

  bgraphInit2 (actgrafptr, domdist, domwght0, domwght1);

  if ((srcgrafptr != NULL) &&                     /* If target architecture needs external gains and */
      (indgrafptr->vertnbr != srcgrafptr->vertnbr)) {  /* If induced subgraph is not original graph  */
    if (bgraphInit3 (actgrafptr, srcgrafptr, mapptr, domsubtab) != 0) { /* Add external loads        */
      bgraphExit (actgrafptr);
      return     (1);
    }
  }

#ifdef SCOTCH_DEBUG_BGRAPH2
  if (bgraphCheck (actgrafptr) != 0) {
    errorPrint ("bgraphInit: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

  return (0);
}

void
bgraphInit2 (
Bgraph * restrict const         actgrafptr,       /* Active graph                       */
const Anum                      domdist,          /* Distance between both subdomains   */
const Anum                      domwght0,         /* Processor workforce in each domain */
const Anum                      domwght1)
{
  actgrafptr->fronnbr       = 0;                  /* No frontier since all vertices set to part 0 */
  actgrafptr->compload0min  = 0;                  /* No external constraints on bipartition (yet) */
  actgrafptr->compload0max  = actgrafptr->s.velosum;
  actgrafptr->compload0avg  = (Gnum) (((double) actgrafptr->s.velosum * (double) domwght0) / (double) (domwght0 + domwght1));
  actgrafptr->compload0dlt  = actgrafptr->s.velosum - actgrafptr->compload0avg;
  actgrafptr->compload0     = actgrafptr->s.velosum;
  actgrafptr->compsize0     = actgrafptr->s.vertnbr;
  actgrafptr->commload      = 0;
  actgrafptr->commloadextn0 = 0;
  actgrafptr->commgainextn  = 0;
  actgrafptr->commgainextn0 = 0;
  actgrafptr->domdist       = domdist;
  actgrafptr->domwght[0]    = domwght0;
  actgrafptr->domwght[1]    = domwght1;
  actgrafptr->levlnum       = 0;

  memSet (actgrafptr->parttax + actgrafptr->s.baseval, 0, actgrafptr->s.vertnbr * sizeof (GraphPart)); /* Set all vertices to part 0 */
}

/* This routine adds external gain data to
** the active graph given to it, according
** to the initial source graph, the current
** mapping, and the two subdomains.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
bgraphInit3 (
Bgraph * restrict const         actgrafptr,       /*+ Active graph being built +*/
const Graph * restrict const    srcgrafptr,       /*+ Original source graph    +*/
const Mapping * restrict const  mapptr,           /*+ Partial mapping          +*/
const ArchDom                   domsub[])         /*+ Subdomains               +*/
{
  const Arch * restrict tgtarchptr;               /* Pointer to the target architecture */
  Gnum                  actvertnum;               /* Number of current active vertex    */
  Gnum                  commloadextn0;            /* External communication load        */
  Gnum                  commgainextn0;            /* External communication gain        */
  Gnum * restrict       veextax;                  /* External gain array                */
  Gnum                  veexflagval;              /* Flag set if external array useful  */
  
  tgtarchptr = &mapptr->archdat;                  /* Get target architecture */

  if ((veextax = (Gnum *) memAlloc (actgrafptr->s.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("bgraphInit3: out of memory");
    return     (1);
  }
  veextax -= actgrafptr->s.baseval;

  veexflagval   =                                 /* No useful array entry yet     */
  commloadextn0 =                                 /* No external communication yet */
  commgainextn0 = 0;
  for (actvertnum = actgrafptr->s.baseval;        /* Compute external loads */
       actvertnum < actgrafptr->s.vertnnd; actvertnum ++) {
    Gnum                commgainextn;             /* External communication gain for current vertex */
    Gnum                srcvertnum;               /* Number of current original vertex              */

    commgainextn = 0;                             /* Initialize external loads               */
    srcvertnum   = actgrafptr->s.vnumtax[actvertnum]; /* Get vertex number in original graph */

    if ((srcgrafptr->vendtax[srcvertnum] - srcgrafptr->verttax[srcvertnum]) != /* If vertex has external edges */
        (actgrafptr->s.vendtax[actvertnum] - actgrafptr->s.verttax[actvertnum])) {
      Gnum                commloadextn;           /* External communication load for current vertex */
      Gnum                srcedgenum;
      Gnum                srcedloval;

      commloadextn = 0;
      srcedgenum   = srcgrafptr->verttax[srcvertnum];
      srcedloval   = 1;                           /* Assume no edge loads */
      if (actgrafptr->s.vendtax[actvertnum] > actgrafptr->s.verttax[actvertnum]) { /* If vertex has active edges */
        Gnum                actedgenum;
        Gnum                srcvertend;

        for (actedgenum = actgrafptr->s.verttax[actvertnum], srcvertend = actgrafptr->s.vnumtax[actgrafptr->s.edgetax[actedgenum]]; ;
             srcedgenum ++) {
          ArchDom * restrict  srcdomnptr;         /* Pointer to domain of current source edge vertex */

#ifdef SCOTCH_DEBUG_BGRAPH2
          if (srcedgenum >= srcgrafptr->vendtax[srcvertnum]) {
            errorPrint ("bgraphInit3: internal error");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

          if (srcgrafptr->edgetax[srcedgenum] == srcvertend) { /* If source graph edge is current active edge    */
            if (++ actedgenum >= actgrafptr->s.vendtax[actvertnum]) { /* If all active edges found               */
              srcedgenum ++;                      /* Process next edge as external edge in next loop             */
              break;                              /* All remaining external edges will be processed in next loop */
            }
            srcvertend = actgrafptr->s.vnumtax[actgrafptr->s.edgetax[actedgenum]]; /* Get new end number to look for */
            continue;                             /* Skip internal edge                                              */
          }

          srcdomnptr = mapDomain (mapptr, srcgrafptr->edgetax[srcedgenum]);
          if (srcgrafptr->edlotax != NULL)
            srcedloval = srcgrafptr->edlotax[srcedgenum];

          commloadextn += srcedloval * archDomDist (tgtarchptr, domsub,     srcdomnptr);
          commgainextn += srcedloval * archDomDist (tgtarchptr, domsub + 1, srcdomnptr);
        }
      }
      for ( ; srcedgenum < srcgrafptr->vendtax[srcvertnum]; srcedgenum ++) {
        ArchDom * restrict  srcdomnptr;           /* Pointer to domain of current source edge vertex */

        srcdomnptr = mapDomain (mapptr, srcgrafptr->edgetax[srcedgenum]);
        if (srcgrafptr->edlotax != NULL)
          srcedloval = srcgrafptr->edlotax[srcedgenum];

        commloadextn += srcedloval * archDomDist (tgtarchptr, domsub,     srcdomnptr);
        commgainextn += srcedloval * archDomDist (tgtarchptr, domsub + 1, srcdomnptr);
      }

      commgainextn  -= commloadextn;              /* Compute vertex gain        */
      commloadextn0 += commloadextn;              /* Account for external edges */
      commgainextn0 += commgainextn;
    }

    veextax[actvertnum] = commgainextn;           /* Record external gain value */
    veexflagval        |= commgainextn;           /* Accumulate non-zero values */
  }

  if (veexflagval == 0) {                         /* If external gain array is useless */
    memFree (veextax + actgrafptr->s.baseval);    /* Forget about it                   */
    return  (0);
  }

  actgrafptr->s.flagval |= BGRAPHFREEVEEX;        /* Keep external gain array */
  actgrafptr->veextax    = veextax;

  actgrafptr->commload      = commloadextn0;
  actgrafptr->commgainextn  = commgainextn0;
  actgrafptr->commloadextn0 = commloadextn0;
  actgrafptr->commgainextn0 = commgainextn0;

  return (0);
}

/* This routine frees the contents
** of the given active graph.
** It returns:
** - VOID  : in all cases.
*/

void
bgraphExit (
Bgraph * restrict const     grafptr)
{
  if ((grafptr->veextax != NULL) &&               /* External gain array is private */
      ((grafptr->s.flagval & BGRAPHFREEVEEX) != 0))
    memFree (grafptr->veextax + grafptr->s.baseval);
  if ((grafptr->frontab != NULL) &&
      ((grafptr->s.flagval & BGRAPHFREEFRON) != 0))
    memFree (grafptr->frontab);
  if ((grafptr->parttax != NULL) &&
      ((grafptr->s.flagval & BGRAPHFREEPART) != 0))
    memFree (grafptr->parttax + grafptr->s.baseval);

  graphExit (&grafptr->s);                        /* Free re-allocated arrays of cloned source graph, if any */

#ifdef SCOTCH_DEBUG_BGRAPH2
  memSet (grafptr, ~0, sizeof (Bgraph));
#endif /* SCOTCH_DEBUG_BGRAPH2 */
}

/* This routine swaps all of the graph
** vertices from one part to another, and
** recomputes the resulting gains.
** It returns:
** - VOID  : in all cases.
*/

void
bgraphSwal (
Bgraph * restrict const     grafptr)
{
  Gnum                vertnum;

  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++)
    grafptr->parttax[vertnum] ^= 1;

  grafptr->compload0    =   grafptr->s.velosum - grafptr->compload0;
  grafptr->compload0dlt =   grafptr->s.velosum - grafptr->compload0dlt - 2 * grafptr->compload0avg;
  grafptr->compsize0    =   grafptr->s.vertnbr - grafptr->compsize0;
  grafptr->commload    +=   grafptr->commgainextn;
  grafptr->commgainextn = - grafptr->commgainextn;
}

/* This routine moves all of the graph
** vertices to the first part, and
** computes the resulting gains.
** It returns:
** - VOID  : in all cases.
*/

void
bgraphZero (
Bgraph * restrict const     grafptr)
{
  memSet (grafptr->parttax + grafptr->s.baseval, 0, grafptr->s.vertnbr * sizeof (GraphPart)); /* Set all vertices to part 0 */

  grafptr->fronnbr      = 0;                      /* No frontier vertices */
  grafptr->compload0    = grafptr->s.velosum;
  grafptr->compload0dlt = grafptr->s.velosum - grafptr->compload0avg;
  grafptr->compsize0    = grafptr->s.vertnbr;
  grafptr->commload     = grafptr->commloadextn0; /* Initialize communication load */
  grafptr->commgainextn = grafptr->commgainextn0;
  grafptr->bbalval      = (double) grafptr->compload0dlt / (double) grafptr->compload0avg;
}
