/* Copyright 2004,2007,2008,2010,2011 ENSEIRB, INRIA & CNRS
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
Kgraph * restrict const         actgrafptr,       /* Active graph */
const Graph * restrict const    srcgrafptr,       /* Source graph */
const Mapping * restrict const  mappptr)          /* Mapping      */
{
  const Arch * restrict archptr;                  /* Pointer to mapping architecture */
  ArchDom               domfrst;                  /* First, largest domain           */
  Gnum                  domfrstload;              /* Load of first domain            */
  Anum                  termnum;

  actgrafptr->s          = *srcgrafptr;           /* Clone source graph   */
  actgrafptr->s.flagval &= ~GRAPHFREETABS;        /* Do not allow freeing */

  if (mappptr != &actgrafptr->m)
    actgrafptr->m = *mappptr;                     /* Clone mapping */

  if ((actgrafptr->comploadavg = (Gnum *) memAlloc (actgrafptr->m.domnmax * sizeof (Gnum) * 2)) == NULL) {
    errorPrint ("kgraphInit: out of memory");
    return     (1);
  }
  actgrafptr->comploaddlt = actgrafptr->comploadavg + actgrafptr->m.domnnbr; /* TRICK: can send both arrays in one piece */

  archptr = &mappptr->archdat;
  archDomFrst (archptr, &domfrst);                /* Get first, largest domain */
  domfrstload = archDomWght (archptr, &domfrst);  /* Get its load              */

  actgrafptr->comploadavg[0] = (archDomWght (archptr, &actgrafptr->m.domntab[0]) * actgrafptr->s.velosum) / domfrstload;
  actgrafptr->comploaddlt[0] = actgrafptr->s.velosum - actgrafptr->comploadavg[0];
  for (termnum = 1; termnum < actgrafptr->m.domnnbr; termnum ++) {
    actgrafptr->comploadavg[termnum] = (archDomWght (archptr, &actgrafptr->m.domntab[termnum]) * actgrafptr->s.velosum) / domfrstload;
    actgrafptr->comploaddlt[termnum] = - actgrafptr->comploadavg[termnum];
  }

  actgrafptr->fronnbr     = 0;                    /* No frontier yet */
  actgrafptr->frontab     = NULL;
  actgrafptr->comploadrat = (double) actgrafptr->s.velosum / (double) domfrstload;
  actgrafptr->commload    = 0;
  actgrafptr->levlnum     = 0;

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
  if ((grafptr->m.parttax != NULL) &&
      ((grafptr->s.flagval & KGRAPHFREEPART) != 0))
    memFree (grafptr->m.parttax + grafptr->m.baseval);
  if (grafptr->frontab != NULL)                   /* Free frontier if it exists */
    memFree (grafptr->frontab);
  if (grafptr->comploadavg != NULL)               /* Free load array if it exists */
    memFree (grafptr->comploadavg);

  graphExit (&grafptr->s);
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
  memSet (grafptr->m.parttax + grafptr->m.baseval, 0, grafptr->m.vertnbr * sizeof (Anum)); /* Set all vertices to subdomain 0 */

  grafptr->m.domntab[0] = grafptr->m.domnorg;     /* Point to first domain */
  grafptr->fronnbr      = 0;                      /* No frontier vertices  */
}
