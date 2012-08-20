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
/**   NAME       : metis_graph_part.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the compatibility        **/
/**                library for the MeTiS partitioning      **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 08 sep 2006     **/
/**                                 to     07 jun 2007     **/
/**                # Version 5.1  : from : 06 jun 2009     **/
/**                                 to     30 jun 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "metis.h"                                /* Our "metis.h" file */

/************************************/
/*                                  */
/* These routines are the C API for */
/* MeTiS graph ordering routines.   */
/*                                  */
/************************************/

/* This routine is the interface between MeTiS
** and Scotch. It computes the partition of a
** weighted or unweighted graph.
** It returns:
** - 0   : if the partition could be computed.
** - !0  : on error.
*/

static
int
_SCOTCH_METIS_PartGraph (
const int * const           n,
const int * const           xadj,
const int * const           adjncy,
const int * const           vwgt,
const int * const           adjwgt,
const int * const           numflag,
const int * const           nparts,
int * const                 part)
{
  SCOTCH_Graph        grafdat;                    /* Scotch graph object to interface with libScotch */
  SCOTCH_Strat        stradat;
  SCOTCH_Num          baseval;
  SCOTCH_Num          vertnbr;
  int                 o;

  if (sizeof (SCOTCH_Num) != sizeof (int)) {
    errorPrint ("METIS_PartGraph* (as of SCOTCH): SCOTCH_Num type should equate to int");
    return (1);
  }

  SCOTCH_graphInit (&grafdat);

  baseval = *numflag;
  vertnbr = *n;

  o = 1;                                          /* Assume something will go wrong */
  if (SCOTCH_graphBuild (&grafdat,
                         baseval, vertnbr, xadj, xadj + 1, vwgt, NULL,
                         xadj[vertnbr] - baseval, adjncy, adjwgt) == 0) {
    SCOTCH_stratInit (&stradat);
#ifdef SCOTCH_DEBUG_ALL
    if (SCOTCH_graphCheck (&grafdat) == 0)        /* TRICK: next instruction called only if graph is consistent */
#endif /* SCOTCH_DEBUG_ALL */
    o = SCOTCH_graphPart (&grafdat, *nparts, &stradat, part);
    SCOTCH_stratExit (&stradat);
  }
  SCOTCH_graphExit (&grafdat);

  if (baseval != 0) {                             /* MeTiS part array is based, Scotch is not */
    SCOTCH_Num          vertnum;

    for (vertnum = 0; vertnum < vertnbr; vertnum ++)
      part[vertnum] += baseval;
  }

  return (o);
}

/*
**
*/

void
METISNAMEU(METIS_PartGraphKway) (
const int * const           n,
const int * const           xadj,
const int * const           adjncy,
const int * const           vwgt,
const int * const           adjwgt,
const int * const           wgtflag,
const int * const           numflag,
const int * const           nparts,
const int * const           options,
int * const                 edgecut,
int * const                 part)
{
  METISNAMEU(METIS_PartGraphRecursive) (n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}

void
METISNAMEU(METIS_PartGraphRecursive) (
const int * const           n,
const int * const           xadj,
const int * const           adjncy,
const int * const           vwgt,
const int * const           adjwgt,
const int * const           wgtflag,
const int * const           numflag,
const int * const           nparts,
const int * const           options,
int * const                 edgecut,
int * const                 part)
{
  const int *           vwgt2;
  const int *           adjwgt2;
  const int * restrict  parttax;
  const int * restrict  verttax;
  const int * restrict  edgetax;
  int                   vertnnd;
  int                   vertnum;
  int                   edgenum;
  int                   commcut;

  vwgt2   = (((*wgtflag & 2) != 0) ? vwgt   : NULL);
  adjwgt2 = (((*wgtflag & 1) != 0) ? adjwgt : NULL);

  if (_SCOTCH_METIS_PartGraph (n, xadj, adjncy, vwgt2, adjwgt2, numflag, nparts, part) != 0)
    return;

  parttax = part   - *numflag;
  verttax = xadj   - *numflag;
  edgetax = adjncy - *numflag;
  edgenum = *numflag;
  vertnum = *numflag;
  vertnnd = *n + vertnum;
  commcut = 0;

  if (adjwgt2 == NULL) {                          /* If graph does not have edge weights */
    for ( ; vertnum < vertnnd; vertnum ++) {
      int                 edgennd;
      int                 partval;

      partval = parttax[vertnum];
      for (edgennd = verttax[vertnum + 1]; edgenum < edgennd; edgenum ++) {
        if (parttax[edgetax[edgenum]] != partval)
          commcut ++;
      }
    }
  }
  else {                                          /* Graph has edge weights */
    const int * restrict  edlotax;

    edlotax = adjwgt2 - *numflag;
    for ( ; vertnum < vertnnd; vertnum ++) {
      int                 edgennd;
      int                 partval;

      partval = parttax[vertnum];
      for (edgennd = verttax[vertnum + 1]; edgenum < edgennd; edgenum ++) {
        int                 vertend;

        vertend = edgetax[edgenum];
        if (parttax[vertend] != partval)
          commcut += edlotax[edgenum];
      }
    }
  }
  *edgecut = commcut / 2;
}

/* Scotch does not directly consider communication volume.
** Instead, wertex communication loads are added to the edge
** loads so as to emulate this behavior : heavily weighted
** edges, connected to heavily communicating vertices, will
** be less likely to be cut.
*/

void
METISNAMEU(METIS_PartGraphVKway) (
const int * const           n,
const int * const           xadj,
const int * const           adjncy,
const int * const           vwgt,
const int * const           vsize,
const int * const           wgtflag,
const int * const           numflag,
const int * const           nparts,
const int * const           options,
int * const                 volume,
int * const                 part)
{
  int                   baseval;
  const int *           vwgt2;
  const int *           vsize2;
  int                   vsizval;                  /* Communication volume of current vertex */
  int                   vertnbr;
  int                   vertnum;
  int                   edgenum;
  const int * restrict  edgetax;
  const int * restrict  parttax;
  int * restrict        nghbtab;
  int                   commvol;

  vsize2  = ((*wgtflag & 1) != 0) ? vsize : NULL;
  vwgt2   = ((*wgtflag & 2) != 0) ? vwgt  : NULL;
  baseval = *numflag;
  vertnbr = *n;
  edgetax = adjncy - baseval;

  if (vsize2 == NULL)                             /* If no communication load data provided */
    _SCOTCH_METIS_PartGraph (n, xadj, adjncy, vwgt2, NULL, numflag, nparts, part);
  else {                                          /* Will have to turn communication volumes into edge loads */
    const int * restrict  vsiztax;
    int                   edgenbr;
    int * restrict        edlotax;
    int                   o;

    edgenbr = xadj[vertnbr] - baseval;
    if ((edlotax = memAlloc (edgenbr * sizeof (int))) == NULL)
      return;
    edlotax -= baseval;                           /* Base access to edlotax */
    vsiztax  = vsize2 - baseval;

    for (vertnum = 0, edgenum = baseval;          /* Un-based scan of vertex array xadj */
         vertnum < vertnbr; vertnum ++) {
      int                 vsizval;                /* Communication size of current vertex */
      int                 edgennd;

      vsizval = vsize2[vertnum];
      for (edgennd = xadj[vertnum + 1]; edgenum < edgennd; edgenum ++) { /* Based traversal of edge array adjncy */
        int                 vertend;              /* Based end vertex number                                     */

        vertend = edgetax[edgenum];
        edlotax[edgenum] = vsizval + vsiztax[vertend];
      }
    }

    o = _SCOTCH_METIS_PartGraph (n, xadj, adjncy, vwgt2, edlotax + baseval, numflag, nparts, part);

    memFree (edlotax + baseval);

    if (o != 0)
      return;
  }

  if ((nghbtab = memAlloc (*nparts * sizeof (int))) == NULL)
    return;
  memSet (nghbtab, ~0, *nparts * sizeof (int));

  parttax = part - baseval;
  vsizval = 1;                                      /* Assume no vertex communication sizes */
  for (vertnum = 0, edgenum = baseval, commvol = 0; /* Un-based scan of vertex array xadj   */
       vertnum < vertnbr; vertnum ++) {
    int                 partval;
    int                 edgennd;

    partval = part[vertnum];
    nghbtab[partval] = vertnum;                   /* Do not count local neighbors in communication volume */
    if (vsize2 != NULL)
      vsizval = vsize2[vertnum];

    for (edgennd = xadj[vertnum + 1]; edgenum < edgennd; edgenum ++) { /* Based traversal of edge array adjncy */
      int                 vertend;                /* Based end vertex number                                   */
      int                 partend;

      vertend = edgetax[edgenum];
      partend = parttax[vertend];
      if (nghbtab[partend] != vertnum) {          /* If first neighbor in this part */
        nghbtab[partend] = vertnum;               /* Set part as accounted for      */
        commvol += vsizval;
      }
    }
  }
  *volume = commvol;

  memFree (nghbtab);
}
