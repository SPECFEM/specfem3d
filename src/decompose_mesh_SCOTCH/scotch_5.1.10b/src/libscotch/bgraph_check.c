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
/**   NAME       : bgraph_check.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the bipartition    **/
/**                graph consistency checking routine.     **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 08 jan 2004     **/
/**                                 to     07 dec 2005     **/
/**                # Version 5.1  : from : 04 oct 2009     **/
/**                                 to     04 oct 2009     **/
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
#include "bgraph.h"

/*************************/
/*                       */
/* These routines handle */
/* bipartition graphs.   */
/*                       */
/*************************/

/* This routine checks the consistency
** of the given bipartition graph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
bgraphCheck (
const Bgraph * restrict const grafptr)
{
  int * restrict      flagtax;                    /* Frontier flag array       */
  Gnum                vertnum;                    /* Number of current vertex  */
  Gnum                fronnum;                    /* Number of frontier vertex */
  Gnum                compload[2];
  Gnum                compsize[2];
  Gnum                commcut[2];
  Gnum                commloadintn;
  Gnum                commloadextn;
  Gnum                commgainextn;
  Gnum                edloval;

  const Gnum * restrict const       verttax = grafptr->s.verttax;
  const Gnum * restrict const       vendtax = grafptr->s.vendtax;
  const Gnum * restrict const       velotax = grafptr->s.velotax;
  const Gnum * restrict const       edgetax = grafptr->s.edgetax;
  const Gnum * restrict const       edlotax = grafptr->s.edlotax;
  const GraphPart * restrict const  parttax = grafptr->parttax;

  if ((flagtax = memAlloc (grafptr->s.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("bgraphCheck: out of memory");
    return     (1);
  }
  memSet (flagtax, ~0, grafptr->s.vertnbr * sizeof (Gnum));
  flagtax -= grafptr->s.baseval;

  if (grafptr->compload0 != (grafptr->compload0avg + grafptr->compload0dlt)) {
    errorPrint ("bgraphCheck: invalid balance");
    return     (1);
  }

  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
    if ((parttax[vertnum] < 0) ||
        (parttax[vertnum] > 1)) {
      errorPrint ("bgraphCheck: invalid part array");
      return     (1);
    }
  }

  if ((grafptr->fronnbr < 0) ||
      (grafptr->fronnbr > grafptr->s.vertnbr)) {
    errorPrint ("bgraphCheck: invalid number of frontier vertices");
    return     (1);
  }
  for (fronnum = 0; fronnum < grafptr->fronnbr; fronnum ++) {
    Gnum                vertnum;
    Gnum                edgenum;
    GraphPart           partval;
    GraphPart           flagval;

    vertnum = grafptr->frontab[fronnum];
    if ((vertnum < grafptr->s.baseval) || (vertnum >= grafptr->s.vertnnd)) {
      errorPrint ("bgraphCheck: invalid vertex index in frontier array");
      return     (1);
    }
    if (flagtax[vertnum] != ~0) {
      errorPrint ("bgraphCheck: duplicate vertex in frontier array");
      return     (1);
    }
    flagtax[vertnum] = 0;
    partval = parttax[vertnum];

    for (edgenum = verttax[vertnum], flagval = 0;
         edgenum < vendtax[vertnum]; edgenum ++)
      flagval |= parttax[edgetax[edgenum]] ^ partval; /* Flag set if neighbor part differs from vertex part */

    if (flagval == 0) {
      errorPrint ("bgraphCheck: invalid vertex in frontier array");
      return     (1);
    }
  }

  compload[0]  =
  compload[1]  = 0;
  compsize[0]  =
  compsize[1]  = 0;
  commloadintn = 0;
  commloadextn = grafptr->commloadextn0;
  commgainextn = 0;
  edloval      = 1;                               /* Assume edges are not weighted */
  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
    Gnum                partval;                  /* Part of current vertex */
    Gnum                edgenum;                  /* Number of current edge */

    partval = (Gnum) parttax[vertnum];
    if (grafptr->veextax != NULL) {
      Gnum                veexval;

      veexval = grafptr->veextax[vertnum];
      commloadextn += veexval * partval;
      commgainextn += veexval * (1 - 2 * partval);
    }

    compload[partval] += (velotax == NULL) ? 1 : velotax[vertnum];
    compsize[partval] ++;

    commcut[0] =
    commcut[1] = 0;
    for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
      int                 partend;
      int                 partdlt;

      if (edlotax != NULL)
        edloval = edlotax[edgenum];
      partend = parttax[edgetax[edgenum]];
      partdlt = partval ^ partend;
      commcut[partend] ++;
      commloadintn += partdlt * edloval * partend; /* Only count loads once, when (partend == 1) */
    }

    if ((commcut[0] != 0) && (commcut[1] != 0) && /* If vertex should be in frontier array */
        (flagtax[vertnum] != 0)) {
      errorPrint ("bgraphCheck: vertex should be in frontier array");
      return     (1);
    }
  }
  if (compsize[0] != grafptr->compsize0) {
    errorPrint ("bgraphCheck: invalid part size");
    return     (1);
  }
  if ((commloadintn * grafptr->domdist + commloadextn) != grafptr->commload) {
    errorPrint ("bgraphCheck: invalid communication loads");
    return     (1);
  }
  if (commgainextn != grafptr->commgainextn) {
    errorPrint ("bgraphCheck: invalid communication gains");
    return     (1);
  }

  memFree (flagtax + grafptr->s.baseval);

  return (0);
}
