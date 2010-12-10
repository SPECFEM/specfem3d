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
/**   NAME       : kgraph_check.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the mapping graph  **/
/**                consistency checking routine.           **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 13 jul 2010     **/
/**                                 to     13 jul 2010     **/
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

/*************************/
/*                       */
/* These routines handle */
/* mapping graphs.       */
/*                       */
/*************************/

/* This routine checks the consistency
** of the given mapping graph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
kgraphCheck (
const Kgraph * restrict const grafptr)
{
  int * restrict      flagtax;                    /* Frontier flag array       */
  Gnum                vertnum;                    /* Number of current vertex  */
  Gnum                fronnum;                    /* Number of frontier vertex */

  const Gnum * restrict const verttax = grafptr->s.verttax;
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  const Anum * restrict const parttax = grafptr->m.parttax;

  if ((flagtax = memAlloc (grafptr->s.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("kgraphCheck: out of memory");
    return     (1);
  }
  memSet (flagtax, ~0, grafptr->s.vertnbr * sizeof (Gnum));
  flagtax -= grafptr->s.baseval;

  if ((grafptr->m.domnmax <= 0) || (grafptr->m.domnnbr <= 0) || (grafptr->m.domnnbr > grafptr->m.domnmax)) { 
    errorPrint ("kgraphCheck: invalid number of domains");
    return     (1);
  }

  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
    if ((parttax[vertnum] < 0) ||
        (parttax[vertnum] >= grafptr->m.domnnbr)) {
      errorPrint ("kgraphCheck: invalid part array");
      return     (1);
    }
  }

  if ((grafptr->fronnbr < 0) ||
      (grafptr->fronnbr > grafptr->s.vertnbr)) {
    errorPrint ("kgraphCheck: invalid number of frontier vertices");
    return     (1);
  }
  for (fronnum = 0; fronnum < grafptr->fronnbr; fronnum ++) {
    Gnum                vertnum;
    Gnum                edgenum;
    Anum                partval;
    Anum                flagval;

    vertnum = grafptr->frontab[fronnum];
    if ((vertnum < grafptr->s.baseval) || (vertnum >= grafptr->s.vertnnd)) {
      errorPrint ("kgraphCheck: invalid vertex index in frontier array");
      return     (1);
    }
    if (flagtax[vertnum] != ~0) {
      errorPrint ("kgraphCheck: duplicate vertex in frontier array");
      return     (1);
    }
    flagtax[vertnum] = 0;
    partval = parttax[vertnum];

    for (edgenum = verttax[vertnum], flagval = 0;
         edgenum < vendtax[vertnum]; edgenum ++)
      flagval |= parttax[edgetax[edgenum]] ^ partval; /* Flag set if neighbor part differs from vertex part */

    if (flagval == 0) {
      errorPrint ("kgraphCheck: invalid vertex in frontier array");
      return     (1);
    }
  }

  memFree (flagtax + grafptr->s.baseval);

  return (0);
}
