/* Copyright 2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : library_graph_coarsen.c                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the graph    **/
/**                coarsening routine of the libSCOTCH     **/
/**                library.                                **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 07 aug 2011     **/
/**                                 to     07 aug 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "graph.h"
#include "graph_coarsen.h"
#include "scotch.h"

/*********************************/
/*                               */
/* This routine is the C API for */
/* the graph coarsening routine. */
/*                               */
/*********************************/

/*+ This routine creates a coarse graph from the
*** given fine graph, unless the coarse graph is
*** smaller than some threshold size or the
*** coarsening ratio is above some other threshold.
*** If the coarse graph is created, a coarse-to-fine
*** vertex array is created, that contains a pair of
*** fine indices for each coarse index. The contents
*** of the Scotch internal array are copied to the
*** array provided by the user.
*** It returns:
*** - 0  : if the graph has been coarsened.
*** - 1  : if the graph could not be coarsened.
*** - 2  : on error.
+*/

int
SCOTCH_graphCoarsen (
const SCOTCH_Graph * restrict const finegrafptr,  /* Fine graph structure to fill      */
SCOTCH_Graph * restrict const       coargrafptr,  /* Coarse graph                      */
SCOTCH_Num * restrict const         coarmulttab,  /* Pointer to multinode array        */
const SCOTCH_Num                    coarnbr,      /* Minimum number of coarse vertices */
const double                        coarrat)      /* Maximum contraction ratio         */
{
  GraphCoarsenMulti * restrict  coarmultptr;
  int                           o;

  o = graphCoarsen ((const Graph * restrict const) finegrafptr, (Graph * restrict const) coargrafptr,
                    &coarmultptr, coarnbr, coarrat, GRAPHCOARHEM);

  if (o == 0) {                                   /* If coarsening succeeded */
    SCOTCH_Num          coarvertnbr;

    SCOTCH_graphSize (coargrafptr, &coarvertnbr, NULL); /* Get number of coarse vertices */
    memCpy (coarmulttab, coarmultptr, coarvertnbr * 2 * sizeof (Gnum));
  }

  return (o);
}
