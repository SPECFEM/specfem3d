/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : graph_coarsen.h                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the source graph coarsening         **/
/**                functions.                              **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 02 dec 1992     **/
/**                                 to     18 may 1993     **/
/**                # Version 1.3  : from : 30 apr 1994     **/
/**                                 to     18 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     18 aug 1994     **/
/**                # Version 3.0  : from : 07 jul 1995     **/
/**                                 to     28 sep 1995     **/
/**                # Version 3.1  : from : 28 nov 1995     **/
/**                                 to     28 nov 1995     **/
/**                # Version 3.2  : from : 07 sep 1996     **/
/**                                 to     17 sep 1998     **/
/**                # Version 4.0  : from : 13 dec 2001     **/
/**                                 to     05 dec 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/** Prime number for cache-friendly perturbations. **/

#define GRAPHCOARPERTPRIME          179           /* Prime number */

/** Prime number for hashing vertex numbers. **/

#define GRAPHCOARHASHPRIME          1049          /* Prime number */

/*
**  The type and structure definitions.
*/

/*+ Here are the edge matching function types for coarsening. +*/

typedef enum GraphCoarsenType_ {
  GRAPHCOARHEM,                                   /*+ Heavy-edge matching         +*/
  GRAPHCOARSCN,                                   /*+ Scanning (first) matching   +*/
  GRAPHCOARCSC,                                   /*+ Crystal scanning matching   +*/
  GRAPHCOARCHE,                                   /*+ Crystal heavy-edge matching +*/
  GRAPHCOARNBR                                    /*+ Number of matching types    +*/
} GraphCoarsenType;

/*+ The multinode table element, which contains
    pairs of based indices of collapsed vertices.
    Both values are equal for uncollapsed vertices.
    As the base values of the fine and coarse graphs
    may be different, the values of the collapsed
    vertices are set with respect to the base value
    of the fine graph.                               +*/

typedef struct GraphCoarsenMulti_  {
  Gnum                      vertnum[2];           /*+ Numbers of the collapsed vertices of a multinode +*/
} GraphCoarsenMulti;

/*+ A table made of such elements is used during
    coarsening to build the edge array of the new
    graph, after the labeling of the vertices.    +*/

typedef struct GraphCoarsenHash_ {
  Gnum                      vertorgnum;           /*+ Origin vertex (i.e. pass) number +*/
  Gnum                      vertendnum;           /*+ Other end vertex number          +*/
  Gnum                      edgenum;              /*+ Number of corresponding edge     +*/
} GraphCoarsenHash;

/*
**  The function prototypes.
*/

#ifndef GRAPH_COARSEN
#define static
#endif

int                         graphCoarsen        (const Graph * restrict const, Graph * restrict const, GraphCoarsenMulti * restrict * const, const Gnum, const double, const GraphCoarsenType);

static void                 graphCoarsenEdgeLl  (const Graph * const, const Gnum * const, const GraphCoarsenMulti * restrict const, Graph * const, GraphCoarsenHash * const, const Gnum);
static void                 graphCoarsenEdgeLu  (const Graph * const, const Gnum * const, const GraphCoarsenMulti * restrict const, Graph * const, GraphCoarsenHash * const, const Gnum);

static Gnum                 graphCoarsenMatchHy (const Graph * const, Gnum *, const Gnum, const Gnum);
static Gnum                 graphCoarsenMatchSc (const Graph * const, Gnum *, const Gnum, const Gnum);
static Gnum                 graphCoarsenMatchCs (const Graph * const, Gnum *, const Gnum, const Gnum);
static Gnum                 graphCoarsenMatchCh (const Graph * const, Gnum *, const Gnum, const Gnum);

#undef static
