/* Copyright 2004,2007,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : bgraph_bipart_fm.h                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for our Improved Fiduccia-Mattheyses    **/
/**                bipartitioning algorithm.               **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 30 sep 1993     **/
/**                                 to     09 oct 1993     **/
/**                # Version 1.1  : from : 15 oct 1993     **/
/**                                 to     15 oct 1993     **/
/**                # Version 1.3  : from : 06 apr 1994     **/
/**                                 to     13 apr 1994     **/
/**                # Version 2.0  : from : 04 jul 1994     **/
/**                                 to     25 nov 1994     **/
/**                # Version 3.0  : from : 06 jul 1995     **/
/**                                 to     06 jul 1995     **/
/**                # Version 3.1  : from : 06 nov 1995     **/
/**                                 to     07 jun 1996     **/
/**                # Version 3.2  : from : 21 sep 1996     **/
/**                                 to     13 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     12 mar 1999     **/
/**                # Version 4.0  : from : 27 aug 2004     **/
/**                                 to     27 aug 2004     **/
/**                # Version 5.1  : from : 27 mar 2011     **/
/**                                 to     27 mar 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Gain table subbits. +*/

#define BGRAPHBIPARTFMSUBBITS       4

/*+ Prime number for hashing vertex numbers. +*/

#define BGRAPHBIPARTFMHASHPRIME     17            /*+ Prime number for hashing +*/

/** Gain table vertex status. **/

#define BGRAPHBIPARTFMSTATEFREE     ((GainLink *) 0) /*+ Vertex in initial state           +*/
#define BGRAPHBIPARTFMSTATEUSED     ((GainLink *) 1) /*+ Swapped vertex                    +*/
#define BGRAPHBIPARTFMSTATELINK     ((GainLink *) 2) /*+ Currently in gain table if higher +*/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct BgraphBipartFmParam_ {
  INT                       movenbr;              /*+ Maximum number of uneffective moves that can be done +*/
  INT                       passnbr;              /*+ Number of passes to be performed (-1 : infinite)     +*/
  double                    deltval;              /*+ Maximum weight imbalance ratio                       +*/
} BgraphBipartFmParam;

/*+ The hash vertex structure. For trick
    reasons, the gain table data structure
    must be the first field of the structure. +*/

typedef struct BgraphBipartFmVertex_ {
  GainLink                  gainlink;             /*+ Gain link: FIRST                        +*/
  Gnum                      vertnum;              /*+ Number of vertex                        +*/
  int                       partval;              /*+ Vertex part                             +*/
  Gnum                      compgain;             /*+ Computation gain                        +*/
  Gnum                      commgain;             /*+ Communication gain                      +*/
  Gnum                      commcut;              /*+ Cut edges                               +*/
  Gnum                      mswpnum;              /*+ Number of move sweep when data recorded +*/
} BgraphBipartFmVertex;

/*+ The move recording structure. +*/

typedef struct BgraphBipartFmSave_ {
  Gnum                      hashnum;              /*+ Number of vertex slot +*/
  int                       partval;              /*+ Vertex part           +*/
  Gnum                      compgain;             /*+ Computation gain      +*/
  Gnum                      commgain;             /*+ Communication gain    +*/
  Gnum                      commcut;              /*+ Cut edges             +*/
} BgraphBipartFmSave;

/*
**  The function prototypes.
*/

#ifndef BGRAPH_BIPART_FM
#define static
#endif

static BgraphBipartFmVertex * bgraphBipartFmTablGet (GainTabl * restrict const, const Gnum, const Gnum, const Gnum);

int                         bgraphBipartFm      (Bgraph * restrict const, const BgraphBipartFmParam * const);

static int                  bgraphBipartFmResize (BgraphBipartFmVertex * restrict *, Gnum * restrict const, Gnum * const, BgraphBipartFmSave * restrict *, const Gnum, GainTabl * const, BgraphBipartFmVertex ** const);
#ifdef SCOTCH_DEBUG_BGRAPH3
static int                  bgraphBipartFmCheck (const Bgraph * restrict const, const BgraphBipartFmVertex * restrict const, const Gnum, const int, const Gnum, const Gnum, const Gnum);
#endif /* SCOTCH_DEBUG_BGRAPH3 */

#undef static
