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
/**   NAME       : bgraph_bipart_gg.h                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Luca SCARANO (v3.1)                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the function       **/
/**                declarations for the greedy graph       **/
/**                growing bipartitioning method.          **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 07 jan 1996     **/
/**                                 to     29 apr 1996     **/
/**                # Version 3.2  : from : 20 sep 1996     **/
/**                                 to     13 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 09 jan 2004     **/
/**                                 to     09 jan 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ System-defined constants. +*/

#define BGRAPHBIPARTGGGAINTABLSUBBITS 1

#define BGRAPHBIPARTGGSTATEFREE     ((GainLink *) 0) /*+ Vertex in initial state (TRICK: must be 0) +*/
#define BGRAPHBIPARTGGSTATEUSED     ((GainLink *) 1) /*+ Swapped vertex                             +*/
#define BGRAPHBIPARTGGSTATELINK     ((GainLink *) 2) /*+ Currently in gain table if higher          +*/

/*
**  The type and structure definitions.
*/

/** Method parameters. **/

typedef struct BgraphBipartGgParam_ {
  INT                       passnbr;              /*+ Number of passes to do +*/
} BgraphBipartGgParam;

/*+ The complementary vertex structure. For
    trick reasons, the gain table data structure
    must be the first field of the structure.    +*/

typedef struct BgraphBipartGgVertex_ {
  GainLink                  gainlink;             /*+ Gain link: FIRST                       +*/
  Gnum                      commgain0;            /*+ Gain if vertex and neighbors in part 0 +*/
  Gnum                      commgain;             /*+ Gain value                             +*/
} BgraphBipartGgVertex;

/*
**  The function prototypes.
*/

#ifndef BGRAPH_BIPART_GG
#define static
#endif

int                         bgraphBipartGg      (Bgraph * restrict const, const BgraphBipartGgParam * const);

#undef static
