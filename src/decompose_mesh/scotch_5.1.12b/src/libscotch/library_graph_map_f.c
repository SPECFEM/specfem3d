/* Copyright 2004,2007,2010,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : library_graph_map_f.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the Fortran API for the  **/
/**                mapping routines of the libSCOTCH       **/
/**                library.                                **/
/**                                                        **/
/**   DATES      : # Version 3.4  : from : 02 dec 1999     **/
/**                                 to     15 nov 2001     **/
/**                # Version 4.0  : from : 12 jan 2004     **/
/**                                 to     12 dec 2005     **/
/**                # Version 5.1  : from : 27 mar 2010     **/
/**                                 to     31 aug 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "scotch.h"

/**************************************/
/*                                    */
/* These routines are the Fortran API */
/* for the mapping routines.          */
/*                                    */
/**************************************/

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHMAPINIT, scotchfgraphmapinit, (     \
const SCOTCH_Graph * const  grafptr,            \
SCOTCH_Mapping * const      mapptr,             \
const SCOTCH_Arch * const   archptr,            \
SCOTCH_Num * const          maptab,             \
int * const                 revaptr),           \
(grafptr, mapptr, archptr, maptab, revaptr))
{
  *revaptr = SCOTCH_graphMapInit (grafptr, mapptr, archptr, maptab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHMAPEXIT, scotchfgraphmapexit, (     \
const SCOTCH_Graph * const  grafptr,            \
SCOTCH_Mapping * const      mapptr),            \
(grafptr, mapptr))
{
  SCOTCH_graphMapExit (grafptr, mapptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHMAPLOAD, scotchfgraphmapload, (     \
const SCOTCH_Graph * const  grafptr,            \
SCOTCH_Mapping * const      mapptr,             \
int * const                 fileptr,            \
int * const                 revaptr),           \
(grafptr, mapptr, fileptr, revaptr))
{
  FILE *              stream;                     /* Stream to build from handle */
  int                 filenum;                    /* Duplicated handle           */
  int                 o;

  if ((filenum = dup (*fileptr)) < 0) {           /* If cannot duplicate file descriptor */
    errorPrint ("SCOTCHFGRAPHMAPLOAD: cannot duplicate handle");

    *revaptr = 1;                                 /* Indicate error */
    return;
  }
  if ((stream = fdopen (filenum, "r")) == NULL) { /* Build stream from handle */
    errorPrint ("SCOTCHFGRAPHMAPLOAD: cannot open input stream");
    close      (filenum);
    *revaptr = 1;
    return;
  }

  o = SCOTCH_graphMapLoad (grafptr, mapptr, stream);

  fclose (stream);                                /* This closes filenum too */

  *revaptr = o;
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHMAPSAVE, scotchfgraphmapsave, (     \
const SCOTCH_Graph * const  grafptr,            \
SCOTCH_Mapping * const      mapptr,             \
int * const                 fileptr,            \
int * const                 revaptr),           \
(grafptr, mapptr, fileptr, revaptr))
{
  FILE *              stream;                     /* Stream to build from handle */
  int                 filenum;                    /* Duplicated handle           */
  int                 o;

  if ((filenum = dup (*fileptr)) < 0) {           /* If cannot duplicate file descriptor */
    errorPrint ("SCOTCHFGRAPHMAPSAVE: cannot duplicate handle");

    *revaptr = 1;                                 /* Indicate error */
    return;
  }
  if ((stream = fdopen (filenum, "w")) == NULL) { /* Build stream from handle */
    errorPrint ("SCOTCHFGRAPHMAPSAVE: cannot open output stream");
    close      (filenum);
    *revaptr = 1;
    return;
  }

  o = SCOTCH_graphMapSave (grafptr, mapptr, stream);

  fclose (stream);                                /* This closes filenum too */

  *revaptr = o;
}

/*
**
*/

FORTRAN (                                         \
SCOTCHFGRAPHMAPCOMPUTE, scotchfgraphmapcompute, ( \
SCOTCH_Graph * const        grafptr,              \
SCOTCH_Mapping * const      mapptr,               \
SCOTCH_Strat * const        stratptr,             \
int * const                 revaptr),             \
(grafptr, mapptr, stratptr, revaptr))
{
  *revaptr = SCOTCH_graphMapCompute (grafptr, mapptr, stratptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHMAP, scotchfgraphmap, (             \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Arch * const   archptr,            \
SCOTCH_Strat * const        stratptr,           \
SCOTCH_Num * const          maptab,             \
int * const                 revaptr),           \
(grafptr, archptr, stratptr, maptab, revaptr))
{
  *revaptr = SCOTCH_graphMap (grafptr, archptr, stratptr, maptab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHPART, scotchfgraphpart, (           \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Num * const    partptr,            \
SCOTCH_Strat * const        stratptr,           \
SCOTCH_Num * const          maptab,             \
int * const                 revaptr),           \
(grafptr, partptr, stratptr, maptab, revaptr))
{
  *revaptr = SCOTCH_graphPart (grafptr, *partptr, stratptr, maptab);
}

/* String lengths are passed at the very
** end of the argument list.
*/

FORTRAN (                                       \
SCOTCHFSTRATGRAPHMAP, scotchfstratgraphmap, (   \
SCOTCH_Strat * const        stratptr,           \
const char * const          string,             \
int * const                 revaptr,            \
const int                   strnbr),            \
(stratptr, string, revaptr, strnbr))
{
  char * restrict     strtab;                     /* Pointer to null-terminated string */

  if ((strtab = (char *) memAlloc (strnbr + 1)) == NULL) { /* Allocate temporary space */
    errorPrint ("SCOTCHFSTRATGRAPHMAP: out of memory (1)");
    *revaptr = 1;
  }
  memCpy (strtab, string, strnbr);                /* Copy string contents */
  strtab[strnbr] = '\0';                          /* Terminate string     */

  *revaptr = SCOTCH_stratGraphMap (stratptr, strtab); /* Call original routine */

  memFree (strtab);
}

/*
**
*/

FORTRAN (                                               \
SCOTCHFSTRATGRAPHMAPBUILD, scotchfstratgraphmapbuild, ( \
SCOTCH_Strat * const        stratptr,                   \
const SCOTCH_Num * const    flagval,                    \
const SCOTCH_Num * const    partnbr,                    \
const double * const        balrat,                     \
int * const                 revaptr),                   \
(stratptr, flagval, partnbr, balrat, revaptr))
{
  *revaptr = SCOTCH_stratGraphMapBuild (stratptr, *flagval, *partnbr, *balrat);
}

/*
**
*/

FORTRAN (                                                       \
SCOTCHFSTRATGRAPHCLUSTERBUILD, scotchfstratgraphclusterbuild, ( \
SCOTCH_Strat * const        stratptr,                           \
const SCOTCH_Num * const    flagval,                            \
const SCOTCH_Num * const    pwgtval,                            \
const double * const        densval,                            \
const double * const        bbalval,                            \
int * const                 revaptr),                           \
(stratptr, flagval, pwgtval, densval, bbalval, revaptr))
{
  *revaptr = SCOTCH_stratGraphClusterBuild (stratptr, *flagval, *pwgtval, *densval, *bbalval);
}
