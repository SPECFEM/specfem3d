/* Copyright 2004,2007-2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : mapping.c                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles (partial) mappings. **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 31 mar 1993     **/
/**                                 to     31 mar 1993     **/
/**                # Version 1.0  : from : 04 oct 1993     **/
/**                                 to     06 oct 1993     **/
/**                # Version 1.1  : from : 15 oct 1993     **/
/**                                 to     15 oct 1993     **/
/**                # Version 1.3  : from : 09 apr 1994     **/
/**                                 to     11 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     17 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     18 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     19 oct 1995     **/
/**                # Version 3.1  : from : 30 oct 1995     **/
/**                                 to     14 jun 1996     **/
/**                # Version 3.2  : from : 23 aug 1996     **/
/**                                 to     07 sep 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to     30 mar 1999     **/
/**                # Version 3.4  : from : 11 sep 2001     **/
/**                                 to     08 nov 2001     **/
/**                # Version 4.0  : from : 16 jan 2004     **/
/**                                 to     05 jan 2005     **/
/**                # Version 5.1  : from : 25 jun 2008     **/
/**                                 to     28 apr 2009     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MAPPING

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"

/***********************************/
/*                                 */
/* These routines handle mappings. */
/*                                 */
/***********************************/

/* This routine builds a mapping.
** It returns:
** - 0   : if mapping successfully initialized.
** - !0  : on error.
*/

int
mapInit2 (
Mapping * restrict const        mappptr,          /*+ Mapping structure to fill      +*/
const Gnum                      baseval,          /*+ Base value                     +*/
const Gnum                      vertnbr,
const Arch * restrict const     archptr,
const ArchDom * restrict const  domnptr)          /*+ Pointer to initial (sub)domain +*/
{
  Anum                domnmax;                    /* Maximum number of domains       */
  Anum * restrict     parttab;                    /* Temporary pointer to part array */

  if (archVar (archptr))                          /* If target architecture is variable-sized */
    domnmax = (vertnbr > 1024) ? 1024 : vertnbr;  /* Pre-set number of domains                */
  else                                            /* Else if fixed architecture               */
    domnmax = archDomSize (archptr, domnptr);     /* Get architecture size                    */

#ifdef SCOTCH_DEBUG_MAP2
  if (domnmax <= 0) {
    errorPrint ("mapInit2: target architecture must have at least one domain");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_MAP2 */

  mappptr->baseval = baseval;
  mappptr->vertnbr = vertnbr;
  mappptr->domnmax = domnmax + 1;                 /* +1 for empty domain in mapLoad */
  mappptr->domnnbr = 1;                           /* One domain in mapping to date  */
  mappptr->archdat = *archptr;
  mappptr->domnorg = *domnptr;

  if ((parttab = (Anum *) memAlloc (vertnbr * sizeof (Anum))) == NULL) { /* Allocate part array first as it will never move */
    errorPrint ("mapInit: out of memory (1)");
    return     (1);
  }
  mappptr->parttax = parttab - baseval;
  memSet (parttab, 0, vertnbr * sizeof (Anum));   /* All vertices mapped to first domain */

  if ((mappptr->domntab = (ArchDom *) memAlloc ((domnmax + 1) * sizeof (ArchDom))) == NULL) { /* Allocate possibly variable-sized domain array */
    errorPrint ("mapInit: out of memory (2)");
    return (1);
  }
  mappptr->domntab[0] = *domnptr;                 /* Set first domain */

  return (0);
}

int
mapInit (
Mapping * restrict const      mappptr,
const Gnum                    baseval,
const Gnum                    vertnbr,
const Arch * restrict const   archptr)
{
  ArchDom             domnorg;                    /* First, largest, domain */

  archDomFrst (archptr, &domnorg);                /* Get first domain of target architecture */
  return (mapInit2 (mappptr, baseval, vertnbr, archptr, &domnorg));
}

/* This routine frees the contents
** of the given mapping.
** It returns:
** - VOID  : in all cases.
*/

void
mapExit (
Mapping * const             mappptr)
{
  if (mappptr->domntab != NULL)
    memFree (mappptr->domntab);
  if (mappptr->parttax != NULL)
    memFree (mappptr->parttax + mappptr->baseval);

#ifdef SCOTCH_DEBUG_MAP2
  memSet (mappptr, ~0, sizeof (Mapping));
#endif /* SCOTCH_DEBUG_MAP2 */
}
