/* Copyright 2004,2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : arch_torus.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the torus graph     **/
/**                target architectures.                   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to   : 24 mar 1993     **/
/**                # Version 1.2  : from : 04 feb 1994     **/
/**                                 to   : 11 feb 1994     **/
/**                # Version 1.3  : from : 20 apr 1994     **/
/**                                 to   : 20 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to   : 23 dec 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to   : 29 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     08 sep 1995     **/
/**                # Version 3.1  : from : 07 may 1996     **/
/**                                 to     22 jul 1996     **/
/**                # Version 3.2  : from : 16 oct 1996     **/
/**                                 to     14 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 05 nov 2003     **/
/**                                 to     10 mar 2005     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     11 aug 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_TORUS

#include "module.h"
#include "common.h"
#include "arch.h"
#include "arch_torus.h"

/***********************************************/
/*                                             */
/* These are the 2-dimensional torus routines. */
/*                                             */
/***********************************************/

/* This routine loads the
** bidimensional torus architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archTorus2ArchLoad (
ArchTorus2 * restrict const archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchTorus2)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTorus2Dom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTorus2ArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((intLoad (stream, &archptr->c[0]) != 1) ||
      (intLoad (stream, &archptr->c[1]) != 1) ||
      (archptr->c[0] < 1) || (archptr->c[1] < 1)) {
    errorPrint ("archTorus2ArchLoad: bad input");
    return     (1);
  }

  return (0);
}

/* This routine saves the
** bidimensional torus architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archTorus2ArchSave (
const ArchTorus2 * const    archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchTorus2)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTorus2Dom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTorus2ArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " ",
               (Anum) archptr->c[0],
               (Anum) archptr->c[1]) == EOF) {
    errorPrint ("archTorus2ArchSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archTorus2DomNum (
const ArchTorus2 * const    archptr,
const ArchTorus2Dom * const domptr)
{
  return ((domptr->c[1][0] * archptr->c[0]) + domptr->c[0][0]); /* Return vertex number */
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archTorus2DomTerm (
const ArchTorus2 * const    archptr,
ArchTorus2Dom * const       domptr,
const ArchDomNum            domnum)
{
  if (domnum < (archptr->c[0] * archptr->c[1])) { /* If valid label */
    domptr->c[0][0] =                             /* Set the domain */
    domptr->c[0][1] = domnum % archptr->c[0];
    domptr->c[1][0] =
    domptr->c[1][1] = domnum / archptr->c[0];

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the rectangular domain.
*/

Anum 
archTorus2DomSize (
const ArchTorus2 * const    archptr,
const ArchTorus2Dom * const domptr)
{
  return ((domptr->c[0][1] - domptr->c[0][0] + 1) *
          (domptr->c[1][1] - domptr->c[1][0] + 1));
}

/* This function returns the average
** distance between two rectangular
** domains (in fact the distance between
** the centers of the domains).
*/

Anum 
archTorus2DomDist (
const ArchTorus2 * const    archptr,
const ArchTorus2Dom * const dom0ptr,
const ArchTorus2Dom * const dom1ptr)
{
  Anum               dc0, dc1;
  Anum               ds0, ds1;

  dc0 = abs (dom0ptr->c[0][0] + dom0ptr->c[0][1] -
             dom1ptr->c[0][0] - dom1ptr->c[0][1]);
  ds0 = (dc0 > archptr->c[0]) ? (2 * archptr->c[0] - dc0) : dc0;

  dc1 = abs (dom0ptr->c[1][0] + dom0ptr->c[1][1] -
             dom1ptr->c[1][0] - dom1ptr->c[1][1]);
  ds1 = (dc1 > archptr->c[1]) ? (2 * archptr->c[1] - dc1) : dc1;

  return ((ds0 + ds1) >> 1);
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archTorus2DomFrst (
const ArchTorus2 * const        archptr,
ArchTorus2Dom * restrict const  domptr)
{
  domptr->c[0][0] =
  domptr->c[1][0] = 0;
  domptr->c[0][1] = archptr->c[0] - 1;
  domptr->c[1][1] = archptr->c[1] - 1;

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archTorus2DomLoad (
const ArchTorus2 * const        archptr,
ArchTorus2Dom * restrict const  domptr,
FILE * restrict const           stream)
{
  if ((intLoad (stream, &domptr->c[0][0]) != 1) ||
      (intLoad (stream, &domptr->c[1][0]) != 1) ||
      (intLoad (stream, &domptr->c[0][1]) != 1) ||
      (intLoad (stream, &domptr->c[1][1]) != 1)) {
    errorPrint ("archTorus2DomLoad: bad input");
    return     (1);
  }

  return (0);
}

/* This routine saves domain information
** to the given stream.
** - 0   : on success.
** - !0  : on error.
*/

int
archTorus2DomSave (
const ArchTorus2 * const    archptr,
const ArchTorus2Dom * const domptr,
FILE * restrict const       stream)
{
  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " " ANUMSTRING " " ANUMSTRING " ",
               (Anum) domptr->c[0][0], (Anum) domptr->c[1][0],
               (Anum) domptr->c[0][1], (Anum) domptr->c[1][1]) == EOF) {
    errorPrint ("archTorus2DomSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function tries to split a rectangular
** domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archTorus2DomBipart (
const ArchTorus2 * const        archptr,
const ArchTorus2Dom * const     domptr,
ArchTorus2Dom * restrict const  dom0ptr,
ArchTorus2Dom * restrict const  dom1ptr)
{
  Anum                dimsiz[2];
  int                 dimval;                     /* Dimension along which to split */

  dimsiz[0] = domptr->c[0][1] - domptr->c[0][0];
  dimsiz[1] = domptr->c[1][1] - domptr->c[1][0];

  if ((dimsiz[0] | dimsiz[1]) == 0)               /* Return if cannot bipartition more */
    return (1);

  dimval = 1;
  if ((dimsiz[0] > dimsiz[1]) ||                  /* Split domain in two along largest dimension */
      ((dimsiz[0] == dimsiz[1]) && (archptr->c[0] > archptr->c[1])))
    dimval = 0;

  if (dimval == 0) {                              /* Split across the X dimension */
    dom0ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = (domptr->c[0][0] + domptr->c[0][1]) / 2;
    dom1ptr->c[0][0] = dom0ptr->c[0][1] + 1;
    dom1ptr->c[0][1] = domptr->c[0][1];
    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];
  }
  else {                                          /* Split across the Y dimension */
    dom0ptr->c[0][0] = dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];
    dom0ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = (domptr->c[1][0] + domptr->c[1][1]) / 2;
    dom1ptr->c[1][0] = dom0ptr->c[1][1] + 1;
    dom1ptr->c[1][1] = domptr->c[1][1];
  }

  return (0);
}

/* This function creates the MPI_Datatype for
** 2D torus domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archTorus2DomMpiType (
const ArchTorus2 * const      archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (4, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */

/***********************************************/
/*                                             */
/* These are the 3-dimensional torus routines. */
/*                                             */
/***********************************************/

/* This routine loads the
** tridimensional torus architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archTorus3ArchLoad (
ArchTorus3 * restrict const archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchTorus3)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTorus3Dom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTorus3ArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((intLoad (stream, &archptr->c[0]) != 1) ||
      (intLoad (stream, &archptr->c[1]) != 1) ||
      (intLoad (stream, &archptr->c[2]) != 1) ||
      (archptr->c[0] < 1) || (archptr->c[1] < 1) || (archptr->c[2] < 1)) {
    errorPrint ("archTorus3ArchLoad: bad input");
    return     (1);
  }

  return (0);
}

/* This routine saves the
** tridimensional torus architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archTorus3ArchSave (
const ArchTorus3 * const    archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchTorus3)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTorus3Dom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTorus3ArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " " ANUMSTRING " ",
               (Anum) archptr->c[0], (Anum) archptr->c[1], (Anum) archptr->c[2]) == EOF) {
    errorPrint ("archTorus3ArchSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archTorus3DomNum (
const ArchTorus3 * const    archptr,
const ArchTorus3Dom * const domptr)
{
  return ((((domptr->c[2][0]  * archptr->c[1]) +  /* Return the vertex number */
             domptr->c[1][0]) * archptr->c[0]) +
             domptr->c[0][0]);
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archTorus3DomTerm (
const ArchTorus3 * const    archptr,
ArchTorus3Dom * const       domptr,
const ArchDomNum            domnum)
{
  if (domnum < (archptr->c[0] * archptr->c[1] * archptr->c[2])) { /* If valid label */
    domptr->c[0][0] =                             /* Set the domain                 */
    domptr->c[0][1] = domnum % archptr->c[0];
    domptr->c[1][0] =
    domptr->c[1][1] = (domnum / archptr->c[0]) % archptr->c[1];
    domptr->c[2][0] =
    domptr->c[2][1] = domnum / (archptr->c[0] * archptr->c[1]);

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the cubic domain.
*/

Anum 
archTorus3DomSize (
const ArchTorus3 * const    archptr,
const ArchTorus3Dom * const domptr)
{
  return ((domptr->c[0][1] - domptr->c[0][0] + 1) *
          (domptr->c[1][1] - domptr->c[1][0] + 1) *
          (domptr->c[2][1] - domptr->c[2][0] + 1));
}

/* This function returns the average distance
** between two cubic domains (in fact the
** distance between the centers of the domains).
*/

Anum 
archTorus3DomDist (
const ArchTorus3 * const    archptr,
const ArchTorus3Dom * const dom0ptr,
const ArchTorus3Dom * const dom1ptr)
{
  Anum               dc0, dc1, dc2;
  Anum               ds0, ds1, ds2;

  dc0 = abs (dom0ptr->c[0][0] + dom0ptr->c[0][1] -
             dom1ptr->c[0][0] - dom1ptr->c[0][1]);
  ds0 = (dc0 > archptr->c[0]) ? (2 * archptr->c[0] - dc0) : dc0;

  dc1 = abs (dom0ptr->c[1][0] + dom0ptr->c[1][1] -
             dom1ptr->c[1][0] - dom1ptr->c[1][1]);
  ds1 = (dc1 > archptr->c[1]) ? (2 * archptr->c[1] - dc1) : dc1;

  dc2 = abs (dom0ptr->c[2][0] + dom0ptr->c[2][1] -
             dom1ptr->c[2][0] - dom1ptr->c[2][1]);
  ds2 = (dc2 > archptr->c[2]) ? (2 * archptr->c[2] - dc2) : dc2;

  return ((ds0 + ds1 + ds2) >> 1);
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archTorus3DomFrst (
const ArchTorus3 * const        archptr,
ArchTorus3Dom * restrict const  domptr)
{
  domptr->c[0][0] =
  domptr->c[1][0] =
  domptr->c[2][0] = 0;
  domptr->c[0][1] = archptr->c[0] - 1;
  domptr->c[1][1] = archptr->c[1] - 1;
  domptr->c[2][1] = archptr->c[2] - 1;

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archTorus3DomLoad (
const ArchTorus3 * const        archptr,
ArchTorus3Dom * restrict const  domptr,
FILE * restrict const           stream)
{
  if ((intLoad (stream, &domptr->c[0][0]) != 1) ||
      (intLoad (stream, &domptr->c[1][0]) != 1) ||
      (intLoad (stream, &domptr->c[2][0]) != 1) ||
      (intLoad (stream, &domptr->c[0][1]) != 1) ||
      (intLoad (stream, &domptr->c[1][1]) != 1) ||
      (intLoad (stream, &domptr->c[2][1]) != 1)) {
    errorPrint ("archTorus3DomLoad: bad input");
    return     (1);
  }

  return (0);
}

/* This routine saves domain information
** to the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archTorus3DomSave (
const ArchTorus3 * const    archptr,
const ArchTorus3Dom * const domptr,
FILE * restrict const       stream)
{
  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " " ANUMSTRING " " ANUMSTRING " " ANUMSTRING " " ANUMSTRING " ",
               (Anum) domptr->c[0][0], (Anum) domptr->c[1][0], (Anum) domptr->c[2][0],
               (Anum) domptr->c[0][1], (Anum) domptr->c[1][1], (Anum) domptr->c[2][1]) == EOF) {
    errorPrint ("archTorus3DomSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function tries to split a cubic
** domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archTorus3DomBipart (
const ArchTorus3 * const        archptr,
const ArchTorus3Dom * const     domptr,
ArchTorus3Dom * restrict const  dom0ptr,
ArchTorus3Dom * restrict const  dom1ptr)
{
  Anum                dimsiz[3];
  int                 dimtmp;
  int                 dimval;

  dimsiz[0] = domptr->c[0][1] - domptr->c[0][0];
  dimsiz[1] = domptr->c[1][1] - domptr->c[1][0];
  dimsiz[2] = domptr->c[2][1] - domptr->c[2][0];

  if ((dimsiz[0] | dimsiz[1] | dimsiz[2]) == 0)   /* Return if cannot bipartition more */
    return (1);

  dimval = (archptr->c[1] > archptr->c[0]) ? 1 : 0; /* Assume all subdomain dimensions are equal */
  if (archptr->c[2] > archptr->c[dimval])         /* Find priviledged dimension                  */
    dimval = 2;

  dimtmp = dimval;                                /* Find best dimension */
  if (dimsiz[(dimtmp + 1) % 3] > dimsiz[dimval])
    dimval = (dimtmp + 1) % 3;
  if (dimsiz[(dimtmp + 2) % 3] > dimsiz[dimval])
    dimval = (dimtmp + 2) % 3;

  if (dimval == 0) {                              /* Split domain in two along largest dimension */
    dom0ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = (domptr->c[0][0] + domptr->c[0][1]) / 2;
    dom1ptr->c[0][0] = dom0ptr->c[0][1] + 1;
    dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];

    dom0ptr->c[2][0] = dom1ptr->c[2][0] = domptr->c[2][0];
    dom0ptr->c[2][1] = dom1ptr->c[2][1] = domptr->c[2][1];
  }
  else if (dimval == 1) {
    dom0ptr->c[0][0] = dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = (domptr->c[1][0] + domptr->c[1][1]) / 2;
    dom1ptr->c[1][0] = dom0ptr->c[1][1] + 1;
    dom1ptr->c[1][1] = domptr->c[1][1];

    dom0ptr->c[2][0] = dom1ptr->c[2][0] = domptr->c[2][0];
    dom0ptr->c[2][1] = dom1ptr->c[2][1] = domptr->c[2][1];
  }
  else {
    dom0ptr->c[0][0] = dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];

    dom0ptr->c[2][0] = domptr->c[2][0];
    dom0ptr->c[2][1] = (domptr->c[2][0] + domptr->c[2][1]) / 2;
    dom1ptr->c[2][0] = dom0ptr->c[2][1] + 1;
    dom1ptr->c[2][1] = domptr->c[2][1];
  }

  return (0);
}

/* This function creates the MPI_Datatype for
** 3D torus domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archTorus3DomMpiType (
const ArchTorus3 * const      archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (6, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */
