/* Copyright 2004,2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : arch_torus.h                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the tori graph target architecture  **/
/**                functions.                              **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to   : 24 mar 1993     **/
/**                # Version 1.2  : from : 04 feb 1994     **/
/**                                 to   : 11 feb 1994     **/
/**                # Version 1.3  : from : 20 apr 1994     **/
/**                                 to   : 20 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to   : 12 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to   : 30 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     17 aug 1995     **/
/**                # Version 3.1  : from : 22 jul 1996     **/
/**                                 to     23 jul 1996     **/
/**                # Version 3.2  : from : 16 oct 1996     **/
/**                                 to     14 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 05 nov 2003     **/
/**                                 to     05 nov 2003     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     21 jan 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ The 2D-torus definitions. +*/

typedef struct ArchTorus2_ {
  Anum                      c[2];                 /*+ Mesh dimensions +*/
} ArchTorus2;

typedef struct ArchTorus2Dom_ {
  Anum                      c[2][2];              /*+ Inclusive X and Y coordinates +*/
} ArchTorus2Dom;

/*+ The 3D-torus definitions. +*/

typedef struct ArchTorus3_ {
  Anum                      c[3];                 /*+ Mesh dimensions +*/
} ArchTorus3;

typedef struct ArchTorus3Dom_ {
  Anum                      c[3][2];              /*+ Inclusive X, Y, and Z coordinates +*/
} ArchTorus3Dom;

/*
**  The function prototypes.
*/

#ifndef ARCH_TORUS
#define static
#endif

int                         archTorus2ArchLoad  (ArchTorus2 * restrict const, FILE * restrict const);
int                         archTorus2ArchSave  (const ArchTorus2 * const, FILE * restrict const);
#define archTorus2ArchFree          NULL
ArchDomNum                  archTorus2DomNum    (const ArchTorus2 * const, const ArchTorus2Dom * const);
int                         archTorus2DomTerm   (const ArchTorus2 * const, ArchTorus2Dom * restrict const, const ArchDomNum);
Anum                        archTorus2DomSize   (const ArchTorus2 * const, const ArchTorus2Dom * const);
#define archTorus2DomWght           archTorus2DomSize
Anum                        archTorus2DomDist   (const ArchTorus2 * const, const ArchTorus2Dom * const, const ArchTorus2Dom * const);
int                         archTorus2DomFrst   (const ArchTorus2 * const, ArchTorus2Dom * const);
int                         archTorus2DomLoad   (const ArchTorus2 * const, ArchTorus2Dom * const, FILE * restrict const);
int                         archTorus2DomSave   (const ArchTorus2 * const, const ArchTorus2Dom * const, FILE * restrict const);
int                         archTorus2DomBipart (const ArchTorus2 * const, const ArchTorus2Dom * const, ArchTorus2Dom * restrict const, ArchTorus2Dom * restrict const);
int                         archTorus2DomBipartO (const ArchTorus2 * const, const ArchTorus2Dom * const, ArchTorus2Dom * restrict const, ArchTorus2Dom * restrict const);
int                         archTorus2DomBipartU (const ArchTorus2 * const, const ArchTorus2Dom * const, ArchTorus2Dom * restrict const, ArchTorus2Dom * restrict const);
#ifdef SCOTCH_PTSCOTCH
int                         archTorus2DomMpiType (const ArchTorus2 * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

int                         archTorus3ArchLoad  (ArchTorus3 * restrict const, FILE * restrict const);
int                         archTorus3ArchSave  (const ArchTorus3 * const, FILE * restrict const);
#define archTorus3ArchFree          NULL
ArchDomNum                  archTorus3DomNum    (const ArchTorus3 * const, const ArchTorus3Dom * const);
int                         archTorus3DomTerm   (const ArchTorus3 * const, ArchTorus3Dom * restrict const, const ArchDomNum);
Anum                        archTorus3DomSize   (const ArchTorus3 * const, const ArchTorus3Dom * const);
#define archTorus3DomWght           archTorus3DomSize
Anum                        archTorus3DomDist   (const ArchTorus3 * const, const ArchTorus3Dom * const, const ArchTorus3Dom * const);
int                         archTorus3DomFrst   (const ArchTorus3 * const, ArchTorus3Dom * const);
int                         archTorus3DomLoad   (const ArchTorus3 * const, ArchTorus3Dom * const, FILE * restrict const);
int                         archTorus3DomSave   (const ArchTorus3 * const, const ArchTorus3Dom * const, FILE * restrict const);
int                         archTorus3DomBipart (const ArchTorus3 * const, const ArchTorus3Dom * const, ArchTorus3Dom * restrict const, ArchTorus3Dom * restrict const);
#ifdef SCOTCH_PTSCOTCH
int                         archTorus3DomMpiType (const ArchTorus3 * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

#undef static
