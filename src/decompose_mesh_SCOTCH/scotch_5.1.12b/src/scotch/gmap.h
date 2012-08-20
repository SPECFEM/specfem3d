/* Copyright 2004,2007,2008,2010,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : gmap.h                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a graph static mapper.          **/
/**                These lines are the data declaration    **/
/**                for the main routine.                   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 05 jan 1993     **/
/**                                 to     12 may 1993     **/
/**                # Version 1.3  : from : 09 apr 1994     **/
/**                                 to     30 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     08 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     09 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     15 aug 1995     **/
/**                # Version 3.1  : from : 07 nov 1995     **/
/**                                 to     10 nov 1995     **/
/**                # Version 3.2  : from : 04 oct 1996     **/
/**                                 to     18 jul 1997     **/
/**                # Version 3.3  : from : 07 oct 1998     **/
/**                                 to   : 31 may 1999     **/
/**                # Version 4.0  : from : 16 jan 2004     **/
/**                                 to   : 16 jan 2004     **/
/**                # Version 5.0  : from : 12 jun 2008     **/
/**                                 to   : 18 jun 2008     **/
/**                # Version 5.1  : from : 28 aug 2010     **/
/**                                 to   : 18 jul 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ File name aliases. +*/

#define C_FILENBR                   4             /* Number of files in list */

#define C_filenamesrcinp            C_fileTab[0].name /* Source graph input file name        */
#define C_filenametgtinp            C_fileTab[1].name /* Target architecture input file name */
#define C_filenamemapout            C_fileTab[2].name /* Mapping result output file name     */
#define C_filenamelogout            C_fileTab[3].name /* Log file name                       */

#define C_filepntrsrcinp            C_fileTab[0].pntr /* Source graph input file        */
#define C_filepntrtgtinp            C_fileTab[1].pntr /* Target architecture input file */
#define C_filepntrmapout            C_fileTab[2].pntr /* Mapping result output file     */
#define C_filepntrlogout            C_fileTab[3].pntr /* Log file                       */

/*+ Process flags. +*/

#define C_FLAGNONE                  0x0000        /* No flags            */
#define C_FLAGPART                  0x0001        /* Partitioning        */
#define C_FLAGVERBSTR               0x0002        /* Verbose flags       */
#define C_FLAGVERBTIM               0x0004
#define C_FLAGVERBMAP               0x0008
#define C_FLAGKBALVAL               0x0010        /* Imbalance tolerance */
#define C_FLAGCLUSTER               0x0020        /* Clustering          */
