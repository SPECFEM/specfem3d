/* Copyright 2008-2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : parmetis_dgraph_part.c                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the compatibility        **/
/**                library for the ParMeTiS ordering       **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 19 jun 2008     **/
/**                                 to     30 jun 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "ptscotch.h"
#include "parmetis.h"                             /* Our "parmetis.h" file */

/************************************/
/*                                  */
/* These routines are the C API for */
/* the ParMeTiS graph ordering      */
/* routine.                         */
/*                                  */
/************************************/

void
METISNAMEU(ParMETIS_V3_PartKway) (
const int * const           vtxdist,
int * const                 xadj,
int * const                 adjncy,
int * const                 vwgt,
int * const                 adjwgt,
const int * const           wgtflag,
const int * const           numflag,
const int * const           ncon,                 /* Not used */
const int * const           nparts,
const float * const         tpwgts,
const float * const         ubvec,                /* Not used */
const int * const           options,              /* Not used */
int * const                 edgecut,
int * const                 part,
MPI_Comm *                  comm)
{
  MPI_Comm            proccomm;
  int                 procglbnbr;
  int                 proclocnum;
  SCOTCH_Num          baseval;
  SCOTCH_Arch         archdat;
  SCOTCH_Dgraph       grafdat;                    /* Scotch distributed graph object to interface with libScotch   */
  SCOTCH_Dmapping     mappdat;                    /* Scotch distributed mapping object to interface with libScotch */
  SCOTCH_Strat        stradat;
  SCOTCH_Num          vertlocnbr;
  SCOTCH_Num *        veloloctab;
  SCOTCH_Num          edgelocnbr;
  SCOTCH_Num *        edloloctab;
  SCOTCH_Num *        velotab;
  double *            vwgttab;
  int                 i;

  if (sizeof (SCOTCH_Num) != sizeof (int)) {
    SCOTCH_errorPrint ("ParMETIS_V3_PartKway (as of SCOTCH): SCOTCH_Num type should equate to int");
    return;
  }

  if ((vwgttab = malloc (*nparts * sizeof (double))) == NULL)
    return;
  if ((velotab = malloc (*nparts * sizeof (SCOTCH_Num))) == NULL) {
    free (vwgttab);
    return;
  }
  for (i = 0; i < *nparts; i ++)
    vwgttab[i] = (double) tpwgts[i] * (double) (*nparts);
  for (i = 0; i < *nparts; i ++) {
    double deltval;
    deltval = fabs (vwgttab[i] - floor (vwgttab[i] + 0.5));
    if (deltval > 0.01) {
      int                 j;

      deltval = 1.0 / deltval;
      for (j = 0; j < *nparts; j ++)
        vwgttab[j] *= deltval;
    }
  }
  for (i = 0; i < *nparts; i ++)
    velotab[i] = (SCOTCH_Num) (vwgttab[i] + 0.5);

  proccomm = *comm;
  if (SCOTCH_dgraphInit (&grafdat, proccomm) != 0)
    return;

  MPI_Comm_size (proccomm, &procglbnbr);
  MPI_Comm_rank (proccomm, &proclocnum);
  baseval    = *numflag;
  vertlocnbr = vtxdist[proclocnum + 1] - vtxdist[proclocnum];
  edgelocnbr = xadj[vertlocnbr] - baseval;
  veloloctab = ((vwgt   != NULL) && ((*wgtflag & 2) != 0)) ? vwgt   : NULL;
  edloloctab = ((adjwgt != NULL) && ((*wgtflag & 1) != 0)) ? adjwgt : NULL;

  if (SCOTCH_dgraphBuild (&grafdat, baseval,
                          vertlocnbr, vertlocnbr, xadj, xadj + 1, veloloctab, NULL,
                          edgelocnbr, edgelocnbr, adjncy, NULL, edloloctab) == 0) {
    SCOTCH_stratInit (&stradat);
#ifdef SCOTCH_DEBUG_ALL
    if (SCOTCH_dgraphCheck (&grafdat) == 0)       /* TRICK: next instruction called only if graph is consistent */
#endif /* SCOTCH_DEBUG_ALL */
    {
      SCOTCH_archInit (&archdat);

      if ((SCOTCH_archCmpltw (&archdat, *nparts, velotab) == 0) &&
          (SCOTCH_dgraphMapInit (&grafdat, &mappdat, &archdat, part) == 0)) {
        SCOTCH_dgraphMapCompute (&grafdat, &mappdat, &stradat);

        SCOTCH_dgraphMapExit (&grafdat, &mappdat);
      }
      SCOTCH_archExit (&archdat);
    }
    SCOTCH_stratExit (&stradat);
  }
  SCOTCH_dgraphExit (&grafdat);

  *edgecut = 0;                                   /* TODO : compute real edge cut for people who might want it */

  free (vwgttab);
  free (velotab);

  if (baseval != 0) {                             /* MeTiS part array is based, Scotch is not */
    SCOTCH_Num          vertlocnum;

    for (vertlocnum = 0; vertlocnum < vertlocnbr; vertlocnum ++)
      part[vertlocnum] += baseval;
  }
}

/*
**
*/

void
METISNAMEU(ParMETIS_V3_PartGeomKway) (
const int * const           vtxdist,
int * const                 xadj,
int * const                 adjncy,
int * const                 vwgt,
int * const                 adjwgt,
const int * const           wgtflag,
const int * const           numflag,
const int * const           ndims,                /* Not used */
const float * const         xyz,                  /* Not used */
const int * const           ncon,                 /* Not used */
const int * const           nparts,
const float * const         tpwgts,
const float * const         ubvec,
const int * const           options,              /* Not used */
int * const                 edgecut,
int * const                 part,
MPI_Comm *                  commptr)
{
  METISNAMEU(ParMETIS_V3_PartKway) (vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ncon, nparts, tpwgts, ubvec, options, edgecut, part, commptr);
}
