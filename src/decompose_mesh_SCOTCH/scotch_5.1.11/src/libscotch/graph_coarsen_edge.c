/* Copyright 2004,2007,2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : graph_coarsen_edge.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This commodity file contains the edge   **/
/**                arrays building subroutine which is     **/
/**                duplicated, with minor modifications,   **/
/**                into graph_coarsen.c.                   **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 17 dec 2001     **/
/**                                 to     25 feb 2004     **/
/**                # Version 5.0  : from : 13 dec 2006     **/
/**                                 to     14 dec 2006     **/
/**                # Version 5.1  : from : 30 oct 2009     **/
/**                                 to     30 oct 2009     **/
/**                                                        **/
/************************************************************/

static
void
GRAPHCOARSENEDGENAME (
const Graph * restrict const              finegrafptr, /*+ Graph to coarsen      +*/
const Gnum * restrict const               finecoartax, /*+ Fine to coarse array  +*/
const GraphCoarsenMulti * restrict const  coarmulttax, /*+ Multinode array       +*/
Graph * restrict const                    coargrafptr, /*+ Coarse graph to build +*/
GraphCoarsenHash * restrict const         coarhashtab, /*+ End vertex hash table +*/
const Gnum                                coarhashmsk) /*+ Hash table mask       +*/
{
  Gnum                coarvertnum;
  Gnum                coaredgenum;
  Gnum                coardegrmax;
  Gnum                coaredlosum;

  Gnum * restrict const       coarverttax = coargrafptr->verttax;
  Gnum * restrict const       coaredgetax = coargrafptr->edgetax;
  Gnum * restrict const       coaredlotax = coargrafptr->edlotax;
  const Gnum * restrict const fineverttax = finegrafptr->verttax;
  const Gnum * restrict const finevendtax = finegrafptr->vendtax;
  const Gnum * restrict const fineedgetax = finegrafptr->edgetax;
  GRAPHCOARSENEDGEINIT;

  coaredlosum = 0;
  for (coarvertnum = coaredgenum = coargrafptr->baseval, coardegrmax = 0;
       coarvertnum < coargrafptr->vertnnd; coarvertnum ++) {
    Gnum                finevertnum;
    int                 i;

    coarverttax[coarvertnum] = coaredgenum;       /* Set vertex edge index */
    i = 0;
    do {                                          /* For all fine edges of multinode vertices */
      Gnum                fineedgenum;

      finevertnum = coarmulttax[coarvertnum].vertnum[i];
      for (fineedgenum = fineverttax[finevertnum];
           fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
        Gnum                coarvertend;          /* Number of coarse vertex which is end of fine edge */
        Gnum                h;

        coarvertend = finecoartax[fineedgetax[fineedgenum]];
        if (coarvertend != coarvertnum) {         /* If not end of collapsed edge */
          for (h = (coarvertend * GRAPHCOARHASHPRIME) & coarhashmsk; ; h = (h + 1) & coarhashmsk) {
            if (coarhashtab[h].vertorgnum != coarvertnum) { /* If old slot           */
              coarhashtab[h].vertorgnum = coarvertnum; /* Mark it in reference array */
              coarhashtab[h].vertendnum = coarvertend;
              coarhashtab[h].edgenum    = coaredgenum;
              coaredgetax[coaredgenum]  = coarvertend; /* One more edge created */
              GRAPHCOARSENEDGEEDLOINIT;           /* Initialize edge load entry */
              coaredgenum ++;
              break;                              /* Give up hashing */
            }
            if (coarhashtab[h].vertendnum == coarvertend) { /* If coarse edge already exists */
              GRAPHCOARSENEDGEEDLOADD;            /* Accumulate edge load                    */
              break;                              /* Give up hashing                         */
            }
          }
        }
        else
          GRAPHCOARSENEDGEEDLOSUB;
      }
    } while (i ++, finevertnum != coarmulttax[coarvertnum].vertnum[1]); /* Skip to next matched vertex if both vertices not equal */

    if (coardegrmax < (coaredgenum - coarverttax[coarvertnum]))
      coardegrmax = coaredgenum - coarverttax[coarvertnum];
  }
  coarverttax[coarvertnum] = coaredgenum;         /* Mark end of edge array */

  coargrafptr->edlosum = finegrafptr->edlosum + coaredlosum; /* Subtract all matched edges */
  coargrafptr->degrmax = coardegrmax;
}
