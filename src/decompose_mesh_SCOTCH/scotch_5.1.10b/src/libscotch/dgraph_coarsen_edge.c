/* Copyright 2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_coarsen_edge.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This commodity file contains the edge   **/
/**                arrays building subroutine which is     **/
/**                duplicated, with minor modifications,   **/
/**                into dgraph_coarsen.c.                  **/
/**                                                        **/
/**   DATES      : # Version 5.2  : from : 11 dec 2008     **/
/**                                 to   : 11 dec 2008     **/
/**                                                        **/
/************************************************************/

/* This routine performs the coarsening of edges
** with respect to the data structures filled to
** date:
** - the coarmulttax array, which contains the
**   coarse index of each local or ghost vertex,
** - the edgercvtab array, which contains compacted
**   edge data for all of the remote vertices which
**   have been sent to us,
** - the vertloctax array, which contains, for
**   multinodes which have a remote vertex, the
**   index of the degree and coarse end vertex
**   data of the coarse remote vertex.
*/

void
DGRAPHCOARSENEDGENAME (
DgraphCoarsenData * restrict const  coarptr)
{


  for (coarvertlocnum = coaredgelocnum = coargrafptr->baseval; /* For all coarse vertices (that is, multinodes) */
       coarvertlocnum < coarvertlocnnd; coarvertlocnum ++) {

    if (multloctax[multlocnum].vertglbnum[1] < 0) { /* If second multinode vertex is remote    */
      edgedatidx = coarvertloctax[coarvertlocnum]; /* Get index of vertex in remote edge array */
      coarvertloctax[coarvertlocnum] = coaredgelocnum; /* Set beginning of coarse vertex array */
      multloctax[multlocnum].vertglbnum[1] = fineedgeloctax[-2 - multloctax[multlocnum].vertglbnum[1]]; /* Finalize multinode */

      for (fineedgegstnum = edgedataidx + 1, fineedgegstnnd = fineedgegstnum + edgercvtab[edgedataidx];
           fineedgegstnum < fineedgegstnnd; ) {
        Gnum                coarvertglbend;       /* Number of coarse vertex which is end of fine edge */
        Gnum                h;

        coarvertglbend = edgercvtab[fineedgegstnum ++];
        fineedlogstval = edgercvtab[fineedgegstnum ++];

        if (coarvertglbend != coarvertlocnum + coarvertlocadj) { /* If not end of collapsed edge */
          for (h = (coarvertglbend * DGRAPHCOARHASHPRIME) & coarhashmsk; ; h = (h + 1) & coarhashmsk) {
            if (coarhashtab[h].vertorgnum != coarverloctnum) { /* If old slot           */
              coarhashtab[h].vertorgnum = coarvertlocnum; /* Mark it in reference array */
              coarhashtab[h].vertendnum = coarvertglbend;
              coarhashtab[h].edgelocnum = coaredgelocnum;
              coaredgeloctax[coaredgelocnum] = coarvertglbend; /* One more edge created */
              DGRAPHCOARSENEDGEEDLOINIT;          /* Initialize edge load entry         */
              coaredgelocnum ++;
              break;                              /* Give up hashing */
            }
            if (coarhashtab[h].vertendnum == coarvertglbend) { /* If coarse edge already exists */
              DGRAPHCOARSENEDGEEDLOADD;           /* Accumulate edge load                       */
              break;                              /* Give up hashing                            */
            }
          }
        }
        else
          DGRAPHCOARSENEDGEEDLOSUB;
      }
      j = 0;                                      /* Will not explore vertglbnum[1] again */
    }
    else {
      coarvertloctax[coarvertlocnum] = coaredgelocnum;
      j = 1;
    }

    i = 0;
    do {                                          /* For all fine edges of multinode vertices */
      Gnum                fineedgenum;

      finevertnum = coarmulttax[coarvertnum].vertnum[i];
      for (fineedgenum = finegrafptr->verttax[finevertnum];
           fineedgenum < finegrafptr->vendtax[finevertnum]; fineedgenum ++) {
        Gnum                coarvertend;          /* Number of coarse vertex which is end of fine edge */
        Gnum                h;

        coarvertend = finecoartax[finegrafptr->edgetax[fineedgenum]];
        if (coarvertend != coarvertnum) {         /* If not end of collapsed edge */
          for (h = (coarvertend * GRAPHCOARHASHPRIME) & coarhashmsk; ; h = (h + 1) & coarhashmsk) {
            if (coarhashtab[h].vertorgnum != coarvertnum) { /* If old slot           */
              coarhashtab[h].vertorgnum = coarvertnum; /* Mark it in reference array */
              coarhashtab[h].vertendnum = coarvertend;
              coarhashtab[h].edgenum    = coaredgenum;
              coargrafptr->edgetax[coaredgenum] = coarvertend; /* One more edge created */
              GRAPHCOARSENEDGEEDLOINIT;           /* Initialize edge load entry         */
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

    if (coardegrmax < (coaredgenum - coargrafptr->verttax[coarvertnum]))
      coardegrmax = coaredgenum - coargrafptr->verttax[coarvertnum];
  }
}
