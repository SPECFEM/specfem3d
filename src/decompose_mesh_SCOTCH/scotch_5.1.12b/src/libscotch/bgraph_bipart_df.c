/* Copyright 2004,2007,2008,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : bgraph_bipart_df.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a bipartition of   **/
/**                a bipartition graph by using a          **/
/**                diffusion scheme.                       **/
/**                                                        **/
/**   NOTES      : # This algorithm has been designed to   **/
/**                  work on band graphs only, for which   **/
/**                  the two anchor vertices are the two   **/
/**                  last vertices, the before-last as     **/
/**                  anchor of part 0, and the last as     **/
/**                  anchor of part 1.                     **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 09 jan 2007     **/
/**                                 to     10 sep 2007     **/
/**                # Version 5.1  : from : 29 oct 2007     **/
/**                                 to     27 mar 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BGRAPH_BIPART_DF

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "bgraph.h"
#include "bgraph_bipart_df.h"

/*
**  The static variables.
*/

static const Gnum           bgraphbipartdfloadzero = 0;
static const Gnum           bgraphbipartdfloadone  = 1;

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the bipartitioning.
** It returns:
** - 0 : if bipartitioning could be computed.
** - 1 : on error.
*/

int
bgraphBipartDf (
Bgraph * restrict const           grafptr,        /*+ Active graph      +*/
const BgraphBipartDfParam * const paraptr)        /*+ Method parameters +*/
{
  float * restrict      edlstax;                  /* Edge load sum array       */
  float * restrict      veextax;                  /* External leaks array      */
  float * restrict      difotax;                  /* Old diffusion value array */
  float * restrict      difntax;                  /* New diffusion value array */
  float                 cdifval;
  float                 cremval;
  Gnum                  fronnum;
  Gnum                  compload0;
  Gnum                  compload1;
  Gnum                  compsize1;
  Gnum                  commloadintn;
  Gnum                  commloadextn;
  Gnum                  commgainextn;
  Gnum                  veexnbr;
  float                 veexval;
  const Gnum * restrict veexbax;
  Gnum                  veexmsk;
  float                 veloval;
  const Gnum * restrict velobax;
  Gnum                  velomsk;
  Gnum                  vancval0;                 /* Initial values for both anchors */
  Gnum                  vancval1;
  float                 vanctab[2];               /* Value to add to each anchor */
  Gnum                  vertnum;
  INT                   passnum;

  veexnbr = (grafptr->veextax != NULL) ? grafptr->s.vertnbr : 0;
  if (memAllocGroup ((void **) (void *)
                     &edlstax, (size_t) (grafptr->s.vertnbr * sizeof (float)),
                     &veextax, (size_t) (veexnbr            * sizeof (float)),
                     &difotax, (size_t) (grafptr->s.vertnbr * sizeof (float)),
                     &difntax, (size_t) (grafptr->s.vertnbr * sizeof (float)), NULL) == NULL) {
    errorPrint ("bgraphBipartDf: out of memory (1)");
    return     (1);
  }
  edlstax -= grafptr->s.baseval;                  /* Base access to veextax and diffusion arrays */
  difotax -= grafptr->s.baseval;
  difntax -= grafptr->s.baseval;
  veextax  = (grafptr->veextax != NULL) ? veextax - grafptr->s.baseval : NULL;

  compload0 = (paraptr->typeval == BGRAPHBIPARTDFTYPEBAL) /* If balanced parts wanted */
              ? grafptr->compload0avg             /* Target is average                */
              : ( (grafptr->compload0 < grafptr->compload0min) ? grafptr->compload0min : /* Else keep load if not off balance */
                 ((grafptr->compload0 > grafptr->compload0max) ? grafptr->compload0max : grafptr->compload0));
  vancval0 = (float) - compload0;                 /* Values to be injected to anchor vertices at every iteration */
  vancval1 = (float) (grafptr->s.velosum - compload0);

  if (grafptr->s.edlotax == NULL) {               /* If graph doesn't have edge weights */
    Gnum                vertnum;

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
      difotax[vertnum] = 0.0;
      edlstax[vertnum] = (float) (grafptr->s.vendtax[vertnum] - grafptr->s.verttax[vertnum]);
    }
  }
  else {                                          /* If graph has edge weights */
    Gnum                vertnum;

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
      Gnum                edgenum;
      Gnum                edgennd;
      Gnum                edlosum;

      for (edgenum = grafptr->s.verttax[vertnum], edgennd = grafptr->s.vendtax[vertnum], edlosum = 0;
           edgenum < edgennd; edgenum ++)
        edlosum += grafptr->s.edlotax[edgenum];

      difotax[vertnum] = 0.0;
      edlstax[vertnum] = (float) edlosum;
    }
  }

  if (veextax != NULL) {
    Gnum                vertnum;
    Gnum                veexsum;
    Gnum                veexsum0;
    float               idodist;                  /* Inverse of domdist */

    idodist = 1.0 / (float) grafptr->domdist;
    for (vertnum = grafptr->s.baseval, veexsum = veexsum0 = 0;
         vertnum < grafptr->s.vertnnd; vertnum ++) {
      Gnum                veexval;

      veexval = grafptr->veextax[vertnum];
      veexsum  += veexval;                        /* Sum all external gains, positive and negative           */
      veexsum0 += BGRAPHBIPARTDFGNUMSGNMSK (veexval); /* Sum all negative external gains; superscalar update */
      veextax[vertnum] = (float) veexval / idodist;
    }
    vancval0 += (float) veexsum0;
    vancval1 += (float) (veexsum - veexsum0);
  }
  vancval1 -= BGRAPHBIPARTDFEPSILON;              /* Slightly tilt value to add to part 1 */

  difotax[grafptr->s.vertnnd - 2] = vancval0 / edlstax[grafptr->s.vertnnd - 2]; /* Load anchor vertices for first pass */
  difotax[grafptr->s.vertnnd - 1] = vancval1 / edlstax[grafptr->s.vertnnd - 1];

  veexval = 0.0F;                                 /* Assume no external gains */
  veloval = 1.0F;                                 /* Assume no vertex loads   */
  cdifval = paraptr->cdifval;
  cremval = paraptr->cremval;
  for (passnum = 0; passnum < paraptr->passnbr; passnum ++) { /* For all passes */
    Gnum                vertnum;
    Gnum                vertnnd;
    float               vancval;                  /* Value to load vertex with if anchor */
    float *             difttax;                  /* Temporary swap value                */

    vanctab[0] = vancval0;                        /* Copy anchor values to injection variables */
    vanctab[1] = vancval1;
    vancval    = 0.0F;                            /* At first vertices are not anchors */
    vertnum    = grafptr->s.baseval;
    vertnnd    = grafptr->s.vertnnd - 2;
    while (1) {
      for ( ; vertnum < vertnnd; vertnum ++) {
        Gnum                edgenum;
        Gnum                edgennd;
        float               diffval;

        edgenum = grafptr->s.verttax[vertnum];
        edgennd = grafptr->s.vendtax[vertnum];
        diffval = 0.0F;
        if (grafptr->s.edlotax != NULL)
          for ( ; edgenum < edgennd; edgenum ++)
            diffval += difotax[grafptr->s.edgetax[edgenum]] * (float) grafptr->s.edlotax[edgenum];
        else
          for ( ; edgenum < edgennd; edgenum ++)
            diffval += difotax[grafptr->s.edgetax[edgenum]];

        diffval *= cdifval;
        diffval += difotax[vertnum] * cremval * edlstax[vertnum] + vancval;
        if (grafptr->s.velotax != NULL)
          veloval = (float) grafptr->s.velotax[vertnum];
        if (veextax != NULL) {
          veexval = veextax[vertnum] * diffval;
          if (veexval >= 0.0F)                    /* If vertex is already in right part */
            veexval = 0.0F;                       /* Then it will not be impacted       */
          else {
            int                 partval;

            partval = (diffval < 0.0F) ? 0 : 1;   /* Load anchor with load removed from vertex */
            vanctab[partval] += (float) (2 * partval - 1) * veexval;
          }
        }
        if (diffval >= 0.0F) {
          diffval -= veloval + veexval;
          if (diffval <= 0.0F)
            diffval = +BGRAPHBIPARTDFEPSILON;
        }
        else {
          diffval += veloval + veexval;
          if (diffval >= 0.0F)
            diffval = -BGRAPHBIPARTDFEPSILON;
        }
        if (isnan (diffval))                      /* If overflow occured (because of avalanche process)                        */
          goto abort;                             /* Exit main loop without swapping arrays so as to keep last valid iteration */

        difntax[vertnum] = diffval / edlstax[vertnum];
      }
      if (vertnum == grafptr->s.vertnnd)          /* If all vertices processed, exit intermediate infinite loop */
        break;

      vertnnd ++;                                 /* Prepare to go only for one more run      */
      vancval = vanctab[vertnum - grafptr->s.vertnnd + 2]; /* Load variable with anchor value */
    }

    difttax = (float *) difntax;                  /* Swap old and new diffusion arrays          */
    difntax = (float *) difotax;                  /* Casts to prevent IBM compiler from yelling */
    difotax = (float *) difttax;
  }
abort :                                           /* If overflow occured, resume here */

  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) /* Update part according to diffusion state */
    grafptr->parttax[vertnum] = (difotax[vertnum] <= 0.0F) ? 0 : 1;

  if (grafptr->veextax != NULL) {
    veexbax = grafptr->veextax;
    veexmsk = ~((Gnum) 0);
  }
  else {
    veexbax = &bgraphbipartdfloadzero;
    veexmsk = 0;
  }
  if (grafptr->s.velotax != NULL) {
    velobax = grafptr->s.velotax;
    velomsk = ~((Gnum) 0);
  }
  else {
    velobax = &bgraphbipartdfloadone;
    velomsk = 0;
  }

  for (vertnum = grafptr->s.baseval, fronnum = 0, commloadextn = grafptr->commloadextn0, commgainextn = commloadintn = compload1 = compsize1 = 0;
       vertnum < grafptr->s.vertnnd; vertnum ++) {
    Gnum                edgenum;
    Gnum                partval;
    Gnum                commload;                 /* Vertex internal communication load */

    partval = (Gnum) grafptr->parttax[vertnum];
    compsize1    += partval;
    compload1    += partval * velobax[vertnum & velomsk];
    commloadextn += partval * veexbax[vertnum & veexmsk];
    commgainextn += (1 - 2 * partval) * veexbax[vertnum & veexmsk];
    commload      = 0;
    if (grafptr->s.edlotax != NULL) {
      for (edgenum = grafptr->s.verttax[vertnum]; edgenum < grafptr->s.vendtax[vertnum]; edgenum ++) {
        Gnum                partend;

        partend   = (Gnum) grafptr->parttax[grafptr->s.edgetax[edgenum]];
        commload += (partval ^ partend) * grafptr->s.edlotax[edgenum];
      }
    }
    else {
      for (edgenum = grafptr->s.verttax[vertnum]; edgenum < grafptr->s.vendtax[vertnum]; edgenum ++)
        commload += partval ^ (Gnum) grafptr->parttax[grafptr->s.edgetax[edgenum]];
    }
    commloadintn += commload;                     /* Internal loads will be added twice */
    if (commload != 0)                            /* If end vertex is in the other part */
      grafptr->frontab[fronnum ++] = vertnum;     /* Then it belongs to the frontier    */
  }
  grafptr->fronnbr      = fronnum;
  grafptr->compload0    = grafptr->s.velosum - compload1;
  grafptr->compload0dlt = grafptr->compload0 - grafptr->compload0avg;
  grafptr->compsize0    = grafptr->s.vertnbr - compsize1;
  grafptr->commload     = commloadextn + (commloadintn / 2) * grafptr->domdist;
  grafptr->commgainextn = commgainextn;
  grafptr->bbalval      = (double) ((grafptr->compload0dlt < 0) ? (- grafptr->compload0dlt) : grafptr->compload0dlt) / (double) grafptr->compload0avg;

  memFree (edlstax + grafptr->s.baseval);

#ifdef SCOTCH_DEBUG_BGRAPH2
  if (bgraphCheck (grafptr) != 0) {
    errorPrint ("bgraphBipartDf: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

  return (0);
}
