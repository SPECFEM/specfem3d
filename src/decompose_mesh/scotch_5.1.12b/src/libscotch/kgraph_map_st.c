/* Copyright 2004,2007,2009-2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : kgraph_map_st.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the strategy and   **/
/**                method tables for the multi-way static  **/
/**                mapping routines.                       **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 15 oct 1996     **/
/**                                 to     26 may 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to     17 may 1999     **/
/**                # Version 3.4  : from : 12 sep 2001     **/
/**                                 to     12 sep 2001     **/
/**                # Version 4.0  : from : 12 jan 2004     **/
/**                                 to     05 jan 2005     **/
/**                # Version 5.1  : from : 04 oct 2009     **/
/**                                 to     29 mar 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_ST

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "graph_coarsen.h"
#include "arch.h"
#include "mapping.h"
#include "bgraph.h"
#include "bgraph_bipart_st.h"
#include "kgraph.h"
#include "kgraph_map_ml.h"
#include "kgraph_map_rb.h"
#include "kgraph_map_st.h"

/*
**  The static and global variables.
*/

static union {
  KgraphMapMlParam          param;
  StratNodeMethodData       padding;
} kgraphmapstdefaultml = { { 100, 0.8L, GRAPHCOARHEM, &stratdummy, &stratdummy } };

static union {
  KgraphMapRbParam          param;
  StratNodeMethodData       padding;
} kgraphmapstdefaultrb = { { 1, 1, KGRAPHMAPRBPOLINGSIZE, &stratdummy, 0.05 } };

static StratMethodTab       kgraphmapstmethtab[] = { /* Mapping methods array */
                              { KGRAPHMAPSTMETHML, "m",  kgraphMapMl, &kgraphmapstdefaultml },
                              { KGRAPHMAPSTMETHRB, "r",  kgraphMapRb, &kgraphmapstdefaultrb },
                              { -1,                NULL, NULL,        NULL } };

static StratParamTab        kgraphmapstparatab[] = { /* Method parameter list */
                              { KGRAPHMAPSTMETHML,  STRATPARAMSTRAT,  "asc",
                                (byte *) &kgraphmapstdefaultml.param,
                                (byte *) &kgraphmapstdefaultml.param.stratasc,
                                (void *) &kgraphmapststratab },
                              { KGRAPHMAPSTMETHML,  STRATPARAMSTRAT,  "low",
                                (byte *) &kgraphmapstdefaultml.param,
                                (byte *) &kgraphmapstdefaultml.param.stratlow,
                                (void *) &kgraphmapststratab },
                              { KGRAPHMAPSTMETHML,  STRATPARAMCASE,   "type",
                                (byte *) &kgraphmapstdefaultml.param,
                                (byte *) &kgraphmapstdefaultml.param.coartype,
                                (void *) "hscd" },
                              { KGRAPHMAPSTMETHML,  STRATPARAMINT,    "vert",
                                (byte *) &kgraphmapstdefaultml.param,
                                (byte *) &kgraphmapstdefaultml.param.coarnbr,
                                NULL },
                              { KGRAPHMAPSTMETHML,  STRATPARAMDOUBLE, "rat",
                                (byte *) &kgraphmapstdefaultml.param,
                                (byte *) &kgraphmapstdefaultml.param.coarrat,
                                NULL },
                              { KGRAPHMAPSTMETHRB,  STRATPARAMCASE,   "job",
                                (byte *) &kgraphmapstdefaultrb.param,
                                (byte *) &kgraphmapstdefaultrb.param.flagjobtie,
                                (void *) "ut" },
                              { KGRAPHMAPSTMETHRB,  STRATPARAMDOUBLE, "bal",
                                (byte *) &kgraphmapstdefaultrb.param,
                                (byte *) &kgraphmapstdefaultrb.param.kbalval,
                                NULL },
                              { KGRAPHMAPSTMETHRB,  STRATPARAMCASE,   "map",
                                (byte *) &kgraphmapstdefaultrb.param,
                                (byte *) &kgraphmapstdefaultrb.param.flagmaptie,
                                (void *) "ut" },
                              { KGRAPHMAPSTMETHRB,  STRATPARAMCASE,   "poli",
                                (byte *) &kgraphmapstdefaultrb.param,
                                (byte *) &kgraphmapstdefaultrb.param.polival,
                                (void *) "rls LS" },
                              { KGRAPHMAPSTMETHRB,  STRATPARAMSTRAT,  "sep",
                                (byte *) &kgraphmapstdefaultrb.param,
                                (byte *) &kgraphmapstdefaultrb.param.strat,
                                (void *) &bgraphbipartststratab },
                              { KGRAPHMAPSTMETHNBR, STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

static StratParamTab        kgraphmapstcondtab[] = { /* Graph condition parameter table */
                              { STRATNODENBR,        STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

StratTab                    kgraphmapststratab = { /* Strategy tables for graph mapping methods */
                              kgraphmapstmethtab,
                              kgraphmapstparatab,
                              kgraphmapstcondtab };

/****************************************/
/*                                      */
/* This is the generic mapping routine. */
/*                                      */
/****************************************/

/* This routine computes the given
** mapping according to the given
** strategy.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kgraphMapSt (
Kgraph * restrict const       grafptr,            /*+ Mapping graph    +*/
const Strat * restrict const  strat)              /*+ Mapping strategy +*/
{
  StratTest           val;
  int                 o;

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (sizeof (Gnum) != sizeof (INT)) {
    errorPrint ("kgraphMapSt: invalid type specification for parser variables");
    return     (1);
  }
  if ((sizeof (KgraphMapRbParam) > sizeof (StratNodeMethodData))) {
    errorPrint ("kgraphMapSt: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
#ifdef SCOTCH_DEBUG_KGRAPH1
  if ((strat->tabl != &kgraphmapststratab) &&
      (strat       != &stratdummy)) {
    errorPrint ("kgraphMapSt: invalid parameter (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH1 */

  o = 0;
  switch (strat->type) {
    case STRATNODECONCAT :
      o = kgraphMapSt (grafptr, strat->data.concat.strat[0]); /* Apply first strategy          */
      if (o == 0)                                 /* If it worked all right                    */
        o |= kgraphMapSt (grafptr, strat->data.concat.strat[1]); /* Then apply second strategy */
      break;
    case STRATNODECOND :
      o = stratTestEval (strat->data.cond.test, &val, (void *) grafptr); /* Evaluate expression */
      if (o == 0) {                               /* If evaluation was correct                  */
#ifdef SCOTCH_DEBUG_KGRAPH2
        if ((val.typetest != STRATTESTVAL) ||
            (val.typenode != STRATPARAMLOG)) {
          errorPrint ("kgraphMapSt: invalid test result");
          o = 1;
          break;
        }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
        if (val.data.val.vallog == 1)             /* If expression is true                 */
          o = kgraphMapSt (grafptr, strat->data.cond.strat[0]); /* Apply first strategy    */
        else {                                    /* Else if expression is false           */
          if (strat->data.cond.strat[1] != NULL)  /* And if there is an else statement     */
            o = kgraphMapSt (grafptr, strat->data.cond.strat[1]); /* Apply second strategy */
        }
      }
      break;
    case STRATNODEEMPTY :
      break;
    case STRATNODESELECT :
      errorPrint ("kgraphMapSt: selection operator not implemented for k-way strategies");
      return      (1);
#ifdef SCOTCH_DEBUG_KGRAPH1
    case STRATNODEMETHOD :
#else /* SCOTCH_DEBUG_KGRAPH1 */
    default :
#endif /* SCOTCH_DEBUG_KGRAPH1 */
      return (strat->tabl->methtab[strat->data.method.meth].func (grafptr, (void *) &strat->data.method.data));
#ifdef SCOTCH_DEBUG_KGRAPH1
    default :
      errorPrint ("kgraphMapSt: invalid parameter (2)");
      return     (1);
#endif /* SCOTCH_DEBUG_KGRAPH1 */
  }
  return (o);
}
