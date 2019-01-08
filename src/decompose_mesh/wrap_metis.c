/*
!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
*/

// METIS partitioning
//
// uses multi-level K-way partitioning scheme

#include "config.h"
#include <stdlib.h>
#include <stdio.h>

#ifdef USE_METIS

#pragma message ("\nCompiling with: USE_METIS enabled\n")

// metis
#include <metis.h>
// macro check
#if METIS_VER_MAJOR >= 5
// METIS version 5.x
#define METIS_SAFE_CALL(call) do {                                       \
        int METIS_err = call;                                         \
        if (METIS_err != METIS_OK) {                                       \
            fprintf (stderr, "Metis Error in file '%s' in line %i\n", \
                     __FILE__, __LINE__ );     \
            exit(EXIT_FAILURE);                                         \
        }                                                               \
    } while (0)
#endif

#endif // USE_METIS

// wrapper for metis call
void
FC_FUNC_(wrap_metis, WRAP_METIS)(int* nspec, int* ncon, int* xadj, int* adjncy,int* vwgt,int* adjwgt,
                                 int* nparts,float* ubvec, int* edgecut,int* part) {

    // PATOH and METIS libs conflict
#ifdef USE_METIS

#if METIS_VER_MAJOR >= 5
    // METIS version 5

    // note: METIS v5.0 is nearly a complete re-write and provides better support for 64-bit architectures.
    //       As a result, the API routines have changed, making it incompatible with earlier versions.
    //
    // see: glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf

    // METIS version 5.x format
    //METIS_API(int) METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj,
    //                                   idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt,
    //                                   idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options,
    //                                   idx_t *edgecut, idx_t *part);
    printf("Metis: Calling METIS version v5 PartGraphKway\n");

    // options
    int options[METIS_NOPTIONS];
    METIS_SAFE_CALL(METIS_SetDefaultOptions(options));

    // to avoid error due to different bit-size of float and real_t (check with #ifdef METIS_USE_DOUBLEPRECISION)
    real_t *vec = (real_t*) malloc(*ncon * sizeof(real_t));
    int i;
    for(i=0; i<*ncon; i++){
      vec[i] = ubvec[i];
      //debug
      //printf("ubvec %i %f %f\n",i,ubvec[i],vec[i]);
    }

    // Kway partitioning scheme
    if (*ncon > 1){
      // multi-constraint
      METIS_SAFE_CALL(METIS_PartGraphKway(nspec,ncon,xadj,adjncy,vwgt,NULL,adjwgt,nparts,NULL,vec,options,edgecut,part));
    }else{
      // single constraint
      METIS_SAFE_CALL(METIS_PartGraphKway(nspec,ncon,xadj,adjncy,vwgt,NULL,NULL,nparts,NULL,NULL,options,edgecut,part));
    }
    free(vec);

#else
    // METIS version 4

    // unfortunately, METIS version 4 has no definition like METIS_VER_MAJOR to identify its source version.
    // we assume version 4.0.3 provided in external_libs/

    /*
    version 4 format:
      wgtflag - Used to indicate if the graph is weighted.
               wgtflag can take the following values:
               0 No weights (vwgts and adjwgt are NULL)
               1 Weights on the edges only (vwgts = NULL)
               2 Weights on the vertices only (adjwgt = NULL)
               3 Weights both on vertices and edges.

      numflag - Used to indicate which numbering scheme is used for the adjacency structure of the graph.
               numflag can take the following two values:
               0 C-style numbering is assumed that starts from 0
               1 Fortran-style numbering is assumed that starts from 1

      options - This is an array of 5 integers that is used to pass parameters for the various phases of the algorithm.
               If options[0]=0 then default values are used.
    */

    if (*ncon > 1){
      // version 4, multi-constraint partitioning
      //void METIS mCPartGraphKway (int *n, int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
      //                            idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *ubvec, int *options,
      //                            int *edgecut, idxtype *part)
      printf("Metis: Calling METIS version v4 mCPartGraphKway\n");

      if (*ncon > 15){printf("Warning: ncon shouldn't be bigger than 15 for version 4 - currently it is %d\n", *ncon);}

      int wgtflag = 1; // 1 == Weights on the edges (adjwgt not NULL).
      int numflag = 0;
      int options[5] = { 0, 0, 0, 0, 0};

      METIS_mCPartGraphKway(nspec,ncon,xadj,adjncy,vwgt,adjwgt,&wgtflag,&numflag,nparts,ubvec,options,edgecut,part);
    }else{
      // version 4, single constraint
      //void METIS_PartGraphKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
      //                         idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
      //                         int *options, int *edgecut, idxtype *part);
      // no return value, but will have edgecut == -1 if failed...
      printf("Metis: Calling METIS version v4 PartGraphKway\n");
      int wgtflag = 2; // Weights on the vertices only (adjwgt = NULL)
      int numflag = 0;
      int options[5] = { 0, 0, 0, 0, 0};

      METIS_PartGraphKway(nspec,xadj,adjncy,vwgt,NULL,&wgtflag,&numflag,nparts,options,edgecut,part);
    }
    // checks for error
    if (*edgecut < 0){
      fprintf(stderr,"Metis Error: version v4 negative edgecut error %d\n", *edgecut);
      exit(EXIT_FAILURE);
    }
#endif

#else
    // compiled without support
    fprintf(stderr,"METIS Error: Please re-compile with -DUSE_METIS. Probably just need to make clean; make xdecompose_mesh.\n");
    exit(1);
#endif
}
