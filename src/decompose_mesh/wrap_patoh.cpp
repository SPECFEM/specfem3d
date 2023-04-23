/*
!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
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

// uses PaToH version 3.2
// (version from March, 2011)
// https://www.cc.gatech.edu/~umit/software.html
//
// original written by: Bora Ucar, 16 March 2013 (modified from an earlier wrapPaToH/)

/* ---------------------------------------------------------------------- */

// uncomment to initialize graph partitioner with random seeds
//#define SETRANDSEED

/* ---------------------------------------------------------------------- */

#include "config.h"

#ifndef USE_PATOH
/*
// note: this can lead to compilation and linking issues with intel compilers, complaining about
//       .. undefined reference to `__gxx_personality_v0' ..
//       when linking with a c-compiler instead of a c++-compiler or if -lstdc++ is missing as linker argument
//
// stub routine (used for compilation without patoh libraries)
extern "C"
void
FC_FUNC_(wrap_patoh, WRAP_PATOH)(int* __c, int* __n,
                                 int* nwghts,
                                 int *xpins, int *pins,
                                 int *cwghts,
                                 int* __k, int* __nconst,
                                 int *partvec,
                                 int* __cuttype,
                                 int* __cut) {
  fprintf(stderr,"PATOH Error: called without PATOH Support. To enable PATOH support, please re-compile with flag -DUSE_PATOH.\n");
  exit(1);
}

*/

#else

/* ---------------------------------------------------------------------- */
#pragma message ("\nCompiling with: USE_PATOH enabled\n")

// PaToH function
#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <iostream>

// libconfig: installed from libconfig-hr version 1.4.9
// see:  http://www.hyperrealm.com/libconfig/
#include <libconfig.h++>

#include "patoh.h"

using namespace std;
using namespace libconfig;

/* ---------------------------------------------------------------------- */

void loadPatohParameters(PaToH_Parameters &args);
void writePatohParameters(PaToH_Parameters &args);

#define LOADCFG(val) do {                                               \
    try {                                                               \
      args.val = cfg.lookup(#val);                                      \
    }                                                                   \
    catch(const SettingNotFoundException &nfex)                         \
    {                                                                   \
        cerr << "No " << #val << " setting in configuration file." << endl; \
        exit(EXIT_FAILURE);                                             \
    }                                                                   \
    catch(const SettingTypeException &nfex)                             \
    {                                                                   \
        cerr << "Setting Type Failed: " << #val << endl;                \
        exit(EXIT_FAILURE);                                             \
    }                                                                   \
  } while (0)

/* ---------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(wrap_patoh, WRAP_PATOH)(int* __c, int* __n,
                                      int* nwghts,
                                      int *xpins, int *pins,
                                      int *cwghts,
                                      int* __k, int* __nconst,
                                      int *partvec,
                                      int* __useFixCells,
                                      int* __cuttype,
                                      int* __cut) {

  PaToH_Parameters args;

  int *partwghts = (int *)NULL;
  double imbal = 0.03;

  int _k = *__k;
  int cuttype = *__cuttype;
  int _c = *__c;
  int _n = *__n;
  int _nconst = *__nconst;
  int useFixCells = *__useFixCells;

  int i,ret;
  int verbose = 1;

// in case METIS and PATOH libraries conflict
//#ifndef USE_METIS

#ifdef SETRANDSEED
  static int patohseed = -1;
  if (patohseed == -1) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    patohseed = tp.tv_sec;
  } else {
    ++patohseed;
  }
#endif

  // version check
#ifndef PATOH_BALANCE_STRICT
  // version 3.0 interface
  printf("PaToH: interface 3.0\n");
  printf("PaToH Error: PaToH implementation requires PaToH version 3.2, please re-compile with new version.\n");
  exit(1);
#endif
  // assumes version 3.2 interface
  printf("PaToH: interface 3.2\n");

/* from manual.pdf:

  int PaToH_Initialize_Parameters(PPaToH_Parameters pargs, int cuttype, int SBProbType);

  options:

    cuttype - must be either PATOH_CUTPART (== 2) for cutnet metric
              or PATOH_CONPART (== 1) for "Connectivity-1" metric

    PATOH_SUGPARAM_QUALITY: if you could afford a little bit more extra time for a little better quality
                            (such as VLSI partitioning), use this value.

*/

  printf("PaToH: initialize parameters\n");
  ret = PaToH_Initialize_Parameters(&args, cuttype, PATOH_SUGPARAM_QUALITY);
  if (ret != 0){ printf("PaToH error: in initializing parameters for PaToH"); exit(12);}

  // sets partitioning parameters
  args._k = _k; // number of desired partitions

  args.MemMul_Pins += 5; // increases memory for pins array
  args.outputdetail = PATOH_OD_LOW; // low verbosity

  // seed: by default, Patoh uses current time as the seed (args.seed == -1) allowing for random initialization.
  //       here, we fix the seed to make results easier to reproduce.
#ifdef SETRANDSEED
  args.seed = patohseed;
#else
  args.seed = 1; // fixing random generator seed
#endif

  // optional parameters to try out:
  //
  // - refinement:
  // by default, the following is used
  //args.crs_alg = PATOH_CRS_ABSHPC; // Absorption Clustering using Pins
  //args.ref_alg = PATOH_REFALG_BFMKL; // One pass BFM followed by one pass BKL is used
  //
  // instead, use first FM before falling back after 500 refinement level
  //args.crs_alg = PATOH_CRS_HCM; // Heavy connectivity matching
  //args.initp_runno = 20; // increases number of initial partititioning runs
  //args.ref_alg = PATOH_REFALG_D_FM; // FM with dynamic locking
  //args.ref_useafter = 500;
  //
  // - balance:
  // by default, Patoh strictly forces balance
  //args.balance = 0;
  //
  // instead, use each bisection with a maximum imbalance of epsilon
  //args.balance = 2;
  //args.init_imbal = 0.06;
  //args.final_imbal = 0.03;
  //args.fast_initbal_mult = 1.5;

  // stores parameters (for backup)
  printf("PaToH: write parameters\n");
  writePatohParameters(args);

  // checks if reading parameters would work
  printf("PaToH: load parameters\n");
  loadPatohParameters(args);

  printf("PaToH: check parameters\n");
  ret = PaToH_Check_User_Parameters(&args, verbose);
  if (ret != 0){ printf("PaToH error: in checking user parameters for PaToH"); exit(12);}

  partwghts = (int *) calloc(_k * _nconst, sizeof(int));

  printf("PaToH: allocating arrays\n");
  ret = PaToH_Alloc(&args, _c, _n, _nconst, cwghts, nwghts, xpins, pins);
  if (ret != 0){ printf("PaToH error: in allocating memory for PaToH"); exit(12);}

  // PaToH interface version 3.2
  printf("PaToH: calling PaToH_Part\n");
  ret = PaToH_Part(&args, _c, _n, _nconst, useFixCells, cwghts, nwghts,
                   xpins, pins, NULL, partvec, partwghts, __cut);
  if (ret != 0){ printf("PaToH error: in PaToH_Part partitioning for PaToH"); exit(12);}

  free(partwghts);

  ret = PaToH_Free();
  if (ret != 0){ printf("PaToH error: in freeing arrays for PaToH"); exit(12);}


// in case METIS and PATOH libraries conflict
//#else
//  printf("PaToH Error: libPATOH conflicts with libMETIS. Recompile without -DUSE_METIS\n");
//  exit(1);
//#endif // USE_METIS

}


/* ---------------------------------------------------------------------- */

void loadPatohParameters(PaToH_Parameters &args) {

  Config cfg;
  // Read the file. If there is an error, report it and exit.
  try {
    cfg.readFile("patoh.cfg");
  }
  catch(const FileIOException &fioex) {
    cerr << "Couldn't read/find configuration file patoh.cfg." << endl;
    exit(EXIT_FAILURE);
  }
  catch(const ParseException &pex) {
    cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
         << " - " << pex.getError() << endl;
    exit(EXIT_FAILURE);
  }

  // PaToH interface version 3.2
  // Get the parameters
  try {
    LOADCFG(cuttype);

    LOADCFG(outputdetail);
    LOADCFG(seed);
    LOADCFG(doinitperm);
    LOADCFG(bisec_fixednetsizetrsh);
    LOADCFG(bisec_netsizetrsh);
    LOADCFG(bisec_partmultnetsizetrsh);
    LOADCFG(bigVcycle);
    LOADCFG(smallVcycle);
    LOADCFG(usesamematchinginVcycles);
    LOADCFG(usebucket);
    LOADCFG(maxcellinheap);
    LOADCFG(heapchk_mul);
    LOADCFG(heapchk_div);
    LOADCFG(MemMul_CellNet);
    LOADCFG(MemMul_Pins);
    LOADCFG(MemMul_General);

    LOADCFG(crs_VisitOrder);
    LOADCFG(crs_alg);
    LOADCFG(crs_coarsento);
    LOADCFG(crs_coarsentokmult);
    LOADCFG(crs_coarsenper);
    LOADCFG(crs_maxallowedcellwmult);
    LOADCFG(crs_idenafter);
    LOADCFG(crs_iden_netsizetrh);
    LOADCFG(crs_useafter);
    LOADCFG(crs_useafteralg);

    LOADCFG(nofinstances);
    LOADCFG(initp_alg);
    LOADCFG(initp_runno);
    LOADCFG(initp_ghg_trybalance);
    LOADCFG(initp_refalg);

    LOADCFG(ref_alg);
    LOADCFG(ref_useafter);
    LOADCFG(ref_useafteralg);
    LOADCFG(ref_passcnt);
    LOADCFG(ref_maxnegmove);
    LOADCFG(ref_maxnegmovemult);
    LOADCFG(ref_dynamiclockcnt);
    LOADCFG(ref_slow_uncoarsening);

    LOADCFG(balance);
    LOADCFG(init_imbal);
    LOADCFG(final_imbal);
    LOADCFG(fast_initbal_mult);
    LOADCFG(init_sol_discard_mult);
    LOADCFG(final_sol_discard_mult);

    // example:
    // args.initp_runno = cfg.lookup("initp_runno"); // 100 helps
  }
  catch(const SettingNotFoundException &nfex) {
    cerr << "No 'name' setting in configuration file." << endl;
    exit(EXIT_FAILURE);
  }

}

/* ---------------------------------------------------------------------- */

void writePatohParameters(PaToH_Parameters &args) {

  Config cfg;

  Setting &root = cfg.getRoot();

  // note: libconfig doesn't support writing comment lines (#.. or //..)
  //       this comment is added as a configuration item.
  root.add("comment", Setting::TypeString) = "#PATOH configuration - backup file";

  // set the parameters
  root.add("cuttype",Setting::TypeInt) = args.cuttype;
  root.add("k",Setting::TypeInt) = args._k;

  // miscellaneous parameters
  root.add("outputdetail",Setting::TypeInt) = args.outputdetail;
  root.add("seed",Setting::TypeInt) = args.seed;
  root.add("doinitperm",Setting::TypeInt) = args.doinitperm;

  root.add("bisec_fixednetsizetrsh",Setting::TypeInt) = args.bisec_fixednetsizetrsh;
  root.add("bisec_netsizetrsh",Setting::TypeFloat) = args.bisec_netsizetrsh;
  root.add("bisec_partmultnetsizetrsh",Setting::TypeInt) = args.bisec_partmultnetsizetrsh;

  root.add("bigVcycle",Setting::TypeInt) = args.bigVcycle;
  root.add("smallVcycle",Setting::TypeInt) = args.smallVcycle;
  root.add("usesamematchinginVcycles",Setting::TypeInt) = args.usesamematchinginVcycles;

  root.add("usebucket",Setting::TypeInt) = args.usebucket;
  root.add("maxcellinheap",Setting::TypeInt) = args.maxcellinheap;
  root.add("heapchk_mul",Setting::TypeInt) = args.heapchk_mul;
  root.add("heapchk_div",Setting::TypeInt) = args.heapchk_div;

  root.add("MemMul_CellNet",Setting::TypeInt) = args.MemMul_CellNet;
  root.add("MemMul_Pins",Setting::TypeInt) = args.MemMul_Pins;
  root.add("MemMul_General",Setting::TypeInt) = args.MemMul_General;

  // coarsening parameters
  root.add("crs_VisitOrder",Setting::TypeInt) = args.crs_VisitOrder;
  root.add("crs_alg",Setting::TypeInt) = args.crs_alg;
  root.add("crs_coarsento",Setting::TypeInt) = args.crs_coarsento;
  root.add("crs_coarsentokmult",Setting::TypeInt) = args.crs_coarsentokmult;
  root.add("crs_coarsenper",Setting::TypeInt) = args.crs_coarsenper;
  root.add("crs_maxallowedcellwmult",Setting::TypeFloat) = args.crs_maxallowedcellwmult;
  root.add("crs_idenafter",Setting::TypeInt) = args.crs_idenafter;
  root.add("crs_iden_netsizetrh",Setting::TypeInt) = args.crs_iden_netsizetrh;
  root.add("crs_useafter",Setting::TypeInt) = args.crs_useafter;
  root.add("crs_useafteralg",Setting::TypeInt) = args.crs_useafteralg;

  // initial partitioning parameters
  root.add("nofinstances",Setting::TypeInt) = args.nofinstances;
  root.add("initp_alg",Setting::TypeInt) = args.initp_alg;
  root.add("initp_runno",Setting::TypeInt) = args.initp_runno;
  root.add("initp_ghg_trybalance",Setting::TypeInt) = args.initp_ghg_trybalance;
  root.add("initp_refalg",Setting::TypeInt) = args.initp_refalg;

  // uncoarsening parameters
  root.add("ref_alg",Setting::TypeInt) = args.ref_alg;
  root.add("ref_useafter",Setting::TypeInt) = args.ref_useafter;
  root.add("ref_useafteralg",Setting::TypeInt) = args.ref_useafteralg;
  root.add("ref_passcnt",Setting::TypeInt) = args.ref_passcnt;
  root.add("ref_maxnegmove",Setting::TypeInt) = args.ref_maxnegmove;
  root.add("ref_maxnegmovemult",Setting::TypeFloat) = args.ref_maxnegmovemult;
  root.add("ref_dynamiclockcnt",Setting::TypeInt) = args.ref_dynamiclockcnt;
  root.add("ref_slow_uncoarsening",Setting::TypeFloat) = args.ref_slow_uncoarsening;

  root.add("balance",Setting::TypeInt) = args.balance;
  root.add("init_imbal",Setting::TypeFloat) = args.init_imbal;
  root.add("final_imbal",Setting::TypeFloat) = args.final_imbal;
  root.add("fast_initbal_mult",Setting::TypeFloat) = args.fast_initbal_mult;
  root.add("init_sol_discard_mult",Setting::TypeFloat) = args.init_sol_discard_mult;
  root.add("final_sol_discard_mult",Setting::TypeFloat) = args.final_sol_discard_mult;

  // write the file. If there is an error, report it and exit.
  try {
    cfg.writeFile("patoh.cfg");
  }
  catch(const FileIOException &fioex) {
    cerr << "Couldn't write configuration file patoh.cfg" << endl;
    exit(EXIT_FAILURE);
  }
}


#endif //USE_PATOH
