dnl                                               -*- Autoconf -*-
dnl
dnl Process this file with autoconf to produce the 'configure' script:
dnl
dnl autoreconf -i
dnl
dnl You must have recent versions of Autoconf and Automake installed.

############################################################

AC_PREREQ(2.61)
AC_INIT([Specfem3D], m4_normalize(m4_include([VERSION])), [see the wiki], [Specfem3D])
AC_CONFIG_SRCDIR([src/specfem3D/specfem3D.F90])

AC_CONFIG_HEADERS([setup/config.h])

AC_CONFIG_MACRO_DIR([m4])
m4_include(m4/cit_backports.m4)

############################################################

AS_BOX([setting up compilation flags])
AC_CANONICAL_HOST

# 'configure' options

###
### precision
###

AC_ARG_ENABLE([double-precision],
    [AS_HELP_STRING([--enable-double-precision],
        [solver in double precision @<:@default=no@:>@])],
    [want_double_precision="$enableval"],
    [want_double_precision=no])
if test x"$want_double_precision" = xno; then
    CUSTOM_REAL=SIZE_REAL
    CUSTOM_MPI_TYPE=MPI_REAL
else
    CUSTOM_REAL=SIZE_DOUBLE
    CUSTOM_MPI_TYPE=MPI_DOUBLE_PRECISION
fi
AC_SUBST([CUSTOM_REAL])
AC_SUBST([CUSTOM_MPI_TYPE])

###
### debugging
###

AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [build with debugging options enabled @<:@default=no@:>@])],
    [want_debug="$enableval"],
    [want_debug=no])
AM_CONDITIONAL([COND_DEBUG], [test x"$want_debug" != xno])
export COND_DEBUG_FALSE
export COND_DEBUG_TRUE

###
### MPI
###

AC_ARG_WITH([mpi],
    [AS_HELP_STRING([--with-mpi],
        [build parallel version @<:@default=yes@:>@])],
    [want_mpi="$withval"],
    [want_mpi=yes])
AM_CONDITIONAL([COND_MPI], [test x"$want_mpi" = xyes])

###
### FORCE_VECTORIZATION
###

## note: vectorization support has been removed by DK due to problems with Intel compilers (2016).
##       enabling will have no effect.

AC_ARG_ENABLE([vectorization],
    [AS_HELP_STRING([--enable-vectorization],
        [build FORCE_VECTORIZATION enabled version @<:@default=auto@:>@])],
    [want_vec="$enableval"],
    [want_vec=no])
if test x"$want_vec" == xauto; then
    if test x"$want_debug" != xno; then
        want_vec=no
    else
        want_vec=yes
    fi
else
    if test x"$want_vec" != xno && test x"$want_debug" != xno; then
        AC_MSG_ERROR([--enable-debug and --enable-vectorization cannot be specified simultaneously.])
    fi
fi
AM_CONDITIONAL([COND_VECTORIZATION], [test x"$want_vec" != xno])

###
### GPU version
###

## CUDA
AC_ARG_WITH([cuda],
    [AS_HELP_STRING([--with-cuda],
        [build CUDA GPU enabled version @<:@default=no@:>@])],
    [want_cuda="$withval"],
    [want_cuda=no])

# sets default: --with-cuda
if test x"$want_cuda" = xyes; then
    want_cuda=cuda4
fi
AM_CONDITIONAL([COND_CUDA], [test x"$want_cuda" != xno])
AM_CONDITIONAL([COND_CUDA4], [test x"$want_cuda" = xcuda4])
AM_CONDITIONAL([COND_CUDA5], [test x"$want_cuda" = xcuda5])
AM_CONDITIONAL([COND_CUDA6], [test x"$want_cuda" = xcuda6])
AM_CONDITIONAL([COND_CUDA7], [test x"$want_cuda" = xcuda7])
AM_CONDITIONAL([COND_CUDA8], [test x"$want_cuda" = xcuda8])
AM_CONDITIONAL([COND_CUDA9], [test x"$want_cuda" = xcuda9])
AM_CONDITIONAL([COND_CUDA10], [test x"$want_cuda" = xcuda10])
AM_CONDITIONAL([COND_CUDA11], [test x"$want_cuda" = xcuda11])

# cuda linking for cuda 5x and 6x and 7x and 8x and ..
AM_CONDITIONAL([COND_CUDA_PLUS],
  [test "$want_cuda" = cuda4 \
     -o "$want_cuda" = cuda5 \
     -o "$want_cuda" = cuda6 \
     -o "$want_cuda" = cuda7 \
     -o "$want_cuda" = cuda8 \
     -o "$want_cuda" = cuda9 \
     -o "$want_cuda" = cuda10 \
     -o "$want_cuda" = cuda11 \]
)

## HIP
AC_ARG_WITH([hip],
 [AS_HELP_STRING([--with-hip], [build HIP GPU enabled version @<:@default=no@:>@])],
    [want_hip="$withval"],
    [want_hip=no])
AM_CONDITIONAL([COND_HIP], [test x"$want_hip" != xno])
AM_CONDITIONAL([COND_HIP_MI8], [test x"$want_hip" = xMI8])
AM_CONDITIONAL([COND_HIP_MI25], [test x"$want_hip" = xMI25])
AM_CONDITIONAL([COND_HIP_MI50], [test x"$want_hip" = xMI50])
AM_CONDITIONAL([COND_HIP_MI100], [test x"$want_hip" = xMI100])
AM_CONDITIONAL([COND_HIP_MI250], [test x"$want_hip" = xMI250])

AM_CONDITIONAL([COND_HIP_CUDA5], [test x"$want_hip" = xcuda5])
AM_CONDITIONAL([COND_HIP_CUDA6], [test x"$want_hip" = xcuda6])
AM_CONDITIONAL([COND_HIP_CUDA7], [test x"$want_hip" = xcuda7])
AM_CONDITIONAL([COND_HIP_CUDA8], [test x"$want_hip" = xcuda8])
AM_CONDITIONAL([COND_HIP_CUDA9], [test x"$want_hip" = xcuda9])
AM_CONDITIONAL([COND_HIP_CUDA10], [test x"$want_hip" = xcuda10])
AM_CONDITIONAL([COND_HIP_CUDA11], [test x"$want_hip" = xcuda11])

###
### OpenMP
###

AC_ARG_WITH([openmp],
    [AS_HELP_STRING([--with-openmp],
        [build OpenMP enabled version @<:@default=no@:>@])],
    [with_omp="$withval"],
    [with_omp=no])
AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--enable-openmp],
        [build OpenMP enabled version @<:@default=no@:>@])],
    [enable_omp="$enableval"],
    [enable_omp=no])
AS_IF([test x"$enable_omp" != xno || test x"$with_omp" != xno],
    [want_omp=yes],[want_omp=no])
AM_CONDITIONAL([COND_OMP], [test x"$want_omp" != xno])

###
### VTK
###

AC_ARG_ENABLE([vtk],
    [AS_HELP_STRING([--enable-vtk],
        [build VTK enabled version @<:@default=no@:>@])],
    [want_vtk="$enableval"],
    [want_vtk=no])
AM_CONDITIONAL([COND_VTK], [test x"$want_vtk" != xno])

###
### ADIOS
###

AC_ARG_WITH([adios],
    [AS_HELP_STRING([--with-adios],
        [build ADIOS enabled version @<:@default=no@:>@])],
    [want_adios="$withval"],
    [want_adios=no])
AM_CONDITIONAL([COND_ADIOS], [test x"$want_adios" != xno])

###
### ADIOS2
###

AC_ARG_WITH([adios2],
    [AS_HELP_STRING([--with-adios2],
        [build ADIOS enabled version @<:@default=no@:>@])],
    [want_adios2="$withval"],
    [want_adios2=no])
AM_CONDITIONAL([COND_ADIOS2], [test x"$want_adios2" != xno])

###
### ASDF
###

AC_ARG_WITH([asdf],
    [AS_HELP_STRING([--with-asdf],
        [build ASDF enabled version @<:@default=no@:>@])],
    [want_asdf="$withval"],
    [want_asdf=no])
AM_CONDITIONAL([COND_ASDF], [test x$"$want_asdf" != xno])

############################################################

# Checks for programs.

# a courtesy to the installed base of users
if test x"$FC" = x && test x"$F90" != x; then
    FC="$F90"
fi
if test x"$MPIFC" = x && test x"$MPIF90" != x; then
    MPIFC="$MPIF90"
fi
# note: if not set, default FCFLAGS (and CFLAGS) will be set by autoconf/configure script (mostly -g -O2)
if test x"$FCFLAGS" = x && test x"$FLAGS_CHECK" != x ; then
    FCFLAGS=""
fi

AC_PROG_FC
export FC
export MPIFC
F77="$FC"
FFLAGS="$FCFLAGS"

AC_PROVIDE([AC_PROG_F77])
AC_SUBST([FCENV])
AC_SUBST(srcdir)

flags_guess="$SHELL $srcdir/flags.guess"
AC_MSG_NOTICE([running $flags_guess])
flags=`$flags_guess` ||
    AC_MSG_ERROR([$flags_guess failed])
eval $flags

AC_FC_WRAPPERS
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)
AC_FC_FREEFORM()

AC_FC_PP_DEFINE()
AC_SUBST([FC_DEFINE])
AC_FC_PP_SRCEXT(F90) dnl Because AC_FC_PP_DEFINE messes with things.

AC_PROG_CC

############################################################

#checks for Scotch
AS_BOX([SCOTCH])

export SCOTCH_DIR
export SCOTCH_LIBDIR
export SCOTCH_INCLUDEDIR
export USE_BUNDLED_SCOTCH

AC_ARG_WITH([scotch-dir],
AS_HELP_STRING([--with-scotch-dir=DIR],[define the root path to Scotch (e.g. /opt/scotch/)]),
[
 ac_scotch_dir="$withval";
])

AC_ARG_WITH([scotch-includedir],
AS_HELP_STRING([--with-scotch-includedir=DIR],[define the path to the Scotch headers (e.g. /opt/scotch/include)]),
[
 ac_scotch_include_dir="$withval";
])

AC_ARG_WITH([scotch-libdir],
AS_HELP_STRING([--with-scotch-libdir=DIR],[define the path to the Scotch libraries (e.g. /opt/scotch/lib)]),
[
 ac_scotch_lib_dir="$withval";
])

if test "${USE_BUNDLED_SCOTCH}" != "1"; then

   if test -z "${ac_scotch_lib_dir}"; then
      if test -n "${ac_scotch_dir}"; then
         ac_scotch_lib_dir="${ac_scotch_dir}/lib";
      else
        ac_scotch_lib_dir="/usr/lib";
      fi
   fi

  if test -z "${ac_scotch_include_dir}"; then
    if test -n "${ac_scotch_dir}"; then
      ac_scotch_include_dir="${ac_scotch_dir}/include";
    else
      ac_scotch_include_dir="/usr/include/scotch";
    fi
  fi

  scotch_lib=""
  ac_save_ldflags=${LDFLAGS}
  if test -n "${ac_scotch_lib_dir}"; then
    LDFLAGS="${LDFLAGS} -L${ac_scotch_lib_dir}";
  fi
  AC_CHECK_LIB(scotch,scotchfarchinit ,[scotch_lib="yes";],[scotch_lib="no";LDFLAGS=${ac_save_ldflags}], -lscotcherr)

  scotch_include=""
  ac_save_cppflags=${CPPFLAGS}
  # scotch only needed for fortran code
  #if test -n "${ac_scotch_include_dir}"; then
  #  CPPFLAGS="${CPPFLAGS} -I${ac_scotch_include_dir}";
  #fi

dnl This does not work because we are in Fortran mode, which does not
dnl handle headers.
dnl
dnl AC_CHECK_HEADER(scotchf.h,[scotch_include="yes";],[scotch_include="no";CPPFLAGS=${ac_save_cppflags}])

  scotch_usable="${scotch_lib}"

fi

AC_MSG_CHECKING([whether Scotch is usable])
if test "x${scotch_usable}" = "xyes"; then
 AC_DEFINE([HAVE_SCOTCH],[1],[defined if Scotch is installed])
 want_scotch=yes
 USE_BUNDLED_SCOTCH=0
 SCOTCH_DIR="${ac_scotch_dir}"
 SCOTCH_LIBDIR="${ac_scotch_lib_dir}"
 SCOTCH_INCLUDEDIR="${ac_scotch_include_dir}"
 AC_MSG_RESULT([yes])
else
 AC_DEFINE([HAVE_SCOTCH],[1],[defined if Scotch is installed])
 AC_MSG_RESULT([no, using bundled scotch instead])
 AC_PROG_LEX
 if test -z "$LEX" || test "X$LEX" = "Xno"; then
   AC_MSG_ERROR([No suitable lex found])
 fi
 AC_PROG_YACC
 if test -z "$YACC" || test "X$YACC" = "Xno"; then
   AC_MSG_ERROR([No suitable yacc or bison found])
 fi
 ACX_PTHREAD(AC_MSG_RESULT([pthread found]), AC_MSG_ERROR([pthread not found]))

 # scotch only needed with mpi support
 if test "$want_mpi" = yes; then
  want_scotch=yes
  #daniel: scotch bundle
  # uses bundled scotch: current version pointed to by symbolic link scotch/
  USE_BUNDLED_SCOTCH=1

  SCOTCH_DIR="$srcdir/external_libs/scotch"
  SCOTCH_LIBDIR="${SCOTCH_DIR}/lib"
  SCOTCH_INCLUDEDIR="${SCOTCH_DIR}/include"
 else
  # no mpi support, scotch not needed
  want_scotch=no
  USE_BUNDLED_SCOTCH=0
 fi # want_mpi

fi
LDFLAGS=${ac_save_ldflags}

AM_CONDITIONAL([COND_SCOTCH], [test "$want_scotch" = yes])
AS_IF([test "$want_scotch" = yes],[
  AC_MSG_RESULT([SCOTCH is enabled])
],[
  AC_MSG_RESULT([SCOTCH is not enabled])
])


############################################################
# influential environment variables

AC_ARG_VAR(USE_BUNDLED_SCOTCH, [Set to 1 to always use the bundled Scotch library])
AC_ARG_VAR(SCOTCH_DIR, [Directory where Scotch is installed])
AC_ARG_VAR(SCOTCH_INCLUDEDIR, [Directory where Scotch headers are installed])
AC_ARG_VAR(SCOTCH_LIBDIR, [Directory where Scotch libraries are installed])

AC_ARG_VAR(MPIFC, [MPI Fortran compiler command])
AC_ARG_VAR(MPILIBS, [extra libraries for linking MPI programs])
AC_ARG_VAR(MPICC, [MPI C compiler command])

AC_ARG_VAR(FLAGS_CHECK, [Fortran compiler flags])

# scratch disks
AC_ARG_VAR(LOCAL_PATH_IS_ALSO_GLOBAL, [files on a local path on each node are also seen as global with same path @<:@default=true@:>@])

############################################################

# tests
AS_IF([test x"$MPIFC" = x],[
    MPIFC=mpif90
])

AS_IF([test x"$MPICC" = x],[
    MPICC=mpicc
])

AS_IF([test x"$LOCAL_PATH_IS_ALSO_GLOBAL" = x],[
    LOCAL_PATH_IS_ALSO_GLOBAL=true
])

# Checks for typedefs, structures, and compiler characteristics.
# check fortran modules flag
AS_BOX([module extensions])

AC_FC_MODULE_EXTENSION
AC_SUBST([FC_MODEXT])
AC_FC_MODULE_FLAG([], [])
AC_FC_MODULE_OUTPUT_FLAG([
	FCFLAGS_f90="$FC_MODOUT./obj $FC_MODINC./obj $FC_MODINC. $FCFLAGS_f90"
	FC_MODDIR=./obj
], [
	FC_MODDIR=.
])
AC_SUBST([FC_MODDIR])

# Checks for header files.

AC_LANG_PUSH(C)
AC_CHECK_HEADER(emmintrin.h,AC_DEFINE([HAVE_EMMINTRIN],[1],[Define if emmintrin.h]))
AC_CHECK_HEADER(xmmintrin.h,AC_DEFINE([HAVE_XMMINTRIN],[1],[Define if xmmintrin.h]))
# AIX doesn't have err.h so we need to check
AC_CHECK_HEADER(err.h,      AC_DEFINE([HAVE_ERR],      [1],[Define if err.h]))
AC_LANG_POP(C)

# Checks for libraries.

###
### MPI
###

AS_IF([test "$want_mpi" = yes], [
    AS_BOX([MPI])
    # checks MPI include directory
    CIT_MPI_INCDIR([$MPIFC])
])

###
### ADIOS
###

AS_IF([test x"$want_adios" != xno], [
    AS_BOX([ADIOS])
    if test "$want_mpi" != yes; then
        AC_MSG_ERROR([Cannot compile with ADIOS without MPI.])
    fi
    CIT_ADIOS_CONFIG
])

###
### ADIOS2
###

AS_IF([test x"$want_adios2" != xno], [
    AS_BOX([ADIOS2])
    if test "$want_mpi" != yes; then
        AC_MSG_ERROR([Cannot compile with ADIOS2 without MPI.])
    fi

    # for fortran compilation
    CIT_ADIOS2_CONFIG

    # for C file compilation
    AC_MSG_CHECKING([for ADIOS2 C flags and libs])
    AS_IF([test x"$ADIOS2_CONFIG" != x], [
      ADIOS2_CFLAGS=`$ADIOS2_CONFIG --c-flags`
      ADIOS2_CLIBS=`$ADIOS2_CONFIG --c-libs`
      AC_MSG_RESULT(done)
    ],[AC_MSG_RESULT(no adios2_config)])
    AC_SUBST([ADIOS2_CFLAGS])
    AC_SUBST([ADIOS2_CLIBS])

    # in case adios2 c file compilation needs mpi
    AC_SUBST([MPICC])
])

###
### ASDF
###

AS_IF([test x"$want_asdf" != xno], [
    AS_BOX([ASDF])
    AC_ARG_VAR(ASDF_LIBS, [ASDF libraries for linking programs])
])

#mpbl: adding flags to config.h.in through autoheader in order to save them into the adios file.
#configure_flags_str="FC=$FC FCFLAGS=$FCFLAGS $FLAGS_CHECK"
#AC_DEFINE_UNQUOTED([CONFIGURE_FLAGS],["${configure_flags_str}"],
#         [configuration flags to be saved in the adios output file.])
CONFIGURE_FLAGS="FC=$FC FCFLAGS=$FCFLAGS"
AC_SUBST([CONFIGURE_FLAGS])

###
### GPU
###

## CUDA
AS_IF([test x"$want_cuda" != xno], [
    AS_BOX([CUDA])
    CIT_CUDA_CONFIG
])

## HIP
AS_IF([test x"$want_hip" != xno], [
    AS_BOX([HIP])
    # environment variables
    AC_ARG_VAR(HIPCC, [AMD HIP compiler command])
    AC_ARG_VAR(HIP_FLAGS, [HIP compiler flags])
    AC_ARG_VAR(HIP_INC, [Location of HIP include files])
    AC_ARG_VAR(HIP_LIB, [Location of HIP library libhip_hcc])

    # tests HIPCC variable
    AS_IF([test x"$HIPCC" = x],[
        HIPCC=hipcc
    ])

    # Check for compiler
    # checks if program in path
    AC_PATH_PROG(HIPCC_PROG, $HIPCC)
    if test -z "$HIPCC_PROG" ; then
        AC_MSG_ERROR([cannot find '$HIPCC' program, please check your PATH.])
    fi

    # sets default HIP path
    AC_PATH_PROG(HIPCONFIG_PROG, hipconfig)
    AS_IF([test x"$HIPCONFIG_PROG" = x],[
      AC_PATH_PROG(HIPCONFIG_PROG, hipconfig, path=/opt/rocm/bin)
    ])
    AS_IF([test x"$HIPCONFIG_PROG" != x],[
        HIP_PATH=`$HIPCONFIG_PROG --path`
    ])
    AC_MSG_NOTICE([hipconfig program: $HIPCONFIG_PROG])
    AC_MSG_NOTICE([HIP path: $HIP_PATH])

    #AS_IF([test x"$HIP_INC" = x],[
    #    HIP_INC="${HIP_PATH}/include"
    #])

    # adds default HIP library
    AS_IF([test x"$HIP_LIB" = x -a x"${HIP_PATH}" != x],[
        HIP_LIB="${HIP_PATH}/lib"
        AC_MSG_NOTICE([adding default HIP library path: ${HIP_LIB}])
    ])

    # checks platform
    if test x"$HIPCONFIG_PROG" != x ; then
      GPU_PLATFORM=`$HIPCONFIG_PROG --platform`
      AS_IF([test x"$GPU_PLATFORM" = xnvidia], [
          AC_MSG_NOTICE([HIP PLATFORM NVIDIA detected.])
          HIP_PLATFORM=nvcc
          FLAG_PLATFORM="${FC_DEFINE}__HIP_PLATFORM_NVCC__"
      ])

      AS_IF([test x"$GPU_PLATFORM" = xamd], [
          AC_MSG_NOTICE([HIP PLATFORM AMD detected.])
          HIP_PLATFORM=hcc
          FLAG_PLATFORM="${FC_DEFINE}__HIP_PLATFORM_HCC__"
      ])
      # compiler might require a platform selection flag
      case "${HIP_FLAGS}" in
        *HIP_PLATFORM*) AC_MSG_NOTICE([hip flags contain platform specifier: ${HIP_FLAGS}]) ;;
        *) AC_MSG_NOTICE([adding flag for platform specifier ${FLAG_PLATFORM}]);
           HIP_FLAGS="${HIP_FLAGS} ${FLAG_PLATFORM}"
           AC_MSG_NOTICE([HIP_FLAGS: ${HIP_FLAGS}]) ;;
      esac
    fi
    # sets default platform to AMD if not set
    AS_IF([test x"$GPU_PLATFORM" = x],[
      GPU_PLATFORM=amd
    ])

    # for compilation errors like:
    #   /usr/bin/ld: obj/assemble_MPI_scalar_gpu.hip.o: relocation R_X86_64_32 against `.rodata.str1.1' can not be used
    #      when making a PIE object; recompile with -fPIC
    # try to add
    #   HIP_FLAGS="-fPIC"
    # or
    #   HIP_FLAGS="-fno-gpu-rdc"
    # compiler might require a -fPIC flag
    case "${HIP_FLAGS}" in
      *fPIC*) AC_MSG_NOTICE([HIP flags contain position independent code flag -fPIC: ${HIP_FLAGS}]);;
      *) AC_MSG_NOTICE([consider adding compilation flag for position independent code: HIP_FLAGS="-fPIC"]);;
    esac

    # test lib & include
    AS_IF([test x"$HIP_LIB" != x],[
        HIP_LDFLAGS="-L$HIP_LIB"
    ])
    AS_IF([test x"$HIP_INC" != x],[
        HIP_CPPFLAGS="-I$HIP_INC"
    ])
    HIP_LIBS="-lamdhip64"

    # checks header file
    AC_LANG_PUSH(C)
    AC_REQUIRE_CPP
    CFLAGS_save="$CFLAGS"
    LDFLAGS_save="$LDFLAGS"
    LIBS_save="$LIBS"

    # uses hipcc compiler
    CFLAGS="$HIP_FLAGS"
    if test "x$HIP_INC" != "x"; then
      HIP_CPPFLAGS="-I$HIP_INC"
      CFLAGS="$CFLAGS $HIP_CPPFLAGS"
    fi

    AC_MSG_CHECKING([for hip/hip_runtime.h])
    ac_compile='$HIPCC -c $CFLAGS conftest.$ac_ext >&5'
    AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([[#include <hip/hip_runtime.h>]],[[void* ptr = 0;]])
    ], [
        AC_MSG_RESULT(yes)
    ], [
        AC_MSG_RESULT(no)
        AC_MSG_ERROR([HIP runtime header not found; try setting HIP_INC.])
    ])

    # checks for HIP library
    if test "x$HIP_LIB" != "x"; then
      HIP_LDFLAGS="-L$HIP_LIB"
      LDFLAGS="$HIP_LDFLAGS $LDFLAGS"
    fi
    HIP_LIBS="-lamdhip64"
    LIBS="$HIP_LIBS $LIBS"

    # runs compilation test with hipcc
    AC_MSG_CHECKING([hipcc compilation with hipMalloc in -lamdhip64])
    ac_compile='$HIPCC -c $CFLAGS conftest.$ac_ext >&5'
    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([[#include <hip/hip_runtime.h>]],[[void* ptr = 0;hipMalloc(&ptr, 1);]])
    ], [
      AC_MSG_RESULT(yes)
    ], [
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([HIP library function with hipcc compilation failed; try setting HIP_INC.])
    ])

    # runs linking test with hipcc
    AC_MSG_CHECKING([hipcc linking with hipMalloc in -lamdhip64])
    ac_link='$HIPCC -o conftest$ac_exeext $CFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
    AC_LINK_IFELSE([
      AC_LANG_PROGRAM([[
#include <stdio.h>
#include <hip/hip_runtime.h>]],[[void* ptr = 0;hipMalloc(&ptr, 1);]])
    ],[
      AC_MSG_RESULT(yes)
    ],[
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([HIP library linking with HIP failed; try setting HIP_LIB.])
    ])

    #HIP_HEADER_H="hip/hip_runtime.h"
    #AC_CHECK_HEADER([${HIP_HEADER_H}], [], [
    #    AC_MSG_ERROR([HIP header not found; try setting HIP_INC.])
    #])

    CFLAGS="$CFLAGS_save"
    LDFLAGS="$LDFLAGS_save"
    LIBS="$LIBS_save"
    AC_LANG_POP(C)

    # export flags
    AC_SUBST([HIP_CPPFLAGS])
    AC_SUBST([HIP_LDFLAGS])
    AC_SUBST([HIP_LIBS])
])
AM_CONDITIONAL([COND_HIP_PLATFORM_AMD], [test x"$GPU_PLATFORM" = xamd])
AM_CONDITIONAL([COND_HIP_PLATFORM_NVIDIA], [test x"$GPU_PLATFORM" = xnvidia])


###
### OpenMP
###

AS_IF([test x"$want_omp" != xno], [
  AS_BOX([OpenMP])
  AC_MSG_NOTICE([OpenMP compilation is enabled])

  AC_ARG_VAR(OMP_FCFLAGS, [OpenMP Fortran compiler flags])
  AC_ARG_VAR(OMP_LIB, [Location of extra OpenMP libraries])

  # openmp checking
  AC_MSG_NOTICE([OpenMP flag $OMP_FCFLAGS])
  CIT_FC_OPENMP_MODULE([$FC], [$OMP_FCFLAGS])
])

###
### FORCE_VECTORIZATION
###

AS_IF([test x"$want_vec" != xno],[
  AS_BOX([FORCE VECTORIZATION])
  AC_MSG_NOTICE([FORCE_VECTORIZATION is enabled])

  FLAGS_CHECK="${FLAGS_CHECK} ${FC_DEFINE}FORCE_VECTORIZATION"
  CPPFLAGS="${CPPFLAGS} -DFORCE_VECTORIZATION"
])

###
### VTK
###

AS_IF([test x"$want_vtk" != xno], [
    AS_BOX([VTK])
    CIT_OPTIONS_VTK
    CIT_PATH_VTK
])

##
## Git version info
##
AS_BOX([git version info])
CIG_PKG_GIT(SPECFEM3D)

AC_DEFINE_UNQUOTED([SPECFEM3D_VERSION], ["$PACKAGE_VERSION"], [Define SPECFEM3D version])
#AC_DEFINE_UNQUOTED([SPECFEM3D_DOI], ["10.5281/zenodo.xxxxx"], [Define SPECFEM3D doi])

AS_IF([test x"$SPECFEM3D_RELEASE_VERSION" = x0],[
  AC_MSG_NOTICE([building from git repository])
  GIT_PACKAGE_VERSION=$SPECFEM3D_GIT_REVISION
  GIT_COMMIT_VERSION=$SPECFEM3D_GIT_HASH
  GIT_DATE_VERSION=$SPECFEM3D_GIT_DATE
],[
  AC_MSG_NOTICE([not a git repository])
  GIT_PACKAGE_VERSION=$PACKAGE_VERSION
  GIT_COMMIT_VERSION="unknown"
  GIT_DATE_VERSION="unknown"
])
AC_SUBST([GIT_PACKAGE_VERSION])
AC_SUBST([GIT_COMMIT_VERSION])
AC_SUBST([GIT_DATE_VERSION])

############################################################

AS_BOX([setting up default simulation setup])

# Output results.
AC_CONFIG_FILES([
    Makefile
    setup/constants.h
    setup/constants_tomography.h
    setup/precision.h
    setup/config.fh
    setup/version.fh
])

# Scotch makefile
if test "$USE_BUNDLED_SCOTCH" = 1; then
  AC_CONFIG_FILES([
	  ${SCOTCH_DIR}/src/Makefile.inc
  ])
fi

AC_CONFIG_FILES([
    DATA/Par_file:DATA/Par_file
    DATA/CMTSOLUTION:DATA/CMTSOLUTION
    DATA/STATIONS:DATA/STATIONS
])

AC_CONFIG_COMMANDS([bin], [AS_MKDIR_P(bin)])
AC_CONFIG_COMMANDS([obj], [AS_MKDIR_P(obj)])
AC_CONFIG_COMMANDS([DATA], [AS_MKDIR_P(DATA)])
AC_CONFIG_COMMANDS([OUTPUT_FILES], [AS_MKDIR_P(OUTPUT_FILES)])
AC_CONFIG_COMMANDS([DATABASES_MPI], [AS_MKDIR_P(OUTPUT_FILES/DATABASES_MPI)])

AC_OUTPUT

#daniel: adding custom lines to config.h.in through autoheader
AH_BOTTOM([
/* Define to select optimized file i/o for regional simulations */
/* map fails when output files are > 4GB, which is often the case for GPU simulations */
#undef USE_MAP_FUNCTION
])


dnl FYI, this is not a "result", but AS_ECHO is not available in 2.61, and it
dnl basically works like a portable echo.

AC_MSG_RESULT([])
AS_BOX([$PACKAGE_NAME $PACKAGE_VERSION])
AC_MSG_RESULT([])
AC_MSG_RESULT([$0 has completed and set up a default configuration to build.])
AC_MSG_RESULT([])
AC_MSG_RESULT([You may wish to modify the following files before running a simulation:])
AC_MSG_RESULT([AS_HELP_STRING([DATA/Par_file],
                              [Set parameters affecting the simulation.])])
AC_MSG_RESULT([AS_HELP_STRING([DATA/CMTSOLUTION],
                              [Set the source parameters before running the solver.])])
AC_MSG_RESULT([AS_HELP_STRING([DATA/STATIONS],
                              [Set the receiver stations before running the solver.])])
AC_MSG_RESULT([])

dnl end of configure.ac
