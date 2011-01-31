dnl                                               -*- Autoconf -*-
dnl
dnl Process this file with autoconf to produce the 'configure' script:
dnl
dnl ACLOCAL="aclocal -I m4" autoreconf
dnl
dnl You must have recent versions of Autoconf and Automake installed.

# $Id: configure.ac,v 1.1 2005/08/24 22:45:51 leif Exp $


AC_PREREQ(2.61)
AC_INIT([Specfem 3D], [2.0.0], [jtromp AT caltech.edu], [Specfem3D])
AC_CONFIG_SRCDIR([README])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])


# 'configure' options

AC_ARG_WITH([pyre],
    [AC_HELP_STRING([--with-pyre],
        [build Pyrized version @<:@default=no@:>@])],
    [want_pyre="$withval"],
    [want_pyre=no])
AM_CONDITIONAL([COND_PYRE], [test "$want_pyre" = yes])

AC_ARG_WITH([mpi],
    [AC_HELP_STRING([--with-mpi],
        [build parallel version @<:@default=yes@:>@])],
    [want_mpi="$withval"],
    [want_mpi=yes])
AM_CONDITIONAL([COND_MPI], [test "$want_mpi" = yes])

AC_ARG_ENABLE([double-precision],
    [AC_HELP_STRING([--enable-double-precision],
        [solver in double precision @<:@default=no@:>@])],
    [want_double_precision="$enableval"],
    [want_double_precision=no])
if test "$want_double_precision" = no; then
    CUSTOM_REAL=SIZE_REAL
    CUSTOM_MPI_TYPE=MPI_REAL
else
    CUSTOM_REAL=SIZE_DOUBLE
    CUSTOM_MPI_TYPE=MPI_DOUBLE_PRECISION
fi
AC_SUBST([CUSTOM_REAL])
AC_SUBST([CUSTOM_MPI_TYPE])


# Checks for programs.

if test "$want_pyre" = yes; then
    AM_PATH_PYTHON([2.3])
    CIT_PYTHON_SYSCONFIG
fi

# a courtesy to the installed base of users
if test x"$FC" = x && test x"$F90" != x; then
    FC="$F90"
fi
if test x"$MPIFC" = x && test x"$MPIF90" != x; then
    MPIFC="$MPIF90"
fi

AC_PROG_FC
export FC
export MPIFC
F77="$FC"
FFLAGS="$FCFLAGS"
AC_PROVIDE([AC_PROG_F77])
AC_SUBST([FCENV])

flags_guess="$SHELL flags.guess"
AC_MSG_NOTICE([running $flags_guess])
flags=`$flags_guess` ||
    AC_MSG_ERROR([$flags_guess failed])
eval $flags

AC_FC_WRAPPERS
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)
AC_FC_FREEFORM()

AC_PROG_CC


#checks for Scotch

export SCOTCH_DIR
export SCOTCH_LIBDIR
export SCOTCH_INCLUDEDIR
export USE_BUNDLED_SCOTCH

AC_ARG_WITH([scotch-dir],
AC_HELP_STRING([--with-scotch-dir],[define the root path to Scotch (e.g. /opt/scotch/)]),
[
 ac_scotch_dir="$withval";
])

AC_ARG_WITH([scotch-includedir],
AC_HELP_STRING([--with-scotch-includedir],[define the path to the Scotch headers (e.g. /opt/scotch/include)]),
[
 ac_scotch_include_dir="$withval";
])

AC_ARG_WITH([scotch-libdir],
AC_HELP_STRING([--with-scotch-libdir],[define the path to the Scotch libraries (e.g. /opt/scotch/lib)]),
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
  if test -n "${ac_scotch_include_dir}"; then
    CPPFLAGS="${CPPFLAGS} -I${ac_scotch_include_dir}";
  fi

dnl This does not work because we are in Fortran mode, which does not
dnl handle headers.
dnl
dnl AC_CHECK_HEADER(scotchf.h,[scotch_include="yes";],[scotch_include="no";CPPFLAGS=${ac_save_cppflags}])

  scotch_usable="${scotch_lib}"

fi

AC_MSG_CHECKING([whether Scotch is usable])
if test "x${scotch_usable}" = "xyes"; then
 AC_DEFINE([HAVE_SCOTCH],[1],[defined if Scotch is installed])
 USE_BUNDLED_SCOTCH=0
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
 USE_BUNDLED_SCOTCH=1
 SCOTCH_LIBDIR="scotch_5.1.10b/lib"
 SCOTCH_INCLUDEDIR="scotch_5.1.10b/include"
fi
LDFLAGS=${ac_save_ldflags}

# influential environment variables

AC_ARG_VAR(USE_BUNDLED_SCOTCH, [Set to 1 to always use the bundled Scotch library])
AC_ARG_VAR(SCOTCH_DIR, [Directory where Scotch is installed])
AC_ARG_VAR(SCOTCH_INCLUDEDIR, [Directory where Scotch headers are installed])
AC_ARG_VAR(SCOTCH_LIBDIR, [Directory where Scotch libraries are installed])
AC_ARG_VAR(MPIFC, [MPI Fortran compiler command])
AC_ARG_VAR(MPILIBS, [extra libraries for linking MPI programs])
AC_ARG_VAR(MPICC, [MPI C compiler command])
AC_ARG_VAR(FLAGS_CHECK, [Fortran compiler flags for non-critical subroutines])
AC_ARG_VAR(FLAGS_NO_CHECK, [Fortran compiler flags for creating fast, production-run code for critical subroutines])
if test x"$MPIFC" = x; then
    MPIFC=mpif90
fi
if test x"$MPICC" = x; then
    MPICC=mpicc
fi

AC_ARG_VAR(LOCAL_PATH_IS_ALSO_GLOBAL, [files on a local path on each node are also seen as global with same path @<:@default=true@:>@])
if test x"$LOCAL_PATH_IS_ALSO_GLOBAL" = x; then
    LOCAL_PATH_IS_ALSO_GLOBAL=true
fi

AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])



# Checks for Python modules and packages.

AC_LANG(C)
if test "$want_pyre" = yes; then
    builddir=`pwd`
    save_PYTHONPATH="$PYTHONPATH"
    PYTHONPATH="$builddir/python:$PYTHONPATH"; export PYTHONPATH
    # Cheap work-around for the inability of merlin/setuptools to control #! line.
    rm -f pyspecfem3D
    $as_ln_s $PYTHON pyspecfem3D
    pyspecfem3D="$builddir/pyspecfem3D"
    cd $srcdir

    AC_MSG_NOTICE([downloading missing Python dependencies])
    AS_IF([AC_TRY_COMMAND([$pyspecfem3D setup.py install_deps -zmxd $builddir/deps >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD])],
          [],
          [AC_MSG_FAILURE([cannot download missing Python dependencies])])

    AC_MSG_NOTICE([building Python dependencies])
    AS_IF([AC_TRY_COMMAND([$pyspecfem3D setup.py develop -H None -f $builddir/deps -s $builddir -d $builddir/python >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD])],
          [],
          [AC_MSG_FAILURE([building Python dependencies])])

    AC_MSG_CHECKING([for egg-related flags])
    AS_IF([AC_TRY_COMMAND([$pyspecfem3D setup.py egg_flags >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD])],
          [AC_MSG_RESULT(ok)
           . egg-flags.sh
           rm -f egg-flags.sh
          ],
          [AC_MSG_RESULT(failed)
          AC_MSG_FAILURE([cannot scan Python eggs for flags])])

    cd $builddir
    rm -f pyspecfem3D
    PYTHONPATH="$save_PYTHONPATH"
    PYTHONPATH="${pythondir}:${pyexecdir}${save_PYTHONPATH:+:${save_PYTHONPATH}}"

    AC_SUBST(PYTHONPATH)
    AC_SUBST(PYTHON_EGG_CFLAGS)
    AC_SUBST(PYTHON_EGG_CPPFLAGS)
    AC_SUBST(PYTHON_EGG_LDFLAGS)
fi


# Checks for libraries.


# Checks for header files.

if test "$want_mpi" = yes; then
    CIT_FC_MPI_HEADER([$MPIFC], [$FLAGS_NO_CHECK])
fi

# Checks for typedefs, structures, and compiler characteristics.

if test "$want_pyre" = yes; then
    CIT_FC_MAIN
fi


# Checks for library functions.


AC_CONFIG_FILES([Makefile in_data_files/Par_file in_data_files/CMTSOLUTION in_data_files/STATIONS src/meshfem3D/Makefile src/meshfem3D/constants.h src/meshfem3D/precision.h  src/decompose_mesh_SCOTCH/Makefile src/specfem3D/Makefile src/generate_databases/Makefile src/shared/constants.h src/shared/precision.h src/decompose_mesh_SCOTCH/scotch_5.1.10b/src/Makefile.inc src/check_mesh_quality_CUBIT_Abaqus/Makefile])
AC_OUTPUT


dnl end of configure.ac
