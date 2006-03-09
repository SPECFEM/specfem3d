dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.
dnl

# $Id: configure.ac,v 1.1 2005/08/24 22:45:51 leif Exp $

AC_PREREQ(2.59)
AC_INIT([Specfem 3D Basin], [1.3.0], [jtromp@caltech.edu], [Specfem3DBasin])
AC_CONFIG_SRCDIR([specfem3D.f90])
AC_CONFIG_HEADER([obj/config.h:config.h.in])
AC_CONFIG_MACRO_DIR([m4])

# 'configure' options
AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])
AC_ARG_WITH([pyre],
    [AC_HELP_STRING([--with-pyre],
        [build Pyre modules @<:@default=yes@:>@])],
    [want_pyre="$withval"],
    [want_pyre=yes])
AM_CONDITIONAL([COND_PYRE], [test "$want_pyre" = yes])

# Checks for programs.
if test "$want_pyre" = yes; then
    AM_PATH_PYTHON([2.3])
    CIT_PYTHON_SYSCONFIG
fi
AC_PROG_FC
F77=$FC
FFLAGS=$FCFLAGS
AC_PROVIDE([AC_PROG_F77])
AC_FC_WRAPPERS
CIT_FC_MAIN
AC_PROG_CC

# Checks for Python modules and packages.
if test "$want_pyre" = yes; then
    PYTHONPATH="${pythondir}${PYTHONPATH:+:${PYTHONPATH}}"; export PYTHONPATH
    CIT_PATH_PYTHIA([0.8], [mpi])
    AC_SUBST([PYTHONPATH])
fi

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([obj/config:config.in])
AC_OUTPUT

chmod +x obj/config

dnl end of configure.ac
