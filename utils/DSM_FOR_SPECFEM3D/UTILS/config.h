/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define to alternate name for `main' routine that is called from a `main' in
   the Fortran libraries. */
/* #undef FC_MAIN */

/* Define if emmintrin.h */
#define HAVE_EMMINTRIN 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define if you have POSIX threads libraries and header files. */
/* #undef HAVE_PTHREAD */

/* defined if Scotch is installed */
#define HAVE_SCOTCH 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define if xmmintrin.h */
#define HAVE_XMMINTRIN 1

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "jtromp AT princeton.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "Specfem 3D"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "Specfem 3D 2.1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "Specfem3D"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.1.0"

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
/* #undef YYTEXT_POINTER */


/* Uncomment and define to select optimized file i/o for regional simulations */
/* map fails when output files are > 4GB, which is often the case for GPU simulations */
// #define USE_MAP_FUNCTION


