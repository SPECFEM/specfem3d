/*
 !========================================================================
 !
 !                         S p e c f e m 3 D
 !                         -----------------
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
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 ! The full text of the license is available in file "LICENSE".
 !
 !========================================================================
 */

/* Dimitri Komatitsch, University of Toulouse, May 2011: */

/* added code to force Flush-To-Zero (FTZ) on Intel processors */
/* otherwise Gradual Underflow can be extremely slow. With Intel */
/* ifort one can use compiler option -ftz, but no such option exists */
/* in gcc and gfortran therefore we call an assembler routine directly here. */
/* Very precise description available at */
/* http://software.intel.com/en-us/articles/x87-and-sse-floating-point-assists-in-ia-32-flush-to-zero-ftz-and-denormals-are-zero-daz/ */
/* and at http://www.rz.uni-karlsruhe.de/rz/docs/VTune/reference/vc148.htm */
/* from http://software.intel.com/en-us/articles/x87-and-sse-floating-point-assists-in-ia-32-flush-to-zero-ftz-and-denormals-are-zero-daz : */

/* Flush-To-Zero (FTZ) mode */

/* FTZ is a method of bypassing IEEE 754 methods of dealing with */
/* invalid floating-point numbers due to underflows. This mode is much */
/* faster. Two conditions must be met for FTZ processing to occur: */

/*   * The FTZ bit (bit 15) in the MXCSR register must be masked (value = 1). */
/*   * The underflow exception (bit 11) needs to be masked (value = 1). */

/* This routine is not strictly necessary for SPECFEM, thus if it does not compile on your system
   (since it calls some low-level system routines) just suppress all the lines below (i.e. make it an empty file)
   and comment out the call to force_ftz() in the main SPECFEM program */

#include "config.h"

#define FTZ_BIT 15
#define UNDERFLOW_EXCEPTION_MASK 11

#ifdef __GNUC__ // only use these features on gnu compiler
#ifdef HAVE_XMMINTRIN
  #define FORCE_FTZ
  #include <xmmintrin.h>
#elif HAVE_EMMINTRIN
  #include <emmintrin.h>
  #define FORCE_FTZ
#endif
#endif // __GNUC__

void
FC_FUNC_(force_ftz,FORCE_FTZ)()
{

// DK DK Nov 2018: uncomment this if you have any problem compiling this file
//#undef FORCE_FTZ

#ifdef __GNUC__
#ifdef FORCE_FTZ
  unsigned int x;

  /* force FTZ by setting bits 11 and 15 to one */
  x = _mm_getcsr();
  x |= (1 << FTZ_BIT);
  x |= (1 << UNDERFLOW_EXCEPTION_MASK);
  _mm_setcsr(x);
#endif
#endif // __GNUC__
}
