/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

#include <stdio.h>
#include <stdlib.h>

#include "config.h"

// for ASDF reader setup

void FC_FUNC_(asdf_setup,ASDF_SETUP)(void) {
  fprintf(stderr,"ERROR: ASDF_FORMAT enabled without ASDF Support. To enable support, reconfigure with --with-asdf flag.\n");
  exit(1);
}
void FC_FUNC_(asdf_cleanup,ASDF_CLEANUP)(void) {}

// for ASDF writer

void FC_FUNC_(init_asdf_data,INIT_ASDF_DATA)(void) {
  fprintf(stderr,"ERROR: ASDF_FORMAT enabled without ASDF Support. To enable support, reconfigure with --with-asdf flag.\n");
  exit(1);
}
void FC_FUNC_(store_asdf_data,STORE_ASDF_DATA)(void) {}
void FC_FUNC_(close_asdf_data,CLOSE_ASDF_DATA)(void) {}
void FC_FUNC_(write_asdf,WRITE_ASDF)(void) {}

// for ASDF reader

void FC_FUNC_(read_adjoint_sources_asdf,READ_ADJOINT_SOURCES_ASDF)(void) {}
void FC_FUNC_(check_adjoint_sources_asdf,CHECK_ADJOINT_SOURCES_ASDF)(void) {}
