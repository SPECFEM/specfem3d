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

#ifndef SMOOTH_CUDA_H
#define SMOOTH_CUDA_H

#include "mesh_constants_gpu.h"


// smoothing data structure
typedef struct Smooth_data_ {
  realw * x_me;
  realw * y_me;
  realw * z_me;

  realw * data_smooth;
  realw * normalisation;

  realw sigma_h2_inv;
  realw sigma_v2_inv;

  realw h_criterion;
  realw v_criterion;

  int nspec_me;
  int nker;
} Smooth_data;


#endif  // SMOOTH_CUDA_H
