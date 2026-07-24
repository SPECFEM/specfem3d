/*
!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!    Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                             CNRS, France
!                      and Princeton University, USA
!                (there are currently many more authors!)
!                          (c) October 2017
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


__global__ void compute_forces_coriolis_device(realw_p accel,
                                               realw_const_p veloc,
                                               realw_const_p displ,
                                               int size,
                                               realw_const_p rmassx,
                                               realw_const_p rmassy,
                                               realw_const_p rmassz,
                                               realw two_omegax,
                                               realw two_omegay,
                                               realw two_omegaz) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  if (id < size) {
    realw omegax = (realw) 0.5 * two_omegax;
    realw omegay = (realw) 0.5 * two_omegay;
    realw omegaz = (realw) 0.5 * two_omegaz;

    // velocity components
    realw vx = veloc[3*id  ];
    realw vy = veloc[3*id+1];
    realw vz = veloc[3*id+2];

    // term 2 \Omega \times \partial_t u (Coriolis force)
    realw f_coriolis_x = two_omegay * vz - two_omegaz * vy;
    realw f_coriolis_y = two_omegaz * vx - two_omegax * vz;
    realw f_coriolis_z = two_omegax * vy - two_omegay * vx;

    // dynamic displacement components
    realw rx = displ[3*id  ];
    realw ry = displ[3*id+1];
    realw rz = displ[3*id+2];

    // first cross product: v1 = Omega x r
    realw facx = omegay * rz - omegaz * ry;
    realw facy = omegaz * rx - omegax * rz;
    realw facz = omegax * ry - omegay * rx;

    // centripetal acceleration: ac = Omega x (Omega x r)
    realw f_centrifugal_x = omegay * facz - omegaz * facy;
    realw f_centrifugal_y = omegaz * facx - omegax * facz;
    realw f_centrifugal_z = omegax * facy - omegay * facx;

    // mass matrix scaling
    realw f_contrib_x = (f_coriolis_x + f_centrifugal_x) / rmassx[id];
    realw f_contrib_y = (f_coriolis_y + f_centrifugal_y) / rmassy[id];
    realw f_contrib_z = (f_coriolis_z + f_centrifugal_z) / rmassz[id];

    // adds contribution to accel (negative sign)
    accel[3*id  ] -= f_contrib_x;
    accel[3*id+1] -= f_contrib_y;
    accel[3*id+2] -= f_contrib_z;
  }
}
