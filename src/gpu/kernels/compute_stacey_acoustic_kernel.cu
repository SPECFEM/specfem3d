/*
!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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


__global__ void compute_stacey_acoustic_kernel(field* potential_dot_acoustic,
                                               field* potential_dot_dot_acoustic,
                                               int* abs_boundary_ispec,
                                               int* abs_boundary_ijk,
                                               realw* abs_boundary_jacobian2Dw,
                                               int* d_ibool,
                                               realw* rhostore,
                                               realw* kappastore,
                                               int* ispec_is_acoustic,
                                               int SIMULATION_TYPE,
                                               int SAVE_FORWARD,
                                               int num_abs_boundary_faces,
                                               field* b_potential_dot_acoustic,
                                               field* b_potential_dot_dot_acoustic,
                                               field* b_absorb_potential,
                                               int gravity) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,k,iglob,ispec;
  realw rhol,kappal,cpl;
  realw jacobianw;
  field vel,absorbl;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if (iface < num_abs_boundary_faces){

  //  if (igll<NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if (ispec_is_acoustic[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

      iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      // determines bulk sound speed
      rhol = rhostore[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

      kappal = kappastore[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

      cpl = sqrt( kappal / rhol );

      // velocity
      if (gravity ){
        // daniel: TODO - check gravity and stacey condition here...
        // uses a potential definition of: s = grad(chi)
        vel = potential_dot_acoustic[iglob] / rhol ;
      }else{
        // uses a potential definition of: s = 1/rho grad(chi)
        vel = potential_dot_acoustic[iglob] / rhol;
      }

      // gets associated, weighted jacobian
      jacobianw = abs_boundary_jacobian2Dw[INDEX2(NGLL2,igll,iface)];

      // Sommerfeld condition
      absorbl = vel * jacobianw / cpl;

      atomicAdd(&potential_dot_dot_acoustic[iglob],-absorbl);

      // adjoint simulations
      if (SIMULATION_TYPE == 3){
        // Sommerfeld condition
        absorbl = b_absorb_potential[INDEX2(NGLL2,igll,iface)];
        atomicAdd(&b_potential_dot_dot_acoustic[iglob],-absorbl);

      }else if (SIMULATION_TYPE == 1 && SAVE_FORWARD ){
        // saves boundary values
        b_absorb_potential[INDEX2(NGLL2,igll,iface)] = absorbl;
      }
    }
//  }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_stacey_acoustic_single_kernel(field* potential_dot_acoustic,
                                                      field* potential_dot_dot_acoustic,
                                                      int* abs_boundary_ispec,
                                                      int* abs_boundary_ijk,
                                                      realw* abs_boundary_jacobian2Dw,
                                                      int* d_ibool,
                                                      realw* rhostore,
                                                      realw* kappastore,
                                                      int* ispec_is_acoustic,
                                                      int FORWARD_OR_ADJOINT,
                                                      int SIMULATION_TYPE,
                                                      int SAVE_FORWARD,
                                                      int num_abs_boundary_faces,
                                                      field* b_absorb_potential,
                                                      int gravity) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,k,iglob,ispec;
  realw rhol,kappal,cpl;
  realw jacobianw;
  field vel,absorbl;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if (iface < num_abs_boundary_faces){

  //  if (igll<NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if (ispec_is_acoustic[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

      iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      // adjoint simulations
      if (FORWARD_OR_ADJOINT == 3){
        // Sommerfeld condition
        absorbl = b_absorb_potential[INDEX2(NGLL2,igll,iface)];
        atomicAdd(&potential_dot_dot_acoustic[iglob],-absorbl);

      }else{
        // determines bulk sound speed
        rhol = rhostore[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];
        kappal = kappastore[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

        cpl = sqrt( kappal / rhol );

        // velocity
        if (gravity ){
          // daniel: TODO - check gravity and stacey condition here...
          // uses a potential definition of: s = grad(chi)
          vel = potential_dot_acoustic[iglob] / rhol ;
        }else{
          // uses a potential definition of: s = 1/rho grad(chi)
          vel = potential_dot_acoustic[iglob] / rhol;
        }

        // gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw[INDEX2(NGLL2,igll,iface)];

        // Sommerfeld condition
        absorbl = vel * jacobianw / cpl;

        atomicAdd(&potential_dot_dot_acoustic[iglob],-absorbl);

        if (SIMULATION_TYPE == 1 && SAVE_FORWARD){
          // saves boundary values
          b_absorb_potential[INDEX2(NGLL2,igll,iface)] = absorbl;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_stacey_acoustic_undoatt_kernel( field* potential_dot_acoustic,
                                                        field* potential_dot_dot_acoustic,
                                                        int* abs_boundary_ispec,
                                                        int* abs_boundary_ijk,
                                                        realw* abs_boundary_jacobian2Dw,
                                                        int* d_ibool,
                                                        realw* rhostore,
                                                        realw* kappastore,
                                                        int* ispec_is_acoustic,
                                                        int num_abs_boundary_faces,
                                                        int gravity) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,k,iglob,ispec;
  realw rhol,kappal,cpl;
  realw jacobianw;
  field vel,absorbl;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if (iface < num_abs_boundary_faces){

  //  if (igll<NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if (ispec_is_acoustic[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

      iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      // determines bulk sound speed
      rhol = rhostore[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];
      kappal = kappastore[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

      cpl = sqrt( kappal / rhol );

      // velocity
      if (gravity ){
        // daniel: TODO - check gravity and stacey condition here...
        // uses a potential definition of: s = grad(chi)
        vel = potential_dot_acoustic[iglob] / rhol ;
      }else{
        // uses a potential definition of: s = 1/rho grad(chi)
        vel = potential_dot_acoustic[iglob] / rhol;
      }

      // gets associated, weighted jacobian
      jacobianw = abs_boundary_jacobian2Dw[INDEX2(NGLL2,igll,iface)];

      // Sommerfeld condition
      absorbl = vel * jacobianw / cpl;

      atomicAdd(&potential_dot_dot_acoustic[iglob],-absorbl);
    }
  }
}


