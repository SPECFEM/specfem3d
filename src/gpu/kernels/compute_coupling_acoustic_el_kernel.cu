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


__global__ void compute_coupling_acoustic_el_kernel(realw* displ,
                                                    field* potential_dot_dot_acoustic,
                                                    int num_coupling_ac_el_faces,
                                                    int* coupling_ac_el_ispec,
                                                    int* coupling_ac_el_ijk,
                                                    realw* coupling_ac_el_normal,
                                                    realw* coupling_ac_el_jacobian2Dw,
                                                    int* d_ibool,
                                                    int simulation_type,
                                                    int backward_simulation) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,k,iglob,ispec;
  realw displ_x,displ_y,displ_z,displ_n;
  realw nx,ny,nz;
  realw jacobianw;

  if (iface < num_coupling_ac_el_faces){

    // don't compute points outside NGLLSQUARE==NGLL2==25
    // way 2: no further check needed since blocksize = 25
    //  if (igll<NGLL2) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = coupling_ac_el_ispec[iface] - 1;

    i = coupling_ac_el_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)] - 1;
    j = coupling_ac_el_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)] - 1;
    k = coupling_ac_el_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)] - 1;

    iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

    // elastic displacement on global point
    displ_x = displ[iglob*3] ; // (1,iglob)
    displ_y = displ[iglob*3+1] ; // (2,iglob)
    displ_z = displ[iglob*3+2] ; // (3,iglob)

    // adjoint wavefield case
    if (simulation_type != 1 && backward_simulation == 0){
      // handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
      // adjoint definition: \partial_t^2 \bfs^\dagger = - \frac{1}{\rho} \bfnabla\phi^\dagger
      displ_x = - displ_x;
      displ_y = - displ_y;
      displ_z = - displ_z;
    }

    // gets associated normal on GLL point
    nx = coupling_ac_el_normal[INDEX3(NDIM,NGLL2,0,igll,iface)]; // (1,igll,iface)
    ny = coupling_ac_el_normal[INDEX3(NDIM,NGLL2,1,igll,iface)]; // (2,igll,iface)
    nz = coupling_ac_el_normal[INDEX3(NDIM,NGLL2,2,igll,iface)]; // (3,igll,iface)

    // calculates displacement component along normal
    // (normal points outwards of acoustic element)
    displ_n = displ_x*nx + displ_y*ny + displ_z*nz;

    // gets associated, weighted jacobian
    jacobianw = coupling_ac_el_jacobian2Dw[INDEX2(NGLL2,igll,iface)];

    // continuity of pressure and normal displacement on global point

    // note: Newmark time scheme together with definition of scalar potential:
    //          pressure = - chi_dot_dot
    //          requires that this coupling term uses the updated displacement at time step [t+delta_t],
    //          which is done at the very beginning of the time loop
    //          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
    //          it also means you have to calculate and update this here first before
    //          calculating the coupling on the elastic side for the acceleration...
    atomicAdd(&potential_dot_dot_acoustic[iglob],+ jacobianw*displ_n);

  //  }
  }
}

