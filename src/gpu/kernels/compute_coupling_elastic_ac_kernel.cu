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


__global__ void compute_coupling_elastic_ac_kernel(field* potential_dot_dot_acoustic,
                                                    realw* accel,
                                                    int num_coupling_ac_el_faces,
                                                    int* coupling_ac_el_ispec,
                                                    int* coupling_ac_el_ijk,
                                                    realw* coupling_ac_el_normal,
                                                    realw* coupling_ac_el_jacobian2Dw,
                                                    int* d_ibool,
                                                    int gravity,
                                                    realw* minus_g,
                                                    realw* rhostore,
                                                    realw* displ,
                                                    int simulation_type,
                                                    int backward_simulation) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,k,iglob,ispec;
  field pressure;
  realw nx,ny,nz;
  realw jacobianw;
  realw rhol;

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

    // gets associated normal on GLL point
    // note: normal points away from acoustic element
    nx = coupling_ac_el_normal[INDEX3(NDIM,NGLL2,0,igll,iface)]; // (1,igll,iface)
    ny = coupling_ac_el_normal[INDEX3(NDIM,NGLL2,1,igll,iface)]; // (2,igll,iface)
    nz = coupling_ac_el_normal[INDEX3(NDIM,NGLL2,2,igll,iface)]; // (3,igll,iface)

    // gets associated, weighted jacobian
    jacobianw = coupling_ac_el_jacobian2Dw[INDEX2(NGLL2,igll,iface)];

    // acoustic pressure on global point
    if (gravity ){
      // takes density (from acoustic? element)
      rhol = rhostore[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

      // note: uses potential chi such that displacement s = grad(chi),
      //         pressure becomes: p = - kappa ( div( s ) ) = rho ( - dot_dot_chi + g * s )
      //  g only acting in negative z-direction

      // daniel: TODO - check gravity and coupling would be displ * nz  correct?
      pressure = rhol*( - potential_dot_dot_acoustic[iglob]
                       + minus_g[iglob] * displ[iglob*3+2] );

      //daniel: TODO - check gravity and coupling
      //pressure = - potential_dot_dot_acoustic[iglob] ;
      //if (iface == 128 && igll == 5){
      //  printf("coupling acoustic: %f %f \n",potential_dot_dot_acoustic[iglob],
      //             minus_g[iglob] * displ[iglob*3+2]);
      //}

    }else{
      // no gravity: uses potential chi such that displacement s = 1/rho grad(chi)
      //                  pressure p = - kappa ( div( s )) then becomes: p = - dot_dot_chi
      //                  ( multiplied with factor 1/kappa due to setup of equation of motion )
      pressure = - potential_dot_dot_acoustic[iglob];
    }

    if (simulation_type != 1 && backward_simulation == 0){
      // handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
      // adjoint definition: pressure^\dagger = potential^\dagger
      pressure = - pressure;
    }

    // continuity of displacement and pressure on global point
    //
    // note: Newmark time scheme together with definition of scalar potential:
    //          pressure = - chi_dot_dot
    //          requires that this coupling term uses the *UPDATED* pressure (chi_dot_dot), i.e.
    //          pressure at time step [t + delta_t]
    //          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
    //          it means you have to calculate and update the acoustic pressure first before
    //          calculating this term...
    atomicAdd(&accel[iglob*3],+ jacobianw*nx*pressure);
    atomicAdd(&accel[iglob*3+1],+ jacobianw*ny*pressure);
    atomicAdd(&accel[iglob*3+2],+ jacobianw*nz*pressure);
    //  }
  }
}


