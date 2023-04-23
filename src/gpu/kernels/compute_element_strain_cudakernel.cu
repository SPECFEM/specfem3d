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


__global__ void compute_element_strain_cudakernel(int* ispec_is_elastic,
                                                  int* d_ibool,
                                                  realw* displ,
                                                  realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                                  realw* epsilondev_xz,realw* epsilondev_yz,
                                                  realw* epsilondev_trace,
                                                  realw* epsilon_trace_over_3,
                                                  realw* d_xix,realw* d_xiy,realw* d_xiz,
                                                  realw* d_etax,realw* d_etay,realw* d_etaz,
                                                  realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                                  int* d_irregular_element_number,
                                                  realw xix_regular,
                                                  realw* d_hprime_xx,
                                                  int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;
  int ijk_ispec = ijk + NGLL3*ispec;
  //int ijk_ispec_padded = ijk + NGLL3_PADDED*ispec;

  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  float tempx1l,tempx2l,tempx3l;
  float tempy1l,tempy2l,tempy3l;
  float tempz1l,tempz2l,tempz3l;
  float xixl,xiyl,xizl;
  float etaxl,etayl,etazl;
  float gammaxl,gammayl,gammazl;
  float duxdxl,duxdyl,duxdzl;
  float duydxl,duydyl,duydzl;
  float duzdxl,duzdyl,duzdzl;
  float templ,fac1,fac2,fac3;

  // handles case when there is 1 extra block (due to rectangular grid)
  if (ispec < NSPEC_AB) {

    // elastic elements only
    if (ispec_is_elastic[ispec]) {
      int iglob = d_ibool[ijk + NGLL3_PADDED*ispec] - 1;

      if (ijk < NGLL3){
        sh_tempx[ijk] = displ[iglob*3];
        sh_tempy[ijk] = displ[iglob*3+1];
        sh_tempz[ijk] = displ[iglob*3+2];
      }
      // synchronizes threads
      __syncthreads();

      int K = ijk/NGLL2;
      int J = (ijk - K*NGLL2)/NGLLX;
      int I = ijk - K*NGLL2 - J*NGLLX;

      tempx1l = 0.0f;
      tempx2l = 0.0f;
      tempx3l = 0.0f;
      tempy1l = 0.0f;
      tempy2l = 0.0f;
      tempy3l = 0.0f;
      tempz1l = 0.0f;
      tempz2l = 0.0f;
      tempz3l = 0.0f;

      for (int l = 0; l <= NGLLX - 1; l += 1) {
        fac1 = d_hprime_xx[l * NGLLX + I];
        tempx1l = tempx1l + sh_tempx[K * NGLL2 + J * NGLLX + l] * fac1;
        tempy1l = tempy1l + sh_tempy[K * NGLL2 + J * NGLLX + l] * fac1;
        tempz1l = tempz1l + sh_tempz[K * NGLL2 + J * NGLLX + l] * fac1;
        fac2 = d_hprime_xx[l * NGLLX + J];
        tempx2l = tempx2l + sh_tempx[K * NGLL2 + l * NGLLX + I] * fac2;
        tempy2l = tempy2l + sh_tempy[K * NGLL2 + l * NGLLX + I] * fac2;
        tempz2l = tempz2l + sh_tempz[K * NGLL2 + l * NGLLX + I] * fac2;
        fac3 = d_hprime_xx[l * NGLLX + K];
        tempx3l = tempx3l + sh_tempx[l * NGLL2 + J * NGLLX + I] * fac3;
        tempy3l = tempy3l + sh_tempy[l * NGLL2 + J * NGLLX + I] * fac3;
        tempz3l = tempz3l + sh_tempz[l * NGLL2 + J * NGLLX + I] * fac3;
      }

      int ispec_irreg = d_irregular_element_number[ispec] - 1;
      if (ispec_irreg >= 0){
        int offset = ispec_irreg*NGLL3_PADDED + ijk;
        xixl = d_xix[offset];
        xiyl = d_xiy[offset];
        xizl = d_xiz[offset];
        etaxl = d_etax[offset];
        etayl = d_etay[offset];
        etazl = d_etaz[offset];
        gammaxl = d_gammax[offset];
        gammayl = d_gammay[offset];
        gammazl = d_gammaz[offset];

        // derivatives
        duxdxl = xixl * tempx1l + etaxl * tempx2l + gammaxl * tempx3l;
        duxdyl = xiyl * tempx1l + etayl * tempx2l + gammayl * tempx3l;
        duxdzl = xizl * tempx1l + etazl * tempx2l + gammazl * tempx3l;

        duydxl = xixl * tempy1l + etaxl * tempy2l + gammaxl * tempy3l;
        duydyl = xiyl * tempy1l + etayl * tempy2l + gammayl * tempy3l;
        duydzl = xizl * tempy1l + etazl * tempy2l + gammazl * tempy3l;

        duzdxl = xixl * tempz1l + etaxl * tempz2l + gammaxl * tempz3l;
        duzdyl = xiyl * tempz1l + etayl * tempz2l + gammayl * tempz3l;
        duzdzl = xizl * tempz1l + etazl * tempz2l + gammazl * tempz3l;      }
      else{
        // derivatives
        duxdxl = xix_regular*tempx1l;
        duxdyl = xix_regular*tempx2l;
        duxdzl = xix_regular*tempx3l;

        duydxl = xix_regular*tempy1l;
        duydyl = xix_regular*tempy2l;
        duydzl = xix_regular*tempy3l;

        duzdxl = xix_regular*tempz1l;
        duzdyl = xix_regular*tempz2l;
        duzdzl = xix_regular*tempz3l;
      }

      // stores strains
      epsilondev_trace[ijk_ispec] = (duxdxl + duydyl + duzdzl);

      templ = (duxdxl + duydyl + duzdzl) * (0.3333333333333333f);
      epsilon_trace_over_3[ijk_ispec] = templ;

      epsilondev_xx[ijk_ispec] = duxdxl - templ;
      epsilondev_yy[ijk_ispec] = duydyl - templ;
      epsilondev_xy[ijk_ispec] = (duxdyl + duydxl) * (0.5f);
      epsilondev_xz[ijk_ispec] = (duzdxl + duxdzl) * (0.5f);
      epsilondev_yz[ijk_ispec] = (duzdyl + duydzl) * (0.5f);
    }
  }
}


