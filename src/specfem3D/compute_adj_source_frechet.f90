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


! compute the integrated derivatives of source parameters (M_jk and X_s)

  subroutine compute_adj_source_frechet(ispec,displ_s,Mxx,Myy,Mzz,Mxy,Mxz,Myz,eps_s,eps_m_s, &
                                        hxir,hetar,hgammar,hpxir,hpetar,hpgammar, &
                                        hprime_xx,hprime_yy,hprime_zz)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  use specfem_par, only: xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,xix_regular, &
                         irregular_element_number

  implicit none

  ! input
  integer,intent(in) :: ispec
  real(kind=CUSTOM_REAL),intent(in) :: displ_s(NDIM,NGLLX,NGLLY,NGLLZ)
  double precision,intent(in) :: Mxx, Myy, Mzz, Mxy, Mxz, Myz
  ! output
  real(kind=CUSTOM_REAL),intent(inout) :: eps_s(NDIM,NDIM), eps_m_s(NDIM)

  ! receiver Lagrange interpolators
  double precision,dimension(NGLLX),intent(in) :: hxir
  double precision,dimension(NGLLY),intent(in) :: hetar
  double precision,dimension(NGLLZ),intent(in) :: hgammar
  double precision,intent(in) :: hpxir(NGLLX),hpetar(NGLLY),hpgammar(NGLLZ)

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY),intent(in) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ),intent(in) :: hprime_zz

  ! local variables
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l, tempy1l,tempy2l,tempy3l, &
             tempz1l,tempz2l,tempz3l, hp1, hp2, hp3, &
             xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl, &
             duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl, &
             xix_s,xiy_s,xiz_s,etax_s,etay_s,etaz_s,gammax_s,gammay_s,gammaz_s, &
             hlagrange_xi, hlagrange_eta, hlagrange_gamma, hlagrange

  real(kind=CUSTOM_REAL) :: eps(NDIM,NDIM), eps_array(NDIM,NDIM,NGLLX,NGLLY,NGLLZ), &
             eps_m_array(NGLLX,NGLLY,NGLLZ)

  integer :: i,j,k,l,ispec_irreg

  ispec_irreg = irregular_element_number(ispec)

  ! first compute the strain at all the GLL points of the source element
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX

        tempx1l = 0._CUSTOM_REAL
        tempx2l = 0._CUSTOM_REAL
        tempx3l = 0._CUSTOM_REAL

        tempy1l = 0._CUSTOM_REAL
        tempy2l = 0._CUSTOM_REAL
        tempy3l = 0._CUSTOM_REAL

        tempz1l = 0._CUSTOM_REAL
        tempz2l = 0._CUSTOM_REAL
        tempz3l = 0._CUSTOM_REAL

        do l = 1,NGLLX    ! assumes NGLLX == NGLLY == NGLLZ
          hp1 = hprime_xx(i,l)
          tempx1l = tempx1l + displ_s(1,l,j,k)*hp1
          tempy1l = tempy1l + displ_s(2,l,j,k)*hp1
          tempz1l = tempz1l + displ_s(3,l,j,k)*hp1

          hp2 = hprime_yy(j,l)
          tempx2l = tempx2l + displ_s(1,i,l,k)*hp2
          tempy2l = tempy2l + displ_s(2,i,l,k)*hp2
          tempz2l = tempz2l + displ_s(3,i,l,k)*hp2

          hp3 = hprime_zz(k,l)
          tempx3l = tempx3l + displ_s(1,i,j,l)*hp3
          tempy3l = tempy3l + displ_s(2,i,j,l)*hp3
          tempz3l = tempz3l + displ_s(3,i,j,l)*hp3
        enddo

        ! derivatives (dudx,..)
        if (ispec_irreg /= 0) then
          !irregular element
          xixl = xixstore(i,j,k,ispec_irreg)
          xiyl = xiystore(i,j,k,ispec_irreg)
          xizl = xizstore(i,j,k,ispec_irreg)
          etaxl = etaxstore(i,j,k,ispec_irreg)
          etayl = etaystore(i,j,k,ispec_irreg)
          etazl = etazstore(i,j,k,ispec_irreg)
          gammaxl = gammaxstore(i,j,k,ispec_irreg)
          gammayl = gammaystore(i,j,k,ispec_irreg)
          gammazl = gammazstore(i,j,k,ispec_irreg)

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l
        else
          !regular element
          duxdxl = xix_regular*tempx1l
          duxdyl = xix_regular*tempx2l
          duxdzl = xix_regular*tempx3l

          duydxl = xix_regular*tempy1l
          duydyl = xix_regular*tempy2l
          duydzl = xix_regular*tempy3l

          duzdxl = xix_regular*tempz1l
          duzdyl = xix_regular*tempz2l
          duzdzl = xix_regular*tempz3l
        endif

        ! symmetric definition of strain: epsilon = 1/2 ( grad(u) + grad(u)^T)

        ! strain eps_jk
        eps(1,1) = duxdxl                     ! eps_xx = duxdx
        eps(1,2) = (duxdyl + duydxl) / 2      ! eps_xy = 1/2 (duxdy + duydz)
        eps(1,3) = (duxdzl + duzdxl) / 2      ! eps_xz = 1/2 (duxdz + duzdx)
        eps(2,2) = duydyl                     ! eps_yy = duydy
        eps(2,3) = (duydzl + duzdyl) / 2      ! eps_yz = 1/2 (duydz + duzdy)
        eps(3,3) = duzdzl                     ! eps_zz = duzdz
        eps(2,1) = eps(1,2)                   ! symmetry: eps_yx = eps_xy
        eps(3,1) = eps(1,3)                   ! symmetry: eps_zx = eps_xz
        eps(3,2) = eps(2,3)                   ! symmetry: eps_zy = eps_yz

        eps_array(:,:,i,j,k) = eps(:,:)

        ! Mjk eps_jk
        eps_m_array(i,j,k) = Mxx * eps(1,1) + Myy * eps(2,2) + Mzz * eps(3,3) + &
                        2 * (Mxy * eps(1,2) + Mxz * eps(1,3) + Myz * eps(2,3))

      enddo
    enddo
  enddo

  ! interpolate the strain eps_s(:,:) from eps_array(:,:,i,j,k)
  eps_s(:,:) = 0.0_CUSTOM_REAL

  xix_s = 0.0_CUSTOM_REAL; xiy_s = 0.0_CUSTOM_REAL; xiz_s = 0.0_CUSTOM_REAL
  etax_s = 0.0_CUSTOM_REAL; etay_s = 0.0_CUSTOM_REAL; etaz_s = 0.0_CUSTOM_REAL
  gammax_s = 0.0_CUSTOM_REAL; gammay_s = 0.0_CUSTOM_REAL; gammaz_s = 0.0_CUSTOM_REAL

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        hlagrange = hxir(i)*hetar(j)*hgammar(k)

        eps_s(1,1) = eps_s(1,1) + eps_array(1,1,i,j,k)*hlagrange
        eps_s(1,2) = eps_s(1,2) + eps_array(1,2,i,j,k)*hlagrange
        eps_s(1,3) = eps_s(1,3) + eps_array(1,3,i,j,k)*hlagrange
        eps_s(2,2) = eps_s(2,2) + eps_array(2,2,i,j,k)*hlagrange
        eps_s(2,3) = eps_s(2,3) + eps_array(2,3,i,j,k)*hlagrange
        eps_s(3,3) = eps_s(3,3) + eps_array(3,3,i,j,k)*hlagrange

        if (ispec_irreg /= 0) then
          !irregular element
          xix_s = xix_s + xixstore(i,j,k,ispec_irreg)*hlagrange
          xiy_s = xiy_s + xiystore(i,j,k,ispec_irreg)*hlagrange
          xiz_s = xiz_s + xizstore(i,j,k,ispec_irreg)*hlagrange
          etax_s = etax_s + etaxstore(i,j,k,ispec_irreg)*hlagrange
          etay_s = etay_s + etaystore(i,j,k,ispec_irreg)*hlagrange
          etaz_s = etaz_s + etazstore(i,j,k,ispec_irreg)*hlagrange
          gammax_s = gammax_s + gammaxstore(i,j,k,ispec_irreg)*hlagrange
          gammay_s = gammay_s + gammaystore(i,j,k,ispec_irreg)*hlagrange
          gammaz_s = gammaz_s + gammazstore(i,j,k,ispec_irreg)*hlagrange
        else
          !regular element
          xix_s = xix_s + xix_regular*hlagrange
          etay_s = xix_s
          gammaz_s = xix_s
        endif

      enddo
    enddo
  enddo

  ! for completion purpose, not used in specfem3D.f90
  eps_s(2,1) = eps_s(1,2)
  eps_s(3,1) = eps_s(1,3)
  eps_s(3,2) = eps_s(2,3)

  ! compute the gradient of M_jk * eps_jk, and then interpolate it
  eps_m_s(:) = 0._CUSTOM_REAL
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        hlagrange_xi    = hpxir(i) * hetar(j)  * hgammar(k)
        hlagrange_eta   = hxir(i)  * hpetar(j) * hgammar(k)
        hlagrange_gamma = hxir(i)  * hetar(j)  * hpgammar(k)

        eps_m_s(1) = eps_m_s(1) + eps_m_array(i,j,k)*(hlagrange_xi * xix_s + hlagrange_eta * etax_s + hlagrange_gamma * gammax_s)
        eps_m_s(2) = eps_m_s(2) + eps_m_array(i,j,k)*(hlagrange_xi * xiy_s + hlagrange_eta * etay_s + hlagrange_gamma * gammay_s)
        eps_m_s(3) = eps_m_s(3) + eps_m_array(i,j,k)*(hlagrange_xi * xiz_s + hlagrange_eta * etaz_s + hlagrange_gamma * gammaz_s)
      enddo
    enddo
  enddo

  end subroutine compute_adj_source_frechet
