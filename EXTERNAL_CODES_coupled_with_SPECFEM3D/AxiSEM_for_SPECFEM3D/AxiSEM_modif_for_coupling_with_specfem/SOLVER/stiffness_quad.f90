!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
module stiffness_quad

  use global_parameters, only: realkind
  use data_matr
  use data_mesh, only: axis_solid, axis_fluid, nsize
  use data_spec
  use data_source

  use unrolled_loops

  implicit none

  public :: glob_stiffness_quad
  public :: glob_anel_stiffness_quad

  private

contains

!-----------------------------------------------------------------------------------------
!> Wrapper routine to avoid if statements in the timeloop
pure subroutine glob_stiffness_quad(glob_stiffness,u)
  use data_mesh, only: npol, nel_solid

  real(kind=realkind), intent(in)  :: u(0:,0:,:,:)
  real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,3)

  if (npol == 4) then
     call glob_stiffness_quad_4(glob_stiffness, u)
  else
     call glob_stiffness_quad_generic(glob_stiffness, u)
  endif

end subroutine glob_stiffness_quad
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_stiffness_quad_generic(glob_stiffness,u)

  use data_mesh, only: npol, nel_solid

  ! I/O for global arrays
  real(kind=realkind), intent(in)  :: u(0:,0:,:,:)
  real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,1:3)

  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_phi
  real(kind=realkind), dimension(0:npol,0:npol) :: us,uz,uphi

  real(kind=realkind), dimension(0:npol,0:npol) :: m_1l, m_2l, m_3l, m_4l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_5l, m_6l, m_7l, m_8l

  real(kind=realkind), dimension(0:npol,0:npol) :: m_w1l, m_w2l, m_w3l, m_w4l, m_w5l

  real(kind=realkind), dimension(0:npol,0:npol) :: m11sl, m21sl, m41sl, m12sl, m22sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m32sl, m42sl, m11zl, m21zl, m41zl
  real(kind=realkind), dimension(0:npol,0:npol) :: m1phil, m2phil, m4phil

  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w1l, m0_w2l, m0_w3l
  real(kind=realkind), dimension(0:npol) :: m0_w4l, m0_w5l, m0_w6l

  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4, X5, X6 ! MxM arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: S1s, S2s, S1phi, S2phi, S1z, S2z ! Sum

  real(kind=realkind), dimension(0:npol,0:npol) :: c1, c2, c3, c4, c5, c6

  real(kind=realkind), dimension(0:npol) :: V1, V2, V3

  integer :: ielem

  do ielem = 1, nel_solid

     us(0:npol,0:npol)   = u(0:npol,0:npol,ielem,1)
     uphi(0:npol,0:npol) = u(0:npol,0:npol,ielem,2)
     uz(0:npol,0:npol)   = u(0:npol,0:npol,ielem,3)

     m_1l(0:npol,0:npol) = M_1(:,:,ielem)
     m_2l(0:npol,0:npol) = M_2(:,:,ielem)
     m_3l(0:npol,0:npol) = M_3(:,:,ielem)
     m_4l(0:npol,0:npol) = M_4(:,:,ielem)
     m_5l(0:npol,0:npol) = M_5(:,:,ielem)
     m_6l(0:npol,0:npol) = M_6(:,:,ielem)
     m_7l(0:npol,0:npol) = M_7(:,:,ielem)
     m_8l(0:npol,0:npol) = M_8(:,:,ielem)

     m_w1l(0:npol,0:npol) = M_w1(:,:,ielem)
     m_w2l(0:npol,0:npol) = M_w2(:,:,ielem)
     m_w3l(0:npol,0:npol) = M_w3(:,:,ielem)
     m_w4l(0:npol,0:npol) = M_w4(:,:,ielem)
     m_w5l(0:npol,0:npol) = M_w5(:,:,ielem)

     m11sl(0:npol,0:npol) = M11s(:,:,ielem)
     m21sl(0:npol,0:npol) = M21s(:,:,ielem)
     m41sl(0:npol,0:npol) = M41s(:,:,ielem)
     m12sl(0:npol,0:npol) = M12s(:,:,ielem)
     m22sl(0:npol,0:npol) = M22s(:,:,ielem)
     m32sl(0:npol,0:npol) = M32s(:,:,ielem)
     m42sl(0:npol,0:npol) = M42s(:,:,ielem)
     m11zl(0:npol,0:npol) = M11z(:,:,ielem)
     m21zl(0:npol,0:npol) = M21z(:,:,ielem)
     m41zl(0:npol,0:npol) = M41z(:,:,ielem)

     m1phil(0:npol,0:npol) = M1phi(:,:,ielem)
     m2phil(0:npol,0:npol) = M2phi(:,:,ielem)
     m4phil(0:npol,0:npol) = M4phi(:,:,ielem)

     ! First MxM
     if (.not. axis_solid(ielem) ) then
        call mxm(G2T, us, X1)
        call mxm(G2T, uphi, X2)
        call mxm(G2T, uz, X3)
     else
        call mxm(G1T, us, X1)
        call mxm(G1T, uphi, X2)
        call mxm(G1T, uz, X3)
     endif

     call mxm(us, G2, X4)
     call mxm(uphi, G2, X5)
     call mxm(uz, G2, X6)

     ! s and phi components
     ! buffering terms that occure in both components
     c1 = m_2l * X4
     c2 = m_1l * X1
     c3 = m_6l * X5
     c4 = m_5l * X2
     c5 = m_4l * X6
     c6 = m_3l * X3

     loc_stiffness_s = c1 + c2 + 2 * (c3 + c4) + c5 + c6 &
                        + m_w1l * us + m_w2l * uphi + 2 * m_w3l * uz
     loc_stiffness_phi = -2 * (c1 + c2 + c5 + c6) - (c3 + c4) &
                            + m_w2l * us +  m_w4l * uphi - m_w3l * uz

     ! z component
     loc_stiffness_z = 2 * (m_8l * X5 + m_7l * X2) + m_w3l * (2 * us - uphi) + m_w5l * uz

     ! s component
     S1s = m11sl * X4 + m21sl * X1 + m12sl * X6 + m22sl * X3 + m_1l * (us - 2 * uphi)
     S2s = m11sl * X1 + m41sl * X4 + m32sl * X3 + m42sl * X6 + m_2l * (us - 2 * uphi)

     ! z component
     S1z = m11zl * X6 + m21zl * X3 + m32sl * X4 + m22sl * X1 + m_3l * (us - 2 * uphi)
     S2z = m11zl * X3 + m41zl * X6 + m12sl * X1 + m42sl * X4 + m_4l * (us - 2 * uphi)

     ! phi component
     S1phi = m1phil * X5 + m2phil * X2 + m_5l * (2 * us - uphi) + 2 * m_7l * uz
     S2phi = m1phil * X2 + m4phil * X5 + m_6l * (2 * us - uphi) + 2 * m_8l * uz

     !Second MxM
     call mxm(S2s, G2T, X2)
     call mxm(S2phi, G2T, X4)
     call mxm(S2z, G2T, X6)

     if (.not. axis_solid(ielem) ) then
        call mxm(G2, S1s, X1)
        call mxm(G2, S1phi, X3)
        call mxm(G2, S1z, X5)
     else
        call mxm(G1, S1s, X1)
        call mxm(G1, S1phi, X3)
        call mxm(G1, S1z, X5)
     endif

     ! Final Sum
     loc_stiffness_s   = loc_stiffness_s   + X1 + X2
     loc_stiffness_phi = loc_stiffness_phi + X3 + X4
     loc_stiffness_z   = loc_stiffness_z   + X5 + X6

     if ( axis_solid(ielem) ) then

        m0_w1l(0:npol) = M0_w1(0:npol,ielem)
        m0_w2l(0:npol) = M0_w2(0:npol,ielem)
        m0_w4l(0:npol) = M0_w4(0:npol,ielem)
        m0_w6l(0:npol) = M0_w6(0:npol,ielem)

        ! additional anisotropic terms
        m0_w3l(0:npol) = M0_w3(0:npol,ielem)
        m0_w5l(0:npol) = M0_w5(0:npol,ielem)

        ! VxM
        call vxm(G0, us, V1)
        call vxm(G0, uphi, V2)
        call vxm(G0, uz, V3)

        ! Collocations, Sums, Tensorization
        S1s = outerprod(G0, m0_w1l * V1 + m0_w2l * V2 + m0_w3l * V3)

        S1phi = outerprod(G0, m0_w2l * V1 + m0_w4l * V2 + m0_w5l * V3)

        S1z = outerprod(G0, m0_w3l * V1 + m0_w5l * V2 + m0_w6l * V3)
        ! end additional anisotropic terms

        ! Final Sum
        loc_stiffness_s   = loc_stiffness_s   + S1s
        loc_stiffness_phi = loc_stiffness_phi + S1phi
        loc_stiffness_z   = loc_stiffness_z   + S1z

     endif

     glob_stiffness(0:npol,0:npol,ielem,1) = loc_stiffness_s
     glob_stiffness(0:npol,0:npol,ielem,2) = loc_stiffness_phi
     glob_stiffness(0:npol,0:npol,ielem,3) = loc_stiffness_z

  enddo

end subroutine glob_stiffness_quad_generic
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_stiffness_quad_4(glob_stiffness,u)

  use data_mesh, only: nel_solid

  integer, parameter               :: npol = 4
  ! I/O for global arrays
  real(kind=realkind), intent(in)  :: u(0:npol,0:npol,nel_solid,3)
  real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,1:3)

  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_phi
  real(kind=realkind), dimension(0:npol,0:npol) :: us,uz,uphi

  real(kind=realkind), dimension(0:npol,0:npol) :: m_1l, m_2l, m_3l, m_4l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_5l, m_6l, m_7l, m_8l

  real(kind=realkind), dimension(0:npol,0:npol) :: m_w1l, m_w2l, m_w3l, m_w4l, m_w5l

  real(kind=realkind), dimension(0:npol,0:npol) :: m11sl, m21sl, m41sl, m12sl, m22sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m32sl, m42sl, m11zl, m21zl, m41zl
  real(kind=realkind), dimension(0:npol,0:npol) :: m1phil, m2phil, m4phil

  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w1l, m0_w2l, m0_w3l
  real(kind=realkind), dimension(0:npol) :: m0_w4l, m0_w5l, m0_w6l

  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4, X5, X6 ! MxM arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: S1s, S2s, S1phi, S2phi, S1z, S2z ! Sum

  real(kind=realkind), dimension(0:npol,0:npol) :: c1, c2, c3, c4, c5, c6

  real(kind=realkind), dimension(0:npol) :: V1, V2, V3

  integer :: ielem

  do ielem = 1, nel_solid

     us(0:npol,0:npol)   = u(0:npol,0:npol,ielem,1)
     uphi(0:npol,0:npol) = u(0:npol,0:npol,ielem,2)
     uz(0:npol,0:npol)   = u(0:npol,0:npol,ielem,3)

     m_1l(0:npol,0:npol) = M_1(:,:,ielem)
     m_2l(0:npol,0:npol) = M_2(:,:,ielem)
     m_3l(0:npol,0:npol) = M_3(:,:,ielem)
     m_4l(0:npol,0:npol) = M_4(:,:,ielem)
     m_5l(0:npol,0:npol) = M_5(:,:,ielem)
     m_6l(0:npol,0:npol) = M_6(:,:,ielem)
     m_7l(0:npol,0:npol) = M_7(:,:,ielem)
     m_8l(0:npol,0:npol) = M_8(:,:,ielem)

     m_w1l(0:npol,0:npol) = M_w1(:,:,ielem)
     m_w2l(0:npol,0:npol) = M_w2(:,:,ielem)
     m_w3l(0:npol,0:npol) = M_w3(:,:,ielem)
     m_w4l(0:npol,0:npol) = M_w4(:,:,ielem)
     m_w5l(0:npol,0:npol) = M_w5(:,:,ielem)

     m11sl(0:npol,0:npol) = M11s(:,:,ielem)
     m21sl(0:npol,0:npol) = M21s(:,:,ielem)
     m41sl(0:npol,0:npol) = M41s(:,:,ielem)
     m12sl(0:npol,0:npol) = M12s(:,:,ielem)
     m22sl(0:npol,0:npol) = M22s(:,:,ielem)
     m32sl(0:npol,0:npol) = M32s(:,:,ielem)
     m42sl(0:npol,0:npol) = M42s(:,:,ielem)
     m11zl(0:npol,0:npol) = M11z(:,:,ielem)
     m21zl(0:npol,0:npol) = M21z(:,:,ielem)
     m41zl(0:npol,0:npol) = M41z(:,:,ielem)

     m1phil(0:npol,0:npol) = M1phi(:,:,ielem)
     m2phil(0:npol,0:npol) = M2phi(:,:,ielem)
     m4phil(0:npol,0:npol) = M4phi(:,:,ielem)

     ! First MxM
     if (.not. axis_solid(ielem) ) then
        call mxm_4(G2T, us, X1)
        call mxm_4(G2T, uphi, X2)
        call mxm_4(G2T, uz, X3)
     else
        call mxm_4(G1T, us, X1)
        call mxm_4(G1T, uphi, X2)
        call mxm_4(G1T, uz, X3)
     endif

     call mxm_4(us, G2, X4)
     call mxm_4(uphi, G2, X5)
     call mxm_4(uz, G2, X6)

     ! s and phi components
     ! buffering terms that occure in both components
     c1 = m_2l * X4
     c2 = m_1l * X1
     c3 = m_6l * X5
     c4 = m_5l * X2
     c5 = m_4l * X6
     c6 = m_3l * X3

     loc_stiffness_s = c1 + c2 + 2 * (c3 + c4) + c5 + c6 &
                        + m_w1l * us + m_w2l * uphi + 2 * m_w3l * uz
     loc_stiffness_phi = -2 * (c1 + c2 + c5 + c6) - (c3 + c4) &
                            + m_w2l * us +  m_w4l * uphi - m_w3l * uz

     ! z component
     loc_stiffness_z = 2 * (m_8l * X5 + m_7l * X2) + m_w3l * (2 * us - uphi) + m_w5l * uz

     ! s component
     S1s = m11sl * X4 + m21sl * X1 + m12sl * X6 + m22sl * X3 + m_1l * (us - 2 * uphi)
     S2s = m11sl * X1 + m41sl * X4 + m32sl * X3 + m42sl * X6 + m_2l * (us - 2 * uphi)

     ! z component
     S1z = m11zl * X6 + m21zl * X3 + m32sl * X4 + m22sl * X1 + m_3l * (us - 2 * uphi)
     S2z = m11zl * X3 + m41zl * X6 + m12sl * X1 + m42sl * X4 + m_4l * (us - 2 * uphi)

     ! phi component
     S1phi = m1phil * X5 + m2phil * X2 + m_5l * (2 * us - uphi) + 2 * m_7l * uz
     S2phi = m1phil * X2 + m4phil * X5 + m_6l * (2 * us - uphi) + 2 * m_8l * uz

     !Second MxM
     call mxm_4(S2s, G2T, X2)
     call mxm_4(S2phi, G2T, X4)
     call mxm_4(S2z, G2T, X6)

     if (.not. axis_solid(ielem) ) then
        call mxm_4(G2, S1s, X1)
        call mxm_4(G2, S1phi, X3)
        call mxm_4(G2, S1z, X5)
     else
        call mxm_4(G1, S1s, X1)
        call mxm_4(G1, S1phi, X3)
        call mxm_4(G1, S1z, X5)
     endif

     ! Final Sum
     loc_stiffness_s   = loc_stiffness_s   + X1 + X2
     loc_stiffness_phi = loc_stiffness_phi + X3 + X4
     loc_stiffness_z   = loc_stiffness_z   + X5 + X6

     if ( axis_solid(ielem) ) then

        m0_w1l(0:npol) = M0_w1(0:npol,ielem)
        m0_w2l(0:npol) = M0_w2(0:npol,ielem)
        m0_w4l(0:npol) = M0_w4(0:npol,ielem)
        m0_w6l(0:npol) = M0_w6(0:npol,ielem)

        ! additional anisotropic terms
        m0_w3l(0:npol) = M0_w3(0:npol,ielem)
        m0_w5l(0:npol) = M0_w5(0:npol,ielem)

        ! VxM
        call vxm_4(G0, us, V1)
        call vxm_4(G0, uphi, V2)
        call vxm_4(G0, uz, V3)

        ! Collocations, Sums, Tensorization
        S1s = outerprod_4(G0, m0_w1l * V1 + m0_w2l * V2 + m0_w3l * V3)

        S1phi = outerprod_4(G0, m0_w2l * V1 + m0_w4l * V2 + m0_w5l * V3)

        S1z = outerprod_4(G0, m0_w3l * V1 + m0_w5l * V2 + m0_w6l * V3)
        ! end additional anisotropic terms

        ! Final Sum
        loc_stiffness_s   = loc_stiffness_s   + S1s
        loc_stiffness_phi = loc_stiffness_phi + S1phi
        loc_stiffness_z   = loc_stiffness_z   + S1z

     endif

     glob_stiffness(0:npol,0:npol,ielem,1) = loc_stiffness_s
     glob_stiffness(0:npol,0:npol,ielem,2) = loc_stiffness_phi
     glob_stiffness(0:npol,0:npol,ielem,3) = loc_stiffness_z

  enddo

end subroutine glob_stiffness_quad_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_quad(glob_stiffness, R, R_cg, cg)
  use data_mesh, only: npol

  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: R(0:,0:,:,:,:)
  real(kind=realkind), intent(in)    :: R_cg(:,:,:,:)
  logical, intent(in)                :: cg

  if (cg) then
     call glob_anel_stiffness_quad_cg4(glob_stiffness, R_cg)
  else
     if (npol == 4) then
        call glob_anel_stiffness_quad_4(glob_stiffness, R)
     else
        call glob_anel_stiffness_quad_generic(glob_stiffness, R)
     endif
  endif

end subroutine glob_anel_stiffness_quad
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_quad_generic(glob_stiffness, R)

  use attenuation, only: n_sls_attenuation
  use data_mesh, only: npol, nel_solid

  ! I/O global arrays
  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:) !0:npol,0:npol,nel_solid,1:3)
  real(kind=realkind), intent(in)    :: R(0:,0:,:,:,:) !(0:npol,0:npol,6,n_sls_attenuation,nel_solid)

  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_p
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z

  real(kind=realkind), dimension(0:npol,0:npol) :: r1, r2, r3, r4, r5, r6

  real(kind=realkind), dimension(0:npol,0:npol) :: yl
  real(kind=realkind), dimension(0:npol,0:npol) :: v_s_etal, v_s_xil
  real(kind=realkind), dimension(0:npol,0:npol) :: v_z_etal, v_z_xil

  real(kind=realkind), dimension(0:npol,0:npol) :: S1s, S2s
  real(kind=realkind), dimension(0:npol,0:npol) :: S1p, S2p
  real(kind=realkind), dimension(0:npol,0:npol) :: S1z, S2z
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4, X5, X6

  real(kind=realkind), dimension(0:npol) :: y0l
  real(kind=realkind), dimension(0:npol) :: v0_s_etal, v0_s_xil
  real(kind=realkind), dimension(0:npol) :: v0_z_etal, v0_z_xil
  real(kind=realkind), dimension(0:npol) :: V1

  integer :: ielem, j

  do ielem = 1, nel_solid

     yl(:,:) = Y(:,:,ielem)
     v_s_etal(:,:) = V_s_eta(:,:,ielem)
     v_s_xil(:,:)  = V_s_xi(:,:,ielem)
     v_z_etal(:,:) = V_z_eta(:,:,ielem)
     v_z_xil(:,:)  = V_z_xi(:,:,ielem)

     r1(:,:) = 0
     r2(:,:) = 0
     r3(:,:) = 0
     r4(:,:) = 0
     r5(:,:) = 0
     r6(:,:) = 0

     ! sum memory variables first, then compute stiffness terms of the sum
     do j=1, n_sls_attenuation
        r1(:,:) = r1(:,:) + R(:,:,1,j,ielem)
        r2(:,:) = r2(:,:) + R(:,:,2,j,ielem)
        r3(:,:) = r3(:,:) + R(:,:,3,j,ielem)
        r4(:,:) = r4(:,:) + R(:,:,4,j,ielem)
        r5(:,:) = r5(:,:) + R(:,:,5,j,ielem)
        r6(:,:) = r6(:,:) + R(:,:,6,j,ielem)
     enddo

     S1s = v_z_etal * r1 + v_s_etal * r5
     S2s = v_z_xil  * r1 + v_s_xil  * r5

     S1p = v_z_etal * r6 + v_s_etal * r4
     S2p = v_z_xil  * r6 + v_s_xil  * r4

     S1z = v_z_etal * r5 + v_s_etal * r3
     S2z = v_z_xil  * r5 + v_s_xil  * r3

     if (.not. axis_solid(ielem) ) then
        call mxm(G2, S1s, X1)
        call mxm(G2, S1p, X3)
        call mxm(G2, S1z, X5)
     else
        call mxm(G1, S1s, X1)
        call mxm(G1, S1p, X3)
        call mxm(G1, S1z, X5)
     endif

     call mxm(S2s, G2T, X2)
     call mxm(S2p, G2T, X4)
     call mxm(S2z, G2T, X6)

     loc_stiffness_s = X1 + X2 + yl * (r2 - 2 * r6)
     loc_stiffness_p = -X3 - X4 + yl * (r6 - 2 * r2)
     loc_stiffness_z = X5 + X6 - yl * r4

     if (axis_solid(ielem)) then
        y0l(:) = Y0(:,ielem)
        v0_s_etal(:) = V0_s_eta(:,ielem)
        v0_s_xil(:)  = V0_s_xi(:,ielem)
        v0_z_etal(:) = V0_z_eta(:,ielem)
        v0_z_xil(:)  = V0_z_xi(:,ielem)

        ! s - component
        V1 = v0_z_etal * r1(0,:) + y0l * (r2(0,:) - 2 * r6(0,:))
        loc_stiffness_s = loc_stiffness_s + outerprod(G0, V1)

        ! p - component
        V1 = - v0_z_etal * r6(0,:) + y0l * (r6(0,:) - 2 * r2(0,:))
        loc_stiffness_p = loc_stiffness_p + outerprod(G0, V1)

        ! z - component
        V1 = v0_s_etal * r3(0,:)
        loc_stiffness_z = loc_stiffness_z + outerprod(G0, V1)
     endif

     ! subtracting (!) from the global stiffness
     glob_stiffness(0:npol,0:npol,ielem,1) = &
            glob_stiffness(0:npol,0:npol,ielem,1) - loc_stiffness_s
     glob_stiffness(0:npol,0:npol,ielem,2) = &
            glob_stiffness(0:npol,0:npol,ielem,2) - loc_stiffness_p
     glob_stiffness(0:npol,0:npol,ielem,3) = &
            glob_stiffness(0:npol,0:npol,ielem,3) - loc_stiffness_z
  enddo

end subroutine glob_anel_stiffness_quad_generic
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_quad_4(glob_stiffness, R)

  use attenuation, only: n_sls_attenuation
  use data_mesh, only: nel_solid

  ! I/O global arrays
  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:) !0:4,0:4,nel_solid,1:3)
  real(kind=realkind), intent(in)    :: R(0:,0:,:,:,:) !(0:4,0:4,6,n_sls_attenuation,nel_solid)

  ! local variables for all elements
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_s
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_p
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_z

  real(kind=realkind), dimension(0:4,0:4) :: r1, r2, r3, r4, r5, r6

  real(kind=realkind), dimension(0:4,0:4) :: yl
  real(kind=realkind), dimension(0:4,0:4) :: v_s_etal, v_s_xil
  real(kind=realkind), dimension(0:4,0:4) :: v_z_etal, v_z_xil

  real(kind=realkind), dimension(0:4,0:4) :: S1s, S2s
  real(kind=realkind), dimension(0:4,0:4) :: S1p, S2p
  real(kind=realkind), dimension(0:4,0:4) :: S1z, S2z
  real(kind=realkind), dimension(0:4,0:4) :: X1, X2, X3, X4, X5, X6

  real(kind=realkind), dimension(0:4) :: y0l
  real(kind=realkind), dimension(0:4) :: v0_s_etal, v0_s_xil
  real(kind=realkind), dimension(0:4) :: v0_z_etal, v0_z_xil
  real(kind=realkind), dimension(0:4) :: V1

  integer :: ielem, j

  do ielem = 1, nel_solid

     yl(:,:) = Y(:,:,ielem)
     v_s_etal(:,:) = V_s_eta(:,:,ielem)
     v_s_xil(:,:)  = V_s_xi(:,:,ielem)
     v_z_etal(:,:) = V_z_eta(:,:,ielem)
     v_z_xil(:,:)  = V_z_xi(:,:,ielem)

     r1(:,:) = 0
     r2(:,:) = 0
     r3(:,:) = 0
     r4(:,:) = 0
     r5(:,:) = 0
     r6(:,:) = 0

     ! sum memory variables first, then compute stiffness terms of the sum
     do j=1, n_sls_attenuation
        r1(:,:) = r1(:,:) + R(:,:,1,j,ielem)
        r2(:,:) = r2(:,:) + R(:,:,2,j,ielem)
        r3(:,:) = r3(:,:) + R(:,:,3,j,ielem)
        r4(:,:) = r4(:,:) + R(:,:,4,j,ielem)
        r5(:,:) = r5(:,:) + R(:,:,5,j,ielem)
        r6(:,:) = r6(:,:) + R(:,:,6,j,ielem)
     enddo

     S1s = v_z_etal * r1 + v_s_etal * r5
     S2s = v_z_xil  * r1 + v_s_xil  * r5

     S1p = v_z_etal * r6 + v_s_etal * r4
     S2p = v_z_xil  * r6 + v_s_xil  * r4

     S1z = v_z_etal * r5 + v_s_etal * r3
     S2z = v_z_xil  * r5 + v_s_xil  * r3

     if (.not. axis_solid(ielem) ) then
        call mxm_4(G2, S1s, X1)
        call mxm_4(G2, S1p, X3)
        call mxm_4(G2, S1z, X5)
     else
        call mxm_4(G1, S1s, X1)
        call mxm_4(G1, S1p, X3)
        call mxm_4(G1, S1z, X5)
     endif

     call mxm_4(S2s, G2T, X2)
     call mxm_4(S2p, G2T, X4)
     call mxm_4(S2z, G2T, X6)

     loc_stiffness_s = X1 + X2 + yl * (r2 - 2 * r6)
     loc_stiffness_p = -X3 - X4 + yl * (r6 - 2 * r2)
     loc_stiffness_z = X5 + X6 - 2 * yl * r4

     if (axis_solid(ielem)) then
        y0l(:) = Y0(:,ielem)
        v0_s_etal(:) = V0_s_eta(:,ielem)
        v0_s_xil(:)  = V0_s_xi(:,ielem)
        v0_z_etal(:) = V0_z_eta(:,ielem)
        v0_z_xil(:)  = V0_z_xi(:,ielem)

        ! s - component
        V1 = v0_z_etal * r1(0,:) + y0l * (r2(0,:) - 2 * r6(0,:))
        loc_stiffness_s = loc_stiffness_s + outerprod_4(G0, V1)

        ! p - component
        V1 = - v0_z_etal * r6(0,:) + y0l * (r6(0,:) - 2 * r2(0,:))
        loc_stiffness_p = loc_stiffness_p + outerprod_4(G0, V1)

        ! z - component
        V1 = v0_s_etal * r3(0,:)
        loc_stiffness_z = loc_stiffness_z + outerprod_4(G0, V1)
     endif

     ! subtracting (!) from the global stiffness
     glob_stiffness(0:4,0:4,ielem,1) = &
            glob_stiffness(0:4,0:4,ielem,1) - loc_stiffness_s
     glob_stiffness(0:4,0:4,ielem,2) = &
            glob_stiffness(0:4,0:4,ielem,2) - loc_stiffness_p
     glob_stiffness(0:4,0:4,ielem,3) = &
            glob_stiffness(0:4,0:4,ielem,3) - loc_stiffness_z
  enddo

end subroutine glob_anel_stiffness_quad_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_quad_cg4(glob_stiffness, R)

  use attenuation, only: n_sls_attenuation
  use data_mesh, only: npol, nel_solid

  ! I/O global arrays
  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: R(:,:,:,:) !(1:4,6,n_sls_attenuation,nel_solid)

  ! local variables for all elements
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_s
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_p
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_z

  real(kind=realkind), dimension(1:4) :: r1, r2, r3, r4, r5, r6

  real(kind=realkind), dimension(1:4) :: yl
  real(kind=realkind), dimension(1:4) :: v_s_etal, v_s_xil
  real(kind=realkind), dimension(1:4) :: v_z_etal, v_z_xil

  real(kind=realkind), dimension(1:4) :: S1s, S2s
  real(kind=realkind), dimension(1:4) :: S1p, S2p
  real(kind=realkind), dimension(1:4) :: S1z, S2z
  real(kind=realkind), dimension(0:4,0:4) :: X1, X2, X3, X4, X5, X6

  integer :: ielem, j

  do ielem = 1, nel_solid

     yl(:) = Y_cg4(:,ielem)
     v_s_etal(:) = V_s_eta_cg4(:,ielem)
     v_s_xil(:)  = V_s_xi_cg4(:,ielem)
     v_z_etal(:) = V_z_eta_cg4(:,ielem)
     v_z_xil(:)  = V_z_xi_cg4(:,ielem)

     r1(:) = 0
     r2(:) = 0
     r3(:) = 0
     r4(:) = 0
     r5(:) = 0
     r6(:) = 0

     ! sum memory variables first, then compute stiffness terms of the sum
     do j=1, n_sls_attenuation
        r1(:) = r1(:) + R(:,1,j,ielem)
        r2(:) = r2(:) + R(:,2,j,ielem)
        r3(:) = r3(:) + R(:,3,j,ielem)
        r4(:) = r4(:) + R(:,4,j,ielem)
        r5(:) = r5(:) + R(:,5,j,ielem)
        r6(:) = r6(:) + R(:,6,j,ielem)
     enddo

     S1s = v_z_etal * r1 + v_s_etal * r5
     S2s = v_z_xil  * r1 + v_s_xil  * r5

     S1p = v_z_etal * r6 + v_s_etal * r4
     S2p = v_z_xil  * r6 + v_s_xil  * r4

     S1z = v_z_etal * r5 + v_s_etal * r3
     S2z = v_z_xil  * r5 + v_s_xil  * r3

     if (.not. axis_solid(ielem) ) then
        call mxm_cg4_sparse_b(G2, S1s, X1)
        call mxm_cg4_sparse_b(G2, S1p, X3)
        call mxm_cg4_sparse_b(G2, S1z, X5)
     else
        call mxm_cg4_sparse_b(G1, S1s, X1)
        call mxm_cg4_sparse_b(G1, S1p, X3)
        call mxm_cg4_sparse_b(G1, S1z, X5)
     endif

     call mxm_cg4_sparse_a(S2s, G2T, X2)
     call mxm_cg4_sparse_a(S2p, G2T, X4)
     call mxm_cg4_sparse_a(S2z, G2T, X6)

     loc_stiffness_s = X1 + X2
     loc_stiffness_s(1,1) = loc_stiffness_s(1,1) + yl(1) * (r2(1) - 2 * r6(1))
     loc_stiffness_s(1,3) = loc_stiffness_s(1,3) + yl(2) * (r2(2) - 2 * r6(2))
     loc_stiffness_s(3,1) = loc_stiffness_s(3,1) + yl(3) * (r2(3) - 2 * r6(3))
     loc_stiffness_s(3,3) = loc_stiffness_s(3,3) + yl(4) * (r2(4) - 2 * r6(4))

     loc_stiffness_p = -X3 - X4
     loc_stiffness_p(1,1) = loc_stiffness_p(1,1) + yl(1) * (r6(1) - 2 * r2(1))
     loc_stiffness_p(1,3) = loc_stiffness_p(1,3) + yl(2) * (r6(2) - 2 * r2(2))
     loc_stiffness_p(3,1) = loc_stiffness_p(3,1) + yl(3) * (r6(3) - 2 * r2(3))
     loc_stiffness_p(3,3) = loc_stiffness_p(3,3) + yl(4) * (r6(4) - 2 * r2(4))

     loc_stiffness_z = X5 + X6
     loc_stiffness_z(1,1) = loc_stiffness_z(1,1) - 2 * yl(1) * r4(1)
     loc_stiffness_z(1,3) = loc_stiffness_z(1,3) - 2 * yl(2) * r4(2)
     loc_stiffness_z(3,1) = loc_stiffness_z(3,1) - 2 * yl(3) * r4(3)
     loc_stiffness_z(3,3) = loc_stiffness_z(3,3) - 2 * yl(4) * r4(4)

     ! subtracting (!) from the global stiffness
     glob_stiffness(0:4,0:4,ielem,1) = &
            glob_stiffness(0:4,0:4,ielem,1) - loc_stiffness_s
     glob_stiffness(0:4,0:4,ielem,2) = &
            glob_stiffness(0:4,0:4,ielem,2) - loc_stiffness_p
     glob_stiffness(0:4,0:4,ielem,3) = &
            glob_stiffness(0:4,0:4,ielem,3) - loc_stiffness_z
  enddo

end subroutine glob_anel_stiffness_quad_cg4
!-----------------------------------------------------------------------------------------

end module stiffness_quad
!=========================================================================================
