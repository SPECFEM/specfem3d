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
module stiffness_mono

  use global_parameters, only: realkind
  use data_matr
  use data_mesh, only: axis_solid, axis_fluid, nsize
  use data_spec
  use data_source

  use unrolled_loops

  implicit none

  public :: glob_stiffness_mono
  public :: glob_anel_stiffness_mono

  private

contains

!-----------------------------------------------------------------------------------------
!> Wrapper routine to avoid if statements in the timeloop
pure subroutine glob_stiffness_mono(glob_stiffness,u)
  use data_mesh, only: npol, nel_solid

  real(kind=realkind), intent(in)  :: u(0:,0:,:,:)
  real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,3)

  if (npol == 4) then
     call glob_stiffness_mono_4(glob_stiffness, u)
  else
     call glob_stiffness_mono_generic(glob_stiffness, u)
  endif

end subroutine glob_stiffness_mono
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_stiffness_mono_4(glob_stiffness,u)

  use data_mesh, only: nel_solid

  integer, parameter               :: npol = 4
  ! I/O global arrays
  real(kind=realkind), intent(in)  :: u(0:npol,0:npol,nel_solid,3)
  real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,3)

  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z
  real(kind=realkind), dimension(0:npol,0:npol) :: us, uz

  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4     ! MxM arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: S1s, S2s, S1z, S2z ! Sum arrays

  real(kind=realkind), dimension(0:npol) :: V1, V2, V3, V4
  real(kind=realkind), dimension(0:npol) :: uz0

  integer :: ielem

  do ielem = 1, nel_solid

     us(:,:) = u(:,:,ielem,1)
     uz(:,:) = u(:,:,ielem,3)

     if (.not. axis_solid(ielem) ) then
        call mxm_4(G2T, us, X1)
        call mxm_4(G2T, uz, X2)
     else
        call mxm_4(G1T, us, X1)
        call mxm_4(G1T, uz, X2)
     endif

     call mxm_4(us, G2, X3)
     call mxm_4(uz, G2, X4)

     ! lower order terms in s
     loc_stiffness_s = m_4(:,:,ielem) * X4 + m_2(:,:,ielem) * X3 &
                     + m_1(:,:,ielem) * X1 + m_3(:,:,ielem) * X2 + us * m_w1(:,:,ielem)

     ! higher order terms + lower order terms with D_xi mxm ()
     S1s = m11s(:,:,ielem) * X3 + m21s(:,:,ielem) * X1 + m12s(:,:,ielem) * X4 &
         + m22s(:,:,ielem) * X2 + m_1(:,:,ielem) * us
     S2s = m11s(:,:,ielem) * X1 + m41s(:,:,ielem) * X3 + m32s(:,:,ielem) * X2 &
         + m42s(:,:,ielem) * X4 + m_2(:,:,ielem) * us
     S1z = m11z(:,:,ielem) * X4 + m21z(:,:,ielem) * X2 + m32s(:,:,ielem) * X3 &
         + m22s(:,:,ielem) * X1 + m_3(:,:,ielem) * us
     S2z = m11z(:,:,ielem) * X2 + m41z(:,:,ielem) * X4 + m12s(:,:,ielem) * X1 &
         + m42s(:,:,ielem) * X3 + m_4(:,:,ielem) * us

     call mxm_4(S2s, G2T, X2)
     call mxm_4(S2z, G2T, X4)

     if (.not. axis_solid(ielem) ) then
        call mxm_4(G2, S1s, X1)
        call mxm_4(G2, S1z, X3)
     else
        call mxm_4(G1, S1s, X1)
        call mxm_4(G1, S1z, X3)
     endif

     loc_stiffness_s = loc_stiffness_s + X1 + X2
     loc_stiffness_z = X3 + X4

     ! additional axis terms
     if (axis_solid(ielem) ) then
        uz0 = uz(0,:)
        X2 = 0

        call vxm_4(G0, us, V1)
        call vxm_4(uz0, G2, V2)

        V4 = m0_w1(:,ielem) * V1 + m0_w3(:,ielem) * V2

        ! additional anisotropic terms
        call vxm_4(G0, uz, V3)

        V4 = V4 + m0_w2(:,ielem) * V3
        X2 = outerprod_4(G0, m0_w2(:,ielem) * V1)
        ! end additional anisotropic terms

        V2 = m0_w3(:,ielem) * V1
        call vxm_4(V2, G2T, V1)
        X2(0,:) = X2(0,:) + V1

        loc_stiffness_s = loc_stiffness_s + outerprod_4(G0, V4)
        loc_stiffness_z = X2 + loc_stiffness_z
     endif

     glob_stiffness(:,:,ielem,1) = loc_stiffness_s
     glob_stiffness(:,:,ielem,3) = loc_stiffness_z

  enddo

end subroutine glob_stiffness_mono_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_stiffness_mono_generic(glob_stiffness,u)

  use data_mesh, only: npol, nel_solid

  ! I/O global arrays
  real(kind=realkind), intent(in)  :: u(0:,0:,:,:)
  real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,3)

  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z
  real(kind=realkind), dimension(0:npol,0:npol) :: us, uz
  real(kind=realkind), dimension(0:npol,0:npol) :: m_w1l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_1l, m_2l, m_3l, m_4l
  real(kind=realkind), dimension(0:npol,0:npol) :: m11sl, m21sl, m41sl, m12sl, m22sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m32sl, m42sl, m11zl, m21zl, m41zl

  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w1l, m0_w2l, m0_w3l

  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4     ! MxM arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: S1s, S2s, S1z, S2z ! Sum arrays

  real(kind=realkind), dimension(0:npol) :: V1, V2, V3, V4
  real(kind=realkind), dimension(0:npol) :: uz0

  integer :: ielem

  do ielem = 1, nel_solid

     us(:,:) = u(:,:,ielem,1)
     uz(:,:) = u(:,:,ielem,3)

     m_1l(:,:) = M_1(:,:,ielem)
     m_2l(:,:) = M_2(:,:,ielem)
     m_3l(:,:) = M_3(:,:,ielem)
     m_4l(:,:) = M_4(:,:,ielem)

     m_w1l(:,:) = M_w1(:,:,ielem)

     m11sl(:,:) = M11s(:,:,ielem)
     m21sl(:,:) = M21s(:,:,ielem)
     m41sl(:,:) = M41s(:,:,ielem)
     m12sl(:,:) = M12s(:,:,ielem)
     m22sl(:,:) = M22s(:,:,ielem)
     m32sl(:,:) = M32s(:,:,ielem)
     m42sl(:,:) = M42s(:,:,ielem)
     m11zl(:,:) = M11z(:,:,ielem)
     m21zl(:,:) = M21z(:,:,ielem)
     m41zl(:,:) = M41z(:,:,ielem)

     if (.not. axis_solid(ielem) ) then
        call mxm(G2T, us, X1)
        call mxm(G2T, uz, X2)
     else
        call mxm(G1T, us, X1)
        call mxm(G1T, uz, X2)
     endif

     call mxm(us, G2, X3)
     call mxm(uz, G2, X4)

     ! lower order terms in s
     loc_stiffness_s = m_4l * X4 + m_2l * X3 + m_1l * X1 + m_3l * X2 + us * m_w1l

     ! higher order terms + lower order terms with D_xi mxm ()
     S1s = m11sl * X3 + m21sl * X1 + m12sl * X4 + m22sl * X2 + m_1l * us
     S2s = m11sl * X1 + m41sl * X3 + m32sl * X2 + m42sl * X4 + m_2l * us
     S1z = m11zl * X4 + m21zl * X2 + m32sl * X3 + m22sl * X1 + m_3l * us
     S2z = m11zl * X2 + m41zl * X4 + m12sl * X1 + m42sl * X3 + m_4l * us

     call mxm(S2s, G2T, X2)
     call mxm(S2z, G2T, X4)

     if (.not. axis_solid(ielem) ) then
        call mxm(G2, S1s, X1)
        call mxm(G2, S1z, X3)
     else
        call mxm(G1, S1s, X1)
        call mxm(G1, S1z, X3)
     endif

     loc_stiffness_s = loc_stiffness_s + X1 + X2
     loc_stiffness_z = X3 + X4

     ! additional axis terms
     if (axis_solid(ielem) ) then
        m0_w1l(:) = M0_w1(:,ielem)
        m0_w2l(:) = M0_w2(:,ielem)
        m0_w3l(:) = M0_w3(:,ielem)

        uz0 = uz(0,:)
        X2 = 0

        call vxm(G0, us, V1)
        call vxm(uz0, G2, V2)

        V4 = m0_w1l * V1 + m0_w3l * V2

        ! additional anisotropic terms
        call vxm(G0, uz, V3)

        V4 = V4 + m0_w2l * V3
        X2 = outerprod(G0, m0_w2l * V1)
        ! end additional anisotropic terms

        V2 = m0_w3l * V1
        call vxm(V2, G2T, V1)
        X2(0,:) = X2(0,:) + V1

        loc_stiffness_s = loc_stiffness_s + outerprod(G0, V4)
        loc_stiffness_z = X2 + loc_stiffness_z
     endif

     glob_stiffness(:,:,ielem,1) = loc_stiffness_s
     glob_stiffness(:,:,ielem,3) = loc_stiffness_z

  enddo

end subroutine glob_stiffness_mono_generic
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_mono(glob_stiffness, R, R_cg, cg)
  use data_mesh, only: npol

  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: R(0:,0:,:,:,:)
  real(kind=realkind), intent(in)    :: R_cg(:,:,:,:)
  logical, intent(in)                :: cg

  if (cg) then
     call glob_anel_stiffness_mono_cg4(glob_stiffness, R_cg)
  else
     if (npol == 4) then
        call glob_anel_stiffness_mono_4(glob_stiffness, R)
     else
        call glob_anel_stiffness_mono_generic(glob_stiffness, R)
     endif
  endif

end subroutine glob_anel_stiffness_mono
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_mono_generic(glob_stiffness, R)

  use attenuation, only: n_sls_attenuation
  use data_mesh, only: npol, nel_solid

  ! I/O global arrays
  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: R(0:,0:,:,:,:)

  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z

  real(kind=realkind), dimension(0:npol,0:npol) :: r1, r2, r3, r5

  real(kind=realkind), dimension(0:npol,0:npol) :: yl
  real(kind=realkind), dimension(0:npol,0:npol) :: v_s_etal, v_s_xil
  real(kind=realkind), dimension(0:npol,0:npol) :: v_z_etal, v_z_xil

  real(kind=realkind), dimension(0:npol,0:npol) :: S1s, S2s
  real(kind=realkind), dimension(0:npol,0:npol) :: S1z, S2z
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4

  real(kind=realkind), dimension(0:npol) :: y0l
  real(kind=realkind), dimension(0:npol) :: v0_s_etal, v0_s_xil
  real(kind=realkind), dimension(0:npol) :: v0_z_etal, v0_z_xil
  real(kind=realkind), dimension(0:npol) :: V1, V2, V3, V4

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
     r5(:,:) = 0

     ! sum memory variables first, then compute stiffness terms of the sum
     do j=1, n_sls_attenuation
        r1(:,:) = r1(:,:) + R(:,:,1,j,ielem)
        r2(:,:) = r2(:,:) + R(:,:,2,j,ielem)
        r3(:,:) = r3(:,:) + R(:,:,3,j,ielem)
        r5(:,:) = r5(:,:) + R(:,:,5,j,ielem)
     enddo

     S1s = v_z_etal * r1 + v_s_etal * r5
     S2s = v_z_xil  * r1 + v_s_xil  * r5

     S1z = v_z_etal * r5 + v_s_etal * r3
     S2z = v_z_xil  * r5 + v_s_xil  * r3

     if (.not. axis_solid(ielem) ) then
        call mxm(G2, S1s, X1)
        call mxm(G2, S1z, X3)
     else
        call mxm(G1, S1s, X1)
        call mxm(G1, S1z, X3)
     endif

     call mxm(S2s, G2T, X2)
     call mxm(S2z, G2T, X4)

     loc_stiffness_s = X1 + X2 + yl * r2
     loc_stiffness_z = X3 + X4

     if (axis_solid(ielem)) then
        y0l(:) = Y0(:,ielem)
        v0_s_etal(:) = V0_s_eta(:,ielem)
        v0_s_xil(:)  = V0_s_xi(:,ielem)
        v0_z_etal(:) = V0_z_eta(:,ielem)
        v0_z_xil(:)  = V0_z_xi(:,ielem)

        ! s - component
        V1 = v0_z_etal * r1(0,:) + v0_s_etal * r5(0,:) + y0l * r2(0,:)
        loc_stiffness_s = loc_stiffness_s + outerprod(G0, V1)

        ! z - component
        V2 = v0_z_etal * r5(0,:) + v0_s_etal * r3(0,:)
        V3 = v0_z_xil  * r5(0,:) + v0_s_xil  * r3(0,:)
        call vxm(V3, G2T, V4)

        loc_stiffness_z = loc_stiffness_z + outerprod(G0, V2)
        loc_stiffness_z(0,:) = loc_stiffness_z(0,:) + V4
     endif

     ! subtracting (!) from the global stiffness
     glob_stiffness(0:npol,0:npol,ielem,1) = &
            glob_stiffness(0:npol,0:npol,ielem,1) - loc_stiffness_s
     glob_stiffness(0:npol,0:npol,ielem,3) = &
            glob_stiffness(0:npol,0:npol,ielem,3) - loc_stiffness_z
  enddo

end subroutine glob_anel_stiffness_mono_generic
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_mono_4(glob_stiffness, R)

  use attenuation, only: n_sls_attenuation
  use data_mesh, only: nel_solid

  ! I/O global arrays
  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: R(0:,0:,:,:,:)

  ! local variables for all elements
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_s
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_z

  real(kind=realkind), dimension(0:4,0:4) :: r1, r2, r3, r5

  real(kind=realkind), dimension(0:4,0:4) :: yl
  real(kind=realkind), dimension(0:4,0:4) :: v_s_etal, v_s_xil
  real(kind=realkind), dimension(0:4,0:4) :: v_z_etal, v_z_xil

  real(kind=realkind), dimension(0:4,0:4) :: S1s, S2s
  real(kind=realkind), dimension(0:4,0:4) :: S1z, S2z
  real(kind=realkind), dimension(0:4,0:4) :: X1, X2, X3, X4

  real(kind=realkind), dimension(0:4) :: y0l
  real(kind=realkind), dimension(0:4) :: v0_s_etal, v0_s_xil
  real(kind=realkind), dimension(0:4) :: v0_z_etal, v0_z_xil
  real(kind=realkind), dimension(0:4) :: V1, V2, V3, V4

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
     r5(:,:) = 0

     ! sum memory variables first, then compute stiffness terms of the sum
     do j=1, n_sls_attenuation
        r1(:,:) = r1(:,:) + R(:,:,1,j,ielem)
        r2(:,:) = r2(:,:) + R(:,:,2,j,ielem)
        r3(:,:) = r3(:,:) + R(:,:,3,j,ielem)
        r5(:,:) = r5(:,:) + R(:,:,5,j,ielem)
     enddo

     S1s = v_z_etal * r1 + v_s_etal * r5
     S2s = v_z_xil  * r1 + v_s_xil  * r5

     S1z = v_z_etal * r5 + v_s_etal * r3
     S2z = v_z_xil  * r5 + v_s_xil  * r3

     if (.not. axis_solid(ielem) ) then
        call mxm_4(G2, S1s, X1)
        call mxm_4(G2, S1z, X3)
     else
        call mxm_4(G1, S1s, X1)
        call mxm_4(G1, S1z, X3)
     endif

     call mxm_4(S2s, G2T, X2)
     call mxm_4(S2z, G2T, X4)

     loc_stiffness_s = X1 + X2 + yl * r2
     loc_stiffness_z = X3 + X4

     if (axis_solid(ielem)) then
        y0l(:) = Y0(:,ielem)
        v0_s_etal(:) = V0_s_eta(:,ielem)
        v0_s_xil(:)  = V0_s_xi(:,ielem)
        v0_z_etal(:) = V0_z_eta(:,ielem)
        v0_z_xil(:)  = V0_z_xi(:,ielem)

        ! s - component
        V1 = v0_z_etal * r1(0,:) + v0_s_etal * r5(0,:) + y0l * r2(0,:)
        loc_stiffness_s = loc_stiffness_s + outerprod_4(G0, V1)

        ! z - component
        V2 = v0_z_etal * r5(0,:) + v0_s_etal * r3(0,:)
        V3 = v0_z_xil  * r5(0,:) + v0_s_xil  * r3(0,:)
        call vxm_4(V3, G2T, V4)

        loc_stiffness_z = loc_stiffness_z + outerprod_4(G0, V2)
        loc_stiffness_z(0,:) = loc_stiffness_z(0,:) + V4
     endif

     ! subtracting (!) from the global stiffness
     glob_stiffness(0:4,0:4,ielem,1) = &
            glob_stiffness(0:4,0:4,ielem,1) - loc_stiffness_s
     glob_stiffness(0:4,0:4,ielem,3) = &
            glob_stiffness(0:4,0:4,ielem,3) - loc_stiffness_z
  enddo

end subroutine glob_anel_stiffness_mono_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_mono_cg4(glob_stiffness, R)

  use attenuation, only: n_sls_attenuation
  use data_mesh, only: nel_solid

  ! I/O global arrays
  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: R(:,:,:,:)

  ! local variables for all elements
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_s
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_z

  real(kind=realkind), dimension(1:4) :: r1, r2, r3, r5

  real(kind=realkind), dimension(1:4) :: yl
  real(kind=realkind), dimension(1:4) :: v_s_etal, v_s_xil
  real(kind=realkind), dimension(1:4) :: v_z_etal, v_z_xil

  real(kind=realkind), dimension(1:4) :: S1s, S2s
  real(kind=realkind), dimension(1:4) :: S1z, S2z
  real(kind=realkind), dimension(0:4,0:4) :: X1, X2, X3, X4

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
     r5(:) = 0

     ! sum memory variables first, then compute stiffness terms of the sum
     do j=1, n_sls_attenuation
        r1(:) = r1(:) + R(:,1,j,ielem)
        r2(:) = r2(:) + R(:,2,j,ielem)
        r3(:) = r3(:) + R(:,3,j,ielem)
        r5(:) = r5(:) + R(:,5,j,ielem)
     enddo

     S1s = v_z_etal * r1 + v_s_etal * r5
     S2s = v_z_xil  * r1 + v_s_xil  * r5

     S1z = v_z_etal * r5 + v_s_etal * r3
     S2z = v_z_xil  * r5 + v_s_xil  * r3

     if (.not. axis_solid(ielem) ) then
        call mxm_cg4_sparse_b(G2,  S1s, X1)
        call mxm_cg4_sparse_b(G2,  S1z, X3)
     else
        call mxm_cg4_sparse_b(G1,  S1s, X1)
        call mxm_cg4_sparse_b(G1,  S1z, X3)
     endif

     call mxm_cg4_sparse_a(S2s, G2T, X2)
     call mxm_cg4_sparse_a(S2z, G2T, X4)

     loc_stiffness_s = X1 + X2
     loc_stiffness_s(1,1) = loc_stiffness_s(1,1) + yl(1) * r2(1)
     loc_stiffness_s(1,3) = loc_stiffness_s(1,3) + yl(2) * r2(2)
     loc_stiffness_s(3,1) = loc_stiffness_s(3,1) + yl(3) * r2(3)
     loc_stiffness_s(3,3) = loc_stiffness_s(3,3) + yl(4) * r2(4)

     loc_stiffness_z = X3 + X4

     ! subtracting (!) from the global stiffness
     glob_stiffness(0:4,0:4,ielem,1) = &
            glob_stiffness(0:4,0:4,ielem,1) - loc_stiffness_s
     glob_stiffness(0:4,0:4,ielem,3) = &
            glob_stiffness(0:4,0:4,ielem,3) - loc_stiffness_z
  enddo

end subroutine glob_anel_stiffness_mono_cg4
!-----------------------------------------------------------------------------------------

end module stiffness_mono
!=========================================================================================
