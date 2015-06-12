!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
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
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!====================
module stiffness
!====================

  use global_parameters, only: realkind
  use data_matr
  use data_mesh,         only: axis_solid, axis_fluid, nsize
  use data_spec
  use data_source
  
  use unrolled_loops
  
  implicit none
  
  public :: glob_stiffness_mono
  public :: glob_anel_stiffness_mono

  public :: glob_stiffness_di
  public :: glob_anel_stiffness_di

  public :: glob_stiffness_quad
  public :: glob_anel_stiffness_quad

  public :: glob_fluid_stiffness
  
  private

contains

!-----------------------------------------------------------------------------
pure function outerprod(a,b) 
  ! outer product (dyadic) from numerical recipes
  
  real(kind=realkind), dimension(:), intent(in)     :: a, b
  real(kind=realkind), dimension(size(a),size(b))   :: outerprod

  outerprod = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a))
end function outerprod
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
pure function outerprod_4(a,b) 
  ! outer product (dyadic) from numerical recipes
  
  real(kind=realkind), dimension(0:4), intent(in)   :: a, b
  real(kind=realkind), dimension(0:4,0:4)           :: outerprod_4

  outerprod_4 = spread(a, dim=2, ncopies=5) * spread(b, dim=1, ncopies=5)
end function outerprod_4
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
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

     if ( .not. axis_solid(ielem) ) then
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

     if ( .not. axis_solid(ielem) ) then
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

  end do

end subroutine glob_stiffness_mono_4
!=============================================================================

!-----------------------------------------------------------------------------
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

     if ( .not. axis_solid(ielem) ) then
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

     if ( .not. axis_solid(ielem) ) then
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
!=============================================================================

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_mono(glob_stiffness, R, R_cg, cg)
  use data_mesh,    only: npol
  
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
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_mono_generic(glob_stiffness, R)

  use attenuation,  only: n_sls_attenuation
  use data_mesh,    only: npol, nel_solid
  
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

     if ( .not. axis_solid(ielem) ) then
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
!=============================================================================

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_mono_4(glob_stiffness, R)

  use attenuation,  only: n_sls_attenuation
  use data_mesh,    only: nel_solid
  
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

     if ( .not. axis_solid(ielem) ) then
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
!=============================================================================

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_mono_cg4(glob_stiffness, R)

  use attenuation,  only: n_sls_attenuation
  use data_mesh,    only: nel_solid
  
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

     if ( .not. axis_solid(ielem) ) then
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
!=============================================================================

!-----------------------------------------------------------------------------
!> Wrapper routine to avoid if statements in the timeloop
pure subroutine glob_stiffness_di(glob_stiffness,u)
  use data_mesh, only: npol, nel_solid
  
  real(kind=realkind), intent(in)  :: u(0:,0:,:,:)
  real(kind=realkind), intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,3)

  if (npol == 4) then
     call glob_stiffness_di_4(glob_stiffness, u)
  else
     call glob_stiffness_di_generic(glob_stiffness, u)
  endif

end subroutine glob_stiffness_di
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
pure subroutine glob_stiffness_di_4(glob_stiffness,u)

  use data_mesh, only: nel_solid

  integer, parameter              :: npol = 4
  ! I/O for global arrays
  real(kind=realkind),intent(in)  :: u(0:npol,0:npol,nel_solid,3)
  real(kind=realkind),intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,3)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_1
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_2
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_3
  real(kind=realkind), dimension(0:npol,0:npol) :: u1,u2,u3
  real(kind=realkind), dimension(0:npol,0:npol) :: m_w1l, m_w2l, m_w3l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_6l, m_2l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_1l, m_5l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_4l, m_8l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_3l, m_7l
  
  real(kind=realkind), dimension(0:npol,0:npol) :: m11sl, m21sl, m41sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m12sl, m22sl, m42sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m13sl, m23sl, m33sl, m43sl
  
  real(kind=realkind), dimension(0:npol,0:npol) :: m11zl, m21zl, m41zl
  
  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w1l, m0_w2l, m0_w3l, m0_w4l
  real(kind=realkind), dimension(0:npol) :: m0_w6l, m0_w7l, m0_w8l, m0_w9l, m0_w10l
  
  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s2, loc_stiffness_s3
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4, X5, X6, X7, X8 ! MxM 
  real(kind=realkind), dimension(0:npol,0:npol) :: S1p, S1m, S2p, S2m, S1z, S2z ! Sum
  
  real(kind=realkind), dimension(0:npol,0:npol) :: c1, c2, c3
  
  real(kind=realkind), dimension(0:npol) :: V1, V2, V3, V4, V5
  real(kind=realkind), dimension(0:npol) :: u10, u20
  
  integer :: ielem

  do ielem = 1, nel_solid

     u1(:,:) = u(:,:,ielem,1)
     u2(:,:) = u(:,:,ielem,2)
     u3(:,:) = u(:,:,ielem,3)
     
     m_1l(:,:) = M_1(:,:,ielem)
     m_2l(:,:) = M_2(:,:,ielem)
     m_3l(:,:) = M_3(:,:,ielem)
     m_4l(:,:) = M_4(:,:,ielem)
     m_5l(:,:) = M_5(:,:,ielem)
     m_6l(:,:) = M_6(:,:,ielem)
     m_7l(:,:) = M_7(:,:,ielem)
     m_8l(:,:) = M_8(:,:,ielem)
     
     m_w1l(:,:) = M_w1(:,:,ielem)
     m_w2l(:,:) = M_w2(:,:,ielem)
     m_w3l(:,:) = M_w3(:,:,ielem)

     m11sl(:,:) = M11s(:,:,ielem)
     m21sl(:,:) = M21s(:,:,ielem)
     m41sl(:,:) = M41s(:,:,ielem)
     m12sl(:,:) = M12s(:,:,ielem)
     m22sl(:,:) = M22s(:,:,ielem)
     m42sl(:,:) = M42s(:,:,ielem)
     m13sl(:,:) = M13s(:,:,ielem)
     m23sl(:,:) = M32s(:,:,ielem) ! correct!! (static memory reasons,
                                  ! reusing static array from
                                  ! monopole)
     m33sl(:,:) = M33s(:,:,ielem)
     m43sl(:,:) = M43s(:,:,ielem)

     m11zl(:,:) = M11z(:,:,ielem)
     m21zl(:,:) = M21z(:,:,ielem)
     m41zl(:,:) = M41z(:,:,ielem)

     ! First MxM
     call mxm_4(u1, G2, X4)
     call mxm_4(u2, G2, X5)
     call mxm_4(u3, G2, X6)

     if ( .not. axis_solid(ielem) ) then
        call mxm_4(G2T, u1, X1)
        call mxm_4(G2T, u2, X2)
        call mxm_4(G2T, u3, X3)
     else
        call mxm_4(G1T, u1, X1)
        call mxm_4(G1T, u2, X2)
        call mxm_4(G1T, u3, X3)
     endif

     ! Sum for the z-component
     X7 = X1 + X2
     X8 = X4 + X5 

     ! Collocations and sums of W_x and W_x^d terms
     ! - component
     loc_stiffness_s2 = m_8l * X6 + m_7l * X3 + m_1l  * X1 + m_5l  * X2 &
                      + m_2l * X4 + m_6l * X5 + m_w1l * u2 + m_w2l * u3

     ! z component
     loc_stiffness_s3 = m_4l  * X4 - m_4l  * X5 + m_3l * X1 - m_3l * X2 &
                      + m_w2l * u2 + m_w3l * u3
        
     ! + and -
     ! buffering reused terms
     c1 = m13sl * X6
     c2 = m23sl * X3
     c3 = m_3l * u3
     
     s1p = c1 + c2 + c3 + m11sl * X4 + m21sl * X1 + m12sl * X5 + m22sl * X2 + m_1l * u2
     s1m = c1 + c2 - c3 + m11sl * X5 + m21sl * X2 + m12sl * X4 + m22sl * X1 + m_5l * u2
     
     c1 = m33sl * X3
     c2 = m43sl * X6
     c3 = m_4l * u3

     s2p = c1 + c2 + c3 + m11sl * X1 + m41sl * X4 + m12sl * X2 + m42sl * X5 + m_2l * u2
     s2m = c1 + c2 - c3 + m11sl * X2 + m41sl * X5 + m12sl * X1 + m42sl * X4 + m_6l * u2

     ! z component
     S1z = m33sl * X8 + m23sl * X7 + m11zl * X6 + m21zl * X3 + m_7l * u2
     S2z = m13sl * X7 + m43sl * X8 + m11zl * X3 + m41zl * X6 + m_8l * u2

     ! Second MxM
     if ( .not. axis_solid(ielem) ) then
        call mxm_4(G2, S1p, X1)
        call mxm_4(G2, S1m, X3)
        call mxm_4(G2, S1z, X5)
     else
        call mxm_4(G1, S1p, X1)
        call mxm_4(G1, S1m, X3)
        call mxm_4(G1, S1z, X5)
     endif
        
     call mxm_4(S2p, G2T, X2)
     call mxm_4(S2m, G2T, X4)
     call mxm_4(S2z, G2T, X6)

     loc_stiffness_1 = X1 + X2
     loc_stiffness_2 = X3 + X4 + loc_stiffness_s2
     loc_stiffness_3 = X5 + X6 + loc_stiffness_s3

     ! Additional terms for the axial elements
     if ( axis_solid(ielem) ) then
        m0_w1l(:)  = M0_w1(:,ielem)
        m0_w2l(:)  = M0_w2(:,ielem)
        m0_w7l(:)  = M0_w7(:,ielem)
        m0_w8l(:)  = M0_w8(:,ielem)
        m0_w9l(:)  = M0_w9(:,ielem)

        u10 = u1(0,:)
        u20 = u2(0,:)

        ! VxM
        call vxm_4(G0, u1, V1)
        call vxm_4(G0, u2, V2)
        call vxm_4(G0, u3, V3)

        call vxm_4(u10, G2, V4)
        call vxm_4(u20, G2, V5)

        ! zero in isotropic case
        m0_w3l(:)  = M0_w3(:,ielem)
        m0_w4l(:)  = M0_w4(:,ielem)
        m0_w6l(:)  = M0_w6(:,ielem)
        m0_w10l(:) = M0_w10(:,ielem)

        S1p = outerprod_4(G0, m0_w1l * V2 + m0_w3l * V3)
        
        S1m = outerprod_4(G0, m0_w1l * V1 + m0_w2l * V5 + m0_w6l  * V4 &
                                        + m0_w9l * V2 + m0_w10l * V3)
        
        S1z = outerprod_4(G0, m0_w3l * V1 + (m0_w4l + m0_w8l) * V4 + m0_w7l * V3 &
                                        + m0_w10l * V2)

        V4 = (m0_w2l + m0_w6l) * V2 + (m0_w4l + m0_w8l) * V3
        ! end additional anisotropic terms

        ! Final VxM in + component
        call vxm_4(V4, G2T, V1)
        S1p(0,:) = S1p(0,:) + V1
        
        loc_stiffness_1 = loc_stiffness_1 + S1p
        loc_stiffness_2 = loc_stiffness_2 + S1m
        loc_stiffness_3 = loc_stiffness_3 + S1z
     endif

     glob_stiffness(0:npol,0:npol,ielem,1) = loc_stiffness_1
     glob_stiffness(0:npol,0:npol,ielem,2) = loc_stiffness_2
     glob_stiffness(0:npol,0:npol,ielem,3) = loc_stiffness_3

  enddo

end subroutine glob_stiffness_di_4
!=============================================================================


!-----------------------------------------------------------------------------
pure subroutine glob_stiffness_di_generic(glob_stiffness,u)

  use data_mesh, only: npol, nel_solid

  ! I/O for global arrays
  real(kind=realkind),intent(in)  :: u(0:,0:,:,:)
  real(kind=realkind),intent(out) :: glob_stiffness(0:npol,0:npol,nel_solid,3)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_1
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_2
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_3
  real(kind=realkind), dimension(0:npol,0:npol) :: u1,u2,u3
  real(kind=realkind), dimension(0:npol,0:npol) :: m_w1l, m_w2l, m_w3l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_6l, m_2l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_1l, m_5l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_4l, m_8l
  real(kind=realkind), dimension(0:npol,0:npol) :: m_3l, m_7l
  
  real(kind=realkind), dimension(0:npol,0:npol) :: m11sl, m21sl, m41sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m12sl, m22sl, m42sl
  real(kind=realkind), dimension(0:npol,0:npol) :: m13sl, m23sl, m33sl, m43sl
  
  real(kind=realkind), dimension(0:npol,0:npol) :: m11zl, m21zl, m41zl
  
  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w1l, m0_w2l, m0_w3l, m0_w4l
  real(kind=realkind), dimension(0:npol) :: m0_w6l, m0_w7l, m0_w8l, m0_w9l, m0_w10l
  
  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_s2, loc_stiffness_s3
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4, X5, X6, X7, X8 ! MxM 
  real(kind=realkind), dimension(0:npol,0:npol) :: S1p, S1m, S2p, S2m, S1z, S2z ! Sum
  
  real(kind=realkind), dimension(0:npol,0:npol) :: c1, c2, c3
  
  real(kind=realkind), dimension(0:npol) :: V1, V2, V3, V4, V5
  real(kind=realkind), dimension(0:npol) :: u10, u20
  
  integer :: ielem

  do ielem = 1, nel_solid

     u1(:,:) = u(:,:,ielem,1)
     u2(:,:) = u(:,:,ielem,2)
     u3(:,:) = u(:,:,ielem,3)
     
     m_1l(:,:) = M_1(:,:,ielem)
     m_2l(:,:) = M_2(:,:,ielem)
     m_3l(:,:) = M_3(:,:,ielem)
     m_4l(:,:) = M_4(:,:,ielem)
     m_5l(:,:) = M_5(:,:,ielem)
     m_6l(:,:) = M_6(:,:,ielem)
     m_7l(:,:) = M_7(:,:,ielem)
     m_8l(:,:) = M_8(:,:,ielem)
     
     m_w1l(:,:) = M_w1(:,:,ielem)
     m_w2l(:,:) = M_w2(:,:,ielem)
     m_w3l(:,:) = M_w3(:,:,ielem)

     m11sl(:,:) = M11s(:,:,ielem)
     m21sl(:,:) = M21s(:,:,ielem)
     m41sl(:,:) = M41s(:,:,ielem)
     m12sl(:,:) = M12s(:,:,ielem)
     m22sl(:,:) = M22s(:,:,ielem)
     m42sl(:,:) = M42s(:,:,ielem)
     m13sl(:,:) = M13s(:,:,ielem)
     m23sl(:,:) = M32s(:,:,ielem) ! correct!! (static memory reasons,
                                  ! reusing static array from
                                  ! monopole)
     m33sl(:,:) = M33s(:,:,ielem)
     m43sl(:,:) = M43s(:,:,ielem)

     m11zl(:,:) = M11z(:,:,ielem)
     m21zl(:,:) = M21z(:,:,ielem)
     m41zl(:,:) = M41z(:,:,ielem)

     ! First MxM
     call mxm(u1, G2, X4)
     call mxm(u2, G2, X5)
     call mxm(u3, G2, X6)

     if ( .not. axis_solid(ielem) ) then
        call mxm(G2T, u1, X1)
        call mxm(G2T, u2, X2)
        call mxm(G2T, u3, X3)
     else
        call mxm(G1T, u1, X1)
        call mxm(G1T, u2, X2)
        call mxm(G1T, u3, X3)
     endif

     ! Sum for the z-component
     X7 = X1 + X2
     X8 = X4 + X5 

     ! Collocations and sums of W_x and W_x^d terms
     ! - component
     loc_stiffness_s2 = m_8l * X6 + m_7l * X3 + m_1l  * X1 + m_5l  * X2 &
                      + m_2l * X4 + m_6l * X5 + m_w1l * u2 + m_w2l * u3

     ! z component
     loc_stiffness_s3 = m_4l  * X4 - m_4l  * X5 + m_3l * X1 - m_3l * X2 &
                      + m_w2l * u2 + m_w3l * u3
        
     ! + and -
     ! buffering reused terms
     c1 = m13sl * X6
     c2 = m23sl * X3
     c3 = m_3l * u3
     
     s1p = c1 + c2 + c3 + m11sl * X4 + m21sl * X1 + m12sl * X5 + m22sl * X2 + m_1l * u2
     s1m = c1 + c2 - c3 + m11sl * X5 + m21sl * X2 + m12sl * X4 + m22sl * X1 + m_5l * u2
     
     c1 = m33sl * X3
     c2 = m43sl * X6
     c3 = m_4l * u3

     s2p = c1 + c2 + c3 + m11sl * X1 + m41sl * X4 + m12sl * X2 + m42sl * X5 + m_2l * u2
     s2m = c1 + c2 - c3 + m11sl * X2 + m41sl * X5 + m12sl * X1 + m42sl * X4 + m_6l * u2

     ! z component
     S1z = m33sl * X8 + m23sl * X7 + m11zl * X6 + m21zl * X3 + m_7l * u2
     S2z = m13sl * X7 + m43sl * X8 + m11zl * X3 + m41zl * X6 + m_8l * u2

     ! Second MxM
     if ( .not. axis_solid(ielem) ) then
        call mxm(G2, S1p, X1)
        call mxm(G2, S1m, X3)
        call mxm(G2, S1z, X5)
     else
        call mxm(G1, S1p, X1)
        call mxm(G1, S1m, X3)
        call mxm(G1, S1z, X5)
     endif
        
     call mxm(S2p, G2T, X2)
     call mxm(S2m, G2T, X4)
     call mxm(S2z, G2T, X6)

     loc_stiffness_1 = X1 + X2
     loc_stiffness_2 = X3 + X4 + loc_stiffness_s2
     loc_stiffness_3 = X5 + X6 + loc_stiffness_s3

     ! Additional terms for the axial elements
     if ( axis_solid(ielem) ) then
        m0_w1l(:)  = M0_w1(:,ielem)
        m0_w2l(:)  = M0_w2(:,ielem)
        m0_w7l(:)  = M0_w7(:,ielem)
        m0_w8l(:)  = M0_w8(:,ielem)
        m0_w9l(:)  = M0_w9(:,ielem)

        u10 = u1(0,:)
        u20 = u2(0,:)

        ! VxM
        call vxm(G0, u1, V1)
        call vxm(G0, u2, V2)
        call vxm(G0, u3, V3)

        call vxm(u10, G2, V4)
        call vxm(u20, G2, V5)

        ! zero in isotropic case
        m0_w3l(:)  = M0_w3(:,ielem)
        m0_w4l(:)  = M0_w4(:,ielem)
        m0_w6l(:)  = M0_w6(:,ielem)
        m0_w10l(:) = M0_w10(:,ielem)

        S1p = outerprod(G0, m0_w1l * V2 + m0_w3l * V3)
        
        S1m = outerprod(G0, m0_w1l * V1 + m0_w2l * V5 + m0_w6l  * V4 &
                                        + m0_w9l * V2 + m0_w10l * V3)
        
        S1z = outerprod(G0, m0_w3l * V1 + (m0_w4l + m0_w8l) * V4 + m0_w7l * V3 &
                                        + m0_w10l * V2)

        V4 = (m0_w2l + m0_w6l) * V2 + (m0_w4l + m0_w8l) * V3
        ! end additional anisotropic terms

        ! Final VxM in + component
        call vxm(V4, G2T, V1)
        S1p(0,:) = S1p(0,:) + V1
        
        loc_stiffness_1 = loc_stiffness_1 + S1p
        loc_stiffness_2 = loc_stiffness_2 + S1m
        loc_stiffness_3 = loc_stiffness_3 + S1z
     endif

     glob_stiffness(0:npol,0:npol,ielem,1) = loc_stiffness_1
     glob_stiffness(0:npol,0:npol,ielem,2) = loc_stiffness_2
     glob_stiffness(0:npol,0:npol,ielem,3) = loc_stiffness_3

  enddo

end subroutine glob_stiffness_di_generic
!=============================================================================

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_di(glob_stiffness, R, R_cg, cg)
  use data_mesh,    only: npol
  
  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: R(0:,0:,:,:,:)
  real(kind=realkind), intent(in)    :: R_cg(:,:,:,:)
  logical, intent(in)                :: cg

  if (cg) then
     call glob_anel_stiffness_di_cg4(glob_stiffness, R_cg)
  else
     if (npol == 4) then
        call glob_anel_stiffness_di_4(glob_stiffness, R)
     else
        call glob_anel_stiffness_di_generic(glob_stiffness, R)
     endif
  endif

end subroutine glob_anel_stiffness_di
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_di_generic(glob_stiffness, R)

  use attenuation,  only: n_sls_attenuation
  use data_mesh,    only: npol, nel_solid
  
  ! I/O global arrays
  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: R(0:,0:,:,:,:)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_p
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_m
  real(kind=realkind), dimension(0:npol,0:npol) :: loc_stiffness_z
  
  real(kind=realkind), dimension(0:npol,0:npol) :: r1, r2, r3, r4, r5, r6

  real(kind=realkind), dimension(0:npol,0:npol) :: yl
  real(kind=realkind), dimension(0:npol,0:npol) :: v_s_etal, v_s_xil
  real(kind=realkind), dimension(0:npol,0:npol) :: v_z_etal, v_z_xil

  real(kind=realkind), dimension(0:npol,0:npol) :: S1p, S2p
  real(kind=realkind), dimension(0:npol,0:npol) :: S1m, S2m
  real(kind=realkind), dimension(0:npol,0:npol) :: S1z, S2z
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2, X3, X4, X5, X6
  
  real(kind=realkind), dimension(0:npol) :: y0l
  real(kind=realkind), dimension(0:npol) :: v0_s_etal, v0_s_xil
  real(kind=realkind), dimension(0:npol) :: v0_z_etal, v0_z_xil
  real(kind=realkind), dimension(0:npol) :: V1, V2, V3
  
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

     S1p = v_z_etal * (r1 - r6) + v_s_etal * (r5 - r4)
     S2p = v_z_xil  * (r1 - r6) + v_s_xil  * (r5 - r4)
     
     S1m = v_z_etal * (r1 + r6) + v_s_etal * (r5 + r4)
     S2m = v_z_xil  * (r1 + r6) + v_s_xil  * (r5 + r4)
     
     S1z = v_z_etal * r5 + v_s_etal * r3
     S2z = v_z_xil  * r5 + v_s_xil  * r3

     if ( .not. axis_solid(ielem) ) then
        call mxm(G2, S1p, X1)
        call mxm(G2, S1m, X3)
        call mxm(G2, S1z, X5)
     else
        call mxm(G1, S1p, X1)
        call mxm(G1, S1m, X3)
        call mxm(G1, S1z, X5)
     endif

     call mxm(S2p, G2T, X2)
     call mxm(S2m, G2T, X4)
     call mxm(S2z, G2T, X6)

     loc_stiffness_p = X1 + X2
     loc_stiffness_m = X3 + X4 + 2 * yl * (r2 - r6)
     loc_stiffness_z = X5 + X6 - yl * r4

     if (axis_solid(ielem)) then
        y0l(:) = Y0(:,ielem)
        v0_s_etal(:) = V0_s_eta(:,ielem)
        v0_s_xil(:)  = V0_s_xi(:,ielem)
        v0_z_etal(:) = V0_z_eta(:,ielem)
        v0_z_xil(:)  = V0_z_xi(:,ielem)

        ! p - component
        V1 = v0_z_etal * (r1(0,:) - r6(0,:)) + v0_s_etal * (r5(0,:) - r4(0,:))
        
        V2 = v0_z_xil  * (r1(0,:) - r6(0,:)) + v0_s_xil  *  (r5(0,:) - r4(0,:))
        call vxm(V2, G2T, V3)
        
        loc_stiffness_p = loc_stiffness_p + outerprod(G0, V1)
        loc_stiffness_p(0,:) = loc_stiffness_p(0,:) + V3
        
        ! m - component
        V1 = v0_z_etal * (r1(0,:) + r6(0,:)) + v0_s_etal * (r5(0,:) + r4(0,:)) &
                + y0l * 2 * (r2(0,:) - r6(0,:))
        loc_stiffness_m = loc_stiffness_m + outerprod(G0, V1)

        ! z - component
        V1 = v0_z_etal * r5(0,:) + v0_s_etal * r3(0,:) - y0l * r4(0,:)
        loc_stiffness_z = loc_stiffness_z + outerprod(G0, V2)
     endif

     ! subtracting (!) from the global stiffness
     glob_stiffness(0:npol,0:npol,ielem,1) = &
            glob_stiffness(0:npol,0:npol,ielem,1) - loc_stiffness_p
     glob_stiffness(0:npol,0:npol,ielem,2) = &
            glob_stiffness(0:npol,0:npol,ielem,2) - loc_stiffness_m
     glob_stiffness(0:npol,0:npol,ielem,3) = &
            glob_stiffness(0:npol,0:npol,ielem,3) - loc_stiffness_z
  enddo

end subroutine glob_anel_stiffness_di_generic
!=============================================================================

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_di_4(glob_stiffness, R)

  use attenuation,  only: n_sls_attenuation
  use data_mesh,    only: nel_solid
  
  ! I/O global arrays
  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: R(0:,0:,:,:,:)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_p
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_m
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_z
  
  real(kind=realkind), dimension(0:4,0:4) :: r1, r2, r3, r4, r5, r6

  real(kind=realkind), dimension(0:4,0:4) :: yl
  real(kind=realkind), dimension(0:4,0:4) :: v_s_etal, v_s_xil
  real(kind=realkind), dimension(0:4,0:4) :: v_z_etal, v_z_xil

  real(kind=realkind), dimension(0:4,0:4) :: S1p, S2p
  real(kind=realkind), dimension(0:4,0:4) :: S1m, S2m
  real(kind=realkind), dimension(0:4,0:4) :: S1z, S2z
  real(kind=realkind), dimension(0:4,0:4) :: X1, X2, X3, X4, X5, X6
  
  real(kind=realkind), dimension(0:4) :: y0l
  real(kind=realkind), dimension(0:4) :: v0_s_etal, v0_s_xil
  real(kind=realkind), dimension(0:4) :: v0_z_etal, v0_z_xil
  real(kind=realkind), dimension(0:4) :: V1, V2, V3
  
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

     S1p = v_z_etal * (r1 - r6) + v_s_etal * (r5 - r4)
     S2p = v_z_xil  * (r1 - r6) + v_s_xil  * (r5 - r4)
     
     S1m = v_z_etal * (r1 + r6) + v_s_etal * (r5 + r4)
     S2m = v_z_xil  * (r1 + r6) + v_s_xil  * (r5 + r4)
     
     S1z = v_z_etal * r5 + v_s_etal * r3
     S2z = v_z_xil  * r5 + v_s_xil  * r3

     if ( .not. axis_solid(ielem) ) then
        call mxm_4(G2, S1p, X1)
        call mxm_4(G2, S1m, X3)
        call mxm_4(G2, S1z, X5)
     else
        call mxm_4(G1, S1p, X1)
        call mxm_4(G1, S1m, X3)
        call mxm_4(G1, S1z, X5)
     endif

     call mxm_4(S2p, G2T, X2)
     call mxm_4(S2m, G2T, X4)
     call mxm_4(S2z, G2T, X6)

     loc_stiffness_p = X1 + X2
     loc_stiffness_m = X3 + X4 + 2 * yl * (r2 - r6)
     loc_stiffness_z = X5 + X6 - yl * r4

     if (axis_solid(ielem)) then
        y0l(:) = Y0(:,ielem)
        v0_s_etal(:) = V0_s_eta(:,ielem)
        v0_s_xil(:)  = V0_s_xi(:,ielem)
        v0_z_etal(:) = V0_z_eta(:,ielem)
        v0_z_xil(:)  = V0_z_xi(:,ielem)

        ! p - component
        V1 = v0_z_etal * (r1(0,:) - r6(0,:)) + v0_s_etal * (r5(0,:) - r4(0,:))
        
        V2 = v0_z_xil  * (r1(0,:) - r6(0,:)) + v0_s_xil  *  (r5(0,:) - r4(0,:))
        call vxm_4(V2, G2T, V3)
        
        loc_stiffness_p = loc_stiffness_p + outerprod_4(G0, V1)
        loc_stiffness_p(0,:) = loc_stiffness_p(0,:) + V3
        
        ! m - component
        V1 = v0_z_etal * (r1(0,:) + r6(0,:)) + v0_s_etal * (r5(0,:) + r4(0,:)) &
                + y0l * 2 * (r2(0,:) - r6(0,:))
        loc_stiffness_m = loc_stiffness_m + outerprod_4(G0, V1)

        ! z - component
        V1 = v0_z_etal * r5(0,:) + v0_s_etal * r3(0,:) - y0l * r4(0,:)
        loc_stiffness_z = loc_stiffness_z + outerprod_4(G0, V2)
     endif

     ! subtracting (!) from the global stiffness
     glob_stiffness(0:4,0:4,ielem,1) = &
            glob_stiffness(0:4,0:4,ielem,1) - loc_stiffness_p
     glob_stiffness(0:4,0:4,ielem,2) = &
            glob_stiffness(0:4,0:4,ielem,2) - loc_stiffness_m
     glob_stiffness(0:4,0:4,ielem,3) = &
            glob_stiffness(0:4,0:4,ielem,3) - loc_stiffness_z
  enddo

end subroutine glob_anel_stiffness_di_4
!=============================================================================

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_di_cg4(glob_stiffness, R)

  use attenuation,  only: n_sls_attenuation
  use data_mesh,    only: npol, nel_solid
  
  ! I/O global arrays
  real(kind=realkind), intent(inout) :: glob_stiffness(0:,0:,:,:)
  real(kind=realkind), intent(in)    :: R(:,:,:,:) !(1:4,6,n_sls_attenuation,nel_solid)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_p
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_m
  real(kind=realkind), dimension(0:4,0:4) :: loc_stiffness_z
  
  real(kind=realkind), dimension(1:4) :: r1, r2, r3, r4, r5, r6

  real(kind=realkind), dimension(1:4) :: yl
  real(kind=realkind), dimension(1:4) :: v_s_etal, v_s_xil
  real(kind=realkind), dimension(1:4) :: v_z_etal, v_z_xil

  real(kind=realkind), dimension(1:4) :: S1p, S2p
  real(kind=realkind), dimension(1:4) :: S1m, S2m
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

     S1p = v_z_etal * (r1 - r6) + v_s_etal * (r5 - r4)
     S2p = v_z_xil  * (r1 - r6) + v_s_xil  * (r5 - r4)
     
     S1m = v_z_etal * (r1 + r6) + v_s_etal * (r5 + r4)
     S2m = v_z_xil  * (r1 + r6) + v_s_xil  * (r5 + r4)
     
     S1z = v_z_etal * r5 + v_s_etal * r3
     S2z = v_z_xil  * r5 + v_s_xil  * r3

     if ( .not. axis_solid(ielem) ) then
        call mxm_cg4_sparse_b(G2, S1p, X1)
        call mxm_cg4_sparse_b(G2, S1m, X3)
        call mxm_cg4_sparse_b(G2, S1z, X5)
     else
        call mxm_cg4_sparse_b(G1, S1p, X1)
        call mxm_cg4_sparse_b(G1, S1m, X3)
        call mxm_cg4_sparse_b(G1, S1z, X5)
     endif

     call mxm_cg4_sparse_a(S2p, G2T, X2)
     call mxm_cg4_sparse_a(S2m, G2T, X4)
     call mxm_cg4_sparse_a(S2z, G2T, X6)

     loc_stiffness_p = X1 + X2

     loc_stiffness_m = X3 + X4
     loc_stiffness_m(1,1) = loc_stiffness_m(1,1) + 2 * yl(1) * (r2(1) - r6(1))
     loc_stiffness_m(1,3) = loc_stiffness_m(1,3) + 2 * yl(2) * (r2(2) - r6(2))
     loc_stiffness_m(3,1) = loc_stiffness_m(3,1) + 2 * yl(3) * (r2(3) - r6(3))
     loc_stiffness_m(3,3) = loc_stiffness_m(3,3) + 2 * yl(4) * (r2(4) - r6(4))

     loc_stiffness_z = X5 + X6
     loc_stiffness_z(1,1) = loc_stiffness_z(1,1) - yl(1) * r4(1)
     loc_stiffness_z(1,3) = loc_stiffness_z(1,3) - yl(2) * r4(2)
     loc_stiffness_z(3,1) = loc_stiffness_z(3,1) - yl(3) * r4(3)
     loc_stiffness_z(3,3) = loc_stiffness_z(3,3) - yl(4) * r4(4)

     ! subtracting (!) from the global stiffness
     glob_stiffness(0:4,0:4,ielem,1) = &
            glob_stiffness(0:4,0:4,ielem,1) - loc_stiffness_p
     glob_stiffness(0:4,0:4,ielem,2) = &
            glob_stiffness(0:4,0:4,ielem,2) - loc_stiffness_m
     glob_stiffness(0:4,0:4,ielem,3) = &
            glob_stiffness(0:4,0:4,ielem,3) - loc_stiffness_z
  enddo

end subroutine glob_anel_stiffness_di_cg4
!=============================================================================

!-----------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
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
     if ( .not. axis_solid(ielem) ) then
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

     if ( .not. axis_solid(ielem) ) then
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
!=============================================================================

!-----------------------------------------------------------------------------
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
     if ( .not. axis_solid(ielem) ) then
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

     if ( .not. axis_solid(ielem) ) then
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
!=============================================================================

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_quad(glob_stiffness, R, R_cg, cg)
  use data_mesh,    only: npol
  
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
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_quad_generic(glob_stiffness, R)

  use attenuation,  only: n_sls_attenuation
  use data_mesh,    only: npol, nel_solid
  
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

     if ( .not. axis_solid(ielem) ) then
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
!=============================================================================

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_quad_4(glob_stiffness, R)

  use attenuation,  only: n_sls_attenuation
  use data_mesh,    only: nel_solid
  
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

     if ( .not. axis_solid(ielem) ) then
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
!=============================================================================

!-----------------------------------------------------------------------------
pure subroutine glob_anel_stiffness_quad_cg4(glob_stiffness, R)

  use attenuation,  only: n_sls_attenuation
  use data_mesh,    only: npol, nel_solid
  
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

     if ( .not. axis_solid(ielem) ) then
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
!=============================================================================

!-----------------------------------------------------------------------------
!> Wrapper routine to avoid if statements in the timeloop
pure subroutine glob_fluid_stiffness(glob_stiffness_fl, chi)
  use data_mesh, only: npol, nel_fluid
  
  real(kind=realkind), intent(in)  :: chi(0:,0:,:)
  real(kind=realkind), intent(out) :: glob_stiffness_fl(0:npol,0:npol,nel_fluid)

  if (npol == 4) then
     call glob_fluid_stiffness_4(glob_stiffness_fl, chi)
  else
     call glob_fluid_stiffness_generic(glob_stiffness_fl, chi)
  endif

end subroutine glob_fluid_stiffness
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
pure subroutine glob_fluid_stiffness_generic(glob_stiffness_fl, chi)

  use data_mesh, only: npol, nel_fluid
  
  ! I/O for global arrays
  real(kind=realkind), intent(in)  :: chi(0:,0:,:)
  real(kind=realkind), intent(out) :: glob_stiffness_fl(0:npol,0:npol,nel_fluid)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: chi_l, loc_stiffness
  real(kind=realkind), dimension(0:npol,0:npol) :: m_w_fl_l
  real(kind=realkind), dimension(0:npol,0:npol) :: m1chil, m2chil, m4chil
  
  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w_fl_l
  
  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2  ! MxM arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: S1, S2  ! Sum
  
  real(kind=realkind), dimension(0:npol) :: V1
  
  integer :: ielem

  do ielem = 1, nel_fluid

     chi_l(0:npol,0:npol) = chi(0:npol,0:npol,ielem)
     m1chil(0:npol,0:npol) = M1chi_fl(:,:,ielem)
     m2chil(0:npol,0:npol) = M2chi_fl(:,:,ielem)
     m4chil(0:npol,0:npol) = M4chi_fl(:,:,ielem)

     ! First MxM
     call mxm(G2T, chi_l, X1)
     call mxm(chi_l, G2, X2)

     ! Collocations and sums of D terms
     S1 = m1chil * X2 + m2chil * X1
     S2 = m1chil * X1 + m4chil * X2

     !Second MxM
     call mxm(G2, S1, X1)
     call mxm(S2, G2T, X2)
     
     ! Final Sum
     loc_stiffness = X1 + X2

     ! dipole and quadrupole cases: additional 2nd order term
     if (src_type(1) .ne. 'monopole') then

        m_w_fl_l(0:npol,0:npol) = M_w_fl(:,:,ielem)

        loc_stiffness = loc_stiffness + m_w_fl_l * chi_l 

        if ( axis_fluid(ielem) ) then
           m0_w_fl_l(0:npol) = M0_w_fl(0:npol,ielem)
           call vxm(G0,chi_l,V1)
           
           chi_l = outerprod(G0, m0_w_fl_l * V1) !chi_l as dummy 

           loc_stiffness = loc_stiffness + chi_l
        endif

     endif

     glob_stiffness_fl(0:npol,0:npol,ielem) = loc_stiffness

  enddo

end subroutine glob_fluid_stiffness_generic
!=============================================================================

!-----------------------------------------------------------------------------
pure subroutine glob_fluid_stiffness_4(glob_stiffness_fl, chi)

  use data_mesh, only: nel_fluid
  
  integer, parameter               :: npol = 4

  ! I/O for global arrays
  real(kind=realkind), intent(in)  :: chi(0:,0:,:)
  real(kind=realkind), intent(out) :: glob_stiffness_fl(0:npol,0:npol,nel_fluid)
  
  ! local variables for all elements
  real(kind=realkind), dimension(0:npol,0:npol) :: chi_l, loc_stiffness
  real(kind=realkind), dimension(0:npol,0:npol) :: m_w_fl_l
  real(kind=realkind), dimension(0:npol,0:npol) :: m1chil, m2chil, m4chil
  
  ! local variables for axial elements
  real(kind=realkind), dimension(0:npol) :: m0_w_fl_l
  
  ! work arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: X1, X2  ! MxM arrays
  real(kind=realkind), dimension(0:npol,0:npol) :: S1, S2  ! Sum
  
  real(kind=realkind), dimension(0:npol) :: V1
  
  integer :: ielem

  do ielem = 1, nel_fluid

     chi_l(0:npol,0:npol) = chi(0:npol,0:npol,ielem)
     m1chil(0:npol,0:npol) = M1chi_fl(:,:,ielem)
     m2chil(0:npol,0:npol) = M2chi_fl(:,:,ielem)
     m4chil(0:npol,0:npol) = M4chi_fl(:,:,ielem)

     ! First MxM
     call mxm_4(G2T, chi_l, X1)
     call mxm_4(chi_l, G2, X2)

     ! Collocations and sums of D terms
     S1 = m1chil * X2 + m2chil * X1
     S2 = m1chil * X1 + m4chil * X2

     !Second MxM
     call mxm_4(G2, S1, X1)
     call mxm_4(S2, G2T, X2)
     
     ! Final Sum
     loc_stiffness = X1 + X2

     ! dipole and quadrupole cases: additional 2nd order term
     if (src_type(1) .ne. 'monopole') then

        m_w_fl_l(0:npol,0:npol) = M_w_fl(:,:,ielem)

        loc_stiffness = loc_stiffness + m_w_fl_l * chi_l 

        if ( axis_fluid(ielem) ) then
           m0_w_fl_l(0:npol) = M0_w_fl(0:npol,ielem)
           call vxm_4(G0,chi_l,V1)
           
           chi_l = outerprod(G0, m0_w_fl_l * V1) !chi_l as dummy 

           loc_stiffness = loc_stiffness + chi_l
        endif

     endif

     glob_stiffness_fl(0:npol,0:npol,ielem) = loc_stiffness

  enddo

end subroutine glob_fluid_stiffness_4
!=============================================================================

!====================
end module stiffness
!====================
