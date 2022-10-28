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
module stiffness_fluid

  use global_parameters, only: realkind
  use data_matr
  use data_mesh, only: axis_solid, axis_fluid, nsize
  use data_spec
  use data_source

  use unrolled_loops

  implicit none

  public :: glob_fluid_stiffness

  private

contains

!-----------------------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
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
     if ( axis_fluid(ielem) ) then
        call mxm(G1T, chi_l, X1)
     else
        call mxm(G2T, chi_l, X1)
     endif
     call mxm(chi_l, G2, X2)

     ! Collocations and sums of D terms
     S1 = m1chil * X2 + m2chil * X1
     S2 = m1chil * X1 + m4chil * X2

     !Second MxM
     if ( axis_fluid(ielem) ) then
        call mxm(G1, S1, X1)
     else
        call mxm(G2, S1, X1)
     endif
     call mxm(S2, G2T, X2)

     ! Final Sum
     loc_stiffness = X1 + X2

     ! dipole and quadrupole cases: additional 2nd order term
     if (src_type(1) /= 'monopole') then

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
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
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
     if ( axis_fluid(ielem) ) then
        call mxm_4(G1T, chi_l, X1)
     else
        call mxm_4(G2T, chi_l, X1)
     endif
     call mxm_4(chi_l, G2, X2)

     ! Collocations and sums of D terms
     S1 = m1chil * X2 + m2chil * X1
     S2 = m1chil * X1 + m4chil * X2

     !Second MxM
     if ( axis_fluid(ielem) ) then
        call mxm_4(G1, S1, X1)
     else
        call mxm_4(G2, S1, X1)
     endif
     call mxm_4(S2, G2T, X2)

     ! Final Sum
     loc_stiffness = X1 + X2

     ! dipole and quadrupole cases: additional 2nd order term
     if (src_type(1) /= 'monopole') then

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
!-----------------------------------------------------------------------------------------

end module stiffness_fluid
!=========================================================================================
