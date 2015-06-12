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

!===================
module data_matr
!===================
  !
  ! Global arrays (i.e. defined on each GLL point) that are
  ! needed for the mass, stiffness and boundary terms of the 
  ! temporal ODE. 
  
  use global_parameters

  implicit none
  public 

  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  !	Mass matrix arrays
  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  real(kind=realkind), protected, dimension(:,:,:), allocatable :: inv_mass_rho
  real(kind=realkind), protected, dimension(:,:,:), allocatable :: inv_mass_fluid
  real(kind=realkind), dimension(:,:,:), allocatable :: unassem_mass_rho_solid
  real(kind=realkind), dimension(:,:,:), allocatable :: unassem_mass_lam_fluid 

  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Precomputed stiffness matrices
  !++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Static solid matrices for all source types:
  real(kind=realkind), dimension(:,:,:), allocatable :: M11s, M21s
  real(kind=realkind), dimension(:,:,:), allocatable :: M41s, M12s
  real(kind=realkind), dimension(:,:,:), allocatable :: M22s, M42s
  real(kind=realkind), dimension(:,:,:), allocatable :: M11z, M21z, M41z
  real(kind=realkind), dimension(:,:,:), allocatable :: M32s
  ! for dipole: 
  real(kind=realkind), dimension(:,:,:), allocatable :: M13s, M33s, M43s
  ! for quadrupole: 
  real(kind=realkind), dimension(:,:,:), allocatable :: M1phi, M2phi, M4phi
  
  
  ! for all source types:
  real(kind=realkind), dimension(:,:,:), allocatable :: M_1, M_2, M_3, M_4
  ! for dipole and quadpole:
  real(kind=realkind), allocatable :: M_5(:,:,:), M_6(:,:,:), M_7(:,:,:), M_8(:,:,:) 
  
  
  ! for all source types:
  real(kind=realkind), dimension(:,:,:), allocatable :: M_w1
  ! for dipole and quadpole:
  real(kind=realkind), allocatable :: M_w2(:,:,:), M_w3(:,:,:)
  ! for quadrupole: 
  real(kind=realkind), allocatable :: M_w4(:,:,:), M_w5(:,:,:)
  

  ! for all source types:
  real(kind=realkind), dimension(:,:), allocatable :: M0_w1, M0_w2, M0_w3
  ! for dipole and quadpole:
  real(kind=realkind), allocatable :: M0_w4(:,:), M0_w5(:,:), M0_w6(:,:)
  ! for dipole
  real(kind=realkind), allocatable :: M0_w7(:,:), M0_w8(:,:), M0_w9(:,:), M0_w10(:,:)

  ! Fluid matrices
  real(kind=realkind), dimension(:,:,:), allocatable :: M1chi_fl, M2chi_fl
  real(kind=realkind), dimension(:,:,:), allocatable :: M4chi_fl
  real(kind=realkind), dimension(:,:,:), allocatable :: M_w_fl
  real(kind=realkind), dimension(:,:), allocatable   :: M0_w_fl

  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Precomputed solid-fluid boundary matrices
  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Solid-fluid boundary matrix
  real(kind=realkind), dimension(:,:,:), allocatable :: bdry_matr
  real(kind=realkind), dimension(:,:,:), allocatable :: bdry_matr_fluid
  real(kind=realkind), dimension(:,:,:), allocatable :: bdry_matr_solid
  real(kind=dp), dimension(:)          , allocatable :: solflubdry_radius
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  !	Attenuation
  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Q_mu and Q_kappa, assumed to be homogeneous within an element
  real(kind=realkind), allocatable :: Q_mu(:), Q_kappa(:)
  real(kind=sp), allocatable       :: points_solid(:,:,:,:) ! for memory variable output
  real(kind=realkind), allocatable :: delta_mu(:,:,:), delta_kappa(:,:,:)
  real(kind=realkind), allocatable :: delta_mu_cg4(:,:), delta_kappa_cg4(:,:)

  ! Anelastic precomputable matrices 
  real(kind=realkind), allocatable :: Y(:,:,:)
  real(kind=realkind), allocatable :: Y_cg4(:,:)
  real(kind=realkind), allocatable :: Y0(:,:)

  real(kind=realkind), allocatable :: V_s_eta(:,:,:), V_s_xi(:,:,:)
  real(kind=realkind), allocatable :: V_z_eta(:,:,:), V_z_xi(:,:,:)

  real(kind=realkind), allocatable :: V_s_eta_cg4(:,:), V_s_xi_cg4(:,:)
  real(kind=realkind), allocatable :: V_z_eta_cg4(:,:), V_z_xi_cg4(:,:)

  real(kind=realkind), allocatable :: V0_s_eta(:,:), V0_s_xi(:,:)
  real(kind=realkind), allocatable :: V0_z_eta(:,:), V0_z_xi(:,:)

contains

subroutine set_mass_matrices(npol, nel_solid, nel_fluid, inv_mass_rho_loc, inv_mass_fluid_loc)

  integer, intent(in)        :: npol, nel_solid, nel_fluid
  real(kind=realkind), intent(in)  :: inv_mass_rho_loc(:,:,:), inv_mass_fluid_loc(:,:,:)

  allocate(inv_mass_rho(0:npol,0:npol,1:nel_solid))
  allocate(inv_mass_fluid(0:npol,0:npol,1:nel_fluid))

  inv_mass_rho   = inv_mass_rho_loc
  inv_mass_fluid = inv_mass_fluid_loc
  
end subroutine set_mass_matrices


!=====================
end module data_matr
!=====================
