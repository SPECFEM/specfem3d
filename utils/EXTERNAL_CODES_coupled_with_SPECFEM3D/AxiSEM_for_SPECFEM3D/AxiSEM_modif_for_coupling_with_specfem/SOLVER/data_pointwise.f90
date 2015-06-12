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

!> This module is only known during the time loop if the strain tensor 
!! is computed on-the-fly. The fluid section is additionally known if global 
!! snapshots are dumped (to compute the displacement in the fluid).
!===================
 module data_pointwise
!===================
  
  use global_parameters, only: realkind
  
  implicit none
  public 
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++
  !  Precomputed matrices for pointwise derivatives
  !+++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !> Inverse density in fluid - only needed when computing pointwise displacement
  real(kind=realkind), allocatable :: inv_rho_fluid(:,:,:)

  !> (s rho)^-1 in fluid - only for phi comp. of fluid displacement
  real(kind=realkind), allocatable :: prefac_inv_s_rho_fluid(:,:,:)
  real(kind=realkind), allocatable :: inv_s_fluid(:,:,:)

  !> (s)^-1 in solid - needed for the strain tensor, if computed on-the-fly
  real(kind=realkind), allocatable :: inv_s_solid(:,:,:)

  ! Note that these matrices may include minus signs where necessary.

  ! Solid region only
  real(kind=realkind), allocatable :: DsDeta_over_J_sol(:,:,:)
  real(kind=realkind), allocatable :: DzDeta_over_J_sol(:,:,:)
  real(kind=realkind), allocatable :: DsDxi_over_J_sol(:,:,:)
  real(kind=realkind), allocatable :: DzDxi_over_J_sol(:,:,:)

  real(kind=realkind), allocatable :: DsDeta_over_J_sol_cg4(:,:)
  real(kind=realkind), allocatable :: DzDeta_over_J_sol_cg4(:,:)
  real(kind=realkind), allocatable :: DsDxi_over_J_sol_cg4(:,:)
  real(kind=realkind), allocatable :: DzDxi_over_J_sol_cg4(:,:)

  ! Fluid region only
  real(kind=realkind), allocatable :: DsDeta_over_J_flu(:,:,:)
  real(kind=realkind), allocatable :: DzDeta_over_J_flu(:,:,:)
  real(kind=realkind), allocatable :: DsDxi_over_J_flu(:,:,:)
  real(kind=realkind), allocatable :: DzDxi_over_J_flu(:,:,:)

!=======================
 end module data_pointwise
!=======================
