!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

! for elastic solver

  subroutine compute_coupling_viscoelastic_ac(NSPEC_AB,NGLOB_AB, &
                                              ibool,accel,potential_dot_dot_acoustic, &
                                              num_coupling_ac_el_faces, &
                                              coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                              coupling_ac_el_normal, &
                                              coupling_ac_el_jacobian2Dw, &
                                              iphase, &
                                              PML_CONDITIONS, &
                                              SIMULATION_TYPE,backward_simulation, &
                                              potential_acoustic,potential_dot_acoustic)

! returns the updated acceleration array: accel

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,NGLLSQUARE
  use pml_par, only: rmemory_coupling_el_ac_potential_dot_dot,is_CPML,spec_to_CPML,NSPEC_CPML

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB,SIMULATION_TYPE
  logical,intent(in) :: backward_simulation

  ! displacement and pressure
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic

  ! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! acoustic-elastic coupling surface
  integer,intent(in) :: num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL),intent(in) :: coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces)
  real(kind=CUSTOM_REAL),intent(in) :: coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces)
  integer,intent(in) :: coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces)
  integer,intent(in) :: coupling_ac_el_ispec(num_coupling_ac_el_faces)

  ! communication overlap
  integer,intent(in) :: iphase

  ! local parameters
  real(kind=CUSTOM_REAL) :: pressure_x,pressure_y,pressure_z
  real(kind=CUSTOM_REAL) :: nx,ny,nz,jacobianw

  integer :: iface,igll,ispec,iglob
  integer :: i,j,k

  ! CPML
  integer :: ispec_CPML
  logical :: PML_CONDITIONS

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! loops on all coupling faces
  do iface = 1,num_coupling_ac_el_faces

    ! gets corresponding spectral element
    ! (note: can be either acoustic or elastic element, no need to specify since
    !           no material properties are needed for this coupling term)
    ispec = coupling_ac_el_ispec(iface)

    ! loops over common GLL points
    do igll = 1, NGLLSQUARE
      i = coupling_ac_el_ijk(1,igll,iface)
      j = coupling_ac_el_ijk(2,igll,iface)
      k = coupling_ac_el_ijk(3,igll,iface)

      ! gets global index of this common GLL point
      ! (note: should be the same as for corresponding i',j',k',ispec_elastic or ispec_elastic )
      iglob = ibool(i,j,k,ispec)

      ! acoustic pressure on global point
      pressure_x = - potential_dot_dot_acoustic(iglob)
      pressure_y = pressure_x
      pressure_z = pressure_x

      ! adjoint wavefield case
      if (SIMULATION_TYPE /= 1 .and. (.not. backward_simulation)) then
        ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
        ! adjoint definition: pressure^\dagger = potential^\dagger
        pressure_x = - pressure_x
        pressure_y = - pressure_y
        pressure_z = - pressure_z
      endif

      ! CPML overwrite cases
      if (PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
        if (is_CPML(ispec)) then
          if (SIMULATION_TYPE == 1) then
            ispec_CPML = spec_to_CPML(ispec)
            call pml_compute_memory_variables_elastic_acoustic(ispec_CPML,iface,iglob,i,j,k, &
                                            pressure_x,pressure_y,pressure_z,potential_acoustic, &
                                            potential_dot_acoustic,potential_dot_dot_acoustic, &
                                            num_coupling_ac_el_faces,rmemory_coupling_el_ac_potential_dot_dot)
            pressure_x = - pressure_x
            pressure_y = - pressure_y
            pressure_z = - pressure_z
          endif

          if (SIMULATION_TYPE == 3) then
            ispec_CPML = spec_to_CPML(ispec)
            call pml_compute_memory_variables_elastic_acoustic(ispec_CPML,iface,iglob,i,j,k, &
                                            pressure_x,pressure_y,pressure_z,potential_acoustic, &
                                            potential_dot_acoustic,potential_dot_dot_acoustic, &
                                            num_coupling_ac_el_faces,rmemory_coupling_el_ac_potential_dot_dot)
          endif
        endif
      endif

      ! gets associated normal on GLL point
      ! (note convention: pointing outwards of acoustic element)
      nx = coupling_ac_el_normal(1,igll,iface)
      ny = coupling_ac_el_normal(2,igll,iface)
      nz = coupling_ac_el_normal(3,igll,iface)

      ! gets associated, weighted 2D jacobian
      ! (note: should be the same for elastic and acoustic element)
      jacobianw = coupling_ac_el_jacobian2Dw(igll,iface)

      ! continuity of displacement and pressure on global point
      !
      ! note: Newmark time scheme together with definition of scalar potential:
      !          pressure = - chi_dot_dot
      !          requires that this coupling term uses the *UPDATED* pressure (chi_dot_dot), i.e.
      !          pressure at time step [t + delta_t]
      !          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
      !          it means you have to calculate and update the acoustic pressure first before
      !          calculating this term...
      accel(1,iglob) = accel(1,iglob) + jacobianw * nx * pressure_x
      accel(2,iglob) = accel(2,iglob) + jacobianw * ny * pressure_y
      accel(3,iglob) = accel(3,iglob) + jacobianw * nz * pressure_z

    enddo ! igll

  enddo ! iface

  end subroutine compute_coupling_viscoelastic_ac

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_ocean(NSPEC_AB,NGLOB_AB, &
                                    ibool,rmassx,rmassy,rmassz, &
                                    rmass_ocean_load,accel, &
                                    free_surface_normal,free_surface_ijk,free_surface_ispec, &
                                    num_free_surface_faces)

! updates acceleration with ocean load term:
! approximates ocean-bottom continuity of pressure & displacement for longer period waves (> ~20s ),
! assuming incompressible fluid column above bathymetry ocean bottom

  use constants

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB),intent(inout) :: accel
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: rmassx,rmassy,rmassz
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: rmass_ocean_load

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! free surface
  integer,intent(in) :: num_free_surface_faces
  real(kind=CUSTOM_REAL),intent(in) :: free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces)
  integer,intent(in) :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer,intent(in) :: free_surface_ispec(num_free_surface_faces)

! local parameters
  real(kind=CUSTOM_REAL) :: nx,ny,nz
  real(kind=CUSTOM_REAL) :: force_normal_comp
  integer :: i,j,k,ispec,iglob
  integer :: igll,iface
  logical,dimension(NGLOB_AB) :: updated_dof_ocean_load

  !   initialize the updates
  updated_dof_ocean_load(:) = .false.

  ! for surface elements exactly at the top of the model (ocean bottom)
  do iface = 1,num_free_surface_faces

    ispec = free_surface_ispec(iface)
    do igll = 1, NGLLSQUARE
      i = free_surface_ijk(1,igll,iface)
      j = free_surface_ijk(2,igll,iface)
      k = free_surface_ijk(3,igll,iface)

      ! get global point number
      iglob = ibool(i,j,k,ispec)

      ! only update once
      if (.not. updated_dof_ocean_load(iglob)) then

        ! get normal
        nx = free_surface_normal(1,igll,iface)
        ny = free_surface_normal(2,igll,iface)
        nz = free_surface_normal(3,igll,iface)

        ! make updated component of right-hand side
        ! we divide by rmass() which is 1 / M
        ! we use the total force which includes the Coriolis term above
        force_normal_comp = accel(1,iglob)*nx / rmassx(iglob) &
                            + accel(2,iglob)*ny / rmassy(iglob) &
                            + accel(3,iglob)*nz / rmassz(iglob)

        accel(1,iglob) = accel(1,iglob) &
          + (rmass_ocean_load(iglob) - rmassx(iglob)) * force_normal_comp * nx
        accel(2,iglob) = accel(2,iglob) &
          + (rmass_ocean_load(iglob) - rmassy(iglob)) * force_normal_comp * ny
        accel(3,iglob) = accel(3,iglob) &
          + (rmass_ocean_load(iglob) - rmassz(iglob)) * force_normal_comp * nz

        ! done with this point
        updated_dof_ocean_load(iglob) = .true.

      endif

    enddo ! igll
  enddo ! iface

  end subroutine compute_coupling_ocean
!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_ocean_backward(NSPEC_AB,NGLOB_AB, &
                                    ibool,rmassx,rmassy,rmassz,rmass_ocean_load, &
                                    free_surface_normal,free_surface_ijk,free_surface_ispec, &
                                    num_free_surface_faces,SIMULATION_TYPE, &
                                    NGLOB_ADJOINT,b_accel)

! updates acceleration with ocean load term:
! approximates ocean-bottom continuity of pressure & displacement for longer period waves (> ~20s ),
! assuming incompressible fluid column above bathymetry ocean bottom

  use constants

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: rmassx,rmassy,rmassz
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: rmass_ocean_load

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! free surface
  integer :: num_free_surface_faces
  real(kind=CUSTOM_REAL) :: free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ispec(num_free_surface_faces)

  ! adjoint simulations
  integer :: SIMULATION_TYPE,NGLOB_ADJOINT
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_ADJOINT):: b_accel

! local parameters
  real(kind=CUSTOM_REAL) :: nx,ny,nz
  integer :: i,j,k,ispec,iglob
  integer :: igll,iface
  logical,dimension(NGLOB_AB) :: updated_dof_ocean_load
  ! adjoint locals
  real(kind=CUSTOM_REAL) :: b_force_normal_comp

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return

  !   initialize the updates
  updated_dof_ocean_load(:) = .false.

  ! for surface elements exactly at the top of the model (ocean bottom)
  do iface = 1,num_free_surface_faces

    ispec = free_surface_ispec(iface)
    do igll = 1, NGLLSQUARE
      i = free_surface_ijk(1,igll,iface)
      j = free_surface_ijk(2,igll,iface)
      k = free_surface_ijk(3,igll,iface)

      ! get global point number
      iglob = ibool(i,j,k,ispec)

      ! only update once
      if (.not. updated_dof_ocean_load(iglob)) then

        ! get normal
        nx = free_surface_normal(1,igll,iface)
        ny = free_surface_normal(2,igll,iface)
        nz = free_surface_normal(3,igll,iface)

        ! make updated component of right-hand side
        ! we divide by rmass() which is 1 / M
        ! we use the total force which includes the Coriolis term above

        ! adjoint simulations
        b_force_normal_comp = b_accel(1,iglob)*nx / rmassx(iglob) &
                              + b_accel(2,iglob)*ny / rmassy(iglob) &
                              + b_accel(3,iglob)*nz / rmassz(iglob)

        b_accel(1,iglob) = b_accel(1,iglob) &
          + (rmass_ocean_load(iglob) - rmassx(iglob)) * b_force_normal_comp * nx
        b_accel(2,iglob) = b_accel(2,iglob) &
          + (rmass_ocean_load(iglob) - rmassy(iglob)) * b_force_normal_comp * ny
        b_accel(3,iglob) = b_accel(3,iglob) &
          + (rmass_ocean_load(iglob) - rmassz(iglob)) * b_force_normal_comp * nz

        ! done with this point
        updated_dof_ocean_load(iglob) = .true.

      endif

    enddo ! igll
  enddo ! iface

  end subroutine compute_coupling_ocean_backward


