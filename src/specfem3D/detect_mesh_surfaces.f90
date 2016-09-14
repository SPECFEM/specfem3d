!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine detect_mesh_surfaces()

  use specfem_par
  use specfem_par_movie
  use specfem_par_acoustic
  use specfem_par_elastic
  implicit none
  integer :: ier

  ! for mesh surface
  allocate(ispec_is_surface_external_mesh(NSPEC_AB), &
          iglob_is_surface_external_mesh(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'error allocating array for mesh surface'

  ! determines model surface
  if (.not. RECEIVERS_CAN_BE_BURIED .or. MOVIE_SURFACE .or. CREATE_SHAKEMAP) then
    ! returns surface points/elements
    ! in ispec_is_surface_external_mesh / iglob_is_surface_external_mesh and
    ! number of faces in nfaces_surface
    call detect_surface(NPROC,NGLOB_AB,NSPEC_AB,ibool, &
                      ispec_is_surface_external_mesh, &
                      iglob_is_surface_external_mesh, &
                      nfaces_surface, &
                      num_interfaces_ext_mesh, &
                      max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh, &
                      my_neighbours_ext_mesh, &
                      ibool_interfaces_ext_mesh)
  endif

  ! takes cross-section surfaces instead
  if (MOVIE_SURFACE .or. CREATE_SHAKEMAP) then
    if (MOVIE_TYPE == 2 .and. PLOT_CROSS_SECTIONS) then
      call detect_surface_cross_section(NPROC,NGLOB_AB,NSPEC_AB,ibool, &
                              ispec_is_surface_external_mesh, &
                              iglob_is_surface_external_mesh, &
                              nfaces_surface, &
                              num_interfaces_ext_mesh, &
                              max_nibool_interfaces_ext_mesh, &
                              nibool_interfaces_ext_mesh, &
                              my_neighbours_ext_mesh, &
                              ibool_interfaces_ext_mesh, &
                              CROSS_SECTION_X,CROSS_SECTION_Y,CROSS_SECTION_Z, &
                              xstore,ystore,zstore,myrank)
    endif
  endif

  ! takes number of faces for top, free surface only
  if (MOVIE_TYPE == 1) then
    nfaces_surface = num_free_surface_faces
  endif

  ! handles movies and shakemaps
  if (MOVIE_SURFACE .or. CREATE_SHAKEMAP) then
    call setup_movie_meshes()
  endif

  ! stores wavefields for whole volume
  if (MOVIE_VOLUME) then
    ! acoustic
    if (ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION) then
      allocate(velocity_x(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
               velocity_y(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
               velocity_z(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if (ier /= 0) stop 'error allocating array movie velocity_x etc.'
    endif
    ! elastic only
    if (ELASTIC_SIMULATION) then
      allocate(div(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
               curl_x(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
               curl_y(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
               curl_z(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if (ier /= 0) stop 'error allocating array movie div and curl'
      div(:,:,:,:) = 0._CUSTOM_REAL
      curl_x(:,:,:,:) = 0._CUSTOM_REAL
      curl_y(:,:,:,:) = 0._CUSTOM_REAL
      curl_z(:,:,:,:) = 0._CUSTOM_REAL
    endif
  endif

  ! initializes cross-section gif image
  if (PNM_IMAGE) then
    call write_PNM_initialize()
  endif

  end subroutine detect_mesh_surfaces

