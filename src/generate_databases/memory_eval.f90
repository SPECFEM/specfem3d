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


! compute the approximate amount of memory needed to run the solver

 subroutine memory_eval(NSPEC_AB,NGLOB_AB,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh, &
                        APPROXIMATE_OCEAN_LOAD,memory_size)

  use constants
  use shared_parameters, only: ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION
  use generate_databases_par, only: PML_CONDITIONS,nspec_cpml
  use create_regions_mesh_ext_par, only: NSPEC_ANISO,NSPEC_PORO,nspec_irregular

  implicit none

  ! input
  integer, intent(in) :: NSPEC_AB,NGLOB_AB
  integer, intent(in) :: max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh
  logical, intent(in) :: APPROXIMATE_OCEAN_LOAD
  ! output
  double precision, intent(out) :: memory_size

  memory_size = 0.d0

! add size of each set of arrays multiplied by the number of such arrays

  ! see: initialize_simulation.f90
  ! ibool
  memory_size = memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_AB*dble(SIZE_INTEGER)

  ! xix,xiy,xiz,
  ! etax,etay,etaz,
  ! gammax,gammay,gammaz,jacobian
  memory_size = memory_size + 10.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*nspec_irregular*dble(CUSTOM_REAL)

  ! irregular_element_number
  memory_size = memory_size + NSPEC_AB*dble(SIZE_INTEGER)

  ! xstore,ystore,zstore
  memory_size = memory_size + 3.d0*NGLOB_AB*dble(CUSTOM_REAL)

  ! kappastore,mustore
  memory_size = memory_size + 2.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_AB*dble(CUSTOM_REAL)

  ! ispec_acoustic,ispec_elastic,ispec_is_poroelastic (logical)
  memory_size = memory_size + 3.d0*NSPEC_AB*dble(SIZE_LOGICAL)

  ! see: read_mesh_databases.f90
  ! acoustic arrays
  if (ACOUSTIC_SIMULATION) then
    ! potential_acoustic, potentical_dot_acoustic, potential_dot_dot_acoustic
    memory_size = memory_size + 3.d0*NGLOB_AB*dble(CUSTOM_REAL)
    ! rmass_acoustic
    memory_size = memory_size + NGLOB_AB*dble(CUSTOM_REAL)
    ! rhostore
    memory_size = memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_AB*dble(CUSTOM_REAL)
  endif

  ! see: read_mesh_databases.f90 and pml_allocate_arrays.f90
  ! C-PML arrays
  ! is_CPML
  memory_size = memory_size + NSPEC_AB*dble(SIZE_LOGICAL)
  if (PML_CONDITIONS) then
     ! CPML_regions,CPML_to_spec,CPML_type
     memory_size = memory_size + 3.d0*nspec_cpml*dble(SIZE_INTEGER)

     ! spec_to_CPML
     memory_size = memory_size + NSPEC_AB*dble(SIZE_INTEGER)

     ! d_store_x,d_store_y,d_store_z,d_store_x,d_store_y,d_store_z,alpha_store_x,alpha_store_y,alpha_store_z
     memory_size = memory_size + 9.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*nspec_cpml*dble(CUSTOM_REAL)

     ! rmemory_dux_dxl_x,rmemory_dux_dyl_x,rmemory_dux_dzl_x,rmemory_duy_dxl_x,
     ! rmemory_duy_dyl_x,rmemory_duz_dxl_x,rmemory_duz_dzl_x,
     ! rmemory_dux_dxl_y,rmemory_dux_dyl_y,rmemory_duy_dxl_y,rmemory_duy_dyl_y,
     ! rmemory_duy_dzl_y,rmemory_duz_dyl_y,rmemory_duz_dzl_y,
     ! rmemory_dux_dxl_z,rmemory_dux_dzl_z,rmemory_duy_dyl_z,rmemory_duy_dzl_z,
     ! rmemory_duz_dxl_z,rmemory_duz_dyl_z,rmemory_duz_dzl_z
     memory_size = memory_size + 21.d0*3.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*nspec_cpml*dble(CUSTOM_REAL)

     ! rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl
     memory_size = memory_size + 3.d0*3.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*nspec_cpml*dble(CUSTOM_REAL)

     ! rmemory_displ_elastic
     memory_size = memory_size + 3.d0*dble(NDIM)*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*nspec_cpml*dble(CUSTOM_REAL)

     ! rmemory_potential_acoustic
     memory_size = memory_size + 3.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*nspec_cpml*dble(CUSTOM_REAL)

     ! second derivative of the potential
     memory_size = memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*nspec_cpml*dble(CUSTOM_REAL)
  endif

  ! elastic arrays
  if (ELASTIC_SIMULATION) then
    ! displacement,velocity,acceleration
    memory_size = memory_size + 3.d0*dble(NDIM)*NGLOB_AB*dble(CUSTOM_REAL)

    ! rmass
    memory_size = memory_size + 3*NGLOB_AB*dble(CUSTOM_REAL)

    ! rho_vp,rho_vs
    memory_size = memory_size + 2.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_AB*dble(CUSTOM_REAL)

    ! qmu_attenaution_store
    memory_size = memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_AB*dble(CUSTOM_REAL)

    ! c11store,...c66store
    memory_size = memory_size + 21.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_ANISO*dble(CUSTOM_REAL)

    if (APPROXIMATE_OCEAN_LOAD) then
      ! rmass_ocean_load
      memory_size = memory_size + NGLOB_AB*dble(CUSTOM_REAL)
      ! updated_dof_ocean_load
      memory_size = memory_size + NGLOB_AB*dble(SIZE_LOGICAL)
    endif
  endif

  ! elastic arrays
  if (POROELASTIC_SIMULATION) then
    ! displs_poroelastic,..
    memory_size = memory_size + 6.d0*dble(NDIM)*NGLOB_AB*dble(CUSTOM_REAL)
    ! rmass_solid_poroelastic,..
    memory_size = memory_size + 2.d0*NGLOB_AB*dble(CUSTOM_REAL)
    ! rhoarraystore,..
    memory_size = memory_size + 17.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_PORO*dble(CUSTOM_REAL)
  endif

  ! skipping boundary surfaces
  ! skipping free surfaces
  ! skipping acoustic-elastic coupling surfaces

  ! MPI interfaces
  ! my_neighbors_ext_mesh,nibool_interfaces_ext_mesh
  memory_size = memory_size + 2.d0*num_interfaces_ext_mesh*dble(SIZE_INTEGER)

  ! ibool_interfaces_ext_mesh
  memory_size = memory_size + max_nibool_interfaces_ext_mesh*num_interfaces_ext_mesh*dble(SIZE_INTEGER)

  ! MPI communications
  if (ACOUSTIC_SIMULATION) then
    ! buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh
    memory_size = memory_size + 2.d0*max_nibool_interfaces_ext_mesh*num_interfaces_ext_mesh*dble(CUSTOM_REAL)
    ! request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh
    memory_size = memory_size + 2.d0*num_interfaces_ext_mesh*dble(SIZE_INTEGER)
  endif

  if (ELASTIC_SIMULATION) then
    ! buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh
    memory_size = memory_size + 2.d0*dble(NDIM)*max_nibool_interfaces_ext_mesh*num_interfaces_ext_mesh*dble(CUSTOM_REAL)
    ! request_send_vector_ext_mesh,request_recv_vector_ext_mesh
    memory_size = memory_size + 2.d0*num_interfaces_ext_mesh*dble(SIZE_INTEGER)
  endif

  if (POROELASTIC_SIMULATION) then
    ! buffer_send_vector_ext_mesh_s,..
    memory_size = memory_size + 4.d0*dble(NDIM)*max_nibool_interfaces_ext_mesh*num_interfaces_ext_mesh*dble(CUSTOM_REAL)
    ! request_send_vector_ext_mesh_s,..
    memory_size = memory_size + 4.d0*num_interfaces_ext_mesh*dble(SIZE_INTEGER)
  endif

  ! ispec_is_inner
  memory_size = memory_size + NSPEC_AB*dble(SIZE_LOGICAL)

  ! skipping phase_ispec_inner_acoustic
  ! skipping phase_ispec_inner_elastic

  ! see: prepare_timerun.f90
  ! skipping attenuation R_xx,..R_yz and epsilondev_xx,...epsilondev_yz  no information yet about NSPEC_ATTENUATION_AB

  ! note: no adjoint array evaluation, since it depends on SIMULATION_TYPE which can vary for each run
  !         and is undependant of mesh databases

  ! note: no dyamic arrays like for seismograms, receivers and sources,
  !         since it depends on number of timesteps and number of stations etc., which can also vary for each run

  end subroutine memory_eval

!
!-------------------------------------------------------------------------------------------------
!

! compute the approximate amount of memory needed to run the mesher

 subroutine memory_eval_mesher(nspec,npointot,nnodes_ext_mesh, &
                               nelmnts_ext_mesh,nmat_ext_mesh,num_interfaces_ext_mesh, &
                               max_interface_size_ext_mesh,nspec2D_xmin,nspec2D_xmax, &
                               nspec2D_ymin,nspec2D_ymax,nspec2D_bottom,nspec2D_top, &
                               memory_size_request)

  use constants
  use generate_databases_par, only: NGNOD,NGNOD2D

  implicit none

  integer,intent(in) :: nspec,npointot,nnodes_ext_mesh,nelmnts_ext_mesh, &
           nmat_ext_mesh,num_interfaces_ext_mesh, &
           max_interface_size_ext_mesh,nspec2D_xmin,nspec2D_xmax, &
           nspec2D_ymin,nspec2D_ymax,nspec2D_bottom,nspec2D_top

  double precision,intent(inout) :: memory_size_request

  ! local parameters
  double precision, parameter :: BYTES_TO_MB = 1. / 1024. / 1024.
  double precision :: memory_size

  ! memory usage, in generate_database() routine so far
  ! due to overflow issues for larger meshes, we divide bytes to get MB
  memory_size = NGLLX*NGLLY*NGLLZ*nspec*dble(4) * BYTES_TO_MB
  memory_size = memory_size + 3*NGLLX*NGLLY*NGLLZ*nspec*dble(8) * BYTES_TO_MB
  memory_size = memory_size + NDIM*nnodes_ext_mesh*dble(8) * BYTES_TO_MB
  memory_size = memory_size + NGNOD*nelmnts_ext_mesh*dble(4) * BYTES_TO_MB
  memory_size = memory_size + 2*nelmnts_ext_mesh*dble(4) * BYTES_TO_MB
  memory_size = memory_size + 5*nmat_ext_mesh*dble(8) * BYTES_TO_MB
  memory_size = memory_size + 3*num_interfaces_ext_mesh*dble(1) * BYTES_TO_MB
  memory_size = memory_size + 6*max_interface_size_ext_mesh*num_interfaces_ext_mesh*dble(4) * BYTES_TO_MB
  memory_size = memory_size + NGLLX*NGLLX*max_interface_size_ext_mesh*num_interfaces_ext_mesh*dble(4) * BYTES_TO_MB
  memory_size = memory_size + nspec2D_xmin*20*dble(1) * BYTES_TO_MB
  memory_size = memory_size + nspec2D_xmax*20*dble(1) * BYTES_TO_MB
  memory_size = memory_size + nspec2D_ymin*20*dble(1) * BYTES_TO_MB
  memory_size = memory_size + nspec2D_ymax*20*dble(1) * BYTES_TO_MB
  memory_size = memory_size + nspec2D_bottom*20*dble(1) * BYTES_TO_MB
  memory_size = memory_size + nspec2D_top*20*dble(1) * BYTES_TO_MB

  ! memory usage, in create_regions_mesh_ext() routine requested approximately
  memory_size_request = 3*NGNOD*dble(8) * BYTES_TO_MB
  memory_size_request = memory_size_request + NGLLX*NGLLY*NGLLZ*nspec*dble(4) * BYTES_TO_MB
  memory_size_request = memory_size_request + 6*nspec*dble(1) * BYTES_TO_MB
  memory_size_request = memory_size_request + 6*NGLLX*dble(8) * BYTES_TO_MB
  memory_size_request = memory_size_request + NGNOD*NGLLX*NGLLY*NGLLZ*dble(8) * BYTES_TO_MB
  memory_size_request = memory_size_request + NDIM*NGNOD*NGLLX*NGLLY*NGLLZ*dble(8) * BYTES_TO_MB
  memory_size_request = memory_size_request + 4*NGNOD2D*NGLLY*NGLLZ*dble(8) * BYTES_TO_MB
  memory_size_request = memory_size_request + 4*NDIM2D*NGNOD2D*NGLLX*NGLLY*dble(8) * BYTES_TO_MB
  memory_size_request = memory_size_request + 17*NGLLX*NGLLY*NGLLY*nspec*dble(CUSTOM_REAL) * BYTES_TO_MB
  memory_size_request = memory_size_request + (1+NDIM)*NGLLY*NGLLZ*nspec2D_xmin*dble(CUSTOM_REAL) * BYTES_TO_MB
  memory_size_request = memory_size_request + (1+NDIM)*NGLLY*NGLLZ*nspec2D_xmax*dble(CUSTOM_REAL) * BYTES_TO_MB
  memory_size_request = memory_size_request + (1+NDIM)*NGLLX*NGLLZ*nspec2D_ymin*dble(CUSTOM_REAL) * BYTES_TO_MB
  memory_size_request = memory_size_request + (1+NDIM)*NGLLX*NGLLZ*nspec2D_ymax*dble(CUSTOM_REAL) * BYTES_TO_MB
  memory_size_request = memory_size_request + (1+NDIM)*NGLLX*NGLLY*NSPEC2D_BOTTOM*dble(CUSTOM_REAL) * BYTES_TO_MB
  memory_size_request = memory_size_request + (1+NDIM)*NGLLX*NGLLY*NSPEC2D_TOP*dble(CUSTOM_REAL) * BYTES_TO_MB
  memory_size_request = memory_size_request + 2*npointot*dble(4) * BYTES_TO_MB
  memory_size_request = memory_size_request + npointot*dble(1) * BYTES_TO_MB
  memory_size_request = memory_size_request + 3*npointot*dble(8) * BYTES_TO_MB

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  minimum memory used so far     : ', sngl(memory_size),'MB per process'
    write(IMAIN,*) '  minimum total memory requested : ', sngl(memory_size + memory_size_request),'MB per process'
    write(IMAIN,*)
  endif

  end subroutine memory_eval_mesher

