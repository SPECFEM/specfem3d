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

! for HDF file I/0

  subroutine save_arrays_solver_mesh_hdf5()

#ifdef USE_HDF5
  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM,NGLLSQUARE,IMAIN,USE_MESH_COLORING_GPU,CUSTOM_REAL, &
                       MAX_STRING_LEN

  use shared_parameters, only: ACOUSTIC_SIMULATION, ELASTIC_SIMULATION, POROELASTIC_SIMULATION, &
    APPROXIMATE_OCEAN_LOAD, ANISOTROPY, &
    COUPLE_WITH_INJECTION_TECHNIQUE, MESH_A_CHUNK_OF_THE_EARTH,NPROC

  ! global indices
  use generate_databases_par, only: NSPEC_AB, ibool, NGLOB_AB

  use generate_databases_par, only: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, &
    NSPEC2D_BOTTOM, NSPEC2D_TOP, &
    ibelm_xmin, ibelm_xmax,ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top, &
    SIMULATION_TYPE, SAVE_FORWARD, &
    STACEY_ABSORBING_CONDITIONS, &
    LOCAL_PATH, myrank

  ! MPI interfaces
  use generate_databases_par, only: num_interfaces_ext_mesh,my_neighbors_ext_mesh, &
    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh

  ! PML
  use generate_databases_par, only: PML_CONDITIONS, &
    nspec_cpml,CPML_width_x,CPML_width_y,CPML_width_z,CPML_to_spec, &
    CPML_regions,is_CPML,min_distance_between_CPML_parameter, &
    d_store_x,d_store_y,d_store_z,k_store_x,k_store_y,k_store_z, &
    alpha_store_x,alpha_store_y,alpha_store_z, &
    nglob_interface_PML_acoustic,points_interface_PML_acoustic, &
    nglob_interface_PML_elastic,points_interface_PML_elastic

  ! mesh surface
  use generate_databases_par, only: ispec_is_surface_external_mesh,iglob_is_surface_external_mesh, &
    nfaces_surface

  use create_regions_mesh_ext_par

  use manager_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! local parameters
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: max_nibool_interfaces_ext_mesh

  integer :: ier,i
  character(len=MAX_STRING_LEN) :: filename

  ! MPI variables
  integer :: info, comm

  ! if collective write
  logical, parameter :: if_col = .true.

  ! hdf5 valiables
  character(len=64) :: dset_name, tempstr

  ! element node connectivity for movie output
  ! the node ids are stored after dividing one NGLL*-th order spectral element into NGLLX*NGLLY*NGLLZ elements
  integer, dimension(9,nspec_ab*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)) :: spec_elm_conn_xdmf

  ! dummy arrays
  integer, dimension(1,1), parameter                    :: i2d_dummy = reshape((/0/),(/1,1/))
  integer, dimension(1,1,1), parameter                  :: i3d_dummy = reshape((/0/),(/1,1,1/))
  real(kind=CUSTOM_REAL), dimension(1,1), parameter     :: r2d_dummy = reshape((/0.0/),(/1,1/))
  real(kind=CUSTOM_REAL), dimension(1,1,1), parameter   :: r3d_dummy = reshape((/0.0/),(/1,1,1/))

  ! offset arrays
  integer, dimension(0:NPROC-1) :: offset_nglob
  integer, dimension(0:NPROC-1) :: offset_nspec
  integer, dimension(0:NPROC-1) :: offset_nspec_irregular
  integer, dimension(0:NPROC-1) :: offset_nglob_ocean
  integer, dimension(0:NPROC-1) :: offset_nspecporo
  integer, dimension(0:NPROC-1) :: offset_nspeccpml
  integer, dimension(0:NPROC-1) :: offset_nglob_interface_PML_acoustic
  integer, dimension(0:NPROC-1) :: offset_nglob_interface_PML_elastic
  integer, dimension(0:NPROC-1) :: offset_num_abs_boundary_faces
  integer, dimension(0:NPROC-1) :: offset_nglob_xy
  integer, dimension(0:NPROC-1) :: offset_nspec2D_xmin
  integer, dimension(0:NPROC-1) :: offset_nspec2D_xmax
  integer, dimension(0:NPROC-1) :: offset_nspec2D_ymin
  integer, dimension(0:NPROC-1) :: offset_nspec2D_ymax
  integer, dimension(0:NPROC-1) :: offset_nspec2D_bottom_ext
  integer, dimension(0:NPROC-1) :: offset_nspec2D_top_ext
  integer, dimension(0:NPROC-1) :: offset_num_free_surface_faces
  integer, dimension(0:NPROC-1) :: offset_num_coupling_ac_el_faces
  integer, dimension(0:NPROC-1) :: offset_num_coupling_ac_po_faces
  integer, dimension(0:NPROC-1) :: offset_num_coupling_el_po_faces
  integer, dimension(0:NPROC-1) :: offset_num_interfaces_ext_mesh
  integer, dimension(0:NPROC-1) :: offset_max_ni_bool_interfaces_ext_mesh
  integer, dimension(0:NPROC-1) :: offset_nspec_aniso
  integer, dimension(0:NPROC-1) :: offset_num_phase_ispec_acoustic
  integer, dimension(0:NPROC-1) :: offset_num_phase_ispec_elastic
  integer, dimension(0:NPROC-1) :: offset_num_phase_ispec_poroelastic
  integer, dimension(0:NPROC-1) :: offset_num_colors_acoustic
  integer, dimension(0:NPROC-1) :: offset_num_colors_elastic
  integer, dimension(0:NPROC-1) :: offset_nspec_ab
  integer, dimension(0:NPROC-1) :: offset_nglob_ab

  ! saves mesh file external_mesh.h5
  tempstr = "/external_mesh.h5"
  filename = LOCAL_PATH(1:len_trim(LOCAL_PATH))//trim(tempstr)

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '     database file: ',trim(filename)
    write(IMAIN,*) '     using HDF5 file format'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  !---------------------------.
  ! Setup the values to write |
  !---------------------------'
  !MPI interfaces
  ! gather interface information
  max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh(:))
  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 650')
  if (ier /= 0) stop 'error allocating array'
  do i = 1, num_interfaces_ext_mesh
    ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,i)
  enddo
  call synchronize_all()

  ! get MPI parameters
  call world_get_comm(comm)
  call world_get_info_null(info)

  ! initialize h5 object
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROC)

  !
  ! prepare offset arrays
  !
  call gather_all_all_singlei(nglob_ab,offset_nglob,NPROC) ! n globs in each proc
  call gather_all_all_singlei(nspec_ab,offset_nspec,NPROC) ! n spec in each proc

  call get_connectivity_for_movie(nspec_ab, ibool, spec_elm_conn_xdmf, sum(offset_nglob(0:myrank-1)))

  call gather_all_all_singlei(nspec_irregular,offset_nspec_irregular,NPROC) ! n spec in each proc

  if (APPROXIMATE_OCEAN_LOAD) &
    call gather_all_all_singlei(NGLOB_OCEAN,offset_nglob_ocean,NPROC)

  if (POROELASTIC_SIMULATION) &
    call gather_all_all_singlei(NSPEC_PORO,offset_nspecporo,NPROC)

  if (PML_CONDITIONS) &
    call gather_all_all_singlei(nspec_cpml,offset_nspeccpml,NPROC)

  if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
    call gather_all_all_singlei(nglob_interface_PML_acoustic, &
                                offset_nglob_interface_PML_acoustic,NPROC)
    call gather_all_all_singlei(nglob_interface_PML_elastic, &
                                offset_nglob_interface_PML_elastic,NPROC)
  endif

  call gather_all_all_singlei(num_abs_boundary_faces,offset_num_abs_boundary_faces,NPROC)
  call gather_all_all_singlei(nglob_xy,offset_nglob_xy,NPROC)
  call gather_all_all_singlei(nspec2D_xmin,offset_nspec2D_xmin,NPROC)
  call gather_all_all_singlei(nspec2D_xmax,offset_nspec2D_xmax,NPROC)
  call gather_all_all_singlei(nspec2D_ymin,offset_nspec2D_ymin,NPROC)
  call gather_all_all_singlei(nspec2D_ymax,offset_nspec2D_ymax,NPROC)
  call gather_all_all_singlei(nspec2D_bottom,offset_nspec2D_bottom_ext,NPROC)
  call gather_all_all_singlei(nspec2D_top,offset_nspec2D_top_ext,NPROC)

  call gather_all_all_singlei(num_free_surface_faces,offset_num_free_surface_faces,NPROC)
  call gather_all_all_singlei(num_coupling_ac_el_faces,offset_num_coupling_ac_el_faces,NPROC)
  call gather_all_all_singlei(num_coupling_ac_po_faces,offset_num_coupling_ac_po_faces,NPROC)
  call gather_all_all_singlei(num_coupling_el_po_faces,offset_num_coupling_el_po_faces,NPROC)
  call gather_all_all_singlei(num_interfaces_ext_mesh,offset_num_interfaces_ext_mesh,NPROC)
  call gather_all_all_singlei(max_nibool_interfaces_ext_mesh,offset_max_ni_bool_interfaces_ext_mesh,NPROC)

  if (ELASTIC_SIMULATION .and. ANISOTROPY) &
    call gather_all_all_singlei(nspec_aniso,offset_nspec_aniso,NPROC)

  if (ACOUSTIC_SIMULATION) &
    call gather_all_all_singlei(num_phase_ispec_acoustic,offset_num_phase_ispec_acoustic,NPROC)

  if (ELASTIC_SIMULATION) &
    call gather_all_all_singlei(num_phase_ispec_elastic,offset_num_phase_ispec_elastic,NPROC)

  if (POROELASTIC_SIMULATION) &
    call gather_all_all_singlei(num_phase_ispec_poroelastic,offset_num_phase_ispec_poroelastic,NPROC)

  if (USE_MESH_COLORING_GPU) then
    if (ACOUSTIC_SIMULATION) &
      call gather_all_all_singlei(num_colors_outer_acoustic+num_colors_inner_acoustic, &
                                  offset_num_colors_acoustic,NPROC)
    if (ELASTIC_SIMULATION) &
      call gather_all_all_singlei(num_colors_outer_elastic+num_colors_inner_elastic, &
                                  offset_num_colors_elastic,NPROC)
  endif

  call gather_all_all_singlei(nspec_ab,offset_nspec_ab,NPROC)
  call gather_all_all_singlei(nglob_ab,offset_nglob_ab,NPROC)

  !
  ! make datasets by main
  !

  ! offset arrays
  if (myrank == 0) then
    call h5_create_file(filename)

    call h5_write_dataset_no_group("offset_nglob",offset_nglob)
    call h5_write_dataset_no_group("offset_nspec",offset_nspec)
    call h5_write_dataset_no_group("offset_nspec_irregular",offset_nspec_irregular)

    if (APPROXIMATE_OCEAN_LOAD) &
      call h5_write_dataset_no_group("offset_nglob_ocean",offset_nglob_ocean)
    if (POROELASTIC_SIMULATION) &
      call h5_write_dataset_no_group("offset_nspecporo",offset_nspecporo)
    if (PML_CONDITIONS) &
      call h5_write_dataset_no_group("offset_nspeccpml",offset_nspeccpml)

    if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
      call h5_write_dataset_no_group("offset_nglob_interface_PML_acoustic",offset_nglob_interface_PML_acoustic)
      call h5_write_dataset_no_group("offset_nglob_interface_PML_elastic",offset_nglob_interface_PML_elastic)
    endif

    call h5_write_dataset_no_group("offset_num_abs_boundary_faces",offset_num_abs_boundary_faces)
    call h5_write_dataset_no_group("offset_nglob_xy",offset_nglob_xy)
    call h5_write_dataset_no_group("offset_nspec2D_xmin",offset_nspec2D_xmin)
    call h5_write_dataset_no_group("offset_nspec2D_xmax",offset_nspec2D_xmax)
    call h5_write_dataset_no_group("offset_nspec2D_ymin",offset_nspec2D_ymin)
    call h5_write_dataset_no_group("offset_nspec2D_ymax",offset_nspec2D_ymax)
    call h5_write_dataset_no_group("offset_nspec2D_bottom_ext",offset_nspec2D_bottom_ext)
    call h5_write_dataset_no_group("offset_nspec2D_top_ext",offset_nspec2D_top_ext)
    call h5_write_dataset_no_group("offset_num_free_surface_faces",offset_num_free_surface_faces)
    call h5_write_dataset_no_group("offset_num_coupling_ac_el_faces",offset_num_coupling_ac_el_faces)
    call h5_write_dataset_no_group("offset_num_coupling_ac_po_faces",offset_num_coupling_ac_po_faces)
    call h5_write_dataset_no_group("offset_num_coupling_el_po_faces",offset_num_coupling_el_po_faces)
    call h5_write_dataset_no_group("offset_num_interfaces_ext_mesh",offset_num_interfaces_ext_mesh)

    if (ELASTIC_SIMULATION .and. ANISOTROPY) &
      call h5_write_dataset_no_group("offset_nspec_aniso",offset_nspec_aniso)

    if (ACOUSTIC_SIMULATION) then
      if (sum(offset_num_phase_ispec_acoustic) > 0) &
        call h5_write_dataset_no_group("offset_num_phase_ispec_acoustic",offset_num_phase_ispec_acoustic)
    endif

    if (ELASTIC_SIMULATION) then
      if (sum(offset_num_phase_ispec_elastic) > 0) &
        call h5_write_dataset_no_group("offset_num_phase_ispec_elastic",offset_num_phase_ispec_elastic)
    endif

    if (POROELASTIC_SIMULATION) then
      if (sum(offset_num_phase_ispec_poroelastic) > 0) &
        call h5_write_dataset_no_group("offset_num_phase_ispec_poroelastic",offset_num_phase_ispec_poroelastic)
    endif

    if (USE_MESH_COLORING_GPU) then
      if (ACOUSTIC_SIMULATION) call h5_write_dataset_no_group("offset_num_colors_acoustic",offset_num_colors_acoustic)
      if (ELASTIC_SIMULATION)  call h5_write_dataset_no_group("offset_num_colors_elastic",offset_num_colors_elastic)
    endif

    call h5_write_dataset_no_group("offset_nspec_ab",offset_nspec_ab)
    call h5_write_dataset_no_group("offset_nglob_ab",offset_nglob_ab)

    ! other datasets

    dset_name = "nspec" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "nglob" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "nspec_irregular" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "ibool" ! 4 i (/0,0,0, offset_nglob/)
    call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nglob(:))/), 4, 1)

    dset_name = "xstore_unique" ! 1 r (/offset_nglob/)
    call h5_create_dataset_gen(dset_name, (/sum(offset_nglob(:))/),1, CUSTOM_REAL)
    dset_name = "ystore_unique" ! 1 r (/offset_nglob/)
    call h5_create_dataset_gen(dset_name, (/sum(offset_nglob(:))/),1, CUSTOM_REAL)
    dset_name = "zstore_unique" ! 1 r (/offset_nglob/)
    call h5_create_dataset_gen(dset_name, (/sum(offset_nglob(:))/),1, CUSTOM_REAL)

    dset_name = "irregular_element_number" ! 1 i (/offset_nspec/)
    call h5_create_dataset_gen(dset_name, (/sum(offset_nspec(:))/),1, 1)

    dset_name = "xix_regular" ! 1 r (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, CUSTOM_REAL)
    dset_name = "jacobian_regular" ! 1 r (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, CUSTOM_REAL)

    if (sum(offset_nspec_irregular) > 0) then
      dset_name = "xixstore" ! 4 r  (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_irregular(:))/), 4, CUSTOM_REAL)
      dset_name = "xiystore" ! 4 r  (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_irregular(:))/), 4, CUSTOM_REAL)
      dset_name = "xizstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_irregular(:))/), 4, CUSTOM_REAL)
      dset_name = "etaxstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_irregular(:))/), 4, CUSTOM_REAL)
      dset_name = "etaystore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_irregular(:))/), 4, CUSTOM_REAL)
      dset_name = "etazstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_irregular(:))/), 4, CUSTOM_REAL)
      dset_name = "gammaxstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_irregular(:))/), 4, CUSTOM_REAL)
      dset_name = "gammaystore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_irregular(:))/), 4, CUSTOM_REAL)
      dset_name = "gammazstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_irregular(:))/), 4, CUSTOM_REAL)
      dset_name = "jacobianstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_irregular(:))/), 4, CUSTOM_REAL)
    else
      ! dummy
      dset_name = "xixstore" ! 4 r  (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/1,1,1,1/), 4, CUSTOM_REAL)
      dset_name = "xiystore" ! 4 r  (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/1,1,1,1/), 4, CUSTOM_REAL)
      dset_name = "xizstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/1,1,1,1/), 4, CUSTOM_REAL)
      dset_name = "etaxstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/1,1,1,1/), 4, CUSTOM_REAL)
      dset_name = "etaystore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/1,1,1,1/), 4, CUSTOM_REAL)
      dset_name = "etazstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/1,1,1,1/), 4, CUSTOM_REAL)
      dset_name = "gammaxstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/1,1,1,1/), 4, CUSTOM_REAL)
      dset_name = "gammaystore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/1,1,1,1/), 4, CUSTOM_REAL)
      dset_name = "gammazstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/1,1,1,1/), 4, CUSTOM_REAL)
      dset_name = "jacobianstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/1,1,1,1/), 4, CUSTOM_REAL)
    endif

    dset_name = "kappastore" ! 4 r (/0,0,0,offset_nspec/)
    call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec(:))/), 4, CUSTOM_REAL)
    dset_name = "mustore" ! 4 r (/0,0,0,offset_nspec/)
    call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec(:))/), 4, CUSTOM_REAL)
    dset_name = "ispec_is_acoustic" ! 1 l (/0,0,0,offset_nspec/)
    call h5_create_dataset_gen(dset_name,(/sum(offset_nspec(:))/), 1, 0)
    dset_name = "ispec_is_elastic" ! 1 l (/0,0,0,offset_nspec/)
    call h5_create_dataset_gen(dset_name,(/sum(offset_nspec(:))/), 1, 0)
    dset_name = "ispec_is_poroelastic" ! 1 l (/0,0,0,offset_nspec/)
    call h5_create_dataset_gen(dset_name,(/sum(offset_nspec(:))/), 1, 0)

    ! acoustic
    if (ACOUSTIC_SIMULATION) then
      dset_name = "rmass_acoustic" ! 1 r (/offset_nglob/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_nglob(:))/), 1, CUSTOM_REAL)
    endif

    ! this array is needed for acoustic simulations but also for elastic simulations with CPML,
    ! thus we allocate it and read it in all cases (whether the simulation is acoustic, elastic, or acoustic/elastic)
    dset_name = "rhostore" ! 4 r (/0,0,0,offset_nspec/)
    call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec(:))/), 4, CUSTOM_REAL)

    ! elastic
    if (ELASTIC_SIMULATION) then
      dset_name = "rmass" ! 1 r (/offset_nglob/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_nglob(:))/), 1, CUSTOM_REAL)

      if (APPROXIMATE_OCEAN_LOAD) then
        dset_name = "rmass_ocean_load" ! 1 r (/offset_nglob_ocean/)
        call h5_create_dataset_gen(dset_name,(/sum(offset_nglob_ocean(:))/), 1, CUSTOM_REAL)
      endif

      !pll Stacey
      dset_name = "rho_vp" ! 4 r (/0,0,0,offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec(:))/), 4, CUSTOM_REAL)
      dset_name = "rho_vs" ! 4 r (/0,0,0,offset_nspec/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec(:))/), 4, CUSTOM_REAL)
    endif

    ! poroelastic
    if (POROELASTIC_SIMULATION) then
      dset_name = "rmass_solid_poroelastic" ! 1 r (/offset_nglob/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_nglob(:))/), 1, CUSTOM_REAL)
      dset_name = "rmass_fluid_poroelastic" ! 1 r (/offset_nglob/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_nglob(:))/), 1, CUSTOM_REAL)
      dset_name = "rhoarraystore" ! 5 r (/0,0,0,0,offset_nspecporo/)
      call h5_create_dataset_gen(dset_name,(/2,NGLLX,NGLLY,NGLLZ,sum(offset_nspecporo(:))/), 5, CUSTOM_REAL)
      dset_name = "kappaarraystore" ! 5 r (/0,0,0,0,offset_nspecporo/)
      call h5_create_dataset_gen(dset_name,(/3,NGLLX,NGLLY,NGLLZ,sum(offset_nspecporo(:))/), 5, CUSTOM_REAL)
      dset_name = "etastore" ! 4 r (/0,0,0,offset_nspecporo/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspecporo(:))/), 4, CUSTOM_REAL)
      dset_name = "tortstore" ! 4 r (/0,0,0,offset_nspecporo/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspecporo(:))/), 4, CUSTOM_REAL)
      dset_name = "permstore" ! 5 r (/0,0,0,0,offset_nspecporo/)
      call h5_create_dataset_gen(dset_name,(/6,NGLLX,NGLLY,NGLLZ,sum(offset_nspecporo(:))/), 5, CUSTOM_REAL)
      dset_name = "phistore" ! 4 r (/0,0,0,offset_nspecporo/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspecporo(:))/), 4, CUSTOM_REAL)
      dset_name = "rho_vpI" ! 4 r (/0,0,0,offset_nspecporo/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspecporo(:))/), 4, CUSTOM_REAL)
      dset_name = "rho_vpII" ! 4 r (/0,0,0,offset_nspecporo/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspecporo(:))/), 4, CUSTOM_REAL)
      dset_name = "rho_vsI" ! 4 r (/0,0,0,offset_nspecporo/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspecporo(:))/), 4, CUSTOM_REAL)
    endif

    ! C-PML absorbing boundary conditions
    if (PML_CONDITIONS) then
      dset_name = "nspec_cpml" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "CPML_width_x" ! 1 r (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, CUSTOM_REAL)
      dset_name = "CPML_width_y" ! 1 r (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, CUSTOM_REAL)
      dset_name = "CPML_width_z" ! 1 r (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, CUSTOM_REAL)
      dset_name = "min_distance_between_CPML_parameter" ! 1 r (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, CUSTOM_REAL)

      if (sum(offset_nspeccpml) > 0) then
        dset_name = "CPML_regions" ! 1 i (/offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/sum(offset_nspeccpml(:))/), 1, 1)
        dset_name = "CPML_to_spec" ! 1 i (/offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/sum(offset_nspeccpml(:))/), 1, 1)
        dset_name = "is_CPML" ! 1 l (/offset_nspecab/)
        call h5_create_dataset_gen(dset_name,(/sum(offset_nspec_ab(:))/), 1, 1)
        dset_name = "d_store_x" ! 4 r (/0,0,0,offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspeccpml(:))/), 4, CUSTOM_REAL)
        dset_name = "d_store_y" ! 4 r (/0,0,0,offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspeccpml(:))/), 4, CUSTOM_REAL)
        dset_name = "d_store_z" ! 4 r (/0,0,0,offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspeccpml(:))/), 4, CUSTOM_REAL)
        dset_name = "k_store_x" ! 4 r (/0,0,0,offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspeccpml(:))/), 4, CUSTOM_REAL)
        dset_name = "k_store_y" ! 4 r (/0,0,0,offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspeccpml(:))/), 4, CUSTOM_REAL)
        dset_name = "k_store_z" ! 4 r (/0,0,0,offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspeccpml(:))/), 4, CUSTOM_REAL)
        dset_name = "alpha_store_x" ! 4 r (/0,0,0,offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspeccpml(:))/), 4, CUSTOM_REAL)
        dset_name = "alpha_store_y" ! 4 r (/0,0,0,offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspeccpml(:))/), 4, CUSTOM_REAL)
        dset_name = "alpha_store_z" ! 4 r (/0,0,0,offset_nspeccpml/)
        call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspeccpml(:))/), 4, CUSTOM_REAL)

        ! --------------------------------------------------------------------------------------------
        ! for adjoint tomography
        ! save the array stored the points on interface between PML and interior computational domain
        ! --------------------------------------------------------------------------------------------
        if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
          dset_name = "nglob_interface_PML_acoustic" ! 1 i (/myrank/)
          call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
          dset_name = "nglob_interface_PML_elastic" ! 1 i (/myrank/)
          call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
          !if (nglob_interface_PML_acoustic > 0) then
          if (sum(offset_nglob_interface_PML_acoustic) > 0) then
            dset_name = "points_interface_PML_acoustic" ! 1 i (/offset_nglob_interface_PML_acoustic/)
            call h5_create_dataset_gen(dset_name,(/sum(offset_nglob_interface_PML_acoustic)/), 1, 1)
          endif
          !if (nglob_interface_PML_elastic > 0) then
          if (sum(offset_nglob_interface_PML_elastic) > 0) then
            dset_name = "points_interface_PML_elastic" ! 1 i (/offset_nglob_interface_PML_elastic/)
            call h5_create_dataset_gen(dset_name,(/sum(offset_nglob_interface_PML_elastic)/), 1, 1)
          endif
        endif
      endif
    endif   ! PML

    ! absorbing boundary surface
    dset_name = "num_abs_boundary_faces" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)

    if (sum(offset_num_abs_boundary_faces) > 0) then
      dset_name = "abs_boundary_ispec" ! 1 i (/offset_num_abs_boundary_faces/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_num_abs_boundary_faces(:))/), 1, 1)
      dset_name = "abs_boundary_ijk" ! 3 i (/0,0,offset_num_abs_boundary_faces/)
      call h5_create_dataset_gen(dset_name,(/3,NGLLSQUARE,sum(offset_num_abs_boundary_faces(:))/), 3, 1)
      dset_name = "abs_boundary_jacobian2Dw" ! 2 r (/0,offset_num_abs_boundary_faces/)
      call h5_create_dataset_gen(dset_name,(/NGLLSQUARE,sum(offset_num_abs_boundary_faces(:))/), 2, CUSTOM_REAL)
      dset_name = "abs_boundary_normal" ! 3 r (/0,0,offset_num_abs_boundary_faces/)
      call h5_create_dataset_gen(dset_name,(/NDIM,NGLLSQUARE,sum(offset_num_abs_boundary_faces(:))/), 3, CUSTOM_REAL)

      if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
        ! store mass matrix contributions
        if (ELASTIC_SIMULATION) then
          dset_name = "rmassx" ! 1 r (/offset_nglob_xy/)
          call h5_create_dataset_gen(dset_name, (/sum(offset_nglob_xy(:))/),1,CUSTOM_REAL)
          dset_name = "rmassy" ! 1 r (/offset_nglob_xy/)
          call h5_create_dataset_gen(dset_name, (/sum(offset_nglob_xy(:))/),1,CUSTOM_REAL)
          dset_name = "rmassz" ! 1 r (/offset_nglob_xy/)
          call h5_create_dataset_gen(dset_name, (/sum(offset_nglob_xy(:))/),1,CUSTOM_REAL)
        endif
        if (ACOUSTIC_SIMULATION) then
          dset_name = "rmassz_acoustic" ! 1 r (/offset_nglob_xy/)
          call h5_create_dataset_gen(dset_name, (/sum(offset_nglob_xy(:))/),1,CUSTOM_REAL)
        endif
      endif
    else
      ! dummy
      dset_name = "abs_boundary_ispec" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "abs_boundary_ijk" ! 3 i (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, 1)
      dset_name = "abs_boundary_jacobian2Dw" ! 2 r (/0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,NPROC/), 2, CUSTOM_REAL)
      dset_name = "abs_boundary_normal" ! 3 r (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, CUSTOM_REAL)

      if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
        ! store mass matrix contributions
        if (ELASTIC_SIMULATION) then
          dset_name = "rmassx" ! 1 r (/myrank/)
          call h5_create_dataset_gen(dset_name,(/NPROC/), 1, CUSTOM_REAL)
          dset_name = "rmassy" ! 1 r (/myrank/)
          call h5_create_dataset_gen(dset_name,(/NPROC/), 1, CUSTOM_REAL)
          dset_name = "rmassz" ! 1 r (/myrank/)
          call h5_create_dataset_gen(dset_name,(/NPROC/), 1, CUSTOM_REAL)
       endif
        if (ACOUSTIC_SIMULATION) then
          dset_name = "rmassz_acoustic" ! 1 r (/myrank/)
          call h5_create_dataset_gen(dset_name,(/NPROC/), 1, CUSTOM_REAL)
       endif
      endif
    endif

    dset_name = "nspec2D_xmin" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "nspec2D_xmax" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "nspec2D_ymin" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "nspec2D_ymax" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "NSPEC2D_BOTTOM" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "NSPEC2D_TOP" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)

    if (sum(offset_nspec2D_xmin) > 0) then
      dset_name = "ibelm_xmin" ! 1 i (/offset_nspec2D_xmin/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_nspec2D_xmin(:))/), 1, 1)
    endif
    if (sum(offset_nspec2D_xmax) > 0) then
      dset_name = "ibelm_xmax" ! 1 i (/offset_nspec2D_xmax/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_nspec2D_xmax(:))/), 1, 1)
    endif
    if (sum(offset_nspec2D_ymin) > 0) then
      dset_name = "ibelm_ymin" ! 1 i (/offset_nspec2D_ymin/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_nspec2D_ymin(:))/), 1, 1)
    endif
    if (sum(offset_nspec2D_ymax) > 0) then
      dset_name = "ibelm_ymax" ! 1 i (/offset_nspec2D_ymax/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_nspec2D_ymax(:))/), 1, 1)
    endif
    if (sum(offset_nspec2D_bottom_ext) > 0) then
      dset_name = "ibelm_bottom" ! 1 i (/offset_nspec2D_bottom_ext/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_nspec2D_bottom_ext(:))/), 1, 1)
    endif
    if (sum(offset_nspec2D_top_ext) > 0) then
      dset_name = "ibelm_top" ! 1 i (/offset_nspec2D_top_ext/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_nspec2D_top_ext(:))/), 1, 1)
    endif

    ! free surface
    dset_name = "num_free_surface_faces" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)

    if (sum(offset_num_free_surface_faces) > 0) then
      dset_name = "free_surface_ispec" ! 1 i (/offset_num_free_surface_faces/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_num_free_surface_faces(:))/), 1, 1)
      dset_name = "free_surface_ijk" ! 3 i (/0,0,offset_num_free_surface_faces/)
      call h5_create_dataset_gen(dset_name,(/3,NGLLSQUARE,sum(offset_num_free_surface_faces(:))/), 3, 1)
      dset_name = "free_surface_jacobian2Dw" ! 2 r (/0,offset_num_free_surface_faces/)
      call h5_create_dataset_gen(dset_name,(/NGLLSQUARE,sum(offset_num_free_surface_faces(:))/), 2, CUSTOM_REAL)
      dset_name = "free_surface_normal" ! 3 r (/0,0,offset_num_free_surface_faces/)
      call h5_create_dataset_gen(dset_name,(/NDIM,NGLLSQUARE,sum(offset_num_free_surface_faces(:))/), 3, CUSTOM_REAL)
    else
      dset_name = "free_surface_ispec" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "free_surface_ijk" ! 3 i (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, 1)
      dset_name = "free_surface_jacobian2Dw" ! 2 r (/0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,NPROC/), 2, CUSTOM_REAL)
      dset_name = "free_surface_normal" ! 3 r (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, CUSTOM_REAL)
    endif

    ! acoustic-elastic coupling surface
    dset_name = "num_coupling_ac_el_faces" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)

    if (sum(offset_num_coupling_ac_el_faces) > 0) then
      dset_name = "coupling_ac_el_ispec" ! 1 i (/offset_num_coupling_ac_el_faces/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_num_coupling_ac_el_faces(:))/), 1, 1)
      dset_name = "coupling_ac_el_ijk" ! 3 i (/0,0,offset_num_coupling_ac_el_faces/)
      call h5_create_dataset_gen(dset_name,(/3,NGLLSQUARE,sum(offset_num_coupling_ac_el_faces(:))/), 3, 1)
      dset_name = "coupling_ac_el_jacobian2Dw" ! 2 r (/0,offset_num_coupling_ac_el_faces/)
      call h5_create_dataset_gen(dset_name,(/NGLLSQUARE,sum(offset_num_coupling_ac_el_faces(:))/), 2, CUSTOM_REAL)
      dset_name = "coupling_ac_el_normal" ! 3 r (/0,0,offset_num_coupling_ac_el_faces/)
      call h5_create_dataset_gen(dset_name,(/NDIM,NGLLSQUARE,sum(offset_num_coupling_ac_el_faces(:))/), 3, CUSTOM_REAL)
    else
      dset_name = "coupling_ac_el_ispec" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "coupling_ac_el_ijk" ! 3 i (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, 1)
      dset_name = "coupling_ac_el_jacobian2Dw" ! 2 r (/0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,NPROC/), 2, CUSTOM_REAL)
      dset_name = "coupling_ac_el_normal" ! 3 r (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, CUSTOM_REAL)
    endif

    ! acoustic-poroelastic coupling surface
    dset_name = "num_coupling_ac_po_faces" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)

    if (sum(offset_num_coupling_ac_po_faces) > 0) then
      dset_name = "coupling_ac_po_ispec" ! 1 i (/offset_num_coupling_ac_po_faces/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_num_coupling_ac_po_faces(:))/), 1, 1)
      dset_name = "coupling_ac_po_ijk" ! 3 i (/0,0,offset_num_coupling_ac_po_faces/)
      call h5_create_dataset_gen(dset_name,(/3,NGLLSQUARE,sum(offset_num_coupling_ac_po_faces(:))/), 3, 1)
      dset_name = "coupling_ac_po_jacobian2Dw" ! 2 r (/0,offset_num_coupling_ac_po_faces/)
      call h5_create_dataset_gen(dset_name,(/NGLLSQUARE,sum(offset_num_coupling_ac_po_faces(:))/), 2, CUSTOM_REAL)
      dset_name = "coupling_ac_po_normal" ! 3 r (/0,0,offset_num_coupling_ac_po_faces/)
      call h5_create_dataset_gen(dset_name,(/NDIM,NGLLSQUARE,sum(offset_num_coupling_ac_po_faces(:))/), 3, CUSTOM_REAL)
    else
      dset_name = "coupling_ac_po_ispec" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "coupling_ac_po_ijk" ! 3 i (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, 1)
      dset_name = "coupling_ac_po_jacobian2Dw" ! 2 r (/0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,NPROC/), 2, CUSTOM_REAL)
      dset_name = "coupling_ac_po_normal" ! 3 r (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, CUSTOM_REAL)
    endif

    ! elastic-poroelastic coupling surface
    dset_name = "num_coupling_el_po_faces" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)

    if (sum(offset_num_coupling_el_po_faces) > 0) then
      dset_name = "coupling_el_po_ispec" ! 1 i (/offset_num_coupling_el_po_faces/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_num_coupling_el_po_faces(:))/), 1, 1)
      dset_name = "coupling_po_el_ispec" ! 1 i (/offset_num_coupling_el_po_faces/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_num_coupling_el_po_faces(:))/), 1, 1)
      dset_name = "coupling_el_po_ijk" ! 3 i (/0,0,offset_num_coupling_el_po_faces/)
      call h5_create_dataset_gen(dset_name,(/3,NGLLSQUARE,sum(offset_num_coupling_el_po_faces(:))/), 3, 1)
      dset_name = "coupling_po_el_ijk" ! 3 i (/0,0,offset_num_coupling_el_po_faces/)
      call h5_create_dataset_gen(dset_name,(/3,NGLLSQUARE,sum(offset_num_coupling_el_po_faces(:))/), 3, 1)
      dset_name = "coupling_el_po_jacobian2Dw" ! 2 r (/0,offset_num_coupling_el_po_faces/)
      call h5_create_dataset_gen(dset_name,(/NGLLSQUARE,sum(offset_num_coupling_el_po_faces(:))/), 2, CUSTOM_REAL)
      dset_name = "coupling_el_po_normal" ! 3 r (/0,0,offset_num_coupling_el_po_faces/)
      call h5_create_dataset_gen(dset_name,(/NDIM,NGLLSQUARE,sum(offset_num_coupling_el_po_faces(:))/), 3, CUSTOM_REAL)
    else
      dset_name = "coupling_el_po_ispec" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "coupling_po_el_ispec" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "coupling_el_po_ijk" ! 3 i (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, 1)
      dset_name = "coupling_po_el_ijk" ! 3 i (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, 1)
      dset_name = "coupling_el_po_jacobian2Dw" ! 2 r (/0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,NPROC/), 2, CUSTOM_REAL)
      dset_name = "coupling_el_po_normal" ! 3 r (/0,0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,1,NPROC/), 3, CUSTOM_REAL)
    endif

    dset_name = "num_interfaces_ext_mesh" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)

    if (sum(offset_num_interfaces_ext_mesh) > 0) then
      dset_name = "max_nibool_interfaces_ext_mesh" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "my_neighbors_ext_mesh" ! 1 i (/offset_num_interfaces_ext_mesh/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_num_interfaces_ext_mesh(:))/), 1, 1)
      dset_name = "nibool_interfaces_ext_mesh" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/sum(offset_num_interfaces_ext_mesh(:))/), 1, 1)
      dset_name = "ibool_interfaces_ext_mesh_dummy" ! 2 i (/offset_max_ni_bool_interfaces_ext_mesh/)
      call h5_create_dataset_gen(dset_name, &
              (/maxval(offset_max_ni_bool_interfaces_ext_mesh),sum(offset_num_interfaces_ext_mesh(:))/), 2, 1)
    else
      dset_name = "max_nibool_interfaces_ext_mesh" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "my_neighbors_ext_mesh" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "nibool_interfaces_ext_mesh" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "ibool_interfaces_ext_mesh_dummy" ! 2 i (/0,myrank/)
      call h5_create_dataset_gen(dset_name,(/1,NPROC/), 2, 1)
    endif

    ! anisotropy
    if (ELASTIC_SIMULATION .and. ANISOTROPY) then
      dset_name = "c11store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c12store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c13store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c14store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c15store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c16store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c22store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c23store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c24store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c25store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c26store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c33store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c34store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c35store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c36store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c44store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c45store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c46store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c55store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c56store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
      dset_name = "c66store" ! 4 r (/0,0,0,offset_nspec_aniso/)
      call h5_create_dataset_gen(dset_name,(/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_aniso(:))/), 4, CUSTOM_REAL)
    endif

    ! inner/outer elements
    dset_name = "ispec_is_inner" ! 1 l (/offset_nspec/)
    call h5_create_dataset_gen(dset_name,(/sum(offset_nspec(:))/), 1, 0)

    if (ACOUSTIC_SIMULATION) then
      dset_name = "nspec_inner_acoustic" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "nspec_outer_acoustic" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "num_phase_ispec_acoustic" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)

      if (sum(offset_num_phase_ispec_acoustic) > 0) then
        dset_name = "phase_ispec_inner_acoustic" ! 2 i (/offset_num_phase_ispec_acoustic, 0/)
        call h5_create_dataset_gen(dset_name,(/sum(offset_num_phase_ispec_acoustic(:)),2/), 2, 1)
      else
        dset_name = "phase_ispec_inner_acoustic" ! 2 i (/myrank,0/)
        call h5_create_dataset_gen(dset_name,(/NPROC,2/), 2, 1)
      endif
    endif

    if (ELASTIC_SIMULATION) then
      dset_name = "nspec_inner_elastic" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "nspec_outer_elastic" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "num_phase_ispec_elastic" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)

      if (sum(offset_num_phase_ispec_elastic) > 0) then
        dset_name = "phase_ispec_inner_elastic" ! 2 i (/offset_num_phase_ispec_elastic,0/)
        call h5_create_dataset_gen(dset_name,(/sum(offset_num_phase_ispec_elastic(:)),2/), 2, 1)
      else
        dset_name = "phase_ispec_inner_elastic" ! 2 i (/myrank, 0/)
        call h5_create_dataset_gen(dset_name,(/NPROC,2/), 2, 1)
      endif
    endif

    if (POROELASTIC_SIMULATION) then
      dset_name = "nspec_inner_poroelastic" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "nspec_outer_poroelastic" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
      dset_name = "num_phase_ispec_poroelastic" ! 1 i (/myrank/)
      call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)

      if (sum(offset_num_phase_ispec_poroelastic) > 0) then
        dset_name = "phase_ispec_inner_poroelastic" ! 2 i (/offset_num_phase_ispec_poroelastic/)
        call h5_create_dataset_gen(dset_name,(/sum(offset_num_phase_ispec_poroelastic(:)),2/), 2, 1)
      else
        dset_name = "phase_ispec_inner_poroelastic" ! 2 i (/myrank,0/)
        call h5_create_dataset_gen(dset_name,(/NPROC,2/), 2, 1)
      endif
    endif

    ! mesh coloring
    if (USE_MESH_COLORING_GPU) then
      if (ACOUSTIC_SIMULATION) then
        dset_name = "num_colors_outer_acoustic" ! 1 i (/myrank/)
        call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
        dset_name = "num_colors_inner_acoustic" ! 1 i (/myrank/)
        call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
        dset_name = "num_elem_colors_acoustic" ! 1 i (/offset_num_colors_outer_acoustic+num_colors_inner_acoustic/)
        call h5_create_dataset_gen(dset_name,(/sum(offset_num_colors_acoustic(:))/), 1, 1)
      endif
      if (ELASTIC_SIMULATION) then
        dset_name = "num_colors_outer_elastic" ! 1 i (/myrank/)
        call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
        dset_name = "num_colors_inner_elastic" ! 1 i (/myrank/)
        call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
        dset_name = "num_elem_colors_elastic" ! 1 i (/offset_num_colors_outer_elastic+num_colors_inner_elastic/)
        call h5_create_dataset_gen(dset_name,(/sum(offset_num_colors_elastic(:))/), 1, 1)
      endif
    endif

    ! surface points
    dset_name = "nfaces_surface" ! 1 i (/myrank/)
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "ispec_is_surface_external_mesh" ! 1 l (/offset_nspec_ab/)
    call h5_create_dataset_gen(dset_name,(/sum(offset_nspec_ab(:))/), 1, 0)
    dset_name = "iglob_is_surface_external_mesh" ! 1 l (/offset_nglob_ab/)
    call h5_create_dataset_gen(dset_name,(/sum(offset_nglob_ab(:))/), 1, 0)
    ! arrays for visualization
    dset_name = "spec_elm_conn_xdmf" ! 2 i (/0,offset_nspec*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)/)
    call h5_create_dataset_gen(dset_name,(/9,sum(offset_nspec(:)*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1))/), 2, 1)

    call h5_close_file()
  endif
  call synchronize_all()

  !
  ! write arrays into datasets collectively
  !
  call h5_open_file_p_collect(filename)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "     start dataset preparation and write in h5 format"
    call flush_IMAIN()
  endif

  ! set dwrite flagif_colto pre_define the dataset on file before write.

  dset_name = "nspec" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name,(/nspec_ab/), (/myrank/),if_col)
  dset_name = "nglob" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name,(/nglob_ab/), (/myrank/),if_col)
  dset_name = "nspec_irregular" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name,(/nspec_irregular/), (/myrank/),if_col)
  dset_name = "ibool" ! 4 i (/0,0,0, offset_nglobs/)
  call h5_write_dataset_collect_hyperslab(dset_name, ibool, (/0,0,0,sum(offset_nglob(0:myrank-1))/), if_col)
  dset_name = "xstore_unique" ! 1 r (/offset_nglobs/)
  call h5_write_dataset_collect_hyperslab(dset_name,xstore_unique,(/sum(offset_nglob(0:myrank-1))/),if_col)
  dset_name = "ystore_unique" ! 1 r (/offset_nglobs/)
  call h5_write_dataset_collect_hyperslab(dset_name,ystore_unique,(/sum(offset_nglob(0:myrank-1))/),if_col)
  dset_name = "zstore_unique" ! 1 r (/offset_nglobs/)
  call h5_write_dataset_collect_hyperslab(dset_name,zstore_unique,(/sum(offset_nglob(0:myrank-1))/),if_col)
  dset_name = "irregular_element_number" ! 1 i (/offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,irregular_element_number,(/sum(offset_nspec(0:myrank-1))/),if_col)
  dset_name = "xix_regular" ! 1 r (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name,(/xix_regular/),(/myrank/),if_col)
  dset_name = "jacobian_regular" ! 1 r (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name,(/jacobian_regular/),(/myrank/),if_col)
  dset_name = "xixstore" ! 4 r  (/0,0,0,offset_nspec_irregular=offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,xixstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
  dset_name = "xiystore" ! 4 r  (/0,0,0,offset_nspec_irregular=offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,xiystore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
  dset_name = "xizstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,xizstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
  dset_name = "etaxstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,etaxstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
  dset_name = "etaystore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,etaystore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
  dset_name = "etazstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,etazstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
  dset_name = "gammaxstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,gammaxstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
  dset_name = "gammaystore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,gammaystore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
  dset_name = "gammazstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,gammazstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
  dset_name = "jacobianstore" ! 4 r (/0,0,0,offset_nspec_irregular=offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,jacobianstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
  dset_name = "kappastore" ! 4 r (/0,0,0,offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,kappastore,(/0,0,0,sum(offset_nspec(0:myrank-1))/),if_col)
  dset_name = "mustore" ! 4 r (/0,0,0,offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name,mustore,(/0,0,0,sum(offset_nspec(0:myrank-1))/),if_col)
  dset_name = "ispec_is_acoustic" ! 1 l (/0,0,0,offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name, ispec_is_acoustic, (/sum(offset_nspec(0:myrank-1))/),if_col)
  dset_name = "ispec_is_elastic" ! 1 l (/0,0,0,offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name, ispec_is_elastic, (/sum(offset_nspec(0:myrank-1))/),if_col)
  dset_name = "ispec_is_poroelastic" ! 1 l (/0,0,0,offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name, ispec_is_poroelastic, (/sum(offset_nspec(0:myrank-1))/),if_col)

  ! acoustic
  if (ACOUSTIC_SIMULATION) then
    dset_name = "rmass_acoustic" ! 1 r (/offset_nglob/)
    call h5_write_dataset_collect_hyperslab(dset_name, rmass_acoustic,(/sum(offset_nglob(0:myrank-1))/),if_col)
  endif

  ! this array is needed for acoustic simulations but also for elastic simulations with CPML,
  ! thus we allocate it and read it in all cases (whether the simulation is acoustic, elastic, or acoustic/elastic)
  dset_name = "rhostore" ! 4 r (/0,0,0,offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name, rhostore, (/0,0,0,sum(offset_nspec(0:myrank-1))/),if_col)

  ! elastic
  if (ELASTIC_SIMULATION) then
    dset_name = "rmass" ! 1 r (/offset_nglob/)
    call h5_write_dataset_collect_hyperslab(dset_name, rmass, (/sum(offset_nglob(0:myrank-1))/), if_col)
    if (APPROXIMATE_OCEAN_LOAD) then
      dset_name = "rmass_ocean_load" ! 1 r (/offset_nglob_ocean/)
      call h5_write_dataset_collect_hyperslab(dset_name, rmass_ocean_load, (/sum(offset_nglob_ocean(0:myrank-1))/),if_col)
    endif
    !pll Stacey
    dset_name = "rho_vp" ! 4 r (/0,0,0,offset_nspec/)
    call h5_write_dataset_collect_hyperslab(dset_name, rho_vp, (/0,0,0,sum(offset_nspec(0:myrank-1))/),if_col)
    dset_name = "rho_vs" ! 4 r (/0,0,0,offset_nspec/)
    call h5_write_dataset_collect_hyperslab(dset_name, rho_vs, (/0,0,0,sum(offset_nspec(0:myrank-1))/),if_col)
  endif

  ! poroelastic
  if (POROELASTIC_SIMULATION) then
    dset_name = "rmass_solid_poroelastic" ! 1 r (/offset_nglob/)
    call h5_write_dataset_collect_hyperslab(dset_name, rmass_solid_poroelastic, (/sum(offset_nglob(0:myrank-1))/),if_col)
    dset_name = "rmass_fluid_poroelastic" ! 1 r (/offset_nglob/)
    call h5_write_dataset_collect_hyperslab(dset_name, rmass_fluid_poroelastic, (/sum(offset_nglob(0:myrank-1))/),if_col)
    dset_name = "rhoarraystore" ! 5 r (/0,0,0,0,offset_nspecporo/)
    call h5_write_dataset_collect_hyperslab(dset_name, rhoarraystore, (/0,0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
    dset_name = "kappaarraystore" ! 5 r (/0,0,0,0,offset_nspecporo/)
    call h5_write_dataset_collect_hyperslab(dset_name, kappaarraystore, (/0,0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
    dset_name = "etastore" ! 4 r (/0,0,0,offset_nspecporo/)
    call h5_write_dataset_collect_hyperslab(dset_name, etastore, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
    dset_name = "tortstore" ! 4 r (/0,0,0,offset_nspecporo/)
    call h5_write_dataset_collect_hyperslab(dset_name, tortstore, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
    dset_name = "permstore" ! 5 r (/0,0,0,0,offset_nspecporo/)
    call h5_write_dataset_collect_hyperslab(dset_name, permstore, (/0,0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
    dset_name = "phistore" ! 4 r (/0,0,0,offset_nspecporo/)
    call h5_write_dataset_collect_hyperslab(dset_name, phistore, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
    dset_name = "rho_vpI" ! 4 r (/0,0,0,offset_nspecporo/)
    call h5_write_dataset_collect_hyperslab(dset_name, rho_vpI, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
    dset_name = "rho_vpII" ! 4 r (/0,0,0,offset_nspecporo/)
    call h5_write_dataset_collect_hyperslab(dset_name, rho_vpII, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
    dset_name = "rho_vsI" ! 4 r (/0,0,0,offset_nspecporo/)
    call h5_write_dataset_collect_hyperslab(dset_name, rho_vsI, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
  endif

  ! C-PML absorbing boundary conditions
  if (PML_CONDITIONS) then
    dset_name = "nspec_cpml" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/nspec_cpml/), (/myrank/),if_col)
    dset_name = "CPML_width_x" ! 1 r (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/CPML_width_x/), (/myrank/),if_col)
    dset_name = "CPML_width_y" ! 1 r (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/CPML_width_y/), (/myrank/),if_col)
    dset_name = "CPML_width_z" ! 1 r (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/CPML_width_z/), (/myrank/),if_col)
    dset_name = "min_distance_between_CPML_parameter" ! 1 r (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/min_distance_between_CPML_parameter/), (/myrank/),if_col)

    if (sum(offset_nspeccpml) > 0) then
      dset_name = "CPML_regions" ! 1 i (/offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, CPML_regions, (/sum(offset_nspeccpml(0:myrank-1))/),if_col)
      dset_name = "CPML_to_spec" ! 1 i (/offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, CPML_to_spec, (/sum(offset_nspeccpml(0:myrank-1))/),if_col)
      dset_name = "is_CPML" ! 1 l (/offset_nspecab/)
      call h5_write_dataset_collect_hyperslab(dset_name, is_CPML, (/sum(offset_nspec_ab(0:myrank-1))/),if_col)
      dset_name = "d_store_x" ! 4 r (/0,0,0,offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, d_store_x, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
      dset_name = "d_store_y" ! 4 r (/0,0,0,offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, d_store_y, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
      dset_name = "d_store_z" ! 4 r (/0,0,0,offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, d_store_z, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
      dset_name = "k_store_x" ! 4 r (/0,0,0,offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, k_store_x, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
      dset_name = "k_store_y" ! 4 r (/0,0,0,offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, k_store_y, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
      dset_name = "k_store_z" ! 4 r (/0,0,0,offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, k_store_z, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
      dset_name = "alpha_store_x" ! 4 r (/0,0,0,offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, alpha_store_x, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
      dset_name = "alpha_store_y" ! 4 r (/0,0,0,offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, alpha_store_y, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
      dset_name = "alpha_store_z" ! 4 r (/0,0,0,offset_nspeccpml/)
      call h5_write_dataset_collect_hyperslab(dset_name, alpha_store_z, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)

      ! --------------------------------------------------------------------------------------------
      ! for adjoint tomography
      ! save the array stored the points on interface between PML and interior computational domain
      ! --------------------------------------------------------------------------------------------
      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        dset_name = "nglob_interface_PML_acoustic" ! 1 i (/myrank/)
        call h5_write_dataset_collect_hyperslab(dset_name, (/nglob_interface_PML_acoustic/), (/myrank/),if_col)
        dset_name = "nglob_interface_PML_elastic" ! 1 i (/myrank/)
        call h5_write_dataset_collect_hyperslab(dset_name, (/nglob_interface_PML_elastic/), (/myrank/),if_col)
        if (sum(offset_nglob_interface_PML_acoustic) > 0) then
          dset_name = "points_interface_PML_acoustic" ! 1 i (/offset_nglob_interface_PML_acoustic/)
          call h5_write_dataset_collect_hyperslab(dset_name, &
              points_interface_PML_acoustic, (/sum(offset_nglob_interface_PML_acoustic(0:myrank-1))/),if_col)
        endif
        if (sum(offset_nglob_interface_PML_elastic) > 0) then
          dset_name = "points_interface_PML_elastic" ! 1 i (/offset_nglob_interface_PML_elastic/)
          call h5_write_dataset_collect_hyperslab(dset_name, &
              points_interface_PML_elastic, (/sum(offset_nglob_interface_PML_elastic(0:myrank-1))/),if_col)
        endif
      endif
    endif
  endif

  ! absorbing boundary surface
  dset_name = "num_abs_boundary_faces" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/num_abs_boundary_faces/), (/myrank/),if_col)

  if (sum(offset_num_abs_boundary_faces) > 0) then
    dset_name = "abs_boundary_ispec" ! 1 i (/offset_num_abs_boundary_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, abs_boundary_ispec, &
            (/sum(offset_num_abs_boundary_faces(0:myrank-1))/),if_col)
    dset_name = "abs_boundary_ijk" ! 3 i (/0,0,offset_num_abs_boundary_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, abs_boundary_ijk, &
            (/0,0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), if_col)
    dset_name = "abs_boundary_jacobian2Dw" ! 2 r (/0,offset_num_abs_boundary_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, abs_boundary_jacobian2Dw, &
            (/0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), if_col)
    dset_name = "abs_boundary_normal" ! 3 r (/0,0,offset_num_abs_boundary_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, abs_boundary_normal, &
            (/0,0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), if_col)

    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! store mass matrix contributions
      if (ELASTIC_SIMULATION) then
        dset_name = "rmassx" ! 1 r (/offset_nglob_xy/)
        call h5_write_dataset_collect_hyperslab(dset_name, rmassx, (/sum(offset_nglob_xy(0:myrank-1))/),if_col)
        dset_name = "rmassy" ! 1 r (/offset_nglob_xy/)
        call h5_write_dataset_collect_hyperslab(dset_name, rmassy, (/sum(offset_nglob_xy(0:myrank-1))/),if_col)
        dset_name = "rmassz" ! 1 r (/offset_nglob_xy/)
        call h5_write_dataset_collect_hyperslab(dset_name, rmassz, (/sum(offset_nglob_xy(0:myrank-1))/),if_col)
     endif
      if (ACOUSTIC_SIMULATION) then
        dset_name = "rmassz_acoustic" ! 1 r (/offset_nglob_xy/)
        call h5_write_dataset_collect_hyperslab(dset_name, rmassz_acoustic, (/sum(offset_nglob_xy(0:myrank-1))/),if_col)
     endif
    endif
  else
    dset_name = "abs_boundary_ispec" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/0/), (/myrank/),if_col)
    dset_name = "abs_boundary_ijk" ! 3 i (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, i3d_dummy, (/0,0,myrank/),if_col)
    dset_name = "abs_boundary_jacobian2Dw" ! 2 r (/0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, r2d_dummy, (/0,myrank/),if_col)
    dset_name = "abs_boundary_normal" ! 3 r (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, r3d_dummy, (/0,0,myrank/),if_col)

    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! store mass matrix contributions
      if (ELASTIC_SIMULATION) then
        dset_name = "rmassx" ! 1 r (/myrank/)
        call h5_write_dataset_collect_hyperslab(dset_name, (/0.0/), (/myrank/),if_col)
        dset_name = "rmassy" ! 1 r (/myrank/)
        call h5_write_dataset_collect_hyperslab(dset_name, (/0.0/), (/myrank/),if_col)
        dset_name = "rmassz" ! 1 r (/myrank/)
        call h5_write_dataset_collect_hyperslab(dset_name, (/0.0/), (/myrank/),if_col)
      endif
      if (ACOUSTIC_SIMULATION) then
        dset_name = "rmassz_acoustic" ! 1 r (/myrank/)
        call h5_write_dataset_collect_hyperslab(dset_name, (/0.0/), (/myrank/),if_col)
      endif
    endif

  endif

  dset_name = "nspec2D_xmin" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/nspec2D_xmin/), (/myrank/),if_col)
  dset_name = "nspec2D_xmax" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/nspec2D_xmax/), (/myrank/),if_col)
  dset_name = "nspec2D_ymin" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/nspec2D_ymin/), (/myrank/),if_col)
  dset_name = "nspec2D_ymax" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/nspec2D_ymax/), (/myrank/),if_col)
  dset_name = "NSPEC2D_BOTTOM" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/nspec2D_bottom/), (/myrank/),if_col)
  dset_name = "NSPEC2D_TOP" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/nspec2D_top/), (/myrank/),if_col)

  if (sum(offset_nspec2D_xmin) > 0) then
    dset_name = "ibelm_xmin" ! 1 i (/offset_nspec2D_xmin/)
    call h5_write_dataset_collect_hyperslab(dset_name, ibelm_xmin, (/sum(offset_nspec2D_xmin(0:myrank-1))/),if_col)
  endif
  if (sum(offset_nspec2D_xmax) > 0) then
    dset_name = "ibelm_xmax" ! 1 i (/offset_nspec2D_xmax/)
    call h5_write_dataset_collect_hyperslab(dset_name, ibelm_xmax, (/sum(offset_nspec2D_xmax(0:myrank-1))/),if_col)
  endif
  if (sum(offset_nspec2D_ymin) > 0) then
    dset_name = "ibelm_ymin" ! 1 i (/offset_nspec2D_ymin/)
    call h5_write_dataset_collect_hyperslab(dset_name, ibelm_ymin, (/sum(offset_nspec2D_ymin(0:myrank-1))/),if_col)
  endif
  if (sum(offset_nspec2D_ymax) > 0) then
    dset_name = "ibelm_ymax" ! 1 i (/offset_nspec2D_ymax/)
    call h5_write_dataset_collect_hyperslab(dset_name, ibelm_ymax, (/sum(offset_nspec2D_ymax(0:myrank-1))/),if_col)
  endif
  if (sum(offset_nspec2D_bottom_ext) > 0) then
    dset_name = "ibelm_bottom" ! 1 i (/offset_nspec2D_bottom_ext/)
    call h5_write_dataset_collect_hyperslab(dset_name, ibelm_bottom, (/sum(offset_nspec2D_bottom_ext(0:myrank-1))/),if_col)
  endif
  if (sum(offset_nspec2D_top_ext) > 0) then
    dset_name = "ibelm_top" ! 1 i (/offset_nspec2D_top_ext/)
    call h5_write_dataset_collect_hyperslab(dset_name, ibelm_top, (/sum(offset_nspec2D_top_ext(0:myrank-1))/),if_col)
  endif

  ! free surface
  dset_name = "num_free_surface_faces" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/num_free_surface_faces/), (/myrank/),if_col)

  if (sum(offset_num_free_surface_faces) > 0) then
    dset_name = "free_surface_ispec" ! 1 i (/offset_num_free_surface_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, free_surface_ispec, &
            (/sum(offset_num_free_surface_faces(0:myrank-1))/),if_col)
    dset_name = "free_surface_ijk" ! 3 i (/0,0,offset_num_free_surface_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, free_surface_ijk, &
            (/0,0,sum(offset_num_free_surface_faces(0:myrank-1))/),if_col)
    dset_name = "free_surface_jacobian2Dw" ! 2 r (/0,offset_num_free_surface_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, free_surface_jacobian2Dw, &
            (/0,sum(offset_num_free_surface_faces(0:myrank-1))/),if_col)
    dset_name = "free_surface_normal" ! 3 r (/0,0,offset_num_free_surface_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, free_surface_normal, &
            (/0,0,sum(offset_num_free_surface_faces(0:myrank-1))/),if_col)
  else
    dset_name = "free_surface_ispec" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/0/), (/myrank/),if_col)
    dset_name = "free_surface_ijk" ! 3 i (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, i3d_dummy, (/0,0,myrank/),if_col)
    dset_name = "free_surface_jacobian2Dw" ! 2 r (/0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, r2d_dummy, (/0,myrank/),if_col)
    dset_name = "free_surface_normal" ! 3 r (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, r3d_dummy, (/0,0,myrank/),if_col)
  endif

  ! acoustic-elastic coupling surface
  dset_name = "num_coupling_ac_el_faces" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/num_coupling_ac_el_faces/), (/myrank/),if_col)

  if (sum(offset_num_coupling_ac_el_faces) > 0) then
    dset_name = "coupling_ac_el_ispec" ! 1 i (/offset_num_coupling_ac_el_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_ac_el_ispec, &
              (/sum(offset_num_coupling_ac_el_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_ac_el_ijk" ! 3 i (/0,0,offset_num_coupling_ac_el_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_ac_el_ijk, &
              (/0,0,sum(offset_num_coupling_ac_el_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_ac_el_jacobian2Dw" ! 2 r (/0,offset_num_coupling_ac_el_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_ac_el_jacobian2Dw, &
             (/0,sum(offset_num_coupling_ac_el_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_ac_el_normal" ! 3 r (/0,0,offset_num_coupling_ac_el_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_ac_el_normal, &
             (/0,0,sum(offset_num_coupling_ac_el_faces(0:myrank-1))/),if_col)
  else
    dset_name = "coupling_ac_el_ispec" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/0/), (/myrank/),if_col)
    dset_name = "coupling_ac_el_ijk" ! 3 i (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, i3d_dummy, (/0,0,myrank/),if_col)
    dset_name = "coupling_ac_el_jacobian2Dw" ! 2 r (/0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, r2d_dummy, (/0,myrank/),if_col)
    dset_name = "coupling_ac_el_normal" ! 3 r (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, r3d_dummy, (/0,0,myrank/),if_col)
  endif

  ! acoustic-poroelastic coupling surface
  dset_name = "num_coupling_ac_po_faces" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/num_coupling_ac_po_faces/), (/myrank/),if_col)

  if (sum(offset_num_coupling_ac_po_faces) > 0) then
    dset_name = "coupling_ac_po_ispec" ! 1 i (/offset_num_coupling_ac_po_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_ac_po_ispec, &
            (/sum(offset_num_coupling_ac_po_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_ac_po_ijk" ! 3 i (/0,0,offset_num_coupling_ac_po_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_ac_po_ijk, &
            (/0,0,sum(offset_num_coupling_ac_po_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_ac_po_jacobian2Dw" ! 2 r (/0,offset_num_coupling_ac_po_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_ac_po_jacobian2Dw, &
            (/0,sum(offset_num_coupling_ac_po_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_ac_po_normal" ! 3 r (/0,0,offset_num_coupling_ac_po_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_ac_po_normal, &
            (/0,0,sum(offset_num_coupling_ac_po_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_ac_po_ispec" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/0/), (/myrank/),if_col)
    dset_name = "coupling_ac_po_ijk" ! 3 i (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, i3d_dummy, (/0,0,myrank/),if_col)
    dset_name = "coupling_ac_po_jacobian2Dw" ! 2 r (/0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, r2d_dummy, (/0,myrank/),if_col)
    dset_name = "coupling_ac_po_normal" ! 3 r (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, r3d_dummy, (/0,0,myrank/),if_col)
  endif

  ! elastic-poroelastic coupling surface
  dset_name = "num_coupling_el_po_faces" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/num_coupling_el_po_faces/), (/myrank/),if_col)

  if (sum(offset_num_coupling_el_po_faces) > 0) then
    dset_name = "coupling_el_po_ispec" ! 1 i (/offset_num_coupling_el_po_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_el_po_ispec, &
            (/sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_po_el_ispec" ! 1 i (/offset_num_coupling_el_po_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_po_el_ispec, &
            (/sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_el_po_ijk" ! 3 i (/0,0,offset_num_coupling_el_po_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_el_po_ijk, &
            (/0,0,sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_po_el_ijk" ! 3 i (/0,0,offset_num_coupling_el_po_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_po_el_ijk, &
            (/0,0,sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_el_po_jacobian2Dw" ! 2 r (/0,offset_num_coupling_el_po_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_el_po_jacobian2Dw, &
            (/0,sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
    dset_name = "coupling_el_po_normal" ! 3 r (/0,0,offset_num_coupling_el_po_faces/)
    call h5_write_dataset_collect_hyperslab(dset_name, coupling_el_po_normal, &
            (/0,0,sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
  else
    dset_name = "coupling_el_po_ispec" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/0/), (/myrank/),if_col)
    dset_name = "coupling_po_el_ispec" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/0/), (/myrank/),if_col)
    dset_name = "coupling_el_po_ijk" ! 3 i (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, i3d_dummy, (/0,0,myrank/),if_col)
    dset_name = "coupling_po_el_ijk" ! 3 i (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, i3d_dummy, (/0,0,myrank/),if_col)
    dset_name = "coupling_el_po_jacobian2Dw" ! 2 r (/0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, r2d_dummy, (/0,myrank/),if_col)
    dset_name = "coupling_el_po_normal" ! 3 r (/0,0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, r3d_dummy, (/0,0,myrank/),if_col)
  endif

  dset_name = "num_interfaces_ext_mesh" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/num_interfaces_ext_mesh/), (/myrank/),if_col)

  if (sum(offset_num_interfaces_ext_mesh) > 0) then
    dset_name = "max_nibool_interfaces_ext_mesh" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/max_nibool_interfaces_ext_mesh/), (/myrank/),if_col)
    dset_name = "my_neighbors_ext_mesh" ! 1 i (/offset_num_interfaces_ext_mesh/)
    call h5_write_dataset_collect_hyperslab(dset_name, my_neighbors_ext_mesh, &
            (/sum(offset_num_interfaces_ext_mesh(0:myrank-1))/),if_col)
    dset_name = "nibool_interfaces_ext_mesh" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, nibool_interfaces_ext_mesh, &
            (/sum(offset_num_interfaces_ext_mesh(0:myrank-1))/),if_col)
    dset_name = "ibool_interfaces_ext_mesh_dummy" ! 2 i (/offset_max_ni_bool_interfaces_ext_mesh/)
    call h5_write_dataset_collect_hyperslab(dset_name, ibool_interfaces_ext_mesh_dummy, &
            (/0,sum(offset_num_interfaces_ext_mesh(0:myrank-1))/),if_col)
  else
    dset_name = "max_nibool_interfaces_ext_mesh" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/0/), (/myrank/),if_col)
    dset_name = "my_neighbors_ext_mesh" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/0/), (/myrank/),if_col)
    dset_name = "nibool_interfaces_ext_mesh" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/0/), (/myrank/),if_col)
    dset_name = "ibool_interfaces_ext_mesh_dummy" ! 2 i (/0,myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, i2d_dummy, (/myrank/),if_col)
  endif

  ! anisotropy
  if (ELASTIC_SIMULATION .and. ANISOTROPY) then
    dset_name = "c11store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c11store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c12store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c12store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c13store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c13store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c14store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c14store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c15store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c15store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c16store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c16store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c22store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c22store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c23store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c23store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c24store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c24store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c25store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c25store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c26store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c26store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c33store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c33store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c34store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c34store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c35store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c35store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c36store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c36store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c44store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c44store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c45store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c45store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c46store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c46store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c55store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c55store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c56store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c56store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    dset_name = "c66store" ! 4 r (/0,0,0,offset_nspec_aniso/)
    call h5_write_dataset_collect_hyperslab(dset_name, c66store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
  endif

  ! inner/outer elements
  dset_name = "ispec_is_inner" ! 1 l (/offset_nspec/)
  call h5_write_dataset_collect_hyperslab(dset_name, ispec_is_inner,(/sum(offset_nspec(0:myrank-1))/),if_col)

  if (ACOUSTIC_SIMULATION) then
     dset_name = "nspec_inner_acoustic" ! 1 i (/myrank/)
     call h5_write_dataset_collect_hyperslab(dset_name, (/nspec_inner_acoustic/),(/myrank/),if_col)
     dset_name = "nspec_outer_acoustic" ! 1 i (/myrank/)
     call h5_write_dataset_collect_hyperslab(dset_name, (/nspec_outer_acoustic/),(/myrank/),if_col)
     dset_name = "num_phase_ispec_acoustic" ! 1 i (/myrank/)
     call h5_write_dataset_collect_hyperslab(dset_name, (/num_phase_ispec_acoustic/),(/myrank/),if_col)
    if (sum(offset_num_phase_ispec_acoustic) > 0) then
      dset_name = "phase_ispec_inner_acoustic" ! 2 i (/offset_num_phase_ispec_acoustic, 0/)
      call h5_write_dataset_collect_hyperslab(dset_name, phase_ispec_inner_acoustic, &
              (/sum(offset_num_phase_ispec_acoustic(0:myrank-1)),0/),if_col)
    else
      dset_name = "phase_ispec_inner_acoustic" ! 2 i (/myrank,0/)
      call h5_write_dataset_collect_hyperslab(dset_name, i2d_dummy,(/myrank,0/),if_col)
    endif
  endif

  if (ELASTIC_SIMULATION) then
    dset_name = "nspec_inner_elastic" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/nspec_inner_elastic/),(/myrank/),if_col)
    dset_name = "nspec_outer_elastic" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/nspec_outer_elastic/),(/myrank/),if_col)
    dset_name = "num_phase_ispec_elastic" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/num_phase_ispec_elastic/),(/myrank/),if_col)
    if (sum(offset_num_phase_ispec_elastic) > 0) then
      dset_name = "phase_ispec_inner_elastic" ! 2 i (/offset_num_phase_ispec_elastic,0/)
      call h5_write_dataset_collect_hyperslab(dset_name, phase_ispec_inner_elastic, &
              (/sum(offset_num_phase_ispec_elastic(0:myrank-1)),0/),if_col)
    else
      dset_name = "phase_ispec_inner_elastic" ! 2 i (/myrank, 0/)
      call h5_write_dataset_collect_hyperslab(dset_name, i2d_dummy,(/myrank,0/),if_col)
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    dset_name = "nspec_inner_poroelastic" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/nspec_inner_poroelastic/),(/myrank/),if_col)
    dset_name = "nspec_outer_poroelastic" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/nspec_outer_poroelastic/),(/myrank/),if_col)
    dset_name = "num_phase_ispec_poroelastic" ! 1 i (/myrank/)
    call h5_write_dataset_collect_hyperslab(dset_name, (/num_phase_ispec_poroelastic/), &
            (/myrank/),if_col)
    if (sum(offset_num_phase_ispec_poroelastic) > 0) then
      dset_name = "phase_ispec_inner_poroelastic" ! 2 i (/offset_num_phase_ispec_poroelastic/)
      call h5_write_dataset_collect_hyperslab(dset_name, phase_ispec_inner_poroelastic, &
            (/sum(offset_num_phase_ispec_poroelastic(0:myrank-1)),0/),if_col)
      dset_name = "phase_ispec_inner_poroelastic" ! 2 i (/myrank,0/)
      call h5_write_dataset_collect_hyperslab(dset_name, i2d_dummy,(/myrank,0/),if_col)
    endif
  endif

  ! mesh coloring
  if (USE_MESH_COLORING_GPU) then
    if (ACOUSTIC_SIMULATION) then
      dset_name = "num_colors_outer_acoustic" ! 1 i (/myrank/)
      call h5_write_dataset_collect_hyperslab(dset_name, &
              (/num_colors_outer_acoustic/),(/myrank/),if_col)
      dset_name = "num_colors_inner_acoustic" ! 1 i (/myrank/)
      call h5_write_dataset_collect_hyperslab(dset_name, &
              (/num_colors_inner_acoustic/),(/myrank/),if_col)
      dset_name = "num_elem_colors_acoustic" ! 1 i (/offset_num_colors_outer_acoustic+num_colors_inner_acoustic/)
      call h5_write_dataset_collect_hyperslab(dset_name, &
              num_elem_colors_acoustic,(/sum(offset_num_colors_acoustic(0:myrank-1))/),if_col)
    endif
    if (ELASTIC_SIMULATION) then
      dset_name = "num_colors_outer_elastic" ! 1 i (/myrank/)
      call h5_write_dataset_collect_hyperslab(dset_name, &
              (/num_colors_outer_elastic/),(/myrank/),if_col)
      dset_name = "num_colors_inner_elastic" ! 1 i (/myrank/)
      call h5_write_dataset_collect_hyperslab(dset_name, &
              (/num_colors_inner_elastic/),(/myrank/),if_col)
      dset_name = "num_elem_colors_elastic" ! 1 i (/offset_num_colors_outer_elastic+num_colors_inner_elastic/)
      call h5_write_dataset_collect_hyperslab(dset_name, num_elem_colors_elastic, &
              (/sum(offset_num_colors_elastic(0:myrank-1))/),if_col)
    endif
  endif

  ! surface points
  dset_name = "nfaces_surface" ! 1 i (/myrank/)
  call h5_write_dataset_collect_hyperslab(dset_name, (/nfaces_surface/), (/myrank/), if_col)
  dset_name = "ispec_is_surface_external_mesh" ! 1 l (/offset_nspec_ab/)
  call h5_write_dataset_collect_hyperslab(dset_name, ispec_is_surface_external_mesh, &
          (/sum(offset_nspec_ab(0:myrank-1))/), if_col)
  dset_name = "iglob_is_surface_external_mesh" ! 1 l (/offset_nglob_ab/)
  call h5_write_dataset_collect_hyperslab(dset_name, iglob_is_surface_external_mesh, &
          (/sum(offset_nglob_ab(0:myrank-1))/), if_col)

  ! arrays for visualization
  dset_name = "spec_elm_conn_xdmf" ! 2 i (/0,offset_nspec*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)/)
  call h5_write_dataset_collect_hyperslab(dset_name, spec_elm_conn_xdmf, &
          (/0,sum(offset_nspec(0:myrank-1))*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)/), if_col)

  ! stores arrays in binary files
  !if (SAVE_MESH_FILES) then
  !  call save_arrays_solver_files()
  !endif

  ! if SAVE_MESH_FILES is true then the files have already been saved, no need to save them again
  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
    call save_arrays_solver_injection_boundary()
  endif

  ! synchronizes processes
  call synchronize_all()

  ! cleanup
  if (allocated(ibool_interfaces_ext_mesh_dummy)) then
    deallocate(ibool_interfaces_ext_mesh_dummy,stat=ier)
    if (ier /= 0) stop 'error deallocating array ibool_interfaces_ext_mesh_dummy'
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "     write mesh dataset in h5 format finished"
    call flush_IMAIN()
  endif

  call h5_close_file_p()
  call h5_finalize()

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) "     all h5 objects closed"
    call flush_IMAIN()
  endif

#else
  ! no HDF5 compilation support
  ! user output
  print *
  print *, "Error: HDF5 routine save_arrays_solver_mesh_hdf5() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *
  stop 'Error HDF5 save_arrays_solver_mesh_hdf5(): called without compilation support'

#endif

  end subroutine save_arrays_solver_mesh_hdf5

!
!-------------------------------------------------------------------------------
!

  subroutine get_connectivity_for_movie(nspec,ibool,elm_conn,o)

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ

  implicit none

  integer, intent(in)                                          :: nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in)       :: ibool
  integer, intent(in) :: o ! node id offset (starting global element id of each proc)
  integer, dimension(9,nspec*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)), intent(inout) :: elm_conn

  integer :: ispec,ii,icub,jcub,kcub
  integer, parameter :: cell_type = 9

  do ispec=1, nspec
    ! extract information from full GLL grid
    ! node order follows vtk format

    do icub=0,NGLLX-2
      do jcub=0,NGLLY-2
        do kcub=0,NGLLZ-2
          ii = 1+(ispec-1)*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1) + (icub*(NGLLY-1)*(NGLLZ-1)+jcub*(NGLLZ-1)+kcub)
          elm_conn(1, ii)  = cell_type
          elm_conn(2, ii)  = ibool(icub+1,jcub+1,kcub+1,ispec)-1+o ! node id starts 0 in xdmf rule
          elm_conn(3, ii)  = ibool(icub+2,jcub+1,kcub+1,ispec)-1+o
          elm_conn(4, ii)  = ibool(icub+2,jcub+2,kcub+1,ispec)-1+o
          elm_conn(5, ii)  = ibool(icub+1,jcub+2,kcub+1,ispec)-1+o
          elm_conn(6, ii)  = ibool(icub+1,jcub+1,kcub+2,ispec)-1+o
          elm_conn(7, ii)  = ibool(icub+2,jcub+1,kcub+2,ispec)-1+o
          elm_conn(8, ii)  = ibool(icub+2,jcub+2,kcub+2,ispec)-1+o
          elm_conn(9, ii)  = ibool(icub+1,jcub+2,kcub+2,ispec)-1+o
        enddo
      enddo
    enddo
  enddo

  end subroutine get_connectivity_for_movie


