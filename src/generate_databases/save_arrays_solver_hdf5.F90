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

! for external mesh

  subroutine save_arrays_solver_ext_mesh_h5(nspec,nglob,APPROXIMATE_OCEAN_LOAD,ibool, &
                    num_interfaces_ext_mesh,my_neighbors_ext_mesh,nibool_interfaces_ext_mesh, &
                    max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                    SAVE_MESH_FILES,ANISOTROPY)

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,IOUT, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
    SIMULATION_TYPE,SAVE_FORWARD,mask_ibool_interior_domain, &
    STACEY_ABSORBING_CONDITIONS,USE_MESH_COLORING_GPU, &
    myrank

  ! PML
  use generate_databases_par, only: PML_CONDITIONS, &
    nspec_cpml,CPML_width_x,CPML_width_y,CPML_width_z,CPML_to_spec, &
    CPML_regions,is_CPML,min_distance_between_CPML_parameter,nspec_cpml_tot, &
    d_store_x,d_store_y,d_store_z,k_store_x,k_store_y,k_store_z, &
    alpha_store_x,alpha_store_y,alpha_store_z, &
    nglob_interface_PML_acoustic,points_interface_PML_acoustic, &
    nglob_interface_PML_elastic,points_interface_PML_elastic

  ! mesh surface
  use generate_databases_par, only: ispec_is_surface_external_mesh,iglob_is_surface_external_mesh, &
    nfaces_surface,nspec_irregular

  use create_regions_mesh_ext_par

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH,LOCAL_PATH,NPROC

  use constants, only: CUSTOM_REAL

  use phdf5_utils

  implicit none

  integer,intent(in) :: nspec,nglob
  ! ocean load
  logical,intent(in) :: APPROXIMATE_OCEAN_LOAD
  ! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  ! MPI interfaces
  integer,intent(in) :: num_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(in) :: my_neighbors_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer,intent(in) :: max_interface_size_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh

  logical,intent(in) :: SAVE_MESH_FILES
  logical,intent(in) :: ANISOTROPY

  ! local parameters
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: max_nibool_interfaces_ext_mesh

  integer :: ier,i,iproc
  character(len=MAX_STRING_LEN) :: filename

  ! mpi variables
  integer :: info, comm

  ! hdf5 valiables
  character(len=64) :: group_name
  character(len=64) :: dset_name
  character(len=10) :: tempstr
  character(len=5)  :: gname_proc_head = "proc_"
  type(h5io)        :: h5

  ! dump dataset size
  integer, dimension(NPROC,4) :: dsize_dump


  ! saves mesh file external_mesh.h5
  prname = "/external_mesh.h5"
  filename = LOCAL_PATH(1:len_trim(LOCAL_PATH))//prname

  h5 = h5io()

  ! get mpi parameters
  call world_get_comm(comm)
  call get_info_null(info)
  ! initialize h5 object
  call h5_init(h5, filename)
  call h5_set_mpi_info(h5, comm, info, myrank, NPROC)


  ! create hdf5 file and write datasets independently
  ! create file, groups and datasets of all mpi rank
  if (myrank == 0) then
    call h5_create_file(h5)
    ! create group
    do iproc = 0, NPROC-1
      write(tempstr, "(i6.6)") iproc
      group_name = gname_proc_head // trim(tempstr)

      call h5_create_group(h5, group_name)
    enddo
    call h5_close_file(h5)
  endif
  call synchronize_all()


  ! create dataset
  if (myrank == 0) then
     call h5_open_file(h5)
     write(tempstr, "(i6.6)") myrank
     group_name = gname_proc_head // trim(tempstr)
     call h5_open_group(h5, group_name)
  endif

  ! create datasets
  dset_name = "nspec"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  dset_name = "nglob"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  dset_name = "nspec_irregular"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  dset_name = "ibool"
  call h5_create_dataset_setter(h5, dset_name, shape(ibool), 4, 1)

  dset_name = "xstore_dummy"
  call h5_create_dataset_setter(h5, dset_name, shape(xstore_dummy), 1, CUSTOM_REAL)
  dset_name = "ystore_dummy"
  call h5_create_dataset_setter(h5, dset_name, shape(ystore_dummy), 1, CUSTOM_REAL)
  dset_name = "zstore_dummy"
  call h5_create_dataset_setter(h5, dset_name, shape(zstore_dummy), 1, CUSTOM_REAL)

  dset_name = "irregular_element_number"
  call h5_create_dataset_setter(h5, dset_name, shape(irregular_element_number), 1, 1)
  dset_name = "xix_regular"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, CUSTOM_REAL)
  dset_name = "jacobian_regular"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, CUSTOM_REAL)

  dset_name = "xixstore"
  call h5_create_dataset_setter(h5, dset_name, shape(xixstore), 4, CUSTOM_REAL)
  dset_name = "xiystore"
  call h5_create_dataset_setter(h5, dset_name, shape(xiystore), 4, CUSTOM_REAL)
  dset_name = "xizstore"
  call h5_create_dataset_setter(h5, dset_name, shape(xizstore), 4, CUSTOM_REAL)
  dset_name = "etaxstore"
  call h5_create_dataset_setter(h5, dset_name, shape(etaxstore), 4, CUSTOM_REAL)
  dset_name = "etaystore"
  call h5_create_dataset_setter(h5, dset_name, shape(etaystore), 4, CUSTOM_REAL)
  dset_name = "etazstore"
  call h5_create_dataset_setter(h5, dset_name, shape(etazstore), 4, CUSTOM_REAL)
  dset_name = "gammaxstore"
  call h5_create_dataset_setter(h5, dset_name, shape(gammaxstore), 4, CUSTOM_REAL)
  dset_name = "gammaystore"
  call h5_create_dataset_setter(h5, dset_name, shape(gammaystore), 4, CUSTOM_REAL)
  dset_name = "gammazstore"
  call h5_create_dataset_setter(h5, dset_name, shape(gammazstore), 4, CUSTOM_REAL)
  dset_name = "jacobianstore"
  call h5_create_dataset_setter(h5, dset_name, shape(jacobianstore), 4, CUSTOM_REAL)

  dset_name = "kappastore"
  call h5_create_dataset_setter(h5, dset_name, shape(kappastore), 4, CUSTOM_REAL)
  dset_name = "mustore"
  call h5_create_dataset_setter(h5, dset_name, shape(mustore), 4, CUSTOM_REAL)

  dset_name = "ispec_is_acoustic"
  call h5_create_dataset_setter(h5, dset_name, shape(ispec_is_acoustic), 1, 1)
  dset_name = "ispec_is_elastic"
  call h5_create_dataset_setter(h5, dset_name, shape(ispec_is_elastic), 1, 1)
  dset_name = "ispec_is_poroelastic"
  call h5_create_dataset_setter(h5, dset_name, shape(ispec_is_poroelastic), 1, 1)
  
  ! acoustic
  if (ACOUSTIC_SIMULATION) then
    dset_name = "rmass_acoustic"
    call h5_create_dataset_setter(h5, dset_name, shape(rmass_acoustic), 1, CUSTOM_REAL)
  endif
  
  ! this array is needed for acoustic simulations but also for elastic simulations with CPML,
  ! thus we allocate it and read it in all cases (whether the simulation is acoustic, elastic, or acoustic/elastic)
  dset_name = "rhostore"
  call h5_create_dataset_setter(h5, dset_name, shape(rhostore), 4, CUSTOM_REAL)
  
  
  ! elastic
  if (ELASTIC_SIMULATION) then
    dset_name = "rmass"
    call h5_create_dataset_setter(h5, dset_name, shape(rmass), 1, CUSTOM_REAL)
    if (APPROXIMATE_OCEAN_LOAD) then
      dset_name = "rmass_ocean_load"
      call h5_create_dataset_setter(h5, dset_name, shape(rmass_ocean_load), 1, CUSTOM_REAL)
    endif
    !pll Stacey
    dset_name = "rho_vp"
    call h5_create_dataset_setter(h5, dset_name, shape(rho_vp), 4, CUSTOM_REAL)
    dset_name = "rho_vs"
    call h5_create_dataset_setter(h5, dset_name, shape(rho_vs), 4, CUSTOM_REAL)
  endif
  
  ! poroelastic
  if (POROELASTIC_SIMULATION) then
    dset_name = "rmass_solid_poroelastic"
    call h5_create_dataset_setter(h5, dset_name, shape(rmass_solid_poroelastic), 1, CUSTOM_REAL)
    dset_name = "rmass_fluid_poroelastic"
    call h5_create_dataset_setter(h5, dset_name, shape(rmass_fluid_poroelastic), 1, CUSTOM_REAL)
    dset_name = "rhoarraystore"
    call h5_create_dataset_setter(h5, dset_name, shape(rhoarraystore), 5, CUSTOM_REAL)
    dset_name = "kappaarraystore"
    call h5_create_dataset_setter(h5, dset_name, shape(kappaarraystore), 5, CUSTOM_REAL)
    dset_name = "etastore"
    call h5_create_dataset_setter(h5, dset_name, shape(etastore), 4, CUSTOM_REAL)
    dset_name = "tortstore"
    call h5_create_dataset_setter(h5, dset_name, shape(tortstore), 4, CUSTOM_REAL)
    dset_name = "permstore"
    call h5_create_dataset_setter(h5, dset_name, shape(permstore), 5, CUSTOM_REAL)
    dset_name = "phistore"
    call h5_create_dataset_setter(h5, dset_name, shape(phistore), 4, CUSTOM_REAL)
    dset_name = "rho_vpI"
    call h5_create_dataset_setter(h5, dset_name, shape(rho_vpI), 4, CUSTOM_REAL)
    dset_name = "rho_vpII"
    call h5_create_dataset_setter(h5, dset_name, shape(rho_vpII), 4, CUSTOM_REAL)
    dset_name = "rho_vsI"
    call h5_create_dataset_setter(h5, dset_name, shape(rho_vsI), 4, CUSTOM_REAL)
  endif

  ! C-PML absorbing boundary conditions
  if (PML_CONDITIONS) then
    dset_name = "nspec_cpml"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
    dset_name = "CPML_width_x"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, CUSTOM_REAL)
    dset_name = "CPML_width_y"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, CUSTOM_REAL)
    dset_name = "CPML_width_z"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, CUSTOM_REAL)
    dset_name = "min_distance_between_CPML_parameter"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, CUSTOM_REAL)
    if (nspec_cpml > 0) then
      dset_name = "CPML_regions"
      call h5_create_dataset_setter(h5, dset_name, shape(CPML_regions), 1, 1)
      dset_name = "CPML_to_spec"
      call h5_create_dataset_setter(h5, dset_name, shape(CPML_to_spec), 1, 1)
      dset_name = "is_CPML"
      call h5_create_dataset_setter(h5, dset_name, shape(is_CPML), 1, 1)
      dset_name = "d_store_x"
      call h5_create_dataset_setter(h5, dset_name, shape(d_store_x), 4, CUSTOM_REAL)
      dset_name = "d_store_y"
      call h5_create_dataset_setter(h5, dset_name, shape(d_store_y), 4, CUSTOM_REAL)
      dset_name = "d_store_z"
      call h5_create_dataset_setter(h5, dset_name, shape(d_store_z), 4, CUSTOM_REAL)
      dset_name = "k_store_x"
      call h5_create_dataset_setter(h5, dset_name, shape(k_store_x), 4, CUSTOM_REAL)
      dset_name = "k_store_y"
      call h5_create_dataset_setter(h5, dset_name, shape(k_store_y), 4, CUSTOM_REAL)
      dset_name = "k_store_z"
      call h5_create_dataset_setter(h5, dset_name, shape(k_store_z), 4, CUSTOM_REAL)
      dset_name = "alpha_store_x"
      call h5_create_dataset_setter(h5, dset_name, shape(alpha_store_x), 4, CUSTOM_REAL)
      dset_name = "alpha_store_y"
      call h5_create_dataset_setter(h5, dset_name, shape(alpha_store_y), 4, CUSTOM_REAL)
      dset_name = "alpha_store_z"
      call h5_create_dataset_setter(h5, dset_name, shape(alpha_store_z), 4, CUSTOM_REAL)
      ! --------------------------------------------------------------------------------------------
      ! for adjoint tomography
      ! save the array stored the points on interface between PML and interior computational domain
      ! --------------------------------------------------------------------------------------------
      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        dset_name = "nglob_interface_PML_acoustic"
        call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
        dset_name = "nglob_interface_PML_elastic"
        call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
        if (nglob_interface_PML_acoustic > 0) then
          dset_name = "points_interface_PML_acoustic"
          call h5_create_dataset_setter(h5, dset_name, shape(points_interface_PML_acoustic), 1, 1)
        endif
        if (nglob_interface_PML_elastic > 0) then
          dset_name = "points_interface_PML_elastic"
          call h5_create_dataset_setter(h5, dset_name, shape(points_interface_PML_elastic), 1, 1)
        endif
      endif
    endif
  endif

  ! absorbing boundary surface
  write(IOUT) num_abs_boundary_faces
  if (num_abs_boundary_faces > 0) then
    dset_name = "abs_boundary_ispec"
    call h5_create_dataset_setter(h5, dset_name, shape(abs_boundary_ispec), 1, 1)
    dset_name = "abs_boundary_ijk"
    call h5_create_dataset_setter(h5, dset_name, shape(abs_boundary_ijk), 3, 1)
    dset_name = "abs_boundary_jacobian2Dw"
    call h5_create_dataset_setter(h5, dset_name, shape(abs_boundary_jacobian2Dw), 2, CUSTOM_REAL)
    dset_name = "abs_boundary_normal"
    call h5_create_dataset_setter(h5, dset_name, shape(abs_boundary_normal), 3, CUSTOM_REAL)

    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! store mass matrix contributions
      if (ELASTIC_SIMULATION) then
        dset_name = "rmassx"
        call h5_create_dataset_setter(h5, dset_name, shape(rmassx), 1, CUSTOM_REAL)
        dset_name = "rmassy"
        call h5_create_dataset_setter(h5, dset_name, shape(rmassy), 1, CUSTOM_REAL)
        dset_name = "rmassz"
        call h5_create_dataset_setter(h5, dset_name, shape(rmassz), 1, CUSTOM_REAL)
     endif
      if (ACOUSTIC_SIMULATION) then
        dset_name = "rmassz_acoustic"
        call h5_create_dataset_setter(h5, dset_name, shape(rmassz_acoustic), 1, CUSTOM_REAL)
     endif
    endif
  endif

  dset_name = "nspec2D_xmin"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  dset_name = "nspec2D_xmax"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  dset_name = "nspec2D_ymin"
  call h5_create_dataset_setter(h5, dset_name, (/1/),1 , 1)
  dset_name = "nspec2D_ymax"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  dset_name = "NSPEC2D_BOTTOM"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  dset_name = "NSPEC2D_TOP"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  dset_name = "ibelm_xmin"
  call h5_create_dataset_setter(h5, dset_name, shape(ibelm_xmin), 1, 1)
  dset_name = "ibelm_xmax"
  call h5_create_dataset_setter(h5, dset_name, shape(ibelm_xmax), 1, 1)
  dset_name = "ibelm_ymin"
  call h5_create_dataset_setter(h5, dset_name, shape(ibelm_ymin), 1, 1)
  dset_name = "ibelm_ymax"
  call h5_create_dataset_setter(h5, dset_name, shape(ibelm_ymax), 1, 1)
  dset_name = "ibelm_bottom"
  call h5_create_dataset_setter(h5, dset_name, shape(ibelm_bottom), 1, 1)
  dset_name = "ibelm_top"
  call h5_create_dataset_setter(h5, dset_name, shape(ibelm_top), 1, 1)

  ! free surface
  write(IOUT) num_free_surface_faces
  if (num_free_surface_faces > 0) then
    dset_name = "free_surface_ispec"
    call h5_create_dataset_setter(h5, dset_name, shape(free_surface_ispec), 1, 1)
    dset_name = "free_surface_ijk"
    call h5_create_dataset_setter(h5, dset_name, shape(free_surface_ijk), 3, 1)
    dset_name = "free_surface_jacobian2Dw"
    call h5_create_dataset_setter(h5, dset_name, shape(free_surface_jacobian2Dw), 2, CUSTOM_REAL)
    dset_name = "free_surface_normal"
    call h5_create_dataset_setter(h5, dset_name, shape(free_surface_normal), 3, CUSTOM_REAL)
  endif

  ! acoustic-elastic coupling surface
  write(IOUT) num_coupling_ac_el_faces
  if (num_coupling_ac_el_faces > 0) then
    dset_name = "coupling_ac_el_ispec"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_ac_el_ispec), 1, 1)
    dset_name = "coupling_ac_el_ijk"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_ac_el_ijk), 3, 1)
    dset_name = "coupling_ac_el_jacobian2Dw"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_ac_el_jacobian2Dw), 2, CUSTOM_REAL)
    dset_name = "coupling_ac_el_normal"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_ac_el_normal), 3, CUSTOM_REAL)
  endif
  
  ! acoustic-poroelastic coupling surface
  dset_name = "num_coupling_ac_po_faces"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  if (num_coupling_ac_po_faces > 0) then
    dset_name = "coupling_ac_po_ispec"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_ac_po_ispec), 1, 1)
    dset_name = "coupling_ac_po_ijk"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_ac_po_ijk), 3, 1)
    dset_name = "coupling_ac_po_jacobian2Dw"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_ac_po_jacobian2Dw), 2, CUSTOM_REAL)
    dset_name = "coupling_ac_po_normal"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_ac_po_normal), 3, CUSTOM_REAL)
  endif

  ! elastic-poroelastic coupling surface
  dset_name = "num_coupling_el_po_faces"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  if (num_coupling_el_po_faces > 0) then
    dset_name = "coupling_el_po_ispec"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_el_po_ispec), 1, 1)
    dset_name = "coupling_po_el_ispec"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_po_el_ispec), 1, 1)
    dset_name = "coupling_el_po_ijk"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_el_po_ijk), 3, 1)
    dset_name = "coupling_po_el_ijk"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_po_el_ijk), 3, 1)
    dset_name = "coupling_el_po_jacobian2Dw"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_el_po_jacobian2Dw), 2, CUSTOM_REAL)
    dset_name = "coupling_el_po_normal"
    call h5_create_dataset_setter(h5, dset_name, shape(coupling_el_po_normal), 3, CUSTOM_REAL)
  endif

  !MPI interfaces
  max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh(:))

  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 650')
  if (ier /= 0) stop 'error allocating array'
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,i)
  enddo

  dset_name = "num_interfaces_ext_mesh"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  if (num_interfaces_ext_mesh > 0) then
    dset_name = "max_nibool_interfaces_ext_mesh"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
    dset_name = "my_neighbors_ext_mesh"
    call h5_create_dataset_setter(h5, dset_name, shape(my_neighbors_ext_mesh), 1, 1)
    dset_name = "nibool_interfaces_ext_mesh"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
    dset_name = "ibool_interfaces_ext_mesh_dummy"
    call h5_create_dataset_setter(h5, dset_name, shape(ibool_interfaces_ext_mesh_dummy), 2, 1)
  endif
  
  ! anisotropy
  if (ELASTIC_SIMULATION .and. ANISOTROPY) then
    dset_name = "c11store"
    call h5_create_dataset_setter(h5, dset_name, shape(c11store), 4, CUSTOM_REAL)
    dset_name = "c12store"
    call h5_create_dataset_setter(h5, dset_name, shape(c12store), 4, CUSTOM_REAL)
    dset_name = "c13store"
    call h5_create_dataset_setter(h5, dset_name, shape(c13store), 4, CUSTOM_REAL)
    dset_name = "c14store"
    call h5_create_dataset_setter(h5, dset_name, shape(c14store), 4, CUSTOM_REAL)
    dset_name = "c15store"
    call h5_create_dataset_setter(h5, dset_name, shape(c15store), 4, CUSTOM_REAL)
    dset_name = "c16store"
    call h5_create_dataset_setter(h5, dset_name, shape(c16store), 4, CUSTOM_REAL)
    dset_name = "c22store"
    call h5_create_dataset_setter(h5, dset_name, shape(c22store), 4, CUSTOM_REAL)
    dset_name = "c23store"
    call h5_create_dataset_setter(h5, dset_name, shape(c23store), 4, CUSTOM_REAL)
    dset_name = "c24store"
    call h5_create_dataset_setter(h5, dset_name, shape(c24store), 4, CUSTOM_REAL)
    dset_name = "c25store"
    call h5_create_dataset_setter(h5, dset_name, shape(c25store), 4, CUSTOM_REAL)
    dset_name = "c26store"
    call h5_create_dataset_setter(h5, dset_name, shape(c26store), 4, CUSTOM_REAL)
    dset_name = "c33store"
    call h5_create_dataset_setter(h5, dset_name, shape(c33store), 4, CUSTOM_REAL)
    dset_name = "c34store"
    call h5_create_dataset_setter(h5, dset_name, shape(c34store), 4, CUSTOM_REAL)
    dset_name = "c35store"
    call h5_create_dataset_setter(h5, dset_name, shape(c35store), 4, CUSTOM_REAL)
    dset_name = "c36store"
    call h5_create_dataset_setter(h5, dset_name, shape(c36store), 4, CUSTOM_REAL)

    dset_name = "c44store"
    call h5_create_dataset_setter(h5, dset_name, shape(c44store), 4, CUSTOM_REAL)
    dset_name = "c45store"
    call h5_create_dataset_setter(h5, dset_name, shape(c45store), 4, CUSTOM_REAL)
    dset_name = "c46store"
    call h5_create_dataset_setter(h5, dset_name, shape(c46store), 4, CUSTOM_REAL)
    dset_name = "c55store"
    call h5_create_dataset_setter(h5, dset_name, shape(c55store), 4, CUSTOM_REAL)
    dset_name = "c56store"
    call h5_create_dataset_setter(h5, dset_name, shape(c56store), 4, CUSTOM_REAL)
    dset_name = "c66store"
    call h5_create_dataset_setter(h5, dset_name, shape(c66store), 4, CUSTOM_REAL)
  endif

  ! inner/outer elements
  dset_name = "ispec_is_inner"
  call h5_create_dataset_setter(h5, dset_name, shape(ispec_is_inner), 1, 1)

  if (ACOUSTIC_SIMULATION) then
     dset_name = "nspec_inner_acoustic"
     call h5_create_dataset_setter(h5, dset_name, (/1/), 1 ,1)
     dset_name = "nspec_outer_acoustic"
     call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
     dset_name = "num_phase_ispec_acoustic"
     call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
    if (num_phase_ispec_acoustic > 0) then
      dset_name = "phase_ispec_inner_acoustic"
      call h5_create_dataset_setter(h5, dset_name, shape(phase_ispec_inner_acoustic), 2, 1)
    endif
  endif

  if (ELASTIC_SIMULATION) then
    dset_name = "nspec_inner_elastic"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
    dset_name = "nspec_outer_elastic"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
    dset_name = "num_phase_ispec_elastic"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
    if (num_phase_ispec_elastic > 0) then
      dset_name = "phase_ispec_inner_elastic"
      call h5_create_dataset_setter(h5, dset_name, shape(phase_ispec_inner_elastic), 2, 1)
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    dset_name = "nspec_inner_poroelastic"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
    dset_name = "nspec_outer_poroelastic"
    call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
    dset_name = "num_phase_ispec_poroelastic"
    if (num_phase_ispec_poroelastic > 0) then
      dset_name = "phase_ispec_inner_poroelastic"
      call h5_create_dataset_setter(h5, dset_name, shape(phase_ispec_inner_poroelastic), 2, 1)
    endif
  endif

  ! mesh coloring
  if (USE_MESH_COLORING_GPU) then
    if (ACOUSTIC_SIMULATION) then
      dset_name = "num_colors_outer_acoustic"
      call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
      dset_name = "num_colors_inner_acoustic"
      call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
      dset_name = "num_elem_colors_acoustic"
      call h5_create_dataset_setter(h5, dset_name, shape(num_elem_colors_acoustic), 1, 1)
    endif
    if (ELASTIC_SIMULATION) then
      dset_name = "num_colors_outer_elastic"
      call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
      dset_name = "num_colors_inner_elastic"
      call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
      dset_name = "num_elem_colors_elastic"
      call h5_create_dataset_setter(h5, dset_name, shape(num_elem_colors_elastic), 1, 1)
    endif
  endif

  ! surface points
  dset_name = "nfaces_surface"
  call h5_create_dataset_setter(h5, dset_name, (/1/), 1, 1)
  dset_name = "ispec_is_surface_external_mesh"
  call h5_create_dataset_setter(h5, dset_name, shape(ispec_is_surface_external_mesh), 1, 1)
  dset_name = "iglob_is_surface_external_mesh"
  call h5_create_dataset_setter(h5, dset_name, shape(iglob_is_surface_external_mesh), 1, 1)

  print *, "dataset pre-setting finished"

  call synchronize_all()
 

!
! write dataset data
!

   !call synchronize_all()
   if (myrank == 0) then
    call h5_close_group(h5)
    call h5_close_file(h5)
   endif
   call synchronize_all()
 
!
!
  print*, "start writing datasets"
  call h5_open_file_p(h5)
  write(tempstr, "(i6.6)") myrank
  group_name = gname_proc_head // trim(tempstr)
  print *, "opening group: ", group_name
  call h5_open_group(h5, group_name)
!
!
  dset_name = "nspec"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec/))
  dset_name = "nglob"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/nglob/))

  dset_name = "nspec_irregular"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec_irregular/))
!
  dset_name = "ibool"
  call h5_write_dataset_p_4d_i(h5, dset_name, ibool)

  dset_name = "xstore_dummy"
  call h5_write_dataset_p_1d_r(h5, dset_name, xstore_dummy)
  dset_name = "ystore_dummy"
  call h5_write_dataset_p_1d_r(h5, dset_name, ystore_dummy)
  dset_name = "zstore_dummy"
  call h5_write_dataset_p_1d_r(h5, dset_name, zstore_dummy)

  dset_name = "irregular_element_number"
  call h5_write_dataset_p_1d_i(h5, dset_name, irregular_element_number)
  dset_name = "xix_regular"
  call h5_write_dataset_p_1d_r(h5, dset_name, (/xix_regular/))
  dset_name = "jacobian_regular"
  call h5_write_dataset_p_1d_r(h5, dset_name, (/jacobian_regular/))

  dset_name = "xixstore"
  call h5_write_dataset_p_4d_r(h5, dset_name, xixstore)
  dset_name = "xiystore"
  call h5_write_dataset_p_4d_r(h5, dset_name, xiystore)
  dset_name = "xizstore"
  call h5_write_dataset_p_4d_r(h5, dset_name, xizstore)
  dset_name = "etaxstore"
  call h5_write_dataset_p_4d_r(h5, dset_name, etaxstore)
  dset_name = "etaystore"
  call h5_write_dataset_p_4d_r(h5, dset_name, etaystore)
  dset_name = "etazstore"
  call h5_write_dataset_p_4d_r(h5, dset_name, etazstore)
  dset_name = "gammaxstore"
  call h5_write_dataset_p_4d_r(h5, dset_name, gammaxstore)
  dset_name = "gammaystore"
  call h5_write_dataset_p_4d_r(h5, dset_name, gammaystore)
  dset_name = "gammazstore"
  call h5_write_dataset_p_4d_r(h5, dset_name, gammazstore)
  dset_name = "jacobianstore"
  call h5_write_dataset_p_4d_r(h5, dset_name, jacobianstore)

  dset_name = "kappastore"
  call h5_write_dataset_p_4d_r(h5, dset_name, kappastore)
  dset_name = "mustore"
  call h5_write_dataset_p_4d_r(h5, dset_name, mustore)

  dset_name = "ispec_is_acoustic"
  call h5_write_dataset_p_1d_l(h5, dset_name, ispec_is_acoustic)
  dset_name = "ispec_is_elastic"
  call h5_write_dataset_p_1d_l(h5, dset_name, ispec_is_elastic)
  dset_name = "ispec_is_poroelastic"
  call h5_write_dataset_p_1d_l(h5, dset_name, ispec_is_poroelastic)
  
  ! acoustic
  if (ACOUSTIC_SIMULATION) then
    dset_name = "rmass_acoustic"
    call h5_write_dataset_p_1d_r(h5, dset_name, rmass_acoustic)
  endif
  
  ! this array is needed for acoustic simulations but also for elastic simulations with CPML,
  ! thus we allocate it and read it in all cases (whether the simulation is acoustic, elastic, or acoustic/elastic)
  dset_name = "rhostore"
  call h5_write_dataset_p_4d_r(h5, dset_name, rhostore)
  
  
  ! elastic
  if (ELASTIC_SIMULATION) then
    dset_name = "rmass"
    call h5_write_dataset_p_1d_r(h5, dset_name, rmass)
    if (APPROXIMATE_OCEAN_LOAD) then
      dset_name = "rmass_ocean_load"
      call h5_write_dataset_p_1d_r(h5, dset_name, rmass_ocean_load)
    endif
    !pll Stacey
    dset_name = "rho_vp"
    call h5_write_dataset_p_4d_r(h5, dset_name, rho_vp)
    dset_name = "rho_vs"
    call h5_write_dataset_p_4d_r(h5, dset_name, rho_vs)
  endif
  
  ! poroelastic
  if (POROELASTIC_SIMULATION) then
    dset_name = "rmass_solid_poroelastic"
    call h5_write_dataset_p_1d_r(h5, dset_name, rmass_solid_poroelastic)
    dset_name = "rmass_fluid_poroelastic"
    call h5_write_dataset_p_1d_r(h5, dset_name, rmass_fluid_poroelastic)
    dset_name = "rhoarraystore"
    call h5_write_dataset_p_5d_r(h5, dset_name, rhoarraystore)
    dset_name = "kappaarraystore"
    call h5_write_dataset_p_5d_r(h5, dset_name, kappaarraystore)
    dset_name = "etastore"
    call h5_write_dataset_p_4d_r(h5, dset_name, etastore)
    dset_name = "tortstore"
    call h5_write_dataset_p_4d_r(h5, dset_name, tortstore)
    dset_name = "permstore"
    call h5_write_dataset_p_5d_r(h5, dset_name, permstore)
    dset_name = "phistore"
    call h5_write_dataset_p_4d_r(h5, dset_name, phistore)
    dset_name = "rho_vpI"
    call h5_write_dataset_p_4d_r(h5, dset_name, rho_vpI)
    dset_name = "rho_vpII"
    call h5_write_dataset_p_4d_r(h5, dset_name, rho_vpII)
    dset_name = "rho_vsI"
    call h5_write_dataset_p_4d_r(h5, dset_name, rho_vsI)
  endif

  ! C-PML absorbing boundary conditions
  if (PML_CONDITIONS) then
    dset_name = "nspec_cpml"
    call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec_cpml/))
    dset_name = "CPML_width_x"
    call h5_write_dataset_p_1d_r(h5, dset_name, (/CPML_width_x/))
    dset_name = "CPML_width_y"
    call h5_write_dataset_p_1d_r(h5, dset_name, (/CPML_width_y/))
    dset_name = "CPML_width_z"
    call h5_write_dataset_p_1d_r(h5, dset_name, (/CPML_width_z/))
    dset_name = "min_distance_between_CPML_parameter"
    call h5_write_dataset_p_1d_r(h5, dset_name, (/min_distance_between_CPML_parameter/))
    if (nspec_cpml > 0) then
      dset_name = "CPML_regions"
      call h5_write_dataset_p_1d_i(h5, dset_name, CPML_regions)
      dset_name = "CPML_to_spec"
      call h5_write_dataset_p_1d_i(h5, dset_name, CPML_to_spec)
      dset_name = "is_CPML"
      call h5_write_dataset_p_1d_l(h5, dset_name, is_CPML)
      dset_name = "d_store_x"
      call h5_write_dataset_p_4d_r(h5, dset_name, d_store_x)
      dset_name = "d_store_y"
      call h5_write_dataset_p_4d_r(h5, dset_name, d_store_y)
      dset_name = "d_store_z"
      call h5_write_dataset_p_4d_r(h5, dset_name, d_store_z)
      dset_name = "k_store_x"
      call h5_write_dataset_p_4d_r(h5, dset_name, k_store_x)
      dset_name = "k_store_y"
      call h5_write_dataset_p_4d_r(h5, dset_name, k_store_y)
      dset_name = "k_store_z"
      call h5_write_dataset_p_4d_r(h5, dset_name, k_store_z)
      dset_name = "alpha_store_x"
      call h5_write_dataset_p_4d_r(h5, dset_name, alpha_store_x)
      dset_name = "alpha_store_y"
      call h5_write_dataset_p_4d_r(h5, dset_name, alpha_store_y)
      dset_name = "alpha_store_z"
      call h5_write_dataset_p_4d_r(h5, dset_name, alpha_store_z)
      ! --------------------------------------------------------------------------------------------
      ! for adjoint tomography
      ! save the array stored the points on interface between PML and interior computational domain
      ! --------------------------------------------------------------------------------------------
      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        dset_name = "nglob_interface_PML_acoustic"
        call h5_write_dataset_p_1d_i(h5, dset_name, (/nglob_interface_PML_acoustic/))
        dset_name = "nglob_interface_PML_elastic"
        call h5_write_dataset_p_1d_i(h5, dset_name, (/nglob_interface_PML_elastic/))
        if (nglob_interface_PML_acoustic > 0) then
          dset_name = "points_interface_PML_acoustic"
          call h5_write_dataset_p_1d_i(h5, dset_name, points_interface_PML_acoustic)
        endif
        if (nglob_interface_PML_elastic > 0) then
          dset_name = "points_interface_PML_elastic"
          call h5_write_dataset_p_1d_i(h5, dset_name, points_interface_PML_elastic)
        endif
      endif
    endif
  endif

  ! absorbing boundary surface
  write(IOUT) num_abs_boundary_faces
  if (num_abs_boundary_faces > 0) then
    dset_name = "abs_boundary_ispec"
    call h5_write_dataset_p_1d_i(h5, dset_name, abs_boundary_ispec)
    dset_name = "abs_boundary_ijk"
    call h5_write_dataset_p_3d_i(h5, dset_name, abs_boundary_ijk)
    dset_name = "abs_boundary_jacobian2Dw"
    call h5_write_dataset_p_2d_r(h5, dset_name, abs_boundary_jacobian2Dw)
    dset_name = "abs_boundary_normal"
    call h5_write_dataset_p_3d_r(h5, dset_name, abs_boundary_normal)

    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! store mass matrix contributions
      if (ELASTIC_SIMULATION) then
        dset_name = "rmassx"
        call h5_write_dataset_p_1d_r(h5, dset_name, rmassx)
        dset_name = "rmassy"
        call h5_write_dataset_p_1d_r(h5, dset_name, rmassy)
        dset_name = "rmassz"
        call h5_write_dataset_p_1d_r(h5, dset_name, rmassz)
     endif
      if (ACOUSTIC_SIMULATION) then
        dset_name = "rmassz_acoustic"
        call h5_write_dataset_p_1d_r(h5, dset_name, rmassz_acoustic)
     endif
    endif
  endif

  dset_name = "nspec2D_xmin"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec2D_xmin/))
  dset_name = "nspec2D_xmax"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec2D_xmax/))
  dset_name = "nspec2D_ymin"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec2D_ymin/))
  dset_name = "nspec2D_ymax"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec2D_ymax/))
  dset_name = "NSPEC2D_BOTTOM"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/NSPEC2D_BOTTOM/))
  dset_name = "NSPEC2D_TOP"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/NSPEC2D_TOP/))
  dset_name = "ibelm_xmin"
  call h5_write_dataset_p_1d_i(h5, dset_name, ibelm_xmin)
  dset_name = "ibelm_xmax"
  call h5_write_dataset_p_1d_i(h5, dset_name, ibelm_xmax)
  dset_name = "ibelm_ymin"
  call h5_write_dataset_p_1d_i(h5, dset_name, ibelm_ymin)
  dset_name = "ibelm_ymax"
  call h5_write_dataset_p_1d_i(h5, dset_name, ibelm_ymax)
  dset_name = "ibelm_bottom"
  call h5_write_dataset_p_1d_i(h5, dset_name, ibelm_bottom)
  dset_name = "ibelm_top"
  call h5_write_dataset_p_1d_i(h5, dset_name, ibelm_top)

  ! free surface
  write(IOUT) num_free_surface_faces
  if (num_free_surface_faces > 0) then
    dset_name = "free_surface_ispec"
    call h5_write_dataset_p_1d_i(h5, dset_name, free_surface_ispec)
    dset_name = "free_surface_ijk"
    call h5_write_dataset_p_3d_i(h5, dset_name, free_surface_ijk)
    dset_name = "free_surface_jacobian2Dw"
    call h5_write_dataset_p_2d_r(h5, dset_name, free_surface_jacobian2Dw)
    dset_name = "free_surface_normal"
    call h5_write_dataset_p_3d_r(h5, dset_name, free_surface_normal)
  endif

  ! acoustic-elastic coupling surface
  write(IOUT) num_coupling_ac_el_faces
  if (num_coupling_ac_el_faces > 0) then
    dset_name = "coupling_ac_el_ispec"
    call h5_write_dataset_p_1d_i(h5, dset_name, coupling_ac_el_ispec)
    dset_name = "coupling_ac_el_ijk"
    call h5_write_dataset_p_3d_i(h5, dset_name, coupling_ac_el_ijk)
    dset_name = "coupling_ac_el_jacobian2Dw"
    call h5_write_dataset_p_2d_r(h5, dset_name, coupling_ac_el_jacobian2Dw)
    dset_name = "coupling_ac_el_normal"
    call h5_write_dataset_p_3d_r(h5, dset_name, coupling_ac_el_normal)
  endif
  
  ! acoustic-poroelastic coupling surface
  dset_name = "num_coupling_ac_po_faces"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/num_coupling_ac_po_faces/))
  if (num_coupling_ac_po_faces > 0) then
    dset_name = "coupling_ac_po_ispec"
    call h5_write_dataset_p_1d_i(h5, dset_name, coupling_ac_po_ispec)
    dset_name = "coupling_ac_po_ijk"
    call h5_write_dataset_p_3d_i(h5, dset_name, coupling_ac_po_ijk)
    dset_name = "coupling_ac_po_jacobian2Dw"
    call h5_write_dataset_p_2d_r(h5, dset_name, coupling_ac_po_jacobian2Dw)
    dset_name = "coupling_ac_po_normal"
    call h5_write_dataset_p_3d_r(h5, dset_name, coupling_ac_po_normal)
  endif

  ! elastic-poroelastic coupling surface
  dset_name = "num_coupling_el_po_faces"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/num_coupling_el_po_faces/))
  if (num_coupling_el_po_faces > 0) then
    dset_name = "coupling_el_po_ispec"
    call h5_write_dataset_p_1d_i(h5, dset_name, coupling_el_po_ispec)
    dset_name = "coupling_po_el_ispec"
    call h5_write_dataset_p_1d_i(h5, dset_name, coupling_po_el_ispec)
    dset_name = "coupling_el_po_ijk"
    call h5_write_dataset_p_3d_i(h5, dset_name, coupling_el_po_ijk)
    dset_name = "coupling_po_el_ijk"
    call h5_write_dataset_p_3d_i(h5, dset_name, coupling_po_el_ijk)
    dset_name = "coupling_el_po_jacobian2Dw"
    call h5_write_dataset_p_2d_r(h5, dset_name, coupling_el_po_jacobian2Dw)
    dset_name = "coupling_el_po_normal"
    call h5_write_dataset_p_3d_r(h5, dset_name, coupling_el_po_normal)
  endif

  !MPI interfaces

  dset_name = "num_interfaces_ext_mesh"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/num_interfaces_ext_mesh/))
  if (num_interfaces_ext_mesh > 0) then
    dset_name = "max_nibool_interfaces_ext_mesh"
    call h5_write_dataset_p_1d_i(h5, dset_name, (/max_nibool_interfaces_ext_mesh/))
    dset_name = "my_neighbors_ext_mesh"
    call h5_write_dataset_p_1d_i(h5, dset_name, my_neighbors_ext_mesh)
    dset_name = "nibool_interfaces_ext_mesh"
    call h5_write_dataset_p_1d_i(h5, dset_name, (/nibool_interfaces_ext_mesh/))
    dset_name = "ibool_interfaces_ext_mesh_dummy"
    call h5_write_dataset_p_2d_i(h5, dset_name, ibool_interfaces_ext_mesh_dummy)
  endif
  
  ! anisotropy
  if (ELASTIC_SIMULATION .and. ANISOTROPY) then
    dset_name = "c11store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c11store)
    dset_name = "c12store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c12store)
    dset_name = "c13store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c13store)
    dset_name = "c14store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c14store)
    dset_name = "c15store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c15store)
    dset_name = "c16store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c16store)
    dset_name = "c22store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c22store)
    dset_name = "c23store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c23store)
    dset_name = "c24store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c24store)
    dset_name = "c25store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c25store)
    dset_name = "c26store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c26store)
    dset_name = "c33store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c33store)
    dset_name = "c34store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c34store)
    dset_name = "c35store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c35store)
    dset_name = "c36store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c36store)
    dset_name = "c44store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c44store)
    dset_name = "c45store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c45store)
    dset_name = "c46store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c46store)
    dset_name = "c55store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c55store)
    dset_name = "c56store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c56store)
    dset_name = "c66store"
    call h5_write_dataset_p_4d_r(h5, dset_name, c66store)
  endif

  ! inner/outer elements
  dset_name = "ispec_is_inner"
  call h5_write_dataset_p_1d_l(h5, dset_name, ispec_is_inner)

  if (ACOUSTIC_SIMULATION) then
     dset_name = "nspec_inner_acoustic"
     call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec_inner_acoustic/))
     dset_name = "nspec_outer_acoustic"
     call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec_outer_acoustic/))
     dset_name = "num_phase_ispec_acoustic"
     call h5_write_dataset_p_1d_i(h5, dset_name, (/num_phase_ispec_acoustic/))
    if (num_phase_ispec_acoustic > 0) then
      dset_name = "phase_ispec_inner_acoustic"
      call h5_write_dataset_p_2d_i(h5, dset_name, phase_ispec_inner_acoustic)
    endif
  endif

  if (ELASTIC_SIMULATION) then
    dset_name = "nspec_inner_elastic"
    call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec_inner_elastic/))
    dset_name = "nspec_outer_elastic"
    call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec_outer_elastic/))
    dset_name = "num_phase_ispec_elastic"
    call h5_write_dataset_p_1d_i(h5, dset_name, (/num_phase_ispec_elastic/))
    if (num_phase_ispec_elastic > 0) then
      dset_name = "phase_ispec_inner_elastic"
      call h5_write_dataset_p_2d_i(h5, dset_name, phase_ispec_inner_elastic)
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    dset_name = "nspec_inner_poroelastic"
    call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec_inner_poroelastic/))
    dset_name = "nspec_outer_poroelastic"
    call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec_outer_poroelastic/))
    dset_name = "num_phase_ispec_poroelastic"
    if (num_phase_ispec_poroelastic > 0) then
      dset_name = "phase_ispec_inner_poroelastic"
      call h5_write_dataset_p_2d_i(h5, dset_name, phase_ispec_inner_poroelastic)
    endif
  endif

  ! mesh coloring
  if (USE_MESH_COLORING_GPU) then
    if (ACOUSTIC_SIMULATION) then
      dset_name = "num_colors_outer_acoustic"
      call h5_write_dataset_p_1d_i(h5, dset_name, (/num_colors_outer_acoustic/))
      dset_name = "num_colors_inner_acoustic"
      call h5_write_dataset_p_1d_i(h5, dset_name, (/num_colors_inner_acoustic/))
      dset_name = "num_elem_colors_acoustic"
      call h5_write_dataset_p_1d_i(h5, dset_name, num_elem_colors_acoustic)
    endif
    if (ELASTIC_SIMULATION) then
      dset_name = "num_colors_outer_elastic"
      call h5_write_dataset_p_1d_i(h5, dset_name, (/num_colors_outer_elastic/))
      dset_name = "num_colors_inner_elastic"
      call h5_write_dataset_p_1d_i(h5, dset_name, (/num_colors_inner_acoustic/))
      dset_name = "num_elem_colors_elastic"
      call h5_write_dataset_p_1d_i(h5, dset_name, num_elem_colors_elastic)
    endif
  endif

  ! surface points
  dset_name = "nfaces_surface"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/nfaces_surface/))
  dset_name = "ispec_is_surface_external_mesh"
  call h5_write_dataset_p_1d_l(h5, dset_name, ispec_is_surface_external_mesh)
  dset_name = "iglob_is_surface_external_mesh"
  call h5_write_dataset_p_1d_l(h5, dset_name, iglob_is_surface_external_mesh)


  ! stores arrays in binary files
  if (SAVE_MESH_FILES) call save_arrays_solver_files(nspec,nglob,ibool)

  ! if SAVE_MESH_FILES is true then the files have already been saved, no need to save them again
  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
    call save_arrays_solver_injection_boundary(nspec,ibool)
  endif

  ! cleanup
  deallocate(ibool_interfaces_ext_mesh_dummy,stat=ier)
  if (ier /= 0) stop 'error deallocating array ibool_interfaces_ext_mesh_dummy'

  ! PML
  deallocate(is_CPML,stat=ier); if (ier /= 0) stop 'error deallocating array is_CPML'
  if (nspec_cpml_tot > 0) then
     deallocate(CPML_to_spec,stat=ier); if (ier /= 0) stop 'error deallocating array CPML_to_spec'
     deallocate(CPML_regions,stat=ier); if (ier /= 0) stop 'error deallocating array CPML_regions'
  endif

  if (PML_CONDITIONS) then
     deallocate(d_store_x,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_x'
     deallocate(d_store_y,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_y'
     deallocate(d_store_z,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_z'
     deallocate(k_store_x,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_x'
     deallocate(k_store_y,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_y'
     deallocate(k_store_z,stat=ier); if (ier /= 0) stop 'error deallocating array d_store_z'
     deallocate(alpha_store_x,stat=ier); if (ier /= 0) stop 'error deallocating array alpha_store_x'
     deallocate(alpha_store_y,stat=ier); if (ier /= 0) stop 'error deallocating array alpha_store_y'
     deallocate(alpha_store_z,stat=ier); if (ier /= 0) stop 'error deallocating array alpha_store_z'
     if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
       deallocate(mask_ibool_interior_domain,stat=ier)
       if (ier /= 0) stop 'error deallocating array mask_ibool_interior_domain'

       if (nglob_interface_PML_acoustic > 0) then
         deallocate(points_interface_PML_acoustic,stat=ier)
         if (ier /= 0) stop 'error deallocating array points_interface_PML_acoustic'
       endif

       if (nglob_interface_PML_elastic > 0) then
         deallocate(points_interface_PML_elastic,stat=ier)
         if (ier /= 0) stop 'error deallocating array points_interface_PML_elastic'
       endif
     endif
  endif
  print *, "write mesh dataset finished"
  call h5_close_group(h5)
  call h5_close_file(h5)
  call h5_destructor(h5)
  print *, "all h5 object closed"
  
  end subroutine save_arrays_solver_ext_mesh_h5


!! below all the functions of original save_arrays are kept for the further dev.

!!
!!-------------------------------------------------------------------------------------------------
!!
!
!  subroutine save_arrays_solver_files(nspec,nglob,ibool)
!
!! outputs binary files for single mesh parameters (for example vp, vs, rho, ..)
!
!  use generate_databases_par, only: myrank,NGLLX,NGLLY,NGLLZ,NGLLSQUARE,IMAIN,IOUT,FOUR_THIRDS
!
!  ! MPI interfaces
!  use generate_databases_par, only: nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,num_interfaces_ext_mesh
!
!  use create_regions_mesh_ext_par
!
!  use shared_parameters, only: NPROC
!
!  implicit none
!
!  integer,intent(in) :: nspec,nglob
!  ! mesh coordinates
!  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
!
!  ! local parameters
!  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: v_tmp
!  integer,dimension(:),allocatable :: v_tmp_i
!  integer :: ier,i,j
!  integer, dimension(:), allocatable :: iglob_tmp
!  integer :: inum, num_points
!  character(len=MAX_STRING_LEN) :: filename
!
!  !----------------------------------------------------------------------
!  ! outputs mesh files in vtk-format for visualization
!  ! (mostly for free-surface and acoustic/elastic coupling surfaces)
!  logical,parameter :: SAVE_MESH_FILES_ADDITIONAL = .true.
!
!  !----------------------------------------------------------------------
!
!  if (myrank == 0) then
!    write(IMAIN,*) '     saving mesh files for AVS, OpenDX, Paraview'
!    call flush_IMAIN()
!  endif
!
!  ! mesh arrays used for example in combine_vol_data.f90
!  !--- x coordinate
!  open(unit=IOUT,file=prname(1:len_trim(prname))//'x.bin',status='unknown',form='unformatted',iostat=ier)
!  if (ier /= 0) stop 'error opening file x.bin'
!  write(IOUT) xstore_dummy
!  close(IOUT)
!
!  !--- y coordinate
!  open(unit=IOUT,file=prname(1:len_trim(prname))//'y.bin',status='unknown',form='unformatted',iostat=ier)
!  if (ier /= 0) stop 'error opening file y.bin'
!  write(IOUT) ystore_dummy
!  close(IOUT)
!
!  !--- z coordinate
!  open(unit=IOUT,file=prname(1:len_trim(prname))//'z.bin',status='unknown',form='unformatted',iostat=ier)
!  if (ier /= 0) stop 'error opening file z.bin'
!  write(IOUT) zstore_dummy
!  close(IOUT)
!
!  ! ibool
!  open(unit=IOUT,file=prname(1:len_trim(prname))//'ibool.bin',status='unknown',form='unformatted',iostat=ier)
!  if (ier /= 0) stop 'error opening file ibool.bin'
!  write(IOUT) ibool
!  close(IOUT)
!
!  allocate(v_tmp(NGLLX,NGLLY,NGLLZ,nspec), stat=ier)
!  if (ier /= 0) call exit_MPI_without_rank('error allocating array 651')
!  if (ier /= 0) call exit_MPI_without_rank('error allocating array')
!
!  ! vp (for checking the mesh and model)
!  !minimum = minval( abs(rho_vp) )
!  !if (minimum(1) /= 0.0) then
!  !  v_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp
!  !else
!  !  v_tmp = 0.0
!  !endif
!  v_tmp = 0.0
!  where( rho_vp /= 0._CUSTOM_REAL ) v_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp
!  open(unit=IOUT,file=prname(1:len_trim(prname))//'vp.bin',status='unknown',form='unformatted',iostat=ier)
!  if (ier /= 0) stop 'error opening file vp.bin'
!  write(IOUT) v_tmp
!  close(IOUT)
!
!  ! vp values - VTK file output
!  filename = prname(1:len_trim(prname))//'vp'
!  call write_VTK_data_gll_cr(nspec,nglob, &
!                      xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!                      v_tmp,filename)
!
!
!  ! vs (for checking the mesh and model)
!  !minimum = minval( abs(rho_vs) )
!  !if (minimum(1) /= 0.0) then
!  !  v_tmp = mustore / rho_vs
!  !else
!  !  v_tmp = 0.0
!  !endif
!  v_tmp = 0.0
!  where( rho_vs /= 0._CUSTOM_REAL )  v_tmp = mustore / rho_vs
!  open(unit=IOUT,file=prname(1:len_trim(prname))//'vs.bin',status='unknown',form='unformatted',iostat=ier)
!  if (ier /= 0) stop 'error opening file vs.bin'
!  write(IOUT) v_tmp
!  close(IOUT)
!
!  ! vs values - VTK file output
!  filename = prname(1:len_trim(prname))//'vs'
!  call write_VTK_data_gll_cr(nspec,nglob, &
!                      xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!                      v_tmp,filename)
!
!  ! outputs density model for check
!  v_tmp = 0.0
!  where( rho_vp /= 0._CUSTOM_REAL ) v_tmp = rho_vp**2 / (FOUR_THIRDS * mustore + kappastore)
!  open(unit=IOUT,file=prname(1:len_trim(prname))//'rho.bin',status='unknown',form='unformatted',iostat=ier)
!  if (ier /= 0) stop 'error opening file rho.bin'
!  write(IOUT) v_tmp
!  close(IOUT)
!
!  ! attenuation
!  ! shear attenuation Qmu
!  open(unit=IOUT,file=prname(1:len_trim(prname))//'qmu.bin',status='unknown',form='unformatted',iostat=ier)
!  if (ier /= 0) stop 'error opening file qmu.bin'
!  write(IOUT) qmu_attenuation_store
!  close(IOUT)
!
!  ! shear attenuation - VTK file output
!  filename = prname(1:len_trim(prname))//'qmu'
!  call write_VTK_data_gll_cr(nspec,nglob, &
!                      xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!                      qmu_attenuation_store,filename)
!
!  ! bulk attenuation Qkappa
!  open(unit=IOUT,file=prname(1:len_trim(prname))//'qkappa.bin',status='unknown',form='unformatted',iostat=ier)
!  if (ier /= 0) stop 'error opening file qkappa.bin'
!  write(IOUT) qkappa_attenuation_store
!  close(IOUT)
!
!  ! bulk attenuation - VTK file output
!  filename = prname(1:len_trim(prname))//'qkappa'
!  call write_VTK_data_gll_cr(nspec,nglob, &
!                      xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!                      qkappa_attenuation_store,filename)
!
!  ! frees temporary array
!  deallocate(v_tmp)
!
!  ! additional VTK file output
!  if (SAVE_MESH_FILES_ADDITIONAL) then
!    ! user output
!    call synchronize_all()
!    if (myrank == 0) then
!      write(IMAIN,*) '     saving additonal mesh files with surface/coupling points'
!      call flush_IMAIN()
!    endif
!
!    ! saves free surface points
!    if (num_free_surface_faces > 0) then
!      ! saves free surface interface points
!      allocate( iglob_tmp(NGLLSQUARE*num_free_surface_faces),stat=ier)
!      if (ier /= 0) call exit_MPI_without_rank('error allocating array 652')
!      if (ier /= 0) stop 'error allocating array iglob_tmp'
!      inum = 0
!      iglob_tmp(:) = 0
!      do i=1,num_free_surface_faces
!        do j=1,NGLLSQUARE
!          inum = inum+1
!          iglob_tmp(inum) = ibool(free_surface_ijk(1,j,i), &
!                                  free_surface_ijk(2,j,i), &
!                                  free_surface_ijk(3,j,i), &
!                                  free_surface_ispec(i) )
!        enddo
!      enddo
!      filename = prname(1:len_trim(prname))//'free_surface'
!      call write_VTK_data_points(nglob, &
!                        xstore_dummy,ystore_dummy,zstore_dummy, &
!                        iglob_tmp,NGLLSQUARE*num_free_surface_faces, &
!                        filename)
!
!      deallocate(iglob_tmp)
!    endif
!
!    ! acoustic-elastic domains
!    if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION) then
!      ! saves points on acoustic-elastic coupling interface
!      num_points = NGLLSQUARE*num_coupling_ac_el_faces
!      allocate( iglob_tmp(num_points),stat=ier)
!      if (ier /= 0) call exit_MPI_without_rank('error allocating array 653')
!      if (ier /= 0) stop 'error allocating array iglob_tmp'
!      inum = 0
!      iglob_tmp(:) = 0
!      do i = 1,num_coupling_ac_el_faces
!        do j = 1,NGLLSQUARE
!          inum = inum+1
!          iglob_tmp(inum) = ibool(coupling_ac_el_ijk(1,j,i), &
!                                  coupling_ac_el_ijk(2,j,i), &
!                                  coupling_ac_el_ijk(3,j,i), &
!                                  coupling_ac_el_ispec(i) )
!        enddo
!      enddo
!      filename = prname(1:len_trim(prname))//'coupling_acoustic_elastic'
!      call write_VTK_data_points(nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
!                                 iglob_tmp,num_points,filename)
!
!      ! saves acoustic/elastic flag
!      allocate(v_tmp_i(nspec),stat=ier)
!      if (ier /= 0) call exit_MPI_without_rank('error allocating array 654')
!      if (ier /= 0) stop 'error allocating array v_tmp_i'
!      do i = 1,nspec
!        if (ispec_is_acoustic(i)) then
!          v_tmp_i(i) = 1
!        else if (ispec_is_elastic(i)) then
!          v_tmp_i(i) = 2
!        else
!          v_tmp_i(i) = 0
!        endif
!      enddo
!      filename = prname(1:len_trim(prname))//'acoustic_elastic_flag'
!      call write_VTK_data_elem_i(nspec,nglob,xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!                                 v_tmp_i,filename)
!
!      deallocate(iglob_tmp,v_tmp_i)
!    endif !if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION )
!
!    ! acoustic-poroelastic domains
!    if (ACOUSTIC_SIMULATION .and. POROELASTIC_SIMULATION) then
!      ! saves points on acoustic-poroelastic coupling interface
!      num_points = NGLLSQUARE*num_coupling_ac_po_faces
!      allocate( iglob_tmp(num_points),stat=ier)
!      if (ier /= 0) call exit_MPI_without_rank('error allocating array 655')
!      if (ier /= 0) stop 'error allocating array iglob_tmp'
!      inum = 0
!      iglob_tmp(:) = 0
!      do i = 1,num_coupling_ac_po_faces
!        do j = 1,NGLLSQUARE
!          inum = inum+1
!          iglob_tmp(inum) = ibool(coupling_ac_po_ijk(1,j,i), &
!                                  coupling_ac_po_ijk(2,j,i), &
!                                  coupling_ac_po_ijk(3,j,i), &
!                                  coupling_ac_po_ispec(i) )
!        enddo
!      enddo
!      filename = prname(1:len_trim(prname))//'coupling_acoustic_poroelastic'
!      call write_VTK_data_points(nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
!                                 iglob_tmp,num_points,filename)
!
!      ! saves acoustic/poroelastic flag
!      allocate(v_tmp_i(nspec),stat=ier)
!      if (ier /= 0) call exit_MPI_without_rank('error allocating array 656')
!      if (ier /= 0) stop 'error allocating array v_tmp_i'
!      do i = 1,nspec
!        if (ispec_is_acoustic(i)) then
!          v_tmp_i(i) = 1
!        else if (ispec_is_poroelastic(i)) then
!          v_tmp_i(i) = 2
!        else
!          v_tmp_i(i) = 0
!        endif
!      enddo
!      filename = prname(1:len_trim(prname))//'acoustic_poroelastic_flag'
!      call write_VTK_data_elem_i(nspec,nglob,xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!                                 v_tmp_i,filename)
!
!      deallocate(v_tmp_i,iglob_tmp)
!    endif !if (ACOUSTIC_SIMULATION .and. POROELASTIC_SIMULATION )
!
!    ! elastic-poroelastic domains
!    if (ELASTIC_SIMULATION .and. POROELASTIC_SIMULATION) then
!      ! saves points on elastic-poroelastic coupling interface
!      num_points = NGLLSQUARE*num_coupling_el_po_faces
!      allocate( iglob_tmp(num_points),stat=ier)
!      if (ier /= 0) call exit_MPI_without_rank('error allocating array 657')
!      if (ier /= 0) stop 'error allocating array iglob_tmp'
!      inum = 0
!      iglob_tmp(:) = 0
!      do i = 1,num_coupling_el_po_faces
!        do j = 1,NGLLSQUARE
!          inum = inum+1
!          iglob_tmp(inum) = ibool(coupling_el_po_ijk(1,j,i), &
!                                  coupling_el_po_ijk(2,j,i), &
!                                  coupling_el_po_ijk(3,j,i), &
!                                  coupling_el_po_ispec(i) )
!        enddo
!      enddo
!      filename = prname(1:len_trim(prname))//'coupling_elastic_poroelastic'
!      call write_VTK_data_points(nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
!                                 iglob_tmp,num_points,filename)
!
!      ! saves elastic/poroelastic flag
!      allocate(v_tmp_i(nspec),stat=ier)
!      if (ier /= 0) call exit_MPI_without_rank('error allocating array 658')
!      if (ier /= 0) stop 'error allocating array v_tmp_i'
!      do i=1,nspec
!        if (ispec_is_elastic(i)) then
!          v_tmp_i(i) = 1
!        else if (ispec_is_poroelastic(i)) then
!          v_tmp_i(i) = 2
!        else
!          v_tmp_i(i) = 0
!        endif
!      enddo
!      filename = prname(1:len_trim(prname))//'elastic_poroelastic_flag'
!      call write_VTK_data_elem_i(nspec,nglob,xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!                                 v_tmp_i,filename)
!
!      deallocate(v_tmp_i,iglob_tmp)
!    endif !if (ACOUSTIC_SIMULATION .and. POROELASTIC_SIMULATION
!
!    ! MPI
!    if (NPROC > 1) then
!      ! saves MPI interface points
!      num_points = sum(nibool_interfaces_ext_mesh(1:num_interfaces_ext_mesh))
!      allocate( iglob_tmp(num_points),stat=ier)
!      if (ier /= 0) call exit_MPI_without_rank('error allocating array 659')
!      if (ier /= 0) stop 'error allocating array iglob_tmp'
!      inum = 0
!      iglob_tmp(:) = 0
!      do i = 1,num_interfaces_ext_mesh
!        do j = 1, nibool_interfaces_ext_mesh(i)
!          inum = inum + 1
!          iglob_tmp(inum) = ibool_interfaces_ext_mesh(j,i)
!        enddo
!      enddo
!
!      filename = prname(1:len_trim(prname))//'MPI_points'
!      call write_VTK_data_points(nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
!                                 iglob_tmp,num_points,filename)
!      deallocate(iglob_tmp)
!    endif ! NPROC > 1
!  endif  !if (SAVE_MESH_FILES_ADDITIONAL)
!
!  end subroutine save_arrays_solver_files
!
!!
!!-------------------------------------------------------------------------------------------------
!!
!
!  subroutine save_arrays_solver_injection_boundary(nspec,ibool)
!
!  use generate_databases_par, only: myrank,NGLLX,NGLLY,NGLLZ,NGLLSQUARE,IMAIN,IOUT
!
!  use create_regions_mesh_ext_par
!
!  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH
!
!  implicit none
!
!  integer,intent(in) :: nspec
!  ! mesh coordinates
!  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
!
!  ! local parameters
!  integer :: ier,i,j,k
!  integer :: iface, ispec, iglob, igll
!  real(kind=CUSTOM_REAL) :: nx,ny,nz
!  character(len=MAX_STRING_LEN) :: filename
!
!  if (myrank == 0) then
!    write(IMAIN,*) '     saving mesh files for coupled injection boundary'
!    call flush_IMAIN()
!  endif
!
!  ! checks if anything to do
!  if (.not. (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH)) return
!
!  filename = prname(1:len_trim(prname))//'absorb_dsm'
!  open(IOUT,file=filename(1:len_trim(filename)),status='unknown',form='unformatted',iostat=ier)
!  if (ier /= 0) stop 'error opening file absorb_dsm'
!  write(IOUT) num_abs_boundary_faces
!  write(IOUT) abs_boundary_ispec
!  write(IOUT) abs_boundary_ijk
!  write(IOUT) abs_boundary_jacobian2Dw
!  write(IOUT) abs_boundary_normal
!  close(IOUT)
!
!  filename = prname(1:len_trim(prname))//'inner'
!  open(IOUT,file=filename(1:len_trim(filename)),status='unknown',form='unformatted',iostat=ier)
!  write(IOUT) ispec_is_inner
!  write(IOUT) ispec_is_elastic
!  close(IOUT)
!
!  !! VM VM write an ascii file for instaseis input
!  filename = prname(1:len_trim(prname))//'normal.txt'
!  open(IOUT,file=filename(1:len_trim(filename)),status='unknown',iostat=ier)
!  write(IOUT, *) ' number of points :', num_abs_boundary_faces*NGLLSQUARE
!
!  do iface = 1,num_abs_boundary_faces
!     ispec = abs_boundary_ispec(iface)
!     if (ispec_is_elastic(ispec)) then
!        do igll = 1,NGLLSQUARE
!
!           ! gets local indices for GLL point
!           i = abs_boundary_ijk(1,igll,iface)
!           j = abs_boundary_ijk(2,igll,iface)
!           k = abs_boundary_ijk(3,igll,iface)
!
!           iglob = ibool(i,j,k,ispec)
!
!           nx = abs_boundary_normal(1,igll,iface)
!           ny = abs_boundary_normal(2,igll,iface)
!           nz = abs_boundary_normal(3,igll,iface)
!
!           write(IOUT,'(6f25.10)') xstore_dummy(iglob), ystore_dummy(iglob), zstore_dummy(iglob), nx, ny, nz
!
!        enddo
!     endif
!  enddo
!  close(IOUT)
!
!  end subroutine save_arrays_solver_injection_boundary
