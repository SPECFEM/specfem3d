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

! HDF5 file i/o

  subroutine read_mesh_for_init_hdf5()

! reads in the value of NSPEC_AB and NGLOB_AB

#ifdef USE_HDF5
  use specfem_par
  use manager_hdf5

  implicit none
  ! Local variables
  character(len=MAX_STRING_LEN) :: database_hdf5
  character(len=30) :: dsetname
  character(len=64) :: tempstr

  ! MPI variables
  integer :: info, comm

  ! if collective read
  logical :: if_col = .true.

  if (I_should_read_the_database) then
    ! set file name
    tempstr = "/external_mesh.h5"
    database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//trim(tempstr)

    ! initialze hdf5
    call world_get_comm(comm)
    call world_get_info_null(info)

    call h5_initialize()
    call h5_set_mpi_info(comm, info, myrank, NPROC)

    ! opens file
    call h5_open_file_p_collect(database_hdf5)

    ! read datasets
    dsetname = "nspec"
    !call h5_read_dataset_p_scalar(dsetname, NSPEC_AB)
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, NSPEC_AB, (/myrank/),if_col)
    dsetname = "nglob"
    !call h5_read_dataset_p_scalar(dsetname, NGLOB_AB)
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, NGLOB_AB, (/myrank/),if_col)
    dsetname = "nspec_irregular"
    !call h5_read_dataset_p_scalar(dsetname, NSPEC_IRREGULAR)
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, NSPEC_IRREGULAR, (/myrank/),if_col)

    ! close hdf5
    call h5_close_file()
    call h5_finalize()
  endif

  call bcast_all_i_for_database(NSPEC_AB, 1)
  call bcast_all_i_for_database(NGLOB_AB, 1)
  call bcast_all_i_for_database(NSPEC_IRREGULAR, 1)

#else
  ! no HDF5 compilation support
  ! user output
  print *
  print *, "Error: HDF5 routine read_mesh_for_init_hdf5() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *
  stop 'Error HDF5 read_mesh_for_init_hdf5(): called without compilation support'

#endif

  end subroutine read_mesh_for_init_hdf5


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_hdf5()

! reads in full database

#ifdef USE_HDF5
  use pml_par

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use specfem_par_movie

  use manager_hdf5

  implicit none

  ! Local variables
  integer :: itest,ier
  character(len=MAX_STRING_LEN) :: database_hdf5
  character(len=40) :: dsetname
  character(len=64) :: tempstr

  ! MPI variables
  integer :: info, comm

  ! if collective read
  logical :: if_col = .true.

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
  !integer, dimension(0:NPROC-1) :: offset_max_ni_bool_interfaces_ext_mesh
  integer, dimension(0:NPROC-1) :: offset_nspec_aniso
  integer, dimension(0:NPROC-1) :: offset_num_phase_ispec_acoustic
  integer, dimension(0:NPROC-1) :: offset_num_phase_ispec_elastic
  integer, dimension(0:NPROC-1) :: offset_num_phase_ispec_poroelastic
  integer, dimension(0:NPROC-1) :: offset_num_colors_acoustic
  integer, dimension(0:NPROC-1) :: offset_num_colors_elastic
  integer, dimension(0:NPROC-1) :: offset_nspec_ab
  integer, dimension(0:NPROC-1) :: offset_nglob_ab

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "Reading mesh databases..."
    write(IMAIN,*) "  database file          : external_mesh.h5"
    write(IMAIN,*) "  from directory         : ",trim(LOCAL_PATH)
    write(IMAIN,*) "  using HDF5 file format"
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! start reading the databases
  if (I_should_read_the_database) then
    ! set file name
    tempstr = "/external_mesh.h5"
    database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//trim(tempstr)

    ! initialze hdf5
    call world_get_comm(comm)
    call world_get_info_null(info)

    call h5_initialize()
    call h5_set_mpi_info(comm, info, myrank, NPROC)

    ! open file
    call h5_open_file_p_collect(database_hdf5)

    ! gets offsets
    call h5_read_dataset_collect_hyperslab("offset_nglob",offset_nglob, (/0/), if_col)
    call h5_read_dataset_collect_hyperslab("offset_nspec",offset_nspec, (/0/), if_col)
    call h5_read_dataset_collect_hyperslab("offset_nspec_irregular",offset_nspec_irregular, (/0/), if_col)

    call h5_read_dataset_collect_hyperslab("offset_nspec2D_xmin",offset_nspec2D_xmin, (/0/), if_col)
    call h5_read_dataset_collect_hyperslab("offset_nspec2D_xmax",offset_nspec2D_xmax, (/0/), if_col)
    call h5_read_dataset_collect_hyperslab("offset_nspec2D_ymin",offset_nspec2D_ymin, (/0/), if_col)
    call h5_read_dataset_collect_hyperslab("offset_nspec2D_ymax",offset_nspec2D_ymax, (/0/), if_col)
    call h5_read_dataset_collect_hyperslab("offset_nspec2D_bottom_ext",offset_nspec2D_bottom_ext, (/0/), if_col)
    call h5_read_dataset_collect_hyperslab("offset_nspec2D_top_ext",offset_nspec2D_top_ext, (/0/), if_col)

    call h5_read_dataset_collect_hyperslab("offset_nspec_ab",offset_nspec_ab, (/0/), if_col)
    call h5_read_dataset_collect_hyperslab("offset_nglob_ab",offset_nglob_ab, (/0/), if_col)

    ! info about external mesh simulation
    dsetname = "nspec"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, NSPEC_AB, (/myrank/),if_col)
    dsetname = "nglob"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, NGLOB_AB, (/myrank/),if_col)
    dsetname = "nspec_irregular"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, NSPEC_IRREGULAR, (/myrank/),if_col)

    dsetname = "ibool"
    call h5_read_dataset_collect_hyperslab(dsetname, ibool, (/0,0,0,sum(offset_nglob(0:myrank-1))/), if_col)
    dsetname = "xstore_unique"
    call h5_read_dataset_collect_hyperslab(dsetname,xstore,(/sum(offset_nglob(0:myrank-1))/),if_col)
    dsetname = "ystore_unique"
    call h5_read_dataset_collect_hyperslab(dsetname,ystore,(/sum(offset_nglob(0:myrank-1))/),if_col)
    dsetname = "zstore_unique"
    call h5_read_dataset_collect_hyperslab(dsetname,zstore,(/sum(offset_nglob(0:myrank-1))/),if_col)

    dsetname = "irregular_element_number"
    call h5_read_dataset_collect_hyperslab(dsetname,irregular_element_number,(/sum(offset_nspec(0:myrank-1))/),if_col)
    dsetname = "xix_regular"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname,xix_regular,(/myrank/),if_col)
    dsetname = "jacobian_regular"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname,jacobian_regular,(/myrank/),if_col)

    dsetname = "xixstore"
    call h5_read_dataset_collect_hyperslab(dsetname,xixstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
    dsetname = "xiystore"
    call h5_read_dataset_collect_hyperslab(dsetname,xiystore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
    dsetname = "xizstore"
    call h5_read_dataset_collect_hyperslab(dsetname,xizstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
    dsetname = "etaxstore"
    call h5_read_dataset_collect_hyperslab(dsetname,etaxstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
    dsetname = "etaystore"
    call h5_read_dataset_collect_hyperslab(dsetname,etaystore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
    dsetname = "etazstore"
    call h5_read_dataset_collect_hyperslab(dsetname,etazstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
    dsetname = "gammaxstore"
    call h5_read_dataset_collect_hyperslab(dsetname,gammaxstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
    dsetname = "gammaystore"
    call h5_read_dataset_collect_hyperslab(dsetname,gammaystore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
    dsetname = "gammazstore"
    call h5_read_dataset_collect_hyperslab(dsetname,gammazstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)
    dsetname = "jacobianstore"
    call h5_read_dataset_collect_hyperslab(dsetname,jacobianstore,(/0,0,0,sum(offset_nspec_irregular(0:myrank-1))/),if_col)

    dsetname = "ispec_is_acoustic"
    call h5_read_dataset_collect_hyperslab(dsetname, ispec_is_acoustic, (/sum(offset_nspec(0:myrank-1))/),if_col)
    dsetname = "ispec_is_elastic"
    call h5_read_dataset_collect_hyperslab(dsetname, ispec_is_elastic, (/sum(offset_nspec(0:myrank-1))/),if_col)
    dsetname = "ispec_is_poroelastic"
    call h5_read_dataset_collect_hyperslab(dsetname, ispec_is_poroelastic, (/sum(offset_nspec(0:myrank-1))/),if_col)

    dsetname = "kappastore"
    call h5_read_dataset_collect_hyperslab(dsetname,kappastore,(/0,0,0,sum(offset_nspec(0:myrank-1))/),if_col)
    dsetname = "mustore"
    call h5_read_dataset_collect_hyperslab(dsetname,mustore,(/0,0,0,sum(offset_nspec(0:myrank-1))/),if_col)
  endif

  call bcast_all_i_for_database(NSPEC_AB, 1)
  call bcast_all_i_for_database(NGLOB_AB, 1)
  call bcast_all_i_for_database(NSPEC_IRREGULAR, 1)
  call bcast_all_i_for_database(ibool(1,1,1,1), size(ibool))
  call bcast_all_i_for_database(irregular_element_number(1), size(irregular_element_number))
  call bcast_all_cr_for_database(xix_regular, 1)
  call bcast_all_cr_for_database(jacobian_regular, 1)

  if (size(xstore) > 0) call bcast_all_cr_for_database(xstore(1), size(xstore))
  if (size(ystore) > 0) call bcast_all_cr_for_database(ystore(1), size(ystore))
  if (size(zstore) > 0) call bcast_all_cr_for_database(zstore(1), size(zstore))

  call bcast_all_cr_for_database(xixstore(1,1,1,1), size(xixstore))
  call bcast_all_cr_for_database(xiystore(1,1,1,1), size(xiystore))
  call bcast_all_cr_for_database(xizstore(1,1,1,1), size(xizstore))
  call bcast_all_cr_for_database(etaxstore(1,1,1,1), size(etaxstore))
  call bcast_all_cr_for_database(etaystore(1,1,1,1), size(etaystore))
  call bcast_all_cr_for_database(etazstore(1,1,1,1), size(etazstore))
  call bcast_all_cr_for_database(gammaxstore(1,1,1,1), size(gammaxstore))
  call bcast_all_cr_for_database(gammaystore(1,1,1,1), size(gammaystore))
  call bcast_all_cr_for_database(gammazstore(1,1,1,1), size(gammazstore))
  call bcast_all_cr_for_database(jacobianstore(1,1,1,1), size(jacobianstore))

  call bcast_all_cr_for_database(kappastore(1,1,1,1), size(kappastore))
  call bcast_all_cr_for_database(mustore(1,1,1,1), size(mustore))

  if (size(ispec_is_acoustic) > 0) &
    call bcast_all_l_for_database(ispec_is_acoustic(1), size(ispec_is_acoustic))
  if (size(ispec_is_elastic) > 0) &
    call bcast_all_l_for_database(ispec_is_elastic(1), size(ispec_is_elastic))
  if (size(ispec_is_poroelastic) > 0) &
    call bcast_all_l_for_database(ispec_is_poroelastic(1), size(ispec_is_poroelastic))

  ! all processes will have acoustic_simulation set if any flag is .true.
  call any_all_l( ANY(ispec_is_acoustic), ACOUSTIC_SIMULATION )
  call any_all_l( ANY(ispec_is_elastic), ELASTIC_SIMULATION )
  call any_all_l( ANY(ispec_is_poroelastic), POROELASTIC_SIMULATION )

  ! number of acoustic elements in this partition
  nspec_acoustic = count(ispec_is_acoustic(:))
  ! number of elastic elements in this partition
  nspec_elastic = count(ispec_is_elastic(:))
  ! number of elastic elements in this partition
  nspec_poroelastic = count(ispec_is_poroelastic(:))

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  simulation w/ acoustic    domain: ',ACOUSTIC_SIMULATION
    write(IMAIN,*) '  simulation w/ elastic     domain: ',ELASTIC_SIMULATION
    write(IMAIN,*) '  simulation w/ poroelastic domain: ',POROELASTIC_SIMULATION
    write(IMAIN,*)
    write(IMAIN,*) '  slice 0 has:'
    write(IMAIN,*) '  number of elements acoustic   :',nspec_acoustic
    write(IMAIN,*) '  number of elements elastic    :',nspec_elastic
    write(IMAIN,*) '  number of elements poroelastic:',nspec_poroelastic
    call flush_IMAIN()
  endif

  ! acoustic
  if (ACOUSTIC_SIMULATION) then
    ! potentials
    ! NB_RUNS_ACOUSTIC_GPU is set to 1 by default in constants.h
    allocate(potential_acoustic(NGLOB_AB*NB_RUNS_ACOUSTIC_GPU),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1420')
    if (ier /= 0) stop 'Error allocating array potential_acoustic'
    allocate(potential_dot_acoustic(NGLOB_AB*NB_RUNS_ACOUSTIC_GPU),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1421')
    if (ier /= 0) stop 'Error allocating array potential_dot_acoustic'
    allocate(potential_dot_dot_acoustic(NGLOB_AB*NB_RUNS_ACOUSTIC_GPU),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1422')
    if (ier /= 0) stop 'Error allocating array potential_dot_dot_acoustic'
    potential_acoustic(:) = 0.0_CUSTOM_REAL; potential_dot_acoustic(:) = 0.0_CUSTOM_REAL
    potential_dot_dot_acoustic(:) = 0.0_CUSTOM_REAL

    !if (SIMULATION_TYPE /= 1) then
    !  allocate(potential_acoustic_adj_coupling(NGLOB_AB),stat=ier)
    !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1423')
    !  if (ier /= 0) stop 'Error allocating array potential_acoustic_adj_coupling'
    !endif

    ! mass matrix, density
    allocate(rmass_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1424')
    if (ier /= 0) stop 'Error allocating array rmass_acoustic'
    rmass_acoustic(:) = 0.0_CUSTOM_REAL

    if (I_should_read_the_database) then
      dsetname = "rmass_acoustic"
      call h5_read_dataset_collect_hyperslab(dsetname, rmass_acoustic,(/sum(offset_nglob(0:myrank-1))/),if_col)
    endif
    if (size(rmass_acoustic) > 0) call bcast_all_cr_for_database(rmass_acoustic(1), size(rmass_acoustic))

    ! initializes mass matrix contribution
    allocate(rmassz_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1425')
    if (ier /= 0) stop 'Error allocating array rmassz_acoustic'
    rmassz_acoustic(:) = 0.0_CUSTOM_REAL
  endif

! rho array is needed for acoustic simulations but also for elastic simulations with CPML,
! read it in all cases (whether the simulation is acoustic, elastic, or acoustic/elastic)
  if (I_should_read_the_database) then
    dsetname = "rhostore"
    !call h5_read_dataset_p(dsetname, rhostore)
    call h5_read_dataset_collect_hyperslab(dsetname, rhostore, (/0,0,0,sum(offset_nspec(0:myrank-1))/),if_col)
  endif
  call bcast_all_cr_for_database(rhostore(1,1,1,1), size(rhostore))

  ! elastic simulation
  if (ELASTIC_SIMULATION) then
    ! checks
    if (NB_RUNS_ACOUSTIC_GPU > 1) &
      stop 'NB_RUNS_ACOUSTIC_GPU > 1 not compatible with elastic or coupled simulations'

    ! displacement,velocity,acceleration
    allocate(displ(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1427')
    if (ier /= 0) stop 'Error allocating array displ'
    allocate(veloc(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1428')
    if (ier /= 0) stop 'Error allocating array veloc'
    allocate(accel(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1429')
    if (ier /= 0) stop 'Error allocating array accel'
    displ(:,:) = 0.0_CUSTOM_REAL; veloc(:,:) = 0.0_CUSTOM_REAL; accel(:,:) = 0.0_CUSTOM_REAL

    if (SIMULATION_TYPE /= 1) then
      allocate(accel_adj_coupling(NDIM,NGLOB_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1430')
      if (ier /= 0) stop 'Error allocating array accel_adj_coupling'
      accel_adj_coupling(:,:) = 0.0_CUSTOM_REAL
    endif

    ! allocates mass matrix
    allocate(rmass(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1431')
    rmass(:) = 0.0_CUSTOM_REAL

    if (ier /= 0) stop 'Error allocating array rmass'
    ! initializes mass matrix contributions
    allocate(rmassx(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1432')
    allocate(rmassy(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1433')
    allocate(rmassz(NGLOB_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1434')
    if (ier /= 0) stop 'Error allocating array rmassx,rmassy,rmassz'
    rmassx(:) = 0._CUSTOM_REAL
    rmassy(:) = 0._CUSTOM_REAL
    rmassz(:) = 0._CUSTOM_REAL

    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1435')
    if (ier /= 0) stop 'Error allocating array rho_vp'
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1436')
    if (ier /= 0) stop 'Error allocating array rho_vs'
    rho_vp(:,:,:,:) = 0.0_CUSTOM_REAL; rho_vs(:,:,:,:) = 0.0_CUSTOM_REAL

    allocate(c11store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1437')
    allocate(c12store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1438')
    allocate(c13store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1439')
    allocate(c14store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1440')
    allocate(c15store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1441')
    allocate(c16store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1442')
    allocate(c22store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1443')
    allocate(c23store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1444')
    allocate(c24store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1445')
    allocate(c25store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1446')
    allocate(c26store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1447')
    allocate(c33store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1448')
    allocate(c34store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1449')
    allocate(c35store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1450')
    allocate(c36store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1451')
    allocate(c44store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1452')
    allocate(c45store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1453')
    allocate(c46store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1454')
    allocate(c55store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1455')
    allocate(c56store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1456')
    allocate(c66store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1457')
    if (ier /= 0) stop 'Error allocating array c11store etc.'
    c11store(:,:,:,:) = 0.0_CUSTOM_REAL; c12store(:,:,:,:) = 0.0_CUSTOM_REAL; c13store(:,:,:,:) = 0.0_CUSTOM_REAL
    c14store(:,:,:,:) = 0.0_CUSTOM_REAL; c15store(:,:,:,:) = 0.0_CUSTOM_REAL; c16store(:,:,:,:) = 0.0_CUSTOM_REAL
    c22store(:,:,:,:) = 0.0_CUSTOM_REAL; c23store(:,:,:,:) = 0.0_CUSTOM_REAL; c24store(:,:,:,:) = 0.0_CUSTOM_REAL
    c25store(:,:,:,:) = 0.0_CUSTOM_REAL; c26store(:,:,:,:) = 0.0_CUSTOM_REAL; c33store(:,:,:,:) = 0.0_CUSTOM_REAL
    c34store(:,:,:,:) = 0.0_CUSTOM_REAL; c35store(:,:,:,:) = 0.0_CUSTOM_REAL; c36store(:,:,:,:) = 0.0_CUSTOM_REAL
    c44store(:,:,:,:) = 0.0_CUSTOM_REAL; c45store(:,:,:,:) = 0.0_CUSTOM_REAL; c46store(:,:,:,:) = 0.0_CUSTOM_REAL
    c55store(:,:,:,:) = 0.0_CUSTOM_REAL; c56store(:,:,:,:) = 0.0_CUSTOM_REAL; c66store(:,:,:,:) = 0.0_CUSTOM_REAL

    ! note: currently, they need to be defined, as they are used in some subroutine arguments
    allocate(R_xx(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1458')
    allocate(R_yy(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1459')
    allocate(R_xy(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1460')
    allocate(R_xz(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1461')
    allocate(R_yz(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1462')
    if (ier /= 0) stop 'Error allocating array R_xx etc.'
    R_xx(:,:,:,:,:) = 0.0_CUSTOM_REAL; R_yy(:,:,:,:,:) = 0.0_CUSTOM_REAL; R_xy(:,:,:,:,:) = 0.0_CUSTOM_REAL
    R_xz(:,:,:,:,:) = 0.0_CUSTOM_REAL; R_yz(:,:,:,:,:) = 0.0_CUSTOM_REAL

    ! needed for attenuation and/or kernel computations
    allocate(epsilondev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1463')
    allocate(epsilondev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1464')
    allocate(epsilondev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1465')
    allocate(epsilondev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1466')
    allocate(epsilondev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1467')
    allocate(epsilondev_trace(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1468')
    if (ier /= 0) stop 'Error allocating array epsilondev_xx etc.'
    epsilondev_xx(:,:,:,:) = 0.0_CUSTOM_REAL; epsilondev_yy(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilondev_xy(:,:,:,:) = 0.0_CUSTOM_REAL; epsilondev_xz(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilondev_yz(:,:,:,:) = 0.0_CUSTOM_REAL; epsilondev_trace(:,:,:,:) = 0.0_CUSTOM_REAL

    allocate(R_trace(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1469')
    if (ier /= 0) stop 'Error allocating array R_trace etc.'
    R_trace(:,:,:,:,:) = 0.0_CUSTOM_REAL

    ! note: needed for some subroutine arguments
    allocate(epsilon_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1470')
    if (ier /= 0) stop 'Error allocating array epsilon_trace_over_3'
    epsilon_trace_over_3(:,:,:,:) = 0.0_CUSTOM_REAL

    ! needed for attenuation
    allocate(factor_common(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1471')
    if (ier /= 0) stop 'Error allocating array factor_common etc.'
    factor_common(:,:,:,:,:) = 0.0_CUSTOM_REAL

    allocate(factor_common_kappa(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1472')
    if (ier /= 0) stop 'Error allocating array factor_common_kappa etc.'
    factor_common_kappa(:,:,:,:,:) = 0.0_CUSTOM_REAL

    ! reads mass matrices
    if (I_should_read_the_database) then
      dsetname = "rmass"
      call h5_read_dataset_collect_hyperslab(dsetname, rmass, (/sum(offset_nglob(0:myrank-1))/), if_col)
    endif
    call bcast_all_cr_for_database(rmass(1), size(rmass))
    if (ier /= 0) stop 'Error reading in array rmass'

    if (APPROXIMATE_OCEAN_LOAD) then
      ! ocean mass matrix
      allocate(rmass_ocean_load(NGLOB_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1473')
      if (ier /= 0) stop 'Error allocating array rmass_ocean_load'
      rmass_ocean_load(:) = 0.0_CUSTOM_REAL

      if (I_should_read_the_database) then
        call h5_read_dataset_collect_hyperslab("offset_nglob_ocean",offset_nglob_ocean, (/0/), if_col)
        dsetname = "rmass_ocean_load"
        call h5_read_dataset_collect_hyperslab(dsetname, rmass_ocean_load, (/sum(offset_nglob_ocean(0:myrank-1))/),if_col)
      endif
      if (size(rmass_ocean_load) > 0) call bcast_all_cr_for_database(rmass_ocean_load(1), size(rmass_ocean_load))
    else
      ! dummy allocation
      allocate(rmass_ocean_load(1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1474')
      if (ier /= 0) stop 'Error allocating dummy array rmass_ocean_load'
    endif

    !pll material parameters for stacey conditions
    if (I_should_read_the_database) then
      dsetname = "rho_vp"
      call h5_read_dataset_collect_hyperslab(dsetname, rho_vp, (/0,0,0,sum(offset_nspec(0:myrank-1))/),if_col)
    endif
    if (size(rho_vp) > 0) call bcast_all_cr_for_database(rho_vp(1,1,1,1), size(rho_vp))
    if (I_should_read_the_database) then
      dsetname = "rho_vs"
      call h5_read_dataset_collect_hyperslab(dsetname, rho_vs, (/0,0,0,sum(offset_nspec(0:myrank-1))/),if_col)
    endif
    if (size(rho_vs) > 0) call bcast_all_cr_for_database(rho_vs(1,1,1,1), size(rho_vs))

  else
    ! no elastic attenuation & anisotropy
    ATTENUATION = .false.
    ANISOTROPY = .false.
  endif

  ! poroelastic
  if (POROELASTIC_SIMULATION) then
    ! checks
    if (GPU_MODE) &
      call exit_mpi(myrank,'POROELASTICITY not supported by GPU mode yet...')
    if (NB_RUNS_ACOUSTIC_GPU > 1) &
      call exit_mpi(myrank,'NB_RUNS_ACOUSTIC_GPU > 1 not compatible with poroelastic or coupled simulations')

    ! displacement,velocity,acceleration for the solid (s) & fluid (w) phases
    allocate(displs_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1475')
    if (ier /= 0) stop 'Error allocating array displs_poroelastic'
    allocate(velocs_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1476')
    if (ier /= 0) stop 'Error allocating array velocs_poroelastic'
    allocate(accels_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1477')
    if (ier /= 0) stop 'Error allocating array accels_poroelastic'
    allocate(displw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1478')
    if (ier /= 0) stop 'Error allocating array displw_poroelastic'
    allocate(velocw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1479')
    if (ier /= 0) stop 'Error allocating array velocw_poroelastic'
    allocate(accelw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1480')
    if (ier /= 0) stop 'Error allocating array accelw_poroelastic'
    displs_poroelastic(:,:) = 0.0_CUSTOM_REAL; velocs_poroelastic(:,:) = 0.0_CUSTOM_REAL
    accels_poroelastic(:,:) = 0.0_CUSTOM_REAL
    displw_poroelastic(:,:) = 0.0_CUSTOM_REAL; velocw_poroelastic(:,:) = 0.0_CUSTOM_REAL
    accelw_poroelastic(:,:) = 0.0_CUSTOM_REAL

    allocate(rmass_solid_poroelastic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1481')
    if (ier /= 0) stop 'Error allocating array rmass_solid_poroelastic'
    allocate(rmass_fluid_poroelastic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1482')
    if (ier /= 0) stop 'Error allocating array rmass_fluid_poroelastic'
    rmass_solid_poroelastic(:) = 0.0_CUSTOM_REAL; rmass_fluid_poroelastic(:) = 0.0_CUSTOM_REAL

    allocate(rhoarraystore(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1483')
    allocate(kappaarraystore(3,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1484')
    allocate(etastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1485')
    allocate(tortstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1486')
    allocate(phistore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1487')
    allocate(permstore(6,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1488')
    allocate(rho_vpI(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1489')
    allocate(rho_vpII(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1490')
    allocate(rho_vsI(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1491')
    if (ier /= 0) stop 'Error allocating array poroelastic properties'
    rhoarraystore(:,:,:,:,:) = 0.0_CUSTOM_REAL; kappaarraystore(:,:,:,:,:) = 0.0_CUSTOM_REAL
    etastore(:,:,:,:) = 0.0_CUSTOM_REAL; tortstore(:,:,:,:) = 0.0_CUSTOM_REAL
    phistore(:,:,:,:) = 0.0_CUSTOM_REAL; permstore(:,:,:,:,:) = 0.0_CUSTOM_REAL
    rho_vpI(:,:,:,:) = 0.0_CUSTOM_REAL; rho_vpII(:,:,:,:) = 0.0_CUSTOM_REAL
    rho_vsI(:,:,:,:) = 0.0_CUSTOM_REAL

    ! needed for kernel computations
    allocate(epsilonsdev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1492')
    allocate(epsilonsdev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1493')
    allocate(epsilonsdev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1494')
    allocate(epsilonsdev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1495')
    allocate(epsilonsdev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1496')
    allocate(epsilonwdev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1497')
    allocate(epsilonwdev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1498')
    allocate(epsilonwdev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1499')
    allocate(epsilonwdev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1500')
    allocate(epsilonwdev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1501')
    if (ier /= 0) stop 'Error allocating array epsilonsdev_xx etc.'
    epsilonsdev_xx(:,:,:,:) = 0.0_CUSTOM_REAL; epsilonsdev_yy(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilonsdev_xy(:,:,:,:) = 0.0_CUSTOM_REAL; epsilonsdev_xz(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilonsdev_yz(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilonwdev_xx(:,:,:,:) = 0.0_CUSTOM_REAL; epsilonwdev_yy(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilonwdev_xy(:,:,:,:) = 0.0_CUSTOM_REAL; epsilonwdev_xz(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilonwdev_yz(:,:,:,:) = 0.0_CUSTOM_REAL

    allocate(epsilons_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1502')
    allocate(epsilonw_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1503')
    if (ier /= 0) stop 'Error allocating array epsilons_trace_over_3 etc.'
    epsilons_trace_over_3(:,:,:,:) = 0.0_CUSTOM_REAL; epsilonw_trace_over_3(:,:,:,:) = 0.0_CUSTOM_REAL

    if (I_should_read_the_database) then
      call h5_read_dataset_collect_hyperslab("offset_nspecporo",offset_nspecporo, (/0/), if_col)
      dsetname = "rmass_solid_poroelastic"
      call h5_read_dataset_collect_hyperslab(dsetname, rmass_solid_poroelastic, (/sum(offset_nglob(0:myrank-1))/),if_col)
      dsetname = "rmass_fluid_poroelastic"
      call h5_read_dataset_collect_hyperslab(dsetname, rmass_fluid_poroelastic, (/sum(offset_nglob(0:myrank-1))/),if_col)
      dsetname = "rhoarraystore"
      call h5_read_dataset_collect_hyperslab(dsetname, rhoarraystore, (/0,0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
      dsetname = "kappaarraystore"
      call h5_read_dataset_collect_hyperslab(dsetname, kappaarraystore, &
                                                  (/0,0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
      dsetname = "etastore"
      call h5_read_dataset_collect_hyperslab(dsetname, etastore, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
      dsetname = "tortstore"
      call h5_read_dataset_collect_hyperslab(dsetname, tortstore, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
      dsetname = "permstore"
      call h5_read_dataset_collect_hyperslab(dsetname, permstore, (/0,0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
      dsetname = "phistore"
      call h5_read_dataset_collect_hyperslab(dsetname, phistore, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
      dsetname = "rho_vpI"
      call h5_read_dataset_collect_hyperslab(dsetname, rho_vpI, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
      dsetname = "rho_vpII"
      call h5_read_dataset_collect_hyperslab(dsetname, rho_vpII, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
      dsetname = "rho_vsI"
      call h5_read_dataset_collect_hyperslab(dsetname, rho_vsI, (/0,0,0,sum(offset_nspecporo(0:myrank-1))/),if_col)
    endif
    if (size(rmass_solid_poroelastic) > 0) &
      call bcast_all_cr_for_database(rmass_solid_poroelastic(1), size(rmass_solid_poroelastic))
    if (size(rmass_fluid_poroelastic) > 0) &
      call bcast_all_cr_for_database(rmass_fluid_poroelastic(1), size(rmass_fluid_poroelastic))
    if (size(rhoarraystore) > 0) &
      call bcast_all_cr_for_database(rhoarraystore(1,1,1,1,1), size(rhoarraystore))
    if (size(kappaarraystore) > 0) &
      call bcast_all_cr_for_database(kappaarraystore(1,1,1,1,1), size(kappaarraystore))
    if (size(etastore) > 0) &
      call bcast_all_cr_for_database(etastore(1,1,1,1), size(etastore))
    if (size(tortstore) > 0) &
      call bcast_all_cr_for_database(tortstore(1,1,1,1), size(tortstore))
    if (size(permstore) > 0) &
      call bcast_all_cr_for_database(permstore(1,1,1,1,1), size(permstore))
    if (size(phistore) > 0) &
      call bcast_all_cr_for_database(phistore(1,1,1,1), size(phistore))
    if (size(rho_vpI) > 0) &
      call bcast_all_cr_for_database(rho_vpI(1,1,1,1), size(rho_vpI))
    if (size(rho_vpII) > 0) &
      call bcast_all_cr_for_database(rho_vpII(1,1,1,1), size(rho_vpII))
    if (size(rho_vsI) > 0) &
      call bcast_all_cr_for_database(rho_vsI(1,1,1,1), size(rho_vsI))
  else
    ! dummy allocations (needed for subroutine arguments)
    allocate(rhoarraystore(2,1,1,1,1), &
             kappaarraystore(3,1,1,1,1), &
             etastore(1,1,1,1), &
             tortstore(1,1,1,1), &
             phistore(1,1,1,1), &
             permstore(6,1,1,1,1), &
             rho_vpI(1,1,1,1), &
             rho_vpII(1,1,1,1), &
             rho_vsI(1,1,1,1))
  endif

  ! checks simulation types are valid
  if ((.not. ACOUSTIC_SIMULATION) .and. &
      (.not. ELASTIC_SIMULATION) .and. &
      (.not. POROELASTIC_SIMULATION)) then
    call exit_mpi(myrank,'Error no simulation type defined')
  endif

  ! C-PML absorbing boundary conditions
  ! we allocate this array even when PMLs are absent because we need it in logical tests in "if" statements
  allocate(is_CPML(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1504')
  if (ier /= 0) stop 'Error allocating array is_CPML'
! make sure there are no PMLs by default,
! and then below if NSPEC_CPML > 0 we will read the real flags for this mesh from the disk
  is_CPML(:) = .false.
  NSPEC_CPML = 0

  if (PML_CONDITIONS) then
    if (I_should_read_the_database) then
      call h5_read_dataset_collect_hyperslab("offset_nspeccpml",offset_nspeccpml, (/0/), if_col)
      dsetname = "nspec_cpml"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, NSPEC_CPML, (/myrank/),if_col)
      dsetname = "CPML_width_x"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, CPML_width_x, (/myrank/),if_col)
      dsetname = "CPML_width_y"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, CPML_width_y, (/myrank/),if_col)
      dsetname = "CPML_width_z"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, CPML_width_z, (/myrank/),if_col)
      dsetname = "min_distance_between_CPML_parameter"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, min_distance_between_CPML_parameter, (/myrank/),if_col)
    endif
    call bcast_all_i_for_database(NSPEC_CPML, 1)
    call bcast_all_cr_for_database(CPML_width_x, 1)
    call bcast_all_cr_for_database(CPML_width_y, 1)
    call bcast_all_cr_for_database(CPML_width_z, 1)
    call bcast_all_cr_for_database(min_distance_between_CPML_parameter, 1)

    if (sum(offset_nspeccpml) > 0) then
      allocate(CPML_regions(NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1505')
      if (ier /= 0) stop 'Error allocating array CPML_regions'
      allocate(CPML_to_spec(NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1506')
      if (ier /= 0) stop 'Error allocating array CPML_to_spec'
      CPML_regions(:) = 0; CPML_to_spec(:) = 0

      allocate(d_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1507')
      if (ier /= 0) stop 'Error allocating array d_store_x'
      allocate(d_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1508')
      if (ier /= 0) stop 'Error allocating array d_store_y'
      allocate(d_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1509')
      if (ier /= 0) stop 'Error allocating array d_store_z'
      allocate(K_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1510')
      if (ier /= 0) stop 'Error allocating array K_store_x'
      allocate(K_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1511')
      if (ier /= 0) stop 'Error allocating array K_store_y'
      allocate(K_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1512')
      if (ier /= 0) stop 'Error allocating array K_store_z'
      allocate(alpha_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1513')
      if (ier /= 0) stop 'Error allocating array alpha_store'
      allocate(alpha_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1514')
      if (ier /= 0) stop 'Error allocating array alpha_store'
      allocate(alpha_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1515')
      if (ier /= 0) stop 'Error allocating array alpha_store'
      d_store_x(:,:,:,:) = 0.0_CUSTOM_REAL; d_store_y(:,:,:,:) = 0.0_CUSTOM_REAL; d_store_z(:,:,:,:) = 0.0_CUSTOM_REAL
      K_store_x(:,:,:,:) = 0.0_CUSTOM_REAL; K_store_y(:,:,:,:) = 0.0_CUSTOM_REAL; K_store_z(:,:,:,:) = 0.0_CUSTOM_REAL
      alpha_store_x(:,:,:,:) = 0.0_CUSTOM_REAL
      alpha_store_y(:,:,:,:) = 0.0_CUSTOM_REAL
      alpha_store_z(:,:,:,:) = 0.0_CUSTOM_REAL

      if (I_should_read_the_database) then
        dsetname = "CPML_regions"
        call h5_read_dataset_collect_hyperslab(dsetname, CPML_regions, (/sum(offset_nspeccpml(0:myrank-1))/),if_col)
        dsetname = "CPML_to_spec"
        call h5_read_dataset_collect_hyperslab(dsetname, CPML_to_spec, (/sum(offset_nspeccpml(0:myrank-1))/),if_col)
        dsetname = "is_CPML"
        call h5_read_dataset_collect_hyperslab(dsetname, is_CPML, (/sum(offset_nspec_ab(0:myrank-1))/),if_col)
        dsetname = "d_store_x"
        call h5_read_dataset_collect_hyperslab(dsetname, d_store_x, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
        dsetname = "d_store_y"
        call h5_read_dataset_collect_hyperslab(dsetname, d_store_y, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
        dsetname = "d_store_z"
        call h5_read_dataset_collect_hyperslab(dsetname, d_store_z, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
        dsetname = "k_store_x"
        call h5_read_dataset_collect_hyperslab(dsetname, k_store_x, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
        dsetname = "k_store_y"
        call h5_read_dataset_collect_hyperslab(dsetname, k_store_y, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
        dsetname = "k_store_z"
        call h5_read_dataset_collect_hyperslab(dsetname, k_store_z, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
        dsetname = "alpha_store_x"
        call h5_read_dataset_collect_hyperslab(dsetname, alpha_store_x, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
        dsetname = "alpha_store_y"
        call h5_read_dataset_collect_hyperslab(dsetname, alpha_store_y, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
        dsetname = "alpha_store_z"
        call h5_read_dataset_collect_hyperslab(dsetname, alpha_store_z, (/0,0,0,sum(offset_nspeccpml(0:myrank-1))/),if_col)
      endif

      if (size(CPML_regions) > 0) call bcast_all_i_for_database(CPML_regions(1), size(CPML_regions))
      if (size(CPML_to_spec) > 0) call bcast_all_i_for_database(CPML_to_spec(1), size(CPML_to_spec))
      if (size(is_CPML) > 0) call bcast_all_l_for_database(is_CPML(1), size(is_CPML))
      if (size(d_store_x) > 0) call bcast_all_cr_for_database(d_store_x(1,1,1,1), size(d_store_x))
      if (size(d_store_y) > 0) call bcast_all_cr_for_database(d_store_y(1,1,1,1), size(d_store_y))
      if (size(d_store_z) > 0) call bcast_all_cr_for_database(d_store_z(1,1,1,1), size(d_store_z))
      if (size(k_store_x) > 0) call bcast_all_cr_for_database(k_store_x(1,1,1,1), size(k_store_x))
      if (size(k_store_y) > 0) call bcast_all_cr_for_database(k_store_y(1,1,1,1), size(k_store_y))
      if (size(k_store_z) > 0) call bcast_all_cr_for_database(k_store_z(1,1,1,1), size(k_store_z))
      if (size(alpha_store_x) > 0) call bcast_all_cr_for_database(alpha_store_x(1,1,1,1), size(alpha_store_x))
      if (size(alpha_store_y) > 0) call bcast_all_cr_for_database(alpha_store_y(1,1,1,1), size(alpha_store_y))
      if (size(alpha_store_z) > 0) call bcast_all_cr_for_database(alpha_store_z(1,1,1,1), size(alpha_store_z))

      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        ! acoustic
        if (I_should_read_the_database) then
          call h5_read_dataset_collect_hyperslab( &
                           "offset_nglob_interface_PML_acoustic",offset_nglob_interface_PML_acoustic, (/0/), if_col)
          dsetname = "nglob_interface_PML_acoustic"
          call h5_read_dataset_scalar_collect_hyperslab(dsetname, nglob_interface_PML_acoustic, (/myrank/),if_col)
        endif
        call bcast_all_i_for_database(nglob_interface_PML_acoustic, 1)
        ! elastic
        if (I_should_read_the_database) then
          call h5_read_dataset_collect_hyperslab( &
                            "offset_nglob_interface_PML_elastic",offset_nglob_interface_PML_elastic, (/0/), if_col)
          dsetname = "nglob_interface_PML_elastic"
          call h5_read_dataset_scalar_collect_hyperslab(dsetname, nglob_interface_PML_elastic, (/myrank/),if_col)
        endif
        call bcast_all_i_for_database(nglob_interface_PML_elastic, 1)
        ! acoustic allocation
        if (sum(offset_nglob_interface_PML_acoustic) > 0) then
          allocate(points_interface_PML_acoustic(nglob_interface_PML_acoustic),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 1516')
          if (ier /= 0) stop 'Error allocating array points_interface_PML_acoustic'
          points_interface_PML_acoustic(:) = 0

          if (I_should_read_the_database) then
            dsetname = "points_interface_PML_acoustic"
            call h5_read_dataset_collect_hyperslab(dsetname, &
                points_interface_PML_acoustic, (/sum(offset_nglob_interface_PML_acoustic(0:myrank-1))/),if_col)
          endif
          if (size(points_interface_PML_acoustic) > 0) &
            call bcast_all_i_for_database(points_interface_PML_acoustic(1), size(points_interface_PML_acoustic))
        endif
        ! elastic allocation
        if (sum(offset_nglob_interface_PML_elastic) > 0) then
          allocate(points_interface_PML_elastic(nglob_interface_PML_elastic),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 1517')
          if (ier /= 0) stop 'Error allocating array points_interface_PML_elastic'
          points_interface_PML_elastic(:) = 0

          if (I_should_read_the_database) then
            dsetname = "points_interface_PML_elastic"
            call h5_read_dataset_collect_hyperslab(dsetname, &
                points_interface_PML_acoustic, (/sum(offset_nglob_interface_PML_elastic(0:myrank-1))/),if_col)
          endif
          if (size(points_interface_PML_elastic) > 0) &
            call bcast_all_i_for_database(points_interface_PML_elastic(1), size(points_interface_PML_elastic))
        endif
      endif
    endif
  endif

  ! absorbing boundary surface
  if (I_should_read_the_database) then
    dsetname = "num_abs_boundary_faces"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, num_abs_boundary_faces, (/myrank/),if_col)
  endif
  call bcast_all_i_for_database(num_abs_boundary_faces, 1)

  ! checks
  if (num_abs_boundary_faces < 0) then
    print *,'read_mesh_databases: reading in negative num_abs_boundary_faces ',num_abs_boundary_faces,'...resetting to zero'
    num_abs_boundary_faces = 0
  endif
  allocate(abs_boundary_ispec(num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1518')
  allocate(abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1519')
  allocate(abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1520')
  allocate(abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1521')
  if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'
  abs_boundary_ispec(:) = 0; abs_boundary_ijk(:,:,:) = 0
  abs_boundary_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL; abs_boundary_normal(:,:,:) = 0.0_CUSTOM_REAL

  if (I_should_read_the_database) then
    call h5_read_dataset_collect_hyperslab("offset_num_abs_boundary_faces",offset_num_abs_boundary_faces, (/0/), if_col)
  endif
  call bcast_all_i_array_for_database(offset_num_abs_boundary_faces, size(offset_num_abs_boundary_faces))

  if (sum(offset_num_abs_boundary_faces) > 0) then
    if (I_should_read_the_database) then
      call h5_read_dataset_collect_hyperslab("offset_nglob_xy",offset_nglob_xy, (/0/), if_col)
    endif
    call bcast_all_i_array_for_database(offset_nglob_xy, size(offset_nglob_xy))

    if (I_should_read_the_database) then
      dsetname = "abs_boundary_ispec"
      call h5_read_dataset_collect_hyperslab(dsetname, abs_boundary_ispec, &
              (/sum(offset_num_abs_boundary_faces(0:myrank-1))/),if_col)
      dsetname = "abs_boundary_ijk"
      call h5_read_dataset_collect_hyperslab(dsetname, abs_boundary_ijk, &
              (/0,0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), if_col)
      dsetname = "abs_boundary_jacobian2Dw"
      call h5_read_dataset_collect_hyperslab(dsetname, abs_boundary_jacobian2Dw, &
              (/0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), if_col)
      dsetname = "abs_boundary_normal"
      call h5_read_dataset_collect_hyperslab(dsetname, abs_boundary_normal, &
              (/0,0,sum(offset_num_abs_boundary_faces(0:myrank-1))/), if_col)
    endif
    if (size(abs_boundary_ispec) > 0) &
      call bcast_all_i_for_database(abs_boundary_ispec(1), size(abs_boundary_ispec))
    if (size(abs_boundary_ijk) > 0) &
      call bcast_all_i_for_database(abs_boundary_ijk(1,1,1), size(abs_boundary_ijk))
    if (size(abs_boundary_jacobian2Dw) > 0) &
      call bcast_all_cr_for_database(abs_boundary_jacobian2Dw(1,1), size(abs_boundary_jacobian2Dw))
    if (size(abs_boundary_normal) > 0) &
      call bcast_all_cr_for_database(abs_boundary_normal(1,1,1), size(abs_boundary_normal))

    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! store mass matrix contributions
      if (ELASTIC_SIMULATION) then
        if (I_should_read_the_database) then
          dsetname = "rmassx"
          call h5_read_dataset_collect_hyperslab(dsetname, rmassx(1:offset_nglob_xy(myrank)), &
                                                          (/sum(offset_nglob_xy(0:myrank-1))/),if_col)
          dsetname = "rmassy"
          call h5_read_dataset_collect_hyperslab(dsetname, rmassy(1:offset_nglob_xy(myrank)), &
                                                          (/sum(offset_nglob_xy(0:myrank-1))/),if_col)
          dsetname = "rmassz"
          call h5_read_dataset_collect_hyperslab(dsetname, rmassz(1:offset_nglob_xy(myrank)), &
                                                          (/sum(offset_nglob_xy(0:myrank-1))/),if_col)
        endif
        if (size(rmassx) > 0) call bcast_all_cr_for_database(rmassx(1), size(rmassx))
        if (size(rmassy) > 0) call bcast_all_cr_for_database(rmassy(1), size(rmassy))
        if (size(rmassz) > 0) call bcast_all_cr_for_database(rmassz(1), size(rmassz))
      endif
      if (ACOUSTIC_SIMULATION) then
        if (I_should_read_the_database) then
          dsetname = "rmassz_acoustic"
          call h5_read_dataset_collect_hyperslab(dsetname, rmassz_acoustic(1:offset_nglob_xy(myrank)), &
                                                          (/sum(offset_nglob_xy(0:myrank-1))/),if_col)
        endif
        if (size(rmassz_acoustic) > 0) call bcast_all_cr_for_database(rmassz_acoustic(1), size(rmassz_acoustic))
      endif
    endif
  endif

  ! boundaries
  if (I_should_read_the_database) then
    dsetname = "nspec2D_xmin"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec2D_xmin, (/myrank/),if_col)
    dsetname = "nspec2D_xmax"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec2D_xmax, (/myrank/),if_col)
    dsetname = "nspec2D_ymin"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec2D_ymin, (/myrank/),if_col)
    dsetname = "nspec2D_ymax"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec2D_ymax, (/myrank/),if_col)
    dsetname = "NSPEC2D_BOTTOM"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec2D_bottom, (/myrank/),if_col)
    dsetname = "NSPEC2D_TOP"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec2D_top, (/myrank/),if_col)
  endif
  call bcast_all_i_for_database(nspec2D_xmin, 1)
  call bcast_all_i_for_database(nspec2D_xmax, 1)
  call bcast_all_i_for_database(nspec2D_ymin, 1)
  call bcast_all_i_for_database(nspec2D_ymax, 1)
  call bcast_all_i_for_database(NSPEC2D_BOTTOM, 1)
  call bcast_all_i_for_database(NSPEC2D_TOP, 1)

  allocate(ibelm_xmin(nspec2D_xmin),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1522')
  allocate(ibelm_xmax(nspec2D_xmax),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1523')
  allocate(ibelm_ymin(nspec2D_ymin),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1524')
  allocate(ibelm_ymax(nspec2D_ymax),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1525')
  allocate(ibelm_bottom(NSPEC2D_BOTTOM),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1526')
  allocate(ibelm_top(NSPEC2D_TOP),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1527')
  if (ier /= 0) stop 'Error allocating arrays ibelm_xmin,ibelm_xmax etc.'
  ibelm_xmin(:) = 0; ibelm_xmax(:) = 0
  ibelm_ymin(:) = 0; ibelm_ymax(:) = 0
  ibelm_bottom(:) = 0; ibelm_top(:) = 0

  if (I_should_read_the_database) then
    if (sum(offset_nspec2D_xmin) > 0) then
      dsetname = "ibelm_xmin"
      call h5_read_dataset_collect_hyperslab(dsetname, ibelm_xmin, (/sum(offset_nspec2D_xmin(0:myrank-1))/),if_col)
    endif
    if (sum(offset_nspec2D_xmax) > 0) then
      dsetname = "ibelm_xmax"
      call h5_read_dataset_collect_hyperslab(dsetname, ibelm_xmax, (/sum(offset_nspec2D_xmax(0:myrank-1))/),if_col)
    endif
    if (sum(offset_nspec2D_ymin) > 0) then
      dsetname = "ibelm_ymin"
      call h5_read_dataset_collect_hyperslab(dsetname, ibelm_ymin, (/sum(offset_nspec2D_ymin(0:myrank-1))/),if_col)
    endif
    if (sum(offset_nspec2D_ymax) > 0) then
      dsetname = "ibelm_ymax"
      call h5_read_dataset_collect_hyperslab(dsetname, ibelm_ymax, (/sum(offset_nspec2D_ymax(0:myrank-1))/),if_col)
    endif
    if (sum(offset_nspec2D_bottom_ext) > 0) then
      dsetname = "ibelm_bottom"
      call h5_read_dataset_collect_hyperslab(dsetname, ibelm_bottom, (/sum(offset_nspec2D_bottom_ext(0:myrank-1))/),if_col)
    endif
    if (sum(offset_nspec2D_top_ext) > 0) then
      dsetname = "ibelm_top"
      call h5_read_dataset_collect_hyperslab(dsetname, ibelm_top, (/sum(offset_nspec2D_top_ext(0:myrank-1))/),if_col)
    endif
  endif
  if (size(ibelm_xmin) > 0) call bcast_all_i_for_database(ibelm_xmin(1), size(ibelm_xmin))
  if (size(ibelm_xmax) > 0) call bcast_all_i_for_database(ibelm_xmax(1), size(ibelm_xmax))
  if (size(ibelm_ymin) > 0) call bcast_all_i_for_database(ibelm_ymin(1), size(ibelm_ymin))
  if (size(ibelm_ymax) > 0) call bcast_all_i_for_database(ibelm_ymax(1), size(ibelm_ymax))
  if (size(ibelm_bottom) > 0) call bcast_all_i_for_database(ibelm_bottom(1), size(ibelm_bottom))
  if (size(ibelm_top) > 0) call bcast_all_i_for_database(ibelm_top(1), size(ibelm_top))

  ! free surface
  if (I_should_read_the_database) then
    dsetname = "num_free_surface_faces"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, num_free_surface_faces, (/myrank/),if_col)
  endif
  call bcast_all_i_for_database(num_free_surface_faces, 1)

  allocate(free_surface_ispec(num_free_surface_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1528')
  allocate(free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1529')
  allocate(free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1530')
  allocate(free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1531')
  if (ier /= 0) stop 'Error allocating arrays free_surface_ispec etc.'
  free_surface_ispec(:) = 0; free_surface_ijk(:,:,:) = 0
  free_surface_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL; free_surface_normal(:,:,:) = 0.0_CUSTOM_REAL

  if (I_should_read_the_database) then
    call h5_read_dataset_collect_hyperslab("offset_num_free_surface_faces",offset_num_free_surface_faces, (/0/), if_col)
  endif
  call bcast_all_i_array_for_database(offset_num_free_surface_faces, size(offset_num_free_surface_faces))

  if (sum(offset_num_free_surface_faces) > 0) then
    if (I_should_read_the_database) then
      dsetname = "free_surface_ispec"
      call h5_read_dataset_collect_hyperslab(dsetname, free_surface_ispec, &
              (/sum(offset_num_free_surface_faces(0:myrank-1))/),if_col)
      dsetname = "free_surface_ijk"
      call h5_read_dataset_collect_hyperslab(dsetname, free_surface_ijk, &
              (/0,0,sum(offset_num_free_surface_faces(0:myrank-1))/),if_col)
      dsetname = "free_surface_jacobian2Dw"
      call h5_read_dataset_collect_hyperslab(dsetname, free_surface_jacobian2Dw, &
              (/0,sum(offset_num_free_surface_faces(0:myrank-1))/),if_col)
      dsetname = "free_surface_normal"
      call h5_read_dataset_collect_hyperslab(dsetname, free_surface_normal, &
              (/0,0,sum(offset_num_free_surface_faces(0:myrank-1))/),if_col)
    endif
    if (size(free_surface_ispec) > 0) &
      call bcast_all_i_for_database(free_surface_ispec(1), size(free_surface_ispec))
    if (size(free_surface_ijk) > 0) &
      call bcast_all_i_for_database(free_surface_ijk(1,1,1), size(free_surface_ijk))
    if (size(free_surface_jacobian2Dw) > 0) &
      call bcast_all_cr_for_database(free_surface_jacobian2Dw(1,1), size(free_surface_jacobian2Dw))
    if (size(free_surface_normal) > 0) &
      call bcast_all_cr_for_database(free_surface_normal(1,1,1), size(free_surface_normal))
  endif

  ! acoustic-elastic coupling surface
  if (I_should_read_the_database) then
    dsetname = "num_coupling_ac_el_faces"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, num_coupling_ac_el_faces, (/myrank/),if_col)
  endif
  call bcast_all_i_for_database(num_coupling_ac_el_faces, 1)

  allocate(coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1532')
  allocate(coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1533')
  allocate(coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1534')
  allocate(coupling_ac_el_ispec(num_coupling_ac_el_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1535')
  if (ier /= 0) stop 'Error allocating array coupling_ac_el_normal etc.'
  coupling_ac_el_ispec(:) = 0; coupling_ac_el_ijk(:,:,:) = 0
  coupling_ac_el_normal(:,:,:) = 0.0_CUSTOM_REAL; coupling_ac_el_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL

  if (I_should_read_the_database) then
    call h5_read_dataset_collect_hyperslab("offset_num_coupling_ac_el_faces",offset_num_coupling_ac_el_faces,(/0/),if_col)
  endif
  call bcast_all_i_array_for_database(offset_num_coupling_ac_el_faces, size(offset_num_coupling_ac_el_faces))

  if (sum(offset_num_coupling_ac_el_faces) > 0) then
    if (I_should_read_the_database) then
      dsetname = "coupling_ac_el_ispec"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_ac_el_ispec, &
            (/sum(offset_num_coupling_ac_el_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_ac_el_ijk"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_ac_el_ijk, &
                (/0,0,sum(offset_num_coupling_ac_el_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_ac_el_jacobian2Dw"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_ac_el_jacobian2Dw, &
               (/0,sum(offset_num_coupling_ac_el_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_ac_el_normal"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_ac_el_normal, &
               (/0,0,sum(offset_num_coupling_ac_el_faces(0:myrank-1))/),if_col)
    endif
    if (size(coupling_ac_el_ispec) > 0) &
      call bcast_all_i_for_database(coupling_ac_el_ispec(1), size(coupling_ac_el_ispec))
    if (size(coupling_ac_el_ijk) > 0) &
      call bcast_all_i_for_database(coupling_ac_el_ijk(1,1,1), size(coupling_ac_el_ijk))
    if (size(coupling_ac_el_jacobian2Dw) > 0) &
      call bcast_all_cr_for_database(coupling_ac_el_jacobian2Dw(1,1), size(coupling_ac_el_jacobian2Dw))
    if (size(coupling_ac_el_normal) > 0) &
      call bcast_all_cr_for_database(coupling_ac_el_normal(1,1,1), size(coupling_ac_el_normal))
  endif

  ! acoustic-poroelastic coupling surface
  if (I_should_read_the_database) then
    dsetname = "num_coupling_ac_po_faces"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, num_coupling_ac_po_faces, (/myrank/),if_col)
  endif
  call bcast_all_i_for_database(num_coupling_ac_po_faces, 1)

  allocate(coupling_ac_po_normal(NDIM,NGLLSQUARE,num_coupling_ac_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1536')
  allocate(coupling_ac_po_jacobian2Dw(NGLLSQUARE,num_coupling_ac_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1537')
  allocate(coupling_ac_po_ijk(3,NGLLSQUARE,num_coupling_ac_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1538')
  allocate(coupling_ac_po_ispec(num_coupling_ac_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1539')
  if (ier /= 0) stop 'Error allocating array coupling_ac_po_normal etc.'
  coupling_ac_po_ispec(:) = 0; coupling_ac_po_ijk(:,:,:) = 0
  coupling_ac_po_normal(:,:,:) = 0.0_CUSTOM_REAL; coupling_ac_po_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL

  if (I_should_read_the_database) then
    call h5_read_dataset_collect_hyperslab("offset_num_coupling_ac_po_faces",offset_num_coupling_ac_po_faces, (/0/), if_col)
  endif
  call bcast_all_i_array_for_database(offset_num_coupling_ac_po_faces,size(offset_num_coupling_ac_po_faces))

  if (sum(offset_num_coupling_ac_po_faces) > 0) then
    if (I_should_read_the_database) then
      dsetname = "coupling_ac_po_ispec"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_ac_po_ispec, &
              (/sum(offset_num_coupling_ac_po_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_ac_po_ijk"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_ac_po_ijk, &
              (/0,0,sum(offset_num_coupling_ac_po_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_ac_po_jacobian2Dw"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_ac_po_jacobian2Dw, &
              (/0,sum(offset_num_coupling_ac_po_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_ac_po_normal"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_ac_po_normal, &
              (/0,0,sum(offset_num_coupling_ac_po_faces(0:myrank-1))/),if_col)
    endif
    if (size(coupling_ac_po_ispec) > 0) &
      call bcast_all_i_for_database(coupling_ac_po_ispec(1), size(coupling_ac_po_ispec))
    if (size(coupling_ac_po_ijk) > 0) &
      call bcast_all_i_for_database(coupling_ac_po_ijk(1,1,1), size(coupling_ac_po_ijk))
    if (size(coupling_ac_po_jacobian2Dw) > 0) &
      call bcast_all_cr_for_database(coupling_ac_po_jacobian2Dw(1,1), size(coupling_ac_po_jacobian2Dw))
    if (size(coupling_ac_po_normal) > 0) &
      call bcast_all_cr_for_database(coupling_ac_po_normal(1,1,1), size(coupling_ac_po_normal))
  endif

  ! elastic-poroelastic coupling surface
  if (I_should_read_the_database) then
    dsetname = "num_coupling_el_po_faces"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, num_coupling_el_po_faces, (/myrank/),if_col)
  endif
  call bcast_all_i_for_database(num_coupling_el_po_faces, 1)

  allocate(coupling_el_po_normal(NDIM,NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1540')
  allocate(coupling_el_po_jacobian2Dw(NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1541')
  allocate(coupling_el_po_ijk(3,NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1542')
  allocate(coupling_po_el_ijk(3,NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1543')
  allocate(coupling_el_po_ispec(num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1544')
  allocate(coupling_po_el_ispec(num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1545')
  if (ier /= 0) stop 'Error allocating array coupling_el_po_normal etc.'
  coupling_el_po_ispec(:) = 0; coupling_el_po_ijk(:,:,:) = 0
  coupling_po_el_ispec(:) = 0; coupling_po_el_ijk(:,:,:) = 0
  coupling_el_po_normal(:,:,:) = 0.0_CUSTOM_REAL; coupling_el_po_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL

  if (I_should_read_the_database) then
    call h5_read_dataset_collect_hyperslab("offset_num_coupling_el_po_faces",offset_num_coupling_el_po_faces, (/0/), if_col)
  endif
  call bcast_all_i_array_for_database(offset_num_coupling_el_po_faces,size(offset_num_coupling_el_po_faces))

  if (sum(offset_num_coupling_el_po_faces) > 0) then
    if (I_should_read_the_database) then
      dsetname = "coupling_el_po_ispec"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_el_po_ispec, &
              (/sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_po_el_ispec"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_po_el_ispec, &
              (/sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_el_po_ijk"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_el_po_ijk, &
              (/0,0,sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_po_el_ijk"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_po_el_ijk, &
              (/0,0,sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_el_po_jacobian2Dw"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_el_po_jacobian2Dw, &
              (/0,sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
      dsetname = "coupling_el_po_normal"
      call h5_read_dataset_collect_hyperslab(dsetname, coupling_el_po_normal, &
              (/0,0,sum(offset_num_coupling_el_po_faces(0:myrank-1))/),if_col)
    endif
    if (size(coupling_el_po_ispec) > 0) &
      call bcast_all_i_for_database(coupling_el_po_ispec(1), size(coupling_el_po_ispec))
    if (size(coupling_po_el_ispec) > 0) &
      call bcast_all_i_for_database(coupling_po_el_ispec(1), size(coupling_po_el_ispec))
    if (size(coupling_el_po_ijk) > 0) &
      call bcast_all_i_for_database(coupling_el_po_ijk(1,1,1), size(coupling_el_po_ijk))
    if (size(coupling_po_el_ijk) > 0) &
      call bcast_all_i_for_database(coupling_po_el_ijk(1,1,1), size(coupling_po_el_ijk))
    if (size(coupling_el_po_jacobian2Dw) > 0) &
      call bcast_all_cr_for_database(coupling_el_po_jacobian2Dw(1,1), size(coupling_el_po_jacobian2Dw))
    if (size(coupling_el_po_normal) > 0) &
      call bcast_all_cr_for_database(coupling_el_po_normal(1,1,1), size(coupling_el_po_normal))
  endif

  ! MPI interfaces
  if (I_should_read_the_database) then
     dsetname = "num_interfaces_ext_mesh"
     call h5_read_dataset_scalar_collect_hyperslab(dsetname, num_interfaces_ext_mesh, (/myrank/),if_col)
  endif
  call bcast_all_i_for_database(num_interfaces_ext_mesh, 1)

  allocate(my_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1546')
  allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1547')
  if (ier /= 0) stop 'Error allocating array my_neighbors_ext_mesh etc.'
  my_neighbors_ext_mesh(:) = -1; nibool_interfaces_ext_mesh(:) = 0

  if (I_should_read_the_database) then
    call h5_read_dataset_collect_hyperslab("offset_num_interfaces_ext_mesh",offset_num_interfaces_ext_mesh, (/0/), if_col)
  endif
  call bcast_all_i_array_for_database(offset_num_interfaces_ext_mesh,size(offset_num_interfaces_ext_mesh))

  if (sum(offset_num_interfaces_ext_mesh) > 0) then
    if (I_should_read_the_database) then
      dsetname = "max_nibool_interfaces_ext_mesh"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, max_nibool_interfaces_ext_mesh, (/myrank/),if_col)
    endif
    call bcast_all_i_for_database(max_nibool_interfaces_ext_mesh, 1)
    allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1548')
    if (ier /= 0) stop 'Error allocating array ibool_interfaces_ext_mesh'
    ibool_interfaces_ext_mesh(:,:) = 0

    if (I_should_read_the_database) then
      dsetname = "my_neighbors_ext_mesh"
      call h5_read_dataset_collect_hyperslab(dsetname, my_neighbors_ext_mesh, &
              (/sum(offset_num_interfaces_ext_mesh(0:myrank-1))/),if_col)
      dsetname = "nibool_interfaces_ext_mesh"
      call h5_read_dataset_collect_hyperslab(dsetname, nibool_interfaces_ext_mesh, &
              (/sum(offset_num_interfaces_ext_mesh(0:myrank-1))/),if_col)
      dsetname = "ibool_interfaces_ext_mesh_dummy"
      call h5_read_dataset_collect_hyperslab(dsetname, ibool_interfaces_ext_mesh, &
              (/0,sum(offset_num_interfaces_ext_mesh(0:myrank-1))/),if_col)
    endif
    if (size(my_neighbors_ext_mesh) > 0) &
      call bcast_all_i_for_database(my_neighbors_ext_mesh(1), size(my_neighbors_ext_mesh))
    if (size(nibool_interfaces_ext_mesh) > 0) &
      call bcast_all_i_for_database(nibool_interfaces_ext_mesh(1), size(nibool_interfaces_ext_mesh))
    if (size(ibool_interfaces_ext_mesh) > 0) &
      call bcast_all_i_for_database(ibool_interfaces_ext_mesh(1,1), size(ibool_interfaces_ext_mesh))
  else
    ! no interfaces
    max_nibool_interfaces_ext_mesh = 0
    ! dummy allocations
    allocate(ibool_interfaces_ext_mesh(0,0),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1549')
  endif

  ! material properties
  if (ELASTIC_SIMULATION .and. ANISOTROPY) then
    if (I_should_read_the_database) then
      call h5_read_dataset_collect_hyperslab("offset_nspec_aniso",offset_nspec_aniso, (/0/), if_col)
      dsetname = "c11store"
      call h5_read_dataset_collect_hyperslab(dsetname, c11store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c12store"
      call h5_read_dataset_collect_hyperslab(dsetname, c12store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c13store"
      call h5_read_dataset_collect_hyperslab(dsetname, c13store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c14store"
      call h5_read_dataset_collect_hyperslab(dsetname, c14store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c15store"
      call h5_read_dataset_collect_hyperslab(dsetname, c15store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c16store"
      call h5_read_dataset_collect_hyperslab(dsetname, c16store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c22store"
      call h5_read_dataset_collect_hyperslab(dsetname, c22store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c23store"
      call h5_read_dataset_collect_hyperslab(dsetname, c23store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c24store"
      call h5_read_dataset_collect_hyperslab(dsetname, c24store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c25store"
      call h5_read_dataset_collect_hyperslab(dsetname, c25store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c26store"
      call h5_read_dataset_collect_hyperslab(dsetname, c26store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c33store"
      call h5_read_dataset_collect_hyperslab(dsetname, c33store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c34store"
      call h5_read_dataset_collect_hyperslab(dsetname, c34store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c35store"
      call h5_read_dataset_collect_hyperslab(dsetname, c35store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c36store"
      call h5_read_dataset_collect_hyperslab(dsetname, c36store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c44store"
      call h5_read_dataset_collect_hyperslab(dsetname, c44store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c45store"
      call h5_read_dataset_collect_hyperslab(dsetname, c45store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c46store"
      call h5_read_dataset_collect_hyperslab(dsetname, c46store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c55store"
      call h5_read_dataset_collect_hyperslab(dsetname, c55store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c56store"
      call h5_read_dataset_collect_hyperslab(dsetname, c56store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
      dsetname = "c66store"
      call h5_read_dataset_collect_hyperslab(dsetname, c66store, (/0,0,0,sum(offset_nspec_aniso(0:myrank-1))/), if_col)
    endif
    if (size(c11store) > 0) call bcast_all_cr_for_database(c11store(1,1,1,1), size(c11store))
    if (size(c12store) > 0) call bcast_all_cr_for_database(c12store(1,1,1,1), size(c12store))
    if (size(c13store) > 0) call bcast_all_cr_for_database(c13store(1,1,1,1), size(c13store))
    if (size(c14store) > 0) call bcast_all_cr_for_database(c14store(1,1,1,1), size(c14store))
    if (size(c15store) > 0) call bcast_all_cr_for_database(c15store(1,1,1,1), size(c15store))
    if (size(c16store) > 0) call bcast_all_cr_for_database(c16store(1,1,1,1), size(c16store))
    if (size(c22store) > 0) call bcast_all_cr_for_database(c22store(1,1,1,1), size(c22store))
    if (size(c23store) > 0) call bcast_all_cr_for_database(c23store(1,1,1,1), size(c23store))
    if (size(c24store) > 0) call bcast_all_cr_for_database(c24store(1,1,1,1), size(c24store))
    if (size(c25store) > 0) call bcast_all_cr_for_database(c25store(1,1,1,1), size(c25store))
    if (size(c26store) > 0) call bcast_all_cr_for_database(c26store(1,1,1,1), size(c26store))
    if (size(c33store) > 0) call bcast_all_cr_for_database(c33store(1,1,1,1), size(c33store))
    if (size(c34store) > 0) call bcast_all_cr_for_database(c34store(1,1,1,1), size(c34store))
    if (size(c35store) > 0) call bcast_all_cr_for_database(c35store(1,1,1,1), size(c35store))
    if (size(c36store) > 0) call bcast_all_cr_for_database(c36store(1,1,1,1), size(c36store))
    if (size(c44store) > 0) call bcast_all_cr_for_database(c44store(1,1,1,1), size(c44store))
    if (size(c45store) > 0) call bcast_all_cr_for_database(c45store(1,1,1,1), size(c45store))
    if (size(c46store) > 0) call bcast_all_cr_for_database(c46store(1,1,1,1), size(c46store))
    if (size(c55store) > 0) call bcast_all_cr_for_database(c55store(1,1,1,1), size(c55store))
    if (size(c56store) > 0) call bcast_all_cr_for_database(c56store(1,1,1,1), size(c56store))
    if (size(c66store) > 0) call bcast_all_cr_for_database(c66store(1,1,1,1), size(c66store))
  endif

  ! inner / outer elements
  allocate(ispec_is_inner(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1550')
  if (ier /= 0) stop 'Error allocating array ispec_is_inner'
  ispec_is_inner(:) = .false.

  if (I_should_read_the_database) then
    dsetname = "ispec_is_inner"
    call h5_read_dataset_collect_hyperslab(dsetname, ispec_is_inner,(/sum(offset_nspec(0:myrank-1))/),if_col)
  endif
  if (size(ispec_is_inner) > 0) call bcast_all_l_for_database(ispec_is_inner(1), size(ispec_is_inner))

  if (ACOUSTIC_SIMULATION) then
    if (I_should_read_the_database) then
      dsetname = "nspec_inner_acoustic"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec_inner_acoustic,(/myrank/),if_col)
      dsetname = "nspec_outer_acoustic"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec_outer_acoustic,(/myrank/),if_col)
      dsetname = "num_phase_ispec_acoustic"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, num_phase_ispec_acoustic,(/myrank/),if_col)
    endif
    call bcast_all_i_for_database(nspec_inner_acoustic, 1)
    call bcast_all_i_for_database(nspec_outer_acoustic, 1)
    call bcast_all_i_for_database(num_phase_ispec_acoustic, 1)
    if (num_phase_ispec_acoustic < 0) stop 'Error acoustic simulation: num_phase_ispec_acoustic is < zero'

    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1551')
    if (ier /= 0) stop 'Error allocating array phase_ispec_inner_acoustic'
    phase_ispec_inner_acoustic(:,:) = 0

    call h5_read_dataset_collect_hyperslab("offset_num_phase_ispec_acoustic", &
            offset_num_phase_ispec_acoustic, (/0/), if_col)
    if (sum(offset_num_phase_ispec_acoustic) > 0) then
     if (I_should_read_the_database) then
        dsetname = "phase_ispec_inner_acoustic"
        call h5_read_dataset_collect_hyperslab(dsetname, phase_ispec_inner_acoustic, &
                (/sum(offset_num_phase_ispec_acoustic(0:myrank-1)),0/),if_col)
      endif
      if (size(phase_ispec_inner_acoustic) > 0) &
            call bcast_all_i_for_database(phase_ispec_inner_acoustic(1,1), size(phase_ispec_inner_acoustic))
    endif
  endif

  if (ELASTIC_SIMULATION) then
    if (I_should_read_the_database) then
      dsetname = "nspec_inner_elastic"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec_inner_elastic,(/myrank/),if_col)
      dsetname = "nspec_outer_elastic"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec_outer_elastic,(/myrank/),if_col)
      dsetname = "num_phase_ispec_elastic"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, num_phase_ispec_elastic,(/myrank/),if_col)
    endif
    call bcast_all_i_for_database(nspec_inner_elastic, 1)
    call bcast_all_i_for_database(nspec_outer_elastic, 1)
    call bcast_all_i_for_database(num_phase_ispec_elastic, 1)
    if (num_phase_ispec_elastic < 0) stop 'Error elastic simulation: num_phase_ispec_elastic is < zero'

    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1552')
    if (ier /= 0) stop 'Error allocating array phase_ispec_inner_elastic'
    phase_ispec_inner_elastic(:,:) = 0

    if (I_should_read_the_database) then
      call h5_read_dataset_collect_hyperslab("offset_num_phase_ispec_elastic",offset_num_phase_ispec_elastic, (/0/), if_col)
    endif
    call bcast_all_i_array_for_database(offset_num_phase_ispec_elastic,size(offset_num_phase_ispec_elastic))

    if (sum(offset_num_phase_ispec_elastic) > 0) then
      if (I_should_read_the_database) then
         dsetname = "phase_ispec_inner_elastic"
        call h5_read_dataset_collect_hyperslab(dsetname, phase_ispec_inner_elastic, &
                (/sum(offset_num_phase_ispec_elastic(0:myrank-1)),0/),if_col)
      endif
      if (size(phase_ispec_inner_elastic) > 0) &
        call bcast_all_i_for_database(phase_ispec_inner_elastic(1,1), size(phase_ispec_inner_elastic))
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    if (I_should_read_the_database) then
      dsetname = "nspec_inner_poroelastic"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec_inner_poroelastic,(/myrank/),if_col)
      dsetname = "nspec_outer_poroelastic"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, nspec_outer_poroelastic,(/myrank/),if_col)
      dsetname = "num_phase_ispec_poroelastic"
      call h5_read_dataset_scalar_collect_hyperslab(dsetname, num_phase_ispec_poroelastic, (/myrank/),if_col)
    endif
    call bcast_all_i_for_database(nspec_inner_poroelastic, 1)
    call bcast_all_i_for_database(nspec_outer_poroelastic, 1)
    call bcast_all_i_for_database(num_phase_ispec_poroelastic, 1)
    if (num_phase_ispec_poroelastic < 0) stop 'Error poroelastic simulation: num_phase_ispec_poroelastic is < zero'

    allocate( phase_ispec_inner_poroelastic(num_phase_ispec_poroelastic,2),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1553')
    if (ier /= 0) stop 'Error allocating array phase_ispec_inner_poroelastic'
    phase_ispec_inner_poroelastic(:,:) = 0

    if (I_should_read_the_database) then
      call h5_read_dataset_collect_hyperslab("offset_num_phase_ispec_poroelastic", &
            offset_num_phase_ispec_poroelastic, (/0/), if_col)
    endif
    call bcast_all_i_array_for_database(offset_num_phase_ispec_poroelastic,size(offset_num_phase_ispec_poroelastic))
    if (sum(offset_num_phase_ispec_poroelastic) > 0) then
      if (I_should_read_the_database) then
        dsetname = "phase_ispec_inner_poroelastic"
        call h5_read_dataset_collect_hyperslab(dsetname, phase_ispec_inner_poroelastic, &
              (/sum(offset_num_phase_ispec_poroelastic(0:myrank-1)),0/),if_col)
      endif
      if (size(phase_ispec_inner_poroelastic) > 0) &
        call bcast_all_i_for_database(phase_ispec_inner_poroelastic(1,1), size(phase_ispec_inner_poroelastic))
    endif
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! acoustic domain colors
    if (ACOUSTIC_SIMULATION) then
      if (I_should_read_the_database) then
        dsetname = "num_colors_outer_acoustic"
        call h5_read_dataset_scalar_collect_hyperslab(dsetname,num_colors_outer_acoustic,(/myrank/),if_col)
        dsetname = "num_colors_inner_acoustic"
        call h5_read_dataset_scalar_collect_hyperslab(dsetname,num_colors_inner_acoustic,(/myrank/),if_col)
      endif
      call bcast_all_i_for_database(num_colors_outer_acoustic, 1)
      call bcast_all_i_for_database(num_colors_inner_acoustic, 1)

      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1554')
      if (ier /= 0) stop 'Error allocating num_elem_colors_acoustic array'
      num_elem_colors_acoustic(:) = 0

      if (I_should_read_the_database) then
        call h5_read_dataset_collect_hyperslab( &
                         "offset_num_colors_acoustic",offset_num_colors_acoustic, (/0/), if_col)
        dsetname = "num_elem_colors_acoustic"
        call h5_read_dataset_collect_hyperslab(dsetname, &
                num_elem_colors_acoustic,(/sum(offset_num_colors_acoustic(0:myrank-1))/),if_col)
      endif
      if (size(num_elem_colors_acoustic) > 0) &
        call bcast_all_i_for_database(num_elem_colors_acoustic(1), size(num_elem_colors_acoustic))
    endif
    ! elastic domain colors
    if (ELASTIC_SIMULATION) then
      if (I_should_read_the_database) then
        dsetname = "num_colors_outer_elastic"
        call h5_read_dataset_scalar_collect_hyperslab(dsetname,num_colors_outer_elastic,(/myrank/),if_col)
        dsetname = "num_colors_inner_elastic"
        call h5_read_dataset_scalar_collect_hyperslab(dsetname,num_colors_inner_elastic,(/myrank/),if_col)
      endif
      call bcast_all_i_for_database(num_colors_outer_elastic, 1)
      call bcast_all_i_for_database(num_colors_inner_elastic, 1)

      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1555')
      if (ier /= 0) stop 'Error allocating num_elem_colors_elastic array'
      num_elem_colors_elastic(:) = 0

      if (I_should_read_the_database) then
        call h5_read_dataset_collect_hyperslab( &
                  "offset_num_colors_elastic",offset_num_colors_elastic, (/0/), if_col)
        dsetname = "num_elem_colors_elastic"
        call h5_read_dataset_collect_hyperslab(dsetname, num_elem_colors_elastic, &
                (/sum(offset_num_colors_elastic(0:myrank-1))/),if_col)
      endif
      if (size(num_elem_colors_elastic) > 0) &
        call bcast_all_i_for_database(num_elem_colors_elastic(1), size(num_elem_colors_elastic))
    endif
  else
    ! allocates dummy arrays
    if (ACOUSTIC_SIMULATION) then
      num_colors_outer_acoustic = 0
      num_colors_inner_acoustic = 0
      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1556')
      if (ier /= 0) stop 'Error allocating num_elem_colors_acoustic array'
    endif
    if (ELASTIC_SIMULATION) then
      num_colors_outer_elastic = 0
      num_colors_inner_elastic = 0
      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1557')
      if (ier /= 0) stop 'Error allocating num_elem_colors_elastic array'
    endif
  endif

  ! for mesh surface
  allocate(ispec_is_surface_external_mesh(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1558')
  allocate(iglob_is_surface_external_mesh(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1559')
  if (ier /= 0) stop 'error allocating array for mesh surface'
  ispec_is_surface_external_mesh(:) = .false.; iglob_is_surface_external_mesh(:) = .false.

  ! determines model surface
  ! returns surface points/elements in ispec_is_surface_external_mesh / iglob_is_surface_external_mesh
  ! and number of faces in nfaces_surface
  ! used for receiver detection, movie files and shakemaps
  if (I_should_read_the_database) then
    dsetname = "nfaces_surface"
    call h5_read_dataset_scalar_collect_hyperslab(dsetname, nfaces_surface, (/myrank/), if_col)
    dsetname = "ispec_is_surface_external_mesh"
    call h5_read_dataset_collect_hyperslab(dsetname, ispec_is_surface_external_mesh, &
            (/sum(offset_nspec_ab(0:myrank-1))/), if_col)
    dsetname = "iglob_is_surface_external_mesh"
    call h5_read_dataset_collect_hyperslab(dsetname, iglob_is_surface_external_mesh, &
            (/sum(offset_nglob_ab(0:myrank-1))/), if_col)
  endif
  call bcast_all_i_for_database(nfaces_surface, 1)
  call bcast_all_l_for_database(ispec_is_surface_external_mesh(1), size(ispec_is_surface_external_mesh))
  call bcast_all_l_for_database(iglob_is_surface_external_mesh(1), size(iglob_is_surface_external_mesh))

  ! done reading database
  if (I_should_read_the_database) then
    ! close hdf5
    call h5_close_file()
    call h5_finalize()
  endif

  ! checks MPI interface arrays
  if (NPROC > 1) then
    ! checks if valid neighbors array
    if (any(my_neighbors_ext_mesh(:) < 0) .or. any(my_neighbors_ext_mesh(:) >= NPROC)) then
      print *,'Error: rank ',myrank,'has invalid MPI interface my_neighbors: ',my_neighbors_ext_mesh(:)
      stop 'Invalid MPI interface my_neighbors'
    endif
    ! checks if valid nibool_interfaces array
    if (any(nibool_interfaces_ext_mesh(:) <= 0) .or. any(nibool_interfaces_ext_mesh(:) >= NGLOB_AB)) then
      print *,'Error: rank ',myrank,'has invalid MPI interface nibool_interfaces: ',nibool_interfaces_ext_mesh(:)
      stop 'Invalid MPI interface nibool_interfaces'
    endif
    ! checks if valid ibool entries
    do itest = 1,num_interfaces_ext_mesh
      if (any(ibool_interfaces_ext_mesh(1:nibool_interfaces_ext_mesh(itest),itest) <= 0) .or. &
          any(ibool_interfaces_ext_mesh(1:nibool_interfaces_ext_mesh(itest),itest) > NGLOB_AB)) then
        print *,'Error: rank ',myrank,'has invalid ibool_interfaces at interface ',itest, &
                'entries: ', ibool_interfaces_ext_mesh(1:nibool_interfaces_ext_mesh(itest),itest)
        stop 'Invalid MPI interface ibool_interfaces'
      endif
    enddo
  endif

  ! MPI communications
  if (ACOUSTIC_SIMULATION) then
    allocate(buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh*NB_RUNS_ACOUSTIC_GPU,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1560')
    allocate(buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh*NB_RUNS_ACOUSTIC_GPU,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1561')
    buffer_send_scalar_ext_mesh(:,:) = 0.0_CUSTOM_REAL; buffer_recv_scalar_ext_mesh(:,:) = 0.0_CUSTOM_REAL

    allocate(request_send_scalar_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1562')
    allocate(request_recv_scalar_ext_mesh(num_interfaces_ext_mesh), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1563')
    if (ier /= 0) stop 'Error allocating array buffer_send_scalar_ext_mesh,.. for acoustic simulations'
    request_send_scalar_ext_mesh(:) = 0; request_recv_scalar_ext_mesh(:) = 0
  endif

  if (ELASTIC_SIMULATION) then
    allocate(buffer_send_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1564')
    allocate(buffer_recv_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1565')
    buffer_send_vector_ext_mesh(:,:,:) = 0.0_CUSTOM_REAL; buffer_recv_vector_ext_mesh(:,:,:) = 0.0_CUSTOM_REAL

    allocate(request_send_vector_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1566')
    allocate(request_recv_vector_ext_mesh(num_interfaces_ext_mesh), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1567')
    if (ier /= 0) stop 'Error allocating array buffer_send_vector_ext_mesh,.. for elastic simulations'
    request_send_vector_ext_mesh(:) = 0; request_recv_vector_ext_mesh(:) = 0
  endif

  if (POROELASTIC_SIMULATION) then
    allocate(buffer_send_vector_ext_mesh_s(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1568')
    allocate(buffer_recv_vector_ext_mesh_s(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1569')
    allocate(buffer_send_vector_ext_mesh_w(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1570')
    allocate(buffer_recv_vector_ext_mesh_w(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1571')
    buffer_send_vector_ext_mesh_s(:,:,:) = 0.0_CUSTOM_REAL; buffer_recv_vector_ext_mesh_s(:,:,:) = 0.0_CUSTOM_REAL
    buffer_send_vector_ext_mesh_w(:,:,:) = 0.0_CUSTOM_REAL; buffer_recv_vector_ext_mesh_w(:,:,:) = 0.0_CUSTOM_REAL

    allocate(request_send_vector_ext_mesh_s(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1572')
    allocate(request_recv_vector_ext_mesh_s(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1573')
    allocate(request_send_vector_ext_mesh_w(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1574')
    allocate(request_recv_vector_ext_mesh_w(num_interfaces_ext_mesh), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1575')
    if (ier /= 0) stop 'Error allocating array buffer_send_vector_ext_mesh_s,.. for poroelastic simulations'
    request_send_vector_ext_mesh_s(:) = 0; request_recv_vector_ext_mesh_s(:) = 0
    request_send_vector_ext_mesh_w(:) = 0; request_recv_vector_ext_mesh_w(:) = 0
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  done"
    write(IMAIN,*)
    call flush_IMAIN()
  endif

#else
  ! no HDF5 compilation support
  ! user output
  print *
  print *, "Error: HDF5 routine read_mesh_databases_hdf5() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *
  stop 'Error HDF5 read_mesh_databases_hdf5(): called without compilation support'

#endif

  end subroutine read_mesh_databases_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_moho_hdf5()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none

  !#TODO: HDF5 moho version not implemented yet
  stop 'HDF5_ENABLED version of routine read_mesh_databases_moho_hdf5() not implemented yet'

  end subroutine read_mesh_databases_moho_hdf5

