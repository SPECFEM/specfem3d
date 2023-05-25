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


  subroutine write_mesh_databases_hdf5()

! writes out new Databases files for each partition in HDF5-format

#ifdef USE_HDF5
  use decompose_mesh_par
  use fault_scotch, only: ANY_FAULT,nodes_coords_open,write_fault_database

  ! HDF5 file i/o
  use manager_hdf5
  use part_decompose_mesh_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! local parameters
  integer :: ier
  integer :: my_ninterface,max_interface_size

  ! processor dependent group names
  character(len=64) :: dset_name

  ! offset arrays
  integer, dimension(0:nparts-1, 6) :: n_elmet_on_bound_all_part
  integer, dimension(0:nparts-1) :: my_ninterface_all
  integer, dimension(0:nparts-1,0:ninterfaces-1) :: my_interfaces_all, my_nb_interfaces_all
  integer, dimension(0:nparts-1) :: offset_nnodes           ! number of nodes in each partition
  integer, dimension(0:nparts-1) :: offset_nelems           ! number of elements in each partition
  integer, dimension(0:nparts-1) :: offset_nelems_cpml      ! number of cpml elements
  integer, dimension(0:nparts-1) :: offset_n_elms_bounds    ! number of elements contacting to the bounds
  integer, dimension(0:nparts-1) :: offset_n_elms_interface ! number of elements on the interfaces
  integer, dimension(0:nparts-1) :: offset_nb_interfaces    ! number of interfaces

  ! user output
  print *, 'HDF5 Database output'
  print *

  ! initializes
  my_ninterface = 0
  max_interface_size = 0

  ! allocates MPI interfaces
  if (ninterfaces > 0) then
    allocate(my_interfaces(0:ninterfaces-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 115')
    if (ier /= 0) stop 'Error allocating array my_interfaces'
    allocate(my_nb_interfaces(0:ninterfaces-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 116')
    if (ier /= 0) stop 'Error allocating array my_nb_interfaces'
  else
    ! dummy allocation
    allocate(my_interfaces(1),my_nb_interfaces(1))
  endif
  my_interfaces(:) = 0; my_nb_interfaces(:) = 0

  ! user output
  print *, 'partitions: '
  print *, '  num         = ',nparts
  print *, '  ninterfaces = ',ninterfaces
  print *

  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) open(124,file='Numglob2loc_elmn.txt')

  ! opens output file
  prname = "Database.h5"
  name_database_hdf5 = outputpath_name(1:len_trim(outputpath_name))//'/'//prname

  ! initialize hdf5 io
  call h5_initialize()

  ! create a hdf5 file
  call h5_create_file(name_database_hdf5)

  ! continue opening hdf5 file till the end of write process
  call h5_open_file(name_database_hdf5)

  ! prepare number of nodes and  offset arrays
  do ipart = 0, nparts-1
    ! gets number of nodes
     call write_glob2loc_nodes_database_h5(ipart, nodes_coords, &
                                           glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                           glob2loc_nodes, nnodes, 1, offset_nnodes)
    ! gets number of elements
     call write_partition_database_h5(ipart, nspec, elmnts, &
                                      glob2loc_elmnts, glob2loc_nodes_nparts, &
                                      glob2loc_nodes_parts, glob2loc_nodes, part, mat, NGNOD, 1, offset_nelems, offset_nnodes)

    ! gets number absorbing/free-surface boundaries
    call write_boundaries_database_h5(ipart, nspec, &
                                      nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, &
                                      nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                                      ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                                      ibelm_ymax, ibelm_bottom, ibelm_top, &
                                      nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                                      nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
                                      glob2loc_elmnts, glob2loc_nodes_nparts, &
                                      glob2loc_nodes_parts, glob2loc_nodes, part, &
                                      NGNOD2D, n_elmet_on_bound_all_part(ipart,:), 1, &
                                      offset_n_elms_bounds)

    ! count C-PML elements indices, CPML-regions and thickness of C-PML layer
    call write_cpml_database_h5(ipart, nspec, nspec_cpml, CPML_to_spec, &
                                CPML_regions, is_CPML, glob2loc_elmnts, part, offset_nelems_cpml, offset_nelems, 1)

    ! count number of MPI interfaces
    call write_interfaces_database_h5(tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                      my_ninterface_all(ipart), my_interfaces_all(ipart,:), my_nb_interfaces_all(ipart,:), &
                                      glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                      glob2loc_nodes, 1, &
                                      nparts, offset_n_elms_interface, offset_nb_interfaces)
  enddo

  ! offset arrays
  dset_name = "offset_nnodes"
  call h5_write_dataset_no_group(dset_name, offset_nnodes)
  dset_name = "offset_nelems"
  call h5_write_dataset_no_group(dset_name, offset_nelems)
  dset_name = "offset_nelems_cpml"
  call h5_write_dataset_no_group(dset_name, offset_nelems_cpml)
  dset_name = "offset_n_elms_bounds"
  call h5_write_dataset_no_group(dset_name, offset_n_elms_bounds)
  dset_name = "offset_n_elms_interface"
  call h5_write_dataset_no_group(dset_name, offset_n_elms_interface)
  dset_name = "offset_nb_interfaces"
  call h5_write_dataset_no_group(dset_name, offset_nb_interfaces)

  ! database arrays
  dset_name = "nodes_coords"
  call h5_create_dataset_gen(dset_name,(/3,sum(offset_nnodes(:))/), 2, 8)
  dset_name = "glob2loc_nodes"
  call h5_create_dataset_gen(dset_name,(/sum(offset_nnodes(:))/), 1, 1)
  dset_name = "nnodes_loc"
  call h5_create_dataset_gen(dset_name,(/nparts/), 1, 1)
  dset_name = "elm_conn"
  call h5_create_dataset_gen(dset_name,(/NGNOD,sum(offset_nelems(:))/), 2, 1)
  dset_name = "elm_conn_xdmf"
  call h5_create_dataset_gen(dset_name,(/NGNOD+1,sum(offset_nelems(:))/), 2, 1)
  dset_name = "nspec_local"
  call h5_create_dataset_gen(dset_name,(/nparts/), 1, 1)
  dset_name = "mat_mesh"
  call h5_create_dataset_gen(dset_name,(/2,sum(offset_nelems(:))/), 2, 1)
  dset_name = "ispec_local"
  call h5_create_dataset_gen(dset_name,(/sum(offset_nelems(:))/), 1, 1)
  dset_name = "n_elms_on_bound"
  call h5_create_dataset_gen(dset_name,(/6*nparts/), 1, 1)
  dset_name = "glob2loc_elms"
  call h5_create_dataset_gen(dset_name, (/1+NGNOD2D,sum(offset_n_elms_bounds(:))/),2,1)
  dset_name = "elements_cpml"
  call h5_create_dataset_gen(dset_name, (/2,sum(offset_nelems_cpml(:))/),2,1)
  dset_name = "if_cpml"
  call h5_create_dataset_gen(dset_name, (/sum(offset_nelems(:))/),1,1)
  dset_name = "nspec_cpml_globloc"
  call h5_create_dataset_gen(dset_name, (/2*nparts/),1,1)
  dset_name = "my_nb_interfaces"
  call h5_create_dataset_gen(dset_name,(/2, sum(offset_nb_interfaces(:))/), 2, 1)
  dset_name = "my_interfaces"
  call h5_create_dataset_gen(dset_name,(/6, sum(offset_n_elms_interface(:))/), 2, 1)
  dset_name = "my_ninterface_and_max"
  call h5_create_dataset_gen(dset_name,(/nparts*2/), 1, 1)

  ! writes out Database file for each partition
  do ipart = 0, nparts-1
     ! writes out node coordinate locations
     call write_glob2loc_nodes_database_h5(ipart, nodes_coords, &
                                           glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                           glob2loc_nodes, nnodes, 2, offset_nnodes)

     ! writes out spectral element indices
     !write(IIN_database) nspec_local
     call write_partition_database_h5(ipart, nspec, elmnts, &
                                      glob2loc_elmnts, glob2loc_nodes_nparts, &
                                      glob2loc_nodes_parts, glob2loc_nodes, part, mat, NGNOD, 2, offset_nelems, offset_nnodes)

     ! writes out absorbing/free-surface boundaries
     call write_boundaries_database_h5(ipart, nspec, &
                                       nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, &
                                       nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                                       ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                                       ibelm_ymax, ibelm_bottom, ibelm_top, &
                                       nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                                       nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
                                       glob2loc_elmnts, glob2loc_nodes_nparts, &
                                       glob2loc_nodes_parts, glob2loc_nodes, part, &
                                       NGNOD2D, n_elmet_on_bound_all_part(ipart,:), 2, &
                                       offset_n_elms_bounds)

     ! writes out C-PML elements indices, CPML-regions and thickness of C-PML layer
     call write_cpml_database_h5(ipart, nspec, nspec_cpml, CPML_to_spec, &
                                 CPML_regions, is_CPML, glob2loc_elmnts, part, offset_nelems_cpml, offset_nelems, 2)

     ! writes out MPI interfaces elements
     call write_interfaces_database_h5(tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                       my_ninterface_all(ipart), my_interfaces_all(ipart,:), my_nb_interfaces_all(ipart,:), &
                                       glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                       glob2loc_nodes, 2, &
                                       nparts, offset_n_elms_interface, offset_nb_interfaces)

     ! writes out moho surface (optional)
     call write_moho_surface_database_h5(ipart, nspec, &
                                         glob2loc_elmnts, glob2loc_nodes_nparts, &
                                         glob2loc_nodes_parts, glob2loc_nodes, part, &
                                         nspec2D_moho, ibelm_moho, nodes_ibelm_moho, NGNOD2D)
  enddo ! end iproc loop

  ! write material properties
  call write_material_props_database_h5(count_def_mat, count_undef_mat, &
                                        mat_prop, undef_mat_prop)

  call h5_close_file()
  call h5_finalize()

  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) close(124)

  ! #TODO: HDF5 support for fault simulation not implemented yet
  ! writes fault database
  if (ANY_FAULT) then
    do ipart = 0, nparts-1
      ! opens file
      write(prname, "(i6.6,'_Database_fault')") ipart

      open(unit=16,file=outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname, &
           status='replace', action='write', form='unformatted', iostat = ier)
      if (ier /= 0) then
        print *,'Error file open:',outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname
        print *
        print *,'check if path exists:',outputpath_name(1:len_trim(outputpath_name))
        stop
      endif

      ! writes out fault element indices
      call write_fault_database(16, ipart, nspec, &
                                glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                glob2loc_nodes, part)

      ! gets number of local nodes nnodes_loc in this partition
      call write_glob2loc_nodes_database(16, ipart, nnodes_loc, nodes_coords_open, &
                                         glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                         glob2loc_nodes, nnodes, 1)

      ! ascii file format
      !write(16,*) nnodes_loc
      ! binary file format
      write(16) nnodes_loc

      ! writes out coordinates nodes_coords_open of local nodes
      ! note that we output here the "open" split node positions for the fault nodal points
      call write_glob2loc_nodes_database(16, ipart, nnodes_loc, nodes_coords_open, &
                                         glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                         glob2loc_nodes, nnodes, 2)
      close(16)
    enddo
  endif

  ! cleanup
  deallocate(CPML_to_spec,stat=ier); if (ier /= 0) stop 'Error deallocating array CPML_to_spec'
  deallocate(CPML_regions,stat=ier); if (ier /= 0) stop 'Error deallocating array CPML_regions'
  deallocate(is_CPML,stat=ier); if (ier /= 0) stop 'Error deallocating array is_CPML'

  ! user output
  print *
  print *, 'Databases files in directory: ',outputpath_name(1:len_trim(outputpath_name))
  print *

#else
  ! no HDF5 compilation support
  ! user output
  print *
  print *, "Error: HDF5 routine write_mesh_databases_hdf5() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *
  stop 'Error HDF5 write_mesh_databases_hdf5(): called without compilation support'

#endif

  end subroutine write_mesh_databases_hdf5

