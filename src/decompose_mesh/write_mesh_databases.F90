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


  subroutine write_mesh_databases()

! writes out new Databases files for each partition

  use decompose_mesh_par

  implicit none

  ! local parameters
  integer :: ier
  integer :: my_ninterface,max_interface_size

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

  ! writes out Database file for each partition
  do ipart = 0, nparts-1

    ! opens output file
    write(prname, "(i6.6,'_Database')") ipart

    open(unit=IIN_database,file=outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname, &
         status='unknown', action='write', form='unformatted', iostat = ier)
    if (ier /= 0) then
      print *,'Error file open:',outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname
      print *
      print *,'check if path exists:',outputpath_name(1:len_trim(outputpath_name))
      stop 'Error file open Database'
    endif

    ! gets number of nodes
    call write_glob2loc_nodes_database(IIN_database, ipart, nnodes_loc, nodes_coords, &
                                       glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                       glob2loc_nodes, nnodes, 1)

    call write_partition_database(IIN_database, ipart, nspec_local, nspec, elmnts, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, mat, NGNOD, 1)

    !debug
    !print *, ipart,": nspec_local=",nspec_local, " nnodes_local=", nnodes_loc

    ! writes out node coordinate locations
    write(IIN_database) nnodes_loc

    call write_glob2loc_nodes_database(IIN_database, ipart, nnodes_loc, nodes_coords, &
                                       glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                       glob2loc_nodes, nnodes, 2)

    call write_material_props_database(IIN_database,count_def_mat,count_undef_mat, &
                                       mat_prop, undef_mat_prop)

    ! writes out spectral element indices
    write(IIN_database) nspec_local

    call write_partition_database(IIN_database, ipart, nspec_local, nspec, elmnts, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, mat, NGNOD, 2)

    ! writes out absorbing/free-surface boundaries
    call write_boundaries_database(IIN_database, ipart, nspec, nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, &
                                  nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                                  ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                                  ibelm_ymax, ibelm_bottom, ibelm_top, &
                                  nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                                  nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, &
                                  glob2loc_nodes_parts, glob2loc_nodes, part, NGNOD2D)

    ! writes out C-PML elements indices, CPML-regions and thickness of C-PML layer
    call write_cpml_database(IIN_database, ipart, nspec, nspec_cpml, CPML_to_spec, &
                             CPML_regions, is_CPML, glob2loc_elmnts, part)

    ! gets number of MPI interfaces
    call write_interfaces_database(IIN_database, tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                   my_ninterface, my_interfaces, my_nb_interfaces, &
                                   glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                   glob2loc_nodes, 1, nparts)

    ! writes out MPI interfaces elements
    if (my_ninterface == 0) then
      max_interface_size = 0
    else
      max_interface_size = maxval(my_nb_interfaces)
    endif
    ! user output
    print *,'  partition ',ipart,'has number of MPI interfaces: ',my_ninterface,'maximum size',max_interface_size

    ! format: #number_of_MPI_interfaces  #maximum_number_of_elements_on_each_interface
    write(IIN_database) my_ninterface, max_interface_size

    call write_interfaces_database(IIN_database, tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                   my_ninterface, my_interfaces, my_nb_interfaces, &
                                   glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                   glob2loc_nodes, 2, nparts)

    ! writes out moho surface (optional)
    call write_moho_surface_database(IIN_database, ipart, nspec, &
                                     glob2loc_elmnts, glob2loc_nodes_nparts, &
                                     glob2loc_nodes_parts, glob2loc_nodes, part, &
                                     nspec2D_moho, ibelm_moho, nodes_ibelm_moho, NGNOD2D)

    close(IIN_database)


     ! write fault database
     if (ANY_FAULT) then
        write(prname, "(i6.6,'_Database_fault')") ipart
        open(unit=16,file=outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname, &
             status='replace', action='write', form='unformatted', iostat = ier)
        if (ier /= 0) then
          print *,'Error file open:',outputpath_name(1:len_trim(outputpath_name))//'/proc'//prname
          print *
          print *,'check if path exists:',outputpath_name(1:len_trim(outputpath_name))
          stop
        endif
        call write_fault_database(16, ipart, nspec, &
                                  glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                  glob2loc_nodes, part)
        !write(16,*) nnodes_loc
        write(16) nnodes_loc

        call write_glob2loc_nodes_database(16, ipart, nnodes_loc, nodes_coords_open, &
                                           glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                           glob2loc_nodes, nnodes, 2)
        close(16)
     endif


  enddo

  ! cleanup
  deallocate(CPML_to_spec,stat=ier); if (ier /= 0) stop 'Error deallocating array CPML_to_spec'
  deallocate(CPML_regions,stat=ier); if (ier /= 0) stop 'Error deallocating array CPML_regions'
  deallocate(is_CPML,stat=ier); if (ier /= 0) stop 'Error deallocating array is_CPML'

  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) close(124)

  ! user output
  print *
  print *, 'Databases files in directory: ',outputpath_name(1:len_trim(outputpath_name))
  print *

  end subroutine write_mesh_databases

