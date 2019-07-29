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


  subroutine write_mesh_databases_hdf5()

! writes out new Databases files for each partition

  use hdf5

  use decompose_mesh_par

  implicit none

  ! local parameters
  integer :: ier
  integer(HID_T) :: file_id
  
  ! processor dependent group names
  character(len=5)  :: gname_proc_head = "proc_"
  character(len=11) :: gname_proc
  character(len=10) :: tempstr

  allocate(my_interfaces(0:ninterfaces-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 115')
  if (ier /= 0) stop 'Error allocating array my_interfaces'
  allocate(my_nb_interfaces(0:ninterfaces-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 116')
  if (ier /= 0) stop 'Error allocating array my_nb_interfaces'

  ! user output
  print *, 'partitions: '
  print *, '  num = ',nparts
  print *

  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) open(124,file='Numglob2loc_elmn.txt')

  ! opens output file
  prname = "Database.h5"
  IIN_database_hdf5 = outputpath_name(1:len_trim(outputpath_name))//prname

  ! initialize hdf5 file
  call h5open_f(ier)
  call hdf5_initialize(IIN_database_hdf5, nparts)

  ! continue opening hdf5 file till the end of write process
  call h5fopen_f(IIN_database_hdf5, H5F_ACC_RDWR_F, file_id, ier) ! file open !!!! error
 
  ! writes out Database file for each partition
  do ipart = 0, nparts-1

     ! set hdf5 group name
     write(tempstr, "(i6.6)") ipart
     gname_proc = gname_proc_head // trim(tempstr)
 
     ! gets number of nodes
     call write_glob2loc_nodes_database_h5(IIN_database_hdf5, gname_proc, file_id, ipart, nnodes_loc, nodes_coords, &
                                glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                glob2loc_nodes, nnodes, 1)

     call write_partition_database_h5(IIN_database_hdf5, gname_proc, file_id, ipart, nspec_local, nspec, elmnts, &
                                glob2loc_elmnts, glob2loc_nodes_nparts, &
                                glob2loc_nodes_parts, glob2loc_nodes, part, mat, NGNOD, 1)

     !debug
     !print *, ipart,": nspec_local=",nspec_local, " nnodes_local=", nnodes_loc

     ! writes out node coordinate locations
     ! write(IIN_database) nnodes_loc !!!!!!!!!!!!! nnodes_loc will be containded as an attribute of glob2loc_nodes group

     call write_glob2loc_nodes_database_h5(IIN_database_hdf5, gname_proc, file_id, ipart, nnodes_loc, nodes_coords, &
                                glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                glob2loc_nodes, nnodes, 2)

     ! writes out spectral element indices
     !write(IIN_database) nspec_local
     call write_partition_database_h5(IIN_database_hdf5, gname_proc, file_id, ipart, nspec_local, nspec, elmnts, &
                                glob2loc_elmnts, glob2loc_nodes_nparts, &
                                glob2loc_nodes_parts, glob2loc_nodes, part, mat, NGNOD, 2)

     ! writes out absorbing/free-surface boundaries
     call write_boundaries_database_h5(IIN_database_hdf5, gname_proc, file_id, ipart, nspec, &
                                nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, &
                                nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                                ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                                ibelm_ymax, ibelm_bottom, ibelm_top, &
                                nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                                nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
                                glob2loc_elmnts, glob2loc_nodes_nparts, &
                                glob2loc_nodes_parts, glob2loc_nodes, part, NGNOD2D)

     ! writes out C-PML elements indices, CPML-regions and thickness of C-PML layer
     call write_cpml_database_h5(IIN_database_hdf5, gname_proc, file_id, ipart, nspec, nspec_cpml, CPML_to_spec, &
          CPML_regions, is_CPML, glob2loc_elmnts, part)

     ! gets number of MPI interfaces
     call write_interfaces_database_h5(IIN_database_hdf5, gname_proc, file_id, &
                                tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                my_ninterface, my_interfaces, my_nb_interfaces, &
                                glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                glob2loc_nodes, 1, nparts)

     ! writes out MPI interfaces elements
     !print *,' my interfaces:',my_ninterface,maxval(my_nb_interfaces)
     !if (my_ninterface == 0) then
     ! write(IIN_database) my_ninterface, 0       ! avoids problem with maxval for empty array my_nb_interfaces
     !else
     ! write(IIN_database) my_ninterface, maxval(my_nb_interfaces)
     !endif

     call write_interfaces_database_h5(IIN_database_hdf5, gname_proc, file_id, &
                                tab_interfaces, tab_size_interfaces, ipart, ninterfaces, &
                                my_ninterface, my_interfaces, my_nb_interfaces, &
                                glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                glob2loc_nodes, 2, nparts)

     ! writes out moho surface (optional)
     call write_moho_surface_database_h5(IIN_database_hdf5, gname_proc, file_id, ipart, nspec, &
                                glob2loc_elmnts, glob2loc_nodes_nparts, &
                                glob2loc_nodes_parts, glob2loc_nodes, part, &
                                nspec2D_moho, ibelm_moho, nodes_ibelm_moho, NGNOD2D)

     !close(IIN_database)


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

  enddo ! end iproc loop

  ! write material properties
  call write_material_props_database_h5(IIN_database_hdf5,file_id,count_def_mat,count_undef_mat, &
                                     mat_prop, undef_mat_prop)

  call h5fclose_f(file_id, ier)
  call h5close_f(ier)

  ! cleanup
  deallocate(CPML_to_spec,stat=ier); if (ier /= 0) stop 'Error deallocating array CPML_to_spec'
  deallocate(CPML_regions,stat=ier); if (ier /= 0) stop 'Error deallocating array CPML_regions'
  deallocate(is_CPML,stat=ier); if (ier /= 0) stop 'Error deallocating array is_CPML'

  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) close(124)

  ! user output
  print *, 'Databases files in directory: ',outputpath_name(1:len_trim(outputpath_name))

  end subroutine write_mesh_databases_hdf5

