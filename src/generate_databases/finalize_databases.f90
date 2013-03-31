!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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


  subroutine finalize_databases

! checks user input parameters

  use generate_databases_par
  implicit none

  integer :: i

! print number of points and elements in the mesh
  call sum_all_i(NGLOB_AB,nglob_total)
  call sum_all_i(NSPEC_AB,nspec_total)
  call sync_all()
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Repartition of elements:'
    write(IMAIN,*) '-----------------------'
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in each slice: ',NSPEC_AB
    write(IMAIN,*) 'total number of points in each slice: ',NGLOB_AB
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in entire mesh: ',nspec_total
    write(IMAIN,*) 'total number of points in entire mesh: ',nglob_total
    write(IMAIN,*) 'total number of DOFs in entire mesh: ',nglob_total*NDIM
    write(IMAIN,*)
    write(IMAIN,*) 'total number of time steps in the solver will be: ',NSTEP
    write(IMAIN,*)
    ! write information about precision used for floating-point operations
    if(CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    write(IMAIN,*)
  endif

! gets number of surface elements (for movie outputs)
  allocate( ispec_is_surface_external_mesh(NSPEC_AB), &
           iglob_is_surface_external_mesh(NGLOB_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating array'
  max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh)
  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array'
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,:) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,:)
  enddo
  call sync_all()
  call detect_surface(NPROC,NGLOB_AB,NSPEC_AB,ibool, &
                        ispec_is_surface_external_mesh, &
                        iglob_is_surface_external_mesh, &
                        nfaces_surface_ext_mesh, &
                        num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh, &
                        ibool_interfaces_ext_mesh_dummy )

  deallocate(ibool)
  deallocate(ispec_is_surface_external_mesh)
  deallocate(iglob_is_surface_external_mesh)
  deallocate(ibool_interfaces_ext_mesh_dummy)

  ! takes number of faces for top, free surface only
  if( MOVIE_TYPE == 1 ) then
    nfaces_surface_ext_mesh = NSPEC2D_TOP
  endif

! number of surface faces for all partitions together
  call sum_all_i(nfaces_surface_ext_mesh,nfaces_surface_glob_ext_mesh)

! copy number of elements and points in an include file for the solver
  if( myrank == 0 ) then

    call save_header_file(NSPEC_AB,NGLOB_AB,NPROC, &
               ATTENUATION,ANISOTROPY,NSTEP,DT,STACEY_INSTEAD_OF_FREE_SURFACE, &
               SIMULATION_TYPE,max_memory_size,nfaces_surface_glob_ext_mesh)
  endif

! elapsed time since beginning of mesh generation
  if(myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for mesh generation and buffer creation in seconds = ',tCPU
    write(IMAIN,*) 'End of mesh generation'
    write(IMAIN,*)
  endif

! close main output file
  if(myrank == 0) then
    write(IMAIN,*) 'done'
    write(IMAIN,*)
    close(IMAIN)
  endif

! synchronize all the processes to make sure everybody has finished
  call sync_all()

  end subroutine finalize_databases

