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
!


  subroutine setup_mesh

! mesh creation for solver

  use generate_databases_par

  implicit none

  ! local parameters
  integer :: iface,icorner,inode,ier

  ! compute maximum number of points
  npointot = NSPEC_AB * NGLLX * NGLLY * NGLLZ

  ! use dynamic allocation to allocate memory for arrays
  ! local to global indices array
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 605')
  if (ier /= 0) stop 'error allocating array ibool'
  ibool(:,:,:,:) = 0

  ! node coordinates defined on local level
  allocate(xstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 606')
  if (ier /= 0) stop 'error allocating array xstore'
  xstore(:,:,:,:) = 0.d0
  allocate(ystore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 607')
  if (ier /= 0) stop 'error allocating array ystore'
  ystore(:,:,:,:) = 0.d0
  allocate(zstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 608')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  zstore(:,:,:,:) = 0.d0

  ! estimates memory requirement
  call memory_eval_mesher(NSPEC_AB,npointot,nnodes_ext_mesh, &
                          nelmnts_ext_mesh,nmat_ext_mesh,num_interfaces_ext_mesh, &
                          max_interface_size_ext_mesh,nspec2D_xmin,nspec2D_xmax, &
                          nspec2D_ymin,nspec2D_ymax,nspec2D_bottom,nspec2D_top, &
                          max_memory_size_request)

  max_memory_size = max_memory_size_request

  ! make sure everybody is synchronized
  call synchronize_all()

  ! main working routine to create all the regions of the mesh
  if (myrank == 0) then
    write(IMAIN,*) 'create regions:'
  endif
  call create_regions_mesh()

  ! print min and max of topography included
  min_elevation = HUGEVAL
  max_elevation = -HUGEVAL
  do iface = 1,nspec2D_top_ext
     do icorner = 1,NGNOD2D
        inode = nodes_ibelm_top(icorner,iface)
        if (nodes_coords_ext_mesh(3,inode) < min_elevation) then
           min_elevation = nodes_coords_ext_mesh(3,inode)
        endif
        if (nodes_coords_ext_mesh(3,inode) > max_elevation) then
           max_elevation = nodes_coords_ext_mesh(3,inode)
        endif
     enddo
  enddo

  ! compute the maximum of the maxima for all the slices using an MPI reduction
  call min_all_dp(min_elevation,min_elevation_all)
  call max_all_dp(max_elevation,max_elevation_all)

  if (myrank == 0) then
    if (min_elevation_all /= HUGEVAL .and. max_elevation_all /= -HUGEVAL) then
      write(IMAIN,*)
      write(IMAIN,*) 'min and max of elevation (i.e. height of the upper surface of the mesh) included in mesh in m is ', &
                           min_elevation_all,' ',max_elevation_all
      write(IMAIN,*)
    else
      write(IMAIN,*)
      write(IMAIN,*) 'no upper surface of the mesh detected (no "topography" included in the mesh), there is something wrong'
      call exit_MPI(myrank,'wrong or empty definition of upper surface of the mesh in file free_or_absorbing_surface_file_zmax')
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! make sure everybody is synchronized
  call synchronize_all()

  ! clean-up
  deallocate(xstore,ystore,zstore)
  deallocate(ibool)
  deallocate(ispec_is_surface_external_mesh)
  deallocate(iglob_is_surface_external_mesh)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'done mesh setup'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  end subroutine setup_mesh
