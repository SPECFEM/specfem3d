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


  subroutine setup_mesh

! mesh creation for solver

  use generate_databases_par

  implicit none

  ! compute maximum number of points
  npointot = NSPEC_AB * NGLLCUBE

! use dynamic allocation to allocate memory for arrays
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibool'
  allocate(xstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating array xstore'
  allocate(ystore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ystore'
  allocate(zstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  call memory_eval_mesher(myrank,NSPEC_AB,npointot,nnodes_ext_mesh, &
                        nelmnts_ext_mesh,nmat_ext_mesh,num_interfaces_ext_mesh, &
                        max_interface_size_ext_mesh,nspec2D_xmin,nspec2D_xmax, &
                        nspec2D_ymin,nspec2D_ymax,nspec2D_bottom,nspec2D_top, &
                        max_memory_size_request)

  max_memory_size = max_memory_size_request

! make sure everybody is synchronized
  call sync_all()

! main working routine to create all the regions of the mesh
  if(myrank == 0) then
    write(IMAIN,*) 'create regions: '
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

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'min and max of topography included in mesh in m is ',min_elevation_all,' ',max_elevation_all
    write(IMAIN,*)
    call flush_IMAIN()
  endif

! clean-up
  deallocate(xstore,ystore,zstore)

! make sure everybody is synchronized
  call sync_all()

  end subroutine setup_mesh
