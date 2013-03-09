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
! United States and French Government Sponsorship Acknowledged.

subroutine pml_output_VTKs()

  ! outputs informations about C-PML elements in VTK-file format

  use pml_par
  use specfem_par, only: NGLOB_AB,NSPEC_AB,myrank,prname,xstore,ystore,zstore,ibool
  use constants, only: NGLLX,NGLLY,NGLLZ,IMAIN

  implicit none

  ! local parameters
  integer :: ispec,ispec_CPML,ier
  integer, dimension(:), allocatable :: temp_CPML_regions
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable:: temp_d_store_x,temp_d_store_y,temp_d_store_z
  character(len=256) :: vtkfilename

  if(myrank == 0) write(IMAIN,*) 'Writing informations about C-PML elements in VTK-file format'

  ! C-PML regions
  allocate(temp_CPML_regions(NSPEC_AB),stat=ier)
  if(ier /= 0) stop 'error allocating array temp_CPML_regions'

  temp_CPML_regions(:) = 0

  do ispec_CPML=1,nspec_cpml
     ispec = CPML_to_spec(ispec_CPML)

     temp_CPML_regions(ispec) = CPML_regions(ispec_CPML)
  enddo

  if(myrank == 0) write(IMAIN,*) 'Generating CPML_regions VTK file'

  vtkfilename = prname(1:len_trim(prname))//'CPML_regions'
  call write_VTK_data_elem_i(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,temp_CPML_regions,vtkfilename)

  deallocate(temp_CPML_regions)

  ! C-PML damping profile arrays
  allocate(temp_d_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if(ier /= 0) stop 'error allocating array temp_d_store_x'
  allocate(temp_d_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if(ier /= 0) stop 'error allocating array temp_d_store_y'
  allocate(temp_d_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if(ier /= 0) stop 'error allocating array temp_d_store_z'

  temp_d_store_x(:,:,:,:) = 0._CUSTOM_REAL
  temp_d_store_y(:,:,:,:) = 0._CUSTOM_REAL
  temp_d_store_z(:,:,:,:) = 0._CUSTOM_REAL

  do ispec_CPML=1,nspec_cpml
     ispec = CPML_to_spec(ispec_CPML)

     temp_d_store_x(:,:,:,ispec) = d_store_x(:,:,:,ispec_CPML)
     temp_d_store_y(:,:,:,ispec) = d_store_y(:,:,:,ispec_CPML)
     temp_d_store_z(:,:,:,ispec) = d_store_z(:,:,:,ispec_CPML)
  enddo

  if(myrank == 0) write(IMAIN,*) 'Generating CPML_damping_dx, CPML_damping_dy and CPML_damping_dz VTK files'

  vtkfilename = prname(1:len_trim(prname))//'CPML_damping_dx'
  call write_VTK_data_gll_cr(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,temp_d_store_x,vtkfilename)

  vtkfilename = prname(1:len_trim(prname))//'CPML_damping_dy'
  call write_VTK_data_gll_cr(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,temp_d_store_y,vtkfilename)

  vtkfilename = prname(1:len_trim(prname))//'CPML_damping_dz'
  call write_VTK_data_gll_cr(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,temp_d_store_z,vtkfilename)

  deallocate(temp_d_store_x)
  deallocate(temp_d_store_y)
  deallocate(temp_d_store_z)

end subroutine pml_output_VTKs
