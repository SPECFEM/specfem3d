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

  subroutine create_CPML_regions(nspec,nglob,nodes_coords)

  use meshfem3D_par, only: myrank,ibool, &
    nspec_CPML,is_CPML,CPML_to_spec,CPML_regions, &
    THICKNESS_OF_X_PML,THICKNESS_OF_Y_PML,THICKNESS_OF_Z_PML

  ! create the different regions of the mesh
  use constants, only: IMAIN,CUSTOM_REAL,SMALL_PERCENTAGE_TOLERANCE, &
                      CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  ! CPML
  use shared_parameters, only: PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE

  implicit none

  integer,intent(in):: nspec,nglob
  double precision, dimension(nglob,3), intent(in) :: nodes_coords

  ! local parameters
  real(kind=CUSTOM_REAL) :: xmin,xmax,ymin,ymax,zmin,zmax,limit
  real(kind=CUSTOM_REAL) :: xmin_all,xmax_all,ymin_all,ymax_all,zmin_all,zmax_all

  logical, dimension(:), allocatable :: is_X_CPML,is_Y_CPML,is_Z_CPML

  integer :: nspec_CPML_total
  integer :: i1,i2,i3,i4,i5,i6,i7,i8
  integer :: ispec,ispec_CPML
  integer :: ier

  ! CPML allocation
  allocate(is_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1310')
  if (ier /= 0) stop 'Error allocating is_CPML array'

  ! initializes CPML elements
  is_CPML(:) = .false.
  nspec_CPML = 0

  ! checks if anything to do
  if (.not. PML_CONDITIONS) then
    ! dummy allocation
    allocate(CPML_to_spec(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1311')
    allocate(CPML_regions(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1312')
    if (ier /= 0) stop 'Error allocating dummy CPML arrays'

    ! nothing to do anymore
    return
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'creating PML region'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! compute the min and max values of each coordinate
  xmin = minval(nodes_coords(:,1))
  xmax = maxval(nodes_coords(:,1))

  ymin = minval(nodes_coords(:,2))
  ymax = maxval(nodes_coords(:,2))

  zmin = minval(nodes_coords(:,3))
  zmax = maxval(nodes_coords(:,3))

  call min_all_all_cr(xmin,xmin_all)
  call min_all_all_cr(ymin,ymin_all)
  call min_all_all_cr(zmin,zmin_all)

  call max_all_all_cr(xmax,xmax_all)
  call max_all_all_cr(ymax,ymax_all)
  call max_all_all_cr(zmax,zmax_all)

  allocate(is_X_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1313')
  allocate(is_Y_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1314')
  allocate(is_Z_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1315')

  is_X_CPML(:) = .false.
  is_Y_CPML(:) = .false.
  is_Z_CPML(:) = .false.

  do ispec = 1,nspec

    i1=ibool(1,1,1,ispec)
    i2=ibool(2,1,1,ispec)
    i3=ibool(2,2,1,ispec)
    i4=ibool(1,2,1,ispec)
    i5=ibool(1,1,2,ispec)
    i6=ibool(2,1,2,ispec)
    i7=ibool(2,2,2,ispec)
    i8=ibool(1,2,2,ispec)

    ! Xmin CPML
    limit=xmin_all+THICKNESS_OF_X_PML * SMALL_PERCENTAGE_TOLERANCE
    if (  nodes_coords(i1,1) < limit .and. nodes_coords(i2,1) < limit .and. &
         nodes_coords(i3,1) < limit .and. nodes_coords(i4,1) < limit .and. &
         nodes_coords(i5,1) < limit .and. nodes_coords(i6,1) < limit .and. &
         nodes_coords(i7,1) < limit .and. nodes_coords(i8,1) < limit) then
       is_X_CPML(ispec) = .true.
    endif

    ! Xmax CPML
    limit = xmax_all - THICKNESS_OF_X_PML * SMALL_PERCENTAGE_TOLERANCE
    if (  nodes_coords(i1,1) > limit .and. nodes_coords(i2,1) > limit .and. &
         nodes_coords(i3,1) > limit .and. nodes_coords(i4,1) > limit .and. &
         nodes_coords(i5,1) > limit .and. nodes_coords(i6,1) > limit .and. &
         nodes_coords(i7,1) > limit .and. nodes_coords(i8,1) > limit) then
       is_X_CPML(ispec) = .true.
    endif

    ! Ymin CPML
    limit=ymin_all+THICKNESS_OF_Y_PML * SMALL_PERCENTAGE_TOLERANCE
    if (  nodes_coords(i1,2) < limit .and. nodes_coords(i2,2) < limit .and. &
         nodes_coords(i3,2) < limit .and. nodes_coords(i4,2) < limit .and. &
         nodes_coords(i5,2) < limit .and. nodes_coords(i6,2) < limit .and. &
         nodes_coords(i7,2) < limit .and. nodes_coords(i8,2) < limit) then
       is_Y_CPML(ispec) = .true.
    endif

    ! Ymax CPML
    limit = ymax_all - THICKNESS_OF_Y_PML * SMALL_PERCENTAGE_TOLERANCE
    if (  nodes_coords(i1,2) > limit .and. nodes_coords(i2,2) > limit .and. &
         nodes_coords(i3,2) > limit .and. nodes_coords(i4,2) > limit .and. &
         nodes_coords(i5,2) > limit .and. nodes_coords(i6,2) > limit .and. &
         nodes_coords(i7,2) > limit .and. nodes_coords(i8,2) > limit) then
       is_Y_CPML(ispec) = .true.
    endif

    ! Zmin CPML
    limit=zmin_all+THICKNESS_OF_Z_PML * SMALL_PERCENTAGE_TOLERANCE
    if (  nodes_coords(i1,3) < limit .and. nodes_coords(i2,3) < limit .and. &
         nodes_coords(i3,3) < limit .and. nodes_coords(i4,3) < limit .and. &
         nodes_coords(i5,3) < limit .and. nodes_coords(i6,3) < limit .and. &
         nodes_coords(i7,3) < limit .and. nodes_coords(i8,3) < limit) then
       is_Z_CPML(ispec) = .true.
    endif

    if (PML_INSTEAD_OF_FREE_SURFACE) then
       ! Zmax CPML
       limit = zmax_all - THICKNESS_OF_Z_PML * SMALL_PERCENTAGE_TOLERANCE
       if (  nodes_coords(i1,3) > limit .and. nodes_coords(i2,3) > limit .and. &
            nodes_coords(i3,3) > limit .and. nodes_coords(i4,3) > limit .and. &
            nodes_coords(i5,3) > limit .and. nodes_coords(i6,3) > limit .and. &
            nodes_coords(i7,3) > limit .and. nodes_coords(i8,3) > limit) then
          is_Z_CPML(ispec) = .true.
       endif
    endif

    if (is_X_CPML(ispec) .or. is_Y_CPML(ispec) .or. is_Z_CPML(ispec)) &
         nspec_CPML = nspec_CPML + 1

  enddo

  ! outputs total number of CPML elements
  nspec_CPML_total = 0
  call sum_all_i(nspec_CPML,nspec_CPML_total)
  call bcast_all_singlei(nspec_CPML_total)

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Created a total of ',nspec_CPML_total,' unique CPML elements'
    write(IMAIN,*)'   (i.e., ',100.*nspec_CPML/real(nspec),'% of the mesh)'
    write(IMAIN,*)
  endif

  ! allocates arrays
  allocate(CPML_to_spec(nspec_CPML),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1316')
  allocate(CPML_regions(nspec_CPML),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1317')
  if (ier /= 0) stop 'Error allocating CPML arrays'

  ispec_CPML=0
  do ispec=1,nspec
    if (is_X_CPML(ispec) .and. is_Y_CPML(ispec) .and. is_Z_CPML(ispec)) then
       ispec_CPML=ispec_CPML+1
       CPML_to_spec(ispec_CPML)=ispec
       CPML_regions(ispec_CPML)=CPML_XYZ
       is_CPML(ispec)=.true.
    else if (is_Y_CPML(ispec) .and. is_Z_CPML(ispec)) then
       ispec_CPML=ispec_CPML+1
       CPML_to_spec(ispec_CPML)=ispec
       CPML_regions(ispec_CPML)=CPML_YZ_ONLY
       is_CPML(ispec)=.true.

    else if (is_X_CPML(ispec) .and. is_Z_CPML(ispec)) then
       ispec_CPML=ispec_CPML+1
       CPML_to_spec(ispec_CPML)=ispec
       CPML_regions(ispec_CPML)=CPML_XZ_ONLY
       is_CPML(ispec)=.true.

    else if (is_X_CPML(ispec) .and. is_Y_CPML(ispec)) then
       ispec_CPML=ispec_CPML+1
       CPML_to_spec(ispec_CPML)=ispec
       CPML_regions(ispec_CPML)=CPML_XY_ONLY
       is_CPML(ispec)=.true.

    else if (is_Z_CPML(ispec)) then
       ispec_CPML=ispec_CPML+1
       CPML_to_spec(ispec_CPML)=ispec
       CPML_regions(ispec_CPML)=CPML_Z_ONLY
       is_CPML(ispec)=.true.

    else if (is_Y_CPML(ispec)) then
       ispec_CPML=ispec_CPML+1
       CPML_to_spec(ispec_CPML)=ispec
       CPML_regions(ispec_CPML)=CPML_Y_ONLY
       is_CPML(ispec)=.true.

    else if (is_X_CPML(ispec)) then
       ispec_CPML=ispec_CPML+1
       CPML_to_spec(ispec_CPML)=ispec
       CPML_regions(ispec_CPML)=CPML_X_ONLY
       is_CPML(ispec)=.true.

   endif

  enddo

  ! checks
  if (ispec_CPML /= nspec_CPML) stop 'Error number of CPML element is not consistent'

  end subroutine create_CPML_regions
