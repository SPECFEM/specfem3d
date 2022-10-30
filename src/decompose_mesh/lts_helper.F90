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


  subroutine lts_add_overlap_elements(nglob,nspec,iglob_p_refine,ispec_p_refine)

  use decompose_mesh_par, only: NGNOD, elmnts

  implicit none

  integer,intent(in) :: nspec,nglob
  integer,dimension(nglob),intent(inout) :: iglob_p_refine
  integer,dimension(nspec),intent(inout) :: ispec_p_refine

  ! local parameters
  integer :: ispec,iglob,i,p

  ! user output
  print *,'  adding overlap region'
  print *

  ! adds overlap by one element
  ! 1. overlap
  do ispec = 1,NSPEC
    ! sets refinement on all points
    ! this enlarges the finer region by one element
    p = ispec_p_refine(ispec)
    do i = 1,NGNOD
      iglob = elmnts(i,ispec)
      if ( p > iglob_p_refine(iglob) ) iglob_p_refine(iglob) = p
    enddo
  enddo
  ! re-sets p refinement for element due to overlap addition above
  ! (adds also elements touching global points with higher refinement)
  ispec_p_refine(:) = 0
  do ispec = 1,NSPEC
    do i = 1,NGNOD
      iglob = elmnts(i,ispec)
      p = iglob_p_refine(iglob)
      if ( p > ispec_p_refine(ispec)) ispec_p_refine(ispec) = p
    enddo
  enddo
  ! checks that all elements have a valid p-value assigned
  if ( minval(ispec_p_refine(:)) == 0 ) stop 'error ispec_p_refine has zero entry'

  ! 2. overlap
  if ( .false. ) then
    do ispec = 1,NSPEC
      ! sets refinement on all points
      ! this enlarges the finer region by one element
      p = ispec_p_refine(ispec)
      do i = 1,NGNOD
        iglob = elmnts(i,ispec)
        if ( p > iglob_p_refine(iglob) ) iglob_p_refine(iglob) = p
      enddo
    enddo
    ! re-sets p refinement for element due to overlap addition above
    ! (adds also elements touching global points with higher refinement)
    ispec_p_refine(:) = 0
    do ispec = 1,NSPEC
      do i = 1,NGNOD
        iglob = elmnts(i,ispec)
        p = iglob_p_refine(iglob)
        if ( p > ispec_p_refine(ispec)) ispec_p_refine(ispec) = p
      enddo
    enddo
    ! checks that all elements have a valid p-value assigned
    if ( minval(ispec_p_refine(:)) == 0 ) stop 'error ispec_p_refine has zero entry'
  endif

  end subroutine lts_add_overlap_elements

!
!------------------------------------------------------------------------------------
!

  subroutine lts_save_p_level_partitions()

  use constants, only: MAX_STRING_LEN

  use decompose_mesh_par, only: part, ispec_p_refine, nspec, nnodes, nodes_coords, elmnts, NGNOD, LOCAL_PATH

  implicit none

  ! local parameters
  ! global coordinates
  double precision, dimension(:),allocatable :: xstore_dummy,ystore_dummy,zstore_dummy
  ! element flag arrays
  integer, dimension(:,:),allocatable :: tmp_elmnts
  integer, dimension(:),allocatable :: tmp_elem_flag
  integer :: ispec,ipart,inum
  integer :: num_parts,nglob
  character(len=MAX_STRING_LEN) :: filename

  nglob = nnodes
  num_parts = maxval(part)+1

  print *,"  Saving p-values for ", num_parts, "partitions"

  ! dummy arrays for vtk output
  allocate(tmp_elmnts(NGNOD,nspec))
  tmp_elmnts(:,:) = 0

  allocate(tmp_elem_flag(nspec))
  tmp_elem_flag(:) = 0

  allocate(xstore_dummy(nglob),ystore_dummy(nglob),zstore_dummy(nglob))
  xstore_dummy(:) = nodes_coords(1,:)
  ystore_dummy(:) = nodes_coords(2,:)
  zstore_dummy(:) = nodes_coords(3,:)

  do ipart = 0,num_parts-1
    ! counts number of element in this partition and sets tmp values for this partition
    inum = 0
    tmp_elmnts(:,:) = 0
    tmp_elem_flag(:) = 0
    do ispec = 1,NSPEC
      if (part(ispec) == ipart) then
        inum = inum + 1
        tmp_elmnts(:,inum) = elmnts(:,ispec)+1
        tmp_elem_flag(inum) = ispec_p_refine(ispec)
      endif
    enddo

    ! filename for current partition
    write(filename,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',ipart,'_lts_p_value'

    ! saves as vtk file
    call write_VTK_data_ngnod_elem_i(inum,nglob,NGNOD,xstore_dummy,ystore_dummy,zstore_dummy, &
                                     tmp_elmnts,tmp_elem_flag,filename)

  enddo
  print *,'  written file(s): ',trim(LOCAL_PATH)//'/proc******_lts_p_value'//'.vtk'
  print *

  deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
  deallocate(tmp_elmnts,tmp_elem_flag)

  end subroutine lts_save_p_level_partitions

!
!-----------------------------------------------------------------------------
!

  subroutine lts_communication_load(ispec_p_refine, nspec)

  use decompose_mesh_par, only: xadj, adjncy, part, p_level, num_p_level

  implicit none

  integer :: nspec
  integer :: ispec, ier,iconn, load_amt
  integer :: num_neighbors, neighbor

  integer, dimension(nspec) :: ispec_p_refine
  integer :: total_communication_load
  integer, dimension(:), allocatable :: communication_load_level

  integer :: num_p,p,ilevel
  integer :: edgecut
  real :: comm_ratio

  ! user output
  print *, 'communication distribution (LTS):'

  num_p = maxval(ispec_p_refine)

  total_communication_load = 0
  edgecut = 0

  allocate(communication_load_level(num_p),stat=ier)
  if (ier /= 0) stop 'Error allocating array communication_load_level'

  communication_load_level(:) = 0

!!!BORA: update so that different levels and different procs get diffent count
  do ispec = 1,nspec
    num_neighbors = xadj(ispec+1)-xadj(ispec)
    do iconn = 1,num_neighbors
      ! get neighbor (+1 for C->Fortran)
      neighbor = adjncy(xadj(ispec) + iconn)+1

      ! neighbor is *not* in our partition
      if (part(ispec) /= part(neighbor)) then
        edgecut = edgecut + 1
        ! neighbor same p-level as ispec
        load_amt = 0
        if (ispec_p_refine(ispec) == ispec_p_refine(neighbor)) then
          load_amt = ispec_p_refine(ispec)
        else if (ispec_p_refine(ispec) > ispec_p_refine(neighbor)) then
          ! communication only required at coarser updates
          load_amt = ispec_p_refine(neighbor)
        else if (ispec_p_refine(ispec) < ispec_p_refine(neighbor)) then
          ! communication only required at coarser updates
          load_amt = ispec_p_refine(ispec)
        endif

        total_communication_load = total_communication_load + load_amt

        communication_load_level(ispec_p_refine(ispec)) = communication_load_level(ispec_p_refine(ispec))&
                                                          + load_amt
      endif
    enddo
  enddo
  if (edgecut > 0) then
    ! imbalance in percent
    comm_ratio = real(total_communication_load)/edgecut
  else
    comm_ratio = 1.0
  endif

  ! user output
  do ilevel = 1,num_p_level
    p = p_level(ilevel)
    print *,'  p-level: ',ilevel,' p = ',p,' communication load = ',communication_load_level(p)
  enddo
  print *,'  total communication load: ', total_communication_load
  print *,'               vs. edgecut: ', edgecut
  print *,'                          : ratio = ',comm_ratio
  print *,'                            ( > 1  being p-level refinement driven, 1 being default)'
  print *

  end subroutine lts_communication_load
