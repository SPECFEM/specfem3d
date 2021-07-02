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

module module_partition

  use constants, only: NDIM
  use shared_parameters, only: NGNOD
  use module_qsort

  integer                                                  :: nE
  integer,             dimension(:),  allocatable          :: ipart

contains

  !--------------------------------------------
  ! Heuristic mesh decomposition
  !---------------------------------------------

  !
  !  decomposition is made using a distance criterion from
  !  a given face of a computational box domain taking
  !  into account load of each element.
  !
  subroutine partition_mesh(elmnts, nodes_coords, load_elmnts,  nspec, nnodes, npart_1, npart_2, npart_3)

    implicit none

    ! input mesh
    integer,                                     intent(in) :: nspec, nnodes
    double precision,    dimension(NDIM,nnodes),    intent(in) :: nodes_coords
    double precision,    dimension(nspec),       intent(in) :: load_elmnts
    integer,             dimension(NGNOD,nspec), intent(in) :: elmnts

    ! partition
    integer,                                     intent(in) :: npart_1, npart_2, npart_3

    ! local
    double precision                                        :: xmin, xmax, ymin, ymax, zmin, zmax
    double precision,    dimension(3)                       :: ref_point
    double precision,    dimension(:,:), allocatable        :: elmnts_center
    double precision,    dimension(:,:), allocatable        :: elmnts_center_1, elmnts_center_2, elmnts_center_3
    double precision,    dimension(:),   allocatable        :: cri_load_1, cri_load_2, cri_load_3
    double precision,    dimension(:),   allocatable        :: sum_load_1, sum_load_2, sum_load_3
    double precision,    dimension(:),   allocatable        :: load_elmnts_1, load_elmnts_2, load_elmnts_3
    integer,             dimension(:),   allocatable        :: ipart_1, ipart_2, ipart_3
    integer,             dimension(:),   allocatable        :: iperm_1, iperm_2, iperm_3
    integer,             dimension(:),   allocatable        :: old_num_1, old_num_2, old_num_3
    integer,             dimension(:),   allocatable        :: nEipart_1, nEipart_2, nEipart_3
    integer                                                 :: nE_1, nE_2, nE_3
    integer                                                 :: kpart_2, kpart_3, p1, p2, p3
    integer                                                 :: i, iE, idir, num_original_element, ier
    !

    nE = nspec
    allocate(ipart(nE),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 67')
    ipart(:) = -1

    xmin = minval(nodes_coords(1,:))
    xmax = maxval(nodes_coords(1,:))

    ymin = minval(nodes_coords(2,:))
    ymax = maxval(nodes_coords(2,:))

    zmin = minval(nodes_coords(3,:))
    zmax = maxval(nodes_coords(3,:))

    ref_point(1) = xmin
    ref_point(2) = ymin
    ref_point(3) = zmin
    write(27,*)
    WRITE(27,*) ' xmin, ymin, zmin ', xmin, ymin, zmin
    write(27,*) ' nspec, nnodes ',  nspec, nnodes
    allocate(elmnts_center(3,nE),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 68')
    write(27,*)
    call compute_elmnts_center(elmnts_center, elmnts, nodes_coords, nspec, nnodes)

    ! partition in direction 1 on the whole mesh
    idir = 1
    nE_1 = nE
    allocate(sum_load_1(nE_1),cri_load_1(nE_1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 69')
    allocate(ipart_1(nE_1), nEipart_1(nE_1), iperm_1(nE_1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 70')
    allocate(load_elmnts_1(nE_1), elmnts_center_1(3,nE_1), old_num_1(nE_1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 71')
    elmnts_center_1(:,:)=elmnts_center(:,:)
    load_elmnts_1(:)=load_elmnts(:)
    do i = 1,nE
      old_num_1(i)=i
    enddo

    call compute_partition(ipart_1, nEipart_1, npart_1, sum_load_1, cri_load_1, &
                           load_elmnts_1, elmnts_center_1, iperm_1, nE_1, ref_point, idir)

    write(27,*) ' made partion in 1st direction ', npart_1

    ! partition of all obtained slice in direction 2
    do kpart_2 = 1, npart_1

       idir = 2
       nE_2 = nEipart_1(kpart_2)
       allocate(sum_load_2(nE_2), cri_load_2(nE_2),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 72')
       allocate(load_elmnts_2(nE_2), elmnts_center_2(3,nE_2),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 73')
       allocate(ipart_2(nE_2), nEipart_2(npart_2), iperm_2(nE_2), old_num_2(nE_2),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 74')

       call extract_partition(load_elmnts_2, elmnts_center_2, old_num_2, nE_2, &
            ipart_1, load_elmnts_1, elmnts_center_1, old_num_1, kpart_2, nE_1)

       call compute_partition(ipart_2, nEipart_2, npart_2, sum_load_2, cri_load_2, &
            load_elmnts_2, elmnts_center_2, iperm_2, nE_2, ref_point, idir)

       ! partition of all remained slice in direction 3
       do kpart_3 = 1, npart_2
          idir = 3
          nE_3 = nEipart_2(kpart_3)

          allocate(sum_load_3(nE_3), cri_load_3(nE_3),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 75')
          allocate(load_elmnts_3(nE_3), elmnts_center_3(3,nE_3),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 76')
          allocate(ipart_3(nE_3), nEipart_3(npart_3), iperm_3(nE_3), old_num_3(nE_3),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 77')

          call extract_partition(load_elmnts_3, elmnts_center_3, old_num_3, nE_3, &
               ipart_2, load_elmnts_2, elmnts_center_2, old_num_2, kpart_3, nE_2)

          call compute_partition(ipart_3, nEipart_3, npart_3, sum_load_3, cri_load_3, &
               load_elmnts_3, elmnts_center_3, iperm_3, nE_3, ref_point, idir)

          do iE=1, nE_3
             p1 = kpart_2
             p2 = kpart_3
             p3 = ipart_3(iE)
             num_original_element = old_num_3(iE)
             ipart(num_original_element) = p1 + npart_1*(p2-1) + npart_1*npart_2*(p3-1)
          enddo

          deallocate(load_elmnts_3, elmnts_center_3)
          deallocate(sum_load_3, cri_load_3)
          deallocate(ipart_3,nEipart_3, iperm_3, old_num_3)

       enddo

       deallocate(load_elmnts_2, elmnts_center_2)
       deallocate(sum_load_2, cri_load_2)
       deallocate(ipart_2, nEipart_2, iperm_2, old_num_2)
    enddo

    write(27,*) 'made partition in 2nd and 3th direction', npart_2, npart_3

    deallocate(sum_load_1, cri_load_1)
    deallocate(ipart_1, nEipart_1)
    deallocate(iperm_1, old_num_1)

  end subroutine partition_mesh


  !--------------------------------------------
  ! Geometry mesh decomposition
  !---------------------------------------------

  !
  !  decomposition is made only using a distance criterion from
  !  a given face of a computational box domain
  !

  subroutine partition_mesh_distance(elmnts, nodes_coords, nspec, nnodes, npart_1, npart_2, npart_3)

    implicit none

    ! input mesh
    integer,                                     intent(in) :: nspec, nnodes
    double precision,    dimension(NDIM,nnodes),    intent(in) :: nodes_coords
    integer,             dimension(NGNOD,nspec), intent(in) :: elmnts

    ! partition
    integer,                                     intent(in) :: npart_1, npart_2, npart_3

    ! local
    double precision                                        :: xmin, xmax, ymin,ymax, zmin, zmax
    double precision                                        :: x_bin, y_bin, z_bin
    double precision,    dimension(3)                       :: ref_point
    double precision,    dimension(:,:), allocatable        :: elmnts_center
    integer    :: p1, p2, p3, iE, ier

    nE = nspec
    allocate(ipart(nE),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 67')
    ipart(:) = -1

    xmin = minval(nodes_coords(1,:))
    xmax = maxval(nodes_coords(1,:))
    ymin = minval(nodes_coords(2,:))
    ymax = maxval(nodes_coords(2,:))
    zmin = minval(nodes_coords(3,:))
    zmax = maxval(nodes_coords(3,:))

    ref_point(1) = xmin
    ref_point(2) = ymin
    ref_point(3) = zmin
    write(27,*)
    WRITE(27,*) ' xmin, ymin, zmin ', xmin, ymin, zmin
    write(27,*) ' nspec, nnodes ',  nspec, nnodes
    allocate(elmnts_center(3,nE),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 68')
    write(27,*)
    call compute_elmnts_center(elmnts_center, elmnts, nodes_coords, nspec, nnodes)


    x_bin = (xmax - xmin)/npart_1
    y_bin = (ymax - ymin)/npart_2
    z_bin = (zmax - zmin)/npart_3
    do iE=1, nE
       p1 = floor((elmnts_center(1,iE)-xmin)/x_bin) + 1
       p2 = floor((elmnts_center(2,iE)-ymin)/y_bin) + 1
       p3 = floor((elmnts_center(3,iE)-zmin)/z_bin) + 1
       ipart(iE) = p1 + npart_1*(p2-1) + npart_1*npart_2*(p3-1)
    enddo

    deallocate(elmnts_center)

end subroutine partition_mesh_distance

  !---------
  ! compute partition in one direction
  !---------
  subroutine compute_partition(ipart_tmp, nEipart_tmp, npart_tmp, sum_load_tmp, cri_load_perm, &
                               load_elmnts_tmp, elmnts_center_tmp, iperm_tmp, nE_tmp, ref_point, idir)

    implicit none

    integer,                                  intent(in)    :: nE_tmp, npart_tmp, idir
    double precision,    dimension(3) ,       intent(in)    :: ref_point
    double precision,    dimension(nE_tmp),   intent(in)    :: load_elmnts_tmp
    double precision,    dimension(3,nE_tmp), intent(in)    :: elmnts_center_tmp
    integer,             dimension(npart_tmp),intent(inout) :: nEipart_tmp
    integer,             dimension(nE_tmp),   intent(inout) :: ipart_tmp, iperm_tmp
    double precision,    dimension(nE_tmp),   intent(inout) :: sum_load_tmp, cri_load_perm

    integer                                                 :: i, k
    double precision                                        :: Load_by_part

    ! initialise permutation
    do i = 1,nE_tmp
      iperm_tmp(i) = i
    enddo

    call compute_criteria(cri_load_perm, elmnts_center_tmp, nE_tmp, ref_point, idir)

    call QsortC(cri_load_perm,  iperm_tmp)

    call compute_sum_load(sum_load_tmp, load_elmnts_tmp, iperm_tmp, nE_tmp)

    nEipart_tmp(:)=0
!! DK DK Oct 2018: added this safety test
    if (nE_tmp <= 0) stop 'Error: cannot use an array that has been declared with a size of zero'
    load_by_part = floor(sum_load_tmp(nE_tmp) / real(npart_tmp,8)) + 1

    write(27,*) ' Load value by partition  ', Load_by_part, sum_load_tmp(nE_tmp), real(npart_tmp,8)

    do i = 1,nE_tmp
       k = iperm_tmp(i)
       ipart_tmp(k) = int(floor(sum_load_tmp(i) / Load_by_part)) + 1  !! must have the smaller closest integer + 1
       nEipart_tmp(ipart_tmp(k)) = nEipart_tmp(ipart_tmp(k)) + 1      !! (sum_load_tmp is not sorted)
    enddo

  end subroutine compute_partition

  !----------
  ! extract sub arrays for each partition
  !----------
  subroutine extract_partition(load_elmnts_1, elmnts_center_1, old_num_1, nE_1, &
                               ipart_0, load_elmnts_0, elmnts_center_0, old_num_0, kpart_0, nE_0)

    implicit none

    integer,                                   intent(in)     :: nE_0, nE_1, kpart_0
    integer,              dimension(nE_0),     intent(in)     :: ipart_0, old_num_0
    double precision,     dimension(3,nE_0),   intent(in)     :: elmnts_center_0
    double precision,     dimension(nE_0),     intent(in)     :: load_elmnts_0
    integer,              dimension(nE_1),     intent(inout)  :: old_num_1
    double precision,     dimension(3,nE_1),   intent(inout)  :: elmnts_center_1
    double precision,     dimension(nE_1),     intent(inout)  :: load_elmnts_1

    integer                                                   :: i, k

    k = 0
    do i = 1, nE_0
       if (ipart_0(i) == kpart_0) then
          k = k + 1
          elmnts_center_1(1,k) = elmnts_center_0(1,i)
          elmnts_center_1(2,k) = elmnts_center_0(2,i)
          elmnts_center_1(3,k) = elmnts_center_0(3,i)
          load_elmnts_1(k) = load_elmnts_0(i)
          old_num_1(k) = old_num_0(i)
       endif
    enddo

  end subroutine extract_partition


  !---------
  ! compute criteria for partition : distance from one face
  !---------
  subroutine compute_criteria(cri_load_perm, elmnts_center_tmp, nE_tmp, ref_point, idir)

    implicit none

    integer,                                  intent(in)    :: nE_tmp, idir
    double precision,    dimension(3),        intent(in)    :: ref_point
    double precision,    dimension(3,nE_tmp), intent(in)    :: elmnts_center_tmp
    double precision,    dimension(nE_tmp),   intent(inout) :: cri_load_perm

    integer                                                 :: i

    do i=1, nE_tmp
       cri_load_perm(i) = abs( elmnts_center_tmp(idir,i) - ref_point(idir) )
    enddo

  end subroutine compute_criteria


  !---------
  ! compute sum load of partition
  !---------
  subroutine compute_sum_load(sum_load, load_elmnts, iperm_tmp, nE_tmp)

    implicit none

    integer,                                  intent(in)    :: nE_tmp
    integer,             dimension(nE_tmp),   intent(in)    :: iperm_tmp
    double precision,    dimension(nE_tmp),   intent(in)    :: load_elmnts
    double precision,    dimension(nE_tmp),   intent(inout) :: sum_load

    integer                                                 :: i, k

    sum_load(:) = 0

!! DK DK Oct 2018: added these two safety tests
    if (nE_tmp < 0) stop 'error: negative nE_tmp in compute_sum_load(), this should not happen'
    if (nE_tmp == 0) stop 'error: null nE_tmp in compute_sum_load(), this should not happen'

    k = iperm_tmp(1)
    sum_load(1) = load_elmnts(k)
    do i = 2, nE_tmp
       k = iperm_tmp(i)
       sum_load(i) = sum_load(i-1) + load_elmnts(k)
    enddo

  end subroutine compute_sum_load

  !-----------
  ! compute the center of each elements
  !-----------
  subroutine compute_elmnts_center(elmnts_center, elmnts, nodes_coords, nspec, nnodes)

    implicit none

    ! input mesh
    integer,                                     intent(in)    :: nspec, nnodes
    double precision,    dimension(NDIM,nnodes),    intent(in)    :: nodes_coords
    integer,             dimension(NGNOD,nspec), intent(in)    :: elmnts
    double precision,    dimension(NDIM,nspec),     intent(inout) :: elmnts_center
    ! local parameters
    integer                                                    :: iE, i

    elmnts_center(:,:) = 0.d0
    do iE =1, nE
       do i = 1, NGNOD
          elmnts_center(1,iE) = elmnts_center(1,iE) + nodes_coords(1,elmnts(i,iE))
          elmnts_center(2,iE) = elmnts_center(2,iE) + nodes_coords(2,elmnts(i,iE))
          elmnts_center(3,iE) = elmnts_center(3,iE) + nodes_coords(3,elmnts(i,iE))
       enddo
       elmnts_center(1,iE) = elmnts_center(1,iE) / real(NGNOD,8)
       elmnts_center(2,iE) = elmnts_center(2,iE) / real(NGNOD,8)
       elmnts_center(3,iE) = elmnts_center(3,iE) / real(NGNOD,8)
    enddo


  end subroutine compute_elmnts_center

end module module_partition

