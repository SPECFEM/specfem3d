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

!----
!---- assemble the contributions between slices using non-blocking MPI
!----

  subroutine assemble_MPI_vector_blocking(NPROC,NGLOB_AB,array_val, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          my_neighbors_ext_mesh)

! assembles vector field in blocking way, only returns after values have been received and assembled

  use constants, only: NDIM,CUSTOM_REAL,itag

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: array_val

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh

  ! local parameters

  ! send/receive temporary buffers
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector

  ! requests
  integer, dimension(:), allocatable :: request_send_vector
  integer, dimension(:), allocatable :: request_recv_vector

  integer :: ipoin,iinterface,ier

! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    allocate(buffer_send_vector(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1985')
    if (ier /= 0) stop 'error allocating array buffer_send_vector'
    allocate(buffer_recv_vector(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1986')
    if (ier /= 0) stop 'error allocating array buffer_recv_vector'
    allocate(request_send_vector(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1987')
    if (ier /= 0) stop 'error allocating array request_send_vector'
    allocate(request_recv_vector(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1988')
    if (ier /= 0) stop 'error allocating array request_recv_vector'

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        buffer_send_vector(:,ipoin,iinterface) = array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      call isend_cr(buffer_send_vector(1,1,iinterface), &
                     NDIM*nibool_interfaces_ext_mesh(iinterface), &
                     my_neighbors_ext_mesh(iinterface), &
                     itag, &
                     request_send_vector(iinterface))
      call irecv_cr(buffer_recv_vector(1,1,iinterface), &
                     NDIM*nibool_interfaces_ext_mesh(iinterface), &
                     my_neighbors_ext_mesh(iinterface), &
                     itag, &
                     request_recv_vector(iinterface))
    enddo

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_vector(iinterface))
    enddo

    ! adding contributions of neighbors
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
             array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector(:,ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_vector(iinterface))
    enddo

    deallocate(buffer_send_vector)
    deallocate(buffer_recv_vector)
    deallocate(request_send_vector)
    deallocate(request_recv_vector)

  endif

  end subroutine assemble_MPI_vector_blocking

!
!-------------------------------------------------------------------------------------------------
!
  subroutine synchronize_MPI_vector_blocking_ord(NPROC,NGLOB_AB,array_val, &
                                                 num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                                 nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                                 my_neighbors_ext_mesh)

! kbai added this subroutine to synchronize a vector field
! to ensure that its values at nodes on MPI interfaces stay equal on all processors that share the node.
!
! Synchronize by setting the value to that of the processor with highest rank
! (it doesn't really matter which value we take, as long as all procs end up with exactly the same value).
!
! We assume that the interfaces are ordered by increasing rank of the neighbor.
! Uses blocking communication: only returns after values have been received and assembled

  use constants, only: NDIM,CUSTOM_REAL,itag,myrank

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: array_val

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  ! local parameters

  ! send/receive temporary buffers
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector
  ! requests
  integer, dimension(:), allocatable :: request_send_vector
  integer, dimension(:), allocatable :: request_recv_vector

  integer :: ipoin,iinterface,ier,iglob

! setting the value to that of the processor with highest rank

  if (NPROC > 1) then
    ! arrays for MPI halo exchange
    allocate(buffer_send_vector(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1989')
    if (ier /= 0) stop 'error allocating array buffer_send_vector'
    allocate(buffer_recv_vector(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1990')
    if (ier /= 0) stop 'error allocating array buffer_recv_vector'
    allocate(request_send_vector(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1991')
    if (ier /= 0) stop 'error allocating array request_send_vector'
    allocate(request_recv_vector(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1992')
    if (ier /= 0) stop 'error allocating array request_recv_vector'

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
        buffer_send_vector(:,ipoin,iinterface) = array_val(:,iglob)
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      call isend_cr(buffer_send_vector(1,1,iinterface), &
                     NDIM*nibool_interfaces_ext_mesh(iinterface), &
                     my_neighbors_ext_mesh(iinterface), &
                     itag, &
                     request_send_vector(iinterface))
      call irecv_cr(buffer_recv_vector(1,1,iinterface), &
                     NDIM*nibool_interfaces_ext_mesh(iinterface), &
                     my_neighbors_ext_mesh(iinterface), &
                     itag, &
                     request_recv_vector(iinterface))
    enddo

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_vector(iinterface))
    enddo

    ! set the value to that of the highest-rank processor
    do iinterface = 1, num_interfaces_ext_mesh
      if (myrank < my_neighbors_ext_mesh(iinterface)) then
        do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
          iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
          array_val(:,iglob) = buffer_recv_vector(:,ipoin,iinterface)
        enddo
      endif
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_vector(iinterface))
    enddo

    deallocate(buffer_send_vector)
    deallocate(buffer_recv_vector)
    deallocate(request_send_vector)
    deallocate(request_recv_vector)
  endif

  end subroutine synchronize_MPI_vector_blocking_ord

!
!------------------------------------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_async_send(NPROC,NGLOB_AB,array_val, &
                                            buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            my_neighbors_ext_mesh, &
                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

! sends data

  use constants, only: NDIM,CUSTOM_REAL,itag

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: array_val

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
       buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  ! local parameters
  integer ipoin,iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        buffer_send_vector_ext_mesh(:,ipoin,iinterface) = &
             array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      call isend_cr(buffer_send_vector_ext_mesh(1,1,iinterface), &
                    NDIM*nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_send_vector_ext_mesh(iinterface))
      call irecv_cr(buffer_recv_vector_ext_mesh(1,1,iinterface), &
                    NDIM*nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_recv_vector_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_async_send

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_async_recv(NPROC,NGLOB_AB,array_val, &
                                            buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                            my_neighbors_ext_mesh)

! waits for send/receiver to be completed and assembles contributions

  use constants, only: NDIM,CUSTOM_REAL
  use specfem_par, only: FAULT_SIMULATION

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: array_val

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
    buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(in) :: my_neighbors_ext_mesh

  ! local parameters
  integer :: ipoin,iinterface,iglob

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC == 1) return

  ! fault ruptures
  if (FAULT_SIMULATION) then
    ! receives MPI buffers with ordered assembly
    call assemble_MPI_vector_async_w_ord(NPROC,NGLOB_AB,array_val, &
                                         buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                         max_nibool_interfaces_ext_mesh, &
                                         nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                         request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                         my_neighbors_ext_mesh)
    ! all done
    return
  endif

  ! wait for communications completion (recv)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_recv_vector_ext_mesh(iinterface))
  enddo

  ! adding all neighbor contributions
  do iinterface = 1, num_interfaces_ext_mesh
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
      array_val(:,iglob) = array_val(:,iglob) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
    enddo
  enddo

  ! wait for communications completion (send)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_send_vector_ext_mesh(iinterface))
  enddo

  end subroutine assemble_MPI_vector_async_recv


!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_async_w_ord(NPROC,NGLOB_AB,array_val, &
                                             buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                             max_nibool_interfaces_ext_mesh, &
                                             nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                             request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                             my_neighbors_ext_mesh)

! waits for data to receive and assembles

! The goal of this version is to avoid different round-off errors in different processors.
! The contribution of each processor is added following the order of its rank.
! This guarantees that the sums are done in the same order on all processors.
!
! NOTE: this version assumes that the interfaces are ordered by increasing rank of the neighbor.
! That is currently done so in subroutine write_interfaces_database in decompose_mesh_SCOTCH/part_decompose_mesh_SCOTCH.f90
! A safety test could be added here.
!
! October 2012 - Surendra Somala and Jean-Paul Ampuero - Caltech Seismolab

  use constants, only: NDIM,CUSTOM_REAL,myrank

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: array_val

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
       buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(in) :: my_neighbors_ext_mesh

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: mybuffer
  integer :: ipoin,iinterface,iglob
  logical :: need_add_my_contrib

! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC == 1) return

  ! move interface values of array_val to local buffers
  do iinterface = 1, num_interfaces_ext_mesh
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
      mybuffer(:,ipoin,iinterface) = array_val(:,iglob)
      ! set them to zero right away to avoid counting it more than once during assembly:
      ! buffers of higher rank get zeros on nodes shared with current buffer
      array_val(:,iglob) = 0._CUSTOM_REAL
    enddo
  enddo

  ! wait for communications completion (recv)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_recv_vector_ext_mesh(iinterface))
  enddo

  ! adding all contributions in order of processor rank
  need_add_my_contrib = .true.
  do iinterface = 1, num_interfaces_ext_mesh
    if (need_add_my_contrib .and. myrank < my_neighbors_ext_mesh(iinterface)) call add_my_contrib()
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
      array_val(:,iglob) = array_val(:,iglob) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
    enddo
  enddo
  if (need_add_my_contrib) call add_my_contrib()

  ! wait for communications completion (send)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_send_vector_ext_mesh(iinterface))
  enddo

  contains

    subroutine add_my_contrib()

    integer :: my_iinterface,my_ipoin

    do my_iinterface = 1, num_interfaces_ext_mesh
      do my_ipoin = 1, nibool_interfaces_ext_mesh(my_iinterface)
        iglob = ibool_interfaces_ext_mesh(my_ipoin,my_iinterface)
        array_val(:,iglob) = array_val(:,iglob) + mybuffer(:,my_ipoin,my_iinterface)
      enddo
    enddo
    need_add_my_contrib = .false.

    end subroutine add_my_contrib

  end subroutine assemble_MPI_vector_async_w_ord


!
!--------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_poro_s(NPROC,NGLOB_AB,array_val1,array_val2, &
                                        buffer_send_vector_ext_mesh_s,buffer_recv_vector_ext_mesh_s, &
                                        buffer_send_vector_ext_mesh_w,buffer_recv_vector_ext_mesh_w, &
                                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,my_neighbors_ext_mesh, &
                                        request_send_vector_ext_mesh_s,request_recv_vector_ext_mesh_s, &
                                        request_send_vector_ext_mesh_w,request_recv_vector_ext_mesh_w)

! sends data

  use constants, only: NDIM,CUSTOM_REAL,itag

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB
  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: array_val1,array_val2

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
       buffer_send_vector_ext_mesh_s,buffer_recv_vector_ext_mesh_s, &
       buffer_send_vector_ext_mesh_w,buffer_recv_vector_ext_mesh_w

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh_s,request_recv_vector_ext_mesh_s
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh_w,request_recv_vector_ext_mesh_w

  ! local parameters
  integer :: ipoin,iinterface

! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        buffer_send_vector_ext_mesh_s(:,ipoin,iinterface) = array_val1(:,ibool_interfaces_ext_mesh(ipoin,iinterface))
        buffer_send_vector_ext_mesh_w(:,ipoin,iinterface) = array_val2(:,ibool_interfaces_ext_mesh(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      !val1
      call isend_cr(buffer_send_vector_ext_mesh_s(1,1,iinterface), &
                    NDIM*nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_send_vector_ext_mesh_s(iinterface))
      call irecv_cr(buffer_recv_vector_ext_mesh_s(1,1,iinterface), &
                    NDIM*nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_recv_vector_ext_mesh_s(iinterface))
      !val2
      call isend_cr(buffer_send_vector_ext_mesh_w(1,1,iinterface), &
                    NDIM*nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_send_vector_ext_mesh_w(iinterface))
      call irecv_cr(buffer_recv_vector_ext_mesh_w(1,1,iinterface), &
                    NDIM*nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_recv_vector_ext_mesh_w(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_poro_s

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_poro_w(NPROC,NGLOB_AB,array_val1,array_val2, &
                                        buffer_recv_vector_ext_mesh_s,buffer_recv_vector_ext_mesh_w, &
                                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                        request_send_vector_ext_mesh_s,request_recv_vector_ext_mesh_s, &
                                        request_send_vector_ext_mesh_w,request_recv_vector_ext_mesh_w)

! waits for data to receive and assembles

  use constants, only: NDIM,CUSTOM_REAL

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: array_val1,array_val2

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
       buffer_recv_vector_ext_mesh_s,buffer_recv_vector_ext_mesh_w

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh_s,request_recv_vector_ext_mesh_s
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh_w,request_recv_vector_ext_mesh_w

  ! local parameters
  integer :: ipoin,iinterface

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_vector_ext_mesh_s(iinterface))
      call wait_req(request_recv_vector_ext_mesh_w(iinterface))
    enddo

    ! adding contributions of neighbors
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        array_val1(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
             array_val1(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh_s(:,ipoin,iinterface)
        array_val2(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
             array_val2(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh_w(:,ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_vector_ext_mesh_s(iinterface))
      call wait_req(request_send_vector_ext_mesh_w(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_poro_w



!
!-------------------------------------------------------------------------------------------------
!

! with cuda functions...

  subroutine transfer_boundary_to_device(NPROC, Mesh_pointer, &
                                         buffer_recv_vector_ext_mesh, &
                                         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                         request_recv_vector_ext_mesh)

  use constants, only: NDIM,CUSTOM_REAL

  implicit none

  integer,intent(in) :: NPROC
  integer(kind=8),intent(in) :: Mesh_pointer

  ! array to assemble
  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
       buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_recv_vector_ext_mesh

  ! local parameters
  integer :: iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    !write(IMAIN,*) "sending MPI_wait"
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_vector_ext_mesh(iinterface))
    enddo

    ! send contributions to GPU
    call transfer_boundary_to_device_a(Mesh_pointer, buffer_recv_vector_ext_mesh, &
                                       num_interfaces_ext_mesh, max_nibool_interfaces_ext_mesh)
  endif

  ! This step is done via previous function transfer_and_assemble...
  ! do iinterface = 1, num_interfaces_ext_mesh
  !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
  !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
  !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
  !   enddo
  ! enddo

  end subroutine transfer_boundary_to_device

!
!-------------------------------------------------------------------------------------------------
!


  subroutine assemble_MPI_vector_send_cuda(NPROC, &
                                           buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

! sends data
! note: array to assemble already filled into buffer_send_vector_ext_mesh array

  use constants, only: NDIM,CUSTOM_REAL,itag

  implicit none

  integer,intent(in) :: NPROC

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
       buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  ! local parameters
  integer :: iinterface

  ! note: preparation of the contribution between partitions using MPI
  !          already done in transfer_boun_accel routine

  ! send only if more than one partition
  if (NPROC > 1) then

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      call isend_cr(buffer_send_vector_ext_mesh(1,1,iinterface), &
                     NDIM*nibool_interfaces_ext_mesh(iinterface), &
                     my_neighbors_ext_mesh(iinterface), &
                     itag, &
                     request_send_vector_ext_mesh(iinterface))
      call irecv_cr(buffer_recv_vector_ext_mesh(1,1,iinterface), &
                     NDIM*nibool_interfaces_ext_mesh(iinterface), &
                     my_neighbors_ext_mesh(iinterface), &
                     itag, &
                     request_recv_vector_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_send_cuda


!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_write_cuda(NPROC,NGLOB_AB,array_val, Mesh_pointer, &
                                            buffer_recv_vector_ext_mesh, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                            FORWARD_OR_ADJOINT)

! waits for data to receive and assembles

  use constants, only: NDIM,CUSTOM_REAL

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB
  integer(kind=8),intent(in) :: Mesh_pointer

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: array_val

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
       buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer,intent(in) :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface
  real(kind=CUSTOM_REAL) :: dummy_val ! to avoid compiler warning

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_vector_ext_mesh(iinterface))
    enddo

    ! adding contributions of neighbors
    call transfer_asmbl_accel_to_device(Mesh_pointer, &
                                        buffer_recv_vector_ext_mesh, &
                                        num_interfaces_ext_mesh, max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh, &
                                        ibool_interfaces_ext_mesh,FORWARD_OR_ADJOINT)

    ! to avoid compiler warning
    dummy_val = array_val(1,1)

    ! This step is done via previous function transfer_and_assemble...
    ! do iinterface = 1, num_interfaces_ext_mesh
    !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
    !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
    !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
    !   enddo
    ! enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_vector_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_write_cuda


!
!-------------------------------------------------------------------------------------------------
!


  subroutine synchronize_MPI_vector_write_cuda(NPROC,NGLOB_AB,array_val, Mesh_pointer, &
                                               buffer_recv_vector_ext_mesh, &
                                               num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                               nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                               request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                               FORWARD_OR_ADJOINT)

! waits for data to receive and assembles

  use constants, only: NDIM,CUSTOM_REAL

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB
  integer(kind=8),intent(in) :: Mesh_pointer

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: array_val

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
       buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer,intent(in) :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface
  real(kind=CUSTOM_REAL) :: dummy_val ! to avoid compiler warnings

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_vector_ext_mesh(iinterface))
    enddo

    ! adding contributions of neighbors
    call transfer_sync_accel_to_device(Mesh_pointer, &
                                       buffer_recv_vector_ext_mesh, &
                                       num_interfaces_ext_mesh, max_nibool_interfaces_ext_mesh, &
                                       nibool_interfaces_ext_mesh, &
                                       ibool_interfaces_ext_mesh,FORWARD_OR_ADJOINT)

    ! to avoid compiler warning
    dummy_val = array_val(1,1)

    ! This step is done via previous function transfer_and_assemble...
    ! do iinterface = 1, num_interfaces_ext_mesh
    !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
    !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
    !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
    !   enddo
    ! enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_vector_ext_mesh(iinterface))
    enddo

  endif

  end subroutine synchronize_MPI_vector_write_cuda


!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_send_cuda(NPROC, &
                                           buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

! non-blocking MPI send

  ! sends data
  ! note: assembling data already filled into buffer_send_scalar_ext_mesh array

  use constants, only: CUSTOM_REAL,NB_RUNS_ACOUSTIC_GPU,itag

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh*NB_RUNS_ACOUSTIC_GPU,num_interfaces_ext_mesh),intent(inout) :: &
    buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  ! local parameters
  integer :: iinterface

  ! sends only if more than one partition
  if (NPROC > 1) then

    ! note: partition border copy into the buffer has already been done
    !          by routine transfer_boun_pot_from_device()

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      ! note: passing arguments:
      !          **array**(1:nibool_interfaces_ext_mesh(iinterface),iinterface)
      !       might lead to additional copy memory operations for certain compilers (slows down performance);
      !       to avoid this in Fortran, one might just pass the pointer to the first array value:
      !          **array**(1,iinterface)

      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar_ext_mesh(1,iinterface), &
                    nibool_interfaces_ext_mesh(iinterface)*NB_RUNS_ACOUSTIC_GPU, &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_send_scalar_ext_mesh(iinterface))
      ! receive request
      call irecv_cr(buffer_recv_scalar_ext_mesh(1,iinterface), &
                    nibool_interfaces_ext_mesh(iinterface)*NB_RUNS_ACOUSTIC_GPU, &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_recv_scalar_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_send_cuda

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,array_val, &
                                            Mesh_pointer, &
                                            buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
!                                             nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                                            FORWARD_OR_ADJOINT)

! waits for send/receiver to be completed and assembles contributions

  use constants, only: CUSTOM_REAL,NB_RUNS_ACOUSTIC_GPU

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB
  integer(kind=8),intent(in) :: Mesh_pointer

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: array_val

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh*NB_RUNS_ACOUSTIC_GPU,num_interfaces_ext_mesh),intent(inout) :: &
       buffer_recv_scalar_ext_mesh

! integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
! integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  integer :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface ! ipoin
  real(kind=CUSTOM_REAL) :: dummy_val ! to avoid compiler warning

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_scalar_ext_mesh(iinterface))
    enddo

    ! adding contributions of neighbors
    call transfer_asmbl_pot_to_device(Mesh_pointer,buffer_recv_scalar_ext_mesh,FORWARD_OR_ADJOINT)

    ! to avoid compiler warning
    dummy_val = array_val(1)

    ! note: adding contributions of neighbors has been done just above for cuda
    !do iinterface = 1, num_interfaces_ext_mesh
    !  do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
    !    array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
    !         array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) &
    !         + buffer_recv_scalar_ext_mesh(ipoin,iinterface)
    !  enddo
    !enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_scalar_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_write_cuda

