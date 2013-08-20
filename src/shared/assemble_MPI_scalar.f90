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

!----
!---- assemble the contributions between slices using non-blocking MPI
!----

  subroutine assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,array_val, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)

! assembles scalar field in a blocking way, returns only after values have been assembled

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_scalar
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_recv_scalar
  integer, dimension(:), allocatable :: request_send_scalar
  integer, dimension(:), allocatable :: request_recv_scalar


  integer ipoin,iinterface,ier

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if(NPROC > 1) then

    allocate(buffer_send_scalar(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array buffer_send_scalar'
    allocate(buffer_recv_scalar(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array buffer_recv_scalar'
    allocate(request_send_scalar(num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array request_send_scalar'
    allocate(request_recv_scalar(num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array request_recv_scalar'

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        buffer_send_scalar(ipoin,iinterface) = array_val(ibool_interfaces_ext_mesh(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
                     nibool_interfaces_ext_mesh(iinterface), &
                     my_neighbours_ext_mesh(iinterface), &
                     itag, &
                     request_send_scalar(iinterface) )
      ! receive request
      call irecv_cr(buffer_recv_scalar(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
                     nibool_interfaces_ext_mesh(iinterface), &
                     my_neighbours_ext_mesh(iinterface), &
                     itag, &
                     request_recv_scalar(iinterface) )
    enddo

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_scalar(iinterface))
    enddo

    ! adding contributions of neighbours
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
             array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_scalar(ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_scalar(iinterface))
    enddo

    deallocate(buffer_send_scalar)
    deallocate(buffer_recv_scalar)
    deallocate(request_send_scalar)
    deallocate(request_recv_scalar)

  endif

  end subroutine assemble_MPI_scalar_blocking

!
!----
!

  subroutine assemble_MPI_scalar_i_blocking(NPROC,NGLOB_AB,array_val, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB

! array to assemble
  integer, dimension(NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

  integer, dimension(:,:), allocatable :: buffer_send_scalar
  integer, dimension(:,:), allocatable :: buffer_recv_scalar
  integer, dimension(:), allocatable :: request_send_scalar
  integer, dimension(:), allocatable :: request_recv_scalar

  integer :: ipoin,iinterface,ier

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if(NPROC > 1) then

    allocate(buffer_send_scalar(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array buffer_send_scalar'
    allocate(buffer_recv_scalar(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array buffer_recv_scalar'
    allocate(request_send_scalar(num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array request_send_scalar'
    allocate(request_recv_scalar(num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array request_recv_scalar'

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        buffer_send_scalar(ipoin,iinterface) = array_val(ibool_interfaces_ext_mesh(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      ! non-blocking synchronous send request
      call isend_i(buffer_send_scalar(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
                   nibool_interfaces_ext_mesh(iinterface), &
                   my_neighbours_ext_mesh(iinterface), &
                   itag, &
                   request_send_scalar(iinterface) )
      ! receive request
      call irecv_i(buffer_recv_scalar(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
                   nibool_interfaces_ext_mesh(iinterface), &
                   my_neighbours_ext_mesh(iinterface), &
                   itag, &
                   request_recv_scalar(iinterface) )
    enddo

    ! wait for communications completion
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_scalar(iinterface))
    enddo

    ! adding contributions of neighbours
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
             array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_scalar(ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_scalar(iinterface))
    enddo

    deallocate(buffer_send_scalar)
    deallocate(buffer_recv_scalar)
    deallocate(request_send_scalar)
    deallocate(request_recv_scalar)

  endif

  end subroutine assemble_MPI_scalar_i_blocking

!
!----
!

  subroutine assemble_MPI_scalar_async_send(NPROC,NGLOB_AB,array_val, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

! non-blocking MPI send

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

! array to send
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val


  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  integer ipoin,iinterface

! sends only if more than one partition
  if(NPROC > 1) then

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        buffer_send_scalar_ext_mesh(ipoin,iinterface) = &
          array_val(ibool_interfaces_ext_mesh(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
           nibool_interfaces_ext_mesh(iinterface), &
           my_neighbours_ext_mesh(iinterface), &
           itag, &
           request_send_scalar_ext_mesh(iinterface) &
           )
      ! receive request
      call irecv_cr(buffer_recv_scalar_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
           nibool_interfaces_ext_mesh(iinterface), &
           my_neighbours_ext_mesh(iinterface), &
           itag, &
           request_recv_scalar_ext_mesh(iinterface) &
           )

    enddo

  endif

  end subroutine assemble_MPI_scalar_async_send

!
!----
!

  subroutine assemble_MPI_scalar_async_recv(NPROC,NGLOB_AB,array_val, &
                        buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

! waits for send/receiver to be completed and assembles contributions

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val


  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_scalar_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  integer ipoin,iinterface

! assemble only if more than one partition
  if(NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_scalar_ext_mesh(iinterface))
    enddo

    ! adding contributions of neighbours
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
             array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) &
             + buffer_recv_scalar_ext_mesh(ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_scalar_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_async_recv

