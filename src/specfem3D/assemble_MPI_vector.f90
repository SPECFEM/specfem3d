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

  subroutine assemble_MPI_vector_blocking(NPROC,NGLOB_AB,array_val, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)

! assembles vector field in blocking way, only returns after values have been received and assembled

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

  ! local parameters

  ! send/receive temporary buffers
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector

  ! requests
  integer, dimension(:), allocatable :: request_send_vector
  integer, dimension(:), allocatable :: request_recv_vector

  integer ipoin,iinterface,ier


! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if(NPROC > 1) then

    allocate(buffer_send_vector(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array buffer_send_vector'
    allocate(buffer_recv_vector(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array buffer_recv_vector'
    allocate(request_send_vector(num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array request_send_vector'
    allocate(request_recv_vector(num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array request_recv_vector'

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
                     my_neighbours_ext_mesh(iinterface), &
                     itag, &
                     request_send_vector(iinterface) )
      call irecv_cr(buffer_recv_vector(1,1,iinterface), &
                     NDIM*nibool_interfaces_ext_mesh(iinterface), &
                     my_neighbours_ext_mesh(iinterface), &
                     itag, &
                     request_recv_vector(iinterface) )
    enddo

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_vector(iinterface))
    enddo

    ! adding contributions of neighbours
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

  subroutine assemble_MPI_vector_async_send(NPROC,NGLOB_AB,array_val, &
                                           buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                           my_neighbours_ext_mesh, &
                                           request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

    ! sends data

    implicit none

    include "constants.h"

    integer :: NPROC
    integer :: NGLOB_AB

    ! array to assemble
    real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val

    integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

    real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
         buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh

    integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
    integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
    integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

    integer ipoin,iinterface

    ! here we have to assemble all the contributions between partitions using MPI

    ! assemble only if more than one partition
    if(NPROC > 1) then

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
               my_neighbours_ext_mesh(iinterface), &
               itag, &
               request_send_vector_ext_mesh(iinterface) &
               )
          call irecv_cr(buffer_recv_vector_ext_mesh(1,1,iinterface), &
               NDIM*nibool_interfaces_ext_mesh(iinterface), &
               my_neighbours_ext_mesh(iinterface), &
               itag, &
               request_recv_vector_ext_mesh(iinterface) &
               )
       enddo

    endif

  end subroutine assemble_MPI_vector_async_send

!
!-------------------------------------------------------------------------------------------------
!

! unused, new ordered routine is used now...
!
!  subroutine assemble_MPI_vector_async_recv(NPROC,NGLOB_AB,array_val, &
!                                            buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
!                                            max_nibool_interfaces_ext_mesh, &
!                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
!                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
!
!! waits for data to receive and assembles
!
!  implicit none
!
!  include "constants.h"
!
!  integer :: NPROC
!  integer :: NGLOB_AB
!
!! array to assemble
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val
!
!  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
!
!  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
!       buffer_recv_vector_ext_mesh
!
!  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
!  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
!  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh
!
!  integer ipoin,iinterface
!
!! here we have to assemble all the contributions between partitions using MPI
!
!! assemble only if more than one partition
!  if(NPROC > 1) then
!
!! wait for communications completion (recv)
!  do iinterface = 1, num_interfaces_ext_mesh
!    call wait_req(request_recv_vector_ext_mesh(iinterface))
!  enddo
!
!! adding contributions of neighbours
!  do iinterface = 1, num_interfaces_ext_mesh
!    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
!      array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
!           array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
!    enddo
!  enddo
!
!! wait for communications completion (send)
!  do iinterface = 1, num_interfaces_ext_mesh
!    call wait_req(request_send_vector_ext_mesh(iinterface))
!  enddo
!
!  endif
!
!  end subroutine assemble_MPI_vector_async_recv
!
!

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_async_w_ord(NPROC,NGLOB_AB,array_val, &
                                            buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                            my_neighbours_ext_mesh,myrank)

! waits for data to receive and assembles

! The goal of this version is to avoid different round-off errors in different processors.
! The contribution of each processor is added following the order of its rank.
! This guarantees that the sums are done in the same order on all processors.
!
! NOTE: this version assumes that the interfaces are ordered by increasing rank of the neighbour.
! That is currently done so in subroutine write_interfaces_database in decompose_mesh_SCOTCH/part_decompose_mesh_SCOTCH.f90
! A safety test could be added here.
!
! October 2012 - Surendra Somala and Jean-Paul Ampuero - Caltech Seismolab

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,myrank

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh

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
    if (need_add_my_contrib .and. myrank < my_neighbours_ext_mesh(iinterface)) call add_my_contrib()
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
                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,my_neighbours_ext_mesh, &
                      request_send_vector_ext_mesh_s,request_recv_vector_ext_mesh_s, &
                      request_send_vector_ext_mesh_w,request_recv_vector_ext_mesh_w &
                      )

! sends data

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val1,array_val2

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_vector_ext_mesh_s,buffer_recv_vector_ext_mesh_s,&
       buffer_send_vector_ext_mesh_w,buffer_recv_vector_ext_mesh_w

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh_s,request_recv_vector_ext_mesh_s
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh_w,request_recv_vector_ext_mesh_w

  integer ipoin,iinterface

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if(NPROC > 1) then

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
         my_neighbours_ext_mesh(iinterface), &
         itag, &
         request_send_vector_ext_mesh_s(iinterface) &
         )
    call irecv_cr(buffer_recv_vector_ext_mesh_s(1,1,iinterface), &
         NDIM*nibool_interfaces_ext_mesh(iinterface), &
         my_neighbours_ext_mesh(iinterface), &
         itag, &
         request_recv_vector_ext_mesh_s(iinterface) &
         )
!val2
    call isend_cr(buffer_send_vector_ext_mesh_w(1,1,iinterface), &
         NDIM*nibool_interfaces_ext_mesh(iinterface), &
         my_neighbours_ext_mesh(iinterface), &
         itag, &
         request_send_vector_ext_mesh_w(iinterface) &
         )
    call irecv_cr(buffer_recv_vector_ext_mesh_w(1,1,iinterface), &
         NDIM*nibool_interfaces_ext_mesh(iinterface), &
         my_neighbours_ext_mesh(iinterface), &
         itag, &
         request_recv_vector_ext_mesh_w(iinterface) &
         )
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

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val1,array_val2

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh_s,buffer_recv_vector_ext_mesh_w

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh_s,request_recv_vector_ext_mesh_s
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh_w,request_recv_vector_ext_mesh_w

  integer ipoin,iinterface

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if(NPROC > 1) then

! wait for communications completion (recv)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_recv_vector_ext_mesh_s(iinterface))
    call wait_req(request_recv_vector_ext_mesh_w(iinterface))
  enddo

! adding contributions of neighbours
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
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,&
                                            request_recv_vector_ext_mesh)

  implicit none

  include "constants.h"

  integer :: NPROC
  integer(kind=8) :: Mesh_pointer

  ! array to assemble
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: request_recv_vector_ext_mesh

  ! local parameters
  integer :: iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then

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

! not used...
!  subroutine assemble_MPI_vector_write_cuda_no_transfer(NPROC,NGLOB_AB,array_val, Mesh_pointer, &
!                                            buffer_recv_vector_ext_mesh, &
!                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
!                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
!                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
!                                            FORWARD_OR_ADJOINT )
!
!  implicit none
!
!  include "constants.h"
!
!  integer :: NPROC
!  integer :: NGLOB_AB
!  integer(kind=8) :: Mesh_pointer
!  ! array to assemble
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val
!
!  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
!
!  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
!       buffer_recv_vector_ext_mesh
!
!  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
!  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
!  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh
!  !integer, dimension(num_interfaces_ext_mesh) :: request_recv_vector_ext_mesh
!
!  integer :: FORWARD_OR_ADJOINT
!
!  ! local parameters
!  integer :: iinterface
!
!  ! here we have to assemble all the contributions between partitions using MPI
!
!  ! assemble only if more than one partition
!  if(NPROC > 1) then
!
!     ! adding contributions of neighbours
!     call assemble_accel_on_device(Mesh_pointer, array_val, buffer_recv_vector_ext_mesh, &
!          num_interfaces_ext_mesh, max_nibool_interfaces_ext_mesh, &
!          nibool_interfaces_ext_mesh,&
!          ibool_interfaces_ext_mesh,FORWARD_OR_ADJOINT)
!
!     ! This step is done via previous function transfer_and_assemble...
!     ! do iinterface = 1, num_interfaces_ext_mesh
!     !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
!     !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
!     !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
!     !   enddo
!     ! enddo
!
!     ! wait for communications completion (send)
!     do iinterface = 1, num_interfaces_ext_mesh
!        call wait_req(request_send_vector_ext_mesh(iinterface))
!     enddo
!  endif
!
!  end subroutine assemble_MPI_vector_write_cuda_no_transfer

!
!-------------------------------------------------------------------------------------------------
!


  subroutine assemble_MPI_vector_send_cuda(NPROC, &
                                          buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh, &
                                          my_neighbours_ext_mesh, &
                                          request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

! sends data
! note: array to assemble already filled into buffer_send_vector_ext_mesh array

  implicit none

  include "constants.h"

  integer :: NPROC

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  ! local parameters
  integer :: iinterface

  ! note: preparation of the contribution between partitions using MPI
  !          already done in transfer_boun_accel routine

  ! send only if more than one partition
  if(NPROC > 1) then

     ! send messages
     do iinterface = 1, num_interfaces_ext_mesh
        call isend_cr(buffer_send_vector_ext_mesh(1,1,iinterface), &
                       NDIM*nibool_interfaces_ext_mesh(iinterface), &
                       my_neighbours_ext_mesh(iinterface), &
                       itag, &
                       request_send_vector_ext_mesh(iinterface) )
        call irecv_cr(buffer_recv_vector_ext_mesh(1,1,iinterface), &
                       NDIM*nibool_interfaces_ext_mesh(iinterface), &
                       my_neighbours_ext_mesh(iinterface), &
                       itag, &
                       request_recv_vector_ext_mesh(iinterface) )
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
                                            FORWARD_OR_ADJOINT )

! waits for data to receive and assembles

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB
  integer(kind=8) :: Mesh_pointer

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_vector_ext_mesh(iinterface))
    enddo

    ! adding contributions of neighbours
    call transfer_asmbl_accel_to_device(Mesh_pointer, array_val, &
                                        buffer_recv_vector_ext_mesh, &
                                        num_interfaces_ext_mesh, max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh,&
                                        ibool_interfaces_ext_mesh,FORWARD_OR_ADJOINT)

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

  subroutine assemble_MPI_scalar_send_cuda(NPROC, &
                                           buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbours_ext_mesh, &
                                           request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

! non-blocking MPI send

  ! sends data
  ! note: assembling data already filled into buffer_send_scalar_ext_mesh array
  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  ! local parameters
  integer :: iinterface

  ! sends only if more than one partition
  if(NPROC > 1) then

    ! note: partition border copy into the buffer has already been done
    !          by routine transfer_boun_pot_from_device()

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      ! note: passing arguments:
      !          **array**(1:nibool_interfaces_ext_mesh(iinterface),iinterface)
      !       might lead to additional copy memory operations for certain compilers (slows down performance);
      !       to avoid this in fortran, one might just pass the pointer to the first array value:
      !          **array**(1,iinterface)

      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
                    nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbours_ext_mesh(iinterface), &
                    itag, &
                    request_send_scalar_ext_mesh(iinterface) )
      ! receive request
      call irecv_cr(buffer_recv_scalar_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
                    nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbours_ext_mesh(iinterface), &
                    itag, &
                    request_recv_scalar_ext_mesh(iinterface) )
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
!                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                        FORWARD_OR_ADJOINT)

! waits for send/receiver to be completed and assembles contributions

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB
  integer(kind=8) :: Mesh_pointer

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val


  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_scalar_ext_mesh

! integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
! integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  integer :: FORWARD_OR_ADJOINT

  integer :: iinterface ! ipoin

  ! assemble only if more than one partition
  if(NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_scalar_ext_mesh(iinterface))
    enddo

    ! adding contributions of neighbours
    call transfer_asmbl_pot_to_device(Mesh_pointer,array_val,buffer_recv_scalar_ext_mesh,FORWARD_OR_ADJOINT)

    ! note: adding contributions of neighbours has been done just above for cuda
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

