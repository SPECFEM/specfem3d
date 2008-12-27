!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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
!---- assemble the contributions between slices and chunks using MPI
!----

  subroutine assemble_MPI_scalar(array_val,iproc_xi,iproc_eta,addressing, &
            iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_xi,npoin2D_eta, &
            NPROC_XI,NPROC_ETA,NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NPOIN2DMAX_XY)

  implicit none

  include "constants.h"

! include values created by the mesher
  include "OUTPUT_FILES/values_from_mesher.h"

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB_VAL) :: array_val

  integer iproc_xi,iproc_eta
  integer npoin2D_xi,npoin2D_eta

  integer NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NPOIN2DMAX_XY
  integer NPROC_XI,NPROC_ETA

! for addressing of the slices
  integer, dimension(0:NPROC_XI-1,0:NPROC_ETA-1) :: addressing

! 2-D addressing and buffers for summation between slices
  integer, dimension(NPOIN2DMAX_XMIN_XMAX) :: iboolleft_xi,iboolright_xi
  integer, dimension(NPOIN2DMAX_YMIN_YMAX) :: iboolleft_eta,iboolright_eta

  real(kind=CUSTOM_REAL), dimension(NPOIN2DMAX_XY) :: buffer_send_faces_scalar,buffer_received_faces_scalar

  integer ipoin
  integer sender,receiver

  integer, external :: proc_null

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! here we have to assemble all the contributions between slices using MPI

!----
!---- first assemble along xi using the 2-D topology
!----

! assemble along xi only if more than one slice
  if(NPROC_XI > 1) then

! slices copy the right face into the buffer
  do ipoin=1,npoin2D_xi
    buffer_send_faces_scalar(ipoin) = array_val(iboolright_xi(ipoin))
  enddo

! send messages forward along each row
  if(iproc_xi == 0) then
    sender = proc_null()
  else
    sender = addressing(iproc_xi - 1,iproc_eta)
  endif
  if(iproc_xi == NPROC_XI-1) then
    receiver = proc_null()
  else
    receiver = addressing(iproc_xi + 1,iproc_eta)
  endif
  call sendrecv_all_cr(buffer_send_faces_scalar,npoin2D_xi,receiver,itag2, &
                       buffer_received_faces_scalar,npoin2D_xi,sender,itag)

! all slices add the buffer received to the contributions on the left face
  if(iproc_xi > 0) then
  do ipoin=1,npoin2D_xi
      array_val(iboolleft_xi(ipoin)) = array_val(iboolleft_xi(ipoin)) + &
                                buffer_received_faces_scalar(ipoin)
  enddo
  endif

! the contributions are correctly assembled on the left side of each slice
! now we have to send the result back to the sender
! all slices copy the left face into the buffer
  do ipoin=1,npoin2D_xi
    buffer_send_faces_scalar(ipoin) = array_val(iboolleft_xi(ipoin))
  enddo

! send messages backward along each row
  if(iproc_xi == NPROC_XI-1) then
    sender = proc_null()
  else
    sender = addressing(iproc_xi + 1,iproc_eta)
  endif
  if(iproc_xi == 0) then
    receiver = proc_null()
  else
    receiver = addressing(iproc_xi - 1,iproc_eta)
  endif
  call sendrecv_all_cr(buffer_send_faces_scalar,npoin2D_xi,receiver,itag2, &
                       buffer_received_faces_scalar,npoin2D_xi,sender,itag)

! all slices copy the buffer received to the contributions on the right face
  if(iproc_xi < NPROC_XI-1) then
  do ipoin=1,npoin2D_xi
    array_val(iboolright_xi(ipoin)) = buffer_received_faces_scalar(ipoin)
  enddo
  endif

  endif

!----
!---- then assemble along eta using the 2-D topology
!----

! assemble along eta only if more than one slice
  if(NPROC_ETA > 1) then

! slices copy the right face into the buffer
  do ipoin=1,npoin2D_eta
    buffer_send_faces_scalar(ipoin) = array_val(iboolright_eta(ipoin))
  enddo

! send messages forward along each row
  if(iproc_eta == 0) then
    sender = proc_null()
  else
    sender = addressing(iproc_xi,iproc_eta - 1)
  endif
  if(iproc_eta == NPROC_ETA-1) then
    receiver = proc_null()
  else
    receiver = addressing(iproc_xi,iproc_eta + 1)
  endif
  call sendrecv_all_cr(buffer_send_faces_scalar,npoin2D_eta,receiver,itag2, &
                       buffer_received_faces_scalar,npoin2D_eta,sender,itag)

! all slices add the buffer received to the contributions on the left face
  if(iproc_eta > 0) then
  do ipoin=1,npoin2D_eta
      array_val(iboolleft_eta(ipoin)) = array_val(iboolleft_eta(ipoin)) + &
                                buffer_received_faces_scalar(ipoin)
  enddo
  endif

! the contributions are correctly assembled on the left side of each slice
! now we have to send the result back to the sender
! all slices copy the left face into the buffer
  do ipoin=1,npoin2D_eta
    buffer_send_faces_scalar(ipoin) = array_val(iboolleft_eta(ipoin))
  enddo

! send messages backward along each row
  if(iproc_eta == NPROC_ETA-1) then
    sender = proc_null()
  else
    sender = addressing(iproc_xi,iproc_eta + 1)
  endif
  if(iproc_eta == 0) then
    receiver = proc_null()
  else
    receiver = addressing(iproc_xi,iproc_eta - 1)
  endif
  call sendrecv_all_cr(buffer_send_faces_scalar,npoin2D_eta,receiver,itag2, &
                       buffer_received_faces_scalar,npoin2D_eta,sender,itag)

! all slices copy the buffer received to the contributions on the right face
  if(iproc_eta < NPROC_ETA-1) then
  do ipoin=1,npoin2D_eta
    array_val(iboolright_eta(ipoin)) = buffer_received_faces_scalar(ipoin)
  enddo
  endif

  endif

  end subroutine assemble_MPI_scalar

!
!----
!

  subroutine assemble_MPI_scalar_ext_mesh(NPROC,NGLOB_AB,array_val, &
            buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
            ninterfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,my_neighbours_ext_mesh, &
            request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh &
            )

  implicit none

  include "constants.h"

! include values created by the mesher
  include "OUTPUT_FILES/values_from_mesher.h"

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val

  integer :: NPROC
  integer :: NGLOB_AB

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,ninterfaces_ext_mesh) :: &
       buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh

  integer :: ninterfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  integer, dimension(ninterfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,ninterfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(ninterfaces_ext_mesh) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  integer ipoin,iinterface
  integer sender,receiver

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if(NPROC > 1) then

! partition border copy into the buffer
  do iinterface = 1, ninterfaces_ext_mesh
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      buffer_send_scalar_ext_mesh(ipoin,iinterface) = array_val(ibool_interfaces_ext_mesh(ipoin,iinterface))
    enddo
  enddo

! send messages
  do iinterface = 1, ninterfaces_ext_mesh
    call issend_cr(buffer_send_scalar_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
         nibool_interfaces_ext_mesh(iinterface), &
         my_neighbours_ext_mesh(iinterface), &
         itag, &
         request_send_scalar_ext_mesh(iinterface) &
         )
    call irecv_cr(buffer_recv_scalar_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
         nibool_interfaces_ext_mesh(iinterface), &
         my_neighbours_ext_mesh(iinterface), &
         itag, &
         request_recv_scalar_ext_mesh(iinterface) &
         )
  enddo

! wait for communications completion
  do iinterface = 1, ninterfaces_ext_mesh
    call wait_req(request_recv_scalar_ext_mesh(iinterface))
  enddo

! adding contributions of neighbours
  do iinterface = 1, ninterfaces_ext_mesh
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
           array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_scalar_ext_mesh(ipoin,iinterface)
    enddo
  enddo

  endif

  end subroutine assemble_MPI_scalar_ext_mesh

!
!----
!

  subroutine assemble_MPI_scalar_i_ext_mesh(NPROC,NGLOB_AB,array_val, &
            buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
            ninterfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,my_neighbours_ext_mesh, &
            request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh &
            )

  implicit none

  include "constants.h"

! include values created by the mesher
  include "OUTPUT_FILES/values_from_mesher.h"

! array to assemble
  integer, dimension(NGLOB_AB) :: array_val

  integer :: NPROC
  integer :: NGLOB_AB

  integer, dimension(max_nibool_interfaces_ext_mesh,ninterfaces_ext_mesh) :: &
       buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh

  integer :: ninterfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  integer, dimension(ninterfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,ninterfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(ninterfaces_ext_mesh) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  integer ipoin,iinterface
  integer sender,receiver

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if(NPROC > 1) then

! partition border copy into the buffer
  do iinterface = 1, ninterfaces_ext_mesh
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      buffer_send_scalar_ext_mesh(ipoin,iinterface) = array_val(ibool_interfaces_ext_mesh(ipoin,iinterface))
    enddo
  enddo

! send messages
  do iinterface = 1, ninterfaces_ext_mesh
    call issend_i(buffer_send_scalar_ext_mesh(1,iinterface), &
         nibool_interfaces_ext_mesh(iinterface), &
         my_neighbours_ext_mesh(iinterface), &
         itag, &
         request_send_scalar_ext_mesh(iinterface) &
         )
    call irecv_i(buffer_recv_scalar_ext_mesh(1,iinterface), &
         nibool_interfaces_ext_mesh(iinterface), &
         my_neighbours_ext_mesh(iinterface), &
         itag, &
         request_recv_scalar_ext_mesh(iinterface) &
         )
  enddo

! wait for communications completion
  do iinterface = 1, ninterfaces_ext_mesh
    call wait_req(request_recv_scalar_ext_mesh(iinterface))
  enddo

! adding contributions of neighbours
  do iinterface = 1, ninterfaces_ext_mesh
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
           array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_scalar_ext_mesh(ipoin,iinterface)
    enddo
  enddo

! wait for communications completion (send)
  do iinterface = 1, ninterfaces_ext_mesh
    call wait_req(request_send_scalar_ext_mesh(iinterface))
  enddo


  endif

  end subroutine assemble_MPI_scalar_i_ext_mesh
