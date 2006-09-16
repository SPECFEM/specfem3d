!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!----
!---- assemble the contributions between slices and chunks using MPI
!----

  subroutine assemble_MPI_vector(array_val,iproc_xi,iproc_eta,addressing, &
            iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta, &
            buffer_send_faces_vector,buffer_received_faces_vector,npoin2D_xi,npoin2D_eta, &
            NPROC_XI,NPROC_ETA,NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NPOIN2DMAX_XY)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

! include values created by the mesher
  include "OUTPUT_FILES/values_from_mesher.h"

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val

  integer iproc_xi,iproc_eta
  integer npoin2D_xi,npoin2D_eta

  integer NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NPOIN2DMAX_XY
  integer NPROC_XI,NPROC_ETA

! for addressing of the slices
  integer, dimension(0:NPROC_XI-1,0:NPROC_ETA-1) :: addressing

! 2-D addressing and buffers for summation between slices
  integer, dimension(NPOIN2DMAX_XMIN_XMAX) :: iboolleft_xi,iboolright_xi
  integer, dimension(NPOIN2DMAX_YMIN_YMAX) :: iboolleft_eta,iboolright_eta

  real(kind=CUSTOM_REAL), dimension(NDIM,NPOIN2DMAX_XY) :: buffer_send_faces_vector,buffer_received_faces_vector

! MPI status of messages to be received
  integer msg_status(MPI_STATUS_SIZE)

  integer ipoin
  integer sender,receiver,ier

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! here we have to assemble all the contributions between slices using MPI

!----
!---- first assemble along xi using the 2-D topology
!----

! assemble along xi only if more than one slice
  if(NPROC_XI > 1) then

! slices copy the right face into the buffer
  do ipoin=1,npoin2D_xi
    buffer_send_faces_vector(:,ipoin) = array_val(:,iboolright_xi(ipoin))
  enddo

! send messages forward along each row
  if(iproc_xi == 0) then
    sender = MPI_PROC_NULL
  else
    sender = addressing(iproc_xi - 1,iproc_eta)
  endif
  if(iproc_xi == NPROC_XI-1) then
    receiver = MPI_PROC_NULL
  else
    receiver = addressing(iproc_xi + 1,iproc_eta)
  endif
  call MPI_SENDRECV(buffer_send_faces_vector,NDIM*npoin2D_xi,CUSTOM_MPI_TYPE,receiver, &
        itag2,buffer_received_faces_vector,NDIM*npoin2D_xi,CUSTOM_MPI_TYPE,sender, &
        itag,MPI_COMM_WORLD,msg_status,ier)

! all slices add the buffer received to the contributions on the left face
  if(iproc_xi > 0) then
  do ipoin=1,npoin2D_xi
      array_val(:,iboolleft_xi(ipoin)) = array_val(:,iboolleft_xi(ipoin)) + &
                                buffer_received_faces_vector(:,ipoin)
  enddo
  endif

! the contributions are correctly assembled on the left side of each slice
! now we have to send the result back to the sender
! all slices copy the left face into the buffer
  do ipoin=1,npoin2D_xi
    buffer_send_faces_vector(:,ipoin) = array_val(:,iboolleft_xi(ipoin))
  enddo

! send messages backward along each row
  if(iproc_xi == NPROC_XI-1) then
    sender = MPI_PROC_NULL
  else
    sender = addressing(iproc_xi + 1,iproc_eta)
  endif
  if(iproc_xi == 0) then
    receiver = MPI_PROC_NULL
  else
    receiver = addressing(iproc_xi - 1,iproc_eta)
  endif
  call MPI_SENDRECV(buffer_send_faces_vector,NDIM*npoin2D_xi,CUSTOM_MPI_TYPE,receiver, &
        itag2,buffer_received_faces_vector,NDIM*npoin2D_xi,CUSTOM_MPI_TYPE,sender, &
        itag,MPI_COMM_WORLD,msg_status,ier)

! all slices copy the buffer received to the contributions on the right face
  if(iproc_xi < NPROC_XI-1) then
  do ipoin=1,npoin2D_xi
    array_val(:,iboolright_xi(ipoin)) = buffer_received_faces_vector(:,ipoin)
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
    buffer_send_faces_vector(:,ipoin) = array_val(:,iboolright_eta(ipoin))
  enddo

! send messages forward along each row
  if(iproc_eta == 0) then
    sender = MPI_PROC_NULL
  else
    sender = addressing(iproc_xi,iproc_eta - 1)
  endif
  if(iproc_eta == NPROC_ETA-1) then
    receiver = MPI_PROC_NULL
  else
    receiver = addressing(iproc_xi,iproc_eta + 1)
  endif
  call MPI_SENDRECV(buffer_send_faces_vector,NDIM*npoin2D_eta,CUSTOM_MPI_TYPE,receiver, &
    itag2,buffer_received_faces_vector,NDIM*npoin2D_eta,CUSTOM_MPI_TYPE,sender, &
    itag,MPI_COMM_WORLD,msg_status,ier)

! all slices add the buffer received to the contributions on the left face
  if(iproc_eta > 0) then
  do ipoin=1,npoin2D_eta
      array_val(:,iboolleft_eta(ipoin)) = array_val(:,iboolleft_eta(ipoin)) + &
                                buffer_received_faces_vector(:,ipoin)
  enddo
  endif

! the contributions are correctly assembled on the left side of each slice
! now we have to send the result back to the sender
! all slices copy the left face into the buffer
  do ipoin=1,npoin2D_eta
    buffer_send_faces_vector(:,ipoin) = array_val(:,iboolleft_eta(ipoin))
  enddo

! send messages backward along each row
  if(iproc_eta == NPROC_ETA-1) then
    sender = MPI_PROC_NULL
  else
    sender = addressing(iproc_xi,iproc_eta + 1)
  endif
  if(iproc_eta == 0) then
    receiver = MPI_PROC_NULL
  else
    receiver = addressing(iproc_xi,iproc_eta - 1)
  endif
  call MPI_SENDRECV(buffer_send_faces_vector,NDIM*npoin2D_eta,CUSTOM_MPI_TYPE,receiver, &
    itag2,buffer_received_faces_vector,NDIM*npoin2D_eta,CUSTOM_MPI_TYPE,sender, &
    itag,MPI_COMM_WORLD,msg_status,ier)

! all slices copy the buffer received to the contributions on the right face
  if(iproc_eta < NPROC_ETA-1) then
  do ipoin=1,npoin2D_eta
    array_val(:,iboolright_eta(ipoin)) = buffer_received_faces_vector(:,ipoin)
  enddo
  endif

  endif

  end subroutine assemble_MPI_vector

