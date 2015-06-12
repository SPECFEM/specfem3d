!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!===================
!> Message-passing communication variables
!! Note: For the easy serialization of the code (see commun.f90 and commpi.f90),
!! one still needs this module as these quantities are read in from the mesher.
module data_comm
!===================

  use global_parameters
  use linked_list
  
  implicit none
  public 

  integer                              :: comm
  integer                              :: mpi_realkind

  ! Solid domain message passing

  integer                              :: num_comm_gll_solid
  integer                              :: sizemsgrecvmax_solid
  integer                              :: sizemsgsendmax_solid
  integer                              :: sizemsgmax_solid
  integer, allocatable                 :: glob2el_solid(:,:)
  real(kind=realkind), allocatable     :: buffs_solid(:,:), buffr_solid(:,:)
  type(list)                           :: buffs_all_solid, buffr_all_solid

  class(link), pointer                 :: buffs, buffr

  integer                              :: sizerecv_solid, sizesend_solid
  integer, dimension(:),   allocatable :: listrecv_solid, sizemsgrecv_solid
  integer, dimension(:),   allocatable :: listsend_solid, sizemsgsend_solid
  integer, dimension(:,:), allocatable :: glocal_index_msg_recv_solid
  integer, dimension(:,:), allocatable :: glocal_index_msg_send_solid

  integer, dimension(:), allocatable   :: recv_request_solid, send_request_solid

  ! Fluid domain message passing
  integer                              :: num_comm_gll_fluid
  integer                              :: sizemsgrecvmax_fluid
  integer                              :: sizemsgsendmax_fluid
  integer                              :: sizemsgmax_fluid
  integer, allocatable                 :: glob2el_fluid(:,:)
  real(kind=realkind), allocatable     :: buffs_fluid(:,:), buffr_fluid(:,:)
  type(list)                           :: buffs_all_fluid, buffr_all_fluid

  integer                              :: sizerecv_fluid, sizesend_fluid
  integer, dimension(:),   allocatable :: listrecv_fluid, sizemsgrecv_fluid
  integer, dimension(:),   allocatable :: listsend_fluid, sizemsgsend_fluid
  integer, dimension(:,:), allocatable :: glocal_index_msg_recv_fluid
  integer, dimension(:,:), allocatable :: glocal_index_msg_send_fluid

  integer, dimension(:), allocatable   :: recv_request_fluid, send_request_fluid

!=======================
end module data_comm
!=======================
