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

!==================
module data_pdb
!==================

  use global_parameters, only                : sp, dp
  implicit none
  public

  integer                                   :: nproc
  integer, dimension(:), allocatable        :: nel
  integer, dimension(:), allocatable        :: nel_fluid
  integer, dimension(:), allocatable        :: nel_solid
  integer, dimension(:,:), allocatable      :: procel
  integer, dimension(:), allocatable        :: inv_procel
  integer, dimension(:,:), allocatable      :: procel_fluid, procel_fluidp
  integer, dimension(:,:), allocatable      :: procel_solid, procel_solidp
  integer, dimension(:,:), allocatable      :: inv_procel_solidp
  integer, dimension(:,:), allocatable      :: inv_procel_fluidp
  integer, dimension(:), allocatable        :: nbelong
  integer, dimension(:), allocatable        :: nprocb
  integer, dimension(:,:), allocatable      :: lprocb
  integer, dimension(:), allocatable        :: el2proc 

  ! Glocal message passing...redundant eventually!
  integer, dimension(:), allocatable        :: sizerecvp, sizesendp
  integer, dimension(:,:), allocatable      :: listrecvp, listsendp
  integer, dimension(:,:), allocatable      :: sizemsgrecvp, sizemsgsendp
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_recvp
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_sendp

  ! Slocal - Solid message passing, these all go into the database
  integer, dimension(:), allocatable        :: sizerecvp_solid, sizesendp_solid
  integer, dimension(:,:), allocatable      :: listrecvp_solid, listsendp_solid
  integer, dimension(:,:), allocatable      :: sizemsgrecvp_solid, sizemsgsendp_solid
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_recvp_solid
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_sendp_solid

  ! Flocal - Fluid message passing, these all go into the database
  integer, dimension(:), allocatable        :: sizerecvp_fluid, sizesendp_fluid
  integer, dimension(:,:), allocatable      :: listrecvp_fluid, listsendp_fluid
  integer, dimension(:,:), allocatable      :: sizemsgrecvp_fluid, sizemsgsendp_fluid
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_recvp_fluid
  integer, dimension(:,:,:), allocatable    :: glocal_index_msg_sendp_fluid

  ! Solid-fluid
  integer, dimension(:), allocatable        :: nbelong_solid, nbelong_fluid
  integer, dimension(:), allocatable        :: nprocb_solid, nprocb_fluid
  integer, dimension(:,:), allocatable      :: lprocb_solid, lprocb_fluid

  integer, dimension(:,:), allocatable      :: igloc_solid, igloc_fluid
  integer, dimension(:), allocatable        :: nglobp_solid, nglobp_fluid

  integer, dimension(:), allocatable        :: slob2sloc, flob2floc

  character(len=6), dimension(:,:), allocatable :: eltypep_solid, eltypep_fluid


  ! Solid-fluid boundary
  integer, allocatable                      :: nbdry_el(:)
  integer, allocatable                      :: belemp(:,:)
  integer, allocatable                      :: bdry_solid_elp(:,:), bdry_fluid_elp(:,:)
  integer, dimension(:,:), allocatable      :: bdry_jpol_solidp, bdry_jpol_fluidp
  logical, allocatable                      :: have_bdry_elemp(:)

  ! glocal arrays
  integer, dimension(:,:), allocatable      :: igloc
  integer, dimension(:), allocatable        :: nglobp

  ! Serendipity arrays
  integer, dimension(:), allocatable         :: nglobmeshp
  real(kind=dp), dimension(:,:), allocatable :: scpp, zcpp
  integer, dimension(:,:), allocatable       :: iglobcp
  integer, dimension(:,:,:), allocatable     :: lnodescp
  character(len=6), dimension(:,:), allocatable :: eltypep
  logical, dimension(:,:), allocatable       :: coarsingp

  real(kind=dp), dimension(:), allocatable   :: theta_min_proc, theta_max_proc

  ! global to glocal mapping
  integer, dimension(:,:), allocatable       :: glob2gloc

!======================
end module data_pdb
!======================
