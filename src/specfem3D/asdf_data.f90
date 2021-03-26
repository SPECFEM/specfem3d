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

!> Module defining the data structure for ASDF
!!     * Waveforms defined per event
!!     * Waveform attributes defined per seismogram per event

module asdf_data

  implicit none

  private

  type asdf_record
    real, pointer :: record(:)
  end type asdf_record

  type asdf_event
    !scalars
    character(len=13) :: event

    !size info
    integer :: nreceivers
    integer :: nrec_local

    !seismic record info
    real, pointer :: receiver_lat(:), receiver_lo(:)
    real, pointer :: receiver_el(:), receiver_dpt(:)

    character(len=32),pointer :: receiver_name_array(:)
    character(len=8),pointer :: network_array(:)
    character(len=3),pointer :: component_array(:)

    !seismograms
    type (asdf_record), pointer :: records(:)
  end type asdf_event

  ! ASDF
  type(asdf_event) :: asdf_container

  public :: asdf_container

end module asdf_data
