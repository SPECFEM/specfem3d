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


  subroutine read_adjoint_sources_ASDF(adj_source_name, adj_source, index_start, index_end)

  use specfem_par, only: CUSTOM_REAL, myrank, current_asdf_handle
  use iso_c_binding, only: C_NULL_CHAR

  implicit none

  character(len=*), intent(in) :: adj_source_name
  real(kind=CUSTOM_REAL), dimension(*),intent(out) :: adj_source ! NSTEP block size
  integer,intent(in) :: index_start, index_end

  ! local parameters
  integer :: itime, offset, nsamples
  !--- Error variable
  integer ier

  ! Fortran/C index convension. Meaning Fortran starts at 1 and C starts
  ! at 0 so
  ! we need to subtract 1 for the C subroutine read_partial_waveform_f
  offset = index_start - 1 ! the value to start reading from
  nsamples = index_end - index_start + 1 ! this is how many points we want to read in from the adjoint source

  ! print *, myrank, " myrank ", trim(adj_source_name)
  ! print *, " offset, ", offset
  ! print *, " nsamples ", nsamples
  ! print *, adj_source_name, " reading"

  call ASDF_read_partial_waveform_f(current_asdf_handle, &
                                    "AuxiliaryData/AdjointSources/" // trim(adj_source_name) // C_NULL_CHAR, &
                                    offset, nsamples, adj_source, ier)

  if (ier /= 0) then
    print *,'Error reading adjoint source: ',trim(adj_source_name)
    print *,'rank ',myrank,' - time step: ',itime,' index_start:',index_start,' index_end: ',index_end
    print *,'  ',trim(adj_source_name)//'has wrong length, please check with your simulation duration'
    call exit_MPI(myrank,'Adjoint source '//trim(adj_source_name)//' has wrong length, please check with your simulation duration')
  endif

  end subroutine read_adjoint_sources_ASDF

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_adjoint_sources_ASDF(irec, nadj_files_found)

  use specfem_par

  use iso_c_binding, only: C_NULL_CHAR

  implicit none

  integer,intent(in) :: irec
  integer,intent(inout) :: nadj_files_found

  ! local parameters
  integer :: nsamples_infered
  integer :: icomp,ier
  integer :: adjoint_source_exists
  character(len=MAX_STRING_LEN) :: adj_filename,adj_source_file
  character(len=3),dimension(NDIM) :: comp

  adj_source_file = trim(network_name(irec))//'_'//trim(station_name(irec))

  ! prepares channel names
  do icomp = 1,NDIM
    call write_channel_name(icomp,comp(icomp))
  enddo

  ! loops over file components E/N/Z
  do icomp = 1,NDIM

    ! name of adjoint source file for this component
    adj_filename = trim(adj_source_file) // '_'// comp(icomp)

    ! checks if adjoint source exists in ASDF file
    call ASDF_adjoint_source_exists_f(current_asdf_handle, trim(adj_filename) // C_NULL_CHAR, adjoint_source_exists)

    if (adjoint_source_exists == 0) then
      !adjoint source not found
      !stops simulation
      call exit_MPI(myrank,'adjoint source '//trim(adj_filename)//' not found, please check STATIONS_ADJOINT file and ASDF file')
    endif

    ! checks length of file
    call ASDF_get_num_elements_from_path_f(current_asdf_handle, &
                                           "AuxiliaryData/AdjointSources/" // trim(adj_filename) // C_NULL_CHAR, &
                                           nsamples_infered, ier)

    ! checks length
    if (nsamples_infered /= NSTEP) then
      print *,'ASDF adjoint source error: ',"AuxiliaryData/AdjointSources/" // trim(adj_filename), &
              ' has length',nsamples_infered,' but should be',NSTEP
      call exit_MPI(myrank, &
        'file ' // "AuxiliaryData/AdjointSources/" // trim(adj_filename) // &
        ' length is wrong, please check your adjoint sources and your simulation duration')
    endif

    ! updates counter for found files
    nadj_files_found = nadj_files_found + 1

  enddo

  end subroutine check_adjoint_sources_ASDF
