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

  subroutine write_output_ASCII_or_binary(one_seismogram, &
                                          NSTEP,it,SIMULATION_TYPE,DT,t0, &
                                          iorientation,sisname,final_LOCAL_PATH)

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function

  use constants

  use specfem_par, only: myrank,USE_BINARY_FOR_SEISMOGRAMS,SAVE_ALL_SEISMOS_IN_ONE_FILE,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
    seismo_offset,seismo_current,NTSTEP_BETWEEN_OUTPUT_SAMPLE

  implicit none

  integer,intent(in) :: NSTEP,it,SIMULATION_TYPE

  real(kind=CUSTOM_REAL), dimension(NDIM,NTSTEP_BETWEEN_OUTPUT_SEISMOS/NTSTEP_BETWEEN_OUTPUT_SAMPLE),intent(in) :: one_seismogram

  double precision,intent(in) :: t0,DT

  integer,intent(in) :: iorientation

  character(len=MAX_STRING_LEN),intent(in) :: sisname,final_LOCAL_PATH

  ! local parameter
  integer :: isample,ier,it_current
  double precision :: time_t_db
  real(kind=CUSTOM_REAL) :: value,time_t

  real, dimension(1:seismo_current) :: tr

  ! opens seismogram file
  if (USE_BINARY_FOR_SEISMOGRAMS) then

    ! binary format case
    if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
      ! trace name
      write(IOUT) sisname
    else
      if (seismo_offset == 0) then
        open(unit=IOUT, file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)), &
             form='unformatted', action='write', access='stream', status='replace', iostat=ier)
      else
        open(unit=IOUT, file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)), &
             form='unformatted', action='write', access='stream', status='old', position='append', iostat=ier)
      endif
      if (ier /= 0) call exit_mpi(myrank,'Error opening file: '// &
                                  final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)))
    endif

  else

    ! ASCII format case
    if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
      ! trace name
      write(IOUT,*) sisname(1:len_trim(sisname))
    else
      if (seismo_offset == 0) then
        open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)), &
             action='write', status='unknown',iostat=ier)
      else
        open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)), &
             action='write', status='old',position='append',iostat=ier)
      endif
      if (ier /= 0) call exit_mpi(myrank,'Error opening file: '// &
                                  final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)))
    endif

  endif

  ! make sure we never write more than the maximum number of time steps
  ! subtract half duration of the source to make sure travel time is correct
  do isample = 1,seismo_current

    ! seismogram value
    value = one_seismogram(iorientation,isample)

    ! current time increment
    it_current = seismo_offset + isample
    if (it_current > it) stop 'Invalid current time step in writing output ASCII/binary'

    ! time
    if (SIMULATION_TYPE == 1) then
      ! forward simulation
      ! distinguish between single and double precision for reals
      time_t_db = dble( (it_current-1) * NTSTEP_BETWEEN_OUTPUT_SAMPLE) * DT - t0
    else if (SIMULATION_TYPE == 3) then
      ! adjoint simulation: backward/reconstructed wavefields
      ! distinguish between single and double precision for reals
      ! note: compare time_t with time used for source term
      time_t_db = dble(NSTEP - it_current * NTSTEP_BETWEEN_OUTPUT_SAMPLE) * DT - t0
    endif

    ! converts time to CUSTOM_REAL for output
    time_t = real(time_t_db,kind=CUSTOM_REAL)

    if (USE_BINARY_FOR_SEISMOGRAMS) then
      ! binary format case
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
        ! outputs CUSTOM_REAL values
        write(IOUT) time_t,value
      else
        tr(isample) = value
      endif
    else
      ! ASCII format case
      write(IOUT,*) time_t,value
    endif

  enddo

  if (.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) then
    ! binary format case
    ! writes out whole trace into binary file
    if (USE_BINARY_FOR_SEISMOGRAMS) write(IOUT) tr(:)
    ! closes files
    close(IOUT)
  endif

  end subroutine write_output_ASCII_or_binary

