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

  subroutine write_output_ASCII(one_seismogram, &
              NSTEP,it,SIMULATION_TYPE,DT,t0,myrank, &
              iorientation,irecord,sisname,final_LOCAL_PATH)

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function

  implicit none

  include "constants.h"

  integer :: NSTEP,it,SIMULATION_TYPE

  real(kind=CUSTOM_REAL), dimension(NDIM,NSTEP) :: one_seismogram

  double precision t0,DT

  integer myrank
  integer iorientation,irecord

  character(len=256) sisname,final_LOCAL_PATH

  ! local parameter
  integer isample,nt_s,ier
  real, dimension(:), allocatable :: tr
  real(kind=CUSTOM_REAL) :: time_t

  ! opens seismogram file
  if(SEISMOGRAMS_BINARY)then
     ! allocate trace
     nt_s = min(it,NSTEP)
     allocate(tr(nt_s),stat=ier)
     if( ier /= 0 ) stop 'error allocating array tr'

     ! binary format case
     open(unit=IOUT, file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//&
          sisname(1:len_trim(sisname)), form='unformatted', access='direct', recl=4*(nt_s))
     tr(:)=0
  else
     ! ASCII format case
     open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//&
          sisname(1:len_trim(sisname)),status='unknown')
  endif

  ! make sure we never write more than the maximum number of time steps
  ! subtract half duration of the source to make sure travel time is correct
  do isample = 1,min(it,NSTEP)
    if(irecord == 1) then

      ! forward simulation
      if( SIMULATION_TYPE == 1 ) then
        ! distinguish between single and double precision for reals
        if(CUSTOM_REAL == SIZE_REAL) then
          time_t = sngl( dble(isample-1)*DT - t0 )
        else
          time_t = dble(isample-1)*DT - t0
        endif
      endif

      ! adjoint simulation: backward/reconstructed wavefields
      if( SIMULATION_TYPE == 3 ) then
        ! distinguish between single and double precision for reals
        ! note: compare time_t with time used for source term
        if(CUSTOM_REAL == SIZE_REAL) then
          time_t = sngl( dble(NSTEP-isample)*DT - t0 )
        else
          time_t = dble(NSTEP-isample)*DT - t0
        endif
      endif

      if(SEISMOGRAMS_BINARY) then
         ! binary format case
         !tr(isample)=seismograms(iorientation,irec_local,isample)
         tr(isample) = one_seismogram(iorientation,isample)
      else
         ! ASCII format case
         !write(IOUT,*) time_t,' ',seismograms(iorientation,irec_local,isample)
         write(IOUT,*) time_t,' ',one_seismogram(iorientation,isample)
      endif

    else
      call exit_MPI(myrank,'incorrect record label')
    endif
  enddo

  ! binary format case
  if(SEISMOGRAMS_BINARY) then
    ! writes out whole trace into binary file
    write(IOUT,rec=1) tr
    deallocate(tr)
  endif

  close(IOUT)

  end subroutine write_output_ASCII

