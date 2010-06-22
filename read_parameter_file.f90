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

  subroutine read_parameter_file( NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                        OCEANS,ANISOTROPY,ABSORBING_CONDITIONS, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,NTSTEP_BETWEEN_OUTPUT_INFO, &
                        SIMULATION_TYPE,SAVE_FORWARD )

  implicit none

  include "constants.h"

  integer NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,SIMULATION_TYPE
  integer NSOURCES,NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO,UTM_PROJECTION_ZONE

  double precision DT,HDUR_MOVIE

  logical ATTENUATION,USE_OLSEN_ATTENUATION,OCEANS,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT,USE_HIGHRES_FOR_MOVIES
  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,SUPPRESS_UTM_PROJECTION

  character(len=256) LOCAL_PATH,CMTSOLUTION

! local variables
  integer ::ios,icounter,isource,idummy,nproc_eta_old,nproc_xi_old
  double precision :: hdur,minval_hdur
  character(len=256) :: dummystring
  integer, external :: err_occurred

  ! opens file DATA/Par_file
  call open_parameter_file()

  ! reads in parameters
  call read_value_integer(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if(err_occurred() /= 0) return
  call read_value_logical(SAVE_FORWARD, 'solver.SAVE_FORWARD')
  if(err_occurred() /= 0) return
  call read_value_integer(UTM_PROJECTION_ZONE, 'mesher.UTM_PROJECTION_ZONE')
  if(err_occurred() /= 0) return
  call read_value_logical(SUPPRESS_UTM_PROJECTION, 'mesher.SUPPRESS_UTM_PROJECTION')
  if(err_occurred() /= 0) return
  ! total number of processors 
  call read_value_integer(NPROC, 'mesher.NPROC')
  if(err_occurred() /= 0) then
    ! checks if it's using an old Par_file format
    call read_value_integer(nproc_eta_old, 'mesher.NPROC_ETA')
    if( err_occurred() /= 0 ) then
      print*,'please specify the number of processes in Par_file as:'
      print*,'NPROC           =    <my_number_of_desired_processes> '
      return
    endif
    ! checks if it's using an old Par_file format
    call read_value_integer(nproc_xi_old, 'mesher.NPROC_XI')
    if( err_occurred() /= 0 ) then
      print*,'please specify the number of processes in Par_file as:'
      print*,'NPROC           =    <my_number_of_desired_processes> '
      return
    endif
    NPROC = nproc_eta_old * nproc_xi_old    
  endif  
  call read_value_integer(NSTEP, 'solver.NSTEP')
  if(err_occurred() /= 0) return
  call read_value_double_precision(DT, 'solver.DT')
  if(err_occurred() /= 0) return
  call read_value_logical(OCEANS, 'model.OCEANS')
  if(err_occurred() /= 0) return
  call read_value_logical(ATTENUATION, 'model.ATTENUATION')
  if(err_occurred() /= 0) return
  call read_value_logical(USE_OLSEN_ATTENUATION, 'model.USE_OLSEN_ATTENUATION')
  if(err_occurred() /= 0) return
  call read_value_logical(ANISOTROPY, 'model.ANISOTROPY')
  if(err_occurred() /= 0) return
  call read_value_logical(ABSORBING_CONDITIONS, 'solver.ABSORBING_CONDITIONS')
  if(err_occurred() /= 0) return
  call read_value_logical(MOVIE_SURFACE, 'solver.MOVIE_SURFACE')
  if(err_occurred() /= 0) return
  call read_value_logical(MOVIE_VOLUME, 'solver.MOVIE_VOLUME')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'solver.NTSTEP_BETWEEN_FRAMES')
  if(err_occurred() /= 0) return
  call read_value_logical(CREATE_SHAKEMAP, 'solver.CREATE_SHAKEMAP')
  if(err_occurred() /= 0) return
  call read_value_logical(SAVE_DISPLACEMENT, 'solver.SAVE_DISPLACEMENT')
  if(err_occurred() /= 0) return
  call read_value_logical(USE_HIGHRES_FOR_MOVIES, 'solver.USE_HIGHRES_FOR_MOVIES')
  if(err_occurred() /= 0) return
  call read_value_double_precision(HDUR_MOVIE, 'solver.HDUR_MOVIE')
  if(err_occurred() /= 0) return
  call read_value_logical(SAVE_MESH_FILES, 'mesher.SAVE_MESH_FILES')
  if(err_occurred() /= 0) return
  call read_value_string(LOCAL_PATH, 'LOCAL_PATH')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_INFO, 'solver.NTSTEP_BETWEEN_OUTPUT_INFO')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'solver.NTSTEP_BETWEEN_OUTPUT_SEISMOS')
  if(err_occurred() /= 0) return
  call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'solver.PRINT_SOURCE_TIME_FUNCTION')
  if(err_occurred() /= 0) return


  ! compute the total number of sources in the CMTSOLUTION file
  ! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file
  call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION', 'DATA/CMTSOLUTION')
  
  open(unit=1,file=CMTSOLUTION,iostat=ios,status='old',action='read')
  if(ios /= 0) stop 'error opening CMTSOLUTION file'

  icounter = 0
  do while(ios == 0)
    read(1,"(a)",iostat=ios) dummystring
    if(ios == 0) icounter = icounter + 1
  enddo
  close(1)

  if(mod(icounter,NLINES_PER_CMTSOLUTION_SOURCE) /= 0) &
    stop 'total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'

  NSOURCES = icounter / NLINES_PER_CMTSOLUTION_SOURCE
  if(NSOURCES < 1) stop 'need at least one source in CMTSOLUTION file'

  ! compute the minimum value of hdur in CMTSOLUTION file
  open(unit=1,file=CMTSOLUTION,status='old',action='read')
  minval_hdur = HUGEVAL
  do isource = 1,NSOURCES

    ! skip other information
    do idummy = 1,3
      read(1,"(a)") dummystring
    enddo

    ! read half duration and compute minimum
    read(1,"(a)") dummystring
    read(dummystring(15:len_trim(dummystring)),*) hdur
    minval_hdur = min(minval_hdur,hdur)

    ! skip other information
    do idummy = 1,9
      read(1,"(a)") dummystring
    enddo

  enddo
  close(1)

! one cannot use a Heaviside source for the movies
  if((MOVIE_SURFACE .or. MOVIE_VOLUME) .and. sqrt(minval_hdur**2 + HDUR_MOVIE**2) < TINYVAL) &
    stop 'hdur too small for movie creation, movies do not make sense for Heaviside source'

! close parameter file
  call close_parameter_file()
  
  end subroutine read_parameter_file

