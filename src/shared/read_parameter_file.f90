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

  subroutine read_parameter_file(NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT,NGNOD,NGNOD2D, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,TOMOGRAPHY_PATH, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                        APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,ANISOTROPY,STACEY_ABSORBING_CONDITIONS,MOVIE_TYPE, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,NTSTEP_BETWEEN_OUTPUT_INFO, &
                        SIMULATION_TYPE,SAVE_FORWARD,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                        USE_FORCE_POINT_SOURCE,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        USE_RICKER_TIME_FUNCTION,OLSEN_ATTENUATION_RATIO,PML_CONDITIONS, &
                        PML_INSTEAD_OF_FREE_SURFACE,f0_FOR_PML,IMODEL,FULL_ATTENUATION_SOLID,TRAC_PATH)

  implicit none

  include "constants.h"

  integer NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,SIMULATION_TYPE, NTSTEP_BETWEEN_READ_ADJSRC
  integer NSOURCES,NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO,UTM_PROJECTION_ZONE
  integer NOISE_TOMOGRAPHY,NGNOD,NGNOD2D,MOVIE_TYPE
  integer IMODEL

  double precision DT,HDUR_MOVIE,OLSEN_ATTENUATION_RATIO,f0_FOR_PML

  logical ATTENUATION,USE_OLSEN_ATTENUATION,APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,STACEY_ABSORBING_CONDITIONS,SAVE_FORWARD
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT,USE_HIGHRES_FOR_MOVIES
  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,SUPPRESS_UTM_PROJECTION
  logical USE_FORCE_POINT_SOURCE,STACEY_INSTEAD_OF_FREE_SURFACE,USE_RICKER_TIME_FUNCTION
  logical PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE,FULL_ATTENUATION_SOLID

  character(len=256) LOCAL_PATH,TOMOGRAPHY_PATH,CMTSOLUTION,FORCESOLUTION,TRAC_PATH ! VM VM adds TRAC_PATH

! local variables
  integer ::ios,icounter,isource,idummy,nproc_eta_old,nproc_xi_old
  double precision :: hdur,minval_hdur
  character(len=256) :: dummystring
  integer, external :: err_occurred

  character(len=150) MODEL
  integer :: i,irange

  ! opens file Par_file
  call open_parameter_file()

  ! reads in parameters
  call read_value_integer(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if(err_occurred() /= 0) return
  call read_value_integer(NOISE_TOMOGRAPHY, 'solver.NOISE_TOMOGRAPHY')
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

  ! number of nodes for 2D and 3D shape functions for quadrilaterals and hexahedra
  call read_value_integer(NGNOD, 'solver.NGNOD')
  if(err_occurred() /= 0) return

  ! define the velocity model
  call read_value_string(MODEL, 'model.MODEL')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MODEL'

  call read_value_logical(APPROXIMATE_OCEAN_LOAD, 'model.APPROXIMATE_OCEAN_LOAD')
  if(err_occurred() /= 0) return
  call read_value_logical(TOPOGRAPHY, 'model.TOPOGRAPHY')
  if(err_occurred() /= 0) return
  call read_value_logical(ATTENUATION, 'model.ATTENUATION')
  if(err_occurred() /= 0) return
  call read_value_logical(FULL_ATTENUATION_SOLID, 'model.FULL_ATTENUATION_SOLID')
  if(err_occurred() /= 0) return
  call read_value_logical(ANISOTROPY, 'model.ANISOTROPY')
  if(err_occurred() /= 0) return
  call read_value_string(TOMOGRAPHY_PATH, 'TOMOGRAPHY_PATH')
  if(err_occurred() /= 0) return
  call read_value_logical(USE_OLSEN_ATTENUATION, 'model.USE_OLSEN_ATTENUATION')
  if(err_occurred() /= 0) return
  call read_value_double_precision(OLSEN_ATTENUATION_RATIO, 'model.OLSEN_ATTENUATION_RATIO')
  if(err_occurred() /= 0) return
  call read_value_logical(PML_CONDITIONS, 'solver.PML_CONDITIONS')
  if(err_occurred() /= 0) return
  call read_value_logical(PML_INSTEAD_OF_FREE_SURFACE, 'model.PML_INSTEAD_OF_FREE_SURFACE')
  if(err_occurred() /= 0) return
  call read_value_double_precision(f0_FOR_PML, 'model.f0_FOR_PML')
  if(err_occurred() /= 0) return
  !call read_value_logical(ROTATE_PML_ACTIVATE, 'solver.ROTATE_PML_ACTIVATE')
  !if(err_occurred() /= 0) return
  !call read_value_double_precision(ROTATE_PML_ANGLE, 'solver.ROTATE_PML_ANGLE')
  !if(err_occurred() /= 0) return
  call read_value_logical(STACEY_ABSORBING_CONDITIONS, 'solver.STACEY_ABSORBING_CONDITIONS')
  if(err_occurred() /= 0) return
  call read_value_logical(STACEY_INSTEAD_OF_FREE_SURFACE, 'model.STACEY_INSTEAD_OF_FREE_SURFACE')
  if(err_occurred() /= 0) return
  call read_value_logical(CREATE_SHAKEMAP, 'solver.CREATE_SHAKEMAP')
  if(err_occurred() /= 0) return
  call read_value_logical(MOVIE_SURFACE, 'solver.MOVIE_SURFACE')
  if(err_occurred() /= 0) return
  call read_value_integer(MOVIE_TYPE, 'solver.MOVIE_TYPE')
  if(err_occurred() /= 0) return
  call read_value_logical(MOVIE_VOLUME, 'solver.MOVIE_VOLUME')
  if(err_occurred() /= 0) return
  call read_value_logical(SAVE_DISPLACEMENT, 'solver.SAVE_DISPLACEMENT')
  if(err_occurred() /= 0) return
  call read_value_logical(USE_HIGHRES_FOR_MOVIES, 'solver.USE_HIGHRES_FOR_MOVIES')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'solver.NTSTEP_BETWEEN_FRAMES')
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
  call read_value_integer(NTSTEP_BETWEEN_READ_ADJSRC, 'solver.NTSTEP_BETWEEN_READ_ADJSRC')
  if(err_occurred() /= 0) return
  call read_value_logical(USE_FORCE_POINT_SOURCE, 'solver.USE_FORCE_POINT_SOURCE')
  if(err_occurred() /= 0) return
  call read_value_logical(USE_RICKER_TIME_FUNCTION, 'solver.USE_RICKER_TIME_FUNCTION')
  if(err_occurred() /= 0) return
  call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'solver.PRINT_SOURCE_TIME_FUNCTION')
  if(err_occurred() /= 0) return
  !! VM VM read the traction path directory
  if (OLD_TEST_TO_FIX_ONE_DAY) then
     call read_value_string(TRAC_PATH, 'TRAC_PATH')
     if(err_occurred() /= 0) return
  endif

  ! close parameter file
  call close_parameter_file()

  ! noise simulations:
  ! double the number of time steps, if running noise simulations (+/- branches)
  if( NOISE_TOMOGRAPHY /= 0 )   NSTEP = 2*NSTEP-1

  ! for noise simulations, we need to save movies at the surface (where the noise is generated)
  ! and thus we force MOVIE_SURFACE to be .true., in order to use variables defined for surface movies later
  if( NOISE_TOMOGRAPHY /= 0 ) then
     MOVIE_TYPE = 1
     MOVIE_SURFACE = .true.
     USE_HIGHRES_FOR_MOVIES = .true.     ! we need to save surface movie everywhere, i.e. at all GLL points on the surface
  endif

  ! the default value of NTSTEP_BETWEEN_READ_ADJSRC (0) is to read the whole trace at the same time
  if( NTSTEP_BETWEEN_READ_ADJSRC == 0 )  NTSTEP_BETWEEN_READ_ADJSRC = NSTEP

  ! total times steps must be dividable by adjoint source chunks/blocks
  if ( mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) /= 0 ) then
     print*,'When NOISE_TOMOGRAPHY is not equal to zero, ACTUAL_NSTEP=2*NSTEP-1'
     stop 'error: mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) must be zero! Please modify Par_file and recompile solver'
  endif

  ! checks number of nodes for 2D and 3D shape functions for quadrilaterals and hexahedra
  ! curvature (i.e. HEX27 elements) is not handled by our internal mesher, for that use Gmsh (CUBIT does not handle it either)
  if( NGNOD == 8 ) then
     NGNOD2D = 4
  else if( NGNOD == 27 ) then
     NGNOD2D = 9
  else if( NGNOD /= 8 .and. NGNOD /= 27 ) then
     stop 'elements should have 8 or 27 control nodes, please modify NGNOD in Par_file and recompile solver'
  endif

  if( USE_FORCE_POINT_SOURCE ) then
     ! compute the total number of sources in the FORCESOLUTION file
     ! there are NLINES_PER_FORCESOLUTION_SOURCE lines per source in that file
     call get_value_string(FORCESOLUTION, 'solver.FORCESOLUTION',&
          IN_DATA_FILES_PATH(1:len_trim(IN_DATA_FILES_PATH))//'FORCESOLUTION')

     open(unit=21,file=trim(FORCESOLUTION),iostat=ios,status='old',action='read')
     if(ios /= 0) stop 'error opening FORCESOLUTION file'

     icounter = 0
     do while(ios == 0)
        read(21,"(a)",iostat=ios) dummystring
        if(ios == 0) icounter = icounter + 1
     enddo
     close(21)

     if( mod(icounter,NLINES_PER_FORCESOLUTION_SOURCE) /= 0 ) &
          stop 'error: total number of lines in FORCESOLUTION file should be a multiple of NLINES_PER_FORCESOLUTION_SOURCE'

     NSOURCES = icounter / NLINES_PER_FORCESOLUTION_SOURCE
     if(NSOURCES < 1) stop 'error: need at least one source in FORCESOLUTION file'

  else
     ! compute the total number of sources in the CMTSOLUTION file
     ! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file
     call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION',&
          IN_DATA_FILES_PATH(1:len_trim(IN_DATA_FILES_PATH))//'CMTSOLUTION')

     open(unit=21,file=trim(CMTSOLUTION),iostat=ios,status='old',action='read')
     if(ios /= 0) stop 'error opening CMTSOLUTION file'

     icounter = 0
     do while(ios == 0)
        read(21,"(a)",iostat=ios) dummystring
        if(ios == 0) icounter = icounter + 1
     enddo
     close(21)

     if(mod(icounter,NLINES_PER_CMTSOLUTION_SOURCE) /= 0) &
          stop 'error: total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'

     NSOURCES = icounter / NLINES_PER_CMTSOLUTION_SOURCE
     if(NSOURCES < 1) stop 'error: need at least one source in CMTSOLUTION file'

     ! compute the minimum value of hdur in CMTSOLUTION file
     open(unit=21,file=trim(CMTSOLUTION),status='old',action='read')
     minval_hdur = HUGEVAL
     do isource = 1,NSOURCES

        ! skip other information
        do idummy = 1,3
           read(21,"(a)") dummystring
        enddo

        ! read half duration and compute minimum
        read(21,"(a)") dummystring
        read(dummystring(15:len_trim(dummystring)),*) hdur
        minval_hdur = min(minval_hdur,hdur)

        ! skip other information
        do idummy = 1,9
           read(21,"(a)") dummystring
        enddo

     enddo
     close(21)

     ! one cannot use a Heaviside source for the movies
     if( (MOVIE_SURFACE .or. MOVIE_VOLUME) .and. sqrt(minval_hdur**2 + HDUR_MOVIE**2) < TINYVAL ) &
          stop 'error: hdur too small for movie creation, movies do not make sense for Heaviside source'
  endif

  ! converts all string characters to lowercase
  irange = iachar('a') - iachar('A')
  do i = 1,len_trim(MODEL)
    if( lge(MODEL(i:i),'A') .and. lle(MODEL(i:i),'Z') ) then
      MODEL(i:i) = achar( iachar(MODEL(i:i)) + irange )
    endif
  enddo

  ! determines velocity model
  select case( trim(MODEL) )

  ! default mesh model
  case( 'default' )
    IMODEL = IMODEL_DEFAULT

  ! 1-D models
  case( '1d_prem' )
    IMODEL = IMODEL_1D_PREM
  case( '1d_socal' )
    IMODEL = IMODEL_1D_SOCAL
  case( '1d_cascadia')
    IMODEL = IMODEL_1D_CASCADIA

  ! user models
  case( '1d_prem_pb' )
    IMODEL = IMODEL_1D_PREM_PB
  case( 'aniso' )
    IMODEL = IMODEL_DEFAULT
    ANISOTROPY = .true.
  case( 'external' )
    IMODEL = IMODEL_USER_EXTERNAL
  case( 'ipati' )
    IMODEL = IMODEL_IPATI
  case( 'ipati_water' )
    IMODEL = IMODEL_IPATI_WATER
  case( 'gll' )
    IMODEL = IMODEL_GLL
  case( 'salton_trough')
    IMODEL = IMODEL_SALTON_TROUGH
  case( 'tomo' )
    IMODEL = IMODEL_TOMO

  case default
    print*
    print*,'********** model not recognized: ',trim(MODEL),' **************'
    print*,'********** using model: default',' **************'
    print*
    IMODEL = IMODEL_DEFAULT
  end select

  ! check
  if( IMODEL == IMODEL_IPATI .or. IMODEL == IMODEL_IPATI_WATER ) then
    if( USE_RICKER_TIME_FUNCTION .eqv. .false. ) &
         stop 'error: please set USE_RICKER_TIME_FUNCTION to .true. in Par_file and recompile solver'
  endif


  ! absorbing conditions
  ! stacey absorbing free surface must have also stacey turned on
  if( STACEY_INSTEAD_OF_FREE_SURFACE ) STACEY_ABSORBING_CONDITIONS = .true.

  ! only one absorbing condition is possible
  if( PML_CONDITIONS ) then
    if( STACEY_ABSORBING_CONDITIONS )&
      stop 'error: please set STACEY_ABSORBING_CONDITIONS and STACEY_INSTEAD_OF_FREE_SURFACE to .false. in Par_file for PML'
  endif

  end subroutine read_parameter_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_gpu_mode(GPU_MODE,GRAVITY)

  implicit none
  include "constants.h"

  logical :: GPU_MODE
  logical :: GRAVITY

  ! initializes flags
  GPU_MODE = .false.
  GRAVITY = .false.

  ! opens file Par_file
  call open_parameter_file()

  call read_value_logical(GPU_MODE, 'solver.GPU_MODE')
  call read_value_logical(GRAVITY, 'solver.GRAVITY')

  ! close parameter file
  call close_parameter_file()

  end subroutine read_gpu_mode

!===============================================================================
!> Read ADIOS related flags from the Par_file
!! \param ADIOS_ENABLED Main flag to decide if ADIOS is used. If setted to
!!                      false no other parameter is taken into account.
!! \param ADIOS_FOR_DATABASES Flag to indicate if the databases are written
!!                            and read with the help of ADIOS.
!! \param ADIOS_FOR_MESH flag to indicate if the mesh (generate database) is
!!                       written using ADIOS.
!! \param ADIOS_FOR_FORWARD_ARRAYS flag to indicate if the solver forward arrays
!!                                 are written using ADIOS.
!! \param ADIOS_FOR_KERNELS flag to indicate if the kernels are saved using
!!                          adios
!! \author MPBL
subroutine read_adios_parameters(ADIOS_ENABLED, ADIOS_FOR_DATABASES,       &
                                 ADIOS_FOR_MESH, ADIOS_FOR_FORWARD_ARRAYS, &
                                 ADIOS_FOR_KERNELS)
  implicit none
  include "constants.h"

  logical, intent(out) :: ADIOS_ENABLED, ADIOS_FOR_DATABASES,       &
                          ADIOS_FOR_MESH, ADIOS_FOR_FORWARD_ARRAYS, &
                          ADIOS_FOR_KERNELS

  ! initialize flags to false
  ADIOS_ENABLED            = .false.
  ADIOS_FOR_DATABASES      = .false.
  ADIOS_FOR_MESH           = .false.
  ADIOS_FOR_FORWARD_ARRAYS = .false.
  ADIOS_FOR_KERNELS        = .false.
  ! opens file Par_file
  call open_parameter_file()
  call read_value_logical(ADIOS_ENABLED, 'solver.ADIOS_ENABLED')
  if (ADIOS_ENABLED) then
    call read_value_logical(ADIOS_FOR_DATABASES, 'solver.ADIOS_FOR_DATABASES')
    call read_value_logical(ADIOS_FOR_MESH, 'solver.ADIOS_FOR_MESH')
    call read_value_logical(ADIOS_FOR_FORWARD_ARRAYS, &
                           'solver.ADIOS_FOR_FORWARD_ARRAYS')
    call read_value_logical(ADIOS_FOR_KERNELS, 'solver.ADIOS_FOR_KERNELS')
  endif
  call close_parameter_file()

end subroutine read_adios_parameters


