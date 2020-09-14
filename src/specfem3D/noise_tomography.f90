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

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file is first implemented for SPECFEM3D_GLOBE, and therefore it follows variables in GLOBAL package.
! Do NOT worry about strange names containing, cause they are from GLOBAL package.
! However, please be informed that some of those varibles (even input parameters) can be dummy (not used),
!    since they may exist only in GLOBAL package.
! Also, slight modifications are included due to different structures in SPECFEM3D_GLOBE and SPECFEM3D.
! In general, this file in SPECFEM3D should be better than the one in SPECFEM3D_GLOBE,
!    since it is less dependent on some other files (values_from_mesher.h)
! Should you be interested in details, please refer to how those subroutines are called in other files
!    ("grep NOISE_TOMOGRAPHY src/*/*90" should return you where those subroutines are called)
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module user_noise_distribution

!daniel: TODO -- setting USE_PIERO_DISTRIBUTION = .true. will produce errors
!            when using with the default example in "example/noise_tomography/"
!            i left it here so that Max can run his example without changing this every time...
  logical,parameter :: USE_PIERO_DISTRIBUTION = .false.

contains

! wrapper function
! this subroutine must be modified by USERS for their own noise distribution

  subroutine noise_distribution_direction(xcoord_in,ycoord_in,zcoord_in, &
                                          normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                                          mask_noise_out)

  use constants

  implicit none

  ! input parameters
  real(kind=CUSTOM_REAL) :: xcoord_in,ycoord_in,zcoord_in
  ! output parameters
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out

  ! Setup for NOISE_TOMOGRAPHY by Piero Basini
  if (USE_PIERO_DISTRIBUTION) then
    call noise_distribution_dir_non_uni(xcoord_in,ycoord_in,zcoord_in, &
                                        normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                                        mask_noise_out)
  else
    ! DEFAULT routine
    call noise_distribution_direction_d(xcoord_in,ycoord_in,zcoord_in, &
                                        normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                                        mask_noise_out)
  endif

  end subroutine noise_distribution_direction

  !
  !-----------------------------------------------------------------------------------------------
  !

! characterizes noise statistics:
!     for a given point (xcoord,ycoord,zcoord), specify the noise direction "normal_x/y/z_noise"
!     and noise distribution "mask_noise"
!
! USERS: need to modify this subroutine for their own noise characteristics
  subroutine noise_distribution_direction_d(xcoord_in,ycoord_in,zcoord_in, &
                                            normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                                            mask_noise_out)
  use constants

  implicit none

  ! input parameters
  real(kind=CUSTOM_REAL) :: xcoord_in,ycoord_in,zcoord_in
  ! output parameters
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  ! local parameters
  real(kind=CUSTOM_REAL) :: ldummy

  !*****************************************************************************************************************
  !******************************** change your noise characteristics below ****************************************
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! noise direction
  !!!!! here, the noise is assumed to be vertical
  normal_x_noise_out = 0.0
  normal_y_noise_out = 0.0
  normal_z_noise_out = 1.0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  noise distribution
  !!!!! here, the noise is assumed to be uniform
  mask_noise_out = 1.0
  !******************************** change your noise characteristics above ****************************************
  !*****************************************************************************************************************

  ! dummy assign to avoid compiler warnings
  ldummy = xcoord_in
  ldummy = ycoord_in
  ldummy = zcoord_in

  end subroutine noise_distribution_direction_d

  !
  !-----------------------------------------------------------------------------------------------
  !

  subroutine noise_distribution_dir_non_uni(xcoord_in,ycoord_in,zcoord_in, &
                                            normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                                            mask_noise_out)

  use constants

  implicit none

  ! input parameters
  real(kind=CUSTOM_REAL) :: xcoord_in,ycoord_in,zcoord_in
  ! output parameters
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  ! local parameters
  !PB VARIABLES TO DEFINE THE REGION OF NOISE
  real(kind=CUSTOM_REAL) :: xcoord,ycoord,zcoord !,xcoord_center,ycoord_center
  real :: lon,lat,colat,lon_cn,lat_cn,dsigma,d,dmax

  ! coordinates "x/y/zcoord_in" actually contain r theta phi, therefore convert back to x y z
  ! call rthetaphi_2_xyz(xcoord,ycoord,zcoord, xcoord_in,ycoord_in,zcoord_in)
  xcoord = xcoord_in
  ycoord = ycoord_in
  zcoord = zcoord_in

  !PB NOT UNIF DISTRIBUTION OF NOISE ON THE SURFACE OF A SPHERE
  !PB lon lat colat ARE IN RADIANS SINCE ARE OBTAINED FROM Cartesian COORDINATES
  !PB lon_cn lat_cn (cn = CENTER OF NOISE REGION) IF NOT, MUST BE CONVERTED IN RADIANS
  !PB lon_cn lat_cn are inserted directly here for simplicity

  lon_cn = (3.89)*PI/180
  lat_cn = (45.113)*PI/180

  if (abs(xcoord) > 1.e-24 .or. abs(ycoord) > 1.e-24) then
    if (xcoord >= 0) then
      lon = asin(ycoord/(sqrt(xcoord**2+ycoord**2)))
    else
      lon = (PI-(asin(ycoord/(sqrt(xcoord**2+ycoord**2)))))
    endif
  else
    lon = 0.0
  endif

  if (abs(zcoord) > 1.e-24) then
    colat = atan(sqrt(xcoord**2+ycoord**2)/zcoord)
  else
    colat = 0.0
  endif
  lat = (PI/2)-colat

  !PB CALCULATE THE DISTANCE BETWEEN CENTER OF NOISE REGION AND EACH
  ! POINT OF THE MODEL'S FREE SURFACE  !PB dsigma IS THE "3D" ANGLE BETWEEN
  ! THE TWO POINTS, THEN d = R*dsigma
  dsigma = acos(sin(lon)*sin(lon_cn)+cos(lon)*cos(lon_cn)*cos(lat-lat_cn))
  d = sqrt(xcoord**2+ycoord**2+zcoord**2)*dsigma

  !PB IF YOU WANT TO USE A NONUNIFORM DISTRIBUTION OF NOISE IN THE EXAMPLE
  !PROVIDED WITH THE CODE, THEN UNCOMMENT THE FOLLOWING LINES (before definition
  !of dmax)

  !  xcoord_center = 30000
  !  ycoord_center = 30000
  !  d = sqrt((xcoord_center-xcoord)**2+(ycoord_center-ycoord)**2)

  !PB NOTE THAT d IS EXPRESSED IN METERS REMEBER THAT WHEN YOU SET THE PARAMETER dmax
  !PB dmax IS THE RADIUS OF THE AREA IN WHICH masc_noise_out IS 1 (NOISE IS DEFINED)

  dmax = 300000

  ! NOTE that all coordinates are non-dimensionalized in GLOBAL package!
  ! USERS are free to choose which set to use,
  ! either "r theta phi" (xcoord_in,ycoord_in,zcoord_in)
  ! or     "x y z"       (xcoord,ycoord,zcoord)


  !*****************************************************************************************************************
  !******************************** change your noise characteristics below ****************************************
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! noise direction
  !!!!! here, the noise is assumed to be vertical
  normal_x_noise_out = 0.0
  normal_y_noise_out = 0.0
  normal_z_noise_out = 1.0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  noise distribution
  !!!!! here, the noise is assumed to be uniform
  !  mask_noise_out = 1.0

  !HERE IS NOT UNIFORM
  if (d <= dmax) then
    mask_noise_out = 1.0
  else
    mask_noise_out = 0.0
  endif

  !******************************** change your noise characteristics above ****************************************
  !*****************************************************************************************************************

  end subroutine noise_distribution_dir_non_uni

end module user_noise_distribution

!
! =============================================================================================================
!

! read parameters
  subroutine read_parameters_noise(nrec,NSTEP,nmovie_points, &
                                   islice_selected_rec,xi_receiver,eta_receiver,gamma_receiver,nu_rec, &
                                   noise_sourcearray,xigll,yigll,zigll, &
                                   ibool, &
                                   xstore,ystore,zstore, &
                                   irec_main_noise,normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                                   nspec,nglob, &
                                   num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                                   ispec_is_acoustic)

  use constants
  use user_noise_distribution

  implicit none

  ! input parameters
  integer,intent(in) :: nrec, NSTEP, nmovie_points
  integer,intent(in) :: nspec,nglob

  integer, dimension(nrec),intent(in) :: islice_selected_rec

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  double precision, dimension(nrec),intent(in)  :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll
  double precision, dimension(NDIM,NDIM,nrec),intent(inout) :: nu_rec
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore,ystore,zstore

  integer,intent(in) :: num_free_surface_faces
  integer, dimension(num_free_surface_faces),intent(in) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces),intent(in) :: free_surface_ijk

  logical, dimension(nspec),intent(in) :: ispec_is_acoustic

  ! output parameters
  integer,intent(inout) :: irec_main_noise
  real(kind=CUSTOM_REAL),intent(inout) :: noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP)
  real(kind=CUSTOM_REAL), dimension(nmovie_points),intent(inout) :: normal_x_noise,normal_y_noise,normal_z_noise,mask_noise

  ! local parameters
  integer :: ipoin,ispec,i,j,k,iglob,ier,iface,igll
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  character(len=MAX_STRING_LEN) :: filename

  ! read main receiver ID -- the ID in "STATIONS"
  filename = trim(OUTPUT_FILES)//'/..//NOISE_TOMOGRAPHY/irec_main_noise'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0) &
    call exit_MPI(myrank, 'file '//trim(filename)//' does NOT exist! This file contains the ID of the main receiver')

  read(IIN_NOISE,*,iostat=ier) irec_main_noise
  if (ier /= 0) call exit_MPI(myrank,'error reading file irec_main_noise')

  close(IIN_NOISE)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  main station: ID = ',irec_main_noise
    call flush_IMAIN()
  endif

  ! checks value
  if (irec_main_noise <= 0 .or. irec_main_noise > nrec) then
    print *,'Error: irec_main_noise value:',irec_main_noise,'must be positive and less than ',nrec
    call exit_MPI(myrank,'error irec_main_noise value')
  endif

  ! writes out main as file info
  if (myrank == 0) then
    open(unit=IOUT_NOISE,file=trim(OUTPUT_FILES)//'/irec_main_noise', &
            status='unknown',action='write',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'error opening file '//trim(OUTPUT_FILES)//'/irec_main_noise')
    write(IOUT_NOISE,*) 'The main receiver is: (RECEIVER ID)', irec_main_noise
    close(IOUT_NOISE)
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  main station: xi/eta/gamma = ', &
      sngl(xi_receiver(irec_main_noise)),sngl(eta_receiver(irec_main_noise)),sngl(gamma_receiver(irec_main_noise))
    write(IMAIN,*) '  main station: in slice ',islice_selected_rec(irec_main_noise)
    call flush_IMAIN()
  endif

  ! compute source arrays for "ensemble forward source", which is source of "ensemble forward wavefield"
  if (myrank == islice_selected_rec(irec_main_noise) .or. myrank == 0) then ! myrank == 0 is used for output only
    call compute_arrays_source_noise(xi_receiver(irec_main_noise), &
                                     eta_receiver(irec_main_noise), &
                                     gamma_receiver(irec_main_noise), &
                                     nu_rec(:,:,irec_main_noise),noise_sourcearray, &
                                     xigll,yigll,zigll,NSTEP)
  endif

  ! noise distribution and noise direction
  ipoin = 0

  ! loops over surface points
  ! puts noise distrubution and direction onto the surface points
  do iface = 1, num_free_surface_faces

    ispec = free_surface_ispec(iface)

    ! checks if surface element belongs to elastic domain
    if (ispec_is_acoustic(ispec)) then
      print *,'Error noise simulation: element',ispec,'is acoustic'
      call exit_MPI(myrank,'Error: noise for acoustic elements not implemented yet!')
    endif

    do igll = 1, NGLLSQUARE
      i = free_surface_ijk(1,igll,iface)
      j = free_surface_ijk(2,igll,iface)
      k = free_surface_ijk(3,igll,iface)

      ipoin = ipoin + 1
      iglob = ibool(i,j,k,ispec)

      ! this subroutine must be modified by USERS in module user_noise_distribution
      call noise_distribution_direction(xstore(iglob), &
                                        ystore(iglob),zstore(iglob), &
                                        normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                                        mask_noise_out)

      normal_x_noise(ipoin) = normal_x_noise_out
      normal_y_noise(ipoin) = normal_y_noise_out
      normal_z_noise(ipoin) = normal_z_noise_out
      mask_noise(ipoin)     = mask_noise_out
    enddo

  enddo

  end subroutine read_parameters_noise

!
! =============================================================================================================
!

  subroutine check_parameters_noise(NOISE_TOMOGRAPHY,SIMULATION_TYPE,SAVE_FORWARD, &
                                    LOCAL_PATH,NSPEC_TOP,NSTEP)

! check for consistency of the parameters

  use constants, only: CUSTOM_REAL,NDIM,NGLLSQUARE,MAX_STRING_LEN,IOUT_NOISE,OUTPUT_FILES,myrank

  implicit none

  ! input parameters
  integer,intent(in) :: NOISE_TOMOGRAPHY,SIMULATION_TYPE,NSPEC_TOP,NSTEP
  character(len=MAX_STRING_LEN),intent(in) :: LOCAL_PATH
  logical,intent(in) :: SAVE_FORWARD

  ! local parameters
  integer :: reclen
  integer(kind=8) :: filesize
  character(len=MAX_STRING_LEN) :: outputname

  if (myrank == 0) then
    open(unit=IOUT_NOISE,file=trim(OUTPUT_FILES)//'NOISE_SIMULATION', &
         status='unknown',action='write')
    write(IOUT_NOISE,*) '*******************************************************************************'
    write(IOUT_NOISE,*) '*******************************************************************************'
    write(IOUT_NOISE,*) 'WARNING!!!!!!!!!!!!'
    write(IOUT_NOISE,*) 'You are running simulations using NOISE TOMOGRAPHY techniques.'
    write(IOUT_NOISE,*) 'Please make sure you understand the procedures before you have a try.'
    write(IOUT_NOISE,*) 'Displacements everywhere at the free surface are saved every timestep,'
    write(IOUT_NOISE,*) 'so make sure that LOCAL_PATH in Par_file is not global.'
    write(IOUT_NOISE,*) 'Otherwise the disk storage may be a serious issue, as is the speed of I/O.'
    write(IOUT_NOISE,*) 'Also make sure that NO earthquakes are included,'
    write(IOUT_NOISE,*) 'i.e., set moment tensor to be ZERO in CMTSOLUTION'
    write(IOUT_NOISE,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(IOUT_NOISE,*) 'If you just want a regular EARTHQUAKE simulation,'
    write(IOUT_NOISE,*) 'set NOISE_TOMOGRAPHY=0 in Par_file'
    write(IOUT_NOISE,*) '*******************************************************************************'
    write(IOUT_NOISE,*) '*******************************************************************************'
    close(IOUT_NOISE)
  endif

  !no dependancy on movies ...
  !if (.not. USE_HIGHRES_FOR_MOVIES) &
  !  call exit_mpi(myrank,'Please set USE_HIGHRES_FOR_MOVIES in Par_file to be .true.')

  if (NOISE_TOMOGRAPHY == 1) then
    if (SIMULATION_TYPE /= 1) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=1 requires SIMULATION_TYPE=1! check Par_file')

  else if (NOISE_TOMOGRAPHY == 2) then
    if (SIMULATION_TYPE /= 1) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=2 requires SIMULATION_TYPE=1! check Par_file')
    if (.not. SAVE_FORWARD) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=2 requires SAVE_FORWARD=.true.! check Par_file')

  else if (NOISE_TOMOGRAPHY == 3) then
    if (SIMULATION_TYPE /= 3) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=3 requires SIMULATION_TYPE=3! check Par_file')
    if (SAVE_FORWARD) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=3 requires SAVE_FORWARD=.false.! check Par_file')
  endif

  if (NOISE_TOMOGRAPHY /= 0) then
    ! save/read the surface movie using the same c routine as we do for absorbing boundaries (file ID is 2)

    ! size of single record
    reclen = CUSTOM_REAL*NDIM*NGLLSQUARE*NSPEC_TOP

    ! only open files if there are surface faces in this paritition
    if (NSPEC_TOP > 0) then

      ! check integer size limit: size of b_reclen_field must fit onto an 4-byte integer
      if (NSPEC_TOP > int(2147483646.0 / (CUSTOM_REAL * NGLLSQUARE * NDIM))) then
        print *,'reclen of noise surface_movie needed exceeds integer 4-byte limit: ',reclen
        print *,'  ',CUSTOM_REAL, NDIM, NGLLSQUARE, NSPEC_TOP
        print *,'bit size Fortran: ',bit_size(NSPEC_TOP)
        call exit_MPI(myrank,"error NSPEC_TOP integer limit")
      endif

      ! total file size
      filesize = reclen
      filesize = filesize*NSTEP

      write(outputname,"('/proc',i6.6,'_surface_movie')") myrank
      select case(NOISE_TOMOGRAPHY)
      case (1)
        ! write movie files
        call open_file_abs_w(2,trim(LOCAL_PATH)//trim(outputname),len_trim(trim(LOCAL_PATH)//trim(outputname)),filesize)
      case (2)
        ! read movie files
        call open_file_abs_r(2,trim(LOCAL_PATH)//trim(outputname),len_trim(trim(LOCAL_PATH)//trim(outputname)),filesize)
      case (3)
        ! read movie files
        call open_file_abs_r(2,trim(LOCAL_PATH)//trim(outputname),len_trim(trim(LOCAL_PATH)//trim(outputname)),filesize)
      case default
        stop 'Invalid NOISE_TOMOGRAPHY value for noise simulation'
      end select
    endif
  endif

  end subroutine check_parameters_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! read and construct the "source" (source time function based upon noise spectrum)
! for "ensemble forward source"
  subroutine compute_arrays_source_noise(xi_noise,eta_noise,gamma_noise, &
                                         nu_single,noise_sourcearray, &
                                         xigll,yigll,zigll,NSTEP)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,MAX_STRING_LEN, &
    IIN_NOISE,IOUT_NOISE,IMAIN,OUTPUT_FILES,myrank,TINYVAL_SNGL

  implicit none

  ! input parameters
  integer :: NSTEP
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll
  double precision, dimension(NDIM,NDIM) :: nu_single  ! rotation matrix at the main receiver
  ! output parameters
  real(kind=CUSTOM_REAL) :: noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP)

  ! local parameters
  integer :: itime, i, j, k, ier, nlines
  real(kind=CUSTOM_REAL) :: junk
  real(kind=CUSTOM_REAL) :: noise_src(NSTEP),noise_src_u(NDIM,NSTEP)
  double precision, dimension(NDIM) :: nu_main       ! component direction chosen at the main receiver
  double precision :: xi_noise, eta_noise, gamma_noise ! main receiver location

  ! receiver Lagrange interpolators
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar
  double precision :: hpxir(NGLLX), hpetar(NGLLY),hpgammar(NGLLZ)

  character(len=MAX_STRING_LEN) :: filename

  ! main receiver component direction, \nu_main
  filename = trim(OUTPUT_FILES)//'/..//NOISE_TOMOGRAPHY/nu_main'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0 .and. myrank == 0) then
    call exit_MPI(myrank, &
      'file '//trim(filename)//' does NOT exist! nu_main is the component direction (ENZ) for main receiver')
  endif

  do i = 1,3
    read(IIN_NOISE,*,iostat=ier) nu_main(i)
    if (ier /= 0 .and. myrank == 0) &
      call exit_MPI(myrank, &
        'file '//trim(filename)//' has wrong length, the vector should have three components (ENZ)')
  enddo
  close(IIN_NOISE)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  main station: direction vector nu = ',nu_main(:)
    call flush_IMAIN()
  endif

  ! outputs to file for checking
  if (myrank == 0) then
    open(unit=IOUT_NOISE,file=trim(OUTPUT_FILES)//'nu_main',status='unknown',action='write')
    write(IOUT_NOISE,*) 'The direction (ENZ) of selected component of main receiver is', nu_main
    close(IOUT_NOISE)
  endif

  ! initializes
  noise_src(:) = 0._CUSTOM_REAL

  ! noise file (source time function)
  filename = trim(OUTPUT_FILES)//'/..//NOISE_TOMOGRAPHY/S_squared'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0 .and. myrank == 0) then
    call exit_MPI(myrank, 'file '//trim(filename)//' does NOT exist! This file should have been generated using Matlab scripts')
  endif

  ! counts line reads noise source S(t)
  nlines = 0
  do while(ier == 0)
    read(IIN_NOISE,*,iostat=ier) junk, junk
    if (ier == 0)  nlines = nlines + 1
  enddo
  rewind(IIN_NOISE)

  ! checks to be sure that file is matching simulation setup
  if (nlines /= NSTEP) then
    print *,'Error: invalid number of lines ',nlines,' in file NOISE_TOMOGRAPHY/S_squared'
    print *,'       should be equal to NSTEP = ',NSTEP
    print *,'Please check file...'
    call exit_MPI(myrank,'Error invalid number of lines in file S_squared')
  endif

  ! reads noise source S(t)
  do itime = 1,NSTEP
    read(IIN_NOISE,*,iostat=ier) junk, noise_src(itime)
    if (ier /= 0)  call exit_MPI(myrank, &
        'file '//trim(filename)//' has wrong length, please check your simulation duration')
  enddo

  close(IIN_NOISE)

  ! normalizes
  if (maxval(abs(noise_src)) > TINYVAL_SNGL) then
    noise_src(:) = noise_src(:)/maxval(abs(noise_src))
  else
    print *,'Error: noise source S_squared is (almost) zero: absolute max = ',maxval(abs(noise_src))
    print *,'Please check source file NOISE_TOMOGRAPHY/S_squared'
    stop 'Error source S_squared zero'
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  noise source S_squared normalized: min/max = ',minval(noise_src(:)),maxval(noise_src(:))
    call flush_IMAIN()
  endif

  ! rotates to Cartesian
  do itime = 1, NSTEP
    noise_src_u(:,itime) = nu_single(1,:) * noise_src(itime) * nu_main(1) &
                         + nu_single(2,:) * noise_src(itime) * nu_main(2) &
                         + nu_single(3,:) * noise_src(itime) * nu_main(3)
  enddo

  ! receiver interpolators
  call lagrange_any(xi_noise,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta_noise,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma_noise,NGLLZ,zigll,hgammar,hpgammar)

  ! adds interpolated source contribution to all GLL points within this element
  noise_sourcearray(:,:,:,:,:) = 0.0_CUSTOM_REAL
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        do itime = 1, NSTEP
          noise_sourcearray(:,i,j,k,itime) = hxir(i) * hetar(j) * hgammar(k) * noise_src_u(:,itime)
        enddo
      enddo
    enddo
  enddo

  end subroutine compute_arrays_source_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! step 1: calculate the "ensemble forward source"
! add noise spectrum to the location of main receiver
  subroutine add_source_main_rec_noise(nrec,NSTEP,accel,noise_sourcearray, &
                                         ibool,islice_selected_rec,ispec_selected_rec, &
                                         it,irec_main_noise, &
                                         nspec,nglob)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,myrank

  implicit none

  ! input parameters
  integer,intent(in) :: nrec,NSTEP,irec_main_noise
  integer,intent(in) :: nspec,nglob
  integer, dimension(nrec),intent(in) :: islice_selected_rec,ispec_selected_rec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP),intent(in) :: noise_sourcearray

  real(kind=CUSTOM_REAL),dimension(NDIM,nglob),intent(inout) :: accel  ! both input and output

  ! local parameters
  integer :: i,j,k,iglob,ispec, it

  if (irec_main_noise <= 0) then
    print *,'Error rank',myrank,'invalid main id ',irec_main_noise
    stop 'Error invalid irec_main_noise'
  endif

  ! adds noise source (only if this proc carries the noise)
  if (myrank == islice_selected_rec(irec_main_noise)) then

    ispec = ispec_selected_rec(irec_main_noise)

    ! adds nosie source contributions
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          accel(:,iglob) = accel(:,iglob) + noise_sourcearray(:,i,j,k,it)
        enddo
      enddo
    enddo

  endif

  end subroutine add_source_main_rec_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! step 1: calculate the "ensemble forward source"
! save surface movie (displacement) at every time steps, for step 2 & 3.

  subroutine noise_save_surface_movie()

  use constants, only: CUSTOM_REAL,NDIM,NGLLSQUARE

  use specfem_par, only: it,ibool, &
    num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
    Mesh_pointer,GPU_MODE

  use specfem_par_elastic, only: displ
  use specfem_par_noise, only: noise_surface_movie

  implicit none

  ! local parameters
  integer :: ispec,i,j,k,iglob,iface,igll

  ! checks if anything to do
  if (.not. (num_free_surface_faces > 0)) return

  ! writes out wavefield at surface
  if (.not. GPU_MODE) then
    ! on CPU
    ! loops over surface points
    ! get coordinates of surface mesh and surface displacement
    do iface = 1, num_free_surface_faces

      ispec = free_surface_ispec(iface)

      do igll = 1, NGLLSQUARE
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)

        iglob = ibool(i,j,k,ispec)
        noise_surface_movie(:,igll,iface) = displ(:,iglob)
      enddo
    enddo
  else
    ! on GPU
    ! GPU_MODE == 1
    call transfer_surface_to_host(Mesh_pointer,noise_surface_movie)
  endif

  ! save surface motion to disk
  call write_abs(2,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLSQUARE*num_free_surface_faces,it)

  end subroutine noise_save_surface_movie

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! step 2/3: calculate/reconstruct the "ensemble forward wavefield"
! read surface movie (displacement) at every time steps, injected as the source of "ensemble forward wavefield"
! in step 2, call noise_read_add_surface_movie(..., NSTEP-it+1 ,...)
! in step 3, call noise_read_add_surface_movie(..., it ,...)
  subroutine noise_read_add_surface_movie(nmovie_points, &
                                          accel, &
                                          normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                                          ibool,noise_surface_movie,it, &
                                          nspec,nglob, &
                                          num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                                          free_surface_jacobian2Dw)

  use constants

  implicit none

  ! input parameters
  integer,intent(in) :: it,nmovie_points
  integer,intent(in) :: nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel ! both input and output

  real(kind=CUSTOM_REAL), dimension(nmovie_points),intent(in) :: normal_x_noise,normal_y_noise,normal_z_noise, mask_noise

  integer,intent(in) :: num_free_surface_faces
  integer, dimension(num_free_surface_faces),intent(in) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces),intent(in) :: free_surface_ijk
  real(kind=CUSTOM_REAL),intent(in) :: free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces)

  ! local parameters
  integer :: ipoin,ispec,i,j,k,iglob,iface,igll
  real(kind=CUSTOM_REAL) :: eta,val
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_free_surface_faces) :: noise_surface_movie

  ! reads in ensemble noise sources at surface
  if (num_free_surface_faces > 0) then

    ! read surface movie
    call read_abs(2,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLSQUARE*num_free_surface_faces,it)

    ! get coordinates of surface mesh and surface displacement
    ipoin = 0

    ! loops over surface points
    ! puts noise distrubution and direction onto the surface points
    do iface = 1, num_free_surface_faces

      ispec = free_surface_ispec(iface)

      do igll = 1, NGLLSQUARE
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)

        ipoin = ipoin + 1
        iglob = ibool(i,j,k,ispec)

        ! along noise direction
        eta = noise_surface_movie(1,igll,iface) * normal_x_noise(ipoin) + &
              noise_surface_movie(2,igll,iface) * normal_y_noise(ipoin) + &
              noise_surface_movie(3,igll,iface) * normal_z_noise(ipoin)

        val = eta * mask_noise(ipoin) * free_surface_jacobian2Dw(igll,iface) ! wgllwgll_xy(i,j) * jacobian2D_top(i,j,iface)

        accel(1,iglob) = accel(1,iglob) + val * normal_x_noise(ipoin)   ! x-component
        accel(2,iglob) = accel(2,iglob) + val * normal_y_noise(ipoin)
        accel(3,iglob) = accel(3,iglob) + val * normal_z_noise(ipoin)
      enddo
    enddo
  endif

  end subroutine noise_read_add_surface_movie

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================
! On GPU
! step 2/3: calculate/reconstruct the "ensemble forward wavefield"
! read surface movie (displacement) at every time steps, injected as the source of "ensemble forward wavefield"
! in step 2, call noise_read_add_surface_movie_GPU(..., NSTEP-it+1 ,...)
! in step 3, call noise_read_add_surface_movie(..., it ,...)
  subroutine noise_read_add_surface_movie_GPU(noise_surface_movie,it,num_free_surface_faces, &
                                              Mesh_pointer,NOISE_TOMOGRAPHY)

  use constants

  implicit none

  ! input parameters
  integer,intent(in) :: it,num_free_surface_faces

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_free_surface_faces),intent(inout) :: noise_surface_movie

  ! GPU_MODE parameters
  integer(kind=8),intent(in) :: Mesh_pointer
  integer,intent(in) :: NOISE_TOMOGRAPHY

  ! reads in ensemble noise sources at surface
  if (num_free_surface_faces > 0) then

    ! read surface movie
    call read_abs(2,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLSQUARE*num_free_surface_faces,it)

    call noise_read_add_surface_movie_cu(Mesh_pointer,noise_surface_movie,NOISE_TOMOGRAPHY)

  endif

  end subroutine noise_read_add_surface_movie_GPU


! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! step 3: constructing noise source strength kernel

  subroutine compute_kernels_strength_noise(nmovie_points,ibool, &
                                            sigma_kl,displ,deltat,it, &
                                            normal_x_noise,normal_y_noise,normal_z_noise, &
                                            noise_surface_movie, &
                                            nspec,nglob, &
                                            num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                                            GPU_MODE,Mesh_pointer)

  use constants

  implicit none

  ! input parameters
  integer,intent(in) :: it
  integer,intent(in) :: nmovie_points
  integer,intent(in) :: nspec,nglob

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL),intent(in) :: deltat
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ
  real(kind=CUSTOM_REAL), dimension(nmovie_points),intent(in) :: normal_x_noise,normal_y_noise,normal_z_noise

  integer,intent(in) :: num_free_surface_faces
  integer, dimension(num_free_surface_faces),intent(in) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces),intent(in) :: free_surface_ijk

  ! output parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: sigma_kl

  ! local parameters
  integer :: i,j,k,ispec,iglob,ipoin,iface,igll
  real(kind=CUSTOM_REAL) :: eta
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_free_surface_faces) :: noise_surface_movie

  ! GPU_MODE parameters
  integer(kind=8) :: Mesh_pointer
  logical :: GPU_MODE

  ! updates contribution to noise strength kernel
  if (num_free_surface_faces > 0) then

    ! read surface movie, needed for sigma_kl
    call read_abs(2,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLSQUARE*num_free_surface_faces,it)

    if (.not. GPU_MODE) then

      ! noise source strength kernel
      ! to keep similar structure to other kernels, the source strength kernel is saved as a volumetric kernel
      ! but only updated at the surface, because the noise is generated there
      ipoin = 0

      ! loops over surface points
      ! puts noise distrubution and direction onto the surface points
      do iface = 1, num_free_surface_faces

        ispec = free_surface_ispec(iface)

        do igll = 1, NGLLSQUARE
          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)

          ipoin = ipoin + 1
          iglob = ibool(i,j,k,ispec)

          eta = noise_surface_movie(1,igll,iface) * normal_x_noise(ipoin) + &
                noise_surface_movie(2,igll,iface) * normal_y_noise(ipoin) + &
                noise_surface_movie(3,igll,iface) * normal_z_noise(ipoin)

          sigma_kl(i,j,k,ispec) =  sigma_kl(i,j,k,ispec) &
               + deltat * eta * ( normal_x_noise(ipoin) * displ(1,iglob) &
               + normal_y_noise(ipoin) * displ(2,iglob) &
               + normal_z_noise(ipoin) * displ(3,iglob) )
        enddo

      enddo

    else ! GPU_MODE == 1
      call compute_kernels_strgth_noise_cu(Mesh_pointer,noise_surface_movie,deltat)
    endif ! GPU_MODE

  endif

  end subroutine compute_kernels_strength_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! step 3: save noise source strength kernel
  subroutine save_kernels_strength_noise(LOCAL_PATH,sigma_kl,nspec)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,MAX_STRING_LEN,IOUT_NOISE,myrank

  implicit none

  ! input parameters
  integer,intent(in) :: nspec
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: sigma_kl
  character(len=MAX_STRING_LEN),intent(in) :: LOCAL_PATH

  ! local parameters
  character(len=MAX_STRING_LEN) :: prname

  call create_name_database(prname,myrank,LOCAL_PATH)

  open(unit=IOUT_NOISE,file=trim(prname)//'sigma_kernel.bin',status='unknown', &
       form='unformatted',action='write')
  write(IOUT_NOISE) sigma_kl
  close(IOUT_NOISE)

  end subroutine save_kernels_strength_noise
