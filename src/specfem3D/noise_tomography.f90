!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! characterizes noise statistics:
!     for a given point (xcoord,ycoord,zcoord), specify the noise direction "normal_x/y/z_noise"
!     and noise distribution "mask_noise"
!
! USERS: need to modify this subroutine for their own noise characteristics
  subroutine noise_distribution_direction(xcoord_in,ycoord_in,zcoord_in, &
                  normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                  mask_noise_out)
  implicit none
  include "constants.h"
  ! input parameters
  real(kind=CUSTOM_REAL) :: xcoord_in,ycoord_in,zcoord_in
  ! output parameters
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  ! local parameters


  !*****************************************************************************************************************
  !******************************** change your noise characteristics below ****************************************
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! noise direction
  !!!!! here, the noise is assumed to be vertical (SESAME)
  normal_x_noise_out = 0.0
  normal_y_noise_out = 0.0
  normal_z_noise_out = 1.0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  noise distribution
  !!!!! here, the noise is assumed to be uniform
  mask_noise_out = 1.0
  !******************************** change your noise characteristics above ****************************************
  !*****************************************************************************************************************

  end subroutine noise_distribution_direction

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! read parameters
  subroutine read_parameters_noise(myrank,nrec,NSTEP,nmovie_points, &
                                   islice_selected_rec,xi_receiver,eta_receiver,gamma_receiver,nu, &
                                   noise_sourcearray,xigll,yigll,zigll, &
                                   ibool, &                                   
                                   xstore,ystore,zstore, &
                                   irec_master_noise,normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                                   NSPEC_AB_VAL,NGLOB_AB_VAL, &
                                   num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                                   ispec_is_acoustic)
  implicit none
  include "constants.h"
  ! input parameters
  integer :: myrank, nrec, NSTEP, nmovie_points 
  integer :: NSPEC_AB_VAL,NGLOB_AB_VAL
  
  integer, dimension(nrec) :: islice_selected_rec

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB_VAL) :: ibool
  double precision, dimension(nrec)  :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll
  double precision, dimension(NDIM,NDIM,nrec) :: nu
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB_VAL) :: xstore,ystore,zstore

  integer :: num_free_surface_faces 
  integer, dimension(num_free_surface_faces) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces) :: free_surface_ijk
  
  logical, dimension(NSPEC_AB_VAL) :: ispec_is_acoustic

!daniel: from global code...
  !integer, dimension(NSPEC2D_TOP_VAL) :: ibelm_top ! equals free_surface_ispec
  !integer :: NSPEC2D_TOP_VAL ! equals num_free_surface_faces
  !integer :: nspec_top ! equals num_free_surface_faces
  
  ! output parameters
  integer :: irec_master_noise
  real(kind=CUSTOM_REAL) :: noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP)
  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise,mask_noise
  ! local parameters
  integer :: ipoin,ispec,i,j,k,iglob,ios,iface,igll 
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  character(len=256) :: filename

  ! read master receiver ID -- the ID in "STATIONS"
  filename = trim(OUTPUT_FILES_PATH)//'/../OUTPUT_FILES/NOISE_TOMOGRAPHY/irec_master_noise'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ios)
  if( ios /= 0 ) &
    call exit_MPI(myrank, 'file '//trim(filename)//' does NOT exist! This file contains the ID of the master receiver')
  read(IIN_NOISE,*,iostat=ios) irec_master_noise
  if( ios /= 0 ) call exit_MPI(myrank,'error reading file irec_master_noise')
  close(IIN_NOISE)

  ! checks value
  if( irec_master_noise <= 0 ) then
    write(IOUT,*) 'error: irec_master_noise value:',irec_master_noise,'must be positive'
    call exit_MPI(myrank,'error irec_master_noise value')
  endif

  if (myrank == 0) then
    open(unit=IOUT_NOISE,file=trim(OUTPUT_FILES_PATH)//'/irec_master_noise', &
            status='unknown',action='write',iostat=ios)
    if( ios /= 0 ) call exit_MPI(myrank,'error opening file '//trim(OUTPUT_FILES_PATH)//'/irec_master_noise')
    WRITE(IOUT_NOISE,*) 'The master receiver is: (RECEIVER ID)', irec_master_noise
    close(IOUT_NOISE)
  endif

  ! compute source arrays for "ensemble forward source", which is source of "ensemble forward wavefield"
  if(myrank == islice_selected_rec(irec_master_noise) .OR. myrank == 0) then ! myrank == 0 is used for output only
    call compute_arrays_source_noise(myrank, &
              xi_receiver(irec_master_noise),eta_receiver(irec_master_noise),gamma_receiver(irec_master_noise), &
              nu(:,:,irec_master_noise),noise_sourcearray, xigll,yigll,zigll,NSTEP)
  endif

  ! noise distribution and noise direction
  ipoin = 0

  !daniel: from global code, carefull: ngllz must not be face on top...
  !  do ispec2D = 1, nspec_top
  !    ispec = ibelm_top(ispec2D)
  !    k = NGLLZ

  ! loops over surface points
  ! puts noise distrubution and direction onto the surface points
  do iface = 1, num_free_surface_faces

    ispec = free_surface_ispec(iface)

    ! checks if surface element belongs to elastic domain
    if( ispec_is_acoustic(ispec) ) then
      print*,'error noise simulation: element',ispec,'is acoustic'
      stop 'error: noise for acoustic elements not implemented yet!'
    endif

    do igll = 1, NGLLSQUARE
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)

        ipoin = ipoin + 1
        iglob = ibool(i,j,k,ispec)
        ! this subroutine must be modified by USERS
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



!!  !!!BEGIN!!! save mask_noise for check, a file called "mask_noise" is saved in "./OUTPUT_FIELS/"
!!    ipoin = 0
!!      do ispec2D = 1, nspec_top ! NSPEC2D_TOP(IREGION)
!!          ispec = ibelm_top(ispec2D)
!!          k = NGLLZ
!!        ! loop on all the points inside the element
!!          do j = 1,NGLLY,NIT
!!             do i = 1,NGLLX,NIT
!!                ipoin = ipoin + 1
!!                iglob = ibool(i,j,k,ispec)
!!                store_val_x(ipoin) = xstore(iglob)
!!                store_val_y(ipoin) = ystore(iglob)
!!                store_val_z(ipoin) = zstore(iglob)
!!                store_val_ux(ipoin) = mask_noise(ipoin)
!!                store_val_uy(ipoin) = mask_noise(ipoin)
!!                store_val_uz(ipoin) = mask_noise(ipoin)
!!             enddo
!!          enddo
!!      enddo
!!
!!  ! gather info on master proc
!!      ispec = nmovie_points
!!      call MPI_GATHER(store_val_x,ispec,CUSTOM_MPI_TYPE,store_val_x_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
!!      call MPI_GATHER(store_val_y,ispec,CUSTOM_MPI_TYPE,store_val_y_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
!!      call MPI_GATHER(store_val_z,ispec,CUSTOM_MPI_TYPE,store_val_z_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
!!      call MPI_GATHER(store_val_ux,ispec,CUSTOM_MPI_TYPE,store_val_ux_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
!!      call MPI_GATHER(store_val_uy,ispec,CUSTOM_MPI_TYPE,store_val_uy_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
!!      call MPI_GATHER(store_val_uz,ispec,CUSTOM_MPI_TYPE,store_val_uz_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
!!
!!  ! save maks_noise data to disk in home directory
!!  ! this file can be viewed the same way as surface movie data (xcreate_movie_AVS_DX)
!!  ! create_movie_AVS_DX.f90 needs to be modified in order to do that,
!!  ! i.e., instead of showing the normal component, change it to either x, y or z component, or the norm.
!!    if(myrank == 0) then
!!        open(unit=IOUT_NOISE,file='OUTPUT_FILES/mask_noise',status='unknown',form='unformatted',action='write')
!!        write(IOUT_NOISE) store_val_x_all
!!        write(IOUT_NOISE) store_val_y_all
!!        write(IOUT_NOISE) store_val_z_all
!!        write(IOUT_NOISE) store_val_ux_all
!!        write(IOUT_NOISE) store_val_uy_all
!!        write(IOUT_NOISE) store_val_uz_all
!!        close(IOUT_NOISE)
!!     endif
!!  !!!END!!! save mask_noise for check, a file called "mask_noise" is saved in "./OUTPUT_FIELS/"

  end subroutine read_parameters_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! check for consistency of the parameters
  subroutine check_parameters_noise(myrank,NOISE_TOMOGRAPHY,SIMULATION_TYPE,SAVE_FORWARD, &
                                    LOCAL_PATH,NSPEC_TOP,NSTEP)
  implicit none
  include "constants.h"
  ! input parameters
  integer :: myrank,NOISE_TOMOGRAPHY,SIMULATION_TYPE,NSPEC_TOP,NSTEP
  character(len=256) :: LOCAL_PATH
  logical :: SAVE_FORWARD
  ! local parameters
  integer :: reclen
  integer(kind=8) :: filesize
  character(len=256) :: outputname

  if (myrank == 0) then
     open(unit=IOUT_NOISE,file=trim(OUTPUT_FILES_PATH)//'NOISE_SIMULATION', &
          status='unknown',action='write')
     WRITE(IOUT_NOISE,*) '*******************************************************************************'
     WRITE(IOUT_NOISE,*) '*******************************************************************************'
     WRITE(IOUT_NOISE,*) 'WARNING!!!!!!!!!!!!'
     WRITE(IOUT_NOISE,*) 'You are running simulations using NOISE TOMOGRAPHY techniques.'
     WRITE(IOUT_NOISE,*) 'Please make sure you understand the procedures before you have a try.'
     WRITE(IOUT_NOISE,*) 'Displacements everywhere at the free surface are saved every timestep,'
     WRITE(IOUT_NOISE,*) 'so make sure that LOCAL_PATH in Par_file is not global.'
     WRITE(IOUT_NOISE,*) 'Otherwise the disk storage may be a serious issue, as is the speed of I/O.'
     WRITE(IOUT_NOISE,*) 'Also make sure that NO earthquakes are included,'
     WRITE(IOUT_NOISE,*) 'i.e., set moment tensor to be ZERO in CMTSOLUTION'
     WRITE(IOUT_NOISE,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     WRITE(IOUT_NOISE,*) 'If you just want a regular EARTHQUAKE simulation,'
     WRITE(IOUT_NOISE,*) 'set NOISE_TOMOGRAPHY=0 in Par_file'
     WRITE(IOUT_NOISE,*) '*******************************************************************************'
     WRITE(IOUT_NOISE,*) '*******************************************************************************'
     close(IOUT_NOISE)
  endif

  !no dependancy on movies ...
  !if (.not. USE_HIGHRES_FOR_MOVIES) &
  !  call exit_mpi(myrank,'Please set USE_HIGHRES_FOR_MOVIES in Par_file to be .true.')

  if (NOISE_TOMOGRAPHY==1) then
    if (SIMULATION_TYPE/=1) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=1 requires SIMULATION_TYPE=1! check Par_file')

  else if (NOISE_TOMOGRAPHY==2) then
    if (SIMULATION_TYPE/=1) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=2 requires SIMULATION_TYPE=1! check Par_file')
    if (.not. SAVE_FORWARD) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=2 requires SAVE_FORWARD=.true.! check Par_file')

  else if (NOISE_TOMOGRAPHY==3) then
    if (SIMULATION_TYPE/=3) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=3 requires SIMULATION_TYPE=3! check Par_file')
    if (SAVE_FORWARD) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=3 requires SAVE_FORWARD=.false.! check Par_file')
  endif

  if (NOISE_TOMOGRAPHY/=0) then
     ! save/read the surface movie using the same c routine as we do for absorbing boundaries (file ID is 2)

     ! size of single record
     reclen=CUSTOM_REAL*NDIM*NGLLSQUARE*NSPEC_TOP
     ! total file size
     filesize = reclen
     filesize = filesize*NSTEP
     
     write(outputname,"('/proc',i6.6,'_surface_movie')") myrank
     if (NOISE_TOMOGRAPHY==1) call open_file_abs_w(2,trim(LOCAL_PATH)//trim(outputname), &
                                      len_trim(trim(LOCAL_PATH)//trim(outputname)), &
                                      filesize)
     if (NOISE_TOMOGRAPHY==2) call open_file_abs_r(2,trim(LOCAL_PATH)//trim(outputname), &
                                      len_trim(trim(LOCAL_PATH)//trim(outputname)), &
                                      filesize)
     if (NOISE_TOMOGRAPHY==3) call open_file_abs_r(2,trim(LOCAL_PATH)//trim(outputname), &
                                      len_trim(trim(LOCAL_PATH)//trim(outputname)), &
                                      filesize)
  endif

  end subroutine check_parameters_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! read and construct the "source" (source time function based upon noise spectrum)
! for "ensemble forward source"
  subroutine compute_arrays_source_noise(myrank, &
                                         xi_noise,eta_noise,gamma_noise,nu_single,noise_sourcearray, &
                                         xigll,yigll,zigll,NSTEP)
  implicit none
  include 'constants.h'
  ! input parameters
  integer :: myrank, NSTEP
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll
  double precision, dimension(NDIM,NDIM) :: nu_single  ! rotation matrix at the master receiver
  ! output parameters
  real(kind=CUSTOM_REAL) :: noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP)
  ! local parameters
  integer itime, i, j, k, ios
  real(kind=CUSTOM_REAL) :: junk
  real(kind=CUSTOM_REAL) :: noise_src(NSTEP),noise_src_u(NDIM,NSTEP)
  double precision, dimension(NDIM) :: nu_master       ! component direction chosen at the master receiver
  double precision :: xi_noise, eta_noise, gamma_noise ! master receiver location
  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
        hgammar(NGLLZ), hpgammar(NGLLZ)
  character(len=256) :: filename


  noise_src(:) = 0._CUSTOM_REAL
  ! noise file (source time function)
  filename = trim(OUTPUT_FILES_PATH)//'/../OUTPUT_FILES/NOISE_TOMOGRAPHY/S_squared'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ios)
  if( ios /= 0 .and. myrank == 0 )  &
    call exit_MPI(myrank, 'file '//trim(filename)//' does NOT exist! This file should have been generated using Matlab scripts')

  do itime =1,NSTEP
    read(IIN_NOISE,*,iostat=ios) junk, noise_src(itime)
    if( ios /= 0)  call exit_MPI(myrank,&
        'file '//trim(filename)//' has wrong length, please check your simulation duration')
  enddo
  close(IIN_NOISE)
  noise_src(:)=noise_src(:)/maxval(abs(noise_src))


  ! master receiver component direction, \nu_master
  filename = trim(OUTPUT_FILES_PATH)//'/../OUTPUT_FILES/NOISE_TOMOGRAPHY/nu_master'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ios)
  if( ios /= 0 .and. myrank == 0 ) &
    call exit_MPI(myrank,&
      'file '//trim(filename)//' does NOT exist! nu_master is the component direction (ENZ) for master receiver')

  do itime =1,3
    read(IIN_NOISE,*,iostat=ios) nu_master(itime)
    if( ios /= 0 .and. myrank == 0 ) &
      call exit_MPI(myrank,&
        'file '//trim(filename)//' has wrong length, the vector should have three components (ENZ)')
  enddo
  close(IIN_NOISE)

  if (myrank == 0) then
     open(unit=IOUT_NOISE,file=trim(OUTPUT_FILES_PATH)//'nu_master',status='unknown',action='write')
     WRITE(IOUT_NOISE,*) 'The direction (ENZ) of selected component of master receiver is', nu_master
     close(IOUT_NOISE)
  endif

  ! rotates to cartesian
  do itime = 1, NSTEP
    noise_src_u(:,itime) = nu_single(1,:) * noise_src(itime) * nu_master(1) &
                         + nu_single(2,:) * noise_src(itime) * nu_master(2) &
                         + nu_single(3,:) * noise_src(itime) * nu_master(3)
  enddo

  ! receiver interpolators
  call lagrange_any(xi_noise,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta_noise,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma_noise,NGLLZ,zigll,hgammar,hpgammar)

  ! adds interpolated source contribution to all GLL points within this element
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
! add noise spectrum to the location of master receiver
  subroutine add_source_master_rec_noise(myrank,nrec, &
                                NSTEP,accel,noise_sourcearray, &
                                ibool,islice_selected_rec,ispec_selected_rec, &
                                it,irec_master_noise, &
                                NSPEC_AB_VAL,NGLOB_AB_VAL)
  implicit none
  include "constants.h"
  ! input parameters
  integer :: myrank,nrec,NSTEP,irec_master_noise
  integer :: NSPEC_AB_VAL,NGLOB_AB_VAL
  integer, dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB_VAL) :: ibool
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP) :: noise_sourcearray
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB_VAL) :: accel  ! both input and output
  ! output parameters
  ! local parameters
  integer :: i,j,k,iglob,ispec, it

  if( irec_master_noise <= 0 ) then
    print*,'error rank',myrank,irec_master_noise
    stop 'error irec_master_noise'
  endif

  ! adds noise source (only if this proc carries the noise)
  if(myrank == islice_selected_rec(irec_master_noise)) then
    
    ispec = ispec_selected_rec(irec_master_noise)
    
    ! adds nosie source contributions
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          accel(:,iglob) = accel(:,iglob) &
                        + noise_sourcearray(:,i,j,k,it)
        enddo
      enddo
    enddo
  endif

  end subroutine add_source_master_rec_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! step 1: calculate the "ensemble forward source"
! save surface movie (displacement) at every time steps, for step 2 & 3.

  subroutine noise_save_surface_movie(displ, &
                    ibool, &
                    noise_surface_movie,it, &
                    NSPEC_AB_VAL,NGLOB_AB_VAL, &
                    num_free_surface_faces,free_surface_ispec,free_surface_ijk)
  implicit none
  include "constants.h"
  ! input parameters
  integer :: it
  integer :: NSPEC_AB_VAL,NGLOB_AB_VAL
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB_VAL) :: ibool
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB_VAL) ::  displ

  integer :: num_free_surface_faces
  integer, dimension(num_free_surface_faces) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces) :: free_surface_ijk

  !from global code...
  !integer :: nspec_top ! equals num_free_surface_faces
  !integer :: NSPEC2D_TOP_VAL ! equals num_free_surface_faces
  !integer, dimension(NSPEC2D_TOP_VAL) :: ibelm_top ! equals free_surface_ispec
  !integer :: ispec2D ! equals iface
  
  ! local parameters
  integer :: ispec,i,j,k,iglob,iface,igll 
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,num_free_surface_faces) :: noise_surface_movie

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
                  ibool, &
                  noise_surface_movie, &
                  it, &
                  NSPEC_AB_VAL,NGLOB_AB_VAL, &
                  num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                  free_surface_jacobian2Dw)
  implicit none
  include "constants.h"
  ! input parameters
  integer :: it,nmovie_points
  integer :: NSPEC_AB_VAL,NGLOB_AB_VAL
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB_VAL) :: ibool
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB_VAL) :: accel ! both input and output
  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise, mask_noise

  integer :: num_free_surface_faces 
  integer, dimension(num_free_surface_faces) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces) :: free_surface_ijk
  real(kind=CUSTOM_REAL) :: free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces)

  ! from global code...
  !integer :: nspec_top ! equals num_free_surface_faces
  !integer :: NSPEC2D_TOP_VAL ! equal num_free_surface_faces
  !integer, dimension(NSPEC2D_TOP_VAL) :: ibelm_top ! equals free_surface_ispec
  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_VAL) :: jacobian2D_top 
                    ! equals to:                   free_surface_jacobian2Dw including weights wgllwgll
  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  
  ! local parameters
  integer :: ipoin,ispec,i,j,k,iglob,iface,igll
  real(kind=CUSTOM_REAL) :: eta
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_free_surface_faces) :: noise_surface_movie

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

      eta = noise_surface_movie(1,igll,iface) * normal_x_noise(ipoin) + &
            noise_surface_movie(2,igll,iface) * normal_y_noise(ipoin) + &
            noise_surface_movie(3,igll,iface) * normal_z_noise(ipoin)

      accel(1,iglob) = accel(1,iglob) + eta * mask_noise(ipoin) * normal_x_noise(ipoin) &
                                  * free_surface_jacobian2Dw(igll,iface) 
      accel(2,iglob) = accel(2,iglob) + eta * mask_noise(ipoin) * normal_y_noise(ipoin) &
                                  * free_surface_jacobian2Dw(igll,iface)
      accel(3,iglob) = accel(3,iglob) + eta * mask_noise(ipoin) * normal_z_noise(ipoin) &
                                  * free_surface_jacobian2Dw(igll,iface) ! wgllwgll_xy(i,j) * jacobian2D_top(i,j,iface)
    enddo

  enddo

  end subroutine noise_read_add_surface_movie


! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! step 3: constructing noise source strength kernel

  subroutine compute_kernels_strength_noise(nmovie_points,ibool, &
                          Sigma_kl,displ,deltat,it, &
                          normal_x_noise,normal_y_noise,normal_z_noise, &
                          noise_surface_movie, &
                          NSPEC_AB_VAL,NGLOB_AB_VAL, &
                          num_free_surface_faces,free_surface_ispec,free_surface_ijk)
  implicit none
  include "constants.h"
  ! input parameters
  integer :: it
  integer :: nmovie_points
  integer :: NSPEC_AB_VAL,NGLOB_AB_VAL
  
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB_VAL) :: ibool
  real(kind=CUSTOM_REAL) :: deltat
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB_VAL) :: displ
  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise

  integer :: num_free_surface_faces 
  integer, dimension(num_free_surface_faces) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces) :: free_surface_ijk

  ! from global code...
  !integer :: nspec_top ! equals num_free_surface_faces
  !integer :: NSPEC2D_TOP_VAL ! equals num_free_surface_faces
  !integer, dimension(NSPEC2D_TOP_VAL) :: ibelm_top ! equals free_surface_ispec

  ! output parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB_VAL) :: Sigma_kl

  ! local parameters
  integer :: i,j,k,ispec,iglob,ipoin,iface,igll
  real(kind=CUSTOM_REAL) :: eta
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_free_surface_faces) :: noise_surface_movie

  ! read surface movie, needed for Sigma_kl
  call read_abs(2,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLSQUARE*num_free_surface_faces,it)

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

      Sigma_kl(i,j,k,ispec) =  Sigma_kl(i,j,k,ispec) &
         + deltat * eta * ( normal_x_noise(ipoin) * displ(1,iglob) &
                          + normal_y_noise(ipoin) * displ(2,iglob) &
                          + normal_z_noise(ipoin) * displ(3,iglob) )
    enddo

  enddo

  end subroutine compute_kernels_strength_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! step 3: save noise source strength kernel
  subroutine save_kernels_strength_noise(myrank,LOCAL_PATH,Sigma_kl,NSPEC_AB_VAL)
  implicit none
  include "constants.h"
  ! input parameters
  integer myrank
  integer :: NSPEC_AB_VAL
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB_VAL) :: Sigma_kl
  character(len=256) :: LOCAL_PATH
  ! local parameters
  character(len=256) :: prname

  call create_name_database(prname,myrank,LOCAL_PATH)
  
  open(unit=IOUT_NOISE,file=trim(prname)//'sigma_kernel.bin',status='unknown', &
        form='unformatted',action='write')
  write(IOUT_NOISE) Sigma_kl
  close(IOUT_NOISE)

  end subroutine save_kernels_strength_noise
