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

!----
!---- locate_receivers finds the correct position of the receivers
!----
  subroutine locate_receivers(ibool,myrank,NSPEC_AB,NGLOB_AB,NGNOD, &
                 xstore,ystore,zstore,xigll,yigll,zigll,rec_filename, &
                 nrec,islice_selected_rec,ispec_selected_rec, &
                 xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
                 NPROC,utm_x_source,utm_y_source, &
                 UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                 iglob_is_surface_external_mesh,ispec_is_surface_external_mesh, &
                 num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  implicit none

  include "constants.h"

  integer :: myrank
  integer :: NSPEC_AB,NGLOB_AB,NGNOD

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  ! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  ! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

  ! input receiver file name
  character(len=*) rec_filename

  ! receivers
  integer :: nrec
  integer, dimension(nrec),intent(out) :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(nrec),intent(out) :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(3,3,nrec),intent(out) :: nu
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec),intent(out) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec),intent(out) :: network_name

  integer :: NPROC,UTM_PROJECTION_ZONE
  double precision :: utm_x_source,utm_y_source
  logical :: SUPPRESS_UTM_PROJECTION

  ! free surface
  integer :: num_free_surface_faces
  logical, dimension(NGLOB_AB) :: iglob_is_surface_external_mesh
  logical, dimension(NSPEC_AB) :: ispec_is_surface_external_mesh
  integer, dimension(num_free_surface_faces) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces) :: free_surface_ijk

  ! local parameters

  ! for surface locating and normal computing with external mesh
  real(kind=CUSTOM_REAL), dimension(3) :: u_vector,v_vector,w_vector
  integer :: pt0_ix,pt0_iy,pt0_iz,pt1_ix,pt1_iy,pt1_iz,pt2_ix,pt2_iy,pt2_iz

  integer, allocatable, dimension(:) :: ix_initial_guess,iy_initial_guess,iz_initial_guess
  integer iprocloop
  integer ios

  double precision,dimension(1) :: altitude_rec,distmin_ele
  !double precision,dimension(4) :: elevation_node,dist_node
  double precision,dimension(NPROC) :: distmin_ele_all,elevation_all

  real(kind=CUSTOM_REAL) :: xloc,yloc,loc_ele,loc_distmin

  double precision, allocatable, dimension(:) :: x_target,y_target,z_target

  double precision, allocatable, dimension(:) :: horiz_dist
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision dist

  double precision xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma
  double precision x,y,z
  double precision xix,xiy,xiz
  double precision etax,etay,etaz
  double precision gammax,gammay,gammaz

  ! coordinates of the control points of the surface element
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer irec
  integer i,j,k,ispec,iglob
  integer imin,imax,jmin,jmax,kmin,kmax
  !integer iface,inode,igll,jgll,kgll

!  integer iselected,jselected,iface_selected,iadjust,jadjust
  integer iproc(1)

  ! topology of the control points of the surface element
  integer iax,iay,iaz
  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

  integer iter_loop,ispec_iterate
  integer ia

  ! timer MPI
  double precision, external :: wtime
  double precision time_start,tCPU

  ! use dynamic allocation
  double precision, dimension(:), allocatable :: final_distance
  double precision distmin,final_distance_max

  ! receiver information
  ! station information for writing the seismograms
!  integer :: iglob_selected
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur,stutm_x,stutm_y,elevation
  double precision, allocatable, dimension(:) :: x_found_all,y_found_all,z_found_all
  double precision, dimension(:), allocatable :: final_distance_all
  double precision, allocatable, dimension(:) :: xi_receiver_all,eta_receiver_all,gamma_receiver_all
  double precision, allocatable, dimension(:,:,:) :: nu_all
  integer, allocatable, dimension(:) :: ispec_selected_rec_all

  integer :: ier
  character(len=256) OUTPUT_FILES

  real(kind=CUSTOM_REAL) :: xmin,xmax,ymin,ymax,zmin,zmax
  real(kind=CUSTOM_REAL) :: xmin_ELE,xmax_ELE,ymin_ELE,ymax_ELE,zmin_ELE,zmax_ELE
  integer :: imin_temp,imax_temp,jmin_temp,jmax_temp,kmin_temp,kmax_temp
  integer,dimension(NGLLX*NGLLY*NGLLZ) :: iglob_temp

  ! SU_FORMAT parameters
  double precision :: llat,llon,lele,lbur
  logical :: SU_station_file_exists

  ! get MPI starting time
  time_start = wtime()

  ! user output
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    write(IMAIN,'(1x,a,a,a)') 'reading receiver information from ', trim(rec_filename), ' file'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! dimension of model in current proc
  xmin=minval(xstore(:));          xmax=maxval(xstore(:))
  ymin=minval(ystore(:));          ymax=maxval(ystore(:))
  zmin=minval(zstore(:));          zmax=maxval(zstore(:))

  ! loop limits
  if(FASTER_RECEIVERS_POINTS_ONLY) then
    imin_temp = 1; imax_temp = NGLLX
    jmin_temp = 1; jmax_temp = NGLLY
    kmin_temp = 1; kmax_temp = NGLLZ
  else
    imin_temp = 2; imax_temp = NGLLX - 1
    jmin_temp = 2; jmax_temp = NGLLY - 1
    kmin_temp = 2; kmax_temp = NGLLZ - 1
  endif

  ! define topology of the control element
  call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)

  ! opens STATIONS file
  open(unit=IIN,file=trim(rec_filename),status='old',action='read',iostat=ios)
  if (ios /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))

  ! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)))

  ! checks if station locations already available
  if( SU_FORMAT ) then
    ! checks if file with station infos located from previous run exists
    INQUIRE(file=trim(OUTPUT_FILES)//'/SU_stations_info.bin',exist=SU_station_file_exists)
    if ( SU_station_file_exists ) then
      ! all processes read in stations names from STATIONS file
      do irec=1,nrec
        read(IIN,*,iostat=ios) station_name(irec),network_name(irec),llat,llon,lele,lbur
        if (ios /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))
      enddo
      close(IIN)
      ! master reads in available station information
      if( myrank == 0 ) then
        open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/SU_stations_info.bin', &
              status='old',action='read',form='unformatted',iostat=ios)
        if (ios /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))

        write(IMAIN,*) 'station details from SU_stations_info.bin'
        call flush_IMAIN()

        allocate(x_found(nrec),y_found(nrec),z_found(nrec))
        ! reads in station infos
        read(IOUT_SU) islice_selected_rec,ispec_selected_rec
        read(IOUT_SU) xi_receiver,eta_receiver,gamma_receiver
        read(IOUT_SU) x_found,y_found,z_found
        read(IOUT_SU) nu
        close(IOUT_SU)
        ! write the locations of stations, so that we can load them and write them to SU headers later
        open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/output_list_stations.txt', &
              status='unknown',action='write',iostat=ios)
        if( ios /= 0 ) &
          call exit_mpi(myrank,'error opening file '//trim(OUTPUT_FILES)//'/output_list_stations.txt')

        do irec=1,nrec
          write(IOUT_SU,*) x_found(irec),y_found(irec),z_found(irec)
        enddo

        close(IOUT_SU)
        deallocate(x_found,y_found,z_found)
      endif
      ! main process broadcasts the results to all the slices
      call bcast_all_i(islice_selected_rec,nrec)
      call bcast_all_i(ispec_selected_rec,nrec)
      call bcast_all_dp(xi_receiver,nrec)
      call bcast_all_dp(eta_receiver,nrec)
      call bcast_all_dp(gamma_receiver,nrec)
      call bcast_all_dp(nu,NDIM*NDIM*nrec)
      call sync_all()
      ! user output
      if( myrank == 0 ) then
        ! elapsed time since beginning of mesh generation
        tCPU = wtime() - time_start
        write(IMAIN,*)
        write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
        write(IMAIN,*)
        write(IMAIN,*) 'End of receiver detection - done'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
      ! everything done
      return
    endif
  endif

  ! allocate memory for arrays using number of stations
  allocate(stlat(nrec), &
          stlon(nrec), &
          stele(nrec), &
          stbur(nrec), &
          stutm_x(nrec), &
          stutm_y(nrec), &
          horiz_dist(nrec), &
          elevation(nrec), &
          ix_initial_guess(nrec), &
          iy_initial_guess(nrec), &
          iz_initial_guess(nrec), &
          x_target(nrec), &
          y_target(nrec), &
          z_target(nrec), &
          x_found(nrec), &
          y_found(nrec), &
          z_found(nrec), &
          final_distance(nrec), &
          ispec_selected_rec_all(nrec), &
          xi_receiver_all(nrec), &
          eta_receiver_all(nrec), &
          gamma_receiver_all(nrec), &
          x_found_all(nrec), &
          y_found_all(nrec), &
          z_found_all(nrec), &
          final_distance_all(nrec), &
          nu_all(3,3,nrec),stat=ier)
  if( ier /= 0 ) stop 'error allocating arrays for locating receivers'

  ! loop on all the stations
  do irec=1,nrec

    read(IIN,*,iostat=ios) station_name(irec),network_name(irec), &
                          stlat(irec),stlon(irec),stele(irec),stbur(irec)

    if (ios /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))

    ! convert station location to UTM
    call utm_geo(stlon(irec),stlat(irec),stutm_x(irec),stutm_y(irec),&
                UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

    ! compute horizontal distance between source and receiver in km
    horiz_dist(irec) = dsqrt((stutm_y(irec)-utm_y_source)**2 &
                            + (stutm_x(irec)-utm_x_source)**2) / 1000.d0

    ! print some information about stations
    if(myrank == 0) then
      ! limits user output if too many receivers
      if( nrec < 1000 .and. ( .not. SU_FORMAT )  ) then
        write(IMAIN,*) 'Station #',irec,': ',station_name(irec)(1:len_trim(station_name(irec))), &
                       '.',network_name(irec)(1:len_trim(network_name(irec))), &
                       '    horizontal distance:  ',sngl(horiz_dist(irec)),' km'
      endif
    endif

    ! get approximate topography elevation at source long/lat coordinates
    xloc = stutm_x(irec)
    yloc = stutm_y(irec)
    call get_topo_elevation_free(xloc,yloc,loc_ele,loc_distmin, &
                                NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                num_free_surface_faces,free_surface_ispec,free_surface_ijk)
    altitude_rec(1) = loc_ele
    distmin_ele(1) = loc_distmin

    !  MPI communications to determine the best slice
    call gather_all_dp(distmin_ele,1,distmin_ele_all,1,NPROC)
    call gather_all_dp(altitude_rec,1,elevation_all,1,NPROC)

    if(myrank == 0) then
      iproc = minloc(distmin_ele_all)
      altitude_rec(1) = elevation_all(iproc(1))
    endif
    call bcast_all_dp(altitude_rec,1)
    elevation(irec) = altitude_rec(1)


!     get the Cartesian components of n in the model: nu

    ! orientation consistent with the UTM projection
    ! X coordinate - East
    nu(1,1,irec) = 1.d0
    nu(1,2,irec) = 0.d0
    nu(1,3,irec) = 0.d0
    ! Y coordinate - North
    nu(2,1,irec) = 0.d0
    nu(2,2,irec) = 1.d0
    nu(2,3,irec) = 0.d0
    ! Z coordinate - Vertical
    nu(3,1,irec) = 0.d0
    nu(3,2,irec) = 0.d0
    nu(3,3,irec) = 1.d0

    x_target(irec) = stutm_x(irec)
    y_target(irec) = stutm_y(irec)

    ! receiver's Z coordinate
    if( USE_SOURCES_RECEIVERS_Z ) then
      ! alternative: burial depth is given as z value directly
      z_target(irec) = stbur(irec)
    else
      ! burial depth in STATIONS file given in m
      z_target(irec) = elevation(irec) - stbur(irec)
    endif

    ! reset distance to huge initial value
    distmin=HUGEVAL

    if (.not. SU_FORMAT) then
      ! determines closest GLL point
      ispec_selected_rec(irec) = 0
      do ispec=1,NSPEC_AB

        ! define the interval in which we look for points
        if(FASTER_RECEIVERS_POINTS_ONLY) then
          imin = 1
          imax = NGLLX
          jmin = 1
          jmax = NGLLY
          kmin = 1
          kmax = NGLLZ
        else
          ! loop only on points inside the element
          ! exclude edges to ensure this point is not shared with other elements
          imin = 2
          imax = NGLLX - 1
          jmin = 2
          jmax = NGLLY - 1
          kmin = 2
          kmax = NGLLZ - 1
        endif

        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax

              iglob = ibool(i,j,k,ispec)

              if (.not. RECEIVERS_CAN_BE_BURIED) then
                if ((.not. iglob_is_surface_external_mesh(iglob)) &
                  .or. (.not. ispec_is_surface_external_mesh(ispec))) then
                  cycle
                endif
              endif

              dist = (x_target(irec)-dble(xstore(iglob)))**2 &
                   + (y_target(irec)-dble(ystore(iglob)))**2 &
                   + (z_target(irec)-dble(zstore(iglob)))**2

              ! keep this point if it is closer to the receiver
              if(dist < distmin) then
                distmin = dist
                ispec_selected_rec(irec) = ispec
                ix_initial_guess(irec) = i
                iy_initial_guess(irec) = j
                iz_initial_guess(irec) = k

                xi_receiver(irec) = dble(ix_initial_guess(irec))
                eta_receiver(irec) = dble(iy_initial_guess(irec))
                gamma_receiver(irec) = dble(iz_initial_guess(irec))
                x_found(irec) = xstore(iglob)
                y_found(irec) = ystore(iglob)
                z_found(irec) = zstore(iglob)
              endif

            enddo
          enddo
        enddo

        ! compute final distance between asked and found (converted to km)
        final_distance(irec) = dsqrt((x_target(irec)-x_found(irec))**2 &
                                   + (y_target(irec)-y_found(irec))**2 &
                                   + (z_target(irec)-z_found(irec))**2)

      ! end of loop on all the spectral elements in current slice
      enddo
    else
      ! SeismicUnix format
      ispec_selected_rec(irec) = 0
      ix_initial_guess(irec) = 0
      iy_initial_guess(irec) = 0
      iz_initial_guess(irec) = 0
      final_distance(irec) = HUGEVAL
      if ( (x_target(irec)>=xmin .and. x_target(irec)<=xmax) .and. &
          (y_target(irec)>=ymin .and. y_target(irec)<=ymax) .and. &
          (z_target(irec)>=zmin .and. z_target(irec)<=zmax) ) then
        do ispec=1,NSPEC_AB
          iglob_temp=reshape(ibool(:,:,:,ispec),(/NGLLX*NGLLY*NGLLZ/))
          xmin_ELE=minval(xstore(iglob_temp))
          xmax_ELE=maxval(xstore(iglob_temp))
          ymin_ELE=minval(ystore(iglob_temp))
          ymax_ELE=maxval(ystore(iglob_temp))
          zmin_ELE=minval(zstore(iglob_temp))
          zmax_ELE=maxval(zstore(iglob_temp))
          if ( (x_target(irec)>=xmin_ELE .and. x_target(irec)<=xmax_ELE) .and. &
              (y_target(irec)>=ymin_ELE .and. y_target(irec)<=ymax_ELE) .and. &
              (z_target(irec)>=zmin_ELE .and. z_target(irec)<=zmax_ELE) ) then
            ! we find the element (ispec) which "may" contain the receiver (irec)
            ! so we only need to compute distances
            !(which is expensive because of "dsqrt") within those elements
            ispec_selected_rec(irec) = ispec
            do k = kmin_temp,kmax_temp
              do j = jmin_temp,jmax_temp
                do i = imin_temp,imax_temp
                  iglob = ibool(i,j,k,ispec)
                  ! for comparison purpose, we do not have to do "dsqrt", which is expensive
                  dist = ((x_target(irec)-dble(xstore(iglob)))**2 &
                        + (y_target(irec)-dble(ystore(iglob)))**2 &
                        + (z_target(irec)-dble(zstore(iglob)))**2)
                  if(dist < distmin) then
                    distmin = dist
                    ix_initial_guess(irec) = i
                    iy_initial_guess(irec) = j
                    iz_initial_guess(irec) = k
                    xi_receiver(irec) = dble(ix_initial_guess(irec))
                    eta_receiver(irec) = dble(iy_initial_guess(irec))
                    gamma_receiver(irec) = dble(iz_initial_guess(irec))
                    x_found(irec) = xstore(iglob)
                    y_found(irec) = ystore(iglob)
                    z_found(irec) = zstore(iglob)
                  endif
                enddo
              enddo
            enddo
            final_distance(irec) = dsqrt((x_target(irec)-x_found(irec))**2 &
                                       + (y_target(irec)-y_found(irec))**2 &
                                       + (z_target(irec)-z_found(irec))**2)
          endif ! if receiver "may" be within this element
        enddo ! do ispec=1,NSPEC_AB
      endif ! if receiver "may" be within this proc
    endif !if (.not. SU_FORMAT)

    if (ispec_selected_rec(irec) == 0) then
      ! receiver is NOT within this proc, assign trivial values
      ispec_selected_rec(irec) = 1
      ix_initial_guess(irec) = 1
      iy_initial_guess(irec) = 1
      iz_initial_guess(irec) = 1
      final_distance(irec) = HUGEVAL
    endif

    ! get normal to the face of the hexaedra if receiver is on the surface
    if ((.not. RECEIVERS_CAN_BE_BURIED) .and. &
       .not. (ispec_selected_rec(irec) == 1)) then
      pt0_ix = -1
      pt0_iy = -1
      pt0_iz = -1
      pt1_ix = -1
      pt1_iy = -1
      pt1_iz = -1
      pt2_ix = -1
      pt2_iy = -1
      pt2_iz = -1
      ! we get two vectors of the face (three points) to compute the normal
      if (ix_initial_guess(irec) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(1,2,2,ispec_selected_rec(irec)))) then
        pt0_ix = 1
        pt0_iy = NGLLY
        pt0_iz = 1
        pt1_ix = 1
        pt1_iy = 1
        pt1_iz = 1
        pt2_ix = 1
        pt2_iy = NGLLY
        pt2_iz = NGLLZ
      endif
      if (ix_initial_guess(irec) == NGLLX .and. &
         iglob_is_surface_external_mesh(ibool(NGLLX,2,2,ispec_selected_rec(irec)))) then
        pt0_ix = NGLLX
        pt0_iy = 1
        pt0_iz = 1
        pt1_ix = NGLLX
        pt1_iy = NGLLY
        pt1_iz = 1
        pt2_ix = NGLLX
        pt2_iy = 1
        pt2_iz = NGLLZ
      endif
      if (iy_initial_guess(irec) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(2,1,2,ispec_selected_rec(irec)))) then
        pt0_ix = 1
        pt0_iy = 1
        pt0_iz = 1
        pt1_ix = NGLLX
        pt1_iy = 1
        pt1_iz = 1
        pt2_ix = 1
        pt2_iy = 1
        pt2_iz = NGLLZ
      endif
      if (iy_initial_guess(irec) == NGLLY .and. &
         iglob_is_surface_external_mesh(ibool(2,NGLLY,2,ispec_selected_rec(irec)))) then
        pt0_ix = NGLLX
        pt0_iy = NGLLY
        pt0_iz = 1
        pt1_ix = 1
        pt1_iy = NGLLY
        pt1_iz = 1
        pt2_ix = NGLLX
        pt2_iy = NGLLY
        pt2_iz = NGLLZ
      endif
      if (iz_initial_guess(irec) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(2,2,1,ispec_selected_rec(irec)))) then
        pt0_ix = NGLLX
        pt0_iy = 1
        pt0_iz = 1
        pt1_ix = 1
        pt1_iy = 1
        pt1_iz = 1
        pt2_ix = NGLLX
        pt2_iy = NGLLY
        pt2_iz = 1
      endif
      if (iz_initial_guess(irec) == NGLLZ .and. &
         iglob_is_surface_external_mesh(ibool(2,2,NGLLZ,ispec_selected_rec(irec)))) then
        pt0_ix = 1
        pt0_iy = 1
        pt0_iz = NGLLZ
        pt1_ix = NGLLX
        pt1_iy = 1
        pt1_iz = NGLLZ
        pt2_ix = 1
        pt2_iy = NGLLY
        pt2_iz = NGLLZ
      endif

      if (pt0_ix<0 .or.pt0_iy<0 .or. pt0_iz<0 .or. &
         pt1_ix<0 .or. pt1_iy<0 .or. pt1_iz<0 .or. &
         pt2_ix<0 .or. pt2_iy<0 .or. pt2_iz<0) then
        stop 'error in computing normal for receivers.'
      endif

      u_vector(1) = xstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_rec(irec))) &
           - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))
      u_vector(2) = ystore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_rec(irec))) &
           - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))
      u_vector(3) = zstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_rec(irec))) &
           - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))
      v_vector(1) = xstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_rec(irec))) &
           - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))
      v_vector(2) = ystore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_rec(irec))) &
           - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))
      v_vector(3) = zstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_rec(irec))) &
           - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))

      ! cross product
      w_vector(1) = u_vector(2)*v_vector(3) - u_vector(3)*v_vector(2)
      w_vector(2) = u_vector(3)*v_vector(1) - u_vector(1)*v_vector(3)
      w_vector(3) = u_vector(1)*v_vector(2) - u_vector(2)*v_vector(1)

      ! normalize vector w
      w_vector(:) = w_vector(:)/sqrt(w_vector(1)**2+w_vector(2)**2+w_vector(3)**2)

      ! build the two other vectors for a direct base: we normalize u, and v=w^u
      u_vector(:) = u_vector(:)/sqrt(u_vector(1)**2+u_vector(2)**2+u_vector(3)**2)
      v_vector(1) = w_vector(2)*u_vector(3) - w_vector(3)*u_vector(2)
      v_vector(2) = w_vector(3)*u_vector(1) - w_vector(1)*u_vector(3)
      v_vector(3) = w_vector(1)*u_vector(2) - w_vector(2)*u_vector(1)

      ! build rotation matrice nu for seismograms
      if (EXTERNAL_MESH_RECEIVERS_NORMAL) then
        !     East (u)
        nu(1,1,irec) = u_vector(1)
        nu(1,2,irec) = v_vector(1)
        nu(1,3,irec) = w_vector(1)

        !     North (v)
        nu(2,1,irec) = u_vector(2)
        nu(2,2,irec) = v_vector(2)
        nu(2,3,irec) = w_vector(2)

        !     Vertical (w)
        nu(3,1,irec) = u_vector(3)
        nu(3,2,irec) = v_vector(3)
        nu(3,3,irec) = w_vector(3)
      else
        !     East
        nu(1,1,irec) = 1.d0
        nu(1,2,irec) = 0.d0
        nu(1,3,irec) = 0.d0

        !     North
        nu(2,1,irec) = 0.d0
        nu(2,2,irec) = 1.d0
        nu(2,3,irec) = 0.d0

        !     Vertical
        nu(3,1,irec) = 0.d0
        nu(3,2,irec) = 0.d0
        nu(3,3,irec) = 1.d0
      endif

    endif ! of if (.not. RECEIVERS_CAN_BE_BURIED)

  ! end of loop on all the stations
  enddo

  ! close receiver file
  close(IIN)

! ****************************************
! find the best (xi,eta,gamma) for each receiver
! ****************************************

  if(.not. FASTER_RECEIVERS_POINTS_ONLY) then

    ! loop on all the receivers to iterate in that slice
    do irec = 1,nrec

      ispec_iterate = ispec_selected_rec(irec)

      ! use initial guess in xi and eta
      xi = xigll(ix_initial_guess(irec))
      eta = yigll(iy_initial_guess(irec))
      gamma = zigll(iz_initial_guess(irec))

      ! define coordinates of the control points of the element
      do ia=1,NGNOD
        iax = 0
        iay = 0
        iaz = 0
        if(iaddx(ia) == 0) then
          iax = 1
        else if(iaddx(ia) == 1) then
          iax = (NGLLX+1)/2
        else if(iaddx(ia) == 2) then
          iax = NGLLX
        else
          call exit_MPI(myrank,'incorrect value of iaddx')
        endif

        if(iaddy(ia) == 0) then
          iay = 1
        else if(iaddy(ia) == 1) then
          iay = (NGLLY+1)/2
        else if(iaddy(ia) == 2) then
          iay = NGLLY
        else
          call exit_MPI(myrank,'incorrect value of iaddy')
        endif

        if(iaddz(ia) == 0) then
          iaz = 1
        else if(iaddz(ia) == 1) then
          iaz = (NGLLZ+1)/2
        else if(iaddz(ia) == 2) then
          iaz = NGLLZ
        else
          call exit_MPI(myrank,'incorrect value of iaddz')
        endif

        iglob = ibool(iax,iay,iaz,ispec_iterate)
        xelm(ia) = dble(xstore(iglob))
        yelm(ia) = dble(ystore(iglob))
        zelm(ia) = dble(zstore(iglob))

      enddo

      ! iterate to solve the non linear system
      do iter_loop = 1,NUM_ITER

        ! impose receiver exactly at the surface
        !    gamma = 1.d0

        ! recompute jacobian for the new point
        call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,NGNOD)

        ! compute distance to target location
        dx = - (x - x_target(irec))
        dy = - (y - y_target(irec))
        dz = - (z - z_target(irec))

        ! compute increments
        ! gamma does not change since we know the receiver is exactly on the surface
        dxi  = xix*dx + xiy*dy + xiz*dz
        deta = etax*dx + etay*dy + etaz*dz
        dgamma = gammax*dx + gammay*dy + gammaz*dz

        ! update values
        xi = xi + dxi
        eta = eta + deta
        gamma = gamma + dgamma

        ! impose that we stay in that element
        ! (useful if user gives a receiver outside the mesh for instance)
        ! we can go slightly outside the [1,1] segment since with finite elements
        ! the polynomial solution is defined everywhere
        ! this can be useful for convergence of itertive scheme with distorted elements
        if (xi > 1.10d0) xi = 1.10d0
        if (xi < -1.10d0) xi = -1.10d0
        if (eta > 1.10d0) eta = 1.10d0
        if (eta < -1.10d0) eta = -1.10d0
        if (gamma > 1.10d0) gamma = 1.10d0
        if (gamma < -1.10d0) gamma = -1.10d0

      ! end of non linear iterations
      enddo

      ! impose receiver exactly at the surface after final iteration
      !  gamma = 1.d0

      ! compute final coordinates of point found
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,NGNOD)

      ! store xi,eta and x,y,z of point found
      xi_receiver(irec) = xi
      eta_receiver(irec) = eta
      gamma_receiver(irec) = gamma
      x_found(irec) = x
      y_found(irec) = y
      z_found(irec) = z

      ! compute final distance between asked and found (converted to km)
      final_distance(irec) = dsqrt((x_target(irec)-x_found(irec))**2 + &
        (y_target(irec)-y_found(irec))**2 + (z_target(irec)-z_found(irec))**2)

    enddo

  endif ! of if (.not. FASTER_RECEIVERS_POINTS_ONLY)

  ! synchronize all the processes to make sure all the estimates are available
  call sync_all()
  ! for MPI version, gather information from all the nodes
  if (myrank/=0) then ! gather information from other processors (one at a time)
     call send_i(ispec_selected_rec, nrec,0,0)
     call send_dp(xi_receiver,       nrec,0,1)
     call send_dp(eta_receiver,      nrec,0,2)
     call send_dp(gamma_receiver,    nrec,0,3)
     call send_dp(final_distance,    nrec,0,4)
     call send_dp(x_found,           nrec,0,5)
     call send_dp(y_found,           nrec,0,6)
     call send_dp(z_found,           nrec,0,7)
     call send_dp(nu,            3*3*nrec,0,8)
  else
     islice_selected_rec(:) = 0
     do iprocloop=1,NPROC-1
        call recv_i(ispec_selected_rec_all, nrec,iprocloop,0)
        call recv_dp(xi_receiver_all,       nrec,iprocloop,1)
        call recv_dp(eta_receiver_all,      nrec,iprocloop,2)
        call recv_dp(gamma_receiver_all,    nrec,iprocloop,3)
        call recv_dp(final_distance_all,    nrec,iprocloop,4)
        call recv_dp(x_found_all,           nrec,iprocloop,5)
        call recv_dp(y_found_all,           nrec,iprocloop,6)
        call recv_dp(z_found_all,           nrec,iprocloop,7)
        call recv_dp(nu_all,            3*3*nrec,iprocloop,8)
        do irec=1,nrec
           if (final_distance_all(irec) < final_distance(irec)) then
              final_distance(irec) = final_distance_all(irec)
              islice_selected_rec(irec) = iprocloop
              ispec_selected_rec(irec) = ispec_selected_rec_all(irec)
              xi_receiver(irec) = xi_receiver_all(irec)
              eta_receiver(irec) = eta_receiver_all(irec)
              gamma_receiver(irec) = gamma_receiver_all(irec)
              x_found(irec) = x_found_all(irec)
              y_found(irec) = y_found_all(irec)
              z_found(irec) = z_found_all(irec)
              nu(:,:,irec) = nu_all(:,:,irec)
           endif
        enddo
     enddo
  endif

  ! this is executed by main process only
  if(myrank == 0) then

    do irec=1,nrec

      ! checks stations location
      if(final_distance(irec) == HUGEVAL) then
        write(IMAIN,*) 'error locating station # ',irec,'    ',station_name(irec),network_name(irec)
        call exit_MPI(myrank,'error locating receiver')
      endif

      ! limits user output if too many receivers
      if( nrec < 1000 .and. (.not. SU_FORMAT ) ) then

      write(IMAIN,*)
      write(IMAIN,*) 'station # ',irec,'    ',station_name(irec),network_name(irec)

      write(IMAIN,*) '     original latitude: ',sngl(stlat(irec))
      write(IMAIN,*) '     original longitude: ',sngl(stlon(irec))
      if( SUPPRESS_UTM_PROJECTION ) then
        write(IMAIN,*) '     original x: ',sngl(stutm_x(irec))
        write(IMAIN,*) '     original y: ',sngl(stutm_y(irec))
      else
        write(IMAIN,*) '     original UTM x: ',sngl(stutm_x(irec))
        write(IMAIN,*) '     original UTM y: ',sngl(stutm_y(irec))
      endif
      if( USE_SOURCES_RECEIVERS_Z ) then
        write(IMAIN,*) '     original z: ',sngl(stbur(irec))
      else
        write(IMAIN,*) '     original depth: ',sngl(stbur(irec)),' m'
      endif
      write(IMAIN,*) '     horizontal distance: ',sngl(horiz_dist(irec))
      write(IMAIN,*) '     target x, y, z: ',sngl(x_target(irec)),sngl(y_target(irec)),sngl(z_target(irec))

      write(IMAIN,*) '     closest estimate found: ',sngl(final_distance(irec)),' m away'
      write(IMAIN,*) '     in slice ',islice_selected_rec(irec),' in element ',ispec_selected_rec(irec)
      if(FASTER_RECEIVERS_POINTS_ONLY) then
        write(IMAIN,*) '     in point i,j,k = ',nint(xi_receiver(irec)), &
                                       nint(eta_receiver(irec)), &
                                       nint(gamma_receiver(irec))
        write(IMAIN,*) '     nu1 = ',nu(1,:,irec)
        write(IMAIN,*) '     nu2 = ',nu(2,:,irec)
        write(IMAIN,*) '     nu3 = ',nu(3,:,irec)
      else
        write(IMAIN,*) '     at coordinates: '
        write(IMAIN,*) '     xi    = ',xi_receiver(irec)
        write(IMAIN,*) '     eta   = ',eta_receiver(irec)
        write(IMAIN,*) '     gamma = ',gamma_receiver(irec)
      endif
      if( SUPPRESS_UTM_PROJECTION ) then
        write(IMAIN,*) '     x: ',x_found(irec)
        write(IMAIN,*) '     y: ',y_found(irec)
      else
        write(IMAIN,*) '     UTM x: ',x_found(irec)
        write(IMAIN,*) '     UTM y: ',y_found(irec)
      endif
      if( USE_SOURCES_RECEIVERS_Z ) then
        write(IMAIN,*) '     z: ',z_found(irec)
      else
        write(IMAIN,*) '     depth: ',dabs(z_found(irec) - elevation(irec)),' m'
        write(IMAIN,*) '     z: ',z_found(irec)
      endif
      write(IMAIN,*)

      ! add warning if estimate is poor
      ! (usually means receiver outside the mesh given by the user)
      if(final_distance(irec) > 3000.d0) then
        write(IMAIN,*) '*******************************************************'
        write(IMAIN,*) '***** WARNING: receiver location estimate is poor *****'
        write(IMAIN,*) '*******************************************************'
      endif

      write(IMAIN,*)
      endif

    enddo

    ! compute maximal distance for all the receivers
    final_distance_max = maxval(final_distance(:))

    ! display maximum error for all the receivers
    write(IMAIN,*) 'maximum error in location of all the receivers: ',sngl(final_distance_max),' m'

    ! add warning if estimate is poor
    ! (usually means receiver outside the mesh given by the user)
    if(final_distance_max > 1000.d0) then
      write(IMAIN,*)
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '***** WARNING: at least one receiver is poorly located *****'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
    endif

    ! write the locations of stations, so that we can load them and write them to SU headers later
    open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/output_list_stations.txt', &
              status='unknown',action='write',iostat=ios)
    if( ios /= 0 ) &
      call exit_mpi(myrank,'error opening file '//trim(OUTPUT_FILES)//'/output_list_stations.txt')

    do irec=1,nrec
      write(IOUT_SU,*) x_found(irec),y_found(irec),z_found(irec)
    enddo

    close(IOUT_SU)

    ! stores station infos for later runs
    if( SU_FORMAT ) then
      open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/SU_stations_info.bin', &
            status='unknown',action='write',form='unformatted',iostat=ios)
      if( ios == 0 ) then
        write(IOUT_SU) islice_selected_rec,ispec_selected_rec
        write(IOUT_SU) xi_receiver,eta_receiver,gamma_receiver
        write(IOUT_SU) x_found,y_found,z_found
        write(IOUT_SU) nu
        close(IOUT_SU)
      endif
    endif

    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of receiver detection - done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif    ! end of section executed by main process only

  ! main process broadcasts the results to all the slices
  call bcast_all_i(islice_selected_rec,nrec)
  call bcast_all_i(ispec_selected_rec,nrec)
  call bcast_all_dp(xi_receiver,nrec)
  call bcast_all_dp(eta_receiver,nrec)
  call bcast_all_dp(gamma_receiver,nrec)

  ! deallocate arrays
  deallocate(stlat)
  deallocate(stlon)
  deallocate(stele)
  deallocate(stbur)
  deallocate(stutm_x)
  deallocate(stutm_y)
  deallocate(horiz_dist)
  deallocate(ix_initial_guess)
  deallocate(iy_initial_guess)
  deallocate(iz_initial_guess)
  deallocate(x_target)
  deallocate(y_target)
  deallocate(z_target)
  deallocate(x_found)
  deallocate(y_found)
  deallocate(z_found)
  deallocate(final_distance)
  deallocate(ispec_selected_rec_all)
  deallocate(xi_receiver_all)
  deallocate(eta_receiver_all)
  deallocate(gamma_receiver_all)
  deallocate(x_found_all)
  deallocate(y_found_all)
  deallocate(z_found_all)
  deallocate(final_distance_all)

  ! synchronize all the processes to make sure everybody has finished
  call sync_all()

  end subroutine locate_receivers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine station_filter(SUPPRESS_UTM_PROJECTION,UTM_PROJECTION_ZONE, &
                            myrank,filename,filtered_filename,nfilter, &
                            LATITUDE_MIN, LATITUDE_MAX, LONGITUDE_MIN, LONGITUDE_MAX)

  implicit none

  include 'constants.h'

! input
  logical :: SUPPRESS_UTM_PROJECTION
  integer :: UTM_PROJECTION_ZONE
  integer :: myrank
  character(len=*) :: filename,filtered_filename
  double precision :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

! output
  integer :: nfilter

  integer :: nrec, nrec_filtered, ios

  double precision :: stlat,stlon,stele,stbur,stutm_x,stutm_y
  double precision :: minlat,minlon,maxlat,maxlon
  character(len=MAX_LENGTH_STATION_NAME) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name
  character(len=256) :: dummystring

  nrec = 0
  nrec_filtered = 0

  ! counts number of lines in stations file
  open(unit=IIN, file=trim(filename), status = 'old', iostat = ios)
  if (ios /= 0) call exit_mpi(myrank, 'No file '//trim(filename)//', exit')
  do while(ios == 0)
    read(IIN,"(a256)",iostat = ios) dummystring
    if(ios /= 0) exit

    if( len_trim(dummystring) > 0 ) nrec = nrec + 1
  enddo
  close(IIN)

  ! reads in station locations
  open(unit=IIN, file=trim(filename), status = 'old', iostat = ios)
  do while(ios == 0)
    read(IIN,"(a256)",iostat = ios) dummystring
    if( ios /= 0 ) exit

    ! counts number of stations in min/max region
    if( len_trim(dummystring) > 0 ) then
        dummystring = trim(dummystring)
        read(dummystring, *) station_name, network_name, stlat, stlon, stele, stbur

        ! convert station location to UTM
        call utm_geo(stlon,stlat,stutm_x,stutm_y,&
             UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

        ! counts stations within lon/lat region
        if( stutm_y >= LATITUDE_MIN .and. stutm_y <= LATITUDE_MAX .and. &
           stutm_x >= LONGITUDE_MIN .and. stutm_x <= LONGITUDE_MAX) &
          nrec_filtered = nrec_filtered + 1
     endif
  enddo
  close(IIN)

  ! writes out filtered stations file
  if (myrank == 0) then
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ios)
    open(unit=IOUT,file=trim(filtered_filename),status='unknown')
    do while(ios == 0)
      read(IIN,"(a256)",iostat = ios) dummystring
      if( ios /= 0 ) exit

      !read(IIN,*) station_name,network_name,stlat,stlon,stele,stbur
      if( len_trim(dummystring) > 0 ) then
        dummystring = trim(dummystring)
        read(dummystring, *) station_name, network_name, stlat, stlon, stele, stbur

        ! convert station location to UTM
        call utm_geo(stlon,stlat,stutm_x,stutm_y,&
             UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

        if( stutm_y >= LATITUDE_MIN .and. stutm_y <= LATITUDE_MAX .and. &
           stutm_x >= LONGITUDE_MIN .and. stutm_x <= LONGITUDE_MAX) then

          ! with specific format
          write(IOUT,'(a10,1x,a10,4e18.6)') &
                            trim(station_name),trim(network_name), &
                            sngl(stlat),sngl(stlon),sngl(stele),sngl(stbur)

       endif
      endif
    enddo
    close(IIN)
    close(IOUT)

    write(IMAIN,*)
    write(IMAIN,*) 'there are ',nrec,' stations in file ', trim(filename)
    write(IMAIN,*) 'saving ',nrec_filtered,' stations inside the model in file ', trim(filtered_filename)
    write(IMAIN,*) 'excluding ',nrec - nrec_filtered,' stations located outside the model'
    write(IMAIN,*)

    if( nrec_filtered < 1 ) then
      write(IMAIN,*) 'error filtered stations:'
      write(IMAIN,*) '  simulation needs at least 1 station but got ',nrec_filtered
      write(IMAIN,*)
      write(IMAIN,*) '  check that stations in file '//trim(filename)//' are within'

      if( SUPPRESS_UTM_PROJECTION ) then
        write(IMAIN,*) '    latitude min/max : ',LATITUDE_MIN,LATITUDE_MAX
        write(IMAIN,*) '    longitude min/max: ',LONGITUDE_MIN,LONGITUDE_MAX
      else
        ! convert edge locations from UTM back to lat/lon
        call utm_geo(minlon,minlat,LONGITUDE_MIN,LATITUDE_MIN,&
             UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)
        call utm_geo(maxlon,maxlat,LONGITUDE_MAX,LATITUDE_MAX,&
             UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)
        write(IMAIN,*) '    longitude min/max: ',minlon,maxlon
        write(IMAIN,*) '    latitude min/max : ',minlat,maxlat
        write(IMAIN,*) '    UTM x min/max: ',LONGITUDE_MIN,LONGITUDE_MAX
        write(IMAIN,*) '    UTM y min/max : ',LATITUDE_MIN,LATITUDE_MAX
      endif

      write(IMAIN,*)
    endif

  endif

  nfilter = nrec_filtered

  end subroutine station_filter

