!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

! for elastic solver

  subroutine compute_add_sources_viscoelastic(NSPEC_AB,NGLOB_AB,accel, &
                        ibool, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
                        sourcearrays, &
                        ispec_is_elastic,SIMULATION_TYPE,NSTEP, &
                        nrec,islice_selected_rec,ispec_selected_rec, &
                        nadj_rec_local,adj_sourcearrays, &
                        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY)

  use constants
  use specfem_par, only: xigll,yigll,zigll,xi_receiver,eta_receiver,gamma_receiver, &
                        station_name,network_name,adj_source_file, &
                        num_free_surface_faces,free_surface_ispec, &
                        free_surface_ijk,free_surface_jacobian2Dw, &
                        noise_sourcearray,irec_master_noise, &
                        normal_x_noise,normal_y_noise,normal_z_noise, &
                        mask_noise,noise_surface_movie, &
                        nrec_local,number_receiver_global, &
                        nsources_local,tshift_src,dt,t0,SU_FORMAT, &
                        USE_LDDRK,istage, &
                        EXTERNAL_STF,user_source_time_function

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_add_sources_viscoelastic_1.F90"
#endif

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! source
  integer :: NSOURCES,myrank,it
  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
  real(kind=CUSTOM_REAL), dimension(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrays

  logical, dimension(NSPEC_AB) :: ispec_is_elastic

!adjoint simulations
  integer:: SIMULATION_TYPE,NSTEP
  integer:: nrec
  integer,dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  integer:: nadj_rec_local
  logical :: ibool_read_adj_arrays
  integer :: it_sub_adj,itime,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY
  real(kind=CUSTOM_REAL),dimension(nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ):: &
    adj_sourcearrays

! local parameters
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),allocatable:: adj_sourcearray
  real(kind=CUSTOM_REAL) stf_used

  double precision :: stf,time_source_dble
  double precision,external :: get_stf_viscoelastic

  integer :: isource,iglob,i,j,k,ispec
  integer :: irec_local,irec, ier

! adjoint sources in SU format
  integer :: it_start,it_end
  real(kind=CUSTOM_REAL) :: adj_temp(NSTEP)
  real(kind=CUSTOM_REAL) :: adj_src(NTSTEP_BETWEEN_READ_ADJSRC,NDIM)
  character(len=MAX_STRING_LEN) :: procname
  integer,parameter :: nheader=240      ! 240 bytes
  !integer(kind=2) :: i2head(nheader/2)  ! 2-byte-integer
  !integer(kind=4) :: i4head(nheader/4)  ! 4-byte-integer
  real(kind=4)    :: r4head(nheader/4)  ! 4-byte-real
  !equivalence (i2head,i4head,r4head)    ! share the same 240-byte-memory
  double precision :: hxir(NGLLX),hpxir(NGLLX),hetar(NGLLY),hpetar(NGLLY),hgammar(NGLLZ),hpgammar(NGLLZ)

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_add_sources_viscoelastic_2.F90"
#endif

  ! forward simulations
  if (SIMULATION_TYPE == 1 .and. NOISE_TOMOGRAPHY == 0 .and. nsources_local > 0) then

    do isource = 1,NSOURCES

      !   add the source (only if this proc carries the source)
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        ! debug
        !if (it == 1) print *,'debug',myrank, it,'source',isource, ispec,ispec_is_elastic(ispec)

        if (ispec_is_elastic(ispec)) then
          ! current time
          if (USE_LDDRK) then
            time_source_dble = dble(it-1)*DT + dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
          else
            time_source_dble = dble(it-1)*DT - t0 - tshift_src(isource)
          endif

          ! determines source time function value
          stf = get_stf_viscoelastic(time_source_dble,isource)

          !! VM VM add external source time function
          if (EXTERNAL_STF) then
            stf = user_source_time_function(it, isource)
          endif

          ! distinguishes between single and double precision for reals
          stf_used = real(stf,kind=CUSTOM_REAL)

          ! adds source array
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                accel(:,iglob) = accel(:,iglob) + sourcearrays(isource,:,i,j,k)*stf_used
              enddo
            enddo
          enddo

        endif ! ispec_is_elastic
      endif ! myrank
    enddo ! NSOURCES
  endif ! forward

! NOTE: adjoint sources and backward wavefield timing:
!             idea is to start with the backward field b_displ,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_displ( it ) would correspond to (NSTEP - it - 1)*DT - t0
!       if we read in saved wavefields b_displ() before Newmark time scheme
!       (see sources for simulation_type 1 and seismograms)
!       since at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_displ( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0
!
!       b_displ is now read in after Newmark time scheme:
!       we read the backward/reconstructed wavefield at the end of the first time loop,
!       such that b_displ(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!       assuming that until that end the backward/reconstructed wavefield and adjoint fields
!       have a zero contribution to adjoint kernels.
!       thus the correct indexing is NSTEP - it + 1, instead of NSTEP - it
!
! adjoint wavefields:
!       since the adjoint source traces were derived from the seismograms,
!       it follows that for the adjoint wavefield, the time equivalent to ( T - t ) uses the time-reversed
!       adjoint source traces which start at -t0 and end at time (NSTEP-1)*DT - t0
!       for step it=1: (NSTEP -it + 1)*DT - t0 for backward wavefields corresponds to time T

! adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then

    ! adds adjoint source in this partitions
    if (nadj_rec_local > 0) then

      ! read in adjoint sources block by block (for memory consideration)
      ! e.g., in exploration experiments, both the number of receivers (nrec) and
      ! the number of time steps (NSTEP) are huge,
      ! which may cause problems since we have a large array:
      !   adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLY,NGLLZ)

      ! figure out if we need to read in a chunk of the adjoint source at this timestep
      it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )   !chunk_number
      ibool_read_adj_arrays = (((mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0)) .and. (nadj_rec_local > 0))

      ! needs to read in a new chunk/block of the adjoint source
      ! note that for each partition, we divide it into two parts --- boundaries and interior --- indicated by 'iphase'
      ! we first do calculations for the boudaries, and then start communication
      ! with other partitions while calculate for the inner part
      ! this must be done carefully, otherwise the adjoint sources may be added twice
      if (ibool_read_adj_arrays) then

        ! allocates temporary source array
        allocate(adj_sourcearray(NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
        if (ier /= 0) stop 'error allocating array adj_sourcearray'

        if (.not. SU_FORMAT) then
          !!! read ascii adjoint sources
          irec_local = 0
          do irec = 1, nrec
            ! compute source arrays
            if (myrank == islice_selected_rec(irec)) then
              irec_local = irec_local + 1
              ! reads in **net**.**sta**.**BH**.adj files
              adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
              call compute_arrays_adjoint_source(myrank,adj_source_file, &
                                                 xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec), &
                                                 adj_sourcearray, xigll,yigll,zigll, &
                                                 it_sub_adj,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC)

              do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
                adj_sourcearrays(irec_local,itime,:,:,:,:) = adj_sourcearray(itime,:,:,:,:)
              enddo
            endif
          enddo
        else
          !!! read SU adjoint sources
          ! range of the block we need to read
          it_start = NSTEP - it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC + 1
          it_end   = it_start + NTSTEP_BETWEEN_READ_ADJSRC - 1
          write(procname,"(i4)") myrank
          procname = adjustl(procname)
          ! read adjoint sources
          open(unit=IIN_SU1, file=trim(OUTPUT_FILES)//'../SEM/'//trim(procname)//'_dx_SU.adj', &
               status='old', access='direct', recl=240+4*NSTEP, iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'file '//trim(OUTPUT_FILES) &
                                    //'../SEM/'//trim(procname)//'_dx_SU.adj does not exist')
          open(unit=IIN_SU2, file=trim(OUTPUT_FILES)//'../SEM/'//trim(procname)//'_dy_SU.adj', &
               status='old', access='direct', recl=240+4*NSTEP, iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'file '//trim(OUTPUT_FILES) &
                                    //'../SEM/'//trim(procname)//'_dy_SU.adj does not exist')
          open(unit=IIN_SU3, file=trim(OUTPUT_FILES)//'../SEM/'//trim(procname)//'_dz_SU.adj', &
               status='old', access='direct', recl=240+4*NSTEP, iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'file '//trim(OUTPUT_FILES) &
                                    //'../SEM/'//trim(procname)//'_dz_SU.adj does not exist')

          do irec_local = 1,nrec_local
            irec = number_receiver_global(irec_local)
            read(IIN_SU1,rec=irec_local) r4head, adj_temp
            adj_src(:,1)=adj_temp(it_start:it_end)
            read(IIN_SU2,rec=irec_local) r4head, adj_temp
            adj_src(:,2)=adj_temp(it_start:it_end)
            read(IIN_SU3,rec=irec_local) r4head, adj_temp
            adj_src(:,3)=adj_temp(it_start:it_end)
            ! lagrange interpolators for receiver location
            call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
            call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
            call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
            ! interpolates adjoint source onto GLL points within this element
            do k = 1, NGLLZ
              do j = 1, NGLLY
                do i = 1, NGLLX
                  adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
                enddo
              enddo
            enddo
            do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
              adj_sourcearrays(irec_local,itime,:,:,:,:) = adj_sourcearray(itime,:,:,:,:)
            enddo
          enddo
          close(IIN_SU1)
          close(IIN_SU2)
          close(IIN_SU3)
        endif !if (.not. SU_FORMAT)

        deallocate(adj_sourcearray)
      endif ! if (ibool_read_adj_arrays)


      if (it < NSTEP) then
        ! receivers act as sources
        irec_local = 0
        do irec = 1,nrec

          ! add the source (only if this proc carries the source)
          if (myrank == islice_selected_rec(irec)) then
            irec_local = irec_local + 1

            ispec = ispec_selected_rec(irec)
            if (ispec_is_elastic(ispec)) then

              ! adds source array
              do k = 1,NGLLZ
                do j = 1,NGLLY
                  do i = 1,NGLLX
                    iglob = ibool(i,j,k,ispec_selected_rec(irec))

                    accel(:,iglob) = accel(:,iglob)  &
                              + adj_sourcearrays(irec_local, &
                              NTSTEP_BETWEEN_READ_ADJSRC - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC), &
                              :,i,j,k)
                  enddo
                enddo
              enddo
            endif ! ispec_is_elastic
          endif
        enddo ! nrec
      endif ! it
    endif ! nadj_rec_local
  endif !adjoint

  ! for noise simulations
  if (NOISE_TOMOGRAPHY > 0) then
    ! we have two loops indicated by iphase ("inner elements/points" or "boundary elements/points")
    ! here, we add all noise sources once, when we are calculating for boundary points (iphase==1),
    ! because boundary points are calculated first!
    if (NOISE_TOMOGRAPHY == 1) then
      ! the first step of noise tomography is to use |S(\omega)|^2 as a point force source at one of the receivers.
      ! hence, instead of a moment tensor 'sourcearrays', a 'noise_sourcearray' for a point force is needed.
      ! furthermore, the CMTSOLUTION needs to be zero, i.e., no earthquakes.
      ! now this must be manually set in DATA/CMTSOLUTION, by USERS.
      call add_source_master_rec_noise(myrank,nrec, &
              NSTEP,accel,noise_sourcearray, &
              ibool,islice_selected_rec,ispec_selected_rec, &
              it,irec_master_noise, &
              NSPEC_AB,NGLOB_AB)
    else if (NOISE_TOMOGRAPHY == 2) then
      ! second step of noise tomography, i.e., read the surface movie saved at every timestep
      ! use the movie to drive the ensemble forward wavefield
      call noise_read_add_surface_movie(NGLLSQUARE*num_free_surface_faces,accel, &
                             normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                             ibool,noise_surface_movie,NSTEP-it+1,NSPEC_AB,NGLOB_AB, &
                             num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                             free_surface_jacobian2Dw)
      ! be careful, since ensemble forward sources are reversals of generating wavefield "eta"
      ! hence the "NSTEP-it+1", i.e., start reading from the last timestep
      ! note the ensemble forward sources are generally distributed on the surface of the earth
      ! that's to say, the ensemble forward source is kind of a surface force density, not a body force density
      ! therefore, we must add it here, before applying the inverse of mass matrix
    endif
  endif

  end subroutine compute_add_sources_viscoelastic
!
!=====================================================================
! for elastic solver

  subroutine compute_add_sources_viscoelastic_backward( NSPEC_AB,NGLOB_AB, &
                        ibool, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
                        sourcearrays, &
                        ispec_is_elastic,SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                        b_accel,NOISE_TOMOGRAPHY)

  use constants
  use specfem_par, only: num_free_surface_faces,free_surface_ispec, &
                        free_surface_ijk,free_surface_jacobian2Dw, &
                        normal_x_noise,normal_y_noise,normal_z_noise, &
                        mask_noise,noise_surface_movie, &
                        nsources_local,tshift_src,dt,t0, &
                        USE_LDDRK,istage,EXTERNAL_STF,user_source_time_function

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_add_sources_viscoelastic_1.F90"
#endif

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! source
  integer :: NSOURCES,myrank,it
  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
  real(kind=CUSTOM_REAL), dimension(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrays

  logical, dimension(NSPEC_AB) :: ispec_is_elastic

!adjoint simulations
  integer:: SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_ADJOINT):: b_accel
  integer :: NOISE_TOMOGRAPHY

! local parameters
  real(kind=CUSTOM_REAL) stf_used

  double precision :: stf,time_source_dble
  double precision,external :: get_stf_viscoelastic

  integer :: isource,iglob,i,j,k,ispec

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_add_sources_viscoelastic_2.F90"
#endif

! NOTE: adjoint sources and backward wavefield timing:
!             idea is to start with the backward field b_displ,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_displ( it ) would correspond to (NSTEP - it - 1)*DT - t0
!       if we read in saved wavefields b_displ() before Newmark time scheme
!       (see sources for simulation_type 1 and seismograms)
!       since at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_displ( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0
!
!       b_displ is now read in after Newmark time scheme:
!       we read the backward/reconstructed wavefield at the end of the first time loop,
!       such that b_displ(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!       assuming that until that end the backward/reconstructed wavefield and adjoint fields
!       have a zero contribution to adjoint kernels.
!       thus the correct indexing is NSTEP - it + 1, instead of NSTEP - it
!
! adjoint wavefields:
!       since the adjoint source traces were derived from the seismograms,
!       it follows that for the adjoint wavefield, the time equivalent to ( T - t ) uses the time-reversed
!       adjoint source traces which start at -t0 and end at time (NSTEP-1)*DT - t0
!       for step it=1: (NSTEP -it + 1)*DT - t0 for backward wavefields corresponds to time T

! adjoint simulations
  if (SIMULATION_TYPE == 3 .and. NOISE_TOMOGRAPHY == 0 .and. nsources_local > 0) then

    ! backward source reconstruction
    do isource = 1,NSOURCES

      ! add the source (only if this proc carries the source)
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        if (ispec_is_elastic(ispec)) then
          ! current time
          if (USE_LDDRK) then
            time_source_dble = dble(NSTEP-it)*DT - dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
          else
            time_source_dble = dble(NSTEP-it)*DT - t0 - tshift_src(isource)
          endif

          ! determines source time function value
          stf = get_stf_viscoelastic(time_source_dble,isource)

          !! VM VM add external source time function
          if (EXTERNAL_STF) then
            ! time-reversed
            stf = user_source_time_function(NSTEP-it+1, isource)
          endif

          ! distinguishes between single and double precision for reals
          stf_used = real(stf,kind=CUSTOM_REAL)

          ! adds source
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                b_accel(:,iglob) = b_accel(:,iglob) + sourcearrays(isource,:,i,j,k) * stf_used
              enddo
            enddo
          enddo

        endif ! elastic
      endif ! myrank
    enddo ! NSOURCES
  endif ! adjoint

  ! for noise simulations
  if (NOISE_TOMOGRAPHY > 0) then
    ! we have two loops indicated by iphase ("inner elements/points" or "boundary elements/points")
    ! here, we add all noise sources once, when we are calculating for boundary points (iphase==1),
    ! because boundary points are calculated first!
    if (NOISE_TOMOGRAPHY == 3) then
      ! third step of noise tomography, i.e., read the surface movie saved at every timestep
      ! use the movie to reconstruct the ensemble forward wavefield
      ! the ensemble adjoint wavefield is done as usual
      ! note instead of "NSTEP-it+1", now we us "it", since reconstruction is a reversal of reversal
      call noise_read_add_surface_movie(NGLLSQUARE*num_free_surface_faces,b_accel, &
                            normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                            ibool,noise_surface_movie,it,NSPEC_AB,NGLOB_AB, &
                            num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                            free_surface_jacobian2Dw)
    endif
  endif

  end subroutine compute_add_sources_viscoelastic_backward

!
!=====================================================================
! for elastic solver on GPU

  subroutine compute_add_sources_viscoelastic_GPU(NSPEC_AB, &
                                                  NSOURCES,myrank,it, &
                                                  ispec_is_elastic,SIMULATION_TYPE,NSTEP, &
                                                  nrec,islice_selected_rec,ispec_selected_rec, &
                                                  nadj_rec_local,adj_sourcearrays, &
                                                  NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                                                  Mesh_pointer)

  use constants
  use specfem_par, only: xigll,yigll,zigll,xi_receiver,eta_receiver,gamma_receiver, &
                        station_name,network_name,adj_source_file, &
                        num_free_surface_faces, &
                        irec_master_noise,noise_surface_movie, &
                        nrec_local,number_receiver_global, &
                        nsources_local,tshift_src,dt,t0,SU_FORMAT, &
                        USE_LDDRK,istage,EXTERNAL_STF,user_source_time_function

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_add_sources_viscoelastic_1.F90"
#endif

  implicit none

  integer :: NSPEC_AB

! arrays with mesh parameters per slice

! source
  integer :: NSOURCES,myrank,it

  logical, dimension(NSPEC_AB) :: ispec_is_elastic

  ! adjoint simulations
  integer:: SIMULATION_TYPE,NSTEP
  integer(kind=8) :: Mesh_pointer
  integer:: nrec
  integer,dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  integer:: nadj_rec_local
  logical :: ibool_read_adj_arrays
  integer :: it_sub_adj,itime,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY
  real(kind=CUSTOM_REAL),dimension(nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ):: &
    adj_sourcearrays

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),allocatable:: adj_sourcearray

  double precision :: stf,time_source_dble
  double precision,external :: get_stf_viscoelastic

  ! for GPU_MODE
  double precision, dimension(NSOURCES) :: stf_pre_compute

  integer :: isource,i,j,k
  integer :: irec_local,irec, ier

  ! adjoint sources in SU format
  integer :: it_start,it_end
  real(kind=CUSTOM_REAL) :: adj_temp(NSTEP)
  real(kind=CUSTOM_REAL) :: adj_src(NTSTEP_BETWEEN_READ_ADJSRC,NDIM)
  character(len=MAX_STRING_LEN) :: procname
  integer,parameter :: nheader=240      ! 240 bytes
  !integer(kind=2) :: i2head(nheader/2)  ! 2-byte-integer
  !integer(kind=4) :: i4head(nheader/4)  ! 4-byte-integer
  real(kind=4)    :: r4head(nheader/4)  ! 4-byte-real
  !equivalence (i2head,i4head,r4head)    ! share the same 240-byte-memory
  double precision :: hxir(NGLLX),hpxir(NGLLX),hetar(NGLLY),hpetar(NGLLY),hgammar(NGLLZ),hpgammar(NGLLZ)

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_add_sources_viscoelastic_2.F90"
#endif

  ! forward simulations
  if (SIMULATION_TYPE == 1 .and. NOISE_TOMOGRAPHY == 0 .and. nsources_local > 0) then

    if (NSOURCES > 0) then
      do isource = 1,NSOURCES
        ! current time
        if (USE_LDDRK) then
          time_source_dble = dble(it-1)*DT + dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
        else
          time_source_dble = dble(it-1)*DT - t0 - tshift_src(isource)
        endif

        ! determines source time function value
        stf = get_stf_viscoelastic(time_source_dble,isource)

        !! VM VM add external source time function
        if (EXTERNAL_STF) then
           stf = user_source_time_function(it, isource)
        endif

        ! stores precomputed source time function factor
        stf_pre_compute(isource) = stf
      enddo

      ! only implements SIMTYPE=1 and NOISE_TOM=0
      ! write(*,*) "fortran dt = ", dt
      ! change dt -> DT
      call compute_add_sources_el_cuda(Mesh_pointer,stf_pre_compute,NSOURCES)
    endif
  endif ! forward

! NOTE: adjoint sources and backward wavefield timing:
!             idea is to start with the backward field b_displ,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_displ( it ) would correspond to (NSTEP - it - 1)*DT - t0
!       if we read in saved wavefields b_displ() before Newmark time scheme
!       (see sources for simulation_type 1 and seismograms)
!       since at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_displ( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0
!
!       b_displ is now read in after Newmark time scheme:
!       we read the backward/reconstructed wavefield at the end of the first time loop,
!       such that b_displ(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!       assuming that until that end the backward/reconstructed wavefield and adjoint fields
!       have a zero contribution to adjoint kernels.
!       thus the correct indexing is NSTEP - it + 1, instead of NSTEP - it
!
! adjoint wavefields:
!       since the adjoint source traces were derived from the seismograms,
!       it follows that for the adjoint wavefield, the time equivalent to ( T - t ) uses the time-reversed
!       adjoint source traces which start at -t0 and end at time (NSTEP-1)*DT - t0
!       for step it=1: (NSTEP -it + 1)*DT - t0 for backward wavefields corresponds to time T

! adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! adds adjoint source in this partitions
    if (nadj_rec_local > 0) then

      ! read in adjoint sources block by block (for memory consideration)
      ! e.g., in exploration experiments, both the number of receivers (nrec) and
      ! the number of time steps (NSTEP) are huge,
      ! which may cause problems since we have a large array:
      !   adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLY,NGLLZ)

      ! figure out if we need to read in a chunk of the adjoint source at this timestep
      it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )   !chunk_number
      ibool_read_adj_arrays = (((mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0)) .and. (nadj_rec_local > 0))

      ! needs to read in a new chunk/block of the adjoint source
      ! note that for each partition, we divide it into two parts --- boundaries and interior --- indicated by 'iphase'
      ! we first do calculations for the boudaries, and then start communication
      ! with other partitions while calculate for the inner part
      ! this must be done carefully, otherwise the adjoint sources may be added twice
      if (ibool_read_adj_arrays) then

        ! allocates temporary source array
        allocate(adj_sourcearray(NTSTEP_BETWEEN_READ_ADJSRC,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
        if (ier /= 0) stop 'error allocating array adj_sourcearray'

        if (.not. SU_FORMAT) then
          !!! read ascii adjoint sources
          irec_local = 0
          do irec = 1, nrec
            ! compute source arrays
            if (myrank == islice_selected_rec(irec)) then
              irec_local = irec_local + 1
              ! reads in **net**.**sta**.**BH**.adj files
              adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
              call compute_arrays_adjoint_source(myrank,adj_source_file, &
                                                 xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec), &
                                                 adj_sourcearray, xigll,yigll,zigll, &
                                                 it_sub_adj,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC)

              do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
                adj_sourcearrays(irec_local,itime,:,:,:,:) = adj_sourcearray(itime,:,:,:,:)
              enddo
            endif
          enddo
        else
          !!! read SU adjoint sources
          ! range of the block we need to read
          it_start = NSTEP - it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC + 1
          it_end   = it_start + NTSTEP_BETWEEN_READ_ADJSRC - 1
          write(procname,"(i4)") myrank
          procname = adjustl(procname)
          ! read adjoint sources
          open(unit=IIN_SU1, file=trim(OUTPUT_FILES)//'../SEM/'//trim(procname)//'_dx_SU.adj', &
               status='old', access='direct', recl=240+4*NSTEP, iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'file '//trim(OUTPUT_FILES) &
                                    //'../SEM/'//trim(procname)//'_dx_SU.adj does not exist')
          open(unit=IIN_SU2, file=trim(OUTPUT_FILES)//'../SEM/'//trim(procname)//'_dy_SU.adj', &
               status='old', access='direct', recl=240+4*NSTEP, iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'file '//trim(OUTPUT_FILES) &
                                    //'../SEM/'//trim(procname)//'_dy_SU.adj does not exist')
          open(unit=IIN_SU3, file=trim(OUTPUT_FILES)//'../SEM/'//trim(procname)//'_dz_SU.adj', &
               status='old', access='direct', recl=240+4*NSTEP, iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'file '//trim(OUTPUT_FILES) &
                                    //'../SEM/'//trim(procname)//'_dz_SU.adj does not exist')

          do irec_local = 1,nrec_local
            irec = number_receiver_global(irec_local)
            read(IIN_SU1,rec=irec_local) r4head, adj_temp
            adj_src(:,1)=adj_temp(it_start:it_end)
            read(IIN_SU2,rec=irec_local) r4head, adj_temp
            adj_src(:,2)=adj_temp(it_start:it_end)
            read(IIN_SU3,rec=irec_local) r4head, adj_temp
            adj_src(:,3)=adj_temp(it_start:it_end)
            ! lagrange interpolators for receiver location
            call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
            call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
            call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
            ! interpolates adjoint source onto GLL points within this element
            do k = 1, NGLLZ
              do j = 1, NGLLY
                do i = 1, NGLLX
                  adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src(:,:)
                enddo
              enddo
            enddo
            do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
              adj_sourcearrays(irec_local,itime,:,:,:,:) = adj_sourcearray(itime,:,:,:,:)
            enddo
          enddo
          close(IIN_SU1)
          close(IIN_SU2)
          close(IIN_SU3)
        endif !if (.not. SU_FORMAT)

        deallocate(adj_sourcearray)
      endif ! if (ibool_read_adj_arrays)


      if (it < NSTEP) then
        call add_sources_el_sim_type_2_or_3(Mesh_pointer,adj_sourcearrays, &
                                            ispec_is_elastic, &
                                            ispec_selected_rec, &
                                            nrec, &
                                            NTSTEP_BETWEEN_READ_ADJSRC - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC), &
                                            islice_selected_rec,nadj_rec_local, &
                                            NTSTEP_BETWEEN_READ_ADJSRC)
      endif ! it
    endif ! nadj_rec_local
  endif !adjoint

! note:  b_displ() is read in after Newmark time scheme, thus
!           b_displ(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!           thus indexing is NSTEP - it , instead of NSTEP - it - 1

! adjoint simulations
  if (SIMULATION_TYPE == 3 .and. NOISE_TOMOGRAPHY == 0 .and. nsources_local > 0) then

    if (NSOURCES > 0) then
      do isource = 1,NSOURCES
        ! current time
        if (USE_LDDRK) then
          time_source_dble = dble(NSTEP-it)*DT - dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
        else
          time_source_dble = dble(NSTEP-it)*DT - t0 - tshift_src(isource)
        endif

        ! determines source time function value
        stf = get_stf_viscoelastic(time_source_dble,isource)

        !! VM VM add external source time function
        if (EXTERNAL_STF) then
           stf = user_source_time_function(NSTEP-it+1, isource)
        endif

        ! stores precomputed source time function factor
        stf_pre_compute(isource) = stf
      enddo

      ! only implements SIMTYPE=3
      call compute_add_sources_el_s3_cuda(Mesh_pointer,stf_pre_compute,NSOURCES)
    endif
  endif ! adjoint

  ! for noise simulations
  if (NOISE_TOMOGRAPHY > 0) then
    ! we have two loops indicated by iphase ("inner elements/points" or "boundary elements/points")
    ! here, we add all noise sources once, when we are calculating for boundary points (iphase==1),
    ! because boundary points are calculated first!
    if (NOISE_TOMOGRAPHY == 1) then
      ! the first step of noise tomography is to use |S(\omega)|^2 as a point force source at one of the receivers.
      ! hence, instead of a moment tensor 'sourcearrays', a 'noise_sourcearray' for a point force is needed.
      ! furthermore, the CMTSOLUTION needs to be zero, i.e., no earthquakes.
      ! now this must be manually set in DATA/CMTSOLUTION, by USERS.
      call add_source_master_rec_noise_cu(Mesh_pointer,it,irec_master_noise,islice_selected_rec)
    else if (NOISE_TOMOGRAPHY == 2) then
      ! second step of noise tomography, i.e., read the surface movie saved at every timestep
      ! use the movie to drive the ensemble forward wavefield
      call noise_read_add_surface_movie_GPU(noise_surface_movie,NSTEP-it+1,num_free_surface_faces, &
                                            Mesh_pointer,NOISE_TOMOGRAPHY)
      ! be careful, since ensemble forward sources are reversals of generating wavefield "eta"
      ! hence the "NSTEP-it+1", i.e., start reading from the last timestep
      ! note the ensemble forward sources are generally distributed on the surface of the earth
      ! that's to say, the ensemble forward source is kind of a surface force density, not a body force density
      ! therefore, we must add it here, before applying the inverse of mass matrix
    else if (NOISE_TOMOGRAPHY == 3) then
      ! third step of noise tomography, i.e., read the surface movie saved at every timestep
      ! use the movie to reconstruct the ensemble forward wavefield
      ! the ensemble adjoint wavefield is done as usual
      ! note instead of "NSTEP-it+1", now we us "it", since reconstruction is a reversal of reversal
      call noise_read_add_surface_movie_GPU(noise_surface_movie,it,num_free_surface_faces, &
                                            Mesh_pointer,NOISE_TOMOGRAPHY)
    endif
  endif

  end subroutine compute_add_sources_viscoelastic_GPU


!
!=====================================================================
!

  double precision function get_stf_viscoelastic(time_source_dble,isource)

! returns source time function value for specified time

  use specfem_par, only: USE_FORCE_POINT_SOURCE,USE_RICKER_TIME_FUNCTION,hdur,hdur_Gaussian

  implicit none

  double precision,intent(in) :: time_source_dble
  integer,intent(in) :: isource

  ! local parameters
  double precision :: stf

  double precision, external :: comp_source_time_function,comp_source_time_function_rickr, &
    comp_source_time_function_gauss

  ! determines source time function value
  if (USE_FORCE_POINT_SOURCE) then
    ! single point force
    if (USE_RICKER_TIME_FUNCTION) then
      ! Ricker
      ! f0 has been stored in the hdur() array in the case of FORCESOLUTION,
      ! to use the same array as for CMTSOLUTION
      stf = comp_source_time_function_rickr(time_source_dble,hdur(isource))
    else
      ! Gaussian
      ! stf = comp_source_time_function_gauss(time_source_dble,5.d0*DT)
      !! COMMENTED BY FS FS -> do no longer use hard-coded hdur_Gaussian = 5*DT, but actual value of hdur_Gaussian

      stf = comp_source_time_function_gauss(time_source_dble,hdur_Gaussian(isource))
      !! ADDED BY FS FS -> use actual value of hdur_Gaussian as half duration
    endif
  else
    ! moment-tensor
    if (USE_RICKER_TIME_FUNCTION) then
      ! Ricker
      stf = comp_source_time_function_rickr(time_source_dble,hdur(isource))
    else
      ! Heaviside
      stf = comp_source_time_function(time_source_dble,hdur_Gaussian(isource))
    endif

    ! source encoding
    ! not supported yet for viscoelastic elements... sign of moment-tensor needs to be determined prior to running simulation
    !if (USE_SOURCE_ENCODING) stf = stf * pm1_source_encoding(isource)

  endif ! USE_FORCE_POINT_SOURCE

  ! return value
  get_stf_viscoelastic = stf

  end function get_stf_viscoelastic

