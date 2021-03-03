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

! for elastic solver

  subroutine compute_add_sources_viscoelastic()

  use constants
  use specfem_par, only: station_name,network_name,num_free_surface_faces,free_surface_ispec, &
                        free_surface_ijk,free_surface_jacobian2Dw, &
                        nsources_local,tshift_src,dt,t0,SU_FORMAT, &
                        USE_LDDRK,istage,USE_EXTERNAL_SOURCE_FILE,user_source_time_function, &
                        USE_BINARY_FOR_SEISMOGRAMS,NSPEC_AB,NGLOB_AB,ibool,NSOURCES,myrank,it,islice_selected_source, &
                        ispec_selected_source,sourcearrays,SIMULATION_TYPE,NSTEP,READ_ADJSRC_ASDF, &
                        nrec,islice_selected_rec,ispec_selected_rec, &
                        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                        hxir_adjstore,hetar_adjstore,hgammar_adjstore,source_adjoint,number_adjsources_global,nadj_rec_local, &
                        INVERSE_FWI_FULL_PROBLEM

  use specfem_par_elastic, only: accel,ispec_is_elastic

  ! noise
  use specfem_par_noise, only: noise_sourcearray,irec_main_noise, &
    normal_x_noise,normal_y_noise,normal_z_noise,mask_noise,noise_surface_movie

  ! coupling
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE

  ! faults
  use specfem_par, only: FAULT_SIMULATION

  implicit none

! local parameters
  real(kind=CUSTOM_REAL) :: stf_used,hlagrange
  logical :: ibool_read_adj_arrays
  double precision :: stf,time_source_dble,time_t
  double precision,external :: get_stf_viscoelastic

  integer :: isource,iglob,i,j,k,ispec,it_sub_adj
  integer :: irec_local,irec

  character(len=MAX_STRING_LEN) :: adj_source_file

  ! no source inside the mesh if we are coupling with DSM
  ! because the source is precisely the wavefield coming from the DSM traction file
  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE == 1) return

  ! sets current initial time
  if (USE_LDDRK) then
    ! LDDRK
    ! note: the LDDRK scheme updates displacement after the stiffness computations and
    !       after adding boundary/coupling/source terms.
    !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
    !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
    time_t = dble(it-1-1)*DT + dble(C_LDDRK(istage))*DT - t0
  else
    time_t = dble(it-1)*DT - t0
  endif

  ! forward simulations
  if (SIMULATION_TYPE == 1 .and. NOISE_TOMOGRAPHY == 0 .and. nsources_local > 0) then
    ! ignore CMT sources for fault rupture simulations
    if (FAULT_SIMULATION) return

! openmp solver
!$OMP PARALLEL if (NSOURCES > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(isource,time_source_dble,stf_used,stf,iglob,ispec,i,j,k)

    ! adds elastic sources
!$OMP DO
    do isource = 1,NSOURCES

      ! add the source (only if this proc carries the source)
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        if (ispec_is_elastic(ispec)) then
          ! current time
          time_source_dble = time_t - tshift_src(isource)

          !! add external source time function
          if (USE_EXTERNAL_SOURCE_FILE) then
             stf = user_source_time_function(it, isource)
          else
             ! determines source time function value
             stf = get_stf_viscoelastic(time_source_dble,isource)
          endif

          ! distinguishes between single and double precision for reals
          stf_used = real(stf,kind=CUSTOM_REAL)

          ! adds source array
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
!$OMP ATOMIC
                accel(1,iglob) = accel(1,iglob) + sourcearrays(isource,1,i,j,k)*stf_used
!$OMP ATOMIC
                accel(2,iglob) = accel(2,iglob) + sourcearrays(isource,2,i,j,k)*stf_used
!$OMP ATOMIC
                accel(3,iglob) = accel(3,iglob) + sourcearrays(isource,3,i,j,k)*stf_used
              enddo
            enddo
          enddo

        endif ! ispec_is_elastic
      endif ! myrank
    enddo ! NSOURCES
!$OMP ENDDO
!$OMP END PARALLEL

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
      if (ibool_read_adj_arrays .and. .not. INVERSE_FWI_FULL_PROBLEM) then
        ! reads adjoint source files
        if (.not. (SU_FORMAT .or. READ_ADJSRC_ASDF)) then
          ! ASCII formant
          if (USE_BINARY_FOR_SEISMOGRAMS) stop 'Adjoint simulations not supported with .bin format, please use SU format instead'
          !!! read ascii adjoint sources
          do irec_local = 1, nadj_rec_local
            ! reads in **net**.**sta**.**BH**.adj files
            irec = number_adjsources_global(irec_local)
            adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
            call compute_arrays_adjoint_source(adj_source_file,irec_local)
          enddo
        else if (READ_ADJSRC_ASDF) then
          ! ASDF format
          do irec_local = 1, nadj_rec_local
            ! reads in **net**.**sta**.**BH**.adj files
            irec = number_adjsources_global(irec_local)
            adj_source_file = trim(network_name(irec))//'_'//trim(station_name(irec))
            call compute_arrays_adjoint_source(adj_source_file,irec_local)
          enddo
          call compute_arrays_adjoint_source(adj_source_file, irec_local)
        else
          ! SU format
          call compute_arrays_adjoint_source_SU()
        endif !if (.not. SU_FORMAT)
      endif ! if (ibool_read_adj_arrays)

      ! adds source term
      if (it < NSTEP) then
        ! receivers act as sources
        do irec_local = 1, nadj_rec_local
          irec = number_adjsources_global(irec_local)
          ! element index
          ispec = ispec_selected_rec(irec)
          if (ispec_is_elastic(ispec)) then
            ! adds source array
            do k = 1,NGLLZ
              do j = 1,NGLLY
                do i = 1,NGLLX
                  iglob = ibool(i,j,k,ispec)

                  hlagrange = hxir_adjstore(i,irec_local) * hetar_adjstore(j,irec_local) * hgammar_adjstore(k,irec_local)

                  accel(:,iglob) = accel(:,iglob) &
                    + source_adjoint(:,irec_local,NTSTEP_BETWEEN_READ_ADJSRC - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC)) * hlagrange
                enddo
              enddo
            enddo
          endif ! ispec_is_elastic
        enddo
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
      call add_source_main_rec_noise(nrec,NSTEP,accel,noise_sourcearray, &
                                     ibool,islice_selected_rec,ispec_selected_rec, &
                                     it,irec_main_noise, &
                                     NSPEC_AB,NGLOB_AB)
    else if (NOISE_TOMOGRAPHY == 2) then
      ! second step of noise tomography, i.e., read the surface movie saved at every timestep
      ! use the movie to drive the ensemble forward wavefield
      call noise_read_add_surface_movie(NGLLSQUARE*num_free_surface_faces, &
                                        accel, &
                                        normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                                        ibool,noise_surface_movie, &
                                        NSTEP-it+1, &
                                        NSPEC_AB,NGLOB_AB, &
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

  subroutine compute_add_sources_viscoelastic_backward()

  use constants
  use specfem_par, only: num_free_surface_faces,free_surface_ispec, &
                        free_surface_ijk,free_surface_jacobian2Dw, &
                        nsources_local,tshift_src,dt,t0, &
                        USE_LDDRK,istage,USE_EXTERNAL_SOURCE_FILE,user_source_time_function, &
                        NSPEC_AB,NGLOB_AB,ibool, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
                        sourcearrays,SIMULATION_TYPE,NSTEP,NOISE_TOMOGRAPHY

  use specfem_par_elastic, only: b_accel,ispec_is_elastic

  ! noise
  use specfem_par_noise, only: normal_x_noise,normal_y_noise,normal_z_noise, &
    mask_noise,noise_surface_movie

  ! coupling
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE

  ! undo_att
  use specfem_par, only: UNDO_ATTENUATION_AND_OR_PML,NSUBSET_ITERATIONS,NT_DUMP_ATTENUATION, &
                         iteration_on_subset,it_of_this_subset

  ! faults
  use specfem_par, only: FAULT_SIMULATION

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) stf_used

  double precision :: stf,time_source_dble,time_t
  double precision,external :: get_stf_viscoelastic

  integer :: isource,iglob,i,j,k,ispec,it_tmp

  ! no source inside the mesh if we are coupling with DSM
  ! because the source is precisely the wavefield coming from the DSM traction file
  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE == 3) return

  ! ignore CMT sources for fault rupture simulations
  if (FAULT_SIMULATION) return

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return

  ! iteration step
  if (UNDO_ATTENUATION_AND_OR_PML) then
    ! example: NSTEP is a multiple of NT_DUMP_ATTENUATION
    !         NT_DUMP_ATTENUATION = 301, NSTEP = 1204, NSUBSET_ITERATIONS = 4, iteration_on_subset = 1 -> 4,
    !              1. subset, it_temp goes from 301 down to 1
    !              2. subset, it_temp goes from 602 down to 302
    !              3. subset, it_temp goes from 903 down to 603
    !              4. subset, it_temp goes from 1204 down to 904
    !valid for multiples only:
    !it_tmp = iteration_on_subset * NT_DUMP_ATTENUATION - it_of_this_subset + 1
    !
    ! example: NSTEP is **NOT** a multiple of NT_DUMP_ATTENUATION
    !          NT_DUMP_ATTENUATION = 301, NSTEP = 900, NSUBSET_ITERATIONS = 3, iteration_on_subset = 1 -> 3
    !              1. subset, it_temp goes from (900 - 602) = 298 down to 1
    !              2. subset, it_temp goes from (900 - 301) = 599 down to 299
    !              3. subset, it_temp goes from (900 - 0)   = 900 down to 600
    !works always:
    it_tmp = NSTEP - (NSUBSET_ITERATIONS - iteration_on_subset)*NT_DUMP_ATTENUATION - it_of_this_subset + 1
  else
    it_tmp = it
  endif

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

  ! sets current initial time
  if (USE_LDDRK) then
    ! LDDRK
    ! note: the LDDRK scheme updates displacement after the stiffness computations and
    !       after adding boundary/coupling/source terms.
    !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
    !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
    if (UNDO_ATTENUATION_AND_OR_PML) then
      ! stepping moves forward from snapshot position
      time_t = dble(NSTEP-it_tmp-1)*DT + dble(C_LDDRK(istage))*DT - t0
    else
      time_t = dble(NSTEP-it_tmp-1)*DT - dble(C_LDDRK(istage))*DT - t0
    endif
  else
    time_t = dble(NSTEP-it_tmp)*DT - t0
  endif

! adjoint simulations
  if (NOISE_TOMOGRAPHY == 0 .and. nsources_local > 0) then

    ! backward source reconstruction
    do isource = 1,NSOURCES

      ! add the source (only if this proc carries the source)
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        if (ispec_is_elastic(ispec)) then
          ! current time
          time_source_dble = time_t - tshift_src(isource)

          !! add external source time function
          if (USE_EXTERNAL_SOURCE_FILE) then
             ! time-reversed
             stf = user_source_time_function(NSTEP-it+1, isource)
          else
             ! determines source time function value
             stf = get_stf_viscoelastic(time_source_dble,isource)
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
      call noise_read_add_surface_movie(NGLLSQUARE*num_free_surface_faces, &
                                        b_accel, &
                                        normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                                        ibool,noise_surface_movie, &
                                        it, &
                                        NSPEC_AB,NGLOB_AB, &
                                        num_free_surface_faces,free_surface_ispec,free_surface_ijk, &
                                        free_surface_jacobian2Dw)
    endif
  endif

  end subroutine compute_add_sources_viscoelastic_backward

!
!=====================================================================
! for elastic solver on GPU

  subroutine compute_add_sources_viscoelastic_GPU()

  use constants
  use specfem_par, only: station_name,network_name, &
                        num_free_surface_faces, &
                        nsources_local,tshift_src,dt,t0,SU_FORMAT, &
                        USE_LDDRK,istage,USE_EXTERNAL_SOURCE_FILE,user_source_time_function,USE_BINARY_FOR_SEISMOGRAMS, &
                        UNDO_ATTENUATION_AND_OR_PML, &
                        NSOURCES,it,SIMULATION_TYPE,NSTEP, &
                        nrec,islice_selected_rec, &
                        nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                        Mesh_pointer, &
                        source_adjoint,nadj_rec_local,number_adjsources_global, &
                        INVERSE_FWI_FULL_PROBLEM,GPU_MODE

  use specfem_par_noise, only: irec_main_noise,noise_surface_movie

  ! coupling
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE

  ! faults
  use specfem_par, only: FAULT_SIMULATION

  implicit none

  ! local parameters
  double precision :: stf,time_source_dble,time_t
  double precision,external :: get_stf_viscoelastic
  logical ibool_read_adj_arrays
  ! for GPU_MODE
  double precision, dimension(NSOURCES) :: stf_pre_compute

  integer :: isource,it_sub_adj
  integer :: irec_local,irec

  character(len=MAX_STRING_LEN) :: adj_source_file

  ! no source inside the mesh if we are coupling with DSM
  ! because the source is precisely the wavefield coming from the DSM traction file
  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE == 1) return

  ! checks if anything to do
  if (.not. GPU_MODE) return

  ! forward simulations
  if (SIMULATION_TYPE == 1 .and. NOISE_TOMOGRAPHY == 0 .and. nsources_local > 0) then
    ! ignore CMT sources for fault rupture simulations
    if (FAULT_SIMULATION) return

    if (NSOURCES > 0) then
      ! sets current initial time
      if (USE_LDDRK) then
        ! LDDRK
        ! note: the LDDRK scheme updates displacement after the stiffness computations and
        !       after adding boundary/coupling/source terms.
        !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
        !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
        time_t = dble(it-1-1)*DT + dble(C_LDDRK(istage))*DT - t0
      else
        time_t = dble(it-1)*DT - t0
      endif

      do isource = 1,NSOURCES
        ! current time
        time_source_dble = time_t - tshift_src(isource)

        !! add external source time function
        if (USE_EXTERNAL_SOURCE_FILE) then
           stf = user_source_time_function(it, isource)
        else
           ! determines source time function value
           stf = get_stf_viscoelastic(time_source_dble,isource)
        endif

        ! stores precomputed source time function factor
        stf_pre_compute(isource) = stf
      enddo

      ! only implements SIMTYPE=1 and NOISE_TOM=0
      ! write(*,*) "Fortran dt = ", dt
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
      if (ibool_read_adj_arrays .and. .not. INVERSE_FWI_FULL_PROBLEM) then

        if (.not. SU_FORMAT) then
          ! ASCII format
          if (USE_BINARY_FOR_SEISMOGRAMS) stop 'Adjoint simulations not supported with .bin format, please use SU format instead'
          !!! read ascii adjoint sources
          do irec_local = 1, nadj_rec_local
            ! reads in **net**.**sta**.**BH**.adj files
            irec = number_adjsources_global(irec_local)
            adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
            call compute_arrays_adjoint_source(adj_source_file,irec_local)
          enddo
        else
          ! SU format
          call compute_arrays_adjoint_source_SU()
        endif !if (.not. SU_FORMAT)

      endif ! if (ibool_read_adj_arrays)

      if (it < NSTEP) then
        call add_sources_el_sim_type_2_or_3(Mesh_pointer, &
                                            source_adjoint, &
                                            nrec, &
                                            nadj_rec_local, &
                                            NTSTEP_BETWEEN_READ_ADJSRC, &
                                            it)
      endif ! it
    endif ! nadj_rec_local
  endif !adjoint

! note:  b_displ() is read in after Newmark time scheme, thus
!           b_displ(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!           thus indexing is NSTEP - it , instead of NSTEP - it - 1

! adjoint simulations
  if (SIMULATION_TYPE == 3 .and. NOISE_TOMOGRAPHY == 0 .and. nsources_local > 0) then
    ! ignore CMT sources for fault rupture simulations
    if (FAULT_SIMULATION) return

    if (NSOURCES > 0) then
      do isource = 1,NSOURCES
        ! current time
        if (USE_LDDRK) then
          ! LDDRK
          ! note: the LDDRK scheme updates displacement after the stiffness computations and
          !       after adding boundary/coupling/source terms.
          !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
          !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
          if (UNDO_ATTENUATION_AND_OR_PML) then
            ! stepping moves forward from snapshot position
            time_source_dble = dble(NSTEP-it-1)*DT + dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
          else
            time_source_dble = dble(NSTEP-it-1)*DT - dble(C_LDDRK(istage))*DT - t0 - tshift_src(isource)
          endif
        else
          time_source_dble = dble(NSTEP-it)*DT - t0 - tshift_src(isource)
        endif

        !! add external source time function
        if (USE_EXTERNAL_SOURCE_FILE) then
           stf = user_source_time_function(NSTEP-it+1, isource)
        else
           ! determines source time function value
           stf = get_stf_viscoelastic(time_source_dble,isource)
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
      call add_source_main_rec_noise_cu(Mesh_pointer,it,irec_main_noise,islice_selected_rec)
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

!=====================================================================

  subroutine compute_add_sources_viscoelastic_backward_GPU()

  use constants
  use specfem_par, only: nsources_local,tshift_src,dt,t0, &
                        USE_LDDRK,istage,USE_EXTERNAL_SOURCE_FILE,user_source_time_function, &
                        NSOURCES,it,SIMULATION_TYPE,NSTEP, &
                        NOISE_TOMOGRAPHY, &
                        Mesh_pointer,GPU_MODE

  ! coupling
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE

  ! undo_att
  use specfem_par, only: UNDO_ATTENUATION_AND_OR_PML,NSUBSET_ITERATIONS,NT_DUMP_ATTENUATION, &
                         iteration_on_subset,it_of_this_subset

  ! faults
  use specfem_par, only: FAULT_SIMULATION

  implicit none

  ! local parameters
  double precision :: stf,time_source_dble,time_t
  double precision,external :: get_stf_viscoelastic
  ! for GPU_MODE
  double precision, dimension(NSOURCES) :: stf_pre_compute

  integer :: isource,it_tmp

  ! no source inside the mesh if we are coupling with DSM
  ! because the source is precisely the wavefield coming from the DSM traction file
  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE == 1) return

  ! ignore CMT sources for fault rupture simulations
  if (FAULT_SIMULATION) return

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return
  if (.not. GPU_MODE) return

  ! iteration step
  if (UNDO_ATTENUATION_AND_OR_PML) then
    ! example: NSTEP is a multiple of NT_DUMP_ATTENUATION
    !         NT_DUMP_ATTENUATION = 301, NSTEP = 1204, NSUBSET_ITERATIONS = 4, iteration_on_subset = 1 -> 4,
    !              1. subset, it_temp goes from 301 down to 1
    !              2. subset, it_temp goes from 602 down to 302
    !              3. subset, it_temp goes from 903 down to 603
    !              4. subset, it_temp goes from 1204 down to 904
    !valid for multiples only:
    !it_tmp = iteration_on_subset * NT_DUMP_ATTENUATION - it_of_this_subset + 1
    !
    ! example: NSTEP is **NOT** a multiple of NT_DUMP_ATTENUATION
    !          NT_DUMP_ATTENUATION = 301, NSTEP = 900, NSUBSET_ITERATIONS = 3, iteration_on_subset = 1 -> 3
    !              1. subset, it_temp goes from (900 - 602) = 298 down to 1
    !              2. subset, it_temp goes from (900 - 301) = 599 down to 299
    !              3. subset, it_temp goes from (900 - 0)   = 900 down to 600
    !works always:
    it_tmp = NSTEP - (NSUBSET_ITERATIONS - iteration_on_subset)*NT_DUMP_ATTENUATION - it_of_this_subset + 1
  else
    it_tmp = it
  endif

  ! forward simulations
  if (NOISE_TOMOGRAPHY == 0 .and. nsources_local > 0) then
    ! sets current initial time
    if (USE_LDDRK) then
      ! LDDRK
      ! note: the LDDRK scheme updates displacement after the stiffness computations and
      !       after adding boundary/coupling/source terms.
      !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
      !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
      if (UNDO_ATTENUATION_AND_OR_PML) then
        ! stepping moves forward from snapshot position
        time_t = dble(NSTEP-it_tmp-1)*DT + dble(C_LDDRK(istage))*DT - t0
      else
        time_t = dble(NSTEP-it_tmp-1)*DT - dble(C_LDDRK(istage))*DT - t0
      endif
    else
      time_t = dble(NSTEP-it_tmp)*DT - t0
    endif

    do isource = 1,NSOURCES
      ! current time
      time_source_dble = time_t - tshift_src(isource)
      !! add external source time function
      if (USE_EXTERNAL_SOURCE_FILE) then
         stf = user_source_time_function(it, isource)
      else
         ! determines source time function value
         stf = get_stf_viscoelastic(time_source_dble,isource)
      endif
      ! stores precomputed source time function factor
      stf_pre_compute(isource) = stf
    enddo
    ! only implements SIMTYPE=3
    call compute_add_sources_el_s3_cuda(Mesh_pointer,stf_pre_compute,NSOURCES)
  endif

  ! for noise simulations
  if (NOISE_TOMOGRAPHY > 0) then
    stop 'for NOISE simulations, backward GPU routine is not implemented yet'
  endif

  end subroutine compute_add_sources_viscoelastic_backward_GPU


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

