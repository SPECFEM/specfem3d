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

! for acoustic solver

  subroutine compute_add_sources_acoustic()

  use constants
  use specfem_par, only: station_name,network_name, &
                         nsources_local,tshift_src,DT,t0,SU_FORMAT,USE_LDDRK,istage, &
                         hxir_adjstore,hetar_adjstore,hgammar_adjstore,source_adjoint,number_adjsources_global,nadj_rec_local, &
                         USE_BINARY_FOR_SEISMOGRAMS, &
                         ibool,NSOURCES,myrank,it,ispec_selected_source,islice_selected_source, &
                         sourcearrays,kappastore,SIMULATION_TYPE,NSTEP, &
                         ispec_selected_rec, &
                         nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC, INVERSE_FWI_FULL_PROBLEM

  use specfem_par_acoustic, only: potential_dot_dot_acoustic,ispec_is_acoustic

  ! faults
  use specfem_par, only: FAULT_SIMULATION

  implicit none

! local parameters
  real(kind=CUSTOM_REAL) :: stf_used,hlagrange
  logical :: ibool_read_adj_arrays
  double precision :: stf,time_source_dble,time_t
  double precision,external :: get_stf_acoustic

  integer :: isource,iglob,ispec,i,j,k
  integer :: irec_local,irec,it_sub_adj,it_index

  character(len=MAX_STRING_LEN) :: adj_source_file

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
  if (SIMULATION_TYPE == 1 .and. nsources_local > 0) then
    ! ignore pressure sources for fault rupture simulations
    if (FAULT_SIMULATION) return

! openmp solver
!$OMP PARALLEL if (NSOURCES > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(isource,time_source_dble,stf_used,stf,iglob,ispec,i,j,k)

    ! adds acoustic sources
!$OMP DO
    do isource = 1,NSOURCES

      !   add the source (only if this proc carries the source)
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        if (ispec_is_acoustic(ispec)) then
          ! current time
          time_source_dble = time_t - tshift_src(isource)

          ! determines source time function value
          stf = get_stf_acoustic(time_source_dble,isource,it)

          ! distinguishes between single and double precision for reals
          stf_used = real(stf,kind=CUSTOM_REAL)

          ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
          ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
          ! to add minus the source to Chi_dot_dot to get plus the source in pressure

          ! adds source array
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                ! adds source contribution
                ! note: acoustic source for pressure gets divided by kappa
                iglob = ibool(i,j,k,ispec)
!$OMP ATOMIC
                potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                  - sourcearrays(isource,1,i,j,k) * stf_used / kappastore(i,j,k,ispec)
              enddo
            enddo
          enddo

        endif ! ispec_is_acoustic
      endif ! myrank
    enddo ! NSOURCES
!$OMP ENDDO
!$OMP END PARALLEL

  endif

! NOTE: adjoint sources and backward wavefield timing:
!             idea is to start with the backward field b_potential,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_potential( it ) would correspond to (NSTEP - it - 1)*DT - t0
!       if we read in saved wavefields b_potential() before Newmark time scheme
!       (see sources for simulation_type 1 and seismograms)
!       since at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_potential( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0
!
!       b_potential is now read in after Newmark time scheme:
!       we read the backward/reconstructed wavefield at the end of the first time loop,
!       such that b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
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
      !     adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLY,NGLLZ)

      ! figure out if we need to read in a chunk of the adjoint source at this timestep
      it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )   !chunk_number
      ibool_read_adj_arrays = (((mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0)) .and. (nadj_rec_local > 0))

      ! needs to read in a new chunk/block of the adjoint source
      ! note that for each partition, we divide it into two parts --- boundaries and interior --- indicated by 'iphase'
      ! we first do calculations for the boudaries, and then start communication
      ! with other partitions while we calculate for the inner part
      ! this must be done carefully, otherwise the adjoint sources may be added twice
      if (ibool_read_adj_arrays .and. .not. INVERSE_FWI_FULL_PROBLEM) then

        if (.not. SU_FORMAT) then
          ! ASCII format
          if (USE_BINARY_FOR_SEISMOGRAMS) stop 'Adjoint simulations not supported with .bin format, please use SU format instead'

          !!! read ascii adjoint sources
          do irec_local = 1, nadj_rec_local
            irec = number_adjsources_global(irec_local)
            ! reads in **net**.**sta**.**BH**.adj files
            adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
            call compute_arrays_adjoint_source(adj_source_file,irec)
          enddo
        else
          ! SU format
          call compute_arrays_adjoint_source_SU()
        endif !if (.not. SU_FORMAT)

      endif ! if (ibool_read_adj_arrays)

      if (it < NSTEP) then
        ! receivers act as sources
        do irec_local = 1, nadj_rec_local
          irec = number_adjsources_global(irec_local)
          ! element index
          ispec = ispec_selected_rec(irec)
          if (ispec_is_acoustic(ispec)) then
            ! adds source
            do k = 1,NGLLZ
              do j = 1,NGLLY
                do i = 1,NGLLX
                  iglob = ibool(i,j,k,ispec)

                  hlagrange = hxir_adjstore(i,irec_local) * hetar_adjstore(j,irec_local) * hgammar_adjstore(k,irec_local)

                  ! time step index
                  it_index = NTSTEP_BETWEEN_READ_ADJSRC - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC)

                  ! beware, for acoustic medium, a pressure source would be taking the negative
                  ! and divide by Kappa of the fluid;
                  ! this would have to be done when constructing the adjoint source.
                  !
                  ! note: we take the first component of the adj_sourcearrays
                  !          the idea is to have e.g. a pressure source, where all 3 components would be the same
                  !
                  ! adjoint source of Peter et al. (A8):
                  !   f^adj = - sum_i \partial_t^2 (p^syn - p^obs)(T-t) \delta(x - x_i)
                  ! note that using the adjoint source derived from the optimization problem, there is no 1/kappa term applied
                  ! to the adjoint source. the negative sign also is part of the construction of the adjoint source.
                  !
                  ! since we don't know which formulation of adjoint source is used for the input,
                  ! we add the adjoint source as is, without 1/kappa factor, and with a positive sign.
                  potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                    + source_adjoint(1,irec_local,it_index) * hlagrange
                enddo
              enddo
            enddo
          endif
        enddo
      endif ! it
    endif ! nadj_rec_local > 0
  endif

  end subroutine compute_add_sources_acoustic

!
!=====================================================================

! for acoustic solver for back propagation wave field

  subroutine compute_add_sources_acoustic_backward()

  use constants
  use specfem_par, only: nsources_local,tshift_src,DT,t0,USE_LDDRK,istage, &
                         ibool,NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
                         sourcearrays,kappastore,SIMULATION_TYPE,NSTEP

  ! undo_att
  use specfem_par, only: UNDO_ATTENUATION_AND_OR_PML,NSUBSET_ITERATIONS,NT_DUMP_ATTENUATION, &
                         iteration_on_subset,it_of_this_subset

  use specfem_par_acoustic, only: ispec_is_acoustic,b_potential_dot_dot_acoustic

  ! faults
  use specfem_par, only: FAULT_SIMULATION

  implicit none

! local parameters
  real(kind=CUSTOM_REAL) :: stf_used

  double precision :: stf,time_source_dble,time_t
  double precision, external :: get_stf_acoustic

  integer :: isource,iglob,ispec,i,j,k,it_tmp

  ! ignore pressure sources for fault rupture simulations
  if (FAULT_SIMULATION) return

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return

  ! checks if this slice has sources to add
  if (nsources_local == 0) return

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
!             idea is to start with the backward field b_potential,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_potential( it ) would correspond to (NSTEP - it - 1)*DT - t0
!       if we read in saved wavefields b_potential() before Newmark time scheme
!       (see sources for simulation_type 1 and seismograms)
!       since at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_potential( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0
!
!       b_potential is now read in after Newmark time scheme:
!       we read the backward/reconstructed wavefield at the end of the first time loop,
!       such that b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!       assuming that until that end the backward/reconstructed wavefield and adjoint fields
!       have a zero contribution to adjoint kernels.
!       thus the correct indexing is NSTEP - it + 1, instead of NSTEP - it
!
! adjoint wavefields:
!       since the adjoint source traces were derived from the seismograms,
!       it follows that for the adjoint wavefield, the time equivalent to ( T - t ) uses the time-reversed
!       adjoint source traces which start at -t0 and end at time (NSTEP-1)*DT - t0
!       for step it=1: (NSTEP -it + 1)*DT - t0 for backward wavefields corresponds to time T

! note:  b_potential() is read in after Newmark time scheme, thus
!        b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!        thus indexing is NSTEP - it , instead of NSTEP - it - 1
!
!        this leads to the timing (NSTEP-(it-1)-1)*DT-t0-tshift_src for the source time function here

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

  ! adds sources in backward fields
  ! adds acoustic sources
  do isource = 1,NSOURCES

    !   add the source (only if this proc carries the source)
    if (myrank == islice_selected_source(isource)) then

      ispec = ispec_selected_source(isource)

      if (ispec_is_acoustic(ispec)) then
        ! current time
        time_source_dble = time_t - tshift_src(isource)

        ! determines source time function value
        stf = get_stf_acoustic(time_source_dble,isource,NSTEP-it_tmp+1)

        ! distinguishes between single and double precision for reals
        stf_used = real(stf,kind=CUSTOM_REAL)

        ! add source array
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              ! adds source contribution
              ! note: acoustic source for pressure gets divided by kappa
              iglob = ibool(i,j,k,ispec)
              b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                                                  - sourcearrays(isource,1,i,j,k) * stf_used / kappastore(i,j,k,ispec)
            enddo
          enddo
        enddo

      endif ! ispec_is_acoustic
    endif ! myrank
  enddo ! NSOURCES

  end subroutine compute_add_sources_acoustic_backward

!
!=====================================================================

! for acoustic solver on GPU

  subroutine compute_add_sources_acoustic_GPU()

  use constants
  use specfem_par, only: station_name,network_name, &
                         nsources_local,tshift_src,DT,t0,SU_FORMAT,USE_LDDRK,istage, &
                         source_adjoint,nadj_rec_local,number_adjsources_global, &
                         USE_BINARY_FOR_SEISMOGRAMS, &
                         NSOURCES,it,SIMULATION_TYPE,NSTEP,nrec, &
                         NTSTEP_BETWEEN_READ_ADJSRC,Mesh_pointer, &
                         INVERSE_FWI_FULL_PROBLEM,run_number_of_the_source, &
                         GPU_MODE

  ! faults
  use specfem_par, only: FAULT_SIMULATION

  implicit none

  ! local parameters
  logical :: ibool_read_adj_arrays
  integer :: it_sub_adj

  double precision :: stf,time_source_dble,time_t
  double precision, external :: get_stf_acoustic

  double precision, dimension(NSOURCES) :: stf_pre_compute

  integer :: isource
  integer :: irec_local,irec

  character(len=MAX_STRING_LEN) :: adj_source_file

  ! checks if anything to do
  if (.not. GPU_MODE) return

  ! forward simulations
  if (SIMULATION_TYPE == 1 .and. nsources_local > 0) then
    ! ignore pressure sources for fault rupture simulations
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

        ! determines source time function value
        stf = get_stf_acoustic(time_source_dble,isource,it)

        ! stores precomputed source time function factor
        stf_pre_compute(isource) = stf
      enddo

      ! only implements SIMTYPE=1 and NOISE_TOM=0
      ! write(*,*) "Fortran dt = ", dt
      ! change dt -> DT
      call compute_add_sources_ac_cuda(Mesh_pointer,NSOURCES,stf_pre_compute,run_number_of_the_source)
    endif
  endif

! NOTE: adjoint sources and backward wavefield timing:
!             idea is to start with the backward field b_potential,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_potential( it ) would correspond to (NSTEP - it - 1)*DT - t0
!       if we read in saved wavefields b_potential() before Newmark time scheme
!       (see sources for simulation_type 1 and seismograms)
!       since at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_potential( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0
!
!       b_potential is now read in after Newmark time scheme:
!       we read the backward/reconstructed wavefield at the end of the first time loop,
!       such that b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
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
      !     adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLY,NGLLZ)

      ! figure out if we need to read in a chunk of the adjoint source at this timestep
      it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )   !chunk_number
      ibool_read_adj_arrays = (((mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0)) .and. (nadj_rec_local > 0))

      ! needs to read in a new chunk/block of the adjoint source
      ! note that for each partition, we divide it into two parts --- boundaries and interior --- indicated by 'iphase'
      ! we first do calculations for the boudaries, and then start communication
      ! with other partitions while we calculate for the inner part
      ! this must be done carefully, otherwise the adjoint sources may be added twice
      if (ibool_read_adj_arrays .and. .not. INVERSE_FWI_FULL_PROBLEM) then

        if (.not. SU_FORMAT) then
          ! ASCII format
          if (USE_BINARY_FOR_SEISMOGRAMS) stop 'Adjoint simulations not supported with .bin format, please use SU format instead'
          !!! read ascii adjoint sources
          do irec_local = 1, nadj_rec_local
            irec = number_adjsources_global(irec_local)
            ! reads in **net**.**sta**.**BH**.adj files
            adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
            call compute_arrays_adjoint_source(adj_source_file,irec)
          enddo
        else
          ! SU format
          call compute_arrays_adjoint_source_SU()
        endif !if (.not. SU_FORMAT)

      endif ! if (ibool_read_adj_arrays)

      if (it < NSTEP) then
        ! receivers act as sources
        ! on GPU
        call add_sources_ac_sim_2_or_3_cuda(Mesh_pointer,source_adjoint,nrec,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,it)
      endif ! it
    endif ! nadj_rec_local > 0
  endif

  end subroutine compute_add_sources_acoustic_GPU

!
!=====================================================================
!

  subroutine compute_add_sources_acoustic_backward_GPU()

  use constants
  use specfem_par, only: nsources_local,tshift_src,DT,t0,USE_LDDRK,istage, &
                         NSOURCES,myrank,it, &
                         SIMULATION_TYPE,NSTEP, &
                         GPU_MODE,Mesh_pointer,run_number_of_the_source
  ! undo_att
  use specfem_par, only: UNDO_ATTENUATION_AND_OR_PML,NSUBSET_ITERATIONS,NT_DUMP_ATTENUATION, &
                         iteration_on_subset,it_of_this_subset

  ! faults
  use specfem_par, only: FAULT_SIMULATION

  implicit none

  ! local parameters
  double precision :: stf,time_source_dble,time_t
  double precision, external :: get_stf_acoustic

  double precision, dimension(NSOURCES) :: stf_pre_compute

  integer :: isource,it_tmp

  ! ignore pressure sources for fault rupture simulations
  if (FAULT_SIMULATION) return

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return
  if (.not. GPU_MODE) return

  ! checks if this slice has sources to add
  if (nsources_local == 0) return

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
!             idea is to start with the backward field b_potential,.. at time (T)
!             and convolve with the adjoint field at time (T-t)
!
! backward/reconstructed wavefields:
!       time for b_potential( it ) would correspond to (NSTEP - it - 1)*DT - t0
!       if we read in saved wavefields b_potential() before Newmark time scheme
!       (see sources for simulation_type 1 and seismograms)
!       since at the beginning of the time loop, the numerical Newmark time scheme updates
!       the wavefields, that is b_potential( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0
!
!       b_potential is now read in after Newmark time scheme:
!       we read the backward/reconstructed wavefield at the end of the first time loop,
!       such that b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!       assuming that until that end the backward/reconstructed wavefield and adjoint fields
!       have a zero contribution to adjoint kernels.
!       thus the correct indexing is NSTEP - it + 1, instead of NSTEP - it
!
! adjoint wavefields:
!       since the adjoint source traces were derived from the seismograms,
!       it follows that for the adjoint wavefield, the time equivalent to ( T - t ) uses the time-reversed
!       adjoint source traces which start at -t0 and end at time (NSTEP-1)*DT - t0
!       for step it=1: (NSTEP -it + 1)*DT - t0 for backward wavefields corresponds to time T

! note:  b_potential() is read in after Newmark time scheme, thus
!        b_potential(it=1) corresponds to -t0 + (NSTEP-1)*DT.
!        thus indexing is NSTEP - it , instead of NSTEP - it - 1
!
!        this leads to the timing (NSTEP-(it-1)-1)*DT-t0-tshift_src for the source time function here

  ! sets current initial time
  if (USE_LDDRK) then
    ! LDDRK
    ! note: the LDDRK scheme updates displacement after the stiffness computations and
    !       after adding boundary/coupling/source terms.
    !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
    !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
    if (UNDO_ATTENUATION_AND_OR_PML) then
      ! steps forward in time
      time_t = dble(NSTEP-it_tmp-1)*DT + dble(C_LDDRK(istage))*DT - t0
    else
      ! steps backward
      time_t = dble(NSTEP-it_tmp-1)*DT - dble(C_LDDRK(istage))*DT - t0
    endif
  else
    time_t = dble(NSTEP-it_tmp)*DT - t0
  endif

  ! adds sources in backward fields
  do isource = 1,NSOURCES
    ! current time
    time_source_dble = time_t - tshift_src(isource)

    ! determines source time function value
    stf = get_stf_acoustic(time_source_dble,isource,NSTEP-it_tmp+1)

    ! stores precomputed source time function factor
    stf_pre_compute(isource) = stf
  enddo

  ! only implements SIMTYPE=3
  call compute_add_sources_ac_s3_cuda(Mesh_pointer,NSOURCES,stf_pre_compute,run_number_of_the_source)

  end subroutine compute_add_sources_acoustic_backward_GPU

!
!=====================================================================
!

  double precision function get_stf_acoustic(time_source_dble,isource,it_tmp_ext)

! returns source time function value for specified time

  use specfem_par, only: USE_FORCE_POINT_SOURCE,USE_RICKER_TIME_FUNCTION,USE_TRICK_FOR_BETTER_PRESSURE, &
                         USE_SOURCE_ENCODING,pm1_source_encoding, &
                         hdur,hdur_Gaussian,force_stf

  ! for external STFs
  use specfem_par, only: USE_EXTERNAL_SOURCE_FILE

  implicit none

  double precision,intent(in) :: time_source_dble
  integer,intent(in) :: isource
  integer,intent(in) :: it_tmp_ext

  ! local parameters
  double precision :: stf,f0

  double precision, external :: comp_source_time_function,comp_source_time_function_rickr, &
    comp_source_time_function_d2rck,comp_source_time_function_gauss,comp_source_time_function_d2gau, &
    comp_source_time_function_gauss_2,comp_source_time_function_mono,comp_source_time_function_ext

  ! external source time function
  if (USE_EXTERNAL_SOURCE_FILE) then
    ! gets stf value
    stf = comp_source_time_function_ext(it_tmp_ext, isource)

    ! source encoding
    if (USE_SOURCE_ENCODING) stf = stf * pm1_source_encoding(isource)

    ! returns value
    get_stf_acoustic = stf
    return
  endif

  ! determines source time function value
  if (USE_FORCE_POINT_SOURCE) then
    ! single point force
    select case(force_stf(isource))
    case (0)
      ! Gaussian source time function value
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      if (USE_TRICK_FOR_BETTER_PRESSURE) then
        stf = comp_source_time_function_d2gau(time_source_dble,hdur_Gaussian(isource))
      else
        stf = comp_source_time_function_gauss(time_source_dble,hdur_Gaussian(isource))
      endif
    case (1)
      ! Ricker source time function
      f0 = hdur(isource) ! using hdur as a FREQUENCY just to avoid changing FORCESOLUTION file format
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      if (USE_TRICK_FOR_BETTER_PRESSURE) then
        stf = comp_source_time_function_d2rck(time_source_dble,f0)
      else
        stf = comp_source_time_function_rickr(time_source_dble,f0)
      endif
    case (2)
      ! Heaviside (step) source time function
      stf = comp_source_time_function(time_source_dble,hdur_Gaussian(isource))
    case (3)
      ! Monochromatic source time function
      f0 = 1.d0 / hdur(isource) ! using hdur as a PERIOD just to avoid changing FORCESOLUTION file format
      stf = comp_source_time_function_mono(time_source_dble,f0)
    case (4)
      ! Gaussian source time function by Meschede et al. (2011)
      stf = comp_source_time_function_gauss_2(time_source_dble,hdur(isource))
    case default
      stop 'unsupported force_stf value!'
    end select

  else
    ! moment-tensor source in acoustic domain
    if (USE_RICKER_TIME_FUNCTION) then
      ! Ricker
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      if (USE_TRICK_FOR_BETTER_PRESSURE) then
        stf = comp_source_time_function_d2rck(time_source_dble,hdur(isource))
      else
        stf = comp_source_time_function_rickr(time_source_dble,hdur(isource))
      endif
    else
      ! Gaussian source time
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      if (USE_TRICK_FOR_BETTER_PRESSURE) then
        stf = comp_source_time_function_d2gau(time_source_dble,hdur_Gaussian(isource))
      else
        stf = comp_source_time_function_gauss(time_source_dble,hdur_Gaussian(isource))
      endif

      ! quasi-Heaviside for acoustic sources:
      ! careful if using this, as imposing permanent source offset in fluids is non-physical
      ! stf = comp_source_time_function(time_source_dble,hdur_Gaussian(isource))
    endif

  endif ! USE_FORCE_POINT_SOURCE

  ! source encoding
  if (USE_SOURCE_ENCODING) stf = stf * pm1_source_encoding(isource)

  ! return value
  get_stf_acoustic = stf

  end function get_stf_acoustic
