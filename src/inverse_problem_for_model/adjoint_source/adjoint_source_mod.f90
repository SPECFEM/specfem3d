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

module adjoint_source

  !! IMPORT VARIABLES FROM SPECFEM -------------------------------------------------------------------------------------------------
  use specfem_par, only: CUSTOM_REAL, NDIM, MAX_STRING_LEN, OUTPUT_FILES, IIN, network_name, station_name, nrec_local, &
       seismograms_d, seismograms_p, number_receiver_global, &
       ELASTIC_SIMULATION,ACOUSTIC_SIMULATION

  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic

  !!--------------------------------------------------------------------------------------------------------------------------------
  !! IMPORT inverse_problem VARIABLES
  use inverse_problem_par
  use signal_processing
  use interpolation_mod, only: mycorrelation, myconvolution
  use rotations_mod
  implicit none

  !! PRIVATE ATTRIBUTE -------------------------------------------------------------------------------------------------------------
  !real(kind=CUSTOM_REAL), private, dimension(:,:), allocatable      :: elastic_adjoint_source (not used yet)
  real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: raw_residuals, fil_residuals, filfil_residuals
  real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: data_trace_to_use, wkstmp, w_tap
  real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: residuals_for_cost
  !real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: signal (not used yet)

  real(kind=CUSTOM_REAL), private                                   :: fl, fh
  real(kind=CUSTOM_REAL), private                                   :: dt_data
  real(kind=CUSTOM_REAL), private                                   :: prior_data_std
  real(kind=CUSTOM_REAL), private                                   :: nb_traces_tot, window_length
  real(kind=CUSTOM_REAL), private                                   :: cost_function, data_std, nb_data_std

  integer,                private                                   :: nstep_data
  integer,                private                                   :: current_ifrq

  ! bandpass filter
  logical,                private                                   :: use_band_pass_filter
  ! bandpass filter order
  integer,                private, parameter                        :: norder_filter = 4, irek_filter = 1

  ! tapering trace ends should be default
  logical,                private, parameter                        :: taper_ends = .true.
  real(kind=CUSTOM_REAL), private, parameter                        :: taper_width = 0.02  ! 2% percent of trace length
  integer,                private, parameter                        :: taper_onset = 10    ! 10 samples

  ! data weights
  logical,                private                                   :: set_event_trace_weights

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!  Write adjoint sources in SEM directory
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine write_adjoint_sources_for_specfem(acqui_simu,  inversion_param, ievent, myrank)

    implicit none

    integer,                                     intent(in)     :: myrank, ievent
    type(inver),                                 intent(inout)  :: inversion_param
    type(acqui),  dimension(:),                  intent(inout)  :: acqui_simu
    ! local
    real(kind=CUSTOM_REAL)                                      :: tmp_val
    integer                                                     :: irec, irec_local, ispec
    integer                                                     :: lw, i0, i1, i2, i3

    ! initializes
    acqui_simu(ievent)%adjoint_sources(:,:,:) = 0._CUSTOM_REAL

    ! sets current event data infos
    nstep_data = acqui_simu(ievent)%Nt_data
    dt_data    = acqui_simu(ievent)%dt_data

    use_band_pass_filter = inversion_param%use_band_pass_filter
    prior_data_std       = inversion_param%prior_data_std
    nb_traces_tot        = inversion_param%nb_traces_tot

    current_ifrq         = iter_frq

    ! only sets data weights once at beginning of fwi
    if (iter_frq == 1 .and. iter_inverse == 0 .and. iter_wolfe == 0) then
      set_event_trace_weights = .true.
    else
      set_event_trace_weights = .false.
    endif

    !! initialize cost function for each MPI porcess
    cost_function = 0._CUSTOM_REAL
    data_std      = 0._CUSTOM_REAL
    nb_data_std   = 0._CUSTOM_REAL

    ! log output
    if (myrank == 0 .and. (VERBOSE_MODE .or. DEBUG_MODE) ) then
      write(INVERSE_LOG_FILE,*) ' adjoint sources : ',trim(acqui_simu(ievent)%adjoint_source_type)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! allocate temporary arrays
    call allocate_adjoint_source_working_arrays()

    !! define taper on adjoint sources (if not window selected by user)
    !! to do: define user window
    if (taper_ends) then
      ! taper
      lw = int(taper_width * nstep_data)  ! width
      i0 = taper_onset                    ! onset
      i1 = i0 + lw
      i3 = nstep_data - taper_onset
      i2 = i3 - lw
      call taper_window_W(w_tap,i0,i1,i2,i3,nstep_data,1._CUSTOM_REAL)
    else
      ! no tapering
      w_tap(:) = 1._CUSTOM_REAL
    endif

    do irec_local = 1, nrec_local
      ! check: assumes SIMULATION_TYPE == 3, with adjoint sources and receivers being identical
      !        -  nrec_local ==  nadj_rec_local
      !        -  number_receiver_global == number_adjsources_global
      irec = number_receiver_global(irec_local)
      ispec = acqui_simu(ievent)%ispec_selected_rec(irec)

      !! ALLOW TO CHOOSE COMPONENT :
      !! UX UY UZ (check if in solid region)
      !! Pr (check if in fluid region)
      !! VX VY VZ

      ! receiver in elastic domains
      if (ispec_is_elastic(ispec)) then
        !! compute adjoint source according to cost L2 function
        call compute_elastic_adjoint_source_displacement(irec_local, ievent, acqui_simu, inversion_param)
      endif

      ! receiver in acoustic domains
      if (ispec_is_acoustic(ispec)) then
        !! compute adjoint source according to cost L2 function
        call compute_acoustic_adjoint_source_pressure_dot_dot(irec_local, ievent, acqui_simu)
      endif
    enddo

    !! compute cost function   : allreduce cost_function
    tmp_val = cost_function
    call sum_all_all_cr(tmp_val,cost_function)

    !! add the cost function over all sources in group
    inversion_param%total_current_cost = inversion_param%total_current_cost + cost_function

    !! save cost function for the current source (not used yet)
    !inversion_param%current_cost(ievent) = cost_function

    !! standard deviation on data
    tmp_val = data_std
    call sum_all_all_cr(tmp_val,data_std)
    inversion_param%data_std = inversion_param%data_std + data_std

    tmp_val = nb_data_std
    call sum_all_all_cr(tmp_val,nb_data_std)
    inversion_param%nb_data_std = inversion_param%nb_data_std + nb_data_std

    ! log output
    if (myrank == 0 .and. (VERBOSE_MODE .or. DEBUG_MODE) ) then
      write(INVERSE_LOG_FILE,*) '   Cost function for this event : ', cost_function
      write(INVERSE_LOG_FILE,*) '   weight on data    : ', 1._CUSTOM_REAL/prior_data_std
      write(INVERSE_LOG_FILE,*) '   number of traces  : ', nb_traces_tot
      write(INVERSE_LOG_FILE,*) '   total weight      : ', &
                                1._CUSTOM_REAL / sqrt(nb_traces_tot) / prior_data_std / sqrt(nstep_data*dt_data)
      !write(INVERSE_LOG_FILE,*) '   total standard data dev : ', inversion_param%data_std
      !write(INVERSE_LOG_FILE,*) '   total standard data nb  : ', inversion_param%nb_data_std
      !write(INVERSE_LOG_FILE,*) '   total current costs          : ', inversion_param%total_current_cost
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! free temporary arrays
    call deallocate_adjoint_source_working_arrays()

  end subroutine write_adjoint_sources_for_specfem

!----------------------------------------------------------------------------------------------------------------------------------
  subroutine deallocate_adjoint_source_working_arrays()

    deallocate(raw_residuals, fil_residuals,  filfil_residuals, w_tap, residuals_for_cost, &
               data_trace_to_use, wkstmp)
    !deallocate(elastic_adjoint_source, signal) (not used yet)

  end subroutine deallocate_adjoint_source_working_arrays

!----------------------------------------------------------------------------------------------------------------------------------
  subroutine  allocate_adjoint_source_working_arrays()

    integer :: ier

    allocate(raw_residuals(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 280')
    raw_residuals(:) = 0._CUSTOM_REAL

    allocate(fil_residuals(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 281')
    fil_residuals(:) = 0._CUSTOM_REAL

    allocate(filfil_residuals(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 282')
    filfil_residuals(:) = 0._CUSTOM_REAL

    allocate(w_tap(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 283')
    w_tap(:) = 0._CUSTOM_REAL

    allocate(residuals_for_cost(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 285')
    residuals_for_cost(:) = 0._CUSTOM_REAL

    allocate(data_trace_to_use(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 288')
    data_trace_to_use(:) = 0._CUSTOM_REAL

    allocate(wkstmp(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 289')
    wkstmp(:) = 0._CUSTOM_REAL

    !(not used yet)
    !allocate(elastic_adjoint_source(NDIM,nstep_data),stat=ier)
    !if (ier /= 0) call exit_MPI_without_rank('error allocating array 286')
    !elastic_adjoint_source(:,:) = 0._CUSTOM_REAL

    !allocate(signal(nstep_data),stat=ier)
    !if (ier /= 0) call exit_MPI_without_rank('error allocating array 284')
    !signal(:) = 0._CUSTOM_REAL

  end subroutine allocate_adjoint_source_working_arrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
! define adjoint sources
!-----------------------------------------------------------------------------------------------------------------------------------
!
! to implement a new cost function and adjoint source :
!
!    1 /   define character to name your cost function and adjoint source type :  acqui_simu(ievent)%adjoint_source_type
!
!    2/    add a case for your new adjoint source
!
!    3/    compute the adjoint source which is stored in acqui(ievent)%adjoint_sources(NCOM, NREC_LOCAL, NT)
!          note that this arrays will be directly used by specfem as adjoint source thus you need to
!          do any processing here : filter, rotation, ....
!
!    4/    compute the cost function and store it in cost_function variable (you have to perform summation over sources)
!
!
!
  subroutine compute_elastic_adjoint_source_displacement(irec_local, ievent, acqui_simu, inversion_param)

    integer,                                     intent(in)    :: ievent, irec_local
    type(acqui),  dimension(:),                  intent(inout) :: acqui_simu
    type(inver),                                 intent(inout) :: inversion_param
    ! local
    integer                                                    :: icomp, irec_glob
    real(kind=CUSTOM_REAL), dimension(:), allocatable          :: wavelet, filfil_residuals_tmp, tmpl
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable        :: trace_cal_1, trace_cal_2
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable        :: trace_obs_1, trace_obs_2
    real(kind=CUSTOM_REAL)                                     :: cost_value
    double precision                                           :: lat0, lon0, azi0
    integer                                                    :: ier

    !! store residuals and filter

    !! TO DO : convolve seismograms_d(icomp,irec_local,:) by source time function stf (if need)
    !! if (acqui_simu(ievent)%convlove_residuals_by_wavelet) then
    !!     signal(:) = seismograms_d(icomp,irec_local,:)
    !!     call convolution_by_wavelet(wavelet, signal, seismograms_d(icomp,irec_local,:), nstep, nw)
    !! endif

    !! TO DO : for now only L2 fwi is used, we should consider other adjoint sources
    !!
    !!   -- travel time kernel
    !!   -- Enveloppe
    !!   -- rotation to change components
    !!

    select case (trim(adjustl(acqui_simu(ievent)%adjoint_source_type)))

    case ('L2_FWI_TELESEISMIC')

      ! Define temporary trace vector
      if (.not. allocated(trace_cal_1)) then
        allocate(trace_cal_1(3,nstep_data),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 290')
      endif
      if (.not. allocated(trace_cal_2)) then
        allocate(trace_cal_2(3,nstep_data),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 291')
      endif
      if (.not. allocated(trace_obs_1)) then
        allocate(trace_obs_1(3,nstep_data),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 292')
      endif
      if (.not. allocated(trace_obs_2)) then
        allocate(trace_obs_2(3,nstep_data),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 293')
      endif

      lat0 = acqui_simu(ievent)%Origin_chunk_lat
      lon0 = acqui_simu(ievent)%Origin_chunk_lon
      azi0 = acqui_simu(ievent)%Origin_chunk_azi

      ! global receiver index
      irec_glob = acqui_simu(ievent)%number_receiver_global(irec_local)

      ! set traces for rotations
      do icomp = 1, NDIM
        ! synthetics
        trace_cal_1(icomp,:) = seismograms_d(icomp,irec_local,:)
        trace_cal_2(icomp,:) = seismograms_d(icomp,irec_local,:)
        ! observation data
        trace_obs_1(icomp,:) = acqui_simu(ievent)%data_traces(irec_local,:,icomp)
        trace_obs_2(icomp,:) = acqui_simu(ievent)%data_traces(irec_local,:,icomp)

        ! Convolve synthetic data with wavelet
        if (inversion_param%convolution_by_wavelet) then
          if (.not. allocated(wavelet)) then
            allocate(wavelet(nstep_data),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 294')
          endif
          wavelet(:) = acqui_simu(ievent)%user_source_time_function(1,:)
          call myconvolution(trace_cal_2(icomp,:),wavelet,nstep_data,nstep_data,tmpl,0)
          trace_cal_1(icomp,:) = tmpl(:) * dt_data
        endif
      enddo

      ! Do rotation of data (would be easier if lat and lon are split over MPI slices)
      select case (inversion_param%inverted_data_sys)
      case ('xyz')
        ! Do nothing
      case('enz')
        trace_cal_2 = trace_cal_1
        trace_obs_2 = trace_obs_1
        ! Data are in standard coordinate system
        ! Data rotation required to pass in mesh system (zen -> xyz)
        call define_mesh_rotation_matrix(lat0,lon0,azi0)
        ! synthetics
        call rotate_comp_mesh2glob(trace_cal_2(1,:), trace_cal_2(2,:), trace_cal_2(3,:), &
                                   acqui_simu(ievent)%read_station_position(1,irec_glob), &
                                   acqui_simu(ievent)%read_station_position(2,irec_glob), &
                                   nstep_data, 1, trace_cal_1(3,:), trace_cal_1(2,:), trace_cal_1(1,:))
        ! observation data
        call rotate_comp_mesh2glob(trace_obs_2(1,:), trace_obs_2(2,:), trace_obs_2(3,:), &
                                   acqui_simu(ievent)%read_station_position(1,irec_glob), &
                                   acqui_simu(ievent)%read_station_position(2,irec_glob), &
                                   nstep_data, 1, trace_obs_1(3,:), trace_obs_1(2,:), trace_obs_1(1,:))
      case('rtz')
        trace_cal_2 = trace_cal_1
        trace_obs_2 = trace_obs_1
        ! Data are in the souce receiver coordinate system
        ! Data rotation required (baz-azi) (rtz -> zne)
        ! Data rotation required to pass in mesh system (zen -> xyz)
        call define_mesh_rotation_matrix(lat0,lon0,azi0)
        ! synthetics
        call rotate_comp_mesh2glob(trace_cal_2(1,:), trace_cal_2(2,:), trace_cal_2(3,:), &
                                   acqui_simu(ievent)%read_station_position(1,irec_glob), &
                                   acqui_simu(ievent)%read_station_position(2,irec_glob), &
                                   nstep_data, 1, trace_cal_1(3,:), trace_cal_1(2,:), trace_cal_1(1,:))
        ! observation data
        call rotate_comp_mesh2glob(trace_obs_2(1,:), trace_obs_2(2,:), trace_obs_2(3,:), &
                                   acqui_simu(ievent)%read_station_position(1,irec_glob), &
                                   acqui_simu(ievent)%read_station_position(2,irec_glob), &
                                   nstep_data, 1, trace_obs_1(3,:), trace_obs_1(2,:), trace_obs_1(1,:))

        trace_cal_2 = trace_cal_1
        trace_obs_2 = trace_obs_1
        call rotate_ZNE_to_ZRT(trace_cal_2(3,:), trace_cal_2(2,:), trace_cal_2(1,:), &
                               trace_cal_1(3,:), trace_cal_1(1,:), trace_cal_1(2,:), &
                               1,nstep_data,acqui_simu(ievent)%baz(irec_glob))
        call rotate_ZNE_to_ZRT(trace_obs_2(3,:), trace_obs_2(2,:), trace_obs_2(1,:), &
                               trace_obs_1(3,:), trace_obs_1(1,:), trace_obs_1(2,:), &
                               1,nstep_data,acqui_simu(ievent)%baz(irec_glob))
      case('qtl')
        trace_cal_2 = trace_cal_1
        trace_obs_2 = trace_obs_1
        write(*,*)'ERROR : qtl is not implemented yet'
        write(*,*)'NOW STOP'
        stop
      end select

      ! Finally compute residuals, filter and cross-correlate
      do icomp = 1, NDIM
        ! get (rotated) data trace
        wkstmp(:) = trace_obs_1(icomp,:)

        ! Residual (syn - obs)
        raw_residuals(:) =  trace_cal_1(icomp,:) - trace_obs_1(icomp,:)

        ! save residuals because the adjoint source is not residuals
        residuals_for_cost(:) = raw_residuals(:)

        ! tapers (avoids non-zero values at start/end for better filtering results)
        if (taper_ends) then
          residuals_for_cost(:) = residuals_for_cost(:) * w_tap(:)
          raw_residuals(:) = raw_residuals(:) * w_tap(:)
          wkstmp(:) = wkstmp(:) * w_tap(:)
        endif

        ! Filter
        if (use_band_pass_filter) then
          ! band pass filter
          fil_residuals(:) = 0._CUSTOM_REAL
          filfil_residuals(:) = 0._CUSTOM_REAL

          ! teleseismic filter limits
          fl = acqui_simu(ievent)%freqcy_to_invert(icomp,1,irec_local)
          fh = acqui_simu(ievent)%freqcy_to_invert(icomp,2,irec_local)

          ! filters for misfit calculation
          call bwfilt (raw_residuals, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

          ! filters twice for adjoint source
          ! note: filtered adjoint kernels would require to filter both the forward/reconstructed and adjoint wavefields.
          !       to avoid filtering the forward/reconstructed wavefield, it is equivalent to filtering twice only
          !       the adjoint source that produces the adjoint wavefield.
          call bwfilt (fil_residuals, filfil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

          ! filters for misfit calculation
          call bwfilt (residuals_for_cost, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

          ! filter data for storage
          call bwfilt (wkstmp, data_trace_to_use, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
        else
          ! no filter
          fil_residuals(:) = residuals_for_cost(:)
          filfil_residuals(:) = raw_residuals(:)
          data_trace_to_use(:) = wkstmp(:)
        endif

        !! save filtered data
        if (VERBOSE_MODE .or. DEBUG_MODE) then
          acqui_simu(ievent)%synt_traces(icomp, irec_local,:) = data_trace_to_use(:)
        endif

        ! Apply weighting (for teleseismic, weights defined on all stations -> irec_glob)
        fil_residuals(:) = fil_residuals(:) * acqui_simu(ievent)%weight_trace(icomp,irec_glob,:)
        filfil_residuals(:) = filfil_residuals(:) * acqui_simu(ievent)%weight_trace(icomp,irec_glob,:)**2

        !! compute cost function value
        cost_value = sum(fil_residuals(:)**2) * 0.5 * dt_data
        cost_function = cost_function + cost_value

        !! compute raw standard deviation
        data_std = data_std + sum( residuals_for_cost(:)**2 )
        nb_data_std = nb_data_std + size(residuals_for_cost(:))

        ! Finally cross-correlate residuals with wavelet
        if (inversion_param%convolution_by_wavelet) then
          if (.not. allocated(filfil_residuals_tmp)) then
            allocate(filfil_residuals_tmp(nstep_data),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 295')
          endif
          filfil_residuals_tmp(:) = filfil_residuals(:)
          call mycorrelation(filfil_residuals_tmp,wavelet,nstep_data,nstep_data,tmpl,0)
          filfil_residuals(:) = tmpl(:) * dt_data
        endif

        ! tapers adjoint source (again to avoid non-zero onsets)
        if (taper_ends) then
          filfil_residuals(:) = filfil_residuals(:) * w_tap(:)
        endif

        !! store the adjoint source
        acqui_simu(ievent)%adjoint_sources(icomp,irec_local,:) = filfil_residuals(:)

        !(not used yet)
        !elastic_adjoint_source(icomp,:) = filfil_residuals(:)
      enddo

    case ('L2_OIL_INDUSTRY')

      do icomp = 1, NDIM
        !! get data
        wkstmp(:) = acqui_simu(ievent)%data_traces(irec_local,:,icomp)

        ! note: for data comparison, we would want to filter both data and synthetics with the same filter.
        !       given the filter is linear, this is identical to filter the residual: F(syn) - F(obs) = F(syn - obs)

        ! Residual (syn - obs)
        raw_residuals(:) = (seismograms_d(icomp,irec_local,:) - acqui_simu(ievent)%data_traces(irec_local,:,icomp))

        ! save residuals because the adjoint source is not residuals
        residuals_for_cost(:) = raw_residuals(:)

        ! tapering
        if (taper_ends) then
          residuals_for_cost(:) = residuals_for_cost(:) * w_tap(:)
          raw_residuals(:) = raw_residuals(:) * w_tap(:)
          wkstmp(:) = wkstmp(:) * w_tap(:)
        endif

        ! Filter
        if (use_band_pass_filter) then
          ! band pass filter
          fil_residuals(:) = 0._CUSTOM_REAL
          filfil_residuals(:) = 0._CUSTOM_REAL

          ! event filter limits
          fl = acqui_simu(ievent)%fl_event(current_ifrq)
          fh = acqui_simu(ievent)%fh_event(current_ifrq)

          ! filters for misfit calculation
          call bwfilt (raw_residuals, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

          ! filters twice for adjoint source
          ! note: filtered adjoint kernels would require to filter both the forward/reconstructed and adjoint wavefields.
          !       to avoid filtering the forward/reconstructed wavefield, it is equivalent to filtering twice only
          !       the adjoint source that produces the adjoint wavefield.
          call bwfilt (fil_residuals, filfil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

          ! filters for misfit calculation
          call bwfilt (residuals_for_cost, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

          ! filter data for storage
          call bwfilt(wkstmp, data_trace_to_use, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
        else
          ! no filter
          fil_residuals(:) = residuals_for_cost(:)
          filfil_residuals(:) = raw_residuals(:)
          data_trace_to_use(:) = wkstmp(:)
        endif

        !! save filtered data
        if (VERBOSE_MODE .or. DEBUG_MODE) then
           acqui_simu(ievent)%synt_traces(icomp, irec_local,:) = data_trace_to_use(:)
        endif

        !! sets weights
        if (set_event_trace_weights) then
          !! define energy renormalisation / trace energy
          !!$acqui_simu(ievent)%weight_trace(icomp,irec_local,1) = 100._CUSTOM_REAL / &
          !!$                       ((sum( acqui_simu(ievent)%synt_traces(icomp,irec_local,:) )**2) *0.5*dt_data)
          ! by prior, number of traces, trace length
          window_length = nstep_data * dt_data
          acqui_simu(ievent)%weight_trace(icomp,irec_local,1) = acqui_simu(ievent)%weight_trace(icomp,irec_local,1) &
                                   * 1._CUSTOM_REAL / prior_data_std / sqrt(nb_traces_tot) / sqrt(window_length)
        endif

        ! Apply weighting
        fil_residuals(:) = fil_residuals(:) * acqui_simu(ievent)%weight_trace(icomp,irec_local,1)
        filfil_residuals(:) = filfil_residuals(:) * acqui_simu(ievent)%weight_trace(icomp,irec_local,1)**2

        !! compute cost function value
        cost_value = sum(fil_residuals(:)**2) * 0.5 * dt_data
        cost_function = cost_function + cost_value

        !! compute raw standard deviation
        data_std = data_std + sum( residuals_for_cost(:)**2 )
        nb_data_std = nb_data_std + size(residuals_for_cost(:))

        ! tapers adjoint source (again to avoid non-zero onsets)
        if (taper_ends) then
          filfil_residuals(:) = filfil_residuals(:) * w_tap(:)
        endif

        !! store the adjoint source
        acqui_simu(ievent)%adjoint_sources(icomp,irec_local,:) = filfil_residuals(:)

        !(not used yet)
        !elastic_adjoint_source(icomp,:) = filfil_residuals(:)
      enddo

    case default
      write(*,*) 'Warning: unknown adjoint source type :',trim(adjustl(acqui_simu(ievent)%adjoint_source_type))
      stop
    end select

  end subroutine compute_elastic_adjoint_source_displacement

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
! define adjoint sources
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine compute_acoustic_adjoint_source_pressure_dot_dot(irec_local, ievent, acqui_simu)

    integer,                                     intent(in)    :: ievent, irec_local
    type(acqui),  dimension(:),                  intent(inout) :: acqui_simu
    ! local
    integer                                                    :: irec_glob
    real(kind=CUSTOM_REAL)                                     :: cost_value
    !! IN FLUIDS WITH need to save pressure and store the second time derivative of residuals
    !! raw_residuals(:)=seismograms_p(icomp,irec_local,:) - acqui_simu(ievent)%data_traces(irec_local,:,icomp)
    !! in fluid only one component : the pressure
    integer, parameter                                         :: P_ICOMP = 1

    !! TO DO : convolve seismograms_d(P_ICOMP,irec_local,:) by source time function stf (if need)
    !! if (acqui_simu(ievent)%convlove_residuals_by_wavelet) then
    !!     signal(:) = seismograms_d(P_ICOMP,irec_local,:)
    !!     call convolution_by_wavelet(wavelet, signal, seismograms_d(P_ICOMP,irec_local,:), nstep, nw)
    !! endif

    !! for acoustics need - sign (eg : Luo and Tromp Geophysics 2013)

    select case (trim(adjustl(acqui_simu(ievent)%adjoint_source_type)))

    case ('L2_FWI_TELESEISMIC')
      ! get data
      wkstmp(:) = acqui_simu(ievent)%data_traces(irec_local,:,P_ICOMP)

      ! residual (syn - obs)
      raw_residuals(:) = - (seismograms_p(P_ICOMP,irec_local,:) - wkstmp(:))

      ! save residuals because the adjoint source is not residuals
      residuals_for_cost(:) = raw_residuals(:)

      ! tapering
      if (taper_ends) then
        residuals_for_cost(:) = residuals_for_cost(:) * w_tap(:)
        raw_residuals(:) = raw_residuals(:) * w_tap(:)
        wkstmp(:) = wkstmp(:) * w_tap(:)
      endif

      ! choice: depending on acoustic formulation, we can either apply time derivatives or not for the adjoint source.
      !         given the kernel implementations of SPECFEM, we prefer to apply derivatives.
      if (.false.) then
        !! compute second time derivative of raw_residuals
        ! note: SPECFEM uses a potential formulation for fluid domains, as well as for the kernel expressions.
        !       this differs to an acoustic equation using a pressure formulation. the resulting adjoint kernels
        !       are only identical, if in the potential formulation, the adjoint source gets a second-order time derivative
        !       applied to the pressure misfit -(p_syn - p_obs).
        call FD2nd(raw_residuals, dt_data, NSTEP_DATA)
      else
        !! Here: we use the difference between observed pressure and computed pressure,
        !!       not the approach in Luo and Tromp Geophysics 2013,
        !!       which define the adjoint source as " minus second time derivatives of previous residuals "
        !!       We consider that the forward modeling is written in pressure thus
        !!       the adjoint is rho*displacement potential.
        !not needed, as the arrays are already the same: raw_residuals(:) = residuals_for_cost(:)
      endif

      ! Filter
      if (use_band_pass_filter) then
        ! bandpass filter
        fil_residuals(:) = 0._CUSTOM_REAL
        filfil_residuals(:) = 0._CUSTOM_REAL

        ! teleseismic filter limits
        fl = acqui_simu(ievent)%freqcy_to_invert(P_ICOMP,1,irec_local)
        fh = acqui_simu(ievent)%freqcy_to_invert(P_ICOMP,2,irec_local)

        ! filters twice for adjoint source
        ! note: filtered adjoint kernels would require to filter both the forward/reconstructed and adjoint wavefields.
        !       to avoid filtering the forward/reconstructed wavefield, it is equivalent to filtering twice only
        !       the adjoint source that produces the adjoint wavefield.
        call bwfilt (raw_residuals, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
        call bwfilt (fil_residuals, filfil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

        ! filters for misfit calculation
        call bwfilt (residuals_for_cost, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

        ! filter data for storage
        call bwfilt(wkstmp, data_trace_to_use, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
      else
        ! no filter
        fil_residuals(:) = residuals_for_cost(:)
        filfil_residuals(:) = raw_residuals(:)
        data_trace_to_use(:) = wkstmp(:)
      endif

      ! note: computing the derivatives for adjoint source;
      !       on the filtered residuals leads to less numerical artifacts
      if (.false.) then
        call FD2nd(filfil_residuals, dt_data, NSTEP_DATA)
      endif

      !! save filtered data
      if (VERBOSE_MODE .or. DEBUG_MODE) then
        acqui_simu(ievent)%synt_traces(P_ICOMP, irec_local,:) = data_trace_to_use(:)
      endif

      ! global receiver index
      irec_glob = acqui_simu(ievent)%number_receiver_global(irec_local)

      ! Apply weighting (for teleseismic, weights defined on all stations -> irec_glob)
      fil_residuals(:) = fil_residuals(:) * acqui_simu(ievent)%weight_trace(P_ICOMP,irec_glob,:)
      filfil_residuals(:) = filfil_residuals(:) * acqui_simu(ievent)%weight_trace(P_ICOMP,irec_glob,:)**2

      !! compute costs
      !! compute costs
      cost_value = sum(fil_residuals(:)**2) * 0.5 * dt_data
      cost_function = cost_function + cost_value

      !! compute raw standard deviation
      data_std = data_std + sum( residuals_for_cost(:)**2 )
      nb_data_std = nb_data_std + size(residuals_for_cost(:))

      ! tapers adjoint source (again to avoid non-zero onsets)
      if (taper_ends) then
        filfil_residuals(:) = filfil_residuals(:) * w_tap(:)
      endif

      !! store adjoint source
      acqui_simu(ievent)%adjoint_sources(1,irec_local,:) = filfil_residuals(:)

    case ('L2_OIL_INDUSTRY')
      ! get data
      wkstmp(:) = acqui_simu(ievent)%data_traces(irec_local,:,P_ICOMP)

      ! residual (syn - obs)
      raw_residuals(:) =  - (seismograms_p(P_ICOMP,irec_local,:) - wkstmp(:))

      ! save residuals because the adjoint source is not residuals
      residuals_for_cost(:) = raw_residuals(:)

      ! tapering
      if (taper_ends) then
        residuals_for_cost(:) = residuals_for_cost(:) * w_tap(:)
        raw_residuals(:) = raw_residuals(:) * w_tap(:)
        wkstmp(:) = wkstmp(:) * w_tap(:)
      endif

      ! choice: depending on acoustic formulation, we can either apply time derivatives or not for the adjoint source.
      !         given the kernel implementations of SPECFEM, we prefer to apply derivatives.
      if (.false.) then
        !! compute second time derivative of raw_residuals
        ! note: SPECFEM uses a potential formulation for fluid domains, as well as for the kernel expressions.
        !       this differs to an acoustic equation using a pressure formulation. the resulting adjoint kernels
        !       are only identical, if in the potential formulation, the adjoint source gets a second-order time derivative
        !       applied to the pressure misfit -(p_syn - p_obs).
        call FD2nd(raw_residuals, dt_data, NSTEP_DATA)
      else
        !! Here: we use the difference between observed pressure and computed pressure,
        !!       not the approach in Luo and Tromp Geophysics 2013,
        !!       which define the adjoint source as " minus second time derivatives of previous residuals "
        !!       We consider that the forward modeling is written in pressure thus
        !!       the adjoint is rho*displacement potential.
        !not needed, as the arrays are already the same: raw_residuals(:) = residuals_for_cost(:)
      endif

      ! Filter
      if (use_band_pass_filter) then
        ! bandpass filter
        fil_residuals(:) = 0._CUSTOM_REAL
        filfil_residuals(:) = 0._CUSTOM_REAL

        ! event filter limits
        fl = acqui_simu(ievent)%fl_event(current_ifrq)
        fh = acqui_simu(ievent)%fh_event(current_ifrq)

        ! filters twice for adjoint source
        ! note: filtered adjoint kernels would require to filter both the forward/reconstructed and adjoint wavefields.
        !       to avoid filtering the forward/reconstructed wavefield, it is equivalent to filtering twice only
        !       the adjoint source that produces the adjoint wavefield.
        call bwfilt (raw_residuals, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
        call bwfilt (fil_residuals, filfil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

        ! filter residuals for computing costs
        call bwfilt (residuals_for_cost, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

        ! filter data for storage
        call bwfilt (wkstmp, data_trace_to_use, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
      else
        ! no filter
        fil_residuals(:) = residuals_for_cost(:)
        filfil_residuals(:) = raw_residuals(:)
        data_trace_to_use(:) = wkstmp(:)
      endif

      ! note: computing the derivatives for adjoint source;
      !       on the filtered residuals leads to less numerical artifacts
      if (.false.) then
        call FD2nd(filfil_residuals, dt_data, NSTEP_DATA)
      endif

      ! save filtered data
      if (VERBOSE_MODE .or. DEBUG_MODE) then
        acqui_simu(ievent)%synt_traces(P_ICOMP, irec_local,:) = data_trace_to_use(:)
      endif

      !! sets weights
      if (set_event_trace_weights) then
        !! define energy renormalisation / trace energy
        !!$acqui_simu(ievent)%weight_trace(P_ICOMP,irec_local,:) = 100._CUSTOM_REAL / &
        !!$                  ((sum( acqui_simu(ievent)%synt_traces(P_ICOMP,irec_local,:) )**2) *0.5*dt_data)
        ! by prior, number of traces, trace length
        window_length = nstep_data * dt_data
        acqui_simu(ievent)%weight_trace(P_ICOMP,irec_local,1) = acqui_simu(ievent)%weight_trace(P_ICOMP,irec_local,1) &
                                 * 1._CUSTOM_REAL / prior_data_std / sqrt(nb_traces_tot) / sqrt(window_length)
      endif

      ! Apply weighting
      fil_residuals(:) = fil_residuals(:) * acqui_simu(ievent)%weight_trace(P_ICOMP,irec_local,1)
      filfil_residuals(:) = filfil_residuals(:) * acqui_simu(ievent)%weight_trace(P_ICOMP,irec_local,1)**2

      !! compute cost
      cost_value = sum(fil_residuals(:)**2) * 0.5 * dt_data
      cost_function = cost_function + cost_value

      !! compute raw standard deviation
      data_std = data_std + sum( residuals_for_cost(:)**2 )
      nb_data_std = nb_data_std + size(residuals_for_cost(:))

      ! tapers adjoint source (again to avoid non-zero onsets)
      if (taper_ends) then
        filfil_residuals(:) = filfil_residuals(:) * w_tap(:)
      endif

      !! store adjoint source
      acqui_simu(ievent)%adjoint_sources(1,irec_local,:) = filfil_residuals(:)

    case default
      write(*,*) 'Warning: unknown adjoint source type :',trim(adjustl(acqui_simu(ievent)%adjoint_source_type))
      stop
    end select
    !-------------------------------------------------------------------------------------------------

    !! TO DO : cross correlate filfil_residuals by source time function (if need)
    !! if (acqui_simu(ievent)%convlove_residuals_by_wavelet) then
    !!    signal(:) =  filfil_residuals(:);
    !!    call crosscor_by_wavelet(wavelet, signal, filfil_residuals, nstep, nw)
    !! endif

  end subroutine compute_acoustic_adjoint_source_pressure_dot_dot

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                         methods to interface with calling from external module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! unused so far...
!------------------------------------------
!subroutine get_filfil_residuals(sout)
!  real(kind=CUSTOM_REAL)     :: sout(*)
!  integer                    :: i
!  do i=1, NSTEP_DATA
!     sout(i)=filfil_residuals(i)
!  enddo
!end subroutine get_filfil_residuals

!------------------------------------------
!subroutine get_filfil_residuals_in_reverse(sout)
!  real(kind=CUSTOM_REAL)     :: sout(*)
!  integer                    :: i
!  do i = 1, NSTEP_DATA
!     sout(i) = filfil_residuals(NSTEP_DATA-i+1)
!  enddo
!end subroutine get_filfil_residuals_in_reverse

!------------------------------------------
!subroutine get_elastic_adj_src(sout,ir)
!  real(kind=CUSTOM_REAL) , dimension(:,:,:) :: sout
!  integer                    :: i, ic, ir
!  do i = 1, NSTEP_DATA
!     do ic = 1,3
!        sout(ic,ir,i) = elastic_adjoint_source(ic,NSTEP_DATA-i+1)
!     enddo
!  enddo
!end subroutine get_elastic_adj_src

!------------------------------------------
!subroutine get_acoustic_adj_src(sout,ir)
!  real(kind=CUSTOM_REAL), dimension(:,:,:) :: sout
!  integer                    :: i, ic, ir
!  do i = 1, NSTEP_DATA
!     do ic = 1,3
!        sout(ic,ir,i) = filfil_residuals(NSTEP_DATA-i+1)
!     enddo
!  enddo
!end subroutine get_acoustic_adj_src

!--------------------------------------------
!subroutine put_dt_nstep_for_adjoint_sources(delta_t, nb_time_step)
!  integer,                intent(in) :: nb_time_step
!  real(kind=CUSTOM_REAL), intent(in) :: delta_t
!  nstep_data = nb_time_step
!  dt_data = delta_t
!end subroutine put_dt_nstep_for_adjoint_sources

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
end module adjoint_source


