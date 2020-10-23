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
  real(kind=CUSTOM_REAL), private, dimension(:,:), allocatable      :: elastic_adjoint_source, elastic_misfit
  real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: raw_residuals, fil_residuals, filfil_residuals
  real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: residuals, data_trace_to_use, wkstmp
  real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: residuals_for_cost
  real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: signal, w_tap
  real(kind=CUSTOM_REAL), private                                   :: fl, fh
  real(kind=CUSTOM_REAL), private                                   :: dt_data
  real(kind=CUSTOM_REAL), private                                   :: cost_value
  real(kind=CUSTOM_REAL), private                                   :: prior_data_std, data_std, nb_data_std
  real(kind=CUSTOM_REAL), private                                   :: nb_traces_tot, window_lenght
  integer,                private                                   :: norder_filter=4, irek_filter=1
  integer,                private                                   :: nstep_data
  integer,                private                                   :: current_ifrq
  logical,                private                                   :: use_band_pass_filter

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!  Write adjoint sources in SEM directory
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine write_adjoint_sources_for_specfem(acqui_simu,  inversion_param, ievent, myrank)

    implicit none

    integer,                                     intent(in)     :: myrank, ievent
    type(inver),                                 intent(inout)  :: inversion_param
    type(acqui),  dimension(:), allocatable,     intent(inout)  :: acqui_simu

    integer                                                     :: icomp
    integer                                                     :: irec, irec_local, ispec

    real(kind=CUSTOM_REAL)                                      :: cost_function, cost_function_reduced
    real(kind=CUSTOM_REAL)                                      :: cost_function_rec


    integer                                                     :: lw, i0, i1, i2, i3
    integer                                                     :: current_iter

    nstep_data = acqui_simu(ievent)%Nt_data
    dt_data = acqui_simu(ievent)%dt_data

    current_iter=inversion_param%current_iteration
    current_ifrq=inversion_param%current_ifrq
    use_band_pass_filter=inversion_param%use_band_pass_filter
    prior_data_std = inversion_param%prior_data_std
    nb_traces_tot =  inversion_param%nb_traces_tot

    !! initialize cost function for each MPI porcess
    cost_function  = 0._CUSTOM_REAL
    data_std = 0._CUSTOM_REAL
    nb_data_std =  0._CUSTOM_REAL

    call allocate_adjoint_source_working_arrays()

    !! define taper on adjoint sources (if not window selected by user)
    lw=0.02*nstep_data  !!!! WARNGING HARDCODED !!!!!!!!!!!!
    i0=10
    i1=i0 + lw
    i3=nstep_data-10
    i2=i3 - lw
    call taper_window_W(w_tap,i0,i1,i2,i3,nstep_data,1._CUSTOM_REAL)
    !! to do define user window

    do irec_local = 1, nrec_local
       ! check: assumes SIMULATION_TYPE == 3, with adjoint sources and receivers being identical
       !        -  nrec_local ==  nadj_rec_local
       !        -  number_receiver_global == number_adjsources_global
       irec = number_receiver_global(irec_local)
       ispec = acqui_simu(ievent)%ispec_selected_rec(irec)
       cost_function_rec=0.

       !! ALLOW TO CHOOSE COMPONENT :
       !! UX UY UZ (check if in solid region)
       !! Pr (check if in fluid region)
       !! VX VY VZ

       if (ELASTIC_SIMULATION) then
          if (ispec_is_elastic(ispec)) then

             !! ---------------------------------------------------------------------------------------------------------------
             !! compute adjoint source according to cost L2 function
             call compute_elastic_adjoint_source_displacement(irec_local, ievent, current_iter, &
                  acqui_simu, cost_function, inversion_param)

          endif
       endif

       !! IN FLUIDS WITH need to save pressure and store the second time derivative of residuals
       !! raw_residuals(:)=seismograms_p(icomp,irec_local,:) - acqui_simu(ievent)%data_traces(irec_local,:,icomp)
       !! in fluid only one component : the pressure
       icomp=1
       if (ACOUSTIC_SIMULATION) then
          if (ispec_is_acoustic(ispec)) then

             !! ---------------------------------------------------------------------------------------------------------------
             !! compute adjoint source according to cost L2 function
             call compute_acoustic_adjoint_source_pressure_dot_dot(icomp, irec_local, ievent, acqui_simu, cost_function)

          endif

       endif

    enddo


    !! compute cost function   : allreduce cost_function
    cost_function_reduced=0._CUSTOM_REAL
    call sum_all_all_cr(cost_function, cost_function_reduced)
    !! add the cost function over all sources in group
    inversion_param%total_current_cost =inversion_param%total_current_cost + cost_function_reduced
    !! save cost function for the current source
    inversion_param%current_cost(ievent) = cost_function_reduced

    if (myrank == 0 .and. (VERBOSE_MODE .or. DEBUG_MODE) ) then
       write(INVERSE_LOG_FILE,*) '      Cost function for this event : ', cost_function_reduced
       write(INVERSE_LOG_FILE,*) '     weight on data : ', 1._CUSTOM_REAL/prior_data_std
       write(INVERSE_LOG_FILE,*) '     number of traces  : ',  nb_traces_tot
       write(INVERSE_LOG_FILE,*) '     total weight      : ',   1._CUSTOM_REAL/&
            sqrt(nb_traces_tot)/prior_data_std/sqrt(nstep_data*dt_data)
    endif

    !! standard deviation on data
    cost_function_reduced=0._CUSTOM_REAL
    call sum_all_all_cr(data_std, cost_function_reduced)
    inversion_param%data_std = inversion_param%data_std + data_std

    cost_function_reduced=0._CUSTOM_REAL
    call sum_all_all_cr(nb_data_std, cost_function_reduced)
    inversion_param%nb_data_std = inversion_param%nb_data_std + nb_data_std

    call deallocate_adjoint_source_working_arrays()

  end subroutine write_adjoint_sources_for_specfem

!----------------------------------------------------------------------------------------------------------------------------------
  subroutine deallocate_adjoint_source_working_arrays()
    deallocate(residuals, raw_residuals, fil_residuals,  filfil_residuals, w_tap, signal, residuals_for_cost, &
         elastic_adjoint_source, &
         elastic_misfit, data_trace_to_use, wkstmp)
  end subroutine deallocate_adjoint_source_working_arrays

!----------------------------------------------------------------------------------------------------------------------------------
  subroutine  allocate_adjoint_source_working_arrays()

    integer :: ier

    allocate(residuals(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 279')
    allocate(raw_residuals(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 280')
    allocate(fil_residuals(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 281')
    allocate(filfil_residuals(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 282')
    allocate(w_tap(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 283')
    allocate(signal(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 284')
    allocate(residuals_for_cost(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 285')
    allocate(elastic_adjoint_source(NDIM,nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 286')
    allocate(elastic_misfit(NDIM,nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 287')
    allocate(data_trace_to_use(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 288')
    allocate(wkstmp(nstep_data),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 289')

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
!    3/    compute the adjoint source whcih is stored in acqui(ievent)%adjoint_sources(NCOM, NREC_LOCAL, NT)
!          note that this arrays will be directly use be specfem as adjoint source thus you need to
!          do any proccessing here : filter, rotation, ....
!
!    4/    compute the cost function and store it in cost_function variable (you have to perform sommation over sources)
!
!
!
  subroutine compute_elastic_adjoint_source_displacement(irec_local, ievent, current_iter, acqui_simu, cost_function, &
       inversion_param)


    integer,                                     intent(in)    :: ievent, irec_local, current_iter
    type(acqui),  dimension(:), allocatable,     intent(inout) :: acqui_simu
    real(kind=CUSTOM_REAL),                      intent(inout) :: cost_function
    integer                                                    :: icomp, idim, irec_glob !, icomp_tmp
    real(kind=CUSTOM_REAL), dimension(:), allocatable          :: wavelet, filfil_residuals_tmp, tmpl
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable        :: trace_cal_1, trace_cal_2
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable        :: trace_obs_1, trace_obs_2
    double precision                                           :: lat0, lon0, azi0
    type(inver),                                 intent(inout) :: inversion_param

    integer :: ier

    !!----------------------------------------------------------------------------------------------------
    !! store residuals and filter  ---------------------------

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
       irec_glob = acqui_simu(ievent)%number_receiver_global(irec_local)

       do idim = 1, ndim

          ! Get data
          trace_cal_1(idim,:) = seismograms_d(idim,irec_local,:)
          trace_cal_2(idim,:) = seismograms_d(idim,irec_local,:)

          trace_obs_1(idim,:) = acqui_simu(ievent)%data_traces(irec_local,:,idim)
          trace_obs_2(idim,:) = acqui_simu(ievent)%data_traces(irec_local,:,idim)

          ! Convolve synthetic data with wavelet
          if (inversion_param%convolution_by_wavelet) then
             if (.not. allocated(wavelet)) then
               allocate(wavelet(nstep_data),stat=ier)
               if (ier /= 0) call exit_MPI_without_rank('error allocating array 294')
             endif
             wavelet = acqui_simu(ievent)%user_source_time_function(1,:)
             call myconvolution(trace_cal_2(idim,:),wavelet,nstep_data,nstep_data,tmpl,0)
             trace_cal_1(idim,:) = tmpl * dt_data
          endif

       enddo

       ! Do rotation of data (would be easier is lat and lon are split over MPI slices)
       select case (inversion_param%inverted_data_sys)
       case ('xyz')
          ! Do nothing
       case('enz')
          trace_cal_2 = trace_cal_1
          trace_obs_2 = trace_obs_1
          ! Data are in standard coordinate system
          ! Data rotation required to pass in mesh system (zen -> xyz)
          call define_mesh_rotation_matrix(lat0,lon0,azi0)
          call rotate_comp_mesh2glob(trace_cal_2(1,:), trace_cal_2(2,:), trace_cal_2(3,:), &
               acqui_simu(ievent)%read_station_position(1,irec_glob), &
               acqui_simu(ievent)%read_station_position(2,irec_glob), &
               nstep_data, 1, trace_cal_1(3,:), trace_cal_1(2,:), trace_cal_1(1,:))
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
          call rotate_comp_mesh2glob(trace_cal_2(1,:), trace_cal_2(2,:), trace_cal_2(3,:), &
               acqui_simu(ievent)%read_station_position(1,irec_glob), &
               acqui_simu(ievent)%read_station_position(2,irec_glob), &
               nstep_data, 1, trace_cal_1(3,:), trace_cal_1(2,:), trace_cal_1(1,:))
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
          write(*,*)'CATASTROPHIC ERROR'
          write(*,*)'qtl is not implemented yet'
          write(*,*)'NOW STOP'
          stop
       end select

       ! Finally compute residuals, filter and cross-correlate
       do idim = 1, ndim

          ! Resiudal
          raw_residuals(:) =  trace_cal_1(idim,:) - trace_obs_1(idim,:)

          ! Filter
          fil_residuals(:)=0._CUSTOM_REAL
          filfil_residuals(:)=0._CUSTOM_REAL
          fl=acqui_simu(ievent)%freqcy_to_invert(idim,1,irec_local)
          fh=acqui_simu(ievent)%freqcy_to_invert(idim,2,irec_local)
          call bwfilt (raw_residuals, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
          call bwfilt (fil_residuals, filfil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

          ! Apply weighting
          fil_residuals(:) = fil_residuals(:) * acqui_simu(ievent)%weight_trace(idim,irec_glob,:)
          filfil_residuals(:) = filfil_residuals(:) * acqui_simu(ievent)%weight_trace(idim,irec_glob,:)**2

          !! compute cost function value
          cost_value=sum(fil_residuals(:)**2) * 0.5 * dt_data
          cost_function = cost_function + cost_value

          ! Finally cross-correlate residuals with wavelet
          if (inversion_param%convolution_by_wavelet) then
             if (.not. allocated(filfil_residuals_tmp)) then
               allocate(filfil_residuals_tmp(nstep_data),stat=ier)
               if (ier /= 0) call exit_MPI_without_rank('error allocating array 295')
             endif
             filfil_residuals_tmp(:) = filfil_residuals(:)
             call mycorrelation(filfil_residuals_tmp,wavelet,nstep_data,nstep_data,tmpl,0)
             filfil_residuals = tmpl * dt_data
          endif

          !! store the adjoint source
          elastic_adjoint_source(icomp,:) = filfil_residuals(:)
          acqui_simu(ievent)%adjoint_sources(icomp,irec_local,:) = filfil_residuals(:) !*w_tap(:)
       enddo

       !!----------------------------------------------------------------------------------------------------

       case ('L2_OIL_INDUSTRY')

          do icomp = 1, NDIM

             !! get data ------------------------
             if (use_band_pass_filter) then
                !! filter the data
                fil_residuals(:)=0._CUSTOM_REAL
                fl=acqui_simu(ievent)%fl_event(current_ifrq)
                fh=acqui_simu(ievent)%fh_event(current_ifrq)

                wkstmp(:)= acqui_simu(ievent)%data_traces(irec_local,:,icomp)
                call bwfilt(wkstmp, data_trace_to_use, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

                if (VERBOSE_MODE .or. DEBUG_MODE) then
                   !! save filtered data
                   acqui_simu(ievent)%synt_traces(icomp, irec_local,:)=  data_trace_to_use(:)
                endif

             else
                data_trace_to_use(:)=acqui_simu(ievent)%data_traces(irec_local,:,icomp)
             endif

             !! define energy renormalisation
             if (current_iter == 0) then
!!$                do icomp_tmp = 1, NDIM
!!$                   acqui_simu(ievent)%weight_trace(icomp,irec_local)=100._CUSTOM_REAL / &
!!$                        ((sum( acqui_simu(ievent)%synt_traces(irec_local,:,icomp_tmp) )**2) *0.5*dt_data)
!!$                enddo
                window_lenght =  nstep_data * dt_data
                acqui_simu(ievent)%weight_trace(icomp,irec_local,1)=1._CUSTOM_REAL/prior_data_std/&
                     sqrt(nb_traces_tot)/sqrt(window_lenght)
             endif

             !! compute residuals residuals
             residuals(:)= (seismograms_d(icomp,irec_local,:) - data_trace_to_use(:))*&
                  acqui_simu(ievent)%weight_trace(icomp,irec_local,1)

             !! compute cost
             cost_value=sum(residuals(:)**2) * 0.5 * dt_data
             cost_function = cost_function + cost_value

             !! compute raw standard deviation
             data_std = data_std + sum((seismograms_d(icomp,irec_local,:) - data_trace_to_use(:))**2 )
             nb_data_std = nb_data_std + size(residuals(:))

             ! store adjoint source
             acqui_simu(ievent)%adjoint_sources(icomp,irec_local,:)=residuals(:)*w_tap(:)*&
                  acqui_simu(ievent)%weight_trace(icomp,irec_local,1)

          enddo

       case default

       end select



  end subroutine compute_elastic_adjoint_source_displacement

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
! define adjoint sources
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine compute_acoustic_adjoint_source_pressure_dot_dot(icomp, irec_local, ievent, acqui_simu, cost_function)

    integer,                                     intent(in)    :: ievent, icomp, irec_local
    type(acqui),  dimension(:), allocatable,     intent(inout) :: acqui_simu
    real(kind=CUSTOM_REAL),                      intent(inout) :: cost_function
    !!--------------------------------------------------------------------------------------------------
    !! TO DO : convolve seismograms_d(icomp,irec_local,:) by source time function stf (if need)
    !! if (acqui_simu(ievent)%convlove_residuals_by_wavelet) then
    !!     signal(:) = seismograms_d(icomp,irec_local,:)
    !!     call convolution_by_wavelet(wavelet, signal, seismograms_d(icomp,irec_local,:), nstep, nw)
    !! endif
    !! for acoustics need - sign (eg : Luo and Tromp Geophysics 2013)

    select case (trim(adjustl(acqui_simu(ievent)%adjoint_source_type)))

    case ('L2_FWI_TELESEISMIC')
       raw_residuals(:)=-(seismograms_p(icomp,irec_local,:) - acqui_simu(ievent)%data_traces(irec_local,:,icomp))
       residuals_for_cost(:) = raw_residuals(:) !! save residuals because the adjoint source is not residuals

       !! compute second time derivative of raw_residuals
       call FD2nd(raw_residuals, dt_data, NSTEP_DATA)

       fil_residuals(:)=0._CUSTOM_REAL
       filfil_residuals(:)=0._CUSTOM_REAL
       fl=acqui_simu(ievent)%freqcy_to_invert(icomp,1,irec_local)
       fh=acqui_simu(ievent)%freqcy_to_invert(icomp,2,irec_local)
       call bwfilt (raw_residuals, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
       call bwfilt (fil_residuals, filfil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

       !! compute residuals
       call bwfilt (residuals_for_cost, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
       cost_value=sum(fil_residuals(:)**2) * 0.5 * dt_data
       cost_function = cost_function + cost_value

       !! store adjoint source
       acqui_simu(ievent)%adjoint_sources(1,irec_local,:)=filfil_residuals(:)*w_tap(:)

    case ('L2_OIL_INDUSTRY')

       if (use_band_pass_filter) then
          !! filter the data
          fil_residuals(:)=0._CUSTOM_REAL
          fl=acqui_simu(ievent)%fl_event(current_ifrq)
          fh=acqui_simu(ievent)%fh_event(current_ifrq)
          raw_residuals(:)= acqui_simu(ievent)%data_traces(irec_local,:,icomp)
          call bwfilt(raw_residuals, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

          !! save filtered data
          acqui_simu(ievent)%synt_traces(icomp, irec_local,:)= fil_residuals(:)
       else
          fil_residuals(:)=acqui_simu(ievent)%data_traces(irec_local,:,icomp)
       endif

       !! save residuals for adjoint source. Note we use the difference between
       !! obseved pressure and computed pressure, not the approach in Luo and Tromp Gepohysics 2013
       !! which define the adjoint source as " minus second time derivatives of previous residuals "
       !! We consider that the forward modeling is writen in pressure thus
       !! the adjoint is rho*displacement potential.
       residuals_for_cost(:) =  - (seismograms_p(icomp,irec_local,:) - fil_residuals(:))

       !! compute cost
       cost_value=sum(residuals_for_cost(:)**2) * 0.5 * dt_data
       cost_function = cost_function + cost_value

       !! store adjoint source
       acqui_simu(ievent)%adjoint_sources(1,irec_local,:)=residuals_for_cost(:)*w_tap(:)


    case default

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
!------------------------------------------
subroutine get_filfil_residuals(sout)

  real(kind=CUSTOM_REAL)     :: sout(*)
  integer                    :: i
  do i=1, NSTEP_DATA
     sout(i)=filfil_residuals(i)
  enddo

end subroutine get_filfil_residuals
!------------------------------------------
subroutine get_filfil_residuals_in_reverse(sout)

  real(kind=CUSTOM_REAL)     :: sout(*)
  integer                    :: i
  do i=1, NSTEP_DATA
     sout(i)=filfil_residuals(NSTEP_DATA-i+1)
  enddo

end subroutine get_filfil_residuals_in_reverse
!------------------------------------------
subroutine get_elastic_adj_src(sout,ir)

  real(kind=CUSTOM_REAL) , dimension(:,:,:),allocatable     :: sout
  integer                    :: i, ic, ir
  do i=1, NSTEP_DATA
     do ic=1,3
        sout(ic,ir,i)=elastic_adjoint_source(ic,NSTEP_DATA-i+1)
     enddo
  enddo
end subroutine get_elastic_adj_src
!------------------------------------------
subroutine get_acoustic_adj_src(sout,ir)

  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable     :: sout
  integer                    :: i, ic, ir
  do i=1, NSTEP_DATA
     do ic=1,3
        sout(ic,ir,i)=filfil_residuals(NSTEP_DATA-i+1)
     enddo
  enddo
end subroutine get_acoustic_adj_src

!--------------------------------------------
subroutine put_dt_nstep_for_adjoint_sources(delta_t, nb_time_step)
  integer,                intent(in) :: nb_time_step
  real(kind=CUSTOM_REAL), intent(in) :: delta_t
  nstep_data=nb_time_step
  dt_data=delta_t
end subroutine put_dt_nstep_for_adjoint_sources

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
end module adjoint_source


