module adjoint_source

  !! IMPORT VARIABLES FROM SPECFEM -------------------------------------------------------------------------------------------------
  use specfem_par, only: CUSTOM_REAL, NDIM, MAX_STRING_LEN, OUTPUT_FILES, IIN, network_name, station_name, nrec_local, &
       seismograms_d, seismograms_p, number_receiver_global

  use specfem_par_elastic, only: ispec_is_elastic, ELASTIC_SIMULATION
  use specfem_par_acoustic, only: ispec_is_acoustic, ACOUSTIC_SIMULATION

  !!--------------------------------------------------------------------------------------------------------------------------------
  !! IMPORT inverse_problem VARIABLES
  use inverse_problem_par
  use signal_processing

  implicit none

  !! PRIVATE ATTRIBUTE -------------------------------------------------------------------------------------------------------------
  real(kind=CUSTOM_REAL), private, dimension(:,:), allocatable      :: elastic_adjoint_source, elastic_misfit
  real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: raw_residuals, fil_residuals, filfil_residuals
  real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: residuals_for_cost
  real(kind=CUSTOM_REAL), private, dimension(:), allocatable        :: signal, w_tap
  real(kind=CUSTOM_REAL), private                                   :: fl, fh
  real(kind=CUSTOM_REAL), private                                   :: dt_data
  real(kind=CUSTOM_REAL), private                                   :: cost_value
  integer,                private                                   :: norder_filter=4, irek_filter=1
  integer,                private                                   :: nstep_data


contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!  Write adjoint sources in SEM directory
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine write_adjoint_sources_for_specfem(acqui_simu,  inversion_param, isource, myrank)

    implicit none

    integer,                                     intent(in)     :: myrank, isource
    type(inver),                                 intent(inout)  :: inversion_param
    type(acqui),  dimension(:), allocatable,     intent(inout)  :: acqui_simu

    integer                                                     :: icomp
    integer                                                     :: irec, irec_local, ispec

    character(len=3),          dimension(NDIM)                  :: comp
    character(len=3)                                            :: pressure

    real(kind=CUSTOM_REAL)                                      :: cost_function, cost_function_reduced
    real(kind=CUSTOM_REAL)                                      :: cost_function_rec


    !! gets channel names for displacement
    do icomp=1,NDIM
       call write_channel_name(icomp,comp(icomp))
    enddo

    !! gets channel name for pressure
    call write_channel_name(4, pressure)

    nstep_data = acqui_simu(isource)%Nt_data
    dt_data = acqui_simu(isource)%dt_data

    !! initialize cost function for each MPI porcess
    cost_function  = 0._CUSTOM_REAL

    call allocate_adjoint_source_working_arrays()


    call taper_window_W(w_tap,10,30,nstep_data-30,nstep_data-10,nstep_data,1._CUSTOM_REAL)  !!!! WARNGING HARDCODED !!!!!!!!!!!!


    do irec_local = 1, nrec_local

       irec = number_receiver_global(irec_local)
       ispec = acqui_simu(isource)%ispec_selected_rec(irec)
       cost_function_rec=0.



       !! ALLOW TO CHOOSE COMPONENT :
       !! UX UY UZ (check if in solid region)
       !! Pr (check if in fluid region)
       !! VX VY VZ

       if (ELASTIC_SIMULATION) then
          if (ispec_is_elastic(ispec)) then

             !---------------------------------------------------------------------------------------------------------------
             !! compute adjoint source according to cost L2 function
             call compute_elastic_adjoint_source_displacement(irec_local, isource, acqui_simu, cost_function)

          endif
       endif

       !! IN FLUIDS WITH need to save pressure and store the second time derivative of residuals
       !! raw_residuals(:)=seismograms_p(icomp,irec_local,:) - acqui_simu(isource)%data_traces(irec_local,:,icomp)
       !! in fluid only one component : the pressure
       icomp=1
       if (ACOUSTIC_SIMULATION) then
          if (ispec_is_acoustic(ispec)) then

             !! ---------------------------------------------------------------------------------------------------------------
             !! compute adjoint source according to cost L2 function
             call compute_acoustic_adjoint_source_pressure_dot_dot(icomp, irec_local, isource, acqui_simu, cost_function)

          endif

       endif

    enddo


    !! compute cost function   : allreduce cost_function
    cost_function_reduced=0._CUSTOM_REAL
    call sum_all_all_cr(cost_function, cost_function_reduced)

    !! add the cost function over all sources
    inversion_param%total_current_cost =inversion_param%total_current_cost + cost_function_reduced

    !! save cost function for the current source
    inversion_param%current_cost(isource) = cost_function_reduced

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*) '      Cost function for this source : ', cost_function_reduced
    endif

    call deallocate_adjoint_source_working_arrays()

  end subroutine write_adjoint_sources_for_specfem

!----------------------------------------------------------------------------------------------------------------------------------
  subroutine deallocate_adjoint_source_working_arrays()
    deallocate(raw_residuals, fil_residuals,  filfil_residuals, w_tap, signal, residuals_for_cost, elastic_adjoint_source, &
         elastic_misfit)
  end subroutine deallocate_adjoint_source_working_arrays

!----------------------------------------------------------------------------------------------------------------------------------
  subroutine  allocate_adjoint_source_working_arrays()

    allocate(raw_residuals(nstep_data), fil_residuals(nstep_data), filfil_residuals(nstep_data), &
         w_tap(nstep_data), signal(nstep_data), residuals_for_cost(nstep_data), &
         elastic_adjoint_source(NDIM,nstep_data), elastic_misfit(NDIM,nstep_data))

  end subroutine allocate_adjoint_source_working_arrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
! define adjoint sources
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine compute_elastic_adjoint_source_displacement(irec_local, isource, acqui_simu, cost_function)


    integer,                                     intent(in)    :: isource, irec_local
    type(acqui),  dimension(:), allocatable,     intent(inout) :: acqui_simu
    real(kind=CUSTOM_REAL),                      intent(inout) :: cost_function
    integer                                                    :: icomp
    !!----------------------------------------------------------------------------------------------------
    !! store residuals and filter  ---------------------------

    !! TO DO : convolve seismograms_d(icomp,irec_local,:) by source time function stf (if need)
    !! if (acqui_simu(isource)%convlove_residuals_by_wavelet) then
    !!     signal(:) = seismograms_d(icomp,irec_local,:)
    !!     call convolution_by_wavelet(wavelet, signal, seismograms_d(icomp,irec_local,:), nstep, nw)
    !! endif


    !! TO DO : for now only L2 fwi is used, we should consider other adjoint sources
    !!
    !!   -- travel time kernel
    !!   -- Enveloppe
    !!   -- rotation to change components
    !!

    select case (trim(adjustl(acqui_simu(isource)%adjoint_source_type)))

    case ('L2_FWI_TELESEISMIC', 'L2_OIL_INDUSTRY')

       do icomp = 1, NDIM

          raw_residuals(:)=seismograms_d(icomp,irec_local,:) - acqui_simu(isource)%data_traces(irec_local,:,icomp)

          fil_residuals(:)=0._CUSTOM_REAL
          filfil_residuals(:)=0._CUSTOM_REAL
          fl=acqui_simu(isource)%freqcy_to_invert(icomp,1,irec_local)
          fh=acqui_simu(isource)%freqcy_to_invert(icomp,2,irec_local)
          call bwfilt (raw_residuals, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
          call bwfilt (fil_residuals, filfil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)


          !! TO DO : cross correlate filfil_residuals by source time function (if need)
          !! if (acqui_simu(isource)%convlove_residuals_by_wavelet) then
          !!    signal(:) =  filfil_residuals(:);
          !!    call crosscor_by_wavelet(wavelet, signal, filfil_residuals, nstep, nw)
          !! endif


          !! TO DO : choose component to invert

          !! remove component if not used
          if (trim(acqui_simu(isource)%component(icomp)) == '0' .or. &
               trim(acqui_simu(isource)%component(icomp)) == '  ' ) then
             filfil_residuals(:)=0._CUSTOM_REAL
             fil_residuals(:)=0._CUSTOM_REAL
          endif

          !! compute cost function value
          cost_value=sum(fil_residuals(:)**2) * 0.5 * dt_data
          cost_function = cost_function + cost_value

          !! TO DO : cross correlate filfil_residuals by source time function (if need)
          !! if (acqui_simu(isource)%convlove_residuals_by_wavelet) then
          !!    signal(:) =  filfil_residuals(:);
          !!    call crosscor_by_wavelet(wavelet, signal, filfil_residuals, nstep, nw)
          !! endif


          !!----------------------------------------------------------------------------------------------------


          !! store the adjoint source
          elastic_adjoint_source(icomp,:) = filfil_residuals(:)
          acqui_simu(isource)%adjoint_sources(icomp,irec_local,:) = filfil_residuals(:)
       enddo

       !!----------------------------------------------------------------------------------------------------

       !case ('L2_OIL_INDUSTRY')
          !! filter the source wavelet then the synthetics are in the good frequency range
          !! need just to filter the data trace


       case default

       end select



  end subroutine compute_elastic_adjoint_source_displacement

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
! define adjoint sources
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine compute_acoustic_adjoint_source_pressure_dot_dot(icomp, irec_local, isource, acqui_simu, cost_function)

    integer,                                     intent(in)    :: isource, icomp, irec_local
    type(acqui),  dimension(:), allocatable,     intent(inout) :: acqui_simu
    real(kind=CUSTOM_REAL),                      intent(inout) :: cost_function
    !!--------------------------------------------------------------------------------------------------
    !! TO DO : convolve seismograms_d(icomp,irec_local,:) by source time function stf (if need)
    !! if (acqui_simu(isource)%convlove_residuals_by_wavelet) then
    !!     signal(:) = seismograms_d(icomp,irec_local,:)
    !!     call convolution_by_wavelet(wavelet, signal, seismograms_d(icomp,irec_local,:), nstep, nw)
    !! endif
    !! for acoustics need - sign (eg : Luo et al Geophysics 2013)

    select case (trim(adjustl(acqui_simu(isource)%adjoint_source_type)))

    case ('L2_FWI_TELESEISMIC', 'L2_OIL_INDUSTRY' )
       raw_residuals(:)=-(seismograms_p(icomp,irec_local,:) - acqui_simu(isource)%data_traces(irec_local,:,icomp))
       residuals_for_cost(:) = raw_residuals(:) !! save residuals because the adjoint source is not residuals

       !! compute second time derivative of raw_residuals
       call FD2nd(raw_residuals, dt_data, NSTEP_DATA)

       fil_residuals(:)=0._CUSTOM_REAL
       filfil_residuals(:)=0._CUSTOM_REAL
       fl=acqui_simu(isource)%freqcy_to_invert(icomp,1,irec_local)
       fh=acqui_simu(isource)%freqcy_to_invert(icomp,2,irec_local)
       call bwfilt (raw_residuals, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
       call bwfilt (fil_residuals, filfil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)

       !! compute residuals
       call bwfilt (residuals_for_cost, fil_residuals, dt_data, nstep_data, irek_filter, norder_filter, fl, fh)
       cost_value=sum(fil_residuals(:)**2) * 0.5 * dt_data
       cost_function = cost_function + cost_value

       !! store adjoint source
       acqui_simu(isource)%adjoint_sources(1,irec_local,:)=filfil_residuals(:)

    !case ('L2_OIL_INDUSTRY')

    case default

    end select
    !-------------------------------------------------------------------------------------------------

    !! TO DO : cross correlate filfil_residuals by source time function (if need)
    !! if (acqui_simu(isource)%convlove_residuals_by_wavelet) then
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


