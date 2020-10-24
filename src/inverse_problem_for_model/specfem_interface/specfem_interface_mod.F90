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

module specfem_interface

  !! IMPORT SPECFEM ALL VARAIBLES
  use specfem_par
  use shared_parameters
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  !! IMPORT inverse_problem VARIABLES
  use inverse_problem_par

  !! IMPORT inverse problem modules
  use adjoint_source
  use input_output
  use signal_processing

  implicit none

  !! Arrays for saving GPU kernels because at each GPU run the kernels are set to 0 and to perform
  !! summation over events we need to use those temporarry arrays
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable :: rho_ac_kl_GPU, kappa_ac_kl_GPU
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable :: rho_kl_GPU, kappa_kl_GPU, mu_kl_GPU
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable :: hess_rho_ac_kl_GPU, hess_kappa_ac_kl_GPU
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable :: hess_rho_kl_GPU, hess_kappa_kl_GPU, hess_mu_kl_GPU
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable :: cijkl_kl_GPU

contains

!----------------------------------------------------------------------------------------------------------------------------------
!> call specfem solver for forward problem only
!----------------------------------------------------------------------------------------------------------------------------------

  subroutine ComputeSismosPerEvent(ievent, acqui_simu, iter_inverse, inversion_param, myrank)

  implicit none

  integer,                                        intent(in)    ::  ievent, iter_inverse, myrank
  type(acqui),  dimension(:), allocatable,        intent(inout) ::  acqui_simu
  type(inver),                                    intent(inout) ::  inversion_param
  character(len=MAX_LEN_STRING)                                 ::  name_file_tmp

  !! set the parameters to perform the forward simulation
  SIMULATION_TYPE=1
  SAVE_FORWARD=.true.
  COMPUTE_AND_STORE_STRAIN = .false.
  SAVE_MESH_FILES=.false.

  !! forward solver
  call InitSpecfemForOneRun(acqui_simu, ievent, inversion_param, iter_inverse)
  call iterate_time()
  call FinalizeSpecfemForOneRun(acqui_simu, ievent)

  select case(trim(type_input))
  case ('teleseismic')
     name_file_tmp = 'data'
     call write_pif_data_gather(ievent, acqui_simu, inversion_param, seismograms_d, name_file_tmp, myrank)

  case default
     select case ( trim(acqui_simu(ievent)%component(1)) )
     case('UX', 'UY', 'UZ')
        !! array seismogram in displacement
        name_file_tmp = trim(acqui_simu(ievent)%data_file_gather)

        write(INVERSE_LOG_FILE,*) '  ... Writing simulated data gather for event :', ievent

        call write_bin_sismo_on_disk(ievent, acqui_simu, seismograms_d,  name_file_tmp, myrank)

     case('PR')
        !! array seismogram in pressure
        name_file_tmp = trim(acqui_simu(ievent)%data_file_gather)

        write(INVERSE_LOG_FILE,*) '  ... Writing simulated data gather for event :  ', ievent

        call write_bin_sismo_on_disk(ievent, acqui_simu, seismograms_p,  name_file_tmp, myrank)

     case default
        write(*,*) ' ERROR Component not known : ', trim(acqui_simu(ievent)%component(1))
        stop
     end select
  end select

  end subroutine ComputeSismosPerEvent


!-----------------------------------------------------------------------------------------------------------------
!> call specfem solver for computing gradient
!----------------------------------------------------------------------------------------------------------------

  subroutine ComputeGradientPerEvent(ievent, iter_inverse, acqui_simu,  inversion_param)

  implicit none

  integer,                                        intent(in)    ::  ievent, iter_inverse
  type(acqui),  dimension(:), allocatable,        intent(inout) ::  acqui_simu
  type(inver),                                    intent(inout) ::  inversion_param
  character(len=MAX_LEN_STRING)                                 ::  name_file_tmp

  logical                                                       :: save_COUPLE_WITH_INJECTION_TECHNIQUE

  save_COUPLE_WITH_INJECTION_TECHNIQUE = COUPLE_WITH_INJECTION_TECHNIQUE

  if (myrank == 0) write(INVERSE_LOG_FILE,*) '  - > Compute gradient for event :', ievent , ' iteration : ', iter_inverse

  !! choose parameters to perform the forward simulation
  SIMULATION_TYPE=1
  SAVE_FORWARD=.true.
  COMPUTE_AND_STORE_STRAIN = .false.

  !! forward solver
  call InitSpecfemForOneRun(acqui_simu, ievent, inversion_param, iter_inverse)

  call iterate_time()

  call FinalizeSpecfemForOneRun(acqui_simu, ievent)

  !! define adjoint sources
  call write_adjoint_sources_for_specfem(acqui_simu, inversion_param, ievent, myrank)

  select case(trim(type_input))
  case ('teleseismic')
     name_file_tmp = 'adjoint_source'
     call write_pif_data_gather(ievent, acqui_simu, inversion_param, acqui_simu(ievent)%synt_traces, name_file_tmp, myrank)

  case default
     !! dump synthetics and adjoint sources to ckeck
     if (VERBOSE_MODE .or. DEBUG_MODE) then
        call dump_adjoint_sources(iter_inverse, ievent, acqui_simu, myrank)

        select case (trim(acqui_simu(ievent)%component(1)))
        case('UX', 'UY', 'UZ')
           call dump_seismograms(iter_inverse, ievent, seismograms_d, acqui_simu, myrank)
           call dump_filtered_data(iter_inverse,ievent,acqui_simu(ievent)%synt_traces, acqui_simu, myrank)
        case('PR')
           call dump_seismograms(iter_inverse, ievent, seismograms_p, acqui_simu, myrank)
           call dump_filtered_data(iter_inverse,ievent,acqui_simu(ievent)%synt_traces, acqui_simu, myrank)
        end select
     endif
  end select

  !! choose parameters to perform both the forward and adjoint simulation
  SIMULATION_TYPE = 3
  SAVE_FORWARD = .false.
  COMPUTE_AND_STORE_STRAIN = .true.
  APPROXIMATE_HESS_KL = .true.

  !! forward and adjoint runs
  call InitSpecfemForOneRun(acqui_simu, ievent, inversion_param, iter_inverse)

  COUPLE_WITH_INJECTION_TECHNIQUE = .false.  !! do not use coupling since the direct run is runining in backward from boundary

  call iterate_time()
  call FinalizeSpecfemForOneRun(acqui_simu, ievent)

  COUPLE_WITH_INJECTION_TECHNIQUE = save_COUPLE_WITH_INJECTION_TECHNIQUE  !! restore the initial value of variable

  end subroutine ComputeGradientPerEvent


!----------------------------------------------------------------------------------------------------------------------------------
!> initialize specfem before each call for the event ievent
!----------------------------------------------------------------------------------------------------------------------------------

  subroutine InitSpecfemForOneRun(acqui_simu, ievent, inversion_param, iter_inverse)

  implicit none

  integer,                                        intent(in)    ::  ievent, iter_inverse
  type(inver),                                    intent(in)    ::  inversion_param
  type(acqui),  dimension(:), allocatable,        intent(in)    ::  acqui_simu

  integer                                                       :: irec,isrc,ier
  integer                                                       :: icomp, it, irec_local
  real(kind=CUSTOM_REAL)                                        :: DT_cr, lw_tap
  double precision                                              :: DT
  character(len=256)                                            :: name_file
  character(len=MAX_STRING_LEN)                                 :: TRAC_PATH, dsname
  integer(kind=8)                                               :: filesize
  real(kind=CUSTOM_REAL), dimension(:), allocatable             :: raw_stf, filt_stf
  character(len=MAX_LEN_STRING)                                 :: name_file_tmp, ch_to_add

  if (myrank == 0 .and. DEBUG_MODE) write(INVERSE_LOG_FILE,*) ' initialize event number  : ', ievent

  ! time discretization
  NSTEP = acqui_simu(ievent)%Nt_data
  DT_cr = acqui_simu(ievent)%dt_data
  DT = DT_cr
  NSTEP_STF = 1

  deltat = real(DT,kind=CUSTOM_REAL)
  deltatover2 = deltat/2._CUSTOM_REAL
  deltatsqover2 = deltat*deltat/2._CUSTOM_REAL

  b_deltat = - real(DT,kind=CUSTOM_REAL)
  b_deltatover2 = b_deltat/2._CUSTOM_REAL
  b_deltatsqover2 = b_deltat*b_deltat/2._CUSTOM_REAL

  ! prepare source (only one source allowed for now)
  NSOURCES = acqui_simu(ievent)%nsources_tot  !! VM VM replace nsource_loc by nsource_tot since all arrays below
  if (USE_EXTERNAL_SOURCE_FILE) then          !! are using nsource_tot
    NSOURCES_STF = NSOURCES
    NSTEP_STF    = NSTEP
  else
    NSOURCES_STF = 1
    NSTEP_STF    = 1
  endif

  select case (acqui_simu(ievent)%source_type)

  case ('moment','force','shot')
     nsources_local =  acqui_simu(ievent)%nsources_local
     if (allocated(sourcearrays)) deallocate(sourcearrays)
     if (allocated(islice_selected_source)) deallocate(islice_selected_source)
     if (allocated(ispec_selected_source)) deallocate(ispec_selected_source)
     if (allocated(hdur)) deallocate(hdur)
     if (allocated(hdur_Gaussian)) deallocate(hdur_Gaussian)
     if (allocated(tshift_src)) deallocate(tshift_src)

     allocate(sourcearrays(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 490')
     allocate(islice_selected_source(NSOURCES),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 491')
     allocate(ispec_selected_source(NSOURCES),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 492')
     allocate(hdur(NSOURCES),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 493')
     allocate(hdur_Gaussian(NSOURCES),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 494')
     allocate(tshift_src(NSOURCES),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 495')
     sourcearrays(:,:,:,:,:)=acqui_simu(ievent)%sourcearrays(:,:,:,:,:)


     islice_selected_source(:)=acqui_simu(ievent)%islice_selected_source(:)
     ispec_selected_source(:)=acqui_simu(ievent)%ispec_selected_source(:)
     PRINT_SOURCE_TIME_FUNCTION=.true.
     t0 = acqui_simu(ievent)%t0
     hdur(:)=acqui_simu(ievent)%hdur(:)
     hdur_Gaussian(:)=acqui_simu(ievent)%hdur_Gaussian(:)
     tshift_src(:)=acqui_simu(ievent)%tshift(:)

     if (USE_EXTERNAL_SOURCE_FILE) then
        if (allocated(user_source_time_function)) deallocate(user_source_time_function)
        allocate(user_source_time_function(NSTEP_STF, NSOURCES_STF),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 496')
        if (ier /= 0) stop ' error in allocating user_source_time_function'
        if (inversion_param%only_forward) then
           user_source_time_function(:,:)=acqui_simu(ievent)%user_source_time_function(:,:)
        else
           !! filter the user stf
           !! EB EB Warning, filtering may be done each time we are switching events
           if (inversion_param%use_band_pass_filter) then
              allocate(raw_stf(NSTEP), filt_stf(NSTEP),stat=ier)
              if (ier /= 0) call exit_MPI_without_rank('error allocating array 497')
              do isrc=1,NSOURCES
                 raw_stf(:)=acqui_simu(ievent)%user_source_time_function(:,isrc)
                 call bwfilt (raw_stf, filt_stf, &
                              DT_cr, NSTEP, 1, 4, &
                              acqui_simu(ievent)%fl_event(inversion_param%current_ifrq), &
                              acqui_simu(ievent)%fh_event(inversion_param%current_ifrq))
                 lw_tap = 2.5_CUSTOM_REAL
                 call apodise_sig(filt_stf, NSTEP, lw_tap)
                 user_source_time_function(:,isrc)=filt_stf(:)
              enddo
              deallocate(raw_stf, filt_stf)
            else
               user_source_time_function(:,:)=acqui_simu(ievent)%user_source_time_function(:,:)
            endif
           !! write STF used to check
           if (VERBOSE_MODE .and. myrank == 0) then
               write(ch_to_add,'(a10,i4.4,a1,i4.4,a4)') '_stf_used_',ievent,'_',iter_inverse,'.txt'
               name_file_tmp = trim(acqui_simu(ievent)%data_file_gather)//trim(adjustl(ch_to_add))
               open(IINN,file=trim(adjustl(name_file_tmp)))
               do it=1, NSTEP
                  write(IINN,*) (it-1) * DT, user_source_time_function(it,1)
               enddo
               close(IINN)
           endif

        endif
     endif

     if (DEBUG_MODE) then
        write (IIDD , *)
        write (IIDD , *) 'first source islice , ispec : ', islice_selected_source(1), ispec_selected_source(1)
        write (IIDD , *)
     endif
     COUPLE_WITH_INJECTION_TECHNIQUE = .false.
  case('axisem')
     TRAC_PATH=acqui_simu(ievent)%traction_dir
     call create_name_database(dsname,myrank,TRAC_PATH)
     ! open traction file
     open(unit=IIN_veloc_dsm,file=dsname(1:len_trim(dsname))//'sol_axisem',status='old', &
          action='read',form='unformatted',iostat=ier)
     write(*,*) 'OPENING ', dsname(1:len_trim(dsname))//'sol_axisem'
     COUPLE_WITH_INJECTION_TECHNIQUE = .true.
     INJECTION_TECHNIQUE_TYPE = INJECTION_TECHNIQUE_IS_AXISEM
  case('fk')
     FKMODEL_FILE=acqui_simu(ievent)%source_file
     COUPLE_WITH_INJECTION_TECHNIQUE = .true.
     INJECTION_TECHNIQUE_TYPE = INJECTION_TECHNIQUE_IS_FK
  case default
     write(*,*) 'source ', acqui_simu(ievent)%source_type, 'Not yet '
     stop

  end select

  ! prepare receiver for the current event
  !! store current variables
  nrec = acqui_simu(ievent)%nsta_tot
  nrec_local = acqui_simu(ievent)%nsta_slice
  nadj_rec_local = 0                        !! by default dummy 0

  !! de-allocation
  if (allocated(islice_selected_rec)) deallocate(islice_selected_rec)
  if (allocated(ispec_selected_rec))  deallocate(ispec_selected_rec)
  if (allocated(xi_receiver)) deallocate(xi_receiver)
  if (allocated(eta_receiver)) deallocate(eta_receiver)
  if (allocated(gamma_receiver)) deallocate(gamma_receiver)
  if (allocated(station_name)) deallocate(station_name)
  if (allocated(network_name)) deallocate(network_name)
  if (allocated(nu_rec)) deallocate(nu_rec)

  !! re-allocation
  allocate(islice_selected_rec(nrec),ispec_selected_rec(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 498')
  allocate(xi_receiver(nrec),eta_receiver(nrec),gamma_receiver(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 499')
  allocate(station_name(nrec),network_name(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 500')
  allocate(nu_rec(NDIM,NDIM,nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 501')

  !! store current arrays
  islice_selected_rec(:) = acqui_simu(ievent)%islice_selected_rec(:)
  ispec_selected_rec(:) = acqui_simu(ievent)%ispec_selected_rec(:)
  xi_receiver(:) = acqui_simu(ievent)%xi_rec(:)
  eta_receiver(:) = acqui_simu(ievent)%eta_rec(:)
  gamma_receiver(:) = acqui_simu(ievent)%gamma_rec(:)
  station_name(:) = acqui_simu(ievent)%station_name
  network_name(:) = acqui_simu(ievent)%network_name
  nu_rec(:,:,:) = acqui_simu(ievent)%nu_rec(:,:,:)

  if (nrec_local > 0) then

     if (DEBUG_MODE) then
        write(IIDD,*)
        write(IIDD,*)  'init_specfem_for_one_run nrec_local :', nrec_local
        write(IIDD,*)
     endif

     if (allocated(hxir_store)) deallocate(hxir_store)
     if (allocated(hetar_store)) deallocate(hetar_store)
     if (allocated(hgammar_store)) deallocate(hgammar_store)
     if (allocated(hpxir_store)) deallocate(hpxir_store)
     if (allocated(hpetar_store)) deallocate(hpetar_store)
     if (allocated(hpgammar_store)) deallocate(hpgammar_store)
     if (allocated(number_receiver_global)) deallocate(number_receiver_global)

     nadj_rec_local = nrec_local ! assumes SIMULATION_TYPE == 3
     ! checks
     if (SIMULATION_TYPE == 2) stop 'Error invalid simulation type InitSpecfemForOneRun()'

     allocate(hxir_store(NGLLX,nrec_local),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 502')
     allocate(hetar_store(NGLLY,nrec_local),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 503')
     allocate(hgammar_store(NGLLZ,nrec_local),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 504')
     allocate(hpxir_store(NGLLX,nrec_local),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 505')
     allocate(hpetar_store(NGLLY,nrec_local),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 506')
     allocate(hpgammar_store(NGLLZ,nrec_local),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 507')


     allocate(number_receiver_global(nrec_local),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 508')
     number_receiver_global(:)=acqui_simu(ievent)%number_receiver_global(1:nrec_local)

     ! assumes SIMULATION_TYPE == 3
     nullify(number_adjsources_global)
     nullify(hxir_adjstore)
     nullify(hetar_adjstore)
     nullify(hgammar_adjstore)
     number_adjsources_global => number_receiver_global
     hxir_adjstore => hxir_store
     hetar_adjstore => hetar_store
     hgammar_adjstore => hgammar_store

     do irec=1, nrec_local
        hxir_store(:,irec) = acqui_simu(ievent)%hxi(:,irec)
        hetar_store(:,irec) = acqui_simu(ievent)%heta(:,irec)
        hgammar_store(:,irec) = acqui_simu(ievent)%hgamma(:,irec)
        hpxir_store(:,irec) = acqui_simu(ievent)%hpxi(:,irec)
        hpetar_store(:,irec) = acqui_simu(ievent)%hpeta(:,irec)
        hpgammar_store(:,irec) = acqui_simu(ievent)%hpgamma(:,irec)
     enddo
     hxir_adjstore(:,:) = hxir_store(:,:)
     hetar_adjstore(:,:) = hetar_store(:,:)
     hgammar_adjstore(:,:) = hgammar_store(:,:)

     if (allocated(seismograms_d)) deallocate(seismograms_d)
     if (allocated(seismograms_v)) deallocate(seismograms_v)
     if (allocated(seismograms_a)) deallocate(seismograms_a)
     if (allocated(seismograms_p)) deallocate(seismograms_p)

     ! allocate seismogram array
     allocate(seismograms_d(NDIM,nrec_local,nlength_seismogram),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 509')
     if (ier /= 0) stop 'error allocating array seismograms_d'
     allocate(seismograms_v(NDIM,nrec_local,nlength_seismogram),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 510')
     if (ier /= 0) stop 'error allocating array seismograms_v'
     allocate(seismograms_a(NDIM,nrec_local,nlength_seismogram),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 511')
     if (ier /= 0) stop 'error allocating array seismograms_a'
     allocate(seismograms_p(NDIM,nrec_local*NB_RUNS_ACOUSTIC_GPU,nlength_seismogram),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 512')
     if (ier /= 0) stop 'error allocating array seismograms_p'

     ! initialize seismograms
     seismograms_d(:,:,:) = 0._CUSTOM_REAL
     seismograms_v(:,:,:) = 0._CUSTOM_REAL
     seismograms_a(:,:,:) = 0._CUSTOM_REAL
     seismograms_p(:,:,:) = 0._CUSTOM_REAL

  else


     if (allocated(hxir_store)) deallocate(hxir_store)
     if (allocated(hetar_store)) deallocate(hetar_store)
     if (allocated(hgammar_store)) deallocate(hgammar_store)
     if (allocated(hpxir_store)) deallocate(hpxir_store)
     if (allocated(hpetar_store)) deallocate(hpetar_store)
     if (allocated(hpgammar_store)) deallocate(hpgammar_store)
     if (allocated(number_receiver_global)) deallocate(number_receiver_global)

     ! in Fortran it is legal to allocate dummy arrays with a size of zero
     allocate(hxir_store(0,0),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 513')
     allocate(hetar_store(0,0),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 514')
     allocate(hgammar_store(0,0),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 515')
     allocate(hpxir_store(0,0),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 516')
     allocate(hpetar_store(0,0),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 517')
     allocate(hpgammar_store(0,0),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 518')

     allocate(number_receiver_global(0),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 519')

     ! assumes SIMULATION_TYPE == 3
     nullify(number_adjsources_global)
     nullify(hxir_adjstore)
     nullify(hetar_adjstore)
     nullify(hgammar_adjstore)
     number_adjsources_global => number_receiver_global
     hxir_adjstore => hxir_store
     hetar_adjstore => hetar_store
     hgammar_adjstore => hgammar_store

     if (allocated(seismograms_d)) deallocate(seismograms_d)
     if (allocated(seismograms_v)) deallocate(seismograms_v)
     if (allocated(seismograms_a)) deallocate(seismograms_a)
     if (allocated(seismograms_p)) deallocate(seismograms_p)

     ! allocate seismogram array
     allocate(seismograms_d(NDIM,1,1),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 520')
     if (ier /= 0) stop 'error allocating array seismograms_d'
     allocate(seismograms_v(NDIM,1,1),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 521')
     if (ier /= 0) stop 'error allocating array seismograms_v'
     allocate(seismograms_a(NDIM,1,1),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 522')
     if (ier /= 0) stop 'error allocating array seismograms_a'
     allocate(seismograms_p(NDIM,1,1),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 523')
     if (ier /= 0) stop 'error allocating array seismograms_p'

  endif

  !! this is to skip writing seismogram on disk by specfem (both forward and ajoint)
  !! do not change
  NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP
  INVERSE_FWI_FULL_PROBLEM = .true.
  NTSTEP_BETWEEN_READ_ADJSRC = NSTEP

  ! initializes adjoint sources
  if (allocated(source_adjoint)) deallocate(source_adjoint)
  allocate(source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 524')
  if (ier /= 0) stop 'error allocating array adj_sourcearrays'
  source_adjoint(:,:,:) = 0._CUSTOM_REAL
  if (SIMULATION_TYPE == 3) then
     do icomp = 1,NDIM
        do it = 1,NTSTEP_BETWEEN_READ_ADJSRC
           do irec_local = 1, nadj_rec_local
              source_adjoint(icomp, irec_local, it) = acqui_simu(ievent)%adjoint_sources(icomp,irec_local,it)
           enddo
        enddo
     enddo
  endif

  ! manage the storage in memory and in disk the simulaed wavefields
  !! with GPU we have to be carrefull of warning for seismogram not all things are allowed and this make the code crashes
  if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION ) then

     !! todo recuperer inversion_paral%get_synthetics_**
     SAVE_SEISMOGRAMS_DISPLACEMENT =.true.
     SAVE_SEISMOGRAMS_VELOCITY     =.false.
     SAVE_SEISMOGRAMS_ACCELERATION =.false.
     SAVE_SEISMOGRAMS_PRESSURE     =.true.

  else

     if (ACOUSTIC_SIMULATION) then
        SAVE_SEISMOGRAMS_PRESSURE     =.true.
        SAVE_SEISMOGRAMS_DISPLACEMENT =.false.
        SAVE_SEISMOGRAMS_VELOCITY     =.false.
        SAVE_SEISMOGRAMS_ACCELERATION =.false.
     endif

     if (ELASTIC_SIMULATION) then
        SAVE_SEISMOGRAMS_PRESSURE      =.false.
        SAVE_SEISMOGRAMS_DISPLACEMENT  =.true.
        SAVE_SEISMOGRAMS_VELOCITY      =.false.
        SAVE_SEISMOGRAMS_ACCELERATION  =.false.
     endif

  endif

  !! clean arrays
  ! reset all forward wavefields
  call prepare_wavefields()  !! routine from specfem

  ! memory variables if attenuation
  if (ATTENUATION) then
      ! clear memory variables if attenuation
     ! initialize memory variables for attenuation
     epsilondev_trace(:,:,:,:) = 0._CUSTOM_REAL
     epsilondev_xx(:,:,:,:) = 0._CUSTOM_REAL
     epsilondev_yy(:,:,:,:) = 0._CUSTOM_REAL
     epsilondev_xy(:,:,:,:) = 0._CUSTOM_REAL
     epsilondev_xz(:,:,:,:) = 0._CUSTOM_REAL
     epsilondev_yz(:,:,:,:) = 0._CUSTOM_REAL

     R_trace(:,:,:,:,:) = 0._CUSTOM_REAL
     R_xx(:,:,:,:,:) = 0._CUSTOM_REAL
     R_yy(:,:,:,:,:) = 0._CUSTOM_REAL
     R_xy(:,:,:,:,:) = 0._CUSTOM_REAL
     R_xz(:,:,:,:,:) = 0._CUSTOM_REAL
     R_yz(:,:,:,:,:) = 0._CUSTOM_REAL

     if (FIX_UNDERFLOW_PROBLEM) then
        R_trace(:,:,:,:,:) = VERYSMALLVAL
        R_xx(:,:,:,:,:) = VERYSMALLVAL
        R_yy(:,:,:,:,:) = VERYSMALLVAL
        R_xy(:,:,:,:,:) = VERYSMALLVAL
        R_xz(:,:,:,:,:) = VERYSMALLVAL
        R_yz(:,:,:,:,:) = VERYSMALLVAL
     endif
  endif

  ! reaset all adjoint wavefield
  ! elastic domain
  if (ELASTIC_SIMULATION) then
     ! from prepare_timerun_lddrk()
     if (SIMULATION_TYPE == 3) then
        b_R_xx_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_yy_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_xy_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_xz_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_yz_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_trace_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        if (FIX_UNDERFLOW_PROBLEM) then
           b_R_xx_lddrk(:,:,:,:,:) = VERYSMALLVAL
           b_R_yy_lddrk(:,:,:,:,:) = VERYSMALLVAL
           b_R_xy_lddrk(:,:,:,:,:) = VERYSMALLVAL
           b_R_xz_lddrk(:,:,:,:,:) = VERYSMALLVAL
           b_R_yz_lddrk(:,:,:,:,:) = VERYSMALLVAL
           b_R_trace_lddrk(:,:,:,:,:) = VERYSMALLVAL
        endif
     endif

     ! reconstructed/backward elastic wavefields
     b_displ = 0._CUSTOM_REAL
     b_veloc = 0._CUSTOM_REAL
     b_accel = 0._CUSTOM_REAL
     if (FIX_UNDERFLOW_PROBLEM) b_displ = VERYSMALLVAL

     ! memory variables if attenuation
     if (ATTENUATION) then
        b_R_trace = 0._CUSTOM_REAL
        b_R_xx = 0._CUSTOM_REAL
        b_R_yy = 0._CUSTOM_REAL
        b_R_xy = 0._CUSTOM_REAL
        b_R_xz = 0._CUSTOM_REAL
        b_R_yz = 0._CUSTOM_REAL
        b_epsilondev_trace = 0._CUSTOM_REAL
        b_epsilondev_xx = 0._CUSTOM_REAL
        b_epsilondev_yy = 0._CUSTOM_REAL
        b_epsilondev_xy = 0._CUSTOM_REAL
        b_epsilondev_xz = 0._CUSTOM_REAL
        b_epsilondev_yz = 0._CUSTOM_REAL

        if (FIX_UNDERFLOW_PROBLEM) then
           b_R_trace(:,:,:,:,:) = VERYSMALLVAL
           b_R_xx(:,:,:,:,:) = VERYSMALLVAL
           b_R_yy(:,:,:,:,:) = VERYSMALLVAL
           b_R_xy(:,:,:,:,:) = VERYSMALLVAL
           b_R_xz(:,:,:,:,:) = VERYSMALLVAL
           b_R_yz(:,:,:,:,:) = VERYSMALLVAL
        endif

     endif
  endif

  ! acoustic domain
  if (ACOUSTIC_SIMULATION) then
     ! reconstructed/backward acoustic potentials
     b_potential_acoustic = 0._CUSTOM_REAL
     b_potential_dot_acoustic = 0._CUSTOM_REAL
     b_potential_dot_dot_acoustic = 0._CUSTOM_REAL
     if (FIX_UNDERFLOW_PROBLEM) b_potential_acoustic = VERYSMALLVAL
  endif

  ! poroelastic domain
  if (POROELASTIC_SIMULATION) then
     b_displs_poroelastic = 0._CUSTOM_REAL
     b_velocs_poroelastic = 0._CUSTOM_REAL
     b_accels_poroelastic = 0._CUSTOM_REAL
     b_displw_poroelastic = 0._CUSTOM_REAL
     b_velocw_poroelastic = 0._CUSTOM_REAL
     b_accelw_poroelastic = 0._CUSTOM_REAL
     if (FIX_UNDERFLOW_PROBLEM) b_displs_poroelastic = VERYSMALLVAL
     if (FIX_UNDERFLOW_PROBLEM) b_displw_poroelastic = VERYSMALLVAL
  endif


  !! open boundary files
  if (STACEY_ABSORBING_CONDITIONS) then

     if (SIMULATION_TYPE == 3 ) then  !! read only open

        ! opens existing files
        if (ELASTIC_SIMULATION) then

           ! size of single record
           b_reclen_field = CUSTOM_REAL * NDIM * NGLLSQUARE * num_abs_boundary_faces
           ! check integer size limit: size of b_reclen_field must fit onto an 4-byte integer
           if (num_abs_boundary_faces > int(2147483646.0 / (CUSTOM_REAL * NDIM * NGLLSQUARE))) then
              print *,'reclen needed exceeds integer 4-byte limit: ',b_reclen_field
              print *,'  ',CUSTOM_REAL, NDIM, NGLLSQUARE, num_abs_boundary_faces
              print *,'bit size Fortran: ',bit_size(b_reclen_field)
              call exit_MPI(myrank,"error b_reclen_field integer limit")
           endif
           ! total file size
           filesize = b_reclen_field
           filesize = filesize * NSTEP
           call open_file_abs_r(IOABS,trim(prname)//'absorb_field.bin', &
                len_trim(trim(prname)//'absorb_field.bin'), filesize)
        endif

        if (ACOUSTIC_SIMULATION) then

           ! size of single record
           b_reclen_potential = CUSTOM_REAL * NGLLSQUARE * num_abs_boundary_faces
           ! check integer size limit: size of b_reclen_potential must fit onto an 4-byte integer
           if (num_abs_boundary_faces > int(2147483646.0 / (CUSTOM_REAL * NGLLSQUARE))) then
              print *,'reclen needed exceeds integer 4-byte limit: ',b_reclen_potential
              print *,'  ',CUSTOM_REAL, NGLLSQUARE, num_abs_boundary_faces
              print *,'bit size Fortran: ',bit_size(b_reclen_potential)
              call exit_MPI(myrank,"error b_reclen_potential integer limit")
           endif
           ! total file size (two lines to implicitly convert to 8-byte integers)
           filesize = b_reclen_potential
           filesize = filesize*NSTEP
           call open_file_abs_r(IOABS_AC,trim(prname)//'absorb_potential.bin', &
                len_trim(trim(prname)//'absorb_potential.bin'), filesize)
        endif

     else  !!! write open

        if (SAVE_FORWARD) then

           ! opens new file
           if (ELASTIC_SIMULATION) then
              ! size of single record
              b_reclen_field = CUSTOM_REAL * NDIM * NGLLSQUARE * num_abs_boundary_faces
              ! check integer size limit: size of b_reclen_field must fit onto an 4-byte integer
              if (num_abs_boundary_faces > int(2147483646.0 / (CUSTOM_REAL * NDIM * NGLLSQUARE))) then
                 print *,'reclen needed exceeds integer 4-byte limit: ',b_reclen_field
                 print *,'  ',CUSTOM_REAL, NDIM, NGLLSQUARE, num_abs_boundary_faces
                 print *,'bit size Fortran: ',bit_size(b_reclen_field)
                 call exit_MPI(myrank,"error b_reclen_field integer limit")
              endif
              ! total file size
              filesize = b_reclen_field
              filesize = filesize * NSTEP
              call open_file_abs_w(IOABS,trim(prname)//'absorb_field.bin', &
                   len_trim(trim(prname)//'absorb_field.bin'), filesize)
           endif

           if (ACOUSTIC_SIMULATION) then
              ! size of single record
              b_reclen_potential = CUSTOM_REAL * NGLLSQUARE * num_abs_boundary_faces
              ! check integer size limit: size of b_reclen_potential must fit onto an 4-byte integer
              if (num_abs_boundary_faces > int(2147483646.0 / (CUSTOM_REAL * NGLLSQUARE))) then
                 print *,'reclen needed exceeds integer 4-byte limit: ',b_reclen_potential
                 print *,'  ',CUSTOM_REAL, NGLLSQUARE, num_abs_boundary_faces
                 print *,'bit size Fortran: ',bit_size(b_reclen_potential)
                 call exit_MPI(myrank,"error b_reclen_potential integer limit")
              endif
              ! total file size (two lines to implicitly convert to 8-byte integers)
              filesize = b_reclen_potential
              filesize = filesize*NSTEP
              call open_file_abs_w(IOABS_AC,trim(prname)//'absorb_potential.bin', &
                   len_trim(trim(prname)//'absorb_potential.bin'), filesize)
           endif

        endif
     endif
  endif

  !! reallocate all GPU memory according the Fortran arrays
  if (GPU_MODE) call prepare_GPU()

  !! open new log file for specfem
  if (myrank == 0 .and. SIMULATION_TYPE == 1) then
     close(IMAIN)
     write(name_file,'(a15,i5.5,a1,i5.5,a4)') '/output_solver_',ievent,'_',iter_inverse,'.txt'
     open(unit=IMAIN,file=trim(OUTPUT_FILES)//trim(name_file),status='unknown')
  endif

  !! info on mesh and parameters
  if (SIMULATION_TYPE == 1) then
    call check_mesh_resolution(NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                               ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
                               kappastore,mustore,rhostore, &
                               phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                               DT, model_speed_max,min_resolved_period)
  endif

  end subroutine InitSpecfemForOneRun


!---------------------------------------------------------------------------------------------
!> finalize, specfem after running each event
!---------------------------------------------------------------------------------------------

  subroutine  FinalizeSpecfemForOneRun(acqui_simu, ievent)

  implicit none

  integer,                                        intent(in)    :: ievent
  type(acqui),  dimension(:), allocatable,        intent(in)    :: acqui_simu
  integer                                                       :: ier


  !! manage files due to coupling with axisem
  select case (acqui_simu(ievent)%source_type)

  case('axisem')
     close(IIN_veloc_dsm)

  case('fk', 'moment', 'force', 'shot')
     !! nothing to do for the moment

  case default   !! for stopping if any problem with source_type
     write(*,*) 'source ', acqui_simu(ievent)%source_type, 'Not yet '
     stop

  end select

  !! finalize specfem run
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
     open(unit=IOUT,file=prname(1:len_trim(prname))//'save_forward_arrays.bin', &
          status='unknown',form='unformatted',iostat=ier)
     if (ier /= 0) then
        print *,'error: opening save_forward_arrays.bin'
        print *,'path: ',prname(1:len_trim(prname))//'save_forward_arrays.bin'
        call exit_mpi(myrank,'error opening file save_forward_arrays.bin')
     endif

     if (ACOUSTIC_SIMULATION) then
        write(IOUT) potential_acoustic
        write(IOUT) potential_dot_acoustic
        write(IOUT) potential_dot_dot_acoustic
     endif

     if (ELASTIC_SIMULATION) then
        write(IOUT) displ
        write(IOUT) veloc
        write(IOUT) accel

        if (ATTENUATION) then
           write(IOUT) R_trace
           write(IOUT) R_xx
           write(IOUT) R_yy
           write(IOUT) R_xy
           write(IOUT) R_xz
           write(IOUT) R_yz
           write(IOUT) epsilondev_trace
           write(IOUT) epsilondev_xx
           write(IOUT) epsilondev_yy
           write(IOUT) epsilondev_xy
           write(IOUT) epsilondev_xz
           write(IOUT) epsilondev_yz
        endif
     endif

     if (POROELASTIC_SIMULATION) then
        write(IOUT) displs_poroelastic
        write(IOUT) velocs_poroelastic
        write(IOUT) accels_poroelastic
        write(IOUT) displw_poroelastic
        write(IOUT) velocw_poroelastic
        write(IOUT) accelw_poroelastic
     endif

     close(IOUT)
  endif

  if (ELASTIC_SIMULATION .and. (SIMULATION_TYPE == 3 .or. SAVE_FORWARD) .and. STACEY_ABSORBING_CONDITIONS) then
     call close_file_abs(IOABS)
  endif

  if (ACOUSTIC_SIMULATION .and. (SIMULATION_TYPE == 3 .or. SAVE_FORWARD) .and. STACEY_ABSORBING_CONDITIONS) then
     call close_file_abs(IOABS_AC)
  endif

  if (GPU_MODE .and. SIMULATION_TYPE == 3) then

     if (ACOUSTIC_SIMULATION) then
        rho_ac_kl_GPU(:,:,:,:)=rho_ac_kl_GPU(:,:,:,:)+rho_ac_kl(:,:,:,:)
        kappa_ac_kl_GPU(:,:,:,:)=kappa_ac_kl_GPU(:,:,:,:)+kappa_ac_kl(:,:,:,:)
        if (APPROXIMATE_HESS_KL) then
           hess_rho_ac_kl_GPU(:,:,:,:)=hess_rho_ac_kl_GPU(:,:,:,:)+hess_rho_ac_kl(:,:,:,:)
           hess_kappa_ac_kl_GPU(:,:,:,:)=hess_kappa_ac_kl_GPU(:,:,:,:)+hess_kappa_ac_kl(:,:,:,:)
        endif
     endif

     if (ELASTIC_SIMULATION) then
        rho_kl_GPU(:,:,:,:)=rho_kl_GPU(:,:,:,:)+rho_kl(:,:,:,:)
        if (ANISOTROPIC_KL) then
           cijkl_kl_GPU(:,:,:,:,:)=cijkl_kl_GPU(:,:,:,:,:)+cijkl_kl(:,:,:,:,:)
        else
           kappa_kl_GPU(:,:,:,:)=kappa_kl_GPU(:,:,:,:)+kappa_kl(:,:,:,:)
           mu_kl_GPU(:,:,:,:)= mu_kl_GPU(:,:,:,:)+ mu_kl(:,:,:,:)
           if (APPROXIMATE_HESS_KL) then
              hess_rho_kl_GPU(:,:,:,:)=hess_rho_kl_GPU(:,:,:,:)+hess_rho_kl(:,:,:,:)
              hess_kappa_kl_GPU(:,:,:,:)=hess_kappa_kl_GPU(:,:,:,:)+hess_kappa_kl(:,:,:,:)
              hess_mu_kl_GPU(:,:,:,:)= hess_mu_kl_GPU(:,:,:,:)+ hess_mu_kl(:,:,:,:)
           endif
        endif

     endif

  endif

  end subroutine FinalizeSpecfemForOneRun


!------------------------------------------------------------------------------
! Initialize for one step FWI
!------------------------------------------------------------------------------

  subroutine InitForOneStepFWI(inversion_param)

  implicit none

  type(inver),                                    intent(inout) :: inversion_param

  !! reset cost function
  inversion_param%total_current_cost = 0._CUSTOM_REAL
  inversion_param%data_std = 0._CUSTOM_REAL
  inversion_param%nb_data_std = 0._CUSTOM_REAL
  !! reset kenels  ----------------------------

  ! acoustic domain
  if (ACOUSTIC_SIMULATION) then
     rho_ac_kl(:,:,:,:)=0._CUSTOM_REAL
     kappa_ac_kl(:,:,:,:)=0._CUSTOM_REAL

     if (GPU_MODE) then
        rho_ac_kl_GPU(:,:,:,:)=0._CUSTOM_REAL
        kappa_ac_kl_GPU(:,:,:,:)=0._CUSTOM_REAL
     endif

  endif

  ! elastic domain
  if (ELASTIC_SIMULATION) then

     rho_kl(:,:,:,:)   = 0._CUSTOM_REAL
     if (GPU_MODE)  rho_kl_GPU(:,:,:,:)=0._CUSTOM_REAL

     if (ANISOTROPIC_KL) then
        cijkl_kl(:,:,:,:,:) = 0._CUSTOM_REAL
        if (GPU_MODE) cijkl_kl_GPU(:,:,:,:,:)=0._CUSTOM_REAL
     else
        mu_kl(:,:,:,:)    = 0._CUSTOM_REAL
        kappa_kl(:,:,:,:) = 0._CUSTOM_REAL

        if (GPU_MODE) then
           kappa_kl_GPU(:,:,:,:)=0._CUSTOM_REAL
           mu_kl_GPU(:,:,:,:)=0._CUSTOM_REAL
        endif

     endif

     if (APPROXIMATE_HESS_KL) then
        if (ELASTIC_SIMULATION) then

           hess_kl(:,:,:,:)   = 0._CUSTOM_REAL
           hess_kappa_kl(:,:,:,:)   = 0._CUSTOM_REAL
           hess_mu_kl(:,:,:,:)   = 0._CUSTOM_REAL
           hess_rho_kl(:,:,:,:)   = 0._CUSTOM_REAL

           if (GPU_MODE) then
              hess_kappa_kl_GPU(:,:,:,:)   = 0._CUSTOM_REAL
              hess_mu_kl_GPU(:,:,:,:)   = 0._CUSTOM_REAL
              hess_rho_kl_GPU(:,:,:,:)   = 0._CUSTOM_REAL
           endif

        endif
        if (ACOUSTIC_SIMULATION) then

           hess_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL
           hess_kappa_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL
           hess_rho_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL

           if (GPU_MODE) then
               hess_rho_ac_kl_GPU(:,:,:,:)   = 0._CUSTOM_REAL
               hess_kappa_ac_kl_GPU(:,:,:,:)   = 0._CUSTOM_REAL
           endif

        endif
     endif

  endif

  end subroutine InitForOneStepFWI


!-------------------------------------------------------------------
!> General initialization for specfem in order to do FWI.  Directly some specfem subroutines are called to initialize modeling
!-------------------------------------------------------------------

  subroutine InitializeSpecfemForInversion()

  use specfem_par
  use shared_parameters
  use specfem_par_elastic
  !use time_iteration_mod (not yet used)

  implicit none


  !! following subroutines are directly from specfem3D git devel version
  call initialize_simulation()     !! here : we need to initialize with NGLOB_ADJ=NSPEC_ADJ=NSPEC_STRAIN_ONLY=1

  !! creating dummy inputs : STATION, CMTSOLUTION, STATION_ADJOINT and SEM
  !! to be able to run the specfem subroutines that initialize solver.
  !! In anyway, we will use correct parameters for FWI by initializing before each call to specfem solver
  if (myrank == 0) call CreateInitDummyFiles()

  !! thus store right values hat are not been correctly initialized :
  NGLOB_ADJOINT = NGLOB_AB
  NSPEC_ADJOINT = NSPEC_AB
  NSPEC_STRAIN_ONLY = NSPEC_AB

  !! enforce to allocate all arrays for both adjoint and direct simulation
  SIMULATION_TYPE = 3
  APPROXIMATE_HESS_KL = .true. !! test preconditionner
  PRINT_SOURCE_TIME_FUNCTION = .true.

  call read_mesh_databases()
  call read_mesh_databases_moho()
  call read_mesh_databases_adjoint()

  ! safety check
  if (NSPEC_IRREGULAR /= NSPEC_AB) stop 'Please check inverse problem routine for NSPEC_AB /= NSPEC_IRREGULAR'

  call setup_GLL_points()
  call detect_mesh_surfaces()

  call setup_sources_receivers()  !! we have one dummy source and STATION_ADJOINT to set up without crashes

  SIMULATION_TYPE = 1               !! here we need to prepare the fisrt run
  SAVE_FORWARD = .true.             !! which is mandatory direct and need to save forward wavefield

  call prepare_timerun()          !! absord boundary are opened here

  !! we need to clean the GPU memory because we will change arrays
  if (GPU_MODE)  call prepare_cleanup_device(Mesh_pointer,ACOUSTIC_SIMULATION,ELASTIC_SIMULATION, &
                                             STACEY_ABSORBING_CONDITIONS,NOISE_TOMOGRAPHY,COMPUTE_AND_STORE_STRAIN, &
                                             ATTENUATION,ANISOTROPY,APPROXIMATE_OCEAN_LOAD, &
                                             APPROXIMATE_HESS_KL)

  !!--------------------------------------------------------------------
  !! this is specific prepare_timerun version for FWI (subroutine defined below)
  call PrepareTimerunInverseProblem()          !! absord boundary are closed here
  !! TODO
  !!if (USE_UNDO_ATT) call prepare_time_iteration() !! allocate arrays for store saving snapshots and displacement
  !---------------------------------------------------------------------

  end subroutine InitializeSpecfemForInversion


!---------------------------------------------------------------------------------
!> specfific prepare_timerun for FWI (it is actually missing initialization
!! done in prepare_timerun when SIMULATION_TYPE=1. Here we add all
!! actions performed for SIMULATION_TYPE=3. (We cannot call directly
!! use prepare_timerun with  SIMULATION_TYPE=3 at first time.
!---------------------------------------------------------------------------------

  subroutine PrepareTimerunInverseProblem()

  use specfem_par
  use shared_parameters
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  integer  :: ier

  ! from initialize_simulation_adjoint
  NSPEC_ADJOINT = NSPEC_AB
  NGLOB_ADJOINT = NGLOB_AB

  ! from prepare_timerun_constants()
  b_deltat = - real(DT,kind=CUSTOM_REAL)
  b_deltatover2 = b_deltat/2._CUSTOM_REAL
  b_deltatsqover2 = b_deltat*b_deltat/2._CUSTOM_REAL

  ! from prepare_timerun_lddrk()
  if (ELASTIC_SIMULATION) then
     allocate(b_R_xx_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 527')
     allocate(b_R_yy_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 528')
     allocate(b_R_xy_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 529')
     allocate(b_R_xz_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 530')
     allocate(b_R_yz_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 531')
     if (ier /= 0) stop 'Error allocating array R_**_lddrk etc.'

     allocate(b_R_trace_lddrk(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 532')
     if (ier /= 0) stop 'Error allocating array R_**_lddrk etc.'

     if (SIMULATION_TYPE == 3) then
        b_R_xx_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_yy_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_xy_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_xz_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_yz_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        b_R_trace_lddrk(:,:,:,:,:) = 0._CUSTOM_REAL
        if (FIX_UNDERFLOW_PROBLEM) then
           b_R_xx_lddrk(:,:,:,:,:) = VERYSMALLVAL
           b_R_yy_lddrk(:,:,:,:,:) = VERYSMALLVAL
           b_R_xy_lddrk(:,:,:,:,:) = VERYSMALLVAL
           b_R_xz_lddrk(:,:,:,:,:) = VERYSMALLVAL
           b_R_yz_lddrk(:,:,:,:,:) = VERYSMALLVAL
           b_R_trace_lddrk(:,:,:,:,:) = VERYSMALLVAL
        endif
     endif
  endif

  !! from prepare_timerun_adjoint()
  ! attenuation backward memories
  if (ATTENUATION .and. SIMULATION_TYPE == 3) then
     ! precompute Runge-Kutta coefficients if attenuation
     call get_attenuation_memory_values(tau_sigma,b_deltat,b_alphaval,b_betaval,b_gammaval)
  endif

  ! elastic domain
  if (ELASTIC_SIMULATION) then
     rho_kl(:,:,:,:)   = 0._CUSTOM_REAL

     if (ANISOTROPIC_KL) then
        cijkl_kl(:,:,:,:,:) = 0._CUSTOM_REAL
     else
        mu_kl(:,:,:,:)    = 0._CUSTOM_REAL
        kappa_kl(:,:,:,:) = 0._CUSTOM_REAL
     endif

     if (APPROXIMATE_HESS_KL) then
        hess_kl(:,:,:,:)   = 0._CUSTOM_REAL
        !! VM VM SHIN
        !Shin_hess_kappa_kl(:,:,:,:)   = 0._CUSTOM_REAL
        !Shin_hess_mu_kl(:,:,:,:)   = 0._CUSTOM_REAL
        !Shin_hess_rho_kl(:,:,:,:)   = 0._CUSTOM_REAL
     endif

     ! reconstructed/backward elastic wavefields
     b_displ = 0._CUSTOM_REAL
     b_veloc = 0._CUSTOM_REAL
     b_accel = 0._CUSTOM_REAL
     if (FIX_UNDERFLOW_PROBLEM) b_displ = VERYSMALLVAL

     ! memory variables if attenuation
     if (ATTENUATION) then
        b_R_trace = 0._CUSTOM_REAL
        b_R_xx = 0._CUSTOM_REAL
        b_R_yy = 0._CUSTOM_REAL
        b_R_xy = 0._CUSTOM_REAL
        b_R_xz = 0._CUSTOM_REAL
        b_R_yz = 0._CUSTOM_REAL
        b_epsilondev_trace = 0._CUSTOM_REAL
        b_epsilondev_xx = 0._CUSTOM_REAL
        b_epsilondev_yy = 0._CUSTOM_REAL
        b_epsilondev_xy = 0._CUSTOM_REAL
        b_epsilondev_xz = 0._CUSTOM_REAL
        b_epsilondev_yz = 0._CUSTOM_REAL
     endif

     ! moho kernels
     if (SAVE_MOHO_MESH) moho_kl(:,:) = 0._CUSTOM_REAL
  endif

  ! acoustic domain
  if (ACOUSTIC_SIMULATION) then
     rho_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL
     kappa_ac_kl(:,:,:,:) = 0._CUSTOM_REAL

     if (APPROXIMATE_HESS_KL) then
        hess_ac_kl(:,:,:,:)   = 0._CUSTOM_REAL
        !! VM VM SHIN
        !Shin_hess_kappa_kl(:,:,:,:)   = 0._CUSTOM_REAL
        !Shin_hess_mu_kl(:,:,:,:)   = 0._CUSTOM_REAL
        !Shin_hess_rho_kl(:,:,:,:)   = 0._CUSTOM_REAL
     endif

     ! reconstructed/backward acoustic potentials
     b_potential_acoustic = 0._CUSTOM_REAL
     b_potential_dot_acoustic = 0._CUSTOM_REAL
     b_potential_dot_dot_acoustic = 0._CUSTOM_REAL
     if (FIX_UNDERFLOW_PROBLEM) b_potential_acoustic = VERYSMALLVAL

  endif

  ! poroelastic domain
  if (POROELASTIC_SIMULATION) then
     rhot_kl(:,:,:,:)   = 0._CUSTOM_REAL
     rhof_kl(:,:,:,:)   = 0._CUSTOM_REAL
     sm_kl(:,:,:,:)   = 0._CUSTOM_REAL
     eta_kl(:,:,:,:)   = 0._CUSTOM_REAL
     mufr_kl(:,:,:,:)    = 0._CUSTOM_REAL
     B_kl(:,:,:,:) = 0._CUSTOM_REAL
     C_kl(:,:,:,:) = 0._CUSTOM_REAL
     M_kl(:,:,:,:) = 0._CUSTOM_REAL

     !if (APPROXIMATE_HESS_KL) &
     !  hess_kl(:,:,:,:)   = 0._CUSTOM_REAL

     ! reconstructed/backward elastic wavefields
     b_displs_poroelastic = 0._CUSTOM_REAL
     b_velocs_poroelastic = 0._CUSTOM_REAL
     b_accels_poroelastic = 0._CUSTOM_REAL
     b_displw_poroelastic = 0._CUSTOM_REAL
     b_velocw_poroelastic = 0._CUSTOM_REAL
     b_accelw_poroelastic = 0._CUSTOM_REAL
     if (FIX_UNDERFLOW_PROBLEM) b_displs_poroelastic = VERYSMALLVAL
     if (FIX_UNDERFLOW_PROBLEM) b_displw_poroelastic = VERYSMALLVAL
  endif

  ! close absorb files : close here because we will open at each new specfem run
  if (ELASTIC_SIMULATION .and. (SIMULATION_TYPE == 3 .or. SAVE_FORWARD) .and. STACEY_ABSORBING_CONDITIONS) &
                      call close_file_abs(IOABS)
  if (ACOUSTIC_SIMULATION .and. (SIMULATION_TYPE == 3 .or. SAVE_FORWARD) .and. STACEY_ABSORBING_CONDITIONS) &
                      call close_file_abs(IOABS_AC)

  !! allocate arrays for saving the kernel computed by GPU in CPU memory in order to perform summation over events.
  if (GPU_MODE) then
     if (ACOUSTIC_SIMULATION) then
        allocate(rho_ac_kl_GPU(NGLLX,NGLLY,NGLLZ,NSPEC_AB), kappa_ac_kl_GPU(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 533')
        rho_ac_kl_GPU(:,:,:,:)=0._CUSTOM_REAL
        kappa_ac_kl_GPU(:,:,:,:)=0._CUSTOM_REAL
        if (APPROXIMATE_HESS_KL) then
           allocate(hess_rho_ac_kl_GPU(NGLLX,NGLLY,NGLLZ,NSPEC_AB),hess_kappa_ac_kl_GPU(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
           if (ier /= 0) call exit_MPI_without_rank('error allocating array 534')
        endif
     endif

     if (ELASTIC_SIMULATION) then
        allocate(rho_kl_GPU(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 535')
        allocate(kappa_kl_GPU(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 536')
        allocate(mu_kl_GPU(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 537')
        rho_kl_GPU(:,:,:,:)=0._CUSTOM_REAL
        kappa_kl_GPU(:,:,:,:)=0._CUSTOM_REAL
        mu_kl_GPU(:,:,:,:)=0._CUSTOM_REAL
        if (APPROXIMATE_HESS_KL) then
           allocate(hess_rho_kl_GPU(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
           if (ier /= 0) call exit_MPI_without_rank('error allocating array 538')
           allocate(hess_kappa_kl_GPU(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
           if (ier /= 0) call exit_MPI_without_rank('error allocating array 539')
           allocate(hess_mu_kl_GPU(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
           if (ier /= 0) call exit_MPI_without_rank('error allocating array 540')
        endif
     endif

     if (ANISOTROPIC_KL) then
        allocate(cijkl_kl_GPU(21,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 541')
        cijkl_kl_GPU(:,:,:,:,:)=0._CUSTOM_REAL
     endif
  endif

  end subroutine PrepareTimerunInverseProblem


!----------------------------------------------------------------------
!>  copy summed GPU kernel overs events in specfem CPU arrays
!!
!!
!----------------------------------------------------------------------

  subroutine TransfertKernelFromGPUArrays()

  implicit none

  if (ACOUSTIC_SIMULATION) then
     rho_ac_kl(:,:,:,:)=rho_ac_kl_GPU(:,:,:,:)
     kappa_ac_kl(:,:,:,:)=kappa_ac_kl_GPU(:,:,:,:)
     if (APPROXIMATE_HESS_KL) then
        hess_rho_ac_kl(:,:,:,:)=hess_rho_ac_kl_GPU(:,:,:,:)
        hess_kappa_ac_kl(:,:,:,:)=hess_kappa_ac_kl_GPU(:,:,:,:)
     endif
  endif

  if (ELASTIC_SIMULATION) then

     rho_kl(:,:,:,:)=rho_kl_GPU(:,:,:,:)

     if (ANISOTROPIC_KL) then
        cijkl_kl(:,:,:,:,:)=cijkl_kl_GPU(:,:,:,:,:)
     else
        kappa_kl(:,:,:,:)=kappa_kl_GPU(:,:,:,:)
        mu_kl(:,:,:,:)=mu_kl_GPU(:,:,:,:)
        if (APPROXIMATE_HESS_KL) then
           hess_rho_kl(:,:,:,:)=hess_rho_kl_GPU(:,:,:,:)
           hess_kappa_kl(:,:,:,:)=hess_kappa_kl_GPU(:,:,:,:)
           hess_mu_kl(:,:,:,:)=hess_mu_kl_GPU(:,:,:,:)
        endif
     endif

  endif

  end subroutine TransfertKernelFromGPUArrays


!-------------------------------------------------------------------
! test if we can saftely perform a simulation in new model
!-------------------------------------------------------------------

  subroutine  CheckModelSuitabilityForModeling(ModelIsSuitable)

  implicit none

  logical, intent(in out) :: ModelIsSuitable

  ! local parameters
  real(kind=CUSTOM_REAL) :: vpmin,vpmax,vsmin,vsmax,vpmin_glob,vpmax_glob,vsmin_glob,vsmax_glob
  real(kind=CUSTOM_REAL) :: poissonmin,poissonmax,poissonmin_glob,poissonmax_glob
  real(kind=CUSTOM_REAL) :: distance_min,distance_max,distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max,elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min,x_max,x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min,y_max,y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min,z_max,z_min_glob,z_max_glob
  real(kind=CUSTOM_REAL) :: cmax,cmax_glob,pmax,pmax_glob
  real(kind=CUSTOM_REAL) :: dt_suggested,dt_suggested_glob,avg_distance
  real(kind=CUSTOM_REAL) :: vel_min,vel_max
  integer                :: ispec
  logical                :: has_vs_zero

  ModelIsSuitable = .true.

  vpmin_glob = HUGEVAL
  vpmax_glob = -HUGEVAL

  vsmin_glob = HUGEVAL
  vsmax_glob = -HUGEVAL

  poissonmin_glob = HUGEVAL
  poissonmax_glob = -HUGEVAL

  distance_min_glob = HUGEVAL
  distance_max_glob = -HUGEVAL

  x_min_glob = HUGEVAL
  x_max_glob = -HUGEVAL

  y_min_glob = HUGEVAL
  y_max_glob = -HUGEVAL

  z_min_glob = HUGEVAL
  z_max_glob = -HUGEVAL

  elemsize_min_glob = HUGEVAL
  elemsize_max_glob = -HUGEVAL

  cmax_glob = -HUGEVAL
  pmax_glob = -HUGEVAL

  dt_suggested_glob = HUGEVAL

  has_vs_zero = .false.

  do ispec = 1,NSPEC_AB
     ! determines minimum/maximum velocities within this element
     if (ispec_is_elastic(ispec) .or. ispec_is_acoustic(ispec)) then
        ! elastic/acoustic
        !! min max in element
        call  get_vpvs_minmax(vpmin, vpmax, vsmin, vsmax, poissonmin, poissonmax, &
                              ispec, has_vs_zero, &
                              NSPEC_AB, kappastore, mustore, rhostore)
     else
       call exit_MPI(myrank,'Poroelastic elements in specfem_interface_mod not implement yet')
     endif

     !! min/max for whole cpu partition
     vpmin_glob = min(vpmin_glob, vpmin)
     vpmax_glob = max(vpmax_glob, vpmax)

     vsmin_glob = min(vsmin_glob, vsmin)
     vsmax_glob = max(vsmax_glob, vsmax)

     poissonmin_glob = min(poissonmin_glob,poissonmin)
     poissonmax_glob = max(poissonmax_glob,poissonmax)

     ! computes minimum and maximum size of this grid cell
     call get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                              NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

     elemsize_min_glob = min(elemsize_min_glob, elemsize_min)
     elemsize_max_glob = max(elemsize_max_glob, elemsize_max)


     ! estimation of minimum period resolved
     ! based on average GLL distance within element and minimum velocity
     !
     ! rule of thumb (Komatitsch et al. 2005):
     ! "average number of points per minimum wavelength in an element should be around 5."

     ! average distance between GLL points within this element
     avg_distance = elemsize_max / ( NGLLX - 1)  ! since NGLLX = NGLLY = NGLLZ

     ! largest possible minimum period such that number of points per minimum wavelength
     ! npts = ( min(vpmin,vsmin)  * pmax ) / avg_distance  is about ~ NPTS_PER_WAVELENGTH
     !
     ! note: obviously, this estimation depends on the choice of points per wavelength
     !          which is empirical at the moment.
     !          also, keep in mind that the minimum period is just an estimation and
     !          there is no such sharp cut-off period for valid synthetics.
     !          seismograms become just more and more inaccurate for periods shorter than this estimate.
     vel_min = min( vpmin,vsmin)
     pmax = avg_distance / vel_min * NPTS_PER_WAVELENGTH
     pmax_glob = max(pmax_glob,pmax)

     ! computes minimum and maximum distance of neighbor GLL points in this grid cell
     call get_GLL_minmaxdistance(distance_min,distance_max,ispec, &
                                 NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

     distance_min_glob = min(distance_min_glob, distance_min)
     distance_max_glob = max(distance_max_glob, distance_max)

     ! Courant number
     ! based on minimum GLL point distance and maximum velocity
     ! i.e. on the maximum ratio of ( velocity / gridsize )
     cmax = max(vpmax,vsmax) * DT / distance_min
     cmax_glob = max(cmax_glob,cmax)

     ! suggested timestep
     vel_max = max( vpmax,vsmax )
     dt_suggested = COURANT_SUGGESTED * distance_min / vel_max
     dt_suggested_glob = min( dt_suggested_glob, dt_suggested)

  enddo

  ! Vp velocity
  vpmin = vpmin_glob
  vpmax = vpmax_glob
  call min_all_cr(vpmin,vpmin_glob)
  call max_all_cr(vpmax,vpmax_glob)

  ! Vs velocity
  vsmin = vsmin_glob
  if (has_vs_zero) vsmin = 0.0

  vsmax = vsmax_glob
  call min_all_cr(vsmin,vsmin_glob)
  call max_all_cr(vsmax,vsmax_glob)

  ! Poisson's ratio
  poissonmin = poissonmin_glob
  poissonmax = poissonmax_glob
  call min_all_cr(poissonmin,poissonmin_glob)
  call max_all_cr(poissonmax,poissonmax_glob)

  ! GLL point distance
  distance_min = distance_min_glob
  distance_max = distance_max_glob
  call min_all_cr(distance_min,distance_min_glob)
  call max_all_cr(distance_max,distance_max_glob)

  ! element size
  elemsize_min = elemsize_min_glob
  elemsize_max = elemsize_max_glob
  call min_all_cr(elemsize_min,elemsize_min_glob)
  call max_all_cr(elemsize_max,elemsize_max_glob)

  ! model dimensions
  x_min_glob = minval(xstore)
  x_max_glob = maxval(xstore)

  y_min_glob = minval(ystore)
  y_max_glob = maxval(ystore)

  z_min_glob = minval(zstore)
  z_max_glob = maxval(zstore)

  ! min and max dimensions of the model
  x_min = x_min_glob
  x_max = x_max_glob
  call min_all_cr(x_min,x_min_glob)
  call max_all_cr(x_max,x_max_glob)

  y_min = y_min_glob
  y_max = y_max_glob
  call min_all_cr(y_min,y_min_glob)
  call max_all_cr(y_max,y_max_glob)

  z_min = z_min_glob
  z_max = z_max_glob
  call min_all_cr(z_min,z_min_glob)
  call max_all_cr(z_max,z_max_glob)

  ! minimum period
  pmax = pmax_glob
  call max_all_cr(pmax,pmax_glob)

  ! time step
  dt_suggested = dt_suggested_glob
  call min_all_cr(dt_suggested,dt_suggested_glob)

  if (myrank == 0 ) then
     !! CHECK POISSON'S RATIO OF NEW MODEL
     if (poissonmin_glob < -1.0000001d0 .or. poissonmax_glob > 0.50000001d0) then
        ModelIsSuitable=.false.
     endif

     !! CHECK STABILITY FOR NEW MODEL
     if (DT > dt_suggested) then
        ModelIsSuitable=.false.
     endif
  endif
  call bcast_all_singlel(ModelIsSuitable)

  if (ELASTIC_SIMULATION) then
     !! nothing to do
  else if (ACOUSTIC_SIMULATION) then
     deallocate(rho_vp,rho_vs)
  endif

  end subroutine CheckModelSuitabilityForModeling


!-------------------------------------------------------------------
!> create dummy file in order to initialize specfem before setting the right parameters
!! need to improve this by using because that files are not working for every cases.
!! may be use the flag INVERSE_FWI_FULL_PROBLEM
!-------------------------------------------------------------------

  subroutine CreateInitDummyFiles()

  use specfem_par

  implicit none

  integer                                               :: i

  open(666,file=trim(prefix_to_path)//'SEM/XX.S0001.BXX.adj')
  open(667,file=trim(prefix_to_path)//'SEM/XX.S0001.BXY.adj')
  open(668,file=trim(prefix_to_path)//'SEM/XX.S0001.BXZ.adj')

  do i = 1, NSTEP
     write(666,*) (i-1)*dt,0.
     write(667,*) (i-1)*dt,0.
     write(668,*) (i-1)*dt,0.
  enddo

  close(666)
  close(667)
  close(668)

  open(666,file=trim(prefix_to_path)//'SEM/XX.S0001.CXX.adj')
  open(667,file=trim(prefix_to_path)//'SEM/XX.S0001.CXY.adj')
  open(668,file=trim(prefix_to_path)//'SEM/XX.S0001.CXZ.adj')

  do i = 1, NSTEP
     write(666,*) (i-1)*dt,0.
     write(667,*) (i-1)*dt,0.
     write(668,*) (i-1)*dt,0.
  enddo

  close(666)
  close(667)
  close(668)


  open(666,file=trim(prefix_to_path)//'SEM/XX.S0001.HXX.adj')
  open(667,file=trim(prefix_to_path)//'SEM/XX.S0001.HXY.adj')
  open(668,file=trim(prefix_to_path)//'SEM/XX.S0001.HXZ.adj')

  do i = 1, NSTEP
     write(666,*) (i-1)*dt,0.
     write(667,*) (i-1)*dt,0.
     write(668,*) (i-1)*dt,0.
  enddo

  close(666)
  close(667)
  close(668)

  open(666,file=trim(prefix_to_path)//'SEM/XX.S0001.FXX.adj')
  open(667,file=trim(prefix_to_path)//'SEM/XX.S0001.FXY.adj')
  open(668,file=trim(prefix_to_path)//'SEM/XX.S0001.FXZ.adj')

  do i = 1, NSTEP
     write(666,*) (i-1)*dt,0.
     write(667,*) (i-1)*dt,0.
     write(668,*) (i-1)*dt,0.
  enddo

  close(666)
  close(667)
  close(668)

  !! pb il faut ecrire une station qui appartient au maillage pour initialiser le solver specfem
  !! : j'ai mis des points au hasard pour l'instant ...
  open(666,file=trim(prefix_to_path)//'DATA/STATIONS')
  open(667,file=trim(prefix_to_path)//'DATA/STATIONS_ADJOINT')

  write(666,'(a82)') 'S0001  XX 1100 4100 0 0 '
  write(666,'(a82)') 'S0001  XX 4100 1100 0 0 '
  write(666,'(a82)') 'S0001  XX    0    0 0 0 '
  write(666,'(a82)') 'S0001  XX 1000 1000 0 0 '
  write(666,'(a82)') 'S0001  XX 1000 4100 0 0 '
  write(666,'(a82)') 'S0001  XX  100 4100 0 0 '
  write(666,'(a82)') 'S0001  XX 4100  100 0 0 '

  write(667,'(a82)') 'S0001  XX 1100 3100 0 0 '
  write(667,'(a82)') 'S0001  XX 4100 1100 0 0 '
  write(667,'(a82)') 'S0001  XX    0    0 0 0 '
  write(667,'(a82)') 'S0001  XX 1000 1000 0 0 '
  write(667,'(a82)') 'S0001  XX 1000 4100 0 0 '
  write(667,'(a82)') 'S0001  XX  100 4100 0 0 '
  write(667,'(a82)') 'S0001  XX 4100 100 0 0 '


  close(666)
  close(667)

  end subroutine CreateInitDummyFiles

end module specfem_interface
