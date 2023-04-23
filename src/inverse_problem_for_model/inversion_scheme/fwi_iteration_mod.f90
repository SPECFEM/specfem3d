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

module fwi_iteration

  !! IMPORT VARIABLES FROM SPECFEM -------------------------------------------------------------------------------------------------
  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, NSPEC_AB, GPU_MODE
  !---------------------------------------------------------------------------------------------------------------------------------

  use inverse_problem_par

  use inversion_scheme
  use family_parameter
  use specfem_interface
  use regularization_interface
  use precond_mod
  use input_output

  ! arrays
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable :: initial_model, current_model
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable :: prior_model, ref_model
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable :: initial_gradient, current_gradient
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable :: regularization_penalty, gradient_regularization_penalty
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable :: descent_direction
  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable :: fwi_precond, hess_approxim

  !! for line search
  real(kind=CUSTOM_REAL), private                                    :: Q0, Qt, Qp0, Qpt

  public :: OneIterationOptim, InitializeOptimIteration, AllocatememoryForFWI

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-------------------------------------------------------------------------------------------------------------
! perform one step for optimization scheme iterations
!-------------------------------------------------------------------------------------------------------------

  subroutine OneIterationOptim(finished, acqui_simu, inversion_param) !!!!!! , regularization_fd)

    implicit none

!!!!!!    type(regul),  dimension(:), intent(in)    :: regularization_fd
    type(inver),                intent(inout) :: inversion_param
    type(acqui),  dimension(:), intent(inout) :: acqui_simu
    logical,                    intent(inout) :: finished
    !! locals
    real(kind=CUSTOM_REAL)                    :: mwl1, mwl2, td, tg
    real(kind=CUSTOM_REAL)                    :: step_length
    real(kind=CUSTOM_REAL)                    :: NormGrad
    integer                                   :: ievent, Niv, global_iter
    logical                                   :: flag_wolfe
    logical                                   :: ModelIsSuitable
    character(len=MAX_LEN_STRING)             :: prefix_name
    real(kind=CUSTOM_REAL)                    :: tmp_val
    ! timing
    double precision                          :: tCPU,tstart
    integer                                   :: ihours,iminutes,iseconds,int_tCPU
    double precision,                external :: wtime

    ! log output
    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) '*********************************************'
      write(INVERSE_LOG_FILE,*) '***   FWI ITERATION     : ', iter_inverse,'  ***'
      write(INVERSE_LOG_FILE,*) '*********************************************'
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    !! set up parameters for Wolfe's rule ---------------------------------------
    mwl1 =  inversion_param%m1
    mwl2 =  inversion_param%m2
    td = 0._CUSTOM_REAL
    tg = 0._CUSTOM_REAL
    step_length = inversion_param%current_step_length

    iter_wolfe = 0
    flag_wolfe = .false.
    finished = .false.
    ModelIsSuitable = .true.

    Niv = inversion_param%NinvPar

    ! compute and store descent direction
    call ComputeDescentDirection(iter_inverse, descent_direction, fwi_precond)

    ! compute q'(0) = grad . descent_direction for line search
    call Parallel_ComputeInnerProduct(initial_gradient, descent_direction, Niv, Qp0)

    ! save model or gradient or descent direction if asked by user
    call DumpArraysMeshSpecfem(initial_model, initial_gradient, descent_direction, inversion_param)

    !!  Wolfe's sub-iterations  ---------------------------------------------------
    do

      !! Next Wolfe's sub iteration
      iter_wolfe = iter_wolfe + 1

      !! exit if too many failure of Wolfe's sub-iterations
      if (iter_wolfe > inversion_param%Niter_wolfe) then
        finished = .true.
        ! log output
        if (myrank == 0) then
          write(OUTPUT_ITERATION_FILE,*)
          write(OUTPUT_ITERATION_FILE,*) ' STOP : Line search reached maximum sub-iterations allowed  :', &
                                          inversion_param%Niter_wolfe
          write(OUTPUT_ITERATION_FILE,*)
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) ' STOP : Line search reached maximum sub-iterations allowed  :', &
                                    inversion_param%Niter_wolfe
          write(INVERSE_LOG_FILE,*)
        endif
        ! stop, all done
        return
      endif

      ! choose initial step length
      call InitialGuessStep(inversion_param, step_length)

      ! log output
      if (myrank == 0) then
        write(INVERSE_LOG_FILE,*) '*** WOLFE SUB-ITERATION : ', iter_wolfe
        write(INVERSE_LOG_FILE,*) '***          TRIAL STEP : ', step_length
        write(INVERSE_LOG_FILE,*)
        call flush_iunit(INVERSE_LOG_FILE)
      endif

      ! update model for choosen family
      call UpdateModel(inversion_param, step_length, ModelIsSuitable)

      ! if model is not suitable for modeling then try smaller step
      if (.not. ModelIsSuitable) then
        ! write(*,*) 'Model new is not suitable for simulation ', myrank
        step_length = 0.45 * (tg + step_length)
        cycle
      endif

      ! timing
      tstart = wtime()

      ! compute cost function and gradient
      call InitForOneStepFWI(inversion_param)

      do ievent = 1,acqui_simu(1)%nevent_tot
        call ComputeGradientPerEvent(ievent, acqui_simu, inversion_param)
      enddo

      if (GPU_MODE) call TransfertKernelFromGPUArrays()

      ! timing elapsed time
      if (myrank == 0) then
        tCPU = wtime() - tstart
        int_tCPU = int(tCPU)
        ihours = int_tCPU / 3600
        iminutes = (int_tCPU - 3600*ihours) / 60
        iseconds = int_tCPU - 3600*ihours - 60*iminutes
        write(INVERSE_LOG_FILE,"('  all gradients done for frequency band : ',i2,' iteration : ',i3,' Wolfe : ',i3)") &
                                 iter_frq,iter_inverse,iter_wolfe
        write(INVERSE_LOG_FILE,*) ' Elapsed time in seconds  = ',tCPU
        write(INVERSE_LOG_FILE,"('  Elapsed time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
                                 ihours,iminutes,iseconds
        write(INVERSE_LOG_FILE,*)
        call flush_iunit(INVERSE_LOG_FILE)
      endif

      ! store current value of cost function for reduction
      Qt = inversion_param%total_current_cost

      ! communicate gradient and cost function to all simultaneous runs
      call mpi_sum_grad_all_to_all_simultaneous_runs(Qt)

      !! compute standard deviation on data
      tmp_val =  inversion_param%nb_data_std
      call sum_all_all_cr_for_simulatenous_runs(tmp_val, inversion_param%nb_data_std,1)
      tmp_val =  inversion_param%data_std
      call sum_all_all_cr_for_simulatenous_runs(tmp_val, inversion_param%data_std,1)

      if (abs(inversion_param%nb_data_std) > 0._CUSTOM_REAL) then
        inversion_param%data_std = sqrt(inversion_param%data_std / inversion_param%nb_data_std)
      else
        write(*,*) 'Error: no data compared for inversion ',inversion_param%nb_data_std
        stop
      endif

      ! store current gradient in choosen family parameter
      call StoreGradientInfamilyParam(inversion_param, current_gradient, current_model, ref_model, hess_approxim)

      ! compute regularization term and gradient for choosen family
      call AddRegularization(inversion_param, current_model, prior_model, regularization_penalty, &
                             gradient_regularization_penalty, myrank)

      ! add penalty term on cost function
      call Parallel_ComputeL2normSquare(regularization_penalty, inversion_param%NinvPar, inversion_param%cost_penalty)

      Qt = Qt + 0.5 * inversion_param%weight_Tikhonov_normalized * inversion_param%cost_penalty

      current_gradient(:,:,:,:,:) = current_gradient(:,:,:,:,:) + &
                                    inversion_param%weight_Tikhonov_normalized * gradient_regularization_penalty(:,:,:,:,:)

      ! store current cost funtion value
      inversion_param%total_current_cost = Qt

      ! compute q'(t) = grad . descent_direction for line search
      call Parallel_ComputeInnerProduct(current_gradient, descent_direction, Niv, Qpt)

      ! information about costs at current sub-iteration
      if (myrank == 0) then
        ! fwi infos
        write(OUTPUT_FWI_LOG, '(2i10, 2x, 8e18.8 )') iter_inverse, iter_wolfe, step_length, Qt, Q0, Qpt, Qp0, &
                    Qt - 0.5 * inversion_param%weight_Tikhonov_normalized * inversion_param%cost_penalty, &
                    0.5 * inversion_param%weight_Tikhonov_normalized * inversion_param%cost_penalty,  inversion_param%data_std
        call flush_iunit(OUTPUT_FWI_LOG)

        ! log output
        write(INVERSE_LOG_FILE,*) ' current total cost:  ', inversion_param%total_current_cost
        write(INVERSE_LOG_FILE,*) ' current step      :  maximum perturbations      : ', &
                          inversion_param%stats_max_relative_pert(1:inversion_param%Ninvpar),'%'
        call flush_iunit(INVERSE_LOG_FILE)
      endif

      ! write precond and regularization on disk
      if (DEBUG_MODE) then
        if (mygroup <= 0) then
          ! combined index for current frequency stage and inversion iteration
          global_iter = iter_frq * 10000  + iter_inverse
          !prefix_name = 'precond'
          !call DumpArray(fwi_precond, inversion_param, iter_inverse, prefix_name)
          !prefix_name = 'Hess_app'
          !call DumpArray(hess_approxim, inversion_param, iter_inverse, prefix_name)
          prefix_name = 'Grad_Regul'
          call DumpArray(gradient_regularization_penalty, inversion_param, global_iter, prefix_name)
          prefix_name = 'Regul'
          call DumpArray(regularization_penalty, inversion_param, iter_inverse, prefix_name)
        endif
      endif

      ! apply wolfe's rules ---------------------------------------------------------
      call wolfe_rules(mwl1, mwl2, Q0, Qt, Qp0, Qpt, step_length, td, tg, flag_wolfe)

      ! test for exiting ------------------------------------------------------------
      if (flag_wolfe) then
        ! step accepted
        ! define preconditionnner or taper on gradients
        call SetPrecond(inversion_param, current_gradient, hess_approxim, fwi_precond)

        ! store new model and gradient in l-bfgs history for the choosen family parameter
        call StoreModelAndGradientForLBFGS(current_model, current_gradient, iter_inverse+1)

        ! update initial model and gradient
        initial_model(:,:,:,:,:) = current_model(:,:,:,:,:)
        initial_gradient(:,:,:,:,:) = current_gradient(:,:,:,:,:)

        ! update initial costs
        Q0 = Qt
        inversion_param%current_step_length = step_length

        call Parallel_ComputeL2normSquare(current_gradient, Niv, NormGrad)

        !! output information -------------------
        if (myrank == 0) then
          ! iterations output
          write(OUTPUT_ITERATION_FILE,'(i7,"|",e15.8,"|",e15.8,"|",e12.5,"|",e12.5,"|",i3)')  &
              iter_inverse+1, Q0, NormGrad, Q0/inversion_param%Cost_init, NormGrad/inversion_param%Norm_grad_init, iter_wolfe
          call flush_iunit(OUTPUT_ITERATION_FILE)

          ! log output
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) ' Cost function : ', Q0, '   relative cost :', &
                                    100 * Q0/inversion_param%Cost_init,'%'
          write(INVERSE_LOG_FILE,*) ' Gradient Norm : ', NormGrad, '   relative grad :', &
                                    100 * NormGrad/inversion_param%Norm_grad_init, '%'
          write(INVERSE_LOG_FILE,*) ' Data std      : ', inversion_param%data_std
          write(INVERSE_LOG_FILE,*)
          call flush_iunit(INVERSE_LOG_FILE)
        endif

        !! stopping criteria -----------------------------------------------------
        ! zero gradient
        if (inversion_param%Norm_grad_init == 0._CUSTOM_REAL) then
          finished = .true.
          ! log output
          if (myrank == 0) then
            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*) &
              ' FWI STOP : zero gradient :',inversion_param%Norm_grad_init
            write(INVERSE_LOG_FILE,*)
            call flush_iunit(INVERSE_LOG_FILE)
            write(OUTPUT_ITERATION_FILE,*)
            write(OUTPUT_ITERATION_FILE,*) &
              ' FWI STOP : zero gradient :',inversion_param%Norm_grad_init
            write(OUTPUT_ITERATION_FILE,*)
            call flush_iunit(OUTPUT_ITERATION_FILE)
          endif
          ! stop iteration
          return
        endif

        ! sufficient decrease of cost function
        if ( (Q0/inversion_param%Cost_init) <= inversion_param%relat_cost) then
          finished = .true.
          ! log output
          if (myrank == 0) then
            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*) &
              ' FWI STOP : cost function reached maximum relative decrease allowed  :',inversion_param%relat_cost
            write(INVERSE_LOG_FILE,*)
            call flush_iunit(INVERSE_LOG_FILE)
            write(OUTPUT_ITERATION_FILE,*)
            write(OUTPUT_ITERATION_FILE,*) &
              ' FWI STOP : cost function reached maximum relative decrease allowed  :',inversion_param%relat_cost
            write(OUTPUT_ITERATION_FILE,*)
            call flush_iunit(OUTPUT_ITERATION_FILE)
          endif
        endif

        ! sufficient decrease of gradient
        if ( (NormGrad/inversion_param%Norm_grad_init) <= inversion_param%relat_grad) then
          finished = .true.
          ! log output
          if (myrank == 0) then
            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*) &
              ' FWI STOP : gradient of cost function reached maximum relative decrease allowed  :',inversion_param%relat_grad
            write(INVERSE_LOG_FILE,*)
            call flush_iunit(INVERSE_LOG_FILE)
            write(OUTPUT_ITERATION_FILE,*)
            write(OUTPUT_ITERATION_FILE,*) &
              ' FWI STOP : gradient of cost function reached maximum relative decrease allowed  :',inversion_param%relat_grad
            write(OUTPUT_ITERATION_FILE,*)
            call flush_iunit(OUTPUT_ITERATION_FILE)
          endif
        endif

        ! all done
        return !! return to the main iteration loop of FWI

      endif

    enddo

  end subroutine OneIterationOptim

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-------------------------------------------------------------------------------------------------------------
! initialize optimization scheme iterations
!-------------------------------------------------------------------------------------------------------------

  subroutine InitializeOptimIteration(acqui_simu, inversion_param)

    implicit none

    type(acqui),  dimension(:), intent(inout) :: acqui_simu
    type(inver),                intent(inout) :: inversion_param
    ! local
    integer                                   :: ievent,ipar,global_iter
    character(len=MAX_LEN_STRING)             :: prefix_name
    real(kind=CUSTOM_REAL)                    :: tmp_val
    ! timing
    double precision                          :: tCPU,tstart
    integer                                   :: ihours,iminutes,iseconds,int_tCPU
    double precision,                external :: wtime

    ! initializes iterators
    iter_inverse = 0
    iter_wolfe   = 0

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) '*********************************************'
      write(INVERSE_LOG_FILE,*) '***      INITIALIZE FWI ITERATIONS        ***'
      write(INVERSE_LOG_FILE,*) '*********************************************'
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*)  ' frequency group  :',iter_frq
      if (inversion_param%use_band_pass_filter) then
        write(INVERSE_LOG_FILE,*)  ' band pass filter :',acqui_simu(1)%fl_event(iter_frq),'/',acqui_simu(1)%fh_event(iter_frq)
      endif
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! initialize regularization
    call SetUpRegularization(inversion_param, acqui_simu, myrank)

    ! timing
    tstart = wtime()

    ! compute cost function and gradient---------
    call InitForOneStepFWI(inversion_param)

    do ievent = 1,acqui_simu(1)%nevent_tot
      call ComputeGradientPerEvent(ievent, acqui_simu, inversion_param)
    enddo

    if (GPU_MODE) call TransfertKernelFromGPUArrays()

    ! timing elapsed time
    if (myrank == 0) then
      tCPU = wtime() - tstart
      int_tCPU = int(tCPU)
      ihours = int_tCPU / 3600
      iminutes = (int_tCPU - 3600*ihours) / 60
      iseconds = int_tCPU - 3600*ihours - 60*iminutes
      write(INVERSE_LOG_FILE,"('  all gradients done for frequency band : ',i2,' iteration : ',i3)") &
                              iter_frq,iter_inverse
      write(INVERSE_LOG_FILE,*) ' Elapsed time in seconds  = ',tCPU
      write(INVERSE_LOG_FILE,"('  Elapsed time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
                              ihours,iminutes,iseconds
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! store current value of cost function for reduction
    Q0 = inversion_param%total_current_cost

    ! communicate gradient and cost function to all simultaneous runs
    call mpi_sum_grad_all_to_all_simultaneous_runs(Q0)

    !! compute standard deviation on data
    tmp_val =  inversion_param%nb_data_std
    call sum_all_all_cr_for_simulatenous_runs(tmp_val, inversion_param%nb_data_std,1)
    tmp_val =  inversion_param%data_std
    call sum_all_all_cr_for_simulatenous_runs(tmp_val, inversion_param%data_std,1)

    if (abs(inversion_param%nb_data_std) > 0._CUSTOM_REAL) then
      inversion_param%data_std = sqrt(inversion_param%data_std / inversion_param%nb_data_std)
    else
      write(*,*) 'Error: no data compared for inversion ',inversion_param%nb_data_std
      stop
    endif

    ! set model from modeling as reference model
    if (iter_frq == 1) call Modeling2RefInvert(inversion_param, ref_model)

    ! store initial model in choosen family parameter
    call SpecfemParam2Invert(inversion_param, initial_model, ref_model)

    ! store initial gradient in choosen family parameter
    call StoreGradientInFamilyParam(inversion_param, initial_gradient, initial_model, ref_model, hess_approxim)

    if (inversion_param%input_SEM_prior) then
      ! save specific read prior model in family
      write(*,*) ' STOPPING : do not use this option for now (not fully implemented yet): input sem prior'
      stop
      ! (not working yet)
      !call  SpecfemPrior2Invert(inversion_param, prior_model)
    else
      ! save starting model read as prior model for the first frequency group
      if (iter_frq == 1) then
        prior_model(:,:,:,:,:) = initial_model(:,:,:,:,:)
      endif
    endif

    !! compute reference mean model
    if (inversion_param%use_damping_SEM_Tikhonov .and. inversion_param%input_SEM_prior) then
      do ipar = 1, inversion_param%NinvPar
        inversion_param%damp_weight(ipar) = inversion_param%damp_weight(ipar) * &
                                            (size(prior_model(:,:,:,:,ipar)) / sum(prior_model(:,:,:,:,ipar)))**2
      enddo
    endif

    ! compute regularization term and gradient for choosen family
    call AddRegularization(inversion_param, initial_model, prior_model, regularization_penalty, &
                           gradient_regularization_penalty, myrank)

    ! add penalty term on cost function
    call Parallel_ComputeL2normSquare(regularization_penalty, inversion_param%NinvPar, inversion_param%cost_penalty)

    Q0 = Q0 + 0.5 * inversion_param%weight_Tikhonov_normalized * inversion_param%cost_penalty

    initial_gradient(:,:,:,:,:) = initial_gradient(:,:,:,:,:) + &
                                  inversion_param%weight_Tikhonov_normalized * gradient_regularization_penalty(:,:,:,:,:)

    ! store cost function value
    inversion_param%total_current_cost = Q0

    ! store intial values of cost function
    inversion_param%Cost_init = inversion_param%total_current_cost

    ! store initial gradient norm
    call Parallel_ComputeL2normSquare(initial_gradient, inversion_param%NinvPar, inversion_param%Norm_grad_init)

    ! define preconditioner or taper on gradients
    call SetPrecond(inversion_param, initial_gradient, hess_approxim, fwi_precond)

    ! write precond on disk
    if (DEBUG_MODE) then
      if (mygroup <= 0) then
        ! combined index for current frequency stage and inversion iteration
        global_iter = iter_frq * 10000  + iter_inverse
        ! preconditioner
        prefix_name = 'Precond'
        call DumpArray(fwi_precond, inversion_param, global_iter, prefix_name)
        prefix_name = 'Hess_app'
        call DumpArray(hess_approxim, inversion_param, global_iter, prefix_name)
        prefix_name = 'Grad_Regul'
        call DumpArray(gradient_regularization_penalty, inversion_param, global_iter, prefix_name)
        prefix_name = 'Regul'
        call DumpArray(regularization_penalty, inversion_param, global_iter, prefix_name)
        if (inversion_param%use_damping_SEM_Tikhonov &
            .or. inversion_param%use_variable_SEM_damping_sources &
            .or. inversion_param%use_variable_SEM_damping_stations) then
          prefix_name = 'Spatial_damp'
          descent_direction(:,:,:,:,1) = spatial_damping(:,:,:,:) !! use as temporary working array
          call DumpArray(descent_direction, inversion_param, global_iter, prefix_name)
        endif
      endif
    endif

    ! store new model and gradient in l-bfgs history for the choosen family parameter
    call StoreModelAndGradientForLBFGS(initial_model, initial_gradient, 0)

    if (myrank == 0) then
      ! fwi log
      write(OUTPUT_FWI_LOG,*) '#'
      if (inversion_param%use_band_pass_filter) then
        write(OUTPUT_FWI_LOG,*) '# FREQUENCY GROUP :',acqui_simu(1)%fl_event(iter_frq),acqui_simu(1)%fh_event(iter_frq)
      else
        write(OUTPUT_FWI_LOG,*) '# FREQUENCY GROUP :',iter_frq
      endif
      write(OUTPUT_FWI_LOG,*) '#'
      write(OUTPUT_FWI_LOG,'(a100)') &
        '#iter_inverse #iter_wolfe #step_length #Qt #Q0 #Qpt #Qp0 #Qtminuscostpenalty #costpenalty #data_std'
      call flush_iunit(OUTPUT_FWI_LOG)

      ! iteration log
      write(OUTPUT_ITERATION_FILE,*)
      write(OUTPUT_ITERATION_FILE,*) '#'
      if (inversion_param%use_band_pass_filter) then
        write(OUTPUT_ITERATION_FILE,*) '# FREQUENCY GROUP : ',acqui_simu(1)%fl_event(iter_frq),acqui_simu(1)%fh_event(iter_frq)
      else
        write(OUTPUT_ITERATION_FILE,*) '# FREQUENCY GROUP : ',iter_frq
      endif
      write(OUTPUT_ITERATION_FILE,*) '#'
      write(OUTPUT_ITERATION_FILE,'(i7,"|",e15.8,"|",e15.8,"|",f12.5,"|",f12.5,"|")')  &
          0, inversion_param%Cost_init, inversion_param%Norm_grad_init, 1., 1.
      call flush_iunit(OUTPUT_ITERATION_FILE)

      ! inversion log file
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' Initial Cost function : ', inversion_param%Cost_init
      write(INVERSE_LOG_FILE,*) ' Initial Gradient Norm : ', inversion_param%Norm_grad_init
      write(INVERSE_LOG_FILE,*) ' Initial Data std      : ', inversion_param%data_std
      call flush_iunit(INVERSE_LOG_FILE)
    endif

  end subroutine InitializeOptimIteration

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------
! compute current model form initial model and descent direction for a given step length
!-------------------------------------------------------------------------------------------------------------
!
! Since Specfem compute kernel for log of parameter, here we use this convention by default
! BE CAREFUL : the model is stored in parameter and *NOT* log parameter in specfem mesh (rhostore, kappastore, ...)
!
  subroutine UpdateModel(inversion_param, step_length, ModelIsSuitable)

    type(inver),                                    intent(inout) :: inversion_param
    logical,                                        intent(inout) :: ModelIsSuitable
    real(kind=CUSTOM_REAL),                         intent(in)    :: step_length
    integer                                                       :: ipar
    !integer                                                       :: igll, jgll, kgll, ispec

    real(kind=CUSTOM_REAL)                                        :: vmin, vmin_glob
    real(kind=CUSTOM_REAL)                                        :: vmax, vmax_glob
    real(kind=CUSTOM_REAL)                                        :: vmax_rel_ini,vmax_rel_ini_glob
    real(kind=CUSTOM_REAL)                                        :: vmax_rel_pri,vmax_rel_pri_glob

    ! model update
    current_model(:,:,:,:,:) = initial_model(:,:,:,:,:) + step_length * descent_direction(:,:,:,:,:)

    !! store the model on specfem arrays to perform next simulation
    call InvertParam2Specfem(inversion_param, current_model, ref_model)

    call CheckModelSuitabilityForModeling(ModelIsSuitable)

    ! user output
    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*) ' - > update model :  '
      if (.not. ModelIsSuitable) then
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*) ' - > updated model not suitable for simulation trying smaller step  '
        write(INVERSE_LOG_FILE,*)
      endif
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! determines relative perturbations
    do ipar = 1, inversion_param%NinvPar

      select case(inversion_param%parameter_metric)
      case(0)
        !! directly the parameter P
        vmin = minval(current_model(:,:,:,:,ipar))
        vmax = maxval(current_model(:,:,:,:,ipar))
        call min_all_cr(vmin,vmin_glob)
        call max_all_cr(vmax,vmax_glob)

        vmax_rel_ini = maxval( abs(current_model(:,:,:,:,ipar) - initial_model(:,:,:,:,ipar)) / initial_model(:,:,:,:,ipar) )
        vmax_rel_pri = maxval( abs(current_model(:,:,:,:,ipar) - prior_model(:,:,:,:,ipar)) / prior_model(:,:,:,:,ipar) )
        call max_all_cr(vmax_rel_ini,vmax_rel_ini_glob)
        call max_all_cr(vmax_rel_pri,vmax_rel_pri_glob)

      case(1)
        !! use P / Pref
        vmin = minval(current_model(:,:,:,:,ipar) * ref_model(:,:,:,:,ipar))
        vmax = maxval(current_model(:,:,:,:,ipar) * ref_model(:,:,:,:,ipar))
        call min_all_cr(vmin,vmin_glob)
        call max_all_cr(vmax,vmax_glob)

        vmax_rel_ini = maxval( abs(current_model(:,:,:,:,ipar) - initial_model(:,:,:,:,ipar)) / initial_model(:,:,:,:,ipar) )
        vmax_rel_pri = maxval( abs(current_model(:,:,:,:,ipar) - prior_model(:,:,:,:,ipar)) / prior_model(:,:,:,:,ipar) )
        call max_all_cr(vmax_rel_ini,vmax_rel_ini_glob)
        call max_all_cr(vmax_rel_pri,vmax_rel_pri_glob)

      case(2)
        !! use log(P)
        vmin = minval(exp(current_model(:,:,:,:,ipar)))
        vmax = maxval(exp(current_model(:,:,:,:,ipar)))
        call min_all_cr(vmin,vmin_glob)
        call max_all_cr(vmax,vmax_glob)

        vmax_rel_ini = maxval( abs(exp(current_model(:,:,:,:,ipar)) - exp(initial_model(:,:,:,:,ipar))) /  &
                       exp(initial_model(:,:,:,:,ipar)) )
        vmax_rel_pri = maxval( abs(exp(current_model(:,:,:,:,ipar)) - exp(prior_model(:,:,:,:,ipar))) /  &
                       exp(prior_model(:,:,:,:,ipar)) )
        call max_all_cr(vmax_rel_ini,vmax_rel_ini_glob)
        call max_all_cr(vmax_rel_pri,vmax_rel_pri_glob)

      case(3)
        !! ??
        write(*,*) "Error: parameter metric ",inversion_param%parameter_metric," not implemented yet"
        stop
      end select

      ! note: for the first frequency band (iter_frq == 1), the prior model and the initial model are the same.
      !       the initial model refers to the model at the beginning of each iteration step.
      !       the prior model to the model at the beginning of an iteration stage, when multiple frequency bands are used.

      ! user output
      if (myrank == 0) then
        write(INVERSE_LOG_FILE,'( a13, i2, a10, 2(a8, f12.5))') &
            '  Parameter :', ipar, '  -  '//inversion_param%param_inv_name(ipar), &
            '   MIN :',vmin_glob,'   MAX :',vmax_glob
        write(INVERSE_LOG_FILE,'( a33, f10.6, a18, f10.6, a3)') &
             '                 max pert,  prior :', 100.0 * vmax_rel_pri_glob, &
             ' %     previous  :', 100.0 * vmax_rel_ini,' %'
        !write(INVERSE_LOG_FILE,*) ' max pert / starting model : ', 100*vmax_glob0,' %'
        !write(INVERSE_LOG_FILE,*) ' max pert / previous model : ', 100*vmin_glob0,' %'
        call flush_iunit(INVERSE_LOG_FILE)

        ! stores maximum relative perturbations (in %) to previous / initial model
        inversion_param%stats_max_relative_pert(ipar) = 100.0 * vmax_rel_ini_glob
      endif

    enddo

    ! user output
    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

  end subroutine UpdateModel
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------
! compute initial guess for step length to try for line search
!-------------------------------------------------------------------------------------------------------------

  subroutine InitialGuessStep(inversion_param, step_length)

    type(inver),                                    intent(in)    :: inversion_param
    real(kind=CUSTOM_REAL),                         intent(inout) :: step_length
    real(kind=CUSTOM_REAL)                                        :: max_val_tmp, max_val, max_val_model, step_length_guess

    ! gets maximum direction value
    max_val_tmp = maxval(abs(descent_direction(:,:,:,:,:)))
    call  max_all_all_cr(max_val_tmp, max_val)

    max_val_tmp   = maxval(abs(initial_model(:,:,:,:,:)))
    call  max_all_all_cr(max_val_tmp, max_val_model)

    ! safety check
    if (max_val <= 0._CUSTOM_REAL) then
      write(*,*) 'Error: initial guess step has zero descent direction :',max_val
      stop
    endif

    select case(inversion_param%parameter_metric)
    case(0)  !! directly the parameter P
      step_length_guess = inversion_param%max_relative_pert * max_val_model / max_val
    case(1)  !! P / Pref
      step_length_guess = inversion_param%max_relative_pert / max_val
    case(2)  !! log(P)
      step_length_guess = log(1. + inversion_param%max_relative_pert) / max_val
    case(3) !! log(P/Pref)
      step_length_guess = log(1. + inversion_param%max_relative_pert) / max_val
    end select

    ! sets step length for first Wolfe iteration
    if (iter_wolfe == 1) then
      ! for inversion iteration 0 and 1 uses guess, after we use the previous step length
      if (iter_inverse <= 1) then
        step_length = step_length_guess
      else if (iter_inverse <= 5) then
        ! try to use previous step length, unless new guess is (much) smaller or (much) bigger
        if (step_length_guess < 0.5 * step_length) step_length = step_length_guess
        if (step_length_guess > 2.5 * step_length) step_length = step_length_guess
      else
        ! try to use previous step length, unless new guess is (much) smaller
        if (step_length_guess < 0.5 * step_length) step_length = step_length_guess
      endif
    endif

    ! debug output
    if (DEBUG_MODE) then
      write(IIDD,* )
      write(IIDD,* ) 'STEP LENGTH:',  step_length,' GUESS: ',step_length_guess
      write(IIDD,* )
      call flush_iunit(IIDD)
    endif

  end subroutine InitialGuessStep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------
! Prepare and allocate all arrays used in inverse problem
!-------------------------------------------------------------------------------------------------------------
  subroutine AllocatememoryForFWI(inversion_param, nevent)

    type(inver),                                    intent(inout) :: inversion_param
    integer,                                        intent(in)    :: nevent
    ! local
    integer                                                       :: ier, Ninvpar, Mbfgs
    double precision                                              :: memory_size, memory_size_glob

    !! 1/ set the family parameter
    call PrepareArraysfamilyParam(inversion_param)

    Ninvpar = inversion_param%NinvPar
    Mbfgs = inversion_param%max_history_bfgs

    ! memory estimate
    memory_size = 0.d0
    ! bfgs_stored_gradient, bfgs_stored_model
    memory_size = memory_size &
      + 2.d0 * dble(Mbfgs+1) * dble(NGLLX * NGLLY * NGLLZ) * dble(NSPEC_ADJOINT) * dble(Ninvpar) * dble(CUSTOM_REAL)
    ! wks_1, wks_2, wks_1n, wks_2n
    memory_size = memory_size &
      + 4.d0 * dble(NGLLX * NGLLY * NGLLZ) * dble(NSPEC_ADJOINT) * dble(Ninvpar) * dble(CUSTOM_REAL)
    ! initial_model, prior_model, ..
    memory_size = memory_size &
      + 11.d0 * dble(NGLLX * NGLLY * NGLLZ) * dble(NSPEC_ADJOINT) * dble(Ninvpar) * dble(CUSTOM_REAL)
    ! current_cost_prime,.. (not used yet)
    !memory_size = memory_size &
    !  + 4.d0 * dble(NEVENT) * dble(CUSTOM_REAL)

    ! maximum of all processes (may vary e.g. due to different NSPEC_ADJOINT, ..)
    call max_all_dp(memory_size,memory_size_glob)

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' allocate arrays for fwi iterations,  Ninvpar          :', Ninvpar
      write(INVERSE_LOG_FILE,*) '                                      L-BFGS history   :', Mbfgs
      write(INVERSE_LOG_FILE,*) '                                      number of events :', nevent
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' minimum memory requested : ',sngl(memory_size_glob / 1024. / 1024.),'MB per process'
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) '*************************************************************************'
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    !! allocate arrays for L-BFGS inversion scheme
    call AllocateArraysForInversion(inversion_param)

    !! allocate arrays for fwi_iteration
    allocate(initial_model(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 544')
    if (ier /= 0) call exit_MPI(myrank,"error allocation initial_model in AllocatememoryForFWI  subroutine")
    initial_model(:,:,:,:,:) = 0._CUSTOM_REAL

    allocate(prior_model(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 545')
    if (ier /= 0) call exit_MPI(myrank,"error allocation prior_model in AllocatememoryForFWI  subroutine")
    prior_model(:,:,:,:,:) = 0._CUSTOM_REAL

    allocate(ref_model(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 546')
    if (ier /= 0) call exit_MPI(myrank,"error allocation ref_model in AllocatememoryForFWI  subroutine")
    ref_model(:,:,:,:,:) = 0._CUSTOM_REAL

    allocate(current_model(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 547')
    if (ier /= 0) call exit_MPI(myrank,"error allocation current_model in AllocatememoryForFWI  subroutine")
    current_model(:,:,:,:,:) = 0._CUSTOM_REAL

    allocate(initial_gradient(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 548')
    if (ier /= 0) call exit_MPI(myrank,"error allocation initial_gradient in AllocatememoryForFWI  subroutine")
    initial_gradient(:,:,:,:,:) = 0._CUSTOM_REAL

    allocate(current_gradient(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 549')
    if (ier /= 0) call exit_MPI(myrank,"error allocation current_gradient in AllocatememoryForFWI  subroutine")
    current_gradient(:,:,:,:,:) = 0._CUSTOM_REAL

    allocate(descent_direction(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 550')
    if (ier /= 0) call exit_MPI(myrank,"error allocation descent_direction in AllocatememoryForFWI  subroutine")
    descent_direction(:,:,:,:,:) = 0._CUSTOM_REAL

    allocate(fwi_precond(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 551')
    if (ier /= 0) call exit_MPI(myrank,"error allocation fwi_precond in AllocatememoryForFWI subroutine")
    fwi_precond(:,:,:,:,:) = 1._CUSTOM_REAL

    allocate(hess_approxim(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 552')
    if (ier /= 0) call exit_MPI(myrank,"error allocation hess_approxim in AllocatememoryForFWI subroutine")
    hess_approxim(:,:,:,:,:) = 1._CUSTOM_REAL

    allocate(regularization_penalty(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 553')
    if (ier /= 0) call exit_MPI(myrank,"error allocation regularization_penalty in AllocatememoryForFWI subroutine")
    regularization_penalty(:,:,:,:,:) = 0._CUSTOM_REAL

    allocate(gradient_regularization_penalty(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 554')
    if (ier /= 0) call exit_MPI(myrank,"error allocation gradient_regularization_penalty in AllocatememoryForFWI subroutine")
    gradient_regularization_penalty(:,:,:,:,:) = 0._CUSTOM_REAL

    ! costs per event (not used yet)
    !allocate(inversion_param%current_cost(NEVENT),  inversion_param%previous_cost(NEVENT),stat=ier)
    !if (ier /= 0) call exit_MPI_without_rank('error allocating array 556')
    !inversion_param%current_cost(:) = 0._CUSTOM_REAL
    !inversion_param%previous_cost(:) = 0._CUSTOM_REAL

    !allocate(inversion_param%current_cost_prime(NEVENT),  inversion_param%previous_cost_prime(NEVENT),stat=ier)
    !if (ier /= 0) call exit_MPI_without_rank('error allocating array 555')
    !inversion_param%current_cost_prime(:) = 0._CUSTOM_REAL
    !inversion_param%previous_cost_prime(:) = 0._CUSTOM_REAL

    ! model statistics
    allocate(inversion_param%stats_max_relative_pert(NinvPar),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,"error allocation stats_max_relative_pert array")
    inversion_param%stats_max_relative_pert(:) = 0._CUSTOM_REAL

  end subroutine AllocatememoryForFWI


end module fwi_iteration
