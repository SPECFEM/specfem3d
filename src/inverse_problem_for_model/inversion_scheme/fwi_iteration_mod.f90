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

  real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable :: initial_model, current_model, prior_model, ref_model
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

  subroutine OneIterationOptim(iter_inverse, iter_frq, finished, acqui_simu, inversion_param) !!!!!! , regularization_fd)

    implicit none

!!!!!!    type(regul),  dimension(:), allocatable,        intent(in)    :: regularization_fd
    type(inver),                                    intent(inout) :: inversion_param
    type(acqui),  dimension(:), allocatable,        intent(inout) :: acqui_simu
    integer,                                        intent(in)    :: iter_inverse, iter_frq
    logical,                                        intent(inout) :: finished

    !! locals
    real(kind=CUSTOM_REAL)                                        :: mwl1, mwl2, td, tg, step_length
    real(kind=CUSTOM_REAL)                                        :: NormGrad
    integer                                                       :: iter_wolfe, ievent, Niv
    logical                                                       :: flag_wolfe
    logical                                                       :: ModelIsSuitable
    character(len=MAX_LEN_STRING)                                 :: prefix_name
    real(kind=CUSTOM_REAL)                                        :: tmp_val

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '           *********************************************'
       write(INVERSE_LOG_FILE,*) '           ***   FWI ITERATION  : ', iter_inverse,'  ***'
       write(INVERSE_LOG_FILE,*) '           *********************************************'
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
    Niv=inversion_param%NinvPar
    inversion_param%current_iteration=iter_inverse
    inversion_param%current_ifrq=iter_frq

    ! compute and store descent direction
    call ComputeDescentDirection(iter_inverse, descent_direction, fwi_precond)

    ! compute q'(0) = grad . descent_direction for line search
    call Parallel_ComputeInnerProduct(initial_gradient, descent_direction, Niv, Qp0)

    ! save model or gradient or descent direction if asked by user
    call DumpArraysMeshSpecfem(initial_model, initial_gradient, descent_direction, iter_inverse, inversion_param)

    !!  Wolfe's sub-iterations  ---------------------------------------------------
    do

       !! Next Wolfe's sub iteration
       iter_wolfe = iter_wolfe + 1

       !! exit if too many failure of Wolfe's sub-iterations
       if (iter_wolfe > inversion_param%Niter_wolfe) then
          finished=.true.
          if (myrank == 0) then
             write(OUTPUT_ITERATION_FILE,*)
             write(OUTPUT_ITERATION_FILE,*)  ' STOP : Line search reached maximum sub-iterations allowed  :'&
                  ,  inversion_param%Niter_wolfe
             write(OUTPUT_ITERATION_FILE,*)
             write(INVERSE_LOG_FILE,*)
             write(INVERSE_LOG_FILE,*)  ' STOP : Line search reached maximum sub-iterations allowed  :', &
                  inversion_param%Niter_wolfe
             write(INVERSE_LOG_FILE,*)
          endif
          return
       endif

       ! choose initial step length
       ! for iteration 0 and 1 after we use the previous step length
       if (iter_inverse <= 1 .and. iter_wolfe == 1) call InitialGuessStep(inversion_param, step_length)

       if (myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '*** WOLFE SUB-ITERATION FWI ITERATION  : ', iter_wolfe
          write(INVERSE_LOG_FILE,*)'    TRY STEP :', step_length, ' ***'
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
       ! compute cost function and gradient---------
       call InitForOneStepFWI(inversion_param)
       do ievent=1,acqui_simu(1)%nevent_tot
          call ComputeGradientPerEvent(ievent, iter_inverse, acqui_simu, inversion_param)
       enddo
       if (GPU_MODE) call TransfertKernelFromGPUArrays()

       ! store current value of cost function for reduction
       Qt=inversion_param%total_current_cost
       ! communicate gradient and cost function to all simultaneous runs
       call mpi_sum_grad_all_to_all_simultaneous_runs(Qt)

       !! compute standard deviation on data
       tmp_val =  inversion_param%nb_data_std
       call sum_all_all_cr_for_simulatenous_runs(tmp_val, inversion_param%nb_data_std,1)
       tmp_val =  inversion_param%data_std
       call sum_all_all_cr_for_simulatenous_runs(tmp_val, inversion_param%data_std,1)
       inversion_param%data_std = sqrt(inversion_param%data_std / inversion_param%nb_data_std)

       ! store current gradient in choosen family parameter
       call StoreGradientInfamilyParam(inversion_param, current_gradient, current_model, ref_model, hess_approxim)

       ! compute regularization term and gradient for choosen family
       call AddRegularization(inversion_param, current_model, ref_model, prior_model, regularization_penalty, &
            gradient_regularization_penalty, &
            myrank)
       ! add penalty tern on cost function
       call Parallel_ComputeL2normSquare(regularization_penalty, inversion_param%NinvPar, inversion_param%cost_penalty)
       Qt = Qt + 0.5 * inversion_param%weight_Tikonov * inversion_param%cost_penalty
       current_gradient(:,:,:,:,:) = current_gradient(:,:,:,:,:) + &
            inversion_param%weight_Tikonov * gradient_regularization_penalty(:,:,:,:,:)
       ! store current cost funtion value
       inversion_param%total_current_cost=Qt

       ! compute q'(t) = grad . descent_direction for line search
       call Parallel_ComputeInnerProduct(current_gradient, descent_direction, Niv, Qpt)

       !information about costs at current sub-iteration
       if (myrank == 0) then
          write(OUTPUT_FWI_LOG, '(2i10, 2x, 8e20.10 )')  iter_inverse, iter_wolfe, step_length, Qt, Q0, Qpt, Qp0, &
                Qt - 0.5 * inversion_param%weight_Tikonov * inversion_param%cost_penalty, &
                0.5 * inversion_param%weight_Tikonov * inversion_param%cost_penalty,  inversion_param%data_std
          call flush_iunit(OUTPUT_FWI_LOG)
       endif

       ! write precond and regulariztion on disk
       if (mygroup <= 0 .and. VERBOSE_MODE) then
          !prefix_name='precond'
          !call DumpArray(fwi_precond, inversion_param, iter_inverse, prefix_name)
          !prefix_name='Hess_app'
          !call DumpArray(hess_approxim, inversion_param, iter_inverse, prefix_name)
          prefix_name='Grad_Regul'
          call DumpArray(gradient_regularization_penalty, inversion_param, iter_inverse, prefix_name)
          prefix_name='Regul'
          call DumpArray(regularization_penalty, inversion_param, iter_inverse, prefix_name)
       endif


       ! apply wolfe's rules ---------------------------------------------------------
       call wolfe_rules(mwl1, mwl2, Q0, Qt, Qp0, Qpt, step_length, td, tg, flag_wolfe)

       ! test for exiting ------------------------------------------------------------
       if (flag_wolfe) then

          ! define preconditionnner or taper on gradients
          call SetPrecond(iter_inverse, inversion_param, current_gradient, hess_approxim, fwi_precond)

          ! store new model and gradient in l-bfgs history for the choosen family parameter
          call StoreModelAndGradientForLBFGS(current_model, current_gradient, iter_inverse+1)

          initial_model(:,:,:,:,:)=current_model(:,:,:,:,:)
          initial_gradient(:,:,:,:,:)=current_gradient(:,:,:,:,:)
          Q0=Qt
          inversion_param%current_step_length=step_length

          call Parallel_ComputeL2normSquare(current_gradient, Niv, NormGrad)

          !! output information -------------------
          if (myrank == 0) then
             write(OUTPUT_ITERATION_FILE,'(i5,"|",e15.8,"|",e15.8,"|",e12.5,"|",e12.5,"|",i3)')  &
             1+iter_inverse, Q0, NormGrad, Q0/inversion_param%Cost_init, NormGrad/inversion_param%Norm_grad_init, iter_wolfe
             call flush_iunit(OUTPUT_ITERATION_FILE)
             write(INVERSE_LOG_FILE,*)
             write(INVERSE_LOG_FILE,*) ' Cost fuction  : ', Q0,  '   relative cost :', 100*Q0/inversion_param%Cost_init,'%'
             write(INVERSE_LOG_FILE,*) ' Gradient Norm :', NormGrad, '   relative grad :', &
                  100*NormGrad/inversion_param%Norm_grad_init, '%'
             write(INVERSE_LOG_FILE,*) ' Data std :', inversion_param%data_std
          endif



          !! stopping criteria -----------------------------------------------------
          !! sufficient decrease of cost function
          if ( (Q0/inversion_param%Cost_init) <= inversion_param%relat_cost) then
             if (myrank == 0) then
                write(INVERSE_LOG_FILE,*)
                write(INVERSE_LOG_FILE,*) &
                     ' FWI STOP : cost function reached maximum realative decrease allowed  :'&
                     , inversion_param%relat_cost
                write(INVERSE_LOG_FILE,*)
                write(OUTPUT_ITERATION_FILE,*)
                write(OUTPUT_ITERATION_FILE,*) &
                     ' FWI STOP : cost function reached maximum  realative decrease allowed  :'&
                     ,inversion_param%relat_cost
                write(OUTPUT_ITERATION_FILE,*)
             endif
             finished=.true.
          endif
          !! sufficient decrease of gradient
          if ( (NormGrad/inversion_param%Norm_grad_init) <= inversion_param%relat_grad) then
             if (myrank == 0) then
                write(INVERSE_LOG_FILE,*)
                write(INVERSE_LOG_FILE,*) &
                     ' FWI STOP : gradient of cost function reached maximum  realative decrease allowed  :'&
                     ,inversion_param%relat_grad
                write(INVERSE_LOG_FILE,*)
                write(OUTPUT_ITERATION_FILE,*)
                write(OUTPUT_ITERATION_FILE,*) &
                     ' FWI STOP : gradient of cost function reached maximum  realative decrease allowed  :'&
                     ,inversion_param%relat_grad
                write(OUTPUT_ITERATION_FILE,*)
             endif
             finished=.true.

          endif
          !!--------------------------------------------------------------------

          return !! return to the main iteration loop of FWI

       endif

    enddo

  end subroutine OneIterationOptim

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-------------------------------------------------------------------------------------------------------------
! initialize optimization scheme iterations
!-------------------------------------------------------------------------------------------------------------

  subroutine InitializeOptimIteration(iter_frq, acqui_simu, inversion_param)

    implicit none

    type(acqui),  dimension(:), allocatable,        intent(inout) :: acqui_simu
    type(inver),                                    intent(inout) :: inversion_param
    integer,                                        intent(in)    :: iter_frq
    integer                                                       :: ievent,  iter_inverse, ipar
    character(len=MAX_LEN_STRING)                                 :: prefix_name
    real(kind=CUSTOM_REAL)                                        :: tmp_val

    iter_inverse=0
    inversion_param%current_ifrq=iter_frq

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '          *********************************************'
       write(INVERSE_LOG_FILE,*) '          ***      INITIALIZE FWI ITERATIONS        ***'
       write(INVERSE_LOG_FILE,*) '          *********************************************'
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)  '          frequncy band :',iter_frq,  acqui_simu(1)%fl_event(iter_frq), &
            acqui_simu(1)%fh_event(iter_frq)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! initialize regularization
    call SetUpRegularization(inversion_param, acqui_simu, myrank)
    ! compute volume of domain
    regularization_penalty(:,:,:,:,:)=1._CUSTOM_REAL
    call Parallel_ComputeL2normSquare(regularization_penalty, inversion_param%NinvPar, inversion_param%volume_domain)
    inversion_param%volume_domain = inversion_param%volume_domain / inversion_param%NinvPar

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '      Volume of domain : ',  inversion_param%volume_domain
       write(INVERSE_LOG_FILE,*)
       call flush_iunit(INVERSE_LOG_FILE)
    endif


    !! compute normalize coefficient Tikonov weight
    inversion_param%weight_Tikonov = 1. / inversion_param%volume_domain

    ! compute cost function and gradient---------
    call InitForOneStepFWI(inversion_param)

    do ievent=1,acqui_simu(1)%nevent_tot
       call ComputeGradientPerEvent(ievent, iter_inverse, acqui_simu, inversion_param)
    enddo
    if (GPU_MODE) call TransfertKernelFromGPUArrays()

    ! store current value of cost function for reduction
    Q0=inversion_param%total_current_cost
    ! communicate gradient and cost function to all simultaneous runs
    call mpi_sum_grad_all_to_all_simultaneous_runs(Q0)

    !! compute standard deviation on data
    tmp_val =  inversion_param%nb_data_std
    call sum_all_all_cr_for_simulatenous_runs(tmp_val, inversion_param%nb_data_std,1)
    tmp_val =  inversion_param%data_std
    call sum_all_all_cr_for_simulatenous_runs(tmp_val, inversion_param%data_std,1)
    inversion_param%data_std = sqrt(inversion_param%data_std / inversion_param%nb_data_std)

    ! set model from modeling as reference model and
    if (iter_frq == 1) call Modeling2RefInvert(inversion_param, ref_model)

    ! store initial model in choosen family parameter
    call SpecfemParam2Invert(inversion_param, initial_model, ref_model)
    ! store initial gradient in choosen family parameter
    call StoreGradientInFamilyParam(inversion_param, initial_gradient, initial_model, ref_model, hess_approxim)

    if (inversion_param%input_SEM_prior) then
       ! save specific read prior model in family
       write(*,*) ' STOPPING : not use this option for now : input sem prior'
       stop
       call  SpecfemPrior2Invert(inversion_param, prior_model)
    else
       ! save starting model read as prior model for the first frequency group
       if (iter_frq == 1) then
          prior_model(:,:,:,:,:) = initial_model(:,:,:,:,:)
       endif
    endif

    !! compute reference mean model
    if (inversion_param%use_damping_SEM_Tikonov .and. inversion_param%input_SEM_prior) then
       do ipar =1, inversion_param%NinvPar
            inversion_param%damp_weight(ipar) = inversion_param%damp_weight(ipar) * &
                 (size(prior_model(:,:,:,:,ipar)) / sum(prior_model(:,:,:,:,ipar)))**2
       enddo
    endif

    ! compute regularization term and gradient for choosen family
    call AddRegularization(inversion_param, initial_model, ref_model, prior_model, regularization_penalty, &
         gradient_regularization_penalty, &
         myrank)
    ! add penalty tern on cost function
    call Parallel_ComputeL2normSquare(regularization_penalty, inversion_param%NinvPar, inversion_param%cost_penalty)
    Q0 = Q0 + 0.5 *  inversion_param%weight_Tikonov * inversion_param%cost_penalty
    initial_gradient(:,:,:,:,:) = initial_gradient(:,:,:,:,:) + &
         inversion_param%weight_Tikonov * gradient_regularization_penalty(:,:,:,:,:)

   ! store cost function value
    inversion_param%total_current_cost=Q0

    ! store intial values of cost function
    inversion_param%Cost_init=inversion_param%total_current_cost
    call Parallel_ComputeL2normSquare(initial_gradient, inversion_param%NinvPar, inversion_param%Norm_grad_init)

    ! define preconditionnner or taper on gradients
    call SetPrecond(iter_inverse, inversion_param, initial_gradient, hess_approxim, fwi_precond)

    ! write precond on disk
    if (mygroup <= 0 .and. VERBOSE_MODE) then
       prefix_name='precond'
       call DumpArray(fwi_precond, inversion_param, iter_inverse, prefix_name)
       prefix_name='Hess_app'
       call DumpArray(hess_approxim, inversion_param, iter_inverse, prefix_name)
       prefix_name='Grad_Regul'
       call DumpArray(gradient_regularization_penalty, inversion_param, iter_inverse, prefix_name)
       prefix_name='Regul'
       call DumpArray(regularization_penalty, inversion_param, iter_inverse, prefix_name)
       prefix_name='Spatial_damp'
       descent_direction(:,:,:,:,1)=spatial_damping(:,:,:,:) !! use as temporary working array
       call DumpArray(descent_direction, inversion_param, iter_inverse, prefix_name)
    endif

    ! store new model and gradient in l-bfgs history for the choosen family parameter
    call StoreModelAndGradientForLBFGS(initial_model, initial_gradient, 0)

    if (myrank == 0) then
       write(OUTPUT_FWI_LOG,*)
       write(OUTPUT_FWI_LOG,*) '  FREQUENCY GROUP :',  acqui_simu(1)%fl_event(iter_frq), &
            acqui_simu(1)%fh_event(iter_frq)
       write(OUTPUT_FWI_LOG,*)
       write(OUTPUT_ITERATION_FILE,*)
       write(OUTPUT_ITERATION_FILE,*) '  FREQUENCY GROUP :',  acqui_simu(1)%fl_event(iter_frq), &
            acqui_simu(1)%fh_event(iter_frq)
       write(OUTPUT_ITERATION_FILE,*)
       write(OUTPUT_ITERATION_FILE,'(i5,"|",e15.8,"|",e15.8,"|",e12.5,"|",e12.5,"|")')  &
            0, inversion_param%Cost_init, inversion_param%Norm_grad_init, 1., 1.
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) ' Initial Cost function : ', inversion_param%Cost_init
       write(INVERSE_LOG_FILE,*) ' Initial Gradient Norm :', inversion_param%Norm_grad_init
       write(INVERSE_LOG_FILE,*) ' Initial Data std :', inversion_param%data_std
       call flush_iunit(INVERSE_LOG_FILE)
       call flush_iunit(OUTPUT_ITERATION_FILE)
    endif

  end subroutine InitializeOptimIteration

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------
! compute current model form initial model and descent direction for a given step length
!-------------------------------------------------------------------------------------------------------------
!
! Since Specfem compute kernel for log of parameter, here we use this convention by default
! BE CAAREFUL : the model is stored in parameter and *NOT* log parameter in spesfem mesh (rhostore, kappastore, ...)
!
  subroutine UpdateModel(inversion_param, step_length, ModelIsSuitable)

    type(inver),                                    intent(inout) :: inversion_param
    logical,                                        intent(inout) :: ModelIsSuitable
    real(kind=CUSTOM_REAL),                         intent(in)    :: step_length
    integer                                                       :: ipar
    !integer                                                       :: igll, jgll, kgll, ispec

    real(kind=CUSTOM_REAL)                                        :: vmin, vmin_glob
    real(kind=CUSTOM_REAL)                                        :: vmax, vmax_glob
    real(kind=CUSTOM_REAL)                                        :: vmin_glob0
    real(kind=CUSTOM_REAL)                                        :: vmax_glob0

    current_model(:,:,:,:,:) = initial_model(:,:,:,:,:) + step_length * descent_direction(:,:,:,:,:)

    !! store the model on specfem arrays to perform next simulation
    call InvertParam2Specfem(inversion_param, current_model, ref_model)

    call CheckModelSuitabilityForModeling(ModelIsSuitable)

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) '    - > update model :  '
       !write(INVERSE_LOG_FILE,*)
       if (.not. ModelIsSuitable) then
           write(INVERSE_LOG_FILE,*)
           write(INVERSE_LOG_FILE,*) '    - > updated model not suitable for simulation trying smaller step  '
           write(INVERSE_LOG_FILE,*)
        endif
    endif

    do ipar=1, inversion_param%NinvPar

       select case(inversion_param%parameter_metric)

       case(0)  !! directly the parameter P
          vmin =   minval(current_model(:,:,:,:,ipar))
          vmax  =  maxval(current_model(:,:,:,:,ipar))
          call min_all_cr(vmin,vmin_glob)
          call max_all_cr(vmax,vmax_glob)

          vmin =   maxval( abs(current_model(:,:,:,:,ipar) - initial_model(:,:,:,:,ipar)) /   initial_model(:,:,:,:,ipar) )
          vmax  =  maxval( abs(current_model(:,:,:,:,ipar) - prior_model(:,:,:,:,ipar)) /   prior_model(:,:,:,:,ipar) )
          call max_all_cr(vmin,vmin_glob0)
          call max_all_cr(vmax,vmax_glob0)

       case(1)
          vmin =   minval(current_model(:,:,:,:,ipar)*ref_model(:,:,:,:,ipar))
          vmax  =  maxval(current_model(:,:,:,:,ipar)*ref_model(:,:,:,:,ipar))
          call min_all_cr(vmin,vmin_glob)
          call max_all_cr(vmax,vmax_glob)

          vmin =   maxval( abs(current_model(:,:,:,:,ipar) - initial_model(:,:,:,:,ipar)) /   initial_model(:,:,:,:,ipar) )
          vmax  =  maxval( abs(current_model(:,:,:,:,ipar) - prior_model(:,:,:,:,ipar)) /   prior_model(:,:,:,:,ipar) )
          call max_all_cr(vmin,vmin_glob0)
          call max_all_cr(vmax,vmax_glob0)

       case(2)
          vmin =   minval(exp(current_model(:,:,:,:,ipar)))
          vmax  =  maxval(exp(current_model(:,:,:,:,ipar)))
          call min_all_cr(vmin,vmin_glob)
          call max_all_cr(vmax,vmax_glob)

          vmin =   maxval( abs(exp(current_model(:,:,:,:,ipar)) - exp(initial_model(:,:,:,:,ipar))) /  &
               exp(initial_model(:,:,:,:,ipar)) )
          vmax  =  maxval( abs(exp(current_model(:,:,:,:,ipar)) - exp(prior_model(:,:,:,:,ipar))) /  &
               exp(prior_model(:,:,:,:,ipar)) )
          call max_all_cr(vmin,vmin_glob0)
          call max_all_cr(vmax,vmax_glob0)

       case(3)
          !! ??

       end select


       if (myrank == 0) then

          write(INVERSE_LOG_FILE,'( a13, i2, a10, 2(a8, f12.5), a19, f10.6,  a14, f10.6, a3)') &
               '  Parameter :', ipar, inversion_param%param_inv_name(ipar),'   MIN :',vmin_glob, &
               '   MAX :',vmax_glob, &
               ' max pert,  prior :', 100*vmax_glob0, ' % previous  :', 100*vmin_glob0,' %'

          !write(INVERSE_LOG_FILE,*) '            max pert / starting model : ', 100*vmax_glob0,' %'
          !write(INVERSE_LOG_FILE,*) '            max pert / previous model : ', 100*vmin_glob0,' %'

          call flush_iunit(INVERSE_LOG_FILE)
       endif


    enddo
    if (myrank == 0) write(INVERSE_LOG_FILE,*)

  end subroutine UpdateModel
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------
! compute initial guess for step length to try for line search
!-------------------------------------------------------------------------------------------------------------

  subroutine InitialGuessStep(inversion_param, step_length)

    type(inver),                                    intent(in)    :: inversion_param
    real(kind=CUSTOM_REAL),                         intent(inout) :: step_length
    real(kind=CUSTOM_REAL)                                        :: max_val_tmp, max_val

    max_val_tmp = maxval(abs(descent_direction(:,:,:,:,:)))
    call  max_all_all_cr(max_val_tmp, max_val)

    select case(inversion_param%parameter_metric)
    case(0)  !! directly the parameter P
       step_length=inversion_param%max_relative_pert/max_val
    case(1)  !! P / Pref
       step_length=inversion_param%max_relative_pert/max_val
    case(2)  !! log(P)
       step_length= log(1. + inversion_param%max_relative_pert) / max_val
    case(3) !! log(P/Pref)
       step_length= log(1. + inversion_param%max_relative_pert) / max_val
    end select

    if (DEBUG_MODE) then
       write(IIDD,* ) 'STEP  :',  step_length
    endif


  end subroutine InitialGuessStep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------
! Prepare and allocate all arrays used in inverse problem
!-------------------------------------------------------------------------------------------------------------
  subroutine AllocatememoryForFWI(inversion_param, nevent)

    type(inver),                                    intent(inout) :: inversion_param
    integer,                                        intent(in)    :: nevent
    integer                                                       :: ier, Ninvpar

    !! 1/ set the family parameter
    call PrepareArraysfamilyParam(inversion_param)

    Ninvpar =inversion_param%NinvPar


    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*) '  '
       write(INVERSE_LOG_FILE,*) '  allocate arrays for fwi iterations,  Ninvpar:', Ninvpar
       write(INVERSE_LOG_FILE,*) '  '
       write(INVERSE_LOG_FILE,*) '*************************************************************************'
       write(INVERSE_LOG_FILE,*) '  '
       call flush_iunit(INVERSE_LOG_FILE)
    endif

    !! allocate arrays for inversion scheme
    call AllocateArraysForInversion(inversion_param)

    !! allocate arrays for fwi_iteration
    allocate(initial_model(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 544')
    if (ier /= 0) call exit_MPI(myrank,"error allocation initial_model in AllocatememoryForFWI  subroutine")

    allocate(prior_model(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 545')
    if (ier /= 0) call exit_MPI(myrank,"error allocation prior_model in AllocatememoryForFWI  subroutine")

    allocate(ref_model(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 546')
    if (ier /= 0) call exit_MPI(myrank,"error allocation ref_model in AllocatememoryForFWI  subroutine")

    allocate(current_model(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 547')
    if (ier /= 0) call exit_MPI(myrank,"error allocation current_model in AllocatememoryForFWI  subroutine")

    allocate(initial_gradient(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 548')
    if (ier /= 0) call exit_MPI(myrank,"error allocation initial_gradient in AllocatememoryForFWI  subroutine")

    allocate(current_gradient(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 549')
    if (ier /= 0) call exit_MPI(myrank,"error allocation current_gradient in AllocatememoryForFWI  subroutine")

    allocate(descent_direction(NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, Ninvpar),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 550')
    if (ier /= 0) call exit_MPI(myrank,"error allocation descent_direction in AllocatememoryForFWI  subroutine")

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


    allocate(inversion_param%current_cost_prime(NEVENT),  inversion_param%previous_cost_prime(NEVENT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 555')
    allocate(inversion_param%current_cost(NEVENT),  inversion_param%previous_cost(NEVENT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 556')

  end subroutine AllocatememoryForFWI


end module fwi_iteration
