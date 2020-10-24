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

module regularization_interface

  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_AB

  use inverse_problem_par
  use regularization
  use regularization_fd_mod

  implicit none

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: spatial_damping

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine SetUpRegularization(inversion_param, acqui_simu, myrank)

      integer,                                        intent(in)    :: myrank
      type(inver),                                    intent(inout) :: inversion_param
      type(acqui),  dimension(:), allocatable,        intent(inout) :: acqui_simu

      integer :: ier

      if (myrank == 0) then
         write(INVERSE_LOG_FILE,*)
         write(INVERSE_LOG_FILE,*) '          *********************************************'
         write(INVERSE_LOG_FILE,*) '          ***         SET UP REGULARIZATION         ***'
         write(INVERSE_LOG_FILE,*) '          *********************************************'
         write(INVERSE_LOG_FILE,*)
      endif

      if (inversion_param%use_regularization_SEM_Tikonov) then
         if (myrank == 0) then
            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*) '   INITIALIZE TIKONOV SEM BASED regularization  :', inversion_param%weight_Tikonov
         endif
      endif

      if ( inversion_param%use_damping_SEM_Tikonov .or. inversion_param%use_variable_SEM_damping) then
         if (.not. allocated(spatial_damping)) then
           allocate(spatial_damping(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
           if (ier /= 0) call exit_MPI_without_rank('error allocating array 211')
         endif
         spatial_damping(:,:,:,:)=min(1._CUSTOM_REAL, inversion_param%min_damp)
      endif

      if (inversion_param%use_damping_SEM_Tikonov) then
         if (myrank == 0) then

            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*) '   INITIALIZE TIKONOV SEM BASED DAMPING  :', inversion_param%weight_Tikonov

         endif
      endif

      if (inversion_param%use_variable_SEM_damping) then
         if (myrank == 0) then
            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*) '   COMPUTE VARIABLE DAMPING  :', inversion_param%min_damp, inversion_param%max_damp, &
                 inversion_param%distance_from_source
         endif
         call compute_spatial_damping_for_source_singularities(acqui_simu, inversion_param, spatial_damping)
      endif


      !! ----------- FD GRID BASED REGULARIZATION ----------
      if (inversion_param%use_regularization_FD_Tikonov) then
         write(*,*) " ABORT :: BROKEN OPTION : FD tikonov "
         stop
         if (myrank == 0) then
            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*) '   INITIALIZE TIKONOV FD BASED regularization  :', inversion_param%weight_Tikonov
            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*)
         endif
         call setup_FD_regularization(inversion_param%projection_fd, myrank)
      endif


    end subroutine SetUpRegularization

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine AddRegularization(inversion_param, model, ref_model, prior_model, regul_penalty, gradient_regul_penalty, myrank)

    type(inver),                                               intent(inout) :: inversion_param
    integer,                                                   intent(in)    :: myrank
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable, intent(inout) :: model, ref_model, prior_model
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable, intent(inout) :: regul_penalty, gradient_regul_penalty
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),   allocatable                :: model_on_sem, regul_on_sem, gradient_regul_on_sem
    real(kind=CUSTOM_REAL)                                                   :: cost_penalty
    integer                                                                  :: ipar, ier

    inversion_param%cost_penalty= 0._CUSTOM_REAL
    regul_penalty(:,:,:,:,:)= 0._CUSTOM_REAL
    gradient_regul_penalty(:,:,:,:,:)= 0._CUSTOM_REAL

    !! ----------------- SEM BASED REGULARIZATION --------
    if (inversion_param%use_regularization_SEM_Tikonov) then
       if (myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '   USE TIKONOV SEM BASED regularization  :', inversion_param%weight_Tikonov
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
       endif

       allocate(model_on_sem(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 212')
       allocate(regul_on_sem(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 213')
       allocate(gradient_regul_on_sem(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 214')

       do ipar = 1, inversion_param%NinvPar

!!$          if (inversion_param%use_log) then
!!$             model_on_sem(:,:,:,:) = exp(model(:,:,:,:,ipar))
!!$          else
          model_on_sem(:,:,:,:) = model(:,:,:,:,ipar)
!!$          endif
          call compute_bi_laplacian_of_field(model_on_sem, regul_on_sem, gradient_regul_on_sem)

          !! we add penalty on parameter not on log(parameter)
          regul_penalty(:,:,:,:,ipar) = regul_penalty(:,:,:,:,ipar) + &
               sqrt(inversion_param%smooth_weight(ipar))*regul_on_sem(:,:,:,:)

          !! gradient with respect ot the log parameter but with penalty of parameter
          !! thus need to multiply by the model parameter
!!$          if (inversion_param%use_log) then
!!$             gradient_regul_penalty(:,:,:,:,ipar) = gradient_regul_penalty(:,:,:,:,ipar) + &
!!$                  inversion_param%smooth_weight(ipar)* &
!!$                  gradient_regul_on_sem(:,:,:,:)*exp(model(:,:,:,:,ipar))
!!$          else
          gradient_regul_penalty(:,:,:,:,ipar) = gradient_regul_penalty(:,:,:,:,ipar) + &
               inversion_param%smooth_weight(ipar)* &
               gradient_regul_on_sem(:,:,:,:)
!!$          endif

       enddo

       deallocate(model_on_sem, regul_on_sem, gradient_regul_on_sem)

    endif


    if (inversion_param%use_damping_SEM_Tikonov) then
       if (myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '   USE TIKONOV SEM BASED DAMPING  :', inversion_param%weight_Tikonov
          write(INVERSE_LOG_FILE,*)         inversion_param%damp_weight(:)
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
       endif
       do ipar =1, inversion_param%NinvPar

          !! directly damp the log
          !regul_penalty(:,:,:,:,ipar) =  log(model(:,:,:,:,ipar)) - log(prior_model(:,:,:,:,ipar))
          !gradient_regul_penalty(:,:,:,:,ipar) =  log(model(:,:,:,:,ipar)) -  log(prior_model(:,:,:,:,ipar))

!!$          if (inversion_param%use_log) then
!!$
!!$             !! damp the parameter itself (not log)
!!$             regul_penalty(:,:,:,:,ipar) =  regul_penalty(:,:,:,:,ipar) + &
!!$                  sqrt(inversion_param%damp_weight(ipar) * spatial_damping(:,:,:,:) ) * &
!!$                  (exp(model(:,:,:,:,ipar)) - exp(prior_model(:,:,:,:,ipar)))
!!$
!!$             !! gradient is with respect to the log parameter
!!$             gradient_regul_penalty(:,:,:,:,ipar) =   gradient_regul_penalty(:,:,:,:,ipar) + &
!!$                  inversion_param%damp_weight(ipar) * spatial_damping(:,:,:,:) * &
!!$                  ( exp(model(:,:,:,:,ipar)) -  exp(prior_model(:,:,:,:,ipar)) ) * exp(model(:,:,:,:,ipar))
!!$          else
             !! damp the parameter

          select case(inversion_param%parameter_metric)

          case(2) !! for log

             regul_penalty(:,:,:,:,ipar) =  regul_penalty(:,:,:,:,ipar) + &
                  sqrt(inversion_param%damp_weight(ipar) * spatial_damping(:,:,:,:) ) * &
                  (exp(model(:,:,:,:,ipar)) - exp(prior_model(:,:,:,:,ipar)))

             !! gradient is with respect to the log parameter
             gradient_regul_penalty(:,:,:,:,ipar) =   gradient_regul_penalty(:,:,:,:,ipar) + &
                  inversion_param%damp_weight(ipar) * spatial_damping(:,:,:,:) * &
                  ( exp(model(:,:,:,:,ipar)) -  exp(prior_model(:,:,:,:,ipar)) ) * exp(model(:,:,:,:,ipar))

          case default

             regul_penalty(:,:,:,:,ipar) =  regul_penalty(:,:,:,:,ipar) + &
                  sqrt(inversion_param%damp_weight(ipar) * spatial_damping(:,:,:,:) ) * &
                  (model(:,:,:,:,ipar) - prior_model(:,:,:,:,ipar))

             !! gradient is with respect to the parameter
             gradient_regul_penalty(:,:,:,:,ipar) =   gradient_regul_penalty(:,:,:,:,ipar) + &
                  inversion_param%damp_weight(ipar) * spatial_damping(:,:,:,:) * &
                  ( model(:,:,:,:,ipar) -  prior_model(:,:,:,:,ipar) )
          end select
!!$          endif
!!$
       enddo

       !!TO DO : we can add damping on vp and vp/vs to avoid non physical (non numerical) models
       !! need to us family parameter selector
       !regul_penalty(:,:,:,:,2) = ( model(:,:,:,:,2) / model(:,:,:,:,3) - sqrt(3.))
       !gradient_regul_penalty(:,:,:,:,2) = gradient_regul_penalty(:,:,:,:,2) +  (model(:,:,:,:,2) / model(:,:,:,:,3))**2
       !gradient_regul_penalty(:,:,:,:,3) = gradient_regul_penalty(:,:,:,:,3) -  (model(:,:,:,:,2) / model(:,:,:,:,3))**2

    endif


    !! ----------- FD GRID BASED REGULARIZATION ----------
    if (inversion_param%use_regularization_FD_Tikonov) then
       write(*,*) " ABORT :: BROKEN OPTION : FD tikonov "
       stop
       write(*,*) ref_model(1,1,1,1,1)
       if (myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '   USE TIKONOV FD BASED regularization  :', inversion_param%weight_Tikonov
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*)
       endif

       allocate(model_on_sem(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 215')
       allocate(regul_on_sem(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 216')
       allocate(gradient_regul_on_sem(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 217')

       do ipar = 1, inversion_param%NinvPar
          model_on_sem(:,:,:,:) = model(:,:,:,:,ipar)
          call gradient_FD_laplac(model_on_sem, regul_on_sem, gradient_regul_on_sem, cost_penalty, &
                                  inversion_param%projection_fd, myrank, ipar)
          regul_penalty(:,:,:,:,ipar) = regul_on_sem(:,:,:,:)
          gradient_regul_penalty(:,:,:,:,ipar) = gradient_regul_on_sem(:,:,:,:)

       enddo
       deallocate(model_on_sem, regul_on_sem, gradient_regul_on_sem)
    endif

  end subroutine AddRegularization

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module regularization_interface
