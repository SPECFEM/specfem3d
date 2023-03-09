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

module regularization_interface

  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_AB

  use inverse_problem_par
  use regularization
  use regularization_fd_mod
  use inversion_scheme

  implicit none

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: spatial_damping

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine SetUpRegularization(inversion_param, acqui_simu, myrank)

      integer,                    intent(in)                                  :: myrank
      type(inver),                intent(inout)                               :: inversion_param
      type(acqui),  dimension(:), intent(inout)                               :: acqui_simu
      ! local
      integer                                                                 :: ier
      real(kind=CUSTOM_REAL)                                                  :: min_damp_val
      real(kind=CUSTOM_REAL)                                                  :: volume_domain
      real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable               :: array_mesh

      ! user output
      if (myrank == 0) then
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*) '*********************************************'
        write(INVERSE_LOG_FILE,*) '***         SET UP REGULARIZATION         ***'
        write(INVERSE_LOG_FILE,*) '*********************************************'
        write(INVERSE_LOG_FILE,*)
      endif

      ! compute volume of domain (for normalizing Tikhonov weight)
      allocate(array_mesh(NGLLX, NGLLY, NGLLZ, NSPEC_AB, inversion_param%Ninvpar),stat=ier)
      if (ier /= 0) call exit_MPI(myrank,"error allocation array_mesh for volume")
      array_mesh(:,:,:,:,:) = 1._CUSTOM_REAL

      call Parallel_ComputeL2normSquare(array_mesh, inversion_param%NinvPar, volume_domain)

      deallocate(array_mesh)

      ! volume (per parameter)
      volume_domain = volume_domain / inversion_param%NinvPar

      ! log output
      if (myrank == 0) then
        write(INVERSE_LOG_FILE,*) ' Volume of domain : ',  volume_domain
        write(INVERSE_LOG_FILE,*)
        call flush_iunit(INVERSE_LOG_FILE)
      endif

      !! compute normalize coefficient Tikhonov weight
      inversion_param%weight_Tikhonov_normalized = inversion_param%weight_Tikhonov / volume_domain

      ! SEM regularization
      if (inversion_param%use_regularization_SEM_Tikhonov) then
        ! user output
        if (myrank == 0) then
          write(INVERSE_LOG_FILE,*) ' TIKHONOV SEM BASED regularization :', inversion_param%weight_Tikhonov
          write(INVERSE_LOG_FILE,*) '                        normalized :', inversion_param%weight_Tikhonov_normalized
        endif
      endif

      ! damping
      if ( inversion_param%use_damping_SEM_Tikhonov &
          .or. inversion_param%use_variable_SEM_damping_sources &
          .or. inversion_param%use_variable_SEM_damping_stations) then
        ! user output
        if (myrank == 0) then
          write(INVERSE_LOG_FILE,*) ' TIKHONOV SEM BASED DAMPING :', inversion_param%weight_Tikhonov
          write(INVERSE_LOG_FILE,*) '                 normalized :', inversion_param%weight_Tikhonov_normalized
        endif

        ! allocates spatial damping array
        if (.not. allocated(spatial_damping)) then
          allocate(spatial_damping(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 211')
          spatial_damping(:,:,:,:) = 0._CUSTOM_REAL
        endif

        ! initializes with minimum damping value
        min_damp_val = 1._CUSTOM_REAL
        if (inversion_param%use_variable_SEM_damping_sources) min_damp_val = min(min_damp_val,inversion_param%min_damp_src)
        if (inversion_param%use_variable_SEM_damping_stations) min_damp_val = min(min_damp_val,inversion_param%min_damp_sta)
        spatial_damping(:,:,:,:) = min_damp_val

        ! damping source locations
        if (inversion_param%use_variable_SEM_damping_sources) then
          ! user output
          if (myrank == 0) then
            write(INVERSE_LOG_FILE,*) '   COMPUTE VARIABLE SOURCE DAMPING     :', &
                  inversion_param%min_damp_src, inversion_param%max_damp_src, inversion_param%distance_from_source
          endif

          ! adds source damping
          call compute_spatial_damping_for_source_singularities(acqui_simu, inversion_param, spatial_damping)
        endif

        ! damping station location
        if (inversion_param%use_variable_SEM_damping_stations) then
          ! user output
          if (myrank == 0) then
            write(INVERSE_LOG_FILE,*) '   COMPUTE VARIABLE STATION DAMPING    :', &
                  inversion_param%min_damp_sta, inversion_param%max_damp_sta, inversion_param%distance_from_station
          endif

          ! adds station damping
          call compute_spatial_damping_for_station_singularities(acqui_simu, inversion_param, spatial_damping)
        endif
      endif

      !! ----------- FD GRID BASED REGULARIZATION ----------
      if (inversion_param%use_regularization_FD_Tikhonov) then
        ! safety stop
        write(*,*) " ABORT :: BROKEN OPTION : FD Tikhonov "
        stop
        ! to do: needs testing...
        if (myrank == 0) then
          write(INVERSE_LOG_FILE,*) ' INITIALIZE TIKHONOV FD BASED regularization  :', inversion_param%weight_Tikhonov
        endif
        call setup_FD_regularization(inversion_param%projection_fd, myrank)
      endif

      ! user output
      if (myrank == 0) then
        write(INVERSE_LOG_FILE,*)
        call flush_iunit(INVERSE_LOG_FILE)
      endif

    end subroutine SetUpRegularization

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine AddRegularization(inversion_param, model, prior_model, regul_penalty, gradient_regul_penalty, myrank)

    type(inver),                                               intent(inout) :: inversion_param
    integer,                                                   intent(in)    :: myrank
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),              intent(in)    :: model, prior_model
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),              intent(inout) :: regul_penalty, gradient_regul_penalty
    ! local
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),   allocatable                :: model_on_sem, regul_on_sem, gradient_regul_on_sem
    integer                                                                  :: ipar, ier

    ! initializes
    inversion_param%cost_penalty      = 0._CUSTOM_REAL
    regul_penalty(:,:,:,:,:)          = 0._CUSTOM_REAL
    gradient_regul_penalty(:,:,:,:,:) = 0._CUSTOM_REAL

    !! ----------------- SEM BASED REGULARIZATION --------
    if (inversion_param%use_regularization_SEM_Tikhonov) then
      ! log output
      if (myrank == 0) then
        write(INVERSE_LOG_FILE,*) ' USE TIKHONOV SEM BASED regularization  : ', inversion_param%weight_Tikhonov
        write(INVERSE_LOG_FILE,*) '   weighting : ', inversion_param%smooth_weight(:)
        write(INVERSE_LOG_FILE,*)
        call flush_iunit(INVERSE_LOG_FILE)
      endif

      allocate(model_on_sem(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 212')
      model_on_sem(:,:,:,:) = 0._CUSTOM_REAL

      allocate(regul_on_sem(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 213')
      regul_on_sem(:,:,:,:) = 0._CUSTOM_REAL

      allocate(gradient_regul_on_sem(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 214')
      gradient_regul_on_sem(:,:,:,:) = 0._CUSTOM_REAL

      do ipar = 1, inversion_param%NinvPar
        model_on_sem(:,:,:,:) = model(:,:,:,:,ipar)

        call compute_bi_laplacian_of_field(model_on_sem, regul_on_sem, gradient_regul_on_sem)

        !! we add penalty on parameter not on log(parameter)
        regul_penalty(:,:,:,:,ipar) = regul_penalty(:,:,:,:,ipar) + &
                                      sqrt(inversion_param%smooth_weight(ipar)) * regul_on_sem(:,:,:,:)

        !! adds penalty on gradient
        gradient_regul_penalty(:,:,:,:,ipar) = gradient_regul_penalty(:,:,:,:,ipar) + &
                                               inversion_param%smooth_weight(ipar) * gradient_regul_on_sem(:,:,:,:)
      enddo

      deallocate(model_on_sem, regul_on_sem, gradient_regul_on_sem)
    endif

    ! damping
    if (inversion_param%use_damping_SEM_Tikhonov) then
      ! log output
      if (myrank == 0) then
        write(INVERSE_LOG_FILE,*) ' USE TIKHONOV SEM BASED DAMPING  : ', inversion_param%weight_Tikhonov
        write(INVERSE_LOG_FILE,*) '   damping parameters : ', inversion_param%damp_weight(:)
        ! source damping
        if (inversion_param%use_variable_SEM_damping_sources) then
          write(INVERSE_LOG_FILE,*) '   damping sources    : ', inversion_param%min_damp_src,'/',inversion_param%max_damp_src
        endif
        ! station damping
        if (inversion_param%use_variable_SEM_damping_stations) then
          write(INVERSE_LOG_FILE,*) '   damping stations   : ', inversion_param%min_damp_sta,'/',inversion_param%max_damp_sta
        endif
        write(INVERSE_LOG_FILE,*)
        call flush_iunit(INVERSE_LOG_FILE)
      endif

      do ipar = 1, inversion_param%NinvPar
        !! damp the parameter
        select case(inversion_param%parameter_metric)
        case(0,1)
          !! not log
          !! model penalty
          regul_penalty(:,:,:,:,ipar) = regul_penalty(:,:,:,:,ipar) + &
                  sqrt(inversion_param%damp_weight(ipar) * spatial_damping(:,:,:,:) ) * &
                  ( model(:,:,:,:,ipar) - prior_model(:,:,:,:,ipar) )

          !! gradient is with respect to the parameter
          gradient_regul_penalty(:,:,:,:,ipar) = gradient_regul_penalty(:,:,:,:,ipar) + &
                  inversion_param%damp_weight(ipar) * spatial_damping(:,:,:,:) * &
                  ( model(:,:,:,:,ipar) -  prior_model(:,:,:,:,ipar) )
        case(2,3)
          !! for log
          !! model penalty
          regul_penalty(:,:,:,:,ipar) = regul_penalty(:,:,:,:,ipar) + &
                  sqrt(inversion_param%damp_weight(ipar) * spatial_damping(:,:,:,:) ) * &
                  ( exp(model(:,:,:,:,ipar)) - exp(prior_model(:,:,:,:,ipar)) )

          !! gradient is with respect to the log parameter
          gradient_regul_penalty(:,:,:,:,ipar) = gradient_regul_penalty(:,:,:,:,ipar) + &
                  inversion_param%damp_weight(ipar) * spatial_damping(:,:,:,:) * &
                  ( exp(model(:,:,:,:,ipar)) -  exp(prior_model(:,:,:,:,ipar)) ) * exp(model(:,:,:,:,ipar))
        end select
      enddo

      !!TO DO : we can add damping on vp and vp/vs to avoid non physical (non numerical) models
      !! need to use family parameter selector
      !regul_penalty(:,:,:,:,2) = ( model(:,:,:,:,2) / model(:,:,:,:,3) - sqrt(3.))
      !gradient_regul_penalty(:,:,:,:,2) = gradient_regul_penalty(:,:,:,:,2) +  (model(:,:,:,:,2) / model(:,:,:,:,3))**2
      !gradient_regul_penalty(:,:,:,:,3) = gradient_regul_penalty(:,:,:,:,3) -  (model(:,:,:,:,2) / model(:,:,:,:,3))**2
    endif

    !! ----------- FD GRID BASED REGULARIZATION ----------
    if (inversion_param%use_regularization_FD_Tikhonov) then
      ! not working yet
      write(*,*) "ABORT :: BROKEN OPTION : FD Tikhonov "
      stop

      if (myrank == 0) then
        write(INVERSE_LOG_FILE,*) ' USE TIKHONOV FD BASED regularization  :', inversion_param%weight_Tikhonov
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

        call gradient_FD_laplac(model_on_sem, regul_on_sem, gradient_regul_on_sem, &
                                inversion_param%cost_penalty, &
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
