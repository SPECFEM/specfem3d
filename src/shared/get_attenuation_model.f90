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

  module attenuation_model

  implicit none

  ! model_attenuation_storage_var
  type model_attenuation_storage_var
    sequence
    double precision, dimension(:,:), pointer :: tau_eps_storage
    double precision, dimension(:), pointer :: Q_storage
    integer :: Q_resolution
    integer :: Q_max
  end type model_attenuation_storage_var
  type (model_attenuation_storage_var) :: AM_S

  ! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision :: Q  ! Q     = Desired Value of Attenuation or Q
    double precision :: iQ ! iQ    = 1/Q
    double precision, dimension(:), allocatable ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), allocatable :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer :: nf          ! nf    = Number of Frequencies
    integer :: nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables
  type(attenuation_simplex_variables) :: AS_V

  end module attenuation_model

!
!------------------------------------------------------------------------
!

  subroutine get_attenuation_model_olsen(vs_val,Q_mu,OLSEN_ATTENUATION_RATIO)

! uses scaling rule similar to Olsen et al. (2003) to determine attenuation medium
!
! returns: selected (sediment) Q_mu
!
! refers to:
!   K. B. Olsen, S. M. Day and C. R. Bradley, 2003.
!   Estimation of Q for Long-Period (>2 sec) Waves in the Los Angeles Basin
!   BSSA, 93, 2, p. 627-638

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: vs_val
  double precision,intent(out) :: Q_mu

  double precision,intent(in) :: OLSEN_ATTENUATION_RATIO

  !local parameters
  integer :: int_Q_mu

  ! two variations of scaling rule handling
  logical,parameter :: USE_SIMPLE_OLSEN = .false.
  logical,parameter :: USE_DISCRETE_OLSEN = .true.

  ! uses rule Q_mu = constant * v_s
  ! v_s in m/s
  Q_mu = OLSEN_ATTENUATION_RATIO * vs_val

  ! scaling rule
  if (USE_SIMPLE_OLSEN) then
    ! uses a simple, 2-constant model mentioned in Olsen et al. (2003)
    ! vs (in m/s)
    if (vs_val < 2000.0_CUSTOM_REAL) then
      Q_mu = 0.02 * vs_val
    else
      Q_mu = 0.1 * vs_val
    endif
  else if (USE_DISCRETE_OLSEN) then
    ! uses discrete values in sediment range
    int_Q_mu = 10 * nint(Q_mu / 10.)

    ! limits Q to sediment range values
    if (int_Q_mu < 40) int_Q_mu = 40
    if (int_Q_mu > 150) int_Q_mu = 150

    ! converts to double precision value
    Q_mu = dble(int_Q_mu)
  else
    stop 'Error no Olsen model specified, please set rule in get_attenuation_model.f90'
  endif

  end subroutine get_attenuation_model_olsen

!
!------------------------------------------------------------------------
!

  subroutine get_attenuation_model_olsen_qkappa(vp_val,vs_val,Q_mu,Q_kappa, &
                                                USE_ANDERSON_CRITERIA,SCALING_FACTOR_QP_FROM_QS)

! scales Q_kappa from Q_mu and vs/vp ratio

  use constants, only: CUSTOM_REAL,ATTENUATION_COMP_MAXIMUM

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: vp_val,vs_val
  double precision, intent(in) :: Q_mu
  double precision, intent(out) :: Q_kappa

  logical, intent(in) :: USE_ANDERSON_CRITERIA
  double precision, intent(in) :: SCALING_FACTOR_QP_FROM_QS

  ! local parameters
  double precision :: L_val,Q_s,Q_p

  ! Anderson & Hart (1978), Q of the Earth, JGR, 83, No. B12
  ! conversion between (Qp,Qs) and (Qkappa,Qmu)
  ! factor L
  L_val = 4.0d0/3.d0 * (vs_val/vp_val)**2

  ! attenuation Qs (eq.1)
  Q_s = Q_mu

  ! scales Qp from Qs (scaling factor introduced by Zhinan?)
  ! todo: should we scale Qkappa directly? e.g. Q_kappa = 10.d0 * Q_mu
  Q_p = SCALING_FACTOR_QP_FROM_QS * Q_s

  ! Anderson & Hart criteria: Qs/Qp >= L (eq. 4) since Qmu and Qkappa must be positive
  if (USE_ANDERSON_CRITERIA) then
    ! enforces (eq. 4) from Anderson & Hart
    ! note: this might lead to Q_p < Q_s
    if ((Q_s - L_val * Q_p) <= 0.d0 ) then
      ! negligible bulk attenuation (1/Q_kappa -> zero)
      Q_kappa = ATTENUATION_COMP_MAXIMUM
    else
      ! converts to bulk attenuation (eq. 3)
      Q_kappa = (1.0d0 - L_val) * Q_p * Q_s / (Q_s - L_val * Q_p)
    endif
  else
    ! note: this case might lead to: Q_kappa < Q_mu
    !       thus having a solid with stronger bulk attenuation than shear attenuation?

    ! avoids division by zero
    if (abs(Q_s - L_val * Q_p) <= 1.d-5 ) then
      Q_kappa = ATTENUATION_COMP_MAXIMUM
    else
      ! converts to bulk attenuation (eq. 3)
      Q_kappa = (1.0d0 - L_val) * Q_p * Q_s / (Q_s - L_val * Q_p)
    endif
  endif

  end subroutine get_attenuation_model_olsen_qkappa

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_model(nspec,USE_OLSEN_ATTENUATION,OLSEN_ATTENUATION_RATIO, &
                                   mustore,rho_vs,kappastore,rho_vp,qkappa_attenuation_store,qmu_attenuation_store, &
                                   ispec_is_elastic,min_resolved_period,prname,ATTENUATION_f0_REFERENCE)

! precalculates attenuation arrays and stores arrays into files

  use constants, only: myrank,NGLLX,NGLLY,NGLLZ,N_SLS, &
    CUSTOM_REAL,MAX_STRING_LEN,HUGEVAL,IMAIN, &
    ATTENUATION_COMP_MAXIMUM

  use shared_parameters, only: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,COMPUTE_FREQ_BAND_AUTOMATIC,ATT_F_C_SOURCE

  implicit none

  double precision,intent(in) :: OLSEN_ATTENUATION_RATIO,ATTENUATION_f0_REFERENCE
  integer,intent(in) :: nspec
  logical,intent(in) :: USE_OLSEN_ATTENUATION

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: mustore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: rho_vs
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: qkappa_attenuation_store
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: qmu_attenuation_store

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: kappastore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: rho_vp

  logical, dimension(nspec),intent(in) :: ispec_is_elastic
  real(kind=CUSTOM_REAL),intent(in) :: min_resolved_period
  character(len=MAX_STRING_LEN),intent(in) :: prname

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: factor_common
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: scale_factor, scale_factor_kappa
  double precision, dimension(N_SLS) :: tau_sigma_dble,beta_dble,beta_dble_kappa
  double precision factor_scale_dble,one_minus_sum_beta_dble, &
                   factor_scale_dble_kappa,one_minus_sum_beta_dble_kappa
  double precision :: Q_mu,Q_kappa
  double precision :: f_c_source
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tau_sigma
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tauinv
  real(kind=CUSTOM_REAL):: vs_val,vp_val
  integer :: i,j,k,ispec,ier
  double precision :: qmin,qmax,qmin_all,qmax_all
  double precision :: qmin_kappa,qmax_kappa,qmin_kappa_all,qmax_kappa_all

  !-----------------------------------------------------
  ! user parameter

  ! enforces ratio Qs/Qp >= L factor from Anderson & Hart (1978)
  ! IMPORTANT: this flag applies only if USE_OLSEN_ATTENUATION is true
  logical, parameter :: USE_ANDERSON_CRITERIA = .true.

  ! scaling factor to scale Qp from Qs
  double precision,parameter :: SCALING_FACTOR_QP_FROM_QS = 1.5d0

  !-----------------------------------------------------

  ! initializes arrays
  allocate(factor_common(N_SLS,NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1182')
  allocate(scale_factor(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1183')
  if (ier /= 0) call exit_mpi(myrank,'error allocation attenuation arrays')

  allocate(factor_common_kappa(N_SLS,NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1184')
  allocate(scale_factor_kappa(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1185')
  if (ier /= 0) call exit_mpi(myrank,'error allocation attenuation arrays')

  factor_common(:,:,:,:,:) = 1._CUSTOM_REAL
  scale_factor(:,:,:,:) = 1._CUSTOM_REAL

  factor_common_kappa(:,:,:,:,:) = 1._CUSTOM_REAL
  scale_factor_kappa(:,:,:,:) = 1._CUSTOM_REAL

  f_c_source = 0.d0

  ! gets stress relaxation times tau_sigma, i.e.
  ! precalculates tau_sigma depending on period band (constant for all Q_mu), and
  ! determines central frequency f_c_source of attenuation period band
  call get_attenuation_constants(min_resolved_period,tau_sigma_dble, &
                                 f_c_source,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! checks
  if (f_c_source <= 0.d0) call exit_MPI(myrank,'Error: invalid attenuation center frequency, cannot be zero or negative')

  ! stores center frequency as shared parameter
  ATT_F_C_SOURCE = f_c_source

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "Attenuation:"
    write(IMAIN,*) "  The code uses a constant Q quality factor, but approximated"
    write(IMAIN,*) "  based on a series of Zener standard linear solids (SLS)."
    write(IMAIN,*)
    write(IMAIN,*) "  Approximation is performed in the following frequency band:"
    write(IMAIN,*) "  Reference frequency requested by the user (Hz):",sngl(ATTENUATION_f0_REFERENCE), &
                                                        " period (s):",sngl(1.0/ATTENUATION_f0_REFERENCE)
    if (COMPUTE_FREQ_BAND_AUTOMATIC) then
      write(IMAIN,*)
      write(IMAIN,*) "  The following values are computed automatically by the code"
      write(IMAIN,*) "  based on the estimated maximum frequency resolution of your mesh"
      write(IMAIN,*) "  and can thus vary from what you have requested."
    endif

    write(IMAIN,*)
    write(IMAIN,*) "  Frequency band        min/max (Hz):",sngl(1.0/MAX_ATTENUATION_PERIOD),sngl(1.0/MIN_ATTENUATION_PERIOD)
    write(IMAIN,*) "  Period band           min/max (s) :",sngl(MIN_ATTENUATION_PERIOD),sngl(MAX_ATTENUATION_PERIOD)
    write(IMAIN,*) "  Logarithmic central frequency (Hz):",sngl(ATT_F_C_SOURCE)," period (s):",sngl(1.0/ATT_F_C_SOURCE)
    write(IMAIN,*)
    write(IMAIN,*) "  Using full attenuation with both Q_kappa and Q_mu."

    if (USE_OLSEN_ATTENUATION) then
      write(IMAIN,*) "  Using Olsen scaling with attenuation ratio Qmu/vs = ",sngl(OLSEN_ATTENUATION_RATIO)
      if (USE_ANDERSON_CRITERIA) write(IMAIN,*) "  Using Anderson and Hart criteria for ratio Qs/Qp"
    endif

    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! determines inverse of tau_sigma
  tau_sigma(:) = real(tau_sigma_dble(:),kind=CUSTOM_REAL)

  ! precalculates the inverse of tau_sigma
  tauinv(:) = 1._CUSTOM_REAL / tau_sigma(:)

  ! precalculates factors for shear modulus scaling according to attenuation model
  qmin = HUGEVAL
  qmax = 0.0
  qmin_kappa = HUGEVAL
  qmax_kappa = 0.0

  do ispec = 1,nspec

    ! skips non elastic elements
    if (ispec_is_elastic(ispec) .eqv. .false.) cycle

    ! determines attenuation factors for each GLL point
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          ! initializes Q
          Q_mu = 0.d0
          Q_kappa = 0.d0

          ! shear moduli attenuation
          ! gets Q_mu value
          if (USE_OLSEN_ATTENUATION) then
            ! shear attenuation
            ! use scaling rule similar to Olsen et al. (2003)
            vs_val = mustore(i,j,k,ispec) / rho_vs(i,j,k,ispec)
            call get_attenuation_model_olsen(vs_val,Q_mu,OLSEN_ATTENUATION_RATIO)
          else
            ! takes Q set in (CUBIT) mesh
            Q_mu = qmu_attenuation_store(i,j,k,ispec)
          endif

          ! attenuation zero (skips point if zero value)
          if (Q_mu <= 1.e-5) cycle

          ! limits Q
          if (Q_mu < 1.0d0) Q_mu = 1.0d0
          if (Q_mu > ATTENUATION_COMP_MAXIMUM) Q_mu = ATTENUATION_COMP_MAXIMUM

          ! statistics on Q_mu
          if (Q_mu < qmin) qmin = Q_mu
          if (Q_mu > qmax) qmax = Q_mu

          ! bulk moduli attenuation
          ! gets Q_kappa value

          if (USE_OLSEN_ATTENUATION) then
            ! bulk attenuation
            ! compressional wave speed vp
            vp_val = (kappastore(i,j,k,ispec) + 2.0d0 * mustore(i,j,k,ispec) / 3.0d0) / rho_vp(i,j,k,ispec)
            call get_attenuation_model_olsen_qkappa(vp_val,vs_val,Q_mu,Q_kappa, &
                                                    USE_ANDERSON_CRITERIA,SCALING_FACTOR_QP_FROM_QS)
          else
            ! takes Q set in (CUBIT) mesh
            Q_kappa = qkappa_attenuation_store(i,j,k,ispec)
          endif

          ! attenuation zero (means negligible attenuation)
          if (Q_kappa <= 1.e-5) Q_kappa = ATTENUATION_COMP_MAXIMUM

          ! limits Q
          if (Q_kappa < 1.0d0) Q_kappa = 1.0d0
          if (Q_kappa > ATTENUATION_COMP_MAXIMUM) Q_kappa = ATTENUATION_COMP_MAXIMUM

          ! statistics on Q_kappa
          if (Q_kappa < qmin_kappa) qmin_kappa = Q_kappa
          if (Q_kappa > qmax_kappa) qmax_kappa = Q_kappa

          ! gets beta, on_minus_sum_beta and factor_scale
          ! based on calculation of strain relaxation times tau_eps
          call get_attenuation_factors(Q_mu,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
                                       ATT_F_C_SOURCE,tau_sigma_dble, &
                                       beta_dble,one_minus_sum_beta_dble,factor_scale_dble, &
                                       Q_kappa,beta_dble_kappa,one_minus_sum_beta_dble_kappa,factor_scale_dble_kappa, &
                                       ATTENUATION_f0_REFERENCE)

          ! shear attenuation

          ! stores factor for runge-kutta scheme
          ! using factor for modulus defect Delta M_i = - M_relaxed
          ! see e.g. Savage et al. (BSSA, 2010): eq. 11
          !     precomputes factor: 2 ( 1 - tau_eps_i / tau_sigma_i ) / tau_sigma_i
          !EB EB May 2018 : this expression has been corrected, replaced by :
          ! 2 (1 - tau_eps_i / tau_sigma_i ) / tau_sigma_i) / sum(tau_eps / tau_sigma)
          factor_common(:,i,j,k,ispec) = (2._CUSTOM_REAL * beta_dble(:) * tauinv(:)) / one_minus_sum_beta_dble
          scale_factor(i,j,k,ispec) = factor_scale_dble

          ! bulk attenuation
          factor_common_kappa(:,i,j,k,ispec) = (beta_dble_kappa(:) * tauinv(:)) / one_minus_sum_beta_dble_kappa
          scale_factor_kappa(i,j,k,ispec) = factor_scale_dble_kappa

        enddo
      enddo
    enddo
  enddo

  ! statistics
  call min_all_dp(qmin,qmin_all)
  call max_all_dp(qmax,qmax_all)
  call min_all_dp(qmin_kappa,qmin_kappa_all)
  call max_all_dp(qmax_kappa,qmax_kappa_all)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  Q_mu min/max           : ",sngl(qmin_all),sngl(qmax_all)
    write(IMAIN,*) "  Q_kappa min/max        : ",sngl(qmin_kappa_all),sngl(qmax_kappa_all)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! stores attenuation arrays into files
  open(unit=27, file=prname(1:len_trim(prname))//'attenuation.bin', &
        status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error: could not open ',prname(1:len_trim(prname))//'attenuation.bin'
    call exit_mpi(myrank,'error opening attenuation.bin file')
  endif
  write(27) nspec

  ! shear attenuation
  write(27) factor_common
  write(27) scale_factor

  ! bulk attenuation
  write(27) factor_common_kappa
  write(27) scale_factor_kappa

  close(27)

  ! frees memory
  deallocate(factor_common,scale_factor)
  deallocate(factor_common_kappa,scale_factor_kappa)

  end subroutine get_attenuation_model

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_attenuation_memory_values(tau_s, deltat, alphaval,betaval,gammaval)

! returns: runge-kutta scheme terms alphaval, betaval and gammaval

  use constants

  implicit none

  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(in) :: tau_s
  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(out) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL),intent(in) :: deltat

  ! local parameter
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tauinv

  ! inverse of tau_s
  tauinv(:) = - 1._CUSTOM_REAL / tau_s(:)

  ! runge-kutta coefficients
  ! see e.g.: Savage et al. (BSSA, 2010): eq. (11)
  alphaval(:) = 1.0 + deltat*tauinv(:) + deltat**2 * tauinv(:)**2 / 2._CUSTOM_REAL &
                    + deltat**3 * tauinv(:)**3 / 6._CUSTOM_REAL &
                    + deltat**4 * tauinv(:)**4 / 24._CUSTOM_REAL
  betaval(:) = deltat / 2._CUSTOM_REAL + deltat**2 * tauinv(:) / 3._CUSTOM_REAL &
                   + deltat**3 * tauinv(:)**2 / 8._CUSTOM_REAL &
                   + deltat**4 * tauinv(:)**3 / 24._CUSTOM_REAL
  gammaval(:) = deltat / 2._CUSTOM_REAL + deltat**2 * tauinv(:) / 6._CUSTOM_REAL &
                    + deltat**3 * tauinv(:)**2 / 24._CUSTOM_REAL

  end subroutine get_attenuation_memory_values


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_constants(min_resolved_period,tau_sigma, &
                                       f_c_source,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

! returns: period band constants tau_sigma and center frequency f_c_source

  use constants

  use shared_parameters, only: COMPUTE_FREQ_BAND_AUTOMATIC

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: min_resolved_period
  double precision, dimension(N_SLS),intent(inout) :: tau_sigma
  double precision,intent(inout) :: f_c_source
  double precision,intent(inout) :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD

  ! local parameters
  real(kind=CUSTOM_REAL)  :: min_period

  ! determines min/max periods for attenuation band based on minimum resolved period of mesh
  if (COMPUTE_FREQ_BAND_AUTOMATIC) then ! otherwise they were entered as input values by the user in the Par_file
    min_period = 0.99 * min_resolved_period ! uses a small margin
    call get_attenuation_periods(min_period,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)
  endif

  !  sets up stress relaxation times tau_sigma,
  ! equally spaced based on number of standard linear solids and period band
  call get_attenuation_tau_sigma(tau_sigma,N_SLS,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! sets up central frequency
  ! logarithmic mean of frequency interval of absorption band
  call get_attenuation_source_freq(f_c_source,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  end subroutine get_attenuation_constants


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_factors(Q_mu,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
                                     f_c_source,tau_sigma, &
                                     beta,one_minus_sum_beta,factor_scale, &
                                     Q_kappa,beta_kappa,one_minus_sum_beta_kappa,factor_scale_kappa,ATTENUATION_f0_REFERENCE)

! returns: attenuation mechanisms beta,one_minus_sum_beta,factor_scale

! variable frequency range
! variable period range
! variable central logarithmic frequency

! in the future when more memory is available on computers
! it would be more accurate to use four mechanisms instead of three

  use constants

  implicit none

  double precision,intent(in) :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,ATTENUATION_f0_REFERENCE
  double precision,intent(in) :: f_c_source,Q_mu,Q_kappa
  double precision, dimension(N_SLS) :: tau_sigma
  double precision, dimension(N_SLS) :: beta,beta_kappa
  double precision :: one_minus_sum_beta,one_minus_sum_beta_kappa
  double precision :: factor_scale,factor_scale_kappa

  ! local parameters
  double precision, dimension(N_SLS) :: tau_eps,tau_eps_kappa

  ! determines tau_eps for Q_kappa
  call get_attenuation_tau_eps(Q_kappa,tau_sigma,tau_eps_kappa,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! determines one_minus_sum_beta
  call get_attenuation_property_values(tau_sigma,tau_eps_kappa,beta_kappa,one_minus_sum_beta_kappa)

  ! determines the "scale factor"
  call get_attenuation_scale_factor(f_c_source,tau_eps_kappa,tau_sigma,Q_kappa,factor_scale_kappa,ATTENUATION_f0_REFERENCE)
  ! uncomment this to print the constants to use in the 3D viscoelastic analytical code for validation purposes
  ! if (myrank == 0) &
  !   print *,'for Q_Kappa,tau_eps_kappa,tau_sigma,factor_scale_kappa = ',Q_Kappa,tau_eps_kappa(:),tau_sigma(:),factor_scale_kappa

  ! determines tau_eps for Q_mu
  call get_attenuation_tau_eps(Q_mu,tau_sigma,tau_eps,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! determines one_minus_sum_beta
  call get_attenuation_property_values(tau_sigma,tau_eps,beta,one_minus_sum_beta)

  ! determines the "scale factor"
  call get_attenuation_scale_factor(f_c_source,tau_eps,tau_sigma,Q_mu,factor_scale,ATTENUATION_f0_REFERENCE)
  ! uncomment this to print the constants to use in the 3D viscoelastic analytical code for validation purposes
  ! if (myrank == 0) &
  !   print *,'for Q_mu,tau_eps,tau_sigma,factor_scale = ',Q_mu,tau_eps(:),tau_sigma(:),factor_scale

  end subroutine get_attenuation_factors

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_attenuation_property_values(tau_s, tau_eps, beta, one_minus_sum_beta)

! coefficients useful for calculation between relaxed and unrelaxed moduli
!
! returns: coefficients beta, one_minus_sum_beta

  use constants

  implicit none

  double precision,dimension(N_SLS),intent(in) :: tau_s, tau_eps
  double precision,dimension(N_SLS),intent(out) :: beta
  double precision,intent(out):: one_minus_sum_beta

  integer :: i

  ! see e.g. Komatitsch & Tromp 1999, eq. (7)
  ! EB EB May 2018 this equation was wrong and has been corrected here
  ! coefficients beta
  beta(:) = tau_eps(:) / tau_s(:)

  ! sum of coefficients beta
  one_minus_sum_beta = ZERO
  do i = 1,N_SLS
    one_minus_sum_beta = one_minus_sum_beta + beta(i)
    ! this factor will be used later to get the modulus defect
    beta(i) = beta(i)-ONE
  enddo

  end subroutine get_attenuation_property_values

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_scale_factor(f_c_source,tau_eps,tau_sigma,Q_val,scale_factor,ATTENUATION_f0_REFERENCE)

! returns: physical dispersion scaling factor scale_factor

  use constants

  implicit none

  double precision,intent(in) :: Q_val, f_c_source
  ! strain and stress relaxation times
  double precision, dimension(N_SLS),intent(in) :: tau_eps, tau_sigma
  double precision, intent(in) :: ATTENUATION_f0_REFERENCE

  double precision,intent(out) :: scale_factor

  ! local parameters
  double precision :: w_c_source
  double precision :: factor_scale_mu0, factor_scale_mu
  double precision :: xtmp1_nu1,xtmp2_nu1,xtmp_ak_nu1
  integer :: i

  !--- compute central angular frequency of source (non dimensionalized)
  w_c_source = TWO_PI * f_c_source

  !--- quantity by which to scale mu_0 to get mu

  ! this formula can be found for instance in
  ! Liu, H. P., Anderson, D. L. and Kanamori, H., Velocity dispersion due to
  ! anelasticity: implications for seismology and mantle composition,
  ! Geophys. J. R. Astron. Soc., vol. 47, pp. 41-58 (1976)
  ! and in Aki, K. and Richards, P. G., Quantitative seismology, theory and methods,
  ! W. H. Freeman, (1980), second edition, sections 5.5 and 5.5.2, eq. (5.81) p. 170
  factor_scale_mu0 = ONE + TWO * log(f_c_source / ATTENUATION_f0_REFERENCE ) / (PI * Q_val)

  !--- quantity by which to scale mu to get mu_unrelaxed
  xtmp1_nu1 = ONE
  xtmp2_nu1 = ONE

  do i = 1,N_SLS
     xtmp_ak_nu1 = tau_eps(i)/tau_sigma(i) - ONE
     xtmp1_nu1 = xtmp1_nu1 + xtmp_ak_nu1/N_SLS
     xtmp2_nu1 = xtmp2_nu1 + xtmp_ak_nu1/(ONE + ONE/(TWO * PI * f_c_source * tau_sigma(i))**2)/N_SLS
  enddo

  factor_scale_mu = xtmp1_nu1/xtmp2_nu1

  !--- total factor by which to scale mu0 to get mu_unrelaxed
  scale_factor = factor_scale_mu * factor_scale_mu0

  !--- check that the correction factor is close to one
  if (scale_factor < 0.5 .or. scale_factor > 1.5) then
    print *,"Error : in get_attenuation_scale_factor() "
    print *,"  scale factor: ", scale_factor, " should be between 0.5 and 1.5"
    print *,"  factor scale_mu = ",factor_scale_mu," factor scale_mu0 = ",factor_scale_mu0
    print *,"  Q value = ", Q_val, " central frequency = ",f_c_source
    print *,"  ATTENUATION_f0_REFERENCE = ",ATTENUATION_f0_REFERENCE
    print *,"  please check your reference frequency ATTENUATION_f0_REFERENCE"
    call exit_MPI(myrank,'unreliable correction factor in attenuation model')
  endif

  end subroutine get_attenuation_scale_factor

!
!-------------------------------------------------------------------------------------------------
!

!compare: auto_ner.f90, GLOBE package

  subroutine get_attenuation_periods(min_resolved_period,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

! determines min/max periods for attenuation based upon mininum resolved period of mesh

  use constants

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: min_resolved_period
  double precision,intent(out) :: MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD

  ! local parameters
  double precision :: THETA(5)

  ! checks number of standard linear solids
  if (N_SLS < 2 .or. N_SLS > 5) then
     stop 'N_SLS must be greater than 1 or less than 6'
  endif

  ! THETA defines the width of the Attenation Range in Decades
  !   The number defined here were determined by minimizing
  !   the "flatness" of the absoption spectrum.  Each THETA
  !   is defined for a particular N_SLS (constants.h)
  !   THETA(2) is for N_SLS = 2
  THETA(1)           =   0.00d0
  THETA(2)           =   0.75d0
  THETA(3)           =   1.75d0
  THETA(4)           =   2.25d0
  THETA(5)           =   2.85d0
  ! Compute Min Attenuation Period
  !
  ! The Minimum attenuation period = (Grid Spacing in km) / V_min
  !  Grid spacing in km     = Width of an element in km * spacing for GLL point * points per wavelength

  MIN_ATTENUATION_PERIOD = min_resolved_period


  ! Compute Max Attenuation Period
  !
  ! The max attenuation period for 3 SLS is optimally
  !   1.75 decades from the min attenuation period, see THETA above
  !
  ! this uses: theta = log( T_max / T_min ) to calculate T_max for a given T_min

  MAX_ATTENUATION_PERIOD = MIN_ATTENUATION_PERIOD * 10.0d0**THETA(N_SLS)

  end subroutine get_attenuation_periods

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_tau_sigma(tau_s, nsls, min_period, max_period)

! Determines stress relaxation times tau_sigma
! Sets the Tau_sigma (tau_s) to be equally spaced in log10 frequency

  implicit none

  integer :: nsls
  double precision,intent(in) :: min_period, max_period
  double precision,intent(out) :: tau_s(nsls)
  ! local parameters
  double precision :: f1, f2
  double precision :: exp1, exp2
  double precision :: dexpval
  integer :: i
  double precision, parameter :: PI = 3.14159265358979d0

  ! min/max frequencies
  f1 = 1.0d0 / max_period
  f2 = 1.0d0 / min_period

  ! logarithms
  exp1 = log10(f1)
  exp2 = log10(f2)

  ! equally spaced in log10 frequency
  dexpval = (exp2-exp1) / ((nsls*1.0d0) - 1)
  do i = 1,nsls
    tau_s(i) = 1.0 / (PI * 2.0d0 * 10**(exp1 + (i - 1)* 1.0d0 *dexpval))
  enddo

  end subroutine get_attenuation_tau_sigma

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_source_freq(f_c_source,min_period,max_period)

! determines the source frequency

  implicit none

  double precision,intent(out) :: f_c_source
  double precision,intent(in) :: min_period, max_period

  ! local parameters
  double precision :: f1, f2

  ! min/max frequencies
  f1 = 1.0d0 / max_period
  f2 = 1.0d0 / min_period

  ! use the logarithmic central frequency
  f_c_source = 10.0d0**(0.5d0 * (log10(f1) + log10(f2)))

  end subroutine get_attenuation_source_freq

!--------------------------------------------------------------------------------------------------
!
!  This portion of the SPECFEM3D Code was written by:
!  Brian Savage while at
!     California Institute of Technology
!     Department of Terrestrial Magnetism / Carnegie Institute of Washington
!     Univeristy of Rhode Island
!
!   It is based upon formulation in the following references:
!
!   Dahlen and Tromp, 1998
!      Theoretical Global Seismology
!
!   Liu et al. 1976
!      Velocity dispersion due to anelasticity: implications for seismology and mantle composition
!      Geophys, J. R. asts. Soc, Vol 47, pp. 41-58
!
!   The methodology can be found in:
!       Savage, B, D. Komatitsch and J. Tromp, 2010.
!       Effects of 3D Attenuation on Seismic Wave Amplitude and Phase Measurements
!       BSSA, 100, 3, p. 1241-1251.
!
! modifications:
!  - minor modifications by Daniel Peter, november 2010
!--------------------------------------------------------------------------------------------------

  subroutine get_attenuation_tau_eps(Q_in,tau_s,tau_eps,min_period,max_period)

! includes min_period, max_period, and N_SLS
!
! returns: determines strain relaxation times tau_eps

  use constants

  implicit none

  double precision :: Q_in
  double precision, dimension(N_SLS) :: tau_s, tau_eps
  double precision :: min_period,max_period

  ! local parameters
  integer :: rw

  ! note: to speed up this attenuation routine, we will try to compute the attenuation factors only for new Q values.
  !       often, Q values are given for simple Q-models, thus there is no need to recompute the same factors for every GLL point.
  !       we use a storage table AM_S to check/retrieve and store computed tau_eps values for specific Q values.
  !
  ! tries first to READ from storage array
  rw = 1
  call model_attenuation_storage(Q_in, tau_eps, rw)

  ! checks if value was found
  if (rw > 0) return

  ! new Q value, computes tau factors
  call attenuation_invert_by_simplex(min_period, max_period, N_SLS, Q_in, tau_s, tau_eps)

  ! WRITE into storage array, to keep in case for next GLL points
  rw = -1
  call model_attenuation_storage(Q_in, tau_eps, rw)

  end subroutine get_attenuation_tau_eps

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_attenuation_storage(Q_in, tau_eps, rw)

  use constants

  use attenuation_model, only: AM_S

  implicit none

  double precision,intent(in) :: Q_in
  double precision, dimension(N_SLS),intent(inout) :: tau_eps
  integer,intent(inout) :: rw

  ! local parameters
  integer :: Qtmp
  integer :: ier
  !double precision :: Qnew

  double precision, parameter :: ZERO_TOL = 1.e-5

  integer, save :: first_time_called = 1

  ! allocates arrays when first called
  if (first_time_called == 1) then
    first_time_called = 0
    AM_S%Q_resolution = 10**ATTENUATION_COMP_RESOLUTION
    AM_S%Q_max = ATTENUATION_COMP_MAXIMUM
    Qtmp = AM_S%Q_resolution * AM_S%Q_max

    allocate(AM_S%tau_eps_storage(N_SLS, Qtmp), &
             AM_S%Q_storage(Qtmp),stat=ier)
    if (ier /= 0) stop 'error allocating arrays for attenuation storage'
    AM_S%Q_storage(:) = -1
  endif

  if (Q_in < 0.0d0 .or. Q_in > AM_S%Q_max) then
    print *,'Error attenuation_storage()'
    print *,'Attenuation Value out of Range: ', Q_in
    print *,'Attenuation Value out of Range: Min, Max ', 0, AM_S%Q_max
    stop 'Attenuation Value out of Range'
  endif

  ! check for zero Q value
  if (rw > 0 .and. Q_in <= ZERO_TOL) then
    !Q_in = 0.0d0;
    tau_eps(:) = 0.0d0;
    return
  endif

  ! Generate index for Storage Array
  ! and Recast Q using this index
  ! According to Brian, use float
  !Qtmp = Q_in * Q_resolution
  !Q_in = Qtmp / Q_resolution;

  ! by default: resolution is Q_resolution = 10
  ! converts Q to an array integer index:
  ! e.g. Q = 150.31 -> Qtmp = 150.31 * 10 = int( 1503.10 ) = 1503
  Qtmp = int(Q_in * dble(AM_S%Q_resolution))

  ! rounds to corresponding double value:
  ! e.g. Qnew = dble( 1503 ) / dble(10) = 150.30
  ! but Qnew is not used any further...
  !Qnew = dble(Qtmp) / dble(AM_S%Q_resolution)

  if (rw > 0) then
    ! checks
    if (first_time_called == 0) then
      if (.not. associated(AM_S%Q_storage)) &
        stop 'error calling model_attenuation_storage() routine without AM_S array'
    else
      stop 'error calling model_attenuation_storage() routine with first_time_called value invalid'
    endif

    ! READ
    if (AM_S%Q_storage(Qtmp) > 0) then
      ! READ SUCCESSFUL
      tau_eps(:) = AM_S%tau_eps_storage(:,Qtmp)
      ! corresponding Q value would be:
      !Q_in = AM_S%Q_storage(Qtmp)
      rw = 1
    else
      ! READ NOT SUCCESSFUL
      rw = -1
    endif
  else
    ! WRITE SUCCESSFUL
    AM_S%tau_eps_storage(:,Qtmp) = tau_eps(:)
    AM_S%Q_storage(Qtmp) = Q_in
    rw = 1
  endif

  end subroutine model_attenuation_storage

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_invert_by_simplex(t2, t1, n, Q_real, tau_s, tau_eps)

  implicit none

  ! Input / Output
  double precision,intent(in) :: t1, t2
  double precision,intent(in) :: Q_real
  integer,intent(in) :: n
  double precision, dimension(n),intent(in)  :: tau_s
  double precision, dimension(n),intent(out) :: tau_eps

  ! Internal
  integer :: i, iterations, err,prnt
  double precision :: f1, f2, exp1,exp2, min_value !, dexpval

  integer, parameter :: nf = 100
  double precision, dimension(nf) :: f

  !double precision, parameter :: PI = 3.14159265358979d0
  double precision, external :: attenuation_eval

  ! Values to be passed into the simplex minimization routine
  iterations = -1
  min_value  = -1.0e-4
  err        = 0
  prnt       = 0

  !allocate(f(nf))

  ! Determine the min and max frequencies
  f1 = 1.0d0 / t1
  f2 = 1.0d0 / t2

  ! Determine the exponents of the frequencies
  exp1 = log10(f1)
  exp2 = log10(f2)

!  if (f2 < f1 .or. Q_real < 0.0d0 .or. n < 1) then
!     call exit_MPI(0, 'frequencies flipped or Q less than zero or N_SLS < 0')
!  endif

  ! Determine the Source frequency
!  omega_not =  1.0e+03 * 10.0d0**(0.5 * (log10(f1) + log10(f2)))

  ! Determine the Frequencies at which to compare solutions
  !   The frequencies should be equally spaced in log10 frequency
  do i = 1,nf
    f(i) = exp1 + ((i-1)*1.0d0 * (exp2-exp1) / ((nf-1)*1.0d0))
  enddo

  ! Set the Tau_sigma (tau_s) to be equally spaced in log10 frequency
!  dexpval = (exp2-exp1) / ((n*1.0d0) - 1)
!  do i = 1,n
!     tau_s(i) = 1.0 / (PI * 2.0d0 * 10**(exp1 + (i - 1)* 1.0d0 *dexpval))
!  enddo


  ! Shove the parameters into the module
  call attenuation_simplex_setup(nf,n,f,Q_real,tau_s)

  ! Set the Tau_epsilon (tau_eps) to an initial value at omega*tau = 1
  ! tan_delta = 1/Q = (tau_eps - tau_s)/(2 * sqrt(tau e*tau_s))
  !    if we assume tau_eps =~ tau_s
  !    we get the equation below
  do i = 1,n
    tau_eps(i) = tau_s(i) + (tau_s(i) * 2.0d0/Q_real)
  enddo

  ! Run a simplex search to determine the optimum values of tau_eps
  call fminsearch(attenuation_eval, tau_eps, n, iterations, min_value, prnt, err)
  if (err > 0) then
    write(*,*)'Search did not converge for an attenuation of ', Q_real
    write(*,*)'    Iterations: ', iterations
    write(*,*)'    Min Value:  ', min_value
    write(*,*)'    Aborting program'
    call exit_MPI(0,'attenuation_simplex: Search for Strain relaxation times did not converge')
  endif

  !deallocate(f)

  call attenuation_simplex_finish()

  end subroutine attenuation_invert_by_simplex

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_simplex_setup(nf_in,nsls_in,f_in,Q_in,tau_s_in)

!   - Inserts necessary parameters into the module attenuation_simplex_variables
!   - See module for explaination

  use attenuation_model, only: AS_V

  implicit none

  integer,intent(in) :: nf_in, nsls_in
  double precision,intent(in) :: Q_in
  double precision, dimension(nf_in),intent(in)   :: f_in
  double precision, dimension(nsls_in),intent(in) :: tau_s_in

  ! local parameters
  integer :: ier

  allocate(AS_V%f(nf_in),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1188')
  allocate(AS_V%tau_s(nsls_in),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1189')
  if (ier /= 0) stop 'error allocating arrays for attenuation simplex'

  AS_V%nf    = nf_in
  AS_V%nsls  = nsls_in
  AS_V%f     = f_in
  AS_V%Q     = Q_in
  AS_V%iQ    = 1.0d0/AS_V%Q
  AS_V%tau_s = tau_s_in

  end subroutine attenuation_simplex_setup

!
!-------------------------------------------------------------------------------------------------
!

  double precision function attenuation_eval(Xin)

!    - Computes the misfit from a set of relaxation paramters
!          given a set of frequencies and target attenuation
!    - Evaluates only at the given frequencies
!    - Evaluation is done with an L2 norm
!
!    Input
!      Xin = Tau_epsilon, Strain Relaxation Time
!                Note: Tau_sigma the Stress Relaxation Time is loaded
!                      with attenuation_simplex_setup and stored in
!                      attenuation_simplex_variables
!
!    Xi = Sum_i^N sqrt [ (1/Qc_i - 1/Qt_i)^2 / 1/Qt_i^2 ]
!
!     where Qc_i is the computed attenuation at a specific frequency
!           Qt_i is the desired attenuaiton at that frequency
!
!    Uses attenuation_simplex_variables to store constant values
!
!    See attenuation_simplex_setup

  use attenuation_model, only: AS_V

  implicit none

   ! Input
  double precision, dimension(AS_V%nsls),intent(in) :: Xin

  ! local parameters
  double precision, dimension(AS_V%nsls) :: tau_eps
  double precision, dimension(AS_V%nf)   :: A, B, tan_delta

  integer :: i
  double precision :: xi, iQ2

  tau_eps = Xin

  call attenuation_maxwell(AS_V%nf,AS_V%nsls,AS_V%f,AS_V%tau_s,tau_eps,B,A)

  tan_delta = B / A

  attenuation_eval = 0.0d0
  iQ2 = AS_V%iQ**2
  do i = 1,AS_V%nf
    xi = sqrt(( ( (tan_delta(i) - AS_V%iQ) ** 2) / iQ2 ))
    attenuation_eval = attenuation_eval + xi
  enddo

  end function attenuation_eval

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_maxwell(nf,nsls,f,tau_s,tau_eps,B,A)

!   - Computes the Moduli (Maxwell Solid) for a series of
!         Standard Linear Solids
!   - Computes M1 and M2 parameters after Dahlen and Tromp pp.203
!         here called B and A after Liu et al. 1976
!   - Another formulation uses Kelvin-Voigt Solids and computes
!         Compliences J1 and J2 after Dahlen and Tromp pp.203
!
!   Input
!     nf    = Number of Frequencies
!     nsls  = Number of Standard Linear Solids
!     f     = Frequencies (in log10 of frequencies)
!                dimension(nf)
!     tau_s = Tau_sigma  Stress relaxation time (see References)
!                dimension(nsls)
!     tau_eps = Tau_epislon Strain relaxation time (see References)
!                dimension(nsls)!
!   Output
!     B     = Real Moduli      ( M2 Dahlen and Tromp pp.203 )
!                dimension(nf)
!     A     = Imaginary Moduli ( M1 Dahlen and Tromp pp.203 )
!                dimension(nf)
!
!   Dahlen and Tromp, 1998
!      Theoretical Global Seismology
!
!   Liu et al. 1976
!      Velocity dispersion due to anelasticity: implications for seismology and mantle composition
!      Geophys, J. R. asts. Soc, Vol 47, pp. 41-58

!! EB EB May 2018: the expressions of A and B from Liu have been corrected here,
!to incorporate the 1/L
  implicit none

  ! Input
  integer,intent(in) :: nf, nsls
  double precision, dimension(nf),intent(in)   :: f
  double precision, dimension(nsls),intent(in) :: tau_s, tau_eps
  ! Output
  double precision, dimension(nf),intent(out)   :: A,B

  integer :: i,j
  double precision :: w, denom
  double precision,parameter :: PI = 3.14159265358979d0

  A(:) = 0.0d0
  B(:) = 0.0d0

  do i = 1,nf
    w = 2.0d0 * PI * 10**f(i)
    do j = 1,nsls
      !        write(*,*)j,tau_s(j),tau_eps(j)
      denom = 1.0d0 + w**2 * tau_s(j)**2
      A(i) = A(i) + (1.0d0 + (w**2 * tau_eps(j) * tau_s(j)))/ denom
      B(i) = B(i) + w * ( tau_eps(j) - tau_s(j) ) / denom
    enddo
     !     write(*,*)A(i),B(i),10**f(i)
  enddo

  end subroutine attenuation_maxwell

!
!-------------------------------------------------------------------------------------------------
!


  subroutine fminsearch(funk, x, n, itercount, tolf, prnt, err)

! subroutine fminsearch
!   - Computes the minimization of funk(x(n)) using the simplex method
!   - This subroutine is copied from Matlab fminsearch.m
!         and modified to suit my nefarious needs
!   Input
!     funk = double precision function with one input parameter
!                double precision function the_funk(x)
!     x    = Input/Output
!               variables to be minimized
!               dimension(n)
!            Input:  Initial Value
!            Output: Mimimized Value
!     n    = number of variables
!     itercount = Input/Output
!                 Input:  maximum number of iterations
!                         if < 0 default is used (200 * n)
!                 Output: total number of iterations on output
!     tolf      = Input/Output
!                 Input:  minimium tolerance of the function funk(x)
!                 Output: minimium value of funk(x)(i.e. "a" solution)
!     prnt      = Input
!                 3 => report every iteration
!                 4 => report every iteration, total simplex
!     err       = Output
!                 0 => Normal exeecution, converged within desired range
!                 1 => function evaluation exceeded limit
!                 2 => Iterations exceeded limit
!
!     See Matlab fminsearch

  implicit none

  ! Input
  double precision, external :: funk

  integer,intent(in) :: n
  double precision,intent(inout) :: x(n) ! Also Output
  integer,intent(inout) :: itercount
  integer,intent(in) :: prnt
  integer,intent(out) :: err
  double precision,intent(inout) :: tolf

  ! local parameters
  !Internal
  integer :: i,j, how
  integer, parameter :: none             = 0
  integer, parameter :: initial          = 1
  integer, parameter :: expand           = 2
  integer, parameter :: reflect          = 3
  integer, parameter :: contract_outside = 4
  integer, parameter :: contract_inside  = 5
  integer, parameter :: shrink           = 6

  integer :: maxiter, maxfun
  integer :: func_evals
  double precision :: tolx

  double precision :: rho, chi, psi, sigma
  double precision :: xin(n), y(n), v(n,n+1), fv(n+1)
  double precision :: vtmp(n,n+1)
  double precision :: usual_delta, zero_term_delta
  double precision :: xbar(n), xr(n), fxr, xe(n), fxe, xc(n), fxc, fxcc, xcc(n)
  integer :: place(n+1)

  double precision :: max_size_simplex, max_value

  rho   = 1.0d0
  chi   = 2.0d0
  psi   = 0.5d0
  sigma = 0.5d0

  if (itercount > 0) then
     maxiter = itercount
  else
     maxiter = 200 * n
  endif
  itercount = 0
  maxfun  = 200 * n

  if (tolf > 0.0d0) then
     tolx = 1.0e-4
  else
     tolx = 1.0e-4
     tolf = 1.0e-4
  endif

  err = 0

  xin    = x
  v(:,:) = 0.0d0
  fv(:)  = 0.0d0

  v(:,1) = xin
  x      = xin

  fv(1) = funk(xin)

  usual_delta = 0.05
  zero_term_delta = 0.00025

  do j = 1,n
     y = xin
     if (y(j) /= 0.0d0) then
        y(j) = (1.0d0 + usual_delta) * y(j)
     else
        y(j) = zero_term_delta
     endif
     v(:,j+1) = y
     x(:) = y
     fv(j+1) = funk(x)
  enddo

  call qsort_local(fv,n+1,place)

  do i = 1,n+1
     vtmp(:,i) = v(:,place(i))
  enddo
  v = vtmp

  how = initial
  itercount = 1
  func_evals = n+1
  if (prnt == 3) then
     write(*,*)'Iterations   Funk Evals   Value How'
     write(*,*)itercount, func_evals, fv(1), how
  endif
  if (prnt == 4) then
     write(*,*)'How: ',how
     write(*,*)'V: ', v
     write(*,*)'fv: ',fv
     write(*,*)'evals: ',func_evals
  endif

  do while (func_evals < maxfun .and. itercount < maxiter)

     if (max_size_simplex(v,n) <= tolx .and. &
          max_value(fv,n+1) <= tolf) then
        goto 666
     endif
     how = none

     ! xbar = average of the n (NOT n+1) best points
     !     xbar = sum(v(:,1:n), 2)/n
     xbar(:) = 0.0d0
     do i = 1,n
        do j = 1,n
           xbar(i) = xbar(i) + v(i,j)
        enddo
        xbar(i) = xbar(i) / (n*1.0d0)
     enddo
     xr = (1 + rho)*xbar - rho*v(:,n+1)
     x(:) = xr
     fxr = funk(x)
     func_evals = func_evals + 1
     if (fxr < fv(1)) then
        ! Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,n+1)
        x = xe
        fxe = funk(x)
        func_evals = func_evals+1
        if (fxe < fxr) then
           v(:,n+1) = xe
           fv(n+1) = fxe
           how = expand
        else
           v(:,n+1) = xr
           fv(n+1) = fxr
           how = reflect
        endif
     else ! fv(:,1) <= fxr
        if (fxr < fv(n)) then
           v(:,n+1) = xr
           fv(n+1) = fxr
           how = reflect
        else ! fxr >= fv(:,n)
           ! Perform contraction
           if (fxr < fv(n+1)) then
              ! Perform an outside contraction
              xc = (1 + psi*rho)*xbar - psi*rho*v(:,n+1)
              x(:) = xc
              fxc = funk(x)
              func_evals = func_evals+1

              if (fxc <= fxr) then
                 v(:,n+1) = xc
                 fv(n+1) = fxc
                 how = contract_outside
              else
                 ! perform a shrink
                 how = shrink
              endif
           else
              ! Perform an inside contraction
              xcc = (1-psi)*xbar + psi*v(:,n+1)
              x(:) = xcc
              fxcc = funk(x)
              func_evals = func_evals+1

              if (fxcc < fv(n+1)) then
                 v(:,n+1) = xcc
                 fv(n+1) = fxcc
                 how = contract_inside
              else
                 ! perform a shrink
                 how = shrink
              endif
           endif
           if (how == shrink) then
              do j=2,n+1
                 v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1))
                 x(:) = v(:,j)
                 fv(j) = funk(x)
              enddo
              func_evals = func_evals + n
           endif
        endif
     endif

     call qsort_local(fv,n+1,place)
     do i = 1,n+1
        vtmp(:,i) = v(:,place(i))
     enddo
     v = vtmp

     itercount = itercount + 1
     if (prnt == 3) then
        write(*,*)itercount, func_evals, fv(1), how
     else if (prnt == 4) then
        write(*,*)
        write(*,*)'How: ',how
        write(*,*)'v: ',v
        write(*,*)'fv: ',fv
        write(*,*)'evals: ',func_evals
     endif
  enddo

  if (func_evals > maxfun) then
     write(*,*)'function evaluations exceeded prescribed limit', maxfun
     err = 1
  endif
  if (itercount > maxiter) then
     write(*,*)'iterations exceeded prescribed limit', maxiter
     err = 2
  endif

666 continue
  x = v(:,1)
  tolf = fv(1)

  end subroutine fminsearch

!
!-------------------------------------------------------------------------------------------------
!

  double precision function max_value(fv,n)

!    - Finds the maximim value of the difference of between the first
!          value and the remaining values of a vector
!    Input
!      fv = Input
!             Vector
!             dimension(n)
!      n  = Input
!             Length of fv
!
!      returns:
!         Xi = max( || fv(1)- fv(i) || ) for i=2:n
!

  implicit none
  integer,intent(in) :: n
  double precision,intent(in) :: fv(n)

  ! local parameters
  integer :: i
  double precision :: m, z

  m = 0.0d0
  do i = 2,n
     z = abs(fv(1) - fv(i))
     if (z > m) then
        m = z
     endif
  enddo

  max_value = m

  end function max_value

!
!-------------------------------------------------------------------------------------------------
!

  double precision function max_size_simplex(v,n)

!   - Determines the maximum distance between two point in a simplex
!   Input
!     v  = Input
!            Simplex Verticies
!            dimension(n, n+1)
!     n  = Pseudo Length of n
!
!     returns:
!       Xi = max( max( || v(:,1) - v(:,i) || ) ) for i=2:n+1
!

  implicit none
  integer,intent(in) :: n
  double precision,intent(in) :: v(n,n+1)

  ! local parameters
  integer :: i,j
  double precision :: m, z

  m = 0.0d0
  do i = 1,n
     do j = 2,n+1
        z = abs(v(i,j) - v(i,1))
        if (z > m) then
           m = z
        endif
     enddo
  enddo

  max_size_simplex = m

  end function max_size_simplex

!
!-------------------------------------------------------------------------------------------------
!


  subroutine qsort_local(X,n,I)

!    - Implementation of a Bubble Sort Routine
!    Input
!      X = Input/Output
!         Vector to be sorted
!         dimension(n)
!      n = Input
!         Length of X
!      I = Output
!         Sorted indices of vector X
!
!      Example:
!         X = [ 4 3 1 2 ] on Input
!         I = [ 1 2 3 4 ] Computed Internally (in order)
!
!         X = [ 1 2 3 4 ] on Output
!         I = [ 3 4 2 1 ] on Output
!

  implicit none

  integer,intent(in) :: n
  double precision,intent(inout) :: X(n)
  integer,intent(out) :: I(n)

  ! local parameters
  integer :: j,k
  double precision :: rtmp
  integer :: itmp

  do j = 1,n
     I(j) = j
  enddo

  do j = 1,n
     do k = 1,n-j
        if (X(k+1) < X(k)) then
           rtmp   = X(k)
           X(k)   = X(k+1)
           X(k+1) = rtmp

           itmp   = I(k)
           I(k)   = I(k+1)
           I(k+1) = itmp
        endif
     enddo
  enddo

  end subroutine qsort_local

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_simplex_finish()

  use attenuation_model, only: AS_V

  implicit none

  deallocate(AS_V%f)
  deallocate(AS_V%tau_s)

  end subroutine attenuation_simplex_finish

