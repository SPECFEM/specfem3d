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
    double precision, dimension(:), pointer :: Qmu_storage
    integer Q_resolution
    integer Q_max
  end type model_attenuation_storage_var
  type (model_attenuation_storage_var) AM_S

  ! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables
  type(attenuation_simplex_variables) AS_V

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

  use constants

  implicit none

  real(kind=CUSTOM_REAL) :: vs_val
  double precision :: Q_mu
  double precision :: OLSEN_ATTENUATION_RATIO

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

  ! limits Q_mu value range
  if (Q_mu < 1.0d0) Q_mu = 1.0d0
  if (Q_mu > ATTENUATION_COMP_MAXIMUM) Q_mu = ATTENUATION_COMP_MAXIMUM


  end subroutine get_attenuation_model_olsen

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_model(myrank,nspec,USE_OLSEN_ATTENUATION,OLSEN_ATTENUATION_RATIO, &
                                  mustore,rho_vs,kappastore,rho_vp,qkappa_attenuation_store,qmu_attenuation_store, &
                                  ispec_is_elastic,min_resolved_period,prname,ATTENUATION_f0_REFERENCE)

! precalculates attenuation arrays and stores arrays into files

  use constants

  use shared_parameters, only: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,COMPUTE_FREQ_BAND_AUTOMATIC

  implicit none

  double precision,intent(in) :: OLSEN_ATTENUATION_RATIO,ATTENUATION_f0_REFERENCE
  integer,intent(in) :: myrank,nspec
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
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: factor_common
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: one_minus_sum_beta_kappa
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: scale_factor, scale_factor_kappa
  double precision, dimension(N_SLS) :: tau_sigma_dble,beta_dble,beta_dble_kappa
  double precision factor_scale_dble,one_minus_sum_beta_dble, &
                   factor_scale_dble_kappa,one_minus_sum_beta_dble_kappa
  double precision :: Q_mu,Q_kappa,Q_p,Q_s
  double precision :: L_val
  double precision :: f_c_source
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tau_sigma
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tauinv
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: beta,beta_kappa
  real(kind=CUSTOM_REAL):: vs_val,vp_val
  integer :: i,j,k,ispec,ier
  double precision :: qmin,qmax,qmin_all,qmax_all
  double precision :: qmin_kappa,qmax_kappa,qmin_kappa_all,qmax_kappa_all

  !-----------------------------------------------------
  ! user parameter

  ! enforces ratio Qs/Qp >= L factor from Anderson & Hart (1978)
  ! IMPORTANT: this flag applies only if USE_OLSEN_ATTENUATION is true
  logical, parameter :: USE_ANDERSON_CRITERIA = .true.

  !-----------------------------------------------------

  ! initializes arrays
  allocate(one_minus_sum_beta(NGLLX,NGLLY,NGLLZ,nspec), &
          factor_common(N_SLS,NGLLX,NGLLY,NGLLZ,nspec), &
          scale_factor(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_mpi(myrank,'error allocation attenuation arrays')

  allocate(one_minus_sum_beta_kappa(NGLLX,NGLLY,NGLLZ,nspec), &
          factor_common_kappa(N_SLS,NGLLX,NGLLY,NGLLZ,nspec), &
          scale_factor_kappa(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_mpi(myrank,'error allocation attenuation arrays')

  one_minus_sum_beta(:,:,:,:) = 1._CUSTOM_REAL
  factor_common(:,:,:,:,:) = 1._CUSTOM_REAL
  scale_factor(:,:,:,:) = 1._CUSTOM_REAL

  one_minus_sum_beta_kappa(:,:,:,:) = 1._CUSTOM_REAL
  factor_common_kappa(:,:,:,:,:) = 1._CUSTOM_REAL
  scale_factor_kappa(:,:,:,:) = 1._CUSTOM_REAL

  ! gets stress relaxation times tau_sigma, i.e.
  ! precalculates tau_sigma depending on period band (constant for all Q_mu), and
  ! determines central frequency f_c_source of attenuation period band
  call get_attenuation_constants(min_resolved_period,tau_sigma_dble, &
                                 f_c_source,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "Attenuation:"
    write(IMAIN,*) "The code uses a constant Q quality factor,"
    write(IMAIN,*) "but approximated based on a series of Zener standard linear solids (SLS)."
    write(IMAIN,*) "The approximation is performed in the following frequency band:"
    write(IMAIN,*) "  Reference frequency requested by the user (Hz):",sngl(ATTENUATION_f0_REFERENCE), &
                                            " period (s):",sngl(1.0/ATTENUATION_f0_REFERENCE)

    if (COMPUTE_FREQ_BAND_AUTOMATIC) write(IMAIN,*) "  The following values are computed automatically by the code based on &
         &the estimated maximum frequency resolution of your mesh and can thus vary from what you have requested:"

    write(IMAIN,*) "  Frequency band min/max (Hz):",sngl(1.0/MAX_ATTENUATION_PERIOD),sngl(1.0/MIN_ATTENUATION_PERIOD)
    write(IMAIN,*) "  Period band min/max (s):",sngl(MIN_ATTENUATION_PERIOD),sngl(MAX_ATTENUATION_PERIOD)
    write(IMAIN,*) "  Logarithmic central frequency (Hz):",sngl(f_c_source)," period (s):",sngl(1.0/f_c_source)
    write(IMAIN,*) "  Using full attenuation with both Q_kappa and Q_mu."

    if (USE_OLSEN_ATTENUATION) then
      write(IMAIN,*) "  using Olsen scaling with attenuation ratio Qmu/vs = ",sngl(OLSEN_ATTENUATION_RATIO)
      if (USE_ANDERSON_CRITERIA) write(IMAIN,*) "  using Anderson and Hart criteria for ratio Qs/Qp"
    endif

    call flush_IMAIN()
  endif

  ! determines inverse of tau_sigma
  tau_sigma(:) = real(tau_sigma_dble(:),kind=CUSTOM_REAL)

  ! precalculates the negative inverse of tau_sigma
  tauinv(:) = - 1._CUSTOM_REAL / tau_sigma(:)

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

              ! Anderson & Hart (1978), Q of the Earth, JGR, 83, No. B12
              ! conversion between (Qp,Qs) and (Qkappa,Qmu)
              ! factor L
              L_val = 4.0d0/3.d0 * (vs_val/vp_val)**2

              ! attenuation Qs (eq.1)
              Q_s = Q_mu

              ! scales Qp from Qs (scaling factor introduced by Zhinan?)
              ! todo: should we scale Qkappa directly? e.g. Q_kappa = 10.d0 * Q_mu
              Q_p = 1.5d0 * Q_s

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
          call get_attenuation_factors(myrank,Q_mu,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
                                       f_c_source,tau_sigma_dble, &
                                       beta_dble,one_minus_sum_beta_dble,factor_scale_dble, &
                                       Q_kappa,beta_dble_kappa,one_minus_sum_beta_dble_kappa,factor_scale_dble_kappa, &
                                       ATTENUATION_f0_REFERENCE)

          ! shear attenuation
          ! stores factor for unrelaxed parameter
          one_minus_sum_beta(i,j,k,ispec) = one_minus_sum_beta_dble

          ! stores factor for runge-kutta scheme
          ! using factor for modulus defect Delta M_i = - M_relaxed
          ! see e.g. Savage et al. (BSSA, 2010): eq. 11
          !     precomputes factor: 2 ( 1 - tau_eps_i / tau_sigma_i ) / tau_sigma_i
          beta(:) = beta_dble(:)
          factor_common(:,i,j,k,ispec) = 2._CUSTOM_REAL * beta(:) * tauinv(:)

          ! stores scale factor for mu moduli
          scale_factor(i,j,k,ispec) = factor_scale_dble

          ! bulk attenuation
          one_minus_sum_beta_kappa(i,j,k,ispec) = one_minus_sum_beta_dble_kappa
          beta_kappa(:) = beta_dble_kappa(:)
          factor_common_kappa(:,i,j,k,ispec) = beta_kappa(:) * tauinv(:)

          ! stores scale factor for mu moduli
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
  write(27) one_minus_sum_beta
  write(27) factor_common
  write(27) scale_factor

  ! bulk attenuation
  write(27) one_minus_sum_beta_kappa
  write(27) factor_common_kappa
  write(27) scale_factor_kappa

  close(27)

  ! frees memory
  deallocate(one_minus_sum_beta,factor_common,scale_factor)
  deallocate(one_minus_sum_beta_kappa,factor_common_kappa,scale_factor_kappa)

  end subroutine get_attenuation_model

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_attenuation_memory_values(tau_s, deltat, alphaval,betaval,gammaval)

! returns: runge-kutta scheme terms alphaval, betaval and gammaval

  use constants

  implicit none

  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tau_s, alphaval, betaval,gammaval
  real(kind=CUSTOM_REAL) :: deltat

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

  real(kind=CUSTOM_REAL) :: min_resolved_period
  double precision, dimension(N_SLS) :: tau_sigma
  double precision :: f_c_source
  double precision :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD

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

  subroutine get_attenuation_factors(myrank,Q_mu,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
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

  integer:: myrank
  double precision :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,ATTENUATION_f0_REFERENCE
  double precision :: f_c_source,Q_mu,Q_kappa
  double precision, dimension(N_SLS) :: tau_sigma
  double precision, dimension(N_SLS) :: beta,beta_kappa
  double precision :: one_minus_sum_beta,one_minus_sum_beta_kappa
  double precision :: factor_scale,factor_scale_kappa

  ! local parameters
  double precision, dimension(N_SLS) :: tau_eps,tau_eps_kappa


  ! determines tau_eps for Q_mu
  call get_attenuation_tau_eps(Q_mu,tau_sigma,tau_eps, &
                               MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! determines one_minus_sum_beta
  call get_attenuation_property_values(tau_sigma,tau_eps,beta,one_minus_sum_beta)

  ! determines the "scale factor"
  call get_attenuation_scale_factor(myrank,f_c_source,tau_eps,tau_sigma,Q_mu,factor_scale,ATTENUATION_f0_REFERENCE)

  ! determines tau_eps for Q_kappa
  call get_attenuation_tau_eps(Q_kappa,tau_sigma,tau_eps_kappa,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! determines one_minus_sum_beta
  call get_attenuation_property_values(tau_sigma,tau_eps_kappa,beta_kappa,one_minus_sum_beta_kappa)

  ! determines the "scale factor"
  call get_attenuation_scale_factor(myrank,f_c_source,tau_eps_kappa,tau_sigma,Q_kappa,factor_scale_kappa,ATTENUATION_f0_REFERENCE)

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

  ! local parameters
  double precision,dimension(N_SLS) :: tauinv
  integer :: i

  ! inverse of stress relaxation times
  tauinv(:) = -1.0d0 / tau_s(:)

  ! see e.g. Komatitsch & Tromp 1999, eq. (7)

  ! coefficients beta
  beta(:) = 1.0d0 - tau_eps(:) / tau_s(:)

  ! sum of coefficients beta
  one_minus_sum_beta = 1.0d0
  do i = 1,N_SLS
    one_minus_sum_beta = one_minus_sum_beta - beta(i)
  enddo

  end subroutine get_attenuation_property_values

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_scale_factor(myrank,f_c_source,tau_eps,tau_sigma,Q_val,scale_factor,ATTENUATION_f0_REFERENCE)

! returns: physical dispersion scaling factor scale_factor

  use constants

  implicit none

  integer :: myrank
  double precision, intent(in) :: ATTENUATION_f0_REFERENCE
  double precision :: scale_factor, Q_val, f_c_source
  ! strain and stress relaxation times
  double precision, dimension(N_SLS) :: tau_eps, tau_sigma

  ! local parameters
  double precision w_c_source
  double precision factor_scale_mu0, factor_scale_mu
  double precision a_val, b_val
  double precision big_omega
  integer i


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

  !--- compute a, b and Omega parameters
  ! see e.g.:
  !   Liu et al. (1976): eq. 25
  !   using complex modulus Mc = M_R / ( A - i B )
  !   or
  !   Savage et al. (BSSA, 2010): eq. (5) and (6)
  !   complex modulus: M(t) = M_1(t) + i M_2(t)
  a_val = ONE
  b_val = ZERO
  do i = 1,N_SLS
    ! real part M_1 of complex modulus
    a_val = a_val - w_c_source * w_c_source * tau_eps(i) * &
      (tau_eps(i) - tau_sigma(i)) / (1.d0 + w_c_source * w_c_source * tau_eps(i) * tau_eps(i))
    ! imaginary part M_2 of complex modulus
    b_val = b_val + w_c_source * (tau_eps(i) - tau_sigma(i)) / &
      (1.d0 + w_c_source * w_c_source * tau_eps(i) * tau_eps(i))
  enddo

  ! see e.g. Liu et al. (1976): Omega used in equation (20)
  big_omega = a_val * ( sqrt(1.d0 + b_val*b_val/(a_val*a_val)) - 1.d0 )

  !--- quantity by which to scale mu to get mu_relaxed
  factor_scale_mu = b_val * b_val / (TWO * big_omega)

  !--- total factor by which to scale mu0
  scale_factor = factor_scale_mu * factor_scale_mu0

  !--- check that the correction factor is close to one
  if (scale_factor < 0.7 .or. scale_factor > 1.3) then
    write(*,*) "error : in get_attenuation_scale_factor() "
    write(*,*) "  scale factor: ", scale_factor, " should be between 0.7 and 1.3"
    write(*,*) "  Q value = ", Q_val, " central frequency = ",f_c_source
    write(*,*) "  please check your reference frequency ATTENUATION_f0_REFERENCE in constants.h"
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
  f_c_source = 10.0d0**(0.5 * (log10(f1) + log10(f2)))

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

! model_attenuation_variables
!...

  double precision :: Q_in
  double precision, dimension(N_SLS) :: tau_s, tau_eps
  double precision :: min_period,max_period
  double precision :: f0_attenuation,f_min_attenuation,f_max_attenuation

  ! local parameters
  integer :: rw

  ! READ
  rw = 1
  call model_attenuation_storage(Q_in, tau_eps, rw)
  if (rw > 0) return

  call attenuation_invert_by_simplex(min_period, max_period, N_SLS, Q_in, tau_s, tau_eps)

!! DK DK dec 2017: added call to SolvOpt(), which is better (see https://github.com/geodynamics/specfem3d/issues/742 )
  f_min_attenuation = 1.d0 / max_period
  f_max_attenuation = 1.d0 / min_period
  if (N_SLS >= 6) then
    call get_attenuation_source_freq(f0_attenuation,min_period,max_period)
    call compute_attenuation_coeffs(N_SLS,Q_in,f0_attenuation,f_min_attenuation,f_max_attenuation,tau_eps,tau_s)
  endif

  ! WRITE
  rw = -1
  call model_attenuation_storage(Q_in, tau_eps, rw)

  end subroutine get_attenuation_tau_eps

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_attenuation_storage(Qmu, tau_eps, rw)

  use constants

  use attenuation_model, only: AM_S

  implicit none

  double precision :: Qmu, Qmu_new
  double precision, dimension(N_SLS) :: tau_eps
  integer :: rw

  integer :: Qtmp
  integer :: ier

  double precision, parameter :: ZERO_TOL = 1.e-5

  integer, save :: first_time_called = 1

  ! allocates arrays when first called
  if (first_time_called == 1) then
    first_time_called = 0
    AM_S%Q_resolution = 10**ATTENUATION_COMP_RESOLUTION
    AM_S%Q_max = ATTENUATION_COMP_MAXIMUM
    Qtmp = AM_S%Q_resolution * AM_S%Q_max

    allocate(AM_S%tau_eps_storage(N_SLS, Qtmp), &
             AM_S%Qmu_storage(Qtmp),stat=ier)
    if (ier /= 0) stop 'error allocating arrays for attenuation storage'
    AM_S%Qmu_storage(:) = -1
  endif

  if (Qmu < 0.0d0 .or. Qmu > AM_S%Q_max) then
    print *,'Error attenuation_storage()'
    print *,'Attenuation Value out of Range: ', Qmu
    print *,'Attenuation Value out of Range: Min, Max ', 0, AM_S%Q_max
    stop 'Attenuation Value out of Range'
  endif

  if (rw > 0 .and. Qmu <= ZERO_TOL) then
    Qmu = 0.0d0;
    tau_eps(:) = 0.0d0;
    return
  endif
  ! Generate index for Storage Array
  ! and Recast Qmu using this index
  ! Accroding to Brian, use float
  !Qtmp = Qmu * Q_resolution
  !Qmu = Qtmp / Q_resolution;

  ! by default: resolution is Q_resolution = 10
  ! converts Qmu to an array integer index:
  ! e.g. Qmu = 150.31 -> Qtmp = 150.31 * 10 = int( 1503.10 ) = 1503
  Qtmp = int(Qmu * dble(AM_S%Q_resolution))

  ! rounds to corresponding double value:
  ! e.g. Qmu_new = dble( 1503 ) / dble(10) = 150.30
  ! but Qmu_new is not used any further...
  Qmu_new = dble(Qtmp) / dble(AM_S%Q_resolution)

  if (rw > 0) then
    ! checks
    if (first_time_called == 0) then
      if (.not. associated(AM_S%Qmu_storage)) &
        stop 'error calling model_attenuation_storage() routine without AM_S array'
    else
      stop 'error calling model_attenuation_storage() routine with first_time_called value invalid'
    endif

    ! READ
    if (AM_S%Qmu_storage(Qtmp) > 0) then
      ! READ SUCCESSFUL
      tau_eps(:) = AM_S%tau_eps_storage(:,Qtmp)
      Qmu = AM_S%Qmu_storage(Qtmp)
      rw = 1
    else
      ! READ NOT SUCCESSFUL
      rw = -1
    endif
  else
    ! WRITE SUCCESSFUL
    AM_S%tau_eps_storage(:,Qtmp) = tau_eps(:)
    AM_S%Qmu_storage(Qtmp) = Qmu
    rw = 1
  endif

  end subroutine model_attenuation_storage

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_invert_by_simplex(t2, t1, n, Q_real, tau_s, tau_eps)

  implicit none

  ! Input / Output
  double precision  t1, t2
  double precision  Q_real
!  double precision  omega_not
  integer  n
  double precision, dimension(n)   :: tau_s, tau_eps

  ! Internal
  integer i, iterations, err,prnt
  double precision f1, f2, exp1,exp2, min_value !, dexpval
  integer, parameter :: nf = 100
  double precision, dimension(nf) :: f
  double precision, parameter :: PI = 3.14159265358979d0
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


  ! Shove the paramters into the module
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

  integer nf_in, nsls_in
  double precision Q_in
  double precision, dimension(nf_in)   :: f_in
  double precision, dimension(nsls_in) :: tau_s_in
  integer ier

  allocate(AS_V%f(nf_in), &
           AS_V%tau_s(nsls_in),stat=ier)
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
  double precision, dimension(AS_V%nsls) :: Xin
  double precision, dimension(AS_V%nsls) :: tau_eps

  double precision, dimension(AS_V%nf)   :: A, B, tan_delta

  integer i
  double precision xi, iQ2

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

  implicit none

  ! Input
  integer nf, nsls
  double precision, dimension(nf)   :: f
  double precision, dimension(nsls) :: tau_s, tau_eps
  ! Output
  double precision, dimension(nf)   :: A,B

  integer i,j
  double precision w, pi, demon

  PI = 3.14159265358979d0

  A(:) = 1.0d0 -  nsls*1.0d0
  B(:) = 0.0d0
  do i = 1,nf
    w = 2.0d0 * PI * 10**f(i)
    do j = 1,nsls
      !        write(*,*)j,tau_s(j),tau_eps(j)
      demon = 1.0d0 + w**2 * tau_s(j)**2
      A(i) = A(i) + ((1.0d0 + (w**2 * tau_eps(j) * tau_s(j)))/ demon)
      B(i) = B(i) + ((w * (tau_eps(j) - tau_s(j))) / demon)
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

  integer n
  double precision x(n) ! Also Output
  integer itercount, prnt, err
  double precision tolf

  !Internal
  integer i,j, how
  integer, parameter :: none             = 0
  integer, parameter :: initial          = 1
  integer, parameter :: expand           = 2
  integer, parameter :: reflect          = 3
  integer, parameter :: contract_outside = 4
  integer, parameter :: contract_inside  = 5
  integer, parameter :: shrink           = 6

  integer maxiter, maxfun
  integer func_evals
  double precision tolx

  double precision rho, chi, psi, sigma
  double precision xin(n), y(n), v(n,n+1), fv(n+1)
  double precision vtmp(n,n+1)
  double precision usual_delta, zero_term_delta
  double precision xbar(n), xr(n), fxr, xe(n), fxe, xc(n), fxc, fxcc, xcc(n)
  integer place(n+1)

  double precision max_size_simplex, max_value

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
  integer n
  double precision fv(n)

  integer i
  double precision m, z

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
  integer n
  double precision v(n,n+1)

  integer i,j
  double precision m, z

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

  integer n
  double precision X(n)
  integer I(n)

  integer j,k
  double precision rtmp
  integer itmp

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

!
!--------------------------------------------------------------------------------
!

! use of SolvOpt to compute attenuation relaxation mechanisms,
! from Emilie Blanc, Bruno Lombard and Dimitri Komatitsch, CNRS Marseille, France, for a Generalized Zener body model.

! If you use this code for your own research, please cite some (or all) of these articles:
!
! @Article{BlKoChLoXi16,
! Title   = {Highly accurate stability-preserving optimization of the {Z}ener viscoelastic model,
!            with application to wave propagation in the presence of strong attenuation},
! Author  = {\'Emilie Blanc and Dimitri Komatitsch and Emmanuel Chaljub and Bruno Lombard and Zhinan Xie},
! Journal = {Geophysical Journal International},
! Year    = {2016},
! Number  = {1},
! Pages   = {427-439},
! Volume  = {205},
! Doi     = {10.1093/gji/ggw024}}

! The SolvOpt algorithm was developed by Franz Kappel and Alexei V. Kuntsevich
! and is available open source at http://www.uni-graz.at/imawww/kuntsevich/solvopt
!
! It is described in Kappel and Kuntsevich, An Implementation of Shor's r-Algorithm,
! Computational Optimization and Applications, vol. 15, p. 193-205 (2000).

!-------------------------------------------------------------------------

! From Bruno Lombard, May 2014:

! En interne dans le code ci-dessous on travaille en (Theta, Kappa).
! Les Theta sont les points et les Kappa sont les poids.
! Pour repasser en (Tau_Sigma, Tau_Epsilon), on doit appliquer les formules:
!
! Tau_Sigma = 1 / Theta
! Tau_Epsilon = (1 / Theta) * (1 + Nrelax * Kappa) = Tau_Sigma * (1 + Nrelax * Kappa)

! The system to solve can be found in equation (7) of:
! Lombard and Piraux, Numerical modeling of transient two-dimensional viscoelastic waves,
! Journal of Computational Physics, Volume 230, Issue 15, Pages 6099-6114 (2011)

! Suivant les compilateurs et les options de compilation utilisees,
! il peut y avoir des differences au 4eme chiffre significatif. C'est sans consequences sur la precision du calcul :
! l'erreur est de 0.015 % avec optimization non-lineaire, a comparer a 1.47 % avec Emmerich and Korn (1987).
! Si je relance le calcul en initialisant avec le resultat precedent, ce chiffre varie a nouveau tres legerement.

!-------------------------------------------------------------------------

! From Bruno Lombard, June 2014:

! j'ai relu en detail :

! [1] Carcione, Kosslof, Kosslof, "Viscoacoustic wave propagation simulation in the Earth",
!            Geophysics 53-6 (1988), 769-777
!
! [2] Carcione, Kosslof, Kosslof, "Wave propagation simulation in a linear viscoelastic medium",
!            Geophysical Journal International 95 (1988), 597-611
!
! [3] Moczo, Kristek, "On the rheological models used for time-domain methods of seismic wave propagation",
!            Geophysical Research Letters 32 (2005).

! Le probleme provient probablement d'une erreur recurrente dans [1,2] et datant de Liu et al 1976 :
! l'oubli du facteur 1/N dans la fonction de relaxation d'un modele de Zener a N elements.
! Il est effectivement facile de faire l'erreur. Voir l'equation (12) de [3], et le paragraphe qui suit.

! Du coup le module de viscoelasticite est faux dans [1,2], et donc le facteur de qualite,
! et donc les temps de relaxation tau_sigma...

! Apres, [2] calcule une solution analytique juste, mais avec des coefficients sans sens physique.
! Et quand SPECFEM2D obtient un bon accord avec cette solution analytique, ca valide SPECFEM, mais pas le calcul des coefficients.

! Il y a donc une erreur dans [1,2], et [3] a raison.

! Sa solution analytique decoule d'un travail sur ses fonctions de relaxation (A4),
! qu'il injecte ensuite dans la relation de comportement (A1) et travaille en Fourier.

! Le probleme est que sa fonction de relaxation (A4) est fausse : il manque 1/N.
! De ce fait, sa solution analytique est coherente avec sa solution numerique.
! Dans les deux cas, ce sont les memes temps de relaxation qui sont utilises. Mais ces temps sont calcules de facon erronee.

!-------------------------------------------------------------------------

! From Dimitri Komatitsch, June 2014:

! In [2] Carcione, Kosslof, Kosslof, "Wave propagation simulation in a linear viscoelastic medium",
!            Geophysical Journal International 95 (1988), 597-611
! there is another mistake: in Appendix B page 611 Carcione writes omega/(r*v),
! but that is not correct, it should be omega*r/v instead.

!---------------------------------------------------

! From Emilie Blanc, April 2014:

! le programme SolvOpt d'optimization non-lineaire
! avec contrainte. Ce programme prend quatre fonctions en entree :

! - fun() est la fonction a minimiser

! - grad() est le gradient de la fonction a minimiser par rapport a chaque parametre

! - func() est le maximum des residus (= 0 si toutes les contraintes sont satisfaites)

! - gradc() est le gradient du maximum des residus (= 0 si toutes les
! contraintes sont satisfaites)

! Ce programme a ete developpe par Kappel et Kuntsevich. Leur article decrit l'algorithme.

! J'ai utilise ce code pour la poroelasticite haute-frequence, et aussi en
! viscoelasticite fractionnaire (modele d'Andrade, avec Bruno Lombard et
! Cedric Bellis). Nous pouvons interagir sur l'algorithme d'optimization
! pour votre modele visco, et etudier l'effet des coefficients ainsi obtenus.

!---------------------------------------------------

! From Emilie Blanc, March 2014:

! Les entrees du programme principal sont le nombre de variables
! diffusives, le facteur de qualite voulu Qref et la frequence centrale f0.

! Cependant, pour l'optimization non-lineaire, j'ai mis theta_max=100*f0
! et non pas theta_max=2*pi*100*f0. En effet, dans le programme, on
! travaille sur les frequences, et non pas sur les frequences angulaires.
! Cela dit, dans les deux cas j'obtiens les memes coefficients...

!---------------------------------------------------

  subroutine compute_attenuation_coeffs(N,Qref,f0,f_min,f_max,tau_epsilon,tau_sigma)

  implicit none

  integer, intent(in) :: N
  double precision, intent(in) :: Qref,f_min,f_max,f0
  double precision, dimension(1:N), intent(out) :: tau_epsilon,tau_sigma

  integer :: i
  double precision, dimension(1:N) :: point,weight

! nonlinear optimization with constraints
  call nonlinear_optimization(N,Qref,f0,point,weight,f_min,f_max)

  do i = 1,N
    tau_sigma(i) = 1.d0 / point(i)
    tau_epsilon(i) = tau_sigma(i) * (1.d0 + N * weight(i))
  enddo

! print *,'points = '
! do i = 1,N
!   print *,point(i)
! enddo
! print *

! print *,'weights = '
! do i = 1,N
!   print *,weight(i)
! enddo
! print *

! print *,'tau_epsilon = '
! do i = 1,N
!   print *,tau_epsilon(i)
! enddo
! print *

! print *,'tau_sigma = '
! do i = 1,N
!   print *,tau_sigma(i)
! enddo
! print *

  end subroutine compute_attenuation_coeffs

!---------------------------------------------------

! classical calculation of the coefficients based on linear least squares

  subroutine decomposition_LU(a,i_min,n,indx,d)

  implicit none

  integer, intent(in) :: i_min,n
  double precision, intent(out) :: d
  integer, dimension(i_min:n), intent(inout) :: indx
  double precision, dimension(i_min:n,i_min:n), intent(inout) :: a

  integer i,imax,j,k
  double precision big,dum,somme,eps
  double precision, dimension(i_min:n) :: vv

  imax = 0
  d = 1.
  eps = 1.e-20

  do i = i_min,n
    big = 0.
    do j = i_min,n
      if (abs(a(i,j)) > big) then
        big = abs(a(i,j))
      endif
    enddo
    if (big == 0.) then
      print *,'Singular matrix in routine decomposition_LU'
    endif
    vv(i) = 1./big
  enddo

  do j = i_min,n
    do i = i_min,j-1
      somme = a(i,j)
      do k = i_min,i-1
        somme = somme - a(i,k)*a(k,j)
      enddo
      a(i,j) = somme
    enddo

    big = 0.

    do i = j,n
      somme = a(i,j)
      do k = i_min,j-1
        somme = somme - a(i,k)*a(k,j)
      enddo
      a(i,j) = somme
      dum = vv(i)*abs(somme)
      if (dum >= big) then
        big = dum
        imax = i
      endif
    enddo

    if (j /= imax) then
      do k = i_min,n
        dum = a(imax,k)
        a(imax,k) = a(j,k)
        a(j,k) = dum
      enddo
      d = -d
      vv(imax) = vv(j)
    endif

    indx(j) = imax
    if (a(j,j) == 0.) then
      a(j,j) = eps
    endif
    if (j /= n) then
      dum = 1./a(j,j)
      do i = j+1,n
        a(i,j) = a(i,j)*dum
      enddo
    endif
  enddo

  end subroutine decomposition_LU

  subroutine LUbksb(a,i_min,n,indx,b,m)

  implicit none

  integer, intent(in) :: i_min,n,m
  integer, dimension(i_min:n), intent(in) :: indx
  double precision, dimension(i_min:n,i_min:n), intent(in) :: a
  double precision, dimension(i_min:n,i_min:m), intent(inout) :: b

  integer i,ip,j,ii,k
  double precision somme

  do k = i_min,m

    ii = -1

    do i = i_min,n
      ip = indx(i)
      somme = b(ip,k)
      b(ip,k) = b(i,k)
      if (ii /= -1) then
        do j = ii,i-1
          somme = somme - a(i,j)*b(j,k)
        enddo
      else if (somme /= 0.) then
        ii = i
      endif
      b(i,k) = somme
    enddo

    do i = n,i_min,-1
      somme = b(i,k)
      do j = i+1,n
        somme = somme - a(i,j)*b(j,k)
      enddo
      b(i,k) = somme/a(i,i)
    enddo
  enddo

  end subroutine LUbksb

  subroutine syst_LU(a,i_min,n,b,m)

  implicit none

  integer, intent(in) :: i_min,n,m
  double precision, dimension(i_min:n,i_min:n), intent(in) :: a
  double precision, dimension(i_min:n,i_min:m), intent(inout) :: b

  integer i,j
  integer, dimension(i_min:n) :: indx
  double precision d
  double precision, dimension(i_min:n,i_min:n) :: aux

  do j = i_min,n
    indx(j) = 0
    do i = i_min,n
      aux(i,j) = a(i,j)
    enddo
  enddo

  call decomposition_LU(aux,i_min,n,indx,d)
  call LUbksb(aux,i_min,n,indx,b,m)

  end subroutine syst_LU

  subroutine lfit_zener(x,y,sig,ndat,poids,ia,covar,chisq,ma,Qref,point)
! ma = nombre de variable diffusive
! ndat = m = K nombre d'abcisse freq_k

  implicit none

  integer, intent(in) :: ndat,ma
  logical, dimension(1:ma), intent(in) :: ia
  double precision, intent(in) :: Qref
  double precision, intent(out) :: chisq
  double precision, dimension(1:ndat), intent(in) :: x,y,sig
  double precision, dimension(1:ma), intent(in) :: point
  double precision, dimension(1:ma), intent(out) :: poids
  double precision, dimension(1:ma,1:ma), intent(out) :: covar

  integer i,j,k,l,mfit
  double precision ym,wt,sig2i
  double precision, dimension(1:ma) :: afunc
  double precision, dimension(1:ma,1:1) :: beta

  mfit = 0

  do j = 1,ma
    if (ia(j)) then
      mfit = mfit + 1
    endif
  enddo
  if (mfit == 0) then
    print *,'lfit: no parameters to be fitted'
  endif

  do j = 1,mfit
    beta(j,1) = 0.
    do k = 1,mfit
      covar(j,k) = 0.
    enddo
  enddo

  do i = 1,ndat
    call func_zener(x(i),afunc,ma,Qref,point)
    ym = y(i)
    if (mfit < ma) then
      do j = 1,ma
        if (.not. ia(j)) then
          ym = ym - poids(j) * afunc(j)
        endif
      enddo
    endif
    sig2i = 1. / (sig(i) * sig(i))
    j = 0
    do l= 1,ma
      if (ia(l)) then
        j = j+1
        wt = afunc(l) * sig2i
        k = count(ia(1:l))
        covar(j,1:k) = covar(j,1:k) + wt * pack(afunc(1:l),ia(1:l))
        beta(j,1) = beta(j,1) + ym * wt
      endif
    enddo
  enddo

  do j = 2,mfit,1
  do k = 1,j-1,1
    covar(k,j) = covar(j,k)
  enddo
  enddo

  if (ma == 1) then
    poids(1) = beta(1,1)/covar(1,1)
  else if (ma > 1) then
    call syst_LU(covar,1,mfit,beta,1)
    poids(1:ma) = unpack(beta(1:ma,1),ia,poids(1:ma))
  endif

  chisq = 0.
  do i = 1,ndat
    call func_zener(x(i),afunc,ma,Qref,point)
    chisq=chisq+((y(i)-dot_product(poids(1:ma),afunc(1:ma)))/sig(i))**2
  enddo

  end subroutine lfit_zener

  subroutine func_zener(x,afunc,N,Qref,point)

  implicit none

  integer, intent(in) :: N
  double precision, intent(in) :: x,Qref
  double precision, dimension(1:N), intent(in) :: point
  double precision, dimension(1:N), intent(out) :: afunc

  integer k
  double precision num,deno

  do k = 1,N
    num  = x * (point(k) - x / Qref)
    deno = point(k) * point(k) + x * x
    afunc(k) = num / deno
  enddo

  end subroutine func_zener

  subroutine remplit_point(fmin,fmax,N,point)

  use constants, only: TWO_PI

  implicit none

  integer, intent(in) :: N
  double precision, intent(in) :: fmin,fmax
  double precision, dimension(1:N), intent(out) :: point

  integer l

  if (N == 1) then
    point(1) = sqrt(fmin * fmax)
  ELSE
    do l = 1, N, 1
      point(l) = (fmax/fmin) ** ((l-1.)/(N-1.))
! we work in angular frequency, not frequency
      point(l) = TWO_PI * point(l) * fmin
    enddo
  endif

  end subroutine remplit_point

  subroutine classical_linear_least_squares(Qref,poids,point,N,fmin,fmax)

  use constants, only: TWO_PI

  implicit none

  integer, intent(in) :: N
  double precision, intent(in) :: Qref,fmin,fmax
  double precision, dimension(1:N), intent(out) :: point,poids

  integer k,m
  logical, dimension(1:N) :: ia
  double precision ref,freq,chi2
  double precision, dimension(1:N,1:N) :: covar
  double precision, dimension(1:2*N-1) :: x,y_ref,sig

  m = 2*N-1

  call remplit_point(fmin,fmax,N,point)

  ref = 1.0 / Qref

  do k = 1,m
    freq = (fmax/fmin) ** ((k - 1.)/(m - 1.))
! we work in angular frequency, not frequency
    freq = TWO_PI * fmin * freq
    x(k) = freq
    y_ref(k) = ref
    sig(k) = 1.
  enddo

  do k = 1,N
    ia(k) = .true.
  enddo

  call lfit_zener(x,y_ref,sig,m,poids,ia,covar,chi2,N,Qref,point)

  end subroutine classical_linear_least_squares

! Calcul des coefficients par optimization non-lineaire avec contraintes

  subroutine solvopt(n,x,f,fun,flg,grad,options,flfc,func,flgc,gradc,Qref,Kopt,theta_min,theta_max,f_min,f_max)
!-----------------------------------------------------------------------------
! The subroutine SOLVOPT performs a modified version of Shor's r-algorithm in
! order to find a local minimum resp. maximum of a nonlinear function
! defined on the n-dimensional Euclidean space
! or
! a local minimum for a nonlinear constrained problem:
! min { f(x): g(x) ( < )= 0, g(x) in R(m), x in R(n) }.
! Arguments:
! n       is the space dimension (integer*4),
! x       is the n-vector, the coordinates of the starting point
!         at a call to the subroutine and the optimizer at regular return
!         (double precision),
! f       returns the optimum function value
!         (double precision),
! fun     is the entry name of a subroutine which computes the value
!         of the function fun at a point x, should be declared as external
!         in a calling routine,
!        synopsis: fun(x,f)
! grad    is the entry name of a subroutine which computes the gradient
!         vector of the function fun at a point x, should be declared as
!         external in a calling routine,
!         synopsis: grad(x,g)
! func    is the entry name of a subroutine which computes the MAXIMAL
!         RESIDIAL!!! (a scalar) for a set of constraints at a point x,
!         should be declared as external in a calling routine,
!         synopsis: func(x,fc)
! gradc   is the entry name of a subroutine which computes the gradient
!         vector for a constraint with the MAXIMAL RESIDUAL at a point x,
!         should be declared as external in a calling routine,
!        synopsis: gradc(x,gc)
! flg,    (logical) is a flag for the use of a subroutine grad:
!         .true. means gradients are calculated by the user-supplied routine.
! flfc,   (logical) is a flag for a constrained problem:
!         .true. means the maximal residual for a set of constraints
!         is calculated by func.
! flgc,   (logical) is a flag for the use of a subroutine gradc:
!         .true. means gradients of the constraints are calculated
!         by the user-supplied routine.
! options is a vector of optional parameters (double precision):
!     options(1)= H, where sign(H)=-1 resp. sign(H)=+1 means minimize resp.
!         maximize fun (valid only for an unconstrained problem) and
!         H itself is a factor for the initial trial step size
!         (options(1)=-1.d0 by default),
!     options(2)= relative error for the argument in terms of the infinity-norm
!         (1.d-4 by default),
!     options(3)= relative error for the function value (1.d-6 by default),
!     options(4)= limit for the number of iterations (1.5d4 by default),
!     options(5)= control of the display of intermediate results and error
!         resp. warning messages (default value is 0.d0, i.e., no intermediate
!         output but error and warning messages, see the manual for more),
!     options(6)= maximal admissible residual for a set of constraints
!         (options(6)=1.d-8 by default, see the manual for more),
!    *options(7)= the coefficient of space dilation (2.5d0 by default),
!    *options(8)= lower bound for the stepsize used for the difference
!        approximation of gradients (1.d-11 by default,see the manual for more).
!   (* ... changes should be done with care)
! returned optional values:
!     options(9),  the number of iterations, if positive,
!         or an abnormal stop code, if negative (see manual for more),
!                -1: allocation error,
!                -2: improper space dimension,
!                -3: fun returns an improper value,
!                -4: grad returns a zero vector or improper value at the
!                    starting point,
!                -5: func returns an improper value,
!                -6: gradc returns an improper value,
!                -7: function is unbounded,
!                -8: gradient is zero at the point,
!                    but stopping criteria are not fulfilled,
!                -9: iterations limit exceeded,
!               -11: Premature stop is possible,
!               -12: Result may not provide the true optimum,
!               -13: function is flat: result may be inaccurate
!                   in view of a point.
!               -14: function is steep: result may be inaccurate
!                    in view of a function value,
!       options(10), the number of objective function evaluations, and
!       options(11), the number of gradient evaluations.
!       options(12), the number of constraint function evaluations, and
!       options(13), the number of constraint gradient evaluations.
! ____________________________________________________________________________
!
      implicit none
      !include 'messages.inc'

      integer, intent(in) :: Kopt
      double precision, intent(in) :: Qref,theta_min,theta_max,f_min,f_max

      logical flg,flgc,flfc, constr, app, appconstr
      logical FsbPnt, FsbPnt1, termflag, stopf
      logical stopping, dispwarn, Reset, ksm,knan,obj
      integer n, kstore, ajp,ajpp,knorms, k, kcheck, numelem
      integer dispdata, ld, mxtc, termx, limxterm, nzero, krerun
      integer warnno, kflat, stepvanish, i,j,ni,ii, kd,kj,kc,ip
      integer iterlimit, kg,k1,k2, kless,   allocerr
      double precision options(13),doptions(13)
      double precision x(n),f
      double precision nsteps(3), gnorms(10), kk, nx
      double precision ajb,ajs, des, dq,du20,du10,du03
      double precision n_float, cnteps
      double precision low_bound, ZeroGrad, ddx, y
      double precision lowxbound, lowfbound, detfr, detxr, grbnd
      double precision fp,fp1,fc,f1,f2,fm,fopt,frec,fst, fp_rate
      double precision PenCoef, PenCoefNew
      double precision gamma,w,wdef,h1,h,hp
      double precision dx,ng,ngc,nng,ngt,nrmz,ng1,d,dd, laststep
      double precision zero,one,two,three,four,five,six,seven
      double precision eight,nine,ten,hundr
      double precision infty, epsnorm,epsnorm2,powerm12
      double precision, dimension(:,:), allocatable :: B
      double precision, dimension(:), allocatable :: g
      double precision, dimension(:), allocatable :: g0
      double precision, dimension(:), allocatable :: g1
      double precision, dimension(:), allocatable :: gt
      double precision, dimension(:), allocatable :: gc
      double precision, dimension(:), allocatable :: z
      double precision, dimension(:), allocatable :: x1
      double precision, dimension(:), allocatable :: xopt
      double precision, dimension(:), allocatable :: xrec
      double precision, dimension(:), allocatable :: grec
      double precision, dimension(:), allocatable :: xx
      double precision, dimension(:), allocatable :: deltax
      integer, dimension(:), allocatable :: idx
      character(len=100) :: endwarn
      character(len=19) :: allocerrstr
      external fun,grad,func,gradc

      data zero/0.d0/, one/1.d0/, two/2.d0/, three/3.d0/, four/4.d0/, &
         five/5.d0/, six/6.d0/, seven/7.d0/, eight/8.d0/, nine/9.d0/, &
         ten/1.d1/,  hundr/1.d2/, powerm12/1.d-12/, &
         infty /1.d100/, epsnorm /1.d-15/,  epsnorm2 /1.d-30/, &
         allocerrstr/'Allocation Error = '/
! Check the dimension:
      if (n < 2) then
          print *, 'SolvOpt error:'
          print *, 'Improper space dimension.'
         stop 'error in allocate statement in SolvOpt'
        options(9)=-one
        goto 999
      endif
      n_float=dble(n)
! allocate working arrays:
      allocate (B(n,n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (g(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (g0(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (g1(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (gt(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (gc(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (z(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (x1(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (xopt(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (xrec(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (grec(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (xx(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (deltax(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif
      allocate (idx(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         stop 'error in allocate statement in SolvOpt'
      endif

! store flags:
      app= .not. flg
      constr=flfc
      appconstr= .not. flgc
! Default values for options:
      call soptions(doptions)
      do i =1,8
            if (options(i) == zero) then
               options(i)=doptions(i)
            else if (i == 2 .or. i == 3 .or. i == 6) then
               options(i)=dmax1(options(i),powerm12)
               options(i)=dmin1(options(i),one)
               if (i == 2)options(i)=dmax1(options(i),options(8)*hundr)
            else if (i == 7) then
               options(7)=dmax1(options(i),1.5d0)
            endif
      enddo

! WORKING CONSTANTS AND COUNTERS ----{

      options(10)=zero    !! counter for function calculations
      options(11)=zero    !! counter for gradient calculations
      options(12)=zero    !! counter for constraint function calculations
      options(13)=zero    !! counter for constraint gradient calculations
      iterlimit=idint(options(4))
      if (constr) then
        h1=-one           !! NLP: restricted to minimization
        cnteps=options(6)
      else
        h1=dsign(one,options(1))  !! Minimize resp. maximize a function
      endif
      k=0                         !! Iteration counter
      wdef=one/options(7)-one     !! Default space transf. coeff.

! Gamma control ---{
      ajb=one+1.d-1/n_float**2    !! Base I
      ajp=20
      ajpp=ajp                    !! Start value for the power
      ajs=1.15d0                  !! Base II
      knorms=0
      do i =1,10
       gnorms(i)=zero
      enddo
!---}
! Display control ---{
      if (options(5) <= zero) then
         dispdata=0
         if (options(5) == -one) then
            dispwarn=.false.
         else
            dispwarn=.true.
         endif
      else
         dispdata=idnint(options(5))
         dispwarn=.true.
      endif
      ld=dispdata
!---}

! Stepsize control ---{
      dq=5.1d0           !! Step divider (at f_{i+1} > gamma*f_{i})
      du20=two
      du10=1.5d0
      du03=1.05d0        !! Step multipliers (at certain steps made)
      kstore=3
      do i = 1,kstore
       nsteps(i)=zero    !! Steps made at the last 'kstore' iterations
      enddo
      if (app) then
        des=6.3d0        !! Desired number of steps per 1-D search
      else
        des=3.3d0
      endif
      mxtc=3             !! Number of trial cycles (steep wall detect)
!---}
      termx=0
      limxterm=50        !! Counter and limit for x-criterion
! stepsize for gradient approximation
      ddx=dmax1(1.d-11,options(8))

      low_bound=-one+1.d-4     !! Lower bound cosine used to detect a ravine
      ZeroGrad=n_float*1.d-16  !! Lower bound for a gradient norm
      nzero=0                  !! Zero-gradient events counter
! Low bound for the values of variables to take into account
      lowxbound=dmax1(options(2),1.d-3)
! Lower bound for function values to be considered as making difference
      lowfbound=options(3)**2
      krerun=0                 !! Re-run events counter
      detfr=options(3)*hundr   !! Relative error for f/f_{record}
      detxr=options(2)*ten     !! Relative error for norm(x)/norm(x_{record})
      warnno=0                 !! the number of warn.mess. to end with
      kflat=0                  !! counter for points of flatness
      stepvanish=0             !! counter for vanished steps
      stopf=.false.
! ----}  End of setting constants
! ----}  End of the preamble
!--------------------------------------------------------------------
! Compute the function  ( first time ) ----{
      call fun(x,f,Qref,n/2,n,Kopt,f_min,f_max)
      options(10)=options(10)+one
      if (dabs(f) >= infty) then
         if (dispwarn) then
            print *,'SolvOpt error: function equals infinity at the point.'
            print *,'Choose another starting point.'
         endif
         options(9)=-three
         goto 999
      endif
      do i = 1,n
        xrec(i)=x(i)
      enddo
      frec=f     !! record point and function value
! Constrained problem
      if (constr) then
          kless=0
          fp=f
          call func(x,fc,n/2,n,theta_min,theta_max)
          options(12)=options(12)+one
          if (dabs(fc) >= infty) then
             if (dispwarn) then
                print *,'SolvOpt error: FUNC returns infinite value at the point.'
                print *,'Choose another starting point.'
             endif
             options(9)=-five
             goto 999
          endif
        PenCoef=one          !! first rough approximation
        if (fc <= cnteps) then
         FsbPnt=.true.       !! feasible point
         fc=zero
        else
         FsbPnt=.false.
        endif
        f=f+PenCoef*fc
      endif
! ----}
! COMPUTE THE GRADIENT ( FIRST TIME ) ----{
      if (app) then
        do i = 1,n
         deltax(i)=h1*ddx
        enddo
        obj=.true.
        !if (constr) then
           !call apprgrdn()
        !else
           !call apprgrdn()
        !endif
        options(10)=options(10)+n_float
      else
        call grad(x,g,Qref,n/2,n,Kopt,f_min,f_max)
        options(11)=options(11)+one
      endif
      ng=zero
      do i = 1,n
         ng=ng+g(i)*g(i)
      enddo
      ng=dsqrt(ng)
      if (ng >= infty) then
         if (dispwarn) then
            print *,'SolvOpt error:'
            print *,'Gradient equals infinity at the starting point.'
            print *,'Choose another starting point.'
         endif
         options(9)=-four
         goto 999
      else if (ng < ZeroGrad) then
         if (dispwarn) then
            print *,'SolvOpt error:'
            print *,'Gradient equals zero at the starting point.'
            print *,'Choose another starting point.'
         endif
         options(9)=-four
         goto 999
      endif
      if (constr) then
       if (.not. FsbPnt) then
         !if (appconstr) then
            !do j = 1,n
              !if (x(j) >= zero) then
                 !deltax(j)=ddx
              !else
                 !deltax(j)=-ddx
              !endif
            !enddo
            !obj=.false.
            !call apprgrdn()
         if (.not. appconstr) then
            call gradc(x,gc,n/2,n,theta_min,theta_max)
         endif
         ngc=zero
         do i = 1,n
           ngc=ngc+gc(i)*gc(i)
         enddo
         ngc=dsqrt(ngc)
         if (ng >= infty) then
            if (dispwarn) then
               print *,'SolvOpt error: GRADC returns infinite vector at the point.'
               print *,'Choose another starting point.'
            endif
            options(9)=-six
            goto 999
         else if (ng < ZeroGrad) then
            if (dispwarn) then
               print *,'SolvOpt error: GRADC returns zero vector at an infeasible point.'
            endif
            options(9)=-six
            goto 999
         endif
         do i = 1,n
           g(i)=g(i)+PenCoef*gc(i)
         enddo
         ng=zero
         do i = 1,n
           ng=ng+g(i)*g(i)
           grec(i)=g(i)
         enddo
         ng=dsqrt(ng)
       endif
      endif
      do i = 1,n
        grec(i)=g(i)
      enddo
      nng=ng
! ----}
! INITIAL STEP SIZE
      d=zero
      do i = 1,n
        if (d < dabs(x(i))) d=dabs(x(i))
      enddo
      h=h1*dsqrt(options(2))*d                  !! smallest possible stepsize
      if (dabs(options(1)) /= one) then
        h=h1*dmax1(dabs(options(1)),dabs(h))    !! user-supplied stepsize
      else
          h=h1*dmax1(one/dlog(ng+1.1d0),dabs(h)) !! calculated stepsize
      endif

! RESETTING LOOP ----{
      do while (.true.)
        kcheck=0                       !! Set checkpoint counter.
        kg=0                           !! stepsizes stored
        kj=0                           !! ravine jump counter
        do i = 1,n
          do j = 1,n
            B(i,j)=zero
          enddo
          B(i,i)=one                   !! re-set transf. matrix to identity
          g1(i)=g(i)
        enddo
        fst=f
        dx=0
! ----}

! MAIN ITERATIONS ----{

        do while (.true.)
          k=k+1
          kcheck=kcheck+1
          laststep=dx
! ADJUST GAMMA --{
           gamma=one+dmax1(ajb**((ajp-kcheck)*n),two*options(3))
           gamma=dmin1 ( gamma,ajs**dmax1(one,dlog10(nng+one)) )
! --}
       ngt=zero
       ng1=zero
       dd=zero
       do i = 1,n
         d=zero
         do j = 1,n
            d=d+B(j,i)*g(j)
         enddo
         gt(i)=d
         dd=dd+d*g1(i)
         ngt=ngt+d*d
         ng1=ng1+g1(i)*g1(i)
       enddo
       ngt=dsqrt(ngt)
       ng1=dsqrt(ng1)
       dd=dd/ngt/ng1

       w=wdef
! JUMPING OVER A RAVINE ----{
       if (dd < low_bound) then
        if (kj == 2) then
          do i = 1,n
           xx(i)=x(i)
          enddo
        endif
        if (kj == 0) kd=4
        kj=kj+1
        w=-.9d0              !! use large coef. of space dilation
        h=h*two
        if (kj > 2*kd) then
          kd=kd+1
          warnno=1
          endwarn='Premature stop is possible. Try to re-run the routine from the obtained point.'
          do i = 1,n
            if (dabs(x(i)-xx(i)) < epsnorm*dabs(x(i))) then
             if (dispwarn) then
                print *,'SolvOpt warning:'
                print *,'Ravine with a flat bottom is detected.'
             endif
            endif
          enddo
        endif
       else
        kj=0
       endif
! ----}
! DILATION ----{
       nrmz=zero
       do i = 1,n
         z(i)=gt(i)-g1(i)
         nrmz=nrmz+z(i)*z(i)
       enddo
       nrmz=dsqrt(nrmz)
       if (nrmz > epsnorm*ngt) then
        do i = 1,n
         z(i)=z(i)/nrmz
        enddo
! New direction in the transformed space: g1=gt+w*(z*gt')*z and
! new inverse matrix: B = B ( I + (1/alpha -1)zz' )
        d = zero
        do i = 1,n
          d=d+z(i)*gt(i)
        enddo
        ng1=zero
        d = d*w
        do i = 1,n
          dd=zero
          g1(i)=gt(i)+d*z(i)
          ng1=ng1+g1(i)*g1(i)
          do j = 1,n
             dd=dd+B(i,j)*z(j)
          enddo
          dd=w*dd
          do j = 1,n
            B(i,j)=B(i,j)+dd*z(j)
          enddo
        enddo
        ng1=dsqrt(ng1)
       else
        do i = 1,n
         z(i)=zero
         g1(i)=gt(i)
        enddo
        nrmz=zero
       endif
       do i = 1,n
           gt(i)=g1(i)/ng1
       enddo
        do i = 1,n
          d=zero
            do j = 1,n
               d=d+B(i,j)*gt(j)
            enddo
          g0(i)=d
        enddo
! ----}
! RESETTING ----{
        if (kcheck > 1) then
           numelem=0
           do i = 1,n
              if (dabs(g(i)) > ZeroGrad) then
                 numelem=numelem+1
                 idx(numelem)=i
              endif
           enddo
           if (numelem > 0) then
              grbnd=epsnorm*dble(numelem**2)
              ii=0
              do i = 1,numelem
                 j=idx(i)
                 if (dabs(g1(j)) <= dabs(g(j))*grbnd) ii=ii+1
              enddo
              if (ii == n .or. nrmz == zero) then
                if (dispwarn) then
                  print *,'SolvOpt warning: Normal re-setting of a transformation matrix'
                endif
                if (dabs(fst-f) < dabs(f)*1.d-2) then
                   ajp=ajp-10*n
                else
                   ajp=ajpp
                endif
                h=h1*dx/three
                k=k-1
                exit
              endif
           endif
        endif
! ----}
! STORE THE CURRENT VALUES AND SET THE COUNTERS FOR 1-D SEARCH
        do i = 1,n
         xopt(i)=x(i)
        enddo
        fopt=f
        k1=0
        k2=0
        ksm=.false.
        kc=0
        knan=.false.
        hp=h
        if (constr) Reset=.false.
! 1-D SEARCH ----{
        do while (.true.)
         do i = 1,n
          x1(i)=x(i)
         enddo
         f1=f
         if (constr) then
           FsbPnt1=FsbPnt
           fp1=fp
         endif
! NEW POINT
         do i = 1,n
            x(i)=x(i)+hp*g0(i)
         enddo
           ii=0
           do i = 1,n
            if (dabs(x(i)-x1(i)) < dabs(x(i))*epsnorm) ii=ii+1
           enddo
! function value
         call fun(x,f,Qref,n/2,n,Kopt,f_min,f_max)
         options(10)=options(10)+one
         if (h1*f >= infty) then
            if (dispwarn) then
              print *,'SolvOpt error: function is unbounded.'
            endif
            options(9)=-seven
            goto 999
         endif
         if (constr) then
           fp=f
           call func(x,fc,n/2,n,theta_min,theta_max)
           options(12)=options(12)+one
           if (dabs(fc) >= infty) then
               if (dispwarn) then
                  print *,'SolvOpt error: FUNC returns infinite value at the point.'
                  print *,'Choose another starting point.'
               endif
               options(9)=-five
               goto 999
           endif
           if (fc <= cnteps) then
              FsbPnt=.true.
              fc=zero
           else
              FsbPnt=.false.
              fp_rate=fp-fp1
              if (fp_rate < -epsnorm) then
               if (.not. FsbPnt1) then
                d=zero
                do i = 1,n
                  d=d+(x(i)-x1(i))**2
                enddo
                d=dsqrt(d)
                PenCoefNew=-1.5d1*fp_rate/d
                if (PenCoefNew > 1.2d0*PenCoef) then
                  PenCoef=PenCoefNew
                  Reset=.true.
                  kless=0
                  f=f+PenCoef*fc
                  exit
                endif
               endif
              endif
           endif
           f=f+PenCoef*fc
         endif
         if (dabs(f) >= infty) then
             if (dispwarn) then
               print *,'SolvOpt warning: function equals infinity at the point.'
             endif
             if (ksm .or. kc >= mxtc) then
                options(9)=-three
                goto 999
             else
                k2=k2+1
                k1=0
                hp=hp/dq
                do i = 1,n
                 x(i)=x1(i)
                enddo
                f=f1
                knan=.true.
                if (constr) then
                  FsbPnt=FsbPnt1
                  fp=fp1
                endif
             endif
! STEP SIZE IS ZERO TO THE EXTENT OF EPSNORM
         else if (ii == n) then
                stepvanish=stepvanish+1
                if (stepvanish >= 5) then
                    options(9)=-ten-four
                    if (dispwarn) then
                       print *,'SolvOpt: Termination warning:'
                       print *,'Stopping criteria are not fulfilled. The function is very steep at the solution.'
                    endif
                    goto 999
                else
                    do i = 1,n
                     x(i)=x1(i)
                    enddo
                    f=f1
                    hp=hp*ten
                    ksm=.true.
                    if (constr) then
                       FsbPnt=FsbPnt1
                       fp=fp1
                    endif
                endif
! USE SMALLER STEP
         else if (h1*f < h1*gamma**idint(dsign(one,f1))*f1) then
             if (ksm) exit
             k2=k2+1
             k1=0
             hp=hp/dq
             do i = 1,n
              x(i)=x1(i)
             enddo
             f=f1
             if (constr) then
                FsbPnt=FsbPnt1
                fp=fp1
             endif
             if (kc >= mxtc) exit
! 1-D OPTIMIZER IS LEFT BEHIND
         else
             if (h1*f <= h1*f1) exit
! USE LARGER STEP
             k1=k1+1
             if (k2 > 0) kc=kc+1
             k2=0
             if (k1 >= 20) then
                 hp=du20*hp
             else if (k1 >= 10) then
                 hp=du10*hp
             else if (k1 >= 3) then
                 hp=du03*hp
             endif
         endif
        enddo
! ----}  End of 1-D search
! ADJUST THE TRIAL STEP SIZE ----{
        dx=zero
        do i = 1,n
           dx=dx+(xopt(i)-x(i))**2
        enddo
        dx=dsqrt(dx)
        if (kg < kstore)  kg=kg+1
        if (kg >= 2) then
           do i =kg,2,-1
             nsteps(i)=nsteps(i-1)
           enddo
        endif
        d=zero
        do i = 1,n
           d=d+g0(i)*g0(i)
        enddo
        d=dsqrt(d)
        nsteps(1)=dx/(dabs(h)*d)
        kk=zero
        d=zero
        do i = 1,kg
           dd=dble(kg-i+1)
           d=d+dd
           kk=kk+nsteps(i)*dd
        enddo
        kk=kk/d
        if (kk > des) then
             if (kg == 1) then
                h=h*(kk-des+one)
             else
                h=h*dsqrt(kk-des+one)
             endif
        else if (kk < des) then
             h=h*dsqrt(kk/des)
        endif

        if (ksm) stepvanish=stepvanish+1
! ----}
! COMPUTE THE GRADIENT ----{
        if (app) then
          do j = 1,n
            if (g0(j) >= zero) then
               deltax(j)=h1*ddx
            else
               deltax(j)=-h1*ddx
            endif
          enddo
          obj=.true.
          !if (constr) then
             !call apprgrdn()
          !else
             !call apprgrdn()
          !endif
          !options(10)=options(10)+n_float
        else
          call grad(x,g,Qref,n/2,n,Kopt,f_min,f_max)
          options(11)=options(11)+one
        endif
        ng=zero
        do i = 1,n
          ng=ng+g(i)*g(i)
        enddo
        ng=dsqrt(ng)
        if (ng >= infty) then
         if (dispwarn) then
           print *,'SolvOpt error:'
           print *,'Gradient equals infinity at the starting point.'
         endif
         options(9)=-four
         goto 999
        else if (ng < ZeroGrad) then
         if (dispwarn) then
           print *,'SolvOpt warning:'
           print *,'Gradient is zero, but stopping criteria are not fulfilled.'
         endif
         ng=ZeroGrad
        endif
! Constraints:
        if (constr) then
         if (.not. FsbPnt) then
           if (ng < 1.d-2*PenCoef) then
              kless=kless+1
              if (kless >= 20) then
                 PenCoef=PenCoef/ten
                 Reset=.true.
                 kless=0
              endif
           else
              kless=0
           endif
           !if (appconstr) then
                 !do j = 1,n
                   !if (x(j) >= zero) then
                      !deltax(j)=ddx
                   !else
                      !deltax(j)=-ddx
                   !endif
                 !enddo
                 !obj=.false.
                 !call apprgrdn()
                 !options(12)=options(12)+n_float
           if (.not. appconstr) then
                 call gradc(x,gc,n/2,n,theta_min,theta_max)
                 options(13)=options(13)+one
           endif
           ngc=zero
           do i = 1,n
              ngc=ngc+gc(i)*gc(i)
           enddo
           ngc=dsqrt(ngc)
           if (ngc >= infty) then
                  if (dispwarn) then
                     print *,'SolvOpt error: GRADC returns infinite vector at the point.'
                  endif
                  options(9)=-six
                  goto 999
           else if (ngc < ZeroGrad .and. .not. appconstr) then
                  if (dispwarn) then
                     print *,'SolvOpt error: GRADC returns zero vector at an infeasible point.'
                  endif
                  options(9)=-six
                  goto 999
           endif
           do i = 1,n
             g(i)=g(i)+PenCoef*gc(i)
           enddo
           ng=zero
           do i = 1,n
              ng=ng+g(i)*g(i)
           enddo
           ng=dsqrt(ng)
           if (Reset) then
              if (dispwarn) then
                 print *,'SolvOpt warning:'
                 print *,'Re-setting due to the use of a new penalty coefficient.'
              endif
              h=h1*dx/three
              k=k-1
              nng=ng
              exit
           endif
         endif
        endif
        if (h1*f > h1*frec) then
          frec=f
          do i = 1,n
            xrec(i)=x(i)
            grec(i)=g(i)
          enddo
        endif
! ----}
       if (ng > ZeroGrad) then
        if (knorms < 10)  knorms=knorms+1
        if (knorms >= 2) then
          do i =knorms,2,-1
           gnorms(i)=gnorms(i-1)
          enddo
        endif
        gnorms(1)=ng
        nng=one
          do i = 1,knorms
            nng=nng*gnorms(i)
          enddo
        nng=nng**(one/dble(knorms))
       endif
! Norm X:
       nx=zero
       do i = 1,n
        nx=nx+x(i)*x(i)
       enddo
       nx=dsqrt(nx)

! DISPLAY THE CURRENT VALUES ----{
       if (k == ld) then
         print *, &
             'Iteration # ..... function value ..... ', &
             'Step Value ..... Gradient Norm'
         print '(5x,i5,7x,g13.5,6x,g13.5,7x,g13.5)', k,f,dx,ng
         ld=k+dispdata
       endif
!----}
! CHECK THE STOPPING CRITERIA ----{
      termflag=.true.
      if (constr) then
        if (.not. FsbPnt) termflag=.false.
      endif
      if (kcheck <= 5 .or. kcheck <= 12 .and. ng > one)termflag=.false.
      if (kc >= mxtc .or. knan)termflag=.false.
! ARGUMENT
       if (termflag) then
           ii=0
           stopping=.true.
           do i = 1,n
             if (dabs(x(i)) >= lowxbound) then
                ii=ii+1
                idx(ii)=i
                if (dabs(xopt(i)-x(i)) > options(2)*dabs(x(i))) then
                  stopping=.false.
                endif
             endif
           enddo
           if (ii == 0 .or. stopping) then
                stopping=.true.
                termx=termx+1
                d=zero
                do i = 1,n
                  d=d+(x(i)-xrec(i))**2
                enddo
                d=dsqrt(d)
! function
                if (dabs(f-frec) > detfr*dabs(f) .and. &
                  dabs(f-fopt) <= options(3)*dabs(f) .and. &
                  krerun <= 3 .and. .not. constr) then
                   stopping=.false.
                   if (ii > 0) then
                    do i = 1,ii
                     j=idx(i)
                     if (dabs(xrec(j)-x(j)) > detxr*dabs(x(j))) then
                       stopping=.true.
                       exit
                     endif
                    enddo
                   endif
                   if (stopping) then
                      if (dispwarn) then
                        print *,'SolvOpt warning:'
                        print *,'Re-run from recorded point.'
                      endif
                      ng=zero
                      do i = 1,n
                       x(i)=xrec(i)
                       g(i)=grec(i)
                       ng=ng+g(i)*g(i)
                      enddo
                      ng=dsqrt(ng)
                      f=frec
                      krerun=krerun+1
                      h=h1*dmax1(dx,detxr*nx)/dble(krerun)
                      warnno=2
                      endwarn='Result may not provide the optimum. The function apparently has many extremum points.'
                      exit
                   else
                      h=h*ten
                   endif
                else if (dabs(f-frec) > options(3)*dabs(f) .and. &
                  d < options(2)*nx .and. constr) then
                   continue
                else if (dabs(f-fopt) <= options(3)*dabs(f) .or. &
                   dabs(f) <= lowfbound .or. &
                   (dabs(f-fopt) <= options(3) .and. &
                    termx >= limxterm )) then
                  if (stopf) then
                   if (dx <= laststep) then
                    if (warnno == 1 .and. ng < dsqrt(options(3))) then
                       warnno=0
                    endif
                    if (.not. app) then
                      do i = 1,n
                       if (dabs(g(i)) <= epsnorm2) then
                         warnno=3
                         endwarn='Result may be inaccurate in the coordinates. The function is flat at the solution.'
                         exit
                       endif
                      enddo
                    endif
                    if (warnno /= 0) then
                       options(9)=-dble(warnno)-ten
                       if (dispwarn) then
                         print *,'SolvOpt: Termination warning:'
                         print *,endwarn
                         if (app) print *,'The above warning may be reasoned by inaccurate gradient approximation'
                       endif
                    else
                       options(9)=dble(k)
!! DK DK               if (dispwarn) print *,'SolvOpt: Normal termination.'
                    endif
                    goto 999
                   endif
                  else
                   stopf=.true.
                  endif
                else if (dx < powerm12*dmax1(nx,one) .and. &
                       termx >= limxterm) then
                     options(9)=-four-ten
                     if (dispwarn) then
                       print *,'SolvOpt: Termination warning:'
                       print *,'Stopping criteria are not fulfilled. The function is very steep at the solution.'
                       if (app) print *,'The above warning may be reasoned by inaccurate gradient approximation'
                       f=frec
                       do i = 1,n
                        x(i)=xrec(i)
                       enddo
                     endif
                     goto 999
                endif
           endif
       endif
! ITERATIONS LIMIT
            if (k == iterlimit) then
                options(9)=-nine
                if (dispwarn) then
                  print *,'SolvOpt warning:'
                  print *,'Iterations limit exceeded.'
                endif
                goto 999
            endif
! ----}
! ZERO GRADIENT ----{
          if (constr) then
            if (ng <= ZeroGrad) then
                if (dispwarn) then
                  print *,'SolvOpt: Termination warning:'
                  print *,'Gradient is zero, but stopping criteria are not fulfilled.'
                endif
                options(9)=-eight
                goto 999
            endif
          else
            if (ng <= ZeroGrad) then
             nzero=nzero+1
             if (dispwarn) then
               print *,'SolvOpt warning:'
               print *,'Gradient is zero, but stopping criteria are not fulfilled.'
             endif
             if (nzero >= 3) then
               options(9)=-eight
               goto 999
             endif
             do i = 1,n
               g0(i)=-h*g0(i)/two
             enddo
             do i =1,10
               do j = 1,n
                x(j)=x(j)+g0(j)
               enddo
               call fun(x,f,Qref,n/2,n,Kopt,f_min,f_max)
               options(10)=options(10)+one
               if (dabs(f) >= infty) then
                 if (dispwarn) then
                   print *,'SolvOpt error: function equals infinity at the point.'
                 endif
                 options(9)=-three
                 goto 999
               endif
               !if (app) then
                   !do j = 1,n
                     !if (g0(j) >= zero) then
                        !deltax(j)=h1*ddx
                     !else
                        !deltax(j)=-h1*ddx
                     !endif
                   !enddo
                   !obj=.true.
                   !call apprgrdn()
                   !options(10)=options(10)+n_float
               if (.not. app) then
                   call grad(x,g,Qref,n/2,n,Kopt,f_min,f_max)
                   options(11)=options(11)+one
               endif
               ng=zero
               do j = 1,n
                  ng=ng+g(j)*g(j)
               enddo
               ng=dsqrt(ng)
               if (ng >= infty) then
                    if (dispwarn) then
                      print *,'SolvOpt error:'
                      print *,'Gradient equals infinity at the starting point.'
                    endif
                    options(9)=-four
                    goto 999
               endif
               if (ng > ZeroGrad) exit
             enddo
             if (ng <= ZeroGrad) then
                if (dispwarn) then
                  print *,'SolvOpt: Termination warning:'
                  print *,'Gradient is zero, but stopping criteria are not fulfilled.'
                endif
                options(9)=-eight
                goto 999
             endif
             h=h1*dx
             exit
            endif
          endif
! ----}
! function is flat at the point ----{
          if (.not. constr .and. &
             dabs(f-fopt) < dabs(fopt)*options(3) .and. kcheck > 5 .and. ng < one) then

           ni=0
           do i = 1,n
             if (dabs(g(i)) <= epsnorm2) then
               ni=ni+1
               idx(ni)=i
             endif
           enddo
           if (ni >= 1 .and. ni <= n/2 .and. kflat <= 3) then
             kflat=kflat+1
             if (dispwarn) then
                print *,'SolvOpt warning:'
                print *,'The function is flat in certain directions.'
             endif
             warnno=1
             endwarn='Premature stop is possible. Try to re-run the routine from the obtained point.'
             do i = 1,n
               x1(i)=x(i)
             enddo
             fm=f
             do i = 1,ni
              j=idx(i)
              f2=fm
              y=x(j)
              if (y == zero) then
                x1(j)=one
              else if (dabs(y) < one) then
                x1(j)=dsign(one,y)
              else
                x1(j)=y
              endif
              do ip=1,20
               x1(j)=x1(j)/1.15d0
               call fun(x1,f1,Qref,n/2,n,Kopt,f_min,f_max)
               options(10)=options(10)+one
               if (dabs(f1) < infty) then
                 if (h1*f1 > h1*fm) then
                   y=x1(j)
                   fm=f1
                 else if (h1*f2 > h1*f1) then
                   exit
                 else if (f2 == f1) then
                   x1(j)=x1(j)/1.5d0
                 endif
                 f2=f1
               endif
              enddo
              x1(j)=y
             enddo
             if (h1*fm > h1*f) then
              !if (app) then
                !do j = 1,n
                  !deltax(j)=h1*ddx
                !enddo
                !obj=.true.
                !call apprgrdn()
                !options(10)=options(10)+n_float
              if (.not. app) then
                call grad(x1,gt,Qref,n/2,n,Kopt,f_min,f_max)
                options(11)=options(11)+one
              endif
              ngt=zero
              do i = 1,n
                ngt=ngt+gt(i)*gt(i)
              enddo
              if (ngt > epsnorm2 .and. ngt < infty) then
                if (dispwarn) print *,'Trying to recover by shifting insensitive variables.'
                do i = 1,n
                 x(i)=x1(i)
                 g(i)=gt(i)
                enddo
                ng=ngt
                f=fm
                h=h1*dx/three
                options(3)=options(3)/five
                exit
              endif   !! regular gradient
             endif   !! a better value has been found
           endif   !! function is flat
          endif   !! pre-conditions are fulfilled
! ----}
       enddo   !! iterations
      enddo   !! restart

999   continue

! deallocate working arrays:
      deallocate (idx,deltax,xx,grec,xrec,xopt,x1,z,gc,gt,g1,g0,g,B)

  end subroutine solvopt

  subroutine soptions(default)
! SOPTIONS returns the default values for the optional parameters
! used by SolvOpt.

  implicit none

  double precision default(13)

  default(1)  = -1.d0
  default(2)  = 1.d-4
  default(3)  = 1.d-6
  default(4)  = 15.d3
  default(5)  = 0.d0
  default(6)  = 1.d-8
  default(7)  = 2.5d0
  default(8)  = 1.d-12
  default(9)  = 0.d0
  default(10) = 0.d0
  default(11) = 0.d0
  default(12) = 0.d0
  default(13) = 0.d0

  end subroutine soptions

  subroutine func_objective(x,res,freq,Qref,N,Nopt)

  implicit none

  integer, intent(in) :: N,Nopt
  double precision, intent(in) :: freq,Qref
  double precision, intent(out) :: res
  double precision, dimension(1:Nopt), intent(in) :: x

  integer i
  double precision num,deno

  res = 0.d0
  do i = 1,N
    num = x(N+i)*x(N+i)*freq*Qref*(x(i)*x(i) - freq/qref)
    deno = (x(i) ** 4.) + freq*freq
    res = res + num/deno
  enddo

  end subroutine func_objective

  subroutine func_mini(x,res,Qref,N,Nopt,K,f_min,f_max)

! Nopt = 2*N : nombre de coefficients a optimiser

  use constants, only: TWO_PI

  implicit none

  integer, intent(in) :: N,Nopt,K
  double precision, intent(in) :: Qref,f_min,f_max
  double precision, intent(out) :: res
  double precision, dimension(1:Nopt), intent(in) :: x

  integer i
  double precision d,freq,aux

  res = 0.d0
  do i = 1,K
! we work in angular frequency, not frequency
    freq = TWO_PI * f_min*((f_max/f_min)**((i-1.d0)/(K-1.d0)))
    call func_objective(x,aux,freq,Qref,N,Nopt)
    d = aux - 1.d0
    res = res + d*d
  enddo

  end subroutine func_mini

  subroutine grad_func_mini(x,grad,Qref,N,Nopt,K,f_min,f_max)

  use constants, only: TWO_PI

  implicit none

  integer, intent(in) :: N,Nopt,K
  double precision, intent(in) :: Qref,f_min,f_max
  double precision, dimension(1:Nopt), intent(in) :: x
  double precision, dimension(1:Nopt), intent(out) :: grad

  integer i,l
  double precision R,temp0,temp1,temp2,temp3,tamp,aux1,aux2,aux3,aux4
  double precision, dimension(1:N) :: point,poids
  double precision, dimension(1:K) :: freq

  do i = 1,K
! we work in angular frequency, not frequency
    freq(i) = TWO_PI * f_min*((f_max/f_min)**((i-1.d0)/(K-1.d0)))
  enddo

  do l= 1,N
    point(l) = x(l)
    poids(l) = x(N+l)
  enddo

  do l= 1,N
    grad(l) = 0.d0
    grad(N+l) = 0.d0

    do i = 1,K
      call func_objective(x,R,freq(i),Qref,N,Nopt)
      temp3 = R - 1.d0
      temp0 = freq(i)*Qref

      ! derivee par rapport aux poids
      temp1 = temp0*(point(l)*point(l) - freq(i)/qref)
      temp1 = temp1*2.d0*poids(l)
      temp2 = (point(l)**4.d0) + freq(i)*freq(i)
      temp1 = temp1/temp2
      tamp = 2.d0*temp3*temp1
      grad(N+l) = grad(N+l) + tamp

      ! derivee par rapport aux points
      aux1 = -2.d0*(point(l)**5.d0) + 2.d0*point(l)*freq(i)*freq(i) + 4.d0*(point(l)**3.d0)*freq(i)/Qref
      aux3 = temp2*temp2
      aux4 = aux1/aux3
      aux4 = aux4*temp0
      aux2 = aux4*poids(l)*poids(l)
      tamp = 2.d0*temp3*aux2
      grad(l) = grad(l) + tamp
    enddo
  enddo

  end subroutine grad_func_mini

  subroutine max_residu(x,res,N,Nopt,theta_min,theta_max)

  implicit none

  integer, intent(in) :: N,Nopt
  double precision, intent(in) :: theta_min,theta_max
  double precision, intent(out) :: res
  double precision, dimension(1:Nopt), intent(in) :: x

  integer l
  double precision temp,aux

  temp = 0.d0
  res = 0.d0

  do l= 1,N
    aux = res
    temp = max(0.d0,x(l)*x(l)-(theta_max-theta_min))
    res = max(temp,aux)
  enddo

  end subroutine max_residu

  subroutine grad_max_residu(x,grad,N,Nopt,theta_min,theta_max)

  implicit none

  integer, intent(in) :: N,Nopt
  double precision, intent(in) :: theta_min,theta_max
  double precision, dimension(1:Nopt), intent(in) :: x
  double precision, dimension(1:Nopt), intent(out) :: grad

  integer l,l0
  double precision temp,res,aux,temp2
  double precision, dimension(1:N) :: point

  temp = 0.d0
  res = 0.d0

  do l= 1,N
    point(l) = x(l)
  enddo

  l0 = 1
  do l= 1,N
    aux = res
    temp = max(0.d0,point(l)*point(l) - (theta_max-theta_min))
    res = max(temp,aux)
    if (temp > aux) then
      l0 = l
    endif
  enddo

  do l= 1,N
    grad(N+l) = 0.d0
    if (l /= l0) then
      grad(l) = 0.d0
    else
      call max_residu(x,temp2,N,Nopt,theta_min,theta_max)
      if (temp2 == 0.d0) then
        grad(l0) = 0.d0
      else
        grad(l0) = 2.d0*point(l0)
      endif
    endif
  enddo

  end subroutine grad_max_residu

  subroutine nonlinear_optimization(N,Qref,f0,point,poids,f_min,f_max)

  implicit none

  logical, parameter :: USE_SOLVOPT = .true.

  integer, intent(in) :: N
  double precision, intent(in) :: Qref,f0,f_min,f_max
  double precision, dimension(1:N), intent(out) :: point,poids

  external func_mini,grad_func_mini,max_residu,grad_max_residu

  integer K,i
  logical flg,flfc,flgc
  double precision theta_min,theta_max,res
  double precision, dimension(1:2*N) :: x
  double precision, dimension(1:13) :: options

  flg = .true.
  flgc = .true.
  flfc = .true.

  K = 4*N
  theta_min = 0.d0        ! arbitrary lower limit from Bruno Lombard to make sure points never become negative
  theta_max = 1000.d0*f0  ! arbitrary upper limit from Bruno Lombard to make sure points never tend to infinity

  ! this is used as a first guess
  call classical_linear_least_squares(Qref,poids,point,N,f_min,f_max)
  if (.not. USE_SOLVOPT) return

  ! what follows is the nonlinear optimization part

  do i = 1,N
    x(i)   = sqrt(abs(point(i)) - theta_min)
    x(N+i) = sqrt(abs(poids(i)))
  enddo

  call soptions(options)
  call solvopt(2*N,x,res,func_mini,flg,grad_func_mini,options,flfc, &
      max_residu,flgc,grad_max_residu,Qref,K,theta_min,theta_max,f_min,f_max)

  do i = 1,N
    point(i) = theta_min + x(i)*x(i)
    poids(i) = x(N+i)*x(N+i)
  enddo

  end subroutine nonlinear_optimization

