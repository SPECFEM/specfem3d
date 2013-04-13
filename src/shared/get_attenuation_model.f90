!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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


  subroutine get_attenuation_model_olsen(vs_val,Q_mu,OLSEN_ATTENUATION_RATIO)

! uses scaling rule similar to Olsen et al. (2003) to determine attenuation medium
!
! returns: selected (sediment) Q_mu
!
! refers to:
!   K. B. Olsen, S. M. Day and C. R. Bradley, 2003.
!   Estimation of Q for Long-Period (>2 sec) Waves in the Los Angeles Basin
!   BSSA, 93, 2, p. 627-638

  implicit none

  include "constants.h"

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

  ! uses a simple, 2-constant model mentioned in Olsen et al. (2003)
  if( USE_SIMPLE_OLSEN ) then
    ! vs (in m/s)
    if( vs_val < 2000.0_CUSTOM_REAL ) then
      Q_mu = 0.02 * vs_val
    else
      Q_mu = 0.1 * vs_val
    endif
  endif

  ! uses discrete values in sediment range
  if( USE_DISCRETE_OLSEN ) then
    int_Q_mu = 10 * nint(Q_mu / 10.)

    if(int_Q_mu < 40) int_Q_mu = 40
    if(int_Q_mu > 150) int_Q_mu = 150

    if(int_Q_mu == 40) then
      Q_mu = 40.0d0
    else if(int_Q_mu == 50) then
      Q_mu = 50.0d0
    else if(int_Q_mu == 60) then
      Q_mu = 60.0d0
    else if(int_Q_mu == 70) then
      Q_mu = 70.0d0
    else if(int_Q_mu == 80) then
      Q_mu = 80.0d0
    else if(int_Q_mu == 90) then
      Q_mu = 90.0d0
    else if(int_Q_mu == 100) then
      Q_mu = 100.0d0
    else if(int_Q_mu == 110) then
      Q_mu = 110.0d0
    else if(int_Q_mu == 120) then
      Q_mu = 120.0d0
    else if(int_Q_mu == 130) then
      Q_mu = 130.0d0
    else if(int_Q_mu == 140) then
      Q_mu = 140.0d0
    else if(int_Q_mu == 150) then
      Q_mu = 150.0d0
    else
      stop 'incorrect attenuation coefficient'
    endif
  endif

  ! limits Q_mu value range
  if( Q_mu < 1.0d0 ) Q_mu = 1.0d0
  if( Q_mu > ATTENUATION_COMP_MAXIMUM ) Q_mu = ATTENUATION_COMP_MAXIMUM


  end subroutine get_attenuation_model_olsen

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_model(myrank,nspec,USE_OLSEN_ATTENUATION,OLSEN_ATTENUATION_RATIO, &
                                  mustore,rho_vs,kappastore,rho_vp,qkappa_attenuation_store,qmu_attenuation_store, &
                                  ispec_is_elastic,min_resolved_period,prname,FULL_ATTENUATION_SOLID)

! precalculates attenuation arrays and stores arrays into files

  implicit none

  include "constants.h"

  double precision :: OLSEN_ATTENUATION_RATIO
  integer :: myrank,nspec
  logical :: USE_OLSEN_ATTENUATION

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: mustore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rho_vs
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: qkappa_attenuation_store
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: qmu_attenuation_store

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: kappastore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rho_vp

  logical, dimension(nspec) :: ispec_is_elastic
  real(kind=CUSTOM_REAL) :: min_resolved_period
  character(len=256) :: prname

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: factor_common
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: one_minus_sum_beta_kappa
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: scale_factor, scale_factor_kappa
  double precision, dimension(N_SLS) :: tau_sigma_dble,beta_dble,beta_dble_kappa
  double precision factor_scale_dble,one_minus_sum_beta_dble,&
                   factor_scale_dble_kappa,one_minus_sum_beta_dble_kappa
  double precision :: Q_mu,Q_kappa,Q_p,Q_s
  double precision :: f_c_source
  double precision :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tau_sigma
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tauinv
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: beta,beta_kappa
  real(kind=CUSTOM_REAL):: vs_val,vp_val
  integer :: i,j,k,ispec,ier
  double precision :: qmin,qmax,qmin_all,qmax_all
  logical :: FULL_ATTENUATION_SOLID

  ! initializes arrays
  allocate(one_minus_sum_beta(NGLLX,NGLLY,NGLLZ,nspec), &
          factor_common(N_SLS,NGLLX,NGLLY,NGLLZ,nspec), &
          scale_factor(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocation attenuation arrays')

  if(FULL_ATTENUATION_SOLID)then
    allocate(one_minus_sum_beta_kappa(NGLLX,NGLLY,NGLLZ,nspec), &
            factor_common_kappa(N_SLS,NGLLX,NGLLY,NGLLZ,nspec), &
            scale_factor_kappa(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocation attenuation arrays')
  else
    allocate(one_minus_sum_beta_kappa(NGLLX,NGLLY,NGLLZ,1), &
            factor_common_kappa(N_SLS,NGLLX,NGLLY,NGLLZ,1), &
            scale_factor_kappa(NGLLX,NGLLY,NGLLZ,1),stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocation attenuation arrays')
  endif

  one_minus_sum_beta(:,:,:,:) = 1._CUSTOM_REAL
  factor_common(:,:,:,:,:) = 1._CUSTOM_REAL
  scale_factor(:,:,:,:) = 1._CUSTOM_REAL

  one_minus_sum_beta_kappa(:,:,:,:) = 1._CUSTOM_REAL
  factor_common_kappa(:,:,:,:,:) = 1._CUSTOM_REAL
  scale_factor_kappa(:,:,:,:) = 1._CUSTOM_REAL

  ! gets stress relaxation times tau_sigma, i.e.
  ! precalculates tau_sigma depending on period band (constant for all Q_mu), and
  ! determines central frequency f_c_source of attenuation period band
  call get_attenuation_constants(min_resolved_period,tau_sigma_dble,&
              f_c_source,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! user output
  if( myrank == 0 ) then
    write(IMAIN,*)
    write(IMAIN,*) "attenuation: "
    write(IMAIN,*) "  reference period (s)   : ",sngl(1.0/ATTENUATION_f0_REFERENCE), &
                  " frequency: ",sngl(ATTENUATION_f0_REFERENCE)
    write(IMAIN,*) "  period band min/max (s): ",sngl(MIN_ATTENUATION_PERIOD),sngl(MAX_ATTENUATION_PERIOD)
    write(IMAIN,*) "  central period (s)     : ",sngl(1.0/f_c_source), &
                  " frequency: ",sngl(f_c_source)
  endif

  ! determines inverse of tau_sigma
  if(CUSTOM_REAL == SIZE_REAL) then
    tau_sigma(:) = sngl(tau_sigma_dble(:))
  else
    tau_sigma(:) = tau_sigma_dble(:)
  endif
  ! precalculates the negative inverse of tau_sigma
  tauinv(:) = - 1._CUSTOM_REAL / tau_sigma(:)

  ! precalculates factors for shear modulus scaling according to attenuation model
  qmin = HUGEVAL
  qmax = 0.0
  do ispec = 1,nspec

    ! skips non elastic elements
    if( ispec_is_elastic(ispec) .eqv. .false. ) cycle

    ! determines attenuation factors for each GLL point
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          ! shear moduli attenuation
          ! gets Q_mu value
          if(USE_OLSEN_ATTENUATION) then
            ! use scaling rule similar to Olsen et al. (2003)
            vs_val = mustore(i,j,k,ispec) / rho_vs(i,j,k,ispec)
            call get_attenuation_model_olsen(vs_val,Q_mu,OLSEN_ATTENUATION_RATIO)
            if(FULL_ATTENUATION_SOLID)then
              vp_val = (kappastore(i,j,k,ispec) + 2.0d0 * mustore(i,j,k,ispec) / 3.0d0) / rho_vp(i,j,k,ispec)
              Q_s = Q_mu
              Q_p = 1.5d0 * Q_s
              Q_kappa = 1.0d0 / ((1.0/Q_p - 4.0d0/3.0d0*(vp_val/vs_val)**2*(1.d0/Q_mu)) /(1.0d0 - 4.0d0/3.0d0*(vp_val/vs_val)**2))
              if( Q_kappa < 1.0d0 ) Q_kappa = 1.0d0
              if( Q_kappa > ATTENUATION_COMP_MAXIMUM ) Q_kappa = ATTENUATION_COMP_MAXIMUM
            endif
          else
            ! takes Q set in (CUBIT) mesh
            Q_mu = qmu_attenuation_store(i,j,k,ispec)
            Q_kappa = qkappa_attenuation_store(i,j,k,ispec)

            ! attenuation zero
            if( Q_mu <= 1.e-5 ) cycle
            if( Q_kappa <= 1.e-5 ) cycle

            ! limits Q
            if( Q_mu < 1.0d0 ) Q_mu = 1.0d0
            if( Q_mu > ATTENUATION_COMP_MAXIMUM ) Q_mu = ATTENUATION_COMP_MAXIMUM

            if( Q_kappa < 1.0d0 ) Q_kappa = 1.0d0
            if( Q_kappa > ATTENUATION_COMP_MAXIMUM ) Q_kappa = ATTENUATION_COMP_MAXIMUM

          endif

          ! statistics on Q_mu
          if( Q_mu < qmin ) qmin = Q_mu
          if( Q_mu > qmax ) qmax = Q_mu

          ! gets beta, on_minus_sum_beta and factor_scale
          ! based on calculation of strain relaxation times tau_eps
          call get_attenuation_factors(myrank,Q_mu,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
                            f_c_source,tau_sigma_dble, &
                            beta_dble,one_minus_sum_beta_dble,factor_scale_dble, &
                            Q_kappa,beta_dble_kappa,one_minus_sum_beta_dble_kappa,factor_scale_dble_kappa,FULL_ATTENUATION_SOLID)

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

          if(FULL_ATTENUATION_SOLID)then
            one_minus_sum_beta_kappa(i,j,k,ispec) = one_minus_sum_beta_dble_kappa
            beta_kappa(:) = beta_dble_kappa(:)
            factor_common_kappa(:,i,j,k,ispec) = beta_kappa(:) * tauinv(:)

            ! stores scale factor for mu moduli
            scale_factor_kappa(i,j,k,ispec) = factor_scale_dble_kappa
          endif

        enddo
      enddo
    enddo
  enddo

  ! stores attenuation arrays into files
  open(unit=27, file=prname(1:len_trim(prname))//'attenuation.bin', &
        status='unknown',action='write',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error: could not open ',prname(1:len_trim(prname))//'attenuation.bin'
    call exit_mpi(myrank,'error opening attenuation.bin file')
  endif
  write(27) nspec
  write(27) one_minus_sum_beta
  write(27) factor_common
  write(27) scale_factor

  if(FULL_ATTENUATION_SOLID)then
    write(27) one_minus_sum_beta_kappa
    write(27) factor_common_kappa
    write(27) scale_factor_kappa
  endif

  close(27)

  deallocate(one_minus_sum_beta,factor_common,scale_factor)

  if(FULL_ATTENUATION_SOLID)then
    deallocate(one_minus_sum_beta_kappa,factor_common_kappa,scale_factor_kappa)
  endif

  ! statistics
  call min_all_dp(qmin,qmin_all)
  call max_all_dp(qmax,qmax_all)
  ! user output
  if( myrank == 0 ) then
    write(IMAIN,*) "  Q_mu min/max           : ",sngl(qmin_all),sngl(qmax_all)
    write(IMAIN,*)
  endif

  end subroutine get_attenuation_model

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_attenuation_memory_values(tau_s, deltat, alphaval,betaval,gammaval)

! returns: runge-kutta scheme terms alphaval, betaval and gammaval

  implicit none

  include 'constants.h'

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

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) :: min_resolved_period
  double precision, dimension(N_SLS) :: tau_sigma
  double precision :: f_c_source
  double precision :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD

  ! local parameters
  real(kind=CUSTOM_REAL)  :: min_period

  ! determines min/max periods for attenuation band based on minimum resolved period of mesh
  min_period = 0.99 * min_resolved_period ! uses a small margin
  ! debug for comparison with fix values from above
  !min_resolved_period = 0.943
  call get_attenuation_periods(min_period,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  !  sets up stress relaxation times tau_sigma,
  ! equally spaced based on number of standard linear solids and period band
  call get_attenuation_tau_sigma(tau_sigma,N_SLS,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! sets up central frequency
  ! logarithmic mean of frequency interval of absorption band
  call get_attenuation_source_freq(f_c_source,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! debug
  !f_c_source = 0.141421d0
  !tau_sigma(1) =  7.957747154594766669788441504352d0
  !tau_sigma(2) =  1.125395395196382652969191440206d0
  !tau_sigma(3) =  0.159154943091895345608222100964d0

  end subroutine get_attenuation_constants


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_factors(myrank,Q_mu,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
                              f_c_source,tau_sigma, &
                              beta,one_minus_sum_beta,factor_scale, &
                              Q_kappa,beta_kappa,one_minus_sum_beta_kappa,factor_scale_kappa,FULL_ATTENUATION_SOLID)

! returns: attenuation mechanisms beta,one_minus_sum_beta,factor_scale

! variable frequency range
! variable period range
! variable central logarithmic frequency

! in the future when more memory is available on computers
! it would be more accurate to use four mechanisms instead of three


  implicit none

  include "constants.h"

  integer:: myrank
  double precision :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD
  double precision :: f_c_source,Q_mu,Q_kappa
  double precision, dimension(N_SLS) :: tau_sigma
  double precision, dimension(N_SLS) :: beta,beta_kappa
  double precision :: one_minus_sum_beta,one_minus_sum_beta_kappa
  double precision :: factor_scale,factor_scale_kappa
  logical :: FULL_ATTENUATION_SOLID
  ! local parameters
  double precision, dimension(N_SLS) :: tau_eps,tau_eps_kappa


  ! determines tau_eps for Q_mu
  call get_attenuation_tau_eps(Q_mu,tau_sigma,tau_eps, &
                              MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! determines one_minus_sum_beta
  call get_attenuation_property_values(tau_sigma,tau_eps,beta,one_minus_sum_beta)

  ! determines the "scale factor"
  call get_attenuation_scale_factor(myrank,f_c_source,tau_eps,tau_sigma,Q_mu,factor_scale)

  if(FULL_ATTENUATION_SOLID)then
    ! determines tau_eps for Q_kappa
    call get_attenuation_tau_eps(Q_kappa,tau_sigma,tau_eps_kappa, &
                                MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

    ! determines one_minus_sum_beta
    call get_attenuation_property_values(tau_sigma,tau_eps_kappa,beta_kappa,one_minus_sum_beta_kappa)

    ! determines the "scale factor"
    call get_attenuation_scale_factor(myrank,f_c_source,tau_eps_kappa,tau_sigma,Q_kappa,factor_scale_kappa)
  endif

  end subroutine get_attenuation_factors

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_attenuation_property_values(tau_s, tau_eps, beta, one_minus_sum_beta)

! coefficients useful for calculation between relaxed and unrelaxed moduli
!
! returns: coefficients beta, one_minus_sum_beta

  implicit none

  include "constants.h"

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

  subroutine get_attenuation_scale_factor(myrank, f_c_source, tau_eps, tau_sigma, Q_mu, scale_factor)

! returns: physical dispersion scaling factor scale_factor

  implicit none

  include "constants.h"

  integer :: myrank
  double precision :: scale_factor, Q_mu, f_c_source
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
  factor_scale_mu0 = ONE + TWO * log(f_c_source / ATTENUATION_f0_REFERENCE ) / (PI * Q_mu)

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
  if(scale_factor < 0.7 .or. scale_factor > 1.3) then
     write(*,*) "error : in get_attenuation_scale_factor() "
     write(*,*) "  scale factor: ", scale_factor, " should be between 0.7 and 1.3"
     write(*,*) "  please check your reference frequency ATTENUATION_f0_REFERENCE in constants.h"
     call exit_MPI(myrank,'incorrect correction factor in attenuation model')
  endif

  end subroutine get_attenuation_scale_factor



!
!-------------------------------------------------------------------------------------------------
!

!compare: auto_ner.f90, GLOBE package

  subroutine get_attenuation_periods(min_resolved_period,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

! determines min/max periods for attenuation based upon mininum resolved period of mesh

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL),intent(in) :: min_resolved_period
  double precision,intent(out) :: MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD

  ! local parameters
  double precision :: THETA(5)

  ! checks number of standard linear solids
  if(N_SLS < 2 .OR. N_SLS > 5) then
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
  double precision :: dexp
  integer :: i
  double precision, parameter :: PI = 3.14159265358979d0

  ! min/max frequencies
  f1 = 1.0d0 / max_period
  f2 = 1.0d0 / min_period

  ! logarithms
  exp1 = log10(f1)
  exp2 = log10(f2)

  ! equally spaced in log10 frequency
  dexp = (exp2-exp1) / ((nsls*1.0d0) - 1)
  do i = 1,nsls
     tau_s(i) = 1.0 / (PI * 2.0d0 * 10**(exp1 + (i - 1)* 1.0d0 *dexp))
  enddo

  end subroutine get_attenuation_tau_sigma

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_attenuation_source_freq(f_c_source,min_period,max_period)

! Determines the Source Frequency

  implicit none

  double precision,intent(out) :: f_c_source
  double precision,intent(in) :: min_period, max_period

  ! local parameters
  double precision f1, f2,T_c_source

  ! min/max frequencies
  f1 = 1.0d0 / max_period
  f2 = 1.0d0 / min_period

  T_c_source =  1.0e+03 * 10.0d0**(0.5 * (log10(f1) + log10(f2)))
  ! central frequency
  f_c_source = T_c_source / 1000.0d0

  end subroutine get_attenuation_source_freq

!--------------------------------------------------------------------------------------------------
!
!  This portion of the SPECFEM3D Code was written by:
!  Brian Savage while at
!     California Institute of Technology
!     Department of Terrestrial Magnetism / Carnegie Institute of Washington
!     Univeristy of Rhode Island
!
!  <savage@uri.edu>.
!  <savage13@gps.caltech.edu>
!  <savage13@dtm.ciw.edu>
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

  subroutine get_attenuation_tau_eps(Qmu_in,tau_s,tau_eps,min_period,max_period)

! includes min_period, max_period, and N_SLS
!
! returns: determines strain relaxation times tau_eps

  implicit none

  include 'constants.h'

! model_attenuation_variables
!...

  double precision :: Qmu_in
  double precision, dimension(N_SLS) :: tau_s, tau_eps
  double precision :: min_period,max_period

  ! local parameters
  integer :: rw
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

  ! READ
  rw = 1
  call model_attenuation_storage(Qmu_in, tau_eps, rw, AM_S)
  if(rw > 0) return

  call attenuation_invert_by_simplex(min_period, max_period, N_SLS, Qmu_in, tau_s, tau_eps, AS_V)

  ! WRITE
  rw = -1
  call model_attenuation_storage(Qmu_in, tau_eps, rw, AM_S)

  end subroutine get_attenuation_tau_eps

!
!-------------------------------------------------------------------------------------------------
!


  subroutine model_attenuation_storage(Qmu, tau_eps, rw, AM_S)

  implicit none
  include 'constants.h'

! model_attenuation_storage_var
  type model_attenuation_storage_var
    sequence
    double precision, dimension(:,:), pointer :: tau_eps_storage
    double precision, dimension(:), pointer :: Qmu_storage
    integer Q_resolution
    integer Q_max
  end type model_attenuation_storage_var

  type (model_attenuation_storage_var) AM_S
! model_attenuation_storage_var

  double precision Qmu, Qmu_new
  double precision, dimension(N_SLS) :: tau_eps
  integer rw

  integer Qtmp
  integer, save :: first_time_called = 1
  double precision, parameter :: ZERO_TOL = 1.e-5
  integer ier

  if(first_time_called == 1) then
     first_time_called       = 0
     AM_S%Q_resolution = 10**ATTENUATION_COMP_RESOLUTION
     AM_S%Q_max        = ATTENUATION_COMP_MAXIMUM
     Qtmp         = AM_S%Q_resolution * AM_S%Q_max

     allocate(AM_S%tau_eps_storage(N_SLS, Qtmp), &
             AM_S%Qmu_storage(Qtmp),stat=ier)
     if( ier /= 0 ) stop 'error allocating arrays for attenuation storage'
     AM_S%Qmu_storage(:) = -1
  endif

  if(Qmu < 0.0d0 .OR. Qmu > AM_S%Q_max) then
     write(IMAIN,*) 'Error attenuation_storage()'
     write(IMAIN,*) 'Attenuation Value out of Range: ', Qmu
     write(IMAIN,*) 'Attenuation Value out of Range: Min, Max ', 0, AM_S%Q_max
     call exit_MPI(0, 'Attenuation Value out of Range')
  endif

  if(rw > 0 .AND. Qmu <= ZERO_TOL) then
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
  Qtmp    = Qmu * dble(AM_S%Q_resolution)

  ! rounds to corresponding double value:
  ! e.g. Qmu_new = dble( 1503 ) / dble(10) = 150.30
  ! but Qmu_new is not used any further...
  Qmu_new = dble(Qtmp) / dble(AM_S%Q_resolution)

  if(rw > 0) then
     ! READ
     if(AM_S%Qmu_storage(Qtmp) > 0) then
        ! READ SUCCESSFUL
        tau_eps(:)   = AM_S%tau_eps_storage(:, Qtmp)
        Qmu        = AM_S%Qmu_storage(Qtmp)
        rw = 1
     else
        ! READ NOT SUCCESSFUL
        rw = -1
     endif
  else
     ! WRITE SUCCESSFUL
     AM_S%tau_eps_storage(:,Qtmp)    = tau_eps(:)
     AM_S%Qmu_storage(Qtmp)        = Qmu
     rw = 1
  endif

  end subroutine model_attenuation_storage

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_invert_by_simplex(t2, t1, n, Q_real, tau_s, tau_eps, AS_V)

  implicit none

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
  ! attenuation_simplex_variables

  ! Input / Output
  double precision  t1, t2
  double precision  Q_real
!  double precision  omega_not
  integer  n
  double precision, dimension(n)   :: tau_s, tau_eps

  ! Internal
  integer i, iterations, err,prnt
  double precision f1, f2, exp1,exp2, min_value !, dexp
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

!  if(f2 < f1 .OR. Q_real < 0.0d0 .OR. n < 1) then
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
!  dexp = (exp2-exp1) / ((n*1.0d0) - 1)
!  do i = 1,n
!     tau_s(i) = 1.0 / (PI * 2.0d0 * 10**(exp1 + (i - 1)* 1.0d0 *dexp))
!  enddo


  ! Shove the paramters into the module
  call attenuation_simplex_setup(nf,n,f,Q_real,tau_s,AS_V)

  ! Set the Tau_epsilon (tau_eps) to an initial value at omega*tau = 1
  ! tan_delta = 1/Q = (tau_eps - tau_s)/(2 * sqrt(tau e*tau_s))
  !    if we assume tau_eps =~ tau_s
  !    we get the equation below
  do i = 1,n
     tau_eps(i) = tau_s(i) + (tau_s(i) * 2.0d0/Q_real)
  enddo

  ! Run a simplex search to determine the optimum values of tau_eps
  call fminsearch(attenuation_eval, tau_eps, n, iterations, min_value, prnt, err,AS_V)
  if(err > 0) then
     write(*,*)'Search did not converge for an attenuation of ', Q_real
     write(*,*)'    Iterations: ', iterations
     write(*,*)'    Min Value:  ', min_value
     write(*,*)'    Aborting program'
     call exit_MPI(0,'attenuation_simplex: Search for Strain relaxation times did not converge')
  endif

  !deallocate(f)

  call attenuation_simplex_finish(AS_V)

  end subroutine attenuation_invert_by_simplex

!
!-------------------------------------------------------------------------------------------------
!


  subroutine attenuation_simplex_setup(nf_in,nsls_in,f_in,Q_in,tau_s_in,AS_V)

!   - Inserts necessary parameters into the module attenuation_simplex_variables
!   - See module for explaination

  implicit none

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
  ! attenuation_simplex_variables

  integer nf_in, nsls_in
  double precision Q_in
  double precision, dimension(nf_in)   :: f_in
  double precision, dimension(nsls_in) :: tau_s_in
  integer ier

  allocate(AS_V%f(nf_in), &
          AS_V%tau_s(nsls_in),stat=ier)
  if( ier /= 0 ) stop 'error allocating arrays for attenuation simplex'

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


  double precision function attenuation_eval(Xin,AS_V)

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
!    See atteunation_simplex_setup
!

  implicit none

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
  ! attenuation_simplex_variables

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
     xi = sqrt(( ( (tan_delta(i) - AS_V%iQ) ** 2 ) / iQ2 ))
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


  subroutine fminsearch(funk, x, n, itercount, tolf, prnt, err, AS_V)

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
!                 1 => Function Evaluation exceeded limit
!                 2 => Iterations exceeded limit
!
!     See Matlab fminsearch

  implicit none

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
! attenuation_simplex_variables

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


  if(itercount > 0) then
     maxiter = itercount
  else
     maxiter = 200 * n
  endif
  itercount = 0
  maxfun  = 200 * n

  if(tolf > 0.0d0) then
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

  fv(1) = funk(xin,AS_V)

  usual_delta = 0.05
  zero_term_delta = 0.00025

  do j = 1,n
     y = xin
     if(y(j) /= 0.0d0) then
        y(j) = (1.0d0 + usual_delta) * y(j)
     else
        y(j) = zero_term_delta
     endif
     v(:,j+1) = y
     x(:) = y
     fv(j+1) = funk(x,AS_V)
  enddo

  call qsort_local(fv,n+1,place)

  do i = 1,n+1
     vtmp(:,i) = v(:,place(i))
  enddo
  v = vtmp

  how = initial
  itercount = 1
  func_evals = n+1
  if(prnt == 3) then
     write(*,*)'Iterations   Funk Evals   Value How'
     write(*,*)itercount, func_evals, fv(1), how
  endif
  if(prnt == 4) then
     write(*,*)'How: ',how
     write(*,*)'V: ', v
     write(*,*)'fv: ',fv
     write(*,*)'evals: ',func_evals
  endif

  do while (func_evals < maxfun .AND. itercount < maxiter)

     if(max_size_simplex(v,n) <= tolx .AND. &
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
     fxr = funk(x,AS_V)
     func_evals = func_evals + 1
     if (fxr < fv(1)) then
        ! Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,n+1)
        x = xe
        fxe = funk(x,AS_V)
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
              fxc = funk(x,AS_V)
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
              fxcc = funk(x,AS_V)
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
                 fv(j) = funk(x,AS_V)
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

  if(func_evals > maxfun) then
     write(*,*)'function evaluations exceeded prescribed limit', maxfun
     err = 1
  endif
  if(itercount > maxiter) then
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
!      Returns:
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
     if(z > m) then
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
!     Returns:
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
        if(z > m) then
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
!         Sorted Indicies of vecotr X
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
        if(X(k+1) < X(k)) then
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

  subroutine attenuation_simplex_finish(AS_V)

  implicit none

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
  ! attenuation_simplex_variables

  deallocate(AS_V%f)
  deallocate(AS_V%tau_s)

  end subroutine attenuation_simplex_finish
