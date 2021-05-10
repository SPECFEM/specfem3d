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


  subroutine prepare_attenuation()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  use shared_parameters, only: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD

  implicit none

  ! local parameters
  double precision, dimension(N_SLS) :: tau_sigma_dble
  double precision :: f_c_source
  real(kind=CUSTOM_REAL):: scale_factorl
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: scale_factor,scale_factor_kappa

  ! checks if attenuation is on and anything to do
  if (.not. ATTENUATION) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing attenuation"
    write(IMAIN,*) "  The code uses a constant Q quality factor, but approximated"
    write(IMAIN,*) "  based on a series of Zener standard linear solids (SLS)."
    write(IMAIN,*) "  Approximation is performed in the following frequency band:"
    write(IMAIN,*)
    write(IMAIN,*) "  number of SLS bodies: ",N_SLS
    write(IMAIN,*)
    write(IMAIN,*) "  Reference frequency of anelastic model (Hz): ",sngl(ATTENUATION_f0_REFERENCE)
    write(IMAIN,*) "                                   period (s): ",sngl(1.0/ATTENUATION_f0_REFERENCE)
    call flush_IMAIN()
  endif

  ! if attenuation is on, shift shear moduli to center frequency of absorption period band, i.e.
  ! rescale mu to average (central) frequency for attenuation

  ! initializes arrays
  factor_common(:,:,:,:,:) = 1._CUSTOM_REAL

  allocate( scale_factor(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2232')
  if (ier /= 0) call exit_mpi(myrank,'error allocation scale_factor')
  scale_factor(:,:,:,:) = 1._CUSTOM_REAL

  factor_common_kappa(:,:,:,:,:) = 1._CUSTOM_REAL
  allocate( scale_factor_kappa(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2233')
  if (ier /= 0) call exit_mpi(myrank,'error allocation scale_factor_kappa')
  scale_factor_kappa(:,:,:,:) = 1._CUSTOM_REAL

  ! reads in attenuation arrays
  call create_name_database(prname,myrank,LOCAL_PATH)
  if (I_should_read_the_database) then
      open(unit=27, file=prname(1:len_trim(prname))//'attenuation.bin', status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
          print *,'error: could not open ',prname(1:len_trim(prname))//'attenuation.bin'
          call exit_mpi(myrank,'error opening attenuation.bin file')
      endif
  endif

  if (I_should_read_the_database) then
      read(27) ispec
      if (ispec /= NSPEC_ATTENUATION_AB) then
          close(27)
          print *,'error: attenuation file array ',ispec,'should be ',NSPEC_ATTENUATION_AB
          call exit_mpi(myrank,'error attenuation array dimensions, please recompile and rerun generate_databases')
      endif
      read(27) factor_common
      read(27) scale_factor

      read(27) factor_common_kappa
      read(27) scale_factor_kappa

      close(27)
  endif

  ! broadcasts
  call bcast_all_i_for_database(ispec, 1)
  if (size(factor_common) > 0) &
    call bcast_all_cr_for_database(factor_common(1,1,1,1,1), size(factor_common))
  if (size(scale_factor) > 0) &
    call bcast_all_cr_for_database(scale_factor(1,1,1,1), size(scale_factor))
  call bcast_all_cr_for_database(factor_common_kappa(1,1,1,1,1), size(factor_common_kappa))
  call bcast_all_cr_for_database(scale_factor_kappa(1,1,1,1), size(scale_factor_kappa))

  ! gets stress relaxation times tau_sigma, i.e.
  ! precalculates tau_sigma depending on period band (constant for all Q_mu), and
  ! determines central frequency f_c_source of attenuation period band
  call get_attenuation_constants(min_resolved_period,tau_sigma_dble, &
                                 f_c_source,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

  ! checks
  if (f_c_source <= 0.d0) call exit_MPI(myrank,'Error: invalid attenuation center frequency, cannot be zero or negative')

  ! stores center frequency as shared parameter
  ATT_F_C_SOURCE = f_c_source

  ! determines alphaval,betaval,gammaval for runge-kutta scheme
  tau_sigma(:) = real(tau_sigma_dble(:),kind=CUSTOM_REAL)

  ! shifts shear moduli
  do ispec = 1,NSPEC_AB

    ! skips non elastic elements
    if (ispec_is_elastic(ispec) .eqv. .false.) cycle

    ! determines attenuation factors for each GLL point
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          ! scales only mu moduli
          scale_factorl = scale_factor(i,j,k,ispec)
          mustore(i,j,k,ispec) = mustore(i,j,k,ispec) * scale_factorl

          ! scales kappa moduli
          scale_factorl = scale_factor_kappa(i,j,k,ispec)
          kappastore(i,j,k,ispec) = kappastore(i,j,k,ispec) * scale_factorl

        enddo
      enddo
    enddo
  enddo

  deallocate(scale_factor)
  deallocate(scale_factor_kappa)

  ! precompute Runge-Kutta coefficients
  call get_attenuation_memory_values(tau_sigma,deltat,alphaval,betaval,gammaval)

  ! attenuation backward memories
  if (SIMULATION_TYPE == 3) then
    ! precompute Runge-Kutta coefficients if attenuation
    call get_attenuation_memory_values(tau_sigma,b_deltat,b_alphaval,b_betaval,b_gammaval)
  endif

  ! just to be sure to use the forward alpha/beta/gamma-val arrays
  ! (in principle, should be already the same for undo_att as b_deltat is set to delta_t)
  if (UNDO_ATTENUATION_AND_OR_PML) then
    b_alphaval = alphaval
    b_betaval = betaval
    b_gammaval = gammaval
  endif

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

  if (SIMULATION_TYPE == 3) then
    ! memory variables if attenuation
    if (ELASTIC_SIMULATION) then
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
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  Attenuation frequency band min/max (Hz):",sngl(1.0/MAX_ATTENUATION_PERIOD), &
                                                            '/',sngl(1.0/MIN_ATTENUATION_PERIOD)
    write(IMAIN,*) "              period band    min/max (s) :",sngl(MIN_ATTENUATION_PERIOD), &
                                                            '/',sngl(MAX_ATTENUATION_PERIOD)
    write(IMAIN,*) "  Logarithmic center frequency (Hz):",sngl(ATT_F_C_SOURCE)
    write(IMAIN,*) "                     period     (s):",sngl(1.0/ATT_F_C_SOURCE)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_attenuation
