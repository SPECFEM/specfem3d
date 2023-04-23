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


!==============================================================================
! \file save_adjoint_kernels
!
! TODO
! * Better doxygen documentation.
!==============================================================================


!> Save kernels.

  subroutine save_adjoint_kernels()

  use constants, only: SAVE_WEIGHTS
  use shared_parameters, only: ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION

  use specfem_par, only: LOCAL_PATH, NSPEC_AB, ADIOS_FOR_KERNELS, NOISE_TOMOGRAPHY, &
                         APPROXIMATE_HESS_KL, ANISOTROPIC_KL, SAVE_MOHO_MESH

  use specfem_par_noise, only: sigma_kl

  implicit none

  if (ADIOS_FOR_KERNELS) then
    call define_kernel_adios_variables()
  endif

  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    call save_kernels_acoustic()
  endif

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    if (ANISOTROPIC_KL) then
      call save_kernels_elastic_aniso()
    else
      call save_kernels_elastic_iso()
    endif
  endif

  ! poroelastic domains
  if (POROELASTIC_SIMULATION) then
    call save_kernels_poroelastic()
  endif

  ! save weights for volume integration,
  ! in order to benchmark the kernels with analytical expressions
  if (SAVE_WEIGHTS) then
    call save_weights_kernel()
  endif

  ! moho boundary kernels
  if (SAVE_MOHO_MESH) then
    call save_kernels_moho()
  endif

  ! for noise simulations --- noise strength kernel
  if (NOISE_TOMOGRAPHY == 3) then
    call save_kernels_strength_noise(LOCAL_PATH,sigma_kl,NSPEC_AB)
  endif

  ! for preconditioner
  if (APPROXIMATE_HESS_KL) then
    call save_kernels_Hessian()
  endif

  if (ADIOS_FOR_KERNELS) then
    call perform_write_adios_kernels()
  endif

  end subroutine save_adjoint_kernels

!
!-------------------------------------------------------------------------------------------------
!

!> Save weights for volume integration,
!! in order to benchmark the kernels with analytical expressions.

  subroutine save_weights_kernel()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer:: ispec,i,j,k,ier,ispec_irreg
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: weights_kernel
  real(kind=CUSTOM_REAL) :: jacobianl

  allocate(weights_kernel(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2253')
  if (ier /= 0) stop 'error allocating array weights_kernel'
  weights_kernel(:,:,:,:) = 0.0_CUSTOM_REAL

  do ispec = 1, NSPEC_AB
    ispec_irreg = irregular_element_number(ispec)
    if (ispec_irreg == 0) jacobianl = jacobian_regular
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)
          weights_kernel(i,j,k,ispec) = wxgll(i) * wygll(j) * wzgll(k) * jacobianl
        enddo ! i
      enddo ! j
    enddo ! k
  enddo ! ispec

  open(unit=IOUT,file=prname(1:len_trim(prname))//'weights_kernel.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file weights_kernel.bin'
  write(IOUT) weights_kernel
  close(IOUT)

  deallocate(weights_kernel,stat=ier)
  if (ier /= 0) stop 'error allocating array weights_kernel'

  end subroutine save_weights_kernel

!
!-------------------------------------------------------------------------------------------------
!

!> Save acoustic related kernels

  subroutine save_kernels_acoustic()

  use specfem_par
  use specfem_par_acoustic

  implicit none

  ! local parameters
  integer:: ispec,i,j,k,ier
  real(kind=CUSTOM_REAL) :: kappa_invl,rho_invl
  real(kind=CUSTOM_REAL) :: rho_ac_max,rhop_ac_max,kappa_ac_max,alpha_ac_max

  ! note: from Luo et al. (2013), kernels are given for absolute perturbations d(1/rho) and d(1/kappa)
  !       and only change from Peter et al. (2011) for acoustic kernels:
  !         K_kappa = - int_0^T [ phi^adj \partial_t^2 phi ] dt     see (A-27)
  !         K_rho   = - int_0^T [ grad(phi^adj) * grad(phi) ] dt        (A-28)
  !
  !       since we output relative perturbation kernels for elastic domains, we make here also use of the variation
  !         d(1/rho) = - 1/rho^2 d(rho) = - 1/rho d(rho)/rho = - 1/rho dln(rho)
  !
  !       to obtain relative kernels, we start from eq. (A-24)
  !         dChi = int_V [ d(1/kappa) K_kappa + d(1/rho) K_rho ] d^3x (+ elastic kernels)
  !              = int_V [ (-1/kappa K_kappa) dln(kappa) + (- 1/rho K_rho) dln(rho)
  !
  !       and see that relative perturbation kernels are given by
  !          \tilde{K}_kappa = - 1/kappa K_kappa
  !                          = + 1/kappa int_0^T [ phi^adj \partial_t^2 phi ] dt
  !                          = + 1/kappa int_0^T [ \partial_t^2 phi^adj phi ] dt              (equivalent)
  !
  !          \tilde{K}_rho   = - 1/rho   K_rho
  !                          = + 1/rho   int_0^T [ grad(phi^adj) * grad(phi) ] dt

  ! finalizes calculation of acoustic kernels
  do ispec = 1, NSPEC_AB
    ! acoustic simulations
    if (ispec_is_acoustic(ispec)) then
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            ! note: the contributions for acoustic kernels have been all kept positive.
            !       we add here the minus sign which is in front of the kernel definitions for K_rho and K_kappa,
            !       but together with the material factors for relative kernel definitions (dlnrho and dlnkappa)
            !       this becomes positive again with 1/rho and 1/kappa factors
            rho_invl = 1._CUSTOM_REAL / rhostore(i,j,k,ispec)
            kappa_invl = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)

            ! primary kernels
            ! relative kernel definitions (for dlnrho and dlnkappa perturbations):
            !   \tilde{K}_rho = - 1/rho K_rho
            !                 = - 1/rho { - int_0^T [ grad(phi^adj) * grad(phi) ] dt }
            !                 = 1/rho int_0^T [ grad(phi^adj) * grad(phi) ] dt
            rho_ac_kl(i,j,k,ispec) = rho_invl * rho_ac_kl(i,j,k,ispec)
            !   \tilde{K}_kappa = - 1/kappa K_kappa
            !                   = - 1/kappa { - int_0^T [ phi^adj \partial_t^2 phi ] dt }
            !                   = 1/kappa int_0^T [ phi^adj \partial_t^2 phi ] dt
            kappa_ac_kl(i,j,k,ispec) = kappa_invl * kappa_ac_kl(i,j,k,ispec)

            ! secondary, derived kernels
            ! rho prime kernel
            rhop_ac_kl(i,j,k,ispec) = rho_ac_kl(i,j,k,ispec) + kappa_ac_kl(i,j,k,ispec)
            ! vp kernel
            alpha_ac_kl(i,j,k,ispec) = 2._CUSTOM_REAL *  kappa_ac_kl(i,j,k,ispec)
          enddo
        enddo
      enddo
    endif ! acoustic
  enddo

  ! overall min/max value
  call max_all_cr(maxval(rho_ac_kl),rho_ac_max)
  call max_all_cr(maxval(rhop_ac_kl),rhop_ac_max)
  call max_all_cr(maxval(kappa_ac_kl),kappa_ac_max)
  call max_all_cr(maxval(alpha_ac_kl),alpha_ac_max)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Acoustic kernels:'
    write(IMAIN,*) '  maximum value of rho kernel       = ',rho_ac_max
    write(IMAIN,*) '  maximum value of kappa kernel     = ',kappa_ac_max
    write(IMAIN,*)
    write(IMAIN,*) '  maximum value of rho prime kernel = ',rhop_ac_max
    write(IMAIN,*) '  maximum value of alpha kernel     = ',alpha_ac_max
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (ADIOS_FOR_KERNELS) then
    call save_kernels_acoustic_adios()
  else
    ! save kernels to binary files
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rho_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file rho_acoustic_kernel.bin'
    write(IOUT) rho_ac_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'kappa_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file kappa_acoustic_kernel.bin'
    write(IOUT) kappa_ac_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhop_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file rhop_acoustic_kernel.bin'
    write(IOUT) rhop_ac_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'alpha_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file alpha_acoustic_kernel.bin'
    write(IOUT) alpha_ac_kl
    close(IOUT)

  endif

  end subroutine save_kernels_acoustic

!
!-------------------------------------------------------------------------------------------------
!

!> Save elastic isotropic kernels

  subroutine save_kernels_elastic_iso()

  use specfem_par, only: CUSTOM_REAL,NSPEC_AB,NSPEC_ADJOINT,ibool,mustore,kappastore, &
                         ADIOS_FOR_KERNELS,IOUT,prname, &
                         myrank,IMAIN
  use specfem_par_elastic

  implicit none

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhop_kl, alpha_kl, beta_kl

  ! local parameters
  integer:: ispec,i,j,k,iglob,ier
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal

  ! stats
  real(kind=CUSTOM_REAL) :: rho_max,rhop_max,kappa_max,mu_max,alpha_max,beta_max

  ! derived kernels
  ! vp,vs,density prime kernel
  allocate(alpha_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
           beta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
           rhop_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2252')
  if (ier /= 0) stop 'error allocating array alpha_kl, beta_kl, rhop_kl'
  alpha_kl(:,:,:,:) = 0.0_CUSTOM_REAL; beta_kl(:,:,:,:) = 0.0_CUSTOM_REAL
  rhop_kl(:,:,:,:) = 0.0_CUSTOM_REAL

  ! isotropic adjoint kernels (see e.g. Tromp et al. 2005)
  ! for a parameterization: (rho,mu,kappa) "primary" kernels
  ! density kernel
  ! multiplies with rho

  ! finalizes calculation of rhop, beta, alpha kernels
  do ispec = 1, NSPEC_AB

    ! elastic simulations
    if (ispec_is_elastic(ispec)) then

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)

            ! Store local material values
            rhol = rho_vs(i,j,k,ispec)*rho_vs(i,j,k,ispec) / mustore(i,j,k,ispec)
            mul = mustore(i,j,k,ispec)
            kappal = kappastore(i,j,k,ispec)

            rho_kl(i,j,k,ispec) = - rhol * rho_kl(i,j,k,ispec)

            ! shear modulus kernel
            mu_kl(i,j,k,ispec) = - 2.0 * mul * mu_kl(i,j,k,ispec)

            ! bulk modulus kernel
            kappa_kl(i,j,k,ispec) = - kappal * kappa_kl(i,j,k,ispec)

            ! for a parameterization: (rho,alpha,beta)
            ! density prime kernel
            rhop_kl(i,j,k,ispec) = rho_kl(i,j,k,ispec) + kappa_kl(i,j,k,ispec) + mu_kl(i,j,k,ispec)

            ! vs kernel
            beta_kl(i,j,k,ispec) = &
            2.0 * ( mu_kl(i,j,k,ispec) - 4.0 / 3.0 * mul / kappal * kappa_kl(i,j,k,ispec) )

            ! vp kernel
            alpha_kl(i,j,k,ispec) = &
            2.0  * ( 1.0 + 4.0 / 3.0 * mul / kappal ) * kappa_kl(i,j,k,ispec)

          enddo
        enddo
      enddo

    endif ! elastic

  enddo

  ! overall min/max value
  call max_all_cr(maxval(rho_kl),rho_max)
  call max_all_cr(maxval(rhop_kl),rhop_max)
  call max_all_cr(maxval(kappa_kl),kappa_max)
  call max_all_cr(maxval(mu_kl),mu_max)
  call max_all_cr(maxval(alpha_kl),alpha_max)
  call max_all_cr(maxval(beta_kl),beta_max)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Elastic kernels:'
    write(IMAIN,*) '  maximum value of rho  kernel      = ',rho_max
    write(IMAIN,*) '  maximum value of kappa kernel     = ',kappa_max
    write(IMAIN,*) '  maximum value of mu kernel        = ',mu_max
    write(IMAIN,*)
    write(IMAIN,*) '  maximum value of rho prime kernel = ',rhop_max
    write(IMAIN,*) '  maximum value of alpha kernel     = ',alpha_max
    write(IMAIN,*) '  maximum value of beta kernel      = ',beta_max
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (ADIOS_FOR_KERNELS) then
    call save_kernels_elastic_iso_adios(rhop_kl, alpha_kl, beta_kl)
  else
    ! save kernels to binary files
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rho_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file rho_kernel.bin'
    write(IOUT) rho_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'mu_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file mu_kernel.bin'
    write(IOUT) mu_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'kappa_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file kappa_kernel.bin'
    write(IOUT) kappa_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhop_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file rhop_kernel.bin'
    write(IOUT) rhop_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'beta_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file beta_kernel.bin'
    write(IOUT) beta_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'alpha_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file alpha_kernel.bin'
    write(IOUT) alpha_kl
    close(IOUT)
  endif

  deallocate(rhop_kl,alpha_kl,beta_kl)

  end subroutine save_kernels_elastic_iso

!
!-------------------------------------------------------------------------------------------------
!

  !> Save elastic anisotropic kernels

  subroutine save_kernels_elastic_aniso()

  use specfem_par, only: CUSTOM_REAL,NSPEC_AB,NSPEC_ADJOINT,ibool,mustore,kappastore, &
                         SAVE_TRANSVERSE_KL, &
                         ADIOS_FOR_KERNELS,IOUT,prname, &
                         myrank,IMAIN, ANISOTROPY
  use specfem_par_elastic

  implicit none

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: alphav_kl,alphah_kl,betav_kl,betah_kl, eta_kl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: c11_kl,c12_kl,c13_kl,c14_kl,c15_kl,c16_kl, &
                                                             c22_kl,c23_kl,c24_kl,c25_kl,c26_kl, &
                                                             c33_kl,c34_kl,c35_kl,c36_kl, &
                                                             c44_kl,c45_kl,c46_kl, &
                                                             c55_kl,c56_kl, &
                                                             c66_kl

  ! local parameters
  integer:: ispec,i,j,k,iglob,ier
  real(kind=CUSTOM_REAL) :: mul,kappal

  ! Transverse isotropic paramters
  real(kind=CUSTOM_REAL) :: A,N,C,L,F,eta
  real(kind=CUSTOM_REAL), dimension(21) :: cijkl_kl_local
  real(kind=CUSTOM_REAL), dimension(5) :: an_kl

  ! stats
  real(kind=CUSTOM_REAL) :: rho_max
  real(kind=CUSTOM_REAL) :: alphav_max,alphah_max,betav_max,betah_max,eta_max,cijkl_max
  real(kind=CUSTOM_REAL) :: c11_kl_max,c12_kl_max,c13_kl_max,c14_kl_max,c15_kl_max,c16_kl_max, &
                            c22_kl_max,c23_kl_max,c24_kl_max,c25_kl_max,c26_kl_max, &
                            c33_kl_max,c34_kl_max,c35_kl_max,c36_kl_max, &
                            c44_kl_max,c45_kl_max,c46_kl_max, &
                            c55_kl_max,c56_kl_max, &
                            c66_kl_max

  ! allocates temporary transversely isotropic kernels
  if (SAVE_TRANSVERSE_KL) then
    allocate(alphav_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             alphah_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             betav_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             betah_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             eta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2247')
    if (ier /= 0) stop 'error allocating arrays alphav_kl,...'
    alphav_kl(:,:,:,:) = 0.0_CUSTOM_REAL; alphah_kl(:,:,:,:) = 0.0_CUSTOM_REAL
    betav_kl(:,:,:,:) = 0.0_CUSTOM_REAL; betah_kl(:,:,:,:) = 0.0_CUSTOM_REAL
    eta_kl(:,:,:,:) = 0.0_CUSTOM_REAL
  else
    allocate(c11_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c12_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c13_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c14_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c15_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c16_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c22_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c23_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c24_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c25_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c26_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c33_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c34_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c35_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c36_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c44_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c45_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c46_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c55_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c56_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             c66_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array XXXX')
      if (ier /= 0) stop 'error allocating arrays c11_kl,c22_kl,...'
      c11_kl(:,:,:,:) = 0.0_CUSTOM_REAL; c12_kl(:,:,:,:) = 0.0_CUSTOM_REAL
      c13_kl(:,:,:,:) = 0.0_CUSTOM_REAL; c14_kl(:,:,:,:) = 0.0_CUSTOM_REAL
      c15_kl(:,:,:,:) = 0.0_CUSTOM_REAL; c16_kl(:,:,:,:) = 0.0_CUSTOM_REAL
      c22_kl(:,:,:,:) = 0.0_CUSTOM_REAL; c23_kl(:,:,:,:) = 0.0_CUSTOM_REAL
      c24_kl(:,:,:,:) = 0.0_CUSTOM_REAL; c25_kl(:,:,:,:) = 0.0_CUSTOM_REAL
      c26_kl(:,:,:,:) = 0.0_CUSTOM_REAL; c33_kl(:,:,:,:) = 0.0_CUSTOM_REAL
      c34_kl(:,:,:,:) = 0.0_CUSTOM_REAL; c35_kl(:,:,:,:) = 0.0_CUSTOM_REAL
      c36_kl(:,:,:,:) = 0.0_CUSTOM_REAL; c44_kl(:,:,:,:) = 0.0_CUSTOM_REAL
      c45_kl(:,:,:,:) = 0.0_CUSTOM_REAL; c46_kl(:,:,:,:) = 0.0_CUSTOM_REAL
      c55_kl(:,:,:,:) = 0.0_CUSTOM_REAL; c56_kl(:,:,:,:) = 0.0_CUSTOM_REAL
      c66_kl(:,:,:,:) = 0.0_CUSTOM_REAL
  endif

  ! finalizes calculation of rhop, beta, alpha kernels
  do ispec = 1, NSPEC_AB

    ! elastic simulations
    if (ispec_is_elastic(ispec)) then

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)

            cijkl_kl_local(:) = - cijkl_kl(:,i,j,k,ispec)

            if (SAVE_TRANSVERSE_KL) then
              ! SAVE_TRANSVERSE_KL

              if (ANISOTROPY) then
                ! Computes parameters for an anisotropic model
                A = c11store(i,j,k,ispec)
                C = c33store(i,j,k,ispec)
                L = c44store(i,j,k,ispec)
                N = c66store(i,j,k,ispec)
                F = c13store(i,j,k,ispec)
                eta = F / (A - 2.0 * L)
              else
                ! Store local material values
                mul = mustore(i,j,k,ispec)
                kappal = kappastore(i,j,k,ispec)

                ! Computes parameters for an isotropic model
                A = kappal + 4.0 / 3.0 * mul
                C = A
                L = mul
                N = mul
                F = kappal - 2.0 / 3.0 * mul
                eta = 1.0
              endif

              ! note: cijkl_kl_local() is fully anisotropic C_ij kernel
              ! components (non-dimensionalized)
              !          for GLL point at (i,j,k,ispec)

              ! Purpose : compute the kernels for the An coeffs (an_kl)
              ! from the kernels for Cij (cijkl_kl_local)

              ! Definition of the input array cij_kl :
              ! cij_kl(1) = C11 ; cij_kl(2) = C12 ; cij_kl(3) = C13
              ! cij_kl(4) = C14 ; cij_kl(5) = C15 ; cij_kl(6) = C16
              ! cij_kl(7) = C22 ; cij_kl(8) = C23 ; cij_kl(9) = C24
              ! cij_kl(10) = C25 ; cij_kl(11) = C26 ; cij_kl(12) = C33
              ! cij_kl(13) = C34 ; cij_kl(14) = C35 ; cij_kl(15) = C36
              ! cij_kl(16) = C44 ; cij_kl(17) = C45 ; cij_kl(18) = C46
              ! cij_kl(19) = C55 ; cij_kl(20) = C56 ; cij_kl(21) = C66
              ! where the Cij (Voigt's notation) are defined as function of
              ! the components of the elastic tensor in spherical coordinates
              ! by eq. (A.1) of Chen & Tromp, GJI 168 (2007)

              ! From the relations giving Cij in function of An
              ! Checked with Min Chen's results (routine build_cij)

              an_kl(1) = cijkl_kl_local(1) + cijkl_kl_local(2) + cijkl_kl_local(7)    !A
              an_kl(2) = cijkl_kl_local(12)                                           !C
              an_kl(3) = - 2.0 * cijkl_kl_local(2) + cijkl_kl_local(21)               !N
              an_kl(4) = cijkl_kl_local(16) + cijkl_kl_local(19)                      !L
              an_kl(5) = cijkl_kl_local(3) + cijkl_kl_local(8)                        !F

              ! for parameterization: ( alpha_v, alpha_h, beta_v, beta_h, eta, rho )

              ! K_alpha_v
              alphav_kl(i,j,k,ispec) = 2.0 * C * an_kl(2)
              ! K_alpha_h
              alphah_kl(i,j,k,ispec) = 2.0 * A * an_kl(1) + 2.0 * A * eta * an_kl(5)
              ! K_beta_v
              betav_kl(i,j,k,ispec) = 2.0 * L * an_kl(4) - 4.0 * L * eta * an_kl(5)
              ! K_beta_h
              betah_kl(i,j,k,ispec) = 2.0 * N * an_kl(3)
              ! K_eta
              eta_kl(i,j,k,ispec) = F * an_kl(5)
            else
              !Cij_kernels
              c11_kl(i,j,k,ispec) = cijkl_kl_local(1)
              c12_kl(i,j,k,ispec) = cijkl_kl_local(2)
              c13_kl(i,j,k,ispec) = cijkl_kl_local(3)
              c14_kl(i,j,k,ispec) = cijkl_kl_local(4)
              c15_kl(i,j,k,ispec) = cijkl_kl_local(5)
              c16_kl(i,j,k,ispec) = cijkl_kl_local(6)
              c22_kl(i,j,k,ispec) = cijkl_kl_local(7)
              c23_kl(i,j,k,ispec) = cijkl_kl_local(8)
              c24_kl(i,j,k,ispec) = cijkl_kl_local(9)
              c25_kl(i,j,k,ispec) = cijkl_kl_local(10)
              c26_kl(i,j,k,ispec) = cijkl_kl_local(11)
              c33_kl(i,j,k,ispec) = cijkl_kl_local(12)
              c34_kl(i,j,k,ispec) = cijkl_kl_local(13)
              c35_kl(i,j,k,ispec) = cijkl_kl_local(14)
              c36_kl(i,j,k,ispec) = cijkl_kl_local(15)
              c44_kl(i,j,k,ispec) = cijkl_kl_local(16)
              c45_kl(i,j,k,ispec) = cijkl_kl_local(17)
              c46_kl(i,j,k,ispec) = cijkl_kl_local(18)
              c55_kl(i,j,k,ispec) = cijkl_kl_local(19)
              c56_kl(i,j,k,ispec) = cijkl_kl_local(20)
              c66_kl(i,j,k,ispec) = cijkl_kl_local(21)
            endif
          enddo
        enddo
      enddo

    endif ! elastic

  enddo

  ! overall min/max value

  if (SAVE_TRANSVERSE_KL) then
    call max_all_cr(maxval(alphav_kl),alphav_max)
    call max_all_cr(maxval(alphah_kl),alphah_max)
    call max_all_cr(maxval(betav_kl),betav_max)
    call max_all_cr(maxval(betah_kl),betah_max)
    call max_all_cr(maxval(eta_kl),eta_max)
  else
    call max_all_cr(maxval(-rho_kl),rho_max)
    call max_all_cr(maxval(-cijkl_kl),cijkl_max)
    call max_all_cr(maxval(c11_kl),c11_kl_max)
    call max_all_cr(maxval(c12_kl),c12_kl_max)
    call max_all_cr(maxval(c13_kl),c13_kl_max)
    call max_all_cr(maxval(c14_kl),c14_kl_max)
    call max_all_cr(maxval(c15_kl),c15_kl_max)
    call max_all_cr(maxval(c16_kl),c16_kl_max)
    call max_all_cr(maxval(c22_kl),c22_kl_max)
    call max_all_cr(maxval(c23_kl),c23_kl_max)
    call max_all_cr(maxval(c24_kl),c24_kl_max)
    call max_all_cr(maxval(c25_kl),c25_kl_max)
    call max_all_cr(maxval(c26_kl),c26_kl_max)
    call max_all_cr(maxval(c33_kl),c33_kl_max)
    call max_all_cr(maxval(c34_kl),c34_kl_max)
    call max_all_cr(maxval(c35_kl),c35_kl_max)
    call max_all_cr(maxval(c36_kl),c36_kl_max)
    call max_all_cr(maxval(c44_kl),c44_kl_max)
    call max_all_cr(maxval(c45_kl),c45_kl_max)
    call max_all_cr(maxval(c46_kl),c46_kl_max)
    call max_all_cr(maxval(c55_kl),c55_kl_max)
    call max_all_cr(maxval(c56_kl),c56_kl_max)
    call max_all_cr(maxval(c66_kl),c66_kl_max)
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Elastic kernels:'

    if (SAVE_TRANSVERSE_KL) then
      ! tranverse isotropic
      write(IMAIN,*) '  maximum value of alphav kernel     = ',alphav_max
      write(IMAIN,*) '  maximum value of alphah kernel     = ',alphah_max
      write(IMAIN,*) '  maximum value of betav kernel      = ',betav_max
      write(IMAIN,*) '  maximum value of betah kernel      = ',betah_max
      write(IMAIN,*) '  maximum value of eta kernel        = ',eta_max
    else
      ! fully anisotropic
      write(IMAIN,*) '  maximum value of rho  kernel      = ',rho_max
      write(IMAIN,*) '  maximum value of cijkl kernel     = ',cijkl_max
      write(IMAIN,*) '  maximum value of c11 kernel     = ',c11_kl_max
      write(IMAIN,*) '  maximum value of c12 kernel     = ',c12_kl_max
      write(IMAIN,*) '  maximum value of c13 kernel     = ',c13_kl_max
      write(IMAIN,*) '  maximum value of c14 kernel     = ',c14_kl_max
      write(IMAIN,*) '  maximum value of c15 kernel     = ',c15_kl_max
      write(IMAIN,*) '  maximum value of c16 kernel     = ',c16_kl_max
      write(IMAIN,*) '  maximum value of c22 kernel     = ',c22_kl_max
      write(IMAIN,*) '  maximum value of c23 kernel     = ',c23_kl_max
      write(IMAIN,*) '  maximum value of c24 kernel     = ',c24_kl_max
      write(IMAIN,*) '  maximum value of c25 kernel     = ',c25_kl_max
      write(IMAIN,*) '  maximum value of c26 kernel     = ',c26_kl_max
      write(IMAIN,*) '  maximum value of c33 kernel     = ',c33_kl_max
      write(IMAIN,*) '  maximum value of c34 kernel     = ',c34_kl_max
      write(IMAIN,*) '  maximum value of c35 kernel     = ',c35_kl_max
      write(IMAIN,*) '  maximum value of c36 kernel     = ',c36_kl_max
      write(IMAIN,*) '  maximum value of c44 kernel     = ',c44_kl_max
      write(IMAIN,*) '  maximum value of c45 kernel     = ',c45_kl_max
      write(IMAIN,*) '  maximum value of c46 kernel     = ',c46_kl_max
      write(IMAIN,*) '  maximum value of c55 kernel     = ',c55_kl_max
      write(IMAIN,*) '  maximum value of c56 kernel     = ',c56_kl_max
      write(IMAIN,*) '  maximum value of c66 kernel     = ',c66_kl_max
    endif

    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (ADIOS_FOR_KERNELS) then
    call save_kernels_elastic_aniso_adios(alphav_kl, alphah_kl, betav_kl, betah_kl, eta_kl, &
                                          c11_kl,c12_kl,c13_kl,c14_kl,c15_kl,c16_kl, &
                                          c22_kl,c23_kl,c24_kl,c25_kl,c26_kl, &
                                          c33_kl,c34_kl,c35_kl,c36_kl, &
                                          c44_kl,c45_kl,c46_kl, &
                                          c55_kl,c56_kl, &
                                          c66_kl)
  else
    ! outputs transverse isotropic kernels only
    if (SAVE_TRANSVERSE_KL) then
      ! transverse isotropic kernels
      ! (alpha_v, alpha_h, beta_v, beta_h, eta, rho ) parameterization
      open(unit=IOUT,file=trim(prname)//'alphav_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) alphav_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'alphah_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) alphah_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'betav_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) betav_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'betah_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) betah_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'eta_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) eta_kl
      close(IOUT)
    else
      ! fully anisotropic kernels
      ! note: the C_ij and density kernels are not for relative perturbations
      ! (delta ln( m_i) = delta m_i / m_i),
      !          but absolute perturbations (delta m_i = m_i - m_0).
      open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) - rho_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c11_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c11_kernel.bin'
      write(IOUT) c11_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c12_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c12_kernel.bin'
      write(IOUT) c12_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c13_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c13_kernel.bin'
      write(IOUT) c13_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c14_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c14_kernel.bin'
      write(IOUT) c14_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c15_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c15_kernel.bin'
      write(IOUT) c15_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c16_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c16_kernel.bin'
      write(IOUT) c16_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c22_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c22_kernel.bin'
      write(IOUT) c22_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c23_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c23_kernel.bin'
      write(IOUT) c23_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c24_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c24_kernel.bin'
      write(IOUT) c24_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c25_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c25_kernel.bin'
      write(IOUT) c25_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c26_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c26_kernel.bin'
      write(IOUT) c26_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c33_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c33_kernel.bin'
      write(IOUT) c33_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c34_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c34_kernel.bin'
      write(IOUT) c34_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c35_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c35_kernel.bin'
      write(IOUT) c35_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c36_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c36_kernel.bin'
      write(IOUT) c36_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c44_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c44_kernel.bin'
      write(IOUT) c44_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c45_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c45_kernel.bin'
      write(IOUT) c45_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c46_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c46_kernel.bin'
      write(IOUT) c46_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c55_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c55_kernel.bin'
      write(IOUT) c55_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c56_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c56_kernel.bin'
      write(IOUT) c56_kl
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'c66_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file c66_kernel.bin'
      write(IOUT) c66_kl
      close(IOUT)
    endif

  endif

  if (SAVE_TRANSVERSE_KL) then
    deallocate(alphav_kl,alphah_kl,betav_kl,betah_kl,eta_kl)
  else
    deallocate(c11_kl,c12_kl,c13_kl,c14_kl,c15_kl,c16_kl, &
               c22_kl,c23_kl,c24_kl,c25_kl,c26_kl, &
               c33_kl,c34_kl,c35_kl,c36_kl, &
               c44_kl,c45_kl,c46_kl, &
               c55_kl,c56_kl, &
               c66_kl)
  endif

  end subroutine save_kernels_elastic_aniso

!
!-------------------------------------------------------------------------------------------------
!

  !> Save poroelastic related kernels

  subroutine save_kernels_poroelastic()

  use specfem_par
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer:: ispec,i,j,k,ier
  real(kind=CUSTOM_REAL) :: rhol_s,rhol_f,rhol_bar,phil,tortl
  real(kind=CUSTOM_REAL) :: kappal_s ! mul_s
  real(kind=CUSTOM_REAL) :: kappal_f,etal_f
  real(kind=CUSTOM_REAL) :: mul_fr,kappal_fr
  real(kind=CUSTOM_REAL) :: permlxx,permlxy,permlxz,permlyz,permlyy,permlzz
  real(kind=CUSTOM_REAL) :: D_biot,H_biot,C_biot,M_biot,B_biot
  real(kind=CUSTOM_REAL) :: cpIsquare,cpIIsquare,cssquare
  real(kind=CUSTOM_REAL) :: rholb,ratio,dd1,gamma1,gamma2,gamma3,gamma4
  real(kind=CUSTOM_REAL) :: afactor,bfactor,cfactor

  ! finalizes calculation of rhop, beta, alpha kernels
  do ispec = 1, NSPEC_AB

    ! poroelastic simulations
    if (ispec_is_poroelastic(ispec)) then

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

            ! isotropic adjoint kernels (see e.g. Morency et al. 2009)

            ! get poroelastic parameters of current local GLL
            phil = phistore(i,j,k,ispec)
            tortl = tortstore(i,j,k,ispec)
            rhol_s = rhoarraystore(1,i,j,k,ispec)
            rhol_f = rhoarraystore(2,i,j,k,ispec)
            rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f
            kappal_s = kappaarraystore(1,i,j,k,ispec)
            kappal_f = kappaarraystore(2,i,j,k,ispec)
            kappal_fr = kappaarraystore(3,i,j,k,ispec)
            mul_fr = mustore(i,j,k,ispec)
            etal_f = etastore(i,j,k,ispec)
            permlxx = permstore(1,i,j,k,ispec)
            permlxy = permstore(2,i,j,k,ispec)
            permlxz = permstore(3,i,j,k,ispec)
            permlyy = permstore(4,i,j,k,ispec)
            permlyz = permstore(5,i,j,k,ispec)
            permlzz = permstore(6,i,j,k,ispec)

            ! Biot coef
            D_biot = kappal_s*(1._CUSTOM_REAL + phil*(kappal_s/kappal_f - 1._CUSTOM_REAL))
            H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + &
                      kappal_fr + 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
            B_biot = H_biot - 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
            C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
            M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)

            ! Approximated velocities (no viscous dissipation)
            afactor = rhol_bar - phil/tortl*rhol_f
            bfactor = H_biot + phil*rhol_bar/(tortl*rhol_f)*M_biot - 2._CUSTOM_REAL*phil/tortl*C_biot
            cfactor = phil/(tortl*rhol_f)*(H_biot*M_biot - C_biot*C_biot)
            cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4._CUSTOM_REAL*afactor*cfactor))/(2._CUSTOM_REAL*afactor)
            cpIIsquare = (bfactor - sqrt(bfactor*bfactor - 4._CUSTOM_REAL*afactor*cfactor))/(2._CUSTOM_REAL*afactor)
            cssquare = mul_fr/afactor

            ! extras needed
            ! Approximated ratio r = amplitude "w" field/amplitude "s" field (no viscous
            ! dissipation)
            gamma1 = H_biot - phil/tortl*C_biot
            gamma2 = C_biot - phil/tortl*M_biot
            gamma3 = phil/tortl*( M_biot*(afactor/rhol_f + phil/tortl) - C_biot)
            gamma4 = phil/tortl*( C_biot*(afactor/rhol_f + phil/tortl) - H_biot)
            ratio = 0.5_CUSTOM_REAL*(gamma1 - gamma3)/gamma4 + &
                    0.5_CUSTOM_REAL*sqrt((gamma1-gamma3)**2/gamma4**2 + 4._CUSTOM_REAL * gamma2/gamma4)
            rholb = rhol_bar - phil*rhol_f/tortl
            dd1 = (1._CUSTOM_REAL+rholb/rhol_f)*ratio**2 + 2._CUSTOM_REAL*ratio + tortl/phil

            ! primary kernels
            rhot_kl(i,j,k,ispec) = - rhol_bar * rhot_kl(i,j,k,ispec)
            rhof_kl(i,j,k,ispec) = - rhol_f * rhof_kl(i,j,k,ispec)
            sm_kl(i,j,k,ispec) = - rhol_f*tortl/phil * sm_kl(i,j,k,ispec)
            !at the moment suitable for constant permeability
            eta_kl(i,j,k,ispec) = - etal_f/permlxx * eta_kl(i,j,k,ispec)
            mufr_kl(i,j,k,ispec) = - 2._CUSTOM_REAL * mul_fr * mufr_kl(i,j,k,ispec)
            B_kl(i,j,k,ispec) = - B_biot * B_kl(i,j,k,ispec)
            C_kl(i,j,k,ispec) = - C_biot * C_kl(i,j,k,ispec)
            M_kl(i,j,k,ispec) = - M_biot * M_kl(i,j,k,ispec)

            ! density kernels
            rhob_kl(i,j,k,ispec) = rhot_kl(i,j,k,ispec) + B_kl(i,j,k,ispec) + mufr_kl(i,j,k,ispec)
            rhofb_kl(i,j,k,ispec) = rhof_kl(i,j,k,ispec) + C_kl(i,j,k,ispec) + M_kl(i,j,k,ispec) + sm_kl(i,j,k,ispec)
            Bb_kl(i,j,k,ispec) = B_kl(i,j,k,ispec)
            Cb_kl(i,j,k,ispec) = C_kl(i,j,k,ispec)
            Mb_kl(i,j,k,ispec) = M_kl(i,j,k,ispec)
            mufrb_kl(i,j,k,ispec) = mufr_kl(i,j,k,ispec)
            phi_kl(i,j,k,ispec) = - sm_kl(i,j,k,ispec) - M_kl(i,j,k,ispec)

            ! wavespeed kernels
            rhobb_kl(i,j,k,ispec) = rhob_kl(i,j,k,ispec) - &
                      phil*rhol_f/(tortl*B_biot) * &
                      (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil / &
                      tortl*ratio +1._CUSTOM_REAL)/dd1 + &
                      (rhol_bar**2*ratio**2/rhol_f**2*(phil / &
                      tortl*ratio+1)*(phil/tortl*ratio + &
                      phil/tortl * &
                      (1+rhol_f/rhol_bar)-1))/dd1**2) - &
                      4._CUSTOM_REAL/3._CUSTOM_REAL*cssquare )*Bb_kl(i,j,k,ispec) - &
                      rhol_bar*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                      (phil/tortl*ratio + &
                      1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,k,ispec) + &
                      rhol_bar*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                      (phil/tortl*ratio+1._CUSTOM_REAL)/dd1 - &
                      phil*ratio/tortl*(phil / &
                      tortl*ratio+1._CUSTOM_REAL)*&
                      (1+rhol_bar*ratio/rhol_f)/dd1**2)*Cb_kl(i,j,k,ispec)+ &
                      phil*rhol_f*cssquare / &
                      (tortl*mul_fr)*mufrb_kl(i,j,k,ispec)
            rhofbb_kl(i,j,k,ispec) = rhofb_kl(i,j,k,ispec) + &
                       phil*rhol_f/(tortl*B_biot) * &
                       (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil/ &
                       tortl*ratio +1._CUSTOM_REAL)/dd1+&
                       (rhol_bar**2*ratio**2/rhol_f**2*(phil/ &
                       tortl*ratio+1)*(phil/tortl*ratio+ &
                       phil/tortl*&
                       (1+rhol_f/rhol_bar)-1))/dd1**2)- &
                       4._CUSTOM_REAL/3._CUSTOM_REAL*cssquare )*Bb_kl(i,j,k,ispec) + &
                       rhol_bar*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                       (phil/tortl*ratio + &
                       1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,k,ispec) - &
                       rhol_bar*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                       (phil/tortl*ratio+1._CUSTOM_REAL)/dd1 - &
                       phil*ratio/tortl*(phil/ &
                       tortl*ratio+1._CUSTOM_REAL)*&
                       (1+rhol_bar*ratio/rhol_f)/dd1**2)*Cb_kl(i,j,k,ispec)- &
                       phil*rhol_f*cssquare/ &
                       (tortl*mul_fr)*mufrb_kl(i,j,k,ispec)
            phib_kl(i,j,k,ispec) = phi_kl(i,j,k,ispec) - &
                       phil*rhol_bar/(tortl*B_biot) &
                       * ( cpIsquare - rhol_f/rhol_bar*cpIIsquare- &
                       (cpIsquare-cpIIsquare)*( (2._CUSTOM_REAL*ratio**2*phil/ &
                       tortl + (1._CUSTOM_REAL+&
                       rhol_f/rhol_bar)* &
                       (2._CUSTOM_REAL*ratio*phil/tortl+&
                       1._CUSTOM_REAL))/dd1 + (phil/tortl*ratio+ &
                       1._CUSTOM_REAL)*(phil*&
                       ratio/tortl+phil/tortl* &
                       (1._CUSTOM_REAL+rhol_f/&
                       rhol_bar)-1._CUSTOM_REAL)*((1._CUSTOM_REAL+ &
                       rhol_bar/rhol_f-&
                       2._CUSTOM_REAL*phil/tortl)*ratio**2+2._CUSTOM_REAL*ratio)/dd1**2) - &
                       4._CUSTOM_REAL/3._CUSTOM_REAL*rhol_f*cssquare/rhol_bar)*Bb_kl(i,j,k,ispec) + &
                       rhol_f/M_biot * (cpIsquare-cpIIsquare)*(&
                       2._CUSTOM_REAL*ratio*(phil/tortl*ratio+1._CUSTOM_REAL)/dd1 - &
                       (phil/tortl*ratio+1._CUSTOM_REAL)**2*( &
                       (1._CUSTOM_REAL+rhol_bar/&
                       rhol_f-2._CUSTOM_REAL*phil/tortl)*ratio**2+2._CUSTOM_REAL*ratio)/dd1**2 &
                       )*Mb_kl(i,j,k,ispec) + &
                       phil*rhol_f/(tortl*C_biot)* &
                       (cpIsquare-cpIIsquare)*ratio* (&
                       (1._CUSTOM_REAL+rhol_f/rhol_bar*ratio)/dd1 - &
                       (phil/tortl*ratio+1._CUSTOM_REAL)* &
                       (1._CUSTOM_REAL+rhol_bar/&
                       rhol_f*ratio)*((1._CUSTOM_REAL+rhol_bar/rhol_f-2._CUSTOM_REAL*&
                       phil/tortl)*ratio+2._CUSTOM_REAL)/dd1**2&
                        )*Cb_kl(i,j,k,ispec) -&
                       phil*rhol_f*cssquare &
                       /(tortl*mul_fr)*mufrb_kl(i,j,k,ispec)
            cpI_kl(i,j,k,ispec) = 2._CUSTOM_REAL*cpIsquare/B_biot*rhol_bar*( &
                       1._CUSTOM_REAL-phil/tortl + &
                       (phil/tortl*ratio+ &
                       1._CUSTOM_REAL)*(phil/tortl*&
                       ratio+phil/tortl* &
                       (1._CUSTOM_REAL+rhol_f/rhol_bar)-&
                       1._CUSTOM_REAL)/dd1 &
                        )* Bb_kl(i,j,k,ispec) +&
                       2._CUSTOM_REAL*cpIsquare*rhol_f*tortl/(phil*M_biot) *&
                       (phil/tortl*ratio+1._CUSTOM_REAL)**2/dd1*Mb_kl(i,j,k,ispec)+&
                       2._CUSTOM_REAL*cpIsquare*rhol_f/C_biot * &
                       (phil/tortl*ratio+1._CUSTOM_REAL)* &
                       (1._CUSTOM_REAL+rhol_bar/&
                       rhol_f*ratio)/dd1*Cb_kl(i,j,k,ispec)
            cpII_kl(i,j,k,ispec) = 2._CUSTOM_REAL*cpIIsquare*rhol_bar/B_biot * (&
                       phil*rhol_f/(tortl*rhol_bar) - &
                       (phil/tortl*ratio+ &
                       1._CUSTOM_REAL)*(phil/tortl*&
                       ratio+phil/tortl* &
                       (1._CUSTOM_REAL+rhol_f/rhol_bar)-&
                       1._CUSTOM_REAL)/dd1  ) * Bb_kl(i,j,k,ispec) +&
                       2._CUSTOM_REAL*cpIIsquare*rhol_f*tortl/(phil*M_biot) * (&
                       1._CUSTOM_REAL - (phil/tortl*ratio+ &
                       1._CUSTOM_REAL)**2/dd1  )*Mb_kl(i,j,k,ispec) + &
                       2._CUSTOM_REAL*cpIIsquare*rhol_f/C_biot * (&
                       1._CUSTOM_REAL - (phil/tortl*ratio+ &
                       1._CUSTOM_REAL)*(1._CUSTOM_REAL+&
                       rhol_bar/rhol_f*ratio)/dd1)*Cb_kl(i,j,k,ispec)
            cs_kl(i,j,k,ispec) = - 8._CUSTOM_REAL/3._CUSTOM_REAL*cssquare* &
                       rhol_bar/B_biot*(1._CUSTOM_REAL-&
                       phil*rhol_f/(tortl* &
                       rhol_bar))*Bb_kl(i,j,k,ispec) + &
                       2._CUSTOM_REAL*(rhol_bar-rhol_f*&
                       phil/tortl)/&
                       mul_fr*cssquare*mufrb_kl(i,j,k,ispec)
            ratio_kl(i,j,k,ispec) = ratio*rhol_bar*phil/(tortl*B_biot) * &
                       (cpIsquare-cpIIsquare) * ( &
                       phil/tortl*(2._CUSTOM_REAL*ratio+1._CUSTOM_REAL+rhol_f/ &
                       rhol_bar)/dd1 - (phil/tortl*ratio+1._CUSTOM_REAL)*&
                       (phil/tortl*ratio+phil/tortl*(&
                       1._CUSTOM_REAL+rhol_f/rhol_bar)-1._CUSTOM_REAL)*(2._CUSTOM_REAL*ratio*(&
                       1._CUSTOM_REAL+rhol_bar/rhol_f-phil/tortl) +&
                       2._CUSTOM_REAL)/dd1**2  )*Bb_kl(i,j,k,ispec) + &
                       ratio*rhol_f*tortl/(phil*M_biot)*(cpIsquare-cpIIsquare) * &
                       2._CUSTOM_REAL*phil/tortl * (&
                       (phil/tortl*ratio+1._CUSTOM_REAL)/dd1 - &
                       (phil/tortl*ratio+1._CUSTOM_REAL)**2*( &
                       (1._CUSTOM_REAL+rhol_bar/&
                       rhol_f-phil/tortl)*ratio+ &
                       1._CUSTOM_REAL)/dd1**2)*Mb_kl(i,j,k,ispec) +&
                       ratio*rhol_f/C_biot*(cpIsquare-cpIIsquare) * (&
                       (2._CUSTOM_REAL*phil*rhol_bar* &
                       ratio/(tortl*rhol_f)+&
                       phil/tortl+rhol_bar/rhol_f)/dd1 - &
                       2._CUSTOM_REAL*phil/tortl*(phil/tortl*ratio+&
                       1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar/rhol_f*ratio)*((1._CUSTOM_REAL+&
                       rhol_bar/rhol_f- &
                       phil/tortl)*ratio+1._CUSTOM_REAL)/&
                       dd1**2)*Cb_kl(i,j,k,ispec)
          enddo
        enddo
      enddo

    endif ! poroelastic

  enddo

  ! save kernels to binary files
  if (ADIOS_FOR_KERNELS) then
    call save_kernels_poroelastic_adios()
  else
    ! primary kernels
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhot_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file rhot_primeporo_kernel.bin'
    write(IOUT) rhot_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhof_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file rhof_primeporo_kernel.bin'
    write(IOUT) rhof_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'sm_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file sm_primeporo_kernel.bin'
    write(IOUT) sm_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'eta_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file eta_primeporo_kernel.bin'
    write(IOUT) eta_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'mufr_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file mufr_primeporo_kernel.bin'
    write(IOUT) mufr_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'B_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file B_primeporo_kernel.bin'
    write(IOUT) B_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'C_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file C_primeporo_kernel.bin'
    write(IOUT) C_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'M_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file M_primeporo_kernel.bin'
    write(IOUT) M_kl
    close(IOUT)

    ! density kernels
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhob_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file rhob_densityporo_kernel.bin'
    write(IOUT) rhob_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhofb_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file rhofb_densityporo_kernel.bin'
    write(IOUT) rhofb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'phi_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file phi_densityporo_kernel.bin'
    write(IOUT) phi_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'mufrb_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file mufrb_densityporo_kernel.bin'
    write(IOUT) mufrb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'Bb_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file Bb_densityporo_kernel.bin'
    write(IOUT) Bb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'Cb_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file Cb_densityporo_kernel.bin'
    write(IOUT) Cb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'Mb_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file Mb_densityporo_kernel.bin'
    write(IOUT) Mb_kl
    close(IOUT)

    ! wavespeed kernels
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhobb_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file rhobb_waveporo_kernel.bin'
    write(IOUT) rhobb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhofbb_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file rhofbb_waveporo_kernel.bin'
    write(IOUT) rhofbb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'phib_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file phib_waveporo_kernel.bin'
    write(IOUT) phib_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'cs_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file cs_waveporo_kernel.bin'
    write(IOUT) cs_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'cpI_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file cpI_waveporo_kernel.bin'
    write(IOUT) cpI_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'cpII_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file cpII_waveporo_kernel.bin'
    write(IOUT) cpII_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'ratio_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file ratio_waveporo_kernel.bin'
    write(IOUT) ratio_kl
    close(IOUT)

  endif

  end subroutine save_kernels_poroelastic

!
!-------------------------------------------------------------------------------------------------
!

!> Save moho boundary kernels

  subroutine save_kernels_moho()

  use specfem_par, only: CUSTOM_REAL,ADIOS_FOR_KERNELS, &
                         IOUT,prname,myrank,IMAIN,SAVE_MOHO_MESH
  use specfem_par_elastic

  implicit none

  ! local parameters
  integer :: ier
  ! stats
  real(kind=CUSTOM_REAL) :: moho_max

  ! safety check
  if (.not. SAVE_MOHO_MESH) return

  ! overall min/max value
  call max_all_cr(maxval(moho_kl),moho_max)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Moho boundary kernels:'
    write(IMAIN,*) '  maximum value of moho kernel      = ',moho_max
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (ADIOS_FOR_KERNELS) then
    call save_kernels_moho_adios()
  else
    ! binary file output
    open(unit=IOUT,file=prname(1:len_trim(prname))//'moho_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening file moho_kernel.bin'
    write(IOUT) moho_kl
    close(IOUT)
  endif

  end subroutine save_kernels_moho

!
!-------------------------------------------------------------------------------------------------
!

!> Save Hessians

  subroutine save_kernels_Hessian()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  integer :: ier

  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    ! scales approximate Hessian
    hess_ac_kl(:,:,:,:) = 2._CUSTOM_REAL * hess_ac_kl(:,:,:,:)
  endif

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    ! scales approximate Hessian
    hess_kl(:,:,:,:) = 2._CUSTOM_REAL * hess_kl(:,:,:,:)
  endif

  if (ADIOS_FOR_KERNELS) then
    call save_kernels_Hessian_adios()
  else
    ! acoustic domains
    if (ACOUSTIC_SIMULATION) then
      ! stores into file
      open(unit=IOUT,file=trim(prname)//'hess_acoustic_kernel.bin', &
           status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0) stop 'error opening file hess_acoustic_kernel.bin'
      write(IOUT) hess_ac_kl
      close(IOUT)
    endif

    ! elastic domains
    if (ELASTIC_SIMULATION) then
      ! stores into file
      open(unit=IOUT,file=trim(prname)//'hess_kernel.bin', &
           status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0) stop 'error opening file hess_kernel.bin'
      write(IOUT) hess_kl
      close(IOUT)
    endif
  endif

  end subroutine save_kernels_Hessian

!
!-------------------------------------------------------------------------------------------------
!


  subroutine save_kernels_source_derivatives()

  use specfem_par

  implicit none

  ! local parameters
  integer :: irec_local,ier
  character(len=MAX_STRING_LEN) :: outputname

  ! checks
  if (ADIOS_FOR_KERNELS ) stop 'Source derivative kernels not implemented yet for ADIOS'

  ! writes out derivative kernels
  do irec_local = 1, nrec_local
    write(outputname,'(a,i6.6)') OUTPUT_FILES(1:len_trim(OUTPUT_FILES)) // &
        '/src_frechet.',number_receiver_global(irec_local)

    open(unit=IOUT,file=trim(outputname),status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,'error opening file: ',trim(outputname)
      call exit_mpi(myrank,'error opening file src_frechet.**')
    endif

    !
    ! r -> z, theta -> -y, phi -> x
    !
    !  Mrr =  Mzz
    !  Mtt =  Myy
    !  Mpp =  Mxx
    !  Mrt = -Myz
    !  Mrp =  Mxz
    !  Mtp = -Mxy
    write(IOUT,*) Mzz_der(irec_local)
    write(IOUT,*) Myy_der(irec_local)
    write(IOUT,*) Mxx_der(irec_local)
    write(IOUT,*) -Myz_der(irec_local)
    write(IOUT,*) Mxz_der(irec_local)
    write(IOUT,*) -Mxy_der(irec_local)
    write(IOUT,*) sloc_der(1,irec_local)
    write(IOUT,*) sloc_der(2,irec_local)
    write(IOUT,*) sloc_der(3,irec_local)

    close(IOUT)
  enddo

  end subroutine save_kernels_source_derivatives

