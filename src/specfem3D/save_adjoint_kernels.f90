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
!
! United States and French Government Sponsorship Acknowledged.


!==============================================================================
! \file save_adjoint_kernels
!
! TODO
! * Better doxygen documentation.
!==============================================================================


!==============================================================================
!> Save kernels.
subroutine save_adjoint_kernels()

  use constants, only: CUSTOM_REAL, SAVE_TRANSVERSE_KL, ANISOTROPIC_KL, &
                       APPROXIMATE_HESS_KL, NGLLX, NGLLY, NGLLZ
  use specfem_par, only: LOCAL_PATH, myrank, sigma_kl, NSPEC_AB, &
                         ADIOS_FOR_KERNELS, NOISE_TOMOGRAPHY, NSPEC_ADJOINT
  use specfem_par_acoustic, only: ACOUSTIC_SIMULATION
  use specfem_par_elastic, only: ELASTIC_SIMULATION
  use specfem_par_poroelastic, only: POROELASTIC_SIMULATION

  implicit none

  interface
    subroutine save_kernels_elastic(adios_handle, alphav_kl, alphah_kl, &
                                    betav_kl, betah_kl, eta_kl,         &
                                    rhop_kl, alpha_kl, beta_kl)

      integer(kind=8) :: adios_handle
      ! FIXME
      ! Break the CUSTOM_REAL stuff.
      ! put all this file in a module so interface is implicit
      ! OR
      ! redo what was done before SVN revision 22718
      !
      ! see other FIXME below (same than see one)
      real(kind=4), dimension(:,:,:,:), allocatable :: &
          alphav_kl,alphah_kl,betav_kl,betah_kl, &
          eta_kl, rhop_kl, alpha_kl, beta_kl
    end subroutine save_kernels_elastic
  end interface

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: alphav_kl, &
                                                            alphah_kl, &
                                                            betav_kl,  &
                                                            betah_kl,  &
                                                            eta_kl

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: rhop_kl,  &
                                                            alpha_kl, &
                                                            beta_kl

  integer(kind=8) :: adios_handle
  integer :: ier

  ! flag to save GLL weights
  logical,parameter :: SAVE_WEIGHTS = .false.

  if (ADIOS_FOR_KERNELS) then
    call define_kernel_adios_variables(adios_handle, SAVE_WEIGHTS)
  endif

  ! acoustic domains
  if( ACOUSTIC_SIMULATION ) then
    call save_kernels_acoustic(adios_handle)
  endif

  ! elastic domains
  if( ELASTIC_SIMULATION ) then
    ! allocates temporary transversely isotropic kernels
    if( ANISOTROPIC_KL ) then
      if( SAVE_TRANSVERSE_KL ) then
        allocate(alphav_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 alphah_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 betav_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 betah_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 eta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 stat=ier)
        if( ier /=0 ) stop 'error allocating arrays alphav_kl,...'

        ! derived kernels
        ! vp kernel
        allocate(alpha_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
        if( ier /= 0 ) stop 'error allocating array alpha_kl'
        ! vs kernel
        allocate(beta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
        if( ier /= 0 ) stop 'error allocating array beta_kl'
      endif
    else
      ! derived kernels
      ! vp kernel
      allocate(alpha_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if( ier /= 0 ) stop 'error allocating array alpha_kl'
      ! vs kernel
      allocate(beta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if( ier /= 0 ) stop 'error allocating array beta_kl'
      ! density prime kernel
      allocate(rhop_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if( ier /= 0 ) stop 'error allocating array rhop_kl'
    endif

    call save_kernels_elastic(adios_handle, alphav_kl, alphah_kl, &
                              betav_kl, betah_kl, eta_kl,         &
                              rhop_kl, alpha_kl, beta_kl)
  endif

  if( POROELASTIC_SIMULATION ) then
    call save_kernels_poroelastic(adios_handle)
  endif

  ! save weights for volume integration,
  ! in order to benchmark the kernels with analytical expressions
  if( SAVE_WEIGHTS ) then
    call save_weights_kernel()
  endif

  ! for noise simulations --- noise strength kernel
  if (NOISE_TOMOGRAPHY == 3) then
    call save_kernels_strength_noise(myrank,LOCAL_PATH,sigma_kl,NSPEC_AB)
  endif

  ! for preconditioner
  if ( APPROXIMATE_HESS_KL ) then
    call save_kernels_hessian(adios_handle)
  endif

  if (ADIOS_FOR_KERNELS) then
    call perform_write_adios_kernels(adios_handle)
  endif

  if (ELASTIC_SIMULATION) then
    ! frees temporary arrays
    if( ANISOTROPIC_KL ) then
      if( SAVE_TRANSVERSE_KL ) then
        deallocate(alphav_kl,alphah_kl,betav_kl,betah_kl,eta_kl)
        deallocate(alpha_kl,beta_kl)
      endif
    else
      deallocate(rhop_kl,alpha_kl,beta_kl)
    endif
  endif

end subroutine save_adjoint_kernels

!==============================================================================
!> Save weights for volume integration,
!! in order to benchmark the kernels with analytical expressions.
subroutine save_weights_kernel()
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer:: ispec,i,j,k,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: weights_kernel

  allocate(weights_kernel(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating array weights_kernel'
  do ispec = 1, NSPEC_AB
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            weights_kernel(i,j,k,ispec) = wxgll(i) * wygll(j) * wzgll(k) * jacobian(i,j,k,ispec)
          enddo ! i
        enddo ! j
      enddo ! k
  enddo ! ispec

  open(unit=IOUT,file=prname(1:len_trim(prname))//'weights_kernel.bin',status='unknown',form='unformatted',iostat=ier)
  if( ier /= 0 ) stop 'error opening file weights_kernel.bin'
  write(IOUT) weights_kernel
  close(IOUT)

  deallocate(weights_kernel,stat=ier)
  if( ier /= 0 ) stop 'error allocating array weights_kernel'
end subroutine save_weights_kernel

!==============================================================================
!> Save acoustic related kernels
subroutine save_kernels_acoustic(adios_handle)

  use specfem_par
  use specfem_par_acoustic

  implicit none

  integer(kind=8) :: adios_handle

  ! local parameters
  integer:: ispec,i,j,k,ier

  ! finalizes calculation of rhop, beta, alpha kernels
  do ispec = 1, NSPEC_AB

    ! acoustic simulations
    if( ispec_is_acoustic(ispec) ) then

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            ! rho prime kernel
            rhop_ac_kl(i,j,k,ispec) = rho_ac_kl(i,j,k,ispec) + kappa_ac_kl(i,j,k,ispec)

            ! vp kernel
            alpha_ac_kl(i,j,k,ispec) = 2._CUSTOM_REAL *  kappa_ac_kl(i,j,k,ispec)

          enddo
        enddo
      enddo

    endif ! acoustic

  enddo

  if (ADIOS_FOR_KERNELS) then
    call save_kernels_acoustic_adios(adios_handle)
  else
    ! save kernels to binary files
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rho_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rho_acoustic_kernel.bin'
    write(IOUT) rho_ac_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'kappa_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file kappa_acoustic_kernel.bin'
    write(IOUT) kappa_ac_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhop_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rhop_acoustic_kernel.bin'
    write(IOUT) rhop_ac_kl
    close(IOUT)

    open(unit=IOUT,file=prname(1:len_trim(prname))//'alpha_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file alpha_acoustic_kernel.bin'
    write(IOUT) alpha_ac_kl
    close(IOUT)

  endif
end subroutine save_kernels_acoustic

!==============================================================================
!> Save elastic related kernels
subroutine save_kernels_elastic(adios_handle, alphav_kl, alphah_kl, &
                                betav_kl, betah_kl, eta_kl,         &
                                rhop_kl, alpha_kl, beta_kl)

  use specfem_par
  use specfem_par_elastic

  implicit none

  interface
    subroutine save_kernels_elastic_adios(adios_handle, alphav_kl, alphah_kl, &
                                          betav_kl, betah_kl, eta_kl,         &
                                          rhop_kl, alpha_kl, beta_kl)

      integer(kind=8) :: adios_handle
      ! FIXME
      ! see other FIXME above.
      real(kind=4), dimension(:,:,:,:), allocatable :: &
          alphav_kl,alphah_kl,betav_kl,betah_kl, &
          eta_kl, rhop_kl, alpha_kl, beta_kl
    end subroutine save_kernels_elastic_adios
  end interface

  integer(kind=8) :: adios_handle
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    alphav_kl,alphah_kl,betav_kl,betah_kl, &
    eta_kl, rhop_kl, alpha_kl, beta_kl

  ! local parameters
  integer:: ispec,i,j,k,iglob,ier
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal

  ! Transverse isotropic paramters
  real(kind=CUSTOM_REAL) :: A,N,C,L,F,eta
  real(kind=CUSTOM_REAL), dimension(21) :: cijkl_kl_local
  real(kind=CUSTOM_REAL), dimension(5) :: an_kl


  ! finalizes calculation of rhop, beta, alpha kernels
  do ispec = 1, NSPEC_AB

    ! elastic simulations
    if( ispec_is_elastic(ispec) ) then

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)

            ! Store local material values
            rhol = rho_vs(i,j,k,ispec)*rho_vs(i,j,k,ispec) / mustore(i,j,k,ispec)
            mul = mustore(i,j,k,ispec)
            kappal = kappastore(i,j,k,ispec)

            if( ANISOTROPIC_KL ) then
              if( SAVE_TRANSVERSE_KL ) then
                cijkl_kl_local(:) = - cijkl_kl(:,i,j,k,ispec)

                ! Computes parameters for an isotropic model
                A = kappal + FOUR_THIRDS * mul
                C = A
                L = mul
                N = mul
                F = kappal - 2._CUSTOM_REAL/3._CUSTOM_REAL * mul
                eta = 1._CUSTOM_REAL

                ! note: cijkl_kl_local() is fully anisotropic C_ij kernel components (non-dimensionalized)
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

                an_kl(1) = cijkl_kl_local(1)+cijkl_kl_local(2)+cijkl_kl_local(7)  !A
                an_kl(2) = cijkl_kl_local(12)                                     !C
                an_kl(3) = -2*cijkl_kl_local(2)+cijkl_kl_local(21)                !N
                an_kl(4) = cijkl_kl_local(16)+cijkl_kl_local(19)                  !L
                an_kl(5) = cijkl_kl_local(3)+cijkl_kl_local(8)                    !F

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

                ! to check: isotropic kernels from transverse isotropic ones
                alpha_kl(i,j,k,ispec) = alphav_kl(i,j,k,ispec) &
                                                  + alphah_kl(i,j,k,ispec)
                beta_kl(i,j,k,ispec) = betav_kl(i,j,k,ispec) &
                                                  + betah_kl(i,j,k,ispec)
              endif ! SAVE_TRANSVERSE_KL

            else

              ! isotropic kernels

              ! isotropic adjoint kernels (see e.g. Tromp et al. 2005)
              ! for a parameterization: (rho,mu,kappa) "primary" kernels
              ! density kernel
              ! multiplies with rho
              rho_kl(i,j,k,ispec) = - rhol * rho_kl(i,j,k,ispec)

              ! shear modulus kernel
              mu_kl(i,j,k,ispec) = - 2._CUSTOM_REAL * mul * mu_kl(i,j,k,ispec)

              ! bulk modulus kernel
              kappa_kl(i,j,k,ispec) = - kappal * kappa_kl(i,j,k,ispec)

              ! for a parameterization: (rho,alpha,beta)
              ! density prime kernel
              rhop_kl(i,j,k,ispec) = rho_kl(i,j,k,ispec) + kappa_kl(i,j,k,ispec) + mu_kl(i,j,k,ispec)

              ! vs kernel
              beta_kl(i,j,k,ispec) = 2._CUSTOM_REAL * (mu_kl(i,j,k,ispec) &
                    - 4._CUSTOM_REAL * mul / (3._CUSTOM_REAL * kappal) * kappa_kl(i,j,k,ispec))

              ! vp kernel
              alpha_kl(i,j,k,ispec) = 2._CUSTOM_REAL * (1._CUSTOM_REAL &
                    + 4._CUSTOM_REAL * mul / (3._CUSTOM_REAL * kappal) ) * kappa_kl(i,j,k,ispec)

            endif

          enddo
        enddo
      enddo

    endif ! elastic

  enddo

  if (ADIOS_FOR_KERNELS) then
    call save_kernels_elastic_adios(adios_handle, alphav_kl, alphah_kl, &
                                      betav_kl, betah_kl, eta_kl,       &
                                      rhop_kl, alpha_kl, beta_kl)
  else
    if (ANISOTROPIC_KL) then

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

         ! transverse isotropic test kernels
         open(unit=IOUT,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
         write(IOUT)  alpha_kl
         close(IOUT)
         open(unit=IOUT,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
         write(IOUT)  beta_kl
         close(IOUT)

       else
         ! fully anisotropic kernels
         ! note: the C_ij and density kernels are not for relative perturbations (delta ln( m_i) = delta m_i / m_i),
         !          but absolute perturbations (delta m_i = m_i - m_0).
         ! Kappa and mu are for absolute perturbations, can be used to check with purely isotropic versions.
         open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
         write(IOUT)  - rho_kl
         close(IOUT)
         open(unit=IOUT,file=trim(prname)//'cijkl_kernel.bin',status='unknown',form='unformatted',action='write')
         write(IOUT) - cijkl_kl
         close(IOUT)

        endif

    else

      ! save kernels to binary files
      open(unit=IOUT,file=prname(1:len_trim(prname))//'rho_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if( ier /= 0 ) stop 'error opening file rho_kernel.bin'
      write(IOUT) rho_kl
      close(IOUT)

      open(unit=IOUT,file=prname(1:len_trim(prname))//'mu_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if( ier /= 0 ) stop 'error opening file mu_kernel.bin'
      write(IOUT) mu_kl
      close(IOUT)

      open(unit=IOUT,file=prname(1:len_trim(prname))//'kappa_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if( ier /= 0 ) stop 'error opening file kappa_kernel.bin'
      write(IOUT) kappa_kl
      close(IOUT)

      open(unit=IOUT,file=prname(1:len_trim(prname))//'rhop_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if( ier /= 0 ) stop 'error opening file rhop_kernel.bin'
      write(IOUT) rhop_kl
      close(IOUT)

      open(unit=IOUT,file=prname(1:len_trim(prname))//'beta_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if( ier /= 0 ) stop 'error opening file beta_kernel.bin'
      write(IOUT) beta_kl
      close(IOUT)

      open(unit=IOUT,file=prname(1:len_trim(prname))//'alpha_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if( ier /= 0 ) stop 'error opening file alpha_kernel.bin'
      write(IOUT) alpha_kl
      close(IOUT)
    endif

    if (SAVE_MOHO_MESH) then
      open(unit=IOUT,file=prname(1:len_trim(prname))//'moho_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if( ier /= 0 ) stop 'error opening file moho_kernel.bin'
      write(IOUT) moho_kl
      close(IOUT)
    endif
  endif

end subroutine save_kernels_elastic

!==============================================================================
!> Save poroelastic related kernels
subroutine save_kernels_poroelastic(adios_handle)

  use specfem_par
  use specfem_par_poroelastic

  implicit none

  integer(kind=8) :: adios_handle

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
    if( ispec_is_poroelastic(ispec) ) then

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
                      (1+rhol_f/rhol_bar)-1))/dd1**2 ) - &
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
                       (1+rhol_f/rhol_bar)-1))/dd1**2 )- &
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
                       1._CUSTOM_REAL)/dd1**2 )*Mb_kl(i,j,k,ispec) +&
                       ratio*rhol_f/C_biot*(cpIsquare-cpIIsquare) * (&
                       (2._CUSTOM_REAL*phil*rhol_bar* &
                       ratio/(tortl*rhol_f)+&
                       phil/tortl+rhol_bar/rhol_f)/dd1 - &
                       2._CUSTOM_REAL*phil/tortl*(phil/tortl*ratio+&
                       1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar/rhol_f*ratio)*((1._CUSTOM_REAL+&
                       rhol_bar/rhol_f- &
                       phil/tortl)*ratio+1._CUSTOM_REAL)/&
                       dd1**2 )*Cb_kl(i,j,k,ispec)
          enddo
        enddo
      enddo

    endif ! poroelastic

  enddo

  ! save kernels to binary files
  if (ADIOS_FOR_KERNELS) then
    call save_kernels_poroelastic_adios(adios_handle)
  else
    ! primary kernels
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhot_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rhot_primeporo_kernel.bin'
    write(IOUT) rhot_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhof_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rhof_primeporo_kernel.bin'
    write(IOUT) rhof_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'sm_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file sm_primeporo_kernel.bin'
    write(IOUT) sm_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'eta_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file eta_primeporo_kernel.bin'
    write(IOUT) eta_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'mufr_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file mufr_primeporo_kernel.bin'
    write(IOUT) mufr_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'B_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file B_primeporo_kernel.bin'
    write(IOUT) B_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'C_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file C_primeporo_kernel.bin'
    write(IOUT) C_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'M_primeporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file M_primeporo_kernel.bin'
    write(IOUT) M_kl
    close(IOUT)

    ! density kernels
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhob_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rhob_densityporo_kernel.bin'
    write(IOUT) rhob_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhofb_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rhofb_densityporo_kernel.bin'
    write(IOUT) rhofb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'phi_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file phi_densityporo_kernel.bin'
    write(IOUT) phi_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'mufrb_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file mufrb_densityporo_kernel.bin'
    write(IOUT) mufrb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'Bb_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file Bb_densityporo_kernel.bin'
    write(IOUT) Bb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'Cb_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file Cb_densityporo_kernel.bin'
    write(IOUT) Cb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'Mb_densityporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file Mb_densityporo_kernel.bin'
    write(IOUT) Mb_kl
    close(IOUT)

    ! wavespeed kernels
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhobb_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rhobb_waveporo_kernel.bin'
    write(IOUT) rhobb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rhofbb_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rhofbb_waveporo_kernel.bin'
    write(IOUT) rhofbb_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'phib_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file phib_waveporo_kernel.bin'
    write(IOUT) phib_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'cs_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file cs_waveporo_kernel.bin'
    write(IOUT) cs_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'cpI_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file cpI_waveporo_kernel.bin'
    write(IOUT) cpI_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'cpII_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file cpII_waveporo_kernel.bin'
    write(IOUT) cpII_kl
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'ratio_waveporo_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file ratio_waveporo_kernel.bin'
    write(IOUT) ratio_kl
    close(IOUT)

  endif
end subroutine save_kernels_poroelastic

!==============================================================================
!> Save hessians
subroutine save_kernels_hessian(adios_handle)

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  integer(kind=8) :: adios_handle

  integer :: ier

  ! acoustic domains
  if( ACOUSTIC_SIMULATION ) then
    ! scales approximate hessian
    hess_ac_kl(:,:,:,:) = 2._CUSTOM_REAL * hess_ac_kl(:,:,:,:)
  endif

  ! elastic domains
  if( ELASTIC_SIMULATION ) then
    ! scales approximate hessian
    hess_kl(:,:,:,:) = 2._CUSTOM_REAL * hess_kl(:,:,:,:)
  endif

  if (ADIOS_FOR_KERNELS) then
    call save_kernels_hessian_adios(adios_handle)
  else
    ! acoustic domains
    if( ACOUSTIC_SIMULATION ) then
      ! stores into file
      open(unit=IOUT,file=trim(prname)//'hess_acoustic_kernel.bin', &
            status='unknown',form='unformatted',action='write',iostat=ier)
      if( ier /= 0 ) stop 'error opening file hess_acoustic_kernel.bin'
      write(IOUT) hess_ac_kl
      close(IOUT)
    endif

    ! elastic domains
    if( ELASTIC_SIMULATION ) then
      ! stores into file
      open(unit=IOUT,file=trim(prname)//'hess_kernel.bin', &
            status='unknown',form='unformatted',action='write',iostat=ier)
      if( ier /= 0 ) stop 'error opening file hess_kernel.bin'
      write(IOUT) hess_kl
      close(IOUT)
    endif
  endif
end subroutine save_kernels_hessian

