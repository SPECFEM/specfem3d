!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  subroutine save_adjoint_kernels()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none
  integer:: ispec,i,j,k,iglob,ier
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: weights_kernel
  ! flag to save GLL weights
  logical,parameter :: SAVE_WEIGHTS = .true.

  ! finalizes calculation of rhop, beta, alpha kernels
  do ispec = 1, NSPEC_AB

    ! elastic simulations
    if( ispec_is_elastic(ispec) ) then

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)

            ! isotropic adjoint kernels (see e.g. Tromp et al. 2005)
            rhol = rho_vs(i,j,k,ispec)**2 / mustore(i,j,k,ispec)
            mul = mustore(i,j,k,ispec)
            kappal = kappastore(i,j,k,ispec)

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

            ! for a parameterization: (rho,bulk, beta)
            ! where bulk wave speed is c = sqrt( kappa / rho)
            ! note: rhoprime is the same as for (rho,alpha,beta) parameterization
            !bulk_c_kl_crust_mantle(i,j,k,ispec) = 2._CUSTOM_REAL * kappa_kl(i,j,k,ispec)
            !bulk_beta_kl_crust_mantle(i,j,k,ispec ) = 2._CUSTOM_REAL * mu_kl(i,j,k,ispec)

          enddo
        enddo
      enddo

    endif ! elastic

    ! acoustic simulations
    if( ispec_is_acoustic(ispec) ) then

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            ! rho prime kernel
            rhop_ac_kl(i,j,k,ispec) = rho_ac_kl(i,j,k,ispec) + kappa_ac_kl(i,j,k,ispec)

            ! kappa kernel
            alpha_ac_kl(i,j,k,ispec) = TWO *  kappa_ac_kl(i,j,k,ispec)
          enddo
        enddo
      enddo

    endif ! acoustic


  enddo

  ! save kernels to binary files
  if( ELASTIC_SIMULATION ) then
    open(unit=27,file=prname(1:len_trim(prname))//'rho_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rho_kernel.bin'
    write(27) rho_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'mu_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file mu_kernel.bin'
    write(27) mu_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'kappa_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file kappa_kernel.bin'
    write(27) kappa_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'rhop_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rhop_kernel.bin'
    write(27) rhop_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'beta_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file beta_kernel.bin'
    write(27) beta_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'alpha_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file alpha_kernel.bin'
    write(27) alpha_kl
    close(27)

    if (SAVE_MOHO_MESH) then
      open(unit=27,file=prname(1:len_trim(prname))//'moho_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if( ier /= 0 ) stop 'error opening file moho_kernel.bin'
      write(27) moho_kl
      close(27)
    endif

  endif


  ! save kernels to binary files
  if( ACOUSTIC_SIMULATION ) then
    open(unit=27,file=prname(1:len_trim(prname))//'rho_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rho_acoustic_kernel.bin'
    write(27) rho_ac_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'kappa_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file kappa_acoustic_kernel.bin'
    write(27) kappa_ac_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'rhop_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file rhop_acoustic_kernel.bin'
    write(27) rhop_ac_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'alpha_acoustic_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file alpha_acoustic_kernel.bin'
    write(27) alpha_ac_kl
    close(27)

  endif

  ! save weights for volume integration, in order to benchmark the kernels with analytical expressions
  if( SAVE_WEIGHTS ) then
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
    open(unit=27,file=prname(1:len_trim(prname))//'weights_kernel.bin',status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening file weights_kernel.bin'
    write(27) weights_kernel
    close(27)
  endif

  ! for noise simulations --- noise strength kernel
  if (NOISE_TOMOGRAPHY == 3) then
    call save_kernels_strength_noise(myrank,LOCAL_PATH,sigma_kl,NSPEC_AB)
  endif

  ! for preconditioner
  if ( APPROXIMATE_HESS_KL ) then
    call save_kernels_hessian()
  endif
  
  end subroutine save_adjoint_kernels
  
!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_hessian()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none
  integer :: ier
  
  ! acoustic domains
  if( ACOUSTIC_SIMULATION ) then
    ! scales approximate hessian
    hess_ac_kl(:,:,:,:) = 2._CUSTOM_REAL * hess_ac_kl(:,:,:,:)

    ! stores into file
    open(unit=27,file=trim(prname)//'hess_acoustic_kernel.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) stop 'error opening file hess_acoustic_kernel.bin'
    write(27) hess_ac_kl
    close(27)
  endif  

  ! elastic domains
  if( ELASTIC_SIMULATION ) then
    ! scales approximate hessian
    hess_kl(:,:,:,:) = 2._CUSTOM_REAL * hess_kl(:,:,:,:)

    ! stores into file
    open(unit=27,file=trim(prname)//'hess_kernel.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) stop 'error opening file hess_kernel.bin'
    write(27) hess_kl
    close(27)
  endif  
  
  end subroutine save_kernels_hessian
  
