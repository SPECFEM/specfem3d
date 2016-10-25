!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

  subroutine compute_kernel_integral_iso()

! computes volume element associated with points and calculates kernel integral

  use tomography_kernels_iso
  use tomography_model_iso
  implicit none

  ! jacobian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: jacobian
  real(kind=CUSTOM_REAL) :: volumel
  ! integration values
  real(kind=CUSTOM_REAL) :: kernel_integral_alpha,kernel_integral_beta,kernel_integral_rho
  real(kind=CUSTOM_REAL) :: integral_alpha_sum,integral_beta_sum,integral_rho_sum

  real(kind=CUSTOM_REAL) :: volume_glob,volume_glob_sum

  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_beta,norm_rho
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_beta_sum,norm_rho_sum

  ! root-mean square values
  real(kind=CUSTOM_REAL) :: rms_vp,rms_vs,rms_rho
  real(kind=CUSTOM_REAL) :: rms_vp_sum,rms_vs_sum,rms_rho_sum
  real(kind=CUSTOM_REAL) :: dvp,dvs,drho

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll
  ! array with all the weights in the cube
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  integer :: iglob
  integer :: i,j,k,ispec,ier

  ! user output
  if (myrank == 0) then
    print *
    print *,'***********'
    print *,'statistics:'
    print *,'***********'
    print *
  endif

  ! allocates array
  allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating jacobian array'

  ! GLL points
  wgll_cube = 0.0
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        wgll_cube(i,j,k) = wxgll(i)*wygll(j)*wzgll(k)
      enddo
    enddo
  enddo

  ! builds jacobian
  call compute_jacobian(jacobian)

  ! volume associated with global point
  volume_glob = 0._CUSTOM_REAL
  kernel_integral_alpha = 0._CUSTOM_REAL
  kernel_integral_beta = 0._CUSTOM_REAL
  kernel_integral_rho = 0._CUSTOM_REAL
  norm_bulk = 0._CUSTOM_REAL
  norm_beta = 0._CUSTOM_REAL
  norm_rho = 0._CUSTOM_REAL
  rms_vp = 0._CUSTOM_REAL
  rms_vs = 0._CUSTOM_REAL
  rms_rho = 0._CUSTOM_REAL

  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          if (iglob == 0) then
            print *,'iglob zero',i,j,k,ispec
            print *
            print *,'ibool:',ispec
            print *,ibool(:,:,:,ispec)
            print *
            call exit_MPI(myrank,'Error ibool')
          endif

          ! volume associated with GLL point
          volumel = jacobian(i,j,k,ispec) * wgll_cube(i,j,k)
          volume_glob = volume_glob + volumel

          ! kernel integration: for each element
          kernel_integral_alpha = kernel_integral_alpha &
                                 + volumel * kernel_bulk(i,j,k,ispec)

          kernel_integral_beta = kernel_integral_beta &
                                 + volumel * kernel_beta(i,j,k,ispec)

          kernel_integral_rho = kernel_integral_rho &
                                 + volumel * kernel_rho(i,j,k,ispec)

          ! gradient vector norm sqrt(  v^T * v )
          norm_bulk = norm_bulk + kernel_bulk(i,j,k,ispec)**2
          norm_beta = norm_beta + kernel_beta(i,j,k,ispec)**2
          norm_rho = norm_rho + kernel_rho(i,j,k,ispec)**2

          ! checks number (isNaN)
          if (kernel_integral_alpha /= kernel_integral_alpha) then
            print *,'Error NaN: ',kernel_integral_alpha
            print *,'rank:',myrank
            print *,'i,j,k,ispec:',i,j,k,ispec
            print *,'volumel: ',volumel,'kernel_bulk:',kernel_bulk(i,j,k,ispec)
            call exit_MPI(myrank,'Error NaN')
          endif

          ! root-mean square
          ! integrates relative perturbations ( dv / v  using logarithm ) squared
          dvp = log( model_vp_new(i,j,k,ispec) / model_vp(i,j,k,ispec) ) ! alphav
          rms_vp = rms_vp + volumel * dvp*dvp

          dvs = log( model_vs_new(i,j,k,ispec) / model_vs(i,j,k,ispec) ) ! betav
          rms_vs = rms_vs + volumel * dvs*dvs

          drho = log( model_rho_new(i,j,k,ispec) / model_rho(i,j,k,ispec) ) ! rho
          rms_rho = rms_rho + volumel * drho*drho

        enddo
      enddo
    enddo
  enddo

  ! statistics
  ! (note: sum_all_cr() will only return valid results to master process)
  ! kernel integration: for whole volume
  call sum_all_cr(kernel_integral_alpha,integral_alpha_sum)
  call sum_all_cr(kernel_integral_beta,integral_beta_sum)
  call sum_all_cr(kernel_integral_rho,integral_rho_sum)
  call sum_all_cr(volume_glob,volume_glob_sum)

  if (myrank == 0) then
    print *,'integral kernels:'
    print *,'  a   : ',integral_alpha_sum
    print *,'  beta: ',integral_beta_sum
    print *,'  rho : ',integral_rho_sum
    print *
    print *,'  total volume:',volume_glob_sum
    print *
    if (volume_glob_sum < 1.e-25) stop 'Error zero total volume'
  endif

  ! norms: for whole volume
  call sum_all_cr(norm_bulk,norm_bulk_sum)
  call sum_all_cr(norm_beta,norm_beta_sum)
  call sum_all_cr(norm_rho,norm_rho_sum)

  if (myrank == 0) then
    norm_bulk = sqrt(norm_bulk_sum)
    norm_beta = sqrt(norm_beta_sum)
    norm_rho = sqrt(norm_rho_sum)

    print *,'norm kernels:'
    print *,'  a   : ',norm_bulk
    print *,'  beta: ',norm_beta
    print *,'  rho : ',norm_rho
    print *
  endif

  ! root-mean square
  call sum_all_cr(rms_vp,rms_vp_sum)
  call sum_all_cr(rms_vs,rms_vs_sum)
  call sum_all_cr(rms_rho,rms_rho_sum)

  if (myrank == 0) then
    rms_vp  = sqrt( rms_vp_sum / volume_glob_sum )
    rms_vs  = sqrt( rms_vs_sum / volume_glob_sum )
    rms_rho = sqrt( rms_rho_sum / volume_glob_sum )

    print *,'root-mean square of perturbations:'
    print *,'  vp : ',rms_vp
    print *,'  vs : ',rms_vs
    print *,'  rho: ',rms_rho
    print *
  endif
  call synchronize_all()

  ! frees memory
  deallocate(jacobian)

  end subroutine compute_kernel_integral_iso

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernel_integral_tiso()

! computes volume element associated with points

  use tomography_kernels_tiso
  use tomography_model_tiso
  implicit none

  ! jacobian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: jacobian
  real(kind=CUSTOM_REAL) :: volumel

  ! integration values
  real(kind=CUSTOM_REAL) :: integral_bulk_sum,integral_betav_sum, &
    integral_betah_sum,integral_eta_sum
  real(kind=CUSTOM_REAL) :: integral_bulk,integral_betav, &
    integral_betah,integral_eta
  real(kind=CUSTOM_REAL) :: volume_glob,volume_glob_sum

  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_betav,norm_betah,norm_eta
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_betav_sum, &
    norm_betah_sum,norm_eta_sum

  ! root-mean square values
  real(kind=CUSTOM_REAL) :: rms_vpv,rms_vph,rms_vsv,rms_vsh,rms_eta,rms_rho
  real(kind=CUSTOM_REAL) :: rms_vpv_sum,rms_vph_sum,rms_vsv_sum,rms_vsh_sum, &
    rms_eta_sum,rms_rho_sum
  real(kind=CUSTOM_REAL) :: dvpv,dvph,dvsv,dvsh,deta,drho

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll
  ! array with all the weights in the cube
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  integer :: iglob
  integer :: i,j,k,ispec,ier

  ! user output
  if (myrank == 0) then
    print *
    print *,'***********'
    print *,'statistics:'
    print *,'***********'
    print *
  endif

  ! allocates array
  allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating jacobian array'

  ! GLL points
  wgll_cube = 0.0d0
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        wgll_cube(i,j,k) = wxgll(i)*wygll(j)*wzgll(k)
      enddo
    enddo
  enddo

  ! builds jacobian
  call compute_jacobian(jacobian)

  ! volume associated with global point
  volume_glob = 0._CUSTOM_REAL
  integral_bulk = 0._CUSTOM_REAL
  integral_betav = 0._CUSTOM_REAL
  integral_betah = 0._CUSTOM_REAL
  integral_eta = 0._CUSTOM_REAL
  norm_bulk = 0._CUSTOM_REAL
  norm_betav = 0._CUSTOM_REAL
  norm_betah = 0._CUSTOM_REAL
  norm_eta = 0._CUSTOM_REAL
  rms_vpv = 0._CUSTOM_REAL
  rms_vph = 0._CUSTOM_REAL
  rms_vsv = 0._CUSTOM_REAL
  rms_vsh = 0._CUSTOM_REAL
  rms_eta = 0._CUSTOM_REAL
  rms_rho = 0._CUSTOM_REAL
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          if (iglob == 0) then
            print *,'iglob zero',i,j,k,ispec
            print *
            print *,'ibool:',ispec
            print *,ibool(:,:,:,ispec)
            print *
            call exit_MPI(myrank,'Error ibool')
          endif

          ! volume associated with GLL point
          volumel = jacobian(i,j,k,ispec)*wgll_cube(i,j,k)
          volume_glob = volume_glob + volumel

          ! kernel integration: for each element
          integral_bulk = integral_bulk &
                                 + volumel * kernel_bulk(i,j,k,ispec)

          integral_betav = integral_betav &
                                 + volumel * kernel_betav(i,j,k,ispec)

          integral_betah = integral_betah &
                                 + volumel * kernel_betah(i,j,k,ispec)

          integral_eta = integral_eta &
                                 + volumel * kernel_eta(i,j,k,ispec)

          ! gradient vector norm sqrt(  v^T * v )
          norm_bulk = norm_bulk + kernel_bulk(i,j,k,ispec)*kernel_bulk(i,j,k,ispec)
          norm_betav = norm_betav + kernel_betav(i,j,k,ispec)*kernel_betav(i,j,k,ispec)
          norm_betah = norm_betah + kernel_betah(i,j,k,ispec)*kernel_betah(i,j,k,ispec)
          norm_eta = norm_eta + kernel_eta(i,j,k,ispec)*kernel_eta(i,j,k,ispec)

          ! checks number
          if (integral_bulk /= integral_bulk) then
            print *,'Error NaN: ',integral_bulk
            print *,'rank:',myrank
            print *,'i,j,k,ispec:',i,j,k,ispec
            print *,'volumel: ',volumel,'kernel_bulk:',kernel_bulk(i,j,k,ispec)
            call exit_MPI(myrank,'Error NaN')
          endif

          ! root-mean square
          ! integrates relative perturbations ( dv / v  using logarithm ) squared
          dvpv = log( model_vpv_new(i,j,k,ispec) / model_vpv(i,j,k,ispec) ) ! alphav
          rms_vpv = rms_vpv + volumel * dvpv*dvpv

          dvph = log( model_vph_new(i,j,k,ispec) / model_vph(i,j,k,ispec) ) ! alphah
          rms_vph = rms_vph + volumel * dvph*dvph

          dvsv = log( model_vsv_new(i,j,k,ispec) / model_vsv(i,j,k,ispec) ) ! betav
          rms_vsv = rms_vsv + volumel * dvsv*dvsv

          dvsh = log( model_vsh_new(i,j,k,ispec) / model_vsh(i,j,k,ispec) ) ! betah
          rms_vsh = rms_vsh + volumel * dvsh*dvsh

          deta = log( model_eta_new(i,j,k,ispec) / model_eta(i,j,k,ispec) ) ! eta
          rms_eta = rms_eta + volumel * deta*deta

          drho = log( model_rho_new(i,j,k,ispec) / model_rho(i,j,k,ispec) ) ! rho
          rms_rho = rms_rho + volumel * drho*drho

        enddo
      enddo
    enddo
  enddo

  ! statistics
  ! (note: sum_all_cr() will only return valid results to master process)
  ! kernel integration: for whole volume
  call sum_all_cr(integral_bulk,integral_bulk_sum)
  call sum_all_cr(integral_betav,integral_betav_sum)
  call sum_all_cr(integral_betah,integral_betah_sum)
  call sum_all_cr(integral_eta,integral_eta_sum)
  call sum_all_cr(volume_glob,volume_glob_sum)

  if (myrank == 0) then
    print *,'integral kernels:'
    print *,'  bulk : ',integral_bulk_sum
    print *,'  betav : ',integral_betav_sum
    print *,'  betah : ',integral_betah_sum
    print *,'  eta : ',integral_eta_sum
    print *
    print *,'  total volume:',volume_glob_sum
    print *
    if (volume_glob_sum < 1.e-25) stop 'Error zero total volume'
  endif

  ! norms: for whole volume
  call sum_all_cr(norm_bulk,norm_bulk_sum)
  call sum_all_cr(norm_betav,norm_betav_sum)
  call sum_all_cr(norm_betah,norm_betah_sum)
  call sum_all_cr(norm_eta,norm_eta_sum)

  if (myrank == 0) then
    norm_bulk = sqrt(norm_bulk_sum)
    norm_betav = sqrt(norm_betav_sum)
    norm_betah = sqrt(norm_betah_sum)
    norm_eta = sqrt(norm_eta_sum)

    print *,'norm kernels:'
    print *,'  bulk : ',norm_bulk
    print *,'  betav : ',norm_betav
    print *,'  betah : ',norm_betah
    print *,'  eta : ',norm_eta
    print *
  endif

  ! root-mean square
  call sum_all_cr(rms_vpv,rms_vpv_sum)
  call sum_all_cr(rms_vph,rms_vph_sum)
  call sum_all_cr(rms_vsv,rms_vsv_sum)
  call sum_all_cr(rms_vsh,rms_vsh_sum)
  call sum_all_cr(rms_eta,rms_eta_sum)
  call sum_all_cr(rms_rho,rms_rho_sum)

  if (myrank == 0) then
    rms_vpv = sqrt( rms_vpv_sum / volume_glob_sum )
    rms_vph = sqrt( rms_vph_sum / volume_glob_sum )
    rms_vsv = sqrt( rms_vsv_sum / volume_glob_sum )
    rms_vsh = sqrt( rms_vsh_sum / volume_glob_sum )
    rms_eta = sqrt( rms_eta_sum / volume_glob_sum )
    rms_rho = sqrt( rms_rho_sum / volume_glob_sum )

    print *,'root-mean square of perturbations:'
    print *,'  vpv : ',rms_vpv
    print *,'  vph : ',rms_vph
    print *,'  vsv : ',rms_vsv
    print *,'  vsh : ',rms_vsh
    print *,'  eta : ',rms_eta
    print *,'  rho : ',rms_rho
    print *
  endif
  call synchronize_all()

  ! frees memory
  deallocate(jacobian)

  end subroutine compute_kernel_integral_tiso


!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernel_integral_tiso_iso()

! computes volume element associated with points

  use tomography_kernels_iso
  use tomography_model_tiso
  implicit none
  ! jacobian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: jacobian
  real(kind=CUSTOM_REAL) :: volumel

  ! integration values
  real(kind=CUSTOM_REAL) :: integral_bulk_sum,integral_beta_sum,integral_rho_sum
  real(kind=CUSTOM_REAL) :: integral_bulk,integral_beta,integral_rho
  real(kind=CUSTOM_REAL) :: volume_glob,volume_glob_sum

  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_beta,norm_rho
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_beta_sum,norm_rho_sum

  ! root-mean square values
  real(kind=CUSTOM_REAL) :: rms_vpv,rms_vph,rms_vsv,rms_vsh,rms_eta,rms_rho
  real(kind=CUSTOM_REAL) :: rms_vpv_sum,rms_vph_sum,rms_vsv_sum,rms_vsh_sum, &
    rms_eta_sum,rms_rho_sum
  real(kind=CUSTOM_REAL) :: dvpv,dvph,dvsv,dvsh,deta,drho

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll
  ! array with all the weights in the cube
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  integer :: iglob
  integer :: i,j,k,ispec,ier

  ! user output
  if (myrank == 0) then
    print *
    print *,'***********'
    print *,'statistics:'
    print *,'***********'
    print *
  endif

  ! allocates array
  allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating jacobian array'

  ! GLL points
  wgll_cube = 0.0d0
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        wgll_cube(i,j,k) = wxgll(i)*wygll(j)*wzgll(k)
      enddo
    enddo
  enddo

  ! builds jacobian
  call compute_jacobian(jacobian)

  ! volume associated with global point
  volume_glob = 0._CUSTOM_REAL
  integral_bulk = 0._CUSTOM_REAL
  integral_beta = 0._CUSTOM_REAL
  integral_rho = 0._CUSTOM_REAL
  norm_bulk = 0._CUSTOM_REAL
  norm_beta = 0._CUSTOM_REAL
  norm_rho = 0._CUSTOM_REAL
  rms_vpv = 0._CUSTOM_REAL
  rms_vph = 0._CUSTOM_REAL
  rms_vsv = 0._CUSTOM_REAL
  rms_vsh = 0._CUSTOM_REAL
  rms_eta = 0._CUSTOM_REAL
  rms_rho = 0._CUSTOM_REAL

  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          if (iglob == 0) then
            print *,'iglob zero',i,j,k,ispec
            print *
            print *,'ibool:',ispec
            print *,ibool(:,:,:,ispec)
            print *
            call exit_MPI(myrank,'Error ibool')
          endif

          ! volume associated with GLL point
          volumel = jacobian(i,j,k,ispec)*wgll_cube(i,j,k)
          volume_glob = volume_glob + volumel

          ! kernel integration: for each element
          integral_bulk = integral_bulk &
                                 + volumel * kernel_bulk(i,j,k,ispec)

          integral_beta = integral_beta &
                                 + volumel * kernel_beta(i,j,k,ispec)

          integral_rho = integral_rho &
                                 + volumel * kernel_rho(i,j,k,ispec)

          ! gradient vector norm sqrt(  v^T * v )
          norm_bulk = norm_bulk + kernel_bulk(i,j,k,ispec)**2
          norm_beta = norm_beta + kernel_beta(i,j,k,ispec)**2
          norm_rho = norm_rho + kernel_rho(i,j,k,ispec)**2

          ! checks number
          if (integral_bulk /= integral_bulk) then
            print *,'Error NaN: ',integral_bulk
            print *,'rank:',myrank
            print *,'i,j,k,ispec:',i,j,k,ispec
            print *,'volumel: ',volumel,'kernel_bulk:',kernel_bulk(i,j,k,ispec)
            call exit_MPI(myrank,'Error NaN')
          endif

          ! root-mean square
          ! integrates relative perturbations ( dv / v  using logarithm ) squared
          dvpv = log( model_vpv_new(i,j,k,ispec) / model_vpv(i,j,k,ispec) ) ! alphav
          rms_vpv = rms_vpv + volumel * dvpv*dvpv

          dvph = log( model_vph_new(i,j,k,ispec) / model_vph(i,j,k,ispec) ) ! alphah
          rms_vph = rms_vph + volumel * dvph*dvph

          dvsv = log( model_vsv_new(i,j,k,ispec) / model_vsv(i,j,k,ispec) ) ! betav
          rms_vsv = rms_vsv + volumel * dvsv*dvsv

          dvsh = log( model_vsh_new(i,j,k,ispec) / model_vsh(i,j,k,ispec) ) ! betah
          rms_vsh = rms_vsh + volumel * dvsh*dvsh

          deta = log( model_eta_new(i,j,k,ispec) / model_eta(i,j,k,ispec) ) ! eta
          rms_eta = rms_eta + volumel * deta*deta

          drho = log( model_rho_new(i,j,k,ispec) / model_rho(i,j,k,ispec) ) ! rho
          rms_rho = rms_rho + volumel * drho*drho

        enddo
      enddo
    enddo
  enddo

  ! statistics
  ! (note: sum_all_cr() will only return valid results to master process)
  ! kernel integration: for whole volume
  call sum_all_cr(integral_bulk,integral_bulk_sum)
  call sum_all_cr(integral_beta,integral_beta_sum)
  call sum_all_cr(integral_rho,integral_rho_sum)
  call sum_all_cr(volume_glob,volume_glob_sum)

  if (myrank == 0) then
    print *,'integral kernels:'
    print *,'  a   : ',integral_bulk_sum
    print *,'  beta: ',integral_beta_sum
    print *,'  rho : ',integral_rho_sum
    print *
    print *,'  total volume:',volume_glob_sum
    print *
    if (volume_glob_sum < 1.e-25) stop 'Error zero total volume'
  endif

  ! norms: for whole volume
  call sum_all_cr(norm_bulk,norm_bulk_sum)
  call sum_all_cr(norm_beta,norm_beta_sum)
  call sum_all_cr(norm_rho,norm_rho_sum)

  if (myrank == 0) then
    norm_bulk = sqrt(norm_bulk_sum)
    norm_beta = sqrt(norm_beta_sum)
    norm_rho = sqrt(norm_rho_sum)

    print *,'norm kernels:'
    print *,'  a : ',norm_bulk
    print *,'  beta : ',norm_beta
    print *,'  rho : ',norm_rho
    print *
  endif

  ! root-mean square
  call sum_all_cr(rms_vpv,rms_vpv_sum)
  call sum_all_cr(rms_vph,rms_vph_sum)
  call sum_all_cr(rms_vsv,rms_vsv_sum)
  call sum_all_cr(rms_vsh,rms_vsh_sum)
  call sum_all_cr(rms_eta,rms_eta_sum)
  call sum_all_cr(rms_rho,rms_rho_sum)

  if (myrank == 0) then
    rms_vpv = sqrt( rms_vpv_sum / volume_glob_sum )
    rms_vph = sqrt( rms_vph_sum / volume_glob_sum )
    rms_vsv = sqrt( rms_vsv_sum / volume_glob_sum )
    rms_vsh = sqrt( rms_vsh_sum / volume_glob_sum )
    rms_eta = sqrt( rms_eta_sum / volume_glob_sum )
    rms_rho = sqrt( rms_rho_sum / volume_glob_sum )

    print *,'root-mean square of perturbations:'
    print *,'  vpv : ',rms_vpv
    print *,'  vph : ',rms_vph
    print *,'  vsv : ',rms_vsv
    print *,'  vsh : ',rms_vsh
    print *,'  eta : ',rms_eta
    print *,'  rho : ',rms_rho
    print *
  endif
  call synchronize_all()

  ! frees memory
  deallocate(jacobian)

  end subroutine compute_kernel_integral_tiso_iso

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_jacobian(jacobian)

! computes volume element associated with points

  use tomography_par, only: CUSTOM_REAL,NSPEC,NGLOB,NGLLX,NGLLY,NGLLZ,IIN,myrank,MAX_STRING_LEN,REG

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: jacobian

  ! local parameters
  ! dummy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: dummy_sem
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: dummy
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: dummy_ibool
  integer :: ival,ier
  character(len=MAX_STRING_LEN) :: m_file

  ! initializes
  jacobian(:,:,:,:) = 0.0_CUSTOM_REAL

  ! reads jacobian
  write(m_file,'(a,i6.6,a)') 'topo/proc',myrank,trim(REG)//'external_mesh.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif

  read(IIN) ival !nspec
  if (ival /= nspec) call exit_mpi(myrank,'Error invalid nspec value in external_mesh.bin')
  read(IIN) ival !nglob
  if (ival /= nglob) call exit_mpi(myrank,'Error invalid nspec value in external_mesh.bin')

  read(IIN) dummy_ibool ! ibool

  read(IIN) dummy ! x
  read(IIN) dummy ! y
  read(IIN) dummy ! z

  read(IIN) dummy_sem ! xix
  read(IIN) dummy_sem ! xiy
  read(IIN) dummy_sem ! xiz
  read(IIN) dummy_sem ! etax
  read(IIN) dummy_sem ! etay
  read(IIN) dummy_sem ! etaz
  read(IIN) dummy_sem ! gammax
  read(IIN) dummy_sem ! gammay
  read(IIN) dummy_sem ! gammaz
  read(IIN) jacobian

  close(IIN)

  end subroutine compute_jacobian

