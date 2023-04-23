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

  subroutine get_cg_direction_tiso()

! calculates TI gradient based on a conjugate gradient method
!
! based on: Tarantola, Inverse problem theory, 2005.
!                  section 6.22.7 conjugate directions, page 217.
!                  formula for alpha_n based on Polak & Ribiere (1969)
!
! note: we use a preconditioner F_0 = 1, thus lambda_n = gamma_n in (6.322)
!          and use gamma_n as the smoothed kernel (for bulk_c, bulk_betav,..).
!
!          however, one could see smoothing as preconditioner F_0, thus
!          gamma_n would be un-smoothed kernel and lambda_n would be smoothed one...
!          i'm not sure if this makes a difference.

  use tomography_kernels_tiso
  use tomography_kernels_tiso_cg
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha_bulk,alpha_betav,alpha_betah,alpha_eta,alpha_all
  real(kind=CUSTOM_REAL) :: minmax(4),depthmax(2),depthmax_radius(2),max
  real(kind=CUSTOM_REAL) :: r,rmax_vsv,rmax_vsh,depthmax_depth
  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_betav,norm_betah,norm_eta
  real(kind=CUSTOM_REAL) :: norm_bulk_old,norm_betav_old,norm_betah_old,norm_eta_old
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_betav_sum, &
    norm_betah_sum,norm_eta_sum
  real(kind=CUSTOM_REAL) :: min_vsv,min_vsh,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk
  integer :: maxindex(1)
  real(kind=CUSTOM_REAL) :: ratio_bulk,ratio_betav,ratio_betah,ratio_eta
  integer :: iglob
  integer :: i,j,k,ispec,ier

  ! allocate arrays for storing gradient
  ! transversely isotropic arrays
  allocate(model_dbulk(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1089')
  allocate(model_dbetav(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1090')
  allocate(model_dbetah(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1091')
  allocate(model_deta(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1092')
  if (ier /= 0) stop 'error allocating gradient arrays'

  ! initializes arrays
  model_dbulk = 0.0_CUSTOM_REAL
  model_dbetav = 0.0_CUSTOM_REAL
  model_dbetah = 0.0_CUSTOM_REAL
  model_deta = 0.0_CUSTOM_REAL

  ! old kernel/gradient
  ! length ( gamma_(n-1)^T * lambda_(n-1) )
  norm_bulk_old = sum( kernel_bulk_old * kernel_bulk_old )
  norm_betav_old = sum( kernel_betav_old * kernel_betav_old )
  norm_betah_old = sum( kernel_betah_old * kernel_betah_old )
  norm_eta_old = sum( kernel_eta_old * kernel_eta_old )

  call sum_all_cr(norm_bulk_old,norm_bulk_sum)
  call sum_all_cr(norm_betav_old,norm_betav_sum)
  call sum_all_cr(norm_betah_old,norm_betah_sum)
  call sum_all_cr(norm_eta_old,norm_eta_sum)

  ! don't use square root, just take gamma^T * gamma
  norm_bulk_old = norm_bulk_sum
  norm_betav_old = norm_betav_sum
  norm_betah_old = norm_betah_sum
  norm_eta_old = norm_eta_sum

  if (myrank == 0) then
    print *,'norm squared old gradient:'
    print *,'  bulk : ',norm_bulk_old
    print *,'  betav: ',norm_betav_old
    print *,'  betah: ',norm_betah_old
    print *,'  eta  : ',norm_eta_old
    print *

    ! checks lengths
    if (norm_bulk_old < 1.e-22) call exit_mpi(myrank,'norm old gradient bulk is zero')
    if (norm_betav_old < 1.e-22) call exit_mpi(myrank,'norm old gradient betav is zero')
    if (norm_betah_old < 1.e-22) call exit_mpi(myrank,'norm old gradient betah is zero')
    if (norm_eta_old < 1.e-22) call exit_mpi(myrank,'norm old gradient eta is zero')
  endif

  ! Powell, 1977: checks orthogonality between old and new gradient
  ! gets length of ( gamma_(n-1)^T * gamma_n )
  norm_bulk = sum( kernel_bulk_old * kernel_bulk )
  norm_betav = sum( kernel_betav_old * kernel_betav )
  norm_betah = sum( kernel_betah_old * kernel_betah )
  norm_eta = sum( kernel_eta_old * kernel_eta )

  call sum_all_cr(norm_bulk,norm_bulk_sum)
  call sum_all_cr(norm_betav,norm_betav_sum)
  call sum_all_cr(norm_betah,norm_betah_sum)
  call sum_all_cr(norm_eta,norm_eta_sum)

  if (myrank == 0) then
    ! ratio:  ( g_n * g_n-1) / ( g_n-1 * g_n-1)
    ratio_bulk = norm_bulk_sum / norm_bulk_old
    ratio_betav = norm_betav_sum / norm_betav_old
    ratio_betah = norm_betah_sum / norm_betah_old
    ratio_eta = norm_eta_sum / norm_eta_old

    ! if ratio > 0.2 (empirical threshold value), then one should restart with a steepest descent
    print *,'Powell ratio: ( > 0.2 then restart with steepest descent)'
    print *,'  bulk : ',ratio_bulk
    print *,'  betav: ',ratio_betav
    print *,'  betah: ',ratio_betah
    print *,'  eta  : ',ratio_eta
    print *
    if (ratio_bulk > 0.2 .and. ratio_betav > 0.2 .and. ratio_betah > 0.2 &
      .and. ratio_eta > 0.2) then
      print *,'  critical ratio found!'
      print *
      print *,'****************'
      print *
      print *,'  Please consider doing a steepest descent instead cg...'
      print *
      print *,'****************'
    endif
  endif


  ! difference kernel/gradient
  ! length ( ( gamma_n - gamma_(n-1))^T * lambda_n )
  norm_bulk = sum( (kernel_bulk - kernel_bulk_old) * kernel_bulk )
  norm_betav = sum( (kernel_betav - kernel_betav_old) * kernel_betav )
  norm_betah = sum( (kernel_betah - kernel_betah_old) * kernel_betah )
  norm_eta = sum( (kernel_eta - kernel_eta_old) * kernel_eta )

  call sum_all_cr(norm_bulk,norm_bulk_sum)
  call sum_all_cr(norm_betav,norm_betav_sum)
  call sum_all_cr(norm_betah,norm_betah_sum)
  call sum_all_cr(norm_eta,norm_eta_sum)

  ! don't take square root, since norm_bulk_sum could be negative
  ! just use (gamma_n - gamma_n-1)^T * lambda_n
  norm_bulk = norm_bulk_sum
  norm_betav = norm_betav_sum
  norm_betah = norm_betah_sum
  norm_eta = norm_eta_sum

  if (myrank == 0) then
    print *,'norm squared difference gradient:'
    print *,'  bulk : ',norm_bulk
    print *,'  betav: ',norm_betav
    print *,'  betah: ',norm_betah
    print *,'  eta  : ',norm_eta
    print *
  endif

  ! calculates ratio based on Polak & Ribiere (1969)
  if (myrank == 0) then
    if (USE_SEPARATE_CG_STEPLENGTHS) then
      ! calculates steplength alpha for each parameter
      alpha_bulk = norm_bulk / norm_bulk_old
      alpha_betav = norm_betav / norm_betav_old
      alpha_betah = norm_betah / norm_betah_old
      alpha_eta = norm_eta / norm_eta_old

      ! only if contribution is positive it will be considered, otherwise
      ! we set it to zero so that it becomes a steepest descent update
      if (alpha_bulk < 0.0) then
        alpha_bulk = 0.0
      endif
      if (alpha_betav < 0.0) then
        alpha_betav = 0.0
      endif
      if (alpha_betah < 0.0) then
        alpha_betah = 0.0
      endif
      if (alpha_eta < 0.0) then
        alpha_eta = 0.0
      endif

    else
      ! calculates only a single steplength applied to all
      alpha_all = (norm_bulk + norm_betav + norm_betah + norm_eta) &
                  / (norm_bulk_old + norm_betav_old + norm_betah_old + norm_eta_old)

      ! only if contribution is positive it will be considered, otherwise
      ! we set it to zero so that it becomes a steepest descent update
      if (alpha_all < 0.0) then
        alpha_all = 0.0
      endif

      ! sets each steplength to same single one
      alpha_bulk = alpha_all
      alpha_betav = alpha_all
      alpha_betah = alpha_all
      alpha_eta = alpha_all
    endif
    ! user output
    print *,'alpha gradient:'
    print *,'  bulk : ',alpha_bulk
    print *,'  betav: ',alpha_betav
    print *,'  betah: ',alpha_betah
    print *,'  eta  : ',alpha_eta
    print *
  endif
  ! broadcast values from rank 0 to all others
  call bcast_all_singlecr(alpha_bulk)
  call bcast_all_singlecr(alpha_betav)
  call bcast_all_singlecr(alpha_betah)
  call bcast_all_singlecr(alpha_eta)

  ! initializes kernel maximum
  depthmax(:) = 0._CUSTOM_REAL

  ! gradient in negative direction
  if (USE_OLD_GRADIENT) then
    ! uses old kernel/gradient updates ( phi_n-1)
    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

              ! note: uses old gradient update (phi_(n-1) as model_bulk_old, but
              !       given in negative gradient direction

              ! for bulk
              model_dbulk(i,j,k,ispec) = - kernel_bulk(i,j,k,ispec) &
                                         + alpha_bulk * model_dbulk_old(i,j,k,ispec)

              ! for shear
              model_dbetav(i,j,k,ispec) = - kernel_betav(i,j,k,ispec) &
                                          + alpha_betav * model_dbetav_old(i,j,k,ispec)

              model_dbetah(i,j,k,ispec) = - kernel_betah(i,j,k,ispec) &
                                          + alpha_betah * model_dbetah_old(i,j,k,ispec)

              ! for eta
              model_deta(i,j,k,ispec) = - kernel_eta(i,j,k,ispec) &
                                        + alpha_eta * model_deta_old(i,j,k,ispec)

              ! determines maximum kernel betav value within given radius
              if (USE_DEPTH_RANGE_MAXIMUM) then
                ! get radius of point
                iglob = ibool(i,j,k,ispec)
                r = z(iglob)

                ! stores maximum kernel betav/betah value in this depth slice,
                ! since betav/betah are most likely dominating
                if (r < R_TOP .and. r > R_BOTTOM) then
                  ! kernel betav value
                  max_vsv = abs( model_dbetav(i,j,k,ispec) )
                  if (depthmax(1) < max_vsv) then
                    depthmax(1) = max_vsv
                    depthmax_radius(1) = r
                  endif
                  ! kernel betav value
                  max_vsh = abs( model_dbetah(i,j,k,ispec) )
                  if (depthmax(2) < max_vsh) then
                    depthmax(2) = max_vsh
                    depthmax_radius(2) = r
                  endif
                endif
              endif

          enddo
        enddo
      enddo
    enddo
  else
    ! uses only old kernel/gradient
    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

              ! note: uses old kernels (lambda_(n-1) ) in negative gradient direction

              ! for bulk
              model_dbulk(i,j,k,ispec) = - kernel_bulk(i,j,k,ispec) &
                                         - alpha_bulk * kernel_bulk_old(i,j,k,ispec)

              ! for shear
              model_dbetav(i,j,k,ispec) = - kernel_betav(i,j,k,ispec) &
                                          - alpha_betav * kernel_betav_old(i,j,k,ispec)

              model_dbetah(i,j,k,ispec) = - kernel_betah(i,j,k,ispec) &
                                          - alpha_betah * kernel_betah_old(i,j,k,ispec)

              ! for eta
              model_deta(i,j,k,ispec) = - kernel_eta(i,j,k,ispec) &
                                        - alpha_eta * kernel_eta_old(i,j,k,ispec)


              ! determines maximum kernel betav value within given radius
              if (USE_DEPTH_RANGE_MAXIMUM) then
                ! get radius of point
                iglob = ibool(i,j,k,ispec)
                r = z(iglob)

                ! stores maximum kernel betav/betah value in this depth slice,
                ! since betav/betah are most likely dominating
                if (r < R_TOP .and. r > R_BOTTOM) then
                  ! kernel betav value
                  max_vsv = abs( model_dbetav(i,j,k,ispec) )
                  if (depthmax(1) < max_vsv) then
                    depthmax(1) = max_vsv
                    depthmax_radius(1) = r
                  endif
                  ! kernel betav value
                  max_vsh = abs( model_dbetah(i,j,k,ispec) )
                  if (depthmax(2) < max_vsh) then
                    depthmax(2) = max_vsh
                    depthmax_radius(2) = r
                  endif
                endif
              endif

          enddo
        enddo
      enddo
    enddo
  endif

  ! stores model_dbulk, ... arrays
  ! note: stores these new gradient before we scale them with the step length
  call write_gradient_tiso()

  ! statistics
  call min_all_cr(minval(model_dbulk),min_bulk)
  call max_all_cr(maxval(model_dbulk),max_bulk)

  call min_all_cr(minval(model_dbetav),min_vsv)
  call max_all_cr(maxval(model_dbetav),max_vsv)

  call min_all_cr(minval(model_dbetah),min_vsh)
  call max_all_cr(maxval(model_dbetah),max_vsh)

  call min_all_cr(minval(model_deta),min_eta)
  call max_all_cr(maxval(model_deta),max_eta)

  if (myrank == 0) then
    print *,'initial gradient updates:'
    print *,'  bulk min/max : ',min_bulk,max_bulk
    print *,'  betav min/max: ',min_vsv,max_vsv
    print *,'  betah min/max: ',min_vsh,max_vsh
    print *,'  eta min/max  : ',min_eta,max_eta
    print *
  endif

  ! determines maximum kernel betav value within given radius
  if (USE_DEPTH_RANGE_MAXIMUM) then
    ! maximum of all processes stored in max_vsv
    call max_all_cr(depthmax(1),max_vsv)
    call max_all_cr(depthmax(2),max_vsh)
    call max_all_cr(depthmax_radius(1),rmax_vsv)
    call max_all_cr(depthmax_radius(2),rmax_vsh)
  endif

  ! determines step length
  ! based on maximum gradient value (either vsv or vsh)
  if (myrank == 0) then

    ! determines maximum kernel betav value within given radius
    if (USE_DEPTH_RANGE_MAXIMUM) then
      depthmax(1) = max_vsv
      depthmax(2) = max_vsh
      depthmax_radius(1) = rmax_vsv
      depthmax_radius(2) = rmax_vsh

      max = maxval(depthmax)
      maxindex = maxloc(depthmax)
      depthmax_depth = depthmax_radius(maxindex(1))
      ! maximum in given depth range
      print *,'  using depth maximum: '
      print *,'  between depths (top/bottom)   : ',R_TOP,R_BOTTOM
      print *,'  maximum kernel value          : ',max
      print *,'  depth of maximum kernel value : ',depthmax_depth
      print *
    else
      ! maximum gradient values
      minmax(1) = abs(min_vsv)
      minmax(2) = abs(max_vsv)
      minmax(3) = abs(min_vsh)
      minmax(4) = abs(max_vsh)

      ! maximum value of all kernel maxima
      max = maxval(minmax)
    endif
    print *,'step length:'
    print *,'  using kernel maximum: ',max

    ! checks maximum value
    if (max < 1.e-25) stop 'Error maximum kernel value too small for update'

    ! chooses step length such that it becomes the desired, given step factor as inputted
    step_length = step_fac/max

    print *,'  step length : ',step_length
    print *

  endif
  call bcast_all_singlecr(step_length)


  ! gradient length sqrt( v^T * v )
  norm_bulk = sum( model_dbulk * model_dbulk )
  norm_betav = sum( model_dbetav * model_dbetav )
  norm_betah = sum( model_dbetah * model_dbetah )
  norm_eta = sum( model_deta * model_deta )

  call sum_all_cr(norm_bulk,norm_bulk_sum)
  call sum_all_cr(norm_betav,norm_betav_sum)
  call sum_all_cr(norm_betah,norm_betah_sum)
  call sum_all_cr(norm_eta,norm_eta_sum)

  if (myrank == 0) then
    norm_bulk = sqrt(norm_bulk_sum)
    norm_betav = sqrt(norm_betav_sum)
    norm_betah = sqrt(norm_betah_sum)
    norm_eta = sqrt(norm_eta_sum)

    print *,'norm model updates:'
    print *,'  bulk : ',norm_bulk
    print *,'  betav: ',norm_betav
    print *,'  betah: ',norm_betah
    print *,'  eta  : ',norm_eta
    print *
  endif

  ! multiply model updates by a subjective factor that will change the step
  model_dbulk(:,:,:,:) = step_length * model_dbulk(:,:,:,:)
  model_dbetav(:,:,:,:) = step_length * model_dbetav(:,:,:,:)
  model_dbetah(:,:,:,:) = step_length * model_dbetah(:,:,:,:)
  model_deta(:,:,:,:) = step_length * model_deta(:,:,:,:)

  ! statistics
  call min_all_cr(minval(model_dbulk),min_bulk)
  call max_all_cr(maxval(model_dbulk),max_bulk)

  call min_all_cr(minval(model_dbetav),min_vsv)
  call max_all_cr(maxval(model_dbetav),max_vsv)

  call min_all_cr(minval(model_dbetah),min_vsh)
  call max_all_cr(maxval(model_dbetah),max_vsh)

  call min_all_cr(minval(model_deta),min_eta)
  call max_all_cr(maxval(model_deta),max_eta)

  if (myrank == 0) then
    print *,'scaled gradient:'
    print *,'  bulk min/max : ',min_bulk,max_bulk
    print *,'  betav min/max: ',min_vsv,max_vsv
    print *,'  betah min/max: ',min_vsh,max_vsh
    print *,'  eta min/max  : ',min_eta,max_eta
    print *
  endif
  call synchronize_all()

  end subroutine get_cg_direction_tiso

