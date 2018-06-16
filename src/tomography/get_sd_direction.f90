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


  subroutine get_sd_direction_iso()

! calculates gradient by steepest descent method

  use tomography_kernels_iso

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL):: r,max,depth_max
  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_beta,norm_rho
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_beta_sum,norm_rho_sum
  real(kind=CUSTOM_REAL) :: minmax(4)
  real(kind=CUSTOM_REAL) :: min_beta,min_rho,max_beta,max_rho,min_bulk,max_bulk
  integer :: iglob
  integer :: i,j,k,ispec,ier

  ! allocate arrays for storing gradient
  ! isotropic arrays
  allocate(model_dbulk(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1046')
  allocate(model_dbeta(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1047')
  allocate(model_drho(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1048')
  if (ier /= 0) stop 'error allocating gradient arrays'

  ! initializes arrays
  model_dbulk = 0.0_CUSTOM_REAL
  model_dbeta = 0.0_CUSTOM_REAL
  model_drho = 0.0_CUSTOM_REAL

  ! initializes kernel maximum
  max = 0._CUSTOM_REAL

  ! gradient in negative direction for steepest descent
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

            ! for bulk
            model_dbulk(i,j,k,ispec) = - kernel_bulk(i,j,k,ispec)

            ! for shear
            model_dbeta(i,j,k,ispec) = - kernel_beta(i,j,k,ispec)

            ! for rho
            model_drho(i,j,k,ispec) = - kernel_rho(i,j,k,ispec)

            ! determines maximum kernel beta value within given radius
            if (USE_DEPTH_RANGE_MAXIMUM) then
              ! get depth of point (assuming z in vertical direction, up in positive direction)
              iglob = ibool(i,j,k,ispec)
              r = z(iglob)

              ! stores maximum kernel betav value in this depth slice, since betav is most likely dominating
              if (r < R_TOP .and. r > R_BOTTOM) then
                ! shear kernel value
                max_beta = abs( kernel_beta(i,j,k,ispec) )
                if (max < max_beta) then
                  max = max_beta
                  depth_max = r
                endif
              endif
            endif

        enddo
      enddo
    enddo
  enddo

  ! stores model_dbulk, ... arrays
  call write_gradient_iso()

  ! statistics
  call min_all_cr(minval(model_dbulk),min_bulk)
  call max_all_cr(maxval(model_dbulk),max_bulk)

  call min_all_cr(minval(model_dbeta),min_beta)
  call max_all_cr(maxval(model_dbeta),max_beta)

  call min_all_cr(minval(model_drho),min_rho)
  call max_all_cr(maxval(model_drho),max_rho)

  if (myrank == 0) then
    print *,'initial gradient:'
    print *,'  a min/max   : ',min_bulk,max_bulk
    print *,'  beta min/max: ',min_beta,max_beta
    print *,'  rho min/max : ',min_rho,max_rho
    print *
  endif

  ! statistics output
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_gradient_minmax',status='unknown')
    write(IOUT,*) '#min_beta #max_beta #min_bulk #max_bulk #min_rho #max_rho'
    write(IOUT,'(6e24.12)') min_beta, max_beta, min_bulk, max_bulk, min_rho, max_rho
    close(IOUT)
  endif

  ! determines maximum kernel betav value within given radius
  if (USE_DEPTH_RANGE_MAXIMUM) then
    ! maximum of all processes stored in max_vsv
    call max_all_cr(max,max_beta)
    max = max_beta
  endif

  ! determines step length based on maximum gradient value (either shear or bulk)
  if (myrank == 0) then

    ! determines maximum kernel betav value within given radius
    if (USE_DEPTH_RANGE_MAXIMUM) then
      print *,'  using depth maximum: '
      print *,'  between depths (top/bottom)   : ',R_TOP,R_BOTTOM
      print *,'  maximum kernel value          : ',max
      print *,'  depth of maximum kernel value : ',depth_max
      print *
    else
      ! maximum gradient values
      minmax(1) = abs(min_beta)
      minmax(2) = abs(max_beta)
      minmax(3) = abs(min_bulk)
      minmax(4) = abs(max_bulk)

      ! maximum value of all kernel maxima
      max = maxval(minmax)
    endif
    print *,'step length:'
    print *,'  using kernel maximum: ',max

    ! checks maximum value
    if (max < 1.e-25) stop 'Error maximum kernel value too small for update'

    ! chooses step length such that it becomes the desired, given step factor as inputted
    step_length = step_fac/max

    print *,'  step length value   : ',step_length
    print *
  endif
  call bcast_all_singlecr(step_length)


  ! gradient length sqrt( v^T * v )
  norm_bulk = sum( model_dbulk * model_dbulk )
  norm_beta = sum( model_dbeta * model_dbeta )
  norm_rho = sum( model_drho * model_drho )

  call sum_all_cr(norm_bulk,norm_bulk_sum)
  call sum_all_cr(norm_beta,norm_beta_sum)
  call sum_all_cr(norm_rho,norm_rho_sum)

  if (myrank == 0) then
    norm_bulk = sqrt(norm_bulk_sum)
    norm_beta = sqrt(norm_beta_sum)
    norm_rho = sqrt(norm_rho_sum)

    print *,'norm model updates:'
    print *,'  a   : ',norm_bulk
    print *,'  beta: ',norm_beta
    print *,'  rho : ',norm_rho
    print *
  endif

  ! statistics output
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_sum',status='unknown')
    write(IOUT,*) '#norm_beta #norm_bulk #norm_rho'
    write(IOUT,'(3e24.12)') norm_beta, norm_bulk, norm_rho
    close(IOUT)
  endif

  ! multiply model updates by a subjective factor that will change the step
  model_dbulk(:,:,:,:) = step_length * model_dbulk(:,:,:,:)
  model_dbeta(:,:,:,:) = step_length * model_dbeta(:,:,:,:)
  model_drho(:,:,:,:) = step_length * model_drho(:,:,:,:)

  ! statistics
  call min_all_cr(minval(model_dbulk),min_bulk)
  call max_all_cr(maxval(model_dbulk),max_bulk)

  call min_all_cr(minval(model_dbeta),min_beta)
  call max_all_cr(maxval(model_dbeta),max_beta)

  call min_all_cr(minval(model_drho),min_rho)
  call max_all_cr(maxval(model_drho),max_rho)

  if (myrank == 0) then
    print *,'scaled gradient:'
    print *,'  a min/max   : ',min_bulk,max_bulk
    print *,'  beta min/max: ',min_beta,max_beta
    print *,'  rho min/max : ',min_rho,max_rho
    print *
  endif
  call synchronize_all()

  ! statistics output
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_scaled_gradient',status='unknown')
    write(IOUT,*) '#min_beta #max_beta #min_bulk #max_bulk #min_rho #max_rho'
    write(IOUT,'(6e24.12)') min_beta,max_beta,min_bulk,max_bulk,min_rho,max_rho
    close(IOUT)
  endif

  end subroutine get_sd_direction_iso

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_sd_direction_tiso()

! calculates gradient by steepest descent method

  use tomography_kernels_tiso

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL):: r,max,depth_max
  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_betav,norm_betah,norm_eta
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_betav_sum, &
    norm_betah_sum,norm_eta_sum
  real(kind=CUSTOM_REAL) :: min_vsv,min_vsh,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk
  real(kind=CUSTOM_REAL) :: minmax(4)
  integer :: iglob
  integer :: i,j,k,ispec,ier

  ! allocate arrays for storing gradient
  ! transversely isotropic arrays
  allocate(model_dbulk(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1049')
  allocate(model_dbetav(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1050')
  allocate(model_dbetah(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1051')
  allocate(model_deta(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1052')
  if (ier /= 0) stop 'error allocating gradient arrays'

  ! initializes arrays
  model_dbulk = 0.0_CUSTOM_REAL
  model_dbetav = 0.0_CUSTOM_REAL
  model_dbetah = 0.0_CUSTOM_REAL
  model_deta = 0.0_CUSTOM_REAL

  ! initializes kernel maximum
  max = 0._CUSTOM_REAL

  ! gradient in negative direction for steepest descent
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

            ! for bulk
            model_dbulk(i,j,k,ispec) = - kernel_bulk(i,j,k,ispec)

            ! for shear
            model_dbetav(i,j,k,ispec) = - kernel_betav(i,j,k,ispec)
            model_dbetah(i,j,k,ispec) = - kernel_betah(i,j,k,ispec)

            ! for eta
            model_deta(i,j,k,ispec) = - kernel_eta(i,j,k,ispec)

            ! determines maximum kernel betav value within given radius
            if (USE_DEPTH_RANGE_MAXIMUM) then
              ! get depth of point (assuming z in vertical direction, up in positive direction)
              iglob = ibool(i,j,k,ispec)
              r = z(iglob)

              ! stores maximum kernel betav value in this depth slice, since betav is most likely dominating
              if (r < R_TOP .and. r > R_BOTTOM) then
                ! kernel betav value
                max_vsv = abs( kernel_betav(i,j,k,ispec) )
                if (max < max_vsv) then
                  max = max_vsv
                  depth_max = r
                endif
              endif
            endif

        enddo
      enddo
    enddo
  enddo

  ! stores model_dbulk, ... arrays
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
    print *,'initial gradient:'
    print *,'  bulk min/max : ',min_bulk,max_bulk
    print *,'  betav min/max: ',min_vsv,max_vsv
    print *,'  betah min/max: ',min_vsh,max_vsh
    print *,'  eta min/max  : ',min_eta,max_eta
    print *
  endif

  ! determines maximum kernel betav value within given radius
  if (USE_DEPTH_RANGE_MAXIMUM) then
    ! maximum of all processes stored in max_vsv
    call max_all_cr(max,max_vsv)
    max = max_vsv
    depth_max = 6371.0 *( 1.0 - depth_max )
  endif

  ! determines step length
  ! based on maximum gradient value (either vsv or vsh)
  if (myrank == 0) then

    ! determines maximum kernel betav value within given radius
    if (USE_DEPTH_RANGE_MAXIMUM) then
      print *,'  using depth maximum: '
      print *,'  between depths (top/bottom)   : ',R_TOP,R_BOTTOM
      print *,'  maximum kernel value          : ',max
      print *,'  depth of maximum kernel value : ',depth_max
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

    print *,'  step length value   : ',step_length
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

  end subroutine get_sd_direction_tiso

