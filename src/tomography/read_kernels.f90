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


subroutine read_kernels_iso()

! reads in smoothed kernels: bulk, beta, rho

  use tomography_kernels_iso

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: min_vp,min_vs,max_vp,max_vs,min_rho,max_rho
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'reading kernels...'

  ! allocate arrays for storing kernels and perturbations
  ! isotropic arrays
  allocate(kernel_bulk(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 904')
  allocate(kernel_beta(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 905')
  allocate(kernel_rho(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 906')
  if (ier /= 0) stop 'error allocating kernel arrays'

  ! initializes arrays
  kernel_bulk = 0.0_CUSTOM_REAL
  kernel_beta = 0.0_CUSTOM_REAL
  kernel_rho = 0.0_CUSTOM_REAL

  ! reads in smoothed (& summed) event kernel
  if (USE_ALPHA_BETA_RHO) then
    ! reads in alpha kernel
    fname = 'alpha_kernel_smooth'
  else
    ! reads in bulk_c kernel
    fname = 'bulk_c_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_bulk(:,:,:,1:nspec)
  close(IIN)

  ! beta kernel
  if (USE_ALPHA_BETA_RHO) then
    ! reads in beta kernel
    fname = 'beta_kernel_smooth'
  else
    ! reads in bulk_beta kernel
    fname = 'bulk_beta_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_beta(:,:,:,1:nspec)
  close(IIN)

  ! rho kernel
  if (USE_RHO_SCALING) then
    ! uses scaling relation with shear perturbations
    kernel_rho(:,:,:,:) = RHO_SCALING * kernel_beta(:,:,:,:)
    if (myrank == 0) print *,'  rho kernel uses scaling with shear kernel: scaling value = ',RHO_SCALING
  else
    ! uses rho kernel
    fname = 'rho_kernel_smooth'
    write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) kernel_rho(:,:,:,1:nspec)
    close(IIN)
  endif

  ! statistics
  call min_all_cr(minval(kernel_bulk),min_vp)
  call max_all_cr(maxval(kernel_bulk),max_vp)

  call min_all_cr(minval(kernel_beta),min_vs)
  call max_all_cr(maxval(kernel_beta),max_vs)

  call min_all_cr(minval(kernel_rho),min_rho)
  call max_all_cr(maxval(kernel_rho),max_rho)

  if (myrank == 0) then
    print *
    print *,'initial kernels:'
    if (USE_ALPHA_BETA_RHO) then
      print *,'  alpha min/max    : ',min_vp,max_vp
      print *,'  beta min/max     : ',min_vs,max_vs
    else
      print *,'  bulk_c min/max   : ',min_vp,max_vp
      print *,'  bulk_beta min/max: ',min_vs,max_vs
    endif
    print *,'  rho min/max      : ',min_rho,max_rho
    print *
  endif
  call synchronize_all()

  ! statistics output
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_kernels_minmax',status='unknown')
    write(IOUT,*) '#min_vs #max_vs #min_vp #max_vp #min_rho #max_rho'
    write(IOUT,'(4e24.12)') min_vs, max_vs, min_vp, max_vp, min_rho, max_rho
    close(IOUT)
  endif

end subroutine read_kernels_iso

!
!-------------------------------------------------------------------------------------------------
!


subroutine read_kernels_tiso()

! reads in smoothed kernels: bulk, betav, betah, eta

  use tomography_kernels_tiso

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: min_vsv,min_vsh,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'reading kernels...'

  ! allocate arrays for storing kernels and perturbations
  ! transversely isotropic arrays
  allocate(kernel_bulk(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 907')
  allocate(kernel_betav(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 908')
  allocate(kernel_betah(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 909')
  allocate(kernel_eta(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 910')
  if (ier /= 0) stop 'error allocating kernel arrays'

  ! initializes arrays
  kernel_bulk = 0.0_CUSTOM_REAL
  kernel_betav = 0.0_CUSTOM_REAL
  kernel_betah = 0.0_CUSTOM_REAL
  kernel_eta = 0.0_CUSTOM_REAL

  ! bulk kernel
  fname = 'bulk_c_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_bulk(:,:,:,1:nspec)
  close(IIN)

  ! betav kernel
  fname = 'bulk_betav_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betav(:,:,:,1:nspec)
  close(IIN)

  ! betah kernel
  fname = 'bulk_betah_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betah(:,:,:,1:nspec)
  close(IIN)

  ! eta kernel
  fname = 'eta_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_eta(:,:,:,1:nspec)
  close(IIN)


  ! statistics
  call min_all_cr(minval(kernel_bulk),min_bulk)
  call max_all_cr(maxval(kernel_bulk),max_bulk)

  call min_all_cr(minval(kernel_betah),min_vsh)
  call max_all_cr(maxval(kernel_betah),max_vsh)

  call min_all_cr(minval(kernel_betav),min_vsv)
  call max_all_cr(maxval(kernel_betav),max_vsv)

  call min_all_cr(minval(kernel_eta),min_eta)
  call max_all_cr(maxval(kernel_eta),max_eta)

  if (myrank == 0) then
    print *
    print *,'initial kernels:'
    print *,'  bulk min/max : ',min_bulk,max_bulk
    print *,'  betav min/max: ',min_vsv,max_vsv
    print *,'  betah min/max: ',min_vsh,max_vsh
    print *,'  eta min/max  : ',min_eta,max_eta
    print *
  endif
  call synchronize_all()

  ! statistics output
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_kernels_minmax',status='unknown')
    write(IOUT,*) '#min_bulk #max_bulk #min_vsv #max_vsv #min_vsh #max_vsh #min_eta #max_eta'
    write(IOUT,'(4e24.12)') min_bulk, max_bulk, min_vsv, max_vsv, min_vsh, max_vsh, min_eta, max_eta
    close(IOUT)
  endif

end subroutine read_kernels_tiso

