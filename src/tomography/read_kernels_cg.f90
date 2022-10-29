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


subroutine read_kernels_cg_tiso_old()

! reads in smoothed kernels from former iteration in OUTPUT_SUM.old/ : bulk, betav, betah, eta

  use tomography_kernels_tiso_cg

  implicit none
  real(kind=CUSTOM_REAL) :: min_vsv,min_vsh,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk
  logical:: exist,exist_all,use_old_gradient_all
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'reading kernels...'

  ! allocate arrays for storing kernels and perturbations
  ! transversely isotropic arrays
  allocate(kernel_bulk_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1075')
  allocate(kernel_betav_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1076')
  allocate(kernel_betah_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1077')
  allocate(kernel_eta_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1078')
  if (ier /= 0) stop 'error allocating kernel arrays'

  ! initializes arrays
  kernel_bulk_old = 0.0_CUSTOM_REAL
  kernel_betav_old = 0.0_CUSTOM_REAL
  kernel_betah_old = 0.0_CUSTOM_REAL
  kernel_eta_old = 0.0_CUSTOM_REAL

  ! checks if files are available:
  fname = 'bulk_c_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if (.not. exist) then
    print *,'Error file does not exist: ',trim(m_file)
    call exit_mpi(myrank,'file not exist')
  endif
  ! makes sure all processes have same flag
  call any_all_l(exist,exist_all)
  if (.not. exist_all) then
    print *,'old kernels do not exist: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! bulk kernel
  fname = 'bulk_c_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_bulk_old(:,:,:,1:nspec)
  close(IIN)

  ! betav kernel
  fname = 'bulk_betav_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betav_old(:,:,:,1:nspec)
  close(IIN)

  ! betah kernel
  fname = 'bulk_betah_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betah_old(:,:,:,1:nspec)
  close(IIN)

  ! eta kernel
  fname = 'eta_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_eta_old(:,:,:,1:nspec)
  close(IIN)

  ! statistics
  call min_all_cr(minval(kernel_bulk_old),min_bulk)
  call max_all_cr(maxval(kernel_bulk_old),max_bulk)

  call min_all_cr(minval(kernel_betah_old),min_vsh)
  call max_all_cr(maxval(kernel_betah_old),max_vsh)

  call min_all_cr(minval(kernel_betav_old),min_vsv)
  call max_all_cr(maxval(kernel_betav_old),max_vsv)

  call min_all_cr(minval(kernel_eta_old),min_eta)
  call max_all_cr(maxval(kernel_eta_old),max_eta)

  if (myrank == 0) then
    print *
    print *,'old kernels:'
    print *,'  bulk min/max : ',min_bulk,max_bulk
    print *,'  betav min/max: ',min_vsv,max_vsv
    print *,'  betah min/max: ',min_vsh,max_vsh
    print *,'  eta min/max  : ',min_eta,max_eta
    print *
  endif

  ! statistics output
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_kernels_minmax',status='unknown')
    write(IOUT,*) '#min_bulk #max_bulk #min_vsv #max_vsv #min_vsh #max_vsh #min_eta #max_eta'
    write(IOUT,'(4e24.12)') min_bulk, max_bulk, min_vsv, max_vsv, min_vsh, max_vsh, min_eta, max_eta
    close(IOUT)
  endif

  ! reads in old gradient directions (phi_(n-1))
  USE_OLD_GRADIENT = .true.

  ! checks if files are available:
  fname = 'dbulk_c'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if (.not. exist) then
    print *,'old kernel updates do not exist: ',trim(m_file)
    USE_OLD_GRADIENT = .false.
  endif
  ! makes sure all processes have same flag
  call any_all_l(exist,exist_all)
  if (.not. exist_all) then
    if (myrank == 0) print *,'old kernel updates do not exist for all: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! makes sure all processes have same flag
  use_old_gradient_all = .false.
  call synchronize_all()

  call any_all_l(USE_OLD_GRADIENT,use_old_gradient_all)
  if (.not. use_old_gradient_all) then
    print *,'old kernel updates exists, not consistent for all: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! reads in old gradient
  if (USE_OLD_GRADIENT) then

    ! user output
    if (myrank == 0) print *,'reading old gradient...'

    ! allocate arrays for storing old gradient
    ! transversely isotropic arrays
    allocate(model_dbulk_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1079')
    allocate(model_dbetav_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1080')
    allocate(model_dbetah_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1081')
    allocate(model_deta_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1082')
    if (ier /= 0) stop 'error allocating gradient arrays'

    ! initializes arrays
    model_dbulk_old = 0.0_CUSTOM_REAL
    model_dbetav_old = 0.0_CUSTOM_REAL
    model_dbetah_old = 0.0_CUSTOM_REAL
    model_deta_old = 0.0_CUSTOM_REAL

    ! bulk kernel
    fname = 'dbulk_c'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbulk_old(:,:,:,1:nspec)
    close(IIN)

    ! betav kernel
    fname = 'dbetav'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbetav_old(:,:,:,1:nspec)
    close(IIN)

    ! betah kernel
    fname = 'dbetah'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbetah_old(:,:,:,1:nspec)
    close(IIN)

    ! eta kernel
    fname = 'deta'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_deta_old(:,:,:,1:nspec)
    close(IIN)

    ! statistics
    call min_all_cr(minval(model_dbulk_old),min_bulk)
    call max_all_cr(maxval(model_dbulk_old),max_bulk)

    call min_all_cr(minval(model_dbetah_old),min_vsh)
    call max_all_cr(maxval(model_dbetah_old),max_vsh)

    call min_all_cr(minval(model_dbetav_old),min_vsv)
    call max_all_cr(maxval(model_dbetav_old),max_vsv)

    call min_all_cr(minval(model_deta_old),min_eta)
    call max_all_cr(maxval(model_deta_old),max_eta)

    if (myrank == 0) then
      print *
      print *,'old kernel updates:'
      print *,'  bulk min/max : ',min_bulk,max_bulk
      print *,'  betav min/max: ',min_vsv,max_vsv
      print *,'  betah min/max: ',min_vsh,max_vsh
      print *,'  eta min/max  : ',min_eta,max_eta
      print *
    endif

    ! statistics output
    if (PRINT_STATISTICS_FILES .and. myrank == 0) then
      open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_kernel_updates_minmax',status='unknown')
      write(IOUT,*) '#min_bulk #max_bulk #min_vsv #max_vsv #min_vsh #max_vsh #min_eta #max_eta'
      write(IOUT,'(4e24.12)') min_bulk, max_bulk, min_vsv, max_vsv, min_vsh, max_vsh, min_eta, max_eta
      close(IOUT)
    endif

  endif ! USE_OLD_GRADIENT
  call synchronize_all()

end subroutine read_kernels_cg_tiso_old

