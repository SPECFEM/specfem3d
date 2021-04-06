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

program model_update

  use specfem_par, only: kappastore,mustore
  use specfem_par, only: NPROC,OUTPUT_FILES,LOCAL_PATH
  use specfem_par_elastic, only: rho_vp,rho_vs

  use tomography_model_iso
  use tomography_kernels_iso

  implicit none

  ! ======================================================
  ! USER PARAMETERS

  ! directory where the summed and smoothed input kernels are linked
  character(len=MAX_STRING_LEN) :: INPUT_KERNELS_DIR_NAME = 'sum_smooth_kern/'

  ! directory where the mesh files for the NEW model will be written
  character(len=MAX_STRING_LEN) :: LOCAL_PATH_NEW_NAME = 'mesh_files_m01/'

  ! directory where the output files of model_update will be written
  character(len=MAX_STRING_LEN) :: OUTPUT_STATISTICS_DIR_NAME = 'OUTPUT_FILES_MODEL_UPD/'

  ! threshold the old model ("current model")
  logical, parameter :: MINMAX_THRESHOLD_OLD = .false.
  ! threshold the new model
  logical, parameter :: MINMAX_THRESHOLD_NEW = .false.
  ! threshold values
  ! minmax wavespeed values for southern california simulations
  real(kind=CUSTOM_REAL),parameter :: VS_MIN =  600.0
  real(kind=CUSTOM_REAL),parameter :: VS_MAX = 4700.0
  real(kind=CUSTOM_REAL),parameter :: VP_MIN = 1500.0
  real(kind=CUSTOM_REAL),parameter :: VP_MAX = 8200.0
  ! thresholding density
  logical, parameter :: THRESHOLD_RHO = .false.
  real(kind=CUSTOM_REAL),parameter :: RHO_MIN = 1025.0 ! salt water
  real(kind=CUSTOM_REAL),parameter :: RHO_MAX = 5570.0 ! PREM maximum mantle density

  ! ======================================================

  character(len=MAX_STRING_LEN) :: fname
  character(len=MAX_STRING_LEN*2) :: m_file

  integer :: ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: total_model

  ! statistics:
  !models
  real(kind=CUSTOM_REAL) :: vsmin_before, vsmax_before, vpmin_before, vpmax_before
  real(kind=CUSTOM_REAL) :: vsmin_after, vsmax_after, vpmin_after, vpmax_after
  real(kind=CUSTOM_REAL) :: vsmin_new_before, vsmax_new_before, vpmin_new_before, vpmax_new_before
  real(kind=CUSTOM_REAL) :: vsmin_new_after, vsmax_new_after, vpmin_new_after, vpmax_new_after
  real(kind=CUSTOM_REAL) :: rhomin_before, rhomax_before, rhomin_after, rhomax_after
  real(kind=CUSTOM_REAL) :: rhomin_new_before, rhomax_new_before, rhomin_new_after, rhomax_new_after

  integer :: i,j,k,ispec

  ! ============ program starts here =====================

  ! initializes arrays
  call initialize()

  ! default input/output directories
  ! directory where the summed and smoothed input kernels are linked
  if (OUTPUT_FILES(len_trim(OUTPUT_FILES):len_trim(OUTPUT_FILES)) /= '/') then
    INPUT_KERNELS_DIR = trim(OUTPUT_FILES)//'/'//trim(INPUT_KERNELS_DIR_NAME)
  else
    INPUT_KERNELS_DIR = trim(OUTPUT_FILES)//trim(INPUT_KERNELS_DIR_NAME)
  endif

  ! directory where the mesh files for the NEW model will be written
  if (OUTPUT_FILES(len_trim(OUTPUT_FILES):len_trim(OUTPUT_FILES)) /= '/') then
    OUTPUT_MODEL_DIR = trim(OUTPUT_FILES)//'/'//trim(LOCAL_PATH_NEW_NAME)
  else
    OUTPUT_MODEL_DIR = trim(OUTPUT_FILES)//trim(LOCAL_PATH_NEW_NAME)
  endif

  ! directory where the output files of model_update will be written
  PRINT_STATISTICS_FILES = .true.

  if (OUTPUT_FILES(len_trim(OUTPUT_FILES):len_trim(OUTPUT_FILES)) /= '/') then
    OUTPUT_STATISTICS_DIR = trim(OUTPUT_FILES)//'/'//trim(OUTPUT_STATISTICS_DIR_NAME)
  else
    OUTPUT_STATISTICS_DIR = trim(OUTPUT_FILES)//trim(OUTPUT_STATISTICS_DIR_NAME)
  endif

  ! reads in parameters needed
  call read_parameters_tomo()

  ! user output
  if (myrank == 0) then
    print *
    print *,'***********'
    print *,'program model_update: '
    print *,'  NPROC: ',NPROC
    print *,'  NSPEC: ', NSPEC
    print *,'  NGLOB: ', NGLOB
    print *
    print *,'model update for vs & vp & rho'
    print *,'  step_fac = ',step_fac
    print *
    if (USE_ALPHA_BETA_RHO) then
      print *,'kernel parameterization: (alpha,beta,rho)'
    else
      print *,'kernel parameterization: (bulk,beta,rho)'
    endif
    print *
    if (USE_RHO_SCALING) then
      print *,'scaling rho perturbations'
      print *
    endif
    if (MINMAX_THRESHOLD_OLD) print *,'thresholding current (old) model wavespeed values'
    if (MINMAX_THRESHOLD_NEW) print *,'thresholding new (updated) model wavespeed values'
    if (MINMAX_THRESHOLD_OLD .or. MINMAX_THRESHOLD_NEW) then
      if (THRESHOLD_RHO) then
        print *,'  thresholds for: (vs,vp,rho)'
      else
        print *,'  thresholds for: (vs,vp)'
      endif
    endif
    print *,'***********'
    print *
  endif
  call synchronize_all()

  ! reads in external mesh
  call get_external_mesh()

  !===================================================
  ! MODEL UPDATE
  ! allocation
  ! model and kernel variables
  allocate(model_vp(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 911')
  allocate(model_vs(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 912')
  allocate(model_rho(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 913')

  allocate(total_model(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 914')

  ! initialize variables
  ! old model
  model_vp = 0._CUSTOM_REAL
  model_vs = 0._CUSTOM_REAL
  model_rho = 0._CUSTOM_REAL

  !! set up thresholds for old and new models

  if (MINMAX_THRESHOLD_OLD .or. MINMAX_THRESHOLD_NEW) then
    if (myrank == 0) then
      print *,'threshold wavespeed values:'
      print *,'  VS_MIN, VS_MAX  : ',VS_MIN, VS_MAX
      print *,'  VP_MIN, VP_MAX  : ',VP_MIN, VP_MAX
      if (THRESHOLD_RHO) print *,'  RHO_MIN, RHO_MAX: ',RHO_MIN, RHO_MAX
      print *
    endif
    ! statistics output
    if (PRINT_STATISTICS_FILES .and. myrank == 0) then
      if (THRESHOLD_RHO) then
        open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_VS_VP_RHO_MINMAX',status='unknown')
        write(IOUT,*) '#VS_MIN #VS_MAX #VP_MIN #VP_MAX #RHO_MIN #RHO_MAX'
        write(IOUT,'(6e24.12)') VS_MIN, VS_MAX, VP_MIN, VP_MAX, RHO_MIN, RHO_MAX
        close(IOUT)
      else
        open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_VS_VP_MINMAX',status='unknown')
        write(IOUT,*) '#VS_MIN #VS_MAX #VP_MIN #VP_MAX'
        write(IOUT,'(4e24.12)') VS_MIN, VS_MAX, VP_MIN, VP_MAX
        close(IOUT)
      endif
    endif
  endif


  !---------------------------------------------------------------------------------------------
  ! calculate vp,vs,rho values from kappastore,mustore,rho_vp,rho_vs in OLD external_mesh.bin
  !---------------------------------------------------------------------------------------------

  ! vp model
  where( rho_vp /= 0._CUSTOM_REAL ) model_vp = (FOUR_THIRDS * mustore + kappastore) / rho_vp

  ! vs model
  where( rho_vs /= 0._CUSTOM_REAL ) model_vs = mustore / rho_vs

  ! rho model
  where( rho_vp /= 0._CUSTOM_REAL ) model_rho = rho_vp**2 / (FOUR_THIRDS * mustore + kappastore)


  ! compute minmax values of current model
  ! NOTE: mpi_reduce operates on the values from all procs,
  !       but the reduced value only exists on the root proc.
  call min_all_cr(minval(model_vs(:,:,:,1:nspec)), vsmin_before)
  call max_all_cr(maxval(model_vs(:,:,:,1:nspec)), vsmax_before)
  call min_all_cr(minval(model_vp(:,:,:,1:nspec)), vpmin_before)
  call max_all_cr(maxval(model_vp(:,:,:,1:nspec)), vpmax_before)
  call min_all_cr(minval(model_rho(:,:,:,1:nspec)), rhomin_before)
  call max_all_cr(maxval(model_rho(:,:,:,1:nspec)), rhomax_before)

  if (myrank == 0) then
    if (MINMAX_THRESHOLD_OLD) then
      print *,'current model values before thresholding:'
    else
      print *,'current model values:'
    endif
    print *,'  vs min/max : ',vsmin_before, vsmax_before
    print *,'  vp min/max : ',vpmin_before, vpmax_before
    print *,'  rho min/max: ',rhomin_before, rhomax_before
    print *
  endif

  ! statistics output
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    if (MINMAX_THRESHOLD_OLD) then
      m_file = trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_minmax_before'
    else
      m_file = trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_minmax'
    endif
    open(IOUT,file=trim(m_file),status='unknown')
    write(IOUT,*) '#vsmin_before #vsmax_before #vpmin_before #vpmax_before #rhomin_before #rhomax_before'
    write(IOUT,'(6e24.12)') vsmin_before, vsmax_before, vpmin_before, vpmax_before, rhomin_before, rhomax_before
    close(IOUT)
  endif

  ! threshold current model and write out the modified version
  if (MINMAX_THRESHOLD_OLD) then
    do ispec=1,NSPEC
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            ! vs
            if (model_vs(i,j,k,ispec) < VS_MIN) model_vs(i,j,k,ispec) = VS_MIN
            if (model_vs(i,j,k,ispec) > VS_MAX) model_vs(i,j,k,ispec) = VS_MAX
            ! vp
            if (model_vp(i,j,k,ispec) < VP_MIN) model_vp(i,j,k,ispec) = VP_MIN
            if (model_vp(i,j,k,ispec) > VP_MAX) model_vp(i,j,k,ispec) = VP_MAX
            ! rho
            if (THRESHOLD_RHO) then
              if (model_rho(i,j,k,ispec) < RHO_MIN) model_rho(i,j,k,ispec) = RHO_MIN
              if (model_rho(i,j,k,ispec) > RHO_MAX) model_rho(i,j,k,ispec) = RHO_MAX
            endif
          enddo
        enddo
      enddo
    enddo

    ! compute minmax values of the thresholded current model
    call min_all_cr(minval(model_vs(:,:,:,1:nspec)), vsmin_after)
    call max_all_cr(maxval(model_vs(:,:,:,1:nspec)), vsmax_after)
    call min_all_cr(minval(model_vp(:,:,:,1:nspec)), vpmin_after)
    call max_all_cr(maxval(model_vp(:,:,:,1:nspec)), vpmax_after)
    call min_all_cr(minval(model_rho(:,:,:,1:nspec)), rhomin_after)
    call max_all_cr(maxval(model_rho(:,:,:,1:nspec)), rhomax_after)

    if (myrank == 0) then
      print *,'current model values after thresholding:'
      print *,'  vs min/max: ',vsmin_after, vsmax_after
      print *,'  vp min/max: ',vpmin_after, vpmax_after
      print *,'  rho min/max: ',rhomin_after, rhomax_after
      print *
    endif

    ! statistics output
    if (PRINT_STATISTICS_FILES .and. myrank == 0) then
      open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_minmax_after',status='unknown')
      write(IOUT,*) '#vsmin_after #vsmax_after #vpmin_after #vpmax_after #rhomin_after rhomax_after'
      write(IOUT,'(6e24.12)') vsmin_after, vsmax_after, vpmin_after, vpmax_after, rhomin_after, rhomax_after
      close(IOUT)
    endif

    ! file output
    fname = 'vs'
    write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    open(IOUT,file=trim(m_file),form='unformatted')
    write(IOUT) model_vs(:,:,:,1:nspec)
    close(IOUT)

    fname = 'vp'
    write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    open(IOUT,file=trim(m_file),form='unformatted')
    write(IOUT) model_vp(:,:,:,1:nspec)
    close(IOUT)

    if (THRESHOLD_RHO) then
      fname = 'rho'
      write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
      open(IOUT,file=trim(m_file),form='unformatted')
      write(IOUT) model_rho(:,:,:,1:nspec)
      close(IOUT)
    endif
  endif


  !---------------------------------------------------------------------------------------------
  ! reads in smoothed (& summed) event kernel
  !---------------------------------------------------------------------------------------------
  call read_kernels_iso()

  !---------------------------------------------------------------------------------------------
  ! computes search direction
  !---------------------------------------------------------------------------------------------
  call get_sd_direction_iso()

  !---------------------------------------------------------------------------------------------
  ! new model
  !---------------------------------------------------------------------------------------------

  ! computes new model values for alpha, beta and rho
  ! allocate new model arrays
  allocate(model_vp_new(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 915')
  allocate(model_vs_new(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 916')
  allocate(model_rho_new(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 917')
  if (ier /= 0) stop 'Error allocating model arrays'

  ! initializes arrays
  model_vp_new = 0.0_CUSTOM_REAL
  model_vs_new = 0.0_CUSTOM_REAL
  model_rho_new = 0.0_CUSTOM_REAL

  ! S wavespeed model
  model_vs_new = model_vs * exp( model_dbeta )

  ! P wavespeed model
  if (USE_ALPHA_BETA_RHO) then
    ! new vp values use alpha model update
    model_vp_new = model_vp * exp( model_dbulk )
  else
    ! new vp values use bulk model update:
    ! this is based on vp_new = sqrt( bulk_new**2 + 4/3 vs_new**2)
    model_vp_new = sqrt( model_vp**2 * exp(2.0*model_dbulk) + FOUR_THIRDS * model_vs**2 *( &
                            exp(2.0*model_dbeta) - exp(2.0*model_dbulk) ) )
  endif

  ! Rho density model
  model_rho_new = model_rho * exp( model_drho )

  ! statistics
  ! compute minmax values of new model
  call min_all_cr(minval(model_vs_new(:,:,:,1:nspec)), vsmin_new_before)
  call max_all_cr(maxval(model_vs_new(:,:,:,1:nspec)), vsmax_new_before)
  call min_all_cr(minval(model_vp_new(:,:,:,1:nspec)), vpmin_new_before)
  call max_all_cr(maxval(model_vp_new(:,:,:,1:nspec)), vpmax_new_before)
  call min_all_cr(minval(model_rho_new(:,:,:,1:nspec)), rhomin_new_before)
  call max_all_cr(maxval(model_rho_new(:,:,:,1:nspec)), rhomax_new_before)

  if (myrank == 0) then
    if (MINMAX_THRESHOLD_NEW) then
      print *,'new model values before thresholding:'
    else
      print *,'new model values:'
    endif
    print *,'  vs min/max : ',vsmin_new_before, vsmax_new_before
    print *,'  vp min/max : ',vpmin_new_before, vpmax_new_before
    print *,'  rho min/max: ',rhomin_new_before, rhomax_new_before
    print *
  endif

  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    if (MINMAX_THRESHOLD_NEW) then
      m_file = trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_new_minmax_before'
    else
      m_file = trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_new_minmax'
    endif
    open(IOUT,file=trim(m_file),status='unknown')
    write(IOUT,*) '#vsmin_new_before #vsmax_new_before #vpmin_new_before ' &
               // '#vpmax_new_before #rhomin_new_before #rhomax_new_before'
    write(IOUT,'(6e24.12)') vsmin_new_before, vsmax_new_before, &
                            vpmin_new_before, vpmax_new_before,rhomin_new_before, rhomax_new_before
    close(IOUT)
  endif

  ! threshold model according to minmax values specified above
  if (MINMAX_THRESHOLD_NEW) then
    do ispec=1,NSPEC
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            ! vs
            if (model_vs_new(i,j,k,ispec) < VS_MIN) model_vs_new(i,j,k,ispec) = VS_MIN
            if (model_vs_new(i,j,k,ispec) > VS_MAX) model_vs_new(i,j,k,ispec) = VS_MAX
            ! vp
            if (model_vp_new(i,j,k,ispec) < VP_MIN) model_vp_new(i,j,k,ispec) = VP_MIN
            if (model_vp_new(i,j,k,ispec) > VP_MAX) model_vp_new(i,j,k,ispec) = VP_MAX
            ! rho
            if (THRESHOLD_RHO) then
              if (model_rho_new(i,j,k,ispec) < RHO_MIN) model_rho_new(i,j,k,ispec) = RHO_MIN
              if (model_rho_new(i,j,k,ispec) > RHO_MAX) model_rho_new(i,j,k,ispec) = RHO_MAX
            endif
          enddo
        enddo
      enddo
    enddo

    ! write out new models and their global minmax values
    ! compute minmax values of new model
    call min_all_cr(minval(model_vs_new(:,:,:,1:nspec)), vsmin_new_after)
    call max_all_cr(maxval(model_vs_new(:,:,:,1:nspec)), vsmax_new_after)
    call min_all_cr(minval(model_vp_new(:,:,:,1:nspec)), vpmin_new_after)
    call max_all_cr(maxval(model_vp_new(:,:,:,1:nspec)), vpmax_new_after)
    call min_all_cr(minval(model_rho_new(:,:,:,1:nspec)), rhomin_new_after)
    call max_all_cr(maxval(model_rho_new(:,:,:,1:nspec)), rhomax_new_after)

    if (myrank == 0) then
      print *,'new model values after thresholding:'
      print *,'  vs min/max : ',vsmin_new_after, vsmax_new_after
      print *,'  vp min/max : ',vpmin_new_after, vpmax_new_after
      print *,'  rho min/max: ',rhomin_new_after, rhomax_new_after
      print *
    endif

    ! this should only be different if using MINMAX_THRESHOLD_NEW
    if (PRINT_STATISTICS_FILES .and. myrank == 0) then
      open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_new_minmax_after',status='unknown')
      write(IOUT,*) '#vsmin_new_after #vsmax_new_after #vpmin_new_after ' &
                 // '#vpmax_new_after #rhomin_new_after #rhomax_new_after'
      write(IOUT,'(6e24.12)') vsmin_new_after, vsmax_new_after, &
                              vpmin_new_after, vpmax_new_after, rhomin_new_after, rhomax_new_after
      close(IOUT)
    endif
  endif
  call synchronize_all()

  !---------------------------------------------------------------------------------------------
  ! write NEW external_mesh.bin
  ! and store NEW model files vp_new.bin vs_new.bin rho_new.bin vp_new.vtk vs_new.vtk rho_new.vtk
  ! calling save_external_bin_m_up with SAVE_MESH_FILES=true
  ! and also write NEW attenuation files attenuation.bin and attenuation.vtk (this should be equal to the old one)
  !---------------------------------------------------------------------------------------------
  call save_new_databases()


  !---------------------------------------------------------------------------------------------
  ! ALTERNATIVELY:
  ! SAVE_MESH_FILES=false
  ! directly stores new model files vp.bin vs.bin rho.bin vp.vtk vs.vtk rho.vtk
  ! NOTE: save_external_bin_m_up will create NEW external_mesh.bin anyway
  !---------------------------------------------------------------------------------------------

  ! call write_new_model_iso()

! VTK file output
!   ! vs model
!   fname = 'vs_new'
!   write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'/proc',myrank,trim(REG)//trim(fname)
!   call write_VTK_data_gll_cr(nspec,nglob, &
!                       xstore,ystore,zstore,ibool, &
!                       model_vs_new,m_file)
!
!
!   ! vp model
!   fname = 'vp_new'
!   write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'/proc',myrank,trim(REG)//trim(fname)
!   call write_VTK_data_gll_cr(nspec,nglob, &
!                       xstore,ystore,zstore,ibool, &
!                       model_vp_new,m_file)
!
!
!   ! rho model
!   !remove the comment if we are interested in rho
!   fname = 'rho_new'
!   write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'/proc',myrank,trim(REG)//trim(fname)
!   call write_VTK_data_gll_cr(nspec,nglob, &
!                       xstore,ystore,zstore,ibool, &
!                       model_rho_new,m_file)


  !---------------------------------------------------------------------------------------------
  ! P & S & rho model update
  ! stores the linear approximations of the model updates
  !---------------------------------------------------------------------------------------------
  call write_new_model_perturbations_iso()


  !---------------------------------------------------------------------------------------------
  ! compute Poisson's ratio of the current model and new model
  !---------------------------------------------------------------------------------------------

  ! Poisson's ratio of current model
  total_model = 0.
  total_model = 0.5 * ( model_vp**2 - 2.0*model_vs**2) / ( model_vp**2 - model_vs**2)
  fname = 'poisson'
  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  open(IOUT,file=trim(m_file),form='unformatted')
  write(IOUT) total_model(:,:,:,1:nspec)
  close(IOUT)

  ! Poisson's ratio of new model
  total_model = 0.
  total_model = 0.5 * ( model_vp_new**2 - 2.0*model_vs_new**2) / ( model_vp_new**2 - model_vs_new**2)
  fname = 'poisson_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  open(IOUT,file=trim(m_file),form='unformatted')
  write(IOUT) total_model(:,:,:,1:nspec)
  close(IOUT)


  !---------------------------------------------------------------------------------------------
  ! compute bulk wavespeed of the current model and new model
  !---------------------------------------------------------------------------------------------

  ! bulk wavespeed of current model
  total_model = 0.
  total_model = sqrt( model_vp**2 - (4.0/3.0)*model_vs**2)
  fname = 'vb'
  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  open(IOUT,file=trim(m_file),form='unformatted')
  write(IOUT) total_model(:,:,:,1:nspec)
  close(IOUT)

  ! bulk wavespeed of new model
  total_model = 0.
  total_model = sqrt( model_vp_new**2 - (4.0/3.0)*model_vs_new**2)
  fname = 'vb_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  open(IOUT,file=trim(m_file),form='unformatted')
  write(IOUT) total_model(:,:,:,1:nspec)
  close(IOUT)

!===================================================
!===================================================
  deallocate(model_vp, model_vs, model_vp_new, model_vs_new, &
             model_dbulk, model_dbeta, total_model, &
             kernel_bulk, kernel_beta)
  deallocate(kernel_rho)
  deallocate(model_rho, model_rho_new, model_drho)

  ! stop all the processes, and exit
  call finalize_mpi()

end program model_update

!
!-------------------------------------------------------------------------------------------------
!

subroutine initialize()

! initializes arrays

  use tomography_par, only: myrank_tomo => myrank, sizeprocs, NSPEC, NGLOB, USE_ALPHA_BETA_RHO

  use specfem_par, only: NSPEC_AB,NGLOB_AB,NPROC,myrank,ADIOS_ENABLED,ATTENUATION

  implicit none

  logical :: BROADCAST_AFTER_READ

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! reads the parameter file
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(BROADCAST_AFTER_READ)

  if (ADIOS_ENABLED) stop 'Flag ADIOS_ENABLED set to .true. not supported yet for xmodel_update, please rerun program...'

  ! security check
  if (ATTENUATION) then
    print *,'Sorry using ATTENUATION, this routine has qkappa not implemented yet...'
    stop 'Error ATTENUATION flag invalid'
  endif

  if (USE_ALPHA_BETA_RHO .neqv. .true.) then
    print *,'Sorry using USE_ALPHA_BETA_RHO must be set to .true. in constants_tomography.h file; Please recompile...'
    stop 'Error USE_ALPHA_BETA_RHO flag invalid'
  endif

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *, 'Error number of processors supposed to run on: ',NPROC
      print *, 'Error number of MPI processors actually run on: ',sizeprocs
      print *
      print *, 'please rerun with: mpirun -np ',NPROC,' bin/xmodel_update .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif

  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some array sizes below
  call read_mesh_for_init()

  ! sets tomography array dimensions
  NSPEC = NSPEC_AB
  NGLOB = NGLOB_AB
  myrank_tomo = myrank

end subroutine initialize

!
!-------------------------------------------------------------------------------------------------
!

subroutine get_external_mesh()

  use specfem_par, only: CUSTOM_REAL,NSPEC_AB,NGLOB_AB,NSPEC_IRREGULAR,NGLLX,NGLLY,NGLLZ, &
    model_speed_max,DT,myrank,IMAIN,ISTANDARD_OUTPUT

  use specfem_par, only: ibool,xstore,ystore,zstore, &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore, &
    kappastore,mustore,rhostore

  use specfem_par_elastic, only: ispec_is_elastic,min_resolved_period
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_poroelastic, only: ispec_is_poroelastic,rho_vpI,rho_vpII,rho_vsI, &
    phistore,tortstore,rhoarraystore

  use tomography_par, only: OUTPUT_MODEL_DIR

  implicit none
  integer :: ier
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  ! allocate arrays for storing the databases
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 918')

  if (NSPEC_IRREGULAR > 0) then
    allocate(xixstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 919')
    allocate(xiystore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 920')
    allocate(xizstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 921')
    allocate(etaxstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 922')
    allocate(etaystore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 923')
    allocate(etazstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 924')
    allocate(gammaxstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 925')
    allocate(gammaystore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 926')
    allocate(gammazstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 927')
    allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 928')
  else
    allocate(xixstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 929')
    allocate(xiystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 930')
    allocate(xizstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 931')
    allocate(etaxstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 932')
    allocate(etaystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 933')
    allocate(etazstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 934')
    allocate(gammaxstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 935')
    allocate(gammaystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 936')
    allocate(gammazstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 937')
    allocate(jacobianstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 938')
  endif
  if (ier /= 0) stop 'Error allocating arrays for databases'

  ! mesh node locations
  allocate(xstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 939')
  allocate(ystore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 940')
  allocate(zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 941')
  if (ier /= 0) stop 'Error allocating arrays for mesh nodes'

  ! material properties
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 942')
  allocate(mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 943')
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating rho array 943')
  if (ier /= 0) stop 'Error allocating arrays for material properties'

  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 944')
  allocate(ispec_is_elastic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 945')
  allocate(ispec_is_poroelastic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 946')
  if (ier /= 0) stop 'Error allocating arrays for material flags'
  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! reads in external mesh
  call read_mesh_databases()

  ! gets mesh dimensions
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                            x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            distance_min_glob,distance_max_glob)

  ! outputs infos
  if (myrank == 0) then
    print *,'mesh dimensions:'
    print *,'  Xmin and Xmax of the model = ',x_min_glob,x_max_glob
    print *,'  Ymin and Ymax of the model = ',y_min_glob,y_max_glob
    print *,'  Zmin and Zmax of the model = ',z_min_glob,z_max_glob
    print *
    print *,'  Max GLL point distance = ',distance_max_glob
    print *,'  Min GLL point distance = ',distance_min_glob
    print *,'  Max/min ratio = ',distance_max_glob/distance_min_glob
    print *
    print *,'  Max element size = ',elemsize_max_glob
    print *,'  Min element size = ',elemsize_min_glob
    print *,'  Max/min ratio = ',elemsize_max_glob/elemsize_min_glob
    print *
  endif
  call synchronize_all()

  ! resolution check
  ! open main output file, only written to by process 0
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_MODEL_DIR)//'/output_mesh_resolution_initial.txt',status='unknown')

  ! mesh resolution
  call check_mesh_resolution(NSPEC_AB,NGLOB_AB, &
                             ibool,xstore,ystore,zstore, &
                             ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
                             kappastore,mustore,rhostore, &
                             phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                             DT,model_speed_max,min_resolved_period)

  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

end subroutine get_external_mesh

!
!-------------------------------------------------------------------------------------------------
!

subroutine save_new_databases()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_poroelastic, only: ispec_is_poroelastic,rho_vpI,rho_vpII,rho_vsI, &
    phistore,tortstore,rhoarraystore

  use tomography_model_iso, only: model_vs_new,model_vp_new,model_rho_new,OUTPUT_MODEL_DIR

  implicit none

  character(len=MAX_STRING_LEN) :: prname_new
  character(len=MAX_STRING_LEN*2) :: m_file
  integer :: NGLOB_OCEAN

  ! for attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qmu_attenuation_store,qkappa_attenuation_store  ! attenuation
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dummy_g_1,dummy_g_2,dummy_g_3  !xstore,ystore,zstore
  integer, dimension(:), allocatable :: dummy_l_1,dummy_l_2,dummy_l_3,dummy_l_4,dummy_l_5,dummy_l_6,dummy_l_7,dummy_l_8 !ibool-1
  integer, dimension(:), allocatable :: dummy_num
  character(len=MAX_STRING_LEN) :: string1,string2,string3,string4,string5,string6,string7,string8,string9,string10,string11
  integer :: idummy1,idummy2,idummy3,idummy4,idummy5

  integer :: iglob
  real, dimension(:),allocatable :: flag_val
  logical, dimension(:),allocatable :: mask_ibool

  !MPI variables
  integer :: ier

  ! to update the values in external_mesh.bin
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kappastore_new,mustore_new
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp_new,rho_vs_new
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore_new

  integer :: i,j,k,ispec

  ! calculate rmass
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: jacobian_new
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl
  ! mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_old, rmass_new
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_acoustic_new,rmass_solid_poroelastic_new,rmass_fluid_poroelastic_new

  ! user output
  if (myrank == 0) then
    print *,'saving new databases'
  endif

  NGLOB_OCEAN = NGLOB_AB

  ! saves binary mesh files for the NEW model
  call create_name_database(prname_new,myrank,OUTPUT_MODEL_DIR)

  ! new variables for save_external_bin_m_up
  allocate(kappastore_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 951')
  allocate(mustore_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 952')
  allocate(rhostore_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 953')
  allocate(rho_vp_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 954')
  allocate(rho_vs_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 955')

  ! calculate NEW variables to calculate rmass and then for save_external_bin_m_up
  rhostore_new = 0._CUSTOM_REAL
  kappastore_new = 0._CUSTOM_REAL
  mustore_new = 0._CUSTOM_REAL
  rho_vp_new = 0._CUSTOM_REAL
  rho_vs_new = 0._CUSTOM_REAL

  rhostore_new =  model_rho_new
  kappastore_new = model_rho_new * ( (model_vp_new**2) - FOUR_THIRDS * (model_vs_new**2) )
  mustore_new = model_rho_new * model_vs_new * model_vs_new
  rho_vp_new = model_rho_new * model_vp_new
  rho_vs_new = model_rho_new * model_vs_new

  ! jacobian from read_mesh_databases
  ! safety check
  if (NSPEC_IRREGULAR /= NSPEC_AB) stop 'Please check if model_update in save_new_databases() with NSPEC_AB /= NSPEC_IRREGULAR'

  allocate(jacobian_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 956')
  jacobian_new = 0._CUSTOM_REAL
  jacobian_new = jacobianstore

  ! set up coordinates of the Gauss-Lobatto-Legendre points and calculate weights
  ! from constants.h GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! rmass for the OLD model from read_mesh_databases
  allocate(rmass_old(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 957')
  rmass_old = 0._CUSTOM_REAL
  rmass_old = rmass

  ! create mass matrix ONLY for the elastic case
  allocate(rmass_new(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 958')
  rmass_new(:) = 0._CUSTOM_REAL

  ! user output
  if (myrank == 0) then
    print *, '  ...creating mass matrix '
  endif
  call synchronize_all()

  ! note: this does not update the absorbing boundary contributions to the mass matrix
  ! elastic mass matrix
  do ispec=1,NSPEC_AB
    if (ispec_is_elastic(ispec)) then
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)

            weight = wxgll(i)*wygll(j)*wzgll(k)
            jacobianl = jacobian_new(i,j,k,ispec)

            !debug
            !if (myrank == 0) then
            !  print *, 'weight', weight
            !  print *, 'jacobianl', jacobianl
            !endif

            rmass_new(iglob) = rmass_new(iglob) + &
                      real( dble(jacobianl) * weight * dble(rhostore_new(i,j,k,ispec)) ,kind=CUSTOM_REAL)
          enddo
        enddo
      enddo
    endif
  enddo ! nspec
  call synchronize_all()

  ! dummy allocations, arrays are not needed since the update here only works for elastic models
  allocate(rmass_acoustic_new(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 959')
  allocate(rmass_solid_poroelastic_new(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 960')
  allocate(rmass_fluid_poroelastic_new(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 961')
  rmass_acoustic_new = 0._CUSTOM_REAL
  rmass_solid_poroelastic_new = 0._CUSTOM_REAL
  rmass_fluid_poroelastic_new = 0._CUSTOM_REAL

  ! safety check
  if (POROELASTIC_SIMULATION) &
    stop 'Error saving new databases for POROELASTIC models not implemented yet'
  if (ACOUSTIC_SIMULATION) &
    stop 'Error saving new databases for ACOUSTIC models not implemented yet'

  ! new resolution check
  ! open main output file, only written to by process 0
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_MODEL_DIR)//'/output_mesh_resolution_final.txt',status='unknown')

  ! calculate min_resolved_period (needed for attenuation model)
  call check_mesh_resolution(NSPEC_AB,NGLOB_AB, &
                             ibool,xstore,ystore,zstore, &
                             ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
                             kappastore_new,mustore_new,rhostore_new, &
                             phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                             -1.0d0,model_speed_max,min_resolved_period)

  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  !-------- attenuation -------
  ! store the attenuation flag in qmu_attenuation_store
  allocate(qmu_attenuation_store(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 962')
  allocate(qkappa_attenuation_store(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 963')
  qmu_attenuation_store(:,:,:,:) = 0._CUSTOM_REAL
  qkappa_attenuation_store(:,:,:,:) = 0._CUSTOM_REAL

  if (ATTENUATION) then
    ! read the proc*attenuation.vtk for the old model in LOCAL_PATH and store qmu_attenuation_store
    write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_attenuation.vtk'
    open(IIN,file=trim(m_file),status='old',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'Error file not found')
    endif
    read(IIN,'(a)') string1 !text
    read(IIN,'(a)') string2 !text
    read(IIN,'(a)') string3 !text
    read(IIN,'(a)') string4 !text
    read(IIN,'(a,i12,a)') string5, idummy1, string6 !text

    allocate(dummy_g_1(NGLOB_AB),dummy_g_2(NGLOB_AB),dummy_g_3(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 964')
    if (ier /= 0) stop 'Error allocating array dummy etc.'

    read(IIN,'(3e18.6)') dummy_g_1,dummy_g_2,dummy_g_3 !xstore,ystore,zstore for i=1,nglob
    read(IIN,*) !blank line
    deallocate(dummy_g_1,dummy_g_2,dummy_g_3)

    read(IIN,'(a,i12,i12)') string7, idummy2, idummy3 !text

    allocate(dummy_num(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 965')
    allocate(dummy_l_1(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 966')
    allocate(dummy_l_2(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 967')
    allocate(dummy_l_3(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 968')
    allocate(dummy_l_4(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 969')
    allocate(dummy_l_5(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 970')
    allocate(dummy_l_6(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 971')
    allocate(dummy_l_7(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 972')
    allocate(dummy_l_8(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 973')
    if (ier /= 0) stop 'Error allocating array dummy etc.'

    read(IIN,'(9i12)') dummy_num,dummy_l_1,dummy_l_2,dummy_l_3,dummy_l_4, &
                    dummy_l_5,dummy_l_6,dummy_l_7,dummy_l_8 !8,ibool-1 for ispec=1,nspec
    read(IIN,*) !blank line

    deallocate(dummy_num,dummy_l_1,dummy_l_2,dummy_l_3,dummy_l_4,dummy_l_5,dummy_l_6,dummy_l_7,dummy_l_8)

    read(IIN,'(a,i12)') string8, idummy4 !text

    allocate(dummy_num(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 974')
    if (ier /= 0) stop 'Error allocating array dummy etc.'

    read(IIN,*) dummy_num !12 for ispec=1,nspec
    read(IIN,*) !blank line

    deallocate(dummy_num)

    read(IIN,'(a,i12)') string9, idummy5 !text
    read(IIN,'(a)') string10 !text
    read(IIN,'(a)') string11 !text

    allocate(flag_val(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 975')
    if (ier /= 0) stop 'Error allocating flag_val'

    read(IIN,*) flag_val
    read(IIN,*) !blank line
    close(IIN)

    allocate(mask_ibool(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 976')
    if (ier /= 0) stop 'Error allocating mask'

    mask_ibool = .false.
    do ispec=1,NSPEC_AB
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            if (.not. mask_ibool(iglob)) then
              qmu_attenuation_store(i,j,k,ispec) = flag_val(iglob)
              mask_ibool(iglob) = .true.
            endif
          enddo
        enddo
      enddo
    enddo
    call synchronize_all()

    ! calculates and stores attenuation arrays
    call get_attenuation_model(NSPEC_AB,USE_OLSEN_ATTENUATION,OLSEN_ATTENUATION_RATIO, &
                               mustore_new,rho_vs_new,kappastore_new,rho_vp_new, &
                               qkappa_attenuation_store,qmu_attenuation_store, &
                               ispec_is_elastic,min_resolved_period,prname_new,ATTENUATION_f0_REFERENCE)

    deallocate(flag_val,mask_ibool)
  endif

  !----------------------------

  ! user output
  if (myrank == 0) then
    print *,'  writing new databases to directory: ',trim(OUTPUT_MODEL_DIR)
    print *
  endif
  call synchronize_all()

  call save_external_bin_m_up(NSPEC_AB,NGLOB_AB, &
                        rho_vp_new,rho_vs_new,qmu_attenuation_store, &
                        rhostore_new,kappastore_new,mustore_new, &
                        rmass_new,rmass_acoustic_new,rmass_solid_poroelastic_new,rmass_fluid_poroelastic_new, &
                        APPROXIMATE_OCEAN_LOAD,rmass_ocean_load,NGLOB_OCEAN,ibool,xstore,ystore,zstore, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec,num_abs_boundary_faces, &
                        free_surface_normal,free_surface_jacobian2Dw, &
                        free_surface_ijk,free_surface_ispec,num_free_surface_faces, &
                        num_interfaces_ext_mesh,my_neighbors_ext_mesh,nibool_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        prname_new,SAVE_MESH_FILES,ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store, &
                        c22store,c23store,c24store,c25store,c26store,c33store, &
                        c34store,c35store,c36store,c44store,c45store,c46store, &
                        c55store,c56store,c66store, &
                        ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic)

  deallocate(rhostore_new, kappastore_new, mustore_new, rho_vp_new, rho_vs_new)
  deallocate(jacobian_new)
  deallocate(rmass_old,rmass_new)
  deallocate(rmass_acoustic_new,rmass_solid_poroelastic_new,rmass_fluid_poroelastic_new)
  deallocate(qmu_attenuation_store,qkappa_attenuation_store)

end subroutine save_new_databases

