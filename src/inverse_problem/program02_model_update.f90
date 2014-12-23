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

program model_update

  use specfem_par,only: kappastore,mustore !,ibool
  use specfem_par,only: NPROC,OUTPUT_FILES_PATH,LOCAL_PATH
  use specfem_par_elastic,only: rho_vp,rho_vs

  use inverse_problem_model_iso
  use inverse_problem_kernels_iso

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
  if (OUTPUT_FILES_PATH(len_trim(OUTPUT_FILES_PATH):len_trim(OUTPUT_FILES_PATH)) /= '/') then
    INPUT_KERNELS_DIR = trim(OUTPUT_FILES_PATH)//'/'//trim(INPUT_KERNELS_DIR_NAME)
  else
    INPUT_KERNELS_DIR = trim(OUTPUT_FILES_PATH)//trim(INPUT_KERNELS_DIR_NAME)
  endif

  ! directory where the mesh files for the NEW model will be written
  if (OUTPUT_FILES_PATH(len_trim(OUTPUT_FILES_PATH):len_trim(OUTPUT_FILES_PATH)) /= '/') then
    OUTPUT_MODEL_DIR = trim(OUTPUT_FILES_PATH)//'/'//trim(LOCAL_PATH_NEW_NAME)
  else
    OUTPUT_MODEL_DIR = trim(OUTPUT_FILES_PATH)//trim(LOCAL_PATH_NEW_NAME)
  endif

  ! directory where the output files of model_update will be written
  PRINT_STATISTICS_FILES = .true.

  if (OUTPUT_FILES_PATH(len_trim(OUTPUT_FILES_PATH):len_trim(OUTPUT_FILES_PATH)) /= '/') then
    OUTPUT_STATISTICS_DIR = trim(OUTPUT_FILES_PATH)//'/'//trim(OUTPUT_STATISTICS_DIR_NAME)
  else
    OUTPUT_STATISTICS_DIR = trim(OUTPUT_FILES_PATH)//trim(OUTPUT_STATISTICS_DIR_NAME)
  endif

  ! reads in parameters needed
  call read_parameters_invprob()

  ! user output
  if (myrank == 0) then
    print*
    print*,'***********'
    print*,'program model_update: '
    print*,'  NPROC: ',NPROC
    print*,'  NSPEC: ', NSPEC
    print*,'  NGLOB: ', NGLOB
    print*
    print*,'model update for vs & vp & rho'
    print*,'  step_fac = ',step_fac
    print*
    if (USE_ALPHA_BETA_RHO) then
      print*,'kernel parameterization: (alpha,beta,rho)'
    else
      print*,'kernel parameterization: (bulk,beta,rho)'
    endif
    print*
    if (USE_RHO_SCALING) then
      print*,'scaling rho perturbations'
      print*
    endif
    if (MINMAX_THRESHOLD_OLD) print*,'thresholding current (old) model wavespeed values'
    if (MINMAX_THRESHOLD_NEW) print*,'thresholding new (updated) model wavespeed values'
    if (MINMAX_THRESHOLD_OLD .or. MINMAX_THRESHOLD_NEW) then
      if (THRESHOLD_RHO) then
        print*,'  thresholds for: (vs,vp,rho)'
      else
        print*,'  thresholds for: (vs,vp)'
      endif
    endif
    print*,'***********'
    print*
  endif
  call synchronize_all()

  ! reads in external mesh
  call get_external_mesh()

  !===================================================
  ! MODEL UPDATE
  ! allocation
  ! model and kernel variables
  allocate(model_vp(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_vs(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_rho(NGLLX,NGLLY,NGLLZ,NSPEC))

  allocate(total_model(NGLLX,NGLLY,NGLLZ,NSPEC))

  ! initialize variables
  ! old model
  model_vp = 0._CUSTOM_REAL
  model_vs = 0._CUSTOM_REAL
  model_rho = 0._CUSTOM_REAL

  !! set up thresholds for old and new models

  if (MINMAX_THRESHOLD_OLD .or. MINMAX_THRESHOLD_NEW) then
    if (myrank == 0) then
      print*,'threshold wavespeed values:'
      print*,'  VS_MIN, VS_MAX  : ',VS_MIN, VS_MAX
      print*,'  VP_MIN, VP_MAX  : ',VP_MIN, VP_MAX
      if (THRESHOLD_RHO) print*,'  RHO_MIN, RHO_MAX: ',RHO_MIN, RHO_MAX
      print*
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
      print*,'current model values before thresholding:'
    else
      print*,'current model values:'
    endif
    print*,'  vs min/max : ',vsmin_before, vsmax_before
    print*,'  vp min/max : ',vpmin_before, vpmax_before
    print*,'  rho min/max: ',rhomin_before, rhomax_before
    print*
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
      print*,'current model values after thresholding:'
      print*,'  vs min/max: ',vsmin_after, vsmax_after
      print*,'  vp min/max: ',vpmin_after, vpmax_after
      print*,'  rho min/max: ',rhomin_after, rhomax_after
      print*
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
  ! calculates gradient
  !---------------------------------------------------------------------------------------------
  call get_gradient_steepest_iso()

  !---------------------------------------------------------------------------------------------
  ! new model
  !---------------------------------------------------------------------------------------------

  ! computes new model values for alpha, beta and rho
  ! allocate new model arrays
  allocate(model_vp_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_vs_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_rho_new(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
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
      print*,'new model values before thresholding:'
    else
      print*,'new model values:'
    endif
    print*,'  vs min/max : ',vsmin_new_before, vsmax_new_before
    print*,'  vp min/max : ',vpmin_new_before, vpmax_new_before
    print*,'  rho min/max: ',rhomin_new_before, rhomax_new_before
    print*
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
      print*,'new model values after thresholding:'
      print*,'  vs min/max : ',vsmin_new_after, vsmax_new_after
      print*,'  vp min/max : ',vpmin_new_after, vpmax_new_after
      print*,'  rho min/max: ',rhomin_new_after, rhomax_new_after
      print*
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

! vtk file output
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
  total_model = ( model_vp**2 - 2.0*model_vs**2) / ( 2.0*model_vp**2 - 2.0*model_vs**2)
  fname = 'poisson'
  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

  ! Poisson's ratio of new model
  total_model = 0.
  total_model = ( model_vp_new**2 - 2.0*model_vs_new**2) / &
                ( 2.0*model_vp_new**2 - 2.0*model_vs_new**2)
  fname = 'poisson_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)


  !---------------------------------------------------------------------------------------------
  ! compute bulk wavespeed of the current model and new model
  !---------------------------------------------------------------------------------------------

  ! bulk wavespeed of current model
  total_model = 0.
  total_model = sqrt( model_vp**2 - (4.0/3.0)*model_vs**2)
  fname = 'vb'
  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

  ! bulk wavespeed of new model
  total_model = 0.
  total_model = sqrt( model_vp_new**2 - (4.0/3.0)*model_vs_new**2)
  fname = 'vb_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

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

  use inverse_problem_par,only: myrank_invprob => myrank, sizeprocs, NSPEC, NGLOB, USE_ALPHA_BETA_RHO

  use specfem_par,only: &
    NSPEC_AB,NGLOB_AB,NPROC,myrank, &
    ADIOS_ENABLED,ATTENUATION

  implicit none

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! reads the parameter file
  call read_parameter_file()

  ! reads ADIOS flags
  call read_adios_parameters()

  if (ADIOS_ENABLED) stop 'Flag ADIOS_ENABLED set to .true. not supported yet for xmodel_update, please rerun program...'

  ! security check
  if (ATTENUATION) then
    print*,'Sorry using ATTENUATION, this routine has qkappa not implemented yet...'
    stop 'Error ATTENUATION flag invalid'
  endif

  if (USE_ALPHA_BETA_RHO .neqv. .true.) then
    print*,'Sorry using USE_ALPHA_BETA_RHO must be set to .true. in constants_inverse_problem.h file; Please recompile...'
    stop 'Error USE_ALPHA_BETA_RHO flag invalid'
  endif

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print*, 'Error number of processors supposed to run on: ',NPROC
      print*, 'Error number of MPI processors actually run on: ',sizeprocs
      print*
      print*, 'please rerun with: mpirun -np ',NPROC,' bin/xmodel_update .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif

  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some array sizes below
  call read_mesh_for_init()

  ! sets inverse_problem array dimensions
  NSPEC = NSPEC_AB
  NGLOB = NGLOB_AB
  myrank_invprob = myrank

end subroutine initialize

!
!-------------------------------------------------------------------------------------------------
!

subroutine get_external_mesh()

  use specfem_par,only: CUSTOM_REAL,NSPEC_AB,NGLOB_AB,NGLLX,NGLLY,NGLLZ, &
    LOCAL_PATH,SAVE_MESH_FILES,model_speed_max,DT,myrank,IMAIN,ISTANDARD_OUTPUT

  use specfem_par,only: ibool,xstore,ystore,zstore,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian, &
    kappastore,mustore,rhostore

  use specfem_par_elastic,only: ELASTIC_SIMULATION,ispec_is_elastic,rho_vp,rho_vs,min_resolved_period
  use specfem_par_acoustic,only: ACOUSTIC_SIMULATION,ispec_is_acoustic
  use specfem_par_poroelastic,only: POROELASTIC_SIMULATION,ispec_is_poroelastic,rho_vpI,rho_vpII,rho_vsI, &
    phistore,tortstore,rhoarraystore

  use inverse_problem_par,only: OUTPUT_MODEL_DIR

  implicit none
  integer :: ier
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  ! allocate arrays for storing the databases
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xix(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xiy(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xiz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etax(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etay(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etaz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammax(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammay(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for databases'

  ! mesh node locations
  allocate(xstore(NGLOB_AB), &
           ystore(NGLOB_AB), &
           zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for mesh nodes'

  ! material properties
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for material properties'

  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB), &
           ispec_is_elastic(NSPEC_AB), &
           ispec_is_poroelastic(NSPEC_AB),stat=ier)
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
    print*,'mesh dimensions:'
    print*,'  Xmin and Xmax of the model = ',x_min_glob,x_max_glob
    print*,'  Ymin and Ymax of the model = ',y_min_glob,y_max_glob
    print*,'  Zmin and Zmax of the model = ',z_min_glob,z_max_glob
    print*
    print*,'  Max GLL point distance = ',distance_max_glob
    print*,'  Min GLL point distance = ',distance_min_glob
    print*,'  Max/min ratio = ',distance_max_glob/distance_min_glob
    print*
    print*,'  Max element size = ',elemsize_max_glob
    print*,'  Min element size = ',elemsize_min_glob
    print*,'  Max/min ratio = ',elemsize_max_glob/elemsize_min_glob
    print*
  endif
  call synchronize_all()

  ! resolution check
  ! open main output file, only written to by process 0
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_MODEL_DIR)//'/output_mesh_resolution_initial.txt',status='unknown')

  if (ELASTIC_SIMULATION) then
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                               ibool,xstore,ystore,zstore, &
                               kappastore,mustore,rho_vp,rho_vs, &
                               DT,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)

  else if (POROELASTIC_SIMULATION) then
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    rho_vp = 0.0_CUSTOM_REAL
    rho_vs = 0.0_CUSTOM_REAL
    call check_mesh_resolution_poro(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                    DT,model_speed_max,min_resolved_period, &
                                    phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                                    LOCAL_PATH,SAVE_MESH_FILES)
    deallocate(rho_vp,rho_vs)
  else if (ACOUSTIC_SIMULATION) then
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rho_vp'
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rho_vs'
    rho_vp = sqrt( kappastore / rhostore ) * rhostore
    rho_vs = 0.0_CUSTOM_REAL
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                               ibool,xstore,ystore,zstore, &
                               kappastore,mustore,rho_vp,rho_vs, &
                               DT,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)
    deallocate(rho_vp,rho_vs)
  endif

  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

end subroutine get_external_mesh

!
!-------------------------------------------------------------------------------------------------
!

subroutine save_new_databases()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic,only:ACOUSTIC_SIMULATION,ispec_is_acoustic
  use specfem_par_poroelastic,only:POROELASTIC_SIMULATION,ispec_is_poroelastic

  use inverse_problem_model_iso,only: model_vs_new,model_vp_new,model_rho_new,OUTPUT_MODEL_DIR

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
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: jacobianstore
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl
  ! mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_old, rmass_new
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_acoustic_new,rmass_solid_poroelastic_new,rmass_fluid_poroelastic_new

  ! user output
  if (myrank == 0) then
    print*,'saving new databases'
  endif

  NGLOB_OCEAN = NGLOB_AB

  ! saves binary mesh files for the NEW model
  call create_name_database(prname_new,myrank,OUTPUT_MODEL_DIR)

  ! new variables for save_external_bin_m_up
  allocate(kappastore_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           mustore_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           rhostore_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           rho_vp_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           rho_vs_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB))

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
  allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  jacobianstore = 0._CUSTOM_REAL
  jacobianstore = jacobian

  ! set up coordinates of the Gauss-Lobatto-Legendre points and calculate weights
  ! from constants.h GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! rmass for the OLD model from read_mesh_databases
  allocate(rmass_old(NGLOB_AB))
  rmass_old = 0._CUSTOM_REAL
  rmass_old = rmass

  ! create mass matrix ONLY for the elastic case
  allocate(rmass_new(NGLOB_AB))
  rmass_new(:) = 0._CUSTOM_REAL

  ! user output
  if (myrank == 0) then
    print*, '  ...creating mass matrix '
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
            jacobianl = jacobianstore(i,j,k,ispec)

            !debug
            !if (myrank == 0) then
            !  print*, 'weight', weight
            !  print*, 'jacobianl', jacobianl
            !endif

            if (CUSTOM_REAL == SIZE_REAL) then
              rmass_new(iglob) = rmass_new(iglob) + &
                      sngl( dble(jacobianl) * weight * dble(rhostore_new(i,j,k,ispec)) )
            else
              rmass_new(iglob) = rmass_new(iglob) + &
                      jacobianl * weight * rhostore_new(i,j,k,ispec)
            endif
          enddo
        enddo
      enddo
    endif
  enddo ! nspec
  call synchronize_all()

  ! dummy allocations, arrays are not needed since the update here only works for elastic models
  allocate(rmass_acoustic_new(NGLOB_AB))
  allocate(rmass_solid_poroelastic_new(NGLOB_AB))
  allocate(rmass_fluid_poroelastic_new(NGLOB_AB))
  rmass_acoustic_new = 0._CUSTOM_REAL
  rmass_solid_poroelastic_new = 0._CUSTOM_REAL
  rmass_fluid_poroelastic_new = 0._CUSTOM_REAL


  ! new resolution check
  ! open main output file, only written to by process 0
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_MODEL_DIR)//'/output_mesh_resolution_final.txt',status='unknown')

  ! calculate min_resolved_period (needed for attenuation model)
  if (ELASTIC_SIMULATION) then
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                               ibool,xstore,ystore,zstore, &
                               kappastore_new,mustore_new,rho_vp_new,rho_vs_new, &
                               -1.0d0,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)

  else if (POROELASTIC_SIMULATION) then
    stop 'Error saving new databases for POROELASTIC models not implemented yet'
  else if (ACOUSTIC_SIMULATION) then
    stop 'Error saving new databases for ACOUSTIC models not implemented yet'
  endif

  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  !-------- attenuation -------
  ! store the attenuation flag in qmu_attenuation_store
  allocate(qmu_attenuation_store(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           qkappa_attenuation_store(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  qmu_attenuation_store=0._CUSTOM_REAL
  qkappa_attenuation_store=0._CUSTOM_REAL

  if (ATTENUATION) then
    ! read the proc*attenuation.vtk for the old model in LOCAL_PATH and store qmu_attenuation_store
    write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_attenuation.vtk'
    open(12,file=trim(m_file),status='old',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'Error file not found')
    endif
    read(12,'(a)') string1 !text
    read(12,'(a)') string2 !text
    read(12,'(a)') string3 !text
    read(12,'(a)') string4 !text
    read(12,'(a,i12,a)') string5, idummy1, string6 !text

    allocate(dummy_g_1(NGLOB_AB),dummy_g_2(NGLOB_AB),dummy_g_3(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array dummy etc.'

    read(12,'(3e18.6)') dummy_g_1,dummy_g_2,dummy_g_3 !xstore,ystore,zstore for i=1,nglob
    read(12,*) !blank line
    deallocate(dummy_g_1,dummy_g_2,dummy_g_3)

    read(12,'(a,i12,i12)') string7, idummy2, idummy3 !text

    allocate(dummy_num(NSPEC_AB),dummy_l_1(NSPEC_AB),dummy_l_2(NSPEC_AB),dummy_l_3(NSPEC_AB),dummy_l_4(NSPEC_AB), &
             dummy_l_5(NSPEC_AB),dummy_l_6(NSPEC_AB),dummy_l_7(NSPEC_AB),dummy_l_8(NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array dummy etc.'

    read(12,'(9i12)') dummy_num,dummy_l_1,dummy_l_2,dummy_l_3,dummy_l_4, &
                    dummy_l_5,dummy_l_6,dummy_l_7,dummy_l_8 !8,ibool-1 for ispec=1,nspec
    read(12,*) !blank line

    deallocate(dummy_num,dummy_l_1,dummy_l_2,dummy_l_3,dummy_l_4,dummy_l_5,dummy_l_6,dummy_l_7,dummy_l_8)

    read(12,'(a,i12)') string8, idummy4 !text

    allocate(dummy_num(NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array dummy etc.'

    read(12,*) dummy_num !12 for ispec=1,nspec
    read(12,*) !blank line

    deallocate(dummy_num)

    read(12,'(a,i12)') string9, idummy5 !text
    read(12,'(a)') string10 !text
    read(12,'(a)') string11 !text

    allocate(flag_val(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating flag_val'

    read(12,*) flag_val
    read(12,*) !blank line
    close(12)

    allocate(mask_ibool(NGLOB_AB),stat=ier)
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
    call get_attenuation_model(myrank,NSPEC_AB,USE_OLSEN_ATTENUATION,OLSEN_ATTENUATION_RATIO, &
                               mustore_new,rho_vs_new,kappastore_new,rho_vp_new, &
                               qkappa_attenuation_store,qmu_attenuation_store, &
                               ispec_is_elastic,min_resolved_period,prname_new,FULL_ATTENUATION_SOLID)

    deallocate(flag_val,mask_ibool)
  endif

  !----------------------------

  ! user output
  if (myrank == 0) then
    print*,'  writing new databases to directory: ',trim(OUTPUT_MODEL_DIR)
    print*
  endif
  call synchronize_all()

  call save_external_bin_m_up(NSPEC_AB,NGLOB_AB, &
                        xix,xiy,xiz,etax,etay,etaz,&
                        gammax,gammay,gammaz, &
                        jacobianstore,rho_vp_new,rho_vs_new,qmu_attenuation_store, &
                        rhostore_new,kappastore_new,mustore_new, &
                        rmass_new,rmass_acoustic_new,rmass_solid_poroelastic_new,rmass_fluid_poroelastic_new, &
                        APPROXIMATE_OCEAN_LOAD,rmass_ocean_load,NGLOB_OCEAN,ibool,xstore,ystore,zstore, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec,num_abs_boundary_faces, &
                        free_surface_normal,free_surface_jacobian2Dw, &
                        free_surface_ijk,free_surface_ispec,num_free_surface_faces, &
                        num_interfaces_ext_mesh,my_neighbours_ext_mesh,nibool_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        prname_new,SAVE_MESH_FILES,ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store, &
                        c22store,c23store,c24store,c25store,c26store,c33store, &
                        c34store,c35store,c36store,c44store,c45store,c46store, &
                        c55store,c56store,c66store, &
                        ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic)

  deallocate(rhostore_new, kappastore_new, mustore_new, rho_vp_new, rho_vs_new)
  deallocate(jacobianstore)
  deallocate(rmass_old,rmass_new)
  deallocate(rmass_acoustic_new,rmass_solid_poroelastic_new,rmass_fluid_poroelastic_new)
  deallocate(qmu_attenuation_store,qkappa_attenuation_store)

end subroutine save_new_databases

