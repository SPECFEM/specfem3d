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


program model_update

  use :: mpi
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  implicit none

  include 'precision.h'

  ! ======================================================
  ! USER PARAMETERS

  ! root file directory
  character (len=256), parameter :: &
    ROOT_PATH = '/lscratch/users/magnoni/SPECFEM3D/trunk_update/OUTPUT_FILES'

  ! directory where the mesh files for the NEW model will be written
  character (len=256), parameter :: &
    LOCAL_PATH_NEW = trim(ROOT_PATH)//'/'//'DATABASES_MPI/mesh_files_m01'

  ! directory where the output files of model_update will be written
  character (len=256), parameter :: &
    OUTPUT_MODEL_UPD = trim(ROOT_PATH)//'/'//'OUTPUT_FILES_MODEL_UPD'

  ! directory where the summed and smoothed input kernels are linked
  character (len=256), parameter :: &
    INPUT_KERNELS = trim(ROOT_PATH)//'/'//'DATABASES_MPI/sum_smooth_kern'

  ! by default, this algorithm uses (bulk,bulk_beta,rho) kernels to update vp,vs,rho
  ! if you prefer using (alpha,beta,rho) kernels, set this flag to true
  logical, parameter :: USE_ALPHA_BETA_RHO = .true.

  ! ignore rho kernel, but use density perturbations as a scaling of Vs perturbations
  logical, parameter :: USE_RHO_SCALING = .false.

  ! in case of rho scaling, specifies density scaling factor with shear perturbations
  ! see e.g. Montagner & Anderson (1989), Panning & Romanowicz (2006)
  real(kind=CUSTOM_REAL),parameter :: RHO_SCALING = 0.33_CUSTOM_REAL

  !set to true if you want to print out log files with statistics
  logical, parameter :: PRINT_OUT_FILES = .true.

  logical, parameter :: MINMAX_THRESHOLD_OLD = .false.   ! threshold the old model ("current model")
  logical, parameter :: MINMAX_THRESHOLD_NEW = .false.    ! threshold the new model

  ! ======================================================

  character (len=256) :: prname_new
  character(len=256) :: m_file, fname
  integer :: NGLOB_OCEAN
  integer :: NSPEC, NGLOB

  ! for attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qmu_attenuation_store  ! attenuation
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dummy_g_1,dummy_g_2,dummy_g_3  !xstore,ystore,zstore
  integer, dimension(:), allocatable :: dummy_l_1,dummy_l_2,dummy_l_3,dummy_l_4,dummy_l_5,dummy_l_6,dummy_l_7,dummy_l_8 !ibool-1
  integer, dimension(:), allocatable :: dummy_num
  character (len=80) :: string1,string2,string3,string4,string5,string6,string7,string8,string9,string10,string11
  integer :: idummy1,idummy2,idummy3,idummy4,idummy5

  integer :: iglob
  real, dimension(:),allocatable :: flag_val
  logical, dimension(:),allocatable :: mask_ibool

  !MPI variables
  integer :: ier

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_vp, model_vs, model_vp_new, model_vs_new
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_dA, model_dB, total_model
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_a, kernel_b, kernel_rho
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_rho, model_dR
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_rho_new
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_vp_rel, model_vs_rel, model_rho_rel

  ! to update the values in external_mesh.bin
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kappastore_new,mustore_new
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp_new,rho_vs_new
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore_new

  ! steepest descent lengths
  real(kind=CUSTOM_REAL) :: step_fac,step_length
  character(len=150) :: s_step_fac

  ! thresholds
  real(kind=CUSTOM_REAL) :: VS_MIN, VS_MAX, VP_MIN, VP_MAX
!   real(kind=CUSTOM_REAL) :: RHO_MIN, RHO_MAX

  ! statistics:
  !models
  real(kind=CUSTOM_REAL) :: vsmin_before, vsmax_before, vpmin_before, vpmax_before
  real(kind=CUSTOM_REAL) :: vsmin_after, vsmax_after, vpmin_after, vpmax_after
  real(kind=CUSTOM_REAL) :: vsmin_new_before, vsmax_new_before, vpmin_new_before, vpmax_new_before
  real(kind=CUSTOM_REAL) :: vsmin_new_after, vsmax_new_after, vpmin_new_after, vpmax_new_after
  real(kind=CUSTOM_REAL) :: rhomin_before, rhomax_before !, rhomin_after, rhomax_after
  real(kind=CUSTOM_REAL) :: rhomin_new_before, rhomax_new_before, rhomin_new_after, rhomax_new_after
  !kernels
  real(kind=CUSTOM_REAL) :: min_vs_k, max_vs_k, min_vp_k, max_vp_k, min_rho_k, max_rho_k
  !gradients in direction of steepest descent (= - kernels)
  real(kind=CUSTOM_REAL) :: min_vs_g, max_vs_g, min_vp_g, max_vp_g, min_rho_g, max_rho_g
  real(kind=CUSTOM_REAL) :: min_vp,min_vs,max_vp,max_vs,min_rho,max_rho, &
                            max,minmax(4),vs_sum,vp_sum,rho_sum

  integer :: i,j,k,ispec

  ! calculate rmass
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: jacobianstore
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl
  ! mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_old, rmass_new
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_acoustic_new,rmass_solid_poroelastic_new,rmass_fluid_poroelastic_new



  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)


  ! subjective step length to multiply to the gradient
  ! e.g. step_fac = 0.03

  call get_command_argument(1,s_step_fac)

  if (trim(s_step_fac) == '') then
    call exit_MPI(myrank,'Usage: add_model step_factor')
  endif

  ! read in parameter information
  read(s_step_fac,*) step_fac
  if( abs(step_fac) < 1.e-10) then
    print*,'error: step factor ',step_fac
    call exit_MPI(myrank,'error step factor')
  endif

  call sync_all()

  ! reads in parameters and checks for some inconsistencies
  call initialize_simulation()

  call sync_all()

  ! reads in external mesh
  call read_mesh_databases()

  !===================================================
  !===================================================
  ! MODEL UPDATE

  ! number of spectral el for each processor
  NSPEC=NSPEC_AB
  NGLOB=NGLOB_AB

  if( myrank == 0 ) then
    print*,'NSPEC            ', NSPEC
    print*,'NGLOB            ', NGLOB
  endif

  call sync_all()

  !! allocation
  ! model and kernel variables
  allocate(model_vp(NGLLX,NGLLY,NGLLZ,NSPEC), &
          model_vs(NGLLX,NGLLY,NGLLZ,NSPEC), &
          model_vp_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
          model_vs_new(NGLLX,NGLLY,NGLLZ,NSPEC))
  allocate(model_dA(NGLLX,NGLLY,NGLLZ,NSPEC), &
          model_dB(NGLLX,NGLLY,NGLLZ,NSPEC), &
          total_model(NGLLX,NGLLY,NGLLZ,NSPEC))
  allocate(kernel_a(NGLLX,NGLLY,NGLLZ,NSPEC), &
          kernel_b(NGLLX,NGLLY,NGLLZ,NSPEC))
  allocate(kernel_rho(NGLLX,NGLLY,NGLLZ,NSPEC))
  allocate(model_rho(NGLLX,NGLLY,NGLLZ,NSPEC), &
          model_rho_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
          model_dR(NGLLX,NGLLY,NGLLZ,NSPEC))
  ! new variables for save_external_bin_m_up
  allocate(kappastore_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
          mustore_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
          rhostore_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
          rho_vp_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
          rho_vs_new(NGLLX,NGLLY,NGLLZ,NSPEC))


  !! inizialize variables
  ! old model
  model_vp = 0._CUSTOM_REAL
  model_vs = 0._CUSTOM_REAL
  model_rho = 0._CUSTOM_REAL
  model_dA = 0._CUSTOM_REAL
  model_dB = 0._CUSTOM_REAL
  model_dR = 0._CUSTOM_REAL
  kernel_a = 0._CUSTOM_REAL
  kernel_b = 0._CUSTOM_REAL
  kernel_rho = 0._CUSTOM_REAL


  if (myrank == 0) then
    print*,'defaults'
    print*
    print*,'model update for vs & vp & rho'
    print*,'  step_fac = ',step_fac
    print*
    if( USE_ALPHA_BETA_RHO ) then
      print*,'kernel parameterization: (alpha,beta,rho)'
    else
      print*,'kernel parameterization: (bulk,beta,rho)'
    endif
    print*
    if( USE_RHO_SCALING ) then
      print*,'scaling rho perturbations'
      print*
    endif

    if( PRINT_OUT_FILES ) then
     open(20,file=trim(OUTPUT_MODEL_UPD)//'/step_fac',status='unknown')
     write(20,'(1e24.12)') step_fac
     close(20)
    endif
  endif


  !! set up thresholds for old and new models

  if(MINMAX_THRESHOLD_OLD .or. MINMAX_THRESHOLD_NEW) then
     ! minmax wavespeed values for southern california simulations
     VS_MIN = 600.0
     VS_MAX = 4700.0
     VP_MIN = 1500.0
     VP_MAX = 8200.0
!      RHO_MIN =     !remove the comment if we are interested in rho
!      RHO_MAX =     !remove the comment if we are interested in rho


    if( PRINT_OUT_FILES ) then
     if (myrank == 0) then
        open(19,file=trim(OUTPUT_MODEL_UPD)//'/VS_VP_MINMAX',status='unknown')
        write(19,'(4e24.12)') VS_MIN, VS_MAX, VP_MIN, VP_MAX
        close(19)
     endif
    !or
!     if (myrank == 0) then
!        open(19,file=trim(OUTPUT_MODEL_UPD)//'/VS_VP_RHO_MINMAX',status='unknown')
!        write(19,'(4e24.12)') VS_MIN, VS_MAX, VP_MIN, VP_MAX, RHO_MIN, RHO_MAX
!        close(19)
!     endif
    endif
     if( myrank == 0 ) then
       print*,'thresholds:'
       print*,'  VS_MIN, VS_MAX: ',VS_MIN, VS_MAX
       print*,'  VP_MIN, VP_MAX: ',VP_MIN, VP_MAX
       !or
!        print*,'  RHO_MIN, RHO_MAX: ',RHO_MIN, RHO_MAX
       print*
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
  call mpi_reduce(minval(model_vs(:,:,:,1:nspec)), vsmin_before, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vs(:,:,:,1:nspec)), vsmax_before, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vp(:,:,:,1:nspec)), vpmin_before, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vp(:,:,:,1:nspec)), vpmax_before, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_rho(:,:,:,1:nspec)), rhomin_before,1,CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_rho(:,:,:,1:nspec)), rhomax_before,1,CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)


  if( PRINT_OUT_FILES ) then
!   if (myrank == 0) then
!      open(19,file=trim(OUTPUT_MODEL_UPD)//'/vs_vp_minmax_before',status='unknown')
!      write(19,'(4e24.12)') vsmin_before, vsmax_before, vpmin_before, vpmax_before
!      close(19)
!   endif
  !or
   if (myrank == 0) then
      open(19,file=trim(OUTPUT_MODEL_UPD)//'/vs_vp_rho_minmax_before',status='unknown')
      write(19,'(4e24.12)') vsmin_before, vsmax_before, vpmin_before, vpmax_before, rhomin_before, rhomax_before
      close(19)
   endif
  endif

  if( myrank == 0 ) then
    print*,'models before threshold:'
    print*,'  vs min/max: ',vsmin_before, vsmax_before
    print*,'  vp min/max: ',vpmin_before, vpmax_before
    print*,'  rho min/max: ',rhomin_before, rhomax_before
    print*
  endif


  ! threshold current model and write out the modified version
  if(MINMAX_THRESHOLD_OLD) then
     do ispec=1,NSPEC
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 if(model_vs(i,j,k,ispec) < VS_MIN) model_vs(i,j,k,ispec) = VS_MIN
                 if(model_vs(i,j,k,ispec) > VS_MAX) model_vs(i,j,k,ispec) = VS_MAX
                 if(model_vp(i,j,k,ispec) < VP_MIN) model_vp(i,j,k,ispec) = VP_MIN
                 if(model_vp(i,j,k,ispec) > VP_MAX) model_vp(i,j,k,ispec) = VP_MAX
!                  if(model_rho(i,j,k,ispec) < RHO_MIN) model_rho(i,j,k,ispec) = RHO_MIN
!                  if(model_rho(i,j,k,ispec) > RHO_MAX) model_rho(i,j,k,ispec) = RHO_MAX
              enddo
           enddo
        enddo
     enddo

     ! compute minmax values of the thresholded current model
     call mpi_reduce(minval(model_vs(:,:,:,1:nspec)), vsmin_after, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
     call mpi_reduce(maxval(model_vs(:,:,:,1:nspec)), vsmax_after, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
     call mpi_reduce(minval(model_vp(:,:,:,1:nspec)), vpmin_after, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
     call mpi_reduce(maxval(model_vp(:,:,:,1:nspec)), vpmax_after, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
!      call mpi_reduce(minval(model_rho(:,:,:,1:nspec)), rhomin_after,1,CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
!      call mpi_reduce(maxval(model_rho(:,:,:,1:nspec)), rhomax_after,1,CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)


    if( PRINT_OUT_FILES ) then
     if (myrank == 0) then
         open(19,file=trim(OUTPUT_MODEL_UPD)//'/vs_vp_minmax_after',status='unknown')
         write(19,'(4e24.12)') vsmin_after, vsmax_after, vpmin_after, vpmax_after
         close(19)
      endif
!or
!       if (myrank == 0) then
!          open(19,file=trim(OUTPUT_MODEL_UPD)//'/vs_vp_rho_minmax_after',status='unknown')
!          write(19,'(4e24.12)') vsmin_after, vsmax_after, vpmin_after, vpmax_after, rhomin_after, rhomax_after
!          close(19)
!       endif
    endif

     if( myrank == 0 ) then
       print*,'models after threshold:'
       print*,'  vs min/max: ',vsmin_after, vsmax_after
       print*,'  vp min/max: ',vpmin_after, vpmax_after
   !     print*,'  rho min/max: ',rhomin_after, rhomax_after
       print*
     endif


     fname = 'vs'
     write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//trim(fname)//'.bin'
     open(12,file=trim(m_file),form='unformatted')
     write(12) model_vs(:,:,:,1:nspec)
     close(12)

     fname = 'vp'
     write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//trim(fname)//'.bin'
     open(12,file=trim(m_file),form='unformatted')
     write(12) model_vp(:,:,:,1:nspec)
     close(12)

!      fname = 'rho'
!      write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//trim(fname)//'.bin'
!      open(12,file=trim(m_file),form='unformatted')
!      write(12) model_rho(:,:,:,1:nspec)
!      close(12)

  endif


  !---------------------------------------------------------------------------------------------
  ! reads in smoothed (& summed) event kernel
  !---------------------------------------------------------------------------------------------

  ! alpha kernel
  if( USE_ALPHA_BETA_RHO ) then
    ! reads in alpha kernel
    fname = 'alpha_kernel_smooth'
  else
    ! reads in bulk_c kernel
    fname = 'bulk_c_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS)//'/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_a(:,:,:,1:nspec)
  close(12)

  ! beta kernel
  if( USE_ALPHA_BETA_RHO ) then
    ! reads in beta kernel
    fname = 'beta_kernel_smooth'
  else
    ! reads in bulk_beta kernel
    fname = 'bulk_beta_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS)//'/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_b(:,:,:,1:nspec)
  close(12)

  ! rho kernel
  if( USE_RHO_SCALING ) then

    ! uses scaling relation with shear perturbations
    kernel_rho(:,:,:,:) = RHO_SCALING * kernel_b(:,:,:,:)

  else

    ! uses rho kernel
    write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS)//'/proc',myrank,'_rho_kernel_smooth.bin'
    open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(12) kernel_rho(:,:,:,1:nspec)
    close(12)
  endif

  ! statistics
  call mpi_reduce(minval(kernel_a),min_vp_k,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_a),max_vp_k,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(kernel_b),min_vs_k,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_b),max_vs_k,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(kernel_rho),min_rho_k,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_rho),max_rho_k,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)


  if( PRINT_OUT_FILES ) then
   if (myrank == 0) then
      open(19,file=trim(OUTPUT_MODEL_UPD)//'/kernels_minmax',status='unknown')
      write(19,'(4e24.12)') min_vs_k, max_vs_k, min_vp_k, max_vp_k, min_rho_k, max_rho_k
      close(19)
   endif
  endif

  if( myrank == 0 ) then
    print*,'initial kernels:'
    print*,'  beta min/max: ',min_vs_k,max_vs_k
    print*,'  alpha min/max: ',min_vp_k,max_vp_k
    print*,'  rho min/max: ',min_rho_k,max_rho_k
    print*
  endif


  !---------------------------------------------------------------------------------------------
  ! calculates gradient
  !---------------------------------------------------------------------------------------------

  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

            ! gradient in direction of steepest descent

            ! for vp
            model_dA(i,j,k,ispec) = - kernel_a(i,j,k,ispec)

            ! for shear
            model_dB(i,j,k,ispec) = - kernel_b(i,j,k,ispec)

            ! for rho
            model_dR(i,j,k,ispec) = - kernel_rho(i,j,k,ispec)

        enddo
      enddo
    enddo
  enddo

  ! statistics
  call mpi_reduce(minval(model_dA),min_vp_g,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dA),max_vp_g,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_dB),min_vs_g,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dB),max_vs_g,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_dR),min_rho_g,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dR),max_rho_g,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)


  if( PRINT_OUT_FILES ) then
   if (myrank == 0) then
      open(19,file=trim(OUTPUT_MODEL_UPD)//'/gradients_minmax',status='unknown')
      write(19,'(4e24.12)') min_vs_g, max_vs_g, min_vp_g, max_vp_g, min_rho_g, max_rho_g
      close(19)
   endif
  endif

  if( myrank == 0 ) then
    print*,'initial gradients:'
    print*,'  beta min/max : ',min_vs_g,max_vs_g
    print*,'  a min/max: ',min_vp_g,max_vp_g
    print*,'  rho min/max: ',min_rho_g,max_rho_g
    print*
  endif


  !---------------------------------------------------------------------------------------------
  ! master determines step length based on maximum gradient value (either vp or vs)
  !---------------------------------------------------------------------------------------------

  if( myrank == 0 ) then
    minmax(1) = abs(min_vs_g)
    minmax(2) = abs(max_vs_g)

    minmax(3) = abs(min_vp_g)
    minmax(4) = abs(max_vp_g)

    max = maxval(minmax)
    step_length = step_fac/max

    print*,'  step length : ',step_length,max
    print*
  endif
  call mpi_bcast(step_length,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)


  !---------------------------------------------------------------------------------------------
  ! gradient length
  !---------------------------------------------------------------------------------------------

  max_vp = sum( model_dA * model_dA )
  max_vs = sum( model_dB * model_dB )
  max_rho = sum( model_dR * model_dR )
  max_vp = sqrt(max_vp)
  max_vs = sqrt(max_vs)
  max_rho = sqrt(max_rho)

  ! statistics
  call mpi_reduce(max_vp,vp_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(max_vs,vs_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(max_rho,rho_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)


  if( PRINT_OUT_FILES ) then
   if (myrank == 0) then
      open(19,file=trim(OUTPUT_MODEL_UPD)//'/vs_vp_rho_sum',status='unknown')
      write(19,'(4e24.12)') vs_sum, vp_sum, rho_sum
      close(19)
   endif
  endif

  if( myrank == 0 ) then
    print*,'  initial beta length : ',vs_sum
    print*,'  initial a length: ',vp_sum
    print*,'  initial rho length: ',rho_sum
    print*
  endif


  !---------------------------------------------------------------------------------------------
  ! model updates
  !---------------------------------------------------------------------------------------------

  ! multiply model updates by a subjective factor that will change the step

  model_dA(:,:,:,:) = step_length * model_dA(:,:,:,:)
  model_dB(:,:,:,:) = step_length * model_dB(:,:,:,:)
  model_dR(:,:,:,:) = step_length * model_dR(:,:,:,:)


  ! statistics
  call mpi_reduce(minval(model_dA),min_vp,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dA),max_vp,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_dB),min_vs,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dB),max_vs,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_dR),min_rho,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dR),max_rho,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)


  if( PRINT_OUT_FILES ) then
   if (myrank == 0) then
      open(19,file=trim(OUTPUT_MODEL_UPD)//'/scaled_gradients',status='unknown')
      write(19,'(4e24.12)') min_vs,max_vs,min_vp,max_vp,min_rho,max_rho
      close(19)
   endif
  endif

  if( myrank == 0 ) then
    print*,'scaled gradients:'
    print*,'  beta min/max : ',min_vs,max_vs
    print*,'  a min/max: ',min_vp,max_vp
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif


  !---------------------------------------------------------------------------------------------
  ! new model
  !---------------------------------------------------------------------------------------------

  ! computes new model values for alpha, beta and rho

  ! S wavespeed model
  model_vs_new = 0._CUSTOM_REAL
  model_vs_new = model_vs * exp( model_dB )

  ! P wavespeed model
  model_vp_new = 0._CUSTOM_REAL
  if( USE_ALPHA_BETA_RHO ) then
    ! new vp values use alpha model update
    model_vp_new = model_vp * exp( model_dA )
  else
    ! new vp values use bulk model update:
    ! this is based on vp_new = sqrt( bulk_new**2 + 4/3 vs_new**2 )
    model_vp_new = sqrt( model_vp**2 * exp(2.0*model_dA) + FOUR_THIRDS * model_vs**2 *( &
                            exp(2.0*model_dB) - exp(2.0*model_dA) ) )
  endif

  ! Rho density model
  model_rho_new = 0._CUSTOM_REAL
  model_rho_new = model_rho * exp( model_dR )

  ! statistcs
  ! compute minmax values of new model
  call mpi_reduce(minval(model_vs_new(:,:,:,1:nspec)), vsmin_new_before, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vs_new(:,:,:,1:nspec)), vsmax_new_before, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vp_new(:,:,:,1:nspec)), vpmin_new_before, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vp_new(:,:,:,1:nspec)), vpmax_new_before, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_rho_new(:,:,:,1:nspec)),rhomin_new_before,1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_rho_new(:,:,:,1:nspec)),rhomax_new_before,1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)


  if( PRINT_OUT_FILES ) then
!    if (myrank == 0) then
!       open(19,file=trim(OUTPUT_MODEL_UPD)//'/vs_vp_new_minmax_before',status='unknown')
!       write(19,'(4e24.12)') vsmin_new_before, vsmax_new_before, vpmin_new_before, vpmax_new_before
!       close(19)
!    endif
    !or
   if (myrank == 0) then
      open(19,file=trim(OUTPUT_MODEL_UPD)//'/vs_vp_rho_new_minmax_before',status='unknown')
      write(19,'(4e24.12)') vsmin_new_before, vsmax_new_before, &
                           vpmin_new_before, vpmax_new_before, &
                           rhomin_new_before, rhomax_new_before
      close(19)
   endif
  endif

  if( myrank == 0 ) then
    print*,'new models before threshold:'
    print*,'  vs min/max: ',vsmin_new_before, vsmax_new_before
    print*,'  vp min/max: ',vpmin_new_before, vpmax_new_before
    print*,'  rho min/max: ',rhomin_new_before, rhomax_new_before
    print*
  endif


  ! threshold model according to minmax values specified above
  if(MINMAX_THRESHOLD_NEW) then
     do ispec=1,NSPEC
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 if(model_vs_new(i,j,k,ispec) < VS_MIN) model_vs_new(i,j,k,ispec) = VS_MIN
                 if(model_vs_new(i,j,k,ispec) > VS_MAX) model_vs_new(i,j,k,ispec) = VS_MAX
                 if(model_vp_new(i,j,k,ispec) < VP_MIN) model_vp_new(i,j,k,ispec) = VP_MIN
                 if(model_vp_new(i,j,k,ispec) > VP_MAX) model_vp_new(i,j,k,ispec) = VP_MAX
!                  if(model_rho_new(i,j,k,ispec) < RHO_MIN) model_rho_new(i,j,k,ispec) = RHO_MIN
!                  if(model_rho_new(i,j,k,ispec) > RHO_MAX) model_rho_new(i,j,k,ispec) = RHO_MAX
              enddo
           enddo
        enddo
     enddo
  endif


  ! write out new models and their global minmax values
  ! compute minmax values of new model
  call mpi_reduce(minval(model_vs_new(:,:,:,1:nspec)), vsmin_new_after, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vs_new(:,:,:,1:nspec)), vsmax_new_after, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vp_new(:,:,:,1:nspec)), vpmin_new_after, 1, CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vp_new(:,:,:,1:nspec)), vpmax_new_after, 1, CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_rho_new(:,:,:,1:nspec)), rhomin_new_after,1,CUSTOM_MPI_TYPE, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_rho_new(:,:,:,1:nspec)), rhomax_new_after,1,CUSTOM_MPI_TYPE, MPI_MAX, 0, MPI_COMM_WORLD,ier)


  ! this should only be different if using MINMAX_THRESHOLD_NEW
  if( PRINT_OUT_FILES ) then
!    if (myrank == 0) then
!       open(19,file=trim(OUTPUT_MODEL_UPD)//'/vs_vp_new_minmax_after',status='unknown')
!       write(19,'(4e24.12)') vsmin_new_after, vsmax_new_after, vpmin_new_after, vpmax_new_after
!       close(19)
!    endif
   !or
   if (myrank == 0) then
      open(19,file=trim(OUTPUT_MODEL_UPD)//'/vs_vp_rho_new_minmax_after',status='unknown')
      write(19,'(4e24.12)') vsmin_new_after, vsmax_new_after, vpmin_new_after, vpmax_new_after, rhomin_new_after, rhomax_new_after
      close(19)
   endif
  endif

  if( myrank == 0 ) then
    print*,'new models after threshold:'
    print*,'  vs min/max: ',vsmin_new_after, vsmax_new_after
    print*,'  vp min/max: ',vpmin_new_after, vpmax_new_after
    print*,'  rho min/max: ',rhomin_new_after, rhomax_new_after
    print*
  endif


  !---------------------------------------------------------------------------------------------
  ! write NEW external_mesh.bin
  ! and store NEW model files vp_new.bin vs_new.bin rho_new.bin vp_new.vtk vs_new.vtk rho_new.vtk
  ! calling save_external_bin_m_up with SAVE_MESH_FILES=true
  ! and also write NEW attenuation files attenuation.bin and attenuation.vtk (this should be equal to the old one)
  !---------------------------------------------------------------------------------------------

  if( myrank == 0 ) then
    print*
    print*,'  ...saving new databases'
    print*
  endif


  NGLOB_OCEAN = NGLOB


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
  allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,NSPEC))
  jacobianstore = 0._CUSTOM_REAL
  jacobianstore = jacobian

  ! set up coordinates of the Gauss-Lobatto-Legendre points and calculate weights
  ! from constants.h GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! rmass for the OLD model from read_mesh_databases
  allocate(rmass_old(NGLOB))
  rmass_old = 0._CUSTOM_REAL
  rmass_old = rmass

  call sync_all()

  ! create mass matrix ONLY for the elastic case
  allocate(rmass_new(NGLOB))

  call sync_all()

  if( myrank == 0) then
    print*, '  ...creating mass matrix '
  endif

  rmass_new(:) = 0._CUSTOM_REAL

  ! note: this does not update the absorbing boundary contributions to the mass matrix
  ! elastic mass matrix
  do ispec=1,nspec
    if( ispec_is_elastic(ispec) ) then
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)

            weight = wxgll(i)*wygll(j)*wzgll(k)
            jacobianl = jacobianstore(i,j,k,ispec)

            if( myrank == 0) then
              print*, 'weight', weight
              print*, 'jacobianl', jacobianl
            endif

            if(CUSTOM_REAL == SIZE_REAL) then
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


  call sync_all()

  ! dummy allocations, arrays are not needed since the update here only works for elastic models
  allocate(rmass_acoustic_new(NGLOB))
  allocate(rmass_solid_poroelastic_new(NGLOB))
  allocate(rmass_fluid_poroelastic_new(NGLOB))
  rmass_acoustic_new = 0._CUSTOM_REAL
  rmass_solid_poroelastic_new = 0._CUSTOM_REAL
  rmass_fluid_poroelastic_new = 0._CUSTOM_REAL

  !-------- attenuation -------

  ! read the proc*attenuation.vtk for the old model in LOCAL_PATH and store qmu_attenuation_store

  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_attenuation.vtk'
  open(12,file=trim(m_file),status='old',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12,'(a)') string1 !text
  read(12,'(a)') string2 !text
  read(12,'(a)') string3 !text
  read(12,'(a)') string4 !text
  read(12,'(a,i12,a)') string5, idummy1, string6 !text
  allocate(dummy_g_1(NGLOB),dummy_g_2(NGLOB),dummy_g_3(NGLOB),stat=ier)
  if( ier /= 0 ) stop 'error allocating array dummy etc.'
  read(12,'(3e18.6)') dummy_g_1,dummy_g_2,dummy_g_3 !xstore,ystore,zstore for i=1,nglob
  read(12,*) !blank line
  deallocate(dummy_g_1,dummy_g_2,dummy_g_3)

  read(12,'(a,i12,i12)') string7, idummy2, idummy3 !text
  allocate(dummy_num(NSPEC),dummy_l_1(NSPEC),dummy_l_2(NSPEC),dummy_l_3(NSPEC),dummy_l_4(NSPEC), &
           dummy_l_5(NSPEC),dummy_l_6(NSPEC),dummy_l_7(NSPEC),dummy_l_8(NSPEC),stat=ier)
  if( ier /= 0 ) stop 'error allocating array dummy etc.'
  read(12,'(9i12)') dummy_num,dummy_l_1,dummy_l_2,dummy_l_3,dummy_l_4, &
                  dummy_l_5,dummy_l_6,dummy_l_7,dummy_l_8 !8,ibool-1 for ispec=1,nspec
  read(12,*) !blank line
  deallocate(dummy_num,dummy_l_1,dummy_l_2,dummy_l_3,dummy_l_4,dummy_l_5,dummy_l_6,dummy_l_7,dummy_l_8)

  read(12,'(a,i12)') string8, idummy4 !text
  allocate(dummy_num(NSPEC),stat=ier)
  if( ier /= 0 ) stop 'error allocating array dummy etc.'
  read(12,*) dummy_num !12 for ispec=1,nspec
  read(12,*) !blank line
  deallocate(dummy_num)

  read(12,'(a,i12)') string9, idummy5 !text
  read(12,'(a)') string10 !text
  read(12,'(a)') string11 !text
  allocate(flag_val(NGLOB),stat=ier)
  if( ier /= 0 ) stop 'error allocating flag_val'
  read(12,*) flag_val
  read(12,*) !blank line
  close(12)

  ! store the attenuation flag in qmu_attenuation_store
  allocate(qmu_attenuation_store(NGLLX,NGLLY,NGLLZ,NSPEC))
  qmu_attenuation_store=0._CUSTOM_REAL

  allocate(mask_ibool(NGLOB),stat=ier)
  if( ier /= 0 ) stop 'error allocating mask'

  mask_ibool = .false.
  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          if( .not. mask_ibool(iglob) ) then
            qmu_attenuation_store(i,j,k,ispec) = flag_val(iglob)
            mask_ibool(iglob) = .true.
          endif
        enddo
      enddo
    enddo
  enddo


  call sync_all()

  ! calculate min_resolved_period needed for attenuation model
  call check_mesh_resolution(myrank,NSPEC,NGLOB,ibool,&
                            xstore,ystore,zstore, &
                            kappastore_new,mustore_new,rho_vp_new,rho_vs_new, &
                            -1.0d0, model_speed_max,min_resolved_period, &
                            LOCAL_PATH,SAVE_MESH_FILES )

  call sync_all()

  ! saves binary mesh files for attenuation for the NEW model
  call create_name_database(prname_new,myrank,LOCAL_PATH_NEW)
  if( myrank == 0 ) then
    print*,'attenuation.bin new: ', prname_new
  endif

  if( ATTENUATION ) then
    call get_attenuation_model(myrank,NSPEC,USE_OLSEN_ATTENUATION,OLSEN_ATTENUATION_RATIO, &
                          mustore_new,rho_vs_new,kappastore_new,rho_vp_new,qmu_attenuation_store, &
                          ispec_is_elastic,min_resolved_period,prname_new,FULL_ATTENUATION_SOLID)
  endif

  !----------------------------

  call sync_all()

  if( myrank == 0 ) then
    print*,'external_mesh.bin new: ', prname_new
  endif

  call save_external_bin_m_up(NSPEC,NGLOB, &
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




  !---------------------------------------------------------------------------------------------
  ! ALTERNATIVELY:
  ! SAVE_MESH_FILES=false
  ! directly stores new model files vp.bin vs.bin rho.bin vp.vtk vs.vtk rho.vtk
  ! NOTE: save_external_bin_m_up will create NEW external_mesh.bin anyway
  !---------------------------------------------------------------------------------------------

!   ! vs model
!   fname = 'vs_new'
!   write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH_NEW)//'/proc',myrank,'_'//trim(fname)//'.bin'
!   open(12,file=trim(m_file),form='unformatted')
!   write(12) model_vs_new(:,:,:,1:nspec)
!   close(12)
!
!   m_file = trim(LOCAL_PATH_NEW)//'/proc',myrank,'_'//trim(fname)
!   call write_VTK_data_gll_cr(nspec,nglob, &
!                       xstore,ystore,zstore,ibool, &
!                       model_vs_new,m_file)
!
!
!   ! vp model
!   fname = 'vp_new'
!   write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH_NEW)//'/proc',myrank,'_'//trim(fname)//'.bin'
!   open(12,file=trim(m_file),form='unformatted')
!   write(12) model_vp_new(:,:,:,1:nspec)
!   close(12)
!
!   m_file = trim(LOCAL_PATH_NEW)//'/proc',myrank,'_'//trim(fname)
!   call write_VTK_data_gll_cr(nspec,nglob, &
!                       xstore,ystore,zstore,ibool, &
!                       model_vp_new,m_file)
!
!
!   ! rho model
!   !remove the comment if we are interested in rho
!   fname = 'rho_new'
!   write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH_NEW)//'/proc',myrank,'_'//trim(fname)//'.bin'
!   open(12,file=trim(m_file),form='unformatted')
!   write(12) model_rho_new(:,:,:,1:nspec)
!   close(12)

!   m_file = trim(LOCAL_PATH_NEW)//'/proc',myrank,'_'//trim(fname)
!   call write_VTK_data_gll_cr(nspec,nglob, &
!                       xstore,ystore,zstore,ibool, &
!                       model_rho_new,m_file)


  !---------------------------------------------------------------------------------------------
  ! P & S & rho model update
  ! stores the linear approximations of the model updates
  !---------------------------------------------------------------------------------------------

  allocate(model_vp_rel(NGLLX,NGLLY,NGLLZ,NSPEC), &
          model_vs_rel(NGLLX,NGLLY,NGLLZ,NSPEC), &
          model_rho_rel(NGLLX,NGLLY,NGLLZ,NSPEC))

  ! stores relative Vp model perturbations
  where( model_vp /= 0.0 ) model_vp_rel = ( model_vp_new - model_vp) / model_vp

  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH_NEW)//'/proc',myrank,'_dvpvp.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vp_rel(:,:,:,1:nspec)
  close(12)
  call mpi_reduce(maxval(model_vp_rel),max_vp,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vp_rel),min_vp,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! relative S model perturbations
  where ( model_vs > 1.e-10 ) model_vs_rel = ( model_vs_new - model_vs) / model_vs
  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH_NEW)//'/proc',myrank,'_dvsvs.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vs_rel(:,:,:,1:nspec)
  close(12)
  call mpi_reduce(maxval(model_vs_rel),max_vs,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vs_rel),min_vs,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! relative rho model perturbations
!   where ( model_rho > 1.e-10 ) model_rho_rel = ( model_rho_new - model_rho) / model_rho
!   write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH_NEW)//'/proc',myrank,'_drhorho.bin'
!   open(12,file=trim(m_file),form='unformatted',action='write')
!   write(12) model_rho_rel(:,:,:,1:nspec)
!   close(12)
!   call mpi_reduce(maxval(model_rho_rel),max_rho,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
!   call mpi_reduce(minval(model_rho_rel),min_rho,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)


  if( PRINT_OUT_FILES ) then
   if (myrank == 0) then
      open(19,file=trim(OUTPUT_MODEL_UPD)//'/relative_vs_vp',status='unknown')
      write(19,'(4e24.12)') min_vs,max_vs,min_vp,max_vp
      close(19)
   endif
   !or
!    if (myrank == 0) then
!       open(19,file=trim(OUTPUT_MODEL_UPD)//'/relative_vs_vp_rho',status='unknown')
!       write(19,'(4e24.12)') min_vs,max_vs,min_vp,max_vp,min_rho,max_rho
!       close(19)
!    endif
  endif

  if( myrank == 0 ) then
    print*,'relative Vs & Vp update:'
    print*,'  dvs/vs min/max: ',min_vs,max_vs
    print*,'  dvp/vp min/max: ',min_vp,max_vp
!     print*,'  drho/rho min/max: ',min_rho,max_rho
    print*
  endif

  deallocate(model_vp_rel, model_vs_rel, model_rho_rel)

  !---------------------------------------------------------------------------------------------
  ! compute Poisson's ratio of the current model and new model
  !---------------------------------------------------------------------------------------------

  ! Poisson's ratio of current model
  total_model = 0.
  total_model = ( model_vp**2 - 2.0*model_vs**2 ) / ( 2.0*model_vp**2 - 2.0*model_vs**2 )
  fname = 'poisson'
  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

  ! Poisson's ratio of new model
  total_model = 0.
  total_model = ( model_vp_new**2 - 2.0*model_vs_new**2 ) / &
                ( 2.0*model_vp_new**2 - 2.0*model_vs_new**2 )
  fname = 'poisson_new'
  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH_NEW)//'/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)


  !---------------------------------------------------------------------------------------------
  ! compute bulk wavespeed of the current model and new model
  !---------------------------------------------------------------------------------------------

  ! bulk wavespeed of current model
  total_model = 0.
  total_model = sqrt( model_vp**2 - (4.0/3.0)*model_vs**2 )
  fname = 'vb'
  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

  ! bulk wavespeed of new model
  total_model = 0.
  total_model = sqrt( model_vp_new**2 - (4.0/3.0)*model_vs_new**2 )
  fname = 'vb_new'
  write(m_file,'(a,i6.6,a)') trim(LOCAL_PATH_NEW)//'/proc',myrank,'_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted')
  write(12) total_model(:,:,:,1:nspec)
  close(12)

!===================================================
!===================================================


  deallocate(qmu_attenuation_store)
  deallocate(flag_val,mask_ibool)

  deallocate(model_vp, model_vs, model_vp_new, model_vs_new, &
             model_dA, model_dB, total_model, &
             kernel_a, kernel_b)
  deallocate(kernel_rho)
  deallocate(model_rho, model_rho_new, model_dR)
  deallocate(rhostore_new, kappastore_new, mustore_new, &
             rho_vp_new, rho_vs_new)

  deallocate(jacobianstore)
  deallocate(rmass_old,rmass_new)
  deallocate(rmass_acoustic_new,rmass_solid_poroelastic_new,rmass_fluid_poroelastic_new)


  !-----------------------------------------------------

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program model_update
