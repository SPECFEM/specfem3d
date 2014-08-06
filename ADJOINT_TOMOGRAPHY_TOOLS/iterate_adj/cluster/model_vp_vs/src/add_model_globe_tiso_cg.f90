! add_model_globe_tiso_cg
!
! this program can be used to update TRANSVERSE ISOTROPIC model files
! based on smoothed event kernels.
! the kernels are given for tranverse isotropic parameters (bulk_c,bulk_betav,bulk_betah,eta).
!
! the algorithm uses a conjugate gradient method with a step length
! limited by the given maximum update percentage.
!
! input:
!    - step_fac : step length to update the models, f.e. 0.03 for plusminus 3%
!
! setup:
!
!- INPUT_MODEL/  contains:
!       proc000***_reg1_vsv.bin &
!       proc000***_reg1_vsh.bin &
!       proc000***_reg1_vpv.bin &
!       proc000***_reg1_vph.bin &
!       proc000***_reg1_eta.bin &
!       proc000***_reg1_rho.bin
!
!- INPUT_GRADIENT/ contains:
!       proc000***_reg1_bulk_c_kernel_smooth.bin &
!       proc000***_reg1_bulk_betav_kernel_smooth.bin &
!       proc000***_reg1_bulk_betah_kernel_smooth.bin &
!       proc000***_reg1_eta_kernel_smooth.bin &
!       proc000***_reg1_rho_kernel_smooth.bin
!
!- /tigress-hsm/dpeter/SPECFEM3D_GLOBE/KERNELS/OUTPUT_SUM.old/ contains old gradients:
!       proc000***_reg1_bulk_c_kernel_smooth.bin &
!       proc000***_reg1_bulk_betav_kernel_smooth.bin &
!       proc000***_reg1_bulk_betah_kernel_smooth.bin &
!       proc000***_reg1_eta_kernel_smooth.bin &
!       proc000***_reg1_rho_kernel_smooth.bin
!
!- topo/ contains:
!       proc000***_reg1_solver_data_1.bin
!
! new models are stored in
!- OUTPUT_MODEL/ as
!   proc000***_reg1_vpv_new.bin &
!   proc000***_reg1_vph_new.bin &
!   proc000***_reg1_vsv_new.bin &
!   proc000***_reg1_vsh_new.bin &
!   proc000***_reg1_eta_new.bin &
!   proc000***_reg1_rho_new.bin
!
! usage: ./add_model_globe_tiso_cg 0.3
!
!
! NOTE: this routine uses previous model update in OUTPUT_SUM.old/
!             for a conjugate gradient update
!

module model_update_cg

  include 'mpif.h'
  include 'constants_globe.h'
  include 'precision_globe.h'
  include 'values_from_mesher_globe.h'

  ! ======================================================

  ! density scaling factor with shear perturbations
  ! see e.g. Montagner & Anderson (1989), Panning & Romanowicz (2006)
  real(kind=CUSTOM_REAL),parameter :: RHO_SCALING = 0.33_CUSTOM_REAL

  ! constraint on eta model
  real(kind=CUSTOM_REAL),parameter :: LIMIT_ETA_MIN = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: LIMIT_ETA_MAX = 1.5_CUSTOM_REAL

  ! directory which contains kernels/gradients from former iteration
  character(len=150),parameter :: &
    kernel_old_dir = '/tigress-hsm/dpeter/SPECFEM3D_GLOBE/KERNELS/OUTPUT_SUM.old'

  ! ======================================================

  integer, parameter :: NSPEC = NSPEC_CRUST_MANTLE
  integer, parameter :: NGLOB = NGLOB_CRUST_MANTLE

  ! transverse isotropic model files
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        model_vpv,model_vph,model_vsv,model_vsh,model_eta,model_rho
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        model_vpv_new,model_vph_new,model_vsv_new,model_vsh_new,model_eta_new,model_rho_new

  ! model updates
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        model_dbulk,model_dbetah,model_dbetav,model_deta

  ! old gradient updates
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        model_dbulk_old,model_dbetah_old,model_dbetav_old,model_deta_old

  ! kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        kernel_bulk,kernel_betav,kernel_betah,kernel_eta

  ! old kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        kernel_bulk_old,kernel_betav_old,kernel_betah_old,kernel_eta_old


  ! volume
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: x, y, z
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: ibool
  logical, dimension(NSPEC) :: ispec_is_tiso

  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_betav,norm_betah,norm_eta
  real(kind=CUSTOM_REAL) :: norm_bulk_old,norm_betav_old,norm_betah_old,norm_eta_old
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_betav_sum, &
    norm_betah_sum,norm_eta_sum

  ! model update length
  real(kind=CUSTOM_REAL) :: step_fac,step_length

  real(kind=CUSTOM_REAL) :: min_vpv,min_vph,min_vsv,min_vsh, &
    max_vpv,max_vph,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk, &
    min_rho,max_rho

  real(kind=CUSTOM_REAL) :: betav1,betah1,betav0,betah0,rho1,rho0, &
    betaiso1,betaiso0,eta1,eta0,alphav1,alphav0,alphah1,alphah0
  real(kind=CUSTOM_REAL) :: dbetaiso,dbulk

  integer :: nfile, myrank, sizeprocs,  ier
  integer :: i, j, k,ispec, iglob, ishell, n, it, j1, ib, npts_sem, ios
  character(len=150) :: sline, m_file, fname

  logical :: USE_MODEL_OLD

end module model_update_cg

!
!-------------------------------------------------------------------------------------------------
!

program add_model

  use model_update_cg

  implicit none

  ! ============ program starts here =====================

  ! initializes arrays
  call initialize()

  ! reads in parameters needed
  call read_parameters()

  ! reads in current transverse isotropic model files: vpv.. & vsv.. & eta & rho
  call read_model()

  ! reads in smoothed kernels: bulk, betav, betah, eta
  call read_kernels()

  ! reads in old (former inversion) smoothed kernels: bulk, betav, betah, eta
  call read_kernels_old()

  ! calculates gradient
  ! conjugate gradient method
  call get_gradient_cg()

  ! compute new model in terms of alpha, beta, eta and rho
  ! (see also Carl's Latex notes)

  ! model update:
  !   transverse isotropic update only in layer Moho to 220 (where SPECFEM3D_GLOBE considers TISO)
  !   everywhere else uses an isotropic update
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! initial model values
          eta0 = model_eta(i,j,k,ispec)
          betav0 = model_vsv(i,j,k,ispec)
          betah0 = model_vsh(i,j,k,ispec)
          rho0 = model_rho(i,j,k,ispec)
          alphav0 = model_vpv(i,j,k,ispec)
          alphah0 = model_vph(i,j,k,ispec)

          eta1 = 0._CUSTOM_REAL
          betav1 = 0._CUSTOM_REAL
          betah1 = 0._CUSTOM_REAL
          rho1 = 0._CUSTOM_REAL
          alphav1 = 0._CUSTOM_REAL
          alphah1 = 0._CUSTOM_REAL

          ! do not use transverse isotropy except if element is between d220 and Moho
!          if(.not. ( idoubling(ispec)==IFLAG_220_80 .or. idoubling(ispec)==IFLAG_80_MOHO) ) then
          if( .not. ispec_is_tiso(ispec) ) then
            ! isotropic model update

            ! no eta perturbation, since eta = 1 in isotropic media
            eta1 = eta0

            ! shear values
            ! isotropic kernel K_beta = K_betav + K_betah
            ! with same scaling step_length the model update dbeta_iso = dbetav + dbetah
            ! note:
            !   this step length can be twice as big as that given by the input
            dbetaiso = model_dbetav(i,j,k,ispec) + model_dbetah(i,j,k,ispec)
            betav1 = betav0 * exp( dbetaiso )
            betah1 = betah0 * exp( dbetaiso )
            ! note: betah is probably not really used in isotropic layers
            !         (see SPECFEM3D_GLOBE/get_model.f90)

            ! density: uses scaling relation with isotropic shear perturbations
            !               dln rho = RHO_SCALING * dln betaiso
            rho1 = rho0 * exp( RHO_SCALING * dbetaiso )

            ! alpha values
            dbulk = model_dbulk(i,j,k,ispec)
            alphav1 = sqrt( alphav0**2 * exp(2.0*dbulk) + FOUR_THIRDS * betav0**2 * ( &
                                exp(2.0*dbetaiso) - exp(2.0*dbulk) ) )
            alphah1 = sqrt( alphah0**2 * exp(2.0*dbulk) + FOUR_THIRDS * betah0**2 * ( &
                                exp(2.0*dbetaiso) - exp(2.0*dbulk) ) )
            ! note: alphah probably not used in SPECFEM3D_GLOBE

          else

            ! transverse isotropic model update

            ! eta value : limits updated values for eta range constraint
            eta1 = eta0 * exp( model_deta(i,j,k,ispec) )
            if( eta1 < LIMIT_ETA_MIN ) eta1 = LIMIT_ETA_MIN
            if( eta1 > LIMIT_ETA_MAX ) eta1 = LIMIT_ETA_MAX

            ! shear values
            betav1 = betav0 * exp( model_dbetav(i,j,k,ispec) )
            betah1 = betah0 * exp( model_dbetah(i,j,k,ispec) )

            ! density: uses scaling relation with Voigt average of shear perturbations
            betaiso0 = sqrt(  ( 2.0 * betav0**2 + betah0**2 ) / 3.0 )
            betaiso1 = sqrt(  ( 2.0 * betav1**2 + betah1**2 ) / 3.0 )
            dbetaiso = log( betaiso1 / betaiso0 )
            rho1 = rho0 * exp( RHO_SCALING * dbetaiso )

            ! alpha values
            dbulk = model_dbulk(i,j,k,ispec)
            alphav1 = sqrt( alphav0**2 * exp(2.0*dbulk) &
                            + FOUR_THIRDS * betav0**2 * ( &
                                exp(2.0*model_dbetav(i,j,k,ispec)) - exp(2.0*dbulk) ) )
            alphah1 = sqrt( alphah0**2 * exp(2.0*dbulk) &
                            + FOUR_THIRDS * betah0**2 * ( &
                                exp(2.0*model_dbetah(i,j,k,ispec)) - exp(2.0*dbulk) ) )

          endif


          ! stores new model values
          model_vpv_new(i,j,k,ispec) = alphav1
          model_vph_new(i,j,k,ispec) = alphah1
          model_vsv_new(i,j,k,ispec) = betav1
          model_vsh_new(i,j,k,ispec) = betah1
          model_eta_new(i,j,k,ispec) = eta1
          model_rho_new(i,j,k,ispec) = rho1

        enddo
      enddo
    enddo
  enddo

  ! stores new model in files
  call store_new_model()

  ! stores relative model perturbations
  call store_perturbations()

  ! computes volume element associated with points, calculates kernel integral for statistics
  call compute_volume()

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program add_model

!
!-------------------------------------------------------------------------------------------------
!

subroutine initialize()

! initializes arrays

  use model_update_cg
  implicit none

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if( sizeprocs /= nchunks_val*nproc_xi_val*nproc_eta_val ) then
    print*,'sizeprocs:',sizeprocs,nchunks_val,nproc_xi_val,nproc_eta_val
    call exit_mpi(myrank,'error number sizeprocs')
  endif

  ! model
  model_vpv = 0.0_CUSTOM_REAL
  model_vph = 0.0_CUSTOM_REAL
  model_vsv = 0.0_CUSTOM_REAL
  model_vsh = 0.0_CUSTOM_REAL
  model_eta = 0.0_CUSTOM_REAL
  model_rho = 0.0_CUSTOM_REAL

  model_vpv_new = 0.0_CUSTOM_REAL
  model_vph_new = 0.0_CUSTOM_REAL
  model_vsv_new = 0.0_CUSTOM_REAL
  model_vsh_new = 0.0_CUSTOM_REAL
  model_eta_new = 0.0_CUSTOM_REAL
  model_rho_new = 0.0_CUSTOM_REAL

  ! model updates
  model_dbulk = 0.0_CUSTOM_REAL
  model_dbetah = 0.0_CUSTOM_REAL
  model_dbetav = 0.0_CUSTOM_REAL
  model_deta = 0.0_CUSTOM_REAL

  ! gradients
  kernel_bulk = 0.0_CUSTOM_REAL
  kernel_betav = 0.0_CUSTOM_REAL
  kernel_betah = 0.0_CUSTOM_REAL
  kernel_eta = 0.0_CUSTOM_REAL

  ! old gradients
  kernel_bulk_old = 0.0_CUSTOM_REAL
  kernel_betav_old = 0.0_CUSTOM_REAL
  kernel_betah_old = 0.0_CUSTOM_REAL
  kernel_eta_old = 0.0_CUSTOM_REAL

end subroutine initialize

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_parameters()

! reads in parameters needed

  use model_update_cg
  implicit none
  character(len=150) :: s_step_fac

  ! subjective step length to multiply to the gradient
  !step_fac = 0.03

  call getarg(1,s_step_fac)

  if (trim(s_step_fac) == '') then
    call exit_MPI(myrank,'Usage: add_model_globe_tiso_cg step_factor')
  endif

  ! read in parameter information
  read(s_step_fac,*) step_fac
  !if( abs(step_fac) < 1.e-10) then
  !  print*,'error: step factor ',step_fac
  !  call exit_MPI(myrank,'error step factor')
  !endif

  if (myrank == 0) then
    print*,'defaults'
    print*,'  NPROC_XI , NPROC_ETA: ',nproc_xi_val,nproc_eta_val
    print*,'  NCHUNKS: ',nchunks_val
    print*
    print*,'model update for vsv,vsh,vpv,vph,eta,rho:'
    print*,'  step_fac = ',step_fac
    print*

  endif


end subroutine read_parameters

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_model()

! reads in current transverse isotropic model: vpv.. & vsv.. & eta & rho

  use model_update_cg
  implicit none
  integer, dimension(:), allocatable :: idummy
  logical, dimension(:), allocatable :: ldummy

  ! vpv model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vpv.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_vpv(:,:,:,1:nspec)
  close(12)

  ! vph model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vph.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_vph(:,:,:,1:nspec)
  close(12)

  ! vsv model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vsv.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_vsv(:,:,:,1:nspec)
  close(12)

  ! vsh model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vsh.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_vsh(:,:,:,1:nspec)
  close(12)

  ! eta model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_eta.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_eta(:,:,:,1:nspec)
  close(12)

  ! rho model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_rho.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_rho(:,:,:,1:nspec)
  close(12)

  ! statistics
  call mpi_reduce(minval(model_vpv),min_vpv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vpv),max_vpv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_vph),min_vph,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vph),max_vph,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_vsv),min_vsv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vsv),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_vsh),min_vsh,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vsh),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_eta),min_eta,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_eta),max_eta,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_rho),min_rho,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_rho),max_rho,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'initial models:'
    print*,'  vpv min/max: ',min_vpv,max_vpv
    print*,'  vph min/max: ',min_vph,max_vph
    print*,'  vsv min/max: ',min_vsv,max_vsv
    print*,'  vsh min/max: ',min_vsh,max_vsh
    print*,'  eta min/max: ',min_eta,max_eta
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif

  ! allocates temporary array
  allocate(idummy(NSPEC),ldummy(NSPEC),stat=ier)
  if( ier /= 0 ) stop 'error allocating ldummy'

  ! global addressing
  write(m_file,'(a,i6.6,a)') 'topo/proc',myrank,'_reg1_solver_data_2.bin'
  open(11,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(11) x(1:nglob)
  read(11) y(1:nglob)
  read(11) z(1:nglob)
  read(11) ibool(:,:,:,1:nspec)
  read(11) idummy(1:nspec) ! idoubling
  read(11) ldummy(1:nspec) ! is_on_a_slice_edge
  read(11) ispec_is_tiso(1:nspec)
  close(11)

  deallocate(idummy,ldummy)

end subroutine read_model

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_kernels()

! reads in smoothed kernels: bulk, betav, betah, eta

  use model_update_cg
  implicit none

  ! bulk kernel
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_bulk_c_kernel_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_bulk(:,:,:,1:nspec)
  close(12)

  ! betav kernel
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_bulk_betav_kernel_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_betav(:,:,:,1:nspec)
  close(12)

  ! betah kernel
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_bulk_betah_kernel_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_betah(:,:,:,1:nspec)
  close(12)

  ! eta kernel
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_eta_kernel_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_eta(:,:,:,1:nspec)
  close(12)


  ! statistics
  call mpi_reduce(minval(kernel_bulk),min_bulk,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_bulk),max_bulk,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(kernel_betah),min_vsh,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_betah),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(kernel_betav),min_vsv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_betav),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(kernel_eta),min_eta,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_eta),max_eta,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'initial kernels:'
    print*,'  bulk min/max : ',min_bulk,max_bulk
    print*,'  betav min/max: ',min_vsv,max_vsv
    print*,'  betah min/max: ',min_vsh,max_vsh
    print*,'  eta min/max  : ',min_eta,max_eta
    print*
  endif

end subroutine read_kernels

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_kernels_old()

! reads in smoothed kernels from former iteration in OUTPUT_SUM.old/ : bulk, betav, betah, eta

  use model_update_cg
  implicit none
  logical:: exist,exist_all,use_model_old_all

  ! checks if files are available:
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_bulk_c_kernel_smooth.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if( .not. exist ) then
    print*,'error file does not exist: ',trim(m_file)
    call exit_mpi(myrank,'file not exist')
  endif
  ! makes sure all processes have same flag
  call MPI_ALLREDUCE(exist,exist_all,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ier)
  if( .not. exist_all ) then
    print*,'old kernels do not exist: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! bulk kernel
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_bulk_c_kernel_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_bulk_old(:,:,:,1:nspec)
  close(12)

  ! betav kernel
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_bulk_betav_kernel_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_betav_old(:,:,:,1:nspec)
  close(12)

  ! betah kernel
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_bulk_betah_kernel_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_betah_old(:,:,:,1:nspec)
  close(12)

  ! eta kernel
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_eta_kernel_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_eta_old(:,:,:,1:nspec)
  close(12)


  ! statistics
  call mpi_reduce(minval(kernel_bulk_old),min_bulk,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_bulk_old),max_bulk,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(kernel_betah_old),min_vsh,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_betah_old),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(kernel_betav_old),min_vsv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_betav_old),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(kernel_eta_old),min_eta,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_eta_old),max_eta,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'old kernels:'
    print*,'  bulk min/max : ',min_bulk,max_bulk
    print*,'  betav min/max: ',min_vsv,max_vsv
    print*,'  betah min/max: ',min_vsh,max_vsh
    print*,'  eta min/max  : ',min_eta,max_eta
    print*
  endif

  ! reads in old gradient directions (phi_(n-1))
  USE_MODEL_OLD = .true.

  ! checks if files are available:
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_dbulk_c.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if( .not. exist ) then
    print*,'old kernel updates do not exist: ',trim(m_file)
    USE_MODEL_OLD = .false.
  endif
  ! makes sure all processes have same flag
  call MPI_ALLREDUCE(exist,exist_all,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ier)
  if( .not. exist_all ) then
    if( myrank == 0 ) print*,'old kernel updates do not exist for all: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! makes sure all processes have same flag
  use_model_old_all = .false.
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  call MPI_ALLREDUCE(USE_MODEL_OLD,use_model_old_all,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ier)
  if( .not. use_model_old_all ) then
    print*,'old kernel updates exists, not consistent for all: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! reads in old gradients
  if( USE_MODEL_OLD ) then
    ! bulk kernel
    fname = 'dbulk_c'
    write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
    open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(12) model_dbulk_old(:,:,:,1:nspec)
    close(12)

    ! betav kernel
    fname = 'dbetav'
    write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
    open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(12) model_dbetav_old(:,:,:,1:nspec)
    close(12)

    ! betah kernel
    fname = 'dbetah'
    write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
    open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(12) model_dbetah_old(:,:,:,1:nspec)
    close(12)

    ! eta kernel
    fname = 'deta'
    write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
    open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(12) model_deta_old(:,:,:,1:nspec)
    close(12)


    ! statistics
    call mpi_reduce(minval(model_dbulk_old),min_bulk,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
    call mpi_reduce(maxval(model_dbulk_old),max_bulk,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

    call mpi_reduce(minval(model_dbetah_old),min_vsh,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
    call mpi_reduce(maxval(model_dbetah_old),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

    call mpi_reduce(minval(model_dbetav_old),min_vsv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
    call mpi_reduce(maxval(model_dbetav_old),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

    call mpi_reduce(minval(model_deta_old),min_eta,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
    call mpi_reduce(maxval(model_deta_old),max_eta,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

    if( myrank == 0 ) then
      print*,'old kernel updates:'
      print*,'  bulk min/max : ',min_bulk,max_bulk
      print*,'  betav min/max: ',min_vsv,max_vsv
      print*,'  betah min/max: ',min_vsh,max_vsh
      print*,'  eta min/max  : ',min_eta,max_eta
      print*
    endif
  endif ! USE_MODEL_OLD

end subroutine read_kernels_old


!
!-------------------------------------------------------------------------------------------------
!

subroutine get_gradient_cg()

! calculates gradient based on a conjugate gradient method
!
! based on: Tarantola, inverse problem theory, 2005.
!                  section 6.22.7 conjugate directions, page 217.
!                  formula for alpha_n based on Polak & Ribiere (1969)
!
! note: we use a preconditioner F_0 = 1, thus lambda_n = gamma_n in (6.322)
!          and use gamma_n as the smoothed kernel (for bulk_c, bulk_betav,..).
!
!          however, one could see smoothing as preconditioner F_0, thus
!          gamma_n would be un-smoothed kernel and lambda_n would be smoothed one...
!          i'm not sure if this makes a difference.

  use model_update_cg
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha_bulk,alpha_betav,alpha_betah,alpha_eta,alpha_all
  real(kind=CUSTOM_REAL) :: minmax(4),depthmax(2),depthmax_radius(2),max
  real(kind=CUSTOM_REAL) :: r,rmax_vsv,rmax_vsh,depthmax_depth
  integer :: maxindex(1)
  real(kind=CUSTOM_REAL) :: ratio_bulk,ratio_betav,ratio_betah,ratio_eta  

  ! ------------------------------------------------------------------------

  ! sets maximum update in this depth range
  logical,parameter :: use_depth_maximum = .true.
  ! normalized radii
  real(kind=CUSTOM_REAL),parameter :: R_top = (6371.0 - 50.0 ) / R_EARTH_KM ! shallow depth
  real(kind=CUSTOM_REAL),parameter :: R_bottom = (6371.0 - 100.0 ) / R_EARTH_KM ! deep depth
  ! uses separate scaling for each parameters bulk,betav,betah,eta
  ! (otherwise it will calculate a single steplengths to scale all gradients)
  logical,parameter :: use_separate_steplengths = .false.

  ! ------------------------------------------------------------------------

  ! old kernel/gradient
  ! length ( gamma_(n-1)^T * lambda_(n-1) )
  norm_bulk_old = sum( kernel_bulk_old * kernel_bulk_old )
  norm_betav_old = sum( kernel_betav_old * kernel_betav_old )
  norm_betah_old = sum( kernel_betah_old * kernel_betah_old )
  norm_eta_old = sum( kernel_eta_old * kernel_eta_old )

  call mpi_reduce(norm_bulk_old,norm_bulk_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_betav_old,norm_betav_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_betah_old,norm_betah_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_eta_old,norm_eta_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  ! don't use square root, just take gamma^T * gamma
  norm_bulk_old = norm_bulk_sum
  norm_betav_old = norm_betav_sum
  norm_betah_old = norm_betah_sum
  norm_eta_old = norm_eta_sum

  if( myrank == 0 ) then
    print*,'norm squared old gradients:'
    print*,'  bulk : ',norm_bulk_old
    print*,'  betav: ',norm_betav_old
    print*,'  betah: ',norm_betah_old
    print*,'  eta  : ',norm_eta_old
    print*
  endif

  ! checks lengths
  if( myrank == 0 ) then
    if( norm_bulk_old < 1.e-22 ) call exit_mpi(myrank,'norm old gradient bulk is zero')
    if( norm_betav_old < 1.e-22 ) call exit_mpi(myrank,'norm old gradient betav is zero')
    if( norm_betah_old < 1.e-22 ) call exit_mpi(myrank,'norm old gradient betah is zero')
    if( norm_eta_old < 1.e-22 ) call exit_mpi(myrank,'norm old gradient eta is zero')
  endif

  ! Powell, 1977: checks orthogonality between old and new gradients
  ! gets length of ( gamma_(n-1)^T * gamma_n )
  norm_bulk = sum( kernel_bulk_old * kernel_bulk )
  norm_betav = sum( kernel_betav_old * kernel_betav )
  norm_betah = sum( kernel_betah_old * kernel_betah )
  norm_eta = sum( kernel_eta_old * kernel_eta )

  call mpi_reduce(norm_bulk,norm_bulk_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_betav,norm_betav_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_betah,norm_betah_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_eta,norm_eta_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  ! ratio:  ( g_n * g_n-1 ) / ( g_n-1 * g_n-1)
  ratio_bulk = norm_bulk_sum / norm_bulk_old
  ratio_betav = norm_betav_sum / norm_betav_old
  ratio_betah = norm_betah_sum / norm_betah_old
  ratio_eta = norm_eta_sum / norm_eta_old

  ! if ratio > 0.2 (empirical threshold value), then one should restart with a steepest descent
  if( myrank == 0 ) then
    print*,'Powell ratio: (> 0.2 then restart with steepest descent)'
    print*,'  bulk : ',ratio_bulk
    print*,'  betav: ',ratio_betav
    print*,'  betah: ',ratio_betah
    print*,'  eta  : ',ratio_eta
    print*
    if( ratio_bulk > 0.2 .and. ratio_betav > 0.2 .and. ratio_betah > 0.2 &
      .and. ratio_eta > 0.2 ) then
      print*,'  please consider doing a steepest descent instead cg...'
      print*
    endif
  endif


  ! difference kernel/gradients
  ! length ( ( gamma_n - gamma_(n-1))^T * lambda_n )
  norm_bulk = sum( (kernel_bulk - kernel_bulk_old) * kernel_bulk )
  norm_betav = sum( (kernel_betav - kernel_betav_old) * kernel_betav )
  norm_betah = sum( (kernel_betah - kernel_betah_old) * kernel_betah )
  norm_eta = sum( (kernel_eta - kernel_eta_old) * kernel_eta )

  call mpi_reduce(norm_bulk,norm_bulk_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_betav,norm_betav_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_betah,norm_betah_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_eta,norm_eta_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  ! don't take square root, since norm_bulk_sum could be negative
  ! just use (gamma_n - gamma_n-1)^T * lambda_n
  norm_bulk = norm_bulk_sum
  norm_betav = norm_betav_sum
  norm_betah = norm_betah_sum
  norm_eta = norm_eta_sum

  if( myrank == 0 ) then
    print*,'norm squared difference gradients:'
    print*,'  bulk : ',norm_bulk
    print*,'  betav: ',norm_betav
    print*,'  betah: ',norm_betah
    print*,'  eta  : ',norm_eta
    print*
  endif

  ! calculates ratio based on Polak & Ribiere (1969)
  if( myrank == 0 ) then
    if( use_separate_steplengths ) then
      ! calculates steplength alpha for each parameter
      alpha_bulk = norm_bulk / norm_bulk_old
      alpha_betav = norm_betav / norm_betav_old
      alpha_betah = norm_betah / norm_betah_old
      alpha_eta = norm_eta / norm_eta_old

      ! only if contribution is positive it will be considered, otherwise
      ! we set it to zero so that it becomes a steepest descent update
      if( alpha_bulk < 0.0 ) then
        alpha_bulk = 0.0
      endif
      if( alpha_betav < 0.0 ) then
        alpha_betav = 0.0
      endif
      if( alpha_betah < 0.0 ) then
        alpha_betah = 0.0
      endif
      if( alpha_eta < 0.0 ) then
        alpha_eta = 0.0
      endif

    else
      ! calculates only a single steplength applied to all
      alpha_all = (norm_bulk + norm_betav + norm_betah + norm_eta) &
                  / (norm_bulk_old + norm_betav_old + norm_betah_old + norm_eta_old)

      ! only if contribution is positive it will be considered, otherwise
      ! we set it to zero so that it becomes a steepest descent update
      if( alpha_all < 0.0 ) then
        alpha_all = 0.0
      endif

      ! sets each steplength to same single one
      alpha_bulk = alpha_all
      alpha_betav = alpha_all
      alpha_betah = alpha_all
      alpha_eta = alpha_all
    endif
    ! user output
    print*,'alpha gradients:'
    print*,'  bulk : ',alpha_bulk
    print*,'  betav: ',alpha_betav
    print*,'  betah: ',alpha_betah
    print*,'  eta  : ',alpha_eta
    print*
  endif
  ! broadcast values from rank 0 to all others
  call mpi_bcast(alpha_bulk,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call mpi_bcast(alpha_betav,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call mpi_bcast(alpha_betah,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call mpi_bcast(alpha_eta,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  ! initializes kernel maximum
  depthmax(:) = 0._CUSTOM_REAL

  ! gradient in negative direction
  if( USE_MODEL_OLD ) then
    ! uses old kernel/gradient updates ( phi_n-1 )
    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

              ! note: uses old gradient update (phi_(n-1) as model_bulk_old), but
              !           given in negative gradient direction

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
              if( use_depth_maximum ) then
                ! get radius of point
                iglob = ibool(i,j,k,ispec)
                r = sqrt( x(iglob)*x(iglob) + y(iglob)*y(iglob) + z(iglob)*z(iglob) )

                ! stores maximum kernel betav/betah value in this depth slice,
                ! since betav/betah are most likely dominating
                if( r < R_top .and. r > R_bottom ) then
                  ! kernel betav value
                  max_vsv = abs( model_dbetav(i,j,k,ispec) )
                  if( depthmax(1) < max_vsv ) then
                    depthmax(1) = max_vsv
                    depthmax_radius(1) = r
                  endif
                  ! kernel betav value
                  max_vsh = abs( model_dbetah(i,j,k,ispec) )
                  if( depthmax(2) < max_vsh ) then
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
    ! uses only old kernel/gradients
    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

              ! note: uses old gradients (lambda_(n-1) ) in negative gradient direction

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
              if( use_depth_maximum ) then
                ! get radius of point
                iglob = ibool(i,j,k,ispec)
                r = sqrt( x(iglob)*x(iglob) + y(iglob)*y(iglob) + z(iglob)*z(iglob) )

                ! stores maximum kernel betav/betah value in this depth slice,
                ! since betav/betah are most likely dominating
                if( r < R_top .and. r > R_bottom ) then
                  ! kernel betav value
                  max_vsv = abs( model_dbetav(i,j,k,ispec) )
                  if( depthmax(1) < max_vsv ) then
                    depthmax(1) = max_vsv
                    depthmax_radius(1) = r
                  endif
                  ! kernel betav value
                  max_vsh = abs( model_dbetah(i,j,k,ispec) )
                  if( depthmax(2) < max_vsh ) then
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
  ! note: stores these new gradients before we scale them with the step length
  call store_kernel_updates()

  ! statistics
  call mpi_reduce(minval(model_dbulk),min_bulk,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dbulk),max_bulk,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_dbetav),min_vsv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dbetav),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_dbetah),min_vsh,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dbetah),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_deta),min_eta,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_deta),max_eta,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'initial gradient updates:'
    print*,'  bulk min/max : ',min_bulk,max_bulk
    print*,'  betav min/max: ',min_vsv,max_vsv
    print*,'  betah min/max: ',min_vsh,max_vsh
    print*,'  eta min/max  : ',min_eta,max_eta
    print*
  endif

  ! determines maximum kernel betav value within given radius
  if( use_depth_maximum ) then
    ! maximum of all processes stored in max_vsv
    call mpi_reduce(depthmax(1),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
    call mpi_reduce(depthmax(2),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
    call mpi_reduce(depthmax_radius(1),rmax_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
    call mpi_reduce(depthmax_radius(2),rmax_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  endif

  ! determines step length
  ! based on maximum gradient value (either vsv or vsh)
  if( myrank == 0 ) then

    ! determines maximum kernel betav value within given radius
    if( use_depth_maximum ) then
      depthmax(1) = max_vsv
      depthmax(2) = max_vsh
      depthmax_radius(1) = rmax_vsv
      depthmax_radius(2) = rmax_vsh

      max = maxval(depthmax)
      maxindex = maxloc(depthmax)
      depthmax_depth = depthmax_radius(maxindex(1))
      depthmax_depth = 6371.0 *( 1.0 - depthmax_depth )
      ! maximum in given depth range
      print*,'  using depth maximum between 50km and 100km: ',max
      print*,'  approximate depth maximum: ',depthmax_depth
      print*
    else
      ! maximum gradient values
      minmax(1) = abs(min_vsv)
      minmax(2) = abs(max_vsv)
      minmax(3) = abs(min_vsh)
      minmax(4) = abs(max_vsh)

      ! maximum value of all kernel maxima
      max = maxval(minmax)
      print*,'  using maximum: ',max
      print*
    endif

    ! chooses step length such that it becomes the desired, given step factor as inputted
    step_length = step_fac/max

    print*,'  step length : ',step_length
    print*

  endif
  call mpi_bcast(step_length,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)


  ! gradient length sqrt( v^T * v )
  norm_bulk = sum( model_dbulk * model_dbulk )
  norm_betav = sum( model_dbetav * model_dbetav )
  norm_betah = sum( model_dbetah * model_dbetah )
  norm_eta = sum( model_deta * model_deta )

  call mpi_reduce(norm_bulk,norm_bulk_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_betav,norm_betav_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_betah,norm_betah_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_eta,norm_eta_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  norm_bulk = sqrt(norm_bulk_sum)
  norm_betav = sqrt(norm_betav_sum)
  norm_betah = sqrt(norm_betah_sum)
  norm_eta = sqrt(norm_eta_sum)

  if( myrank == 0 ) then
    print*,'norm model updates:'
    print*,'  bulk : ',norm_bulk
    print*,'  betav: ',norm_betav
    print*,'  betah: ',norm_betah
    print*,'  eta  : ',norm_eta
    print*
  endif

  ! multiply model updates by a subjective factor that will change the step
  model_dbulk(:,:,:,:) = step_length * model_dbulk(:,:,:,:)
  model_dbetav(:,:,:,:) = step_length * model_dbetav(:,:,:,:)
  model_dbetah(:,:,:,:) = step_length * model_dbetah(:,:,:,:)
  model_deta(:,:,:,:) = step_length * model_deta(:,:,:,:)


  ! statistics
  call mpi_reduce(minval(model_dbulk),min_bulk,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dbulk),max_bulk,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_dbetav),min_vsv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dbetav),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_dbetah),min_vsh,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dbetah),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_deta),min_eta,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_deta),max_eta,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'scaled gradients:'
    print*,'  bulk min/max : ',min_bulk,max_bulk
    print*,'  betav min/max: ',min_vsv,max_vsv
    print*,'  betah min/max: ',min_vsh,max_vsh
    print*,'  eta min/max  : ',min_eta,max_eta
    print*
  endif

end subroutine get_gradient_cg

!
!-------------------------------------------------------------------------------------------------
!

subroutine store_kernel_updates()

! file output for new model

  use model_update_cg
  implicit none

  ! kernel updates
  fname = 'dbulk_c'
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_dbulk
  close(12)

  fname = 'dbetav'
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_dbetav
  close(12)

  fname = 'dbetah'
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_dbetah
  close(12)

  fname = 'deta'
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_deta
  close(12)

end subroutine store_kernel_updates

!
!-------------------------------------------------------------------------------------------------
!

subroutine store_new_model()

! file output for new model

  use model_update_cg
  implicit none

  ! vpv model
  call mpi_reduce(maxval(model_vpv_new),max_vpv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vpv_new),min_vpv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'vpv_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vpv_new
  close(12)

  ! vph model
  call mpi_reduce(maxval(model_vph_new),max_vph,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vph_new),min_vph,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'vph_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vph_new
  close(12)

  ! vsv model
  call mpi_reduce(maxval(model_vsv_new),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vsv_new),min_vsv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'vsv_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vsv_new
  close(12)

  ! vsh model
  call mpi_reduce(maxval(model_vsh_new),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vsh_new),min_vsh,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'vsh_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vsh_new
  close(12)

  ! eta model
  call mpi_reduce(maxval(model_eta_new),max_eta,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_eta_new),min_eta,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'eta_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_eta_new
  close(12)

  ! rho model
  call mpi_reduce(maxval(model_rho_new),max_rho,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_rho_new),min_rho,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'rho_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_rho_new
  close(12)


  if( myrank == 0 ) then
    print*,'new models:'
    print*,'  vpv min/max: ',min_vpv,max_vpv
    print*,'  vph min/max: ',min_vph,max_vph
    print*,'  vsv min/max: ',min_vsv,max_vsv
    print*,'  vsh min/max: ',min_vsh,max_vsh
    print*,'  eta min/max: ',min_eta,max_eta
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif


end subroutine store_new_model

!
!-------------------------------------------------------------------------------------------------
!

subroutine store_perturbations()

! file output for new model

  use model_update_cg
  implicit none
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: total_model

  ! vpv relative perturbations
  ! logarithmic perturbation: log( v_new) - log( v_old) = log( v_new / v_old )
  total_model = 0.0_CUSTOM_REAL
  where( model_vpv /= 0.0 ) total_model = log( model_vpv_new / model_vpv)
  ! or
  ! linear approximation: (v_new - v_old) / v_old
  !where( model_vpv /= 0.0 ) total_model = ( model_vpv_new - model_vpv) / model_vpv

  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvpvvpv.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_vpv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_vpv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! vph relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vph /= 0.0 ) total_model = log( model_vph_new / model_vph)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvphvph.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_vph,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_vph,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! vsv relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vsv /= 0.0 ) total_model = log( model_vsv_new / model_vsv)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvsvvsv.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_vsv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! vsh relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vsh /= 0.0 ) total_model = log( model_vsh_new / model_vsh)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvshvsh.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_vsh,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! eta relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_eta /= 0.0 ) total_model = log( model_eta_new / model_eta)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_detaeta.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_eta,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_eta,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! rho relative model perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_rho /= 0.0 ) total_model = log( model_rho_new / model_rho)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_drhorho.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_rho,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_rho,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'relative update:'
    print*,'  dvpv/vpv min/max: ',min_vpv,max_vpv
    print*,'  dvph/vph min/max: ',min_vph,max_vph
    print*,'  dvsv/vsv min/max: ',min_vsv,max_vsv
    print*,'  dvsh/vsh min/max: ',min_vsh,max_vsh
    print*,'  deta/eta min/max: ',min_eta,max_eta
    print*,'  drho/rho min/max: ',min_rho,max_rho
    print*
  endif

end subroutine store_perturbations

!
!-------------------------------------------------------------------------------------------------
!

subroutine compute_volume()

! computes volume element associated with points

  use model_update_cg
  implicit none
  ! jacobian
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: jacobian
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl, &
    jacobianl,volumel
  ! integration values
  real(kind=CUSTOM_REAL) :: integral_bulk_sum,integral_betav_sum, &
    integral_betah_sum,integral_eta_sum
  real(kind=CUSTOM_REAL) :: integral_bulk,integral_betav,&
    integral_betah,integral_eta
  real(kind=CUSTOM_REAL) :: volume_glob,volume_glob_sum
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
  write(m_file,'(a,i6.6,a)') 'topo/proc',myrank,'_reg1_solver_data_1.bin'
  open(11,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(11) xix
  read(11) xiy
  read(11) xiz
  read(11) etax
  read(11) etay
  read(11) etaz
  read(11) gammax
  read(11) gammay
  read(11) gammaz
  close(11)

  jacobian = 0.0
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! gets derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
          ! computes the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))
          jacobian(i,j,k,ispec) = jacobianl

          !if( abs(jacobianl) < 1.e-8 ) then
          !  print*,'rank ',myrank,'jacobian: ',jacobianl,i,j,k,wgll_cube(i,j,k)
          !endif

        enddo
      enddo
    enddo
  enddo

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
          if( iglob == 0 ) then
            print*,'iglob zero',i,j,k,ispec
            print*
            print*,'ibool:',ispec
            print*,ibool(:,:,:,ispec)
            print*
            call exit_MPI(myrank,'error ibool')
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
          if( isNaN(integral_bulk) ) then
            print*,'error NaN: ',integral_bulk
            print*,'rank:',myrank
            print*,'i,j,k,ispec:',i,j,k,ispec
            print*,'volumel: ',volumel,'kernel_bulk:',kernel_bulk(i,j,k,ispec)
            call exit_MPI(myrank,'error NaN')
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
  ! kernel integration: for whole volume
  call mpi_reduce(integral_bulk,integral_bulk_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(integral_betav,integral_betav_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(integral_betah,integral_betah_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(integral_eta,integral_eta_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(volume_glob,volume_glob_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'integral kernels:'
    print*,'  bulk : ',integral_bulk_sum
    print*,'  betav : ',integral_betav_sum
    print*,'  betah : ',integral_betah_sum
    print*,'  eta : ',integral_eta_sum
    print*
    print*,'  total volume:',volume_glob_sum
    print*
  endif

  ! norms: for whole volume
  call mpi_reduce(norm_bulk,norm_bulk_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_betav,norm_betav_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_betah,norm_betah_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_eta,norm_eta_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  norm_bulk = sqrt(norm_bulk_sum)
  norm_betav = sqrt(norm_betav_sum)
  norm_betah = sqrt(norm_betah_sum)
  norm_eta = sqrt(norm_eta_sum)

  if( myrank == 0 ) then
    print*,'norm kernels:'
    print*,'  bulk : ',norm_bulk
    print*,'  betav : ',norm_betav
    print*,'  betah : ',norm_betah
    print*,'  eta : ',norm_eta
    print*
  endif

  ! root-mean square
  call mpi_reduce(rms_vpv,rms_vpv_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(rms_vph,rms_vph_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(rms_vsv,rms_vsv_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(rms_vsh,rms_vsh_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(rms_eta,rms_eta_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(rms_rho,rms_rho_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  rms_vpv = sqrt( rms_vpv_sum / volume_glob_sum )
  rms_vph = sqrt( rms_vph_sum / volume_glob_sum )
  rms_vsv = sqrt( rms_vsv_sum / volume_glob_sum )
  rms_vsh = sqrt( rms_vsh_sum / volume_glob_sum )
  rms_eta = sqrt( rms_eta_sum / volume_glob_sum )
  rms_rho = sqrt( rms_rho_sum / volume_glob_sum )

  if( myrank == 0 ) then
    print*,'root-mean square of perturbations:'
    print*,'  vpv : ',rms_vpv
    print*,'  vph : ',rms_vph
    print*,'  vsv : ',rms_vsv
    print*,'  vsh : ',rms_vsh
    print*,'  eta : ',rms_eta
    print*,'  rho : ',rms_rho
    print*
  endif

end subroutine compute_volume

