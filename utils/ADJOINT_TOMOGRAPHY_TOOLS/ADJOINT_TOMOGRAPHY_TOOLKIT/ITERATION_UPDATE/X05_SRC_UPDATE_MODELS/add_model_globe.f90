! add_model_globe_tiso
!
! this program can be used to update TRANSVERSE ISOTROPIC model files
! based on smoothed event kernels.
! the kernels are given for tranverse isotropic parameters (bulk_c,bulk_betav,bulk_betah,eta).
!
! the algorithm uses a steepest descent method with a step length
! determined by the given maximum update percentage.
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
!       proc000***_reg1_eta_kernel_smooth.bin
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
! USAGE: ./add_model_globe_tiso 0.3

module model_update_tiso

  include 'mpif.h'
  include '../../SHARE_FILES/HEADER_FILES/constants.h'
  include '../../SHARE_FILES/HEADER_FILES/precision.h'
  include '../../SHARE_FILES/HEADER_FILES/values_from_mesher.h'

  ! ======================================================

  ! density scaling factor with shear perturbations
  ! see e.g. Montagner & Anderson (1989), Panning & Romanowicz (2006)
  real(kind=CUSTOM_REAL),parameter :: RHO_SCALING = 0.33_CUSTOM_REAL

  ! constraint on eta model
  real(kind=CUSTOM_REAL),parameter :: LIMIT_ETA_MIN = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: LIMIT_ETA_MAX = 1.5_CUSTOM_REAL

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

  ! kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        kernel_bulk,kernel_betav,kernel_betah,kernel_eta

  ! volume
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: x, y, z
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: ibool
  integer, dimension(NSPEC) :: idoubling

  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_betav,norm_betah,norm_eta
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_betav_sum, &
    norm_betah_sum,norm_eta_sum

  ! model update length
  real(kind=CUSTOM_REAL) :: step_fac,step_length

  real(kind=CUSTOM_REAL) :: min_vpv,min_vph,min_vsv,min_vsh, &
    max_vpv,max_vph,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk, &
    min_rho,max_rho,max,minmax(4)

  real(kind=CUSTOM_REAL) :: betav1,betah1,betav0,betah0,rho1,rho0, &
    betaiso1,betaiso0,eta1,eta0,alphav1,alphav0,alphah1,alphah0
  real(kind=CUSTOM_REAL) :: dbetaiso,dbulk

  integer :: nfile, myrank, sizeprocs,  ier
  integer :: i, j, k,ispec, iglob, ishell, n, it, j1, ib, npts_sem, ios
  character(len=256) :: sline, m_file, fname
  character(len=256) :: input_model,input_kernel,output_model

end module model_update_tiso

!
!-------------------------------------------------------------------------------------------------
!

program add_model

  use model_update_tiso

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

  ! computes volume element associated with points, calculates kernel integral for statistics
  call compute_volume()

  ! calculates gradient
  ! steepest descent method
  call get_gradient()

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
          if(.not. ( idoubling(ispec)== IFLAG_670_220 .or.idoubling(ispec)==IFLAG_220_80 .or. idoubling(ispec)==IFLAG_80_MOHO) ) then

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

  ! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program add_model

!
!-------------------------------------------------------------------------------------------------
!

subroutine initialize()

! initializes arrays

  use model_update_tiso
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

end subroutine initialize

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_parameters()

! reads in parameters needed

  use model_update_tiso
  implicit none
  character(len=150) :: s_step_fac

  ! subjective step length to multiply to the gradient
  !step_fac = 0.03

  call getarg(1,s_step_fac)
!> Hejun Zhu
  call getarg(2,input_model)
  call getarg(3,input_kernel)
  call getarg(4,output_model)
!< Hejun Zhu


!> Hejun Zhu
!  if (trim(s_step_fac) == '') then
!    call exit_MPI(myrank,'Usage: add_model_globe_tiso step_factor')
!  endif
  if (trim(s_step_fac) == '' .or. trim(input_model) == '' &
      .or. trim(input_kernel) == ''.or. trim(output_model) == '') then
      call exit_MPI(myrank, 'Usage: add model_globe_tiso step_factor input_model input_kernel output_model')
  endif
!< Hejun Zhu


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
    print*,' input model dir = ', trim(input_model)
    print*,' input gradient dir=', trim(input_kernel)
    print*,' output model dir= ', trim(output_model)
    print*

  endif


end subroutine read_parameters

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_model()

! reads in current transverse isotropic model: vpv.. & vsv.. & eta & rho

  use model_update_tiso
  implicit none

  ! vpv model
  write(m_file,'(a,i6.6,a)') trim(input_model)//'/proc',myrank,'_reg1_vpv.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_vpv(:,:,:,1:nspec)
  close(12)

  ! vph model
  write(m_file,'(a,i6.6,a)') trim(input_model)//'/proc',myrank,'_reg1_vph.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_vph(:,:,:,1:nspec)
  close(12)

  ! vsv model
  write(m_file,'(a,i6.6,a)') trim(input_model)//'/proc',myrank,'_reg1_vsv.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_vsv(:,:,:,1:nspec)
  close(12)

  ! vsh model
  write(m_file,'(a,i6.6,a)') trim(input_model)//'/proc',myrank,'_reg1_vsh.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_vsh(:,:,:,1:nspec)
  close(12)

  ! eta model
  write(m_file,'(a,i6.6,a)') trim(input_model)//'/proc',myrank,'_reg1_eta.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_eta(:,:,:,1:nspec)
  close(12)

  ! rho model
  write(m_file,'(a,i6.6,a)') trim(input_model)//'/proc',myrank,'_reg1_rho.bin'
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

end subroutine read_model

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_kernels()

! reads in smoothed kernels: bulk, betav, betah, eta

  use model_update_tiso
  implicit none

  ! bulk kernel
  write(m_file,'(a,i6.6,a)') trim(input_kernel)//'/proc',myrank,'_reg1_bulk_c_kernel_precond_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_bulk(:,:,:,1:nspec)
  close(12)

  ! betav kernel
  write(m_file,'(a,i6.6,a)') trim(input_kernel)//'/proc',myrank,'_reg1_bulk_betav_kernel_precond_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_betav(:,:,:,1:nspec)
  close(12)

  ! betah kernel
  write(m_file,'(a,i6.6,a)') trim(input_kernel)//'/proc',myrank,'_reg1_bulk_betah_kernel_precond_smooth.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) kernel_betah(:,:,:,1:nspec)
  close(12)

  ! eta kernel
  write(m_file,'(a,i6.6,a)') trim(input_kernel)//'/proc',myrank,'_reg1_eta_kernel_precond_smooth.bin'
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

subroutine compute_volume()

! computes volume element associated with points

  use model_update_tiso
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

  ! global addressing
  write(m_file,'(a,i6.6,a)') &
      '/tigress-hsm/hejunzhu/2011EUROPE_ITERATION_UPDATE/EUROPE_TOPOLOGY_FILE/proc',myrank,'_reg1_solver_data_2.bin'
  open(11,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(11) x(1:nglob)
  read(11) y(1:nglob)
  read(11) z(1:nglob)
  read(11) ibool(:,:,:,1:nspec)
  read(11) idoubling(1:nspec)
  close(11)

  ! builds jacobian
  write(m_file,'(a,i6.6,a)') &
      '/tigress-hsm/hejunzhu/2011EUROPE_ITERATION_UPDATE/EUROPE_TOPOLOGY_FILE/proc',myrank,'_reg1_solver_data_1.bin'
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
  volume_glob = 0.0
  integral_bulk = 0._CUSTOM_REAL
  integral_betav = 0._CUSTOM_REAL
  integral_betah = 0._CUSTOM_REAL
  integral_eta = 0._CUSTOM_REAL
  norm_bulk = 0._CUSTOM_REAL
  norm_betav = 0._CUSTOM_REAL
  norm_betah = 0._CUSTOM_REAL
  norm_eta = 0._CUSTOM_REAL
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

end subroutine compute_volume

!
!-------------------------------------------------------------------------------------------------
!

subroutine get_gradient()

! calculates gradient by steepest descent method

  use model_update_tiso
  implicit none
  ! local parameters
  ! ------------------------------------------------------------------------
  ! sets maximum update in this depth range
  logical,parameter :: use_depth_maximum = .false.
  ! normalized radii
  real(kind=CUSTOM_REAL),parameter :: R_top = (6371.0 - 50.0 ) / R_EARTH_KM ! shallow depth
  real(kind=CUSTOM_REAL),parameter :: R_bottom = (6371.0 - 600.0 ) / R_EARTH_KM ! deep depth
  real(kind=CUSTOM_REAL):: r,depth_max
  ! ------------------------------------------------------------------------

  ! initializes kernel maximum
  max = 0._CUSTOM_REAL

  ! gradient in negative direction for steepest descent
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

            ! for bulk
            model_dbulk(i,j,k,ispec) =  kernel_bulk(i,j,k,ispec) ! no negative sign, in conjugate direction subroutine

            ! for shear
            model_dbetav(i,j,k,ispec) =  kernel_betav(i,j,k,ispec)
            model_dbetah(i,j,k,ispec) =  kernel_betah(i,j,k,ispec)

            ! for eta
            model_deta(i,j,k,ispec) =  kernel_eta(i,j,k,ispec)

            ! determines maximum kernel betav value within given radius
            if( use_depth_maximum ) then
              ! get radius of point
              iglob = ibool(i,j,k,ispec)
              r = sqrt( x(iglob)*x(iglob) + y(iglob)*y(iglob) + z(iglob)*z(iglob) )

              ! stores maximum kernel betav value in this depth slice, since betav is most likely dominating
              if( r < R_top .and. r > R_bottom ) then
                ! kernel betav value
                max_vsv = abs( kernel_betav(i,j,k,ispec) )
                if( max < max_vsv ) then
                  max = max_vsv
                  depth_max = r
                endif
              endif
            endif

        enddo
      enddo
    enddo
  enddo

!> Hejun Zhu
  ! stores model_dbulk, ... arrays
!  call store_kernel_updates()
!< Hejun Zhu

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
    print*,'initial gradients:'
    print*,'  bulk min/max : ',min_bulk,max_bulk
    print*,'  betav min/max: ',min_vsv,max_vsv
    print*,'  betah min/max: ',min_vsh,max_vsh
    print*,'  eta min/max  : ',min_eta,max_eta
    print*
  endif

  ! determines maximum kernel betav value within given radius
  if( use_depth_maximum ) then
    ! maximum of all processes stored in max_vsv
    call mpi_reduce(max,max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
    max = max_vsv
    depth_max = 6371.0 *( 1.0 - depth_max )
  endif

  ! determines step length
  ! based on maximum gradient value (either vsv or vsh)
  if( myrank == 0 ) then

    ! determines maximum kernel betav value within given radius
    if( use_depth_maximum ) then
        print*,'  using depth maximum between 50km and 100km: ',max
        print*,'  approximate depth maximum: ',depth_max
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

end subroutine get_gradient

!
!-------------------------------------------------------------------------------------------------
!

subroutine store_kernel_updates()

! file output for new model

  use model_update_tiso
  implicit none

  ! kernel updates
  fname = 'dbulk_c'
  write(m_file,'(a,i6.6,a)') trim(input_kernel)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_dbulk
  close(12)

  fname = 'dbetav'
  write(m_file,'(a,i6.6,a)') trim(input_kernel)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_dbetav
  close(12)

  fname = 'dbetah'
  write(m_file,'(a,i6.6,a)') trim(input_kernel)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_dbetah
  close(12)

  fname = 'deta'
  write(m_file,'(a,i6.6,a)') trim(input_kernel)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_deta
  close(12)

end subroutine store_kernel_updates

!
!-------------------------------------------------------------------------------------------------
!

subroutine store_new_model()

! file output for new model

  use model_update_tiso
  implicit none

  ! vpv model
  call mpi_reduce(maxval(model_vpv_new),max_vpv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vpv_new),min_vpv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'vpv'
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vpv_new
  close(12)

  ! vph model
  call mpi_reduce(maxval(model_vph_new),max_vph,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vph_new),min_vph,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'vph'
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vph_new
  close(12)

  ! vsv model
  call mpi_reduce(maxval(model_vsv_new),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vsv_new),min_vsv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'vsv'
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vsv_new
  close(12)

  ! vsh model
  call mpi_reduce(maxval(model_vsh_new),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vsh_new),min_vsh,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'vsh'
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vsh_new
  close(12)

  ! eta model
  call mpi_reduce(maxval(model_eta_new),max_eta,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_eta_new),min_eta,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'eta'
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_eta_new
  close(12)

  ! rho model
  call mpi_reduce(maxval(model_rho_new),max_rho,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_rho_new),min_rho,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'rho'
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
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

  use model_update_tiso
  implicit none
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: total_model

  ! vpv relative perturbations
  ! logarithmic perturbation: log( v_new) - log( v_old) = log( v_new / v_old )
  total_model = 0.0_CUSTOM_REAL
  where( model_vpv /= 0.0 ) total_model = log( model_vpv_new / model_vpv)
  ! or
  ! linear approximation: (v_new - v_old) / v_old
  !where( model_vpv /= 0.0 ) total_model = ( model_vpv_new - model_vpv) / model_vpv

  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_dvpvvpv.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_vpv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_vpv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! vph relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vph /= 0.0 ) total_model = log( model_vph_new / model_vph)
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_dvphvph.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_vph,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_vph,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! vsv relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vsv /= 0.0 ) total_model = log( model_vsv_new / model_vsv)
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_dvsvvsv.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_vsv,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_vsv,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! vsh relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vsh /= 0.0 ) total_model = log( model_vsh_new / model_vsh)
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_dvshvsh.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_vsh,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_vsh,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! eta relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_eta /= 0.0 ) total_model = log( model_eta_new / model_eta)
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_detaeta.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_eta,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_eta,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! rho relative model perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_rho /= 0.0 ) total_model = log( model_rho_new / model_rho)
  write(m_file,'(a,i6.6,a)') trim(output_model)//'/proc',myrank,'_reg1_drhorho.bin'
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
