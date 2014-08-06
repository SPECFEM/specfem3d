! add_model_globe_iso
!
! this program can be used to update ISOTROPIC model files with
! (smoothed & summed) event kernels.
! the kernels are given for isotropic parameters (alpha,beta,rho) or ( bulk_c,beta,rho).
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
!       proc000***_reg1_vs.bin &
!       proc000***_reg1_vp.bin &
!       proc000***_reg1_rho.bin
!
!- INPUT_GRADIENT/ contains:
!       proc000***_reg1_bulk_c_kernel_smooth.bin &
!       proc000***_reg1_bulk_beta_kernel_smooth.bin &
!       proc000***_reg1_rho_kernel_smooth.bin
!     or
!       proc000***_reg1_alpha_kernel_smooth.bin &
!       proc000***_reg1_beta_kernel_smooth.bin &
!       proc000***_reg1_rho_kernel_smooth.bin
!
!- topo/ contains:
!       proc000***_reg1_solver_data_1.bin
!
! new models are stored in
!- OUTPUT_MODEL/ as
!   proc000***_reg1_vp_new.bin and
!   proc000***_reg1_vs_new.bin and
!   proc000***_reg1_rho_new.bin and
!
! USAGE: e.g. ./add_model_globe_iso 0.3

module model_update_iso

  include 'mpif.h'
  include 'constants_globe.h'
  include 'precision_globe.h'
  include 'values_from_mesher_globe.h'

  ! ======================================================

  ! USER PARAMETERS

  ! by default, this algorithm uses (bulk,bulk_beta,rho) kernels to update vp,vs,rho
  ! if you prefer using (alpha,beta,rho) kernels, set this flag to true
  logical, parameter :: USE_ALPHA_BETA_RHO = .false.

  ! ignore rho kernel, but use density perturbations as a scaling of Vs perturbations
  logical, parameter :: USE_RHO_SCALING = .false.

  ! in case of rho scaling, specifies density scaling factor with shear perturbations
  ! see e.g. Montagner & Anderson (1989), Panning & Romanowicz (2006)
  real(kind=CUSTOM_REAL),parameter :: RHO_SCALING = 0.33_CUSTOM_REAL

  ! ======================================================

  integer, parameter :: NSPEC = NSPEC_CRUST_MANTLE
  integer, parameter :: NGLOB = NGLOB_CRUST_MANTLE

  ! isotropic model files
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        model_vp,model_vs,model_rho
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        model_vp_new,model_vs_new,model_rho_new

  ! model updates
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        model_dA,model_dB,model_dR

  ! kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        kernel_a,kernel_b,kernel_rho

  ! volume
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: x, y, z
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: ibool

  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_beta,norm_rho
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_beta_sum,norm_rho_sum

  ! steepest descent lengths
  real(kind=CUSTOM_REAL) :: step_fac,step_length

  ! statistics
  real(kind=CUSTOM_REAL) :: min_vp,min_vs,max_vp,max_vs, &
    min_rho,max_rho,minmax(4),vs_sum,vp_sum,rho_sum

  real(kind=CUSTOM_REAL) :: beta1,beta0,rho1,rho0,alpha1,alpha0
  real(kind=CUSTOM_REAL) :: dbetaiso,dA

  integer :: nfile, myrank, sizeprocs,  ier
  integer :: i, j, k,ispec, iglob, ishell, n, it, j1, ib, npts_sem, ios
  character(len=150) :: sline, m_file, fname    

end module model_update_iso

!
!-------------------------------------------------------------------------------------------------
!

program add_model

  use model_update_iso

  implicit none

  ! ============ program starts here =====================

  ! initializes arrays
  call initialize()

  ! reads in parameters needed
  call read_parameters()

  ! reads in current transverse isotropic model files: vpv.. & vsv.. & eta & rho
  call read_model()

  ! reads in smoothed kernels: bulk, beta, rho
  call read_kernels()

  ! calculates gradient
  ! steepest descent method
  call get_gradient()


  ! computes new model values for alpha, beta and rho
  ! and stores new model files

  ! model update:
  !   isotropic update everywhere
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! initial model values
          beta0 = model_vs(i,j,k,ispec)
          rho0 = model_rho(i,j,k,ispec)
          alpha0 = model_vp(i,j,k,ispec)

          beta1 = 0._CUSTOM_REAL
          rho1 = 0._CUSTOM_REAL
          alpha1 = 0._CUSTOM_REAL

          ! isotropic model update

          ! shear values
          dbetaiso = model_dB(i,j,k,ispec)
          beta1 = beta0 * exp( dbetaiso )

          ! density
          rho1 = rho0 * exp( model_dR(i,j,k,ispec) )

          ! alpha values
          dA = model_dA(i,j,k,ispec)
          if( USE_ALPHA_BETA_RHO ) then
            ! new vp values use alpha model update
            alpha1 = alpha0 * exp( dA )
          else
            ! new vp values use bulk model update:
            ! this is based on vp_new = sqrt( bulk_new**2 + 4/3 vs_new**2 )
            alpha1 = sqrt( alpha0**2 * exp(2.0*dA) + FOUR_THIRDS * beta0**2 * ( &
                              exp(2.0*dbetaiso) - exp(2.0*dA) ) )
          endif

          ! stores new model values
          model_vp_new(i,j,k,ispec) = alpha1
          model_vs_new(i,j,k,ispec) = beta1
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

  use model_update_iso
  implicit none

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if( sizeprocs /= nchunks_val*nproc_xi_val*nproc_eta_val ) then
    print*,'sizeprocs:',sizeprocs,nchunks_val,nproc_xi_val,nproc_eta_val
    call exit_mpi(myrank,'error number sizeprocs')
  endif

  model_vp = 0.0_CUSTOM_REAL
  model_vs = 0.0_CUSTOM_REAL
  model_rho = 0.0_CUSTOM_REAL
  
  model_dA = 0.0_CUSTOM_REAL
  model_dB = 0.0_CUSTOM_REAL
  model_dR = 0.0_CUSTOM_REAL
  
  kernel_a = 0.0_CUSTOM_REAL
  kernel_b = 0.0_CUSTOM_REAL
  kernel_rho = 0.0_CUSTOM_REAL

end subroutine initialize

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_parameters()

! reads in parameters needed

  use model_update_iso
  implicit none
  character(len=150) :: s_step_fac

  !--------------------------------------------------------

  ! subjective step length to multiply to the gradient
  ! e.g. step_fac = 0.03

  call getarg(1,s_step_fac)

  if (trim(s_step_fac) == '') then
    call exit_MPI(myrank,'Usage: add_model_globe_iso step_factor')
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

    !open(20,file='step_fac',status='unknown',action='write')
    !write(20,'(1e24.12)') step_fac
    !close(20)
  endif

end subroutine read_parameters

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_model()

! reads in current transverse isotropic model: vpv.. & vsv.. & eta & rho

  use model_update_iso
  implicit none


  ! reads in current vp & vs & rho model
  ! vp model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vp.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_vp(:,:,:,1:nspec)
  close(12)

  ! vs model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vs.bin'
  open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(12) model_vs(:,:,:,1:nspec)
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
  call mpi_reduce(minval(model_vp),min_vp,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vp),max_vp,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_vs),min_vs,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_vs),max_vs,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_rho),min_rho,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_rho),max_rho,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'initial models:'
    print*,'  vs min/max: ',min_vs,max_vs
    print*,'  vp min/max: ',min_vp,max_vp
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif

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
  close(11)

end subroutine read_model

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_kernels()

! reads in smoothed kernels: bulk, betav, betah, eta

  use model_update_iso
  implicit none



  ! reads in smoothed (& summed) event kernel
  if( USE_ALPHA_BETA_RHO ) then
    ! reads in alpha kernel
    fname = 'alpha_kernel_smooth'
  else
    ! reads in bulk_c kernel
    fname = 'bulk_c_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_'//trim(fname)//'.bin'
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
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_'//trim(fname)//'.bin'
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
    write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_rho_kernel_smooth.bin'
    open(12,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(12) kernel_rho(:,:,:,1:nspec)
    close(12)
  endif

  ! statistics
  call mpi_reduce(minval(kernel_a),min_vp,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_a),max_vp,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(kernel_b),min_vs,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_b),max_vs,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(kernel_rho),min_rho,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(kernel_rho),max_rho,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'initial kernels:'
    print*,'  beta min/max: ',min_vs,max_vs
    print*,'  alpha min/max: ',min_vp,max_vp
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif

end subroutine read_kernels

!
!-------------------------------------------------------------------------------------------------
!

subroutine get_gradient()

! calculates gradient by steepest descent method

  use model_update_iso
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL):: r,max,depth_max

  ! ------------------------------------------------------------------------

  ! sets maximum update in this depth range
  logical,parameter :: use_depth_maximum = .true.
  ! normalized radii
  real(kind=CUSTOM_REAL),parameter :: R_top = (6371.0 - 50.0 ) / R_EARTH_KM ! shallow depth
  real(kind=CUSTOM_REAL),parameter :: R_bottom = (6371.0 - 100.0 ) / R_EARTH_KM ! deep depth

  ! ------------------------------------------------------------------------

  ! initializes kernel maximum
  max = 0._CUSTOM_REAL

  ! calculates gradient
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

            ! determines maximum kernel beta value within given radius
            if( use_depth_maximum ) then
              ! get radius of point
              iglob = ibool(i,j,k,ispec)
              r = sqrt( x(iglob)*x(iglob) + y(iglob)*y(iglob) + z(iglob)*z(iglob) )

              ! stores maximum kernel betav value in this depth slice, since betav is most likely dominating
              if( r < R_top .and. r > R_bottom ) then
                ! shear kernel value
                max_vs = abs( kernel_b(i,j,k,ispec) )
                if( max < max_vs ) then
                  max = max_vs
                  depth_max = r
                endif
              endif
            endif

        enddo
      enddo
    enddo
  enddo

  ! statistics
  call mpi_reduce(minval(model_dA),min_vp,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dA),max_vp,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_dB),min_vs,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dB),max_vs,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  call mpi_reduce(minval(model_dR),min_rho,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(maxval(model_dR),max_rho,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'initial gradients:'
    print*,'  a min/max: ',min_vp,max_vp
    print*,'  beta min/max : ',min_vs,max_vs
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif

  ! determines maximum kernel betav value within given radius
  if( use_depth_maximum ) then
    ! maximum of all processes stored in max_vsv
    call mpi_reduce(max,max_vs,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
    max = max_vs
    depth_max = 6371.0 *( 1.0 - depth_max )
  endif

  ! master determines step length based on maximum gradient value (either vp or vs)
  if( myrank == 0 ) then
  
      ! determines maximum kernel betav value within given radius
    if( use_depth_maximum ) then
      print*,'  using depth maximum between 50km and 100km: ',max
      print*,'  approximate depth maximum: ',depth_max
      print*
    else  
      ! maximum gradient values
      minmax(1) = abs(min_vs)
      minmax(2) = abs(max_vs)
      minmax(3) = abs(min_vp)
      minmax(4) = abs(max_vp)

      ! maximum value of all kernel maxima
      max = maxval(minmax)
      print*,'  using maximum: ',max
      print*
    endif
    
    ! chooses step length such that it becomes the desired, given step factor as inputted
    step_length = step_fac/max

    print*,'  step length : ',step_length,max
    print*
  endif
  call mpi_bcast(step_length,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)


  ! gradient length
  max_vp = sum( model_dA * model_dA )
  max_vs = sum( model_dB * model_dB )
  max_rho = sum( model_dR * model_dR )

  call mpi_reduce(max_vp,vp_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(max_vs,vs_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(max_rho,rho_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  max_vp = sqrt(vp_sum)
  max_vs = sqrt(vs_sum)
  max_rho = sqrt(rho_sum)

  if( myrank == 0 ) then
    print*,'norm model updates:'  
    print*,'  initial a length: ',max_vp
    print*,'  initial beta length : ',max_vs
    print*,'  initial rho length: ',max_rho
    print*
  endif

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

  if( myrank == 0 ) then
    print*,'scaled gradients:'
    print*,'  a min/max: ',min_vp,max_vp
    print*,'  beta min/max : ',min_vs,max_vs
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif

end subroutine get_gradient

!
!-------------------------------------------------------------------------------------------------
!

subroutine store_new_model()

! file output for new model

  use model_update_iso
  implicit none

  ! vp model
  call mpi_reduce(maxval(model_vp_new),max_vp,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vp_new),min_vp,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'vp_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vp_new
  close(12)

  ! vs model
  call mpi_reduce(maxval(model_vs_new),max_vs,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(model_vs_new),min_vs,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)
  fname = 'vs_new'
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) model_vs_new
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
    print*,'  vp min/max: ',min_vp,max_vp
    print*,'  vs min/max: ',min_vs,max_vs
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif


end subroutine store_new_model


!
!-------------------------------------------------------------------------------------------------
!

subroutine store_perturbations()

! file output for new model

  use model_update_iso
  implicit none
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: total_model

  ! vp relative perturbations
  ! logarithmic perturbation: log( v_new) - log( v_old) = log( v_new / v_old )
  total_model = 0.0_CUSTOM_REAL
  where( model_vp /= 0.0 ) total_model = log( model_vp_new / model_vp)
  ! or
  ! linear approximation: (v_new - v_old) / v_old
  !where( model_vp /= 0.0 ) total_model = ( model_vp_new - model_vp) / model_vp

  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvpvp.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_vp,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_vp,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  ! vs relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vs /= 0.0 ) total_model = log( model_vs_new / model_vs)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvsvs.bin'
  open(12,file=trim(m_file),form='unformatted',action='write')
  write(12) total_model
  close(12)
  call mpi_reduce(maxval(total_model),max_vs,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(minval(total_model),min_vs,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

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
    print*,'  dvp/vp min/max: ',min_vp,max_vp
    print*,'  dvs/vs min/max: ',min_vs,max_vs
    print*,'  drho/rho min/max: ',min_rho,max_rho
    print*
  endif

end subroutine store_perturbations

!
!-------------------------------------------------------------------------------------------------
!

subroutine compute_volume()

! computes volume element associated with points

  use model_update_iso
  implicit none
  ! jacobian
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: jacobian
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl, &
    jacobianl,volumel

  ! integration values
  real(kind=CUSTOM_REAL) :: kernel_integral_alpha,kernel_integral_beta,kernel_integral_rho
  real(kind=CUSTOM_REAL) :: integral_alpha_sum,integral_beta_sum,integral_rho_sum
  
  real(kind=CUSTOM_REAL) :: volume_glob,volume_glob_sum
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

        enddo
      enddo
    enddo
  enddo

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
          kernel_integral_alpha = kernel_integral_alpha &
                                 + volumel * kernel_a(i,j,k,ispec)

          kernel_integral_beta = kernel_integral_beta &
                                 + volumel * kernel_b(i,j,k,ispec)

          kernel_integral_rho = kernel_integral_rho &
                                 + volumel * kernel_rho(i,j,k,ispec)

          ! gradient vector norm sqrt(  v^T * v )
          norm_bulk = norm_bulk + kernel_a(i,j,k,ispec)**2
          norm_beta = norm_beta + kernel_b(i,j,k,ispec)**2
          norm_rho = norm_rho + kernel_rho(i,j,k,ispec)**2

          ! checks number
          if( isNaN(kernel_integral_alpha) ) then
            print*,'error NaN: ',kernel_integral_alpha
            print*,'rank:',myrank
            print*,'i,j,k,ispec:',i,j,k,ispec
            print*,'volumel: ',volumel,'kernel_bulk:',kernel_a(i,j,k,ispec)
            call exit_MPI(myrank,'error NaN')
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
  ! kernel integration: for whole volume
  call mpi_reduce(kernel_integral_alpha,integral_alpha_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(kernel_integral_beta,integral_beta_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(kernel_integral_rho,integral_rho_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(volume_glob,volume_glob_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  if( myrank == 0 ) then
    print*,'integral kernels:'
    print*,'  a : ',integral_alpha_sum
    print*,'  beta : ',integral_beta_sum
    print*,'  rho : ',integral_rho_sum
    print*
    print*,'  total volume:',volume_glob_sum
    print*
  endif

  ! norms: for whole volume
  call mpi_reduce(norm_bulk,norm_bulk_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_beta,norm_beta_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(norm_rho,norm_rho_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

  norm_bulk = sqrt(norm_bulk_sum)
  norm_beta = sqrt(norm_beta_sum)
  norm_rho = sqrt(norm_rho_sum)

  if( myrank == 0 ) then
    print*,'norm kernels:'
    print*,'  a : ',norm_bulk
    print*,'  beta : ',norm_beta
    print*,'  rho : ',norm_rho
    print*
  endif

  ! root-mean square
  call mpi_reduce(rms_vp,rms_vp_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(rms_vs,rms_vs_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call mpi_reduce(rms_rho,rms_rho_sum,1,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  rms_vp  = sqrt( rms_vp_sum / volume_glob_sum )
  rms_vs  = sqrt( rms_vs_sum / volume_glob_sum )
  rms_rho = sqrt( rms_rho_sum / volume_glob_sum )

  if( myrank == 0 ) then
    print*,'root-mean square of perturbations:'
    print*,'  vp : ',rms_vp
    print*,'  vs : ',rms_vs
    print*,'  rho : ',rms_rho
    print*
  endif

end subroutine compute_volume
