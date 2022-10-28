program random_model

  implicit none

  integer,parameter :: NGLLX = 5, NGLLY = 5, NGLLZ = 5
  integer,parameter :: IOUT = 20
  integer,parameter :: CUSTOM_REAL = 4

  character(len=512),parameter :: LOCAL_PATH='OUTPUT_FILES/DATABASES_MPI/'

  ! Gaussian pertubation: center point
  real(kind=CUSTOM_REAL),parameter :: CENTER_X = 1320.0
  real(kind=CUSTOM_REAL),parameter :: CENTER_Y = 1320.0
  real(kind=CUSTOM_REAL),parameter :: CENTER_Z = -720.0   ! or slightly closer to top: -500.0

  ! wavelength of perturbation
  ! lambda = 0.5 * 2000.0    ! dominant wavelength: 0.5 s from source, 2000 m/s vp
  ! lambda = 0.1 * 2000.0    ! dominant wavelength: 0.1 s from source, 2000 m/s vp
  real(kind=CUSTOM_REAL),parameter :: WAVELENGTH = 0.1 * 2000.0 * 4.0      ! adding factor 4 for better adjoint source

  ! model perturbation flags
  logical, parameter :: DO_PERTURB_VP = .true.  ! default, only Vp velocities will be perturbed
  logical, parameter :: DO_PERTURB_VS = .false.
  logical, parameter :: DO_PERTURB_RHO = .false.

  integer :: NPROC
  integer :: myrank,nspec,nglob
  double precision :: percent
  character(len=256) :: prname,arg,filename
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read
  real, dimension(:,:,:,:),allocatable :: pert_param

  ! random seed
  integer :: n
  integer,dimension(:),allocatable :: myseed
  !real :: tmpharvest

  ! mesh coordinates
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore, ystore, zstore
  integer, dimension(:,:,:,:),allocatable :: ibool
  real(kind=CUSTOM_REAL) :: x,y,z
  real(kind=CUSTOM_REAL) :: dist_v,dist_h,lambda,sigma_h,sigma_v,sigma_h2,sigma_v2
  !real(kind=CUSTOM_REAL) :: norm_h,norm_v,norm

  integer :: i,j,k,ispec,iglob,idummy,ier

  !! input parameters
  if (iargc() /= 2) &
    stop 'Usage: ./xrandom_model percent NPROC [percent must be small enough (~1d-5) for F*dm=S(m+dm)-S(m) to be valid]'

  call get_command_argument(1, arg); read(arg,*,iostat=ier) percent;   if (ier /= 0) stop 'Error reading percent'
  call get_command_argument(2, arg); read(arg,*,iostat=ier)   NPROC;   if (ier /= 0) stop 'Error reading NPROC'

  print *,'random model generation:'
  print *,'  perturbation = ',percent
  print *,'  NPROC        = ',NPROC
  print *

  if (DO_PERTURB_RHO) print *,'  perturbing: rho values'
  if (DO_PERTURB_VP)  print *,'  perturbing: Vp values'
  if (DO_PERTURB_VS)  print *,'  perturbing: Vs values'
  print *

  ! loop over processes
  do myrank = 0,NPROC-1

    ! processors name
    write(prname,'(a,i6.6,a)') trim(LOCAL_PATH)//'proc',myrank,'_'
    prname = adjustl(prname)

    ! mesh
    filename = trim(prname)//'external_mesh.bin'
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if ( ier /= 0 ) then
      print *,'Error: could not open file ',trim(filename)
      stop 'Error opening database'
    endif

    ! nspec & nglob
    read(IOUT) nspec
    read(IOUT) nglob
    read(IOUT) idummy ! skip nspec_irregular

    ! user output
    if (myrank == 0) then
      print *,'mesh:'
      print *,'  nspec = ',nspec
      print *,'  nglob = ',nglob
      print *
    endif

    ! allocates arrays
    allocate(ibool(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if ( ier /= 0 ) stop 'Error allocating array ibool'
    allocate(xstore(nglob),ystore(nglob),zstore(nglob),stat=ier)
    if ( ier /= 0 ) stop 'Error allocating array xstore etc.'

    ! ibool file
    read(IOUT) ibool
    ! global point arrays
    read(IOUT) xstore
    read(IOUT) ystore
    read(IOUT) zstore
    close(IOUT)

    ! rho
    filename = trim(prname)//'rho.bin'
    allocate( rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array rho_read'
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file rho.bin'
    read(IOUT) rho_read
    close(IOUT)

    ! vp
    filename = trim(prname)//'vp.bin'
    allocate( vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array vp_read'
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file vp.bin'
    read(IOUT) vp_read
    close(IOUT)

    ! vs
    filename = trim(prname)//'vs.bin'
    allocate( vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array vs_read'
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file vs.bin'
    read(IOUT) vs_read
    close(IOUT)

    ! model perturbation
    allocate( pert_param(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array pert_param'

    !daniel: this will randomly perturb every GLL point in the model, thus adding like white noise to it.
    !           question: is the signal sensitive to this perturbation? or is it within numerical noise/artefacts?
    !
    ! perturb model randomly
    if ( .false. ) then
      !daniel: change to re-initialize seed with fixed value to make successive runs repeatable
      !CALL RANDOM_SEED()
      ! or
      CALL RANDOM_SEED(size=n)
      allocate(myseed(n))
      myseed(1:n) = myrank * 75347
      CALL RANDOM_SEED(put=myseed)

      ! this should return the same number for repeated runs...
      !call random_number(tmpharvest)
      !print *,'seed size',n
      !print *,'random number: ',tmpharvest

      if (myrank == 0) then
        print *,'random perturbation:'
        print *,'  seed size: ',n
        print *,'  perturbation size = ',percent
        print *
      endif

      ! rho
      if (DO_PERTURB_RHO) then
        CALL RANDOM_NUMBER(pert_param)
        pert_param = pert_param/maxval(abs(pert_param))*2.0-1.0
        rho_read = rho_read*(1.0+percent*pert_param)
      endif

      ! vp
      if (DO_PERTURB_VP) then
        CALL RANDOM_NUMBER(pert_param)
        pert_param = pert_param/maxval(abs(pert_param))*2.0-1.0 ! scales between [-1,1]
        vp_read = vp_read*(1.0+percent*pert_param)
      endif

      ! vs
      if (DO_PERTURB_VS) then
        CALL RANDOM_NUMBER(pert_param)
        pert_param = pert_param/maxval(abs(pert_param))*2.0-1.0
        vs_read = vs_read*(1.0+percent*pert_param)
      endif
    endif

    ! adds a Gaussian perturbation in the middle of the model
    if ( .true. ) then
      ! initializes perturbations
      pert_param(:,:,:,:) = 0.0

      !standard variance
      lambda = WAVELENGTH

      sigma_h = 0.25 * lambda/sqrt(8.0) ! such that scalelength becomes 1/4 of dominant wavelenght
      sigma_v = 0.25 * lambda/sqrt(8.0)
      ! factor two for Gaussian distribution with standard variance sigma
      sigma_h2 = 2.0 * sigma_h * sigma_h
      sigma_v2 = 2.0 * sigma_v * sigma_v

      ! scalelength: approximately S ~ sigma * sqrt(8.0) for a Gaussian smoothing
      if (myrank == 0) then
        print *,'Gaussian perturbation:'
        print *,'  input wavelength                       (m): ',lambda
        print *,'  using scalelengths horizontal,vertical (m): ',sigma_h*sqrt(8.0),sigma_v*sqrt(8.0)
        print *
        print *,'  pertubation center location : x/y/z = ',CENTER_X,'/',CENTER_Y,'/',CENTER_Z
        print *,'  perturbation size           : ',percent
        if (DO_PERTURB_RHO) print *,'  perturbing                  : rho values'
        if (DO_PERTURB_VP)  print *,'  perturbing                  : Vp values'
        if (DO_PERTURB_VS)  print *,'  perturbing                  : Vs values'
        print *
        print *,'rank ',myrank,':'
        print *,'  model dimension: x min/max = ',minval(xstore),maxval(xstore)
        print *,'                   y min/max = ',minval(ystore),maxval(ystore)
        print *,'                   z min/max = ',minval(zstore),maxval(zstore)
        print *
        print *
      endif
      ! theoretic normal value
      ! (see integral over -inf to +inf of exp[- x*x/(2*sigma) ] = sigma * sqrt(2*pi) )
      ! note: smoothing is using a Gaussian (ellipsoid for sigma_h /= sigma_v),
      !norm_h = 2.0*PI*sigma_h**2
      !norm_v = sqrt(2.0*PI) * sigma_v
      !norm   = norm_h * norm_v

      ! sets Gaussian perturbation into the middle of model:
      ! dimension (Width x Length x Depth) : 2640.0 m x 2640.0 m x 1.44 km
      do ispec = 1,nspec
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              ! GLL point location (given in m: dimension 2640 m x 2640 x x 1440 m)
              iglob = ibool(i,j,k,ispec)
              x = xstore(iglob)
              y = ystore(iglob)
              z = zstore(iglob)

              ! vertical distance to center: at - 500 m depth
              dist_v = sqrt( (CENTER_Z - z)*(CENTER_Z - z) )
              ! horizontal distance to center: at 1320 x 1320 m
              dist_h = sqrt( (CENTER_X - x)*(CENTER_X -x) + (CENTER_Y - y)*(CENTER_Y - y) )
              ! Gaussian function:  values between [0,1]
              pert_param(i,j,k,ispec) = exp( - (dist_h*dist_h) / sigma_h2 - (dist_v*dist_v) / sigma_v2 )

              !if (myrank == 0 )print *,pert_param(i,j,k,ispec),x,y,z,dist_v,dist_h
            enddo
          enddo
        enddo
      enddo

      ! adds positive perturbation to model:
      if (DO_PERTURB_RHO) rho_read = rho_read * (1.0 + percent * pert_param)
      if (DO_PERTURB_VP) vp_read  = vp_read * (1.0 + percent * pert_param)
      if (DO_PERTURB_VS) vs_read  = vs_read * (1.0 + percent * pert_param)
    endif

    ! store perturbed model
    ! rho
    open(unit=IOUT,file=trim(prname)//'rho.bin',status='old',action='write',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file rho.bin for writing'
    write(IOUT) rho_read
    close(IOUT)
    ! vp
    open(unit=IOUT,file=trim(prname)//'vp.bin',status='old',action='write',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file vp.bin for writing'
    write(IOUT) vp_read
    close(IOUT)
    ! vs
    open(unit=IOUT,file=trim(prname)//'vs.bin',status='old',action='write',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file vs.bin for writing'
    write(IOUT) vs_read
    close(IOUT)

    print *,'perturbed model files written to: '
    print *,'  ',trim(prname)//'rho.bin'
    print *,'  ',trim(prname)//'vp.bin'
    print *,'  ',trim(prname)//'vs.bin'
    print *

    deallocate(ibool,xstore,ystore,zstore)
    deallocate(rho_read,vp_read,vs_read,pert_param)

  enddo ! NPROC

  print *,'done'

end program random_model
