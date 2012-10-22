program random_model

  implicit none

  include "mpif.h"

  integer,parameter :: NGLLX=5,NGLLY=5,NGLLZ=5,IOUT=20
  integer,parameter :: CUSTOM_REAL = 4
  character(len=512),parameter :: LOCAL_PATH='../OUTPUT_FILES/DATABASES_MPI/'

  integer :: myrank,ier,nspec,nglob,NPROC,ios
  double precision :: percent
  character(len=256) prname,arg
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read,random

  integer :: n
  integer,dimension(:),allocatable :: myseed
  !real :: tmpharvest

  ! mesh coordinates
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore, ystore, zstore
  integer, dimension(:,:,:,:),allocatable :: ibool
  real(kind=CUSTOM_REAL) :: x,y,z
  real(kind=CUSTOM_REAL) :: dist_v,dist_h,lambda,sigma_h,sigma_v,sigma_h2,sigma_v2
  !real(kind=CUSTOM_REAL) :: norm_h,norm_v,norm
  real(kind=CUSTOM_REAL),parameter :: PI = 3.1415926535897931
  integer :: i,j,k,ispec,iglob

  call MPI_Init(ier)
  call MPI_Comm_Rank(MPI_COMM_WORLD,myrank,ier)
  call MPI_Comm_Size(MPI_COMM_WORLD,NPROC,ier)

  !! input parameters
  if( iargc() .ne. 1 ) stop 'Usage: ./xrandom_model percent [must be small enough (~1d-5) for F*dm=S(m+dm)-S(m) to be valid]'
  j=1;  call getarg(j, arg); read(arg,*,iostat=ios) percent;   if (ios /= 0) stop 'Error reading percent'

  ! processors name
  write(prname,'(a,i6.6,a)') trim(LOCAL_PATH)//'proc',myrank,'_'
  ! nspec & nglob
  open(unit=IOUT,file=trim(adjustl(prname))//'external_mesh.bin',status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) stop 'error opening database proc######_external_mesh.bin'
  read(IOUT) nspec
  read(IOUT) nglob

  ! ibool file
  allocate(ibool(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibool'
  read(IOUT) ibool
  ! global point arrays
  allocate(xstore(nglob),ystore(nglob),zstore(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocating array xstore etc.'
  read(IOUT) xstore
  read(IOUT) ystore
  read(IOUT) zstore

  close(IOUT)

  ! rho
  allocate( rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array rho_read'
  open(unit=IOUT,file=trim(adjustl(prname))//'rho.bin',status='old',action='read',form='unformatted')
  read(IOUT) rho_read
  close(IOUT)

  ! vp
  allocate( vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array vp_read'
  open(unit=IOUT,file=trim(adjustl(prname))//'vp.bin',status='old',action='read',form='unformatted')
  read(IOUT) vp_read
  close(IOUT)
  ! vs
  allocate( vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array vs_read'
  open(unit=IOUT,file=trim(adjustl(prname))//'vs.bin',status='old',action='read',form='unformatted')
  read(IOUT) vs_read
  close(IOUT)

!-------------------------------------------------
!daniel: this will randomly perturb every GLL point in the model, thus adding like white noise to it.
!           question: is the signal sensitive to this perturbation? or is it within numerical noise/artefacts?

  ! perturb model randomly
  allocate( random(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array random'

  if( .false. ) then

    CALL RANDOM_SEED()

    !daniel: re-initialize seed with fixed value to make successive runs repeatable
    CALL RANDOM_SEED(size=n)
    allocate(myseed(n))
    myseed(1:n) = myrank * 75347
    CALL RANDOM_SEED(put=myseed)

    ! this should return the same number for repeated runs...
    !call random_number(tmpharvest)
    !print *,'seed size',n
    !print *,'random number: ',tmpharvest

    CALL RANDOM_NUMBER(random)
    random=random/maxval(abs(random))*2.0-1.0
    rho_read=rho_read*(1.0+percent*random)

    CALL RANDOM_NUMBER(random)
    random=random/maxval(abs(random))*2.0-1.0
    vp_read= vp_read*(1.0+percent*random)

    CALL RANDOM_NUMBER(random)
    random=random/maxval(abs(random))*2.0-1.0
    vs_read= vs_read*(1.0+percent*random)

  endif


! adds a gaussian perturbation in the middle of the model
  if( .true. ) then
    ! initializes perturbations
    random(:,:,:,:) = 0.0

    !standard variance
    lambda = 0.5 * 2000.0    ! dominant wavelength: 0.5 s from source, 2000 m/s vp
    sigma_h = 0.25*lambda/sqrt(8.0) ! such that scalelength becomes 1/4 of dominant wavelenght
    sigma_v = 0.25*lambda/sqrt(8.0)
    ! factor two for gaussian distribution with standard variance sigma
    sigma_h2 = 2.0 * sigma_h * sigma_h
    sigma_v2 = 2.0 * sigma_v * sigma_v
    ! scalelength: approximately S ~ sigma * sqrt(8.0) for a gaussian smoothing
    if(myrank == 0 )print*,"  scalelengths horizontal,vertical (m): ",sigma_h*sqrt(8.0),sigma_v*sqrt(8.0)

    ! theoretic normal value
    ! (see integral over -inf to +inf of exp[- x*x/(2*sigma) ] = sigma * sqrt(2*pi) )
    ! note: smoothing is using a gaussian (ellipsoid for sigma_h /= sigma_v),
    !norm_h = 2.0*PI*sigma_h**2
    !norm_v = sqrt(2.0*PI) * sigma_v
    !norm   = norm_h * norm_v

    ! sets gaussian perturbation into the middle of model:
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
            dist_v = sqrt( (-500.0 - z)*(-500.0 - z) )
            ! horizontal distance to center: at 1320 x 1320 m
            dist_h = sqrt( (1320.0 - x)*(1320.0 -x) + (1320.0 - y)*(1320.0 - y) )
            ! gaussian function:  values between [0,1]
            random(i,j,k,ispec) = exp( - (dist_h*dist_h) / sigma_h2 - (dist_v*dist_v) / sigma_v2 )

            !if(myrank == 0 )print*,random(i,j,k,ispec),x,y,z,dist_v,dist_h
          enddo
        enddo
      enddo
    enddo

    ! adds positive perturbation to model:
    !rho_read = rho_read*(1.0+percent*random)
    vp_read  = vp_read*(1.0+percent*random)
    !vs_read  = vs_read*(1.0+percent*random)

  endif

!-------------------------------------------------

  ! store perturbed model
  ! rho
  open(unit=IOUT,file=trim(adjustl(prname))//'rho.bin',status='old',action='write',form='unformatted')
  write(IOUT) rho_read
  close(IOUT)
  ! vp
  open(unit=IOUT,file=trim(adjustl(prname))//'vp.bin',status='old',action='write',form='unformatted')
  write(IOUT) vp_read
  close(IOUT)
  ! vs
  open(unit=IOUT,file=trim(adjustl(prname))//'vs.bin',status='old',action='write',form='unformatted')
  write(IOUT) vs_read
  close(IOUT)

  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  call MPI_Finalize(ier)

end program random_model
