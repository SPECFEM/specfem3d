program random_model

  implicit none

  integer,parameter :: NGLLX=5,NGLLY=5,NGLLZ=5,IOUT=20

  integer :: myrank,ier,nspec,nglob,NPROC,j,ios,NSTEP,irec,irec_total
  real :: DATASPACE,MODELSPACE,DT
  character(len=256) prname,arg,procname,filename
  real(kind=4),dimension(:,:,:,:),allocatable :: vp,vs,rho,weights
  real(kind=4),dimension(:,:,:,:),allocatable :: vp0,vs0,rho0,krhop,kalpha,kbeta
  real(kind=4) :: r4head(60)

  real(kind=4),dimension(:), allocatable :: adj

  !! input parameters
  if( iargc() .ne. 3 ) stop 'Usage: ./xpostprocessing NSTEP DT NPROC'
  j=1;  call getarg(j, arg); read(arg,*,iostat=ios) NSTEP;   if (ios /= 0) stop 'Error reading NSTEP'
  j=2;  call getarg(j, arg); read(arg,*,iostat=ios)    DT;   if (ios /= 0) stop 'Error reading DT'
  j=3;  call getarg(j, arg); read(arg,*,iostat=ios) NPROC;   if (ios /= 0) stop 'Error reading NPROC'

  print*,'postprocessing:'
  print*,'  NSTEP:',NSTEP
  print*,'  DT:',DT
  print*,'  NPROC:',NPROC
  print*

  MODELSPACE=0.0
  DATASPACE=0.0

  irec_total = 0
  allocate(adj(NSTEP))

  do myrank=0,NPROC-1
    !!! calculate inner product in model space --- <F* F dm, dm>
    ! processors name
    write(prname,'(a,i6.6,a)') 'proc',myrank,'_'

    print*,'  ',trim(prname)

    ! nspec & nglob
    open(unit=IOUT,file='../OUTPUT_FILES/DATABASES_MPI/'//trim(adjustl(prname))//'external_mesh.bin',status='old',action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) stop 'error opening database proc######_external_mesh.bin'
    read(IOUT) nspec
    read(IOUT) nglob
    close(IOUT)

    ! weights
    allocate(weights(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array weights'
    open(unit=IOUT,file='../OUTPUT_FILES/DATABASES_MPI/'//trim(adjustl(prname))//'weights_kernel.bin',status='old',action='read',form='unformatted',iostat=ier)
    read(IOUT) weights
    close(IOUT)

    ! kernels
    allocate(krhop(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array krhop'
    open(unit=IOUT,file='../OUTPUT_FILES/DATABASES_MPI/'//trim(adjustl(prname))//'rhop_acoustic_kernel.bin',status='old',action='read',form='unformatted',iostat=ier)
    read(IOUT) krhop
    close(IOUT)

    allocate(kalpha(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array kalpha'
    open(unit=IOUT,file='../OUTPUT_FILES/DATABASES_MPI/'//trim(adjustl(prname))//'alpha_acoustic_kernel.bin',status='old',action='read',form='unformatted',iostat=ier)
    read(IOUT) kalpha
    close(IOUT)

    allocate(kbeta(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array kbeta'

    ! rho
    allocate(rho(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array rho'
    allocate(rho0(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array rho0'
    open(unit=IOUT,file='../models/target_model/'//trim(adjustl(prname))//'rho.bin',status='old',action='read',form='unformatted')
    read(IOUT) rho
    close(IOUT)
    open(unit=IOUT,file='../models/initial_model/'//trim(adjustl(prname))//'rho.bin',status='old',action='read',form='unformatted')
    read(IOUT) rho0
    close(IOUT)

    ! vp
    allocate(vp(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array vp'
    allocate(vp0(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array vp0'
    open(unit=IOUT,file='../models/target_model/'//trim(adjustl(prname))//'vp.bin',status='old',action='read',form='unformatted')
    read(IOUT) vp
    close(IOUT)
    open(unit=IOUT,file='../models/initial_model/'//trim(adjustl(prname))//'vp.bin',status='old',action='read',form='unformatted')
    read(IOUT) vp0
    close(IOUT)

    ! vs
    allocate(vs(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array vs'
    allocate(vs0(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if( ier /= 0 ) stop 'error allocating array vs0'

    ! F dm = S(m+dm) - S(m), we backpropogate syn-dat (see adj_seismogram.f90) => we have to add a minus sign in front of kernels
    MODELSPACE=MODELSPACE + &
               sum(weights*(-krhop)*(rho-rho0)/rho0)+ &
               sum(weights*(-kalpha)*(vp-vp0)/vp0)

    deallocate(rho,rho0,vp,vp0,vs,vs0,weights,krhop,kalpha,kbeta)


    !!! calculate inner product in data space --- <F dm, F dm>
    write(procname,"(i4)") myrank
    filename=trim(adjustl(procname))//"_dx_SU"
    open(111,file="../OUTPUT_FILES/SEM/"//trim(adjustl(filename))//".adj",access='direct',recl=240+4*NSTEP,iostat = ios)
    if( ios /= 0 ) stop 'error opening adjoint trace'
    print*,'  ',trim(adjustl(filename))//".adj"

    irec=1
    do while(ios==0)
       adj(:) = 0.0

       read(111,rec=irec,iostat=ios) r4head,adj
       if (ios /= 0) exit

       DATASPACE=DATASPACE+sum(adj(:)*adj(:))*DT

       irec=irec+1
    enddo
    close(111)
    irec_total = irec_total + irec

!elastic case
!    filename=trim(adjustl(procname))//"_dy_SU"
!    open(111,file="../OUTPUT_FILES/SEM/"//trim(adjustl(filename))//".adj",access='direct',recl=240+4*NSTEP,iostat = ios)
!    if( ios /= 0 ) stop 'error opening adjoint trace'
!    print*,'  ',trim(adjustl(filename))//".adj"
!
!    irec=1
!    do while(ios==0)
!       adj(:) = 0.0
!
!       read(111,rec=irec,iostat=ios) r4head,adj
!       if (ios /= 0) exit
!
!       DATASPACE=DATASPACE+sum(adj(:)*adj(:))*DT
!
!       irec=irec+1
!    enddo
!    close(111)
!
!    filename=trim(adjustl(procname))//"_dz_SU"
!    open(111,file="../OUTPUT_FILES/SEM/"//trim(adjustl(filename))//".adj",access='direct',recl=240+4*NSTEP,iostat = ios)
!    if( ios /= 0 ) stop 'error opening adjoint trace'
!    print*,'  ',trim(adjustl(filename))//".adj"
!
!    irec=1
!    do while(ios==0)
!       adj(:) = 0.0
!
!       read(111,rec=irec,iostat=ios) r4head,adj
!       if (ios /= 0) exit
!
!       DATASPACE=DATASPACE+sum(adj(:)*adj(:))*DT
!
!       irec=irec+1
!    enddo
!    close(111)

  enddo

  print*, 'DATASPACE=',DATASPACE
  print*, 'MODELSPACE=',MODELSPACE
  print*, 'relative error=',(DATASPACE-MODELSPACE)/DATASPACE,irec_total

end program random_model
