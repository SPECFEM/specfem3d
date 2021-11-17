program random_model

  implicit none

  integer,parameter :: NGLLX = 5,NGLLY = 5,NGLLZ = 5
  integer,parameter :: IOUT = 20

  integer :: NPROC,NSTEP,SIM_TYPE
  integer :: myrank,ier,nspec,nglob
  integer :: irec,irec_total,comp_total,icomp

  double precision :: DT
  double precision :: DATASPACE,MODELSPACE
  double precision :: DS_l,MS_l

  character(len=512) :: prname,arg,procname,filename,str_type,name

  real(kind=4),dimension(:,:,:,:),allocatable :: vp,vs,rho,weights
  real(kind=4),dimension(:,:,:,:),allocatable :: vp0,vs0,rho0,krhop,kalpha,kbeta
  real(kind=4) :: r4head(60)

  real(kind=4),dimension(:), allocatable :: adj

  !! input parameters
  if (iargc() /= 4) &
    stop 'Usage: ./xpostprocessing NSTEP DT NPROC type[acoustic/elastic]'

  call get_command_argument(1, arg); read(arg,*,iostat=ier) NSTEP;   if (ier /= 0) stop 'Error reading NSTEP'
  call get_command_argument(2, arg); read(arg,*,iostat=ier)    DT;   if (ier /= 0) stop 'Error reading DT'
  call get_command_argument(3, arg); read(arg,*,iostat=ier) NPROC;   if (ier /= 0) stop 'Error reading NPROC'
  call get_command_argument(4, arg); read(arg,*,iostat=ier) str_type;   if (ier /= 0) stop 'Error reading str_type'

  str_type = adjustl(str_type)
  if (trim(str_type) == "acoustic") then
    SIM_TYPE = 1
  else if (trim(str_type) == "elastic") then
    SIM_TYPE = 2
  else
    print *,'Error: Invalid type '//trim(str_type)//' must be: acoustic or elastic'
    stop 'Invalid type, must be: acoustic or elastic'
  endif

  print *,'postprocessing:'
  print *,'  NSTEP = ',NSTEP
  print *,'  DT    = ',DT
  print *,'  NPROC = ',NPROC
  print *,'  type  = ',SIM_TYPE,'(acoustic == 1 / elastic == 2)'
  print *

  MODELSPACE = 0.d0
  DATASPACE  = 0.d0

  irec_total = 0

  allocate(adj(NSTEP),stat=ier)
  if (ier /= 0) stop 'Error allocating array adj'

  do myrank = 0,NPROC-1

    print *,'rank: ',myrank

    !!! calculate inner product in model space --- < F* F dm, dm>
    print *,'  calculating inner product in model space: < F* F dm, dm>'

    ! processors name
    write(prname,'(a,i6.6,a)') 'proc',myrank,'_'
    prname = adjustl(prname)
    !print *,'  ',trim(prname)

    ! nspec & nglob
    filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'external_mesh.bin'
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if ( ier /= 0 ) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error opening database proc######_external_mesh.bin'
    endif
    read(IOUT) nspec
    read(IOUT) nglob
    close(IOUT)

    ! weights
    allocate(weights(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'error allocating array weights'
    filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'weights_kernel.bin'
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file weights_kernel.bin'
    read(IOUT) weights
    close(IOUT)

    ! kernels
    ! K_rho_prime
    allocate(krhop(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'error allocating array krhop'
    if (SIM_TYPE == 1) then
      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'rhop_acoustic_kernel.bin'
    else
      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'rhop_kernel.bin'
    endif
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file rhop_kernel.bin'
    read(IOUT) krhop
    close(IOUT)

    ! K_alpha
    allocate(kalpha(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'error allocating array kalpha'
    if (SIM_TYPE == 1) then
      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'alpha_acoustic_kernel.bin'
    else
      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'alpha_kernel.bin'
    endif
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file alpha_kernel.bin'
    read(IOUT) kalpha
    close(IOUT)

    ! K_beta, zero shear for acoustic simulation
    allocate(kbeta(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'error allocating array kbeta'
    if (SIM_TYPE == 1) then
      ! acoustic, no shear kernel
      kbeta = 0.0
    else
      ! elastic
      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'beta_kernel.bin'
      open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'Error opening file beta_kernel.bin'
      read(IOUT) kbeta
      close(IOUT)
    endif

    print *,'  kernel rhop : min/max = ',minval(krhop),maxval(krhop)
    print *,'  kernel alpha: min/max = ',minval(kalpha),maxval(kalpha)
    if (SIM_TYPE /= 1) print *,'  kernel beta : min/max = ',minval(kbeta),maxval(kbeta)
    print *

    ! rho
    allocate(rho(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array rho'
    allocate(rho0(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array rho0'
    filename = 'MODELS/target_model/'//trim(prname)//'rho.bin'
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file rho.bin target_model'
    read(IOUT) rho
    close(IOUT)
    filename = 'MODELS/initial_model/'//trim(prname)//'rho.bin'
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file rho.bin initial_model'
    read(IOUT) rho0
    close(IOUT)

    ! vp
    allocate(vp(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array vp'
    allocate(vp0(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array vp0'
    filename = 'MODELS/target_model/'//trim(prname)//'vp.bin'
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file vp.bin target_model'
    read(IOUT) vp
    close(IOUT)
    filename = 'MODELS/initial_model/'//trim(prname)//'vp.bin'
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file vp.bin initial_model'
    read(IOUT) vp0
    close(IOUT)

    ! vs
    allocate(vs(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array vs'
    allocate(vs0(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array vs0'
    if (SIM_TYPE == 1) then
      ! acoustic, zero shear velocity
      vs = 0.0
      vs0 = 0.0
    else
      filename = 'MODELS/target_model/'//trim(prname)//'vs.bin'
      open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'Error opening file vs.bin target_model'
      read(IOUT) vs
      close(IOUT)
      filename = 'MODELS/initial_model/'//trim(prname)//'vs.bin'
      open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'Error opening file vs.bin initial_model'
      read(IOUT) vs0
      close(IOUT)
    endif

    print *,'  relative model perturbation (rho - rhop)/rho0 : min/max = ',minval((rho-rho0)/rho0),maxval((rho-rho0)/rho0)
    print *,'  relative model perturbation (vp - vp0)/vp0    : min/max = ',minval((vp-vp0)/vp0),maxval((vp-vp0)/vp0)
    if (SIM_TYPE == 2) &
      print *,'  relative model perturbation (vs - vs0)/vs0    : min/max = ',minval((vs-vs0)/vs0),maxval((vs-vs0)/vs0)
    print *


    ! F dm = S(m+dm) - S(m)
    ! note: original statement
    !       we backpropogate syn-dat (see adj_seismogram.f90) => we have to add a minus sign in front of kernels
    !       seems not valid for pressure adjoint sources
    if (SIM_TYPE == 1) then
      ! acoustic
      MS_l = sum(weights*(krhop)*(rho-rho0)/rho0) + sum(weights*(kalpha)*(vp-vp0)/vp0)
      MODELSPACE = MODELSPACE + MS_l
    else
      ! elastic
      MS_l = sum(weights*(-krhop)*(rho-rho0)/rho0) + sum(weights*(-kalpha)*(vp-vp0)/vp0) + sum(weights*(-kbeta)*(vs-vs0)/vs0)
      MODELSPACE = MODELSPACE + MS_l
    endif
    print *,'  model space contribution = ',MS_l
    print *,'         total model space = ',MODELSPACE
    print *


    !!! calculate inner product in data space --- < F dm, F dm>
    print *,'  calculating inner product in data space: < F dm, F dm>'

    write(procname,"(i4)") myrank
    procname = adjustl(procname)

    if (SIM_TYPE == 1) then
      ! acoustic
      comp_total = 1  ! only p-comp adjoint trace
    else
      ! elastic
      comp_total = 3   ! all x/y/z-comp adjoint traces
    endif

    do icomp = 1,comp_total

      if (SIM_TYPE == 1) then
        ! acoustic
        name = trim(procname)//"_p_SU"  ! pressure adjoint source
      else
        ! elastic
        ! file
        select case(icomp)
        case (1)
          name = trim(procname)//"_dx_SU"
        case (2)
          name = trim(procname)//"_dy_SU"
        case (3)
          name = trim(procname)//"_dz_SU"
        case default
          stop 'Invalid icomp'
        end select
      endif

      filename = "SEM/"//trim(name)//".adj"
      open(111,file=trim(filename),status='old',access='direct',recl=240+4*NSTEP,iostat=ier)
      if (ier /= 0) then
        print *,'Error opening file: ',trim(filename)
        stop 'Error opening adjoint trace'
      endif
      print *,'  ',trim(filename)

      DS_l = 0.d0
      irec = 1
      do while(ier == 0)
        adj(:) = 0.0

        read(111,rec=irec,iostat=ier) r4head,adj
        if (ier /= 0) exit


        DS_l = DS_l + sum(adj(:)*adj(:)) * DT

        irec = irec+1
      enddo
      close(111)

      DATASPACE = DATASPACE + DS_l
      print *,'  data space contribution = ',DS_l
      print *,'         total data space = ',DATASPACE
      print *

      ! count total (only once)
      if (icomp == 1) irec_total = irec_total + irec

    enddo ! icomp

    deallocate(rho,rho0,vp,vp0,vs,vs0)
    deallocate(weights,krhop,kalpha,kbeta)

  enddo ! NPROC

  print *
  print *,'total number of receivers considered: ',irec_total
  print *
  print *,'DATASPACE      = ',DATASPACE
  print *,'MODELSPACE     = ',MODELSPACE
  print *,'relative error = ',(DATASPACE-MODELSPACE)/DATASPACE
  print *
  print *,'total relative error should be small enough, < 0.001 should be OK'
  print *

end program random_model
