program random_model

  implicit none

  integer,parameter :: NGLLX = 5,NGLLY = 5,NGLLZ = 5
  integer,parameter :: IOUT = 20
  integer,parameter :: CUSTOM_REAL = 4

  ! data space misfit: compares either seismograms or adjoint sources
  logical,parameter :: use_data_traces = .true.     ! .true. == seismogram misfits, .false. == adjoint sources

  integer :: NPROC,NSTEP,SIM_TYPE
  integer :: myrank,ier,nspec,nglob
  integer :: irec,irec_total,comp_total,icomp

  double precision :: DT
  double precision :: DATASPACE,MODELSPACE
  double precision :: DS_l,MS_l,MS_rhol,MS_vpl,MS_vsl,error

  character(len=512) :: prname,arg,procname,str_type,name
  character(len=512) :: filename,filename_dat,filename_syn

  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: vp,vs,rho,weights
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: vp0,vs0,rho0
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: krhop,kalpha,kbeta,ktmp
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: dlnrho,dlnvp,dlnvs

  real(kind=4) :: r4head(60)

  real(kind=CUSTOM_REAL),dimension(:), allocatable :: dat,syn,adj,diff

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
  else if (trim(str_type) == "acoustic-elastic") then
    SIM_TYPE = 3
  else
    print *,'Error: Invalid type '//trim(str_type)//' must be: acoustic, elastic or acoustic-elastic'
    stop 'Invalid type, must be: acoustic, elastic or acoustic-elastic'
  endif

  print *,'postprocessing:'
  print *,'  NSTEP = ',NSTEP
  print *,'  DT    = ',DT
  print *,'  NPROC = ',NPROC
  print *,'  type  = ',SIM_TYPE,'(acoustic == 1 / elastic == 2 / acoustic-elastic == 3)'
  print *

  MODELSPACE = 0.d0
  DATASPACE  = 0.d0

  irec_total = 0

  if (use_data_traces) then
    allocate(dat(NSTEP),syn(NSTEP),diff(NSTEP),stat=ier)
    if (ier /= 0) stop 'Error allocating array dat,syn,diff'
    dat(:) = 0.0; syn(:) = 0.0; diff(:) = 0.0
  else
    allocate(adj(NSTEP),stat=ier)
    if (ier /= 0) stop 'Error allocating array adj'
    adj(:) = 0.0
  endif

  do myrank = 0,NPROC-1
    print *
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

    ! temporary kernel to combine both acoustic and elastic ones for coupled acoustic-elastic case
    if (SIM_TYPE == 3) then
      allocate(ktmp(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'error allocating array ktmp'
      ktmp(:,:,:,:) = 0.0_CUSTOM_REAL
    endif

    ! kernels
    ! K_rho_prime
    allocate(krhop(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'error allocating array krhop'
    krhop(:,:,:,:) = 0.0_CUSTOM_REAL

    if (SIM_TYPE == 1) then
      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'rhop_acoustic_kernel.bin'
    else
      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'rhop_kernel.bin'
    endif
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file rhop_kernel.bin'
    read(IOUT) krhop
    close(IOUT)

    ! coupled acoustic-elastic
    if (SIM_TYPE == 3) then
      ! adds contribution from acoustic kernel (elastic done above)
      ktmp(:,:,:,:) = 0.0_CUSTOM_REAL

      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'rhop_acoustic_kernel.bin'
      open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'Error opening file rhop_acoustic_kernel.bin'
      read(IOUT) ktmp
      close(IOUT)
      ! combined elastic & acoustic
      krhop = krhop + ktmp
    endif

    ! K_alpha
    allocate(kalpha(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'error allocating array kalpha'
    kalpha(:,:,:,:) = 0.0_CUSTOM_REAL

    if (SIM_TYPE == 1) then
      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'alpha_acoustic_kernel.bin'
    else
      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'alpha_kernel.bin'
    endif
    open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening file alpha_kernel.bin'
    read(IOUT) kalpha
    close(IOUT)

    ! coupled acoustic-elastic
    if (SIM_TYPE == 3) then
      ! adds contribution from acoustic kernel (elastic done above)
      ktmp(:,:,:,:) = 0.0_CUSTOM_REAL

      filename = 'OUTPUT_FILES/DATABASES_MPI/'//trim(prname)//'alpha_acoustic_kernel.bin'
      open(unit=IOUT,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'Error opening file rhop_acoustic_kernel.bin'
      read(IOUT) ktmp
      close(IOUT)
      ! combined elastic & acoustic
      kalpha = kalpha + ktmp
    endif

    ! K_beta, zero shear for acoustic simulation
    allocate(kbeta(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'error allocating array kbeta'
    kbeta(:,:,:,:) = 0.0_CUSTOM_REAL

    if (SIM_TYPE == 1) then
      ! acoustic, no shear kernel
      kbeta = 0.0
    else
      ! elastic or acoustic-elastic (w/ elastic elements in slice)
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
    allocate(rho(NGLLX,NGLLY,NGLLZ,nspec), &
             rho0(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if ( ier /= 0 ) stop 'Error allocating array rho,rho0'
    rho(:,:,:,:) = 0.0_CUSTOM_REAL; rho0(:,:,:,:) = 0.0_CUSTOM_REAL

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
    allocate(vp(NGLLX,NGLLY,NGLLZ,nspec), &
             vp0(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if ( ier /= 0 ) stop 'Error allocating array vp,vp0'
    vp(:,:,:,:) = 0.0_CUSTOM_REAL; vp0(:,:,:,:) = 0.0_CUSTOM_REAL

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
    allocate(vs(NGLLX,NGLLY,NGLLZ,nspec), &
             vs0(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if ( ier /= 0 ) stop 'Error allocating array vs,vs0'
    vs(:,:,:,:) = 0.0_CUSTOM_REAL; vs0(:,:,:,:) = 0.0_CUSTOM_REAL

    if (SIM_TYPE == 1) then
      ! acoustic, zero shear velocity
      ! nothing to read in
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

    ! relative perturbations
    allocate(dlnrho(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array dlnrho'
    dlnrho(:,:,:,:) = 0.0_CUSTOM_REAL

    ! relative perturbation
    where(rho0(:,:,:,:) /= 0.0_CUSTOM_REAL) dlnrho = (rho - rho0) / rho0

    allocate(dlnvp(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array dlnvp'
    dlnvp(:,:,:,:) = 0.0_CUSTOM_REAL

    ! relative perturbation
    where(vp0(:,:,:,:) /= 0.0_CUSTOM_REAL) dlnvp = (vp - vp0) / vp0

    allocate(dlnvs(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if ( ier /= 0 ) stop 'Error allocating array dlnvs'
    dlnvs(:,:,:,:) = 0.0_CUSTOM_REAL

    if (SIM_TYPE == 1) then
      ! acoustic
      ! no shear pertubations
    else
      ! relative perturbation
      where(vs0(:,:,:,:) /= 0.0_CUSTOM_REAL) dlnvs = (vs - vs0) / vs0
    endif

    print *,'  relative model perturbation (rho - rhop)/rho0 : min/max = ',minval(dlnrho),maxval(dlnrho)
    print *,'  relative model perturbation (vp - vp0)/vp0    : min/max = ',minval(dlnvp),maxval(dlnvp)
    if (SIM_TYPE /= 1) &
      print *,'  relative model perturbation (vs - vs0)/vs0    : min/max = ',minval(dlnvs),maxval(dlnvs)
    print *

    ! F dm = S(m+dm) - S(m)
    ! note: original statement
    !       we backpropogate syn-dat (see adj_seismogram.f90) => we have to add a minus sign in front of kernels
    if (SIM_TYPE == 1) then
      ! acoustic (no shear contribution)
      MS_rhol = sum(weights * (-krhop) * dlnrho)
      MS_vpl  = sum(weights * (-kalpha) * dlnvp)
      MS_l = MS_rhol + MS_vpl
      print *,'  model space contribution = ',MS_l,' with Mrho = ',MS_rhol,' Mvp = ',MS_vpl

    else
      ! elastic or acoustic-elastic (w/ elastic elements in slice)
      MS_rhol = sum(weights * (-krhop) * dlnrho)
      MS_vpl  = sum(weights * (-kalpha) * dlnvp)
      MS_vsl  = sum(weights * (-kbeta) * dlnvs)
      MS_l = MS_rhol + MS_vpl + MS_vsl
      print *,'  model space contribution = ',MS_l,' with Mrho = ',MS_rhol,' Mvp = ',MS_vpl,' Mvs = ',MS_vsl
    endif

    MODELSPACE = MODELSPACE + MS_l
    print *,'         total model space = ',MODELSPACE
    print *

    !!! calculate inner product in data space --- < F dm, F dm>
    print *,'  calculating inner product in data space: < F dm, F dm>'

    write(procname,"(i4)") myrank
    procname = adjustl(procname)

    if (SIM_TYPE == 1 .or. SIM_TYPE == 3) then
      ! acoustic or acoustic-elastic (w/ receivers in acoustic domain)
      comp_total = 1  ! only p-comp adjoint trace
    else
      ! elastic
      comp_total = 3   ! all x/y/z-comp adjoint traces
    endif

    do icomp = 1,comp_total

      if (SIM_TYPE == 1 .or. SIM_TYPE == 3) then
        ! acoustic or acoustic-elastic
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

      if (use_data_traces) then
        ! reads in (forward) data traces
        ! data d
        filename_dat = "OUTPUT_FILES.dat.forward/" // trim(name)
        open(111,file=trim(filename_dat),status='old',access='direct',recl=240+4*NSTEP,iostat=ier)
        if (ier /= 0) then
          !print *,'Error opening file: ',trim(filename)
          !stop 'Error opening adjoint trace'
          ! skip if rank has no receivers
          cycle
        endif
        print *,'  ',trim(filename_dat)

        ! synthetics s
        filename_syn = "OUTPUT_FILES.syn.forward/" // trim(name)
        open(112,file=trim(filename_syn),status='old',access='direct',recl=240+4*NSTEP,iostat=ier)
        if (ier /= 0) then
          print *,'Error opening file: ',trim(filename)
          stop 'Error opening adjoint trace'
        endif
        print *,'  ',trim(filename_syn)

        DS_l = 0.d0
        irec = 1
        do while(ier == 0)
          dat(:) = 0.0
          syn(:) = 0.0
          diff(:) = 0.0

          ! data d
          read(111,rec=irec,iostat=ier) r4head,dat
          if (ier /= 0) exit

          ! synthetics s
          read(112,rec=irec,iostat=ier) r4head,syn
          if (ier /= 0) exit

          ! misfit (s - d)
          diff(:) = syn(:) - dat(:)

          ! dataspace contribution
          DS_l = DS_l + sum( diff(:) * diff(:) ) * DT

          irec = irec+1
        enddo
        close(111)
        close(112)

      else
        ! reads in adjoint traces
        filename = "SEM/" // trim(name) // ".adj"
        open(111,file=trim(filename),status='old',access='direct',recl=240+4*NSTEP,iostat=ier)
        if (ier /= 0) then
          !print *,'Error opening file: ',trim(filename)
          !stop 'Error opening adjoint trace'
          ! skip if rank has no receivers
          cycle
        endif
        print *,'  ',trim(filename)

        DS_l = 0.d0
        irec = 1
        do while(ier == 0)
          adj(:) = 0.0

          read(111,rec=irec,iostat=ier) r4head,adj
          if (ier /= 0) exit

          ! dataspace contribution
          DS_l = DS_l + sum( adj(:) * adj(:) ) * DT

          irec = irec+1
        enddo
        close(111)
      endif

      DATASPACE = DATASPACE + DS_l
      print *,'  data space contribution = ',DS_l
      print *,'         total data space = ',DATASPACE
      print *

      ! count total (only once)
      if (icomp == 1) irec_total = irec_total + (irec-1)

    enddo ! icomp

    ! free arrays
    deallocate(rho,rho0,vp,vp0,vs,vs0)
    deallocate(dlnrho,dlnvp,dlnvs)
    deallocate(weights,krhop,kalpha,kbeta)
    if (SIM_TYPE == 3) deallocate(ktmp)

  enddo ! NPROC

  ! model vs. data space
  ! relative error
  if (DATASPACE /= 0.0) then
    error = (DATASPACE - MODELSPACE) / DATASPACE
  else
    error = abs(MODELSPACE)
  endif

  print *
  print *,'total number of receivers considered: ',irec_total
  print *
  print *,'DATASPACE      = ',DATASPACE
  print *,'MODELSPACE     = ',MODELSPACE
  print *
  print *,'relative error (dataspace - modelspace)/dataspace = ',error
  print *
  print *,'total relative error should be small enough, < 0.1 should be OK'
  print *

  if (abs(error) < 0.1) then
    print *,'looks ok'
    print *
  else
    print *,'error seems too large, please check...'
    ! exit with an error
    stop 1
  endif

end program random_model
