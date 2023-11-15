program adj_seismogram

  implicit none

  integer :: myrank,ier
  integer :: irec,i,icomp,comp_total

  integer :: NSTEP,NPROC,SIM_TYPE
  double precision :: DT

  real(kind=4),dimension(:),allocatable :: syn,dat,adj
  real(kind=4) :: r4head(60)
  !integer(kind=4) :: header4(1)
  !integer(kind=2) :: header2(2)
  !equivalence(header2,header4)
  real(kind=4) :: diff_max,diff_max_all
  character(len=512) :: filename_in,filename_out,procname,arg,str_type

  ! derivative
  real(kind=4),dimension(:),allocatable :: adj_new
  double precision :: fac,val

  ! input parameters
  if (iargc() /= 4) &
    stop 'Usage: ./xadj NSTEP DT NPROC type[acoustic/elastic/acoustic-elastic]'

  call get_command_argument(1, arg); read(arg,*,iostat=ier)    NSTEP;   if (ier /= 0) stop 'Error reading NSTEP'
  call get_command_argument(2, arg); read(arg,*,iostat=ier)       DT;   if (ier /= 0) stop 'Error reading DT'
  call get_command_argument(3, arg); read(arg,*,iostat=ier)    NPROC;   if (ier /= 0) stop 'Error reading NPROC'
  call get_command_argument(4, arg); read(arg,*,iostat=ier) str_type;   if (ier /= 0) stop 'Error reading str_type'

  if (NSTEP <= 0) stop 'Invalid NSTEP, must be > 0'
  if (NPROC <= 0) stop 'Invalid NPROC, must be > 0'

  str_type = adjustl(str_type)
  if (trim(str_type) == "acoustic") then
    SIM_TYPE = 1
  else if (trim(str_type) == "elastic") then
    SIM_TYPE = 2
  else if (trim(str_type) == "acoustic-elastic") then
    SIM_TYPE = 3
  else
    print *,'Error: invalid type '//trim(str_type)//' must be: acoustic, elastic or acoustic-elastic'
    stop 'Invalid type, must be: acoustic, elastic or acoustic-elastic'
  endif

  ! user output
  print *,'creating adjoint seismograms:'
  print *,'  NSTEP = ',NSTEP
  print *,'  DT    = ',DT
  print *,'  NPROC = ',NPROC
  print *,'  type  = ',SIM_TYPE,'(acoustic == 1 / elastic == 2 / acoustic-elastic == 3)'
  print *

  allocate(syn(NSTEP),dat(NSTEP),adj(NSTEP),stat=ier)
  if (ier /= 0) stop 'Error allocating syn,dat,adj arrays'
  syn(:) = 0.0; dat(:) = 0.0; adj(:) = 0.0

  ! for derivative in acoustic case
  if (SIM_TYPE == 1 .or. SIM_TYPE == 3) then
    print *,"  creating adjoint sources for pressure (taking second-derivative of pressure differences)..."
    print *
    allocate(adj_new(NSTEP),stat=ier)
    if (ier /= 0) stop 'Error allocating adj_new arrays'
    adj_new(:) = 0.0
  endif

  do myrank = 0,NPROC - 1

    write(procname,"(i4)") myrank
    procname = adjustl(procname)

    !!!! read NSTEP from seismograms
    !!!filename_in=trim(procname)//"_dx_SU"
    !!!open(111,file="../OUTPUT_FILES/"//trim(filename_in),access='direct',recl=240,iostat = ier)
    !!!read(111,rec=1,iostat=ier) r4head
    !!!close(111)
    !!!header4=r4head(29)
    !!!NSTEP=header2(2)
    !!!header4=r4head(30)
    !!!DT=header2(1)*1.0d-6
    !!!print *, 'irec=',r4head(1)
    !!!print *, 'xs=',r4head(19)
    !!!print *, 'zs=',r4head(20)
    !!!print *, 'xr=',r4head(21)
    !!!print *, 'zr=',r4head(22)
    !!!print *, 'NSTEP=',NSTEP
    !!!print *, "DT=",DT

    ! read 'syn', 'dat' and then write 'adj'

    ! components
    select case(SIM_TYPE)
    case (1,3)
      ! acoustic / acoustic-elastic (w/ receivers in acoustic domain)
      comp_total = 1 ! single trace
    case (2)
      ! elastic
      comp_total = 3
    case default
      stop 'Invalid SIM_TYPE for component'
    end select

    ! user output
    print *,'rank',myrank,':'

    do icomp = 1,comp_total

      if (SIM_TYPE == 1 .or. SIM_TYPE == 3) then
        ! acoustic / acoustic-elastic
        ! pressure-component
        filename_in  = trim(procname)//"_p_SU"  ! from acceleration
        filename_out = trim(procname)//"_p_SU"  ! adjoint trace assumes _p_SU name for acoustic simulations
      else
        ! elastic
        select case(icomp)
        case(1)
          ! x-component
          filename_in  = trim(procname)//"_dx_SU" ! from displacement
          filename_out = trim(procname)//"_dx_SU" ! adjoint traces assumes d* name
        case (2)
          ! y-component
          filename_in  = trim(procname)//"_dy_SU"
          filename_out = trim(procname)//"_dy_SU" ! adjoint traces assumes d* name
        case (3)
          ! z-component
          filename_in  = trim(procname)//"_dz_SU"
          filename_out = trim(procname)//"_dz_SU" ! adjoint traces assumes d* name
        case default
          stop 'Invalid component'
        end select
      endif

      open(111,file="SEM/syn/"//trim(filename_in), &
                status='old',access='direct',action='read',recl=240+4*NSTEP,iostat=ier)
      if (ier /= 0) then
        !print *,'Error opening file syn: '//"SEM/syn/"//trim(filename_in)
        !stop 'Error opening file syn'
        ! skip if rank has no receivers
        cycle
      endif

      print *,'  input syn : ',"SEM/syn/" // trim(filename_in)
      print *,'  input data: ',"SEM/dat/" // trim(filename_in)

      open(112,file="SEM/dat/"//trim(filename_in), &
                status='old',access='direct',action='read',recl=240+4*NSTEP,iostat=ier)
      if (ier /= 0) then
        print *,'Error opening file dat: '//"SEM/dat/"//trim(filename_in)
        stop 'Error opening file dat'
      endif
      open(113,file="SEM/"//trim(filename_out)//".adj", &
                status='unknown',access='direct',action='write',recl=240+4*NSTEP,iostat=ier)
      if (ier /= 0) then
        print *,'Error opening file .adj: '//"SEM/"//trim(filename_out)//".adj"
        stop 'Error opening file .adj'
      endif

      ! statistics
      diff_max_all = -1.e20

      irec = 1
      do while(ier == 0)
        syn(:) = 0.0
        dat(:) = 0.0
        adj(:) = 0.0

        ! syn
        read(111,rec=irec,iostat=ier) r4head,syn
        if (ier /= 0) exit

        ! dat
        read(112,rec=irec,iostat=ier) r4head,dat
        if (ier /= 0) exit

        ! adjoint source: (s - d)
        adj(:) = syn(:) - dat(:)

        ! acoustic adjoint source
        if (SIM_TYPE == 1 .or. SIM_TYPE == 3) then
          ! for acoustic FWI, L2 adjoint source is the second derivative of pressure difference
          ! (e.g., see Peter et al. 2011, GJI, eq. (A8))
          !
          ! assuming a pressure output for syn and dat, the adjoint source expression is given by (A8) in Peter et al. (2011)
          ! note the negative sign in the definition.
          ! the adjoint source for pressure is: f^adj = - \partial_t^2 ( p_syn - p_obs )
          !
          ! takes second-derivative using central-difference scheme
          adj_new(:) = 0.0

          ! takes second-order derivative

          ! central finite difference (2nd-order scheme)
          !fac = 1.0 / DT**2
          !do i = 2,NSTEP-1
          !  ! central finite difference 2nd-order
          !  val = ( adj(i+1) - 2.d0 * adj(i) + adj(i-1) ) * fac
          !  ! adding negative sign
          !  adj_new(i) = - val
          !enddo

          ! central finite difference (4th-order scheme)
          ! (leads to slightly smoother derivative)
          fac = 1.0 / DT**2
          do i = 3,NSTEP-2
            ! central finite difference 4th-order
            val = ( - 1.d0/12.d0 * adj(i+2) + 4.d0/3.d0 * adj(i+1) - 5.d0/2.d0 * adj(i) &
                    + 4.d0/3.d0 * adj(i-1) - 1.d0/12.d0 * adj(i-2) ) * fac
            ! adding negative sign
            adj_new(i) = - val
          enddo

          ! sets as new adjoint source trace
          adj(:) = adj_new(:)
        endif

        write(113,rec=irec,iostat=ier) r4head,adj
        if (ier /= 0) exit

        ! statistics
        diff_max = maxval(abs(adj))
        if (diff_max > diff_max_all ) diff_max_all = diff_max

        !debug
        !daniel: outputs ascii trace
        if ( myrank == 0 .and. irec == 196 ) then
          open(221,file="SEM/syn/"//trim(filename_in)//".ascii",status='unknown')
          do i=1,NSTEP
            write(221,*) i,syn(i)
          enddo
          close(221)
          open(222,file="SEM/dat/"//trim(filename_in)//".ascii",status='unknown')
          do i=1,NSTEP
            write(222,*) i,dat(i)
          enddo
          close(222)
          open(223,file="SEM/"//trim(filename_out)//".adj.ascii",status='unknown')
          do i=1,NSTEP
            write(223,*) i,adj(i)
          enddo
          close(223)
        endif

        irec = irec+1
      enddo

      ! user output
      print *,'  receivers ',irec-1
      print *,'  adjoint sources written to: ',"SEM/"//trim(filename_out)//".adj"
      if (SIM_TYPE == 1 .or. SIM_TYPE == 3) then
        print *,'  maximum amplitude |f^adj| = |-\partial_t^2(p_syn-p_obs)| = ',diff_max_all
      else
        print *,'  maximum waveform difference (syn - dat) = ',diff_max_all
      endif
      print *

      close(111)
      close(112)
      close(113)

    enddo ! icomp

  enddo ! NPROC

  deallocate(syn,dat)

  print *,'done'

end program adj_seismogram
