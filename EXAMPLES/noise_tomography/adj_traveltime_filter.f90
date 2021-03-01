!!!!! this subroutine is used for preparing adjoint sources in this example ONLY
!!!!! generally, you should make your own measuremnts and calculate corresponding adjoint sources
!!!!! FLEXWIN could be of help, or you can modify this code a little bit

program adj_traveltime

  implicit none

  !--------------------------------------------------------------
  ! USER PARAMETERS

  ! adjoint station
  character(len=10),parameter :: station='X2'
  character(len=10),parameter :: network='DB'
  character(len=3),parameter :: comp = 'MXZ'

  ! measurement window length (in s)
  double precision, parameter :: length_time_window = 70.0d0

  ! taper type (0 = boxcar / 1 = cosine)
  integer, parameter :: taper_type = 1

  ! filters traces (0 = off / 1 = on)
  integer,parameter :: filter_flag = 1

  ! cross-correlation branch for measurement (0 = negative / 1 = positive branch)
  integer,parameter :: branch_type = 1

  !---------------------------------------------------------------

  double precision, parameter :: pi = 3.141592653589793d0

  double precision,dimension(:),allocatable :: data_trace,data_origin,data_temp,adj
  double precision,dimension(:),allocatable :: t_trace
  double precision :: trace_data_max

  character(len=150) :: file_data
  character(len=150) :: file_misfit

  !double precision :: Norm_data_crit,AMP,Norm_data,Norm
  double precision :: Norm_adj_temp
  double precision :: misfit_traveltime, traveltime_delay
  integer :: l,ier

  ! taper
  double precision,dimension(:),allocatable :: taper

  ! window
  double precision,dimension(:),allocatable :: window
  double precision,dimension(:),allocatable :: adj_temp,data_trace_temp
  integer :: itime, length_window, i_start_window, i_end_window

  ! time range
  integer :: nstep
  double precision :: dt

  double precision :: dummy_t,dummy_val
  double precision :: t0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! filter frequencies
  integer :: ifreq, nfreq
  real :: F1,F2,D(8),G,DELT
  integer :: ISW

  ! frequency band
  !real, parameter :: freq_low(2) = (/ 0.0001d , 0.0001d /)
  !real, parameter :: freq_high(2) = (/ 0.5d , 0.05d /)
  ! single frequency band
  real, parameter :: freq_low(1) = 0.01d0    ! 1 / 100s
  real, parameter :: freq_high(1) = 0.2d0    ! 1 / 5s

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! user output
  print *,'noise adjoint source'
  print *
  print *,'station: ',trim(network)//'.'//trim(station)// '.' //comp
  print *,'  using filter flag : ',filter_flag,'(0 = no filter / 1 = bandpass)'
  print *,'  using taper type  : ',taper_type ,'(0 = boxcar / 1 = cosine)'
  print *,'  using measurement window length: ',sngl(length_time_window),'s'
  print *,'  using cross-correlation branch : ',branch_type,'(0 = negative / 1 = positive branch)'
  print *

  ! reads in number of steps (based on first trace)
  ! trace
  file_data = './SEM/'//trim(network)//'.'//trim(station)//'.'//comp//'.semd'
  open(unit=1001,file=trim(file_data),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(file_data)
    stop 'Error opening station trace file in SEM/ '
  endif
  nstep = 0
  dt = 0.d0
  t0 = 0.d0
  do while (ier == 0)
    read(1001,*,iostat=ier) dummy_t,dummy_val
    ! counts time steps
    if (ier == 0) then
      nstep = nstep + 1
    else
      exit
    endif
    ! gets size of time step
    if (nstep == 1) then
      t0 = dummy_t
    else
      dt = abs(dummy_t - t0) / dble(nstep-1)
    endif
  enddo
  close(1001)

  ! user output
  print *,'data file: ',trim(file_data)
  print *,'  number of time steps = ',nstep
  print *,'  time step size       = ',sngl(dt),'s'
  print *,'  trace length         = ',sngl(nstep * dt),'s'
  print *

  ! checks
  if (dt <= 0.d0) stop 'Error time step is zero, please check station trace'
  if (nstep <= 0) stop 'Error number of time steps is zero, please check station trace'

  ! time step (in milliseconds)
  DELT = dt*1.0d3

  ! allocates trace arrays
  allocate(data_temp(nstep), &
           data_trace(nstep), &
           data_origin(nstep), &
           t_trace(nstep), &
           adj(nstep), &
           window(nstep), &
           stat=ier)
  if (ier /= 0) stop 'Error allocating arrays'

  t_trace(:) = 0.d0
  data_trace(:) = 0.d0
  adj(:) = 0.d0

  !!!! loading data !!!!
  ! trace - vertical component
  file_data = './SEM/'//trim(network)//'.'//trim(station)//'.'//comp//'.semd'
  open(unit=1001,file=trim(file_data),status='old',action='read',iostat=ier)
  if (ier /= 0) stop 'Error opening station trace file in SEM/ '

  ! reads in data
  do itime = 1,nstep
    ! original
    !read(1001,*) t_trace(itime),data_trace(itime)
    ! reversed
    ! the reversed seismogram involves $C^\alpha\beta(t)=C^\beta\alpha(-t)$
    read(1001,*) t_trace(itime),data_trace(nstep-itime+1)
  enddo

  close(1001)

  ! stores original data
  data_origin(:) = data_trace(:)

  ! takes time derivative (by central differences)
  data_temp(1) = 0.d0
  data_temp(nstep) = 0.d0
  do itime = 2,nstep-1
    ! displacement --> particle velocity (elastic  ".semd")
    data_temp(itime) = ( data_trace(itime+1) - data_trace(itime-1) )/ (2.d0*dt)
  enddo
  data_trace(:) = data_temp(:)

  ! maximum value
  trace_data_max = maxval(abs(data_trace(:)))
  print *,'  vertical component: maximum value = ',sngl(trace_data_max)
  print *

  ! shifts time line (to start at zero)
  !t0 = t_trace(1)
  !do itime = 1,nstep
  !  t_trace(itime) = t_trace(itime) - t0
  !enddo

  !!!! taper !!!!
  length_window = floor(length_time_window/dt)+1
  if (length_window > (NSTEP+1)/2) length_window = (NSTEP+1)/2

  allocate(taper(length_window), &
           adj_temp(length_window), &
           data_trace_temp(length_window), &
           stat=ier)
  if (ier /= 0) stop 'Error allocating taper/window arrays'

  taper(:) = 1.0d0
  adj_temp(:) = 0.0d0
  data_trace_temp(:) = 0.0d0

  if (taper_type == 1) then
    ! cosine taper, otherwise using a constant (1.0) instead
    do l = 1,length_window
      taper(l) = (1.0-cos(pi*2.0*(l-1)/(length_window-1)))/2.0
    enddo
  endif

  ! misfit file output
  file_misfit = './SEM/misfit_traveltime_delay'
  open(unit=1111,file=trim(file_misfit),status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening misfit output file in SEM/'

  ! number of filter frequencies
  if (filter_flag == 0) then
    nfreq = 1
  else
    nfreq = size(freq_low)
  endif

  !!!! computing adj sources !!!!
  misfit_traveltime = 0.0d0

  ! loops over frequencies
  do ifreq = 1,nfreq
    !!!!!!!!!!!!!!!!!! to be implemented !!!!!!!!!!!!!!!!!
    if (filter_flag == 1) then
      data_temp(:) = data_trace(:)
      ! THIS SECTION CALCULATES THE FILTER AND MUST BE CALLED BEFORE
      ! FILTER IS CALLED
      F1 = freq_low(ifreq)
      F2 = freq_high(ifreq)

      print *,'filtering:'
      print *,'  frequency band: ',ifreq
      print *,'  f_min = ',sngl(F1),'Hz , f_max = ',sngl(F2), 'Hz'
      print *,'  T_min = ',sngl(1.d0/F2),'(s) , T_max = ',sngl(1.d0/F1), '(s)'

      call BNDPAS(F1,F2,DELT,D,G,nstep,ISW)
      !    F1 = LOW FREQUENCY CUTOFF (6 DB DOWN)
      !    F2 = HIGH FREQUENCY CUTOFF (6 DB DOWN)
      !    DELT = SAMPLE INTERVAL IN MILLISECONDS
      !    D = WILL CONTAIN 8 Z DOMAIN COEFFICIENTS OF RECURSIVE FILTER
      !    G = WILL CONTAIN THE GAIN OF THE FILTER,
      call FILTER(data_temp(:),nstep,D,G,2,ISW)
      !     X = DATA VECTOR OF LENGTH N CONTAINING DATA TO BE FILTERED
      !     D = FILTER COEFFICIENTS CALCULATED BY BNDPAS
      !     G = FILTER GAIN
      !     IG = 1  one pass
      !     ig = 2  two passes
      data_trace(:) = data_temp(:)
    endif

    !!!! cross-correlation !!!!
    window(:) = 0.d0

    ! window starting index
    i_start_window = 0
    if (branch_type == 0) then
      !the negative branch
      i_start_window = (nstep+1)/2 - length_window
    else if (branch_type == 1) then
      !the positive branch
      i_start_window = (nstep+1)/2
    else
      stop 'Error invalid cross-correlation branch type, please use either 0 (negative) or 1 (positive)'
    endif

    ! window end index
    i_end_window = i_start_window + length_window - 1
    window(i_start_window : i_end_window) = taper(:)

    ! norms
    !Norm_data_crit = sqrt(DOT_PRODUCT(data_trace(:),data_trace(:)))
    !AMP = Norm_data_crit

    !data_trace_temp(:) = data_trace(i_start_window:i_end_window)
    !Norm_data = sqrt(DOT_PRODUCT( data_trace_temp(:),data_trace_temp(:) ))
    !Norm = Norm_data*Norm_data

    ! tapered/windowed data
    data_trace(:) = data_trace(:) * window(:)

    !!!! normal adjoint sources !!!!
    adj_temp(:) = data_trace(i_start_window : i_end_window)

    ! minus sign comes from integration by part
    Norm_adj_temp = - DOT_PRODUCT(adj_temp(:),adj_temp(:)) * dt

    !!!! choose which one to use !!!!
    ! normalized adjoint sources
    if (abs(Norm_adj_temp) > 0.d0) then
      adj_temp(:) = adj_temp(:) / Norm_adj_temp
    else
      ! zero trace
      adj_temp(:) = 0.d0
    endif
    ! ray density map DeltaT=+1
    !adj_temp(:) = adj_temp(:) / Norm_adj_temp

    ! sets adjoint source
    adj(i_start_window : i_end_window) = adj_temp(:)

    ! traveltime delay
    traveltime_delay = dt
    write(1111,*) traveltime_delay

    ! total misfit
    misfit_traveltime = misfit_traveltime + traveltime_delay * traveltime_delay / 2.d0

    print *,'adjoint source norm     = ',sngl(Norm_adj_temp)
    print *,'traveltime delay        = ',sngl(traveltime_delay)
    print *,'misfit traveltime       = ',sngl(misfit_traveltime)

    ! data window output
    write(file_misfit,'(a,i1,a)') './SEM/misfit_measurement_window',ifreq,'.dat'
    open(unit=1001,file=trim(file_misfit),status='unknown',iostat=ier)
    if (ier /= 0) stop 'Error opening misfit output file in SEM/'
    write(1001,*) '#time   #data_origin    #data_trace    #adjoint_source    #window'
    do itime = 1,nstep
      write(1001,*) t_trace(itime),data_origin(itime),data_trace(itime),adj(itime),window(itime)
    enddo
    close(1001)

  enddo  !do ifreq=1,nfreq

  ! closes misfit file
  close(1111)

  !!!! output !!!!
  ! note: only vertical component has non-zero adjoint source for noise sources

  ! 1st contribution
  file_data = './SEM/adj_sources_contribution1'
  open(unit=1002,file=trim(file_data),status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening output adjoint file in SEM/ '
  do itime = 1,nstep
    write(1002,*) t_trace(itime), adj(nstep-itime+1)
  enddo
  close(1002)
  ! time-reversed
  file_data = './SEM/adj_sources_contribution2'
  open(unit=1002,file=trim(file_data),status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening output adjoint file in SEM/ '
  do itime = 1,nstep
    write(1002,*) t_trace(itime), adj(itime)
  enddo
  close(1002)

  ! traveltime misfit output
  file_misfit = './SEM/misfit_traveltime'
  open(unit=1001,file=trim(file_misfit),status='unknown')
  write(1001,*) misfit_traveltime, traveltime_delay
  close(1001)

  print *
  print *,'done, see adjoint trace contributons in directory: SEM/'
  print *,'  adj_sources_contribution1 -- adjoint source file containing 1st contribution'
  print *,'  adj_sources_contribution2 -- adjoint source file containing time-reversed, 2nd contribution'
  print *

end program adj_traveltime


!
!-----------------------------------------------------------------------------
!
! subroutine for bandpass filter
!
! slighlty modified version from reference in:
!
! E.R. Kanasewich, Time Sequence Analysis in Geophysics, 3rd edition, The University of Alberta Press, 1981
!
! see: (access July 2015)
! https://books.google.com.sa/books?id=k8SSLy-FYagC&pg=PA274&lpg=PA274&dq=bndpas.for&source=bl&ots=gsjFXkN1FZ&sig=W-qvA2kamMr5xIkEIzY_f2yciOI&hl=en&sa=X&redir_esc=y#v=onepage&q=bndpas.for&f=false
!
! or see: website (access ~2012)
! http://www-lgit.obs.ujf-grenoble.fr/users/jrevilla/seiscomp/patch/pack/plugins/seisan/LIB/bndpas.for


  subroutine BNDPAS(F1,F2,DELT,D,G,N,ISW)

! RECURSIVE BUTTERWORTH BAND PASS FILTER (KANASEWICH, TIME SERIES
! ANALYSIS IN GEOPHYSICS, UNIVERSITY OF ALBERTA PRESS, 1975; SHANKS,
! JOHN L, RECURSION FILTERS FOR DIGITAL PROCESSING, GEOPHYSICS, V32,
! FILTER.  THE FILTER WILL HAVE 8 POLES IN THE S PLANE AND IS
! APPLIED IN FORWARD AND REVERSE DIRECTIONS SO AS TO HAVE ZERO
! PHASE SHIFT.  THE GAIN AT THE TWO FREQUENCIES SPECIFIED AS
! CUTOFF FREQUENCIES WILL BE -6DB AND THE ROLLOFF WILL BE ABOUT
! THE FILTER TO PREVENT ALIASING PROBLEMS.

  implicit none

  real, intent(in) :: F1,F2,DELT
  integer, intent(in) :: N

  real, intent(out) :: D(8),G
  integer, intent(out) :: ISW

  ! local parameters
  COMPLEX :: P(4),S(8),Z1,Z2

  real :: DT,TDT,FDT,W1,W2,WW,HWID,A,B,C
  integer :: I

  double precision, parameter :: TWOPI = 2.d0 * 3.141592653589793d0

! THIS SECTION CALCULATES THE FILTER AND MUST BE CALLED BEFORE
! FILTER IS CALLED

!    F1 = LOW FREQUENCY CUTOFF (6 DB DOWN)
!    F2 = HIGH FREQUENCY CUTOFF (6 DB DOWN)
!    DELT = SAMPLE INTERVAL IN MILLISECONDS
!    D = WILL CONTAIN 8 Z DOMAIN COEFFICIENTS OF RECURSIVE FILTER
!    G = WILL CONTAIN THE GAIN OF THE FILTER,

  DT = DELT/1000.0
  TDT = 2.0/DT
  FDT = 4.0/DT

  ISW = 1

  P(1) = CMPLX(-.3826834,.9238795)
  P(2) = CMPLX(-.3826834,-.9238795)
  P(3) = CMPLX(-.9238795,.3826834)
  P(4) = CMPLX(-.9238795,-.3826834)

  W1 = TWOPI*F1
  W2 = TWOPI*F2
  W1 = TDT*TAN(W1/TDT)
  W2 = TDT*TAN(W2/TDT)
  HWID = (W2-W1)/2.0
  WW = W1*W2

  do I = 1,4
    Z1 = P(I)*HWID
    Z2 = Z1*Z1-WW
    Z2 = CSQRT(Z2)
    S(I) = Z1+Z2
    S(I+4) = Z1-Z2
  enddo

  G = 0.5/HWID
  G = G*G
  G = G*G

  do I = 1,7,2
    B = -2.0*REAL(S(I))
    Z1 = S(I)*S(I+1)
    C = REAL(Z1)
    A = TDT+B+C/TDT
    G = G*A
    D(I) = (C*DT-FDT)/A
    D(I+1) = (A-2.0*B)/A
  enddo

  G = G*G

  !debug
  !print *,'debug: FILTER GAIN IS ',G

  end subroutine BNDPAS

  subroutine FILTER(X,N,D,G,IG,ISW)

!     X = DATA VECTOR OF LENGTH N CONTAINING DATA TO BE FILTERED
!     D = FILTER COEFFICIENTS CALCULATED BY BNDPAS
!     G = FILTER GAIN
!     IG = 1  one pass
!     IG = 2  two passes

  implicit none
  integer,intent(in) :: N
  double precision,intent(inout) :: X(N)

  real,intent(in) :: D(N),G
  integer,intent(in) :: IG

  integer, intent(in) :: ISW

  ! local parameters
  real :: XC(3),XD(3),XE(3)
  real :: XM,XM1,XM2,GG
  integer :: K,I,J,M,M1,M2

  ! check
  if (ISW /= 1) then
    print *,'1BNDPAS MUST BE CALLED BEFORE FILTER'
    stop 'Invalid ISW in FILTER() routine'
  endif

  ! APPLY FILTER IN FORWARD DIRECTION
  XM2 = X(1)
  XM1 = X(2)
  XM = X(3)

  XC(1) = XM2
  XC(2) = XM1-D(1)*XC(1)
  XC(3) = XM-XM2-D(1)*XC(2)-D(2)*XC(1)

  XD(1) = XC(1)
  XD(2) = XC(2)-D(3)*XD(1)
  XD(3) = XC(3)-XC(1)-D(3)*XD(2)-D(4)*XD(1)

  XE(1) = XD(1)
  XE(2) = XD(2)-D(5)*XE(1)
  XE(3) = XD(3)-XD(1)-D(5)*XE(2)-D(6)*XE(1)

  X(1) = XE(1)
  X(2) = XE(2)-D(7)*X(1)
  X(3) = XE(3)-XE(1)-D(7)*X(2)-D(8)*X(1)

  do I = 4,N
    XM2 = XM1
    XM1 = XM
    XM = X(I)

    K = I-((I-1)/3)*3

    select case(K)
    case (1)
      M = 1
      M1 = 3
      M2 = 2
    case (2)
      M = 2
      M1 = 1
      M2 = 3
    case (3)
      M = 3
      M1 = 2
      M2 = 1
    case default
      stop 'Invalid K value in FILTER'
    end select

    XC(M) = XM-XM2-D(1)*XC(M1)-D(2)*XC(M2)
    XD(M) = XC(M)-XC(M2)-D(3)*XD(M1)-D(4)*XD(M2)
    XE(M) = XD(M)-XD(M2)-D(5)*XE(M1)-D(6)*XE(M2)

    X(I) = XE(M)-XE(M2)-D(7)*X(I-1)-D(8)*X(I-2)
  enddo

  if (IG /= 1) then
    XM2 = X(N)
    XM1 = X(N-1)
    XM = X(N-2)

    XC(1) = XM2
    XC(2) = XM1-D(1)*XC(1)
    XC(3) = XM-XM2-D(1)*XC(2)-D(2)*XC(1)

    XD(1) = XC(1)
    XD(2) = XC(2)-D(3)*XD(1)
    XD(3) = XC(3)-XC(1)-D(3)*XD(2)-D(4)*XD(1)

    XE(1) = XD(1)
    XE(2) = XD(2)-D(5)*XE(1)
    XE(3) = XD(3)-XD(1)-D(5)*XE(2)-D(6)*XE(1)

    X(N) = XE(1)
    X(N-1) = XE(2)-D(7)*X(1)
    X(N-2) = XE(3)-XE(1)-D(7)*X(2)-D(8)*X(1)

    DO I = 4,N
      XM2 = XM1
      XM1 = XM

      J = N-I+1
      XM = X(J)

      K = I-((I-1)/3)*3

      select case (K)
      case (1)
        M = 1
        M1 = 3
        M2 = 2
      case (2)
        M = 2
        M1 = 1
        M2 = 3
      case (3)
        M = 3
        M1 = 2
        M2 = 1
      case default
        stop 'Invalid K value in FILTER'
      end select

      XC(M) = XM-XM2-D(1)*XC(M1)-D(2)*XC(M2)
      XD(M) = XC(M)-XC(M2)-D(3)*XD(M1)-D(4)*XD(M2)
      XE(M) = XD(M)-XD(M2)-D(5)*XE(M1)-D(6)*XE(M2)
      X(J) = XE(M)-XE(M2)-D(7)*X(J+1)-D(8)*X(J+2)
    enddo
  endif

  if (IG == 1) then
    GG = sqrt(G)   ! if only pass once, modify gain
  else
    GG = G
  endif

  do I = 1,N
    X(I) = X(I)/GG
  enddo

  end subroutine FILTER


