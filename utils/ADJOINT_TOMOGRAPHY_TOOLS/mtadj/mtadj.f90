! this program computes the cross-correlation/transfer function measurements between
! data and synthetics in a WINDOWING file
! and output the corresponding adjoint source for tomography inversions

program mtadj

  use mtadj_constants
  use mtadj_variables
  use mtadj_sub
  !  use mtadj_sub2

  implicit none

  ! sac header information
  character(len=10) :: net,sta,chan

  ! file prefix
  character(len=150) :: datafile, synfile
  character(len=150) :: file_prefix, adj_prefix, meas_prefix

  ! npairs of data and syn, number of windows
  integer :: npairs, nwin, nwin_total, nadj_src
  integer ::  ios, ipair, j, nerr

  ! data and syn measurements/adjoint outpu
  real :: tstart, tend
  logical :: use_adj_trace, use_window
  real*8 :: fstart_dble,fend_dble,dt_dble,dt_adj_dble

  real,dimension(NPT) :: dt_adj_src,amp_adj_src
  real, dimension(NDIM) :: dt_adj_src_win, amp_adj_src_win, &
       dt_adj_src_all, amp_adj_src_all
  real :: dt_chi,amp_chi


  ! --- PROGRAM STARTS HERE ---
  ! read parameter file
  call read_mtadj_par('MTADJ.PAR')
  fstart_dble=fstart
  fend_dble=fend_dble
  dt_adj_dble=dt_adj

  ! loop over measurement file
  open(11,file='MEASUREMENT.WINDOWS',status='old',iostat=ios)
  if (ios /= 0) stop 'Error opening input file: MEASUREMENT WINDOWS'
  read(11,*,iostat=ios) npairs
  if (ios /= 0) stop 'Error reading number of pairs of data/syn'


  nwin_total=0
  do ipair = 1, npairs

     ! read data and syn pair
     read(11,'(a)',iostat=ios) datafile
     if (ios /= 0) stop 'Error reading datafile'
     read(11,'(a)',iostat=ios) synfile
     if (ios /= 0) stop 'Error reading synfile'
     if (DEBUG) print *, trim(datafile), ' ', trim(synfile)

     call read_data_syn(datafile,synfile,sta,net,chan)

     ! output file prefix
     file_prefix=trim(sta)//'.'//trim(net)//'.'//trim(chan)

     dt_dble=dt
     ! filter data and synthetics (check xapiir() usage in sac lib)
     if (BANDPASS) then
        call xapiir(data,npts,'BU',TRBDNDW,APARM,IORD,'BP',fstart_dble,fend_dble,dt_dble,PASSES)
        call xapiir(syn,npts,'BU',TRBDNDW,APARM,IORD,'BP',fstart_dble,fend_dble,dt_dble,PASSES)
     endif

     ! loop over nwin
     read(11,*) nwin
     if (nwin < 0) stop 'Check nwin '
     nwin_total = nwin_total + nwin

     ! initialize
     use_adj_trace = .false.
     dt_adj_src_all = 0.; amp_adj_src_all = 0.; nadj_src = 0

     do j = 1, nwin
        read(11,*,iostat=ios) tstart, tend
        if (ios /= 0) stop 'Error reading tstart and tend'
        print *, ' Measurement No. ', j, '...', tstart, tend
        tstart = max(tstart,b)
        tend = min(tend, b+(npts-1)*dt)

        if (iker == IKER_CC .or. iker == IKER_FD) then
           call cc_fd_measure(file_prefix, tstart, tend)
           if (SELECT_WINDOW) then
              call select_cc_fd_measure(tstart, tend, use_window)
           else
              use_window = .true.
           endif
        endif

        if (OUTPUT_ADJSRC .and. use_window) then
           call mt_adj_src(file_prefix,j,tstart,dt_adj_src,amp_adj_src,dt_chi,amp_chi)

           ! taper and interpolate into adj_src_all(t0,dt,npts)
           call adjust_adj_src(dt_adj_src,amp_adj_src,nlen,tstart,dt, &
                dt_adj_src_win,amp_adj_src_win,npts_adj,b_adj,dt_adj)
        endif

        if (use_window) use_adj_trace = .true.

     enddo ! nwin

     if (OUTPUT_ADJSRC .and. use_adj_trace) then
        dt_adj_src_all = dt_adj_src_all + dt_adj_src_win
        amp_adj_src_all = amp_adj_src_all + amp_adj_src_win
        if (BANDPASS_ADJ) then
           call xapiir(dt_adj_src_all,npts,'BU',TRBDNDW,APARM,IORD,'BP', &
                fstart_dble,fend_dble,dt_adj_dble,PASSES)
        endif
        if (iker /= IKER_FD) then
           adj_prefix=trim(CKER(iker+1))//'.adj'
        else
           adj_prefix=trim(CTAP(itap+1))//'.adj'
        endif
        adj_prefix=trim(adj_dir)//'/'//trim(file_prefix)//'.'//trim(adj_prefix)

        call wsac1(trim(adj_prefix),dt_adj_src_all,npts_adj,b_adj,dt_adj,nerr)
        ! write amplitude adjoint sources here
        nadj_src = nadj_src + 1
        if (nerr > 0) stop 'Error writing adjoint source'
     endif

  enddo ! nfiles

  print *, 'Total number of windows processed ', nwin_total
  print *, 'Total number of adjoint sources ', nadj_src

end program mtadj
