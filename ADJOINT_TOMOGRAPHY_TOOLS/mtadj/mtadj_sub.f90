module mtadj_sub

  use mtadj_constants
  use mtadj_variables
  use mtadj_sub2
  use mtadj_sub3

  implicit none

contains

! =========================================================

  subroutine read_mtadj_par(par_file)

    character(len=*) par_file
    integer :: ios

    open(10,file=trim(par_file),status='old',iostat=ios)
    if (ios /= 0) stop 'Error opening parameter file'

    ! kernel type, WF=0 (waveform), CC=1 (cross-corr), FD=2 (freq-dep)
    read(10,*) iker
    ! taper type for FD meas: BC=0 (boxcar), CS=1 (cosine-taper), MT=2 (multi-taper)
    read(10,*) itap
    read(10,'(a)') meas_dir
    read(10,'(a)') adj_dir

    ! bandpass data and syn
    read(10,*) BANDPASS
    read(10,*) tshort, tlong

    ! cross-correlation time shift maximum
    read(10,*) before_shift

    ! multi-taper setting
    read(10,*) npi,wtr,wtr_mtm

    ! select windows
    read(10,*) SELECT_WINDOW
    ! FD meas window length, window shorter than this length
    ! will resort to CC measurements and kernels
    read(10,*) ncycle_in_window
    read(10,*) dt_fac, err_fac, dt_max_scale ! ???
    !  read(10,*) after_quality, after_tshift, dlna_sigma_min

    ! interpolation of adjoint source
    read(10,*) INCLUDE_ERROR
    read(10,*) MIN_DT_SIGMA,MIN_DlnA_SIGMA
    read(10,*) b_adj, dt_adj, npts_adj
    read(10,*) BANDPASS_ADJ

    close(10)

    ! list input parameters:
    write(*,*) '========= INPUTS FROM MEASUREMENT.PAR ============'
    if (iker == IKER_WF) then
       write(*,*) 'Adjoint source type: Waveform '
    else if (iker == IKER_CC) then
       write(*,*) 'Adjoint source type: Cross-correlation (T & A)'
    else if (iker == IKER_FD) then
       write(*,*) 'Adjoint source type: Frequency-dependent (T & A)'
    else
       stop 'Only allow iker = IKER_WF/CC/FD'
    endif

    if (iker == IKER_FD) then
       if (itap == ITAP_BC) then
          write(*,*) 'Freq dep measurements use taper: boxcar'
          ntaper = 1
       else if (itap == ITAP_CS) then
          write(*,*) 'Freq dep measurements use taper: cosine taper'
          ntaper = 1
       else if (itap == ITAP_MT) then
          write(*,*) 'Freq dep measurements use taper: multi-taper'
          ntaper = 2*npi
          write(*,*) '  only for windows with at least ', ncycle_in_window, ' number of cycles'
          write(*,*) '  number of multi-tapers = ', ntaper
       endif
       write(*,*) ' only output results above water level ', wtr
       if (ntaper > NMAX_TAPER) stop 'ntaper should be <= NMAX_TAPER'
    endif

    if (BANDPASS) then
       write(*,*) 'Run bandpass before making measurements between',tshort, tlong
       fstart = 1./tlong; fend = 1./ tshort
    else
       write(*,*) 'No filtering before making measurements'
    endif
    write(*,*) 'cross-correlation time shift only allowed up to ', before_shift, ' secs'

    if (SELECT_WINDOW) then
       write(*,*) 'Only output adjoint source for selected windows based on:'
       write(*,*) dt_fac,dt_max_scale,err_fac
    else
       write(*,*) 'Adjoint sources for all windows will be used'
    endif

    if (npts_adj > NDIM) stop 'Npts for the adjoint source exceeds the array limit'
    if (NPT > NDIM) stop 'Check NPT > NDIM '
    write(*,*) 'Interpolate output adjoint source onto: ', b_adj, dt_adj, npts_adj
    if (INCLUDE_ERROR) then
       write(*,*) 'Errors are included in the adjoint sources'
    else
       write(*,*) 'Errors are not included in the adjoint sources -- test only !!'
    endif
    if (iker == IKER_CC) then
       write(*,*) '  CC adjoint source minimum error: ', MIN_DT_SIGMA,MIN_DlnA_SIGMA
    endif
    write(*,*) ' '

    ! frequency-domain taper (power of cosine) for adjoint source calculation
    ipwr_w = 4
    ! time domain taper for adjoint source
    ipwr_t = 10

  end subroutine read_mtadj_par

  !========================================================

  subroutine read_data_syn(datafile,synfile,sta,chan,net)

    character(len=*),intent(in) :: datafile,synfile
    character(len=*),intent(out) :: sta,chan,net

    integer :: npts1, npts2, nerr, j
    real :: b1, b2, dt1, dt2

    ! read data and syn
    call rsac1(datafile,data,npts1,b1,dt1,NDIM,nerr)
    if (nerr > 0) stop ' Error reading synthetic file'
    call rsac1(synfile,syn,npts2,b2,dt2,NDIM,nerr)
    if (nerr > 0) stop ' Error reading data file '
    if (npts1 /= npts2 .or. abs(b1-b2) > dt1 .or. abs(dt1-dt2) > EPS_dt) &
         stop 'check if npts1=npts2, b1=b2 and dt1=dt2'

    npts=npts1; b=b1; dt=dt1

    ! sac header sta.net.chan, and file prefixes
    call getkhv('kstnm', sta, nerr)
    call getkhv('kcmpnm', chan, nerr)
    call getkhv('knetwk', net, nerr)

    ! add comparsion for date and time to make sure they correspond
    ! to each other

  end subroutine read_data_syn

!============================================================

  subroutine cc_fd_measure(file_prefix,tstart,tend)

    character(len=*) file_prefix
    real :: tstart, tend
    real, dimension(NPT) :: dataw_save, synw_rc_cc,  synw_rc_fd
    real*8, dimension(NPT) :: ey1, ey2
    integer :: ishift, i, nerr
    real :: f0, df, df_new, ampmax_syn, wtr_amp_syn, ampmax_bot, wtr_amp_bot
    integer :: nf, idf_new, iampmax_syn, i_right_stop, iampmax_bot, ictaper
    integer :: i_ampmax_syn(1)
    complex*16, dimension(NPT) :: cdataw, csynw, cdatawt,csynwt
    complex,dimension(NPT) :: top_fdm, bot_fdm, trans_fdm,csynw_sngl
    character(len=150) :: meas_prefix
    real*8,dimension(NPT,NMAX_TAPER) :: tas_dp
    real*8 :: dt_dble


    ! set measurement file prefix
     if (iker /= IKER_FD) then
        meas_prefix=trim(CKER(iker+1))
     else
        meas_prefix=trim(CTAP(itap+1))
     endif
     meas_prefix=trim(meas_dir)//'/'//trim(file_prefix)//'.'//trim(meas_prefix)

    ! at least 10 points within the window
    if (tend-tstart < 10*dt) stop 'Check if tend > tstart+10*dt'
    nstart = floor((tstart-b)/dt)+1
    nend = floor((tend-b)/dt)
    nlen = nend-nstart+1
    if (nlen > NPT) stop 'Check if the window length is over NPT'

    ! window data and syn (data -> dataw; syn -> synw)
    dataw(1:nlen) = data(nstart:nend)
    synw(1:nlen) = syn(nstart:nend)

    ! for BC or CS tapers, rmean and rtrend
    if (iker == IKER_FD .and. itap == ITAP_BC) then
       call rmean(dataw,nlen); call rmean(synw,nlen)
       call rtrend(dataw,nlen); call rtrend(synw,nlen)
    endif

    if (DEBUG) then
       ! write windowed data and synthetics as sac files
       call wsac1(trim(meas_prefix)//'.obs.sac',dataw,nlen,tstart,dt,nerr)
       if (nerr > 0) stop 'Error writing obs data'
       call wsac1(trim(meas_prefix)//'.syn.sac',synw,nlen,tstart,dt,nerr)
       if (nerr > 0) stop 'Error writing synthetics'
    endif

    ! do not need t-domain tapers for windowed data and syn at this point

    ! ============== CC measurements =======================
    ! save a copy of observed data
    dataw_save(1:nlen) = dataw(1:nlen)

    ! tshift_cc: synw(t) * dataw(t+tshift_cc)
    call compute_time_shift(synw,dataw,nlen,dt,ishift,tshift_cc)
    if (abs(tshift_cc) > BEFORE_SHIFT) &
         stop 'Check if BEFORE_SHIFT is too small for tshift'

    ! align data and synthetics according to CC
    do i = 1, nlen
      if ( (i+nstart-1+ishift) >= 1 .and. (i+nstart-1+ishift) <= npts) then
         dataw(i) = data(i+ishift-1+nstart)
      else
         dataw(i) = 0. ! may be should be replaced by 1st and last value of data(:)
      endif
    enddo
    if (DEBUG) then
       call wsac1(trim(meas_prefix)//'.obs_shifted.sac',dataw,nlen,tstart,dt,nerr)
       if (nerr > 0) stop 'Error writing obs shifted data'
    endif

    ! CC measurements - amplitude dlnA
    dlnA = 0.5 * log( sum(dataw(1:nlen)**2) / sum(synw(1:nlen)**2) )

    ! reconstruct best-fitting syn according to tshift_cc and dlnA
    call reconstruct_syn_cc(synw,nlen,dt,ishift,dlnA,synw_rc_cc)
    if (DEBUG) then
       call wsac1(trim(meas_prefix)//'.syn_rc_cc.sac',synw_rc_cc,nlen,tstart,dt,nerr)
       if (nerr > 0) stop 'Error writing CC reconstructed syn'
    endif

    ! -------
    ! frequency vector from fft
    f0=0; df=1./(NPT*dt); nf=floor(NPT/2.)+1

    ! true independent frequency spacing for nlen
    df_new=1./(tend-tstart); idf_new = int(df_new/df)

    ! FFT windowed data (shifted) and syn
    cdataw = cmplx(0.,0.); csynw = cmplx(0.,0.)
    cdataw(1:nlen)=cmplx(dataw(1:nlen)); csynw(1:nlen)=cmplx(synw(1:nlen))

    call fft(LNPT,cdataw,FORWARD_FFT,dble(dt))
    call fft(LNPT,csynw,FORWARD_FFT,dble(dt))

    ! check the highest trustable frequency according to synthetics water level
    ampmax_syn = maxval(abs(csynw(1:nf)))
    i_ampmax_syn = maxloc(abs(csynw(1:nf)))
    wtr_amp_syn = cmplx(ampmax_syn * wtr, 0.)
    i_pmax=i_ampmax_syn(1)

    ! estimate tshift, dlnA uncertainties
    ! according to misfit between shifted data and reconstructed syn
    call compute_cc_error(dataw_save,synw_rc_cc,nlen,dt,i_pmax,dlnA,&
         sigma_tshift_cc,sigma_dlnA_cc,MIN_DT_SIGMA,MIN_DlnA_SIGMA)

    ! write measurement file
    open(30,file=trim(meas_prefix)//'.dt_dlnA',status='unknown')
    write(30,*) tshift_cc, sigma_tshift_cc
    write(30,*) dlnA, sigma_dlnA_cc
    close(30)

    ! DONE if just cross-correlation measurements
    if (iker /= IKER_FD) then
       dataw(1:nlen) = dataw_save(1:nlen)
       return
    endif

    ! ============= FD measurements ====================

    ! corresponding index to wtr_use_unw -> i_right
    i_right = nf;  i_right_stop = 0
    do i = 1,nf
       if (abs(csynw(i))<=abs(wtr_amp_syn) .and. i_right_stop==0 .and. i>i_pmax ) then
          i_right_stop = 1; i_right = i
       endif
       if (abs(csynw(i))>=10.*abs(wtr_amp_syn) .and. i_right_stop==1 .and. i>i_pmax) then
          i_right_stop = 0; i_right = i
       endif
    enddo
    f_right = df*i_right

    if (itap == ITAP_MT) then
       df_fd = npi*df_new; idf_fd = floor(df_fd/df); ; i_left=idf_fd
    else ! ignore first few points for cosine and boxcar tapers
       df_fd = df_new; idf_fd = floor(df_fd/df); i_left= floor(npi*df_fd/df)
    endif
    f_left=df_fd

    if (DEBUG) then
       print *, 'Frequency of max power in windowed synthetic:'
       print *, 'f_pmax = ', i_pmax * df, 'Hz, T_pmax = ', 1./(i_pmax*df), 'secs'
       print *, 'Frequency spacing df_new = ', df_new, ' Hz'
       print *, '  Longest T = ', 1/df_new, 's;  shortest T = ', 1/(i_right * df),'secs'
       print *, '  i_right = ', i_right
       print *, 'Independent MTM spacing df_meas = ', df_new * npi, ' Hz'
    endif

    ! define tapers for FD measurements: tas(1:nlen,1:ntaper)
    if (itap == ITAP_MT) then
      call staper(nlen, dble(NPI), ntaper, tas_dp, NPT, ey1, ey2)
    else if (itap == ITAP_CS) then
      call costaper(nlen, NPT, tas_dp)
    else if (itap == ITAP_BC) then
      call boxcar(nlen, NPT, tas_dp)
    endif
    tas=tas_dp

    ! compute transfer function for freq-dep measurements
    top_fdm(:)   = cmplx(0.,0.)
    bot_fdm(:)   = cmplx(0.,0.)
    trans_fdm(:) = cmplx(0.,0.)

    do ictaper = 1, ntaper

       ! tapered data (win+shifted) and syn (win)
       cdatawt(:)=cmplx(0.,0.); csynwt(:) = cmplx(0.,0.)
       cdatawt(1:nlen)=cmplx(dataw(1:nlen)*tas(1:nlen,ictaper))
       csynwt(1:nlen)=cmplx(synw(1:nlen)*tas(1:nlen,ictaper))

       ! complex sepctra
       call fft(LNPT,cdatawt,FORWARD_FFT,dble(dt))    ! syn
       call fft(LNPT,csynwt,FORWARD_FFT,dble(dt))    ! data

       ! top and bottom of transfer function
       top_fdm(1:nf) = top_fdm(1:nf) + cdatawt(1:nf) * conjg(csynwt(1:nf))
       bot_fdm(1:nf) = bot_fdm(1:nf) + csynwt(1:nf) * conjg(csynwt(1:nf))

    enddo

    ! water level of bottom
    ampmax_bot = maxval(abs(bot_fdm(1:nf)))
    if (itap == ITAP_MT) then
       wtr_amp_bot = ampmax_bot * (wtr_mtm ** 2)
    else
       wtr_amp_bot = ampmax_bot * (wtr ** 2)
    endif
    do i = 1, nf
      if(abs(bot_fdm(i)) > abs(wtr_amp_bot)) then
         trans_fdm(i) = top_fdm(i) /  bot_fdm(i)
      else if(abs(bot_fdm(i)) < abs(wtr_amp_bot)) then
         trans_fdm(i) = top_fdm(i) / (bot_fdm(i)+wtr_amp_bot)
      endif
    enddo

    ! compute dtau_fdm and dlnA_fdm
    call compute_dtau_dlnA(trans_fdm,dt,tshift_cc,dtau_fdm,dlnA_fdm,i_right)

    ! reconstruct syn with transfer function
    csynw_sngl=csynw
    call reconstruct_syn_fd(csynw_sngl,dtau_fdm,dlnA_fdm,i_right,synw_rc_fd,dt,nlen)

    if (DEBUG) then
       call wsac1(trim(meas_prefix)//'.syn_rc_fd.sac',synw_rc_fd,nlen,tstart,dt,nerr)
       if (nerr > 0) stop 'Error writing FD reconstructed syn'
    endif

    ! estimate transfer function error
    if (itap == ITAP_MT) then
       ! for multi-taper
       call compute_mt_error(ntaper,dataw,synw,tas,&
            nlen,dt,wtr_mtm,i_right,tshift_cc, &
            dtau_fdm,dlnA_fdm,sigma_dtau_fdm, sigma_dlnA_fdm)
    else
       ! for single tapers (BC or CS), use generic error estimation technique
       call compute_fd_error(npi,nlen,i_right,dt,dtau_fdm,dlnA_fdm,&
            sigma_dtau_fdm,sigma_dlnA_fdm)
    endif

    ! write transfer function measurements
    open(30,file=trim(meas_prefix)//'.dtau',status='unknown')
    open(40,file=trim(meas_prefix)//'.dlnA',status='unknown')

    do i = i_left, i_right, idf_new ! not all measurements are indep.
       write(30,'(3f14.4,i6)') df*(i-1), dtau_fdm(i), sigma_dtau_fdm(i), (i-i_left)/idf_fd
       write(40,'(3f14.4,i6)') df*(i-1), dlnA_fdm(i), sigma_dlnA_fdm(i), (i-i_left)/idf_fd
    enddo
    close(30)
    close(40)

    ! adjust errors according to min_dt_sigma, min_dlnA_sigma
    do i = 1, i_right
       sigma_dtau_fdm(i) = max(sigma_dtau_fdm(i),min_dt_sigma)
       sigma_dlnA_fdm(i) = max(sigma_dlnA_fdm(i),min_dlnA_sigma)
    enddo

    dataw(1:nlen) = dataw_save(1:nlen)

  end subroutine cc_fd_measure

! ======================================================================

  subroutine select_cc_fd_measure(tstart, tend, use_window)

    real :: tstart, tend
    logical :: use_window
    real :: df, f_pmax, T_pmax, wlen
    integer :: nf, j
    real, dimension(NPT) :: fvec


    use_window=.true.
    df = 1./(dt*NPT); nf=floor(NPT/2.)+1
    do j = 1, nf
       fvec(j) = df*(j-1)
    enddo

    ! at least N freq points between maximum power and f=0
    f_pmax = df * i_pmax
    T_pmax = 1./ f_pmax
    wlen = dt*nlen
    if( ncycle_in_window * T_pmax > wlen ) then
       use_window = .false.
       print *, 'rejecting window [', tstart, tend, ']'
       if (DEBUG) print *, ' wlen > ncycle * t_pmax= ', wlen, T_pmax
       return
    endif

    ! at least 2 independent measurements between f_left and f_right
    if ( (f_right-f_left) < df_fd) then
       use_window = .false.
       print *, 'rejecting window [',tstart, tend, ']'
       print *, '  --> considering CC adjoint source for this window'
       if (DEBUG) print *, 'f_left,f_right = ', f_left, f_right
       return
    endif

    ! require all FD measurements to be sensible
    if (iker == IKER_FD) then
       do j = i_left, i_right, idf_fd
          if ( abs(dtau_fdm(j)) > 1./(dt_fac*fvec(j)) .or. &
               abs(dtau_fdm(j)) > dt_max_scale*abs(tshift_cc) .or. &
               sigma_dtau_fdm(j) > 1./(err_fac*fvec(j)) ) then
             if (DEBUG) then
                write(*,*) 'Rejecting window ...'
                write(*,*) j, dtau_fdm(j), 1./(dt_fac*fvec(j)), &
                     dt_max_scale*abs(tshift_cc)
                write(*,*) '    ', sigma_dtau_fdm(j), 1./(err_fac*fvec(j))
             endif
             use_window = .false.; return
          endif
       enddo
    endif

  end subroutine select_cc_fd_measure


! =============================================================================


  subroutine mt_adj_src(file_prefix,iwin,tstart,dt_adj_src,amp_adj_src,dt_chi,amp_chi)

    character(len=*) :: file_prefix
    integer, intent(in) :: iwin
    real,intent(in) :: tstart
    real, dimension(NPT) :: dt_adj_src, amp_adj_src
    real :: dt_chi,amp_chi

    real, dimension(NPT) :: synw_veloc, synwt, synwt_veloc,&
         ft_bar, fa_bar, wf_taper, wp_taper, wq_taper, &
         d_bot_mtm, v_bot_mtm
    real*8,dimension(NPT) :: dtau_pj_t, dlnA_qj_t
    real :: ffac, dtau_wtr, dlnA_wtr, err_t, err_a, wtr_v_bot, wtr_d_bot,df
    complex*16, dimension(NPT) :: csynwt, csynwt_veloc, dtau_pj, dlnA_qj
    complex,dimension(NPT,NMAX_TAPER) :: pwc_adj, qwc_adj
    character(len=150) :: file
    integer :: ictaper, i, nerr

     if (iker /= IKER_FD) then
        file=trim(CKER(iker+1))
     else
        file=trim(CTAP(itap+1))
     endif
     file=trim(adj_dir)//'/'//trim(file_prefix)//'.'//trim(file)// '.'//char(iwin+48)

    ! IKER_WF
    if (iker == IKER_WF) then
       dt_adj_src(1:nlen) = synw(1:nlen)-dataw(1:nlen)
       dt_chi = 0.5 * sum(synw(1:nlen)-dataw(1:nlen) ** 2) * dt

    ! IKER_CC
    else if (iker == IKER_CC) then
       call compute_veloc_from_displ(synw,nlen,dt,synw_veloc)
       ! we could have taken away - sign from the measurements, and put them here
       ! just to be consistent with the FD measurements
       ft_bar = 0.; fa_bar = 0.
       ft_bar(1:nlen) = - synw_veloc(1:nlen) / ( sum(synw_veloc(1:nlen)**2) * dt )
       fa_bar(1:nlen) = synw(1:nlen) / ( sum(synw(1:nlen)**2) * dt )
       if (INCLUDE_ERROR) then
          err_t = sigma_tshift_cc; err_a = sigma_dlnA_cc
       else
          err_t = 1; err_a = 1
       endif
       ! include CC measurements (and error)
       dt_adj_src(1:nlen) = - (tshift_cc / err_t**2) * ft_bar(1:nlen)
       amp_adj_src(1:nlen) = - (dlnA / err_a**2 ) * fa_bar(1:nlen)

       dt_chi = 0.5 * (tshift_cc/err_t) ** 2
       amp_chi = 0.5 * (dlnA / err_a) ** 2

    ! IKER_FD
    else if (iker == IKER_FD) then
       df = 1./(dt*NPT)
       ! define frequency-domain taper W(f)
       wf_taper(:) = 0.
       i_left = 1 !! a better choice??
       do i = i_left, i_right
          !wf_taper(i) = 1.                                       ! boxcar
          !wf_taper(i) = 1. - (2.0/nw)**2 * ((i-1) - nw/2.0)**2     ! welch
          wf_taper(i) = 1. - cos(PI*(i-i_left)/(i_right-i_left))**ipwr_w    ! cosine
       enddo
       ! normalize freq taper
       ffac=2*df*sum(wf_taper(i_left:i_right))
       wf_taper(i_left:i_right) = wf_taper(i_left:i_right) / ffac

       ! water level of FD measurements based on average (LQY:can be smaller for mtm!!)
       dtau_wtr = wtr * sum(abs(dtau_fdm(i_left:i_right)))/(i_right-i_left)
       dlnA_wtr = wtr * sum(abs(dlnA_fdm(i_left:i_right)))/(i_right-i_left)

       ! include errors in the W(f) taper
       wp_taper(:) = 0.; wq_taper(:) = 0.
       do i = i_left, i_right
          wp_taper(i) = wf_taper(i); wq_taper(i) = wf_taper(i)
          if (INCLUDE_ERROR) then
             err_t = max(sigma_dtau_fdm(i),dtau_wtr)
             err_a = max(sigma_dlnA_fdm(i),dlnA_wtr)
             wp_taper(i) = wp_taper(i) / (err_t ** 2)
             wq_taper(i) = wq_taper(i) / (err_A ** 2)
          endif
       enddo

       d_bot_mtm = 0. ! should take only 1:i_right values, all positive numbers
       v_bot_mtm = 0.

       pwc_adj = cmplx(0.,0.); qwc_adj = cmplx(0., 0.)

       ! compute bottom (f) of p_j(f) and q_j(f)
       do ictaper = 1, ntaper

          synwt(1:nlen) = synw(1:nlen) * tas(1:nlen,ictaper)
          call compute_veloc_from_displ(synwt,nlen,dt,synwt_veloc)
          ! FFT
          csynwt(:) = 0.; csynwt_veloc = 0.
          csynwt(1:nlen) = cmplx(synwt(1:nlen))
          csynwt_veloc(1:nlen) = cmplx(synwt_veloc(1:nlen))

          call fft(LNPT,csynwt,FORWARD_FFT,dble(dt))
          call fft(LNPT,csynwt_veloc,FORWARD_FFT,dble(dt))

          d_bot_mtm = d_bot_mtm + real(csynwt * conjg(csynwt))
          v_bot_mtm = v_bot_mtm + real(csynwt_veloc * conjg(csynwt_veloc))

          pwc_adj(:,ictaper) =  csynwt_veloc
          qwc_adj(:,ictaper) = - csynwt

       enddo

       ! add water level to bottom of p_j(f) and q_j(f)
       wtr_v_bot = wtr * sum(abs(v_bot_mtm(1:i_right)))/i_right
       wtr_d_bot = wtr * sum(abs(d_bot_mtm(1:i_right)))/i_right
       do i = 1, i_right
          if (v_bot_mtm(i) < wtr_v_bot) v_bot_mtm(i) = v_bot_mtm(i)+wtr_v_bot
          if (d_bot_mtm(i) < wtr_d_bot) d_bot_mtm(i) = d_bot_mtm(i)+wtr_d_bot
       enddo

       ! initialize adjoint source
       dt_adj_src = 0. ; amp_adj_src = 0.
       do ictaper = 1, ntaper
          dtau_pj =  pwc_adj(:,ictaper)
          dlnA_qj =  qwc_adj(:,ictaper)
          dtau_pj(1:i_right) = dtau_pj(1:i_right)/ v_bot_mtm(1:i_right) &
               * cmplx(dtau_fdm(1:i_right), 0.) * cmplx(wp_taper(1:i_right),0.)
          dlnA_qj(1:i_right) = dlnA_qj(1:i_right) / d_bot_mtm(1:i_right) &
               * cmplx(dlnA_fdm(1:i_right), 0.) * cmplx(wp_taper(1:i_right),0.)

          ! IFFT into the time domain
          call fftinv(LNPT,dtau_pj,REVERSE_FFT,dble(dt),dtau_pj_t)
          call fftinv(LNPT,dlnA_qj,REVERSE_FFT,dble(dt),dlnA_qj_t)

          ! add up adjoint source
          dt_adj_src(:) = dt_adj_src(:) + tas(:,ictaper) * dtau_pj_t(:)
          amp_adj_src(:) = amp_adj_src(:) + tas(:,ictaper) * dlnA_qj_t(:)

       enddo

       dt_chi = 0.5 * 2.0 * df * sum( (dtau_fdm(1:i_right))**2 * wp_taper(1:i_right) )
       amp_chi = 0.5 * 2.0 * df * sum( (dlnA_fdm(1:i_right))**2 * wq_taper(1:i_right) )

    endif

    if (DEBUG) then
       print *, 'Writing adjoint file: ', trim(file)//'.adj.sac'
       call wsac1(trim(file)//'.adj.sac',dt_adj_src,nlen,tstart,dt,nerr)
       if (nerr > 0) stop 'Error writing windowed adjoint file'
    endif
    open(40,file=trim(file)//'.chi',status='unknown')
    write(40,*) dt_chi
    write(40,*) amp_chi
    close(40)

  end subroutine mt_adj_src


! ================================================

  subroutine adjust_adj_src(dt_adj_src,amp_adj_src,nlen,tstart,dt,&
                dt_adj_src_win,amp_adj_src_win,npts_adj,b_adj,dt_adj)

    real,dimension(:) :: dt_adj_src_win, amp_adj_src_win
    real, dimension(:) :: dt_adj_src, amp_adj_src
    real :: tstart, dt, b_adj, dt_adj
    integer :: nlen, npts_adj

    real,dimension(NPT) :: time_window
    real :: fac
    integer :: i

    ! time-domain taper
    time_window(1:nlen) = 0.
    do i = 1,nlen
      fac = 1.                                           ! boxcar window
      !fac = 1 - sfac2*((i-1) - dble(nlen1)/2.0)**2       ! welch window
      !fac = 1. - cos(PI*(i-1)/(nlen-1))**ipwr_t          ! cosine window
      time_window(i) = fac
    enddo
    dt_adj_src = dt_adj_src * time_window
    dt_adj_src = dt_adj_src * time_window

    ! interpolate adjoint source
    call interp_adj_src(dt_adj_src,nlen,tstart,dt, &
         dt_adj_src_win,npts_adj,b_adj,dt_adj)

    call interp_adj_src(amp_adj_src,nlen,tstart,dt, &
         amp_adj_src_win,npts_adj,b_adj,dt_adj)


  end subroutine adjust_adj_src


! ================================================
end module mtadj_sub
