module wave2d_sub3

  use wave2d_constants
  use wave2d_variables

  implicit none

! this module contains subroutines pertaining to filtering gridpoints
! and converting amoung UTM, mesh, and index coordinates

contains

  !------------------------------------------

  subroutine mtm_adj(ipick, ievent, nrec, syn, tstart, tend, adj_syn, data)

  ! PARAMETERS ARE INCLUDED IN wave2d_constants.f90
  !parameter(NDIM=8000*4,NDIM2=16000*4,lnpt=14,npt=2**lnpt)
  !parameter(MAXTAPER=5,wtr=0.02,ZZIGN=-1.0)
  !implicit real(a-h,o-z)

  double precision :: tas(NDIM,MAXTAPER), ey1(NDIM), ey2(NDIM), wpi
  !double precision, dimension(NDIM,MAXTAPER) :: phi_mul, dt_mul, abs_mul, abs2_mul
  double precision, dimension(NDIM) :: dzr, dzr2, phi, fr
  double precision, dimension(NDIM) :: dzr_win, dzr2_win, dzr3_win, dzr0_win, dzr30_win, tseis_rec
  !double precision, dimension(NDIM) :: dt_mtm, phi_mtm, abs_mtm, abs2_mtm
  !double precision, dimension(NDIM) :: err_phi, err_dt, err_abs, err_abs2

  ! variables for adjoint sources (p(w) and q(w), etc)
  complex*16, dimension(NDIM,MAXTAPER) :: top_p_ntaper
  complex*16, dimension(NDIM) :: pwc_adj, qwc_adj, top_p, bot_p, ctemp, dtemp, dtau_pj_wc, dlnA_qj_wc
  double precision, dimension(NDIM,MAXTAPER) :: pw_adj, qw_adj, pt_adj, qt_adj, dtau_pj_t, dlnA_qj_t
  double precision, dimension(NDIM) :: fvec, tvec, dlnA_w, dtau_w, fp, fq, w_taper
  double precision :: nw

  ! cross-correlation sources
  double precision, dimension(NDIM) :: datt, synt, syn_veloc, syn_accel, ft_bar_t, fa_bar_t
  double precision :: Nnorm, Mnorm, dlna

  complex*16, dimension(NDIM) :: wseis1, wseis, wseis3, wseis_rec
  complex*16, dimension(NDIM) :: trans, trans_mtm, top_mtm, bot_mtm
  complex*16 :: wtr_use, wtr_use_unw, cci
  character*20 hfmt

  ! -----------------
  integer, intent(in) :: nrec
  double precision, dimension(NSTEP,NCOMP,nrec),intent(in) :: syn
  double precision, dimension(nrec), intent(in) :: tstart, tend
  double precision, dimension(NSTEP,NCOMP,nrec),intent(out) :: adj_syn

  double precision, dimension(NSTEP,NCOMP,nrec),intent(in),optional :: data

  integer itime, icomp, istart, iend, nlen
  double precision, dimension(NSTEP) :: time_window
  double precision :: norm, junk, ntemp
  ! -----------------

  double precision :: twopi, wtr_mtm, fac, cr_shift, cc, cc_max, sfac
  double precision :: tshift_xc, ampmax, ampmax_unw, omega
  double precision :: t_stop, t_start, df, df_new, temp
  integer :: ipwr, ishift, n_left, n_right, ntaper, fnum, ipick, i1
  integer :: i_right, i_right_stop, i_amp_max, i_amp_max_unw, idf_new

  ! error associated with measurement
  double precision :: meas_pert, rand1, ppert

  integer :: i,j,ik,n,ictaper,irec,ievent

!=======================================================

  twopi = 2.*PI
  cci = cmplx(0.,1.)
  wtr_mtm = 1.e-10

  adj_syn(:,:,:) = 0.
  icomp = 1

  print *, ' (mtm_adj.f90) ipick = ', ipick

  ! loop over receivers -- tstart and tend are nrec x 1 arrays
  do irec = 1,nrec

    istart = max(floor(tstart(irec) / DT), 1)
    iend = min(floor(tend(irec) / DT), NSTEP)
    !print *, istart, iend
    if (istart >= iend) stop 'Check if istart < istop'
    nlen = iend - istart

    ! create the time-windowed records for which the measurement is made
    do itime=1,NSTEP
      dzr(itime) = syn(itime,icomp,irec)      ! synthetic
      dzr2(itime) = data(itime,icomp,irec)    ! data
    enddo

    t_start = tstart(irec)
    t_stop = tend(irec)

    ! initialise windowed arrays
    do i = 1, npt
      dzr_win(i) = 0
      dzr2_win(i) = 0
    enddo

    sfac = (2./dble(nlen))**2  ! for Welch window
    ipwr = 10                  ! for cosine window
    time_window(:) = 0.        ! post-processing taper

    do i = 1, nlen
      !fac = 1.                                       ! boxcar window
      !fac = 1 - sfac*((i-1) - dble(nlen)/2.)**2      ! welch window
      fac = 1. - cos(PI*(i-1)/(nlen+1))**ipwr        ! cosine window
      !print *, fac

      time_window(i) = fac

      fac = 1.
      dzr_win(i)  = dzr(i+istart) * fac
      dzr0_win(i) = dzr(i+istart) * fac
      dzr2_win(i) = dzr2(i+istart)* fac
    enddo

    !------------------------------------------------------------------
    ! Shift observed seismogram to maximise alignment prior to measure
    !------------------------------------------------------------------

    cr_shift = 300.0

    n_left = ceiling( (-1.0) * cr_shift / DT )
    n_right = floor( cr_shift / DT )
    ishift = 0
    cc_max = 0
    do i = n_left, n_right, 1
      cc = 0
      do j = 1, nlen
        if ((j+i) > 1 .and. (j+i) < nlen) cc = cc + dzr_win(j) * dzr2_win(j+i)
      enddo
      if ( cc > cc_max) then
        cc_max = cc
        ishift = i
        endif
    enddo
    tshift_xc = ishift*DT  ! KEY: cross-correlation time shift

    !===================================================
    ! if you want a MTM measurement, then go here
    if (IKER == 3 .or. IKER == 4) then

    ! apply time shift to observed seismogram
    write(*,*) 'shift obs seismogram by ', tshift_xc, 'seconds, irec = ', irec
    do i = 1, nlen
      dzr3_win(i) = 0
      if ( (ishift+i) > 1 .and. (ishift+i) < nlen ) dzr3_win(i) = dzr2_win(i+ishift)
      dzr30_win(i) = dzr3_win(i)
    enddo

    ! create complex synthetic seismogram
    do i = 1,npt
      if (i <= nlen) then
        wseis1(i) = cmplx(dzr0_win(i),0)
      else
        wseis1(i) = 0
      endif
    enddo

    ! apply f.t. to get spectrum
    call clogc(lnpt,wseis1,ZZIGN,DT)

    ! multi-taper measurement
    ntaper = 5
    wpi = 2.5
    call staper(nlen, wpi, ntaper, tas, NDIM, ey1, ey2)

    ! might be a problem here
    write(hfmt,'(a,i2.2,a)') '(', ntaper,'e18.6)'

    !------------------------------------------------------------------

    ! initialize transfer function terms
    top_mtm(:) = cmplx(0,0)
    bot_mtm(:) = cmplx(0,0)
    trans_mtm(:) = cmplx(0,0)

    do ictaper = 1, ntaper  ! loop over tapers

      ! apply taper ictaper to synth and obs windowed seismograms
      do i = 1, nlen
        dzr_win(i) = dzr0_win(i) * tas(i,ictaper)
        dzr3_win(i) = dzr30_win(i) * tas(i,ictaper)
      enddo

      ! create complex seismograms
      do i = 1,npt
        if (i <= nlen) then
          wseis(i) = cmplx(dzr_win(i),0)
          wseis3(i) = cmplx(dzr3_win(i),0)
        else
          wseis(i) = 0
          wseis3(i) = 0
        endif
      enddo

      ! apply f.t. to get complex spectra
      call clogc(lnpt,wseis,ZZIGN,DT)
      call clogc(lnpt,wseis3,ZZIGN,DT)

      ! calculate frequency step and number of frequencies
      df = 1/(npt*DT)
      fnum = npt/2 + 1

      ! find water level for single taper measurement
      ampmax = 0
      ampmax_unw = 0
      do i = 1, fnum
        if ( abs(wseis(i)) > ampmax) then
          ampmax =  abs(wseis(i))
          i_amp_max = i
        endif
        if ( abs(wseis1(i)) > ampmax_unw) then
          ampmax_unw =  abs(wseis1(i))
          i_amp_max_unw = i
        endif
      enddo
      wtr_use = cmplx(ampmax * wtr, 0)
      wtr_use_unw = cmplx(ampmax_unw * wtr, 0)

      ! these variables define maximum frequency for measurement
      ! i_right_stop = 1 --> stop at frequency i_right, not fnum
      i_right = fnum
      i_right_stop = 0

      ! loop over frequencies
      do i = 1, fnum

        ! calculate top and bottom of transfer function for multitapers
        ! NOTE THAT THESE QUANTITIES ARE SUMMED OVER THE TAPERS AS WELL
        top_mtm(i) = top_mtm(i) +  wseis3(i) * conjg(wseis(i))
        bot_mtm(i) = bot_mtm(i) +  wseis(i) * conjg(wseis(i))

        ! calculate transfer function for single taper measurement using water level
        if (abs(wseis(i)) > abs(wtr_use))  trans(i) = wseis3(i) / wseis(i)
        if (abs(wseis(i)) <= abs(wtr_use))  trans(i) = wseis3(i) / (wseis(i)+wtr_use)

        ! determine i_right values using the power in the un-tapered synthetic
        if (abs(wseis1(i)) <= abs(wtr_use_unw) .and. i_right_stop == 0 .and. i > i_amp_max_unw) then
          i_right_stop = 1
          i_right = i
        endif
        if (abs(wseis1(i)) >= 10*abs(wtr_use_unw) .and. i_right_stop == 1 .and. i > i_amp_max_unw) then
          i_right_stop = 0
          i_right = i
        endif

      enddo  ! loop over frequencies (i=1,fnum)

      print *, ' taper number ', ictaper, ' out of ', ntaper

    enddo  ! end loop over tapers (ictaper=1,ntaper)

    !------------------------------------------------------------------

    ! calculate frequency spacing of sampling points
    df_new = 1.0 / (t_stop-t_start)
    idf_new = df_new / df

    print *, ' measurement spacing df_new is ', sngl(df_new)
    print *, ' fnum is', fnum

    !------------------------------------------------------------------
    ! Multi-taper measurement, calculation and output
    !------------------------------------------------------------------

    ! find water level for multi taper measurement
    ampmax = 0
    do i = 1, fnum
      if ( abs(bot_mtm(i)) > ampmax) then
        ampmax =  abs(bot_mtm(i))
        i_amp_max = i
      endif
    enddo
    wtr_use = cmplx(ampmax * wtr_mtm**2, 0)

    ! calculate transfer function using water level
    do i = 1, fnum
      if (abs(bot_mtm(i)) > abs(wtr_use)) trans_mtm(i) = top_mtm(i) / bot_mtm(i)
      if (abs(bot_mtm(i)) <= abs(wtr_use)) trans_mtm(i) = top_mtm(i) / (bot_mtm(i)+wtr_use)
    enddo

    !=======================================================
    ! construct time series : tau(omega), dlnA(omega)

    ! taper function for the freq domain
    nw = dble(i_right - 1)

    ! loop to calculate phase and amplitude
    ! NOTE: here we include the time shift
    dtau_w(:) = 0.
    dlnA_w(:) = 0.
    w_taper(:) = 0.
    open(91,file='transfer_freq.dat',status='unknown')
    do i = 1, i_right

      fvec(i)     = df*i     ! do not divide by zero
      !dtau_w(i)   = -atan2(aimag(trans_mtm(i)), real(trans_mtm(i))) - ZZIGN*tshift_xc
      dtau_w(i)   = -(1./fvec(i))*atan2(aimag(trans_mtm(i)), real(trans_mtm(i))) - ZZIGN*tshift_xc
      dlnA_w(i)   = abs(trans_mtm(i))-1

      ! type of filter in the freq domain : boxcar, welch, cosine
      !w_taper(i) = 1.
      !w_taper(i) = 1. - (2./nw)**2 * ((i-1) - nw/2.)**2
      w_taper(i) = 1. - cos(PI*(i-1)/(nw+1))**ipwr

      write(91,'(4e18.8)') fvec(i), dtau_w(i), dlnA_w(i), w_taper(i)
    enddo
    close(91)

    ! Reconstruct mtm fit seismograms : syn*tran
    ! d(w) = s(w) T(w) exp[-i w dT]
    ! trans_mtm is for syn --> SHIFTED data
    wseis_rec(:) = cmplx(0,0)
    do i = 1,i_right
      omega = twopi * (i-1) * df
      !wseis_rec(i) = wseis1(i)*exp(-cci*omega*tshift_xc)*trans_mtm(i)  ! sign
      wseis_rec(i) = wseis1(i) * exp(-cci*omega*dtau_w(i)) * (1.+ dlnA_w(i))
      !wseis_rec(i) = wseis1(i) * exp(-cci*omega*dtau_w(i)*w_taper(i)) * (1.+ dlnA_w(i)*w_taper(i))
    enddo

    ! inverse FFT into time domain
    call ftinv(lnpt,wseis_rec,ZZIGN,DT,tseis_rec)

    open(17,file='recovered_seis.dat',status='unknown')
    do i = 1,nlen
      write(17,'(2e18.8)') DT*i, tseis_rec(i)
    enddo
    close(17)

    endif  ! IKER == 3,4

    !==================================================================
    !     compute cross-correlation adjoint sources
    !==================================================================

    synt(:) = dzr0_win(:)   ! synthetics
    datt(:) = dzr2_win(:)   ! data, unshifted

    ! calculate velocity and acceleration from syn
    do i = 2, nlen-1
      syn_veloc(i) =  (synt(i+1) - synt(i-1)) / (2 * DT)
    enddo
    syn_veloc(1) = (synt(2) - synt(1)) / DT
    syn_veloc(nlen) = (synt(nlen) - synt(nlen-1)) /DT

    do i = 2, nlen-1
      syn_accel(i) =  (syn_veloc(i+1) - syn_veloc(i-1)) / (2 * DT)
    enddo
    syn_accel(1) = (syn_veloc(2) - syn_veloc(1)) / DT
    syn_accel(nlen) = (syn_veloc(nlen) - syn_veloc(nlen-1)) /DT

    ! cross-correlation traveltime adjoint source (normalized velocity)
    ! NOTE sign convention (Qinya, not GJI paper)
    Nnorm = -DT * sum( synt(:) * syn_accel(:) )
    ft_bar_t(:) = -syn_veloc(:) / Nnorm

    ! cross-correlation amplitude adjoint source (normalized displacement)
    Mnorm = DT * sum( synt(:) * synt(:) )
    fa_bar_t(:) = synt(:) / Mnorm

    ! definition of Dahlen and Baig (2002), Eq. 3,17,18 : dlnA = Aobs/Asyn - 1
    dlna = sqrt( (DT * sum( datt(:) * datt(:) )) / (DT * sum( synt(:) * synt(:) )) ) - 1.

    if (0 == 1) then
      print *
      print *, 'cross-correlation measurments:'
      print *, '   dT = ', tshift_xc
      print *, '   dA = ', dlna
    else
      write(*,'(a8,i5,a8,1f18.10,a8,1f18.10)') &
         'irec = ', irec, ', dT = ', tshift_xc, ', dA = ', dlna
    endif

    ! additional files for checking (measure_socal_adj.m)
    if (0 == 1) then
      ! time domain : time, data-disp, syn-disp, syn-vel, syn-accel
      open(29,file='syn_time.dat')
      do i = 1,nlen
        write(29,'(7e18.8)') i*DT, datt(i), synt(i), syn_veloc(i), syn_accel(i)
      enddo
      close(29)

      ! xcorr adjoint sources : traveltime (sampling), traveltime (misfit), amplitude (sampling), amplitude (misfit)
      open(39,file='xcorr_time.dat')
      do i = 1,nlen
        write(39,'(4e18.8)') ft_bar_t(i), -tshift_xc*ft_bar_t(i), fa_bar_t(i), -dlna*fa_bar_t(i)
      enddo
      close(39)
    endif

    !==================================================================
    !     create MTM adjoint sources
    !==================================================================

    if (IKER == 3 .or. IKER == 4) then

       pw_adj(:,:) = 0.      ;  qw_adj(:,:) = 0.
       pt_adj(:,:) = 0.      ;  qt_adj(:,:) = 0.
       dtau_pj_t(:,:) = 0.   ;  dlnA_qj_t(:,:) = 0.
       fp(:) = 0.            ;  fq(:) = 0.
       bot_p(:) = cmplx(0,0)
       top_p_ntaper(:,:) = cmplx(0,0)

       ! This loops to get the denominator of pj(w).
       ! It also stores the displacement field sj(w), which is in the numerator of p(w).
       do ictaper = 1,ntaper

          ! apply TAPER ictaper to windowed synthetics
          do i = 1,nlen
             dzr_win(i)  = dzr0_win(i) * tas(i,ictaper)
          enddo

          ! create complex seismograms
          wseis(:) = cmplx(dzr_win(:),0)

          ! apply f.t. to get complex spectra
          call clogc(lnpt,wseis,ZZIGN,DT)

          ! bottom of p function for multitapers
          do i = 1, i_right
             bot_p(i) = bot_p(i) +  wseis(i) * conjg(wseis(i))
          enddo

          ! term in numerator (sj)
          top_p_ntaper(:,ictaper) = wseis(:)

       enddo  ! loop over tapers

       do ictaper = 1,ntaper   ! loop over tapers
          print *, ' taper number ', ictaper, ' out of ', ntaper

          top_p(:) = top_p_ntaper(:,ictaper)  ! top of p function for multitapers
          pwc_adj(:) = cmplx(0,0)
          qwc_adj(:) = cmplx(0,0)

          ! compute pj(w) and qj(w) -- j term is in top_p
          do i=1,i_right
             omega = twopi*df*i   ! omega not =0

             ! NOTE SIGN
             pwc_adj(i) = cmplx(0,1./omega) * top_p(i) / bot_p(i)  ! (1/w)i = (-iw)/(-w^2)
             qwc_adj(i) = cmplx(0,omega) * pwc_adj(i)              ! d/dt < -- > iw

          enddo

          ! NOTE THAT ftinv CHANGES THE INPUT, SO WE SAVE A COPY HERE.
          ctemp(:) = pwc_adj(:)
          dtemp(:) = qwc_adj(:)

          ! EXTRA : ifft into the time domain : pj(w) --> pj(t) and qj(w) --> qj(t)
          call ftinv(lnpt,pwc_adj,ZZIGN,DT,pt_adj(:,ictaper))
          call ftinv(lnpt,qwc_adj,ZZIGN,DT,qt_adj(:,ictaper))

          ! create dtau(w) pj(w) W(w) and save dtau(t) * pj(t) * W(t)
          dtau_pj_wc(:) = ctemp(:) * cmplx(dtau_w(:),0.) * cmplx(w_taper(:),0.)
          call ftinv(lnpt,dtau_pj_wc,ZZIGN,DT,dtau_pj_t(:,ictaper))

          ! create dlnA(w) qj(w) W(w) and save dlnA(t) * qj(t) * W(t)
          dlnA_qj_wc(:) = dtemp(:) * cmplx(dlnA_w(:),0.) * cmplx(w_taper(:),0.)
          call ftinv(lnpt,dlnA_qj_wc,ZZIGN,DT,dlnA_qj_t(:,ictaper))

          ! create adjoint source
          fp(:) = fp(:) + tas(:,ictaper) * dtau_pj_t(:,ictaper)
          fq(:) = fq(:) + tas(:,ictaper) * dlnA_qj_t(:,ictaper)

       enddo

       ! write sampling adjoint sources to file
       open(21,file='test_fadj_t.dat',status='unknown')
       do i = 1,nlen
          write(21,*) sngl(fp(i)), sngl(fq(i))
       enddo
       close(21)

       ! time domain : tapers and other time series
       open(18,file='test_hj_t.dat',status='unknown')
       open(19,file='test_pj_t.dat',status='unknown')
       open(20,file='test_Pj_t.dat',status='unknown')
       open(21,file='test_Pj_hj_t.dat',status='unknown')
       open(22,file='test_qj_t.dat',status='unknown')
       open(23,file='test_Qj_t.dat',status='unknown')
       open(24,file='test_Qj_hj_t.dat',status='unknown')

       do i = 1,nlen  ! loop over time points
          write(18,hfmt) ( sngl(tas(i,ictaper)), ictaper=1,ntaper )                         ! hj(t)
          write(19,hfmt) ( sngl(pt_adj(i,ictaper)), ictaper=1,ntaper )                      ! pj(t)
          write(20,hfmt) ( sngl(dtau_pj_t(i,ictaper)), ictaper=1,ntaper )                   ! Pj(t)
          write(21,hfmt) ( sngl(dtau_pj_t(i,ictaper) * tas(i,ictaper)), ictaper=1,ntaper )  ! hj(t) Pj(t)

          write(22,hfmt) ( sngl(qt_adj(i,ictaper)), ictaper=1,ntaper )                      ! qj(t)
          write(23,hfmt) ( sngl(dlnA_qj_t(i,ictaper)), ictaper=1,ntaper )                   ! Qj(t)
          write(24,hfmt) ( sngl(dlnA_qj_t(i,ictaper) * tas(i,ictaper)), ictaper=1,ntaper )  ! hj(t) Qj(t)
       enddo

       close(18) ; close(19) ; close(20) ; close(21) ; close(22) ; close(23) ; close(24)

    endif  ! if IKER == 3,4

    ! generate a random number to simulate error in the measurement
    ! ppert determines the range over which the perturbed measurement will be
    ppert = 0.
    call random_number(rand1)
    meas_pert = 1. + ppert*(2.*rand1 - 1.)
    !meas_pert = 1. + ppert

    !print *, tshift_xc, meas_pert, tshift_xc * meas_pert

    ! CREATE THE ADJOINT SOURCES
    ! HERE WE APPLY A TAPER TO FIX THE ENDPOINTS
    do i = 1,nlen

      i1 = istart - 1 + i

      if (ipick == 0) then
        adj_syn(i1,icomp,irec) = ( syn(i1,icomp,irec) -  data(i1,icomp,irec) ) * time_window(i) * meas_pert

      else if (ipick == 1) then
        ! meas_pert = 1.0 for most runs
        adj_syn(i1,icomp,irec) = -tshift_xc * ft_bar_t(i) * time_window(i) * meas_pert

      else if (ipick == 2) then
        adj_syn(i1,icomp,irec) = -dlna * fa_bar_t(i) * time_window(i) * meas_pert

      else if (ipick == 3) then
        adj_syn(i1,icomp,irec) = fp(i) * time_window(i)

      else if (ipick == 4) then
        adj_syn(i1,icomp,irec) = fq(i) * time_window(i)

      else if (ipick == 5) then
        adj_syn(i1,icomp,irec) = ft_bar_t(i) * time_window(i)

      else if (ipick == 6) then
        adj_syn(i1,icomp,irec) = fa_bar_t(i) * time_window(i)
      endif

    enddo

    ! (1) COMPUTE MISFIT function (currently only for waveform, xc-tt, xc-lnA)
    ! (2) COMPUTE MEASUREMENT VECTOR

    imeasure = imeasure + 1    ! global counter variable

      if (ipick == 0) then

        chi(ievent,irec,icomp,1) = 0.5 * sum( adj_syn(:,icomp,irec)**2 ) * DT
        measure_vec(imeasure)    = chi(ievent,irec,icomp,1)

      else if (ipick == 1) then
        chi(ievent,irec,icomp,1) = 0.5 * (tshift_xc * meas_pert)**2
        measure_vec(imeasure)    = tshift_xc * meas_pert

      else if (ipick == 2) then
        chi(ievent,irec,icomp,1) = 0.5 * (dlna * meas_pert)**2
        measure_vec(imeasure)    = dlna * meas_pert

      else if (ipick == 3) then
        chi(ievent,irec,icomp,1) = 0.
        measure_vec(imeasure)    = 0.

      else if (ipick == 4) then
        chi(ievent,irec,icomp,1) = 0.
        measure_vec(imeasure)    = 0.

      else if (ipick == 5) then
        chi(ievent,irec,icomp,1) = 0.
        measure_vec(imeasure)    = 0.

      else if (ipick == 6) then
        chi(ievent,irec,icomp,1) = 0.
        measure_vec(imeasure)    = 0.

      endif


  enddo  ! loop over receivers (irec=1,nrec)

!------------------------------------------------------------------
  end subroutine mtm_adj
!------------------------------------------------------------------
! END MAIN PROGRAM
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine clogc(n,xi,zzign,dt)
!------------------------------------------------------------------
      complex*16, dimension(NDIM) :: xi
      integer :: n
      double precision :: dt,zzign

      complex*16 :: wk, hold, q
      double precision :: m(25)
      double precision :: zign,flx,v
      integer :: lblock,k,fk,jh,ii,istart
      integer :: l,iblock,nblock,i,lbhalf,j,lx

      zign=zzign
      if (zign >= 0.) then
        zign=1.
      else
        zign=-1.
      endif
      lx=2**n
      do 1 i=1,n
    1 m(i)=2**(n-i)
      do 4 l=1,n
      nblock=2**(l-1)
      lblock=lx/nblock
      lbhalf=lblock/2
      k=0
      do 4 iblock=1,nblock
      fk=k
      flx=lx
      v=zign*2.*PI*fk/flx
      wk=cmplx(cos(v),sin(v))
      istart=lblock*(iblock-1)
      do 2 i=1,lbhalf
      j=istart+i
      jh=j+lbhalf
      q=xi(jh)*wk
      xi(jh)=xi(j)-q
      xi(j)=xi(j)+q
    2 continue
      do 3 i=2,n
      ii=i
      if (k < m(i)) goto 4
    3 k=k-m(i)
    4 k=k+m(ii)
      k=0
      do 7 j=1,lx
      if (k < j) goto 5
      hold=xi(j)
      xi(j)=xi(k+1)
      xi(k+1)=hold
    5 do 6 i=1,n
      ii=i
      if (k < m(i)) goto 7
    6 k=k-m(i)
    7 k=k+m(ii)
      if (zign > 0.) goto 9
      flx=flx*dt
      do 8 i=1,lx
    8 xi(i)=xi(i)/flx
      return
    9 do 10 i=1,lx
   10 xi(i)=xi(i)*dt
      return

  end subroutine clogc

!------------------------------------------------------------------
  subroutine ftinv(npow,s,zzign,dt,r)
!------------------------------------------------------------------

      !implicit real*8(a-h,o-z)
      !dimension r(4096*4)
      !complex s(4096*4)

      complex*16 :: s(NDIM)
      double precision :: r(NDIM)
      double precision :: dt,zzign,zign
      integer :: npow, nsmp, nhalf, i

      nsmp=2**npow
      nhalf=nsmp/2
      call rspec(s,nhalf)
      zign = -1.*zzign
      call clogc(npow,s,zign,dt)
      do 10 i=1,nsmp
10    r(i)=real(s(i))
      return

  end subroutine ftinv

!------------------------------------------------------------------
  subroutine rspec(s,np2)
!------------------------------------------------------------------

      !implicit real*8(a-h,o-z)
      !complex s(4096*4)

      complex*16 :: s(NDIM)
      integer :: np2,n,n1,i

      n=2*np2
      n1=np2+1

      s(n1)=0.
!     s(1)=0.
      s(1)= cmplx( real(s(1)),0.)

      do 20 i=1,np2
20    s(np2+i)=conjg(s(np2+2-i))
      return

  end subroutine rspec

!------------------------------------------------------------------
  subroutine remo(ny,nm,nd)
!------------------------------------------------------------------

      !implicit real*8(a-h,o-z)
      !dimension m(12)
      !data m/0,31,59,90,120,151,181,212,243,273,304,334/

      integer :: m(12)
      data m/0,31,59,90,120,151,181,212,243,273,304,334/
      integer :: ny, nm, nd, mm

      if (.not. (nm == 1))goto 23220
      nm=0
23220 continue
      mm=nm
      if (.not. (mm == 0))goto 23222
      return
23222 continue
      nm=0
      nd=nd+m(mm)
      if (.not. (mod(ny,4) == 0 .and. mm > 2))goto 23224
      nd=nd+1
23224 continue
      return

  end subroutine remo

!------------------------------------------------------------------
  subroutine staper(nt, fw, nev, v, ndim, a, w)
!------------------------------------------------------------------
!$$$$ calls tsturm, root
!  Slepian - Thomson multi-taper procedure
!  Slepian, D.     1978  Bell Sys Tech J v57 n5 1371-1430
!  Thomson, D. J.  1982  Proc IEEE v70 n9 1055-1096
!    nt    the number of points in the series
!    fw    the time-bandwidth product (number of Rayleigh bins)
!    nev   the desired number of tapers
!    v     the eigenvectors (tapers) are returned in v(.,nev)
!    a, w  work arrays dimensioned at least nt long (nt+1, nt odd)
!    a(1..nev) contains bandwidth retention factors on output.
!  The tapers are the eigenvectors of the tridiagonal matrix sigma(i,j)
!  [see Slepian(1978) eq 14 and 25.] They are also the eigenvectors of
!  the Toeplitz matrix eq. 18. We solve the tridiagonal system in
!  tsturm for the tapers and use them in Slepians eq 18 to get the
!  bandwidth retention factors (i.e. the eigenvalues) Thomson's
!  normalisation is used with no attention to sign.
      !implicit real*8(a-h,o-z)
      !dimension a(*),w(*),v(ndim,*)
      !parameter (pi=3.14159265358979d0,r2=1.414213562373095d0)

      integer :: nt, nev, ndim
      double precision :: fw
      double precision :: v(NDIM,MAXTAPER), a(NDIM), w(NDIM)

      integer :: i,j,k,m
      integer :: nxi, lh, lp1, neven, nodd, ntot, kk, kmax, nlow, nup
      double precision :: r2,om,com,hn,asav,rbd,dc,sm,s,sn,vmax

      !-------------------------

      r2 = sqrt(2.)

      if (nt < 2) return
      nxi=mod(nt,2)
      lh=(nt/2)+nxi
      lp1=nt+1
      om=2.*PI*fw/nt
      com=cos(om)
      hn=0.5*dble(lp1)
      do 10 i=1,lh
        a(i)=com*(i-hn)**2
   10   w(i)=0.5*dble(i*(nt-i))
      if (nxi == 0) then
        asav=a(lh)-w(lh)
        a(lh)=a(lh)+w(lh)
        rbd=1./(a(lh)+w(lh-1))
      else
        asav=w(lh-1)
        rbd=1./(w(lh)+w(lh-1))
        w(lh-1)=r2*w(lh-1)
      endif
      do 15 i=1,lh
        a(i+lh)=w(i)*rbd
        w(i)=a(i+lh)**2
   15   a(i)=a(i)*rbd
      neven=max0((nev+1)/2,1)
      nodd=nev-neven
!  Do the even tapers
      call tsturm(nt,lh,a,a(lh+1),w,neven,v,ndim,w(lh+1),0)
      do 20 i=1,neven
        k=2*i-1
        if (nxi == 1) v(lh,k)=r2*v(lh,k)
          do 20 j=1,lh
   20     v(lp1-j,k)=v(j,k)
      if (nodd <= 0) goto 34
!  Do the odd tapers
      if (nxi == 0) then
        a(lh)=asav*rbd
      else
        a(nt)=asav*rbd
        w(lh-1)=asav*asav
      endif
      call tsturm(nt,lh-nxi,a,a(lh+1),w,nodd,v,ndim,w(lh+1),1)
      do 30 i=1,nodd
        k=2*i
        if (nxi == 1) v(lh,k)=0.
          do 30 j=1,lh
   30     v(lp1-j,k)=-v(j,k)
   34 ntot=neven+nodd
!  Calculate bandwidth retention parameters
      dc=2.*com
      sm=0.
      s=sin(om)
      w(1)=om/PI
      w(2)=s/PI
      do 35 j=3,nt
        sn=dc*s-sm
        sm=s
        s=sn
   35   w(j)=s/(PI*(j-1))
      do 55 m=1,ntot
        vmax=abs(v(1,m))
        kmax=1
        do 40 kk=2,lh
          if (abs(v(kk,m)) <= vmax) goto 40
          kmax=kk
          vmax=abs(v(kk,m))
   40     continue
        a(m)=0.
        nlow=kmax-1
          do 45 j=1,nlow
   45     a(m)=a(m)+w(j+1)*v(nlow+1-j,m)
        nup=nt-nlow
          do 50 j=1,nup
   50     a(m)=a(m)+w(j)*v(nlow+j,m)
   55 a(m)=a(m)/v(kmax,m)
      return

  end subroutine staper

!------------------------------------------------------------------
  subroutine tsturm(nt,n,a,b,w,nev,r,ndim,ev,ipar)
!------------------------------------------------------------------
!$$$$ calls root
!  Uses bisection and Sturm counting to isolate the eigenvalues of the
!  symmetric tridiagonal matrix with main diagonal a(.) and sub/super
!  diagonal b(.).  Newton's method is used to refine the eigenvalue in
!  subroutine root then direct recursion is used to get the eigenvector
!  as this is always stable.  Note  ipar=0 for even tapers   =1 for odd
!  tapers
      !implicit real*8(a-h,o-z)
      !parameter (epsi=1.d-15,epsi1=5.d-15)
      !dimension a(*),b(*),ev(*),w(*),r(ndim,*)

      double precision, parameter :: epsi = 1.d-15, epsi1 = 5.d-15

      double precision, dimension(NDIM) :: a, b, w, ev
      double precision, dimension(NDIM,MAXTAPER) :: r
      integer :: nt,n,ndim,nev,ipar

      double precision, dimension(NDIM) :: bb
      double precision :: q,el,elam,u,umeps,x,ddot,rnorm
      integer :: i,j,ik,iag,m,jk,jm1

      !-------------------------

      if (n <= 0 .or. nev <= 0) return
      umeps=1.-epsi
      do 5 i=1,nev
    5 ev(i)=-1.
      u=1.
      do 1000 ik=1,nev
      if (ik > 1) u=ev(ik-1)*umeps
      el=min(ev(ik),u)
   10 elam=0.5*(u+el)
      if (abs(u-el) <= epsi1) goto 35
      iag=0
      q=a(1)-elam
      if (q >= 0.) iag=iag+1
      do 15 i=2,n
      if (q == 0.) x=abs(b(i-1))/epsi
      if (q /= 0.) x=w(i-1)/q
      q=a(i)-elam-x
      if (q >= 0.) iag=iag+1
      if (iag > nev) goto 20
   15 continue
      if (iag >= ik) goto 20
      u=elam
      goto 10
   20 if (iag == ik) goto 30
      m=ik+1
      do 25 i=m,iag
   25 ev(i)=elam
      el=elam
      goto 10
   30 el=elam
      call root(u,el,elam,a,b,w,n,ik)
   35 ev(ik)=elam
      jk=2*ik+ipar-1
      r(1,jk)=1.
      r(2,jk)=-(a(1)-ev(ik))/b(1)
      ddot=1.+r(2,jk)*r(2,jk)
      jm1=2
      do 45 j=3,n
      r(j,jk)=-((a(jm1)-ev(ik))*r(jm1,jk)+b(j-2)*r(j-2,jk))/b(jm1)
      ddot=ddot+r(j,jk)*r(j,jk)
   45 jm1=j
      rnorm=sqrt(nt/(2.*ddot))
      do 50 j=1,n
   50 r(j,jk)=r(j,jk)*rnorm
 1000 continue
      return

  end subroutine tsturm

!------------------------------------------------------------------
  subroutine root(u,el,elam,a,bb,w,n,ik)
!------------------------------------------------------------------

      !implicit real*8(a-h,o-z)
      !parameter (epsi = 1.d-15, epsi1 = 5.d-15)
      !dimension a(*),bb(*),w(*)

      double precision, parameter :: epsi = 1.d-15, epsi1 = 5.d-15
      double precision :: u,el,elam
      double precision, dimension(NDIM) :: a,bb,w
      integer :: n,ik

      double precision :: an,b,bm,bn,del,x
      integer :: i,iag

      !----------------------

    5 elam=0.5*(u+el)
   10 if (abs(u-el) <= 1.5*epsi1) return
      an=a(1)-elam
      b=0.
      bn=-1./an
      iag=0
      if (an >= 0.) iag=iag+1
      do 20 i=2,n
      if (an == 0.) x=abs(bb(i-1))/epsi
      if (an /= 0.) x=w(i-1)/an
      an=a(i)-elam-x
      if (an == 0.) an=epsi
      bm=b
      b=bn
      bn=((a(i)-elam)*b-bm*x-1.)/an
      if (an >= 0.) iag=iag+1
   20 continue
      if (iag == ik) goto 25
      u=elam
      goto 30
   25 el=elam
   30 del=1./bn
      if (abs(del) <= epsi1) del=sign(epsi1,del)
      elam=elam-del
      if (elam >= u .or. elam <= el) goto 5
      goto 10

  end subroutine root
!-------------------------------------------

end module wave2d_sub3
