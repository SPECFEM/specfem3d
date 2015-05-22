module ma_sub

  use ma_constants
  use ma_variables
  use ma_sub2
  use ascii_rw

  implicit none

contains

  subroutine mt_measure(datafile,filename,dat_dt,syn_dt,syn_dt_phydisp,t0,dt,npts,tstart,tend, &
         istart,dat_dtw,syn_dtw,syn_dtw_phydisp,nlen,tshift,sigma_dt_cc,dlnA,sigma_dlnA_cc,cc_max,syn_dtw_cc, &
         i_pmax_dat,i_pmax_syn,i_right,trans_w,dtau_w,dlnA_w,sigma_dt,sigma_dlnA,syn_dtw_mt, &
         err_dt,err_dlnA)

  ! ======================================================================
  ! subroutine mt_measure(): making measurements on windowed data and synthetics
  ! Boxcar/Cosine/Multitaper estimates of the transfer function between data and synthetics
  !
  !  Input:
  !        datafile --- note it is dummy right now (LQY)
  !        filename ---  output file prefix (e.g., OUT_DIR/PAS.CI.BHZ.02.mtm)
  !        dat_dt(:), syn_dt(:), t0, dt, npts  -- original data and synthetics array
  !        tstart, tend -- start and end of the measurement window (from FLEXWIN)
  !  Output:
  !        istart -- starting index of the windowed portion in the original trace
  !        dat_dtw(:), syn_dtw(:), nlen -- windowed and shifted data, windowed synthetics
  !        tshift, dlnA, cc_max -- time shift and amplitude cross-correlation measurements
  !        sigma_dt_cc, sigma_dlnA_cc -- corresponding cc error estimates
  !        syn_dtw_cc(:) -- reconstructed synthetics after Dt/DlnA corrections (LQY)
  !        i_pmax_dat, i_pmax_syn: the indices of max amplitude of the (freq.) power for data/syn
  !        i_right -- the maximum reliable frequency estimate index
  !        trans_w(:) -- estimates of transfer function  (MT only ???)
  !        dtau_w(:), dlnA_w(:) -- estimates of travel-time and amplitude anomaly
  !        sigma_dt, sigma_dlnA -- average travel time and amplitude error for mt (same as sigma_dt_cc)
  !        err_dt(:), err_dlnA(:) -- error bar of the travel-time and amplitude estimates (OPTIONAL)
  !
  !  other module variables used:
  !        is_mtm -- taper type: 1 for multitaper, 2 for cosine taper, and 3 or 0 for boxcar taper
  !  original coding in Fortran 77 by Ying Zhou
  !  upgraded to Fortran 90 by Alessia Maggi
  !  organized into package form by Qinya Liu
  !  modifications by Carl Tape and Vala Hjorleifsdottir
  !
  ! ====================================================================

    implicit none

    character(len=150), intent(in) :: filename,datafile
    double precision, dimension(:), intent(in) :: dat_dt,syn_dt, syn_dt_phydisp
    double precision, intent(in) ::  t0,dt,tstart,tend
    integer, intent(in) :: npts ! rarely used

    double precision, intent(out) :: tshift,sigma_dt_cc,dlnA,sigma_dlnA_cc,cc_max,sigma_dt,sigma_dlnA
    integer, intent(out) :: istart,nlen,i_pmax_dat,i_pmax_syn,i_right

    double precision, dimension(:), intent(out) :: syn_dtw,dat_dtw,syn_dtw_phydisp, syn_dtw_cc,syn_dtw_mt,dtau_w,dlnA_w
    complex*16, dimension(:), intent(out) :: trans_w
    double precision, dimension(:), intent(out), optional :: err_dt,err_dlnA


    ! local variables
    double precision, dimension(NPT) :: syn_dtw_mt_dt, &
         syn_dtw_cc_dt, dat_dtw_cc, syn_dtw_h, dat_dtw_h
    double precision :: sfac1,fac,f0,df,df_new,dw, &
         ampmax_unw,wtr_use_unw,ampmax,wtr_use,wtr_mtm,dtau_wa,dlnA_wa !omega
    integer :: ishift,i,ictaper,j,fnum,i_amp_max_unw,i_amp_max,i_right_stop,idf_new,iom,is_mtm_av

    complex*16, dimension(NPT) :: syn_dtwo, dat_dtwo, syn_dtw_ho, dat_dtw_ho,  &
                                  top_mtm, bot_mtm, trans_mtm
    double precision, dimension(NPT) :: wvec, ey1, ey2, dtau_mtm, dlnA_mtm, &
         phi_w, abs_w, err_phi, err_abs, phi_mtm, abs_mtm
    double precision :: eph_ave,edt_ave,eabs_ave,eabs2_ave,eph_iom,edt_iom,eabs_iom,eabs2_iom
    double precision, dimension(:,:),allocatable :: tas,phi_mul,abs_mul,dtau_mul,dlnA_mul

    !-------------------------------------------------------------
    ! check the window [tstart, tend] (the first two unnecessary)
    if ( tstart < t0 .or. tend > t0+(npts-1)*dt .or. tstart >= tend) then
       print *, 'tstart, t0, tend, t0+(npts-1)*dt:'
       print *, tstart, t0, tend, t0+(npts-1)*dt
       stop 'Check tstart and tend'
    else
       print *, '   start and end time of window: ' ,sngl(tstart),sngl(tend)
    endif

    ! initializes i_right
    i_right = 0

    ! water level in spectral divisions for mtm measurements
    wtr_mtm = 1.e-10 ! LQY -- move to par file? too small?

    if (DISPLAY_DETAILS) then
       print *, '   original data/syn length', npts
       call dwascii(trim(filename)//'.obs',dat_dt,npts,t0,dt)
       call dwascii(trim(filename)//'.syn',syn_dt,npts,t0,dt)
    endif
    !--------------------------------------------------------------------------
    ! window and interpolate data and synthetics
    !--------------------------------------------------------------------------

    call interpolate_dat_and_syn(dat_dt,syn_dt,syn_dt_phydisp,tstart,tend,t0,dt,NPT,dat_dtw,syn_dtw,syn_dtw_phydisp,nlen,istart)

    ! some constants
    sfac1 = (2./dble(nlen))**2   ! for Welch window
    ipwr_t = 10                  ! for time-domain cosine taper: 1 - [cos(t)]^(ipwr)

    ! pre-processing time-domain taper
    do i = 1,nlen
      !fac = 1.                                         ! boxcar window
      !fac = 1 - sfac1*((i-1) - dble(nlen)/2.)**2       ! welch window
      fac = 1. - cos(PI*(i-1)/(nlen-1))**ipwr_t        ! cosine window

      syn_dtw(i)  = syn_dtw(i) * fac  ! syn, windowed
      dat_dtw(i) = dat_dtw(i) * fac  ! dat, windowed
      if (USE_PHYSICAL_DISPERSION) then
          syn_dtw_phydisp(i)=syn_dtw_phydisp(i)*fac
      endif
    enddo

    if (DISPLAY_DETAILS) then
       print *, '   new nlen of dat/syn_dtw = ', nlen, '(',sngl(tstart),sngl(tend),')'
       call dwsac1(trim(filename)//'.obs.sac',dat_dtw,nlen,tstart,dt)
       call dwsac1(trim(filename)//'.syn.sac',syn_dtw,nlen,tstart,dt)
    endif

    !------------------------------------------------------------------
    ! cross-correlation traveltime and amplitude measurements
    !------------------------------------------------------------------

    ! compute cross-correlation time shift and also amplitude measurmement
    ! NOTE: records have already been windowed, so no information outside windows is considered
    ! LQY: Ying suggested to align them at relatively long periods first
    call compute_cc(syn_dtw, dat_dtw, nlen, dt, ishift, tshift, dlnA, cc_max)

    ! no need to compute velocity and acceleration of syn in mt_measure()

    ! deconstruct data using (negative) cross-correlation measurments
    call deconstruct_dat_cc(filename,dat_dtw,tstart,dt,nlen, &
        ishift,tshift,dlnA,dat_dtw_cc)

    ! reconstruct synthetics using cross-correlation measurements (plotting purposes only)
    ! syn_dtw_cc (reconstructed syn with both dt/dlnA), syn_dtw_cc_dt (with only dt)
    call reconstruct_syn_cc(syn_dtw,tstart,dt,nlen, &
         ishift,tshift,dlnA,syn_dtw_cc,syn_dtw_cc_dt)

    if (OUTPUT_MEASUREMENT_FILES) then
       call dwsac1(trim(filename)//'.recon_syn_cc.sac',syn_dtw_cc,nlen,tstart,dt)
       call dwsac1(trim(filename)//'.recon_syn_cc_dt.sac',syn_dtw_cc_dt,nlen,tstart,dt)
    endif

    ! compute the estimated uncertainty for the cross-correlation measurment
    sigma_dt_cc = 1.
    sigma_dlnA_cc = 1.
    call compute_average_error(dat_dtw,syn_dtw_cc,syn_dtw_cc_dt,nlen,dt,sigma_dt_cc,sigma_dlnA_cc)

    ! write cross-correlation measurement to file
    is_mtm_av=2
    call write_average_meas(filename,is_mtm_av,tshift,dlnA,sigma_dt_cc,sigma_dlnA_cc)

    !========================================

    ! CHT: if you want a simple waveform difference, then return
    if (is_mtm == 0) return

    !-----------------------------------------------------------------------------
    !  set up FFT for the frequency domain
    !-----------------------------------------------------------------------------

    ! calculate frequency step and number of frequencies used for FFT
    f0 = 0.
    df = 1./(NPT*dt)
    dw = TWOPI * df
    fnum = NPT/2 + 1

    ! calculate frequency spacing for the actual time window
    df_new = 1.0 / (tend-tstart)
    idf_new = df_new / df

    ! assemble omega vector (NPT is the FFT length)
    wvec(:) = 0.
    do j = 1,NPT
      if(j > NPT/2+1) then
        wvec(j) = dw*(j-NPT-1)   ! negative frequencies in second half
      else
        wvec(j) = dw*(j-1)       ! positive frequencies in first half
      endif
    enddo

    ! create complex synthetic seismogram and CC-deconstructed data seismogram
    syn_dtwo = cmplx(0.,0.)
    dat_dtwo = cmplx(0.,0.)
    syn_dtwo(1:nlen) = cmplx(syn_dtw(1:nlen))
    dat_dtwo(1:nlen) = cmplx(dat_dtw_cc(1:nlen))

    call fft(LNPT,syn_dtwo,FORWARD_FFT,dt)
    call fft(LNPT,dat_dtwo,FORWARD_FFT,dt)

    ! index of the freq of the max power in the windowed data/syn
    ampmax_unw = 0.
    i_pmax_dat = 1
    do i = 1, fnum   ! loop over frequencies
      if( abs(dat_dtwo(i)) > ampmax_unw) then
        ampmax_unw =  abs(dat_dtwo(i))
        i_pmax_dat = i
      endif
    enddo

    ampmax_unw = 0.
    i_amp_max_unw = 1
    do i = 1, fnum
      if( abs(syn_dtwo(i)) > ampmax_unw) then
        ampmax_unw =  abs(syn_dtwo(i))
        i_amp_max_unw = i
      endif
    enddo
    i_pmax_syn = i_amp_max_unw

    ! water level based untapered synthetics
    ! used to determine the i_right values (maximum frequency for measurement)
    wtr_use_unw = cmplx(ampmax_unw * WTR, 0.)

    i_right = fnum
    i_right_stop = 0
    do i = 1,fnum
      if( abs(syn_dtwo(i)) <= abs(wtr_use_unw) .and. i_right_stop==0 .and. i > i_amp_max_unw ) then
        i_right_stop = 1  ! power dips below water-level
        i_right = i
      endif
      if( abs(syn_dtwo(i)) >= 10.*abs(wtr_use_unw) .and. i_right_stop==1 .and. i > i_amp_max_unw) then
        i_right_stop = 0  ! power goes above 10*water-level
        i_right = i
      endif
    enddo

    if (DISPLAY_DETAILS) then
      print *, '   frequency of max power in windowed synthetic (Hz):'
      print *, '     i_pmax_syn = ', i_pmax_syn, ', f_pmax = ', sngl(i_pmax_syn * df), ', T_pmax = ', sngl(1./(i_pmax_syn*df))
      print *, '     FFT freq spacing df = ', sngl(df)
      print *, '     measurement freq spacing df_new = ', sngl(df_new)
      print *, '     i_right = ', i_right, ', stopping freq/period = ', sngl(i_right * df),'/', &
           sngl(1./(df*i_right))

      ! write out power for each signal
       call dwascii(trim(filename)//'.obs.power',abs(dat_dtwo(1:i_right)),i_right,df,df)
       call dwascii(trim(filename)//'.syn.power',abs(syn_dtwo(1:i_right)),i_right,df,df)
    endif

    !-------------------------------------------------------------------------------
    ! single-taper estimation of transfer function
    !-------------------------------------------------------------------------------

    ! assign number of tapers (is_mtm can only be 1, 2 or 3 at this point)
    if (is_mtm == 1) then
      ntaper = int(NPI * 2.0)
    else
      ntaper = 1
    endif
    allocate(tas(NPT,ntaper))

    ! calculate the tapers
    if (is_mtm == 1) then
      call staper(nlen, NPI, NTAPER, tas, NPT, ey1, ey2)
    else if (is_mtm == 2) then
      call costaper(nlen, NPT, tas)
    else if (is_mtm == 3) then
      call boxcar(nlen, NPT, tas)
    endif

    ! initialize transfer function terms
    top_mtm(:)   = cmplx(0.,0.)
    bot_mtm(:)   = cmplx(0.,0.)
    trans_mtm(:) = cmplx(0.,0.)

    do ictaper = 1, ntaper

      syn_dtw_ho(:) = cmplx(0.,0.) ! note: this has to be initialized inside the loop
      dat_dtw_ho(:) = cmplx(0.,0.)

      ! apply time-domain taper
      do i = 1, nlen
        syn_dtw_h(i) = syn_dtw(i) * tas(i,ictaper)     ! single-tapered, windowed syn
        dat_dtw_h(i) = dat_dtw_cc(i) * tas(i,ictaper)  ! single-tapered, windowed, shifted data
      enddo

      syn_dtw_ho(1:nlen) = cmplx(syn_dtw_h(1:nlen),0.)
      dat_dtw_ho(1:nlen) = cmplx(dat_dtw_h(1:nlen),0.)

      ! apply FFT to get complex spectra
      call fft(LNPT,syn_dtw_ho,FORWARD_FFT,dt)
      call fft(LNPT,dat_dtw_ho,FORWARD_FFT,dt)

      ! compute water level for single taper measurement by finding max spectral power
      ! in the single-tapered synthetics record
      ampmax = 0.
      i_amp_max = 1
      do i = 1, fnum
        if( abs(syn_dtw_ho(i)) > ampmax) then
          ampmax = abs(syn_dtw_ho(i))
          i_amp_max = i
        endif
      enddo
      wtr_use = cmplx(ampmax * WTR, 0.)

      ! calculate top and bottom of MT transfer function
      do i = 1, fnum
        top_mtm(i) = top_mtm(i) + dat_dtw_ho(i) * conjg(syn_dtw_ho(i))   ! uses data and syn
        bot_mtm(i) = bot_mtm(i) + syn_dtw_ho(i) * conjg(syn_dtw_ho(i))   ! uses syn only

        ! calculate transfer function for single taper measurement using water level
        if (is_mtm /= 1) then
          if (abs(syn_dtw_ho(i)) >  abs(wtr_use)) trans_w(i) = dat_dtw_ho(i) / syn_dtw_ho(i)
          if (abs(syn_dtw_ho(i)) <= abs(wtr_use)) trans_w(i) = dat_dtw_ho(i) / (syn_dtw_ho(i)+wtr_use)
        endif
      enddo

    enddo  ! ictapers

    ! for cosine or boxcar tapers only -- SEE COMMENTS BELOW for the multitaper case
    ! NOTE 1: here we are using trans_w, not trans_mtm
    ! NOTE 2: The single-taper transfer function should give you a perfect fit,
    !         but it is not relevant from the perspective of obtaining a measurement.
    if (is_mtm /= 1) then
       ! phase, abs(trans), travel-time and amplitude as a func of freq for single-tapered measurements
       call write_trans(filename,trans_w,wvec,fnum,i_right,idf_new,df,tshift,dlnA, &
            phi_w,abs_w,dtau_w,dlnA_w,dtau_wa,dlnA_wa)
       call reconstruct_syn(filename,syn_dtwo,wvec,dtau_w,dlnA_w, &
            i_right,tstart,dt,nlen,syn_dtw_mt, syn_dtw_mt_dt)
       !call check_recon_quality(filename,dat_dtw_cc,syn_dtw,dat_dtw,syn_dtw_mt,nlen,dt,tshift,tshift_f1f2,cc_max_f1f2,cc_max)

       return  ! over for single taper (no further error estimates)
    endif

    !-------------------------------------------------------------------------------
    ! multitaper estimation of transfer function
    !-------------------------------------------------------------------------------

    ! water level for multitaper measurements
    ampmax = 0.
    do i = 1, fnum
       if( abs(bot_mtm(i)) > ampmax) then
          ampmax =  abs(bot_mtm(i))
          i_amp_max = i
       endif
    enddo
    wtr_use = cmplx(ampmax * wtr_mtm**2, 0.)
    !wtr_use = cmplx(ampmax * WTR, 0.)

    ! calculate MT transfer function using water level
    do i = 1, fnum
       if(abs(bot_mtm(i)) > abs(wtr_use)) trans_mtm(i) = top_mtm(i) /  bot_mtm(i)
       if(abs(bot_mtm(i)) < abs(wtr_use)) trans_mtm(i) = top_mtm(i) / (bot_mtm(i)+wtr_use)
    enddo

    ! multitaper phase, abs, tt, and amp (freq)
    call write_trans(filename,trans_mtm,wvec,fnum,i_right,idf_new,df,tshift,dlnA, &
        phi_mtm,abs_mtm,dtau_mtm,dlnA_mtm,dtau_wa,dlnA_wa)

    ! apply transfer function to the syn
    call reconstruct_syn(filename,syn_dtwo,wvec,dtau_mtm,dlnA_mtm, &
        i_right,tstart,dt,nlen,syn_dtw_mt,syn_dtw_mt_dt)

    ! check quality
    !call check_recon_quality(filename,dat_dtw_cc,syn_dtw,dat_dtw,syn_dtw_mt,nlen,dt,tshift, tshift_f1f2, cc_max_f1f2,cc_max)

    ! CHT: estimate error using the same procedure as for the cross-correlation error estimate
    sigma_dt = sigma_dt_cc  ;  sigma_dlnA = sigma_dlnA_cc

    ! write average multitaper measurement to file
    is_mtm_av=1
    call write_average_meas(filename, is_mtm_av, dtau_wa, dlnA_wa, sigma_dt, sigma_dlnA)

    !-------------------------------------------------------------------------------
    ! multitaper error estimation
    !-------------------------------------------------------------------------------

    if (ntaper > 1) then  ! this should always be true

      ! no need to save a copy of the control logicals any more (output/display)

      ! allocate Jacknife MT estimates
      allocate(phi_mul(NPT,ntaper))
      allocate(abs_mul(NPT,ntaper))
      allocate(dtau_mul(NPT,ntaper))
      allocate(dlnA_mul(NPT,ntaper))

      do iom = 1, ntaper

        top_mtm(:) = cmplx(0.,0.)
        bot_mtm(:) = cmplx(0.,0.)

        do ictaper = 1, ntaper
          if(ictaper==iom) cycle

          ! apply ictaper-th taper
          syn_dtw_h(1:nlen) = syn_dtw(1:nlen) * tas(1:nlen,ictaper)
          dat_dtw_h(1:nlen) = dat_dtw_cc(1:nlen) * tas(1:nlen,ictaper)

          ! complex tapered series
          syn_dtw_ho(:) = cmplx(0.,0.)
          dat_dtw_ho(:) = cmplx(0.,0.)
          syn_dtw_ho(1:nlen) = cmplx(syn_dtw_h(1:nlen),0.)
          dat_dtw_ho(1:nlen) = cmplx(dat_dtw_h(1:nlen),0.)

          ! apply f.t. to get complex spectra
          call fft(LNPT,syn_dtw_ho,FORWARD_FFT,dt)
          call fft(LNPT,dat_dtw_ho,FORWARD_FFT,dt)

          ! calculate top and bottom of Jacknife transfer function
          do i = 1, fnum
            top_mtm(i) = top_mtm(i) + dat_dtw_ho(i) * conjg(syn_dtw_ho(i))
            bot_mtm(i) = bot_mtm(i) + syn_dtw_ho(i) * conjg(syn_dtw_ho(i))
          enddo
        enddo ! ictaper

        ! water level
        ampmax = 0.
        do i = 1, fnum
          if( abs(bot_mtm(i))>ampmax) then
            ampmax =  abs(bot_mtm(i))
            i_amp_max = i
          endif
        enddo
        wtr_use = cmplx(ampmax * wtr_mtm ** 2, 0.)  ! again very small wtr_mtm is used

        !  calculate transfer function using water level
        do i = 1, fnum
          if(abs(bot_mtm(i))>abs(wtr_use)) trans_mtm(i) = top_mtm(i) / bot_mtm(i)
          if(abs(bot_mtm(i))<=abs(wtr_use)) trans_mtm(i) = top_mtm(i) /(bot_mtm(i)+wtr_use)
        enddo

        call write_trans(filename,trans_mtm,wvec,fnum,i_right,idf_new,df,tshift,dlnA, &
            phi_mul(:,iom),abs_mul(:,iom),dtau_mul(:,iom),dlnA_mul(:,iom))

      enddo ! ntaper

      !----------------------

      if (OUTPUT_MEASUREMENT_FILES) then
         open(10,file=trim(filename)//'.err_ph')
         open(20,file=trim(filename)//'.err_dt')
         open(30,file=trim(filename)//'.err_abs')
         open(40,file=trim(filename)//'.err_dlnA')

      ! CHT: Since all freq. domain points are used in constructing the
      !      adjoint source, we also want to show the entire sigma(f) functions,
      !      not just the sub-sampled version.
         open(50,file=trim(filename)//'.err_dt_full')
         open(60,file=trim(filename)//'.err_dlnA_full')
      endif

      err_phi  = 0.
      err_dt   = 0.
      err_abs  = 0.
      err_dlnA = 0.

      do i = 1, i_right

          eph_ave   = 0.
          edt_ave   = 0.
          eabs_ave  = 0.
          eabs2_ave = 0.

          do iom = 1, ntaper
            eph_iom = ntaper*phi_mtm(i) - (ntaper-1)*phi_mul(i,iom)
            eph_ave = eph_ave + eph_iom

            edt_iom = ntaper*dtau_mtm(i) - (ntaper-1)*dtau_mul(i,iom)
            edt_ave = edt_ave + edt_iom

            eabs_iom = ntaper*abs_mtm(i) - (ntaper-1)*abs_mul(i,iom)
            eabs_ave = eabs_ave + eabs_iom

            eabs2_iom = ntaper*dlnA_mtm(i) - (ntaper-1)*dlnA_mul(i,iom)
            eabs2_ave = eabs2_ave + eabs2_iom
          enddo

          eph_ave   = eph_ave   / (ntaper)
          edt_ave   = edt_ave   / (ntaper)
          eabs_ave  = eabs_ave  / (ntaper)
          eabs2_ave = eabs2_ave / (ntaper)

          do iom = 1, ntaper
            err_phi(i)  = err_phi(i) + ( phi_mul(i,iom) - eph_ave)**2
            err_dt(i)   = err_dt(i)  + (dtau_mul(i,iom) - edt_ave)**2
            err_abs(i)  = err_abs(i) + ( abs_mul(i,iom) - eabs_ave)**2
            err_dlnA(i) = err_dlnA(i)+ (dlnA_mul(i,iom) - eabs2_ave)**2
          enddo

          err_phi(i)  =  sqrt( err_phi(i) / (ntaper * (ntaper-1) ) )
          err_dt(i)   =  sqrt( err_dt(i) / (ntaper * (ntaper-1) ) )
        ! set the error bar for the first point corresponding to
        ! static offset to be large, which makes no contribution to
        ! the adjoint source
          if (i == 1) err_dt(i) = LARGE_VAL
          err_abs(i)  =  sqrt( err_abs(i) / (ntaper * (ntaper-1) ) )
          err_dlnA(i) =  sqrt( err_dlnA(i) / (ntaper * (ntaper-1) ) )

        ! only write out the errors for the 'independent' freq-domain sampling points
        if (mod(i,idf_new) == 0 .and. OUTPUT_MEASUREMENT_FILES) then
          write(10,*) df*(i-1), phi_mtm(i), err_phi(i)
          if (i > 1) write(20,*) df*(i-1), dtau_mtm(i), err_dt(i)
          write(30,*) df*(i-1), abs_mtm(i), err_abs(i)
          write(40,*) df*(i-1), dlnA_mtm(i), err_dlnA(i)
        endif

        ! CHT: write out the entire dt(f) and dlnA(f) for adjoint sources
        if (OUTPUT_MEASUREMENT_FILES) then
           write(50,*) df*(i-1), dtau_mtm(i), err_dt(i)
           write(60,*) df*(i-1), dlnA_mtm(i), err_dlnA(i)
        endif

      enddo ! i_right

      close(10)
      close(20)
      close(30)
      close(40)
      close(50)
      close(60)

      ! pass the MT transfer funnction
      trans_w = trans_mtm
      dtau_w = dtau_mtm
      dlnA_w = dlnA_mtm

    endif

    !     ------------------------------------------------------------------
    !     End error calculation
    !     ------------------------------------------------------------------

  end subroutine mt_measure


  ! =====================================================================================================
  ! subroutine mt_adj()
  ! Compute cross-correlation travel-time/amplitude/banana-donut travel-time/banana-donut amplitude
  ! adjoint sources by assimulate the measurements passed from mt_measure()
  !
  !    Input:
  !      istart -- starting index of the windowed portion of original trace, used to assign
  !                tr/amp_adj_src(:) corresponding to the original synthetics
  !      dat_dtw(:), syn_dtw(:), nlen, dt -- windowed data and synthetics
  !                                           with length nlen and sampling rate dt
  !      tshift, dlnA, sigma_dt_cc,sigma_dlnA_cc -- CC traveltime/amplitude measurements/error estimates
  !      dtau_w(:), dlnA_w(:), err_dtau(:), err_dlnA(:) -- frequency-dependent tt/amp measurements/error
  !      i_left, i_right -- range of indices of valid freqs for frequency-dependent tt/amp measurements
  !
  !    Output:
  !      window_chi(:) -- all available scalar measurement values and chi values (20 in total)
  !      tr_adj_src(:), tr_chi -- travel-time adjoint source and chi value
  !      am_adj_src(:), am_chi -- amplitude adjoint source and chi value (optional)

  !    other global vars used:
  !      imeas -- adjoint source type: 1/2 for waveform, 3/4 for cc/banana-doughnut, 5/6 for cc, 7/8 for mtm
  !      this subroutine is only called when COMPUTE_ADJOINT_SOURCE = true.
  !
  !    original coding by Carl Tape and finalized by Qinya Liu
  ! ======================================================================================================

  subroutine mt_adj(istart,dat_dtw,syn_dtw,syn_dtw_phydisp,nlen,dt,tshift,dlnA,sigma_dt_cc,sigma_dlnA_cc, &
         dtau_w,dlnA_w,err_dtau,err_dlnA,sigma_dt,sigma_dlnA,i_left,i_right, &
         window_chi,tr_adj_src,tr_chi,am_adj_src,am_chi)

    implicit none
    integer, intent(in) :: istart, nlen, i_left, i_right
    double precision, dimension(:), intent(in) :: dat_dtw,syn_dtw,syn_dtw_phydisp
    double precision, intent(in) :: dt, tshift, dlnA, sigma_dt_cc, sigma_dlnA_cc, sigma_dt, sigma_dlnA
    double precision, dimension(:), intent(in) :: dtau_w, dlnA_w, err_dtau, err_dlnA

    double precision, dimension(NCHI), intent(inout) :: window_chi
    double precision, dimension(:), intent(out) :: tr_adj_src, am_adj_src
    double precision, intent(out) :: tr_chi, am_chi
    !double precision, dimension(:), intent(out), optional :: am_adj_src
    !double precision, intent(out), optional :: am_chi

    double precision, dimension(NPT) :: syn_vtw, syn_vtw_h, syn_dtw_h, ey1, ey2
    double precision, dimension(NPT) :: ft_bar_t, fa_bar_t, fp, fq, wp_taper, wq_taper
    complex*16, dimension(NPT) :: d_bot_mtm, v_bot_mtm
    integer :: i, i1, ictaper, ntaper
    double precision :: df,Nnorm,Mnorm,fac,ffac,w_taper(NPT), time_window(NPT)
    double precision, dimension(:,:), allocatable :: tas
    complex*16, dimension(:,:),allocatable :: syn_dtw_ho_all, syn_vtw_ho_all
    complex*16, dimension(NPT) :: pwc_adj,qwc_adj
    double precision, dimension(NPT) :: dtau_pj_t, dlnA_qj_t
    double precision :: dtau_wtr, dlnA_wtr, err_t, err_A
    double precision :: waveform_chi, waveform_d2, waveform_s2, waveform_temp1, waveform_temp2, waveform_temp3

    ! waveform adjoint source is passed by tr_adj_src and tr_chi
    !if (imeas == 0 .and. (present(am_adj_src) .or. present(am_chi))) stop  &
    !   'am_adj_src and am_chi are not needed for imeas = 0 (waveform adjoint source case)'

    ! check the window length
    if (istart + nlen > NDIM) stop 'Check istart + nlen and NPT'

    ! waveform
    if(imeas==1 .or. imeas==2) then
       print *, '     computing waveform adjoint source'
    else if(imeas==3 .or. imeas==4) then
       print *, '     computing banana-doughtnut adjoint source'
    else if(imeas==5 .or. imeas==6) then
       print *, '     computing cross-correlation adjoint source'
    else if(imeas==7 .or. imeas==8) then
       print *, '     computing multitaper adjoint source'
    endif

    ! define post-processing time-domain taper
    ! NOTE: If the adjoint sources will be band-pass filtered at the end,
    !       then perhaps time_window is not necessary (i.e., use boxcar).
    !       However, if you are using a waveform difference, then you want
    !       to make sure that the endpoints of the windows are at zero, since
    !       you would NOT apply the post-processing band-pass filter.
    ! Note also that the dat/syn_dtw(:) has already been time-windowed
    time_window(:) = 0.
    ipwr_t = 10
    do i = 1,nlen
      fac = 1.                                           ! boxcar window
      !fac = 1 - sfac2*((i-1) - dble(nlen1)/2.0)**2       ! welch window
      !fac = 1. - cos(PI*(i-1)/(nlen-1))**ipwr_t          ! cosine window
      time_window(i) = fac
    enddo

    ! ----------------------------------
    ! CROSS CORRELATION ADJOINT SOURCES LQY: does time_window needs to be applied here?
    ! ----------------------------------
    if( (imeas >= 3).and.(imeas <= 6) ) then

      ! compute synthetic velocity
      if (USE_PHYSICAL_DISPERSION) then
        do i = 2, nlen-1
                syn_vtw(i) = (syn_dtw_phydisp(i+1) - syn_dtw_phydisp(i-1)) / (2.0*dt)
        enddo
        syn_vtw(1)    = (syn_dtw_phydisp(2) - syn_dtw_phydisp(1)) / dt
        syn_vtw(nlen) = (syn_dtw_phydisp(nlen) - syn_dtw_phydisp(nlen-1)) / dt

        ! compute CC traveltime and amplitude banana-dougnut kernels
        ft_bar_t = 0.
        Nnorm = dt * sum( syn_vtw(1:nlen) * syn_vtw(1:nlen) )
        ft_bar_t(1:nlen) = -syn_vtw(1:nlen) / Nnorm

        fa_bar_t = 0.
        Mnorm = dt * sum( syn_dtw_phydisp(1:nlen) * syn_dtw_phydisp(1:nlen) )
        fa_bar_t(1:nlen) = syn_dtw_phydisp(1:nlen) / Mnorm
      else
        do i = 2, nlen-1
                syn_vtw(i) = (syn_dtw(i+1) - syn_dtw(i-1)) / (2.0*dt)
        enddo
        syn_vtw(1)    = (syn_dtw(2) - syn_dtw(1)) / dt
        syn_vtw(nlen) = (syn_dtw(nlen) - syn_dtw(nlen-1)) / dt

        ! compute CC traveltime and amplitude banana-dougnut kernels
        ft_bar_t = 0.
        Nnorm = dt * sum( syn_vtw(1:nlen) * syn_vtw(1:nlen) )
        ft_bar_t(1:nlen) = -syn_vtw(1:nlen) / Nnorm

        fa_bar_t = 0.
        Mnorm = dt * sum( syn_dtw(1:nlen) * syn_dtw(1:nlen) )
        fa_bar_t(1:nlen) = syn_dtw(1:nlen) / Mnorm
      endif
    endif

    ! ----------------------------------------------
    ! FREQUENCY-DOMAIN TAPERS FOR MT ADJOINT SOURCES
    ! ----------------------------------------------
    if( is_mtm == 1 ) then

      ! initialize water levels for err_dtau/dlnA division
      dtau_wtr = WTR * sum(abs(dtau_w(i_left:i_right)))/(i_right-i_left)  ! CHT i_left
      dlnA_wtr = WTR * sum(abs(dlnA_w(i_left:i_right)))/(i_right-i_left)  ! CHT i_left

      ! frequency-domain tapers
      ! THIS CHOICE WILL HAVE AN EFFECT ON THE ADJOINT SOURCES
      ipwr_w = 10
      w_taper(:) = 0.
      do i = i_left, i_right    ! CHT: 1 --> i_left
        ! type of filter in the freq domain
        !w_taper(i) = 1.                                       ! boxcar
        !w_taper(i) = 1. - (2.0/nw)**2 * ((i-1) - nw/2.0)**2     ! welch
        w_taper(i) = 1. - cos(PI*(i-i_left)/(i_right-i_left))**ipwr_w    ! cosine
      enddo

      ! compute normalization factor for w_taper
      ! note: 2 is needed for the integration from -inf to inf
      df = 1. /(NPT*dt)
      ffac = 2.0 * df * sum(w_taper(i_left:i_right) )   ! CHT: 1 --> i_left
      if (DISPLAY_DETAILS) print *, '     f-dom taper normalization factor, ffac = ', sngl(ffac)

      ! wp_taper and wq_taper are modified frequency-domain tapers
      ! Notice the option to include the frequency-dependent error.
      wp_taper(:) = 0.
      wq_taper(:) = 0.

      do i = i_left, i_right    ! CHT: i_left = 1

        if (ERROR_TYPE == 0 .or. DO_RAY_DENSITY_SOURCE ) then
          ! no error estimate
          ! only adds normalization factor
          wp_taper(i) = w_taper(i) / ffac
          wq_taper(i) = w_taper(i) / ffac

        else if (ERROR_TYPE == 1) then  ! CC error as a constant for all freqs.
          ! MT error estimate is assigned the CC error estimate
          wp_taper(i) = w_taper(i) / ffac / (sigma_dt ** 2)
          wq_taper(i) = w_taper(i) / ffac / (sigma_dlnA ** 2)

        else if (ERROR_TYPE == 2) then  ! MT jack-knife error estimate
          err_t = err_dtau(i)
          if (err_dtau(i) < dtau_wtr)  err_t = err_t + dtau_wtr
          err_A = err_dlnA(i)
          if (err_dlnA(i) < dlnA_wtr)  err_A = err_A + dlnA_wtr
          wp_taper(i) = w_taper(i) / ffac / (err_t ** 2)
          wq_taper(i) = w_taper(i) / ffac / (err_A ** 2)
        endif
      enddo

!!$    open(88,file='ftaper.dat')
!!$    do i = 1,i_right
!!$       write(88,'(5e18.6)') df*i, w_taper(i), dtau_w(i), dtau_w(i)*w_taper(i), dtau_w(i)*wp_taper(i)
!!$    enddo
!!$    close(88)

      ! allocate MT variables
      ntaper = int(NPI * 2.0)
      allocate(tas(NPT,ntaper))
      allocate(syn_dtw_ho_all(NPT,ntaper))
      allocate(syn_vtw_ho_all(NPT,ntaper))

      ! get the MT tapers
      call staper(nlen, NPI, NTAPER, tas, NPT, ey1, ey2)

      d_bot_mtm = 0.
      v_bot_mtm = 0.

      ! compute the bot required to compute p_j's and q_j's
      do ictaper = 1,ntaper

        ! tapered synthetic displacement
        if (USE_PHYSICAL_DISPERSION) then
                syn_dtw_h(1:nlen) = syn_dtw_phydisp(1:nlen) * tas(1:nlen,ictaper)
        else
                syn_dtw_h(1:nlen) = syn_dtw(1:nlen) * tas(1:nlen,ictaper)
        endif
        ! compute velocity of tapered syn
        do i = 2, nlen-1
          syn_vtw_h(i) = (syn_dtw_h(i+1) - syn_dtw_h(i-1)) / (2.0*dt)
        enddo
        syn_vtw_h(1)    = (syn_dtw_h(2) - syn_dtw_h(1)) / dt
        syn_vtw_h(nlen) = (syn_dtw_h(nlen) - syn_dtw_h(nlen-1)) /dt

        ! single-tapered complex synthetic displacement and velocity
        syn_dtw_ho_all(:,ictaper) = 0.
        syn_vtw_ho_all(:,ictaper) = 0.
        syn_dtw_ho_all(1:nlen,ictaper) = cmplx(syn_dtw_h(1:nlen),0.)
        syn_vtw_ho_all(1:nlen,ictaper) = cmplx(syn_vtw_h(1:nlen),0.)

        ! apply FFT get complex spectra LQY: can't we use disp->velo to save one fft operation?
        call fft(LNPT,syn_dtw_ho_all(:,ictaper),FORWARD_FFT,DT)
        call fft(LNPT,syn_vtw_ho_all(:,ictaper),FORWARD_FFT,DT)

        d_bot_mtm(:) = d_bot_mtm(:) + syn_dtw_ho_all(:,ictaper) * conjg(syn_dtw_ho_all(:,ictaper))
        v_bot_mtm(:) = v_bot_mtm(:) + syn_vtw_ho_all(:,ictaper) * conjg(syn_vtw_ho_all(:,ictaper))

      enddo ! ictaper

      ! compute p_j, q_j, P_j, Q_j and adjoint source fp, fq
      fp = 0.
      fq = 0.
      do ictaper = 1,ntaper

        ! compute p_j(w) and q_j(w)
        pwc_adj(:) = cmplx(0.,0.)
        qwc_adj(:) = cmplx(0.,0.)

        do i = 1, i_right ! LQY: i_left = 1
          pwc_adj(i) =  syn_vtw_ho_all(i,ictaper) / v_bot_mtm(i)
          qwc_adj(i) = -syn_dtw_ho_all(i,ictaper) / d_bot_mtm(i)
        enddo

        ! compute P_j(w) and Q_j(w)
        ! NOTE: the MT measurement is incorporated here
        ! also note that wp_taper and wq_taper can contain uncertainty estimations
        if( DO_RAY_DENSITY_SOURCE ) then
          ! uses a misfit measurement dtau, dlnA  = 1
          pwc_adj(:) = pwc_adj(:) * cmplx(1.0,0.) * cmplx(wp_taper(:),0.)
          qwc_adj(:) = qwc_adj(:) * cmplx(1.0,0.) * cmplx(wq_taper(:),0.)
        else
          ! adds misfit measurement dtau, dlnA
          pwc_adj(:) = pwc_adj(:) * cmplx(dtau_w(:),0.) * cmplx(wp_taper(:),0.)
          qwc_adj(:) = qwc_adj(:) * cmplx(dlnA_w(:),0.) * cmplx(wq_taper(:),0.)
        endif

        ! IFFT into the time domain
        call fftinv(LNPT,pwc_adj,REVERSE_FFT,dt,dtau_pj_t)
        call fftinv(LNPT,qwc_adj,REVERSE_FFT,dt,dlnA_qj_t)

        ! create adjoint source: applies taper to time signal
        fp(:) = fp(:) + tas(:,ictaper) * dtau_pj_t(:)
        fq(:) = fq(:) + tas(:,ictaper) * dlnA_qj_t(:)

      enddo

    endif ! MT adjoint source

    ! -------------------------------------
    !  ASSEMBLE VARIOUS ADJOINT SOURCES
    ! -------------------------------------

    tr_adj_src = 0.
    am_adj_src = 0.

    ! integrated waveform difference squared
    waveform_temp1 = 0. ; waveform_temp2 = 0. ; waveform_temp3 = 0.
    do i = 1,nlen
       waveform_temp1 = waveform_temp1 + ( dat_dtw(i) * time_window(i) )**2
       waveform_temp2 = waveform_temp2 + ( syn_dtw(i) * time_window(i) )**2
       waveform_temp3 = waveform_temp3 + (( dat_dtw(i) - syn_dtw(i) ) * time_window(i) )**2
    enddo
    ! NOTE: does not include DT factor or normalization by duration of window
    waveform_d2  = waveform_temp1
    waveform_s2  = waveform_temp2
    waveform_chi = waveform_temp3

    ! compute traveltime and amplitude adjoint sources for imeas
    do i = 1,nlen
       i1 = istart + i -1 ! start index in the full adjoint source array(1:npts)

       ! waveform
       if(imeas==1 .or. imeas==2) then
          tr_adj_src(i1) = -dat_dtw(i)/waveform_d2 * time_window(i) ! imeas=1
          ! consider normalizing this by waveform_d2
          am_adj_src(i1) = ( syn_dtw(i) - dat_dtw(i) ) * time_window(i) ! imeas=2

          ! use pure data waveform in time window
          if( NO_WAVEFORM_DIFFERENCE ) then
             tr_adj_src(i1) = dat_dtw(i) * time_window(i) ! waveform misfit (imeas=1)
          endif

       ! banana-doughnut kernel adjoint source (no measurement)
       else if(imeas==3 .or. imeas==4) then
          tr_adj_src(i1) = ft_bar_t(i) * time_window(i)  ! imeas=3
          am_adj_src(i1) = fa_bar_t(i) * time_window(i)  ! imreas=4

       ! cross-correlation
       else if(imeas==5 .or. imeas==6) then
          tr_adj_src(i1) = -(tshift / sigma_dt_cc**2 ) * ft_bar_t(i) * time_window(i)
          am_adj_src(i1) = -(dlnA / sigma_dlnA_cc**2 ) * fa_bar_t(i) * time_window(i)

          ! ray density
          if( DO_RAY_DENSITY_SOURCE ) then
             ! uses a misfit measurement of 1
             tr_adj_src(i1) = - (1.0) * ft_bar_t(i) * time_window(i)
             am_adj_src(i1) = - (1.0) * fa_bar_t(i) * time_window(i)
          endif

       ! multitaper
       else if(imeas==7 .or. imeas==8) then
          tr_adj_src(i1) = fp(i) * time_window(i)
          am_adj_src(i1) = fq(i) * time_window(i)
       endif
    enddo

    ! -------------------------------------
    !  COMPUTE MISFIT FUNCTION VALUE
    ! -------------------------------------

    ! CHT: compute misfit function value and measurement value
    ! Note: The taper functions for MT may include error estimates.
    ! 1: multitaper, TT
    ! 2: multitaper, dlnA
    ! 3: cross-correlation, TT
    ! 4: cross-correlation, dlnA
    !window_chi(:) = 0.

    ! misfit function value
    if(is_mtm==1) window_chi(1) = 0.5 * 2.0 * df * sum( (dtau_w(1:i_right))**2 * wp_taper(1:i_right) )
    if(is_mtm==1) window_chi(2) = 0.5 * 2.0 * df * sum( (dlnA_w(1:i_right))**2 * wq_taper(1:i_right) )
    window_chi(3) = 0.5 * (tshift/sigma_dt_cc)**2
    window_chi(4) = 0.5 * (dlnA/sigma_dlnA_cc)**2

    ! cc/averaged mt measurement (no uncertainty estimates)
    if(is_mtm==1) window_chi(5)  = sum( dtau_w(1:i_right) * w_taper(1:i_right) ) / sum(w_taper(1:i_right) )
    if(is_mtm==1) window_chi(6)  = sum( dlnA_w(1:i_right) * w_taper(1:i_right) ) / sum(w_taper(1:i_right) )
    window_chi(7) = tshift
    window_chi(8) = dlnA

    ! replaces misfit function values
    if( DO_RAY_DENSITY_SOURCE ) then
      ! uses misfit measurements equal to 1
      if(is_mtm==1) window_chi(1) = 0.5 * 2.0 * df * sum( (1.0)**2 * wp_taper(1:i_right) )
      if(is_mtm==1) window_chi(2) = 0.5 * 2.0 * df * sum( (1.0)**2 * wq_taper(1:i_right) )
      window_chi(3) = 0.5 * (1.0)**2
      window_chi(4) = 0.5 * (1.0)**2

      if(is_mtm==1) window_chi(5)  = sum( 1.0 * w_taper(1:i_right) ) / sum(w_taper(1:i_right) )
      if(is_mtm==1) window_chi(6)  = sum( 1.0 * w_taper(1:i_right) ) / sum(w_taper(1:i_right) )
      window_chi(7) = 1.0
      window_chi(8) = 1.0
    endif

    ! estimated measurement uncertainties
    if(is_mtm==1) window_chi(9) = sigma_dt
    if(is_mtm==1) window_chi(10) = sigma_dlnA
    window_chi(11) = sigma_dt_cc
    window_chi(12) = sigma_dlnA_cc

    ! for normalization, divide by duration of window
    window_chi(13) = 0.5 * waveform_d2
    window_chi(14) = 0.5 * waveform_s2
    window_chi(15) = 0.5 * waveform_chi
    window_chi(16) = nlen*dt

    if(imeas <= 2) then           ! waveform
      tr_chi = 0.5 * waveform_chi
      am_chi = 0.5 * waveform_chi

    else if( (imeas >= 3).and.(imeas <= 6) ) then  ! cross_correlation
      tr_chi = window_chi(3)
      am_chi = window_chi(4)

    else if( (imeas==7).or.(imeas==8) ) then       ! multitaper
      tr_chi = window_chi(1)
      am_chi = window_chi(2)

    endif

  end subroutine mt_adj

  !==============================================================================
  ! OTHER MAJOR SURBROUTINES
  !==============================================================================

  subroutine bandpass(x,n,delta_t,f1,f2)
    ! a double-precision wrapper around sac xapiir()
    ! modified from FLEXWIN subroutines on 26-July-2009

    implicit none
    integer, intent(in) :: n
    double precision, intent(inout),  dimension(*) :: x
    double precision, intent(in) :: delta_t,f1,f2
    real, dimension(:), allocatable :: x_sngl

    allocate(x_sngl(n))

    x_sngl(1:n) = sngl(x(1:n))
    !  delta_t_sngl = sngl(delta_t)

    ! old version - uses old SacLib
    ! does band-pass filter
    !call xapiir(x_sngl,n,'BU',sngl(TRBDNDW),sngl(APARM),IORD,'BP',sngl(FSTART),sngl(FEND),delta_t_sngl,PASSES)

    ! new version, uses subroutines in libsac.a
    ! does band-pass filter
    ! BU - butterworth
    ! BP - bandpass
    ! LQY: Shouldn't  delta_t_sngl = sngl(delta_t) still be done? same for f1,f2?
    call xapiir(x_sngl,n,'BU',TRBDNDW,APARM,IORD,'BP',f1,f2,delta_t,PASSES)

    x(1:n) = dble(x_sngl(1:n))

    deallocate(x_sngl)

  end subroutine bandpass

  !-----------------------------------------------------------------------------

  subroutine drsac1(datafile,data,npt1,b1,dt1)
    ! read sac file and convert to double precision

    implicit none
    character(len=*),intent(in) :: datafile
    real, dimension(NDIM) :: dat_sngl
    double precision, dimension(NDIM), intent(out) :: data
    integer :: npt1, nerr
    real :: b1_sngl,dt1_sngl
    double precision :: b1,dt1

    ! read file as single precision
    call rsac1(datafile,dat_sngl,npt1,b1_sngl,dt1_sngl,NDIM,nerr)
    if (nerr > 0) then
       print *, 'Error reading sac file', trim(datafile)
       stop
    endif

    ! return double precision quantities
    b1 = dble(b1_sngl)
    dt1 = dble(dt1_sngl)
    data = dble(dat_sngl)

  end subroutine drsac1

  !-----------------------------------------------------------------------------

  subroutine dwsac1(datafile,data,npt1,b1,dt1)
    ! convert to single precision, then write sac file
    ! --> includes an option to add minmax values to sac file,
    !     which are used in the plotting scripts

    implicit none
    character(len=*),intent(in) :: datafile
    integer, intent(in) :: npt1
    double precision, dimension(npt1), intent(in) :: data
    double precision, intent(in) :: b1,dt1
    logical, parameter :: minmax_header = .true.

    real, dimension(npt1) :: dat_sngl,ti_sngl
    real :: b1_sngl,dt1_sngl,xmin_sngl,xmax_sngl
    integer :: nerr,i

    ! convert to single precision
    b1_sngl = real(b1)
    dt1_sngl = real(dt1)
    dat_sngl = real(data)

    if (minmax_header) then
       ! get time vector
       ti_sngl = 0.
       do i = 1,npt1
          ti_sngl(i) = b1_sngl + (i-1)*dt1_sngl
       enddo

       !call newhdr()  ! create a new header

       ! set minmax values in sac file
       xmin_sngl = minval(dat_sngl)
       xmax_sngl = maxval(dat_sngl)
       call setfhv('depmin',xmin_sngl,nerr)
       call setfhv('depmax',xmax_sngl,nerr)

       call setnhv('npts',npt1,nerr)          ! sets number of points
       !call setfhv('b',ti_sngl(1),nerr)       ! sets begin
       !call setfhv('e',ti_sngl(npt1),nerr)    ! sets end
       !call setlhv('leven',.false.,nerr)        ! sets un-even sampling
       !call setihv('iftype','itime',nerr)          ! sets file type: time file

       ! write file with headers (LQY: bug with b in wsac0())
       ! call wsac0(datafile,ti_sngl,dat_sngl,nerr)
       call wsac1(datafile,dat_sngl,npt1,b1_sngl,dt1_sngl,nerr)
    else
       call wsac1(datafile,dat_sngl,npt1,b1_sngl,dt1_sngl,nerr)
    endif
    if (nerr > 0) then
        print *, 'Error writing sac file', trim(datafile)
        stop
    endif

  end subroutine dwsac1

  !-----------------------------------------------------------------------------

  subroutine cc_measure_select(tshift,dlnA,cc_max)

    ! CHT: If the CC timeshift is for some reason larger than the allowable max,
    !      then effectively eliminate the window by zeroing the
    !      cross-correlation traveltime and amplitude measurements.
    ! See subroutine compute_cc in mt_sub.f90.

    implicit none
    double precision, intent(inout) :: tshift, dlnA, cc_max

    if( (cc_max < CC_MIN) .or. (tshift < TSHIFT_MIN) .or. (tshift > TSHIFT_MAX) &
                          .or. (dlnA < DLNA_MIN) .or. (dlnA > DLNA_MAX) ) then
       ! zero the CC measurments
       if (DISPLAY_DETAILS) then
          print *, 'CC measurements: failed because ONE of these is true :'
          print *, ' cc_max      : ', cc_max, CC_MIN, cc_max < CC_MIN
          print *, ' tshift      : ', tshift, TSHIFT_MIN, tshift < TSHIFT_MIN
          print *, ' tshift      : ', tshift, TSHIFT_MAX, tshift > TSHIFT_MAX
          print *, ' dlnA        : ', dlnA, DLNA_MIN, dlnA < DLNA_MIN
          print *, ' dlnA        : ', dlnA, DLNA_MAX, dlnA > DLNA_MAX
       endif

       ! zero the CC measurments
       tshift = 0.0
       dlnA = 0.0
    endif

  end subroutine cc_measure_select

  !-----------------------------------------------------------------------------

  subroutine mt_measure_select(nlen,tshift,i_pmax_syn,dtau_w,err_dt, &
                                dt,i_left,i_right,fstart,fend,use_trace)

    ! an important subroutine to determine whether an MT measurement should be rejected,
    ! in which case a CC measurement is used -- several choices are made here

    implicit none
    integer, intent(in) :: nlen, i_pmax_syn
    double precision, intent(in) :: tshift, dt
    double precision, dimension(:), intent(inout) :: dtau_w, err_dt
    double precision, intent(inout) :: fstart, fend
    integer,intent(inout) :: i_left, i_right
    logical,intent(out) :: use_trace

    double precision :: df, fvec(NPT), f_pmax, T_pmax, Wlen
    integer :: i_right_old, i_left_old
    integer :: j,ntaper
    !logical :: stop_freq

    use_trace = .true.
    df = 1./(dt*NPT)
    f_pmax = df * i_pmax_syn
    T_pmax = 1./ f_pmax
    Wlen = dt*nlen

    if( NCYCLE_IN_WINDOW * T_pmax > Wlen ) then
       print *, '   MTM: rejecting for too few cycles within time window:'
       print *, '   T_pmax : ', sngl(T_pmax)
       print *, '   Wlen : ', sngl(Wlen)
       print *, '   NCYCLE_IN_WINDOW : ', NCYCLE_IN_WINDOW
       print *, '   REJECTION: ', sngl(NCYCLE_IN_WINDOW*T_pmax), &
            sngl(Wlen), NCYCLE_IN_WINDOW * T_pmax < Wlen
       use_trace = .false.
    endif

    !write(*,'(a8,4f12.6)') 'fstart :', fstart, NCYCLE_IN_WINDOW/(Wlen), NCYCLE_IN_WINDOW, Wlen
    !write(*,'(a8,4f12.6)') 'fend :', fend, 1./(2.0*dt), dt

    ! DECREASE the frequency range of the measurement (and adjoint source)
    ! --> note NCYCLE_IN_WINDOW and window length
    ! We subjectively state that we want at least ntaper number of
    ! frequency points between [fstart, fend] for the multitaper measurement.
    fstart = max(fstart, NCYCLE_IN_WINDOW/Wlen)
    fend = min(fend, 1./(2.0*dt))

    ! number of tapers (slepian tapers, type = 1)
    ntaper = int(NPI * 2.0)
    if( ntaper > 10 ) ntaper = 10
    if( ntaper < 1 ) ntaper = 10
    if( use_trace .and. fstart >= fend - ntaper*df ) then
       print *, '   MTM: rejecting for frequency range (NCYCLE_IN_WINDOW/Wlen):'
       print *, '     fstart, fend, df, ntaper : ', sngl(fstart),sngl(fend),sngl(df),ntaper
       print *, '     NCYCLE_IN_WINDOW, Wlen : ', NCYCLE_IN_WINDOW,sngl(Wlen), &
            sngl(NCYCLE_IN_WINDOW/Wlen)
       print *, '     REJECTION fstart >= fend - ntaper*df : ', sngl(fstart), &
            sngl(fend - ntaper*df), fstart >= fend - ntaper*df
       use_trace = .false.
    endif

    ! assemble frequency vector (NPT is the FFT length)
    fvec(:) = 0.
    do j = 1,NPT
      if(j > NPT/2+1) then
        fvec(j) = df*(j-NPT-1)   ! negative frequencies in second half
      else
        fvec(j) = df*(j-1)       ! positive frequencies in first half
      endif
    enddo

!!$    stop_freq = .false.
!!$    do j = 1, i_right
!!$      if (stop_freq) exit
!!$      print *, j, dtau_w(j),stop_freq
!!$      if (abs(dtau_w(j)) > 3 * abs(tshift)) then
!!$        dtau_w(j) = 0
!!$      else if (j /= 1) then
!!$        stop_freq = .true.
!!$      endif
!!$    enddo

    ! determine the indices that denote the new frequency range (CHT)
    ! IT SEEMS LIKE THERE SHOULD BE NO NEED FOR THIS, SINCE THE SIGNAL HAS ALREADY
    ! BEEN BAND-PASSED PRIOR TO MAKING THE MULTITAPER MEASUREMENT.
    ! LQY: this should be useful to constraint the range of which use_trace is checked

    i_left_old = i_left
    i_right_old = i_right
    do j = i_left_old, i_right_old
       if (fvec(j) > fstart) then
          i_left = j-1
          exit
       endif
    enddo
    do j = i_left_old, i_right_old
       if (fvec(j) > fend) then
          i_right = j-1
          exit
       endif
    enddo

    if (DISPLAY_DETAILS) then
       print *, '   determine the validity of MTM over new f bounds:'
       write(*,'(a24,2i6,2f14.8)') '      Old frequency bounds :', i_left_old, i_right_old, &
            sngl(df*(i_left_old-1)), sngl(df*(i_right_old-1))
       write(*,'(a24,2i6,2f14.8)') '      New frequency bounds :', i_left, i_right, &
            sngl(df*(i_left-1)), sngl(df*(i_right-1))
    endif

    ! update the frequency limits
    fstart = (i_left-1)*df
    fend = (i_right-1)*df

    ! if the cross-correlation time-shift is <= a time-step, set dtau(w) to zero
    ! NOTE: this should probably be a user parameter
    if ( abs(tshift) <= 1.01*dt ) then
       dtau_w(:) = 0.
       use_trace = .false.
       if (DISPLAY_DETAILS) then
          print *, 'MTM: rejecting for too small a time shift:'
          print *, '         dt = ', sngl(dt)
          print *, '  tshift = ', sngl(tshift)
       endif
    endif

    ! within the frequency range of interest, check various criteria
    ! CHT: dtau_w(j) --> abs(dtau_w(j)) for the first criterion
    do j = i_left, i_right
       if (use_trace .and. (abs(dtau_w(j)) > 1./(DT_FAC*fvec(j)) .or. err_dt(j) > 1./(ERR_FAC*fvec(j)) &
            .or. abs(dtau_w(j)) > DT_MAX_SCALE*abs(tshift))) then
          use_trace = .false.
          if (DISPLAY_DETAILS) then
             print *, '     MTM: rejecting trace on dtau/err at specific frequency:'
             print *, '       f = ', sngl(fvec(j)), j, sngl(1/fvec(j))
             print *, '       DT_FAC (lower) : ', sngl(abs(dtau_w(j))), &
                  sngl(1./(DT_FAC * fvec(j))), abs(dtau_w(j)) > 1./(DT_FAC * fvec(j))
             print *, '       ERR_FAC (lower) : ', sngl(err_dt(j)), &
                  sngl(1./(ERR_FAC * fvec(j))), err_dt(j) > 1./(ERR_FAC * fvec(j))
             print *, '       DT_MAX_SCALE (lower) : ', sngl(abs(dtau_w(j))), &
                  sngl(DT_MAX_SCALE*abs(tshift)), abs(dtau_w(j)) > DT_MAX_SCALE*abs(tshift)
          endif
       endif
    enddo

  end subroutine mt_measure_select

  !==============================================================================
  !        subroutines used in mtm_measure() and mtm_adj()
  !==============================================================================

  subroutine interpolate_dat_and_syn(data, syn, syn_phydisp, tstart, tend, t0, dt, NPT, &
       dat_win, syn_win, syn_win_phydisp, nlen, istart)

    ! extract windowed portion from original data/syn traces

    implicit none

    double precision, dimension(NPT), intent(in) :: data, syn, syn_phydisp
    double precision, intent(in) :: tstart
    double precision, intent(in) ::  tend, t0, dt
    integer, intent(in) :: NPT

    double precision, dimension(NPT), intent(out) :: dat_win, syn_win, syn_win_phydisp
    integer, intent(out) :: nlen, istart

    integer :: ii, i
    double precision :: time, t1

    nlen = floor((tend-tstart)/dt) + 1

    istart = floor((tstart-t0)/dt) + 1
    ! tstart = t0+(istart-1)*dt ! minor adjustments
    !print *, '*** diff tstart = ', t0+(istart-1)*dt - tstart

    ! limits array bounds
    if( nlen > NPT ) nlen = NPT

    ! move checking inside subroutine
    if (nlen <= 1) stop 'Check the length of the data and syn arrays'

    do i = 1, nlen
      time = tstart + (i-1) * dt
      ii = floor((time-t0)/dt) + 1

      ! checks out-of-bounds (very unlikely event!)
      if( ii >= NPT ) cycle

      t1 = floor((time-t0)/dt) * dt + t0

      dat_win(i) = data(ii) + (data(ii+1)-data(ii)) * (time-t1) / dt
      syn_win(i) = syn(ii) + (syn(ii+1)-syn(ii)) * (time-t1) /dt
      if (USE_PHYSICAL_DISPERSION) then
        syn_win_phydisp(i) = syn_phydisp(ii) + ( syn_phydisp(ii+1) - syn_phydisp(ii) ) * (time-t1)/dt
      endif
    enddo

  end subroutine interpolate_dat_and_syn

  !-----------------------------------------------------------------------------

  subroutine compute_cc(syn, data, nlen, dt, ishift, tshift, dlnA, cc_max)

    ! time shift MEASUREMENT between data (data) and synthetics (syn)
    ! CHT: modified the subroutine to resemble the one used in FLEXWIN

    implicit none
    double precision, dimension(*), intent(in) :: syn, data
    integer, intent(in) :: nlen
    double precision, intent(in) :: dt
    double precision, intent(out) :: tshift, dlnA, cc_max
    integer, intent(out) :: ishift

    double precision :: cc, norm_s, norm ! cr_shift
    integer i1, i2, i, j, i_left, i_right, id_left, id_right

!!$    ! these choices will slide the entire windowed record past the other
!!$    cr_shift = nlen*dt
!!$    i_left  = ceiling( -1.0 * cr_shift / dt )
!!$    i_right = floor( cr_shift / dt )
!!$
!!$    ! cross-correlation
!!$    ishift = 0
!!$    do i = i_left, i_right, 1
!!$
!!$      cc = 0.
!!$      do j = 1, nlen
!!$        if((j+i) > 1 .and. (j+i) < nlen) cc = cc + syn(j) * data(j+i)
!!$      enddo
!!$
!!$      !if(cc > cc_max) then
!!$      ! CHT, 07-Sept-2008: Do not allow time shifts larger than the specified input
!!$      if(cc > cc_max .and. abs(i*dt) <= BEFORE_TSHIFT ) then
!!$        cc_max = cc
!!$        ishift = i
!!$      endif
!!$
!!$    enddo
!!$    tshift = ishift*dt

    ! initialise shift and cross correlation to zero
    ishift = 0
    cc_max = 0.0

    ! index of window limits
    i1 = 1
    i2 = nlen

    ! length of window (number of points, including ends)
    !nlen = i2 - i1 + 1

    ! power of synthetic signal in window
    norm_s = sqrt(sum(syn(i1:i2)*syn(i1:i2)))

    ! left and right limits of index (time) shift search
    ! NOTE: This looks OUTSIDE the time window of interest to compute TSHIFT and CC.
    !       How far to look outside, in theory, should be another parameter.
    !       However, it does not matter as much if the data and synthetics are
    !          zeroed outside the windows.
    i_left = -1*int(nlen/2.0)
    i_right = int(nlen/2.0)

    ! i is the index to shift to be applied to DATA (data)
    do i = i_left, i_right

       ! normalization factor varies as you take different windows of data
       id_left = max(1,i1+i)      ! left index for data window
       id_right = min(nlen,i2+i)  ! right index for data window
       norm = norm_s * sqrt(sum(data(id_left:id_right)*(data(id_left:id_right))))

       ! cc as a function of i
       cc = 0.
       do j = i1, i2   ! loop over full window length
          if((j+i)>=1 .and. (j+i)<=nlen) cc = cc + syn(j)*data(j+i)  ! d is shifted by i
       enddo
       cc = cc/norm

       if (cc > cc_max) then
          ! CHT: do not allow time shifts larger than the specified input range
          ! This is an important criterion, since it may pick TSHIFT_MIN or TSHIFT_MAX
          ! if cc_max within the interval occurs on the boundary.
          if( (i*dt >= TSHIFT_MIN).and.(i*dt <= TSHIFT_MAX) ) then
             cc_max = cc
             ishift = i
          endif
       endif

    enddo
    tshift = ishift*dt

    ! The previously used expression for dlnA of Dahlen and Baig (2002),
    ! is a first-order perturbation of ln(A1/A2) = (A1-A2)/A2 .
    ! The new expression is better suited to getting Gaussian-distributed
    ! values between -1 and 1, with dlnA = 0 indicating perfect fit, as before.
    dlnA = 0.5 * log( sum(data(i1:i2) * data(i1:i2)) / sum(syn(i1:i2) * syn(i1:i2)) )

  end subroutine compute_cc

  ! ---------------------------------------------------------------------------

  subroutine compute_average_error(data_dtw,syn_dtw_cc,syn_dtw_cc_dt,nlen,dt,sigma_dt,sigma_dlnA)

  ! CHT: Estimate the uncertainty in the CC measurement
  !      based on the integrated waveform difference between the data
  !      and the reconstructed synthetics.
  ! NOTE: We implement the exact equations that are in the Latex notes.

    implicit none
    double precision, dimension(*), intent(in) :: data_dtw, syn_dtw_cc, syn_dtw_cc_dt
    integer, intent(in) :: nlen
    double precision, intent(in) :: dt
    double precision, intent(inout) :: sigma_dt, sigma_dlnA

    double precision, dimension(nlen) :: syn_vtw_cc
    double precision :: sigma_dt_top, sigma_dlnA_top, sigma_dt_bot, sigma_dlnA_bot
    integer i

    ! compute synthetic velocity (shifted and stretched)
    do i = 2, nlen-1
      syn_vtw_cc(i) = (syn_dtw_cc(i+1) - syn_dtw_cc(i-1)) / (2.0*dt)
    enddo
    syn_vtw_cc(1)    = (syn_dtw_cc(2) - syn_dtw_cc(1)) / dt
    syn_vtw_cc(nlen) = (syn_dtw_cc(nlen) - syn_dtw_cc(nlen-1)) / dt

    ! estimated uncertainty in cross-correlation traveltime and amplitude
    sigma_dt_top   = sum( (data_dtw(1:nlen) - syn_dtw_cc(1:nlen) )**2 )
    sigma_dlnA_top = sigma_dt_top
    sigma_dt_bot   = sum( syn_vtw_cc(1:nlen)**2 )
    sigma_dlnA_bot = sum( (syn_dtw_cc_dt(1:nlen))**2 )
    sigma_dt       = sqrt( sigma_dt_top / sigma_dt_bot )
    sigma_dlnA     = sqrt( sigma_dlnA_top / sigma_dlnA_bot )

    ! testing
!!$    print *, ' sigma_dt   : ', sigma_dt
!!$    print *, ' sigma_dlnA : ', sigma_dlnA
!!$    open(88,file='tshift.dat')
!!$    do i = 1,nlen
!!$       write(88,'(5e18.6)') (i-1)*dt, data_dtw(i), syn_dtw_cc(i), syn_dtw_cc_dt(i), syn_vtw_cc(i)
!!$    enddo
!!$    close(88)

    ! make final adjustments to uncertainty estimate
    if (ERROR_TYPE == 0) then
       ! set uncertainty factors to 1 if you do not want to incorporate them
       ! into the adjoint sources and the misfit function values
       sigma_dt = 1.0
       sigma_dlnA = 1.0

    else
       ! make sure that the uncertainty estimates are not below the water level;
       ! otherwise, the adjoint sources will blow up unreasonably
       if( sigma_dt < DT_SIGMA_MIN) sigma_dt = DT_SIGMA_MIN
       if( sigma_dlnA < DLNA_SIGMA_MIN) sigma_dlnA = DLNA_SIGMA_MIN

    endif

  end subroutine compute_average_error

  ! ---------------------------------------------------------------------------

  subroutine write_average_meas(filename, is_mtm, dtau_meas, dlnA_meas, dtau_sigma, dlnA_sigma)

    implicit none
    character(len=*), intent(in) :: filename
    double precision, intent(in) :: dtau_meas, dlnA_meas, dtau_sigma, dlnA_sigma
    integer, intent(in) :: is_mtm

    character(len=40) :: stlab, suffix

    if ( is_mtm == 1 ) then
       stlab = 'multitaper averaged' ; suffix = 'average'
    else
       stlab = 'cross-correlation' ; suffix = 'cc'
    endif

    if ( is_mtm > 0) then ! CC or MT
       if (DISPLAY_DETAILS) then
          print *, '   write ', trim(stlab)//'  measurements:'
          print *, '     traveltime :', sngl(dtau_meas), ' +/- ', sngl(dtau_sigma)
          print *, '     amplitude  :', sngl(dlnA_meas), ' +/- ', sngl(dlnA_sigma)
       endif

       ! write average error estimates to file
       if (OUTPUT_MEASUREMENT_FILES) then
          open(71,file=trim(filename)//'.dt_'//trim(suffix))
          write(71,*) dtau_meas, dtau_sigma
          close(71)
          open(72,file=trim(filename)//'.dlnA_'//trim(suffix))
          write(72,*) dlnA_meas, dlnA_sigma
          close(72)
       endif
    endif

  end subroutine write_average_meas

  ! ---------------------------------------------------------------------------

  subroutine write_trans(filename, trans, wvec, fnum, i_right, idf_new, df, tshift, dlnA, &
       phi_wt, abs_wt, dtau_wt, dlnA_wt, dtau_wa, dlnA_wa)

    ! The transfer function maps the synthetics to the CC-deconstructed data;
    ! the CC measurements then need to be applied to match the original data.
    ! dummy: fnum

    implicit none
    character(len=*), intent(in) :: filename
    complex*16, intent(in) :: trans(:)
    double precision, intent(in) :: wvec(:), df, tshift, dlnA
    integer, intent(in) :: fnum, i_right, idf_new

    double precision, dimension(:), intent(out) :: phi_wt, abs_wt, dtau_wt, dlnA_wt
    double precision, intent(out), optional :: dtau_wa, dlnA_wa

    integer :: i, j
    double precision, dimension(NPT) :: fr
    double precision :: smth, smth1, smth2 ! f0
    logical :: ioactive

    ioactive = .false.  ! IO off
    if (present(dtau_wa) .and. present(dlnA_wa)) ioactive = .true.

    abs_wt(:) = 0.
    phi_wt(:) = 0.

    ! note that with the idf_new value, these files are SUB-SAMPLED
    if (OUTPUT_MEASUREMENT_FILES .and. ioactive) then
      open(10,file=trim(filename)//'.ph')
      open(20,file=trim(filename)//'.abs')
      open(30,file=trim(filename)//'.dlnA')
      open(40,file=trim(filename)//'.ph_cor')
      open(50,file=trim(filename)//'.dt')
    endif

    ! loop to calculate phase and amplitude
    do i = 1, i_right
      phi_wt(i) = atan2( aimag(trans(i)) , real(trans(i)) )
      abs_wt(i) = abs(trans(i))
      fr(i) = df*(i-1)
      if (mod(i,idf_new)==0 .and. OUTPUT_MEASUREMENT_FILES .and. ioactive) then
        write(10,*) fr(i), phi_wt(i)
        write(20,*) fr(i), abs_wt(i)
! LQY: should not we apply dlnA from cc before this?
!         write(30,*) fr(i), log(abs_wt(i))
      endif
    enddo

    ! NOTE: the CC measurements dT (tshift) and dlnA are BOTH included
    dtau_wt(1) = tshift
    do i = 1, i_right

      if (i > 1 .and. i < i_right) then
        ! check the smoothness (2nd-order derivative) by 2*pi changes
        smth  =  phi_wt(i+1) + phi_wt(i-1) - 2.0 * phi_wt(i)
        smth1 = (phi_wt(i+1) + TWOPI) + phi_wt(i-1) - 2.0 * phi_wt(i)
        smth2 = (phi_wt(i+1) - TWOPI) + phi_wt(i-1) - 2.0 * phi_wt(i)
        if(abs(smth1)<abs(smth).and.abs(smth1)<abs(smth2).and. abs(phi_wt(i) - phi_wt(i+1)) > PHASE_STEP)then
          if (DISPLAY_DETAILS .and. ioactive) &
               print *, '     phase correction : 2 pi', sngl(fr(i)), sngl(phi_wt(i) - phi_wt(i+1))
          do j = i+1, i_right
            phi_wt(j) = phi_wt(j) + TWOPI
          enddo
        endif
        if(abs(smth2)<abs(smth).and.abs(smth2)<abs(smth1).and. abs(phi_wt(i) - phi_wt(i+1)) > PHASE_STEP)then
          if (DISPLAY_DETAILS .and. ioactive) &
               print *, '     phase correction : - 2 pi', sngl(fr(i)), sngl(phi_wt(i) - phi_wt(i+1))
          do j = i+1, i_right
            phi_wt(j) = phi_wt(j) - TWOPI
          enddo
        endif
      endif

      ! add the CC measurements to the transfer function
      if (i > 1) dtau_wt(i) = (-1./wvec(i)) * phi_wt(i) + tshift
      dlnA_wt(i) = log(abs_wt(i)) + dlnA
      !!dlnA_wt(i) = abs_wt(i) - 1. + dlnA

      if(mod(i,idf_new)==0 .and. OUTPUT_MEASUREMENT_FILES .and. ioactive) then
        write(30,*) fr(i), dlnA_wt(i)
        write(40,*) fr(i), phi_wt(i)
        write(50,*) fr(i), dtau_wt(i)
      endif

    enddo

    if (OUTPUT_MEASUREMENT_FILES .and. ioactive) then
      close(10)
      close(20)
      close(30)
      close(40)
      close(50)
    endif

    ! average values of the transfer functions (optional output argument)
    if (present(dtau_wa) .and. present(dlnA_wa)) then
       dtau_wa = sum( dtau_wt(1:i_right) ) / i_right
       dlnA_wa = sum( dlnA_wt(1:i_right) ) / i_right
    endif

!!$    if (DISPLAY_DETAILS) then
!!$      print *, ' Taper traveltime measurement average : ', sngl(dtau_wa)
!!$      print *, ' Taper amplitude measurement average : ', sngl(dlnA_wa)
!!$      print *, ' i_right : ', i_right
!!$      !f0 = 0.
!!$      !call dwrite_ascfile_f(trim(filename)//'.dt_full',f0,df,i_right,dtau_wt(1:i_right))
!!$      !call dwrite_ascfile_f(trim(filename)//'.dlnA_full',f0,df,i_right,dlnA_wt(1:i_right))
!!$      !call dwrite_ascfile_f(trim(filename)//'.transfer_full',f0,df,i_right,abs(trans(1:i_right)))
!!$    endif

  end subroutine write_trans

  ! --------------------------------------------------------------------

  subroutine deconstruct_dat_cc(filename,dat_dtw,tstart,dt,nlen,&
       ishift,tshift,dlnA,dat_dtw_cc)

    ! Using CC measurements, map the data to the synthetics;
    ! because the windows are picked based on the synthetics,
    ! we apply the transfer function from the synthetics to the
    ! CC-deconstructed data.
    ! dummy inputs: tshift, filename, dt

    implicit none
    character(len=*), intent(in) :: filename
    double precision, dimension(NPT), intent(in) :: dat_dtw
    integer, intent(in) :: ishift, nlen
    double precision, intent(in) :: tshift, dlnA, tstart, dt
    double precision, dimension(NPT), intent(out) :: dat_dtw_cc
    integer i

    ! apply time shift (-dT) to OBSERVED seismogram
    dat_dtw_cc(:) = dat_dtw(:)
    do i = 1, nlen
      if ((ishift+i) > 1 .and. (ishift+i) < nlen ) dat_dtw_cc(i) = dat_dtw(i+ishift)
    enddo
    ! fill the missing time window with the endpoint value
    if (ishift < 0) dat_dtw_cc(1:-ishift+1) = dat_dtw_cc(-ishift+2)
    if (ishift > 0) dat_dtw_cc(nlen-ishift:nlen) = dat_dtw_cc(nlen-ishift-1)

    ! apply cross-correlation amplitude measurement (-DlnA) to OBSERVED seismogram
    dat_dtw_cc(:) = dat_dtw_cc(:) * exp( -dlnA )

    !if (DISPLAY_DETAILS) then
    !   call dwrite_sacfile_f(datafile,'windowed_shifted_data.sac',tstart,nlen,dat_dtw0)
    !endif

  end subroutine deconstruct_dat_cc

  ! --------------------------------------------------------------------

  subroutine reconstruct_syn_cc(syn_dtw,tstart,dt,nlen,ishift,tshift,dlnA,syn_dtw_cc,syn_dtw_cc_dt)

    ! reconstruct the synthetics using cross-correlation measurements:
    !    (1) apply dT to get syn_dtw_cc_dt
    !    (2) apply dT and dlnA to get syn_dtw_cc
    ! dt, tshift are dummy
    implicit none
    double precision, dimension(NPT), intent(in) :: syn_dtw
    integer, intent(in) :: ishift, nlen
    double precision, intent(in) :: tshift, dlnA, tstart, dt
    double precision, dimension(NPT), intent(out) :: syn_dtw_cc, syn_dtw_cc_dt
    integer i

    ! shift synthetics by tshift (in the main program, we shift the data instead)
    ! ishift = tshift * dt
    syn_dtw_cc_dt(:) = syn_dtw(:)
    do i = 1, nlen
      if ((i-ishift) > 1 .and. (i-ishift) < nlen ) syn_dtw_cc_dt(i) = syn_dtw(i-ishift)
    enddo
    ! fill the missing time window with the endpoint value
    if (ishift > 0) syn_dtw_cc_dt(1:ishift+1) = syn_dtw_cc_dt(ishift+2)
    if (ishift < 0) syn_dtw_cc_dt(nlen+ishift:nlen) = syn_dtw_cc_dt(nlen+ishift-1)

    ! apply cross-correlation amplitude measurement
    syn_dtw_cc(:) = 0.
    syn_dtw_cc(:) = syn_dtw_cc_dt * exp( dlnA )    ! based on dlnA = ln(Aobs/Asyn)
    !syn_dtw_cc(:) = syn_dtw_cc_dt * (1. + dlnA)   ! based on first-order approximation of dlnA

  end subroutine reconstruct_syn_cc

  ! -----------------------------------------------------------------------

  subroutine reconstruct_syn(filename, syn_dtwo, wvec, dtau_wt, dlnA_wt, &
       i_right, tstart, dt, nlen, syn_dtw_mt, syn_dtw_mt_dt)

    ! reconstruct the synthetics using multitaper measurements:
    !    (1) apply dtau(w) and dlnA(w) to get syn_dtw_mt0
    implicit none
    character(len=*), intent(in) :: filename
    complex*16, dimension(:), intent(in) ::  syn_dtwo
    double precision, dimension(:), intent(in) :: wvec, dtau_wt, dlnA_wt
    integer, intent(in) :: i_right, nlen
    double precision, intent(in) :: tstart, dt

    double precision, dimension(:), intent(out) :: syn_dtw_mt, syn_dtw_mt_dt

    complex*16, dimension(NPT) :: wseis_recon
    integer i
    double precision omega

    ! apply transfer function to synthetics (phase and amplitude)
    syn_dtw_mt(:) = 0.
    wseis_recon(:) = cmplx(0.,0.)
    do i = 1,i_right
      omega = wvec(i)
      wseis_recon(i) = syn_dtwo(i) * exp(dlnA_wt(i)) * exp(-CCI*omega*dtau_wt(i))
      !wseis_recon(i) = syn_dtwo(i) * (1.+ dlnA_wt(i)) * exp(-CCI*omega*dtau_wt(i))
      !wseis_recon(i) = syn_dtwo(i) * trans_mtm(i) * exp(-CCI*omega*tshift)
    enddo
    call fftinv(LNPT,wseis_recon,REVERSE_FFT,dt,syn_dtw_mt)

    ! apply phase shifts only
    syn_dtw_mt_dt(:) = 0.
    wseis_recon(:) = cmplx(0.,0.)
    do i = 1,i_right
      omega = wvec(i)
      wseis_recon(i) = syn_dtwo(i) * exp(-CCI*omega*dtau_wt(i))
    enddo
    call fftinv(LNPT,wseis_recon,REVERSE_FFT,dt,syn_dtw_mt_dt)

    if (OUTPUT_MEASUREMENT_FILES) then
       call dwsac1(trim(filename)//'.recon_syn.sac',syn_dtw_mt,nlen,tstart,dt)
       call dwsac1(trim(filename)//'.recon_syn_dt.sac',syn_dtw_mt_dt,nlen,tstart,dt)
    endif

  end subroutine reconstruct_syn

  ! -----------------------------------------------------------------------

!!$  subroutine check_recon_quality(filename,dat_dtw_cc,syn_dtw,dat_dtw,syn_dtw_mt,nlen,dt,tshift,tshift_f1f2,cc_max_f1f2,cc_max)
!!$
!!$    character(len=*), intent(in) :: filename
!!$    double precision, dimension(:), intent(in) :: dat_dtw_cc, syn_dtw, dat_dtw, syn_dtw_mt
!!$    double precision, intent(in) :: dt, tshift
!!$    integer, intent(in) :: nlen
!!$    double precision, intent(out) :: tshift_f1f2, cc_max_f1f2, cc_max
!!$
!!$    double precision :: f1,f2, dlnA_f1f2
!!$
!!$    ! Using Alessia's subroutine
!!$    !     First the shifted_obs_win vs the synthetic
!!$    call f1f2_calc(dat_dtw_cc,syn_dtw,nlen,1,nlen,dt, f1,f2,tshift_f1f2,cc_max_f1f2,dlnA_f1f2)
!!$
!!$    cc_max = cc_max_f1f2
!!$    if (OUTPUT_MEASUREMENT_FILES) then
!!$      open(10,file=trim(filename)//'.quality')
!!$      write(10,*) '<--------- F1 ------ F2 ---- tshift -- cc_max --- dlnA -->'
!!$      write(10,"(a,5F10.5)") 'Before',sngl(F1),sngl(F2),sngl(tshift),sngl(cc_max_f1f2),sngl(dlnA_f1f2)
!!$    endif
!!$    if (DISPLAY_DETAILS) then
!!$      write(*,*) '<--------- F1 ------ F2 ---- tshift -- cc_max --- dlnA -->'
!!$      write(*,"(a,5F10.5)") 'Before',sngl(F1),sngl(F2),sngl(tshift),sngl(cc_max_f1f2),sngl(dlnA_f1f2)
!!$    endif
!!$
!!$    !     Then the obs_win vs the reconstructed obs
!!$    call f1f2_calc(dat_dtw,syn_dtw_mt,nlen,1,nlen,dt, f1,f2,tshift_f1f2,cc_max_f1f2,dlnA_f1f2)
!!$
!!$    if (OUTPUT_MEASUREMENT_FILES) then
!!$      write(10,"(a,5F10.5)") 'After ',sngl(F1),sngl(F2),sngl(tshift_f1f2),sngl(cc_max_f1f2),sngl(dlnA_f1f2)
!!$      close(10)
!!$    endif
!!$
!!$    if (DISPLAY_DETAILS) write(*,"(a,5F10.5)") 'After ',sngl(F1),sngl(F2),sngl(tshift_f1f2),sngl(cc_max_f1f2),sngl(dlnA_f1f2)
!!$
!!$  end subroutine check_recon_quality

!-------------------------------------------------------------------

   subroutine interpolate_syn(syn,t1,dt1,npt1,t2,dt2,npt2)

     implicit none
     double precision, dimension(:),intent(inout) :: syn
     integer,intent(in) :: npt1,npt2
     double precision,intent(in) :: t1,dt1,t2,dt2

     double precision :: syn1(NDIM), time, tt
     integer i, ii

     ! initializes trace holding interpolated values
     syn1(1:npt2) = 0.

     ! loops over number of time steps in complete trace
     do i = 1, npt2

       ! sets time (in s) at this time step:
       ! t2 : start time of trace
       ! dt2: delta_t of a single time step
       time = t2 + (i-1) * dt2

       ! checks if time is within measurement window
       ! t1: start time of measurement window
       ! npt1: number of time steps in measurement window
       ! dt1: delta_t of a single time step in measurement window
       if (time > t1 .and. time < t1 + (npt1-1)*dt1) then

         ! sets index of time steps within this window: is 1 at the beginning of window
         ii = floor((time-t1)/dt1) + 1

         ! time increment within this single time step to match the exact value of time
         tt = time - ((ii-1)*dt1 + t1)

         ! interpolates value of trace for the exact time
         syn1(i) = (syn(ii+1)-syn(ii)) * tt/dt1 + syn(ii)
       endif
     enddo

     ! saves interpolated values to output trace
     syn(1:npt2) = syn1(1:npt2)
     ! LQY: zero out any thing beyond npts
     if (npt1 > npt2) syn(npt2+1:npt1)=0.

   end subroutine interpolate_syn

!-------------------------------------------------------------------

   subroutine taper_start(syn,npt,itmax)

     implicit none
     double precision, dimension(:),intent(inout) :: syn
     integer,intent(in) :: npt, itmax
     double precision :: Wt
     integer :: i !,imax

     if (2*itmax > npt) stop 'Check taper_start of adjoint source'

     !imax = maxloc(abs(syn),dim=1)   ! index of the max value
     !Wt = TWOPI / (2.0*(imax-1))    ! period of the taper

     Wt = TWOPI / (2.0*(itmax-1))    ! period of the taper

     if(DISPLAY_DETAILS) print *, 'tapering start of adjoint source from index 1 to index ', itmax

     ! apply a cosine taper from the start to the max value,
     ! such that the starting point is exactly zero
     do i = 1, itmax
        syn(i) = syn(i) * ( 0.5*(1 - cos(Wt*(i-1))) )
     enddo

   end subroutine taper_start

!-------------------------------------------------------------------


   subroutine read_par_file(fstart0,fend0,tt,dtt,nn,chan)

     implicit none
     double precision, intent(out) :: fstart0,fend0,tt,dtt
     integer, intent(out) :: nn
     character(len=10), intent(out) :: chan
     integer :: ios

     ! input file MEASUREMENT.PAR -- see write_par_file.pl for details

     OUT_DIR = 'OUTPUT_FILES'   ! default

     open(10,file='MEASUREMENT.PAR',status='old',iostat=ios)
     if ( ios /= 0) stop 'Error opening MEASUREMENT.PAR file'
     read(10,*) tt,dtt,nn
     read(10,*) imeas0
     read(10,*) chan
     read(10,*) TLONG, TSHORT
     read(10,*) RUN_BANDPASS
     read(10,*) DISPLAY_DETAILS
     read(10,*) OUTPUT_MEASUREMENT_FILES
     read(10,*) COMPUTE_ADJOINT_SOURCE
     read(10,*) TSHIFT_MIN, TSHIFT_MAX
     read(10,*) DLNA_MIN, DLNA_MAX
     read(10,*) CC_MIN
     read(10,*) ERROR_TYPE
     read(10,*) DT_SIGMA_MIN
     read(10,*) DLNA_SIGMA_MIN
     read(10,*) ITAPER
     read(10,*) WTR,NPI
     read(10,*) DT_FAC
     read(10,*) ERR_FAC
     read(10,*) DT_MAX_SCALE
     read(10,*) NCYCLE_IN_WINDOW
     read(10,*) USE_PHYSICAL_DISPERSION
     close(10)

     imeas = imeas0

     ! check the read-in values
     print *, 'INPUTS FROM MEASUREMENT.PAR :'
     print *, '  tt, dtt, nn : ',sngl(tt),sngl(dtt),nn
     print *, '  imeas : ',imeas
     print *, '  chan : ',chan
     print *, '  TLONG, TSHORT : ',sngl(TLONG), sngl(TSHORT)
     fstart0 = 1./TLONG ; fend0 = 1./TSHORT
     print *, '  fstart, fend : ', sngl(fstart0), sngl(fend0)
     print *, '  RUN_BANDPASS : ',RUN_BANDPASS
     print *, '  DISPLAY_DETAILS : ',DISPLAY_DETAILS
     print *, '  OUTPUT_MEASUREMENT_FILES : ',OUTPUT_MEASUREMENT_FILES
     print *, '  COMPUTE_ADJOINT_SOURCE : ',COMPUTE_ADJOINT_SOURCE
     print *, '  TSHIFT_MIN, TSHIFT_MAX : ',sngl(TSHIFT_MIN), sngl(TSHIFT_MAX)
     print *, '  DLNA_MIN, DLNA_MAX : ',sngl(DLNA_MIN), sngl(DLNA_MAX)
     print *, '  CC_MIN : ',sngl(CC_MIN)
     print *, '  ERROR_TYPE : ',ERROR_TYPE
     print *, '  DT_SIGMA_MIN : ',sngl(DT_SIGMA_MIN)
     print *, '  DLNA_SIGMA_MIN : ',sngl(DLNA_SIGMA_MIN)
     print *, '  ITAPER : ',ITAPER
     print *, '  WTR, NPI : ',sngl(WTR),NPI
     print *, '  DT_FAC : ',sngl(DT_FAC)
     print *, '  ERR_FAC : ',sngl(ERR_FAC)
     print *, '  DT_MAX_SCALE : ',sngl(DT_MAX_SCALE)
     print *, '  NCYCLE_IN_WINDOW : ',NCYCLE_IN_WINDOW
     !stop 'checking PAR file input'

    ! old format way..
    !  open(10,file='MEASUREMENT.PAR',status='old',iostat=ios)
    !  read(10,'(a)') out_dir
    !  read(10,*) is_mtm0
    !  read(10,*) wtr,npi
    !  read(10,*) iker0
    !  read(10,*) RUN_BANDPASS
    !  read(10,*) TLONG, TSHORT
    !  read(10,*) tt,dtt,nn
    !  read(10,*) DISPLAY_DETAILS
    !  read(10,*) OUTPUT_MEASUREMENT_FILES
    !  read(10,*) INCLUDE_ERROR
    !  read(10,*) DT_FAC
    !  read(10,*) ERR_FAC
    !  read(10,*) DT_MAX_SCALE
    !  read(10,*) NCYCLE_IN_WINDOW
    !  read(10,*) BEFORE_QUALITY, AFTER_QUALITY
    !  read(10,*) BEFORE_TSHIFT, AFTER_TSHIFT
    !  read(10,*) DT_SIGMA_MIN, DLNA_SIGMA_MIN
    !  close(10)
    !
    !  out_dir = adjustl(out_dir)
    !  iker = iker0
    !  is_mtm = is_mtm0
    !
    !  ! check the read-in values
    !  print *, 'INPUTS FROM MEASUREMENT.PAR :'
    !  print *, '  is_mtm : ',is_mtm
    !  print *, '  wtr, npi : ',wtr,npi
    !  print *, '  iker : ',iker
    !  print *, '  RUN_BANDPASS :',RUN_BANDPASS
    !  print *, '  TLONG, TSHORT : ',TLONG, TSHORT
    !  fstart0 = 1./TLONG ; fend0 = 1./TSHORT
    !  print *, '  fstart, fend :', fstart0, fend0
    !  print *, '  tt, dtt, nn : ',tt,dtt,nn
    !  print *, '  out_dir : ',trim(out_dir)
    !  print *, '  DISPLAY_DETAILS :',DISPLAY_DETAILS
    !  print *, '  OUTPUT_MEASUREMENT_FILES :',OUTPUT_MEASUREMENT_FILES
    !  print *, '  INCLUDE_ERROR :',INCLUDE_ERROR
    !  print *, '  DT_FAC :',DT_FAC
    !  print *, '  ERR_FAC :',ERR_FAC
    !  print *, '  DT_MAX_SCALE :',DT_MAX_SCALE
    !  print *, '  NCYCLE_IN_WINDOW :',NCYCLE_IN_WINDOW
    !  print *, '  BEFORE_QUALITY, AFTER_QUALITY :',BEFORE_QUALITY, AFTER_QUALITY
    !  print *, '  BEFORE_TSHIFT, AFTER_TSHIFT :',BEFORE_TSHIFT, AFTER_TSHIFT
    !  print *, '  DT_SIGMA_MIN, DLNA_SIGMA_MIN :',DT_SIGMA_MIN, DLNA_SIGMA_MIN
    !  !stop 'checking PAR file input'
    ! apply filter (this should EXACTLY match the filter used in the windowing code)
    !trbdndw = 0.3
    !a = 30.
    !iord = 4
    !passes = 2

     ! ray density
     if( DO_RAY_DENSITY_SOURCE ) ERROR_TYPE = 0

     ! assign additional parameters and stop for certain inconsistencies
     if (fstart0 >= fend0) &
          stop 'Check input frequency range of the signal'

     if (nn > NDIM) &
          stop 'Error: Change interpolation nn or NDIM'

     ! LQY: what happens if imeas = 7/8, and itaper = 2,3, is that permitted?
     ! LQY: how about imeas = 1/2, itaper = 1,2,3 matters?

     if ( imeas == 1 .or. imeas == 2 ) then ! waveforms
        is_mtm0 = 0
     else if ( imeas >= 3 .and. imeas <= 6 ) then ! CC Dt/DlnA
        ! for CC kernels, ITAPER must be a single taper (2 or 3)
        if ( ITAPER == 1 ) stop  'Error: Change ITAPER to 2/3 for CC measurements'
        is_mtm0 = ITAPER     ! 2 or 3 for CC adjoint sources
     else if ( imeas==7 .or. imeas==8 ) then
        is_mtm0 = 1          ! multitaper required for MT adjoint source
     else
        stop 'Error: imeas must by 1-8'
     endif

     is_mtm = is_mtm0
     print *, '  is_mtm :',is_mtm
     print *, ' '

   end subroutine read_par_file

!-------------------------------------------------------------------

  subroutine get_sacfile_header(data_file,yr,jda,ho,mi,sec,ntw,sta, &
                                comp,dist,az,baz,slat,slon)

    implicit none
    character(len=*),intent(in) :: data_file

    integer,intent(out):: yr,jda,ho,mi
    double precision,intent(out):: sec,dist,az,baz,slat,slon
    character(len=*),intent(out) :: ntw,sta,comp
    real :: tmp

    integer :: nsec,msec
    integer :: nerr

    ! note here data_file is a dummy argument, and we rely on the earlier
    ! call to rsac1() to retain header and waveform info in computer memory
    ! for data_file trace.


    !  integer header variables  (these values are not really used!)
    call getnhv('nzyear',yr,nerr)
    call getnhv('nzjday',jda,nerr)
    call getnhv('nzhour',ho,nerr)
    call getnhv('nzhour',mi,nerr)
    call getnhv('nzmin',nsec,nerr)
    call getnhv('nzmsec',msec,nerr)

    sec=nsec+msec/1000.0

    ! string headers
    call getkhv('knetwk',ntw,nerr)
    if(nerr /= 0) then
      write(*,*)'Error reading variable: knetwk'
      call exit(-1)
    endif

    call getkhv('kstnm',sta,nerr)
    if(nerr /= 0) then
      write(*,*)'Error reading variable: kstnm'
      call exit(-1)
    endif

    call getkhv('kcmpnm',comp,nerr)
    if(nerr /= 0) then
      write(*,*)'Error reading variable: kcmpnm'
      call exit(-1)
    endif

    ! decimal headers
    call getfhv('dist',tmp,nerr)
    if(nerr /= 0) then
      write(*,*)'Error reading variable: dist'
      call exit(-1)
    endif
    dist = tmp

    call getfhv('az',tmp,nerr)
    if(nerr /= 0) then
      write(*,*)'Error reading variable: az'
      call exit(-1)
    endif
    az = tmp

    call getfhv('baz',tmp,nerr)
    if(nerr /= 0) then
      write(*,*)'Error reading variable: baz'
      call exit(-1)
    endif
    baz = tmp

    call getfhv('stlo',tmp,nerr)
    if(nerr /= 0) then
      write(*,*)'Error reading variable: stlo'
      call exit(-1)
    endif
    slon = tmp

    call getfhv('stla',tmp,nerr)
    if(nerr /= 0) then
      write(*,*)'Error reading variable: stla'
      call exit(-1)
    endif
    slat = tmp

!!$    !  integer header variables
!!$    call saclst_iheader_f(data_file,'nzyear', yr)
!!$    call saclst_iheader_f(data_file,'nzjday', jda)
!!$    call saclst_iheader_f(data_file,'nzhour', ho)
!!$    call saclst_iheader_f(data_file,'nzmin',  mi)
!!$    call saclst_iheader_f(data_file,'nzsec',  nsec)
!!$    call saclst_iheader_f(data_file,'nzmsec', msec)
!!$
!!$    sec=nsec+msec/1000.0
!!$
!!$    call saclst_kheader_f(data_file,'knetwk',ntw,klen)
!!$    call saclst_kheader_f(data_file,'kstnm', sta,klen)
!!$    call saclst_kheader_f(data_file,'kcmpnm',comp,klen)
!!$
!!$    call dsaclst_fheader_f(data_file,'dist',dist)
!!$    call dsaclst_fheader_f(data_file,'az',  az)
!!$    call dsaclst_fheader_f(data_file,'baz', baz)
!!$    call dsaclst_fheader_f(data_file,'stlo',slon)
!!$    call dsaclst_fheader_f(data_file,'stla',slat)

  end subroutine get_sacfile_header

!!==================================================

subroutine setup_weighting(chan_syn)
  !
  ! determines weights based on number of window picks on Z/R/T components
  !
  use ma_weighting

  use ma_constants,only: NDIM
!  use ma_sub,only: get_sacfile_header,drsac1
  use ma_sub2,only: TOL

  implicit none
  character(len=10),intent(in) :: chan_syn

  ! local parameters
  integer :: npairs,ios,ipair,iposition,ipicks
  character(len=150) :: datafile,synfile !,dummy
  character(len=4) :: comp_T,comp_Z,comp_R
  integer :: picks_T, picks_Z, picks_R,npicks
  ! sac header information
  integer :: yr,jda,ho,mi
  double precision :: sec,dist,az,baz,slat,slon,T_surfacewaves
  character(len=10) :: net,sta,chan_dat,chan,cmp
  double precision :: t01, dt1, t02, dt2, t0, dt, tstart, tend
  integer :: npt1, npt2, npts
  double precision, dimension(NDIM) :: data, syn

  ! initializes
  picks_R = 0
  picks_Z = 0
  picks_T = 0

  num_P_SV_V = 0.d0
  num_P_SV_R = 0.d0
  num_SH_T = 0.d0

  num_Rayleigh_V = 0.d0
  num_Rayleigh_R = 0.d0
  num_Love_T = 0.d0

  ! substrings (synthetics components)
  comp_T = trim(chan_syn)//"T."
  comp_R = trim(chan_syn)//"R."
  comp_Z = trim(chan_syn)//"Z."

  ! opens measurement windows
  open(21,file='MEASUREMENT.WINDOWS',status='old',iostat=ios)
  if (ios /= 0) stop 'Error opening input file: MEASUREMENT WINDOWS'
  read(21,*,iostat=ios) npairs
  if (ios /= 0) stop 'Error reading number of pairs of data/syn'
  ! loops through windows
  do ipair=1,npairs

    ! reads in file names
    read(21,'(a)',iostat=ios) datafile
    if (ios /= 0) stop 'Error reading windows datafile'
    read(21,'(a)',iostat=ios) synfile
    if (ios /= 0) stop 'Error reading windows synfile'

    ! read data and syn (read datafile last to take its header later)
    call drsac1(synfile,syn,npt2,t02,dt2)
    call drsac1(datafile,data,npt1,t01,dt1)

    if (max(npt1,npt2) > NDIM) &
        stop 'Error: Too many npts in data or syn'

    ! check if t0 and dt match
    if (abs(dt1-dt2) > TOL) stop 'Error: check if dt match'
    dt = dt1
    npts = min(npt1,npt2)
    if (abs(t01-t02) > dt) then
      print *,'data t0: ',t01
      print *,'syn  t0: ',t02
      stop 'Check if t0 match'
    endif
    t0 = t01

    ! figure out station/network/comp names, etc
    call get_sacfile_header(trim(datafile),yr,jda,ho,mi,sec,net,sta, &
                            chan_dat,dist,az,baz,slat,slon)
    chan = chan_dat
    cmp = chan_dat(3:3)

    ! theoretical surface wave arrival time
    T_surfacewaves = dist / surface_vel

    ! debug output
    !if (DISPLAY_DETAILS) then
    !print *,'debug: '
    !print *,'  yr,jda,ho,mi,sec : ',yr,jda,ho,mi,sec
    !print *,'  net,sta,chan_dat : ',net,sta,chan_dat
    !print *,'  dist,az,baz,slat,slon : ',dist,az,baz,slat,slon
    !print *,'  cmp          = ',cmp
    !print *,'  dist           = ',dist
    !print *,'  T_surfacewaves = ',T_surfacewaves
    !print *
    !endif

    ! reads in window picks
    read(21,*,iostat=ios) npicks
    if (ios /= 0) stop 'Error reading windows npicks'

    ! loops/skips over picks (start/end times)
    do ipicks=1,npicks

      read(21,*,iostat=ios) tstart, tend
      if (ios /= 0) stop 'Error reading window pick: tstart and tend'

      tstart = max(tstart,t0)
      tend = min(tend, t0+(npts-1)*dt)

      ! body wave picks
      if( tend <= T_surfacewaves ) then
        if( cmp(1:1) == "Z" ) num_P_SV_V = num_P_SV_V + 1.d0
        if( cmp(1:1) == "R" ) num_P_SV_R = num_P_SV_R + 1.d0
        if( cmp(1:1) == "T" ) num_SH_T = num_SH_T + 1.d0
      else
      ! surface wave picks
        if( cmp(1:1) == "Z" ) num_Rayleigh_V = num_Rayleigh_V + 1.d0
        if( cmp(1:1) == "R" ) num_Rayleigh_R = num_Rayleigh_R + 1.d0
        if( cmp(1:1) == "T" ) num_Love_T = num_Love_T + 1.d0
      endif

    enddo

    ! determines all picks on a trace component
    ! (also cross-check comp name in filename)
    ! transverse
    iposition = INDEX( trim(synfile), comp_T, .false. )
    if( iposition > 3 .and. iposition < len_trim( synfile) ) then
      if( cmp(1:1) /= "T" ) stop 'error T component pick'
      picks_T = picks_T + npicks
    else
      ! radial
      iposition = INDEX( trim(synfile), comp_R, .false. )
      if( iposition > 3 .and. iposition < len_trim( synfile) ) then
        if( cmp(1:1) /= "R" ) stop 'error R component pick'
        picks_R = picks_R + npicks
      else
        ! vertical
        iposition = INDEX( trim(synfile), comp_Z, .false. )
        if( iposition > 3 .and. iposition < len_trim( synfile) ) then
          if( cmp(1:1) /= "Z" ) stop 'error Z component pick'
          picks_Z = picks_Z + npicks
        endif
      endif
    endif

  enddo
  close(21)


  ! check with total number of picks per component
  if( nint( num_P_SV_R + num_Rayleigh_R ) /= picks_R ) stop 'error R picks'
  if( nint( num_P_SV_V + num_Rayleigh_V ) /= picks_Z ) stop 'error Z picks'
  if( nint( num_SH_T + num_Love_T ) /= picks_T ) stop 'error T picks'

  if( DISPLAY_DETAILS ) then
    print *
    print *,'weighting measurements: '
    print *,'  picks T:',picks_T
    print *,'  picks R:',picks_R
    print *,'  picks Z:',picks_Z
    print *
    print *,'  picks P_SV_R: ',nint(num_P_SV_R)
    print *,'  picks P_SV_V: ',nint(num_P_SV_V)
    print *,'  picks SH_T  : ',nint(num_SH_T)
    print *,'  picks Rayleigh_R: ',nint(num_Rayleigh_R)
    print *,'  picks Rayleigh_V: ',nint(num_Rayleigh_V)
    print *,'  picks Love_T    : ',nint(num_Love_T)
    print *
  endif


  ! sets up weights based on picks
  weight_T = 1.0d0
  weight_R = 1.0d0
  weight_Z = 1.0d0

  ! weighting tries to balance love waves (tranverse) versus rayleigh waves (radial + vertical)
  !if( picks_T > 0 ) then
  !  if( picks_R + picks_Z > 0 ) weight_T = dble(picks_R + picks_Z)/dble(picks_T)
  !endif

  ! use normalization as weights
  if( picks_T > 0 ) weight_T = 1.d0 / picks_T
  if( picks_R > 0 ) weight_R = 1.d0 / picks_R
  if( picks_Z > 0 ) weight_Z = 1.d0 / picks_Z

  ! use normalization (no traces means zero weights)
  if( num_P_SV_R > 0. ) num_P_SV_R = 1.d0 / num_P_SV_R
  if( num_P_SV_V > 0. ) num_P_SV_V = 1.d0 / num_P_SV_V
  if( num_SH_T > 0. ) num_SH_T = 1.d0 / num_SH_T
  if( num_Rayleigh_R > 0. ) num_Rayleigh_R = 1.d0 / num_Rayleigh_R
  if( num_Rayleigh_V > 0. ) num_Rayleigh_V = 1.d0 / num_Rayleigh_V
  if( num_Love_T > 0. ) num_Love_T = 1.d0 / num_Love_T

  print *,'  weight of P_SV_R:',num_P_SV_R
  print *,'  weight of P_SV_V:',num_P_SV_V
  print *,'  weight of SH_T  :',num_SH_T
  print *,'  weight of Rayleigh_R:',num_Rayleigh_R
  print *,'  weight of Rayleigh_V:',num_Rayleigh_V
  print *,'  weight of Love_T    :',num_Love_T
  print *

end subroutine setup_weighting



end module ma_sub
