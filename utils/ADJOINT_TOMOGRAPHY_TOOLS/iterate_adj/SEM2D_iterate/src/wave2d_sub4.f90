module wave2d_sub4

  use wave2d_constants
  use wave2d_variables

  implicit none

  ! This module contains subroutines pertaining to making measurements
  ! between data and synthetics, and then constructing the corresponding adjoint source.
  !
  ! AS OF DEC 2006, THIS MODULE SHOULD BE REPLACED BY THE VERSION THAT QINYA HAS
  ! TESTED AND CHECKED INTO THE SVN SERVER (mt_measure_adj.f90).
  ! AS OF JAN 2008, THE MULTI-TAPER CAPABILITIES HAVE BEEN COMMENTED OUT,
  ! SINCE THE MEASUREMENT TOOLS ARE NO LONGER TESTED HERE.

contains

  !------------------------------------------

!!$  subroutine mtm_test(ievent, nrec, syn, tstart, tend, adj_syn, data, data_recon)
!!$
!!$    integer, intent(in) :: nrec, ievent
!!$    double precision, dimension(NSTEP,NCOMP,nrec),intent(in) :: syn
!!$    double precision, dimension(nrec), intent(in) :: tstart, tend
!!$    double precision, dimension(NSTEP,NCOMP,nrec),intent(out) :: adj_syn
!!$
!!$    double precision, dimension(NSTEP,NCOMP,nrec),intent(in),optional :: data
!!$    double precision, dimension(NSTEP,NCOMP,nrec),optional :: data_recon
!!$
!!$    complex*16, dimension(NDIM) :: wseis_syn, wseis_dat, wseis_dat2
!!$    double precision, dimension(NDIM) :: fp, fq, fq2, wvec
!!$    double precision, dimension(NSTEP) :: syn_displ, syn_veloc
!!$    double precision :: twopi, df, dw, w0, a0, f0
!!$
!!$    complex*16 :: cci
!!$    integer :: i,j,nlen,icomp,irec
!!$
!!$    !=======================================================
!!$
!!$    twopi = 2.*PI
!!$    df = 1./(npt*DT)
!!$    dw = twopi / (npt*DT)
!!$    cci = cmplx(0.,1.)
!!$    w0 = twopi/20.
!!$    f0 = 1./20.
!!$    a0 = 0.04
!!$
!!$    !icomp = 1
!!$
!!$    do icomp = 1,NCOMP
!!$
!!$       ! real time series for analysis
!!$       syn_displ(:) = 0.
!!$       syn_displ(1:NSTEP) = data(1:NSTEP,icomp,1)
!!$       syn_veloc(:) = 0.
!!$       syn_veloc(1:NSTEP) = data(1:NSTEP,icomp,2)
!!$
!!$  !open(29,file='junk1.dat',status='unknown')
!!$  !do i=1,NSTEP
!!$  !   write(29,*) i*DT, syn_displ(i), syn_veloc(i)
!!$  !enddo
!!$  !close(29)
!!$
!!$       nlen = NSTEP
!!$
!!$       ! FFT to obtain complex frequency-domain version
!!$       wseis_syn(:) = cmplx(0.,0.)
!!$       !do j = 1,nlen
!!$       !   wseis_syn(j)       = cmplx(syn_displ(j),0.)
!!$       !   wseis_syn(npt-j+1) = cmplx(syn_displ(j),0.)
!!$       !enddo
!!$       wseis_syn(1:nlen) = cmplx(syn_displ(1:nlen),0.)
!!$       call fft(lnpt,wseis_syn,FORWARD_FFT,DT)   ! displacement
!!$
!!$       ! KEY: assemble omega vector
!!$       wvec(:) = 0.
!!$       do j = 1,npt
!!$          if (j > npt/2) then
!!$             wvec(j) = dw*(j-npt-1)   ! negative frequencies in second half
!!$          else
!!$             wvec(j) = dw*(j-1)       ! positive frequencies in first half
!!$          endif
!!$       enddo
!!$
!!$       ! analytical expression (Harris and Stocker, 1998)
!!$       ! With this IFFT routine, only the positive frequencies need to be filled, in order
!!$       ! to recover the original REAL signal when we IFFT.
!!$       wseis_dat(:) = cmplx(0.,0.)
!!$       do j = 1,npt
!!$          !do j = 1,npt/2
!!$          wseis_dat(j)  = cmplx(w0,0.) / ( (cmplx(wvec(j),0.)*cci + cmplx(a0,0.))**2 + cmplx(w0**2,0.) )
!!$          wseis_dat2(j) = wseis_dat(j) * cmplx(0.,wvec(j))
!!$       enddo
!!$
!!$       open(91,file='test1.dat',status='unknown')
!!$       do j = 1,npt
!!$          write(91,'(4e18.8)') wvec(j)/twopi, abs(wseis_syn(j)), real(wseis_syn(j)), aimag(wseis_syn(j))
!!$       enddo
!!$       close(91)
!!$
!!$       open(91,file='test2.dat',status='unknown')
!!$       do j = 1,npt
!!$          write(91,'(4e18.8)') wvec(j)/twopi, abs(wseis_dat(j)), real(wseis_dat(j)), aimag(wseis_dat(j))
!!$       enddo
!!$       close(91)
!!$
!!$       open(91,file='test3.dat',status='unknown')
!!$       do j = 1,npt
!!$          write(91,'(4e18.8)') wvec(j)/twopi, abs(wseis_dat2(j)), real(wseis_dat2(j)), aimag(wseis_dat2(j))
!!$       enddo
!!$       close(91)
!!$
!!$       !wseis_dat(:) = 0.
!!$
!!$       ! back to real time series
!!$       call fftinv(lnpt,wseis_syn,REVERSE_FFT,DT,fp)    ! displacement -- from NUMERICAL
!!$       call fftinv(lnpt,wseis_dat,REVERSE_FFT,DT,fq)    ! displacement -- from ANALYTICAL
!!$       call fftinv(lnpt,wseis_dat2,REVERSE_FFT,DT,fq2)    ! should be velocity
!!$
!!$       open(91,file='test0.dat',status='unknown')
!!$       do j = 1,nlen
!!$          write(91,'(6e18.8)') j*DT, syn_displ(j), fp(j), fq(j), syn_veloc(j), fq2(j)
!!$       enddo
!!$       close(91)
!!$
!!$       adj_syn(:,:,:) = 0.
!!$
!!$    enddo  ! do icomp=1,NCOMP
!!$
!!$    !------------------------------------------------------------------
!!$  end subroutine mtm_test
!!$  !------------------------------------------------------------------
!!$  ! END TEST PROGRAM
!!$  !------------------------------------------------------------------

  subroutine mtm_adj(ievent, nrec, syn, tstart, tend, adj_syn, data, data_recon)

    double precision :: nw
    double precision, dimension(NSTEP) :: syn_displ, dat_displ, syn_veloc, dat_veloc, syn_accel
    double precision, dimension(NSTEP) :: ft_bar_t, ft_t, fa_bar_t, fa_t
    double precision :: Nnorm, Mnorm, dlnA, dlnAd, dlnAv
    character*100 filename

    integer, intent(in) :: nrec, ievent
    double precision, dimension(NSTEP,NCOMP,nrec),intent(in) :: syn
    double precision, dimension(nrec), intent(in) :: tstart, tend
    double precision, dimension(NSTEP,NCOMP,nrec),intent(out) :: adj_syn

    double precision, dimension(NSTEP,NCOMP,nrec),intent(in),optional :: data
    double precision, dimension(NSTEP,NCOMP,nrec),optional :: data_recon

    double precision, dimension(NDIM) :: dzr, dzr2, dzr_win, dzr2_win, dzr0_win

    integer :: itime, icomp, istart1, iend1, istart2, iend2, nlen
    integer :: n_left, n_right, ishift, i1, ipwr_t
    double precision :: cc, cc_max, tshift_xc, cr_shift, fac, sfac1, t_start, t_stop
    double precision, dimension(NSTEP) :: time_window

    ! error associated with measurement
    double precision :: dlnA_pert, tshift_xc_pert, rand1, ppert
    double precision :: twopi

    integer :: i,j,ik,n,ictaper,irec,imeasure

!====================================================================

!!$    ! PARAMETERS ARE INCLUDED IN wave2d_constants.f90
!!$    !parameter(NDIM=8000*4,NDIM2=16000*4,lnpt=14,npt=2**lnpt)
!!$    !parameter(NTAPER=5,WTR=0.02,ZZIGN=-1.0)
!!$    !implicit real(a-h,o-z)
!!$
!!$    double precision :: tas(NDIM,NTAPER), ey1(NDIM), ey2(NDIM)
!!$    !double precision, dimension(NDIM,NTAPER) :: phi_mul, dt_mul, abs_mul, abs2_mul
!!$    double precision, dimension(NDIM) :: dzr, dzr2, phi, fr
!!$    double precision, dimension(NDIM) :: dzr_win, dzr2_win, dzr3_win, dzr0_win, dzr30_win, tseis_recon
!!$    !double precision, dimension(NDIM) :: dt_mtm, phi_mtm, abs_mtm, abs2_mtm
!!$    !double precision, dimension(NDIM) :: err_phi, err_dt, err_abs, err_abs2
!!$
!!$    ! variables for adjoint sources (p(w) and q(w), etc)
!!$    complex*16, dimension(NDIM,NTAPER) :: top_p_ntaper
!!$    complex*16, dimension(NDIM) :: pwc_adj, qwc_adj, top_p, bot_p, ctemp, dtemp, dtau_pj_wc, dlnA_qj_wc
!!$    double precision, dimension(NDIM,NTAPER) :: pw_adj, qw_adj, pt_adj, qt_adj, dtau_pj_t, dlnA_qj_t
!!$    double precision, dimension(NDIM) :: wvec, tvec, dlnA_w, dtau_w, fp, fq, w_taper, wp_taper, wq_taper
!!$    double precision :: nw
!!$
!!$    ! cross-correlation sources
!!$    double precision, dimension(NSTEP) :: syn_displ, dat_displ, syn_veloc, dat_veloc, syn_accel
!!$    double precision, dimension(NSTEP) :: ft_bar_t, ft_t, fa_bar_t, fa_t
!!$    double precision :: Nnorm, Mnorm, dlnA, dlnAd, dlnAv
!!$
!!$    complex*16, dimension(NDIM) :: wseis_syn, wseis_dat, wseis, wseis3, wseis_rec
!!$    complex*16, dimension(NDIM) :: trans_mtm, top_mtm, bot_mtm   ! removed trans
!!$    complex*16 :: wtr_use, wtr_use_unw, cci, junkc
!!$    character*20 hfmt
!!$    character*100 filename
!!$
!!$    ! -----------------
!!$    integer, intent(in) :: nrec, ievent
!!$    double precision, dimension(NSTEP,NCOMP,nrec),intent(in) :: syn
!!$    double precision, dimension(nrec), intent(in) :: tstart, tend
!!$    double precision, dimension(NSTEP,NCOMP,nrec),intent(out) :: adj_syn
!!$
!!$    double precision, dimension(NSTEP,NCOMP,nrec),intent(in),optional :: data
!!$    double precision, dimension(NSTEP,NCOMP,nrec),optional :: data_recon
!!$
!!$    integer itime, icomp, istart1, iend1, istart2, iend2, nlen
!!$    double precision, dimension(NSTEP) :: time_window
!!$    double precision :: norm, junk, ntemp, Ffac, win_extend, planch_fac
!!$    double precision :: dtau_mtm, dlnA_mtm
!!$    ! -----------------
!!$
!!$    double precision :: twopi, wtr_mtm, fac, cr_shift, cc, cc_max, sfac1, sfac2
!!$    double precision :: tshift_xc, ampmax, ampmax_unw, omega
!!$    double precision :: t_stop, t_start, df, df_new, dw, temp
!!$    integer :: ipwr_t, ipwr_w, ishift, n_left, n_right, fnum, i1
!!$    integer :: i_right, i_right_stop, i_amp_max, i_amp_max_unw, idf_new
!!$
!!$    ! error associated with measurement
!!$    double precision :: meas_pert, rand1, ppert
!!$
!!$    integer :: i,j,ik,n,ictaper,irec

    !=======================================================

    twopi = 2.*PI
    !cci = cmplx(0.,1.)
    !wtr_mtm = 1.e-10

    ! parameters for tapers
    ipwr_t = 10                  ! for 1 - [cos(t)]^(ipwr)
    !ipwr_w = 10                  ! for 1 - [cos(w)]^(ipwr)

    adj_syn(:,:,:) = 0.
    !icomp = 1

    print *, ' (mtm_adj.f90) IKER = ', IKER

    ! KEY: global counter variable
    !imeasure = imeasure + 1

    !===============

    ! loop over receivers -- tstart and tend are nrec x 1 arrays
    do irec = 1,nrec

       ! loop over components
       do icomp = 1,NCOMP

          imeasure = index_data(ievent,irec,icomp)

          print *
          write(*,'(a30,4i10)') ' ievent, irec, icomp, imeasure: ', ievent, irec, icomp, imeasure

          ! create the NON-windowed records for which the measurement is made
          do itime = 1,NSTEP
             dzr(itime)  = syn(itime,icomp,irec)     ! synthetic
             dzr2(itime) = data(itime,icomp,irec)    ! data
          enddo

          ! window -- based on HWIN
          istart1 = max(floor( (tstart(irec)) / DT), 1)
          iend1   = min(floor( (tend(irec)  ) / DT), NSTEP)
          nlen = iend1 - istart1
          print *, istart1, iend1, nlen
          if (istart1 >= iend1) stop 'Check if istart1 < iend1'

          ! start and stop times -- assume t0 = 0
          !t_start = tstart(irec)
          !t_stop  = tend(irec)
          t_start = istart1*DT
          t_stop  = iend1*DT

          ! parameters for time window
          sfac1 = (2./dble(nlen))**2  ! for Welch window
          !ipwr_t = 2                  ! for cosine window

          ! initialise windowed arrays
          dzr_win(:) = 0.             ! length npt (npt is in wave2d_constants.f90)
          dzr0_win(:) = 0.            ! length npt
          dzr2_win(:) = 0.            ! length npt

          ! PRE-PROCESSING -- the type of window does not seem to matter too much
          ! (It is not necessary to fix the endpoints, I think.)
          ! Note that dzr_win is length npt, but dzr is length NSTEP
          do i = 1,nlen
             !fac = 1.                                        ! boxcar window
             !fac = 1 - sfac1*((i-1) - dble(nlen)/2.)**2     ! welch window
             fac = 1. - cos(PI*(i-1)/(nlen+1))**ipwr_t        ! cosine window
             !print *, fac

             dzr_win(i)  =  dzr(i+istart1) * fac  ! syn, windowed
             dzr2_win(i) = dzr2(i+istart1) * fac  ! dat, windowed
          enddo
          dzr0_win(:) = dzr_win(:)  ! syn, windowed

          ! post-processing taper, length NSTEP
          time_window(:) = 0.
          do i = 1,nlen
             !fac = 1.                                        ! boxcar window
             !fac = 1 - sfac2*((i-1) - dble(nlen)/2.)**2     ! welch window
             fac = 1. - cos(PI*(i-1)/(nlen+1))**ipwr_t        ! cosine window

             time_window(i) = fac
          enddo

          !==================================================================
          ! CROSS-CORRELATION MEASUREMENTS
          ! The number of time-steps we deal with is the length of the measurement window.

          ! TRAVELTIME: Check that the output is consistent with your convention.

          !cr_shift = 300.0
          cr_shift = 2.*nlen*DT        ! 14-Nov-2006

          n_left  = ceiling( -1.0 * cr_shift / DT )
          n_right = floor( cr_shift / DT )
          ishift = 0
          cc_max = 0.

          do i = n_left, n_right, 1
             cc = 0.
             do j = 1, nlen
                if ((j+i) > 1 .and. (j+i) < nlen) cc = cc + dzr_win(j) * dzr2_win(j+i)   ! cross-correlation
             enddo
             if ( cc > cc_max) then
                cc_max = cc
                ishift = i
             endif
          enddo
          tshift_xc = ishift*DT  ! KEY: cross-correlation time shift

          ! AMPLITUDE: Check that the output is consistent with your convention.

          !synt(:) = dzr0_win(:)
          !datt(:) = dzr2_win(:)
          syn_displ(:) = 0. ; dat_displ(:) = 0.
          syn_veloc(:) = 0. ; dat_veloc(:) = 0.
          syn_displ(1:nlen) = dzr0_win(1:nlen)     ! windowed synthetic displacement
          dat_displ(1:nlen) = dzr2_win(1:nlen)     ! windowed data displacement, UNSHIFTED

          ! calculate velocity and acceleration from synthetic displacement
          do i = 2, nlen-1
             syn_veloc(i) = (syn_displ(i+1) - syn_displ(i-1)) / (2.*DT)
             dat_veloc(i) = (dat_displ(i+1) - dat_displ(i-1)) / (2.*DT)
          enddo
          syn_veloc(1)    = (syn_displ(2) - syn_displ(1)) / DT
          syn_veloc(nlen) = (syn_displ(nlen) - syn_displ(nlen-1)) /DT
          dat_veloc(1)    = (dat_displ(2) - dat_displ(1)) / DT
          dat_veloc(nlen) = (dat_displ(nlen) - dat_displ(nlen-1)) /DT

          do i = 2, nlen-1
             syn_accel(i) = (syn_veloc(i+1) - syn_veloc(i-1)) / (2.*DT)
          enddo
          syn_accel(1)    = (syn_veloc(2) - syn_veloc(1)) / DT
          syn_accel(nlen) = (syn_veloc(nlen) - syn_veloc(nlen-1)) /DT

          ! cross-correlation amplitude MEASUREMENT
          ! definition of Dahlen and Baig (2002), Eq. 3,17,18 : dlnA = Aobs/Asyn - 1
          !dlnAd = sqrt( (DT * sum( dat_displ(:) * dat_displ(:) )) / (DT * sum( syn_displ(:) * syn_displ(:) )) ) - 1.
          !dlnAv = sqrt( (DT * sum( dat_veloc(:) * dat_veloc(:) )) / (DT * sum( syn_veloc(:) * syn_veloc(:) )) ) - 1.

          ! modified to avoid first-order approximation
          dlnAd = 0.5 * log( sum(dat_displ(:)*dat_displ(:)) / sum(syn_displ(:)*syn_displ(:)) )
          dlnAv = 0.5 * log( sum(dat_veloc(:)*dat_veloc(:)) / sum(syn_veloc(:)*syn_veloc(:)) )

          !----------------------------------------------
          ! ADD MEASUREMENT ERRORS

          ! meas_pert and measure_pert_vec are GLOBAL variables
          if ( ADD_DATA_ERRORS ) then
             meas_pert = measure_pert_vec(imeasure)
          else
             meas_pert = 0.0
          endif
          tshift_xc_pert = tshift_xc + meas_pert

          !----------------------------------------------
          ! CROSS-CORRELATION ADJOINT SOURCES

          ! normalization factors
          ! NOTE sign convention for N (differs from 2005 GJI paper) so that Nnorm >= 0
          Nnorm = -DT * sum( syn_displ(:) * syn_accel(:) )
          Mnorm =  DT * sum( syn_displ(:) * syn_displ(:) )

          ! cross-correlation traveltime adjoint source for banana-doughnut kernel
          ft_bar_t(:) = -syn_veloc(:) / Nnorm

          ! cross-correlation traveltime adjoint source for misfit kernel
          ! NOTE 1: sign convention
          ! NOTE 2: weighted by measurement and (diagonal) data covariance matrix
          ft_t(:)     = -( tshift_xc_pert / cov_data(imeasure) ) * ft_bar_t(:)

          ! cross-correlation amplitude adjoint source
          ! You have TWO OPTIONS: measure the amplitudes based on DISPLACEMENT or VELOCITY
          ! Default option has been IAMP_VEL = 0
          if (IAMP_VEL == 0) then     ! DISPLACEMENT
             fa_bar_t(:) = syn_displ(:) / Mnorm
             dlnA = dlnAd
          else                        ! VELOCITY
             fa_bar_t(:) = -syn_accel(:) / Nnorm
             dlnA = dlnAv
          endif
          fa_t(:) = -(dlnA / cov_data(imeasure) ) * fa_bar_t(:)    ! misfit kernel

          ! for now, we do not allow perturbations for the amplitude measurement
          dlnA_pert = dlnA

          if (0 == 1) then
             !print *
             print *, 'cross-correlation measurments:'
             print *, '       dT = ', tshift_xc
             print *, '   dlnA-d = ', dlnAd
             print *, '   dlnA-v = ', dlnAv
             print *, '        N = ', Nnorm
             print *, '        M = ', Mnorm
          else
             write(*,'(a8,i5,a8,1f18.10,a8,1f18.10)') &
                  'irec = ', irec, ', dT = ', tshift_xc, ', dA = ', dlnA
          endif

          ! additional files for checking (measure_socal_adj.m)
          if (WRITE_SEISMO_RECONSTRUCT) then
             ! time domain : time, data-disp, syn-disp, syn-vel, syn-accel
             write(filename,'(a,i5.5,a)') 'syn_time_', irec, '.dat'
             open(29,file=filename,status='unknown')
             do i = 1,nlen
                write(29,'(5e18.8)') i*DT, dat_displ(i), syn_displ(i), syn_veloc(i), syn_accel(i)
             enddo
             close(29)

             ! xcorr adjoint sources : traveltime (banana-doughnut), traveltime (misfit), amplitude (banana-doughnut), amplitude (misfit)
             write(filename,'(a,i5.5,a)') 'xcorr_time_', irec, '.dat'
             open(39,file=filename,status='unknown')
             do i = 1,nlen
                if (IAMP_VEL == 0) then
                   write(39,'(6e18.8)') ft_bar_t(i), ft_t(i), fa_bar_t(i), fa_t(i), &
                      -syn_accel(i)/Nnorm, -dlnAv*(-syn_accel(i)/Nnorm)
                else
                   write(39,'(6e18.8)') ft_bar_t(i), ft_t(i), fa_bar_t(i), fa_t(i), &
                      syn_displ(i)/Mnorm, -dlnAd*(syn_displ(i)/Mnorm)
                endif
             enddo
             close(39)
          endif

          !if (irec==4) stop 'testing'

!!$          !===================================================
!!$          ! MULTITAPER MEASUREMENTS
!!$
!!$          ! initialize multitaper measurements
!!$          dtau_mtm = 0.
!!$          dlnA_mtm = 0.
!!$
!!$          if (IKER==3 .or. IKER==4) then     ! multitaper measurements
!!$
!!$             ! calculate frequency step and number of frequencies
!!$             df = 1./(npt*DT)
!!$             dw = twopi / (npt*DT)
!!$             fnum = npt/2 + 1
!!$
!!$             ! KEY: assemble omega vector
!!$             wvec(:) = 0.
!!$             do j = 1,npt
!!$                if (j > npt/2) then
!!$                   wvec(j) = dw*(j-npt-1)   ! negative frequencies in second half
!!$                else
!!$                   wvec(j) = dw*(j-1)       ! positive frequencies in first half
!!$                endif
!!$             enddo
!!$
!!$             ! numerical factor for Plancherel theorem
!!$             ! THIS COST US QUITE A BIT OF TIME TO FIGURE OUT
!!$             planch_fac = dble(npt * DT * DT)
!!$
!!$             !----------------
!!$
!!$             ! apply time shift to DATA (observed seismogram)
!!$             print *, ' shift observed seismogram by (s) : ', tshift_xc
!!$             do i = 1, nlen
!!$                dzr3_win(i) = 0.
!!$                if ( (ishift+i) > 1 .and. (ishift+i) < nlen ) dzr3_win(i) = dzr2_win(i+ishift)
!!$                dzr30_win(i) = dzr3_win(i)
!!$             enddo
!!$
!!$             ! create complex synthetic seismogram and complex data seismogram
!!$             ! wseis_syn -- windowed synthetic record in freq domain
!!$             ! wseis_dat -- windowed data record, shifted by xcorr-tt, in freq domain
!!$             wseis_syn(:) = cmplx(0.,0.)
!!$             wseis_dat(:) = cmplx(0.,0.)
!!$             wseis_syn(1:nlen) =  dzr0_win(1:nlen)
!!$             wseis_dat(1:nlen) = dzr30_win(1:nlen)
!!$
!!$             ! apply FFT to get frequency domain versions (syn and data)
!!$             call fft(lnpt,wseis_syn,FORWARD_FFT,DT)
!!$             call fft(lnpt,wseis_dat,FORWARD_FFT,DT)
!!$
!!$             ! get the tapers
!!$             call staper(nlen, NPI, NTAPER, tas, NDIM, ey1, ey2)
!!$
!!$             ! format statement -- might be a problem with some Fortran packages
!!$             write(hfmt,'(a,i2.2,a)') '(', NTAPER,'e18.6)'
!!$
!!$             !------------------------------------------------------------------
!!$
!!$             ! initialize transfer function terms
!!$             top_mtm(:)   = cmplx(0.,0.)
!!$             bot_mtm(:)   = cmplx(0.,0.)
!!$             trans_mtm(:) = cmplx(0.,0.)
!!$
!!$             do ictaper = 1, NTAPER  ! loop over tapers
!!$
!!$                ! time domain: apply taper ictaper to synth and obs windowed seismograms
!!$                ! these get written over at each loop iteration of ictaper
!!$                do i = 1, nlen
!!$                   dzr_win(i)  =  dzr0_win(i) * tas(i,ictaper)   ! syn(t), single-tapered and windowed
!!$                   dzr3_win(i) = dzr30_win(i) * tas(i,ictaper)   ! dat(t), single-tapered and windowed
!!$                enddo
!!$
!!$                ! create complex seismograms
!!$                wseis(:)  = cmplx(0.,0.)
!!$                wseis3(:) = cmplx(0.,0.)
!!$                wseis(1:nlen)  = cmplx(dzr_win(1:nlen),0.)  ! syn(t), single-tapered and windowed
!!$                wseis3(1:nlen) = cmplx(dzr3_win(1:nlen),0.) ! dat(t), single-tapered and windowed
!!$
!!$                ! apply FFT to get complex spectra
!!$                call fft(lnpt,wseis, FORWARD_FFT,DT)   ! syn
!!$                call fft(lnpt,wseis3,FORWARD_FFT,DT)   ! dat
!!$
!!$                ! find max spectral power for single taper
!!$                ! --> could also use maxval and maxloc
!!$                ampmax = 0.
!!$                ampmax_unw = 0.
!!$                do i = 1, fnum   ! loop over frequencies
!!$                   if ( abs(wseis(i)) > ampmax) then              ! syn, single-tapered
!!$                      ampmax = abs(wseis(i))
!!$                      i_amp_max = i
!!$                   endif
!!$                   if ( abs(wseis_syn(i)) > ampmax_unw) then      ! syn
!!$                      ampmax_unw =  abs(wseis_syn(i))
!!$                      i_amp_max_unw = i
!!$                   endif
!!$                enddo
!!$
!!$                ! compute water level for single taper measurement
!!$                wtr_use     = cmplx(ampmax * WTR, 0.)       ! syn
!!$                wtr_use_unw = cmplx(ampmax_unw * WTR, 0.)   ! syn, single-tapered
!!$
!!$                ! determine i_right values using the power in the (untapered) synthetic
!!$                ! these variables define maximum frequency for measurement
!!$                ! i_right_stop = 1 --> stop at frequency i_right, not fnum
!!$                i_right = fnum
!!$                i_right_stop = 0
!!$                do i = 1,fnum             ! loop over frequencies
!!$                   if (i > i_amp_max_unw .and. abs(wseis_syn(i)) <= abs(wtr_use_unw) .and. i_right_stop == 0) then
!!$                      i_right_stop = 1
!!$                      i_right = i
!!$                   endif
!!$                   if (i > i_amp_max_unw .and. abs(wseis_syn(i)) >= 10.*abs(wtr_use_unw) .and. i_right_stop == 1) then
!!$                      i_right_stop = 0
!!$                      i_right = i
!!$                   endif
!!$                enddo  ! frequencies: i = 1,fnum
!!$
!!$                ! loop over frequencies
!!$                do i = 1,fnum
!!$
!!$                   ! calculate top and bottom of transfer function for multitapers
!!$                   ! NOTE THAT THESE QUANTITIES ARE SUMMED OVER THE TAPERS AS WELL
!!$                   top_mtm(i) = top_mtm(i) +  wseis3(i) * conjg(wseis(i))   ! uses data and syn
!!$                   bot_mtm(i) = bot_mtm(i) +  wseis(i)  * conjg(wseis(i))   ! uses syn only
!!$
!!$                   ! calculate transfer function for single taper measurement using water level
!!$                   ! CHT: trans IS NEVER USED HERE
!!$                   !if (abs(wseis(i)) >  abs(wtr_use)) trans(i) = wseis3(i) / wseis(i)
!!$                   !if (abs(wseis(i)) <= abs(wtr_use)) trans(i) = wseis3(i) / (wseis(i)+wtr_use)
!!$
!!$                enddo  ! frequencies: i = 1,fnum
!!$
!!$                !print *, ' taper number ', ictaper, ' out of ', NTAPER
!!$
!!$             enddo  ! tapers: ictaper = 1,NTAPER
!!$
!!$             ! it appears that this has a negligible effect on the recovery of seismograms
!!$             !i_right = fnum
!!$
!!$             !------------------------------------------------------------------
!!$             ! Multi-taper measurement, calculation, and output
!!$             !------------------------------------------------------------------
!!$
!!$             ! find water level for multi-taper measurement
!!$             ampmax = 0.
!!$             do i = 1, fnum
!!$                if ( abs(bot_mtm(i)) > ampmax) then
!!$                   ampmax =  abs(bot_mtm(i))
!!$                   i_amp_max = i
!!$                endif
!!$             enddo
!!$             wtr_use = cmplx(ampmax * wtr_mtm**2, 0.)       ! original expression
!!$             !wtr_use = cmplx(ampmax * 0.01, 0.)            ! Qinya test value
!!$
!!$             ! calculate transfer function using water level
!!$             !do i = 1, fnum
!!$             !  if (abs(bot_mtm(i)) >  abs(wtr_use)) trans_mtm(i) = top_mtm(i) / bot_mtm(i)
!!$             !  if (abs(bot_mtm(i)) <= abs(wtr_use)) trans_mtm(i) = top_mtm(i) / (bot_mtm(i)+wtr_use)
!!$             !enddo
!!$             do i = 1, fnum
!!$                if (abs(bot_mtm(i)) <= abs(wtr_use)) bot_mtm(i) = bot_mtm(i) + wtr_use
!!$             enddo
!!$             trans_mtm(1:fnum) = top_mtm(1:fnum) / bot_mtm(1:fnum)
!!$
!!$             !=======================================================
!!$             ! construct time series : tau(omega), dlnA(omega)
!!$
!!$             ! taper function for the FREQUENCY domain
!!$             nw = dble(i_right - 1)
!!$
!!$             ! loop to calculate phase and amplitude
!!$             ! NOTE: here we include the (cross-correlation) time shift
!!$             dtau_w(:) = 0.
!!$             dlnA_w(:) = 0.
!!$             w_taper(:) = 0.
!!$             do i = 2, i_right   ! start from 1 to avoid dividing by 0
!!$
!!$                !wvec(i)     = dw*i     ! do not divide by zero
!!$                !dtau_w(i)   = -atan2(aimag(trans_mtm(i)), real(trans_mtm(i))) + tshift_xc
!!$                !dtau_w(i)   = -(1./wvec(i)) * atan2(aimag(trans_mtm(i)), real(trans_mtm(i))) + tshift_xc
!!$
!!$                dtau_w(i)   = (-1./wvec(i)) * atan2(aimag(trans_mtm(i)), real(trans_mtm(i))) + tshift_xc
!!$                dlnA_w(i)   = abs(trans_mtm(i)) - 1.
!!$
!!$                ! type of filter in the freq domain : boxcar, welch, cosine
!!$                !w_taper(i) = 1.
!!$                !w_taper(i) = 1. - (2./nw)**2 * ((i-1) - nw/2.)**2
!!$                w_taper(i) = 1. - cos(PI*(i-1)/(nw+1))**ipwr_w
!!$             enddo
!!$
!!$             ! TESTING: create a flat transfer function to compare with adjoint sources
!!$             !dtau_w(:) = tshift_xc
!!$             !dlnA_w(:) = dlnA
!!$
!!$             ! compute normalization factor for W(w)
!!$             ! factor of 2 takes into account integration from -infty to +infty
!!$             Ffac = 2. * dw*sum(w_taper(:) )     ! crude integration
!!$
!!$             ! modified frequency-domain tapers (see Latex notes)
!!$             ! In theory, these would incorporate the error functions sigma_p and sigma_q.
!!$             wp_taper(:) = w_taper(:) / Ffac
!!$             wq_taper(:) = w_taper(:) / Ffac
!!$
!!$             if (WRITE_MTM_FILES) then
!!$                ! write transfer function to file
!!$                write(filename,'(a,i5.5,a)') 'transfer_freq_', irec, '.dat'
!!$                open(91,file=filename,status='unknown')
!!$                do i = 1, i_right
!!$                   write(91,'(4e18.8)') wvec(i)/twopi, dtau_w(i), dlnA_w(i), w_taper(i)
!!$                enddo
!!$                close(91)
!!$
!!$                write(filename,'(a,i5.5,a)') 'transfer_freq_int_', irec, '.dat'
!!$                open(91,file=filename,status='unknown')
!!$                write(91,*) Ffac
!!$                close(91)
!!$             endif
!!$
!!$             ! compute multitaper measurements (crude integration)
!!$             ! These formulas are equivalent to using a BOXCAR window for the taper,
!!$             ! and the values should be close to the cross-correlation values.
!!$             !dtau_mtm = 1. / (i_right*dw) * dw*sum( dtau_w(1:i_right) )
!!$             !dlnA_mtm = 1. / (i_right*dw) * dw*sum( dlnA_w(1:i_right) )
!!$             dtau_mtm = sum( dtau_w(1:i_right) ) / i_right
!!$             dlnA_mtm = sum( dlnA_w(1:i_right) ) / i_right
!!$
!!$             ! reconstruct data (wseis_rec) from synthetics (wseis_syn) using the transfer function (trans_mtm)
!!$             if (WRITE_SEISMO_RECONSTRUCT) then
!!$
!!$                ! Reconstruct mtm fit seismograms : syn*tran
!!$                ! d(w) = s(w) T(w) exp[-i w dT]
!!$                ! trans_mtm is for transferring syn --> SHIFTED data
!!$                wseis_rec(:) = cmplx(0.,0.)
!!$                do i = 1,i_right
!!$                   omega = wvec(i)
!!$
!!$                   ! mathematically, these are identical
!!$                   ! numerically, they produce nearly identical results
!!$                   wseis_rec(i) = wseis_syn(i)  * (1.+ dlnA_w(i)) * exp(-cci*omega*dtau_w(i))
!!$                   !wseis_rec(i) = wseis_syn(i) * trans_mtm(i) * exp(-cci*omega*tshift_xc)
!!$
!!$                   !wseis_rec(i) = wseis_syn(i)*exp(-cci*omega*tshift_xc)*trans_mtm(i)  ! sign
!!$                   !wseis_rec(i) = wseis_syn(i) * exp(-cci*omega*dtau_w(i)*w_taper(i)) * (1.+ dlnA_w(i)*w_taper(i))
!!$                enddo
!!$
!!$                ! inverse FFT into time domain
!!$                call fftinv(lnpt,wseis_rec,REVERSE_FFT,DT,tseis_recon)
!!$
!!$                write(filename,'(a,i5.5,a)') 'recovered_seis_time_', irec, '.dat'
!!$                open(17,file=filename,status='unknown')
!!$                do i = 1,nlen
!!$                   write(17,'(2e18.8)') DT*i, tseis_recon(i)
!!$                enddo
!!$                close(17)
!!$
!!$                ! write power spectra (a.k.a. spectral amplitude) to files
!!$                ! Note that the spectral power for a uniform time shift will not be apparent;
!!$                ! only the amplitudes will appear.
!!$                ! Note that the transfer function (trans_mtm) deals with the SHIFTED synthetics.
!!$                write(filename,'(a,i5.5,a)') 'recovered_seis_freq_', irec, '.dat'
!!$                open(91,file=filename,status='unknown')
!!$                do i = 1,i_right
!!$
!!$                   omega = wvec(i)
!!$                   wseis_rec(i) = wseis_syn(i) * (1.+ dlnA_w(i)) * exp(-cci*omega*dtau_w(i))
!!$                   !wseis_rec(i) = wseis_syn(i) * (1. + dlnA_mtm) * exp(-cci*omega*tshift_xc)  ! xcor shift
!!$                   !wseis_rec(i) = wseis_syn(i) * trans_mtm(i) * exp(-cci*omega*tshift_xc)
!!$
!!$                   ! w, syn, dat, dat-recon
!!$                   write(91,'(4e18.8)') wvec(i), abs(wseis_syn(i)), abs(wseis_dat(i)), abs(wseis_rec(i))
!!$                enddo
!!$                close(91)
!!$
!!$             endif
!!$
!!$             !endif  ! IKER==3,4
!!$
!!$             !===================================================
!!$             ! FFT TESTING: to determine what conventions we are using!
!!$             ! We do this by taking a function, s(t) and its time derivative,
!!$             ! and then computing fourier transforms.
!!$
!!$             if (0==1) then
!!$
!!$                ! create complex synthetic seismogram and complex data seismogram
!!$                wseis_syn(:) = cmplx(0.,0.)
!!$                wseis_dat(:) = cmplx(0.,0.)
!!$                wseis_syn(1:nlen) = cmplx(syn_displ(1:nlen),0.)  ! displacement
!!$                wseis_dat(1:nlen) = cmplx(syn_veloc(1:nlen),0.)  ! velocity
!!$
!!$                print *, '-----------------------------'
!!$                print *, FORWARD_FFT
!!$                print *, wseis_syn(100), wseis_dat(100)
!!$
!!$                ! FFT. to get frequency domain versions
!!$                call fft(lnpt,wseis_syn,FORWARD_FFT,DT)   ! displacement
!!$                call fft(lnpt,wseis_dat,FORWARD_FFT,DT)   ! velocity
!!$
!!$                print *, wseis_syn(100), wseis_dat(100)
!!$                print *, FORWARD_FFT
!!$                print *, '-----------------------------'
!!$
!!$                if (1==1) then    ! check Fourier convention
!!$
!!$                   ! check convention -- s_d(w)*iw should give velocity
!!$                   do i = 1,i_right   ! KEY: do not go too high frequency
!!$                      dtemp(i) = wseis_syn(i) * cmplx(0.,wvec(i))
!!$                      !dtemp(npt-i+1) = wseis_syn(npt-i+1) * cmplx(0.,wvec(npt-i+1))
!!$                   enddo
!!$
!!$                   ! back to real time series
!!$                   ! if our convention is consistent, then wp_taper should match syn_veloc
!!$                   call fftinv(lnpt,wseis_syn,REVERSE_FFT,DT,fp)    ! displacement
!!$                   call fftinv(lnpt,wseis_dat,REVERSE_FFT,DT,fq)    ! velocity
!!$                   call fftinv(lnpt,dtemp,REVERSE_FFT,DT,wp_taper)  ! should be velocity
!!$
!!$                   open(91,file='test.dat',status='unknown')
!!$                   do i = 1,nlen*2
!!$                      write(91,'(5e18.8)') syn_displ(i),fp(i),syn_veloc(i),fq(i),wp_taper(i)
!!$                   enddo
!!$                   close(91)
!!$
!!$                else            ! check Plancherel's theorem
!!$
!!$                   ! discrete form of Plancherel's theorem
!!$                   junk  = sum(syn_displ(:)*syn_veloc(:))
!!$                   junkc = (1./planch_fac) * sum( wseis_syn(:) * conjg(wseis_dat(:)) )
!!$
!!$                   print *, 'checking Plancherel theorem...'
!!$                   print *, junk
!!$                   print *, junkc
!!$                   print *, real(junkc)
!!$                   print *, planch_fac, 1./planch_fac
!!$
!!$                   print *, ' this should be 1 : ', (junk / real(junkc))
!!$                   print *
!!$                   print *, '   Nt = ', nlen
!!$                   print *, '   dt = ', DT
!!$                   print *, '  wNy = ', twopi/(2.*DT)
!!$                   print *, '  wRy = ', twopi/(nlen*DT)
!!$                   print *
!!$                   print *, '   Nw = ', npt
!!$                   print *, '   dw = ', dw
!!$
!!$                endif
!!$
!!$                stop 'testing the FFT from mtm_adj.f90'
!!$
!!$             endif
!!$
!!$             !==================================================================
!!$             ! MULTITAPER ADJOINT SOURCES
!!$
!!$             !if (IKER==3 .or. IKER==4) then
!!$
!!$             pw_adj(:,:) = 0.      ;  qw_adj(:,:) = 0.
!!$             pt_adj(:,:) = 0.      ;  qt_adj(:,:) = 0.
!!$             dtau_pj_t(:,:) = 0.   ;  dlnA_qj_t(:,:) = 0.
!!$             fp(:) = 0.            ;  fq(:) = 0.
!!$
!!$             bot_p(:) = cmplx(0.,0.)
!!$             top_p_ntaper(:,:) = cmplx(0.,0.)
!!$
!!$             ! This loops over the tapers to get the DENOMINATOR term of pj(w).
!!$             ! It also stores the displacement field sj(w), which is in the NUMERATOR term of p(w).
!!$             do ictaper = 1,NTAPER
!!$
!!$                ! apply TAPER ictaper to windowed synthetics
!!$                dzr_win(:) = 0.
!!$                do i = 1,nlen   ! time domain
!!$                   dzr_win(i) = dzr0_win(i) * tas(i,ictaper)
!!$                enddo
!!$
!!$                ! create complex seismograms (tapered synthetics, sj(w))
!!$                wseis(:) = cmplx(dzr_win(:),0.)
!!$
!!$                ! apply f.t. to get complex spectra
!!$                call fft(lnpt,wseis,FORWARD_FFT,DT)
!!$
!!$                ! bottom of p function for multitapers (complex)
!!$                do i = 1, i_right
!!$                   bot_p(i) = bot_p(i) + wseis(i)*conjg(wseis(i))
!!$                enddo
!!$
!!$                ! term in numerator (sj) (complex)
!!$                top_p_ntaper(:,ictaper) = wseis(:)
!!$
!!$             enddo  ! loop over tapers
!!$
!!$             do ictaper = 1,NTAPER   ! loop over tapers
!!$                !print *, ' taper number ', ictaper, ' out of ', NTAPER
!!$
!!$                top_p(:)   = top_p_ntaper(:,ictaper)  ! top of p function for multitapers
!!$                pwc_adj(:) = cmplx(0.,0.)
!!$                qwc_adj(:) = cmplx(0.,0.)
!!$
!!$                ! compute pj(w) and qj(w) -- j term is in top_p
!!$                ! (these get over-written each loop)
!!$                do i = 2, i_right  ! start from index 2 to avoid division by w=0
!!$                   !omega = dw*i   ! omega should not =0
!!$                   omega = wvec(i)
!!$
!!$                   ! formulas for Fourier convention FFT --> e^(-iwt)
!!$                   pwc_adj(i) = cmplx(0.,1./omega) * top_p(i) / bot_p(i)  ! (1/w)i = (-iw)/(-w^2) = -1/(iw)
!!$                   qwc_adj(i) = cmplx(0.,omega) * pwc_adj(i)              ! d/dt <--> -iw
!!$
!!$                   ! formulas for Fourier convention FFT --> e^(iwt)
!!$                   !pwc_adj(i) = cmplx(0.,-1./omega) * top_p(i) / bot_p(i)  ! (-1/w)i = (iw)/(-w^2) = 1/(iw)
!!$                   !qwc_adj(i) = cmplx(0.,-omega) * pwc_adj(i)              ! d/dt <--> -iw
!!$
!!$                enddo
!!$
!!$                ! SAVE A COPY HERE -- fftinv
!!$                ctemp(:) = pwc_adj(:)
!!$                dtemp(:) = qwc_adj(:)
!!$
!!$                ! EXTRA OUTPUT : IFFT into the time domain : pj(w) --> pj(t) and qj(w) --> qj(t)
!!$                if (WRITE_MTM_FILES) then
!!$                   call fftinv(lnpt,pwc_adj,REVERSE_FFT,DT,pt_adj(:,ictaper))
!!$                   call fftinv(lnpt,qwc_adj,REVERSE_FFT,DT,qt_adj(:,ictaper))
!!$                endif
!!$
!!$                ! incorporate measurement
!!$                ! create [dtau(w) pj(w) W(w)] and save [dtau(t) * pj(t) * W(t)]
!!$                ! create [dlnA(w) qj(w) W(w)] and save [dlnA(t) * qj(t) * W(t)]
!!$                dtau_pj_wc(:) = cmplx(0.,0.)
!!$                dlnA_qj_wc(:) = cmplx(0.,0.)
!!$                dtau_pj_wc(:) = twopi * ctemp(:) * cmplx(dtau_w(:),0.) * cmplx(wp_taper(:),0.)
!!$                dlnA_qj_wc(:) = twopi * dtemp(:) * cmplx(dlnA_w(:),0.) * cmplx(wq_taper(:),0.)
!!$                !dtau_pj_wc(:) = twopi * ctemp(:) * cmplx(wp_taper(:),0.)   ! no measurement
!!$                !dlnA_qj_wc(:) = twopi * dtemp(:) * cmplx(wq_taper(:),0.)   ! no measurement
!!$
!!$                ! IFFT into the time domain
!!$                call fftinv(lnpt,dtau_pj_wc,REVERSE_FFT,DT,dtau_pj_t(:,ictaper))
!!$                call fftinv(lnpt,dlnA_qj_wc,REVERSE_FFT,DT,dlnA_qj_t(:,ictaper))
!!$
!!$                ! create adjoint source
!!$                fp(:) = fp(:) + tas(:,ictaper) * dtau_pj_t(:,ictaper)
!!$                fq(:) = fq(:) + tas(:,ictaper) * dlnA_qj_t(:,ictaper)
!!$
!!$             enddo
!!$
!!$             if (WRITE_MTM_FILES) then
!!$                ! write banana-doughnut adjoint sources to file
!!$                write(filename,'(a,i5.5,a)') 'test_fadj_t_', irec, '.dat'
!!$                open(21,file=filename,status='unknown')
!!$                do i = 1,nlen
!!$                   write(21,*) sngl(fp(i)), sngl(fq(i))
!!$                enddo
!!$                close(21)
!!$
!!$                ! write banana-doughnut adjoint sources to file
!!$                open(21,file='test_fadj_t.dat',status='unknown')
!!$                do i = 1,nlen
!!$                   write(21,*) sngl(fp(i)), sngl(fq(i))
!!$                enddo
!!$                close(21)
!!$
!!$                ! time domain : tapers and other time series
!!$                open(18,file='test_hj_t.dat',status='unknown')
!!$                open(19,file='test_pj_t.dat',status='unknown')
!!$                open(20,file='test_Pj_t.dat',status='unknown')
!!$                open(21,file='test_Pj_hj_t.dat',status='unknown')
!!$                open(22,file='test_qj_t.dat',status='unknown')
!!$                open(23,file='test_Qj_t.dat',status='unknown')
!!$                open(24,file='test_Qj_hj_t.dat',status='unknown')
!!$
!!$                do i = 1,nlen  ! loop over time points
!!$                   write(18,hfmt) ( sngl(tas(i,ictaper)), ictaper=1,NTAPER )                         ! hj(t)
!!$                   write(19,hfmt) ( sngl(pt_adj(i,ictaper)), ictaper=1,NTAPER )                      ! pj(t)
!!$                   write(20,hfmt) ( sngl(dtau_pj_t(i,ictaper)), ictaper=1,NTAPER )                   ! Pj(t)
!!$                   write(21,hfmt) ( sngl(dtau_pj_t(i,ictaper) * tas(i,ictaper)), ictaper=1,NTAPER )  ! hj(t) Pj(t)
!!$
!!$                   write(22,hfmt) ( sngl(qt_adj(i,ictaper)), ictaper=1,NTAPER )                      ! qj(t)
!!$                   write(23,hfmt) ( sngl(dlnA_qj_t(i,ictaper)), ictaper=1,NTAPER )                   ! Qj(t)
!!$                   write(24,hfmt) ( sngl(dlnA_qj_t(i,ictaper) * tas(i,ictaper)), ictaper=1,NTAPER )  ! hj(t) Qj(t)
!!$                enddo
!!$
!!$                close(18) ; close(19) ; close(20) ; close(21) ; close(22) ; close(23) ; close(24)
!!$             endif
!!$
!!$          endif  ! if IKER==3,4

!=====================================================

!!$          ! generate a UNIFORM random number to simulate error in the measurement
!!$          ! ppert determines the range over which the perturbed measurement will be
!!$          ppert = 0.
!!$          call random_number(rand1)
!!$          meas_pert = 1. + ppert*(2.*rand1 - 1.)
!!$          !print *, tshift_xc, meas_pert, tshift_xc * meas_pert

          ! CREATE THE ADJOINT SOURCES
          ! HERE WE APPLY A TAPER TO FIX THE ENDPOINTS OF THE TIME WINDOW
          do i = 1,nlen   ! loop over points WITHIN the (OUTER) time window

             i1 = istart1 - 1 + i

!!$             ! store the reconstructed data: d'(w) = T(w) s(w)
!!$             if (WRITE_SEISMO_RECONSTRUCT) then
!!$                if (IKER==4 .or. IKER==4) data_recon(i1,icomp,irec) = tseis_recon(i)
!!$             endif

             if (IKER == 0) then
                adj_syn(i1,icomp,irec) = ( syn(i1,icomp,irec) -  data(i1,icomp,irec) ) * time_window(i)

             else if (IKER == 1) then
                adj_syn(i1,icomp,irec) = ft_t(i) * time_window(i)

             else if (IKER == 2) then
                adj_syn(i1,icomp,irec) = fa_t(i) * time_window(i)

             else if (IKER == 3) then
                stop 'Multitaper measurements NOT an option'
                !adj_syn(i1,icomp,irec) = fp(i) * time_window(i)

             else if (IKER == 4) then
                stop 'Multitaper measurements NOT an option'
                !adj_syn(i1,icomp,irec) = fq(i) * time_window(i)

             else if (IKER == 5) then
                adj_syn(i1,icomp,irec) = ft_bar_t(i) * time_window(i)

             else if (IKER == 6) then
                adj_syn(i1,icomp,irec) = fa_bar_t(i) * time_window(i)
             endif

          enddo

          ! (1) COMPUTE MEASUREMENT VECTOR
          ! (2) COMPUTE MISFIT function (currently only for waveform, xc-tt, xc-lnA)
          ! IKER: (0) waveform
          !       (1) traveltime, cross-correlation, misfit
          !       (2) amplitude, cross-correlation, misfit
          !       (3) traveltime, multitaper
          !       (4) amplitude, multitaper

!!$          ! the 2 factor is needed to offset the 2 factor in Ffac
!!$          dtau_mtm = 2./Ffac * dw*sum( dtau_w(1:i_right) )
!!$          dlnA_mtm = 2./Ffac * dw*sum( dlnA_w(1:i_right) )

          ! measurement values
          !imeasure = imeasure + 1    ! global counter variable

          measure_vec(imeasure,1) = tshift_xc_pert
          measure_vec(imeasure,2) = tshift_xc
          measure_vec(imeasure,3) = dlnA_pert
          measure_vec(imeasure,4) = dlnA
          measure_vec(imeasure,5) = 0.5 * DT*sum( adj_syn(:,icomp,irec)**2 )

!!$          measure_vec(imeasure,4) = dtau_mtm     ! integrated dtau(w) -- compare with tshift_xc
!!$          measure_vec(imeasure,5) = dlnA_mtm     ! integrated dlnA(w) -- compare with dlnA

         ! NOTE THAT THE FACTOR OF 0.5 IS NOT INCLUDED HERE

          if (IKER == 0) then
             ! crude integration of the waveform difference
             chi_data(ievent,irec,icomp,1) = DT*sum( adj_syn(:,icomp,irec)**2 ) / cov_data(imeasure)

          else if (IKER == 1) then
             chi_data(ievent,irec,icomp,1) = (tshift_xc_pert )**2 / cov_data(imeasure)

          else if (IKER == 2) then
             chi_data(ievent,irec,icomp,1) = (dlnA_pert)**2 / cov_data(imeasure)

!!$          else if (IKER==3) then
!!$             chi_data(ievent,irec,icomp,1) = 2.*dw*sum( wp_taper(1:i_right) * (dtau_w(1:i_right)**2) )
!!$
!!$          else if (IKER==4) then
!!$             chi_data(ievent,irec,icomp,1) = 2.*dw*sum( wq_taper(1:i_right) * (dlnA_w(1:i_right)**2) )

          endif

       enddo  ! loop over components (icomp=1,NCOMP)

    enddo  ! loop over receivers (irec=1,nrec)

    !------------------------------------------------------------------
  end subroutine mtm_adj
  !------------------------------------------------------------------
  ! END MAIN PROGRAM
  !------------------------------------------------------------------

!!$!------------------------------------------------------------------
!!$  subroutine fft(n,xi,zign,dt)
!!$! Fourier transform
!!$! This inputs AND outputs a complex function.
!!$! The convention is FFT --> e^(-iwt)
!!$!------------------------------------------------------------------
!!$      complex*16, dimension(NDIM) :: xi
!!$      integer :: n
!!$      double precision :: dt
!!$
!!$      complex*16 :: wk, hold, q
!!$      double precision :: m(25)
!!$      double precision :: zign,flx,v
!!$      integer :: lblock,k,fk,jh,ii,istart
!!$      integer :: l,iblock,nblock,i,lbhalf,j,lx
!!$
!!$      ! sign must be +1. or -1.
!!$      if (zign >= 0.) then
!!$        zign = 1.
!!$      else
!!$        zign = -1.
!!$      endif
!!$
!!$      lx = 2**n
!!$      do 1 i=1,n
!!$    1 m(i) = 2**(n-i)
!!$      do 4 l=1,n
!!$      nblock = 2**(l-1)
!!$      lblock = lx/nblock
!!$      lbhalf = lblock/2
!!$      k = 0
!!$      do 4 iblock=1,nblock
!!$      fk = k
!!$      flx = lx
!!$
!!$      v = zign*2.*PI*fk/flx         ! Fourier convention
!!$
!!$      wk = cmplx(cos(v),-sin(v))   ! sign change to -sin(v) 17-Nov-2006
!!$      istart = lblock*(iblock-1)
!!$
!!$      do 2 i=1,lbhalf
!!$      j  = istart+i
!!$      jh = j+lbhalf
!!$      q = xi(jh)*wk
!!$      xi(jh) = xi(j)-q
!!$      xi(j)  = xi(j)+q
!!$    2 continue
!!$
!!$      do 3 i=2,n
!!$      ii = i
!!$      if (k < m(i)) goto 4
!!$    3 k = k-m(i)
!!$    4 k = k+m(ii)
!!$      k = 0
!!$      do 7 j=1,lx
!!$      if (k < j) goto 5
!!$      hold = xi(j)
!!$      xi(j) = xi(k+1)
!!$      xi(k+1) = hold
!!$    5 do 6 i=1,n
!!$      ii = i
!!$      if (k < m(i)) goto 7
!!$    6 k = k-m(i)
!!$    7 k = k+m(ii)
!!$
!!$      ! final steps deal with dt factors
!!$      if (zign > 0.) then       ! FORWARD FFT
!!$         do i = 1,lx
!!$            xi(i) = xi(i)*dt   ! multiplication by dt
!!$         enddo
!!$
!!$      else                     ! REVERSE FFT
!!$         flx = flx*dt
!!$         do i = 1,lx
!!$            xi(i) = xi(i)/flx  ! division by dt
!!$         enddo
!!$      endif
!!$
!!$  end subroutine fft
!!$
!!$!------------------------------------------------------------------
!!$  subroutine fftinv(npow,s,zign,dt,r)
!!$! inverse Fourier transform -- calls fft
!!$!------------------------------------------------------------------
!!$
!!$      !implicit real*8(a-h,o-z)
!!$      !dimension r(4096*4)
!!$      !complex s(4096*4)
!!$
!!$      complex*16, intent(in) :: s(NDIM)
!!$      double precision, intent(out) :: r(NDIM)   ! note this is REAL
!!$
!!$      !complex*16 :: stemp(NDIM)
!!$      double precision :: dt,zign
!!$      integer :: npow, nsmp, nhalf, i
!!$
!!$      nsmp = 2**npow
!!$      nhalf = nsmp/2
!!$      call rspec(s,nhalf)   ! re-structuring
!!$
!!$      call fft(npow,s,zign,dt)    ! Fourier transform
!!$
!!$      do i = 1,nsmp
!!$        r(i) = real(s(i))     ! REAL part
!!$      enddo
!!$
!!$  end subroutine fftinv
!!$
!!$!------------------------------------------------------------------
!!$  subroutine rspec(s,np2)
!!$!------------------------------------------------------------------
!!$
!!$      !implicit real*8(a-h,o-z)
!!$      !complex s(4096*4)
!!$
!!$      complex*16 :: s(NDIM)
!!$      integer :: np2,n,n1,i
!!$
!!$      n = 2*np2
!!$      n1 = np2+1
!!$
!!$      s(n1) = 0.
!!$!     s(1)  = 0.
!!$      s(1)  = cmplx( real(s(1)),0.)
!!$
!!$      do i = 1,np2
!!$         s(np2+i) = conjg(s(np2+2-i))
!!$      enddo
!!$
!!$  end subroutine rspec
!!$
!!$!------------------------------------------------------------------
!!$  subroutine staper(nt, fw, nev, v, ndim, a, w)
!!$!------------------------------------------------------------------
!!$!$$$$ calls tsturm, root
!!$!  Slepian - Thomson multi-taper procedure
!!$!  Slepian, D.     1978  Bell Sys Tech J v57 n5 1371-1430
!!$!  Thomson, D. J.  1982  Proc IEEE v70 n9 1055-1096
!!$!    nt    the number of points in the series
!!$!    fw    the time-bandwidth product (number of Rayleigh bins)
!!$!    nev   the desired number of tapers
!!$!    v     the eigenvectors (tapers) are returned in v(.,nev)
!!$!    a, w  work arrays dimensioned at least nt long (nt+1, nt odd)
!!$!    a(1..nev) contains bandwidth retention factors on output.
!!$!  The tapers are the eigenvectors of the tridiagonal matrix sigma(i,j)
!!$!  [see Slepian(1978) eq 14 and 25.] They are also the eigenvectors of
!!$!  the Toeplitz matrix eq. 18. We solve the tridiagonal system in
!!$!  tsturm for the tapers and use them in Slepians eq 18 to get the
!!$!  bandwidth retention factors (i.e. the eigenvalues) Thomson's
!!$!  normalisation is used with no attention to sign.
!!$      !implicit real*8(a-h,o-z)
!!$      !dimension a(*),w(*),v(ndim,*)
!!$      !parameter (pi=3.14159265358979d0,r2=1.414213562373095d0)
!!$
!!$      integer :: nt, nev, ndim
!!$      double precision :: fw
!!$      double precision :: v(NDIM,NTAPER), a(NDIM), w(NDIM)
!!$
!!$      integer :: i,j,k,m
!!$      integer :: nxi, lh, lp1, neven, nodd, ntot, kk, kmax, nlow, nup
!!$      double precision :: r2,om,com,hn,asav,rbd,dc,sm,s,sn,vmax
!!$
!!$      !-------------------------
!!$
!!$      r2 = sqrt(2.)
!!$
!!$      if (nt < 2) return
!!$      nxi=mod(nt,2)
!!$      lh=(nt/2)+nxi
!!$      lp1=nt+1
!!$      om=2.*PI*fw/nt
!!$      com=cos(om)
!!$      hn=0.5*dble(lp1)
!!$      do 10 i=1,lh
!!$        a(i)=com*(i-hn)**2
!!$   10   w(i)=0.5*dble(i*(nt-i))
!!$      if (nxi == 0) then
!!$        asav=a(lh)-w(lh)
!!$        a(lh)=a(lh)+w(lh)
!!$        rbd=1./(a(lh)+w(lh-1))
!!$      else
!!$        asav=w(lh-1)
!!$        rbd=1./(w(lh)+w(lh-1))
!!$        w(lh-1)=r2*w(lh-1)
!!$      endif
!!$      do 15 i=1,lh
!!$        a(i+lh)=w(i)*rbd
!!$        w(i)=a(i+lh)**2
!!$   15   a(i)=a(i)*rbd
!!$      neven=max0((nev+1)/2,1)
!!$      nodd=nev-neven
!!$!  Do the even tapers
!!$      call tsturm(nt,lh,a,a(lh+1),w,neven,v,ndim,w(lh+1),0)
!!$      do 20 i=1,neven
!!$        k=2*i-1
!!$        if (nxi == 1) v(lh,k)=r2*v(lh,k)
!!$          do 20 j=1,lh
!!$   20     v(lp1-j,k)=v(j,k)
!!$      if (nodd <= 0) goto 34
!!$!  Do the odd tapers
!!$      if (nxi == 0) then
!!$        a(lh)=asav*rbd
!!$      else
!!$        a(nt)=asav*rbd
!!$        w(lh-1)=asav*asav
!!$      endif
!!$      call tsturm(nt,lh-nxi,a,a(lh+1),w,nodd,v,ndim,w(lh+1),1)
!!$      do 30 i=1,nodd
!!$        k=2*i
!!$        if (nxi == 1) v(lh,k)=0.
!!$          do 30 j=1,lh
!!$   30     v(lp1-j,k)=-v(j,k)
!!$   34 ntot=neven+nodd
!!$!  Calculate bandwidth retention parameters
!!$      dc=2.*com
!!$      sm=0.
!!$      s=sin(om)
!!$      w(1)=om/PI
!!$      w(2)=s/PI
!!$      do 35 j=3,nt
!!$        sn=dc*s-sm
!!$        sm=s
!!$        s=sn
!!$   35   w(j)=s/(PI*(j-1))
!!$      do 55 m=1,ntot
!!$        vmax=abs(v(1,m))
!!$        kmax=1
!!$        do 40 kk=2,lh
!!$          if (abs(v(kk,m)) <= vmax) goto 40
!!$          kmax=kk
!!$          vmax=abs(v(kk,m))
!!$   40     continue
!!$        a(m)=0.
!!$        nlow=kmax-1
!!$          do 45 j=1,nlow
!!$   45     a(m)=a(m)+w(j+1)*v(nlow+1-j,m)
!!$        nup=nt-nlow
!!$          do 50 j=1,nup
!!$   50     a(m)=a(m)+w(j)*v(nlow+j,m)
!!$   55 a(m)=a(m)/v(kmax,m)
!!$      return
!!$
!!$  end subroutine staper
!!$
!!$!------------------------------------------------------------------
!!$  subroutine tsturm(nt,n,a,b,w,nev,r,ndim,ev,ipar)
!!$!------------------------------------------------------------------
!!$!$$$$ calls root
!!$!  Uses bisection and Sturm counting to isolate the eigenvalues of the
!!$!  symmetric tridiagonal matrix with main diagonal a(.) and sub/super
!!$!  diagonal b(.).  Newton's method is used to refine the eigenvalue in
!!$!  subroutine root then direct recursion is used to get the eigenvector
!!$!  as this is always stable.  Note  ipar=0 for even tapers   =1 for odd
!!$!  tapers
!!$      !implicit real*8(a-h,o-z)
!!$      !parameter (epsi=1.d-15,epsi1=5.d-15)
!!$      !dimension a(*),b(*),ev(*),w(*),r(ndim,*)
!!$
!!$      double precision, parameter :: epsi = 1.d-15, epsi1 = 5.d-15
!!$
!!$      double precision, dimension(NDIM) :: a, b, w, ev
!!$      double precision, dimension(NDIM,NTAPER) :: r
!!$      integer :: nt,n,ndim,nev,ipar
!!$
!!$      double precision, dimension(NDIM) :: bb
!!$      double precision :: q,el,elam,u,umeps,x,ddot,rnorm
!!$      integer :: i,j,ik,iag,m,jk,jm1
!!$
!!$      !-------------------------
!!$
!!$      if (n <= 0 .or. nev <= 0) return
!!$      umeps=1.-epsi
!!$      do 5 i=1,nev
!!$    5 ev(i)=-1.
!!$      u=1.
!!$      do 1000 ik=1,nev
!!$      if (ik > 1) u=ev(ik-1)*umeps
!!$      el=min(ev(ik),u)
!!$   10 elam=0.5*(u+el)
!!$      if (abs(u-el) <= epsi1) goto 35
!!$      iag=0
!!$      q=a(1)-elam
!!$      if (q >= 0.) iag=iag+1
!!$      do 15 i=2,n
!!$      if (q == 0.) x=abs(b(i-1))/epsi
!!$      if (q /= 0.) x=w(i-1)/q
!!$      q=a(i)-elam-x
!!$      if (q >= 0.) iag=iag+1
!!$      if (iag > nev) goto 20
!!$   15 continue
!!$      if (iag >= ik) goto 20
!!$      u=elam
!!$      goto 10
!!$   20 if (iag == ik) goto 30
!!$      m=ik+1
!!$      do 25 i=m,iag
!!$   25 ev(i)=elam
!!$      el=elam
!!$      goto 10
!!$   30 el=elam
!!$      call root(u,el,elam,a,b,w,n,ik)
!!$   35 ev(ik)=elam
!!$      jk=2*ik+ipar-1
!!$      r(1,jk)=1.
!!$      r(2,jk)=-(a(1)-ev(ik))/b(1)
!!$      ddot=1.+r(2,jk)*r(2,jk)
!!$      jm1=2
!!$      do 45 j=3,n
!!$      r(j,jk)=-((a(jm1)-ev(ik))*r(jm1,jk)+b(j-2)*r(j-2,jk))/b(jm1)
!!$      ddot=ddot+r(j,jk)*r(j,jk)
!!$   45 jm1=j
!!$      rnorm=sqrt(nt/(2.*ddot))
!!$      do 50 j=1,n
!!$   50 r(j,jk)=r(j,jk)*rnorm
!!$ 1000 continue
!!$      return
!!$
!!$  end subroutine tsturm
!!$
!!$!------------------------------------------------------------------
!!$  subroutine root(u,el,elam,a,bb,w,n,ik)
!!$!------------------------------------------------------------------
!!$
!!$      !implicit real*8(a-h,o-z)
!!$      !parameter (epsi = 1.d-15, epsi1 = 5.d-15)
!!$      !dimension a(*),bb(*),w(*)
!!$
!!$      double precision, parameter :: epsi = 1.d-15, epsi1 = 5.d-15
!!$      double precision :: u,el,elam
!!$      double precision, dimension(NDIM) :: a,bb,w
!!$      integer :: n,ik
!!$
!!$      double precision :: an,b,bm,bn,del,x
!!$      integer :: i,iag
!!$
!!$      !----------------------
!!$
!!$    5 elam=0.5*(u+el)
!!$   10 if (abs(u-el) <= 1.5*epsi1) return
!!$      an=a(1)-elam
!!$      b=0.
!!$      bn=-1./an
!!$      iag=0
!!$      if (an >= 0.) iag=iag+1
!!$      do 20 i=2,n
!!$      if (an == 0.) x=abs(bb(i-1))/epsi
!!$      if (an /= 0.) x=w(i-1)/an
!!$      an=a(i)-elam-x
!!$      if (an == 0.) an=epsi
!!$      bm=b
!!$      b=bn
!!$      bn=((a(i)-elam)*b-bm*x-1.)/an
!!$      if (an >= 0.) iag=iag+1
!!$   20 continue
!!$      if (iag == ik) goto 25
!!$      u=elam
!!$      goto 30
!!$   25 el=elam
!!$   30 del=1./bn
!!$      if (abs(del) <= epsi1) del=sign(epsi1,del)
!!$      elam=elam-del
!!$      if (elam >= u .or. elam <= el) goto 5
!!$      goto 10
!!$
!!$  end subroutine root
!!$!-------------------------------------------

end module wave2d_sub4
