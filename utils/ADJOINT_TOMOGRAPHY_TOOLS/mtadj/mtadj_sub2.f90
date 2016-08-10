module mtadj_sub2

  use mtadj_constants
  use mtadj_sub3

  implicit none

contains

  subroutine compute_time_shift(synw,dataw,nlen,dt,ishift,tshift_cc)

    real, intent(in),dimension(:) :: synw, dataw
    integer, intent(in) :: nlen
    real, intent(in) :: dt
    real, intent(out) :: tshift_cc
    integer, intent(out) :: ishift

    real :: cc,cc_max
    integer :: i,j

    ! these choices will slide the entire windowed record past the other
    if (nlen <= 0) stop 'Error compute_time_shift(): nlen has to be > 0'

    ! cross-correlation (i should be between -nlen+1 to nlen+1)
    ishift = 0
    do i = -nlen, nlen, 1
       cc = 0.
       do j = 1, nlen
          if ((j+i) >= 1 .and. (j+i) <= nlen) cc = cc + synw(j) * dataw(j+i)
       enddo
       if (cc > cc_max) then
          cc_max = cc
          ishift = i
       endif

    enddo
    tshift_cc = ishift*dt

  end subroutine compute_time_shift

  ! ===========================================================================

  subroutine reconstruct_syn_cc(synw,nlen,dt,ishift,dlnA,synw_rc_cc)

    real,intent(in),dimension(:) :: synw
    integer, intent(in) :: nlen, ishift
    real, intent(in) :: dt, dlnA
    real,intent(out),dimension(:) :: synw_rc_cc

    real,dimension(NPT) :: synw_shifted
    integer :: i

    ! shift synthetics
    synw_shifted(:) = synw(:)
    do i = 1, nlen
       if ((i-ishift) >= 1 .and. (i-ishift) <= nlen ) synw_shifted(i) = synw(i-ishift)
    enddo
    ! fill the missing time window with the endpoint value
    if (ishift > 0) synw_shifted(1:ishift) = synw(ishift+1)
    if (ishift < 0) synw_shifted(nlen+ishift:nlen) = synw(nlen+ishift-1)

    ! multiplying amplitude factor
    synw_rc_cc(1:nlen) = synw_shifted(1:nlen) * exp( dlnA )

  end subroutine reconstruct_syn_cc

  ! ======================================================================

  ! this subroutine computes the cross-correlation error
  ! sigma_tshift_cc and sigma_dlnA_cc to take into account
  ! the data and synthetics misfit

  subroutine compute_cc_error(dataw,synw,nlen,dt,i_pmax,dlnA, &
       sigma_tshift_cc,sigma_dlnA_cc, &
       MIN_SIGMA_TSHIFT_CC,MIN_SIGMA_DLNA_CC)

    real,dimension(:),intent(in) :: dataw, synw
    integer, intent(in) :: nlen
    real, intent(in) :: dt, dlnA, MIN_SIGMA_TSHIFT_CC, MIN_SIGMA_DLNA_CC
    integer, intent(in) :: i_pmax
    real, intent(out) :: sigma_tshift_cc, sigma_dlnA_cc

    ! these two parameters seem to be arbitrary, but I guess it is better
    ! to have some estimates than not at all!
    real, parameter :: SIGMA_TSHIFT_CC_SCALE = 0.001
    real, parameter :: SIGMA_DLNA_CC_SCALE = 2.0
    real :: df,Tmax,misfit

    df=1./(NPT*dt)
    if (i_pmax <= 0) stop 'Check if the maximum f has been identified'
    Tmax=1./(df*i_pmax)

    misfit=sum((dataw(1:nlen)-synw(1:nlen))**2)/sum(dataw(1:nlen)**2)

    sigma_tshift_cc = max(MIN_SIGMA_TSHIFT_CC,Tmax*misfit*SIGMA_TSHIFT_CC_SCALE)
    sigma_dlnA_cc = max(MIN_SIGMA_dlnA_CC,abs(dlnA)*misfit*SIGMA_TSHIFT_CC_SCALE)

    print *,

  end subroutine compute_cc_error

  ! ======================================================================

  subroutine reconstruct_syn_fd(csynw,dtau_fdm,dlnA_fdm,i_right, &
       synw_rc_fd,dt,nlen)

    complex, dimension(:), intent(in) :: csynw
    real, dimension(:), intent(in) :: dtau_fdm, dlnA_fdm
    integer, intent(in) :: i_right, nlen
    real,dimension(:), intent(out) :: synw_rc_fd
    real, intent(in) :: dt

    real :: df,fvec(NPT)
    integer :: nf,j
    complex*16,dimension(NPT) :: csynw_recon
    real*8,dimension(NPT) :: synw_rc_fd_dp

    df = 1./(dt*NPT); nf=floor(NPT/2.)+1
    do j = 1, nf
       fvec(j) = df*(j-1)
    enddo
    csynw_recon = cmplx(0.,0.)
    csynw_recon(1:i_right)=csynw(1:i_right) * (1.+ dlnA_fdm(1:i_right))&
         * exp(-CCI*2*pi*fvec(1:i_right)*dtau_fdm(1:i_right))

    call fftinv(LNPT,csynw_recon,REVERSE_FFT,dble(dt),synw_rc_fd_dp)
    synw_rc_fd(1:nlen) = synw_rc_fd_dp(1:nlen)
    synw_rc_fd(nlen+1:NPT) = 0.


  end subroutine reconstruct_syn_fd

  ! ======================================================================

  subroutine compute_dtau_dlnA(trans_fdm,dt,tshift_cc,dtau_fdm,dlnA_fdm,i_right)

    complex, dimension(:),intent(in) :: trans_fdm
    real,intent(in):: tshift_cc,dt
    real,dimension(:), intent(out) :: dtau_fdm, dlnA_fdm
    integer, intent(in) :: i_right

    real :: df, smth, smth1, smth2
    real, dimension(NPT) :: fvec, phi_wt, abs_wt
    integer :: nf, j, i

    df = 1./(dt*NPT); nf=floor(NPT/2.)+1
    do j = 1, nf
       fvec(j) = df*(j-1)
    enddo
    ! obtain amplitude and phase first
    do i = 1, i_right
      phi_wt(i) = atan2( aimag(trans_fdm(i)) , real(trans_fdm(i)) )
      abs_wt(i) = abs(trans_fdm(i))
   enddo

   ! convert phase to travel time delay
   dtau_fdm(1) = 0.  ! tshift_cc is probably not a good choice
   do i = 1, i_right

      dlnA_fdm(i) = abs_wt(i) - 1.
      !
      if (i > 1) dtau_fdm(i) = (-1./(TWOPI*fvec(i))) * phi_wt(i) + tshift_cc

      ! right now phi_wt [-180,180], adjust abrupt phase discontinuities [-90,90]
      if (i > 1 .and. i < i_right) then
        ! check the smoothness (2nd-order derivative) by 2*pi changes
         smth  =  phi_wt(i+1) + phi_wt(i-1) -  2.0 * phi_wt(i)
         smth1 = (phi_wt(i+1) + TWOPI) + phi_wt(i-1) -  2.0 * phi_wt(i)
         smth2 = (phi_wt(i+1) - TWOPI) + phi_wt(i-1) -  2.0 * phi_wt(i)
         if (abs(smth1) < abs(smth) .and. abs(smth1) < abs(smth2) .and. abs(phi_wt(i) - phi_wt(i+1)) > PHASE_STEP) then
            if (DEBUG) print *, '2 pi phase correction:', fvec(i), phi_wt(i) - phi_wt(i+1)
            do j = i+1, i_right
               phi_wt(j) = phi_wt(j) + TWOPI
            enddo
         endif
         if (abs(smth2) < abs(smth) .and. abs(smth2) < abs(smth1) .and. abs(phi_wt(i) - phi_wt(i+1)) > PHASE_STEP) then
            if (DEBUG) print *, '-2 pi phase correction:', fvec(i), phi_wt(i) - phi_wt(i+1)
            do j = i+1, i_right
               phi_wt(j) = phi_wt(j) - TWOPI
            enddo
         endif
      endif
!      write(21,*) fvec(i), aimag(trans_fdm(i)),real(trans_fdm(i)),phi_wt(i)/PI*180
   enddo

  end subroutine compute_dtau_dlnA

  ! ======================================================================
  subroutine compute_mt_error(ntaper,dataw,synw,tas, &
       nlen,dt,wtr_mtm,i_right,tshift_cc, &
       dtau_fdm,dlnA_fdm, sigma_dtau_fdm,sigma_dlnA_fdm)

    integer, intent(in) :: ntaper,  nlen, i_right
    real,dimension(:),intent(in) :: dataw, synw
    real,dimension(:,:),intent(in) :: tas
    real,intent(in) :: dt,tshift_cc,wtr_mtm
    real,dimension(:),intent(in) :: dtau_fdm, dlnA_fdm
    real,dimension(:),intent(out) :: sigma_dtau_fdm, sigma_dlnA_fdm

    integer :: iom, i, nf, ictaper
    complex,dimension(NPT) :: top_mtm, bot_mtm, trans_mtm
    complex*16,dimension(NPT) :: cdatawt, csynwt
    real,dimension(NPT,NMAX_TAPER) :: dtau_mtm, dlnA_mtm
    real :: ampmax_bot, wtr_amp_bot, edt_ave, eabs2_ave, edt_iom, eabs2_iom
    real,dimension(NPT)::err_dt,err_dlnA ,datawt,synwt


    nf=floor(NPT/2.)+1

    do iom = 1, ntaper

       top_mtm(:) = cmplx(0.,0.)
       bot_mtm(:) = cmplx(0.,0.)

       do ictaper = 1, ntaper
          if (ictaper == iom) cycle

          ! apply ictaper'th taper
          datawt(1:nlen) = dataw(1:nlen) * tas(1:nlen,ictaper)
          synwt(1:nlen) = synw(1:nlen) * tas(1:nlen,ictaper)

          ! complex tapered series
          cdatawt(:) = cmplx(0.,0.)
          csynwt(:) = cmplx(0.,0.)
          cdatawt(1:nlen) = cmplx(datawt(1:nlen),0)
          csynwt(1:nlen) = cmplx(synwt(1:nlen),0)

          ! apply f.t. to get complex spectra
          call fft(LNPT,cdatawt,FORWARD_FFT,dble(dt))
          call fft(LNPT,csynwt,FORWARD_FFT,dble(dt))

          ! calculate top and bottom of Jacknife transfer function
          do i = 1, nf
             top_mtm(i) = top_mtm(i) +  cdatawt(i) * conjg(csynwt(i))
             bot_mtm(i) = bot_mtm(i) +  csynwt(i) * conjg(csynwt(i))
          enddo
       enddo ! ictaper

       ! water level
       ampmax_bot = maxval(abs(bot_mtm(1:nf)))
       wtr_amp_bot = ampmax_bot * wtr_mtm ** 2

       !  calculate transfer function using water level
       do i = 1, i_right
          if (abs(bot_mtm(i)) >= abs(wtr_amp_bot)) trans_mtm(i) = top_mtm(i) / bot_mtm(i)
          if (abs(bot_mtm(i)) < abs(wtr_amp_bot)) trans_mtm(i) = top_mtm(i) /(bot_mtm(i)+wtr_amp_bot)
       enddo
       call compute_dtau_dlnA(trans_mtm,dt,tshift_cc,dtau_mtm(:,iom),dlnA_mtm(:,iom),i_right)
    enddo ! iom

    err_dt   = 0.
    err_dlnA = 0.
    do i = 1, i_right

       edt_ave   = 0.
       eabs2_ave = 0.

       do iom = 1, ntaper
          edt_iom = ntaper*dtau_fdm(i) - (ntaper-1)*dtau_mtm(i,iom)
          edt_ave = edt_ave + edt_iom

          eabs2_iom = ntaper*dlnA_fdm(i) - (ntaper-1)*dlnA_mtm(i,iom)
          eabs2_ave = eabs2_ave + eabs2_iom
       enddo

       edt_ave   = edt_ave   / (ntaper)
       eabs2_ave = eabs2_ave / (ntaper)

       do iom = 1, ntaper
          err_dt(i)   = err_dt(i)  + (dtau_mtm(i,iom) - edt_ave)**2
          err_dlnA(i) = err_dlnA(i)+ (dlnA_mtm(i,iom) - eabs2_ave)**2
       enddo

       err_dt(i)   =  sqrt( err_dt(i) / (ntaper * (ntaper-1) ) )
       ! set the error bar for the first point corresponding to
       ! static offset to be large, which makes no contribution to
       ! the adjoint source
       if (i == 1) err_dt(i) = LARGE_VAL
       err_dlnA(i) =  sqrt( err_dlnA(i) / (ntaper * (ntaper-1) ) )

    enddo ! i_right
    sigma_dtau_fdm(1:i_right) = err_dt(1:i_right)
    sigma_dlnA_fdm(1:i_right) = err_dlnA(1:i_right)

  end subroutine compute_mt_error

  ! ======================================================================
  subroutine compute_fd_error(npi,nlen,i_right,dt,dtau_fdm,dlnA_fdm, &
       sigma_dtau_fdm,sigma_dlnA_fdm)

    real,intent(in) :: npi
    integer,intent(in) :: nlen, i_right
    real,intent(in) :: dt
    real,dimension(:),intent(in) :: dtau_fdm,dlnA_fdm
    real,dimension(:),intent(inout) :: sigma_dtau_fdm, sigma_dlnA_fdm

    real :: df
    integer :: idf_new,i

    df = 1./(dt*NPT) * npi
    idf_new = floor(1./(nlen)/df)
    do i=1, i_right-idf_new
       sigma_dtau_fdm(i) = maxval(abs(dtau_fdm(i)-dtau_fdm(i:i+idf_new)))
       sigma_dlnA_fdm(i) = maxval(abs(dlnA_fdm(i)-dlnA_fdm(i:i+idf_new)))
    enddo
    sigma_dtau_fdm(i_right-idf_new+1:i_right) = sigma_dtau_fdm(i_right-idf_new)
    sigma_dlnA_fdm(i_right-idf_new+1:i_right) = sigma_dlnA_fdm(i_right-idf_new)

  end subroutine compute_fd_error

  ! ======================================================================

  subroutine compute_veloc_from_displ(synw,nlen,dt,synw_veloc)

    real, dimension(:), intent(in) :: synw
    integer, intent(in) :: nlen
    real, intent(in) :: dt
    real,dimension(:),intent(out) :: synw_veloc

    integer i

    do i = 2, nlen-1
       synw_veloc(i) = (synw(i+1) - synw(i-1)) / (2.0*dt)
    enddo
    synw_veloc(1)    = (synw(2) - synw(1)) / dt
    synw_veloc(nlen) = (synw(nlen) - synw(nlen-1)) /dt

  end subroutine compute_veloc_from_displ

  ! ======================================================================

  subroutine interp_adj_src(dt_adj_src,nlen,tstart,dt, &
       dt_adj_src_win,npts_adj,b_adj,dt_adj)

    real,intent(in) :: tstart,dt,b_adj,dt_adj
    real,intent(in),dimension(:) :: dt_adj_src
    integer,intent(in) :: nlen,npts_adj
    real,intent(out),dimension(:) :: dt_adj_src_win

    integer :: i,ii
    real :: time, tt

    ! initialize to zero
    dt_adj_src_win(1:npts_adj) = 0.

    do i = 1, npts_adj
       time = b_adj + (i-1) * dt_adj
       if (time > tstart .and. time < tstart + (nlen-1)*dt) then
          ii = floor((time-tstart)/dt) + 1
          tt = time - ((ii-1)*dt + tstart)
          dt_adj_src_win(i) = (dt_adj_src(ii+1)-dt_adj_src(ii)) * tt/dt + dt_adj_src(ii)
       endif
    enddo

  end subroutine interp_adj_src

  ! ======================================================================

end module mtadj_sub2
