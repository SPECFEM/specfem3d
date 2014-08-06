!
! $Id:$
!
!----------------------------------------------------------------------


  subroutine measure_windows_xcorr
  use xcorr_constants
  use seismo_variables
  use measurement_variables
  implicit none

  double precision, dimension(npt) :: syn_lp_local, obs_lp_local
  double precision :: dtau_local, dlnA_local
  double precision :: CC_after_local, dlnA_after_local, tshift_after_local
  double precision, dimension(npt) :: obs_win_local, syn_win_local, obs_rec_local

  integer :: iwin, nlen
  double precision :: lp_min, lp, win_length
  logical :: quality_ok

  double precision, dimension(npt) :: fp_local, fq_local


  
! Loop over windows
! ---------------------------------------------------------------------
  write(*,*) "Cross-correlation measurements"
  do iwin = 1, num_win

    win_length = win_end(iwin) - win_start(iwin)
    nlen=win_length/dt

!   initialise global measurements 
!   ---------------------------------------------------------------------
    n_freq(iwin) = 0 
    fr(:,iwin) = 0
    dtau_w(:,iwin) = 0; dlnA_w(:,iwin) = 0
    F1_after(iwin) = 0; F2_after(iwin) = 0
    Tshift_after(iwin) = 0; CC_after(iwin) = 0 ; dlnA_after = 0


!   minimum lowpass corner frequency for measurement = WIN_LP_FREQ
!   -------------------------------------------------------------------
    lp_min = 1.0/WIN_LP_PERIOD
    lp=lp_min

    quality_ok=.true.
!   initialise all local measurements and erros to zero
!   ---------------------------------------------------------------------
    dtau_local = 0 ; dlnA_local = 0
    obs_win_local(:) = 0 ; syn_win_local(:) = 0; obs_rec_local(:) = 0
    fp_local(:) = 0 ; fq_local(:) = 0
    obs_lp_local(:)=0; syn_lp_local(:)=0
 
    obs_lp_local(1:npts)=obs_lp(1:npts)
    syn_lp_local(1:npts)=synt_lp(1:npts)

!   make the xcorr measurement
!   -------------------------------------------------------------------
    call xcorr_measure(win_start(iwin),win_end(iwin),&
                       syn_lp_local,obs_lp_local,npts,dt,b,&
                       tshift_after_local,CC_after_local,dlnA_after_local,&
                       dtau_local,dlnA_local,&
                       obs_win_local,syn_win_local,obs_rec_local,&
                       CALC_ADJ,fp_local,fq_local,DEBUG)

!   ---------------------------------------------------------------------

!   measurements
    n_freq(iwin) = 1
    fr(1,iwin) = lp
    dtau_w(1,iwin) = dtau_local
    dlnA_w(1,iwin) = dlnA_local


    Tshift_after(iwin) = tshift_after_local
    CC_after(iwin)=CC_after_local
    dlnA_after(iwin)=dlnA_after_local

    syn_win(1:nlen,iwin) = syn_win_local(1:nlen)
    obs_win(1:nlen,iwin) = obs_win_local(1:nlen)
    obs_rec(1:nlen,iwin) = obs_rec_local(1:nlen)

! version without frequency iteration
    write(*,'("Window ",i2," : dTau = ",f5.2," dlnA = ",f5.3," CC = ",f5.3," lp at ",f6.2," s ")') iwin, &
    dtau_local, dlnA_local, CC_after_local, 1/(lp)


    if(calc_adj) then
      fp_adj_win(1:nlen,iwin) = fp_local(1:nlen)
      fq_adj_win(1:nlen,iwin) = fq_local(1:nlen)
    endif

  enddo


  end subroutine measure_windows_xcorr
