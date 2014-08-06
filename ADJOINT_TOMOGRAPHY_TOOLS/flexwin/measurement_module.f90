  module measurement_variables
  use seismo_variables
  implicit none

  ! measurement parameters / variables
  integer, parameter :: N_FREQ_MAX = 9000
  double precision, dimension (N_FREQ_MAX,NWINDOWS) :: fr, dphi_w, dtau_w, dlnA_w
  double precision, dimension (N_FREQ_MAX,NWINDOWS) :: dphi_w_err, dtau_w_err, dlnA_w_err
  integer, dimension (NWINDOWS) :: n_freq
  double precision, dimension (NDIM,NWINDOWS) :: syn_win, obs_win, obs_rec

  ! adjoint source parameters / variables
  logical, parameter :: CALC_ADJ=.true.
  double precision, dimension (NDIM,NWINDOWS) :: fp_adj_win, fq_adj_win

  end module measurement_variables


  ! functions and subroutines that are general and do not depend on type
  ! of measurement


! -----------------------------------------------------------------------

  subroutine calculate_windows_after_quality
  use seismo_variables
  use measurement_variables
  integer :: npts_win

  do i = 1, num_win
    npts_win=i_end(i)-i_start(i)+1
    call calc_criteria(obs_win(1:npts_win,i),obs_rec(1:npts_win,i),npts_win,&
                       1,npts_win,dt,Tshift_after(i),CC_after(i),dlnA_after(i))
  enddo
  end subroutine

! -----------------------------------------------------------------------

 
  subroutine write_measurements_gmt(basename)
  use seismo_variables
  use measurement_variables

  character*120 :: basename
  character*240 :: file_dtau, file_dlnA, file_seis, file_adj
  character*8 :: c_iwin

  integer :: i, iwin, npts_iwin
  double precision    :: max_freq, max_value_dtau, max_value_dlnA
  double precision    :: min_freq, min_value_dtau, min_value_dlnA
  double precision    :: max_seis, max_adj_p, max_adj_q, max_tmp1, max_tmp2


  do iwin = 1, num_win

    ! set filenames for output (use full path)
    call int2string(iwin,c_iwin)
    file_dtau=trim(basename)//'.dtau.'//trim(c_iwin)
    file_dlnA=trim(basename)//'.dlnA.'//trim(c_iwin)
    file_seis=trim(basename)//'.seis.win.'//trim(c_iwin)
    if(CALC_ADJ) file_adj=trim(basename)//'.adj.win.'//trim(c_iwin)

    ! get the extrema of the data to be written out
    min_freq = fr(1,iwin)
    max_freq = fr(n_freq(iwin),iwin)
    min_value_dtau = minval(dtau_w(1:n_freq(iwin),iwin)-dtau_w_err(1:n_freq(iwin),iwin))
    max_value_dtau = maxval(dtau_w(1:n_freq(iwin),iwin)+dtau_w_err(1:n_freq(iwin),iwin))
    min_value_dlnA = minval(dlnA_w(1:n_freq(iwin),iwin)-dlnA_w_err(1:n_freq(iwin),iwin))
    max_value_dlnA = maxval(dlnA_w(1:n_freq(iwin),iwin)+dlnA_w_err(1:n_freq(iwin),iwin))

    npts_iwin = i_end(iwin)-i_start(iwin)
    max_tmp1 = maxval(abs(syn_win(1:npts_iwin,iwin)))
    max_tmp2 = maxval(abs(obs_win(1:npts_iwin,iwin)))
    max_tmp1 = max(max_tmp1, max_tmp2)
    max_tmp2 = maxval(abs(obs_rec(1:npts_iwin,iwin)))
    max_seis = max(max_tmp1, max_tmp2)

    if(CALC_ADJ) then
      max_adj_p = maxval(abs(fp_adj_win(1:npts_iwin,iwin)))
      max_adj_q = maxval(abs(fq_adj_win(1:npts_iwin,iwin)))
    endif

    

    ! write the seismograms and envelopes and f1f2
    ! open the files
    open(unit=11, file=file_dtau)
    open(unit=12, file=file_dlnA)
    open(unit=13, file=file_seis)
    if(CALC_ADJ) open(unit=14, file=file_adj)
    ! write the header - dtau
    write(11,'("# NPTS   = ",i10)') n_freq(iwin)
    write(11,'("# F_START = ",e12.6)') min_freq
    write(11,'("# F_END = ",e12.6)') max_freq
    write(11,'("# DTAU_MIN = ",e12.6)') min_value_dtau
    write(11,'("# DTAU_MAX = ",e12.6)') max_value_dtau
    write(11,'("# CC = ",f10.4)') CC(iwin)
    ! write the header - dlnA
    write(12,'("# NPTS   = ",i10)') n_freq(iwin)
    write(12,'("# F_START = ",e12.6)') min_freq
    write(12,'("# F_END = ",e12.6)') max_freq
    write(12,'("# DLNA_MIN = ",e12.6)') min_value_dlnA
    write(12,'("# DLNA_MAX = ",e12.6)') max_value_dlnA
    write(12,'("# CC = ",f10.4)') CC(iwin)
    ! write the header - seismograms
    write(13,'("# NPTS   = ",i10)') npts_iwin
    write(13,'("# MIN_PERIOD = ",f12.6)') 1/max_freq
    write(13,'("# PLOT_MAX = ",e12.6)') max_seis
    write(13,'("# T_START = ",f10.2)') win_start(iwin)
    write(13,'("# T_END = ",f10.2)') win_end(iwin)
    write(13,'("# CC = ",f10.4)') CC(iwin)
    write(13,*) '# SYN OBS OBS_REC'
    if(CALC_ADJ) then
    ! write the header - adjoint
      write(14,'("# NPTS   = ",i10)') npts_iwin
      write(14,'("# P_MAX = ",e12.6)') max_adj_p
      write(14,'("# Q_MAX = ",e12.6)') max_adj_q
      write(14,'("# T_START = ",f10.2)') win_start(iwin)
      write(14,'("# T_END = ",f10.2)') win_end(iwin)
      write(14,'("# CC = ",f10.4)') CC(iwin)
      write(14,*) '# FP_ADJ FQ_ADJ'
    endif
 
 
    ! write the measurement
    do i = 1, n_freq(iwin)
      write(11,'(e12.6,2(2x,e12.6))') fr(i,iwin), dtau_w(i,iwin), dtau_w_err(i,iwin)
      write(12,'(e12.6,2(2x,e12.6))') fr(i,iwin), dlnA_w(i,iwin), dlnA_w_err(i,iwin)
    enddo

    ! write the seismograms and ajoint sources
    do i = 1, npts_iwin
      write(13,'(4(2x,e12.6))') b+(i_start(iwin)-1)*dt+(i-1)*dt, syn_win(i,iwin), obs_win(i,iwin), obs_rec(i,iwin)
      if(CALC_ADJ) write(14,'(3(2x,e12.6))') b+(i_start(iwin)-1)*dt+(i-1)*dt, fp_adj_win(i,iwin), fq_adj_win(i,iwin)
    enddo

    ! close the files
    close(11)
    close(12)
    close(13)
    if(CALC_ADJ) close(14)

  end do

  end subroutine write_measurements_gmt

 
