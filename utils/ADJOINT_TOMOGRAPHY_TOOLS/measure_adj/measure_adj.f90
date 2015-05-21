program measure_adj

  !  main program that calls the subroutines to make measurements within input time windows
  !  and then compute the corresponding adjoint sources

  ! input parameter:
  !  1. imeas = 1, normalized waveform difference. Adjoint source is constructed from the data
  !                      only, with the form -d(t)/ || d(t) || 2
  !  2. imeas = 2, waveform difference, s(t) - d(t).
  !  3. imeas = 3, cross-correlation traveltime difference for a (banana-doughtnut) sensitivity ker-
  !                       nel. The measurement between data and synthetics is not used in constructing the adjoint
  !                       source.
  !  4. imeas = 4, amplitude difference for a (banana-doughtnut) sensitivity kernel. The measure-
  !                       ment between data and synthetics is not used in constructing the adjoint source.
  !  5. imeas = 5, cross-correlation traveltime difference for an event kernel. The measurement
  !                       between data and synthetics is used in constructing the adjoint source.
  !  6. imeas = 6, amplitude difference for an event kernel. The measurement between data and
  !                        synthetics is used in constructing the adjoint source.
  !  7. imeas = 7, multitaper traveltime difference for an event kernel. The measurement between
  !                       data and synthetics is used in constructing the adjoint source. See multitaper_notes.pdf.
  !  8. imeas = 8, multitaper amplitude difference for an event kernel. The measurement between
  !                       data and synthetics is used in constructing the adjoint source. See multitaper_notes.pdf.

  use ma_variables
  use ma_constants
  use ascii_rw       ! dwascii()
  use ma_sub2        ! fft(), fftinv()
  use ma_sub         ! mt_measure(), mt_adj()
  use ma_weighting

  implicit none

  character(len=150) :: datafile,synfile,synfile_phydisp,file_prefix,file_prefix0,file_prefix2,measure_file_prefix,adj_file_prefix
  integer :: num_meas, j, ios, npt1, npt2,npt3, npts, nn
  double precision, dimension(NDIM) :: data, syn, syn_phydisp, adj_syn_all, &
                        tr_adj_src, am_adj_src, recon_cc_all, syn_dtw_cc, syn_dtw_mt
  double precision :: t01, dt1, t02, dt2, t03, dt3, t0, dt, tstart, tend, tt, dtt, df
  double precision, dimension(NCHI) :: window_chi
  double precision :: fend0, fstart0, fend, fstart

  ! sac header information
  integer :: yr,jda,ho,mi
  double precision :: sec,dist,az,baz,slat,slon
  character(len=10) :: net,sta,chan_dat,chan,cmp,chan_syn
  double precision :: tshift, sigma_dt_cc, dlnA, sigma_dlnA_cc, sigma_dt, sigma_dlnA
  double precision :: all_chi, tr_chi, am_chi, cc_max, T_pmax_dat, T_pmax_syn
  !double precision :: tshift_f1f2, cc_max_f1f2
  double precision, dimension(NPT) :: dtau_w, dlnA_w, err_dt, err_dlnA, syn_dtw, data_dtw,syn_dtw_phydisp
  complex*16, dimension(NPT) :: trans_mtm
  integer :: nlen, i_left, i_pmax_dat, i_pmax_syn, i_right, i_right0, istart, &
        ipair, npairs, nwin, itmax
  logical :: use_trace
  !double precision :: trbdndw, a
  !integer :: iord, passes
  integer :: ipick_type
  double precision :: T_surfacewaves

  !********* PROGRAM STARTS HERE *********************
  ! read in MEASUREMENT.PAR (see ma_sub.f90 and write_par_file.pl)
  ! most parameters are global (see ma_variables.f90)
  call read_par_file(fstart0,fend0,tt,dtt,nn,chan)

  ! uses weights to balance love and rayleigh measurements
  ! we do a normalization of P_SV, P_SH, Love, Rayleigh with the number of measurement picks
  if( DO_WEIGHTING ) call setup_weighting(chan)

  ! input file: MEASUREMENT.WINDOWS
  open(11,file='MEASUREMENT.WINDOWS',status='old',iostat=ios)
  if (ios /= 0) stop 'Error opening input file: MEASUREMENT WINDOWS'

  read(11,*,iostat=ios) npairs
  if (ios /= 0) stop 'Error reading number of pairs of data/syn'
  print *, 'reading in the data and synthetics...'

  ! output files
  open(12,file='window_index',status='unknown',iostat=ios)
  open(13,file='window_chi',status='unknown',iostat=ios)

  nwin = 0; all_chi=0.
  do ipair = 1, npairs

    data(:) = 0.0
    syn(:)  = 0.0

    adj_syn_all(:) = 0.0
    recon_cc_all(:) = 0.0

    ! reads in file names for data and synthetics
    read(11,'(a)',iostat=ios) datafile
    if (ios /= 0) stop 'Error reading windows file'
    read(11,'(a)',iostat=ios) synfile
    if (ios /= 0) stop 'Error reading windows file'
    if (USE_PHYSICAL_DISPERSION) then
            synfile_phydisp=trim(synfile)//'.phydisp'
    endif


    ! read data and syn (in double precision)
    ! LQY: data is read last to be stored in memory for header retrieval later

    call drsac1(datafile,data,npt1,t01,dt1)
    call drsac1(synfile,syn,npt2,t02,dt2)
    if (USE_PHYSICAL_DISPERSION) then
            call drsac1(synfile_phydisp,syn_phydisp,npt3,t03,dt3)
    endif

    if (DISPLAY_DETAILS) then
       print *
       print *, 'data: ',trim(datafile)
       print *, '  min/max: ',sngl(minval(data(:))),sngl(maxval(data(:)))

       print *, 'syn:   ',trim(synfile)
       print *, '  min/max: ',sngl(minval(syn(:))),sngl(maxval(syn(:)))
    endif

    ! check if npts, dt, t0 matches for data and synthetics
    ! no interpolation is necessary at any point
    if (max(npt1,npt2) > NDIM) &
         stop 'Error: Too many number of points in data or syn'
    npts = min(npt1,npt2)

    if (abs(dt1-dt2) > TOL) stop 'Error: check if dt match'
    dt = dt1

    if (abs(t01-t02) > dt)  stop 'Check if t0 match'
    t0 = t01

    if (DISPLAY_DETAILS) print *,'  time, dt, npts :',sngl(t01), sngl(dt), npts


    ! apply bandpass filter to data and synthetics with saclib, if desired
    ! http://www.iris.washington.edu/pipermail/sac-help/2008-March/000376.html
    ! Access to the kidate, xapiir, and getfil is not simple and not
    ! supported under the current state of the SAC code base.

    if(RUN_BANDPASS) then
       call bandpass(data,npts,dt,fstart0,fend0)
       call bandpass(syn,npts,dt,fstart0,fend0)
       if (USE_PHYSICAL_DISPERSION) then
               call bandpass(syn_phydisp,npts,dt,fstart0,fend0)
       endif
    endif

    ! find out station/network/comp names,etc from synthetics
    call get_sacfile_header(trim(synfile),yr,jda,ho,mi,sec,net,sta, &
                          chan_dat,dist,az,baz,slat,slon)

    ! theoretical surface wave arrival time
    T_surfacewaves = dist / surface_vel

    ! synthetics always have the form BH_ or LH_, but the data may not (HH_, LH_, BL_, etc).
    cmp = chan_dat(3:3)
    chan_syn = trim(chan)//trim(cmp)

    ! example: OUT/PAS.CI.BHZ
    file_prefix0 = trim(sta)//'.'//trim(net)//'.'//trim(chan_syn)
    file_prefix2 = trim(OUT_DIR)//'/'//trim(file_prefix0)
    print *
    print *,  trim(file_prefix2), ' --- '

    ! note: MT measurement could revert to CC, but still keep the MT suffix
    write(adj_file_prefix,'(a,i2.2)') trim(file_prefix2)//'.iker', imeas0

    ! reads number of measurement windows
    read(11,*,iostat=ios) num_meas
    if (ios /= 0) stop 'Error reading num_meas'

    do j = 1, num_meas
      ! reads in start and end time of the measurement window
      read(11,*,iostat=ios) tstart, tend
      if (ios /= 0) stop 'Error reading tstart and tend'

      ! checks start and end times of window compared to trace lengths
      tstart = max(tstart,t0)
      tend = min(tend, t0+(npts-1)*dt)
      nlen = floor((tend-tstart)/dt) + 1  ! dummy, replaced later in mt_measure()

      ! write values to output file
      nwin = nwin + 1       ! overall window counter
      write(12,'(a3,a8,a5,a5,3i5,2f12.3)') net,sta,chan_syn,chan_dat,nwin,ipair,j,tstart,tend

      ! add taper type to file prefix: OUT/PAS.CI.BHZ.01.mtm
      write(file_prefix,'(a,i2.2)') trim(file_prefix2)//'.', j

      if (is_mtm == 1) then
        measure_file_prefix = trim(file_prefix) // '.mtm'  ! multitaper taper
      else if (is_mtm == 2) then
        measure_file_prefix = trim(file_prefix) // '.ctp'  ! cosine taper
      else
        measure_file_prefix = trim(file_prefix) // '.btp'  ! boxcar taper
      endif

      print *
      print *, ' Measurement window No.', j, ' ... '

      ! initialize the measurements
      window_chi(:) = 0.

      ! compute integrated waveform difference, normalized by duration of the record
      ! NOTE: (1) this is for the FULL record, not the windowed record
      !       (2) for comparison with waveform_chi, we include the 0.5 factor
      !       (3) we might want to include dt as an integration factor (also for waveform_chi),
      !           but the ratio (d-s)^2 / d^2 avoids the need for dt, nstep, or length of record
      window_chi(17) = 0.5 * sum( data**2 )
      window_chi(18) = 0.5 * sum( syn**2 )
      window_chi(19) = 0.5 * sum( (data-syn)**2 )
      window_chi(20) = npts*dt

      ! make measurements
      ! also compute reconstructed synthetics for CC (and MT, if specified) measurements
      call mt_measure(datafile,measure_file_prefix,data,syn,syn_phydisp,t0,dt,npts,tstart,tend,&
            istart,data_dtw,syn_dtw,syn_dtw_phydisp,nlen,tshift,sigma_dt_cc,dlnA,sigma_dlnA_cc,cc_max,syn_dtw_cc,&
            i_pmax_dat,i_pmax_syn,i_right0,trans_mtm,dtau_w,dlnA_w,sigma_dt,sigma_dlnA,syn_dtw_mt,err_dt,err_dlnA)
      i_right = i_right0
      i_left = 1  ! LQY: is it feasible that i_left is not 1? mt_adj() inherently assumes it.

      ! period of the max power of the synthetic record
      T_pmax_dat = (dt*NPT) / dble(i_pmax_dat)
      T_pmax_syn = (dt*NPT) / dble(i_pmax_syn)

      ! adjust frequency ranges for MT measurements
      ! fstart is constrained by NCYCLE_IN_WINDOW/tlen, fend constrained by i_right
      if (is_mtm == 1) then
         fstart = fstart0  ; fend = fend0
         call mt_measure_select(nlen,tshift,i_pmax_syn,dtau_w,err_dt, &
                              dt,i_left,i_right,fstart,fend,use_trace)
         print *, '     Tlong/Tshort (input) :', sngl(1/fstart0), sngl(1/fend0)
         print *, '     Tlong/Tshort (adjusted)  :', sngl(1/fstart), sngl(1/fend)
         print *, '     period of max data/syn power    :', sngl(T_pmax_dat), sngl(T_pmax_syn)

         ! if MT measurement window is rejected by mt_measure_select, then use a CC measurement
         if(.not. use_trace) then
            !stop 'Check why this MT measurement was rejected'
            print *, '   reverting from MT measurement to CC measurement...'
            imeas = imeas0 - 2
            is_mtm = 3  ! LQY: WHY not is_mtm = 2?
            call mt_measure(datafile,measure_file_prefix,data,syn,syn_phydisp,t0,dt,npts,tstart,tend,&
                  istart,data_dtw,syn_dtw,syn_dtw_phydisp,nlen,&
                  tshift,sigma_dt_cc,dlnA,sigma_dlnA_cc,cc_max,syn_dtw_cc,&
                  i_pmax_dat,i_pmax_syn,i_right,trans_mtm,dtau_w,dlnA_w,sigma_dt,sigma_dlnA,syn_dtw_mt)
         else
            print *, '     using this MTM. '
         endif
      endif

      ! check that the CC measurements are within the specified input range
      if (imeas >= 5) call cc_measure_select(tshift,dlnA,cc_max)

      ! write frequency limits to file
      if (OUTPUT_MEASUREMENT_FILES) then
        df = 1./(dt*NPT)
        open(71,file=trim(measure_file_prefix)//'.freq_limits')
        write(71,'(6f18.8)') fstart0, fend0, df, i_right0*df, fstart, fend
        close(71)
      endif

      ! compute adjoint sources and misfit function values and also the CC-reconstructed records
      if (COMPUTE_ADJOINT_SOURCE) then
         print *, '   Generating adjoint source and chi value for imeas = ', imeas

        ! banana-doughnut kernel (needs only synthetic trace)
        ! LQY: what is this section intended to do?
        ! reset imeas == 3 for adjoint sources without time shift and uncertainty scaling
        ! (pure cross-correlation adjoint source for banana-doughnuts)
        if(imeas == 5 .and. trim(datafile) == trim(synfile) ) then
           print *,'cross-correlation measurement:'
           print *,'  only synthetic file: ',trim(synfile)
           print *,'    without traveltime difference/uncertainty'
           print *
           imeas = 3
        endif

        tr_chi = 0.0 ; am_chi = 0.0    ! must be initialized
        call mt_adj(istart,data_dtw,syn_dtw,syn_dtw_phydisp,nlen,dt,tshift,dlnA,sigma_dt_cc,sigma_dlnA_cc,&
             dtau_w,dlnA_w,err_dt,err_dlnA,sigma_dt,sigma_dlnA,i_left,i_right,&
             window_chi,tr_adj_src,tr_chi,am_adj_src,am_chi)

        ! KEY: write misfit function values to file (two for each window)
        ! Here are the 20 columns of the vector window_chi
        !  1: MT-TT chi,    2: MT-dlnA chi,    3: XC-TT chi,    4: XC-dlnA chi
        !  5: MT-TT meas,   6: MT-dlnA meas,   7: XC-TT meas,   8: XC-dlnA meas
        !  9: MT-TT error, 10: MT-dlnA error, 11: XC-TT error, 12: XC-dlnA error
        ! WINDOW     : 13: data power, 14: syn power, 15: (data-syn) power, 16: window duration
        ! FULL RECORD: 17: data power, 18: syn power, 19: (data-syn) power, 20: record duration
        ! Example of a reduced file: awk '{print $2,$3,$4,$5,$6,$31,$32}' window_chi > window_chi_sub
        write(13,'(a14,a8,a3,a5,i4,i4,2e14.6,20e14.6,2e14.6,2f14.6)') &
           file_prefix0,sta,net,chan_syn,j,imeas,&
           tstart,tend,window_chi(:),tr_chi,am_chi,T_pmax_dat,T_pmax_syn
        print *, '     tr_chi = ', sngl(tr_chi), '  am_chi = ', sngl(am_chi)

        ! uses weighting to balance love / rayleigh measurements
        if( DO_WEIGHTING ) then
           ipick_type = 0
           if( tend <= T_surfacewaves ) then
              ! body wave picks
              if( cmp(1:1) == "Z" ) ipick_type = P_SV_V
              if( cmp(1:1) == "R" ) ipick_type = P_SV_R
              if( cmp(1:1) == "T" ) ipick_type = SH_T
           else
              ! surface wave picks
              if( cmp(1:1) == "Z" ) ipick_type = Rayleigh_V
              if( cmp(1:1) == "R" ) ipick_type = Rayleigh_R
              if( cmp(1:1) == "T" ) ipick_type = Love_T
           endif

          ! LQY: shouldn't chi values be changed accordingly?????
          ! No total chi value is calculated ...

          ! weights by phase types
          select case(ipick_type)
            case( P_SV_V )
              tr_adj_src(:) = tr_adj_src(:) * num_P_SV_V
              tr_chi = tr_chi * num_P_SV_V
            case( P_SV_R )
              tr_adj_src(:) = tr_adj_src(:) * num_P_SV_R
              tr_chi = tr_chi * num_P_SV_R
            case( SH_T )
              tr_adj_src(:) = tr_adj_src(:) * num_SH_T
              tr_chi = tr_chi * num_SH_T
            case( Rayleigh_V )
              tr_adj_src(:) = tr_adj_src(:) * num_Rayleigh_V
              tr_chi = tr_chi * num_Rayleigh_V
            case( Rayleigh_R )
              tr_adj_src(:) = tr_adj_src(:) * num_Rayleigh_R
              tr_chi = tr_chi * num_Rayleigh_R
            case( Love_T )
              tr_adj_src(:) = tr_adj_src(:) * num_Love_T
              tr_chi = tr_chi * num_Love_T
            case default
              stop 'error ipick_type unknown'
          end select
       endif

        ! combine adjoint sources from different measurement windows
       if (mod(imeas,2)==1) then
          adj_syn_all(:) = adj_syn_all(:) + tr_adj_src(:)   ! imeas = 1,3,5,7
          all_chi = all_chi + tr_chi
       else
          adj_syn_all(:) = adj_syn_all(:) + am_adj_src(:)   ! imeas = 2,4,6,8
          all_chi = all_chi + am_chi
       endif

        ! combine CC-reconstructed records
       if (imeas >= 7) then
          recon_cc_all(istart:istart+nlen-1) = recon_cc_all(istart:istart+nlen-1) + syn_dtw_mt(1:nlen)
       else
          recon_cc_all(istart:istart+nlen-1) = recon_cc_all(istart:istart+nlen-1) + syn_dtw_cc(1:nlen)
       endif

     endif ! COMPUTE_ADJOINT_SOURCE

      ! CHT: (re-)set to multitaper parameters, if originally specified
      if (is_mtm0 == 1) then
         imeas = imeas0
         is_mtm = is_mtm0
      endif

   enddo ! nmeas

    !----------------------------
    ! write out the adjoint source for the trace (STA.NI.CMP) by combining contribution from all wins

    if (COMPUTE_ADJOINT_SOURCE) then

    ! write out the CC-reconstructed data from synthetics
       if (OUTPUT_MEASUREMENT_FILES) &
            call dwsac1(trim(file_prefix2)//'.recon.sac',recon_cc_all,npts,t0,dt)

      ! OPTIONAL: A conservative choice is to filter the adjoint source,
      !   since higher frequencies could enter from the tapering operations.
      ! Note: time_window in mt_adj.f90 tapers the windows.

      ! note also:
      ! measurements are done on filtered synthetics F(s) and filtered data F(d), such that DeltaT
      ! is given for filtered data & synthetics.
      ! then kernels,
      ! i.e. for a traveltime measurement: DeltaT = 1/N * int  F(d/dt s) F(ds)
      ! should contain this filter as well.
      !
      ! when we construct the adjoint source here,it is initially a filtered version
      ! as well F(s_adj) since we use/depend on filtered synthetics F(s).
      ! however, for kernel simulations, we do run with a reconstructed forward wavefield,
      ! which is unfiltered (only filter there is by source half-time), but we want to convolve
      !  K = int F*(s_adj) F(s)
      ! using the same (bandpass) filter F() as used for filtereing data & synthetics in the meausurements
      ! We can write the kernel expression as K = int F*{F* (s_adj)}  s
      ! thus we should apply the filter F() twice on the adjoint source
      !
      ! why is this important? the filter, like bandpassing, is usually acausal, that is, it can
      ! introduce a slight phase-shift to the data. but, phase-shifts is what we are interested in
      ! and invert for. so, filtering might affect our inversions...

      ! we do use a bandpass filter here again on the adjoint source. this is slightly different
      ! to the transfer function filter in SAC used initially to filter data & synthetics.
      ! but this seems to be the best and fairly easy what we can do here...
      call bandpass(adj_syn_all,npts,dt,fstart0,fend0) ! sac butterworth filter

      ! cut and interpolate to match time-stepping for SEM
      ! NOTE: This can leave a non-zero value to start the record,
      !       which is NOT GOOD for the SEM simulation.
      call interpolate_syn(adj_syn_all,t0,dt,npts,tt,dtt,nn)

      ! Taper the start of the adjoint source, since cutting the record
      ! may have left a non-zero value to start the record,
      ! which is not good for the SEM simulation.
      itmax = int(TSHORT/dtt)
      call taper_start(adj_syn_all,nn,itmax)

      ! output the adjoint source (or ray density) as ASCII or SAC format
      print *, 'writing adjoint source to file for the full seismogram'
      if( DO_RAY_DENSITY_SOURCE ) then
        call dwascii(trim(adj_file_prefix)//'.density.adj',adj_syn_all,nn,tt,dtt)
      else
        call dwascii(trim(adj_file_prefix)//'.adj',adj_syn_all,nn,tt,dtt)
        ! LQY add the sum of chi values (total misfit), weights included if DO_WEIGHTING on
        open(14,file='window_chi_sum',status='unknown')
        write(14,*) all_chi
        close(14)
      endif

    endif

  enddo ! npairs

  close(11)  ! read: MEASUREMENT.WINDOWS
  close(12)  ! write: window_index
  close(13)  ! write: window_chi

end program measure_adj

