module wave2d_sub

  use wave2d_constants
  use wave2d_variables

  implicit none

contains

  !-----------------------------------------------------

  subroutine write_parameters(filename)

    character(len=*),intent(in) :: filename

    print *, 'writing out parameters'

    open(unit=12, file=filename, status='unknown', iostat=ios)

    write(12,*) 'hey you'

    close(12)

  end subroutine write_parameters

  !-----------------------------------------------------

!!$  subroutine get_source_function(nsrc,origin_time,f0,samp,ti)
!!$
!!$  ! compute the source time function for each source
!!$  ! in general, we use a point source
!!$  ! note that we do not need to know WHERE the source is located
!!$
!!$    ! input
!!$    integer, intent(in) :: nsrc
!!$    double precision, intent(in) :: f0(NCOMP)
!!$    double precision, intent(in) :: origin_time
!!$
!!$    double precision, intent(inout) :: samp(NSTEP,NCOMP,nsrc)
!!$    double precision, intent(out) :: ti(NSTEP)
!!$    integer :: i,icomp
!!$
!!$  ! fill each source index with the source time function
!!$  ! multiplied by the magnitude and direction of the source
!!$  do i = 1,nsrc
!!$    do icomp = 1,NCOMP
!!$       samp(:, icomp, i) = stf(:) * f0(icomp)
!!$    enddo
!!$  enddo
!!$
!!$  end subroutine get_source_function

  !-----------------------------------------------------

  subroutine get_source_time_function(origin_time,stf_vec,ti)

  ! compute the source time function for each source
  ! in general, we use a point source
  ! note that we do not need to know WHERE the source is located

    double precision, intent(in) :: origin_time
    double precision, intent(out) :: ti(NSTEP), stf_vec(NSTEP)

    ! source decay rate (also change in source spectrum if needed)
    double precision, parameter :: decay_rate = 2.628d0
    !double precision, external :: erf

    integer :: itime,i,cyc,icomp,nsrc_plot
    double precision :: alpha, per, t, t1, t2, amp, stf, fgaus, dgaus, tmp
    double precision :: az,kx,ky,c_source
    double precision, dimension(:), allocatable :: d_vec
    character(len=200) :: filename

!-------------------------------------

  print *
  print *, 'compute the forward source time function in get_source_time_function.f90'
  print *, 'hdur = ',sngl(hdur),' s, ISRC_TIME = ', ISRC_TIME
  print *, 'origin time is = ', sngl(origin_time),' s'

  ! parameters for source time function
  ! the non-zero Gaussian is needed for plotting the source time function (perl)
  alpha = decay_rate/hdur
  fgaus = 1.d-8                  ! fraction of amplitude at edge of Gaussian
  dgaus = sqrt(-log(fgaus)) / alpha

  if (ISRC_TIME == 1) then ! Ricker
     amp = -2.*(alpha**3)/dsqrt(PI)

  else if (ISRC_TIME == 2) then ! Gaussian
     amp = alpha/dsqrt(PI)

  else if (ISRC_TIME == 3) then ! truncated sine
     cyc = 3
     per = 2.*hdur
     !t1 = -0.50*per
     t1 = 0.
     t2 = t1 + per*dble(cyc)
     amp = alpha**2*dsqrt(2./PI)*exp(-0.5d0)

  else if (ISRC_TIME == 4) then ! sine
     per = 2.*hdur
     amp = alpha**2*dsqrt(2./PI)*exp(-0.5d0)

!!$  else if (ISRC_TIME==5) then ! plane wave field
!!$
!!$     amp = alpha**2*dsqrt(2./PI)*exp(-0.5d0)   ! amplitude
!!$     az = 25.*PI/180.                          ! azimuth of vector (from north)
!!$     kx = sin(az) ; ky = cos(az)               ! directional unit vector k
!!$     per = 2.*hdur                             ! period of wavefield (s)
!!$
!!$     ! phase velocity of source wavefield (m/s)
!!$     ! NOTE: this signal will be ALIASED if the wavelength is too short
!!$     ! this is based on the Airy dispersion for gravity waves (Bromirski+Duennebier, 5-3)
!!$     ! probably the speed should be based on the slope of the seafloor
!!$     !c_source = 9.81*per/(2*PI)  ! T=16s, c=25 m/s
!!$     c_source = c0
!!$
!!$     ! projection of each source point vector onto the directional vector k
!!$     allocate(d_vec(nsrc))
!!$     d_vec(:) = 0.
!!$     print *
!!$     print *, 'Plane wave source:'
!!$     print *, '   azimuth        : ', sngl(az*180/PI)
!!$     print *, '   period         : ', sngl(per), ' s'
!!$     print *, '   phase velocity : ', sngl(c_source/1000.), ' km/s'
!!$     print *, '   wavelength     : ', sngl(c_source*per/1000.), ' km'
!!$     print *, 'relative distance from plane wave wavefront to each source point:'
!!$     do i=1,nsrc
!!$        d_vec(i) = kx*x(sglob(i)) + ky*z(sglob(i))
!!$        write(*,'(i6,3f12.3)') i, x(sglob(i))/1000., z(sglob(i))/1000., d_vec(i)/1000.
!!$     enddo
!!$     print *
  endif

    do itime = 1, NSTEP
       ti(itime) = dble(itime-1)*DT

       t = ti(itime) - origin_time  ! time shift

       if (ISRC_TIME == 1) then
          ! d/dt[Gaussian] wavelet
          if (t >= -dgaus .and. t <= dgaus) then
             stf = amp*t*exp(-alpha*alpha*t*t)
          else
             stf = 0.
          endif

       else if (ISRC_TIME == 2) then
          ! Error function
          ! source_time_function = 0.5d0*(1.0d0+erf(decay_rate*t/hdur))

          ! Gaussian (this one causes static offset at stations)
          if (t >= -dgaus .and. t <= dgaus) then
             stf = amp*exp(-alpha*alpha*t*t)
          else
             stf = 0.
          endif

       else if (ISRC_TIME == 3) then
          ! truncated sine function (duration is cyc*per seconds)
          if (t >= t1 .and. t <= t2) then
             stf = amp*sin(2*PI*(t-t1)/per)
          else
             stf = 0.
          endif

       else if (ISRC_TIME == 4) then
          ! sine function
          stf = amp*sin(2*PI*t/per)
          !stf = amp/2.*sin(2*PI*t/per) + amp/2.*sin(2*PI*t/(1.1*per))

       !else if (ISRC_TIME==5) then
       !   ! plane wavefield, dependant on source position
       !   tmp = t - d_vec(i)/c_source
       !   !stf = amp*sin( 2*PI/per*tmp )
       !   stf = amp/2.*sin(2*PI*tmp/per) + amp/2.*sin(2*PI*tmp/(1.1*per))

       endif

       ! fill source time function
       stf_vec(itime) = stf

    enddo


  ! taper time series
  ! DO WE WANT TO SIMPLY DETREND THE TIME SERIES?
  if (SRC_TAPER) call taper_series(stf_vec(:),NSTEP)

  end subroutine get_source_time_function

  !-----------------------------------------------------

  subroutine taper_series(x,nt)

  integer, intent(in) :: nt
  double precision, intent(inout) :: x(nt)

  double precision :: ntemp,jtemp,wtemp,cfac
  integer :: i,pwr

  ntemp = dble(nt)/2.0

  ! KEY COMMAND: power of polynomial taper
  ! higher power means affecting only the ends of the series
  pwr = 10   ! Welch : pwr=2

  ! Welch taper (in time)
  do i = 1,nt

     jtemp = dble(i-1)
     wtemp = (jtemp - ntemp) / ntemp
     cfac = 1. - wtemp**pwr
     !cfac = 1 - (2*(i - 1)/(nt - 1) - 1) ** pwr  ! see Qinya code below

     x(i) = cfac*x(i)
  enddo

  end subroutine taper_series

  !-----------------------------------------------------

  subroutine write_snapshot(disp, filename)

    double precision,intent(in) :: disp(NCOMP,NGLOB)
    character(len=200),intent(in) :: filename

    integer :: icomp, iglob, ios

    open(unit = 11, file = trim(filename), status = 'unknown',iostat=ios)
    if (ios /= 0) stop 'Error writing snapshot to disk'
    do iglob = 1, NGLOB
       if (NCOMP == 3) then
          write(11,'(5e12.3)') x(iglob)/LENGTH, z(iglob)/LENGTH, &
                  sngl(disp(1,iglob)),sngl(disp(2,iglob)),sngl(disp(3,iglob))
       else
          write(11,'(5e12.3)') x(iglob)/LENGTH, z(iglob)/LENGTH, sngl(disp(1,iglob))
       endif
    enddo
    close(11)

  end subroutine write_snapshot

  !-----------------------------------------------------

!!$  subroutine write_source_function(nsrc, ti, seis, sglob, seis_name)
!!$
!!$    integer, intent(in) :: nsrc
!!$    integer, intent(in) :: sglob(nsrc)
!!$    double precision, intent(in) ::  seis(NSTEP,NCOMP,nsrc)
!!$    double precision, intent(in) ::  ti(NSTEP)
!!$    character(len=*),intent(in) :: seis_name
!!$
!!$    character(len=200) :: filename
!!$    integer :: i,icomp,itime
!!$
!!$    print *, 'writing out source time function'
!!$
!!$   ! write out source time functions for each point source
!!$   do i = 1,nsrc
!!$      do icomp = 1,NCOMP
!!$
!!$         write(filename,'(a,a,i5.5,a,i1.1)') trim(seis_name), '_', i, '_', icomp
!!$
!!$         open(12,file=filename,status='unknown')
!!$         write(*,*) 'Source #', i, ' at ', sngl(x_lon(sglob(i))), ', ', sngl(z_lat(sglob(i)))
!!$         do itime=1,NSTEP
!!$            write(12,'(f16.6,e16.6)') ti(itime), seis(itime,icomp,i)
!!$            !write(12,'(f16.6,e16.6)') ti(itime), seis(itime,icomp,nsrc)/maxval(seis(:,icomp,nsrc))
!!$         enddo
!!$         close(12)
!!$      enddo
!!$   enddo
!!$
!!$  end subroutine write_source_function

  !-----------------------------------------------------

  subroutine write_seismogram(seis, nrec, seis_name)

    integer, intent(in) :: nrec
    double precision, intent(in) ::  seis(NSTEP,NCOMP,nrec)
    character(len=*),intent(in) :: seis_name

    character(len=200) :: filename
    integer :: irec, icomp, itime

    print *, 'writing out time series (source time function or seismogram)'

    do irec = 1,nrec
       do icomp = 1,NCOMP
          write(filename,'(a,a,i5.5,a,i1.1)') trim(seis_name), '_', irec, '_', icomp
          open(unit = 12, file = filename, status = 'unknown', iostat=ios)
          if (ios /= 0) stop 'Error opening seismogram to write'

          do itime = 1,NSTEP
             write(12,*) DT*itime, seis(itime, icomp, irec)
          enddo

          close(12)
       enddo
    enddo

  end subroutine write_seismogram

  !-----------------------------------------------------

  subroutine write_spectral_map(seis, nrec, rglob, seis_name, write_spectra)

    integer, intent(in) :: nrec
    integer, intent(in) :: rglob(nrec)
    double precision, intent(in) ::  seis(NSTEP,NCOMP,nrec)
    character(len=*),intent(in) :: seis_name
    logical, intent(in) :: write_spectra

    character(len=200) :: filename,filename1,filename2
    integer :: i, irec, icomp, itime

    ! frequency domain
    integer, parameter :: FFTW_FORWARD=-1, FFTW_ESTIMATE=0
    complex*16 :: out(NOUT)
    double precision :: in(NSTEP),ti(NSTEP)
    integer :: plan
    double precision :: wmax,w,re,im,abs_val,ph_val,dw,abs_int
    double precision :: wmin_win, wmax_win

    !------------------------------------

    wmax = PI/DT        ! Nyquist frequency

    ! window for computing integrated power spectrum
    ! fmin,fmax are parameters (constants.f90)
    wmin_win = 2*PI*fmin
    wmax_win = 2*PI*fmax

    do icomp = 1,NCOMP

       ! open spectral map file
       write(filename1,'(a,a,i1.1)') trim(seis_name), '_map_', icomp
       open(unit=99, file=filename1, status='unknown', iostat=ios)
       if (ios /= 0) stop 'Error opening spectral map file to write'

       do irec = 1, nrec

          call dfftw_plan_dft_r2c_1d(plan,NSTEP,in,out,FFTW_ESTIMATE)

          ! specify input time series
          in(:) = seis(:,icomp,irec)

          if (0 == 1) then
            ! write input data to file
            write(filename,'(a,a,i5.5,a,i1.1)') trim(seis_name), '_in_', irec, '_', icomp
            open(unit=10, file=filename, status='unknown', iostat=ios)
            if (ios /= 0) stop 'Error opening seismogram to write'
            do itime = 1,NSTEP
               write(10,'(2e16.6)') DT*itime , in(itime)
            enddo
            close(10)
          endif

          call dfftw_execute(plan)

          if (write_spectra) then
             write(filename2,'(a,a,i5.5,a,i1.1)') trim(seis_name), '_', irec, '_', icomp
             open(unit=12, file=filename2, status='unknown', iostat=ios)
             if (ios /= 0) stop 'Error opening seismogram spectra to write'
          endif

          dw = wmax/dble(NOUT)
          abs_int = 0.
          do i = 1,NOUT
             !w = (i-1)/dble(NOUT) * wmax
             w = (i-1)*dw
             re = real(out(i))
             im = aimag(out(i))
             abs_val = sqrt(re*re+im*im)
             !ph_val = atan2(im,re)

             ! if within the frequency band
             if (w >= wmin_win .and. w <= wmax_win) abs_int = abs_int + abs_val

             if (write_spectra) write(12,'(2e16.6)') w, abs_val
             !if (write_spectra .and. w /= 0.) write(12,'(2e16.6)') (2*PI)/w, abs_val
          enddo
          if (write_spectra) close(12)

          if (0 == 1) then
            write(*,'(a,3f12.4)') ' T, s     (min/0/max) :', (2*PI)/wmax_win , 2*hdur        , (2*PI)/wmin_win
            write(*,'(a,3f12.4)') ' f, Hz    (min/0/max) :', wmin_win/(2*PI) , 1/(2*hdur)    , wmax_win/(2*PI)
            write(*,'(a,3f12.4)') ' w, rad/s (min/0/max) :', wmin_win        , 2*PI/(2*hdur) , wmax_win
            write(*,'(a,e24.8)')  '     integrated power :', dw*abs_int
            print *
          endif

          call dfftw_destroy_plan(plan)

          ! write the spectral amplitude to file
          write(99,'(3e16.6)') x_lon(rglob(irec)), z_lat(rglob(irec)), dw*abs_int
          !write(99,'(3e16.6)') x(rglob(irec)), z(rglob(irec)), dw*abs_int

       enddo  ! irec
       close(99)  ! spectral map file for one component
    enddo  ! icomp

  end subroutine write_spectral_map

  !-----------------------------------------------------

  subroutine filter(ti, seis, nrec)

    integer, intent(in) :: nrec
    double precision, intent(in) :: ti(NSTEP)
    double precision, intent(inout) :: seis(NSTEP,NCOMP,nrec)

    character(len=200) :: filename
    double precision :: data(NSTEP)
    double precision :: dt
    integer irec,itime,icomp,npts

    npts = NSTEP
    dt = ti(2) - ti(1)

    do irec = 1, nrec
       do icomp = 1,NCOMP
          data(:) = 0.
          data(:) = seis(:,icomp,irec)

!!$          write(filename,'(a,i5.5,a,i1.1)') trim(out_dir)//'bpass_in_', irec, '_', icomp
!!$          open(unit=11, file=filename, status='unknown')
!!$          write(11,'(2e16.6)') (ti(itime), data(itime), itime=1,NSTEP)
!!$          close(11)

          ! bandpass filter the records (see constants.f90 for parameter values)
          call rmean(data,npts)
          call rtrend(data,npts)
          call xapiir(data,npts,'BU',trbdndw,a,iord,'BP',fmin, fmax, dt, passes)

!!$          write(filename,'(a,i5.5,a,i1.1)') trim(out_dir)//'bpass_out_', irec, '_', icomp
!!$          open(unit=12, file=filename, status='unknown')
!!$          write(12,'(2e16.6)') (ti(itime), data(itime), itime=1,NSTEP)
!!$          close(12)

          ! replace the unfiltered record by the filtered record
          seis(:,icomp,irec) = data(:)

       enddo
    enddo

  end subroutine filter

  !----------------------------------------------------

  subroutine make_adjoint_source(nrec, syn, tstart, tend, adj_syn, data)

  integer, intent(in) :: nrec
  double precision, dimension(NSTEP,NCOMP,nrec),intent(in) :: syn
  double precision, dimension(nrec), intent(in) :: tstart, tend
  double precision, dimension(NSTEP,NCOMP,nrec),intent(out) :: adj_syn

  double precision, dimension(NSTEP,NCOMP,nrec),intent(in),optional :: data

  integer itime, icomp, istart, iend, nstart, nend, i, j, irec
  double precision, dimension(NSTEP) :: time_window
  double precision, dimension(NSTEP,NCOMP,nrec) :: syn_veloc, syn_accel
  double precision :: norm_mat(nrec,NCOMP)
  double precision :: norm, junk, ntemp

  !----------------------------------------

  adj_syn(:,:,:) = 0.
  norm_mat(:,:) = 0.

  do irec = 1,nrec

    ! time window function

    istart = max(floor(tstart(irec) / DT), 1)
    iend = min(floor(tend(irec) / DT), NSTEP)
    if (istart >= iend) stop 'Check if istart < iend'

    write(*,'(a,i6,a,f12.6,a,f12.6)') ' Receiver ', irec, ' : time window is from ', tstart(irec), ' to', tend(irec)
    write(*,*) '       index ', istart, ' to ', iend

    time_window(:) = 0.

    ! welch taper
    ntemp = (iend - istart + 1)/2.
    do i = istart, iend
      !time_window(i) = 1 - (2 * (i - istart)/(iend - istart) - 1) ** 2

      time_window(i) = 1 - ( (i - istart + 1 - ntemp) / ntemp  ) ** 2
    enddo

    !---------------------------

    ! calculate velocity and acceleration from syn (traveltime adjoint source only)
    if (IKER >= 1) then
       do itime = 2, NSTEP-1
          syn_veloc(itime,:,irec) =  (syn(itime+1,:,irec) - syn(itime-1,:,irec)) / (2 * DT)
       enddo
       syn_veloc(1,:,irec) = (syn(2,:,irec) - syn(1,:,irec)) / DT
       syn_veloc(NSTEP,:,irec) = (syn(NSTEP,:,irec) - syn(NSTEP-1,:,irec)) /DT

       do itime = 2, NSTEP-1
          syn_accel(itime,:,irec) =  (syn_veloc(itime+1,:,irec) - syn_veloc(itime-1,:,irec)) / (2 * DT)
       enddo
       syn_accel(1,:,irec) = (syn_veloc(2,:,irec) - syn_veloc(1,:,irec)) / DT
       syn_accel(NSTEP,:,irec) = (syn_veloc(NSTEP,:,irec) - syn_veloc(NSTEP-1,:,irec)) /DT
    endif

    ! assign adjoint force

    do i = 1,NCOMP

       if (IKER == 0) then       ! waveform

          adj_syn(:,i,irec) = ( syn(:,i,irec) -  data(:,i,irec) ) * time_window(:)

       else if (IKER == 5) then   ! traveltime

          ! minus sign is shifted from norm to adj_syn, in comparison with Tromp et al (2005)
          ! thus, norm is ensured to be POSITIVE (N > 0)
          norm = -DT * sum( time_window(:) * syn(:,i,irec) * syn_accel(:,i,irec) )
          if (abs(norm) > EPS) adj_syn(:,i,irec) = -syn_veloc(:,i,irec) * time_window(:) / norm

       else if (IKER == 6) then  ! amplitude

          ! norm is ensured to be POSITIVE (M > 0)
          norm = DT * sum( time_window(:) * syn(:,i,irec) * syn(:,i,irec) )
          if (abs(norm) > EPS) adj_syn(:,i,irec) = syn(:,i,irec) * time_window(:) / norm

       endif

       norm_mat(irec,i) = norm  ! could be either nothing or traveltime or amplitude

       !adj_syn(:,i,irec) = time_window(:) ! test
       !adj_syn(:,i,irec) = syn(:,i,irec) ! test
    enddo

  enddo ! loop over all receivers

  ! write normalizations to file
  open(15,file=trim(out_dir)//'Nbanana.dat',status='unknown')
  do irec = 1,nrec
     do i = 1,NCOMP
        write(15,'(2i12,1e18.6)') irec, i, norm_mat(irec,i)
     enddo
  enddo
  close(15)

  end subroutine make_adjoint_source

  !-----------------------------------------------------

!!$  subroutine waveform_adjoint_source(nrec, syn, data, tstart, tend, adj_syn)
!!$
!!$    integer, intent(in) :: nrec
!!$    double precision, dimension(NSTEP,NCOMP,nrec),intent(in) :: syn, data
!!$    double precision, dimension(nrec), intent(in) :: tstart, tend
!!$    double precision, dimension(NSTEP,NCOMP,nrec),intent(out) :: adj_syn
!!$
!!$    integer itime, icomp, istart, iend, nstart, nend, i, j, irec
!!$    double precision, dimension(NSTEP) :: time_window
!!$    double precision :: ntemp
!!$
!!$    !----------------------------------------
!!$
!!$    adj_syn(:,:,:) = 0.
!!$
!!$    do irec = 1,nrec
!!$
!!$       ! time window function
!!$
!!$       istart = max(floor(tstart(irec) / DT), 1)
!!$       iend = min(floor(tend(irec) / DT), NSTEP)
!!$       if (istart >= iend) stop 'Check if istart < iend'
!!$
!!$       write(*,'(a,i6,a,f12.6,a,f12.6)') ' Receiver ', irec, ' : time window is from ', tstart(irec), ' to', tend(irec)
!!$       write(*,*) '       index ', istart, ' to ', iend
!!$
!!$       time_window(:) = 0.
!!$
!!$       ! welch taper
!!$       ntemp = (iend - istart + 1)/2.
!!$       do i = istart, iend
!!$         time_window(i) = 1 - ( (i - istart + 1 - ntemp) / ntemp  ) ** 2
!!$       enddo
!!$
!!$      ! assign adjoint force
!!$      do i = 1,NCOMP
!!$
!!$         adj_syn(:,i,irec) = ( syn(:,i,irec) -  data(:,i,irec) ) * time_window(:)
!!$
!!$      enddo
!!$
!!$    enddo ! loop over all receivers
!!$
!!$  end subroutine waveform_adjoint_source

end module wave2d_sub
