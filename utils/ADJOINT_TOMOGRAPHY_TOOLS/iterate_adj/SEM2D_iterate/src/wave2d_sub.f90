module wave2d_sub

  use wave2d_constants
  use wave2d_variables

  implicit none

contains

  !-----------------------------------------------------

  subroutine write_parameters_plot(filename)

    character(len=*),intent(in) :: filename

    print *, 'writing out parameters to ' // trim(filename)

    open(unit=12, file=trim(filename), status='unknown')

    write(12,'(3i10)')   0,NSTEP,NSAVE
    write(12,'(3i10)')   FOR_X, FOR_Y, FOR_Z
    write(12,'(1f16.8)') DT
    write(12,'(1f16.8)') hdur
    close(12)

  end subroutine write_parameters_plot

  !-----------------------------------------------------

  subroutine write_parameters(filename)

    character(len=*),intent(in) :: filename

    print *, 'writing out parameters to ' // trim(filename)

    open(unit=12, file=trim(filename), status='unknown')

    write(12,*) 'IRUNZ',IRUNZ
    write(12,*) 'NFRAME',NFRAME
    write(12,*) 'NSAVE',NSAVE
    write(12,*) 'NSTEP',NSTEP
    write(12,*) 'DT',DT
    write(12,*) 'ISRC_TIME',ISRC_TIME
    write(12,*) 'hdur',hdur
    write(12,*) 'tshift',tshift
    if (SRC_TAPER) write(12,*) 'SRC_TAPER 1'
    if (.not. SRC_TAPER) write(12,*) 'SRC_TAPER 0'
    write(12,*) 'FNORM',FNORM
    write(12,*) 'FOR_X',FOR_X
    write(12,*) 'FOR_Y',FOR_Y
    write(12,*) 'FOR_Z',FOR_Z
    write(12,*) 'REV_X',REV_X
    write(12,*) 'REV_Y',REV_Y
    write(12,*) 'REV_Z',REV_Z
    write(12,*) 'ISRC_SPACE',ISRC_SPACE
    write(12,*) 'IREC_SPACE',IREC_SPACE
    write(12,*) 'NMESH_REC',NMESH_REC
    write(12,*) 'SOURCE_GRID_BUFFER',SOURCE_GRID_BUFFER
    write(12,*) 'STATION_GRID_BUFFER',STATION_GRID_BUFFER
    write(12,*) 'STATION_COAST_BUFFER',STATION_COAST_BUFFER
    write(12,*) 'LAT_MIN',LAT_MIN
    write(12,*) 'LON_MIN',LON_MIN
    write(12,*) 'UTM_PROJECTION_ZONE',UTM_PROJECTION_ZONE
    write(12,*) 'LENGTH',LENGTH
    write(12,*) 'HEIGHT ',HEIGHT
    write(12,*) 'AREA',AREA
    write(12,*) 'NEX',NEX
    write(12,*) 'NEZ',NEZ
    write(12,*) 'R_BETA_OVER_ALPHA',R_BETA_OVER_ALPHA
    write(12,*) 'PBETA',PBETA
    write(12,*) 'PALPHA',PALPHA
    write(12,*) 'PRHO',PRHO
    write(12,*) 'IMODEL_SYN',IMODEL_SYN
    write(12,*) 'IMODEL_DAT',IMODEL_DAT
    if (M0ISMPRIOR) write(12,*) 'M0ISMPRIOR 1'
    if (.not. M0ISMPRIOR) write(12,*) 'M0ISMPRIOR 0'
    write(12,*) 'ISMOOTH_EVENT_KERNEL',ISMOOTH_EVENT_KERNEL
    write(12,*) 'ISMOOTH_MISFIT_KERNEL',ISMOOTH_MISFIT_KERNEL
    write(12,*) 'ISMOOTH_INITIAL_MODEL',ISMOOTH_INITIAL_MODEL
    write(12,*) 'ISMOOTH_MODEL_UPDATE',ISMOOTH_MODEL_UPDATE
    write(12,*) 'SIGMA_SMOOTH_KERNEL',SIGMA_SMOOTH_KERNEL
    write(12,*) 'SIGMA_SMOOTH_MODEL',SIGMA_SMOOTH_MODEL
    write(12,*) 'GAMMA_SMOOTH_KERNEL',GAMMA_SMOOTH_KERNEL
    write(12,*) 'GAMMA_SMOOTH_MODEL',GAMMA_SMOOTH_MODEL
    if (HIGH_RES_SMOOTHING) write(12,*) 'HIGH_RES_SMOOTHING 1'
    if (.not. HIGH_RES_SMOOTHING) write(12,*) 'HIGH_RES_SMOOTHING 0'
    if (EXAMPLE_Gaussian) write(12,*) 'EXAMPLE_Gaussian 1'
    if (.not. EXAMPLE_Gaussian) write(12,*) 'EXAMPLE_Gaussian 0'
    write(12,*) 'IKER',IKER
    write(12,*) 'IAMP_VEL',IAMP_VEL
    write(12,*) 'ISURFACE',ISURFACE
    write(12,*) 'NCOMP',NCOMP
    write(12,*) 'NABSORB',NABSORB
    if (WRITE_STF_F) write(12,*) 'WRITE_STF_F 1'
    if (.not. WRITE_STF_F) write(12,*) 'WRITE_STF_F 0'
    if (WRITE_SEISMO_F) write(12,*) 'WRITE_SEISMO_F 1'
    if (.not. WRITE_SEISMO_F) write(12,*) 'WRITE_SEISMO_F 0'
    if (WRITE_SEISMO_RECONSTRUCT) write(12,*) 'WRITE_SEISMO_RECONSTRUCT 1'
    if (.not. WRITE_SEISMO_RECONSTRUCT) write(12,*) 'WRITE_SEISMO_RECONSTRUCT 0'
    if (WRITE_STF_A) write(12,*) 'WRITE_STF_A 1'
    if (.not. WRITE_STF_A) write(12,*) 'WRITE_STF_A 0'
    if (WRITE_SEISMO_A) write(12,*) 'WRITE_SEISMO_A 1'
    if (.not. WRITE_SEISMO_A) write(12,*) 'WRITE_SEISMO_A 0'
    if (WRITE_KERNELS) write(12,*) 'WRITE_KERNELS 1'
    if (.not. WRITE_KERNELS) write(12,*) 'WRITE_KERNELS 0'
    if (WRITE_KERNEL_SNAPSHOTS) write(12,*) 'WRITE_KERNEL_SNAPSHOTS 1'
    if (.not. WRITE_KERNEL_SNAPSHOTS) write(12,*) 'WRITE_KERNEL_SNAPSHOTS 0'
    if (WRITE_WAVFIELD_SNAPSHOTS) write(12,*) 'WRITE_WAVFIELD_SNAPSHOTS 1'
    if (.not. WRITE_WAVFIELD_SNAPSHOTS) write(12,*) 'WRITE_WAVFIELD_SNAPSHOTS 0'
    if (COMPUTE_KERNELS) write(12,*) 'COMPUTE_KERNELS 1'
    if (.not. COMPUTE_KERNELS) write(12,*) 'COMPUTE_KERNELS 0'
    if (READ_IN) write(12,*) 'READ_IN 1'
    if (.not. READ_IN) write(12,*) 'READ_IN 0'
    if (READ_SINGLE) write(12,*) 'READ_SINGLE 1'
    if (.not. READ_SINGLE) write(12,*) 'READ_SINGLE 0'
    write(12,*) 'NITERATION',NITERATION
    write(12,*) 'VAR_RED_MIN',VAR_RED_MIN
    write(12,*) 'SIGMA_DT',SIGMA_DT
    write(12,*) 'SIGMA_DLNA',SIGMA_DLNA
    write(12,*) 'SIGMA_WAVEFORM',SIGMA_WAVEFORM
    if (ADD_DATA_ERRORS) write(12,*) 'ADD_DATA_ERRORS 1'
    if (.not. ADD_DATA_ERRORS) write(12,*) 'ADD_DATA_ERRORS 0'
    write(12,*) 'POLY_ORDER',POLY_ORDER
    write(12,*) 'PERT_STRUCT_BETA',PERT_STRUCT_BETA
    write(12,*) 'PERT_SOURCE_T',PERT_SOURCE_T
    write(12,*) 'PERT_SOURCE_X',PERT_SOURCE_X
    write(12,*) 'INV_STRUCT_BETA',INV_STRUCT_BETA
    write(12,*) 'INV_SOURCE_T',INV_SOURCE_T
    write(12,*) 'INV_SOURCE_X',INV_SOURCE_X
    if (INCLUDE_MODEL_NORM) write(12,*) 'INCLUDE_MODEL_NORM 1'
    if (.not. INCLUDE_MODEL_NORM) write(12,*) 'INCLUDE_MODEL_NORM 0'
    if (ISOURCE_LOG) write(12,*) 'ISOURCE_LOG 1'
    if (.not. ISOURCE_LOG) write(12,*) 'ISOURCE_LOG 0'
    write(12,*) 'NVAR_STRUCT',NVAR_STRUCT
    write(12,*) 'NVAR_SOURCE',NVAR_SOURCE
    write(12,*) 'NVAR',NVAR
    write(12,*) 'STRUCTURE_PARAMETER_TYPE',STRUCTURE_PARAMETER_TYPE
    write(12,*) 'DENSITY',DENSITY
    write(12,*) 'INCOMPRESSIBILITY',INCOMPRESSIBILITY
    write(12,*) 'RIGIDITY',RIGIDITY
    write(12,*) 'PWAVESPEED',PWAVESPEED
    write(12,*) 'SWAVESPEED',SWAVESPEED
    write(12,*) 'BWAVESPEED',BWAVESPEED
    write(12,*) 'HWIN1',HWIN1
    write(12,*) 'HWIN2',HWIN2
    if (SUPPRESS_UTM_PROJECTION) write(12,*) 'SUPPRESS_UTM_PROJECTION 1'
    if (.not. SUPPRESS_UTM_PROJECTION) write(12,*) 'SUPPRESS_UTM_PROJECTION 0'
    write(12,*) 'ILONGLAT2UTM',ILONGLAT2UTM
    write(12,*) 'IUTM2LONGLAT',IUTM2LONGLAT
    write(12,*) 'ILONLAT2MESH',ILONLAT2MESH
    write(12,*) 'IMESH2LONLAT',IMESH2LONLAT
    write(12,*) 'MAX_SR_FAKE',MAX_SR_FAKE
    write(12,*) 'MAX_EVENT',MAX_EVENT
    write(12,*) 'MAX_SR',MAX_SR
    write(12,*) 'MAX_PHASE',MAX_PHASE
    write(12,*) 'MAX_COMP',MAX_COMP
    write(12,*) 'NELE',NELE
    write(12,*) 'NSPEC',NSPEC
    write(12,*) 'NGLLX',NGLLX
    write(12,*) 'NGLLZ',NGLLZ
    write(12,*) 'NGLL',NGLL
    write(12,*) 'NGLLSQUARE',NGLLSQUARE
    write(12,*) 'NGLOB',NGLOB
    write(12,*) 'NLOCAL',NLOCAL
    write(12,*) 'NSPEC_CORNER',NSPEC_CORNER
    write(12,*) 'NGNOD',NGNOD
    write(12,*) 'NGNOD2D',NGNOD2D
    write(12,*) 'NUM_ITER',NUM_ITER
    write(12,*) 'HUGEVAL',HUGEVAL
    write(12,*) 'TINYVAL',TINYVAL
    write(12,*) 'GAUSSALPHA',GAUSSALPHA
    write(12,*) 'GAUSSBETA',GAUSSBETA
    write(12,*) 'PI',PI
    write(12,*) 'FOUR_THIRDS',FOUR_THIRDS
    write(12,*) 'ONE_THIRD',ONE_THIRD
    write(12,*) 'ONEOVERTWO',ONEOVERTWO
    write(12,*) 'DEG',DEG
    write(12,*) 'NDIM',NDIM
    close(12)

  end subroutine write_parameters

  !-----------------------------------------------------

  subroutine write_chi(dir,nevent,nrec)

    ! Following the simulations for the event kernel, we have chi_data(nevent,nrec,NCOMP,nsrc).
    ! Here we write out various versions of this.

    integer, intent(in) :: nevent, nrec
    integer :: ievent, irec, icomp
    character(len=200) :: dir, filename

    print *, 'writing out chi_data(nevent,nrec,NCOMP,nsrc) to various files...'
    print *, '  nevent = ', nevent
    print *, '    nrec = ', nrec
    print *, '   NCOMP = ', NCOMP
    print *, '    nsrc = ', 1

    ! all chi_data values
    open(19,file=trim(dir)//'chi_data_all.dat',status='unknown')
    do ievent = 1,nevent
       do irec = 1,nrec
          do icomp = 1,NCOMP
             write(19,'(3i8,1e20.10)') ievent, irec, icomp, chi_data(ievent,irec,icomp,1)
          enddo
       enddo
    enddo
    close(19)

    ! summed chi_data value for each event (nevent by 1)
    open(19,file=trim(dir)//'summed_chi_e.dat',status='unknown')
    do ievent = 1,nevent
       write(19,'(1f20.10)') sum(chi_data(ievent,:,:,1))
    enddo
    close(19)

    ! summed chi_data value for each receiver (nrec by 1)
    open(19,file=trim(dir)//'summed_chi_r.dat',status='unknown')
    do irec = 1,nrec
       write(19,'(1f20.10)') sum(chi_data(:,irec,:,1))
    enddo
    close(19)

    ! summed chi_data value for each component (NCOMP by 1)
    open(19,file=trim(dir)//'summed_chi_c.dat',status='unknown')
    do icomp = 1,NCOMP
       write(19,'(1f20.10)') sum(chi_data(:,:,icomp,1))
    enddo
    close(19)

    ! WRITE OUT GLOBAL VARIABLES

!!$    open(19,file=trim(dir)//'chi_model_norm_target.dat',status='unknown')
!!$    write(19,'(1f20.10)') chi_model_norm_target
!!$    close(19)
!!$
!!$    open(19,file=trim(dir)//'chi_model_stop.dat',status='unknown')
!!$    write(19,'(1f20.10)') chi_model_stop
!!$    close(19)

    open(19,file=trim(dir)//'chi_data_stop.dat',status='unknown')
    write(19,'(1f20.10)') chi_data_stop
    close(19)

    open(19,file=trim(dir)//'model_norm.dat',status='unknown')
    write(19,'(1f20.10)') model_norm
    close(19)

    open(19,file=trim(dir)//'data_norm.dat',status='unknown')
    write(19,'(1f20.10)') data_norm
    close(19)

    open(19,file=trim(dir)//'chi.dat',status='unknown')
    write(19,'(1f20.10)') chi_val
    close(19)

    print *, ' Done writing chi values to files'

  end subroutine write_chi

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
    double precision, parameter :: decay_rate = 2.628
    !double precision, external :: erf

    integer :: itime,i,icomp,nsrc_plot
    double precision :: alpha, per, t, t1, t2, amp, stf, fgaus, dgaus, tmp, cyc
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
  fgaus = 1.0d-8                  ! fraction of amplitude at edge of Gaussian
  dgaus = sqrt(-log(fgaus)) / alpha

  if (ISRC_TIME == 1) then ! Ricker
     amp = -2.0*(alpha**3)/dsqrt(PI)

  else if (ISRC_TIME == 2) then ! Gaussian
     amp = alpha/dsqrt(PI)

  else if (ISRC_TIME == 3) then ! truncated sine
     cyc = 3.0
     per = 2.*hdur
     !t1 = -0.50*per
     t1 = 0.0
     t2 = t1 + per*cyc
     amp = alpha**2.0*dsqrt(2.0/PI)*exp(-0.5)

  else if (ISRC_TIME == 4) then ! sine
     per = 2.0*hdur
     amp = alpha**2.0*dsqrt(2.0/PI)*exp(-0.5)

!!$  else if (ISRC_TIME==5) then ! plane wave field
!!$
!!$     amp = alpha**2*dsqrt(2./PI)*exp(-0.5)   ! amplitude
!!$     az = 25.*PI/180.0                          ! azimuth of vector (from north)
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
!!$     d_vec(:) = 0.0
!!$     print *
!!$     print *, 'Plane wave source:'
!!$     print *, '   azimuth        : ', sngl(az*180/PI)
!!$     print *, '   period         : ', sngl(per), ' s'
!!$     print *, '   phase velocity : ', sngl(c_source/1000.0), ' km/s'
!!$     print *, '   wavelength     : ', sngl(c_source*per/1000.0), ' km'
!!$     print *, 'relative distance from plane wave wavefront to each source point:'
!!$     do i=1,nsrc
!!$        d_vec(i) = kx*x(sglob(i)) + ky*z(sglob(i))
!!$        write(*,'(i6,3f12.3)') i, x(sglob(i))/1000.0, z(sglob(i))/1000.0, d_vec(i)/1000.0
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
             stf = 0.0
          endif

       else if (ISRC_TIME == 2) then
          ! Error function
          ! source_time_function = 0.5*(1.0+erf(decay_rate*t/hdur))

          ! Gaussian (this one causes static offset at stations)
          if (t >= -dgaus .and. t <= dgaus) then
             stf = amp*exp(-alpha*alpha*t*t)
          else
             stf = 0.0
          endif

       else if (ISRC_TIME == 3) then
          ! truncated sine function (duration is cyc*per seconds)
          if (t >= t1 .and. t <= t2) then
             stf = amp*sin(2.0*PI*(t-t1)/per)
          else
             stf = 0.0
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
       ! DEBUG ARRAY SIZE
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
          call write_time_series(seis(:,icomp,irec), NSTEP, filename)
!!$          open(unit = 12, file = filename, status = 'unknown', iostat=ios)
!!$          if (ios /= 0) stop 'Error opening seismogram to write'
!!$          do itime = 1,NSTEP
!!$             write(12,*) DT*itime, seis(itime, icomp, irec)
!!$          enddo
!!$          close(12)

       enddo
    enddo

  end subroutine write_seismogram

  !-----------------------------------------------------

  subroutine write_time_series(fvec, tlen, filename)

    integer, intent(in) :: tlen
    double precision, intent(in) ::  fvec(tlen)
    character(len=*),intent(in) :: filename
    integer :: itime

    open(unit=12, file=filename, status='unknown', iostat=ios)
    if (ios /= 0) stop 'Error opening time series file to write'
    do itime = 1,tlen
       write(12,*) DT*itime, fvec(itime)
    enddo
    close(12)

  end subroutine write_time_series

  !-----------------------------------------------------

!!$  subroutine write_spectral_map(seis, nrec, rglob, seis_name, write_spectra)
!!$
!!$    integer, intent(in) :: nrec
!!$    integer, intent(in) :: rglob(nrec)
!!$    double precision, intent(in) ::  seis(NSTEP,NCOMP,nrec)
!!$    character(len=*),intent(in) :: seis_name
!!$    logical, intent(in) :: write_spectra
!!$
!!$    character(len=200) :: filename,filename1,filename2
!!$    integer :: i, irec, icomp, itime
!!$
!!$    ! frequency domain
!!$    integer, parameter :: FFTW_FORWARD=-1, FFTW_ESTIMATE=0
!!$    complex*16 :: out(NOUT)
!!$    double precision :: in(NSTEP),ti(NSTEP)
!!$    integer :: plan
!!$    double precision :: wmax,w,re,im,abs_val,ph_val,dw,abs_int
!!$    double precision :: wmin_win, wmax_win
!!$
!!$    !------------------------------------
!!$
!!$    wmax = PI/DT        ! Nyquist frequency
!!$
!!$    ! window for computing integrated power spectrum
!!$    ! fmin,fmax are parameters (constants.f90)
!!$    wmin_win = 2*PI*fmin
!!$    wmax_win = 2*PI*fmax
!!$
!!$    do icomp = 1,NCOMP
!!$
!!$       ! open spectral map file
!!$       write(filename1,'(a,a,i1.1)') trim(seis_name), '_map_', icomp
!!$       open(unit=99, file=filename1, status='unknown', iostat=ios)
!!$       if (ios /= 0) stop 'Error opening spectral map file to write'
!!$
!!$       do irec = 1, nrec
!!$
!!$          ! KEY: initialize Fourier transform
!!$          call dfftw_plan_dft_r2c_1d(plan,NSTEP,in,out,FFTW_ESTIMATE)
!!$
!!$          ! specify input time series
!!$          in(:) = seis(:,icomp,irec)
!!$
!!$          if (0==1) then
!!$            ! write input data to file
!!$            write(filename,'(a,a,i5.5,a,i1.1)') trim(seis_name), '_in_', irec, '_', icomp
!!$            open(unit=10, file=filename, status='unknown', iostat=ios)
!!$            if (ios /= 0) stop 'Error opening seismogram to write'
!!$            do itime = 1,NSTEP
!!$               write(10,'(2e16.6)') DT*itime , in(itime)
!!$            enddo
!!$            close(10)
!!$          endif
!!$
!!$          ! KEY: Fourier transform
!!$          call dfftw_execute(plan)
!!$
!!$          if (write_spectra) then
!!$             write(filename2,'(a,a,i5.5,a,i1.1)') trim(seis_name), '_', irec, '_', icomp
!!$             open(unit=12, file=filename2, status='unknown', iostat=ios)
!!$             if (ios /= 0) stop 'Error opening seismogram spectra to write'
!!$          endif
!!$
!!$          dw = wmax/dble(NOUT)
!!$          abs_int = 0.0
!!$          do i = 1,NOUT
!!$             !w = (i-1)/dble(NOUT) * wmax
!!$             w = (i-1)*dw
!!$             re = real(out(i))
!!$             im = aimag(out(i))
!!$             abs_val = sqrt(re*re+im*im)
!!$             !ph_val = atan2(im,re)
!!$
!!$             ! if within the frequency band
!!$             if (w >= wmin_win .and. w <= wmax_win) abs_int = abs_int + abs_val
!!$
!!$             if (write_spectra) write(12,'(2e16.6)') w, abs_val
!!$             !if (write_spectra .and. w /= 0.0) write(12,'(2e16.6)') (2*PI)/w, abs_val
!!$          enddo
!!$          if (write_spectra) close(12)
!!$
!!$          if (0==1) then
!!$            write(*,'(a,3f12.4)') ' T, s     (min/0/max) :', (2*PI)/wmax_win , 2*hdur        , (2*PI)/wmin_win
!!$            write(*,'(a,3f12.4)') ' f, Hz    (min/0/max) :', wmin_win/(2*PI) , 1/(2*hdur)    , wmax_win/(2*PI)
!!$            write(*,'(a,3f12.4)') ' w, rad/s (min/0/max) :', wmin_win        , 2*PI/(2*hdur) , wmax_win
!!$            write(*,'(a,e24.8)')  '     integrated power :', dw*abs_int
!!$            print *
!!$          endif
!!$
!!$          call dfftw_destroy_plan(plan)
!!$
!!$          ! write the spectral amplitude to file
!!$          write(99,'(3e16.6)') x_lon(rglob(irec)), z_lat(rglob(irec)), dw*abs_int
!!$          !write(99,'(3e16.6)') x(rglob(irec)), z(rglob(irec)), dw*abs_int
!!$
!!$       enddo  ! irec
!!$       close(99)  ! spectral map file for one component
!!$    enddo  ! icomp
!!$
!!$  end subroutine write_spectral_map

  !-----------------------------------------------------

!!$  subroutine filter(ti, seis, nrec)
!!$
!!$    integer, intent(in) :: nrec
!!$    double precision, intent(in) :: ti(NSTEP)
!!$    double precision, intent(inout) :: seis(NSTEP,NCOMP,nrec)
!!$
!!$    character(len=200) :: filename
!!$    double precision :: data(NSTEP)
!!$    double precision :: dt
!!$    integer irec,itime,icomp,npts
!!$
!!$    npts = NSTEP
!!$    dt = ti(2) - ti(1)
!!$
!!$    do irec = 1, nrec
!!$       do icomp = 1,NCOMP
!!$          data(:) = 0.0
!!$          data(:) = seis(:,icomp,irec)
!!$
!!$          !write(filename,'(a,i5.5,a,i1.1)') trim(out_dir)//'bpass_in_', irec, '_', icomp
!!$          !open(unit=11, file=filename, status='unknown')
!!$          !write(11,'(2e16.6)') (ti(itime), data(itime), itime=1,NSTEP)
!!$          !close(11)
!!$
!!$          ! bandpass filter the records (see constants.f90 for parameter values)
!!$          call rmean(data,npts)
!!$          call rtrend(data,npts)
!!$          call xapiir(data,npts,'BU',trbdndw,a,iord,'BP',fmin, fmax, dt, passes)
!!$
!!$          !write(filename,'(a,i5.5,a,i1.1)') trim(out_dir)//'bpass_out_', irec, '_', icomp
!!$          !open(unit=12, file=filename, status='unknown')
!!$          !write(12,'(2e16.6)') (ti(itime), data(itime), itime=1,NSTEP)
!!$          !close(12)
!!$
!!$          ! replace the unfiltered record by the filtered record
!!$          seis(:,icomp,irec) = data(:)
!!$
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine filter

  !----------------------------------------------------

!!$  subroutine make_adjoint_source(nrec, syn, tstart, tend, adj_syn, data)
!!$
!!$  integer, intent(in) :: nrec
!!$  double precision, dimension(NSTEP,NCOMP,nrec),intent(in) :: syn
!!$  double precision, dimension(nrec), intent(in) :: tstart, tend
!!$  double precision, dimension(NSTEP,NCOMP,nrec),intent(out) :: adj_syn
!!$
!!$  double precision, dimension(NSTEP,NCOMP,nrec),intent(in),optional :: data
!!$
!!$  integer itime, icomp, istart, iend, nstart, nend, i, j, irec
!!$  double precision, dimension(NSTEP) :: time_window
!!$  double precision, dimension(NSTEP,NCOMP,nrec) :: syn_veloc, syn_accel
!!$  double precision :: norm_mat(nrec,NCOMP)
!!$  double precision :: norm, junk, ntemp
!!$
!!$  !----------------------------------------
!!$
!!$  adj_syn(:,:,:) = 0.0
!!$  norm_mat(:,:) = 0.0
!!$
!!$  do irec = 1,nrec
!!$
!!$    ! time window function
!!$
!!$    istart = max(floor(tstart(irec) / DT), 1)
!!$    iend = min(floor(tend(irec) / DT), NSTEP)
!!$    if (istart >= iend) stop 'Check if istart < iend'
!!$
!!$    write(*,'(a,i6,a,f12.6,a,f12.6)') ' Receiver ', irec, ' : time window is from ', tstart(irec), ' to', tend(irec)
!!$    write(*,*) '       index ', istart, ' to ', iend
!!$
!!$    time_window(:) = 0.0
!!$
!!$    ! welch taper
!!$    ntemp = (iend - istart + 1)/2.
!!$    do i = istart, iend
!!$      !time_window(i) = 1 - (2 * (i - istart)/(iend - istart) - 1) ** 2
!!$
!!$      time_window(i) = 1 - ( (i - istart + 1 - ntemp) / ntemp  ) ** 2
!!$    enddo
!!$
!!$    !---------------------------
!!$
!!$    ! calculate velocity and acceleration from syn (traveltime adjoint source only)
!!$    if (IKER >= 1) then
!!$       do itime = 2, NSTEP-1
!!$          syn_veloc(itime,:,irec) =  (syn(itime+1,:,irec) - syn(itime-1,:,irec)) / (2 * DT)
!!$       enddo
!!$       syn_veloc(1,:,irec) = (syn(2,:,irec) - syn(1,:,irec)) / DT
!!$       syn_veloc(NSTEP,:,irec) = (syn(NSTEP,:,irec) - syn(NSTEP-1,:,irec)) /DT
!!$
!!$       do itime = 2, NSTEP-1
!!$          syn_accel(itime,:,irec) =  (syn_veloc(itime+1,:,irec) - syn_veloc(itime-1,:,irec)) / (2 * DT)
!!$       enddo
!!$       syn_accel(1,:,irec) = (syn_veloc(2,:,irec) - syn_veloc(1,:,irec)) / DT
!!$       syn_accel(NSTEP,:,irec) = (syn_veloc(NSTEP,:,irec) - syn_veloc(NSTEP-1,:,irec)) /DT
!!$    endif
!!$
!!$    ! assign adjoint force
!!$
!!$    do i = 1,NCOMP
!!$
!!$       if (IKER==0) then       ! waveform
!!$
!!$          adj_syn(:,i,irec) = ( syn(:,i,irec) -  data(:,i,irec) ) * time_window(:)
!!$
!!$       else if (IKER==5) then   ! traveltime
!!$
!!$          ! minus sign is shifted from norm to adj_syn, in comparison with Tromp et al (2005)
!!$          ! thus, norm is ensured to be POSITIVE (N > 0)
!!$          norm = -DT * sum( time_window(:) * syn(:,i,irec) * syn_accel(:,i,irec) )
!!$          if (abs(norm) > EPS) adj_syn(:,i,irec) = -syn_veloc(:,i,irec) * time_window(:) / norm
!!$
!!$       else if (IKER==6) then  ! amplitude
!!$
!!$          ! norm is ensured to be POSITIVE (M > 0)
!!$          norm = DT * sum( time_window(:) * syn(:,i,irec) * syn(:,i,irec) )
!!$          if (abs(norm) > EPS) adj_syn(:,i,irec) = syn(:,i,irec) * time_window(:) / norm
!!$
!!$       endif
!!$
!!$       norm_mat(irec,i) = norm  ! could be either nothing or traveltime or amplitude
!!$
!!$       !adj_syn(:,i,irec) = time_window(:) ! test
!!$       !adj_syn(:,i,irec) = syn(:,i,irec) ! test
!!$    enddo
!!$
!!$  enddo ! loop over all receivers
!!$
!!$  ! write normalizations to file
!!$  open(15,file=trim(out_dir)//'Nbanana.dat',status='unknown')
!!$  do irec = 1,nrec
!!$     do i = 1,NCOMP
!!$        write(15,'(2i12,1e18.6)') irec, i, norm_mat(irec,i)
!!$     enddo
!!$  enddo
!!$  close(15)
!!$
!!$  end subroutine make_adjoint_source

  !----------------------------------------------------

  subroutine smooth_function(rough_fun, gamma, smooth_fun)

    ! This smooths a locally defined function via Gaussian convolution.

    double precision, intent(in) :: gamma
    double precision, dimension(NGLLX,NGLLZ,NSPEC),intent(in) :: rough_fun
    double precision, dimension(NGLLX,NGLLZ,NSPEC),intent(out) :: smooth_fun
    double precision, dimension(NGLOB) :: rough_global, smooth_global

    !------------

    ! convert local array to global vector
    call local2global(rough_fun, rough_global)

    ! smooth global function
    call smooth_global_function(rough_global, gamma, smooth_global)

    ! convert global vector to local array
    call global2local(smooth_global, smooth_fun)

  end subroutine smooth_function

  !----------------------------------------------------

  subroutine smooth_global_function(rough_global, gamma, smooth_global)

    ! copied from test_smooth.f90 and wave2d_surf.f90 on 08-Jan-2007
    ! This convolves a Gaussian function with a spatial 2D function to obtain
    ! a smoothed version of the initial function.

    double precision, intent(in) :: gamma
    double precision, dimension(NGLOB),intent(in) :: rough_global
    double precision, dimension(NGLOB),intent(out) :: smooth_global

    double precision, dimension(NGLOB) :: gaus_global_ex, gaus_global
    double precision :: dist2, dtrsh2, xcen, zcen, gaus_int
    double precision :: xtar, ztar, d, dmin
    integer :: i,iglob,igaus,icorner
    !character(len=200) :: file_smooth

  !-----------------------

  ! The parameter gamma controls the scalelength of the smoothing Gaussian;
  ! gamma is the full-width of the Gaussian (at 1/e).

  ! All points outside d^2 (dtrsh2) are set to zero.
  dtrsh2 = (1.5*gamma)**2

  if (EXAMPLE_Gaussian) then

     ! EXAMPLE global Gaussian smoothing function for one point
     ! (1) find the closest gridpoint to the target point
     xtar = 0.25*LENGTH
     ztar = 0.25*HEIGHT
     dmin = sqrt(LENGTH**2+HEIGHT**2)  ! max possible distance
     do iglob = 1,NGLOB
        d = sqrt((xtar-x(iglob))**2+(ztar-z(iglob))**2)
        if (d < dmin) then
           igaus = iglob
           dmin = d
        endif
     enddo
     xcen = x(igaus)
     zcen = z(igaus)

     ! (2) compute the example Gaussian
     gaus_global_ex(:) = 0.0
     do iglob = 1,NGLOB
        dist2 = (xcen - x(iglob))**2 + (zcen - z(iglob))**2
        if (dist2 <= dtrsh2) &
           gaus_global_ex(iglob) = (4./(PI*gamma**2)) * exp(-4.*dist2 / (gamma**2))
     enddo

   endif

   !------------------------------------
   ! Compute the SMOOTHED kernel by convolving a Gaussian with the UNSMOOTHED kernel.
   ! This involves integrating NGLOB products between a Gaussian and the unsmoothed kernel.

   print *, 'convolving the kernel with a Gaussian...'

   smooth_global(:) = 0.0

   ! There are FOUR STEPS to the convolution:
   ! (1) Obtain the center-point for the Gaussian; this could be a GLL point for
   !       high-resolution smoothing, or it could be an element corner (low-resolution).
   ! (2) Compute the Gaussian smoothing function at the local level.
   ! (3) Integrate the Gaussian over the grid; this provides the normalization for the Gaussian
   !       and accounts for Gaussians that are partially outside the grid.
   ! (4) Integrate the product of the Gaussian and the rough function.

   if (HIGH_RES_SMOOTHING) then
      ! loop over every GLL point for high-resolution smoothing
      do iglob = 1,NGLOB
         if (mod(iglob,1000) == 0) write(*,*) iglob, ' out of ', NGLOB

         xcen = x(iglob)
         zcen = z(iglob)
         gaus_global(:) = 0.0
         do i = 1,NGLOB
            dist2 = (xcen - x(i))**2 + (zcen - z(i))**2
            if (dist2 <= dtrsh2) &
                 gaus_global(i) = (4./(PI*gamma**2)) * exp(-4.*dist2 / (gamma**2))
         enddo

         !gaus_int_global(iglob) = sum( gaus_global(:) * da_global(:) )
         gaus_int = sum( gaus_global(:) * da_global(:) )
         smooth_global(iglob) = sum( rough_global(:) * gaus_global(:) * da_global(:) ) / gaus_int
      enddo   ! iglob = 1,NGLOB

   else
      ! loop over every element corner-point for low-resolution smoothing
      do icorner = 1,NSPEC_CORNER
         iglob = ielement_corner(icorner)

         xcen = x(iglob)
         zcen = z(iglob)
         gaus_global(:) = 0.0
         do i = 1,NGLOB
            dist2 = (xcen - x(i))**2 + (zcen - z(i))**2
            if (dist2 <= dtrsh2) &
                 gaus_global(i) = (4./(PI*gamma**2)) * exp(-4.*dist2 / (gamma**2))
         enddo

         !gaus_int_global(iglob) = sum( gaus_global(:) * da_global(:) )
         gaus_int = sum( gaus_global(:) * da_global(:) )
         smooth_global(iglob) = sum( rough_global(:) * gaus_global(:) * da_global(:) ) / gaus_int
      enddo   ! icorner = 1,NSPEC_CORNER

   endif

   ! write smooth-related functions to file
   ! (We can also write gaus_int_global(iglob) to file if desired.)
   if (EXAMPLE_Gaussian) then
      open(unit=19,file='fun_smooth.dat',status='unknown')
      if (HIGH_RES_SMOOTHING) then
         do iglob = 1,NGLOB
            write(19,'(8e16.6)') x(iglob), z(iglob), x_lon(iglob), z_lat(iglob), &
                 rough_global(iglob), gaus_global_ex(iglob), &
                 smooth_global(iglob), rough_global(iglob) - smooth_global(iglob)
         enddo
      else
         do icorner = 1,NSPEC_CORNER
            iglob = ielement_corner(icorner)
            write(19,'(8e16.6)') x(iglob), z(iglob), x_lon(iglob), z_lat(iglob), &
                 rough_global(iglob), gaus_global_ex(iglob), &
                 smooth_global(iglob), rough_global(iglob) - smooth_global(iglob)
         enddo
      endif
      close(19)
   endif

  ! plot rough function, Gaussian filter, smooth function, and residual
  !filename1 = 'get_smooth.csh'
  !filename2 = trim(script_dir)//'plot_smoothed_function.pl'
  !open(19,file=filename1,status='unknown')
  !write(19,'(5a,1e16.6)') trim(filename2),' ', trim(out_dir1),' ', &
  !   trim(file_smooth), gamma
  !close(19)
  !call system('chmod 755 get_smooth.csh ; get_smooth.csh')

  end subroutine smooth_global_function

  !----------------------------------------------------

  subroutine local2global(array_local, vec_global)

    ! This converts from a local array to a global array by averaging
    ! the multi-valued GLL points. This has the effect of smoothing
    ! any discontinuities that are honored by the mesh.

    double precision, dimension(NGLLX,NGLLZ,NSPEC),intent(in) :: array_local
    double precision, dimension(NGLOB),intent(out) :: vec_global
    integer :: i,j,ispec,iglob

    vec_global(:) = 0.0

    do ispec = 1,NSPEC
       do j = 1,NGLLZ
          do i = 1,NGLLX
             iglob = ibool(i,j,ispec)
             vec_global(iglob) = vec_global(iglob) + array_local(i,j,ispec)
          enddo
       enddo
    enddo

    ! now divide by the valence
    vec_global(:) = vec_global(:) / valence(:)

  end subroutine local2global

  !----------------------------------------------------

  subroutine global2local(vec_global, array_local)

    ! This converts from a global vector to a local array.
    ! By definition, the resultant array have NO jumps in values
    ! at the element boundaries.
    ! NOTE: This is useful for plotting purposes.

    double precision, dimension(NGLOB),intent(in) :: vec_global
    double precision, dimension(NGLLX,NGLLZ,NSPEC),intent(out) :: array_local
    integer :: i,j,ispec,iglob

    array_local(:,:,:) = 0.0

    do ispec = 1,NSPEC
       do j = 1,NGLLZ
          do i = 1,NGLLX
             iglob = ibool(i,j,ispec)
             array_local(i,j,ispec) = vec_global(iglob)
          enddo
       enddo
    enddo

  end subroutine global2local

  !----------------------------------------------------

  !  subroutine local2mvec(array1, array2, nmod_src, source_vec, nmod, mvec, array_fac0)
  !subroutine local2mvec(array1, array2, nmod_src, source_vec, nmod, mvec)
  subroutine local2mvec(array1, nmod_src, source_vec, nmod, mvec)

    ! This inputs three arrays and creates a vector that can be used
    ! in the conjugate gradient algorithm.

    double precision, dimension(NGLLX,NGLLZ,NSPEC),intent(in) :: array1
    integer, intent(in) :: nmod_src, nmod
    double precision, dimension(nmod_src),intent(in) :: source_vec
    double precision, dimension(nmod), intent(out) :: mvec

    !double precision, dimension(NGLLX,NGLLZ,NSPEC),intent(in), optional :: array_fac0
    double precision, dimension(NGLLX,NGLLZ,NSPEC) :: array_fac

    integer :: i,j,k,ipar,ispec

    !----------

    mvec(:) = 0.0

!!$    ! if the optional argument is not used, then default the factor of 1
!!$    if (.not. present(array_fac0)) then
!!$       array_fac(:,:,:) = 1.
!!$    else
!!$       array_fac(:,:,:) = array_fac0(:,:,:)
!!$    endif

    ! structure parameters
    k = 0
    do ipar = 1,NVAR_STRUCT
       do ispec = 1,NSPEC
          do j = 1,NGLLZ
             do i = 1,NGLLX
                k = k+1
                !if (ipar==1) mvec(k) = array1(i,j,ispec) * array_fac(i,j,ispec)  ! alpha
                !if (ipar==2) mvec(k) = array2(i,j,ispec) * array_fac(i,j,ispec)  ! beta

                if (ipar == 1) mvec(k) = array1(i,j,ispec)  ! btype (kappa, alpha, c)
                !if (ipar==2) mvec(k) = array2(i,j,ispec)  ! atype (mu, beta, beta
             enddo
          enddo
       enddo
    enddo
    if (k /= NVAR_STRUCT * NLOCAL) stop ' check that k = (NVAR_STRUCT * NLOCAL) '

    ! source parameters
    mvec(NVAR_STRUCT*NLOCAL+1 : nmod) = source_vec(:)

  end subroutine local2mvec

  !----------------------------------------------------

  !subroutine mvec2local(nmod, nmod_src, mvec, array1, array2, source_vec)
  subroutine mvec2local(nmod, nmod_src, mvec, array1, source_vec)

    ! This inputs a model vector and outputs the constituent arrays: alpha, beta, source.

    integer, intent(in) :: nmod, nmod_src
    double precision, dimension(nmod), intent(in) :: mvec
    double precision, dimension(NGLLX,NGLLZ,NSPEC), intent(out) :: array1
    double precision, dimension(nmod_src), intent(out) :: source_vec
    integer :: i,j,k,ipar,ispec

    !----------

    array1(:,:,:) = 0.0
    source_vec(:) = 0.0

    ! structure arrays
    k = 0
    do ipar = 1,NVAR_STRUCT
       do ispec = 1,NSPEC
          do j = 1,NGLLZ
             do i = 1,NGLLX
                k = k+1
                if (ipar == 1) array1(i,j,ispec) = mvec(k)  ! btype (mu, beta, beta)
                !if (ipar==2) array2(i,j,ispec) = mvec(k)  ! atype (kappa, alpha, c)
             enddo
          enddo
       enddo
    enddo
    if (k /= NVAR_STRUCT * NLOCAL) stop ' check that k = (NVAR_STRUCT * NLOCAL) '

    ! source vector
    source_vec(:) = mvec(NVAR_STRUCT*NLOCAL+1 : nmod)

  end subroutine mvec2local

  !----------------------------------------------------

  subroutine compute_norm_sq(filename, imnorm, &
        ievent_min, ievent_max, nevent, index_source, nmod, &
        mvec, mvec_prior, cov_model, norm_sq_parts, norm_parts_weight)

    ! This computes the norm-squared of a model vector using the model covariance.
    ! The dimensions of the input model vector are always the same, but this norm
    ! takes into account the number of events used, which may be less than nevent.
    ! UPDATE: If mvec_prior is present, then subtract from mvec: (m-mprior)^T Cm (m-mprior)
    ! NOTE 1: mprior is only used is imnorm = 1

    character(len=200), intent(in) :: filename
    integer, intent(in) :: imnorm, ievent_min, ievent_max, nevent, nmod
    integer, dimension(NVAR_SOURCE, nevent), intent(in) ::index_source
    double precision, dimension(nmod), intent(in) :: mvec, mvec_prior, cov_model
    double precision, dimension(NVAR), intent(in), optional :: norm_parts_weight

    !double precision, intent(out) :: norm_total, norm_struct, norm_source
    double precision, dimension(NVAR), intent(out) :: norm_sq_parts

    double precision, dimension(nmod) :: mtemp, ctemp
    double precision, dimension(NVAR) :: npw, norm_sq_parts0
    integer :: i, ievent, itemp1, itemp2, itemp3

    !----------

    norm_sq_parts0(:) = 0.0
    ctemp(:) = 0.0
    mtemp(:) = 0.0
    npw(:) = 1.0

    ! NOTE 1: if taking the norm of a gradient, use the inverse covariance matrix
    ! NOTE 2: if taking the norm of a model, use m - mprior
    ! NOTE 3: norm_parts_weight is related to the norm-squared weights (not norm)
    if (imnorm == 0) then       ! norm of gradient
        ctemp(:) = 1.0 / cov_model(:)
        mtemp(:) = mvec(:)
        if (present(norm_parts_weight)) npw(:) = 1.0 / norm_parts_weight(:)
        !if (present(norm_parts_weight)) npw(:) = 1.0 / norm_parts_weight(:)**2

    else if (imnorm == 1) then   ! norm of model
        ctemp(:) = cov_model(:)
        mtemp(:) = mvec(:) - mvec_prior(:)
        if (present(norm_parts_weight)) npw(:) = norm_parts_weight(:)
        !if (present(norm_parts_weight)) npw(:) = norm_parts_weight(:)**2

    else
        stop 'imnorm must = 0 or 1'
    endif

    ! structure part of the norm -- BETA only
    ! NOTE: division corresponds to inversion of a diagonal covariance matrix
    norm_sq_parts0(1) = sum( mtemp(1 : NLOCAL)**2 / ctemp(1 : NLOCAL) )

    ! source part of the norm -- only events that you are inverting for
    do ievent = ievent_min, ievent_max
       itemp1 = NLOCAL + index_source(1,ievent)
       itemp2 = NLOCAL + index_source(2,ievent)
       itemp3 = NLOCAL + index_source(3,ievent)

       norm_sq_parts0(2) = norm_sq_parts0(2) + mtemp(itemp1)**2 / ctemp(itemp1)
       norm_sq_parts0(3) = norm_sq_parts0(3) + mtemp(itemp2)**2 / ctemp(itemp2)
       norm_sq_parts0(4) = norm_sq_parts0(4) + mtemp(itemp3)**2 / ctemp(itemp3)
    enddo

    !norm_struct = norm_sq_parts0(1)
    !norm_source = sum( norm_sq_parts0(2:4) )
    !norm_total  = sum( norm_sq_parts0(1:4) )

    ! weight each part of the norm
    norm_sq_parts(:) = norm_sq_parts0(:) * npw(:)

    !-----------------------------------------------------
    ! write parts to file

    ! norm-squared and norm of each parameter class
    open(unit=19,file=filename,status='unknown')
    do i = 1,NVAR
       write(19,'(1i4,6e16.6)') i, &
          norm_sq_parts0(i), sqrt(norm_sq_parts0(i)), &
          norm_sq_parts(i), sqrt(norm_sq_parts(i)), &
          npw(i), sqrt(npw(i))
    enddo

    ! NOTE: THE SUMS ARE MORE RELEVANT IF BALANCING THE NORM-SQUARED
    ! structure
    write(19,'(1i4,6e16.6)') 1, &
       norm_sq_parts0(1), sqrt(norm_sq_parts0(1)), &
       norm_sq_parts(1), sqrt(norm_sq_parts(1)), 0.0, 0.0
    ! source (sum norm-squared parts only)
    write(19,'(1i4,6e16.6)') 2, &
       sum(norm_sq_parts0(2:4)), sqrt(sum(norm_sq_parts0(2:4))), &
       sum(norm_sq_parts(2:4)), sqrt(sum(norm_sq_parts(2:4))), 0.0, 0.0
    ! total (sum norm-squared parts only)
    write(19,'(1i4,6e16.6)') 1, &
       sum(norm_sq_parts0(1:4)), sqrt(sum(norm_sq_parts0(1:4))), &
       sum(norm_sq_parts(1:4)), sqrt(sum(norm_sq_parts(1:4))), 0.0, 0.0

    close(19)

  end subroutine compute_norm_sq

!!$  !-----------------------------------------------------
!!$
!!$  subroutine write_norm_sq(filename,norm_sq_parts)
!!$
!!$    double precision, dimension(NVAR), intent(in) :: norm_sq_parts
!!$    character(len=200), intent(in) :: filename
!!$    integer :: i
!!$
!!$    ! norm-squared and norm of each parameter class
!!$    open(unit=19,file=filename,status='unknown')
!!$    do i = 1,NVAR
!!$       write(19,'(1i10,2e20.10)') i, norm_sq_parts(i), sqrt(norm_sq_parts(i))
!!$    enddo
!!$
!!$    ! NOTE: THE SUMS ARE MORE RELEVANT IF BALANCING THE NORM-SQUARED
!!$    ! structure
!!$    write(19,'(1i10,2e20.10)') 1, norm_sq_parts(1), sqrt(norm_sq_parts(1))
!!$    ! source (sum norm-squared parts only)
!!$    write(19,'(1i10,2e20.10)') 2, sum(norm_sq_parts(2:4)), &
!!$       sqrt(sum(norm_sq_parts(2:4)))
!!$    ! total (sum norm-squared parts only)
!!$    write(19,'(1i10,2e20.10)') 1, sum(norm_sq_parts(1:4)), &
!!$       sqrt(sum(norm_sq_parts(1:4)))
!!$
!!$    close(19)
!!$
!!$  end subroutine write_norm_sq

end module wave2d_sub
