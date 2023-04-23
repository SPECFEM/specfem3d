!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
module data_all

  implicit none
  public

  integer, parameter, private       :: sp = selected_real_kind(6, 37)
  integer, parameter, private       :: dp = selected_real_kind(15, 307)
  integer, parameter                :: qp = selected_real_kind(33, 4931)
  logical                           :: rot_src_rec_true
  character(len=100)                :: stf_file, rot_rec_file, filename, fname

  character(len=7), allocatable, dimension(:)       :: stf_type
  character(len=040), allocatable, dimension(:)     :: recname
  character(len=100), allocatable, dimension(:,:)   :: outname, outname2, rec_full_name
  character(len=10), allocatable, dimension(:,:)    :: src_type
  integer                                           :: nrec, nt_seis

  character(len=1),dimension(3)       :: reccomp
  character(len=100), allocatable     :: simdir(:)
  character(len=4)                    :: appmynum, appidur
  integer                             :: i, it, nsim, isim, iseis, ii, iii, ishift
  real                                :: junk, rec_loc_tol, Mij(6), tshift
  real, dimension(3)                  :: rloc_xyz, rloc_xyz_tmp, rloc_rtp
  real(kind=dp), dimension(:), allocatable     :: time
  real, dimension(:,:), allocatable   :: seis, seis_fil, seis_sglcomp
  real, dimension(:,:,:), allocatable :: mij_phi
  real, allocatable, dimension(:)     :: magnitude
  real(kind=dp)                       :: dt, period, dt_seis
  real, allocatable, dimension(:)     :: thr_orig, phr_orig, th_orig
  real, dimension(3,3)                :: rot_mat
  character(len=3)                    :: rec_comp_sys
  logical                             :: load_snaps
  logical                             :: negative_time
  logical                             :: detailed_output
  character(len=16)                   :: simtype
  logical                             :: use_netcdf
  character(len=4), parameter         :: seistype = 'disp'

  ! discrete dirac sources
  real                                :: shift_fact

  ! parameters from param_post_processing
  real(kind=dp)                       :: conv_period
  real                                :: srccolat, srclon, src_depth
  character(len=7)                    :: conv_stf
  character(len=100)                  :: outdir

  real, allocatable                   :: period_final(:)
  real, allocatable                   :: trans_rot_mat(:,:)
  real, allocatable                   :: colat(:,:), lon(:,:)

  integer                             :: nelem
  integer                             :: nel_fluid
  integer                             :: nproc_mesh

end module data_all
!=========================================================================================

!=========================================================================================
module global_par

  integer, parameter         :: sp = selected_real_kind(6, 37)
  integer, parameter         :: dp = selected_real_kind(15, 307)
  integer, parameter         :: qp = selected_real_kind(33, 4931)
  integer, parameter         :: realkind = sp  ! < Choose solver precision here
  double precision, parameter :: pi = 3.1415926535898

  double precision, parameter :: zero = 0d0, half = 5d-1, third = 1d0/3d0
  double precision, parameter :: quart = 25d-2, one = 1d0, sixth = 1d0/6d0
  double precision, parameter :: two = 2d0, three=3d0, four=4d0, five=5.d0
  double precision, parameter :: fifth = 2.d-1
  double precision, parameter :: epsi = 1d-30
  real, parameter             :: epsi_real=1.e-10
  real(kind=8), parameter     :: smallval_dble = 1e-11
  real, parameter             :: decay = 3.5d0
  real, parameter             :: shift_fact1 = 1.5d0
end module global_par
!=========================================================================================

!=========================================================================================
program post_processing_seis

  use data_all
  use global_par
  use nc_postroutines, only: ncparamtype, nc_open, nc_read_recnames, nc_read_seismograms
  implicit none
  double precision     :: arg1
  type(ncparamtype)    :: ncparams
  real(4), allocatable :: nc_seis(:,:,:,:)
  real(kind=dp), allocatable :: tmp_write(:,:)

  call read_input
  if (use_netcdf) call nc_open(ncparams, nsim, simdir)

  if (rec_comp_sys == 'sph') then
     reccomp(1) = 'th'
     reccomp(2) = 'ph'
     reccomp(3) = 'r'
  else if (rec_comp_sys == 'enz') then
     reccomp(1) = 'N'
     reccomp(2) = 'E'
     reccomp(3) = 'Z'
  else if (rec_comp_sys == 'cyl') then
     reccomp(1) = 's'
     reccomp(2) = 'ph'
     reccomp(3) = 'z'
  else if (rec_comp_sys == 'xyz') then
     reccomp(1) = 'x'
     reccomp(2) = 'y'
     reccomp(3) = 'z'
  else if (rec_comp_sys == 'src') then
     reccomp(1) = 'R'
     reccomp(2) = 'T'
     reccomp(3) = 'Z'
  endif

  write(*,*) 'receiver components: ', (reccomp(i),i=1,3)
  write(*,*)

  ! input seismogram names
  allocate(recname(nrec))
  allocate(thr_orig(nrec))
  allocate(th_orig(nrec))
  allocate(phr_orig(nrec))

  ! these are the rotated receiver coordinates, having the source at the north pole
  allocate(colat(nrec,nsim))
  allocate(lon(nrec,nsim))


  ! these are the original receiver coordinates, having the source at the actual location
  if (use_netcdf) then
      call nc_read_recnames(ncparams, nrec, recname, th_orig, thr_orig, phr_orig)
      do isim = 1, nsim
          colat(:,isim) = 90 - th_orig(:)
          lon(:, isim)  = phr_orig(:)
      enddo
      ! @TODO
  else
      open(unit=61, file=trim(simdir(1))//'/Data/receiver_names.dat')
      do i=1, nrec
         read(61,*) recname(i) !, thr_orig(i), phr_orig(i)
   !      ! @TODO == on a float? maybe > 360 -> - 360
   !      if (phr_orig(i) >= 360.) phr_orig(i) = phr_orig(i) - 360.
      enddo
      close(61)
  endif
  !thr_orig = thr_orig / 180. * pi
  !phr_orig = phr_orig / 180. * pi

  ! compute exact receiver locations (GLL points they are mapped to) and write
  ! to file, should be the same for all subfolders
  do isim=1, 1 !, nsim
     if (.not. use_netcdf) then
         open(unit=20, file=trim(simdir(isim))//'/Data/receiver_pts.dat')
     endif

     open(unit=21, file=trim(outdir)//'/receiver_gll_locs_'//trim(src_type(isim,2))//'.dat')
     write(21,'(6a15)') ' ', 'colat', 'lat', 'lon', 'colat_unrot', 'lon_unrot'
     do i=1,nrec
        if (.not. use_netcdf) then
            read(20,*) colat(i,isim), lon(i,isim), junk
        endif
        if (abs(lon(i,isim)-360.) < 0.01) lon(i,isim) = 0.

        ! transform to xyz
        rloc_xyz(1) = sin(colat(i,isim) * pi / 180.) * cos(lon(i,isim) * pi / 180.)
        rloc_xyz(2) = sin(colat(i,isim) * pi / 180.) * sin(lon(i,isim) * pi / 180.)
        rloc_xyz(3) = cos(colat(i,isim) * pi / 180.)

        ! Rotate to the original (i.e. real src-rec coordinate-based) u_xyz
        rloc_xyz_tmp = rloc_xyz
        rloc_xyz(1) =   cos(srccolat) * cos(srclon) * rloc_xyz_tmp(1) &
                    & - sin(srclon) * rloc_xyz_tmp(2) &
                    & + sin(srccolat) * cos(srclon) * rloc_xyz_tmp(3)
        rloc_xyz(2) =   cos(srccolat) * sin(srclon) * rloc_xyz_tmp(1) &
                    & + cos(srclon) * rloc_xyz_tmp(2) &
                    & + sin(srccolat) * sin(srclon) * rloc_xyz_tmp(3)
        rloc_xyz(3) = - sin(srccolat) * rloc_xyz_tmp(1) &
                    & + cos(srccolat) * rloc_xyz_tmp(3)

        if (rloc_xyz(1) > 1.) rloc_xyz(1) = 1.
        if (rloc_xyz(1) < -1.) rloc_xyz(1) = -1.
        if (rloc_xyz(2) > 1.) rloc_xyz(2) = 1.
        if (rloc_xyz(2) < -1.) rloc_xyz(2) = -1.
        if (rloc_xyz(3) > 1.) rloc_xyz(3) = 1.
        if (rloc_xyz(3) < -1.) rloc_xyz(3) = -1.

        rloc_xyz = rloc_xyz / sqrt(rloc_xyz(1)**2 + rloc_xyz(2)**2 + rloc_xyz(3)**2)

        ! compute colat and lon
        rloc_rtp(2) = acos(rloc_xyz(3))

        arg1 = (rloc_xyz(1)  + smallval_dble) / &
               (sqrt(rloc_xyz(1)**2 + rloc_xyz(2)**2) + smallval_dble)

        if (arg1 > 1.) arg1 = 1.
        if (arg1 < -1.) arg1 = -1.

        if (rloc_xyz(2) >= 0.) then
           rloc_rtp(3) = acos(arg1)
        else
           rloc_rtp(3) = 2.*pi - acos(arg1)
        endif

        ! write exact receiver locations (gll points) in the earth fixed coordinate system to file
        write(21,'(a15,5f15.8)') recname(i), rloc_rtp(2) / pi * 180., &
                                 90. - rloc_rtp(2) / pi * 180., &
                                 rloc_rtp(3) / pi * 180., &
                                 colat(i,isim), &
                                 lon(i,isim)


        thr_orig(i) = rloc_rtp(2)
        phr_orig(i) = rloc_rtp(3)
     enddo
     if (.not. use_netcdf) close(20)
     close(21)
  enddo

  colat = colat / 180. * pi
  lon = lon / 180. * pi

  ! calculate moment tensor and azimuth prefactors/radiation patterns for each simulation
  allocate(mij_phi(nrec,nsim,3))
  call compute_radiation_prefactor(mij_phi, nrec, nsim, lon)

  ! convolve with source period?
  allocate(seis_fil(nt_seis,3)) ! need this for saving...
  allocate(period_final(nsim))


  ! enough to do once? period(:) is checked after reading to have the same value
  ! for all simulations
  do isim=1, nsim
     if (conv_period > 0. ) then
        if (stf_type(isim) /= 'dirac_0' .and. stf_type(isim) /= 'quheavi') then
           write(*,*) ' ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !'
           write(*,*)'    WARNING! You want to convolve although seismograms are upon a ', &
                    trim(stf_type(isim))
           write(*,*) ' ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !'
        endif

        if (conv_period < period) then
           write(*,*) 'Cannot convolve with a period shorter than allowed by the mesh/simulation!'
           write(*,*) 'convolution period, min. mesh period [s]:', conv_period, period
           call flush(6)
           stop 2
        endif

        write(*,*)' Convolve with source period [s]:', conv_period
        period_final(isim) = conv_period
     else
        period_final(isim) = period
     endif

     period = period_final(isim)
  enddo

  ! define time series
  allocate(time(nt_seis), seis(nt_seis,3))
  do iseis=1, nt_seis
     time(iseis) = (iseis - 1) * dt_seis
  enddo
  allocate(seis_sglcomp(nt_seis,3))

  ! output seismogram names
  allocate(outname(nrec,nsim), outname2(nrec,nsim), rec_full_name(nrec,nsim))
  do i=1, nrec
     do isim=1,nsim
        rec_full_name(i,isim) = trim(recname(i))//'_'//seistype//'_post'
        !if (nsim == 4) rec_full_name(i,isim) = trim(rec_full_name(i,isim))//'_mij'
        ! @TODO just a quick fix for a problem with filenames in
        ! postprocessing.csh:
        rec_full_name(i,isim) = trim(rec_full_name(i,isim))//'_mij'

        call define_io_appendix(appidur, int(conv_period))
        rec_full_name(i,isim) = trim(rec_full_name(i,isim))//'_conv'//appidur//''
        outname(i,isim) = trim(outdir)//'/SEISMOGRAMS/'//trim(rec_full_name(i,isim))

        ! if no summation, put output traces in each of the simulation directories
        outname2(i,isim) = trim(outdir)//'/SEISMOGRAMS/UNPROCESSED/'&
                            //trim(recname(i))//'_'//seistype//'_post_'&
                            //trim(src_type(isim,2))

        !write(*,*) 'outname: ', i, trim(outname(i,isim))
        !write(*,*) 'outname2: ', isim, trim(outname2(i,isim))
     enddo
  enddo

  if (conv_period > 0.) then
    tshift = shift_fact1 * conv_period
  else
     tshift = 0.
  endif
  ishift = int((shift_fact+tshift)/(time(2)-time(1)))
  write(*,*) 'ishift1:', shift_fact, tshift,time(2)-time(1)
  write(*,*) 'ishift:', (shift_fact+tshift)/(time(2)-time(1)),int((shift_fact+tshift)/(time(2)-time(1)))

  ! ----------------- Start actual post processing -------------------------
  if (use_netcdf) then
      allocate(nc_seis(nt_seis, 3, nrec, nsim))
      call nc_read_seismograms(ncparams, nt_seis, nrec, nsim, nc_seis)
  endif

  !=================
  do i=1, nrec
  !=================
     call define_io_appendix(appmynum,i)
     seis = 0.
     seis_sglcomp = 0.
     seis_fil = 0.

     !::::::::::::::::::::::::::::::::
     do isim=1,nsim
     !::::::::::::::::::::::::::::::::

        if (use_netcdf) then
            seis_sglcomp(:,:) = nc_seis(:,:,i,isim)
        else
            ! load seismograms from all directories
            open(unit=60,file = trim(simdir(isim))//'/Data/'//trim(recname(i))//'_disp.dat')
            write(*,*) 'opened ', trim(simdir(isim))//'/Data/'//trim(recname(i))//'_disp.dat'

            if (src_type(isim,1) == 'monopole') then
               do iseis=1, nt_seis
                  read(60,*) seis_sglcomp(iseis,1), seis_sglcomp(iseis,3)
               enddo
            else
               do iseis=1, nt_seis
                  read(60,*) seis_sglcomp(iseis,1), seis_sglcomp(iseis,2), seis_sglcomp(iseis,3)
               enddo
            endif
            close(60)
        endif

        if (nsim > 1 .and. detailed_output) then
           ! output seismograms for individual runs into UNPROCESSED/ folder
           ! In original rotated cylindrical frame, without azimuthal radiation factor
           open(unit=150, file=trim(outname2(i,isim))//'_s0.dat', status='new')
           open(unit=151, file=trim(outname2(i,isim))//'_ph0.dat', status='new')
           open(unit=152, file=trim(outname2(i,isim))//'_z0.dat', status='new')

           if (negative_time) then
              do it=1,nt_seis
                 write(150,*) time(it) - shift_fact - tshift, seis_sglcomp(it,1)
                 write(151,*) time(it) - shift_fact - tshift, seis_sglcomp(it,2)
                 write(152,*) time(it) - shift_fact - tshift, seis_sglcomp(it,3)
              enddo
           else
              do it=ishift+1,nt_seis
                 write(150,*) time(it-ishift-1), seis_sglcomp(it,1)
                 write(151,*) time(it-ishift-1), seis_sglcomp(it,2)
                 write(152,*) time(it-ishift-1), seis_sglcomp(it,3)
              enddo
           endif
              close(150)
              close(151)
              close(152)
        endif

        ! sum seismograms upon separate moment tensors and include azimuthal radiation factor
        ! This needs to be done prior to component rotation such that the system is still &
        call sum_individual_wavefields(seis, seis_sglcomp, nt_seis, mij_phi(i,isim,:))

     !::::::::::::::::::::::::::::::::
     enddo !isim
     !::::::::::::::::::::::::::::::::

     ! convolve with a source time function
     if (conv_period > 0.) then
        write(*,*) 'conv period:', conv_stf
        call convolve_with_stf(conv_period, dt_seis, nt_seis, conv_stf, &
                               outdir, seis, seis_fil)
     else
        seis_fil = seis
     endif

     ! rotate receiver components to coordinate system of choice
     ! doing this after the nsim loop assumes that all simulations have the same coordinates.
     ! this is fine within a Mij set, but NOT for a finite fault: In that case, we will have to
     ! amend this by another loop over fault points.
     call rotate_receiver_comp(1, rec_comp_sys, srccolat, srclon, &
                               colat(i,1), lon(i,1), thr_orig(i), phr_orig(i), &
                               nt_seis, seis_fil)

     ! write processed seismograms into joint outdir
     open(unit=50, file=trim(outname(i,1))//'_'//reccomp(1)//'.dat', status='new')
     open(unit=51, file=trim(outname(i,1))//'_'//reccomp(2)//'.dat', status='new')
     open(unit=52, file=trim(outname(i,1))//'_'//reccomp(3)//'.dat', status='new')
     if (.not. allocated(tmp_write)) allocate(tmp_write(2, nt_seis))
     if (negative_time) then
        write(*,*)' writing seismograms into joint directory: negative time'
        tmp_write(1,:) = time - shift_fact - tshift
        tmp_write(2,:) = seis_fil(:,1)
        write(50,'(2ES16.7)') tmp_write
        tmp_write(2,:) = seis_fil(:,2)
        write(51,'(2ES16.7)') tmp_write
        tmp_write(2,:) = seis_fil(:,3)
        write(52,'(2ES16.7)') tmp_write
     else
        write(*,*)' writing seismograms into joint directory: zero time'
        do it=ishift+1, nt_seis
           write(50,*) time(it-ishift-1), seis_fil(it,1)
           write(51,*) time(it-ishift-1), seis_fil(it,2)
           write(52,*) time(it-ishift-1), seis_fil(it,3)
        enddo
     endif
     close(50)
     close(51)
     close(52)

  !=================
  enddo ! nrec
  !=================

  ! plot original source and receiver locations in google earth kml file with link
  ! to seismogram images
  call save_google_earth_kml(srccolat, srclon, src_depth, Mij,period_final(1), &
                             thr_orig, phr_orig, reccomp, src_type(1,1), &
                             nsim,nrec, outdir, rec_full_name(:,1))

  write(*,*) 'writing matlab input files for record sections...'

  open(unit=99, file=trim(outdir)//'/info_matlab.dat')
  write(99,*) nrec
  if (negative_time) then
     write(99,*) nt_seis
  else
     write(99,*) nt_seis - ishift
  endif

  write(99,*) time(2) - time(1)
  write(99,*) conv_period
  do i=1, nrec
     write(99,*) colat(i,1)
  enddo
  close(99)

  open(unit=99, file=trim(outdir)//'/info_seisnames.dat')
  open(unit=98, file=trim(outdir)//'/info_seisstations.dat')
  do ii=1, 3
     do i=1, nrec
        fname = trim(rec_full_name(i,1))//'_'//reccomp(ii)//".dat"
        write(99,*) trim(fname)
        write(98,*) trim(recname(i))
     enddo
  enddo
  close(99)
  close(98)

  open(unit=99, file=trim(outdir)//'/info_seiscomp.dat')
  do i=1,3
       write(99,*) reccomp(i)
  enddo
  close(99)

  ! compute snapshots of the wavefield in a 3D sphere with given top/bottom surface and
  ! opening cross section (see param_snaps)
  if (load_snaps) call compute_3d_wavefields

  write(*,*)' DONE with seismogram post processing!'
21 format(a100)

end program post_processing_seis
!=========================================================================================

!-----------------------------------------------------------------------------------------
subroutine read_input

  use data_all
  use global_par
  implicit none

  integer, allocatable, dimension(:)  :: ishift_deltat, ishift_seisdt, ishift_straindt
  integer, allocatable, dimension(:)  :: ibeg, iend
  integer, allocatable, dimension(:)  :: nt, nt_strain, nrec_tmp, nt_seis_tmp
  real, allocatable, dimension(:)     :: dt_strain, dt_snap, nt_snap, &
                                         srccolat_tmp, srclon_tmp, src_depth_tmp
  real, allocatable, dimension(:)     :: dt_tmp, period_tmp, dt_seis_tmp
  real, allocatable, dimension(:)     :: shift_fact_tmp

  character(len=100),allocatable, dimension(:)      :: bkgrndmodel, rot_rec

  integer             :: i_param_post=500, ioerr, nval
  character(len=256)  :: line
  character(len=256)  :: keyword, keyvalue


  write(*,'(A)', advance='no') '    Reading inparam_basic...'
  open(unit=i_param_post, file='inparam_basic', status='old', action='read', &
       iostat=ioerr)
  if (ioerr /= 0) then
     write(*,*) 'Check input file ''inparam_basic''! Is it still there?'
     stop 2
  endif

  do
    read(i_param_post, fmt='(a256)', iostat=ioerr) line
    if (ioerr < 0) exit
    if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

    read(line,*) keyword, keyvalue

    select case(trim(keyword))

    case('SIMULATION_TYPE')
       if (keyvalue == 'single') then
          nsim = 1
          allocate(simdir(nsim))
          simdir(1) = "./"
       else if (keyvalue == 'moment') then
          nsim = 4
          allocate(simdir(nsim))
          simdir(1) = "MZZ"
          simdir(2) = "MXX_P_MYY"
          simdir(3) = "MXZ_MYZ"
          simdir(4) = "MXY_MXX_M_MYY"
       else if (keyvalue == 'force') then
          write(*,*) 'postprocessing for "force" simulation needs work!'
          stop 2
       endif
    end select
  enddo
  close(i_param_post)


  ! default values:
  rec_comp_sys    = 'enz'
  conv_period     = 0.
  conv_stf        = 'gauss_0'
  !seistype        = 'disp'
  load_snaps      = .false.
  outdir          = './Data_Postprocessing'
  negative_time   = .true.
  detailed_output = .false.

  write(*,'(A)', advance='no') '    Reading param_post_processing...'
  open(unit=i_param_post, file='param_post_processing', status='old', action='read', &
       iostat=ioerr)
  if (ioerr /= 0) then
     write(*,*) 'Check input file ''param_post_processing''! Is it still there?'
     stop 2
  endif

  do
    read(i_param_post, fmt='(a256)', iostat=ioerr) line
    if (ioerr < 0) exit
    if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

    read(line,*) keyword, keyvalue

    parameter_to_read : select case(trim(keyword))

    case('REC_COMP_SYS')
       rec_comp_sys = keyvalue
    case('CONV_PERIOD')
       read(keyvalue, *) conv_period
    case('CONV_STF')
       read(keyvalue, *) conv_stf
    ! Removing this option since velocity seismograms
    ! are not produced in the solver
    !case('SEISTYPE')
    !   seistype = keyvalue
    case('LOAD_SNAPS')
       read(keyvalue, *) load_snaps
    case('DATA_DIR')
       outdir = keyvalue
    case('NEGATIVE_TIME')
       read(keyvalue, *) negative_time
    case('DETAILED_OUTPUT')
       read(keyvalue, *) detailed_output

    end select parameter_to_read

  enddo
  close(i_param_post)

  tshift = 0.

  allocate(bkgrndmodel(nsim), stf_type(nsim))
  allocate(src_type(nsim,2))
  allocate(dt_tmp(nsim), period_tmp(nsim), magnitude(nsim), dt_seis_tmp(nsim))
  allocate(dt_strain(nsim), dt_snap(nsim))
  allocate(rot_rec(nsim))
  allocate(nt(nsim), nrec_tmp(nsim), nt_seis_tmp(nsim), nt_strain(nsim), nt_snap(nsim))
  allocate(ibeg(nsim), iend(nsim), srccolat_tmp(nsim), srclon_tmp(nsim), src_depth_tmp(nsim))
  allocate(shift_fact_tmp(nsim))
  allocate(ishift_deltat(nsim), ishift_seisdt(nsim), ishift_straindt(nsim))

  do isim = 1,nsim
     open(unit=99,file=trim(simdir(isim))//'/simulation.info')
     read(99,*) bkgrndmodel(isim)
     read(99,*) dt_tmp(isim)
     read(99,*) nt(isim)
     read(99,*) src_type(isim,1)
     read(99,*) src_type(isim,2)
     read(99,*) stf_type(isim)
     read(99,*) simtype
     read(99,*) period_tmp(isim)
     read(99,*) src_depth_tmp(isim)
     read(99,*) srccolat_tmp(isim)
     read(99,*) srclon_tmp(isim)
     read(99,*) magnitude(isim)
     read(99,*) nrec_tmp(isim)
     read(99,*) nt_seis_tmp(isim)
     read(99,*) dt_seis_tmp(isim)
     read(99,*) nt_strain(isim)
     read(99,*) dt_strain(isim)
     read(99,*) nt_snap(isim)
     read(99,*) dt_snap(isim)
     read(99,*) rot_rec(isim)
     read(99,*) ibeg(isim)
     read(99,*) iend(isim)
     read(99,*) shift_fact_tmp(isim)
     read(99,*) ishift_deltat(isim)
     read(99,*) ishift_seisdt(isim)
     read(99,*) ishift_straindt(isim)
     read(99,*)
     read(99,*)
     read(99,*) use_netcdf
     read(99,*) nelem
     read(99,*) nel_fluid
     read(99,*) nproc_mesh
     close(99)

     if (src_type(isim,2) == 'thetaforce' .or. src_type(isim,2) == 'phiforce') then
        write(*,'(a,/,a)') 'ERROR: postprocessing for forces with dipole radiation', &
                           '       pattern not yet implemented'
        stop 2
     endif

     write(*,*) 'Simulations: ',isim,trim(simdir(isim))
     write(*,*) '  source type:',src_type(isim,1),' ',src_type(isim,2)
     write(*,*) '  source time function:',stf_type(isim)
     write(*,*) '  simulation  type:', simtype
     write(*,*) '  magnitude:',magnitude(isim)
     write(*,*) '  receivers:',nrec_tmp(isim)
     write(*,*) '  source period:',period_tmp(isim)
     write(*,*) '  time steps:',nt(isim)
     write(*,*) '  rotate recs?',trim(rot_rec(isim))
     write(*,*) '  shift factor of stf:',shift_fact_tmp(isim)
     write(*,*) '  shift factor in samples (dt, dtseis, dtstrain):', &
                    ishift_deltat(isim), ishift_seisdt(isim), ishift_straindt(isim)
     print *
  enddo

  ! input parameter consistency
  if (minval(dt_tmp) /= maxval(dt_tmp) .or. minval(nt) /= maxval(nt) &
       .or. minval(period_tmp) /= maxval(period_tmp) &
       .or. minval(nrec_tmp) /= maxval(nrec_tmp)  .or. minval(nt_seis_tmp) /= maxval(nt_seis_tmp) &
       .or. minval(dt_seis_tmp) /= maxval(dt_seis_tmp) .or. minval(nt_strain) /= maxval(nt_strain) &
       .or. minval(dt_strain) /= maxval(dt_strain) .or. minval(nt_snap) /= maxval(nt_snap) &
       .or. minval(srccolat_tmp) /= maxval(srccolat_tmp) &
       .or. minval(srclon_tmp) /= maxval(srclon_tmp) &
       .or. minval(src_depth_tmp) /= maxval(src_depth_tmp) &
       .or. minval(dt_snap) /= maxval(dt_snap) .or.  minval(shift_fact_tmp) /= maxval(shift_fact_tmp)) then
     write(*,*)'PROBLEM with simulation.info parameters in the respective directories:'
     write(*,*)' one or more of the supposedly equal parameters differ!'
     call flush(6)
     stop 2
  endif

  do isim=1, nsim
     if (trim(bkgrndmodel(isim)) /= trim(bkgrndmodel(1))) then
        write(*,*) 'PROBLEM with simulation.info parameters in the respective directories:'
        write(*,*) '        backgroundmodels differe: ', trim(bkgrndmodel(isim)), &
                   trim(bkgrndmodel(1))
        stop 2
     endif
  enddo

  if (rec_comp_sys /= 'enz' .and. rec_comp_sys /= 'sph' .and. rec_comp_sys /= 'cyl' &
        .and. rec_comp_sys /= 'xyz' .and. rec_comp_sys /= 'src') then
     write(*,*) 'ERROR: unknown coordinate system ', rec_comp_sys
     stop 2
  endif

  srclon = srclon_tmp(1)
  srccolat = srccolat_tmp(1)
  src_depth = src_depth_tmp(1)
  dt = dt_tmp(1)
  period = period_tmp(1)
  dt_seis = dt_seis_tmp(1)
  nrec = nrec_tmp(1)
  nt_seis = nt_seis_tmp(1)
  shift_fact = shift_fact_tmp(1)

  ! unused stuff:
  deallocate(dt_tmp, period_tmp, dt_seis_tmp)
  deallocate(src_depth_tmp, srccolat_tmp, srclon_tmp)
  deallocate(nrec_tmp, nt_seis_tmp)
  deallocate(dt_strain, dt_snap, nt_snap)
  deallocate(bkgrndmodel, rot_rec)
  deallocate(ibeg, iend)
  deallocate(ishift_deltat, ishift_seisdt, ishift_straindt, shift_fact_tmp)


  ! @TODO: make lokal variables for reading and checking, global vars non array

end subroutine read_input
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_radiation_prefactor(mij_prefact, npts, nsim, longit)

  use data_all, only: src_type, magnitude, mij, rot_mat
  use data_all, only: trans_rot_mat, simtype, srccolat, srclon, outdir
  use global_par

  implicit none

  real, dimension(1:npts,1:nsim,1:3), intent(out)   :: mij_prefact
  real, dimension(1:npts), intent(in)               :: longit
  integer, intent(in)                               :: npts, nsim

  real, dimension(6)    :: Mij_scale
  character(len=100)    :: junk, fmtstring
  character(len=256)    :: keyword, keyvalue, line
  integer               :: isim, i, j, iinparam_source, ioerr
  real                  :: transrottmp(1:3,1:3), Mij_matr(3,3), amplitude

  ! This is the rotation matrix of Nissen-Meyer, Dahlen, Fournier, GJI 2007
  ! to rotate xyz coordinates

  if (.not. allocated(trans_rot_mat)) allocate(trans_rot_mat(3,3))

  rot_mat(1,1) = cos(srccolat) * cos(srclon)
  rot_mat(2,2) = cos(srclon)
  rot_mat(3,3) = cos(srccolat)
  rot_mat(2,1) = cos(srccolat) * sin(srclon)
  rot_mat(3,1) = -sin(srccolat)
  rot_mat(3,2) = 0.d0
  rot_mat(1,2) = -sin(srclon)
  rot_mat(1,3) = sin(srccolat) * cos(srclon)
  rot_mat(2,3) = sin(srccolat) * sin(srclon)

  do i=1,3
     do j=1,3
        if (abs(rot_mat(i,j)) < epsi_real) rot_mat(i,j) = 0.0
     enddo
  enddo

  trans_rot_mat(:,:) = transpose(rot_mat)

  select case(simtype)
  case('moment')
     write(*,*)'  reading CMTSOLUTION file....'
     open(unit=20000,file='CMTSOLUTION',POSITION='REWIND',status='old')
     read(20000,*) junk
     read(20000,*) junk
     read(20000,*) junk
     read(20000,*) junk
     read(20000,*) junk
     read(20000,*) junk
     read(20000,*) junk
     read(20000,*) junk, Mij(1) !Mrr
     read(20000,*) junk, Mij(2) !Mtt
     read(20000,*) junk, Mij(3) !Mpp
     read(20000,*) junk, Mij(4) !Mrt
     read(20000,*) junk, Mij(5) !Mrp
     read(20000,*) junk, Mij(6) !Mtp
     close(20000)

     Mij = Mij / 1.E7 ! CMTSOLUTION given in dyn-cm

  case('single')
     iinparam_source = 1132
     open(unit=iinparam_source, file='inparam_source', status='old', action='read', iostat=ioerr)
     if (ioerr /= 0) stop 'Check input file ''inparam_source''! Is it still there?'

     do
        read(iinparam_source, fmt='(a256)', iostat=ioerr) line
        if (ioerr < 0) exit
        if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

        read(line,*) keyword, keyvalue

        if (trim(keyword) == 'SOURCE_AMPLITUDE') then
            read(keyvalue,*) amplitude
        endif
     enddo
     Mij = 0.0
     select case(src_type(1,2))
     case('mrr')
         Mij(1) =  amplitude
     case('mtt_p_mpp')
         Mij(2) =  amplitude
         Mij(3) =  amplitude
     case('mtr', 'mrt')
         Mij(4) =  amplitude
     case('mpr', 'mrp')
         Mij(5) =  amplitude
     case('mtp', 'mpt')
         Mij(6) =  amplitude
     case('mtt_m_mpp')
         Mij(2) =  amplitude
         Mij(3) = -amplitude
     case('explosion')
         Mij(1) =  amplitude
         Mij(2) =  amplitude
         Mij(3) =  amplitude
     case default
         write(*,*) 'unknown source type: ', src_type(1,2)
     end select

  case default
     write(*,*)'unknown simulation type!', simtype
     stop
  end select

  fmtstring = '(6(E12.5))'
  write(*,*) 'Original moment tensor: (Mrr, Mtt, Mpp, Mrt, Mrp, Mtp)'
  write(*,fmtstring) (Mij(i),i=1,6)
  write(*,*) 'magnitudes of each run:'
  write(*,*) (magnitude(isim),isim=1,nsim)

  do isim=1, nsim

     Mij_scale = Mij / magnitude(isim)

     write(*,*)'Mij scaled:'
     write(*,fmtstring) Mij_scale

     select case(src_type(isim,2))
     case('mrr')
       mij_prefact(:,isim,:) = Mij_scale(1)
       mij_prefact(:,isim,2) = 0.
       write(*,*) isim, 'Simulation is mrr, prefact:', &
                  mij_prefact(1,isim,1), mij_prefact(1,isim,2), &
                  mij_prefact(1,isim,3)

     case('mtt_p_mpp')
        mij_prefact(:,isim,:) = Mij_scale(2) + Mij_scale(3)
        mij_prefact(:,isim,2) = 0.
        write(*,*) isim, 'Simulation is mpp, prefact:', &
                   mij_prefact(1,isim,1), mij_prefact(1,isim,2), &
                   mij_prefact(1,isim,3)

     case('mtr', 'mrt', 'mpr', 'mrp')
        mij_prefact(:,isim,1) =   Mij_scale(4) * cos(longit) &
                                + Mij_scale(5) * sin(longit)
        mij_prefact(:,isim,2) = - Mij_scale(4) * sin(longit) &
                                + Mij_scale(5) * cos(longit)
        mij_prefact(:,isim,3) =   Mij_scale(4) * cos(longit) &
                                + Mij_scale(5) * sin(longit)

        write(*,*) isim, 'Simulation is mtr, prefact:', &
                   mij_prefact(1,isim,1), mij_prefact(1,isim,2), &
                   mij_prefact(1,isim,3)

     case('mtp', 'mpt', 'mtt_m_mpp')
        mij_prefact(:,isim,1) = (Mij_scale(2) - Mij_scale(3)) * cos(2. * longit)  &
                                         + 2. * Mij_scale(6)  * sin(2. * longit)
        mij_prefact(:,isim,2) = (Mij_scale(3) - Mij_scale(2)) * sin(2. * longit) &
                                         + 2. * Mij_scale(6)  * cos(2. * longit)
        mij_prefact(:,isim,3) = (Mij_scale(2) - Mij_scale(3)) * cos(2. * longit)  &
                                         + 2. * Mij_scale(6)  * sin(2. * longit)

        write(*,*) isim, 'Simulation is mtp, prefact:', &
                   mij_prefact(1,isim,1), mij_prefact(1,isim,2), &
                   mij_prefact(1,isim,3)

     case('explosion')
        mij_prefact(:,isim,:) = (Mij_scale(1) + Mij_scale(2) + Mij_scale(3)) / 3.
        write(*,*) isim, 'Simulation is explosion, prefact:', &
                   mij_prefact(1,isim,1), mij_prefact(1,isim,2), &
                   mij_prefact(1,isim,3)

     case default
         write(*,*) 'unknown source type ', src_type(isim,2)
         stop

     end select

     write(*,*) 'Mij phi prefactor:', maxval(mij_prefact(:,isim,1)), &
                maxval(mij_prefact(:,isim,2)), maxval(mij_prefact(:,isim,3))
  enddo ! isim
end subroutine compute_radiation_prefactor
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine sum_individual_wavefields(field_sum, field_in, n, mij_prefact)

  implicit none

  integer, intent(in)                 :: n
  real, dimension(n,3), intent(in)    :: field_in
  real, dimension(3),   intent(in)    :: mij_prefact
  real, dimension(n,3), intent(inout) :: field_sum

  field_sum(:,1) = field_sum(:,1) + mij_prefact(1) * field_in(:,1)
  field_sum(:,2) = field_sum(:,2) + mij_prefact(2) * field_in(:,2)
  field_sum(:,3) = field_sum(:,3) + mij_prefact(3) * field_in(:,3)

end subroutine sum_individual_wavefields
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine rotate_receiver_comp(isim, rec_comp_sys, srccolat, srclon, th_rot, ph_rot, &
                                th_orig, ph_orig, nt, seis)

  use data_all, only: nsim, trans_rot_mat
  use global_par
  implicit none

  character(len=3)      :: rec_comp_sys
  real, intent(in)      :: th_rot, ph_rot   ! coordinates in the rotated (src at pole) system
  real, intent(in)      :: th_orig, ph_orig ! coordinates in the unrotated (actual src) system
  real, intent(in)      :: srccolat, srclon ! orginal source coordinates
  integer, intent(in)   :: nt,isim
  real, intent(inout)   :: seis(nt,3)
  real                  :: seis_tmp(nt,3), seisvec(3), rot(3,3)
  integer               :: i

  write(*,*) 'ROTATIONS'
  write(*,*) th_orig * 180. / pi, ph_orig * 180. / pi
  write(*,*) th_rot * 180. / pi, ph_rot * 180. / pi

  ! Need to consider *local* spherical geometry in the first place,
  ! THEN rotate the source-receiver frame to the north pole in the solver.
  ! E.g., consider the difference between source-projected and spherical
  ! coordinates for a source  away from the north pole: they are not the same,
  ! but in the framework below would  be identified as the same.

  ! Source projected frame: transform to spherical without any further rotations
  if (rec_comp_sys == 'src') then
     seis_tmp(:,1) = cos(th_rot) * seis(:,1) - sin(th_rot) * seis(:,3)
     seis_tmp(:,2) = seis(:,2)
     seis_tmp(:,3) = sin(th_rot) * seis(:,1) + cos(th_rot) * seis(:,3)

  ! Rotate from rotated u_sphiz to rotated u_xyz (both in reference, source-at-pole system)
  else
     seis_tmp(:,1) = cos(ph_rot) * seis(:,1) - sin(ph_rot) * seis(:,2)
     seis_tmp(:,2) = sin(ph_rot) * seis(:,1) + cos(ph_rot) * seis(:,2)
     seis_tmp(:,3) = seis(:,3)

     ! Rotate to the original (i.e. real src-rec coordinate-based) u_xyz
     if (srccolat > epsi_real .or. srclon > epsi_real) then
        rot = transpose(trans_rot_mat(:,:))
        do i=1,nt
           seisvec = seis_tmp(i,:)
           seis_tmp(i,:) = matmul(rot,seisvec)
        enddo
     endif

  endif

  ! Rotate to coordinate system of choice
  if (rec_comp_sys == 'enz') then
     seis(:,1) = - cos(th_orig) * cos(ph_orig) * seis_tmp(:,1) &
               & - cos(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
               & + sin(th_orig) * seis_tmp(:,3)
     seis(:,2) = - sin(ph_orig) * seis_tmp(:,1) &
               & + cos(ph_orig) * seis_tmp(:,2)
     seis(:,3) =   sin(th_orig) * cos(ph_orig) * seis_tmp(:,1) &
               & + sin(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
               & + cos(th_orig) * seis_tmp(:,3)


  else if (rec_comp_sys == 'sph') then
     seis(:,1) =   cos(th_orig) * cos(ph_orig) * seis_tmp(:,1) &
               & + cos(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
               & - sin(th_orig) * seis_tmp(:,3)
     seis(:,2) = - sin(ph_orig) * seis_tmp(:,1) &
               & + cos(ph_orig) * seis_tmp(:,2)
     seis(:,3) =   sin(th_orig) * cos(ph_orig) * seis_tmp(:,1) &
               & + sin(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
               & + cos(th_orig) * seis_tmp(:,3)

  else if (rec_comp_sys == 'cyl') then
     seis(:,1) =   cos(ph_orig) * seis_tmp(:,1) + sin(ph_orig) * seis_tmp(:,2)
     seis(:,2) = - sin(ph_orig) * seis_tmp(:,1) + cos(ph_orig) * seis_tmp(:,2)
     seis(:,3) =   seis_tmp(:,3)

  else if (rec_comp_sys == 'xyz') then
     seis = seis_tmp

  else if (rec_comp_sys == 'src') then
     seis = seis_tmp ! taken from above

  else
     write(*,*)'unknown component system',rec_comp_sys
     stop 2
  endif

end subroutine rotate_receiver_comp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!! convolve seismograms computed for dirac delta with a Gaussian
subroutine convolve_with_stf(t_0, dt, nt, stf, outdir, seis, seis_fil)

  use data_all, only: stf_type
  use global_par, only: pi, decay, shift_fact1, dp
  implicit none

  integer, intent(in)            :: nt
  real(kind=dp), intent(in)      :: t_0, dt
  character(len=100), intent(in) :: outdir
  real                           :: time(nt)
  real                           :: tau_j, source, sqrt_pi_inv
  integer                        :: i, j, N_j, irec, lffile, ishift
  real, intent(in)               :: seis(nt,3)
  real, intent(out)              :: seis_fil(nt,3)
  real                           :: src_array(nt), temp_expo, alpha
  character(len=7), intent(in)   :: stf
  character(len=4)               :: appidur, appirec

  write(*,*)
  write(*,*) 'Convolving with period=', t_0
  write(*,*) 'convolve:', stf, stf_type(1), maxval(seis)

  N_j = int(2. * shift_fact1 * t_0 / dt)
  if (N_j > nt) N_j = nt

  call define_io_appendix(appidur, int(t_0))
  alpha = decay / t_0
  sqrt_pi_inv = 1. / dsqrt(pi)

  do i=1, nt
    time(i) = dt * real(i)
    seis_fil(i,:) = 0.
    do j=1, N_j
       tau_j = dble(j) * dt
       ! define source time function
       if (stf == 'gauss_0') then
          temp_expo = alpha*(tau_j-shift_fact1*t_0)
          if (temp_expo < 50.) then
             source = alpha * exp(-temp_expo**2 ) * sqrt_pi_inv / pi
          else
             source = 0.
          endif
       !else if (stf == 'quheavi') then
       !   source = 0.5 * (1.0 + erf((tau_j - shift_fact1 * t_0) / t_0))
       else if (stf == 'gauss_1') then
          source = -2. * (decay / t_0)**2 * (tau_j - shift_fact1 * t_0) * &
                           exp(-( (decay / t_0 * (tau_j - shift_fact1 * t_0))**2) )
          source = source / ( decay / t_0 * sqrt(2.) * exp(-2.) )
       else
          write(*,*) ' other source time function not implemented yet!', stf
          stop 2
       endif
       ! actual time domain convolution
       if (i > j .and. i-j <= nt) seis_fil(i,:) = seis_fil(i,:) + seis(i-j,:) * source * dt
       ! buffer stf for file output
       if (i == 1 ) src_array(j) = source
    enddo
  enddo

  ! @TODO why factor pi? only for Gaussian maybe?
  seis_fil = seis_fil * pi
  write(*,*) 'convolve:', stf, stf_type(1), maxval(seis_fil)

  ! Output source time function as used here
  open(unit=55, file=trim(outdir)//'/stf_'//trim(stf)//'_'//appidur//'sec.dat')
  do i=1, N_j
    write(55,*) time(i), src_array(i)
  enddo
  close(55)

end subroutine convolve_with_stf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine save_google_earth_kml(srccolat1, srclon1, srcdepth, Mij, per, rcvcolat, &
                                 rcvlon, reccomp, src_type, nsim, &
                                 num_rec_glob, outdir, receiver_name)

  use global_par, only: pi
  implicit none

  integer, intent(in)   :: num_rec_glob,nsim
  real, intent(in)      :: srccolat1, srclon1, srcdepth, rcvcolat(1:num_rec_glob), &
                           rcvlon(1:num_rec_glob)
  real, intent(in)      :: Mij(6),per

  character(len=100), intent(in)            :: receiver_name(1:num_rec_glob)
  character(len=100), intent(in)            :: outdir
  character(len=10), intent(in)             :: src_type
  character(len=1),dimension(3), intent(in) :: reccomp

  real                  :: slon, slat, rlon(1:num_rec_glob), rlat(1:num_rec_glob)
  integer               :: i
  character(len=4)      :: app
  character(len=2)      :: comp(3)
  character(len=100)    :: fname2

  write(*,*) 'writing google earth kml file for plotting earthquake and receiver locations/seismograms...'
  write(*,*) 'Check it out: '//trim(outdir)//'/googleearth_src_rec_seis.kml'

  slat=90.-srccolat1*180./pi
  slon=srclon1*180./pi
  if (slon > 180.) slon=slon-360.

  rlat=90.-rcvcolat*180./pi
  rlon=rcvlon*180./pi
  do i=1,num_rec_glob
     if (rlon(i) > 180.) rlon(i)=rlon(i)-360.
  enddo
  open(unit=88,file=trim(outdir)//'/googleearth_src_rec_seis.kml')

  write(88,14)'<?xml version="1.0" encoding="UTF-8"?>'
  write(88,15)'<kml xmlns="http://earth.google.com/kml/2.0">'
  write(88,16)'<Document>'

  write(88,*)
  write(88,*)'<name > earthquake-receiver configuration < /name>'
  write(88,*)'<LookAt>'
  write(88,12)'<longitude>',slon,'</longitude > < latitude>',slat,'</latitude>'
  write(88,*)'<range > 7000000 < /range > < tilt > 0 < /tilt > < heading > 0 < /heading>'
  write(88,*)'</LookAt>'
  write(88,*)
  write(88,*)'......'
  write(88,*)'<Placemark>'
  write(88,*)'<Style id="earthquake">'
  write(88,*)'<IconStyle>'
  write(88,*)'<scale > 5 < /scale>'
  write(88,*)'<Icon>'
  write(88,*)'<href > http://maps.google.com/mapfiles/kml/shapes/earthquake.png < /href>'
  write(88,*)'</Icon>'
  write(88,*)'</IconStyle>'
  write(88,*)'<LabelStyle>'
  write(88,*)'<scale > 5 < /scale>'
  write(88,*)'</LabelStyle>'
  write(88,*)'</Style>'
  write(88,*)'<name > earthquake < /name>'
  write(88,*) '<description > Event details:'
  write(88,20) ' colat,lon [deg]:',srccolat1*180./pi,srclon1*180./pi
  write(88,21)' source depth [km]',srcdepth
  write(88,23)'Mrr=',Mij(1)
  write(88,23)'Mtt=',Mij(2)
  write(88,23)'Mpp=',Mij(3)
  write(88,23)'Mtr=',Mij(4)
  write(88,23)'Mpr=',Mij(5)
  write(88,23)'Mtp=',Mij(6)
  write(88,21)'source period [s]:',per
  write(88,*)'</description>'
  write(88,13)'<Point > < coordinates>',slon,',',slat,'</coordinates > < /Point>'
  write(88,*)'</Placemark>'

  do i=1,num_rec_glob
     write(88,*)
     write(88,*) '<Placemark>'
     write(88,*) '<Style id="cam">'
     write(88,*) '<IconStyle>'
   write(88,*)'<scale > 2 < /scale>'
     write(88,*) '<Icon>'
        write(88,*)'<href > http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png < /href>'
     write(88,*) '</Icon>'
     write(88,*) '</IconStyle>'
  write(88,*)'<LabelStyle>'
  write(88,*)'<scale > 2 < /scale>'
  write(88,*)'</LabelStyle>'
     write(88,*) '</Style>'
     write(88,17) '<name>',trim(receiver_name(i)),'  # ',i,'</name>'
     call define_io_appendix(app,i)
     write(88,119) '<description > station ',trim(receiver_name(i))
     write(88,20) ' colat,lon [deg]:',rcvcolat(i)*180./pi,rcvlon(i)*180./pi

     fname2='GRAPHICS/'//trim(receiver_name(i))//'_'//reccomp(1)//'.gif'
     write(88,*) '<img src="',trim(fname2),'" > < /img>'

     fname2='GRAPHICS/'//trim(receiver_name(i))//'_'//reccomp(3)//'.gif'
     write(88,*) '<img src="',trim(fname2),'" > < /img>'

     if (nsim > 1 .or. src_type /= 'monopole') then
        fname2='GRAPHICS/'//trim(receiver_name(i))//'_'//reccomp(2)//'.gif'
        write(88,*) '<img src="',trim(fname2),'" > < /img>'
     endif
     write(88,*) '</description>'
     write(88,13) '<Point > < coordinates>',rlon(i),',',rlat(i),'</coordinates > < /Point>'
     write(88,*) '</Placemark>'
  enddo

  write(88,*)'......'
  write(88,*)
  write(88,*)'</Document>'
  write(88,*)'</kml>'

  close(88)

12 format(a16,f14.2,a23,f14.2,a12)
13 format(a23,f14.2,a1,f14.2,a23)
14 format(a39)
15 format(a46)
16 format(a11)
17 format(a7,a9,a10,i4,a7)
18 format(a36,a4,a14)
19 format(a24,a15)
119 format(a24,a30)
20 format(A18,f14.2,f14.2)
21 format(A18,f14.2)
23 format(a5,1pe14.2)

end subroutine save_google_earth_kml
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_io_appendix(app,iproc)
! Defines the 4 digit character string appended to any
! data or io file related to process myid.

  integer, intent(in)           :: iproc
  character(len=4), intent(out) :: app

  write(app,"(I4.4)") iproc

end subroutine define_io_appendix
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_3d_wavefields

  use data_all
  use global_par
  implicit none

  integer                               :: iproc, npts, npts_top, npts_bot, npts_meri, &
                                           nphi, snapskip, snap1, snap2
  real, dimension(:,:), allocatable     :: coord, disp, disptot, disptot_sum, azi_xsect2, &
                                           azi_xsect1
  real, dimension(:,:,:), allocatable   :: azi, fluid_prefact, azi_meri, mij_snap
  real, dimension(:), allocatable       :: vp, x, y, z, vptop, vpbot, vpmeri, xtot, ytot, &
                                           ztot, azi_phi_meri
  double precision, dimension(:), allocatable :: xtop, ytop, ztop, xbot, ybot, zbot, &
                                                 azi_phi_top, azi_phi_bot
  real, dimension(:), allocatable       :: xmeri, ymeri, zmeri, longit_snap
  real                                  :: dphi, phi0, prem, r, theta_meri, smallval_meri, &
                                           dr, theta0, theta_bot, theta_top, phi
  integer, dimension(:), allocatable    :: ind_proc_top_tmp, ind_pts_top_tmp, &
                                           ind_proc_bot_tmp, ind_pts_bot_tmp
  integer, dimension(:), allocatable    :: ind_proc_top, ind_pts_top, ind_proc_bot
  integer, dimension(:), allocatable    :: ind_pts_bot, azi_ind_top, azi_ind_bot, &
                                           azi_ind_meri
  integer, dimension(:), allocatable    :: ind_proc_meri_tmp, ind_pts_meri_tmp, &
                                           ind_proc_meri, ind_pts_meri
  real                                  :: smallval_north_top,  smallval_south_top, &
                                           smallval_north_bot,  smallval_south_bot, r0, &
                                           coord9(9, 2), disp9(9, 3)
  character(len=200)                    :: filename1
  logical                               :: use_meri, use_top, use_bot

  character(len=4)                      :: appmynum2
  integer                               :: nptstot, k1, k2, k3, npts_fluid, j, &
                                           naang, npts_read, npts_fluid_read
  real                                  :: rbot, rtop, phi_incr
  double precision                      :: r8, theta8, phi8

  integer             :: i_param_post=500, ioerr, nval
  character(len=256)  :: line
  character(len=256)  :: keyword, keyvalue

  ! read snap plot parameters
  ! default values
  phi0       = 0.
  dphi       = 85.
  rtop       = 6371
  rbot       = 3190
  theta_meri = 60.
  use_top    = .true.
  use_bot    = .true.
  use_meri   = .false.
  snap1      = 1
  snap2      = 1
  snapskip   = 1

  write(*,'(A)', advance='no') '    Reading param_post_processing...'
  open(unit=i_param_post, file='param_post_processing', status='old', action='read', &
       iostat=ioerr)
  if (ioerr /= 0) then
     write(*,*) 'Check input file ''param_post_processing''! Is it still there?'
     stop 2
  endif

  do
    read(i_param_post, fmt='(a256)', iostat=ioerr) line
    if (ioerr < 0) exit
    if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

    read(line,*) keyword, keyvalue

    parameter_to_read : select case(trim(keyword))

    case('3D_PHI_START')
       read(keyvalue, *) phi0
    case('3D_PHI_END')
       read(keyvalue, *) dphi
    case('3D_RTOP')
       read(keyvalue, *) rtop
    case('3D_RBOT')
       read(keyvalue, *) rbot
    case('3D_MERI_COLAT')
       read(keyvalue, *) theta_meri
    case('3D_PLOT_TOP')
       read(keyvalue, *) use_top
    case('3D_PLOT_BOT')
       read(keyvalue, *) use_bot
    case('3D_PLOT_MERI')
       read(keyvalue, *) use_meri
    case('3D_SNAP_BEG')
       read(keyvalue, *) snap1
    case('3D_SNAP_END')
       read(keyvalue, *) snap2
    case('3D_SNAP_STRIDE')
       read(keyvalue, *) snapskip

    end select parameter_to_read

  enddo
  close(i_param_post)

  write(*,*) 'starting azimuth/phi for cross section on the right [deg]:',phi0
  write(*,*) 'ending azimuth/phi for cross section on the left [deg]:',dphi
  write(*,*) 'top surface [km]:',rtop
  write(*,*) 'bottom surface [km]:',rbot
  write(*,*) 'colatitude of meridional cross section:',theta_meri
  write(*,*) '1st,last snap, skipfactor:',snap1,snap2,snapskip
  write(*,*) 'consider meridional cross section?',use_meri
  write(*,*) 'consider top surface?',use_top
  write(*,*) 'consider bottom surface?',use_bot

  phi0 = phi0 / 180. * pi
  dphi = dphi / 180. * pi
  rtop = rtop * 1000.
  rbot = rbot * 1000.
  theta_meri = theta_meri * pi / 180.

  npts = nelem * 16
  npts_read = nelem * 9

  nptstot = npts * nproc_mesh
  write(*,*) 'number of points per proc, total points:', npts, nptstot

  ! load and construct global mesh (one semi-disk)
  write(*,*) 'reading partitioned mesh...'
  allocate(coord(nptstot,2))
  smallval_north_top = rtop
  smallval_south_top = rtop
  smallval_north_bot = rbot
  smallval_south_bot = rbot

  write(*,*) 'Smallest distance to rtop (North,South) BEFORE [km]:', &
             real(smallval_north_top/1000.), real(smallval_south_top/1000.)
  write(*,*) 'Smallest distance to rbot (North,South) BEFORE [km]:', &
             real(smallval_north_bot/1000.), real(smallval_south_bot/1000.)

  do iproc=0, nproc_mesh-1
     call define_io_appendix(appmynum,iproc)
     open(unit=99,file=trim(simdir(1))//'/Data/glob_grid_'//appmynum//'.dat')
     i = 0
     do ii=1,npts_read/9
        do iii=1,9
            read(99,*) coord9(iii,1), coord9(iii,2)
         enddo
         coord(iproc*npts+i+1:iproc*npts+i+2,:) = coord9(1:2,:)
                                                  ! (0,0), (0,npol/2)
         coord(iproc*npts+i+3,:) = coord9(5,:)    ! (npol/2,npol/2)
         coord(iproc*npts+i+4,:) = coord9(4,:)    ! (npol/2,0)
         coord(iproc*npts+i+5:iproc*npts+i+6,:) = coord9(2:3,:)
                                                  ! (0,npol/2),(0,npol)
         coord(iproc*npts+i+7,:) = coord9(6,:)    ! (npol/2,npol)
         coord(iproc*npts+i+8,:) = coord9(5,:)    ! (npol/2,npol/2)
         coord(iproc*npts+i+9:iproc*npts+i+10,:) = coord9(4:5,:)
                                                  ! (npol/2,0),(npol/2,npol/2)
         coord(iproc*npts+i+11,:) = coord9(8,:)   ! (npol,npol/2)
         coord(iproc*npts+i+12,:) = coord9(7,:)   ! (npol,0)
         coord(iproc*npts+i+13:iproc*npts+i+14,:) = coord9(5:6,:)
                                                  ! (npol/2,npol/2),(npol/2,npol)
         coord(iproc*npts+i+15,:) = coord9(9,:)   ! (npol,npol)
         coord(iproc*npts+i+16,:) = coord9(8,:)   ! (npol,npol/2)

         ! determine minimal distance from rtop and rbot
         do iii=1,16
            r0 = sqrt(coord(iproc*npts+i+iii,1)**2 + coord(iproc*npts+i+iii,2)**2)
            if (coord(iproc*npts+i+iii,2) > 0.0 &
                    .and. abs(r0-rtop) < smallval_north_top) then ! north
               smallval_north_top = abs(r0-rtop)
            else if (coord(iproc*npts+i+iii,2) < 0.0 &
                    .and. abs(r0-rtop) < smallval_south_top) then ! south
               smallval_south_top = abs(r0-rtop)
            endif

            if (coord(iproc*npts+i+iii,2) > 0.0 &
                    .and. abs(r0-rbot) < smallval_north_bot) then ! north
               smallval_north_bot = abs(r0-rbot)
            else if (coord(iproc*npts+i+iii,2) < 0.0 &
                    .and. abs(r0 -rbot) < smallval_south_bot) then ! south
               smallval_south_bot = abs(r0-rbot)
            endif
        enddo
        i=i+16
     enddo
     close(99)
  enddo

  ! add some tolerance (.1% of dist to closest point should be smaller then one element)
  smallval_north_top = smallval_north_top * 1.001
  smallval_south_top = smallval_south_top * 1.001
  smallval_north_bot = smallval_north_bot * 1.001
  smallval_south_bot = smallval_south_bot * 1.001

  write(*,*) 'Smallest distance to rtop (North,South) AFTER [km]:', &
       real(smallval_north_top/1000.),real(smallval_south_top/1000.)
  write(*,*) 'Smallest distance to rbot (North,South) AFTER [km]:', &
       real(smallval_north_bot/1000.),real(smallval_south_bot/1000.)

  if (use_top .or. use_meri) then
     allocate(ind_proc_top_tmp(floor(real(nptstot)/10.)))
     allocate(ind_pts_top_tmp(floor(real(nptstot)/10.)))
  endif

  if (use_bot .or. use_meri) then
     allocate(ind_proc_bot_tmp(floor(real(nptstot)/10.)))
     allocate(ind_pts_bot_tmp(floor(real(nptstot)/10.)))
  endif

  k1 = 0
  k2 = 0

  do iproc=0, nproc_mesh-1
     do i=1, npts
        ! check for top and bottom radii
        r0 = sqrt(coord(iproc*npts+i,1)**2 + coord(iproc*npts+i,2)**2)

        if (use_top .or. use_meri) then
           if ( (coord(iproc*npts+i,2) >= 0.0 .and. abs(r0-rtop) <= smallval_north_top) .or. &
                (coord(iproc*npts+i,2) < 0.0 .and. abs(r0-rtop) <= smallval_south_top)) then
              k1 = k1 + 1
              ind_proc_top_tmp(k1) = iproc
              ind_pts_top_tmp(k1) = i
           endif
        endif

        if (use_bot .or. use_meri) then
           if ( (coord(iproc*npts+i,2) >= 0.0 .and. abs(r0-rbot) <= smallval_north_bot) .or. &
                (coord(iproc*npts+i,2) < 0.0 .and. abs(r0-rbot) <= smallval_south_bot)) then
              k2 = k2 + 1
              ind_proc_bot_tmp(k2) = iproc
              ind_pts_bot_tmp(k2) = i
           endif
        endif
     enddo
  enddo

  npts_top = k1
  npts_bot = k2

  write(*,*) '# points on top,bottom:', npts_top, npts_bot

  write(*,*)'allocating index arrays for surfaces....'
  if (use_top .or. use_meri) then
     allocate(ind_proc_top(npts_top))
     allocate(ind_pts_top(npts_top))
     ind_proc_top = ind_proc_top_tmp(1:npts_top)
     ind_pts_top = ind_pts_top_tmp(1:npts_top)
     deallocate(ind_proc_top_tmp, ind_pts_top_tmp)
  endif

  if (use_bot .or. use_meri) then
     allocate(ind_proc_bot(npts_bot))
     allocate(ind_pts_bot(npts_bot))
     ind_proc_bot = ind_proc_bot_tmp(1:npts_bot)
     ind_pts_bot = ind_pts_bot_tmp(1:npts_bot)
     deallocate(ind_proc_bot_tmp,ind_pts_bot_tmp)
  endif

  !----------------------------------------------------------------------
  ! MERIDIONAL based on rtop and rbottom
  if (use_meri) then

     write(*,*)'computing meridional preparameters....'

     ! find closest theta at rbot
     smallval_meri = 2.d0*pi
     do i=1,npts_bot
        r0 = sqrt(coord(ind_proc_bot(i)*npts+ind_pts_bot(i),1)**2 &
                    + coord(ind_proc_bot(i)*npts+ind_pts_bot(i),2)**2)
        theta0 = atan(coord(ind_proc_bot(i)*npts+ind_pts_bot(i),1) &
                        /coord(ind_proc_bot(i)*npts+ind_pts_bot(i),2)+epsi)
        if ( theta0 < 0. ) theta0 = pi + theta0
        if (theta0 == 0.0 .and. coord(iproc*npts+i,2) < 0.0) theta0=pi
        if (abs(theta_meri-theta0) < smallval_meri) then
           smallval_meri = abs(theta_meri-theta0)
           theta_bot = theta0
        endif
     enddo
     write(*,*)'theta meri,theta closest at rbot:',theta_meri*180./pi,theta_bot*180./pi
     theta_meri = theta_bot

     ! find theta at rtop closest to theta from rbot
     smallval_meri = 2.d0*pi
     do i=1,npts_top
        r0 = sqrt(coord(ind_proc_top(i)*npts+ind_pts_top(i),1)**2 &
                    + coord(ind_proc_top(i)*npts+ind_pts_top(i),2)**2)
        theta0 = atan(coord(ind_proc_top(i)*npts+ind_pts_top(i),1) &
                        / coord(ind_proc_top(i)*npts+ind_pts_top(i),2)+epsi)
        if ( theta0 < 0. ) theta0 = pi + theta0
        if (theta0 == 0.0 .and. coord(iproc*npts+i,2) < 0.0) theta0=pi
        if (abs(theta_bot-theta0) < smallval_meri) then
           smallval_meri = abs(theta_bot-theta0)
           theta_top = theta0
        endif
     enddo

     smallval_meri=abs(theta_top-theta_bot)
     write(*,*)'theta closest at rtop and smallval:',theta_top*180./pi,smallval_meri*180./pi
     if (theta_top > theta_bot) then
        theta_meri = theta_bot+ smallval_meri/2.
     else if (theta_top < theta_bot) then
        theta_meri = theta_bot-smallval_meri/2.
     else
        theta_meri = theta_bot
     endif
     smallval_meri=smallval_meri*2.5

     k3=0
     allocate(ind_proc_meri_tmp(floor(real(nptstot)/10.)))
     allocate(ind_pts_meri_tmp(floor(real(nptstot)/10.)))

     do iproc=0,nproc_mesh-1
        do i=1,npts
           r0 = sqrt(coord(iproc*npts+i,1)**2+coord(iproc*npts+i,2)**2)
           theta0=atan(coord(iproc*npts+i,1)/(coord(iproc*npts+i,2)+epsi))
           if ( theta0 < 0. ) theta0 = pi + theta0
           if (theta0 == 0.0 .and. coord(iproc*npts+i,2) < 0.0) theta0=pi
           if ( r0 >= rbot .and. abs(theta0-theta_meri) <= smallval_meri ) then
              k3=k3+1
              ind_proc_meri_tmp(k3) = iproc
              ind_pts_meri_tmp(k3) = i
           endif
        enddo
     enddo
     npts_meri=k3

     write(*,*)'# points on meridional:',npts_meri
     allocate(ind_proc_meri(npts_meri),ind_pts_meri(npts_meri))
     ind_proc_meri = ind_proc_meri_tmp(1:npts_meri)
     ind_pts_meri = ind_pts_meri_tmp(1:npts_meri)
     deallocate(ind_proc_meri_tmp,ind_pts_meri_tmp)
  endif ! use_meri

  ! xyz coordinates-----------------------------------------------------------------------
  write(*,*)'defining xyz...'
  allocate(x(1:2*nptstot),y(1:2*nptstot),z(1:2*nptstot));x=0.;y=0.;z=0.

  ! left cross section--------------------------------------------------------------------
  call sphi2xy(x(1:nptstot),y(1:nptstot),coord(:,1),phi0,nptstot)
  z(1:nptstot)=coord(1:nptstot,2)
  write(*,*)'max s,x,y:',maxval(coord(:,1)),maxval(x),maxval(y)

  ! right cross section-------------------------------------------------------------------
  call sphi2xy(x(nptstot+1:2*nptstot),y(nptstot+1:2*nptstot),coord(:,1),phi0+dphi,nptstot)
  z(nptstot+1:2*nptstot)=coord(1:nptstot,2)

  write(*,*)'max s,x,y:',maxval(coord(:,1)),maxval(x),maxval(y)

  ! save xyz
  allocate(vp(2*nptstot))
  do i=1,2*nptstot
     r= sqrt(x(i)**2+y(i)**2+z(i)**2)
     vp(i)=prem(r,'v_p')
  enddo

  ! rescale vp
  vp=vp/maxval(vp)*0.002-0.001

  filename1=trim(outdir)//'/SNAPS/mesh_xsect'
  call write_VTK_bin_scal(x,y,z,vp,2*nptstot,0,filename1)
  filename1=trim(outdir)//'/SNAPS/mesh_xsect_cell'
  call write_VTK_bin_scal_topology(x,y,z,vp,2*nptstot/4,filename1)


  ! top surface--------------------------------------------------------------------------
  if (use_top) then
     k1 = 0
     write(*,*) 'defining top surface...'

     naang = floor(real(npts_top)/real(2.))**2*6*4
     write(*,*)'points on top surface:',naang
     allocate(xtop(1:naang),ytop(1:naang),ztop(1:naang))
     allocate(azi_phi_top(1:naang),azi_ind_top(1:naang))
     call construct_surface_cubed_sphere(npts_top,npts,dble(rtop),ind_proc_top, &
                                         ind_pts_top,nptstot,dble(coord),k1, &
                                         dble(dphi),dble(phi0),'outside',naang,xtop, &
                                         ytop,ztop,azi_phi_top,azi_ind_top)

     ! extract vp
     write(*,*)'allocating vptop...',k1
     allocate(vptop(k1))
     vptop(1:k1) = vp(ind_proc_top(1)*npts+ind_pts_top(1))

     write(*,*)'save into cell vtk...',size(xtop),k1
     call flush(6)
     filename1=trim(outdir)//'/SNAPS/mesh_top_cell'
     call write_VTK_bin_scal_topology(real(xtop(1:k1)),real(ytop(1:k1)),real(ztop(1:k1)), &
                                      vptop(1:k1),k1/4,filename1)

     call flush(6)
  endif

  ! bottom surface -----------------------------------------------------------------------
  if (use_bot) then
     k2=0
     write(*,*)'defining bot surface...'

     naang = floor(real(npts_bot)/real(2.))**2*6*4
     write(*,*)'points on bottom surface:',naang
     allocate(xbot(1:naang),ybot(1:naang),zbot(1:naang))
     allocate(azi_phi_bot(1:naang),azi_ind_bot(1:naang))

     call construct_surface_cubed_sphere(npts_bot,npts,dble(rbot),ind_proc_bot, &
                                         ind_pts_bot,nptstot,dble(coord),k2, &
                                         dble(dphi),dble(phi0),'innside',naang,xbot, &
                                         ybot,zbot,azi_phi_bot,azi_ind_bot)

     ! extract vp
     allocate(vpbot(k2))
     vpbot(1:k2) = vp(ind_proc_bot(1)*npts+ind_pts_bot(1))

     ! save into cell vtk
     filename1=trim(outdir)//'/SNAPS/mesh_bot_cell'
     call write_VTK_bin_scal_topology(real(xbot(1:k2)),real(ybot(1:k2)),real(zbot(1:k2)), &
                                      vpbot,k2/4,filename1)

     call flush(6)
  endif

  ! meridional cross section -------------------------------------------------------------
  k3=0
  if (use_meri) then
     write(*,*)'defining meridional cross section...'
     dr=(6371000.-rbot)/npts_meri
     write(*,*)'# pts on rmeri and average spacing [km]:',npts_meri,dr/1000.
     allocate(xmeri(1:7*npts_meri**2))
     allocate(ymeri(1:7*npts_meri**2))
     allocate(zmeri(1:7*npts_meri**2))
     xmeri=0.
     ymeri=0.
     zmeri=0.
     allocate(azi_ind_meri(7*npts_meri**2),azi_phi_meri(7*npts_meri**2))
     allocate(vpmeri(1:7*npts_meri**2))

     do i=1,npts_meri
        r0 = sqrt( coord(ind_proc_meri(i)*npts+ind_pts_meri(i),1)**2 &
                    + coord(ind_proc_meri(i)*npts+ind_pts_meri(i),2)**2)
        nphi=max(floor(r0*pi/dr/1.),1)
        phi_incr = pi/nphi
        write(*,*)i,'r0,nphi,phi_incr [km]:',r0/1000.,nphi,phi_incr*r0/1000.
        do j=1,nphi
           k3=k3+1
           phi=(j-1)*phi_incr
           call rthetaphi2xyz(xmeri(k3),ymeri(k3),zmeri(k3),r0,theta_meri,phi)
           azi_ind_meri(k3) = ind_proc_meri(i)*npts+ind_pts_meri(i)
           azi_phi_meri(k3) = phi
           vpmeri(k3) = vp(ind_proc_meri(i)*npts+ind_pts_meri(i))
        enddo
     enddo

     ! save into AVS
     write(*,*)'writing vtk bin file...'
     filename1=trim(outdir)//'/SNAPS/mesh_meri'
     call write_VTK_bin_scal(xmeri,ymeri,zmeri,vpmeri,k3,0,filename1)

  endif !use_meri

  deallocate(coord)
  deallocate(vp)

  ! assembling everything to one coordinate array-----------------------------------------

  allocate(xtot(2*nptstot+k1+k2+k3),ytot(2*nptstot+k1+k2+k3),ztot(2*nptstot+k1+k2+k3))

  xtot(1:2*nptstot)=x(1:2*nptstot)
  ytot(1:2*nptstot)=y(1:2*nptstot)
  ztot(1:2*nptstot)=z(1:2*nptstot)

  if (use_top) then
     xtot(2*nptstot+1:2*nptstot+k1)=real(xtop(1:k1))
     ytot(2*nptstot+1:2*nptstot+k1)=real(ytop(1:k1))
     ztot(2*nptstot+1:2*nptstot+k1)=real(ztop(1:k1))
  endif
  if (use_bot) then
     xtot(2*nptstot+k1+1:2*nptstot+k1+k2)=real(xbot(1:k2))
     ytot(2*nptstot+k1+1:2*nptstot+k1+k2)=real(ybot(1:k2))
     ztot(2*nptstot+k1+1:2*nptstot+k1+k2)=real(zbot(1:k2))
  endif
  if (use_meri) then
     xtot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3)=xmeri(1:k3)
     ytot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3)=ymeri(1:k3)
     ztot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3)=zmeri(1:k3)
  endif

  ! save mesh for kerner
  open(unit=99,file=trim(outdir)//'/SNAPS/mesh_tot.xyz')
  do i=1,2*nptstot+k1+k2
     write(99,*)xtot(i),ytot(i),ztot(i)
  enddo
  close(99)

  if (use_meri) then
     open(unit=99,file=trim(outdir)//'/SNAPS/mesh_meri.xyz')
     do i=2*nptstot+k1+k2+1,2*nptstot+k1+k2+k3
        write(99,*)xtot(i),ytot(i),ztot(i)
     enddo
     close(99)
  endif

  ! load azimuthal prefactors-------------------------------------------------------------
  if (nsim > 1) then
     allocate(longit_snap(2*nptstot+k1+k2))
     longit_snap(1:nptstot) = phi0
     longit_snap(nptstot+1:2*nptstot) = phi0+dphi
     do i=1,k1
        longit_snap(2*nptstot+i) = azi_phi_top(i)
     enddo
     do i=1,k2
        longit_snap(2*nptstot+k1+i) = azi_phi_bot(i)
     enddo
     do i=1,k3
        longit_snap(2*nptstot+k1+k2+i) = azi_phi_meri(i)
     enddo

     allocate(mij_snap(2*nptstot+k1+k2+k3,nsim,3))
     mij_snap= -1.E30
     call compute_radiation_prefactor(mij_snap,2*nptstot+k1+k2+k3,nsim,longit_snap)
     allocate(disptot_sum(2*nptstot+k1+k2+k3,3))
  endif

  ! FLUID REGION
  npts_fluid=nel_fluid*16
  npts_fluid_read=nel_fluid*9
  write(*,*)'points in fluid:',npts_fluid
  allocate(fluid_prefact(npts_fluid*nproc_mesh,2,nsim))

  do isim = 1,nsim
     ! load fluid prefactors-------------------------------------------------------------
     write(*,*) 'loading fluid prefactors...'
     do iproc=0, nproc_mesh-1
        call define_io_appendix(appmynum,iproc)
        open(unit=190,file=trim(simdir(isim))//'/Data/inv_rho_s_fluid_globsnaps_'//appmynum//'.dat')
        i=0
        do ii=1,npts_fluid_read/9
           do iii=1,9
              read(190,*)coord9(iii,1),coord9(iii,2)
           enddo
           fluid_prefact(iproc*npts_fluid+i+1:iproc*npts_fluid+i+2,:,isim) = coord9(1:2,:)
                                                                       ! (0,0), (0,npol/2)
           fluid_prefact(iproc*npts_fluid+i+3,:,isim) = coord9(5,:)    ! (npol/2,npol/2)
           fluid_prefact(iproc*npts_fluid+i+4,:,isim) = coord9(4,:)    ! (npol/2,0)
           fluid_prefact(iproc*npts_fluid+i+5:iproc*npts_fluid+i+6,:,isim) = coord9(2:3,:)
                                                                       ! (0,npol/2),(0,npol)
           fluid_prefact(iproc*npts_fluid+i+7,:,isim) = coord9(6,:)    ! (npol/2,npol)
           fluid_prefact(iproc*npts_fluid+i+8,:,isim) = coord9(5,:)    ! (npol/2,npol/2)
           fluid_prefact(iproc*npts_fluid+i+9:iproc*npts_fluid+i+10,:,isim) = coord9(4:5,:)
                                                                  ! (npol/2,0),(npol/2,npol/2)
           fluid_prefact(iproc*npts_fluid+i+11,:,isim) = coord9(8,:)   ! (npol,npol/2)
           fluid_prefact(iproc*npts_fluid+i+12,:,isim) = coord9(7,:)   ! (npol,0)
           fluid_prefact(iproc*npts_fluid+i+13:iproc*npts_fluid+i+14,:,isim) = coord9(5:6,:)
                                                                  ! (npol/2,npol/2),(npol/2,npol)
           fluid_prefact(iproc*npts_fluid+i+15,:,isim) = coord9(9,:)   ! (npol,npol)
           fluid_prefact(iproc*npts_fluid+i+16,:,isim) = coord9(8,:)   ! (npol,npol/2)
           i=i+16
        enddo  ! npts_fluid_read
        close(190)
     enddo ! nproc_mesh
  enddo !isim

  ! compute longitude and radiation pattern for multiple wavefields

  ! load snaps ===============================================================
  write(*,*)'loading snaps...'
  allocate(disp(2*nptstot+k1+k2+k3,3))

  do j=snap1,snap2,snapskip

     disp = 0.
     if (nsim > 1) disptot_sum = 0.

     do isim=1,nsim
        do iproc=0,nproc_mesh-1
           call define_io_appendix(appmynum,iproc)
           call define_io_appendix(appmynum2,j)
           open(unit=99,file=trim(simdir(isim))//'/Data/snap_'//appmynum//'_'//appmynum2//'.dat', &
                     FORM="UNFORMATTED",STATUS="OLD",POSITION="REWIND")
           i=0
           do ii=1,npts_read/9
              do iii=1,9
                 read(99) disp9(iii,1), disp9(iii,2), disp9(iii,3)
              enddo
              disp(iproc*npts+i+1:iproc*npts+i+2,:) = disp9(1:2,:)
                                                    ! (0,0), (0,npol/2)
              disp(iproc*npts+i+3,:) = disp9(5,:)   ! (npol/2,npol/2)
              disp(iproc*npts+i+4,:) = disp9(4,:)   ! (npol/2,0)
              disp(iproc*npts+i+5:iproc*npts+i+6,:) = disp9(2:3,:)
                                                    ! (0,npol/2),(0,npol)
              disp(iproc*npts+i+7,:) = disp9(6,:)   ! (npol/2,npol)
              disp(iproc*npts+i+8,:) = disp9(5,:)   ! (npol/2,npol/2)
              disp(iproc*npts+i+9:iproc*npts+i+10,:) = disp9(4:5,:)
                                                    ! (npol/2,0),(npol/2,npol/2)
              disp(iproc*npts+i+11,:) = disp9(8,:)  ! (npol,npol/2)
              disp(iproc*npts+i+12,:) = disp9(7,:)  ! (npol,0)
              disp(iproc*npts+i+13:iproc*npts+i+14,:) = disp9(5:6,:)
                                                    ! (npol/2,npol/2),(npol/2,npol)
              disp(iproc*npts+i+15,:) = disp9(9,:)  ! (npol,npol)
              disp(iproc*npts+i+16,:) = disp9(8,:)  ! (npol,npol/2)
              i=i+16
           enddo
           close(99)

           disp(iproc*npts+1:iproc*npts+npts_fluid,1)= &
                disp(iproc*npts+1:iproc*npts+npts_fluid,1) * fluid_prefact(:,1,isim)
           disp(iproc*npts+1:iproc*npts+npts_fluid,2)= &
                disp(iproc*npts+1:iproc*npts+npts_fluid,2) * fluid_prefact(:,1,isim) &
                    * fluid_prefact(:,2,isim)
           disp(iproc*npts+1:iproc*npts+npts_fluid,3)= &
                disp(iproc*npts+1:iproc*npts+npts_fluid,3) * fluid_prefact(:,1,isim)
        enddo

        disp(nptstot+1:2*nptstot,1)=disp(1:nptstot,1)
        disp(nptstot+1:2*nptstot,2)=disp(1:nptstot,2)
        disp(nptstot+1:2*nptstot,3)=disp(1:nptstot,3)

        do i=1,k1
           disp(2*nptstot+i,1) = disp(azi_ind_top(i),1)
           disp(2*nptstot+i,2) = disp(azi_ind_top(i),2)
           disp(2*nptstot+i,3) = disp(azi_ind_top(i),3)
        enddo
        do i=1,k2
           disp(2*nptstot+k1+i,1) = disp(azi_ind_bot(i),1)
           disp(2*nptstot+k1+i,2) = disp(azi_ind_bot(i),2)
           disp(2*nptstot+k1+i,3) = disp(azi_ind_bot(i),3)
        enddo
        if (use_meri) then
           do i=1,k3
              disp(2*nptstot+k1+k2+i,1) = disp(azi_ind_meri(i),1)
              disp(2*nptstot+k1+k2+i,2) = disp(azi_ind_meri(i),2)
              disp(2*nptstot+k1+k2+i,3) = disp(azi_ind_meri(i),3)
           enddo
        endif

        ! @TODO only plot summed stuff when summing to avoid large data. should
        ! be a parameter in param_post_processing at some point
        if (nsim == 1) then
           filename1 = trim(outdir)//'/SNAPS/snap_cell_'&
                       //trim(src_type(isim,2))//'_'//appmynum2//'_z'
           write(*,*)'filename out vtk :',filename1
           call write_VTK_bin_scal_topology(xtot(1:2*nptstot+k1+k2),ytot(1:2*nptstot+k1+k2), &
                                            ztot(1:2*nptstot+k1+k2), &
                                            disp(1:2*nptstot+k1+k2,3),(2*nptstot+k1+k2)/4, &
                                            filename1)
        endif


        if (use_meri) then
           filename1 = trim(simdir(isim))//'/'//trim(outdir)//'/SNAPS/meri_snap_'&
                        //trim(src_type(isim,2))//'_'//appmynum2//'_z'
           call write_VTK_bin_scal(xtot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                ytot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                ztot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                disp(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3,3),k3,0,filename1)
        endif

        if (nsim > 1) then
           disptot_sum(:,1) = disptot_sum(:,1) + mij_snap(:,isim,1)*disp(:,2)
           disptot_sum(:,2) = disptot_sum(:,2) + mij_snap(:,isim,2)*disp(:,2)
           disptot_sum(:,3) = disptot_sum(:,3) + mij_snap(:,isim,3)*disp(:,3)
        endif

     enddo

     if (nsim > 1) then

        filename1=trim(outdir)//'/SNAPS/snap_mij_cell_'//appmynum2//'_z'
        write(*,*)'filename out vtk :',filename1
        call write_VTK_bin_scal_topology(xtot(1:2*nptstot+k1+k2),ytot(1:2*nptstot+k1+k2), &
                                         ztot(1:2*nptstot+k1+k2), &
                                         disptot_sum(1:2*nptstot+k1+k2,3), &
                                         (2*nptstot+k1+k2)/4,filename1)

        if (use_meri) then
           filename1=trim(outdir)//'/SNAPS/meri_snap_mij_'//appmynum2//'_z'
           call write_VTK_bin_scal(xtot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                ytot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                ztot(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3), &
                disptot_sum(2*nptstot+k1+k2+1:2*nptstot+k1+k2+k3,3),k3,0,filename1)
        endif
     endif
  enddo

end subroutine compute_3d_wavefields
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine construct_surface_cubed_sphere(npts_surf, npts, rsurf, ind_proc_surf, &
                                          ind_pts_surf, nptstot, coord1, kpts, &
                                          dphi, phi0, in_or_out, n, xsurf, ysurf, &
                                          zsurf, azi_phi_surf, azi_ind_surf)

    !!!!! BEG CUBED SPHERE
    ! af2tnm: along a great circle, there's the equivalent of two chunks of the cubed
    ! sphere. According to the convention we use in defining the mapping in the cubed sphere
    ! module, that amounts to 2*nang spectral elements.
    ! Assuming we will be using npol_cs=1 in the following, nang has to be even.
    ! We therefore want nang=.5*(npts_surf-1) if npts_surf is odd
    ! We therefore want nang=.5*(npts_surf) if npts_surf is even

    use global_par
    implicit none

    integer, intent(in) :: npts_surf,npts,nptstot,n
    double precision, intent(in) :: rsurf,dphi,phi0
    character(len=7), intent(in) :: in_or_out
    integer, dimension(1:npts_surf), intent(in) :: ind_proc_surf,ind_pts_surf
    double precision, dimension(nptstot,2), intent(in) :: coord1
    double precision, dimension(1:n), intent(out) :: xsurf,ysurf,zsurf,azi_phi_surf
    integer, dimension(1:n), intent(out) :: azi_ind_surf
    integer, intent(out) :: kpts

    double precision, allocatable,dimension(:) :: xi_el,eta_el,r_el
    integer, allocatable, dimension(:,:,:,:) :: number_el
    double precision, allocatable,dimension(:) :: xi_i,xi_j,xi_k
    double precision, allocatable,dimension(:,:,:,:) :: xcol,ycol,zcol
    double precision, allocatable,dimension(:,:,:,:) :: x_el,y_el,z_el
    double precision :: dang,C,D,re,ri,teta,tksi,tr,Xe,Ye,delta
    integer :: npol_cs,ii,iii,izone,nang,nelt,nr,nrpol,iel,nel_surf,jj,ipol,jpol,kpol,j,i
    double precision ::  dist,r_ref,th_ref,th,dr,r,xc,yc,zc,ph

    write(*,*)'computing cubed sphere for surface at r=',rsurf
    write(*,*)'pts surf,total:',npts_surf,npts

    nang=floor(sqrt(real(n))/6./2.)
    write(*,*)'naang,nang,npts_surf:',n,nang,npts_surf

     nr=1 ! one radial layer only
     write(*,*)'NANG:',nang
     allocate(xi_el(0:nang),eta_el(0:nang),r_el(0:nr))
     xi_el = 0.d0;  eta_el = 0.d0 ; r_el = 0.d0
     re=rsurf
     ri=rsurf-1.d0 ! (in km)
     !
     dang=pi/(2.d0*dble(nang))
     dr=(re-ri)/dble(nr)
     do i=0,nang
      xi_el(i)=-pi*0.25d0+dang*dble(i)
      eta_el(i)=-pi*0.25d0+dang*dble(i)
     enddo
     do i=0,nr
      r_el(i)=ri+dr*dble(i)
     enddo
     allocate(x_el(0:nang,0:nang,0:nr,1:6))
     allocate(y_el(0:nang,0:nang,0:nr,1:6))
     allocate(z_el(0:nang,0:nang,0:nr,1:6))
     x_el = 0.d0; y_el=0.d0; z_el = 0.d0
     allocate(number_el(0:nang,0:nang,0:nr,1:6))
     number_el = 0
     do izone = 1, 6
        do iii=0,nr-1     ! loop over  r
           do ii=0,nang-1   ! loop over eta
              do i=0,nang-1   ! loop over ksi
                 number_el(i,ii,iii,izone) = (izone-1)*(nang*nang*nr)+&
                                             iii*(nang*nang)+((ii*nang)+i+1)
              enddo
           enddo
        enddo
     enddo
     nelt = maxval(number_el)
     nrpol=1
     npol_cs=1
     allocate(xi_i(0:npol_cs),xi_j(0:npol_cs),xi_k(0:nrpol))

     xi_i(0) = -1
     xi_i(1) = 1
     xi_j(0) = -1
     xi_j(1) = 1
     xi_k(0) = -1
     xi_k(1) = 1

     write(*,*)'Now define collocation points for the SEM-cubed sphere grid'
     allocate(xcol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
     allocate(ycol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
     allocate(zcol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
     kpts = 0
     nel_surf = 0

     write(*,*)'ZONE NANG:',nr,nang,npol_cs,nelt
     write(*,*)'ZONE NANG 2 nel_surf:',6*nr*nang**2

     do izone=1,6         ! loop over the chunks
      do iii=0,nr-1       ! loop over r   (the global r)
       do ii=0,nang-1     ! loop over eta (the global eta)
        do i=0,nang-1     ! loop over xi  (the global xi)
         iel = number_el(i,ii,iii,izone)
         nel_surf=nel_surf+1
         do kpol = 0, 0 ! TNM nrpol  ! loop over the elemental k index (r direction)
          do jpol = 0, npol_cs  ! loop over the elemental j index (eta direction)
           do ipol = 0, npol_cs ! loop over the elemental i index (xi direction)
            tksi= xi_el(  i) +(1.d0+ xi_i(ipol))*.5d0*dang
            teta=eta_el( ii) +(1.d0 + xi_j(jpol))*.5d0*dang
            tr=    r_el(iii) +(1.d0 + xi_k(kpol))*.5d0*dr
            Xe=tan(tksi)
            Ye=tan(teta)
            C=(1.d0+Xe**2)**(0.5d0)
            D=(1.d0+Ye**2)**(0.5d0)
            delta=1.d0+Xe**2+Ye**2
            kpts=kpts+1
            if (izone == 1) then
             Xcol(ipol,jpol,kpol,iel)=tr*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
            endif
            if (izone == 2) then
             Xcol(ipol,jpol,kpol,iel)=-tr*Xe*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=tr*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
            endif
            if (izone == 3) then
             Xcol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=-tr*Xe*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
            endif
            if (izone == 4) then
             Xcol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
            endif
            if (izone == 5) then
             Xcol(ipol,jpol,kpol,iel)=-tr*Ye*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=tr*delta**(-0.5d0)
            endif
            if (izone == 6) then
             Xcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
             Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
             Zcol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5d0)
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
     ! At this stage we know Xcol, Ycol, Zcol for the cubed sphere
     ! These arrays are four dimensional (ipol,jpol,kpol,iel)
     ! Their knowledge suffice to define the vtk output that will properly take
     ! into account the cubed sphere topology (see the paraview.f90 module)
     !!!!! END CUBED SPHERE

      write(*,*)'number of surface pts,elems,tot pts:',npts_surf,nel_surf,kpts
      write(*,*)'max ind_proc, ind_pts:',maxval(ind_proc_surf),maxval(ind_pts_surf)
      write(*,*)size(xsurf)
      xsurf = 0
      ysurf = 0
      zsurf = 0
      iii = 0
      write(*,*)'constructing 1d array for surface coordinates...'
      do iel = 1, nel_surf
         if ( mod(iel,floor(nel_surf/10.)) == 0  ) then
            write(*,*)'percentage done:',ceiling(real(iel)/real(nel_surf)*100.)
            call flush(6)
         endif
         xc = sum(xcol(:,:,0,iel)) / 4.d0
         yc = sum(ycol(:,:,0,iel)) / 4.d0
         zc = sum(zcol(:,:,0,iel)) / 4.d0

         call xyz2rthetaphi(r, th, ph, xc, yc, zc)

         if ( (in_or_out == 'innside' .and. ph >= phi0 .and. ph <= phi0+dphi) .or. &
              (in_or_out == 'outside' .and. (ph <= phi0 .or. ph >= phi0+dphi) ) ) then

            xsurf(iii+1) = xcol(0,0,0,iel)
            ysurf(iii+1) = ycol(0,0,0,iel)
            zsurf(iii+1) = zcol(0,0,0,iel)

            xsurf(iii+2) = xcol(npol_cs,0,0,iel)
            ysurf(iii+2) = ycol(npol_cs,0,0,iel)
            zsurf(iii+2) = zcol(npol_cs,0,0,iel)

            xsurf(iii+3) = xcol(npol_cs,npol_cs,0,iel)
            ysurf(iii+3) = ycol(npol_cs,npol_cs,0,iel)
            zsurf(iii+3) = zcol(npol_cs,npol_cs,0,iel)

            xsurf(iii+4) = xcol(0,npol_cs,0,iel)
            ysurf(iii+4) = ycol(0,npol_cs,0,iel)
            zsurf(iii+4) = zcol(0,npol_cs,0,iel)

            ! determine the corresponding point in the D-shape domain
            do jj=1,4
               call xyz2rthetaphi(r,th, azi_phi_surf(iii+jj),xsurf(iii+jj),ysurf(iii+jj),zsurf(iii+jj))
               dist=2.d0*pi
               do i=1,npts_surf
                  call get_r_theta(coord1(ind_proc_surf(i)*npts+ind_pts_surf(i),1), &
                       coord1(ind_proc_surf(i)*npts+ind_pts_surf(i),2),r_ref,th_ref)
                  if (abs(th-th_ref) < dist) then
                     dist=abs(th-th_ref)
                     azi_ind_surf(iii+jj)=ind_proc_surf(i)*npts+ind_pts_surf(i)
                  endif
               enddo
            enddo

            iii = iii + 4

            endif ! in_or_out

         enddo

         kpts=iii
         write(*,*)'total points in surface:',kpts

      write(*,*)'done with cubed sphere for r=',rsurf
      deallocate(xcol,ycol,zcol,x_el,y_el,z_el)

end subroutine construct_surface_cubed_sphere
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine sphi2xy(x,y,s,phi,n)
  integer, intent(in) :: n
  real, dimension(1:n), intent(out) :: x,y
  real, dimension(1:n), intent(in) :: s
  real, intent(in) :: phi

  x(1:n)=s(1:n)*cos(phi)
  y(1:n)=s(1:n)*sin(phi)

end subroutine sphi2xy
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_VTK_bin_scal(x,y,z,u1,rows,nelem_disk,filename1)
  use data_all
  implicit none
  integer :: t,rows,nelem_disk,ioerr,dims
  real, dimension(1:rows), intent(in) :: x,y,z,u1
  integer, dimension(1:rows,2) :: cell
  integer, dimension(1:rows) :: cell_type
  real, dimension(1:rows,3) :: W

  character (len=30) :: celltype;
  character (len=100) :: filename1;
  character (len=50) :: ss; !stream

  W(1:rows,1)=x
  W(1:rows,2)=y
  W(1:rows,3)=z
  !points structure
  do i=1,rows
   cell(i,1)=1
   cell(i,2)=i
   cell_type(i)=1
  enddo

   write(*,*)'computing VTK bin file ',trim(filename1)//'.vtk  ...'

#if defined(__gfortran__)
  open(100,file=trim(filename1)//'.vtk',access='stream', &
                          status='replace',convert='big_endian')
#elif defined(__INTEL_COMPILER)
  open(100,file=trim(filename1)//'.vtk',access='stream', &
                          status='replace',convert='big_endian')
#else
  open(100,file=trim(filename1)//'.vtk',access='stream', &
                          status='replace')
#endif

  write(100) '# vtk DataFile Version 3.0'//char(10)
  write(100) 'Cell Fractions'//char(10)
  write(100) 'BINARY'//char(10)
  write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
  write(ss,fmt='(A8,I12,A10)') 'POINTS',rows,' float'
  write(100) ss//char(10)
  !points
  do i=1,rows
  write(100) W(i,1:3)
  enddo
  write(100) char(10)
  !cell topology
  write(ss,fmt='(A5,2I12)') 'CELLS',rows,rows*2
  write(100) char(10)//ss//char(10)
  do i=1,rows
  write(100) cell(i,1:2)
  enddo
  write(100) char(10)
  !cell type
  write(ss,fmt='(A10,2I12)') 'CELL_TYPES',rows
  write(100) char(10)//ss//char(10)
  do i=1,rows
  write(100) cell_type(i)
  enddo
  write(100) char(10)
  !data
  write(ss,fmt='(A10,I12)') 'CELL_DATA',rows
  write(100) char(10)//ss//char(10)
  write(100) 'SCALARS Displ_u1 float 1'//char(10)
  write(100) 'LOOKUP_TABLE default'//char(10) !color table?
  do i=1,rows
  write(100) u1(i)
  enddo
   close(100)
  write(*,*)'...saved ',trim(outdir)//'/'//trim(filename1)//'.vtk'
end subroutine write_vtk_bin_scal
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine write_VTK_bin_scal_topology(x,y,z,u1,elems,filename)
  implicit none
  integer*4 :: i,t,elems
  real*4, dimension(1:elems*4), intent(in) :: x,y,z,u1
  integer*4, dimension(1:elems*5) :: cell
  integer*4, dimension(1:elems) :: cell_type
  character (len=200) :: filename
  character (len=50) :: ss; !stream
  !points structure

  do i=5,elems*5,5
   cell(i-4)=4;
   enddo
  t=0
  do i=5,elems*5,5
  t=t+4
  cell(i-3)=t-4;
  cell(i-2)=t-3;
  cell(i-1)=t-2;
  cell(i)=t-1;
  enddo

  !do i=1,elems
  ! cell_type(i)=9
  !enddo
  cell_type=9
  ! write(*,*)'computing vtk file ',trim(filename),' ...'
#if defined(__gfortran__)
  open(100,file=trim(filename)//'.vtk',access='stream',status='replace',convert='big_endian')
#elif defined(__INTEL_COMPILER)
  open(100,file=trim(filename)//'.vtk',access='stream',status='replace',convert='big_endian')
#else
  open(100,file=trim(filename)//'.vtk',access='stream',status='replace')
#endif

  write(100) '# vtk DataFile Version 4.0'//char(10)
  write(100) 'mittico'//char(10)
  write(100) 'BINARY'//char(10)
  write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
  write(ss,fmt='(A6,I10,A5)') 'POINTS',elems*4,'float'
  write(100) ss//char(10)
  !points
  write(100) (x(i),y(i),z(i),i=1,elems*4)
  write(100) char(10)
  !cell topology
  write(ss,fmt='(A5,2I10)') 'CELLS',elems,elems*5
  write(100) char(10)//ss//char(10)
  write(100) cell
  write(100) char(10)
  !cell type
  write(ss,fmt='(A10,2I10)') 'CELL_TYPES',elems
  write(100) char(10)//ss//char(10)
  write(100) cell_type
  write(100) char(10)
  !data
  write(ss,fmt='(A10,I10)') 'POINT_DATA',elems*4
  write(100) char(10)//ss//char(10)
  write(100) 'SCALARS data float 1'//char(10)
  write(100) 'LOOKUP_TABLE default'//char(10) !color table?
  write(100) u1
   close(100)
  write(*,*)'...saved ',trim(filename)//'.vtk'
end subroutine write_vtk_bin_scal_topology
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine rthetaphi2xyz(x,y,z,r,theta,phi)
  real, intent(out) :: x,y,z
  real, intent(in) :: r,theta,phi

  x = r * sin(theta) * cos(phi)
  y = r * sin(theta) * sin(phi)
  z = r * cos(theta)

end subroutine rthetaphi2xyz
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
real function prem(r0,param)
  ! prem model in terms of domains separated by discontinuities
  ! MvD: - at discontinuities, upper or lower domain is chosen based on numerical
  !        rounding errors
  !      - only used for mesh plots, nor physical meaning


  real, intent(in) :: r0
  real            :: r,x_prem
  real             :: ro_prem,vp_prem,vs_prem
  character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.

  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  if (r > 6356. .and. r <= 6371.01) then        ! upper crustal layer
     ro_prem=2.6
     vp_prem=5.8
     vs_prem=3.2
  else if (r > 6346.6 .and. r <= 6356.) then
     ro_prem=2.9                       ! lower crustal layer
     vp_prem=6.8
     vs_prem=3.9
  else if (r > 6151. .and. r <= 6346.6) then
     ro_prem=2.691+.6924*x_prem             ! upper mantle
     vp_prem=4.1875+3.9382*x_prem
     vs_prem=2.1519+2.3481*x_prem
  else if (r > 5971. .and. r <= 6151. ) then
     ro_prem=7.1089-3.8045*x_prem
     vp_prem=20.3926-12.2569*x_prem
     vs_prem=8.9496-4.4597*x_prem
  else if (r > 5771. .and. r <= 5971.) then
     ro_prem=11.2494-8.0298*x_prem
     vp_prem=39.7027-32.6166*x_prem
     vs_prem=22.3512-18.5856*x_prem
  else if (r > 5701. .and. r <= 5771. ) then
     ro_prem=5.3197-1.4836*x_prem
     vp_prem=19.0957-9.8672*x_prem
     vs_prem=9.9839-4.9324*x_prem
  else if (r > 5600. .and. r <= 5701. ) then   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
  else if (r > 3630. .and. r <= 5600. ) then
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vs_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
  else if (r > 3480. .and. r <= 3630.) then
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+.9783*x_prem**3
  else if (r > 1221.5 .and. r <= 3480. ) then  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vp_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vs_prem=0.00
  else if (r <= 1221.5) then                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vp_prem=11.2622-6.3640*x_prem**2
     vs_prem=3.6678-4.4475*x_prem**2
  ELSE
     write(*,*)'wrong radius!',r
     stop 2
  endif

  if (param == 'rho') then
     prem=ro_prem*1000.
  else if (param == 'v_p') then
     prem=vp_prem*1000.
  else if (param == 'v_s') then
     prem=vs_prem*1000.
  else
     write(*,*)'ERROR IN PREM function:',param,'NOT AN OPTION'
     stop 2
  endif

end function prem
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine xyz2rthetaphi(r,theta,phi,x,y,z)
  use global_par
  ! Alex2TNM:
  ! This one you might find useful for other purposes
  ! TNM@Alex: Indeed, done.
  double precision, intent(out) :: r,theta,phi
  double precision, intent(in) :: x,y,z

  r = dsqrt(x**2+y**2+z**2)
  theta = .5d0*pi-dasin(z/r)
  if (y > 0) then
    if (x > 0) then
      phi = datan(y/(x+1.e-20))
    else
      phi = pi+datan(y/(x+1.e-20))
    endif
  else
    if (x > 0) then
      phi = 2*pi+datan(y/(x+1.e-20))
    else
      phi = pi+datan(y/(x+1.e-20))
    endif
  endif
   if (abs(x) < 1.e-20) then
     if (y > 0.d0) then
       phi = .5d0*pi
     else
       phi = 1.5d0*pi
    endif
  endif
end subroutine xyz2rthetaphi
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_r_theta(s,z,r,th)
  use global_par
  double precision, intent(in)  :: s, z
  double precision, intent(out) :: r, th

  th = datan(s / (z + epsi))

  if ( 0.d0 > th ) th = pi + th
  if (th == zero .and. z < 0.d0) th = pi

  r = dsqrt(s**2 + z**2)

end subroutine get_r_theta
!-----------------------------------------------------------------------------------------
