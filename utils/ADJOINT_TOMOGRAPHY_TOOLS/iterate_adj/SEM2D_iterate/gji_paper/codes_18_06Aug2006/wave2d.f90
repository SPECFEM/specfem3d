program wave2d

  use wave2d_variables  ! global variables
  use wave2d_solver     !
  use wave2d_sub        ! source time function, write seismograms, etc
  use wave2d_sub2       ! UTM conversion
  use wave2d_sub3       ! measurements, including multitaper

  implicit none

! wavefields and seismograms
  double precision, dimension(:,:,:), allocatable :: samp, samp_dat, data, syn, adj_syn

! source time function
  double precision, dimension(NSTEP) :: stf_dat, stf

! socal coast and shelf data
  integer :: UTM_PROJECTION_ZONE, FINITE_SOURCE
  integer :: nshelf,ncoast
  double precision :: d,dmin
  double precision, dimension(:), allocatable :: shelf_lat,shelf_lon,shelf_z,shelf_x
  double precision, dimension(:), allocatable :: coast_lat,coast_lon,coast_z,coast_x
  character(len=200) :: shelf_file,coast_file

  integer :: itemp,itemp1,itemp2,itemp3

! index vectors of the target points within the designated region (NOT gridpoint indices)
!  integer, dimension(:), allocatable :: ifilter_src

!------------------------
! source and receiver variables (and fake receivers for spectral plots)
  integer :: istf, nsrc, nrec, nrecf, nevent, nevent_dat, nevent_run
  !integer, dimension(:), allocatable :: sglob, sglob_dat
  !integer, dimension(:), allocatable :: sglob_temp, rglob_temp, fglob_temp

  integer, dimension(MAX_EVENT):: ifilter_eve_dat
  double precision, dimension(MAX_EVENT)      :: x_eve_lon0_dat, z_eve_lat0_dat, x_eve0_dat, z_eve0_dat
  double precision, dimension(:), allocatable :: x_eve_lon_dat,  z_eve_lat_dat,  x_eve_dat,  z_eve_dat, otime_dat
  double precision, dimension(:), allocatable :: xi_eve_dat, gamma_eve_dat
  integer, dimension(:), allocatable          :: eglob_dat, ispec_eve_dat

  integer, dimension(MAX_EVENT):: ifilter_eve
  double precision, dimension(MAX_EVENT)      :: x_eve_lon0, z_eve_lat0, x_eve0, z_eve0
  double precision, dimension(:), allocatable :: x_eve_lon,  z_eve_lat,  x_eve,  z_eve, otime
  double precision, dimension(:), allocatable :: xi_eve, gamma_eve
  integer, dimension(:), allocatable          :: eglob, ispec_eve

  integer, dimension(MAX_SR):: ifilter_rec, rglob0
  double precision, dimension(MAX_SR)         :: x_rec_lon0, z_rec_lat0, x_rec0, z_rec0
  double precision,dimension(:), allocatable  :: x_rec_lon,  z_rec_lat,  x_rec,  z_rec
  double precision, dimension(:), allocatable :: xi_rec, gamma_rec
  integer, dimension(:), allocatable          :: rglob, ispec_rec

  integer, dimension(MAX_SR_FAKE):: ifilter_recf
  double precision, dimension(MAX_SR_FAKE)    :: x_recf0, z_recf0
  double precision,dimension(:), allocatable  :: x_recf, z_recf
  integer, dimension(:), allocatable          :: fglob

  integer, dimension(MAX_SR):: ifilter_src
  double precision, dimension(MAX_SR)         :: x_src_lon0, z_src_lat0, x_src0, z_src0
  double precision,dimension(:), allocatable  :: x_src_lon,  z_src_lat,  x_src,  z_src
  double precision, dimension(:), allocatable :: xi_src, gamma_src
  integer, dimension(:), allocatable          :: sglob, ispec_src

  integer, dimension(MAX_SR):: ifilter_src_dat
  double precision, dimension(MAX_SR)         :: x_src_lon0_dat, z_src_lat0_dat, x_src0_dat, z_src0_dat
  double precision,dimension(:), allocatable  :: x_src_lon_dat,  z_src_lat_dat,  x_src_dat,  z_src_dat
  double precision, dimension(:), allocatable :: xi_src_dat, gamma_src_dat
  integer, dimension(:), allocatable          :: sglob_dat, ispec_src_dat
!------------------------

! allocate Lagrange interpolators for sources and receivers
  double precision, dimension(:,:), allocatable :: hxir_store, hxie_store, hxied_store, hxis_store, hxisd_store
  double precision, dimension(:,:), allocatable :: hgammar_store, hgammae_store, hgammaed_store, hgammas_store, hgammasd_store
  double precision, dimension(:), allocatable :: hxir,hpxir,hgammar,hpgammar
  double precision, dimension(:), allocatable :: hxie,hpxie,hgammae,hpgammae
  double precision, dimension(:), allocatable :: hxis,hpxis,hgammas,hpgammas
  double precision, dimension(:), allocatable :: hxied,hpxied,hgammaed,hpgammaed
  double precision, dimension(:), allocatable :: hxisd,hpxisd,hgammasd,hpgammasd

  !double precision, dimension(:), allocatable :: x_src_lon, z_src_lat, x_src, z_src

  double precision :: x_src_lon_i, z_src_lat_i, x_src_i, z_src_i, x_src_f, z_src_f
  double precision :: flen,finc,src_az,xd
  double precision :: s_radius,xcen,zcen,xcen_lon,zcen_lat
  double precision :: dmin_trsh_s,dmin_trsh_r,dmin_trsh_f
  character(len=200) :: quake_file

! source function
  double precision :: f0(NCOMP), ti(NSTEP)

! vewlocity structure variables
  integer :: Nfac, IMODEL
  double precision :: w_scale, afac
  double precision, dimension(NGLOB) :: c_glob_dat, c_glob_syn

! source perturbations (01-Feb-2006)
! presently only set up for 2D surface waves
  logical :: ISOURCE_LOG
  integer :: PERT_SOURCE, PERT_STRUCT, INV_SOURCE, INV_STRUCT
  double precision :: src_pert_time, origin_time, origin_time_dat, ptmin, ptmax, ptime
  double precision :: src_pert_dist, rand1, amin, amax, pangle, rmin, rmax, prad
  double precision :: dx_src, dz_src, dt_src
  double precision :: m_scale(3)
  double precision, dimension(:,:,:,:), allocatable :: three_source_model
  double precision, dimension(:), allocatable :: source_gradient, m_scale_src, gradient

  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: gradient_local

! socal earthquake location
  double precision, dimension(MAX_EVENT) :: socal_events_lon, socal_events_lat, socal_events_mag
  integer, dimension(MAX_EVENT) :: socal_events_lab
  integer nmeas, ievent_min, ievent_max, ifirst

!----------------------
! misfit and conjugate gradient algorithm

! misfit
  double precision, dimension(MAX_EVENT) :: chi_etot

! smoothing
  integer :: igaus
  double precision :: dist2, dtrsh2, xtar, ztar, gamma  ! sigma
  double precision, dimension(NGLOB) :: k_rough_global, k_smooth_global
  double precision, dimension(NGLOB) :: k_gaus_global, k_gaus_global_ex, k_gaus_int_global
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: k_temp
  !double precision, dimension(NLOCAL) :: gk_local, m0_local, mt_local, mt2_local, da_local
  !integer, dimension(NLOCAL) :: ibool2
  character(len=200) :: file_smooth

! conjugate gradient optimization (using line-search)
!  double precision, dimension(NGLOB) :: pk, gk, g0, p0, m0, m0_vel, gt, mt, mt_vel, m_vel
  double precision, dimension(:), allocatable :: pk, gk, g0, p0, m0, gt, mt

  double precision, dimension(:), allocatable :: m0_vec, mt_vec

!  double precision, dimension(:), allocatable :: m_vel, m0_vel, mt_vel
  double precision, dimension(:), allocatable :: m_src, m0_src, mt_src, m_src_dat, m_src_syn
  double precision :: beta_val, lam_0_val, lam_t_val, chi_val, chi_t_val, chi_k_val
  double precision :: Pa,Pb,Pc,Pd,qfac,xx1,xx2,yy1,yy2,g1,g2,xmin
  integer :: itest, nmod, nmod_src, nmod_str

!----------------------

! kernels
  double precision, dimension(:), allocatable :: tstart, tend
  character(len=200) :: last_frame_name
  double precision, dimension(:,:,:,:,:), allocatable :: absorb_field
  double precision, dimension(:), allocatable :: rho_kernel, mu_kernel, kappa_kernel

  !double precision :: t0, ft(NCOMP), pxpf(NCOMP), f1(NCOMP), f2(NCOMP), dxdf(NCOMP)
  !double precision :: p_1, p_2, df, phi
  !double precision :: source_time_function
  !character(len=200) :: snap_name, last_frame_name
  !double precision :: xold(NCOMP), fold, gold(NCOMP), pold(NCOMP), xnew(NCOMP), fnew, gnew(NCOMP), pnew(NCOMP)
  !double precision :: lam, beta, eps, tstart, tend, stf

  double precision :: xtemp,ztemp,dx,dz,xmesh,zmesh
  double precision :: t_target, tcen
  double precision :: junk1,junk2,junk3
  double precision :: temp1,temp2,temp3,temp4,temp5

! corners of the grid
  double precision, dimension(4) :: corners_utm_x, corners_utm_z, corners_lon, corners_lat
  double precision :: sq_south, sq_east, sq_north, sq_west, sq_fac, da_mean, sda_mean, m_scale_str, mfac

  character(len=20)  :: data_tag,syn_tag,stf_tag,stfadj_tag
  character(len=200) :: srcfile,recfile,socal_map
  character(len=200) :: socal_dir, model_dir, script_dir, out_dir3, out_dir2, out_dir1, ev_dir
  character(len=200) :: filename, filename1, filename2, file_dat_c, file_syn_c, file_src_rec, file_stf
  character(len=200) :: command1

  integer :: i, j, k, iq, itime, iglob, irec, isrc, ipick, ievent, icomp, irunz
  integer :: isolver, irun0, irun, idat, iopt, ispec, istep, istep_switch, imod

  !double precision :: meas_pert, rand1,  meas, ppert

  !********* PROGRAM STARTS HERE *********************

  call random_seed()        ! random number generator

  out_dir3   = "OUTPUT/"
  in_dir     = "INPUT/"
  script_dir = "scripts/"

  socal_dir = '/home/carltape/socal_modes/local_modes/'
  model_dir = 'model_files/'

  ! exclude sources or receivers that are this close to THE SOCAL COAST
  dmin_trsh_s = 0.d+03                  ! sources
  !dmin_trsh_r = STATION_GRID_BUFFER     ! receivers
  dmin_trsh_r = 0.d+03
  dmin_trsh_f = 0.d+03                  ! fake receivers

  UTM_PROJECTION_ZONE = 11

  ! idat=1 means we run simulations for both data and synthetics (two models)
  idat = 0
  if (IKER <= 4) idat = 1
  if (ISRC_SPACE == 1) FINITE_SOURCE = 0
  if (ISRC_SPACE /= 1) FINITE_SOURCE = 1

!--------------------------------------
! load socal coast and shelf points

  ! read in data files
  nshelf = 101
  ncoast = 810
  allocate(shelf_lat(nshelf),shelf_lon(nshelf),shelf_z(nshelf),shelf_x(nshelf))
  allocate(coast_lat(ncoast),coast_lon(ncoast),coast_z(ncoast),coast_x(ncoast))
  shelf_file = 'oms_shelf'
  coast_file = 'oms_coast'
  open(77,file=trim(in_dir)//trim(shelf_file),status='unknown')
  read(77,*) (shelf_lat(i),shelf_lon(i),i=1,nshelf)
  close(77)
  open(66,file=trim(in_dir)//trim(coast_file),status='unknown')
  read(66,*) (coast_lat(i),coast_lon(i),i=1,ncoast)
  close(66)

  ! convert lat-lon to mesh coordinates
  call mesh_geo(nshelf,shelf_lon,shelf_lat,shelf_x,shelf_z,UTM_PROJECTION_ZONE,ILONLAT2MESH)
  call mesh_geo(ncoast,coast_lon,coast_lat,coast_x,coast_z,UTM_PROJECTION_ZONE,ILONLAT2MESH)
  !write(*,'(4f16.6)') (shelf_lon(i), shelf_lat(i), shelf_x(i)/1000., shelf_z(i)/1000., i=1,nshelf)
  !write(*,'(4f16.6)') (coast_lon(i), coast_lat(i), coast_x(i)/1000., coast_z(i)/1000., i=1,ncoast)

!--------------------------------------

  ! events to use (5-5 or 1-25)
  ievent_min = 1
  ievent_max = 25
  nevent_run = ievent_max - ievent_min + 1

  ! perturbation and inversion parameters
  PERT_SOURCE = 0
  PERT_STRUCT = 1
   INV_SOURCE = 0
   INV_STRUCT = 1

  ISOURCE_LOG = .false.        ! log file showing source loactions

  ! MAX source perturbations for location (lat-lon) and origin time
  src_pert_dist = 5000.0      ! source spatial perturbation (m)
  src_pert_time = -1.0        ! source temporal perturbation (s)
                              ! dt0 = dt0_dat - dt0_syn ( < 0 for delayed synthetics)
  origin_time_dat = tshift    ! data are unperturbed

  !src_pert_dist = 0.
  !src_pert_time = 0.

  ! reference directory for the run
  irunz = 4300

!============================================
! LOOP 1: different tomographic runs
  do iq = 1,1
!============================================

    ! KEY COMMANDS
    IMODEL = 2            ! (0) het map, (1) homo map, (2) checkerboard, (3) read in
    Nfac   = 1            ! scalelength of checker
    gamma  = 15.0d+03     ! scalelength of smoothing Gaussian

    !mfac = 1. - 0.1*(iq-1)

  irun0 = irunz + 20*(iq-1)

  ! name and create the reference output directory for the optimization run
  write(out_dir2,'(a,i4.4,a)') trim(out_dir3)//'run_',irun0,'/'
  command1 = 'mkdir ' // trim(out_dir2)
  call system(command1)
  !open(19,file='temp.csh',status='unknown')
  !write(19,'(2a)') 'mkdir ',out_dir2
  !close(19)
  !call system('chmod 755 temp.csh ; temp.csh')

!--------------------------------------

  print *
  print *, ' max time-step is ', NSTEP
  print *, '            DT is ', sngl(DT)
  print *, '      max time is ', sngl(NSTEP*DT), ' s,', sngl(NSTEP*DT/60), ' min'
  print *

  ! copy wave2d_constants.f90 into output directory
  command1 = 'cp ' // 'src/wave2d_constants.f90 ' // trim(out_dir2)
  call system(command1)

  !call write_parameters(trim(out_dir2)//'parameters1.log')

!--------------------------------------
! mesher

  ! set up model (gets Jacobian)
  call mesher()

  ! da vector and valence vector
  da(:) = 0.
  valence(:) = 0
  do ispec = 1,NSPEC
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob     = ibool(i,j,ispec)
        !da_local(iglob) = wxgll(i)*wzgll(j)*jacobian(i,j,ispec)

        da(iglob) = da(iglob) + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
        valence(iglob) = valence(iglob) + 1
      enddo
    enddo
  enddo
  !da(:) = da_local(:)*valence(:)
  da_mean = sum(da)/NGLOB
  sda_mean = sqrt(da_mean)

  ! global mesh
  open(unit=15,file = trim(out_dir2)//'global_mesh.dat',status='unknown')
  do iglob = 1,NGLOB
    write(15,'(2i8,3e18.8)') iglob, valence(iglob), x(iglob), z(iglob), da(iglob)
  enddo
  close(15)

  if (0 == 1) then

    ! corner points for each element, and centerpoint (in km)
    open(unit=15,file='elements.dat',status='unknown')
    do ispec = 1,nspec
      xtemp = (x1(ispec) + x2(ispec))/2.
      ztemp = (z1(ispec) + z2(ispec))/2.
      write(15,'(i8,6f14.6)') ispec,xtemp/1000.,ztemp/1000,x1(ispec)/1000.,x2(ispec)/1000.,z1(ispec)/1000.,z2(ispec)/1000.
    enddo
    close(15)

    ! GLL points for one element (in km)
    ispec = 292  ! pick an element
    open(unit=15,file='gll_points.dat',status='unknown')
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        write(15,'(4i10,3e18.8)') ispec, i, j, valence(iglob), da(iglob)/1e6, x(iglob)/1000., z(iglob)/1000.
      enddo
    enddo
    close(15)

    stop 'testing'
  endif

  print *
  write(*,'(a,3f20.10)') '  da(min,mean,max) : ', minval(da), da_mean, maxval(da)
  write(*,*)             '        sum [ da ] : ', sum(da)
  write(*,*)             '   LENGTH * HEIGHT : ', LENGTH * HEIGHT
  print *
  write(*,*)             '             gamma : ', gamma
  write(*,*)             '              Nfac : ', Nfac
  write(*,*)             '             irun0 : ', irun0
  write(*,*)             '              IKER : ', IKER
  print *
  print *, ' GLL weights:'
  do i=1,NGLLX
     print *, wxgll(i)
  enddo
  do i=1,NGLLZ
     print *, wzgll(i)
  enddo

!--------------------------------------
! determine the UTM coordinates of your origin and corners

  call utm_geo(LON_MIN,LAT_MIN,utm_xmin,utm_zmin,UTM_PROJECTION_ZONE,ILONGLAT2UTM)
  corners_utm_x(1) = utm_xmin
  corners_utm_z(1) = utm_zmin
  corners_utm_x(2) = utm_xmin + LENGTH
  corners_utm_z(2) = utm_zmin
  corners_utm_x(3) = utm_xmin + LENGTH
  corners_utm_z(3) = utm_zmin + HEIGHT
  corners_utm_x(4) = utm_xmin
  corners_utm_z(4) = utm_zmin + HEIGHT

  print *
  print *,'origin (LON_MIN, LAT_MIN):'
  print *, LON_MIN, LAT_MIN
  print *, '        lon             lat           utm_x            utm_z'
  do i=1,4
     call utm_geo(xtemp,ztemp,corners_utm_x(i),corners_utm_z(i),UTM_PROJECTION_ZONE,IUTM2LONGLAT)
     corners_lon(i) = xtemp
     corners_lat(i) = ztemp
     write(*,'(2f16.6,2e16.6)') corners_lon(i), corners_lat(i), corners_utm_x(i), corners_utm_z(i)
  enddo

  ! convert global gridpoint mesh coordinates to lat-lon
  x_lon(:) = 0.
  z_lat(:) = 0.
  call mesh_geo(NGLOB,x_lon,z_lat,x,z,UTM_PROJECTION_ZONE,IMESH2LONLAT)
  !write(*,'(2f16.6)') (x_lon(iglob), z_lat(iglob), iglob=1,NGLOB)

  ! compute min/max lat-lon
  sq_east  = maxval(x_lon(:))
  sq_west  = minval(x_lon(:))
  sq_north = maxval(z_lat(:))
  sq_south = minval(z_lat(:))

  print *
  print *, 'grid fits within the spherical square patch with these boundaries:'
  print *, '  west  : ', sq_west
  print *, '  east  : ', sq_east
  print *, '  south : ', sq_south
  print *, '  north : ', sq_north

  ! determine the UTM coordinates of your origin
  print *
  print *, 'UTM check for origin of mesh:'
  print *, LON_MIN,LAT_MIN
  call utm_geo(LON_MIN,LAT_MIN,utm_xmin,utm_zmin,UTM_PROJECTION_ZONE,ILONGLAT2UTM)
  print *, utm_xmin, utm_zmin
  call utm_geo(xtemp,ztemp,utm_xmin,utm_zmin,UTM_PROJECTION_ZONE,IUTM2LONGLAT)
  print *, xtemp,ztemp   ! should match LON_MIN,LAT_MIN
  print *

!--------------------------------------
! phase velocity model

  t_target = 2*hdur     ! target period for phase velocity map

  if (IMODEL <= 1) then

     ! write lat-lon gridpoints to file
     filename = trim(socal_dir) // 'socal_input.dat'
     print *, 'Output data file is ',trim(filename)
     open(unit=15,file=filename,status='unknown')
     write(15,*) t_target
     write(15,*) IMODEL       ! formerly IHOMO
     write(15,*) UTM_PROJECTION_ZONE
     write(15,*) NGLOB
     write(15,*) (x_lon(iglob),z_lat(iglob),iglob=1,NGLOB)
     close(15)

     ! run wave2d_socal.f90 to create output file
     ! KEY NOTE: must run this each time you change the grid, hdur, or IMODEL
     call system(trim(model_dir) // 'get_socal_map.csh')

     ! read in phase velocity map
     !write(infile,'(a,a,i4.4,a)') trim(socal_dir), 'socal_T', int(100*t_target),'.dat'
     socal_map = 'socaltest.dat'
     open(unit=16, file = trim(model_dir)//socal_map, status='unknown')
     read(16,*) (junk1,junk2,c_glob(iglob),junk3,iglob=1,NGLOB)
     close(16)

     ! read in reference phase velocity (km/s)
     filename = trim(model_dir) // 'socaltest2.dat'
     open(unit=17,file=filename,status='unknown')
     read(17,*) c0
     close(17)

     ! convert km/s --> m/s
     c_glob(:) = 1000.*c_glob(:)
     c0 = c0*1000.

     c_glob_syn(:) = c0       ! c-maps for synthetics

  else if (IMODEL == 2) then   ! checkerboard

     c0 = 3500.    ! reference phase velocity (m/s) (should be obtained based on hdur)
     !Nfac = 3      ! scalelength of map = N * (wavelength of surface waves for hdur)
     afac = 10.0   ! max PERCENT variation of synthetic model  (10.0)

     ! spatial frequency corresponding to scalelength
     ! factor of 2 is because the scalelength is a 1/2 cycle
     w_scale = 2*pi/(Nfac*2*c0*t_target)

     do iglob=1,NGLOB
        c_glob(iglob) = c0 + c0*afac/100.*(sin(x(iglob)*w_scale) * sin(z(iglob)*w_scale))
        !c_glob(iglob) = c0 + c0*afac/100.
     enddo

     c_glob_syn(:) = c0       ! c-maps for synthetics

  else if (IMODEL == 3) then   ! read c-map for synthetics and data for a 'middle' iteration

     ! read in phase velocity map for data
     open(unit=16,file=trim(model_dir)//'socal_vel_dat.dat', status='unknown')
     read(16,*) (junk1,junk2,c_glob(iglob),junk3,iglob=1,NGLOB)
     close(16)

!!$     ! read in phase velocity map for synthetics
!!$     open(unit=16,file=trim(model_dir)//'socal_vel_syn.dat', status='unknown')
!!$     read(16,*) (junk1,junk2,c_glob_syn(iglob),junk3,iglob=1,NGLOB)
!!$     close(16)
!!$
!!$     ! read in c0
!!$     open(unit=16,file=trim(model_dir)//'socal_vel_c0.dat', status='unknown')
!!$     read(16,*) c0
!!$     close(16)

     ! compute c0 by finding the average of the phase velocity map
     c0 = 0.
     do iglob = 1,NGLOB
        c0 = c0 + c_glob(iglob) * da(iglob)
     enddo
     c0 = c0 / (LENGTH * HEIGHT)

     c_glob_syn(:) = c0

  else
     stop 'IMODEL must be 0,1,2,3'
  endif

  ! c-maps for data
  c_glob_dat(:) = c_glob(:)

  ! UNIFORM PERTURBATION (01-Feb-2006)
  !afac = 5.
  !c_glob_dat(:) = c_glob_syn(:) * (1. + afac/100.)

  ! unperturbed structure
  if (PERT_STRUCT == 0) c_glob_syn(:) = c_glob_dat(:)

  ! write data phase velocity map to file
  file_dat_c = 'socal_vel_dat.dat'
  open(unit=18,file=trim(out_dir2)//file_dat_c,status='unknown')
  do iglob = 1,NGLOB
     write(18,*) sngl(x_lon(iglob)), sngl(z_lat(iglob)), sngl(c_glob_dat(iglob)), sngl(c_glob_dat(iglob)/c0 - 1.)
  enddo
  close(18)

  ! write syn phase velocity map to file
  ! (reference model is 0% perturbation)
  file_syn_c = 'socal_vel_syn.dat'
  open(unit=18,file=trim(out_dir2)//file_syn_c,status='unknown')
  do iglob = 1,NGLOB
     write(18,*) sngl(x_lon(iglob)), sngl(z_lat(iglob)), sngl(c_glob_syn(iglob)), sngl(c_glob_syn(iglob)/c0 - 1.)
  enddo
  close(18)

  ! write in c0 (and period) to file
  open(unit=16,file=trim(out_dir2)//'socal_vel_c0.dat', status='unknown')
  write(16,*) c0
  write(16,*) 2*hdur
  close(16)

  cmin = minval(c_glob)
  cmax = maxval(c_glob)

!!$  ! plot phase velocity map for synthetics (and data)
!!$  iopt = 1 + idat
!!$  filename1 = 'get_model.csh'
!!$  filename2 = trim(script_dir)//'plot_model.pl'
!!$  open(19,file=filename1,status='unknown')
!!$  write(19,'(4a,i5,3a,2f12.6,2a)') trim(filename2),' ', trim(out_dir2),' ', &
!!$     iopt, ' ', trim(file_syn_c), ' ', sngl(c0), sngl(t_target),' ',trim(file_dat_c)
!!$  close(19)
!!$  call system('chmod 755 get_model.csh ; get_model.csh')

  ! assign model parameters to gridpoints
  ihomo = 0
  call set_model_property()

  print *, 'Velocities in km/s :'
  print *, '  S     :', sngl(sqrt(RIGIDITY/DENSITY)/1000.)
  print *, '  cmin  :', sngl(cmin/1000.)
  print *, '  c0    :', sngl(c0/1000.)
  print *, '  cmax  :', sngl(cmax/1000.)
  print *
  print *, ' actual  rigidity (Pa-s)       : ', sngl(RIGIDITY)
  print *, ' desired rigidity (Pa-s), cmin : ', sngl(DENSITY * cmin**2)
  print *, ' desired rigidity (Pa-s), c0   : ', sngl(DENSITY * c0**2)
  print *, ' desired rigidity (Pa-s), cmax : ', sngl(DENSITY * cmax**2)

!--------------------------------------
! SoCal events

  ! read in target events (socal_quakes.m)
  ! socal events M > 4 (1990-2005), selected to get semi-uniform coverage
  !quake_file = 'socal_quakes_v02.dat'
  quake_file = 'socal_quakes_N025.dat'
  open(55,file=trim(in_dir)//trim(quake_file),status='unknown')
  read(55,*) nevent
  if (nevent > MAX_EVENT) stop 'nevent > MAX_EVENT (so increase MAX_EVENT)'
  print *, nevent, ' target events (for synthetics)'
  read(55,'(i14,3f14.7)') (socal_events_lab(i),socal_events_lon(i),socal_events_lat(i),socal_events_mag(i),i=1,nevent)
  close(55)

  print *
  print *, 'SoCal events (label, lon, lat, mag):'
  do i=1,nevent
     write(*,'(i14,3f14.7)') socal_events_lab(i), socal_events_lon(i), socal_events_lat(i), socal_events_mag(i)
  enddo

!--------------------------------------
! Prior to Feb 2006, we simply assigned the source location to the nearst gridpoint.
! But now we compute the locations such that ANY location is allowed.
!--------------------------------------
! events for data (eglob_dat)

  ! fill a portion of the [1:MAX_EVENT] vector with the SoCal events
  x_eve_lon0_dat(1:nevent) = socal_events_lon(:)
  z_eve_lat0_dat(1:nevent) = socal_events_lat(:)

  ! all vectors should be this length
  nevent = MAX_EVENT

  ! convert from lat-lon to mesh coordinates (in meters) -- update nevent
  call mesh_geo(nevent,x_eve_lon0_dat,z_eve_lat0_dat,x_eve0_dat,z_eve0_dat,UTM_PROJECTION_ZONE,ILONLAT2MESH)

  ! filter target points (utm-mesh) (update nevent)
  print *
  print *, 'events for data:'
  call station_filter(nevent, x_eve0_dat, z_eve0_dat, ifilter_eve_dat, SOURCE_GRID_BUFFER)
  !call station_filter(nevent, x_eve, z_eve, dmin_trsh_s, ncoast, coast_x, coast_z)

  if (nevent < 1) stop 'Must have at least one event'

  ! allocate variables
  allocate(x_eve_dat(nevent),z_eve_dat(nevent),x_eve_lon_dat(nevent),z_eve_lat_dat(nevent))
  allocate(eglob_dat(nevent))
  allocate(ispec_eve_dat(nevent),xi_eve_dat(nevent),gamma_eve_dat(nevent))
  allocate(otime_dat(nevent))

  ! perturbation from origin time for data
  otime_dat(:) = origin_time_dat

  ! assign allowable target events to (x_eve_dat, z_eve_dat)
  do i=1,nevent
    x_eve_dat(i) = x_eve0_dat(ifilter_eve_dat(i))
    z_eve_dat(i) = z_eve0_dat(ifilter_eve_dat(i))
  enddo

  !do i=1,nevent
  !   write(*,'(2e24.16)') x_eve_dat(i), z_eve_dat(i)
  !enddo

  ! determine the (eta, xi) corresponding to the target points
  ! this UPDATES x_eve, z_eve; eglob_dat is the index of the closest gridpoint
  call locate_targets(nevent, x_eve_dat, z_eve_dat, eglob_dat, ispec_eve_dat, xi_eve_dat, gamma_eve_dat)

  !do i=1,nevent
  !   write(*,'(2e24.16)') x_eve_dat(i), z_eve_dat(i)
  !enddo

  print *
  print *, 'ievent, x, z, iglob, ispec, xi, gamma'
  do i=1,nevent
     write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_eve_dat(i), z_eve_dat(i), &
        eglob_dat(i), ispec_eve_dat(i), xi_eve_dat(i), gamma_eve_dat(i)
  enddo

  ! convert from mesh to lat-lon
  call mesh_geo(nevent, x_eve_lon_dat, z_eve_lat_dat, x_eve_dat, z_eve_dat, UTM_PROJECTION_ZONE, IMESH2LONLAT)

  ! display target events and final events -- and the distance between (in meters)
  ! The distance error is due to the UTM conversion, but the lat-lon points are only
  ! used for plotting purposes, so this is fine.
  print *
  print *, 'events [x_eve_lon0_dat, z_eve_lat0_dat, x_eve_lon_dat, x_eve_lat_dat, dist (m)]:'
  do i=1,nevent
     temp1 = x_eve_lon0_dat(ifilter_eve_dat(i))
     temp2 = z_eve_lat0_dat(ifilter_eve_dat(i))
     temp3 = x_eve_lon_dat(i)
     temp4 = z_eve_lat_dat(i)
     temp5 = acos( sin(temp2/DEG)*sin(temp4/DEG) + cos(temp2/DEG)*cos(temp4/DEG)*cos(temp1/DEG - temp3/DEG) )
     write(*,'(i8, 5f17.10)') i, temp1, temp2, temp3, temp4, temp5*6371.*1000.
  enddo

  ! events : convert target points from utm to gridpoint index to lat-lon
  !call set_glob(nsrc, x_eve, z_eve, eglob)  ! get nearest gridpoint (eglob) to (x_eve, z_eve)
  !nevent = ns
  !x_eve_lon(:) = 0.  ; z_eve_lat(:) = 0.  ! re-initialize
  !do i=1,nevent
  !   xtemp = x(eglob(i)) + utm_xmin
  !   ztemp = z(eglob(i)) + utm_zmin
  !   call utm_geo(x_eve_lon(i),z_eve_lat(i),xtemp,ztemp,UTM_PROJECTION_ZONE,IUTM2LONGLAT)
  !enddo
  !!call mesh_geo(nevent, x_eve_lon, z_eve_lat, x_eve, z_eve, UTM_PROJECTION_ZONE, IMESH2LONLAT)

  ! write events to file
  open(19,file=trim(out_dir2)//'events_lonlat_dat.dat',status='unknown')
  do ievent = ievent_min, ievent_max
     !write(19,*) sngl(x_lon(eglob_dat(ievent))), sngl(z_lat(eglob_dat(ievent))), ievent
     write(19,'(2f20.10,1i12)') x_eve_lon_dat(ievent), z_eve_lat_dat(ievent), ievent
  enddo
  close(19)
  open(19,file=trim(out_dir2)//'events_xy_dat.dat',status='unknown')
  do ievent = ievent_min, ievent_max
     write(19,'(3f20.10,1i12)') x_eve_dat(ievent), z_eve_dat(ievent), otime_dat(ievent), ievent
  enddo
  close(19)

  !do ievent=1,nevent
  !   print *, x_eve_lon_dat(ievent),z_eve_lat_dat(ievent)
  !enddo

!stop 'testing'

!--------------------------------------
! events for synthetics (eglob)

  ! allocate variables
  allocate(x_eve(nevent),z_eve(nevent),x_eve_lon(nevent),z_eve_lat(nevent))
  allocate(eglob(nevent))
  allocate(ispec_eve(nevent),xi_eve(nevent),gamma_eve(nevent))
  allocate(otime(nevent))

  if (PERT_SOURCE == 1) then  ! source perturbations

    if (1 == 1) then   ! read in perturbed events from another file

      open(19,file='/home/store2/carltape/'//trim(out_dir3)//'run_2550/events_xy.dat',status='unknown')
      do ievent = 1,25
         read(19,'(3f20.10,1i12)') x_eve(ievent), z_eve(ievent), otime(ievent), itemp1
      enddo
      close(19)

    else

       ! perturb the events to get target events for the synthetics
       ! perturbation is described as (r,theta) polar coordinates
       amin = 0.
       amax = 2.*PI
       rmin = 0.
       rmax = src_pert_dist
       ptmin = -abs(src_pert_time)
       ptmax =  abs(src_pert_time)

       do ievent = 1,nevent

          ! origin time perturbation
          call random_number(rand1)
          ptime = (ptmax - ptmin)*rand1 + ptmin
          !ptime = src_pert_time

          ! azimuthal perturbation (in radians)
          call random_number(rand1)
          pangle = (amax - amin)*rand1 + amin
          !pangle = 30.*(PI/180.)               ! fix for TESTING

          ! radial perturbation (in meters)
          call random_number(rand1)
          prad = (rmax - rmin)*rand1 + rmin
          !prad = rmax

          ! fill vectors
          otime(ievent) = otime_dat(ievent) - ptime                ! NOTE THE SIGN
          x_eve(ievent) = x_eve_dat(ievent) + cos(pangle)*prad
          z_eve(ievent) = z_eve_dat(ievent) + sin(pangle)*prad
          temp1 = sqrt((x_eve_dat(ievent)-x_eve(ievent))**2 + (z_eve_dat(ievent)-z_eve(ievent))**2)

          !write(*,'(5f18.8)') x_eve(ievent)/1000., z_eve(ievent)/1000., &
          !     x_eve_dat(ievent)/1000., z_eve_dat(ievent)/1000., temp1/1000.
       enddo

    endif

    ! determine the (eta, xi) corresponding to the target points
    ! this UPDATES x_eve, z_eve; eglob is the index of the closest gridpoint
    call locate_targets(nevent, x_eve, z_eve, eglob, ispec_eve, xi_eve, gamma_eve)

    print *
    print *, 'ievent, x, z, iglob, ispec, xi, gamma'
    do i=1,nevent
      write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_eve(i), z_eve(i), &
         eglob(i), ispec_eve(i), xi_eve(i), gamma_eve(i)
    enddo

    ! convert from mesh to lat-lon
    call mesh_geo(nevent, x_eve_lon, z_eve_lat, x_eve, z_eve, UTM_PROJECTION_ZONE, IMESH2LONLAT)

    ! display target events and perturbed events -- and the distance between (in meters)
    print *
    print *, 'events [x_eve_lon_dat, z_eve_lat_dat, x_eve_lon, x_eve_lat, dist (m)]:'
    do i=1,nevent
       temp1 = x_eve_lon_dat(i)
       temp2 = z_eve_lat_dat(i)
       temp3 = x_eve_lon(i)
       temp4 = z_eve_lat(i)
       temp5 = acos( sin(temp2/DEG)*sin(temp4/DEG) + cos(temp2/DEG)*cos(temp4/DEG)*cos(temp1/DEG - temp3/DEG) )
       write(*,'(i8, 5f17.10)') i, temp1, temp2, temp3, temp4, temp5*6371.*1000.
    enddo

  else  ! no source perturbation

     eglob(:) = eglob_dat(:)
     x_eve(:) = x_eve_dat(:)
     z_eve(:) = z_eve_dat(:)
     x_eve_lon(:) = x_eve_lon_dat(:)
     z_eve_lat(:) = z_eve_lat_dat(:)
     ispec_eve(:) = ispec_eve_dat(:)
     xi_eve(:) = xi_eve_dat(:)
     gamma_eve(:) = gamma_eve_dat(:)

     otime(:) = otime_dat(:)   ! origin time

  endif   ! PERT_SOURCE

  ! write events for synthetics to file
  open(19,file=trim(out_dir2)//'events_lonlat.dat',status='unknown')
  do ievent = ievent_min, ievent_max
    write(19,'(2f20.10,1i12)') x_eve_lon(ievent), z_eve_lat(ievent), ievent
  enddo
  close(19)
  open(19,file=trim(out_dir2)//'events_xy.dat',status='unknown')
  do ievent = ievent_min, ievent_max
     write(19,'(3f20.10,1i12)') x_eve(ievent), z_eve(ievent), otime(ievent), ievent
  enddo
  close(19)

!stop 'testing'

!--------------------------------------
! receivers

  ! allocation and initialization
  !allocate(rglob_temp(MAX_SR))
  !allocate(x_rec0(MAX_SR),z_rec0(MAX_SR),x_rec_lon0(MAX_SR),z_rec_lat0(MAX_SR),ifilter_rec(MAX_SR))
  x_rec0(:) = 0.      ; z_rec0(:) = 0
  x_rec_lon0(:) = 0.  ; z_rec_lat0(:) = 0.
  !rglob_temp(:) = 0

! get the lat-lon of the TARGET RECEIVERS

if (IREC_SPACE == 1) then  ! individual receivers

  ! target receiver
  !x_rec0(1) = 3 * LENGTH/4     ; z_rec0(1) = HEIGHT/2
  x_rec_lon0(1) = -119.6        ; z_rec_lat0(1) = 34.0
  x_rec_lon0(2) = -116.0        ; z_rec_lat0(2) = 32.5
  x_rec_lon0(3) = -116.5        ; z_rec_lat0(3) = 35.5
  !x_rec_lon0(1) = -116.5        ; z_rec_lat0(1) = 35.5

  !x_rec_lon0(1) = -117.450717    ; z_rec_lat0(1) = 34.542814   ! 001
  !x_rec_lon0(1) = -118.141053    ; z_rec_lat0(1) = 35.547570   ! 126

  !x_rec_lon0(1) = socal_events_lon(ievent) + 1.0   ; z_rec_lat0(1) = socal_events_lat(ievent) + 1.0
  !x_rec_lon0(2) = socal_events_lon(ievent) + 1.0   ; z_rec_lat0(2) = socal_events_lat(ievent) - 1.0
  !x_rec_lon0(3) = socal_events_lon(ievent) - 1.0   ; z_rec_lat0(3) = socal_events_lat(ievent) + 1.0
  !x_rec_lon0(4) = socal_events_lon(ievent) - 1.0   ; z_rec_lat0(4) = socal_events_lat(ievent) - 1.0

  nrec = 3

else if (IREC_SPACE == 2) then  ! actual station locations

  ! read in (target) receivers
  recfile = trim(in_dir)//'STATION_149_full'
  !recfile = trim(in_dir)//'STATION_144_noisland'
  !recfile = trim(in_dir)//'STATION_067_nobasin'

  open(88,file=trim(recfile),status='unknown')
  read(88,*) nrec
  read(88,*) (x_rec_lon0(i),z_rec_lat0(i),i=1,nrec)
  close(88)

else if (IREC_SPACE == 4) then ! 'regular' mesh of receivers

   ! calculate mesh spacing
   dx = LENGTH/NMESH_REC
   dz = dx
   j = 0
   do xmesh = STATION_GRID_BUFFER+1.,LENGTH,dx
      do zmesh = STATION_GRID_BUFFER+1.,HEIGHT,dz
            j = j+1
            x_rec0(j) = xmesh
            z_rec0(j) = zmesh
      enddo
   enddo
   nrec = j

   ! get x_rec_lon0, z_rec_lat0
   call mesh_geo(MAX_SR,x_rec_lon0,z_rec_lat0,x_rec0,z_rec0,UTM_PROJECTION_ZONE,IMESH2LONLAT)

else
   stop 'invalid entry for IREC_SPACE'

endif  ! IREC_SPACE

  ! make sure that there are fewer target points than the max allowed
  if (nrec > MAX_SR) then
     print *
     print *, ' IREC_SPACE = ', IREC_SPACE
     print *, '       nrec = ', nrec
     print *, '     MAX_SR = ', MAX_SR
     stop 'rec > MAX_SR so increase MAX_SR'
  endif

  print *
  print *, 'receivers (adjoint sources):'

  ! all vectors should be length MAX_SR
  nrec = MAX_SR

  ! convert from lat-lon to mesh coordinates (in meters)
  call mesh_geo(nrec,x_rec_lon0,z_rec_lat0,x_rec0,z_rec0,UTM_PROJECTION_ZONE,ILONLAT2MESH)

  ! filter target points (utm-mesh) -- update nrec
  call station_filter(nrec, x_rec0, z_rec0, ifilter_rec, STATION_GRID_BUFFER)
  !call station_filter(nrec, x_rec, z_rec, dmin_trsh_r, ncoast, coast_x, coast_z)
  !call station_filter_2(nrec, x_rec, z_rec, -1)  ! -1 for left, +1 for right

  if (nrec < 1) stop 'Must have at least one receiver'

  ! allocate vectors
  allocate(x_rec(nrec),z_rec(nrec),x_rec_lon(nrec),z_rec_lat(nrec))
  allocate(rglob(nrec))
  allocate(ispec_rec(nrec),xi_rec(nrec),gamma_rec(nrec))

  ! assign allowable target receivers to (x_rec, z_rec)
  do i=1,nrec
    x_rec(i) = x_rec0(ifilter_rec(i))
    z_rec(i) = z_rec0(ifilter_rec(i))
  enddo

  ! determine the (eta, xi) corresponding to the target points
  ! this UPDATES x_rec, z_rec; rglob is the index of the closest gridpoint
  call locate_targets(nrec, x_rec, z_rec, rglob, ispec_rec, xi_rec, gamma_rec)

  print *
  print *, 'irec, x, z, iglob, ispec, xi, gamma'
  do i=1,nrec
     write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_rec(i), z_rec(i), rglob(i), ispec_rec(i), xi_rec(i), gamma_rec(i)
  enddo

  ! convert from mesh to lat-lon
  call mesh_geo(nrec, x_rec_lon, z_rec_lat, x_rec, z_rec, UTM_PROJECTION_ZONE, IMESH2LONLAT)

  ! display target receivers and final receivers -- and the distance between (in meters)
  ! The distance error is due to the UTM conversion, but the lat-lon points are only
  ! used for plotting purposes, so this is fine.
  print *
  print *, ' receivers [x_rec_lon0, z_rec_lat0, x_rec_lon, x_rec_lat, dist (m)] :'
  do i=1,nrec
     temp1 = x_rec_lon0(ifilter_rec(i))
     temp2 = z_rec_lat0(ifilter_rec(i))
     temp3 = x_rec_lon(i)
     temp4 = z_rec_lat(i)
     temp5 = acos( sin(temp2/DEG)*sin(temp4/DEG) + cos(temp2/DEG)*cos(temp4/DEG)*cos(temp1/DEG - temp3/DEG) )
     write(*,'(i8, 5f17.10)') i, temp1, temp2, temp3, temp4, temp5*6371.*1000.
  enddo

  ! receiver indices for writing the misfit chi(rec), summed over all events
  rglob0(1:nrec) = rglob(1:nrec)

  !do i=1,nrec
  !   print *, sngl(x(rglob(i))/1000.),sngl(z(rglob(i))/1000.),sngl(x_rec(i)/1000.),sngl(z_rec(i)/1000.)
  !enddo

  !stop 'testing'

  ! deallocate
  !deallocate(x_rec,z_rec)

  ! write receivers to file
  open(19,file=trim(out_dir2)//'recs_lonlat.dat',status='unknown')
  do irec = 1,nrec
     write(19,'(2f20.10,1i12)') x_rec_lon(irec), z_rec_lat(irec), irec
  enddo
  close(19)

!--------------------------------------
! allocate Lagrange interpolators for sources and receivers for both data and synthetics

  ! allocate 1-D Lagrange interpolators and derivatives
  allocate(hxied(NGLLX), hpxied(NGLLX), hgammaed(NGLLZ), hpgammaed(NGLLZ))
  allocate(hxied_store(nevent,NGLLX), hgammaed_store(nevent,NGLLZ))

  allocate(hxie(NGLLX), hpxie(NGLLX), hgammae(NGLLZ), hpgammae(NGLLZ))
  allocate(hxie_store(nevent,NGLLX), hgammae_store(nevent,NGLLZ))

  allocate(hxir(NGLLX), hpxir(NGLLX), hgammar(NGLLZ), hpgammar(NGLLZ))
  allocate(hxir_store(nrec,NGLLX), hgammar_store(nrec,NGLLZ))

  ! events for data
  do ievent = 1,nevent
    call lagrange_poly(xi_eve_dat(ievent),NGLLX,xigll,hxied,hpxied)
    call lagrange_poly(gamma_eve_dat(ievent),NGLLZ,zigll,hgammaed,hpgammaed)
    hxied_store(ievent,:)    = hxied(:)
    hgammaed_store(ievent,:) = hgammaed(:)
    !write(*,'(i8,5f16.10)') ievent, hxied(:)       ! NGLLX
    !write(*,'(i8,5f16.10)') ievent, hgammaed(:)    ! NGLLZ
  enddo
  deallocate(hxied, hpxied, hgammaed, hpgammaed)

  ! events for synthetics
  do ievent = 1,nevent
    call lagrange_poly(xi_eve(ievent),NGLLX,xigll,hxie,hpxie)
    call lagrange_poly(gamma_eve(ievent),NGLLZ,zigll,hgammae,hpgammae)
    hxie_store(ievent,:)    = hxie(:)
    hgammae_store(ievent,:) = hgammae(:)
    !write(*,'(i8,5f16.10)') ievent, hxie(:)       ! NGLLX
    !write(*,'(i8,5f16.10)') ievent, hgammae(:)    ! NGLLZ
  enddo
  deallocate(hxie, hpxie, hgammae, hpgammae)

  ! receivers
  do irec = 1,nrec
    call lagrange_poly(xi_rec(irec),NGLLX,xigll,hxir,hpxir)
    call lagrange_poly(gamma_rec(irec),NGLLZ,zigll,hgammar,hpgammar)
    hxir_store(irec,:)    = hxir(:)
    hgammar_store(irec,:) = hgammar(:)
    !write(*,'(i8,5f16.10)') irec, hxir(:)       ! NGLLX
    !write(*,'(i8,5f16.10)') irec, hgammar(:)    ! NGLLZ
  enddo
  deallocate(hxir, hpxir, hgammar, hpgammar)

!stop 'testing'

!--------------------------------------
! fake receivers for recording seismograms from the adjoint wavefield
! (note that we do not compute lat-lon for these points)
! THIS PORTION COULD BE REMOVED -- IT WAS FOR COMPUTING THE SPECTRAL MAP
! FOR THE OMS SIMULATIONS.

   ! mesh of target points
   x_recf0(:) = 0.
   z_recf0(:) = 0.
   dx = LENGTH/9.
   dz = dx
   i = 0
   do xmesh = 0.,LENGTH,dx
      do zmesh = 0.,HEIGHT,dz
         i = i+1
         if (i > MAX_SR_FAKE) stop 'i > MAX_SR_FAKE so change dx, dz, or MAX_SR_FAKE'
         x_recf0(i) = xmesh
         z_recf0(i) = zmesh
      enddo
   enddo
   !nrecf = i

  print *
  print *, 'fake receivers (for spectral maps):'

  ! all vectors should be length MAX_SR_FAKE
  nrecf = MAX_SR_FAKE

  ! filter target points (utm-mesh) -- update nrecf
  call station_filter(nrecf, x_recf0, z_recf0, ifilter_recf, STATION_GRID_BUFFER)

  if (nrecf < 1) stop 'Must have at least one fake (adjoint) receiver'

  ! allocate vectors
  allocate(x_recf(nrecf),z_recf(nrecf),fglob(nrecf))

  ! assign allowable target receivers to (x_rec, z_rec)
  do i=1,nrecf
    x_recf(i) = x_recf0(ifilter_recf(i))
    z_recf(i) = z_recf0(ifilter_recf(i))
    !print *, i, x_recf(i), z_recf(i)
  enddo

  ! fglob is the index of the CLOSEST GRIDPOINT
  call locate_targets(nrecf, x_recf, z_recf, fglob)

  !do i=1,nrecf
  !  print *, sngl(x(fglob(i))/1000.),sngl(z(fglob(i))/1000.),sngl(x_recf(i)/1000.),sngl(z_recf(i)/1000.)
  !enddo

!!$  ! display target receivers and final fake receivers -- distances in km
!!$  print *
!!$  print *, ' fake receivers [x_rec0, z_rec0, x_rec, x_rec, dist (km)]:'
!!$  do i=1,nrecf
!!$     temp1 = x_recf0(ifilter_recf(i)) / 1000.
!!$     temp2 = z_recf0(ifilter_recf(i)) / 1000.
!!$     temp3 = x_recf(i) / 1000.
!!$     temp4 = z_recf(i) / 1000.
!!$     temp5 = sqrt( (temp3-temp1)**2 + (temp4-temp2)**2 )
!!$     write(*,'(i8,5f17.10)') i, temp1, temp2, temp3, temp4, temp5
!!$  enddo

  ! deallocate: all we need is fglob
  deallocate(x_recf,z_recf)

!stop 'testing'

!--------------------------------------
! initial model vectors

  nmod_str = NGLOB              ! source parameters
  nmod_src = 3*nevent           ! structure parameters
  nmod = nmod_str + nmod_src    ! total model parameters

  ! allocate model vector
  allocate(m0(nmod),mt(nmod),gradient(nmod))
  allocate(m0_vec(nmod),mt_vec(nmod))
  allocate(g0(nmod),gt(nmod),gk(nmod),p0(nmod),pk(nmod))

  !allocate(m_vel(nmod_str),m0_vel(nmod_str),mt_vel(nmod_str))

  m0(:) = 0.
  mt(:) = 0.
  m0_vec(:) = 0.
  mt_vec(:) = 0.
  gradient(:) = 0.

  !------------------
  ! source-related vectors

  ! scaling for source model coefficients
  allocate(source_gradient(nmod_src),m_scale_src(nmod_src))
  allocate(m_src(nmod_src),mt_src(nmod_src),m0_src(nmod_src))
  allocate(m_src_syn(nmod_src),m_src_dat(nmod_src))

  source_gradient(:) = 0.
  m0_src(:) = 0.           ! obsolete
  mt_src(:) = 0.           ! obsolete
  m_src(:) = 0.
  m_src_dat(:) = 0.
  m_src_syn(:) = 0.

  ! fill source model vector for DATA
  do ievent = 1,nevent
    itemp1 = (ievent-1)*3 + 1
    itemp2 = (ievent-1)*3 + 2
    itemp3 = (ievent-1)*3 + 3

    m_src_dat(itemp1) = x_eve_dat(ievent)
    m_src_dat(itemp2) = z_eve_dat(ievent)
    m_src_dat(itemp3) = otime_dat(ievent)
  enddo

  ! fill source model vector for SYNTHETICS
  !origin_time = origin_time_dat - src_pert_time
  do ievent = 1,nevent
    itemp1 = (ievent-1)*3 + 1
    itemp2 = (ievent-1)*3 + 2
    itemp3 = (ievent-1)*3 + 3

    m_src_syn(itemp1) = x_eve(ievent)
    m_src_syn(itemp2) = z_eve(ievent)
    m_src_syn(itemp3) = otime(ievent)
  enddo

  !m_scale_src(:) = m_src_syn(:)

  m_scale(1) = 2*hdur*c0   ! dxs -- wavelength
  m_scale(2) = 2*hdur*c0   ! dzs -- wavelength
  m_scale(3) = 2*hdur      ! dt0 -- period
  do ievent = 1,nevent
    do i = 1,3
      itemp = (ievent-1)*3 + i
      m_scale_src(itemp) = m_scale(i)
    enddo
  enddo

  open(88,file=trim(out_dir2)//'scale_source.dat',status='unknown')
  do i=1,nmod_src
    write(88,'(1e24.12)') m_scale_src(i)
  enddo
  close(88)

  !stop 'testing'

  !------------------
  ! initial model vector (physical coordinates)

  ! structure portion of initial model vector
  !do iglob = 1,NGLOB
  !  !m0(iglob)     = c_glob_syn(iglob) / c0 - 1.  ! fractional perturbation
  !  m0_vec(iglob) = c_glob_syn(iglob)            ! m/s
  !enddo

  do i = 1,nmod
    if (i <= nmod_str) then
      m0_vec(i) = c_glob_syn(i)            ! m/s
    else
      m0_vec(i) = m_src_syn(i-nmod_str)    ! xs, zs, t0
    endif
  enddo

  open(88,file=trim(out_dir2)//'initial_model_vector.dat',status='unknown')
  do i=1,nmod
    write(88,'(1e24.12)') m0_vec(i)
  enddo
  close(88)
  !stop 'testing'

  !------------------

if (ISOURCE_LOG) open(91,file=trim(out_dir2)//'source_vector.log',status='unknown')

itest = 0
!============================================
! LOOP 2: over the models in the CG optimization
  do istep = 0,2*NITERATION
!============================================

  imod = (istep - mod(istep,2))/2          ! index into model number

!!$  if ( mod(imod,2) == 0 ) then
!!$  !if (imod <= 1 .or. imod==5 .or. imod==9 .or. imod==13) then
!!$    INV_SOURCE = 1
!!$    INV_STRUCT = 0
!!$  else
!!$    INV_SOURCE = 0
!!$    INV_STRUCT = 1
!!$  endif

  irun = irun0 + istep
  print *,' ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  == '
  print *,'  istep, imod, irun : ', istep, imod, irun
  if (INV_STRUCT == 1) print *, '  inverting for structure parameters'
  if (INV_SOURCE == 1) print *, '  inverting for source parameters'
  print *,' ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  == '

!enddo
!stop 'itesting'
!do istep = 0,2*NITERATION

  ! name and create the output directory
  write(out_dir1,'(a,i4.4,a)') trim(out_dir3)//'run_',irun,'/'
  command1 = 'mkdir ' // trim(out_dir1)
  call system(command1)
  !open(19,file='temp.csh',status='unknown')
  !write(19,'(2a)') 'mkdir ',out_dir1
  !close(19)
  !call system('chmod 755 temp.csh ; temp.csh')

  ! use the reference model or test model (source and structure)
  ! structure is the top portion; source is the bottom portion
  if (itest == 0) then

     c_glob_syn(:) = m0_vec(1:nmod_str)
     m_src(:)      = m0_vec(nmod_str+1:nmod)
     !m_vel(:) = m0_vel(:)
     !m_src(:) = m0_src(:)

  else if (itest == 1) then
     c_glob_syn(:) = mt_vec(1:nmod_str)
     m_src(:)      = mt_vec(nmod_str+1:nmod)

     !m_vel(:) = mt_vel(:)
     !m_src(:) = mt_src(:)
  endif

  ! initialize measurement vector
  nmeas = nevent * nrec
  allocate(measure_vec(nmeas))
  measure_vec(:) = 0.
  imeasure = 0
  print *, ' nmeas = nevent * nrec = ', nmeas

  ! initialize source perturbation vector
  source_gradient(:) = 0.

  ! initialize summed kernel (structure gradient)
  kernel_basis_sum(:) = 0.

  !stop 'testing'

! loop over events
!============================================
! LOOP 3: over the events
  ifirst = 1
!  do ievent = ievent_min, ievent_max
  do ievent = 5,5    ! 1,5
!============================================

  print *,'------------------------------------'
  print *,'  EVENT NUMBER ',ievent
  print *,'------------------------------------'

  ! NAME THE OUTPUT DIRECTORY
  write(ev_dir,'(a,i3.3,a)') 'event_',ievent,'/'
  out_dir = trim(out_dir1)//trim(ev_dir)
  command1 = 'mkdir ' // trim(out_dir)
  call system(command1)
  !open(19,file='temp.csh',status='unknown')
  !write(19,'(2a)') 'mkdir ',out_dir
  !close(19)
  !call system('chmod 755 temp.csh ; temp.csh')

!--------------------------------------
! source

  ! get the lat-lon of the TARGET RECEIVERS

  if (ISRC_SPACE <= 5) then  ! if you do NOT want a point source from the event list

     x_src_lon0(:) = 0.  ; z_src_lat0(:) = 0.
     x_src0(:) = 0.      ; z_src0(:) = 0.

     if (ISRC_SPACE == 1) then  ! point source(s)

        !x_src0(1) = LENGTH/4      ; z_src0(1) = HEIGHT/2
        !x_src_lon0(1) = -118.5370     ; z_src_lat0(1) = 34.2130   ! Northridge
        !x_src_lon0(1) = -117.918     ; z_src_lat0(1) = 32.3290   ! 2004/06/15
        !x_src_lon0(1) = -117.776     ; z_src_lat0(1) = 33.917    ! Yorba Linda

        !x_src_lon0(1) = -119.0     ; z_src_lat0(1) = 33.0
        !x_src_lon0(2) = -118.0     ; z_src_lat0(2) = 33.5
        !x_src_lon0(3) = -119.5     ; z_src_lat0(3) = 34.3
        !x_src_lon0(1) = -119.5     ; z_src_lat0(1) = 33.0          ! single kernel

        x_src_lon0(1) = x_eve_lon(ievent)
        z_src_lat0(1) = z_eve_lat(ievent)
        nsrc = 1

     else if (ISRC_SPACE == 2) then  ! finite source segment

        ! specify the target starting point of the fault, the azimuth, and the length
        x_src_lon_i = -119.    ;  z_src_lat_i = 33.   ;   flen   = 100.0d+03  ! short fault
        !x_src_lon_i = -118.0   ;  z_src_lat_i = 32.   ;   flen   = 500.0d+03  ! long fault
        src_az = -45.*(PI/180.)        ! azimuth (clockwise from north)
        finc   = 1.5*dh                ! distance betwen target gridpoints (m)

        ! determine the target endpoint of the fault
        call utm_geo(x_src_lon_i, z_src_lat_i, x_src_i, z_src_i, UTM_PROJECTION_ZONE,ILONGLAT2UTM)
        x_src_i = x_src_i - utm_xmin
        z_src_i = z_src_i - utm_zmin
        x_src_f = x_src_i + sin(src_az)*flen
        z_src_f = z_src_i + cos(src_az)*flen

        ! determine the xz target points of the fault (Cartesian polar coordinates)
        i = 0
        do xd = 0,flen,finc
           i = i+1
           x_src0(i) = x_src_i + sin(src_az)*(i-1)*finc
           z_src0(i) = z_src_i + cos(src_az)*(i-1)*finc
           print *, x_src0(i)/1000., z_src0(i)/1000.
        enddo
        nsrc = i

        ! get target fault points in lat-lon
        call mesh_geo(MAX_SR,x_src_lon0,z_src_lat0,x_src0,z_src0,UTM_PROJECTION_ZONE,IUTM2LONGLAT)

     else if (ISRC_SPACE == 3) then  ! California continental shelf (OMS)

        x_src_lon0(1:nshelf) = shelf_lon(1:nshelf)
        z_src_lat0(1:nshelf) = shelf_lat(1:nshelf)
        nsrc = nshelf

     else if (ISRC_SPACE == 4) then  ! California coastline (OMS)

        x_src_lon0(1:ncoast) = coast_lon(1:ncoast)
        z_src_lat0(1:ncoast) = coast_lat(1:ncoast)
        nsrc = ncoast

     else if (ISRC_SPACE == 5) then  ! finite circular region

        ! lat-lon of the center point
        xcen_lon = -119.0
        zcen_lat = 33.0
        x_src_lon0(1) = xcen_lon
        z_src_lat0(1) = zcen_lat
        nsrc = 1

        call mesh_geo(nsrc,x_src_lon0,z_src_lat0,x_src0,z_src0,UTM_PROJECTION_ZONE,ILONLAT2MESH)
        xcen = x_src0(1)
        zcen = z_src0(1)

        ! target points within a radius of the center point
        s_radius = 30.0d+03
        dx = dh/2.
        dz = dx
        i = 0
        do xmesh = xcen-s_radius, xcen+s_radius, dx
           do zmesh = zcen-s_radius, zcen+s_radius, dx
              d = sqrt((xmesh - xcen)**2+(zmesh - zcen)**2)
              if (d < s_radius) then
                 i = i+1
                 x_src0(i) = xmesh
                 z_src0(i) = zmesh
              endif
           enddo
        enddo
        nsrc = i

        ! get circle points in lat-lon
        call mesh_geo(MAX_SR,x_src_lon0,z_src_lat0,x_src0,z_src0,UTM_PROJECTION_ZONE,IUTM2LONGLAT)

     endif  ! ISRC_SPACE

     do i=1,nsrc
        print *, x_src_lon0(i), z_src_lat0(i)
     enddo

     ! make sure that there are fewer target points than the max allowed
     ! (this is not ideal, since you might have a file of points that extend far outside the grid)
     if (nsrc > MAX_SR) then
        print *
        print *, ' ISRC_SPACE = ', ISRC_SPACE
        print *, '       nsrc = ', nsrc
        print *, '     MAX_SR = ', MAX_SR
        stop 'nsrc > MAX_SR so increase MAX_SR'
     endif

     print *
     print *, 'source(s) for event ', ievent, ' :'

     ! all vectors should be length MAX_SR
     nsrc = MAX_SR

     ! convert from lat-lon to mesh coordinates (in meters) -- get (x_src0,z_src0)
     call mesh_geo(nsrc,x_src_lon0,z_src_lat0,x_src0,z_src0,UTM_PROJECTION_ZONE,ILONLAT2MESH)

     ! filter target points (utm-mesh) -- update nsrc
     call station_filter(nsrc, x_src0, z_src0, ifilter_src, SOURCE_GRID_BUFFER)

     if (nsrc < 1) stop 'Must have at least one source'

     ! allocate vectors
     allocate(x_src(nsrc),z_src(nsrc),x_src_lon(nsrc),z_src_lat(nsrc))
     allocate(sglob(nsrc))
     allocate(ispec_src(nsrc),xi_src(nsrc),gamma_src(nsrc))

     ! assign allowable target sources to (x_src, z_src)
     do i=1,nsrc
        x_src(i) = x_src0(ifilter_src(i))
        z_src(i) = z_src0(ifilter_src(i))
     enddo

     ! determine the (eta, xi) corresponding to the target points
     ! this UPDATES x_src, z_src; sglob is the index of the closest gridpoint
     call locate_targets(nsrc, x_src, z_src, sglob, ispec_src, xi_src, gamma_src)

     print *
     do i=1,nsrc
        write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_src(i), z_src(i), sglob(i), ispec_src(i), xi_src(i), gamma_src(i)
     enddo

     ! convert from mesh to lat-lon
     call mesh_geo(nsrc, x_src_lon, z_src_lat, x_src, z_src, UTM_PROJECTION_ZONE, IMESH2LONLAT)

     ! display target sources and final sources -- and the distance between (in meters)
     ! The distance error is due to the UTM conversion, but the lat-lon points are only
     ! used for plotting purposes, so this is fine.
     print *
     print *, 'sources [x_src_lon0, z_src_lat0, x_src_lon, x_src_lat, dist (m)]:'
     do i=1,nsrc
        temp1 = x_src_lon0(ifilter_src(i))
        temp2 = z_src_lat0(ifilter_src(i))
        temp3 = x_src_lon(i)
        temp4 = z_src_lat(i)
        temp5 = acos( sin(temp2/DEG)*sin(temp4/DEG) + cos(temp2/DEG)*cos(temp4/DEG)*cos(temp1/DEG - temp3/DEG) )
        write(*,'(i8, 5f17.10)') i, temp1, temp2, temp3, temp4, temp5*6371.*1000.
     enddo

  else  ! select POINT SOURCE for array of events

     nsrc = 1

     !-----------------
     ! sources for data (no source perturbations)

     ! allocate vectors
     allocate(x_src_dat(nsrc),z_src_dat(nsrc),x_src_lon_dat(nsrc),z_src_lat_dat(nsrc))
     allocate(sglob_dat(nsrc))
     allocate(ispec_src_dat(nsrc),xi_src_dat(nsrc),gamma_src_dat(nsrc))

     ! allocate 1-D Lagrange interpolators and derivatives
     allocate(hxisd_store(nsrc,NGLLX),hgammasd_store(nsrc,NGLLZ))

     sglob_dat(1)     = eglob_dat(ievent)   ! closest gridpoint
     x_src_dat(1)     = x_eve_dat(ievent)
     z_src_dat(1)     = z_eve_dat(ievent)
     x_src_lon_dat(1) = x_eve_lon_dat(ievent)
     z_src_lat_dat(1) = z_eve_lat_dat(ievent)
     ispec_src_dat(1) = ispec_eve_dat(ievent)
     xi_src_dat(1)    = xi_eve_dat(ievent)
     gamma_src_dat(1) = gamma_eve_dat(ievent)
     hxisd_store(1,:) = hxied_store(ievent,:)
     hgammasd_store(1,:) = hgammaed_store(ievent,:)

     !-----------------
     ! sources for synthetics (allows for source perturbations)

     ! source perturbations for this event
     itemp1 = (ievent-1)*3 + 1
     itemp2 = (ievent-1)*3 + 2
     itemp3 = (ievent-1)*3 + 3

     ! allocate vectors
     allocate(x_src(nsrc),z_src(nsrc),x_src_lon(nsrc),z_src_lat(nsrc))
     allocate(sglob(nsrc))
     allocate(ispec_src(nsrc),xi_src(nsrc),gamma_src(nsrc))

     x_src(1) = m_src(itemp1)        ! new x position
     z_src(1) = m_src(itemp2)        ! new z position
     origin_time = m_src(itemp3)     ! new origin time

!!$        ! perturb source location from the previous model
!!$        ! this only changes the source if INV_SRC = 1
!!$        if (istep==0) then  ! initial source
!!$
!!$           x_src(1) = x_eve(ievent)          ! x position, perturbed event
!!$           z_src(1) = z_eve(ievent)          ! z position, perturbed event
!!$           origin_time = origin_time_dat - src_pert_time
!!$
!!$           ! initial source model
!!$           m_src(itemp1) = x_src(1)          ! x position
!!$           m_src(itemp2) = z_src(1)          ! z position
!!$           m_src(itemp3) = origin_time       ! origin time
!!$
!!$           !m_src(:) = m0_src(:)
!!$           m0_vec(nmod_str+1:nmod) = m_src(:)
!!$
!!$        else
!!$
!!$           !origin_time = origin_time_dat     ! testing
!!$
!!$        endif

!!$     ! write model vector to file
!!$     open(unit=19,file=trim(out_dir1)//'test.dat',status='unknown')
!!$     do i = 1,nmod
!!$       write(19,*) m0_vec(i)
!!$     enddo
!!$     close(19)
!!$     stop 'testing'

     ! filter target points (utm-mesh) -- update nsrc
     call station_filter(nsrc, x_src, z_src, ifilter_src, SOURCE_GRID_BUFFER)

     if (nsrc /= 1) stop 'Must be a single point source'

     ! determine the (eta, xi) corresponding to the target points
     ! this UPDATES x_src, z_src; sglob is the index of the closest gridpoint
     call locate_targets(nsrc, x_src, z_src, sglob, ispec_src, xi_src, gamma_src)

     print *
     do i=1,nsrc
        write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_src(i), z_src(i), sglob(i), ispec_src(i), xi_src(i), gamma_src(i)
        write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_src_dat(i), z_src_dat(i), &
           sglob_dat(i), ispec_src_dat(i), xi_src_dat(i), gamma_src_dat(i)
     enddo

     ! convert from mesh to lat-lon
     call mesh_geo(nsrc, x_src_lon, z_src_lat, x_src, z_src, UTM_PROJECTION_ZONE, IMESH2LONLAT)

     allocate(hxis(NGLLX), hpxis(NGLLX), hgammas(NGLLZ), hpgammas(NGLLZ))
     allocate(hxis_store(nsrc,NGLLX), hgammas_store(nsrc,NGLLZ))

     ! key variables for locating synthetic sources
     do isrc = 1,nsrc
       call lagrange_poly(xi_src(isrc),NGLLX,xigll,hxis,hpxis)
       call lagrange_poly(gamma_src(isrc),NGLLZ,zigll,hgammas,hpgammas)
       hxis_store(isrc,:)    = hxis(:)
       hgammas_store(isrc,:) = hgammas(:)
       !write(*,'(i8,5f16.10)') isrc, hxis(:)       ! NGLLX
       !write(*,'(i8,5f16.10)') isrc, hgammas(:)    ! NGLLZ
     enddo
     deallocate(hxis, hpxis, hgammas, hpgammas)

     !-------------------
     ! source for data and source for synthetic
     print *
     print *, 'sources [x_src_lon_dat, z_src_lat_dat, x_src_lon, z_src_lat, dist (m)]:'
     do i=1,nsrc
        temp1 = x_src_lon_dat(i)
        temp2 = z_src_lat_dat(i)
        temp3 = x_src_lon(i)
        temp4 = z_src_lat(i)
        temp5 = acos( sin(temp2/DEG)*sin(temp4/DEG) + cos(temp2/DEG)*cos(temp4/DEG)*cos(temp1/DEG - temp3/DEG) )
        write(*,'(i8, 5f17.10)') i, temp1, temp2, temp3, temp4, temp5*6371.*1000.
     enddo

  endif  ! ISRC_SPACE

  if (ISOURCE_LOG) then
    ! source log file
    itemp1 = (ievent-1)*3 + 1
    itemp2 = (ievent-1)*3 + 2
    itemp3 = (ievent-1)*3 + 3
    write(91,*)
    write(91,*) '------------------------'
    write(91,*) 'istep, imod, ievent, irun : '
    write(91,*) istep, imod, ievent, irun
    write(91,*) 'SOURCE MODEL (xs, ys, t0) :'
    write(91,*) ' [m_src]:'
    write(91,'(a12,3f18.8)') ' xs (km) : ', m_src(itemp1)/1000., &
       x_eve_dat(ievent)/1000., (m_src(itemp1) - x_eve_dat(ievent))/1000.
    write(91,'(a12,3f18.8)') ' zs (km) : ', m_src(itemp2)/1000., &
       z_eve_dat(ievent)/1000., (m_src(itemp2) - z_eve_dat(ievent))/1000.
    write(91,'(a12,3f18.8)') ' t0 (s) : ', m_src(itemp3), &
       origin_time_dat, m_src(itemp3) - origin_time_dat
  endif

!--------------------------------------
! source time function FOR DATA AND SYNTHETICS

  ! source magnitude (same for data and synthetics)
  if (NCOMP == 3) then
     f0(1) = 0.0    ;  f0(2) = FNORM    ; f0(3) = 0.0
  else if (NCOMP == 1) then
     f0(1) = FNORM
  else
     stop 'NCOMP must be 1 or 3'
  endif

  ! source time function for DATA
  stf_dat(:) = 0.
  call get_source_time_function(origin_time_dat,stf_dat,ti)
  !call get_source_time_function(nsrc,origin_time_dat,f0,samp_dat,ti)

  ! source function for data
  allocate(samp_dat(NSTEP,NCOMP,nsrc))
  samp_dat(:,:,:) = 0.0
  do i = 1,nsrc
    do icomp = 1,NCOMP
       samp_dat(:, icomp, i) = stf_dat(:) * f0(icomp)
    enddo
  enddo

  ! source time function for SYNTHETICS (allows for origin time perturbation)
  stf(:) = 0.
  call get_source_time_function(origin_time,stf,ti)
  !call get_source_time_function(nsrc,origin_time,f0,samp,ti)

  ! source function for synthetics
  allocate(samp(NSTEP,NCOMP,nsrc))
  samp(:,:,:) = 0.0
  do i = 1,nsrc
    do icomp = 1,NCOMP
       samp(:, icomp, i) = stf(:) * f0(icomp)
    enddo
  enddo

  ! write out source time functions for the first point source (for plot_model.pl)
  ! other source time functions will be written in write_seismogram.f90
  do i = 1,nsrc
     do icomp = 1,NCOMP
        write(file_stf,'(a,i5.5,a,i1)') trim(out_dir)//'stf_',i,'_',icomp
        open(12,file=file_stf,status='unknown')
        write(*,*)
        write(*,*) 'Event #', ievent, ', Source #', i
        write(*,*) '  actual location       : ', sngl(x_src_lon(i)), ', ', sngl(z_src_lat(i))
        write(*,*) '  closest GLL gridpoint : ', sngl(x_lon(sglob(i))), ', ', sngl(z_lat(sglob(i)))
        do itime=1,NSTEP
           write(12,'(1f16.6,1e16.6)') ti(itime), samp_dat(itime,icomp,i)
        enddo
        close(12)
     enddo
  enddo

!stop 'testing'

!!$  deallocate(samp, sglob, sglob_dat)
!!$  deallocate(x_src,z_src,x_src_lon,z_src_lat)
!!$  deallocate(ispec_src, xi_src, gamma_src)
!!$  enddo
!!$  stop 'testing'
!!$  do ievent = 1,nevent



!!$  do itime = 1, NSTEP
!!$    ti(itime) = dble(itime-1)*DT
!!$    stf = source_time_function(ti(itime)-tshift, hdur, istf)  ! note time shift (constants)
!!$    do i = 1, nsrc
!!$      samp(itime, :, i) = stf * f0(:)
!!$    enddo
!!$  enddo

!--------------------------------------
! testing filter routines

!!$print *, 'testing the fft subroutine using the source time function'
!!$call write_spectral_map(samp, nsrc, sglob, trim(out_dir)//'test')
!!$
!!$print *, 'testing the bandpass filter subroutine using the source time function'
!!$call filter(ti, samp, nsrc)
!!$
!!$print *, 'testing the fft subroutine using the bandpass-filtered source time function'
!!$call write_spectral_map(samp, nsrc, sglob, trim(out_dir)//'test_f')
!!$
!!$stop

!--------------------------------------

  ! write source and receivers to file
  ! NOTE THAT THE RECEIVER LABELING IS SIMPLY NUMERICAL ORDER (FOR NOW)
  !write(filename,'(a,i3.3,a)') trim(out_dir)//'sr_e',ievent,'.txt'
  file_src_rec = trim(out_dir)//'sr.txt'
  open(12,file=file_src_rec,status='unknown',iostat=ios)
  if (ios /= 0) stop 'Error opening out_dir/sr.txt'
  if (ISRC_SPACE /= 5) then
     write(12,'(a,2f12.6,i10)') ('S ', x_lon(sglob(i)), z_lat(sglob(i)), i, i=1,nsrc)
  else
     ! finite area source
     do i=1,nsrc
        d = sqrt((x(sglob(i)) - xcen)**2+(z(sglob(i)) - zcen)**2)
        if ( d > s_radius-dh) then  ! get outermost sources
           write(12,'(a,2f12.6,i10)') 'S ', x_lon(sglob(i)), z_lat(sglob(i)), i
        endif
     enddo
  endif
  write(12,'(a,2f12.6,i10)') ('R ', x_lon(rglob(i)), z_lat(rglob(i)), i, i=1,nrec)
  close(12)

  ! plot phase velocity map with source-receiver geometry and source time function
  iopt = 3 + idat
  if (ISURFACE == 1) then
     !filename1 = 'get_model.csh'
     !filename2 = trim(script_dir)//'plot_model.pl'
     !open(19,file=filename1,status='unknown')
     !write(19,'(4a,i5,3a,2f12.6,7a,i5)') trim(filename2),' ', trim(out_dir1),' ', &
     !    iopt, ' ', trim(file_syn_c), ' ', sngl(c0), sngl(t_target), &
     !    ' ',trim(file_dat_c),' ',trim(file_stf),' ', trim(file_src_rec),' ',FINITE_SOURCE
     !close(19)
     !call system('chmod 755 get_model.csh ; get_model.csh')
  endif

  ! the following are not needed, since we have sglob/rglob/fglob, the index vectors
  ! into the sources, receivers, and fake receivers

  deallocate(x_src,z_src,x_src_lon,z_src_lat)
  deallocate(x_src_dat,z_src_dat,x_src_lon_dat,z_src_lat_dat)

  !stop 'testing'

  !-----------------------------------------------------
  ! write the current models to file

  if (ifirst == 1) then

     ! write phase velocity maps to file (data is always the same)
     open(unit=18,file=trim(out_dir1)//file_dat_c,status='unknown')
     do iglob = 1,NGLOB
        write(18,*) sngl(x_lon(iglob)), sngl(z_lat(iglob)), sngl(c_glob_dat(iglob)), sngl(c_glob_dat(iglob)/c0 - 1.)
     enddo
     close(18)
     open(unit=18,file=trim(out_dir1)//file_syn_c,status='unknown')
     do iglob = 1,NGLOB
        write(18,*) sngl(x_lon(iglob)), sngl(z_lat(iglob)), sngl(c_glob_syn(iglob)), sngl(c_glob_syn(iglob)/c0 - 1.)
     enddo
     close(18)

     ! write source vectors to file (data is always the same)
     open(unit=19,file=trim(out_dir1)//'socal_src_dat.dat',status='unknown')
     do i = 1,nevent
        write(19,'(i8,3f20.12)') i, m_src_dat((i-1)*3 + 1)/1000., &
                                    m_src_dat((i-1)*3 + 2)/1000., &
                                    m_src_dat((i-1)*3 + 3)
     enddo
     close(19)
     open(unit=19,file=trim(out_dir1)//'socal_src_syn.dat',status='unknown')
     do i = 1,nevent
        write(19,'(i8,3f20.12)') i, m_src((i-1)*3 + 1)/1000., &
                                    m_src((i-1)*3 + 2)/1000., &
                                    m_src((i-1)*3 + 3)
     enddo
     close(19)

     ifirst = 0

  endif  ! ifirst

  !stop 'testing'

!=============================================================================================
! ******************************** WAVE SIMLUTAIONS ******************************************
!=============================================================================================

!=========================
! DATA (forward wavefield)

  if (WRITE_STF_F) call write_seismogram(samp_dat, nsrc, trim(out_dir)//'stffor_dat')

  ! compute data for misfit kernels
  allocate(data(NSTEP,NCOMP,nrec))
  data(:,:,:) = 0.0
  if (IKER <= 4) then

     ! set velocity field for the data
     c_glob(:) = c_glob_dat
     call set_model_property()

     ! DATA (forward wavefield)
     ! nsrc, sglob, samp_dat  : #src, src index, source
     ! nrec, rglob, data :      #rec, rec index, seis
     ! NOTE THAT sglob_dat ALLOWS FOR A PERTURBED SOURCE FROM THE SOURCES USED TO COMPUTE THE SYNTHETICS
     isolver = 1
     call solver(isolver, nsrc, sglob_dat, ispec_src_dat, hxisd_store, hgammasd_store, samp_dat, &
                          nrec, rglob,     ispec_rec,      hxir_store,  hgammar_store, data)

     ! write out seismograms at the receivers
     data_tag = 'dat'
     if (WRITE_SEISMO_F) call write_seismogram(data, nrec, trim(out_dir)//data_tag)
  endif

  !stop 'testing'

!=========================
! SYNTHETICS (forward wavefield)

  if (WRITE_STF_F) call write_seismogram(samp, nsrc, trim(out_dir)//'stffor_syn')

  !stop 'testing'

  ! forward wavefield
  allocate(syn(NSTEP,NCOMP,nrec))
  syn(:,:,:) = 0.0

  ! set velocity field for the synthetics (e.g., homogeneous reference model)
  c_glob(:) = c_glob_syn(:)
  call set_model_property()

  isolver = 1
  last_frame_name = trim(out_dir)//'last_frame.txt'
  allocate(absorb_field(NSTEP, NCOMP, NGLL, NELE, NABSORB))
  call solver(isolver, nsrc, sglob, ispec_src, hxis_store, hgammas_store, samp, &
                       nrec, rglob, ispec_rec, hxir_store, hgammar_store,  syn, trim(last_frame_name), absorb_field)

  ! write out seismograms at the receivers
  syn_tag = 'syn'
  !syn_tag = 'forward'
  if (WRITE_SEISMO_F) call write_seismogram(syn, nrec, trim(out_dir)//syn_tag)

  if (WRITE_SPECTRAL_MAP_F) then
     print *, 'compute and write out forward spectral map '
     call write_spectral_map(syn, nrec, rglob, trim(out_dir)//'spectral_forward',WRITE_SPECTRA_F)
  endif

  ! adjoint source function
  allocate(adj_syn(NSTEP,NCOMP,nrec))
  adj_syn(:,:,:) = 0.

  ! initialize time window for adjoint source to be the entire record
  allocate(tstart(nrec),tend(nrec))
  tstart(:) = 0.
  tend(:)   = NSTEP*DT

  ! compute time windows for adjoint source function
  ! d: distance from the FIRST SOURCE to the receivers (assuming homogeneous model)
  ! window with should be wide enough to capture synthetic and data waveforms,
  ! but narrow enough to exclude suprious boundary reflections
  print *
  print *, 'cut times for adjoint sources'
  do i=1,nrec
     d = sqrt( (x(sglob(1)) - x(rglob(i)))**2 + (z(sglob(1)) - z(rglob(i)))**2  )
     tcen      = tshift + d/c0

     ! wider window for data and synthetics FAR APART
     tstart(i) = tcen - 3*hdur
     tend(i)   = tcen + 3*hdur
     !tstart(i) = tcen - 2*hdur
     !tend(i)   = tcen + 2*hdur

     write(*,'(i8,3f12.4)') i, tstart(i), tcen, tend(i)
  enddo
  print *

  ! construct the adjoint source (IKER in wave2d_constants.f90)
  !call make_adjoint_source(nrec, syn, tstart, tend, adj_syn, data)
  call mtm_adj(IKER, ievent, nrec, syn, tstart, tend, adj_syn, data)

  ! write out adjoint source time function at the receivers
  stfadj_tag = 'stfadj'
  if (WRITE_STF_A) call write_seismogram(adj_syn, nrec, trim(out_dir)//stfadj_tag)

  !stop 'testing'

  ! OUTPUT ASCII FILES --> SAC FILES (make_sac_files.pl)
  ! (1) data, (2) synthetics, (3) adjoint source time function

!!$  if (WRITE_SEISMO_F) then
!!$     filename1 = trim(script_dir)//'make_sac_files.csh'
!!$     filename2 = 'make_sac_files.pl'
!!$     open(19,file=filename1,status='unknown')
!!$         if (IKER <= 4) write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', &
!!$             trim(out_dir),' ', trim(data_tag)  ,' ','1', tshift
!!$                       write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', &
!!$             trim(out_dir),' ', trim(syn_tag)   ,' ','1', tshift
!!$     if (WRITE_STF_A)   write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', &
!!$             trim(out_dir),' ', trim(stfadj_tag),' ','1', tshift
!!$     close(19)
!!$     call system('chmod 755 scripts/make_sac_files.csh ; scripts/make_sac_files.csh')
!!$  endif

  ! misfit for each receiver for one event
  chi_etot(ievent) = sum(chi(ievent,:,1,1))
  print *,'---------------------------'
  print *,'misfit for IKER = ',IKER
  print *,'chi-tot = ', sum(chi(ievent,:,1,1))
  print *,'---------------------------'
  open(19,file=trim(out_dir)//'chi_r.dat',status='unknown')
  do irec=1,nrec
     write(19,'(3e16.6)') x_lon(rglob(irec)), z_lat(rglob(irec)), chi(ievent,irec,1,1)
  enddo
  close(19)

  ! total misfit for one event
  open(19,file=trim(out_dir)//'single_chi_r.dat',status='unknown')
  write(19,*) sngl( sum(chi(ievent,:,1,1)) )
  close(19)

!!$  ! write ALL adjoint sources to file
!!$  do ipick = 0,6
!!$     call mtm_adj(ipick, ievent, nrec, syn, tstart, tend, adj_syn, data)
!!$     write(stfadj_tag,'(a,i1)') 'stfadj-', ipick
!!$     call write_source_function(nrec, ti, adj_syn, rglob, trim(out_dir)//stfadj_tag)
!!$  enddo
!!$
!!$  ! plot ALL adjoint sources
!!$  filename1 = trim(script_dir)//'make_sac_files.csh'
!!$  filename2 = 'make_sac_files.pl'
!!$  open(19,file=filename1,status='unknown')
!!$  write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', trim(out_dir),' ', trim(data_tag)  ,' ',"1",tshift
!!$  write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', trim(out_dir),' ', trim(syn_tag)   ,' ',"1",tshift
!!$  do ipick = 0,6
!!$     write(stfadj_tag,'(a,i1)') 'stfadj-', ipick
!!$     write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', trim(out_dir),' ', trim(stfadj_tag),' ',"1",tshift
!!$  enddo
!!$  close(19)
!!$  call system('chmod 755 scripts/make_sac_files.csh ; scripts/make_sac_files.csh')

!stop 'testing adjoint sources'

  !------------
  ! we always evaluate the kernels for the present model
  ! we only evaluate the kernel for the test model if itest==1 and POLY_ORDER==3

  if (itest == 0 .or. POLY_ORDER == 3) then
    print *
    print *, 'compute the kernel via adjoint wavefield interaction'
    print *

    isolver = 3
    allocate(rho_kernel(NGLOB), mu_kernel(NGLOB), kappa_kernel(NGLOB))

    ! initialize source perturbation time series
    allocate(three_source_model(NSTEP,NCOMP,nsrc,3))
    three_source_model(:,:,:,:) = 0.

    ! kernels
    ! samp enters as the forward source time function,
    ! and returns as the (time-reversed) adjoint seismograms
    call solver(isolver, nsrc, sglob, ispec_src, hxis_store, hgammas_store, samp, &
                         nrec, rglob, ispec_rec, hxir_store, hgammar_store, adj_syn, &
                 trim(last_frame_name), absorb_field, &
                 rho_kernel, mu_kernel, kappa_kernel, three_source_model, stf, f0)

      !call write_seismogram(three_source_model(:,:,:,1), nsrc, trim(out_dir)//'pert_src_m01')
      !call write_seismogram(three_source_model(:,:,:,2), nsrc, trim(out_dir)//'pert_src_m02')
      !call write_seismogram(three_source_model(:,:,:,3), nsrc, trim(out_dir)//'pert_src_m03')

      !call write_seismogram(three_source_model(:,:,:,4), nsrc, trim(out_dir)//'pert_src_m04')
      !call write_seismogram(three_source_model(:,:,:,5), nsrc, trim(out_dir)//'pert_src_m05')
      !call write_seismogram(three_source_model(:,:,:,6), nsrc, trim(out_dir)//'pert_src_m06')
      !call write_seismogram(three_source_model(:,:,:,7), nsrc, trim(out_dir)//'pert_src_m07')
      !call write_seismogram(three_source_model(:,:,:,8), nsrc, trim(out_dir)//'pert_src_m08')

      print *
      print *, ' to make the gradient for the source (unscaled)'
      print *, DT*sum(three_source_model(:,1,1,1))
      print *, DT*sum(three_source_model(:,1,1,2))
      print *, DT*sum(three_source_model(:,1,1,3))

      ! construct gradient vector
      ! ievent controls where the values are placed
      ! g(m) is the integrated time series times the source parameter
      do i = 1,3    ! i = 1,2 to test mislocations only
        itemp = (ievent-1)*3 + i
        !source_gradient(itemp) = DT*sum(three_source_model(:,1,1,i)) * m_scale(i)
        source_gradient(itemp) = DT*sum(three_source_model(:,1,1,i)) * m_scale_src(itemp)
      enddo

      ! summed kernel (global variable)
      kernel_basis_sum(:) = kernel_basis_sum(:) + kernel_basis(:)

    deallocate(three_source_model)
    deallocate(rho_kernel, mu_kernel, kappa_kernel)

    ! write out adjoint seismograms at the original sources
    if (WRITE_SEISMO_A) call write_seismogram(samp, nsrc, trim(out_dir)//'synadj')

  endif

!stop 'testing adjoint sources'

  ! deallocate variables associated with one event

  deallocate(data, adj_syn, absorb_field)
  deallocate(syn, tstart, tend)

  deallocate(ispec_src_dat, xi_src_dat, gamma_src_dat)
  deallocate(hxisd_store, hgammasd_store)
  deallocate(samp_dat, sglob_dat)

  deallocate(ispec_src, xi_src, gamma_src)
  deallocate(hxis_store, hgammas_store)
  deallocate(samp, sglob)

  !=====================
  enddo  ! ievent
  !=====================

  print *
  print *, ' Done with all the events for irun =', irun
  print *, ' Now we update the model, or compute a test model'

  !  stop 'testing'

  !-----------------------------------------------------
  ! misfit values

  ! summed chi value for each event (nevent by 1)
  open(19,file=trim(out_dir1)//'summed_chi_e.dat',status='unknown')
  do ievent=1,nevent
     write(19,*) sngl( sum(chi(ievent,:,1,1)) )
     !write(19,*) chi_etot(ievent)
  enddo
  close(19)

  ! summed chi value for each receiver (nrec by 1)
  open(19,file=trim(out_dir1)//'summed_chi_r.dat',status='unknown')
  do irec=1,nrec
     write(19,'(3e16.6)') x_lon(rglob0(irec)), z_lat(rglob0(irec)), sngl( sum(chi(:,irec,1,1)) )
  enddo
  close(19)

  ! overall misfit for all events and all receivers
  chi_val = sum(chi(:,:,1,1))
  open(19,file=trim(out_dir1)//'summed_chi_all.dat',status='unknown')
  write(19,*) chi_val
  close(19)

  ! write all nevent * nrec measurements to file
  open(unit=18,file=trim(out_dir1)//'measure_vec.dat',status='unknown')
  do i = 1,imeasure
     write(18,'(1e20.10)') measure_vec(i)
  enddo
  close(18)
  deallocate(measure_vec)

  if (ISOURCE_LOG) write(91,'(a12,1f18.8)') ' chi(m) : ', chi_val

  ! stopping criterion needs to be better defined (depends on misfit function --- and DT)
  ! (We needed this for the basic source inversions.)
  !if (chi_val <= 0.1 .and. itest==0) stop 'DONE: you have minimized chi(m) to chi(m) <= 0.1'

  print *, ' Written chi values to file'
  print *, ' Now we compute the gradient of the misfit function (using the misfit kernel)'

  !-----------------------------------------------------
  ! COMPUTE THE GRADIENT OF THE MISFIT function, if the present model is not
  ! a test model or if the CG polynomial is a cubic function
  ! DO NOT smooth kernel for test model if quadratic polynomial is being used

  if (itest == 0 .or. POLY_ORDER == 3) then

     print *, ' Computing the gradient of the misfit function for a given model'
     gradient(:) = 0.

     if (INV_STRUCT == 1) then   ! smooth the kernels to remove spurious src-rec effects

        ! summed kernel for all events (NGLOB by 1)
        open(19,file=trim(out_dir1)//'summed_ker.dat',status='unknown')
        do iglob=1,NGLOB
           write(19,'(3e16.6)') x_lon(iglob), z_lat(iglob), sngl(kernel_basis_sum(iglob))
        enddo
        close(19)

        !====================

        k_rough_global(:) = kernel_basis_sum(:)

        ! Gaussian half-width controlling the smoothing (m)
        !dtrsh2 = (3.*sigma)**2  ! all points outside d^2 are set to zero
        dtrsh2 = (9./8.)*gamma**2

        ! EXAMPLE Gaussian smoothing function for one point
        ! find the closest gridpoint to the target point
        xtar = 0.5*LENGTH
        ztar = 0.5*HEIGHT
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
        k_gaus_global_ex(:) = 0.
        do iglob = 1,NGLOB
           dist2 = (xcen - x(iglob))**2 + (zcen - z(iglob))**2
           if (dist2 <= dtrsh2) &
                                !k_gaus_global_ex(iglob) = (1./(2*PI*sigma**2)) * exp(-dist2 / (2.*sigma**2))
                k_gaus_global_ex(iglob) = (4./(PI*gamma**2)) * exp(-4.*dist2 / (gamma**2))
        enddo

        ! compute the SMOOTHED kernel
        k_smooth_global(:) = 0.
        do iglob = 1,NGLOB

           ! compute a Gaussian centered at the iglob point
           ! (part of the Gaussian may be outside the grid)
           xcen = x(iglob)
           zcen = z(iglob)
           k_gaus_global(:) = 0.
           do i = 1,NGLOB
              dist2 = (xcen - x(i))**2 + (zcen - z(i))**2
              if (dist2 <= dtrsh2) &
                 !k_gaus_global(i) = (1./(2.*PI*sigma**2)) * exp(-dist2 / (2.*sigma**2))
                 k_gaus_global(i) = (4./(PI*gamma**2)) * exp(-4.*dist2 / (gamma**2))
           enddo

!!$     ! compute the product of the Gaussian and the rough function, then sum
!!$     k_temp(:,:,:) = 0.
!!$     do ispec = 1,NSPEC
!!$       do j = 1,NGLLZ
!!$         do i = 1,NGLLX
!!$           itemp = ibool(i,j,ispec)
!!$           k_temp(i,j,ispec) = k_rough_global(itemp) * k_gaus_global(itemp) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
!!$         enddo
!!$       enddo
!!$     enddo
!!$     k_smooth_global(iglob) = sum( k_temp(:,:,:) )

           ! integrate the Gaussian over the grid
           ! (this will be constant only for the interior Gaussians)
           k_temp(:,:,:) = 0.
           do ispec = 1,NSPEC
              do j = 1,NGLLZ
                 do i = 1,NGLLX
                    itemp = ibool(i,j,ispec)
                    k_temp(i,j,ispec) = k_gaus_global(itemp) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
                 enddo
              enddo
           enddo
           k_gaus_int_global(iglob) = sum( k_temp(:,:,:) )

           ! compute the product of the Gaussian and the rough function, then sum
           ! normalization by k_gaus_int_global(iglob) accounts for Gaussians that are partially within the grid
           k_temp(:,:,:) = 0.
           do ispec = 1,NSPEC
              do j = 1,NGLLZ
                 do i = 1,NGLLX
                    itemp = ibool(i,j,ispec)
                    k_temp(i,j,ispec) = k_rough_global(itemp) * k_gaus_global(itemp) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
                 enddo
              enddo
           enddo
           k_smooth_global(iglob) = sum( k_temp(:,:,:) ) / k_gaus_int_global(iglob)

        enddo

        ! write smooth-related functions to file
        if (0 == 1) then
          file_smooth = 'fun_smooth.dat'
          open(unit=19,file=trim(out_dir1)//file_smooth,status='unknown')
          do iglob = 1,NGLOB
             write(19,'(8e16.6)') x_lon(iglob), z_lat(iglob), x(iglob), z(iglob), &
                  k_rough_global(iglob), k_gaus_global_ex(iglob), &
                  k_smooth_global(iglob), k_rough_global(iglob) - k_smooth_global(iglob)
          enddo
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

     endif  ! INV_STRUCT

     ! KEY: scaling parameter for structure for (joint) inversions
     mfac = 1.0
     if (istep == 0) then
        if (INV_SOURCE == 1 .and. INV_STRUCT == 1) then
           ! scale structure parameters according to source
           m_scale_str = mfac * sqrt( sum(source_gradient(:)*source_gradient(:)) )  &
              / sqrt( sum(k_smooth_global(:)*k_smooth_global(:) *sqrt(da(:))*sqrt(da(:)) ))

        else
           ! scale model parameters to order of 0.01 variations
           m_scale_str = sda_mean
        endif
     endif

     print *
     print *, ' Scaling parameter for structure gradient:'
     print *, '   mfac = ', mfac
     print *, '   F    = ', m_scale_str

     if (INV_STRUCT == 1) then
        ! KEY: compute gradient in 'irregular block' basis
        ! (special case of local Lagrange polynomial basis)
        do iglob = 1,NGLOB
           !gradient(iglob) = k_smooth_global(iglob)
           !gradient(iglob) = da(iglob) * k_smooth_global(iglob)
           gradient(iglob) = k_smooth_global(iglob) * sqrt(da(iglob)) * m_scale_str  ! old version
        enddo
     endif  ! INV_STRUCT

     ! fill the bottom of the model vector with source parameters
     ! (dx,dy,dt) * nevent
     if (INV_SOURCE == 1) gradient(nmod_str+1:nmod) = source_gradient(:)

     ! write gradient vector to file
     if (0 == 1) then
       open(unit=19,file=trim(out_dir1)//'gradient_vec.dat',status='unknown')
       do i = 1,nmod
         write(19,'(1e20.10)') gradient(i)
       enddo
       close(19)
     endif

  endif  ! itest == 0 .or. POLY_ORDER == 3

  !====================
  ! update search direction and model

  print *, ' Entering CG algorithm to compute new model or test model'

  if (itest == 0) then      ! if the present kernel is for a REAL model

     chi_k_val = chi_val
     gk(:) = gradient(:)

     ! the source model vector is always (0,0,0), since we are
     ! looking at perturbations relative to THE LATEST source position
     m0(:) = 0.

     if (INV_STRUCT == 1) then

       ! KEY: fractional pert from c0
       do iglob = 1,nmod_str
          !m0(iglob) = c_glob_syn(iglob) / c0 - 1.
          !m0(iglob) = (c_glob_syn(iglob) / c0 - 1.) * da(iglob)
          m0(iglob) = (c_glob_syn(iglob) / c0 - 1.) * sqrt(da(iglob)) / m_scale_str  ! old version
       enddo
     endif

     if (INV_SOURCE == 1) then
       ! scaled perturbation from initial source
       do i = 1,nmod_src
          m0(nmod_str+i) = (m_src(i) - m_src_syn(i)) / m_scale_src(i)
       enddo
     endif

     ! update search direction
     if (istep == 0) then
        pk(:) = -gk(:)     ! initial search direction

     else
        beta_val = dot_product( gk, gk(:) - g0(:) ) / dot_product( g0, g0 )
        pk(:) = -gk(:) + beta_val * p0(:)

     endif

     ! test value for line-search to get test model
     !istep_switch = 6
     !if (istep < istep_switch)  lam_t_val = -2.*chi_k_val / dot_product(gk, pk)    ! quadratic extrapolation
     !if (istep >= istep_switch) lam_t_val =    -chi_k_val / dot_product(gk, pk)    ! linear extrapolation
     lam_t_val = -2.*chi_k_val / dot_product(gk, pk)
     mt(:)     = m0(:) + lam_t_val * pk(:)
     !mt(:)     = m0(:) + lam_t_val * pk(:) / da(:)    ! structure only (g = K da)
     itest     = 1

     do i=1,nmod
       if (i <= nmod_str) then    ! structure

          ! get the new (test) structure model in terms of fractional perturbation
          if (INV_STRUCT == 0) then
             mt_vec(i) = c_glob_syn(i)  ! use same structure always

          else
             !mt_vec(i) = c0 * (1. + mt(i))
             !mt_vec(i) = c0 * (1. + mt(i)/ da(i) )
             mt_vec(i) = c0 * (1. + mt(i)/sqrt(da(i))*m_scale_str )    ! old version
          endif

       else                      ! source
          ! get the new source model in terms of (xs, zs, t0)
          if (INV_SOURCE == 0) then
             !mt_vec(i) = m0_vec(i)   ! use same source always
             mt_vec(i) = m_src_syn(i - nmod_str)

          else
             !mt_vec(i) = m0_vec(i) + mt(i)
             !mt_vec(i) = m0_vec(i) + mt(i)*m_scale_src(i - nmod_str)
             !mt_vec(i) = m_src_syn(i - nmod_str) * (1. + mt(i))

             mt_vec(i) = m_src_syn(i - nmod_str) + mt(i)*m_scale_src(i - nmod_str)

          endif

       endif
     enddo

     print *, 'lam_t_val = ', sngl(lam_t_val)

     ! save for next run
     g0(:) = gk(:)
     p0(:) = pk(:)

  else if (itest == 1) then  ! if present kernel is for a test model

     chi_t_val = chi_val

     ! a cubic or quadratic fit requires at least 5 values
     xx1 = 0.
     xx2 = lam_t_val
     yy1 = chi_k_val
     yy2 = chi_t_val
     g1  = dot_product(g0,pk)

     if (POLY_ORDER == 3) then
        ! use cubic polynomial: six values gives an analytical minimum
        ! see Matlab scripts cubic_min_4.m and cubic_min.m

        ! compute gradient of test function (needed for cubic polynomial)
        gt(:) = gradient(:)

        ! reset source vector to zero
        !m0(nmod_str+1:nmod) = 0.

        g2 = dot_product(gt,pk)

        ! coefficients of the cubic polynomial
        Pa = ( -2.*(yy2-yy1) + (g1+g2)*(xx2-xx1) ) / (xx2-xx1)**3
        Pb = ( 3.*(yy2-yy1) - (2.*g1 + g2)*(xx2-xx1) ) / (xx2-xx1)**2
        Pc = g1
        Pd = yy1

        ! get the analytical minimum
        qfac = Pb**2 - 3.*Pa*Pc;
        if (Pa /= 0 .and. qfac >= 0) then
           xmin = (-Pb + sqrt(qfac)) / (3.*Pa)
        else if (Pa == 0 .and. Pb /= 0) then
           xmin = -Pc/(2.*Pb)
        else
           stop 'check the input polynomial'
        endif

        ! write cubic function parameters to file
        open(unit=19,file=trim(out_dir1)//'cubic_poly.dat',status='unknown')
        write(19,'(11e16.6)') xx1,xx2,yy1,yy2,g1,g2,Pa,Pb,Pc,Pd,xmin
        close(19)

     else
        ! use quadratic polynomial: five values gives an analytical minimum
        ! see Matlab script quad_min_4.m

        ! coefficients of the quadratic polynomial (ax^2 + bx + c)
        Pa = ((yy2 - yy1) - g1*(xx2 - xx1)) / (xx2**2 - xx1**2)
        Pb = g1
        Pc = yy1 - Pa*xx1**2 - Pb*xx1

        ! get the analytical minimum (the vertex)
        if (Pa /= 0) then
           xmin = -Pb / (2.*Pa)
        else
           stop 'check the input polynomial'
        endif

        ! write quadratic function parameters to file
        open(unit=19,file=trim(out_dir1)//'quad_poly.dat',status='unknown')
        write(19,'(9e16.6)') xx1,xx2,yy1,yy2,g1,Pa,Pb,Pc,xmin
        close(19)

     endif  ! POLY_ORDER == 3

     ! update model
     lam_0_val = xmin
     !m0(:) = m0(:) + lam_0_val * pk(:) / da(:)   ! structure only (g = K da)
     m0(:) = m0(:) + lam_0_val * pk(:)
     itest = 0

     do i=1,nmod
       if (i <= nmod_str) then    ! structure

          ! get the new structure model in terms of fractional perturbation
          if (INV_STRUCT == 0) then
             m0_vec(i) = c_glob_syn(i)  ! use same structure always

          else
             !m0_vec(i) = c0 * (1. + m0(i))
             !m0_vec(i) = c0 * (1. + m0(i)/ da(i) )
             m0_vec(i) = c0 * (1. + m0(i)/sqrt(da(i))*m_scale_str )    ! old version
          endif

       else                      ! source
          ! get the new source model in terms of (xs, zs, t0)
          if (INV_SOURCE == 0) then
             !m0_vec(i) = m0_vec(i)     ! use same source always
             m0_vec(i) = m_src_syn(i - nmod_str)

          else
            !m0_vec(i) = m0_vec(i) + m0(i)
            !m0_vec(i) = m0_vec(i) + m0(i)*m_scale_src(i - nmod_str)
            !m0_vec(i) = m_src_syn(i - nmod_str) * (1. + m0(i))

            m0_vec(i) = m_src_syn(i - nmod_str) +  m0(i)*m_scale_src(i - nmod_str)

          endif

       endif
     enddo

     print *, 'lam_0_val = ', sngl(lam_0_val)

  endif

  !====================
  ! write updated/test model to file (source + structure)

  ! write CG vectors to file
  if (1 == 1) then
    open(unit=19,file=trim(out_dir1)//'cg_model_vectors.dat',status='unknown')
    do i = 1,nmod
      write(19,'(4e16.6)') m0(i), mt(i),  m0_vec(i), mt_vec(i)
    enddo
    close(19)
  endif
  if (0 == 1) then
    open(unit=19,file=trim(out_dir1)//'cg_grad_vectors.dat',status='unknown')
    do i = 1,nmod
      write(19,'(5e16.6)') g0(i), gt(i), gk(i), p0(i), pk(i)
    enddo
    close(19)
  endif

  ! exit program if model values are unrealistic
  ! NOTE: model parameters must be scaled appropriately
  if (itest == 1) then
    if ( minval(mt) < -10. .or. maxval(mt) > 10. ) stop 'test model is too extreme'
  else
    if ( minval(m0) < -10. .or. maxval(m0) > 10. ) stop 'updated model is too extreme'
  endif

  !====================

  !call system("xv model.jpg &")

!==================================
enddo  ! istep
!==================================

  if (ISOURCE_LOG) close(91)  ! source log file

  ! deallocate event and receiver variables

  deallocate(rglob,fglob)
  deallocate(ispec_rec,xi_rec,gamma_rec)
  deallocate(x_rec,z_rec,x_rec_lon,z_rec_lat)
  deallocate(hxir_store, hgammar_store)

  deallocate(eglob)
  deallocate(ispec_eve,xi_eve,gamma_eve)
  deallocate(x_eve,z_eve,x_eve_lon,z_eve_lat)
  deallocate(hxie_store,hgammae_store)

  deallocate(eglob_dat)
  deallocate(ispec_eve_dat,xi_eve_dat,gamma_eve_dat)
  deallocate(x_eve_dat,z_eve_dat,x_eve_lon_dat,z_eve_lat_dat)
  deallocate(hxied_store, hgammaed_store)

  deallocate(otime,otime_dat)

  ! CG variables

  deallocate(m0,mt,gradient)
  deallocate(g0,gt,gk,p0,pk)
  deallocate(source_gradient,m_scale_src)
  deallocate(m0_vec,mt_vec)
  deallocate(m_src,mt_src,m0_src,m_src_dat,m_src_syn)

!=======================
! COMMANDS FOR OMS RUNS
!=======================
!!$
!!$  ! forward wavefield
!!$  ! nsrc, sglob, samp : #src, src index, source
!!$  ! nrec, rglob, syn  : #rec, rec index, seis
!!$  isolver = 1
!!$  call solver(isolver, nsrc, sglob, samp, nrec, rglob, syn)
!!$
!!$  ! write out seismograms at the receivers
!!$  if (WRITE_SEISMO_F) call write_seismogram(syn, nrec, trim(out_dir)//'forward')
!!$
!!$  print *, 'compute and write out forward spectral map '
!!$  call write_spectral_map(syn, nrec, rglob, trim(out_dir)//'spectral_forward',WRITE_SPECTRA_F)
!!$
!!$  ! bandpass filter the seismograms to create the adjoint sources
!!$  !call filter(ti, syn, nrec)
!!$
!!$  ! nrecf, the number of (fake) receivers where we record the adjoint wavefield
!!$  ! is different from the number of original sources, nsrc.
!!$  ! this is a big memory burden -- perhaps we should write out the seismograms in blocks
!!$  deallocate(samp)
!!$  allocate(samp(NSTEP,NCOMP,nrecf))
!!$  samp(:,:,:) = 0.0
!!$
!!$  ! adjoint receiver = fake receivers
!!$  ! nrecf, fglob, samp : #a-rec, a-rec index, a-seis
!!$  ! nrec,  rglob, syn  : #a-src, a-src index, source
!!$  isolver = 2
!!$  call solver(isolver, nrecf, fglob, samp, nrec, rglob, syn)
!!$
!!$  ! write out adjoint seismograms at the fake receivers
!!$  if (WRITE_SEISMO_A) call write_seismogram(samp, nrecf, trim(out_dir)//'adjoint')
!!$
!!$  print *, 'compute and write out adjoint spectral map '
!!$  call write_spectral_map(samp, nrecf, fglob, trim(out_dir)//'spectral_adjoint',WRITE_SPECTRA_A)
!!$
!!$  ! deallocate variables
!!$  deallocate(sglob,samp,rglob,syn,fglob)

!==================================

!!$  ! adjoint receiver = original (forward) source
!!$  ! nsrc, sglob, samp : #a-rec, a-rec index, a-seis
!!$  ! nrec, rglob, syn  : #a-src, a-src index, source
!!$  call solver(isolver, nsrc, sglob, samp, nrec, rglob, syn)
!!$
!!$  ! write out adjoint seismograms at the original sources
!!$  if (WRITE_SEISMO_A) call write_seismogram(samp, nsrc, trim(out_dir)//'adjoint')

!==================================

  print *
  print *, ' max time-step is ', NSTEP
  print *, '  time-step DT is ', sngl(DT)
  print *, '      max time is ', sngl(NSTEP*DT), ' s,', sngl(NSTEP*DT/60), ' min'
  print *, ' space-step dh is ', sngl(dh)
  print *
  print *, '            c0 is ', sngl(c0)
  print *, '          hdur is ', sngl(hdur)
  print *, '    time shift is ', sngl(tshift)
  print *

  enddo  ! do iq

end program wave2d

