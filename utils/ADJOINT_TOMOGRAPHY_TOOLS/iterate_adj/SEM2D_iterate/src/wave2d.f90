!---------------------------------------------------------------------
! This program simulates elastic wave propagation on a 2D grid.
! The wavefield is stored as 3-component displacement records.
!
! It has the option of computed sensitivity kernels via adjoint methods.
!
! It also contains an embedded conjugate-gradient scheme for recovering
! "target models" from reference models.
!
! For more information, see the following papers:
!   Tromp, Tape, Liu, 2005, GJI.
!   Tape, Liu, Tromp, 2007, GJI.
!---------------------------------------------------------------------
! wave2d.f90
! Carl Tape, 13-Nov-2009
! Forward solver is based on earlier versions by Jeroen Tromp and Qinya Liu.
!
! The spectral-element method papers are by Dimitri Komatitsch, Jeroen Tromp, and others.
! Dimitri's SPECFEM2D code is the primary 2D SEM code on the SVN server.
!---------------------------------------------------------------------

program wave2d

  use wave2d_variables  ! declare global variables
  use wave2d_solver     ! mesher, set_model_property, solver
  use wave2d_sub        ! source time function, write seismograms, smoothing
  use wave2d_sub2       ! UTM conversion
  use wave2d_sub4       ! cross-correlation measurements and adjoint sources

  implicit none

  ! wavefields and seismograms
  double precision, dimension(:,:,:), allocatable :: samp, samp_dat, data, syn, adj_syn, data_recon

  ! source function
  double precision, dimension(NSTEP) :: stf_dat, stf_syn
  double precision :: f0(NCOMP), ti(NSTEP)

  ! socal coast and shelf data (for membrane surface wave simulations)
  integer :: FINITE_SOURCE
  integer :: nshelf,ncoast
  double precision :: d,dmin
  double precision, dimension(:), allocatable :: shelf_lat,shelf_lon,shelf_z,shelf_x
  double precision, dimension(:), allocatable :: coast_lat,coast_lon,coast_z,coast_x
  character(len=200) :: shelf_file,coast_file,idir

  ! corners of the grid  (for membrane surface wave simulations)
  double precision, dimension(4) :: corners_utm_x, corners_utm_z, corners_lon, corners_lat
  double precision :: sq_south, sq_east, sq_north, sq_west

  !------------------------
  ! source and receiver variables (and fake receivers for spectral plots)
  integer :: istf, nsrc, nrec, nrecf, nevent, nevent_dat, nevent_run, idata
  double precision :: dnevent_run
  !integer, dimension(:), allocatable :: sglob, sglob_dat
  !integer, dimension(:), allocatable :: sglob_temp, rglob_temp, fglob_temp

  integer, dimension(MAX_EVENT):: ifilter_eve_dat
  double precision, dimension(MAX_EVENT)      :: x_eve_lon0_dat, z_eve_lat0_dat, x_eve0_dat, z_eve0_dat
  double precision, dimension(:), allocatable :: x_eve_lon_dat,  z_eve_lat_dat,  x_eve_dat,  z_eve_dat, otime_dat
  double precision, dimension(:), allocatable :: xi_eve_dat, gamma_eve_dat
  integer, dimension(:), allocatable          :: eglob_dat, ispec_eve_dat

  integer, dimension(MAX_EVENT):: ifilter_eve_syn
  double precision, dimension(MAX_EVENT)      :: x_eve_lon0_syn, z_eve_lat0_syn, x_eve0_syn, z_eve0_syn
  double precision, dimension(:), allocatable :: x_eve_lon_syn,  z_eve_lat_syn,  x_eve_syn,  z_eve_syn, otime_syn
  double precision, dimension(:), allocatable :: xi_eve_syn, gamma_eve_syn
  integer, dimension(:), allocatable          :: eglob_syn, ispec_eve_syn

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

  ! velocity structure variables
  integer :: iref

  ! source perturbations (01-Feb-2006)
  ! presently only set up for 2D surface waves
  !logical :: ISOURCE_LOG
  !integer :: PERT_SOURCE, PERT_STRUCT, INV_SOURCE, INV_STRUCT_BETA
  logical :: stop_cg
  double precision :: src_pert_time, origin_time, origin_time_dat, origin_time_syn, ptmin, ptmax, ptime
  double precision :: src_pert_dist, rand1, amin, amax, pangle, rmin, rmax, prad
  double precision :: dx_src, dz_src, dt_src
  double precision, dimension(NVAR_STRUCT) :: m_scale_str
  double precision, dimension(NVAR_SOURCE) :: m_scale_src
  double precision, dimension(:,:,:,:), allocatable :: three_source_model
  double precision, dimension(:), allocatable :: source_gradient, m_scale_src_all   ! source_kernel
  double precision, dimension(:), allocatable :: gradient, gradient_data, gradient_model
  double precision, dimension(:,:), allocatable :: gradient_data_all

  !double precision, dimension(NGLLX,NGLLZ,NSPEC) :: gradient_local

  ! socal earthquake location
  double precision, dimension(MAX_EVENT) :: socal_events_lon, socal_events_lat, socal_events_mag
  double precision, dimension(MAX_EVENT) :: x_eve_syn_pert, z_eve_syn_pert, otime_eve_syn_pert, otime
  double precision, dimension(MAX_EVENT) :: x_eve_dat_pert, z_eve_dat_pert, otime_eve_dat_pert
  integer, dimension(MAX_EVENT) :: socal_events_lab
  integer nmeas, nmeas_run, ievent_min, ievent_max, ifirst_event

  !----------------------
  ! misfit and conjugate gradient algorithm

  ! weighting and smoothing
  double precision :: gamma  ! sigma
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: kbeta_smooth   ! kbulk_smooth
!  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: bulk_pert_syn, beta_pert_syn
  double precision :: scale_bulk, scale_beta, scale_source, scale_struct, scale_struct_gji
  double precision :: sigma_bulk0, sigma_beta0, sigma_bulk, sigma_beta, sigma_checker_scale
  double precision :: sigma_xs, sigma_zs, sigma_ts
  double precision :: cov_fac, cov_src_fac, cov_str_fac
  double precision :: afac_target, coverage_str, coverage_src
  double precision :: joint_str_src_grad_ratio, jsw_str, jsw_src, joint_total_weight, jfac
  double precision :: fac_str, fac_ts, fac_xs, fac_ys, fac_total_m ! fac_total_g
  double precision :: ugsq_str, ugsq_ts, ugsq_xs, ugsq_ys, ugsq_src
  double precision, dimension(:), allocatable :: cov_imetric, cov_model  ! icov_metric, icov_model

  ! conjugate gradient optimization (using line-search)
  double precision, dimension(:), allocatable :: pk, gk, g0, p0, m0, gt, mt, mdiff, mprior
  double precision :: beta_val, beta_val_top, beta_val_bot, lam_0_val, lam_t_val, lam_t_val_bot
  double precision :: chi_t_val, chi_k_val, mu_val, gnorm
  double precision :: Pa,Pb,Pc,Pd,qfac,xx1,xx2,yy1,yy2,g1,g2,xmin

  ! arrays describing the structure and source parameters
  ! OBSOLETE: m_src_syn_vec_initial, m_src_syn_initial, m_src_syn_vec, m_src_dat_vec, m0_vec_initial,
  double precision, dimension(:), allocatable :: m0_vec, mt_vec, mtarget
  double precision, dimension(:), allocatable :: m_src_syn, m_src_dat, m_src_prior

  double precision :: vel_min, vel_max
  integer :: itest, nmod, nmod_src, nmod_str

  !----------------------

  ! kernels
  double precision, dimension(:), allocatable :: tstart, tend
  character(len=200) :: last_frame_name
  double precision, dimension(:,:,:,:,:), allocatable :: absorb_field
  double precision, dimension(:,:,:), allocatable :: atype_kernel, btype_kernel, rtype_kernel
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: btype_kernel_sum
  double precision, dimension(NLOCAL) :: ktemp

  double precision :: xtemp,ztemp,dx,dz,xmesh,zmesh
  double precision :: t_target, tcen
  double precision :: junk1,junk2,junk3
  double precision :: temp1,temp2,temp3,temp4,temp5,temp6,temp7

  double precision :: da_local_mean, sda_local_mean, da_global_mean, sda_global_mean
  double precision :: mfac, norm_source_2, norm_bulk_2, norm_beta_2, norm_struct_2

  character(len=20)  :: data_tag,syn_tag,stf_tag,stfadj_tag
  character(len=200) :: srcfile,recfile,socal_map
  character(len=200) :: socal_dir, model_dir, script_dir, out_dir3, out_dir2, out_dir1, ev_dir
  character(len=200) :: filename, filename1, filename2, filename3, file_structure_dat, file_structure_syn, file_src_rec, file_stf
  character(len=200) :: command1

  integer :: itemp,itemp1,itemp2,itemp3, hmod
  integer :: i, j, k, iq, itime, iglob, irec, isrc, ipick, ievent, icomp  ! irunz
  integer :: isolver, irun0, irun, iopt, ispec, istep, istep_switch, imod, ipar, imnorm

  !double precision :: meas_pert, rand1,  meas, ppert

  integer :: GJI_PAPER, qmax

  !-------------------------------------------------------------------------------------
  !********* PROGRAM STARTS HERE *********************
  !-------------------------------------------------------------------------------------

  print *, 'Starting wave2d.f90 ...'

  !--------------------------------------
  ! stop program if there are certain unallowed paramter combinations

  if (READ_IN) then
      print *, 'For READ_IN runs...'
      if (NITERATION /= 0) stop 'NITERATION = 0'
      if (ISRC_SPACE /= 6) stop 'ISRC_SPACE = 6'
      hmod = 1     ! iteration number for read-in models
  endif

  if ( IKER /= 0 .and. IKER /= 1 .and. IKER /= 2 ) stop 'IKER must = 0,1,2'

  if ( ADD_DATA_ERRORS .and. IKER /= 1 ) stop 'for now, IKER = 1 for adding errors to data'

  ! if you do not want kernels, then you cannot iterate
  if ( (.not. COMPUTE_KERNELS) .and. NITERATION /= 0 ) stop 'NITERATION = 0'

  !--------------------------------------

  ! random number generator (for including measurement errors)
  !call random_seed()

  ! specify input and output directories
  out_dir3   = "OUTPUT/"
  in_dir     = "INPUT/"

  ! NOTE: check that these directories exist
  !model_dir  = 'model_files/'
  !script_dir = "scripts/"
  !socal_dir = '/net/denali/scratch1/carltape/socal_modes/socal_modes/local_modes/'

  ! whether a point source or finite source is used
  if (ISRC_SPACE == 1) FINITE_SOURCE = 0
  if (ISRC_SPACE /= 1) FINITE_SOURCE = 1

  origin_time     = tshift       ! reference origin time (s)
  origin_time_dat = origin_time
  origin_time_syn = origin_time
  otime(:) = origin_time         ! origin time for prior model

  ! STANDARD DEVIATIONS of source perturbations (see wave2d_sigmas.m)
  ! NOTE: THESE ARE OVER-WRITTEN FOR ISURFACE=1
  src_pert_dist = 2000.0       ! source spatial perturbation (m)
  src_pert_time = 0.5            ! source temporal perturbation (s)

  !--------------------------------------
  ! load or specify events, depending on whether you are computing
  ! membrane surface waves (ISURFACE==1) or body waves

  if (ISURFACE == 1) then

     !----------
     ! read in target events (socal_quakes.m)
     ! socal events M > 4 (1990-2005), selected to get semi-uniform coverage

     !quake_file = 'socal_quakes_v02.dat'
     quake_file = 'socal_quakes_N025.dat'
     open(55,file=trim(in_dir)//trim(quake_file),status='unknown')
     read(55,*) nevent
     if (nevent > MAX_EVENT) stop 'nevent > MAX_EVENT (so increase MAX_EVENT)'
     print *, nevent, ' target events (for synthetics)'
     read(55,'(i14,3f14.7)') (socal_events_lab(i), &
        socal_events_lon(i), socal_events_lat(i), socal_events_mag(i), i = 1,nevent)
     close(55)

     print *
     print *, 'SoCal events (label, lon, lat, mag):'
     do i = 1,nevent
        write(*,'(i14,3f14.7)') socal_events_lab(i), socal_events_lon(i), socal_events_lat(i), socal_events_mag(i)
     enddo

     !----------
     ! read in perturbations for initial and target models

     otime_eve_syn_pert(:) = 0.0 ; x_eve_syn_pert(:) = 0.0 ;  z_eve_syn_pert(:) = 0.0
     otime_eve_dat_pert(:) = 0.0 ; x_eve_dat_pert(:) = 0.0 ;  z_eve_dat_pert(:) = 0.0

     if (0 == 1) then       ! NEW perturbations

        ! perturbations generated from wave2d_sigmas.m
        open(19,file='INPUT/events_txy_pert_initial.dat',status='old')
        do ievent = 1,nevent
           read(19,*) temp1,temp2,temp3
           otime_eve_syn_pert(ievent) = temp1
           x_eve_syn_pert(ievent) = temp2
           z_eve_syn_pert(ievent) = temp3
        enddo
        close(19)

        ! perturbations generated from wave2d_sigmas.m
        open(20,file='INPUT/events_txy_pert_target.dat',status='old')
        do ievent = 1,nevent
           read(20,*) temp1,temp2,temp3
           otime_eve_dat_pert(ievent) = temp1
           x_eve_dat_pert(ievent) = temp2
           z_eve_dat_pert(ievent) = temp3
        enddo
        close(20)

        ! read sigma values used to generate the source parameter perturbations
        ! NOTE: these over-write what is defined above
        open(21,file='INPUT/events_txy_pert_sigmas.dat',status='old')
        read(21,*) src_pert_time,src_pert_dist,temp1
        close(21)

     else              ! OLD perturbations

        ! perturbations generated from gji_src_inversion.m
        ! NOTE: OLD ordering scheme for source indexing
        open(19,file='INPUT/events_xyt_pert.dat',status='unknown')
        do ievent = 1,nevent
           read(19,*) temp1,temp2,temp3
           !if (ievent==5) then
           !   temp1 = -1.599978278510d3
           !   temp2 = -6.502537573573d2
           !   temp3 = 7.975610515441d-2
           !endif
           x_eve_dat_pert(ievent) = temp1
           z_eve_dat_pert(ievent) = temp2
           otime_eve_dat_pert(ievent) = temp3
        enddo
        close(19)

        open(20,file='INPUT/events_xyt_pert_sigmas.dat',status='unknown')
        read(20,*) temp1,src_pert_dist,src_pert_time
        close(20)

     endif

     ! no perturbations for synthetics --> m0 = mprior
     if (M0ISMPRIOR) then
        otime_eve_syn_pert(:) = 0.0
        x_eve_syn_pert(:) = 0.0
        z_eve_syn_pert(:) = 0.0
     endif

!!$     ! check
!!$     do ievent = 1,nevent
!!$        write(*,'(a20,i6,3f12.3)') 'Initial model :',ievent, &
!!$             otime_eve_syn_pert(ievent), x_eve_syn_pert(ievent), z_eve_syn_pert(ievent)
!!$        write(*,'(a20,i6,3f12.3)') 'Target model :',ievent, &
!!$             otime_eve_dat_pert(ievent), x_eve_dat_pert(ievent), z_eve_dat_pert(ievent)
!!$     enddo
!!$     stop 'testing perturbations'

     !----------
     ! load socal coast and shelf points
     ! (This is useful for some ocean microseismi time-rerversal experiments.)

     ! exclude sources or receivers that are this close to THE SOCAL COAST
     dmin_trsh_s = 0.0d+03                  ! sources
     !dmin_trsh_r = STATION_GRID_BUFFER     ! receivers
     dmin_trsh_r = 0.0d+03
     dmin_trsh_f = 0.0d+03                  ! fake receivers

     ! read in data files
     nshelf = 101
     ncoast = 810
     allocate(shelf_lat(nshelf),shelf_lon(nshelf),shelf_z(nshelf),shelf_x(nshelf))
     allocate(coast_lat(ncoast),coast_lon(ncoast),coast_z(ncoast),coast_x(ncoast))
     shelf_file = 'oms_shelf'
     coast_file = 'oms_coast'
     open(77,file=trim(in_dir)//'BOUNDARIES/'//trim(shelf_file),status='unknown')
     read(77,*) (shelf_lat(i),shelf_lon(i),i = 1,nshelf)
     close(77)
     open(66,file=trim(in_dir)//'BOUNDARIES/'//trim(coast_file),status='unknown')
     read(66,*) (coast_lat(i),coast_lon(i),i = 1,ncoast)
     close(66)

     ! convert lat-lon to mesh coordinates
     call mesh_geo(nshelf,shelf_lon,shelf_lat,shelf_x,shelf_z,UTM_PROJECTION_ZONE,ILONLAT2MESH)
     call mesh_geo(ncoast,coast_lon,coast_lat,coast_x,coast_z,UTM_PROJECTION_ZONE,ILONLAT2MESH)
     !write(*,'(4f16.6)') (shelf_lon(i), shelf_lat(i), shelf_x(i)/1000.0, shelf_z(i)/1000.0, i = 1,nshelf)
     !write(*,'(4f16.6)') (coast_lon(i), coast_lat(i), coast_x(i)/1000.0, coast_z(i)/1000.0, i = 1,ncoast)

  else

     !----------
     ! events for data

     ! all vectors should be this length
     nevent = MAX_EVENT

     if (0 == 1) then    ! read in target events (socal1D.m)

        quake_file = 'socal_quakes_1D_rand10.dat'   ! 200 km width
        !quake_file = 'socal_quakes_1D_rand15.dat'   ! 400 km width
        open(55,file=trim(in_dir)//trim(quake_file),status='unknown')
        read(55,*) nevent
        if (nevent > MAX_EVENT) stop 'nevent > MAX_EVENT (so increase MAX_EVENT)'
        print *, nevent, ' target events (for synthetics)'
        read(55,'(i10,3f24.12)') (junk1,x_eve0_dat(i),z_eve0_dat(i),junk2,i = 1,nevent)
        close(55)

        print *
        print *, 'Target events (index, x, z):'
        do i = 1,nevent
           write(*,*) i,x_eve0_dat(i),z_eve0_dat(i)
        enddo

     else              ! specify target events here

        ! target events (data)
        x_eve0_dat(:) = -99.0            ; z_eve0_dat(:) = -99.0
        !x_eve0_dat(1) = 50.0d3   ; z_eve0_dat(1) = 0.75 * HEIGHT
        !x_eve0_dat(2) = 50.0d3   ; z_eve0_dat(2) = 0.50 * HEIGHT
        !x_eve0_dat(3) = 50.0d3   ; z_eve0_dat(3) = 0.25 * HEIGHT

        !x_eve0_dat(1) = 50.0d3   ; z_eve0_dat(1) = HEIGHT - 2.5d3
        !x_eve0_dat(2) = 50.0d3   ; z_eve0_dat(2) = HEIGHT - 10.0d3
        !x_eve0_dat(3) = 50.0d3   ; z_eve0_dat(3) = HEIGHT - 20.0d3
        !x_eve0_dat(4) = 50.0d3   ; z_eve0_dat(4) = HEIGHT - 40.0d3

        x_eve0_dat(1) = LENGTH/2.0   ; z_eve0_dat(1) = HEIGHT

     endif

     !----------
     ! define 1D model for body waves

     ! SoCal 1D model (Dreger & Helmberger 1990)
     !z_breaks(1) = 5500.0 ; z_breaks(2) = 16000.0 ; z_breaks(3) = 35000.0 ; nlayer = 4

     ! model to honor discontinuities (32 elements in 80 km depth)
     z_breaks(1) = 5000.0 ; z_breaks(2) = 15000.0 ; z_breaks(3) = 35000.0 ; nlayer = 4

     ! four-layer 1D model
     r_layers(1) = 2400.0 ; r_layers(2) = 2670.0 ; r_layers(3) = 2800.0 ; r_layers(4) = 3000.0
     a_layers(1) = 5500.0 ; a_layers(2) = 6300.0 ; a_layers(3) = 6700.0 ; a_layers(4) = 7800.0
     b_layers(1) = 3180.0 ; b_layers(2) = 3640.0 ; b_layers(3) = 3870.0 ; b_layers(4) = 4500.0

     !alpha_min = a_layers(1)
     !alpha_max = a_layers(nlayer)
     !beta_min  = b_layers(1)
     !beta_max  = b_layers(nlayer)
  endif

  !--------------------------------------

  ! events to use (5-5, 1-5, or 1-25)
  !ievent_min = 1  ; ievent_max = nevent
  !ievent_min = 1 ; ievent_max = 5
  ievent_min = 5  ; ievent_max = ievent_min
  nevent_run = ievent_max - ievent_min + 1

  ! number of sources in inversion
  dnevent_run = dble(nevent_run)

  ! MAX source perturbations for location (lat-lon) and origin time
  !    dt0 = dt0_dat - dt0_syn ( < 0 for delayed synthetics)
  !src_pert_dist = 5000.0      ! source spatial perturbation (m)
  !src_pert_time = 1.0         ! source temporal perturbation (s)

  ! reference directory for the run
  !irunz = 500

  GJI_PAPER = 0

  ! if you are not computing kernels, loop over 25 TSVD models to read in
  qmax = 25
  if (COMPUTE_KERNELS .or. (READ_IN .and. READ_SINGLE) ) qmax = 1

  !============================================
  ! LOOP 1: different tomographic runs
  do iq = 1,qmax
  !do iq = 1,250       ! GAUSS
  !============================================

     ! KEY COMMAND: scalelength of checker for velocity models (1,2,3)
     Nfac = 3   ! use Nfac=3 for one-source examples

     if (READ_IN) then
        irun0 = IRUNZ
     else
        irun0 = IRUNZ + 20*(iq-1)    ! increment is 20
     endif

     !irun0 = IRUNZ + iq   ! GAUSS

     ! name the reference output directory for the optimization run
     if (READ_IN) then

        ! inversion steps done in Matlab
        if (READ_SINGLE) then
           write(out_dir2,'(a,i4.4,a,i4.4,a)') trim(out_dir3)//'run_',irun0,'/READ_IN/model_m',hmod,'/'
        else
           write(out_dir2,'(a,i4.4,a,i4.4,a,i4.4,a)') trim(out_dir3)//'run_',irun0,'/READ_IN/model_m',hmod,'/run_p',iq,'/'
        endif

     else
        ! conjugate-gradient-based inversion done in wave2d.f90
        write(out_dir2,'(a,i4.4,a)') trim(out_dir3)//'run_',irun0,'/'
        command1 = 'mkdir ' // trim(out_dir2)
        call system(command1)
        print *, command1
     endif

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

     ! copy wave2d_constants.f90 and wave2d.f90 into output directory
     command1 = 'cp ' // 'src/wave2d_constants.f90 ' // trim(out_dir2)
     call system(command1)
     command1 = 'cp ' // 'src/wave2d.f90 ' // trim(out_dir2)
     call system(command1)

     ! parameters for plotting scripts
     call write_parameters_plot(trim(out_dir2)//'parameters1.log')

     ! all parameters
     call write_parameters(trim(out_dir2)//'wave2d_constants.dat')

     !--------------------------------------
     ! mesher

     ! set up model (gets Jacobian)
     call mesher()

     ! compute da_global, da_local, valence, ielement_corner
     call mesher_additional()

     da_global_mean  = sum(da_global)/NGLOB
     da_local_mean   = sum(da_local)/NLOCAL
     sda_global_mean = sqrt(da_global_mean)
     sda_local_mean  = sqrt(da_local_mean)

     if (1 == 1) then
        ! global mesh
        open(unit=15,file=trim(out_dir2)//'global_mesh.dat',status='unknown')
        do iglob = 1,NGLOB
           write(15,'(2i8,3e18.8)') iglob, valence(iglob), x(iglob), z(iglob), da_global(iglob)
        enddo
        close(15)

        ! local mesh (see subroutine mesher_additional in wave2d_solver.f90)
        open(unit=15,file=trim(out_dir2)//'local_mesh.dat',status='unknown')
        k = 0
        do ispec = 1,NSPEC
           do j = 1,NGLLZ
              do i = 1,NGLLX
                 k = k+1
                 iglob = ibool(i,j,ispec)
                 write(15,'(6i8,4e18.8)') k, ispec, i, j, iglob, valence(iglob), &
                      x(iglob), z(iglob), da_local(i,j,ispec), da_global(iglob)
              enddo
           enddo
        enddo
        close(15)
     endif

     if (0 == 1) then
        ! corner points for each element, and centerpoint (in km)
        open(unit=15,file='elements.dat',status='unknown')
        do ispec = 1,nspec
           xtemp = (x1(ispec) + x2(ispec))/2.0
           ztemp = (z1(ispec) + z2(ispec))/2.0
           write(15,'(i8,6f14.6)') ispec,xtemp/1000.0,ztemp/1000.0, &
              x1(ispec)/1000.0,x2(ispec)/1000.0,z1(ispec)/1000.0,z2(ispec)/1000.0
        enddo
        close(15)

        ! GLL points for one element (in km)
        ispec = 292  ! pick an element
        open(unit=15,file='gll_points.dat',status='unknown')
        do j = 1,NGLLZ
           do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              write(15,'(4i10,3e18.8)') ispec, i, j, valence(iglob), da_global(iglob)/1e6, x(iglob)/1000.0, z(iglob)/1000.0
           enddo
        enddo
        close(15)

        ! GLL points defining the corners of elements
        open(unit=15,file='element_corners.dat',status='unknown')
        do i = 1,NSPEC_CORNER
           iglob = ielement_corner(i)
           write(15,'(1i10,2e18.8)') iglob, x(iglob)/1000.0, z(iglob)/1000.0
        enddo
        close(15)

        stop 'checking the GLL indexing on the elements'
     endif

     print *
     write(*,*)             ' PROPERTIES OF THE GRID :'
     write(*,*)             '     GLOBAL          :'
     write(*,'(a,1f20.10)') '       da_min  (m^2) : ', minval(da_global)
     write(*,'(a,1f20.10)') '       da_mean (m^2) : ', da_global_mean
     write(*,'(a,1f20.10)') '       da_max  (m^2) : ', maxval(da_global)
     write(*,*)             '        sum [ da ] : ', sum(da_global)
     write(*,*)             '     LOCAL           :'
     write(*,'(a,1f20.10)') '       da_min  (m^2) : ', minval(da_local)
     write(*,'(a,1f20.10)') '       da_mean (m^2) : ', da_local_mean
     write(*,'(a,1f20.10)') '       da_max  (m^2) : ', maxval(da_local)
     write(*,*)             '         sum [ da  ] : ', sum(da_local)
     print *
     write(*,*)             '     LENGTH * HEIGHT : ', LENGTH * HEIGHT
     print *
     write(*,*)             ' GAMMA_SMOOTH_KERNEL : ', GAMMA_SMOOTH_KERNEL
     write(*,*)             '                Nfac : ', Nfac
     write(*,*)             '               irun0 : ', irun0
     write(*,*)             '                IKER : ', IKER
     print *
     print *, ' GLL weights:'
     do i = 1,NGLLX
        print *, wxgll(i)
     enddo
     do i = 1,NGLLZ
        print *, wzgll(i)
     enddo

     !stop 'testing'

     !--------------------------------------

     ! determine the UTM coordinates of your origin and corners
     if (ISURFACE == 1) then

        utm_xmin = 0.0 ; utm_zmin = 0.0
        xtemp = LON_MIN ; ztemp = LAT_MIN
        call utm_geo(xtemp,ztemp,utm_xmin,utm_zmin,UTM_PROJECTION_ZONE,ILONGLAT2UTM)
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
        do i = 1,4
           call utm_geo(xtemp,ztemp,corners_utm_x(i),corners_utm_z(i),UTM_PROJECTION_ZONE,IUTM2LONGLAT)
           corners_lon(i) = xtemp
           corners_lat(i) = ztemp
           write(*,'(2f16.6,2e16.6)') corners_lon(i), corners_lat(i), corners_utm_x(i), corners_utm_z(i)
        enddo

        ! convert global gridpoint mesh coordinates to lat-lon
        x_lon(:) = 0.0
        z_lat(:) = 0.0
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
        xtemp = LON_MIN
        ztemp = LAT_MIN
        call utm_geo(xtemp,ztemp,utm_xmin,utm_zmin,UTM_PROJECTION_ZONE,ILONGLAT2UTM)
        print *, utm_xmin, utm_zmin
        call utm_geo(xtemp,ztemp,utm_xmin,utm_zmin,UTM_PROJECTION_ZONE,IUTM2LONGLAT)
        print *, xtemp,ztemp   ! should match LON_MIN,LAT_MIN
        print *

     endif

     !--------------------------------------
     ! assign values for plotting grid

     if (ISURFACE == 1) then
         x_plot(:) = x_lon(:)
         z_plot(:) = z_lat(:)

     else if (ISURFACE == 0) then
         x_plot(:) = x(:)
         z_plot(:) = z(:)

     endif

     !--------------------------------------
     ! phase velocity model

     t_target = 2*hdur     ! target period for phase velocity map

!!$  ! homogeneous or heterogeneous phase velocity map
!!$
!!$  if (IMODEL_DAT == 1 .or. IMODEL_SYN == 1) then     ! read or compute a heterogeneous map
!!$     ! This calls model_files/get_socal_map.csh, which executes wave2d_socal.f90
!!$     ! from /net/denali/scratch1/carltape/socal_modes/local_modes
!!$
!!$     if (0==1) then    ! compute phase velocity map
!!$
!!$        ! write lat-lon gridpoints to file
!!$        filename = trim(socal_dir) // 'socal_input.dat'
!!$        print *, 'Output data file is ',trim(filename)
!!$        open(unit=15,file=filename,status='unknown')
!!$        write(15,*) t_target
!!$        write(15,*) 1               ! formerly IHOMO, IMODEL
!!$        write(15,*) UTM_PROJECTION_ZONE
!!$        write(15,*) NGLOB
!!$        write(15,*) (x_lon(iglob),z_lat(iglob),iglob=1,NGLOB)
!!$        close(15)
!!$
!!$        ! run wave2d_socal.f90 to create output file
!!$        ! KEY NOTE: must run this each time you change the grid, hdur, or IMODEL_DAT
!!$        call system(trim(model_dir) // 'get_socal_map.csh')
!!$
!!$        ! read in phase velocity map
!!$        !write(infile,'(a,a,i4.4,a)') trim(socal_dir), 'socal_T', int(100*t_target),'.dat'
!!$        socal_map = 'socaltest.dat'
!!$        open(unit=16, file = trim(model_dir)//socal_map, status='unknown')
!!$        read(16,*) (junk1,junk2,c_glob(iglob),junk3,iglob=1,NGLOB)
!!$        close(16)
!!$
!!$        ! read in REFERENCE phase velocity (km/s)
!!$        filename = trim(model_dir) // 'socaltest2.dat'
!!$        open(unit=17,file=filename,status='unknown')
!!$        read(17,*) beta0
!!$        close(17)
!!$
!!$        ! convert km/s --> m/s
!!$        c_glob(:) = 1000.*c_glob(:)
!!$        beta0 = beta0*1000.
!!$
!!$     else             ! read in phase velocity map for data
!!$
!!$        open(unit=16,file=trim(model_dir)//'socal_vel_dat.dat', status='unknown')
!!$        read(16,*) (junk1,junk2,c_glob(iglob),junk3,iglob=1,NGLOB)
!!$        close(16)
!!$
!!$     endif
!!$
!!$     ! compute beta0 by finding the average of the phase velocity map
!!$     print *, 'computing the average of the phase velocity map'
!!$     beta0 = 0.0
!!$     do iglob = 1,NGLOB
!!$        beta0 = beta0 + c_glob(iglob) * da_global(iglob)
!!$     enddo
!!$     beta0 = beta0 / (LENGTH * HEIGHT)
!!$
!!$  endif
!!$
!!$  if (IMODEL_SYN == 1) c_glob_syn(:) = c_glob(:)
!!$  if (IMODEL_DAT == 1) c_glob_dat(:) = c_glob(:)

     ! reference phase velocity (m/s) (should be obtained based on hdur)
     ! (This is only relevant for homogeneous models.)
     if (ISURFACE == 1) then
        alpha0 = 5800.0
        beta0 = 3500.0    ! (3500 for GJI)
        rho0 = DENSITY

     else if (ISURFACE == 0) then
        !alpha0 = a_layers(2)    ! layer 2 of SoCal model (1D)
        !beta0  = b_layers(2)    ! layer 2 of SoCal model (1D)
        alpha0 = (a_layers(2)+a_layers(1))/2
        beta0  = (b_layers(2)+b_layers(1))/2
     endif

     ! parameters for checkerboard pattern
     if (IMODEL_DAT == 2 .or. IMODEL_SYN == 2) then
        !Nfac = 3      ! scalelength of map = N * (wavelength of surface waves for hdur)
        afac = 10.0   ! max PERCENT variation of synthetic model  (10.0)

        ! spatial frequency corresponding to scalelength
        ! factor of 2 in demoninator is because the scalelength is a 1/2 cycle
        w_scale = 2*pi/(Nfac*2*beta0*t_target)
     endif

     !----------
     ! load in heterogeneous reference and target models
     ! NOTE: THESE ARE ONLY FOR A SPECIFIC MESH
     if ((IMODEL_SYN == 3) .and. (IMODEL_DAT == 3)) then

        idir = '/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_INPUT/random_fields/'

        ! read in structure model for synthetics
        write(filename2,'(a)') trim(idir)//'model_0001/structure_syn_in.dat'
        !write(filename2,'(a,i4.4,a)') trim(idir)//'model_set_0001/structure_model_',iq,'.dat'
        open(unit=19,file=filename2,status='old')
        do ispec = 1,NSPEC
           do j = 1,NGLLZ
              do i = 1,NGLLX
                 read(19,'(6e20.10)') temp1, temp2, temp3, temp4, temp5, temp6
                 kappa_syn(i,j,ispec) = temp3
                 mu_syn(i,j,ispec) = temp4
                 rho_syn(i,j,ispec) = temp5
              enddo
           enddo
        enddo
        close(19)

        ! read in structure model for data
        write(filename2,'(a)') trim(idir)//'model_0001/structure_dat_in.dat'
        open(unit=20,file=filename2,status='old')
        do ispec = 1,NSPEC
           do j = 1,NGLLZ
              do i = 1,NGLLX
                 read(20,'(6e20.10)') temp1, temp2, temp3, temp4, temp5, temp6
                 kappa_dat(i,j,ispec) = temp3
                 mu_dat(i,j,ispec) = temp4
                 rho_dat(i,j,ispec) = temp5
              enddo
           enddo
        enddo
        close(20)

        ! read in sigma value describing the distribution
        write(filename,'(a)') trim(idir)//'model_0001/structure_sigma.dat'
        open(unit=21,file=filename,status='old')
        read(21,*) sigma_beta0
        close(21)
        sigma_bulk0 = sigma_beta0

        ! compute beta for initial and target models
        beta_syn  = sqrt( mu_syn / rho_syn )
        beta_dat  = sqrt( mu_dat / rho_dat )

        ! read reference values from file -- these were used to make the perturbed field
        idir = '/data1/cig/seismo/3D/ADJOINT_TOMO/iterate_adj/SEM2D_iterate_INPUT/meshes/mesh_0001/'
        open(unit=16,file=trim(idir)//'reference_values.dat', status='unknown')
        read(16,'(6e20.10)') alpha0, beta0, rho0, bulk0, kappa0, mu0
        close(16)

     endif

     !=========================================
     ! (INITIAL) REFERENCE MODEL (synthetics)

     ! initialize
     !kappa_syn = 0.0 ;   mu_syn = 0.0 ;  rho_syn = 0.0
     !alpha_syn = 0.0 ; beta_syn = 0.0 ; bulk_syn = 0.0

     ! set velocity field for the synthetics (e.g., homogeneous reference model)
     ! assigns local arrays mu_syn, kappa_syn, rho_syn
     iref = 1
     call set_model_property(iref)

     ! smooth the reference model
     ! (This was motivated for the imaging principle kernel tests, Feb 2007.)
     if (ISMOOTH_INITIAL_MODEL == 1) then
        call smooth_function(kappa_syn, GAMMA_SMOOTH_MODEL, temp_local1)  ; kappa_syn = temp_local1;
        call smooth_function(   mu_syn, GAMMA_SMOOTH_MODEL, temp_local1)  ;    mu_syn = temp_local1;
        call smooth_function(  rho_syn, GAMMA_SMOOTH_MODEL, temp_local1)  ;   rho_syn = temp_local1;
        bulk_syn  = sqrt( kappa_syn / rho_syn )
        beta_syn  = sqrt( mu_syn / rho_syn )
        alpha_syn = sqrt( (kappa_syn + FOUR_THIRDS*mu_syn) / rho_syn )
     endif

     ! compute alpha0 and beta0 by finding the average of the phase velocity map
     ! (only if you are not reading in these values)
     if ( (IMODEL_SYN /= 3) .or. (IMODEL_DAT /= 3) ) then
        print *, 'computing the average of the phase velocity map'
        kappa0 = sum(kappa_syn * da_local) / (LENGTH * HEIGHT)
        mu0    = sum(   mu_syn * da_local) / (LENGTH * HEIGHT)
        rho0   = sum(  rho_syn * da_local) / (LENGTH * HEIGHT)
        alpha0 = sqrt( (kappa0 + FOUR_THIRDS*mu0) / rho0 )
        beta0  = sqrt(    mu0 / rho0 )
        bulk0  = sqrt( kappa0 / rho0 )
        !alpha0 = sum(alpha_initial_syn * da_local) / (LENGTH * HEIGHT)
        !beta0  = sum(beta_initial_syn * da_local)  / (LENGTH * HEIGHT)
     endif

!!$     ! initial and reference velocity structure
!!$     ! From this model, we compute fractional perturbations.
!!$     bulk_initial_syn(:,:,:) = bulk_syn(:,:,:)
!!$     beta_initial_syn(:,:,:) = beta_syn(:,:,:)

     !=========================================
     ! TARGET MODEL (data)

     ! initialize
     !kappa_dat = 0.0 ;   mu_dat = 0.0 ;  rho_dat = 0.0
     !alpha_dat = 0.0 ; beta_dat = 0.0 ; bulk_dat = 0.0

     ! set velocity field for the data (perturbed reference model)
     ! assigns local arrays mu_dat, kappa_dat, rho_dat
     iref = 0
     call set_model_property(iref)

     ! if READ_IN option, then READ in the structure file (local level)
     ! NOTE: Assume that the mesh and scaling values are IDENTICAL
     !       to what was used in the base directory for the CG algorithm.
     if ( (READ_IN) .and. (INV_STRUCT_BETA == 1) ) then

        kappa_syn = 0.0 ;   mu_syn = 0.0 ;  rho_syn = 0.0
        alpha_syn = 0.0 ; beta_syn = 0.0 ; bulk_syn = 0.0

        ! read in structure model for synthetics
        write(filename2,'(a,i4.4,a)') trim(out_dir2)//'structure_syn_m',hmod,'.dat'
        open(unit=19,file=filename2,status='unknown')
        do ispec = 1,NSPEC
           do j = 1,NGLLZ
              do i = 1,NGLLX
                 ! CURRENT MODEL (synthetics)
                 !read(19,*) temp1, temp2, &
                 !     kappa_syn(i,j,ispec), mu_syn(i,j,ispec), rho_syn(i,j,ispec), &
                 !     temp3, temp4

                 read(19,'(6e20.10)') temp1, temp2, temp3, temp4, temp5, temp6
                 kappa_syn(i,j,ispec) = temp3
                 mu_syn(i,j,ispec) = temp4
                 rho_syn(i,j,ispec) = temp5

                 !print *, temp4, mu_syn(i,j,ispec)

                 !print *, filename2
                 !stop 'CHECK THE FIRST ENTRY'

              enddo
           enddo
        enddo
        close(19)

        ! smooth the input model
        if (ISMOOTH_MODEL_UPDATE == 1) then
           call smooth_function(kappa_syn, GAMMA_SMOOTH_MODEL, temp_local1)  ; kappa_syn = temp_local1;
           call smooth_function(   mu_syn, GAMMA_SMOOTH_MODEL, temp_local1)  ;    mu_syn = temp_local1;
           call smooth_function(  rho_syn, GAMMA_SMOOTH_MODEL, temp_local1)  ;   rho_syn = temp_local1;
        endif

        alpha_syn = sqrt( (kappa_syn + FOUR_THIRDS*mu_syn) / rho_syn )
        beta_syn  = sqrt( mu_syn / rho_syn )
        bulk_syn  = sqrt( kappa_syn / rho_syn )
     endif

     ! unperturbed structure (set target map to be the map used to compute synthetics)
     if (PERT_STRUCT_BETA == 0) then
        mu_dat    = mu_syn
        beta_dat  = beta_syn
        kappa_dat = kappa_syn
        bulk_dat  = bulk_syn
        rho_dat   = rho_syn
        alpha_dat = alpha_syn
     endif

     ! write velocity structure for data and synthetics to file -- LOCAL LEVEL
     file_structure_dat = 'structure_dat.dat'
     file_structure_syn = 'structure_syn.dat'
     open(unit=18,file=trim(out_dir2)//file_structure_dat,status='unknown')
     open(unit=19,file=trim(out_dir2)//file_structure_syn,status='unknown')
     do ispec = 1,NSPEC
        do j = 1,NGLLZ
           do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              ! TARGET MODEL (data)
              write(18,'(6e20.10)') x_plot(iglob), z_plot(iglob), &
                   kappa_dat(i,j,ispec), mu_dat(i,j,ispec), rho_dat(i,j,ispec), &
                   log( beta_dat(i,j,ispec) / beta0 )

              ! CURRENT MODEL (synthetics)
              write(19,'(6e20.10)') x_plot(iglob), z_plot(iglob), &
                   kappa_syn(i,j,ispec), mu_syn(i,j,ispec), rho_syn(i,j,ispec), &
                   log( beta_syn(i,j,ispec) / beta0 )
           enddo
        enddo
     enddo
     close(18)
     close(19)

     !stop 'TESTING'

!!$     ! TARGET MODEL (data)
!!$     call local2global(alpha_dat, temp_global1)
!!$     call local2global(beta_dat, temp_global2)
!!$     open(unit=18,file=trim(out_dir2)//file_structure_dat,status='unknown')
!!$     do iglob = 1,NGLOB
!!$        write(18,'(6e20.10)') x_plot(iglob), z_plot(iglob), temp_global1(iglob), temp_global2(iglob), &
!!$            log( temp_global1(iglob) / alpha0 )*100.0, log( temp_global2(iglob) / beta0 )*100.0
!!$     enddo
!!$     close(18)
!!$
!!$     ! CURRENT MODEL (synthetics)
!!$     call local2global(alpha_syn, temp_global1)     ! P-wave structure
!!$     call local2global(beta_syn, temp_global2)      ! S-wave structure
!!$     open(unit=19,file=trim(out_dir2)//file_structure_syn,status='unknown')
!!$     do iglob = 1,NGLOB
!!$        write(19,'(6e20.10)') x_plot(iglob), z_plot(iglob), temp_global1(iglob), temp_global2(iglob), &
!!$            log( temp_global1(iglob) / alpha0 )*100.0, log( temp_global2(iglob) / beta0 )*100.0
!!$     enddo
!!$     close(19)

     ! write reference structure values to file
     open(unit=16,file=trim(out_dir2)//'reference_values0.dat', status='unknown')
     write(16,'(6e20.10)') PWAVESPEED, SWAVESPEED, DENSITY, BWAVESPEED, INCOMPRESSIBILITY, RIGIDITY
     close(16)

     ! write reference structure values to file
     open(unit=16,file=trim(out_dir2)//'reference_values.dat', status='unknown')
     write(16,'(6e20.10)') alpha0, beta0, rho0, bulk0, kappa0, mu0
     close(16)

     ! write period to file
     open(unit=16,file=trim(out_dir2)//'reference_period.dat', status='unknown')
     write(16,*) 2*hdur
     close(16)

!!$  ! plot phase velocity map for synthetics (and data)
!!$  iopt = 1 + idat
!!$  filename1 = 'get_model.csh'
!!$  filename2 = trim(script_dir)//'plot_model.pl'
!!$  open(19,file=filename1,status='unknown')
!!$  write(19,'(4a,i5,3a,2f12.6,2a)') trim(filename2),' ', trim(out_dir2),' ', &
!!$     iopt, ' ', trim(file_structure_syn), ' ', sngl(beta0), sngl(t_target),' ',trim(file_structure_dat)
!!$  close(19)
!!$  call system('chmod 755 get_model.csh ; get_model.csh')

     ! min/max velocities and gridpoints per wavelength (locally defined)
     alpha_min = minval(alpha_dat)
     alpha_max = maxval(alpha_dat)
     beta_min  = minval(beta_dat)
     beta_max  = maxval(beta_dat)

     print *, ' Mean spacing between GLL points in m :'
     print *, '   dh    :', sngl(dh)
     print *, ' Velocities in km/s :'
     print *, '   Pmin  :', sngl(alpha_min/1000.0)
     print *, '   Pmax  :', sngl(alpha_max/1000.0)
     print *, '   Smin  :', sngl(beta_min/1000.0)
     print *, '   Smax  :', sngl(beta_max/1000.0)
     print *, ' Gridpoints per wavelength :'
     print *, '   Pmin  :', floor(hdur*alpha_min/dh)
     print *, '   Pmax  :', floor(hdur*alpha_max/dh)
     print *, '   Smin  :', floor(hdur*beta_min/dh)
     print *, '   Smax  :', floor(hdur*beta_max/dh)

    ! estimate the time step
    dh = HEIGHT/dble((NGLLZ-1)*NEZ)
    if (dh > LENGTH/dble((NGLLX-1)*NEX)) dh = LENGTH/dble((NGLLX-1)*NEX)
    print *
    print *,'                       space step (km) :', sngl(dh/1000.0)
    if (ISURFACE == 0) then
       print *,'time step est from courant = 0.2, Pmax : ',sngl(0.2*dh/alpha_max),' seconds'
       print *,'time step est from courant = 0.2, Pmin : ',sngl(0.2*dh/alpha_min),' seconds'
       print *,'time step est from courant = 0.2, Smax : ',sngl(0.2*dh/beta_max),' seconds'
       print *,'time step est from courant = 0.2, Smin : ',sngl(0.2*dh/beta_min),' seconds'
       print *,'                      actual time step : ',sngl(DT),' seconds'
    endif

     !--------------------------------------
     ! Prior to Feb 2006, we simply assigned the source location to the nearst gridpoint.
     ! But now we compute the locations such that ANY location is allowed.
     !--------------------------------------
     ! events for synthetics (eglob_syn)

     if (ISURFACE == 1) then
        ! fill a portion of the [1:MAX_EVENT] vector with the SoCal events
        x_eve_lon0_syn(1:nevent) = socal_events_lon(1:nevent)
        z_eve_lat0_syn(1:nevent) = socal_events_lat(1:nevent)

        ! all vectors should be this length
        nevent = MAX_EVENT

        ! convert from lat-lon to mesh coordinates (in meters)
        call mesh_geo(nevent,x_eve_lon0_syn,z_eve_lat0_syn,x_eve0_syn,z_eve0_syn,UTM_PROJECTION_ZONE,ILONLAT2MESH)
     endif

     ! filter target points (utm-mesh) (update nevent)
     print *
     print *, 'nevent: ', nevent
     print *, 'events for synthetics:'
     call station_filter(nevent, x_eve0_syn, z_eve0_syn, ifilter_eve_syn, SOURCE_GRID_BUFFER)
     print *, 'nevent: ', nevent

     if (nevent < 1) stop 'Must have at least one event'

     ! allocate variables
     if (ISURFACE == 1) allocate(x_eve_lon_syn(nevent),z_eve_lat_syn(nevent))
     allocate(x_eve_syn(nevent),z_eve_syn(nevent))
     allocate(eglob_syn(nevent))
     allocate(ispec_eve_syn(nevent),xi_eve_syn(nevent),gamma_eve_syn(nevent))
     allocate(otime_syn(nevent))

     ! origin time for synthetics
     ! initialize to mprior
     otime_syn(1:nevent) = otime(1:nevent)  ! otime_syn is length nevent, otime is length MAX_EVENT
     if (PERT_SOURCE_T == 1) then
        otime_syn(:) = otime_syn(:) + otime_eve_syn_pert(:)
     endif

     ! assign initial synthetic events to (x_eve_syn, z_eve_syn)
     ! initialize to mprior
     do i = 1,nevent
        x_eve_syn(i) = x_eve0_syn(ifilter_eve_syn(i))
        z_eve_syn(i) = z_eve0_syn(ifilter_eve_syn(i))
     enddo
     if (PERT_SOURCE_X == 1) then
        x_eve_syn(:) = x_eve_syn(:) + x_eve_syn_pert(:)
        z_eve_syn(:) = z_eve_syn(:) + z_eve_syn_pert(:)
     endif

     !write(*,'(2e24.16)') (x_eve_syn(i), z_eve_syn(i), i = 1,nevent)

     ! determine the (eta, xi) corresponding to the (initial synthetic) points
     ! this UPDATES x_eve, z_eve; eglob_syn is the index of the closest gridpoint
     call locate_targets(nevent, x_eve_syn, z_eve_syn, eglob_syn, ispec_eve_syn, xi_eve_syn, gamma_eve_syn)

     !write(*,'(2e24.16)') (x_eve_syn(i), z_eve_syn(i), i = 1,nevent)

     print *
     print *, 'ievent, x, z, iglob, ispec, xi, gamma'
     do i = 1,nevent
        write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_eve_syn(i), z_eve_syn(i), &
             eglob_syn(i), ispec_eve_syn(i), xi_eve_syn(i), gamma_eve_syn(i)
     enddo

     ! display target events and final events -- and the distance between (in meters)
     if (ISURFACE == 1) then
        ! convert from mesh to lat-lon (gets x_eve_lon_syn and z_eve_lat_syn)
        call mesh_geo(nevent, x_eve_lon_syn, z_eve_lat_syn, x_eve_syn, z_eve_syn, UTM_PROJECTION_ZONE, IMESH2LONLAT)

        ! The distance error is due to the UTM conversion, but the lat-lon points are only
        ! used for plotting purposes, so this is fine.
        print *
        print *, 'events [x_eve_lon0_syn, z_eve_lat0_syn, x_eve_lon_syn, x_eve_lat_syn, dist (m)]:'
        do i = 1,nevent
           temp1 = x_eve_lon0_syn(ifilter_eve_syn(i))
           temp2 = z_eve_lat0_syn(ifilter_eve_syn(i))
           temp3 = x_eve_lon_syn(i)
           temp4 = z_eve_lat_syn(i)
           temp5 = acos( sin(temp2/DEG)*sin(temp4/DEG) + cos(temp2/DEG)*cos(temp4/DEG)*cos(temp1/DEG - temp3/DEG) )
           write(*,'(i8, 5f17.10)') i, temp1, temp2, temp3, temp4, temp5*6371.0*1000.0
        enddo

     else if (ISURFACE == 0) then
        print *
        print *, 'events, km [x_eve0_syn, z_eve0_syn, x_eve_syn, x_eve_syn, dist]:'
        do i = 1,nevent
           temp1 = x_eve0_syn(ifilter_eve_syn(i))
           temp2 = z_eve0_syn(ifilter_eve_syn(i))
           temp3 = x_eve_syn(i)
           temp4 = z_eve_syn(i)
           temp5 = sqrt( (temp3-temp1)**2 + (temp4 - temp2)**2 )
           write(*,'(i8, 5f17.10)') i,temp1/1000.0,temp2/1000.0,temp3/1000.0,temp4/1000.0,temp5/1000.0
        enddo
     endif

     ! write events for PRIOR MODEL to file
     if (.not. READ_IN ) then
        open(19,file=trim(out_dir2)//'events_prior_xz.dat',status='unknown')
        do ievent = 1,nevent
           write(19,'(3f20.10,1i12)') x_eve0_syn(ievent), z_eve0_syn(ievent), otime(ievent), ievent
        enddo
        close(19)
        if (ISURFACE == 1) then
           open(18,file=trim(out_dir2)//'events_prior_lonlat.dat',status='unknown')
           do ievent = 1,nevent
              write(18,'(2f20.10,1i12)') x_eve_lon0_syn(ievent), z_eve_lat0_syn(ievent), ievent
           enddo
           close(18)
        endif
     endif

     ! write events for SYNTHETICS to file
     if (.not. READ_IN ) then
        open(19,file=trim(out_dir2)//'events_syn_xz.dat',status='unknown')
        do ievent = 1,nevent
           write(19,'(3f20.10,1i12)') x_eve_syn(ievent), z_eve_syn(ievent), otime_syn(ievent), ievent
        enddo
        close(19)
        if (ISURFACE == 1) then
           open(18,file=trim(out_dir2)//'events_syn_lonlat.dat',status='unknown')
           do ievent = 1,nevent
              !write(18,*) sngl(x_lon(eglob_syn(ievent))), sngl(z_lat(eglob_syn(ievent))), ievent
              write(18,'(2f20.10,1i12)') x_eve_lon_syn(ievent), z_eve_lat_syn(ievent), ievent
           enddo
           close(18)
        endif
     endif

     !write(*,'(2e24.16)') (x_eve_lon_syn(i), z_eve_lat_syn(i), i = 1,nevent)

     !stop 'testing event locations'

     !--------------------------------------
     ! events for data (eglob_dat)

     ! allocate variables
     if (ISURFACE == 1) allocate(x_eve_lon_dat(nevent),z_eve_lat_dat(nevent))
     allocate(x_eve_dat(nevent),z_eve_dat(nevent))
     allocate(eglob_dat(nevent))
     allocate(ispec_eve_dat(nevent),xi_eve_dat(nevent),gamma_eve_dat(nevent))
     allocate(otime_dat(nevent))

     if ( PERT_SOURCE_T == 1 .or. PERT_SOURCE_X == 1 ) then  ! source perturbations

        if (ISURFACE == 1) then   ! read in PERTURBED events from another file (for GJI-2007 simulations)
           ! initialize to no perturbations from synthetics
           !otime_dat(:) = otime_syn(:)
           !x_eve_dat(:) = x_eve_syn(:)
           !z_eve_dat(:) = z_eve_syn(:)

           ! initialize to prior model
           ! NOTE: these arrays have different lengths
           otime_dat(1:nevent) = otime(1:nevent)
           x_eve_dat(1:nevent) = x_eve0_syn(1:nevent)
           z_eve_dat(1:nevent) = z_eve0_syn(1:nevent)

           if (PERT_SOURCE_T == 1) then
              otime_dat(:) = otime_dat(:) + otime_eve_dat_pert(:)
           endif

           if (PERT_SOURCE_X == 1) then
              x_eve_dat(:) = x_eve_dat(:) + x_eve_dat_pert(:)
              z_eve_dat(:) = z_eve_dat(:) + z_eve_dat_pert(:)
           endif

!!$           ! perturbations generated from gji_src_inversion.m
!!$           ! NOTE: OLD ordering scheme for source indexing
!!$           open(19,file='INPUT/events_xyt_pert.dat',status='unknown')
!!$           do ievent = 1,25
!!$              read(19,*) temp1,temp2,temp3
!!$              if ( PERT_SOURCE_T == 0 ) temp3 = 0.0
!!$              if ( PERT_SOURCE_X == 0 ) then
!!$                 temp1 = 0.0  ; temp2 = 0.0
!!$              endif
!!$              x_eve_dat(ievent) = x_eve_syn(ievent) + temp1
!!$              z_eve_dat(ievent) = z_eve_syn(ievent) + temp2
!!$              otime_dat(ievent) = otime_syn(ievent) + temp3
!!$           enddo
!!$           close(19)

!!$           open(19,file='OUTPUT/run_2550/events_xy.dat',status='unknown')
!!$           do ievent = 1,25
!!$              read(19,'(3f20.10,1i12)') x_eve_dat(ievent), z_eve_dat(ievent), otime_dat(ievent), itemp1
!!$           enddo
!!$           close(19)
!!$          otime_dat(5) = otime_syn(5) + 0.8   ! testing for event 5 (0.8)

        else if (ISURFACE == 0) then

           ! perturb the events to get target events for the synthetics
           ! perturbation is described as (r,theta) polar coordinates
           amin = 0.0
           amax = 2.*PI
           rmin = 0.0
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
              !pangle = 30.0*(PI/180.0)               ! fix for TESTING

              ! radial perturbation (in meters)
              call random_number(rand1)
              prad = (rmax - rmin)*rand1 + rmin
              !prad = rmax

              ! fill vectors
              otime_dat(ievent) = otime_syn(ievent) - ptime                ! NOTE THE SIGN
              x_eve_dat(ievent) = x_eve_syn(ievent) + cos(pangle)*prad
              z_eve_dat(ievent) = z_eve_syn(ievent) + sin(pangle)*prad
              temp1 = sqrt((x_eve_syn(ievent)-x_eve_dat(ievent))**2 + (z_eve_syn(ievent)-z_eve_dat(ievent))**2)

              !write(*,'(5f18.8)') x_eve_dat(ievent)/1000.0, z_eve_dat(ievent)/1000.0, &
              !     x_eve_syn(ievent)/1000.0, z_eve_syn(ievent)/1000.0, temp1/1000.0
           enddo

        endif

        ! determine the (eta, xi) corresponding to the target points
        ! this UPDATES x_eve, z_eve; eglob is the index of the closest gridpoint
        call locate_targets(nevent, x_eve_dat, z_eve_dat, eglob_dat, ispec_eve_dat, xi_eve_dat, gamma_eve_dat)

        print *
        print *, 'ievent, x, z, iglob, ispec, xi, gamma'
        do i = 1,nevent
           write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_eve_dat(i), z_eve_dat(i), &
                eglob_dat(i), ispec_eve_dat(i), xi_eve_dat(i), gamma_eve_dat(i)
        enddo

        if (ISURFACE == 1) then
           ! convert from mesh to lat-lon
           call mesh_geo(nevent, x_eve_lon_dat, z_eve_lat_dat, x_eve_dat, z_eve_dat, UTM_PROJECTION_ZONE, IMESH2LONLAT)

           ! display target events and perturbed events -- and the distance between (in meters)
           print *
           print *, 'events [x_eve_lon_syn, z_eve_lat_syn, x_eve_lon_dat, z_eve_lat_dat, dist (m)]:'
           do i = 1,nevent
              temp1 = x_eve_lon_syn(i)
              temp2 = z_eve_lat_syn(i)
              temp3 = x_eve_lon_dat(i)
              temp4 = z_eve_lat_dat(i)
              temp5 = acos( sin(temp2/DEG)*sin(temp4/DEG) + cos(temp2/DEG)*cos(temp4/DEG)*cos(temp1/DEG - temp3/DEG) )
              write(*,'(i8, 5f17.10)') i, temp1, temp2, temp3, temp4, temp5*6371.0*1000.0
           enddo
        endif

     else  ! no source perturbation for target sources

        if (ISURFACE == 1) then
           x_eve_lon_dat(:) = x_eve_lon_syn(:)
           z_eve_lat_dat(:) = z_eve_lat_syn(:)
        endif
        eglob_dat(:)     = eglob_syn(:)
        x_eve_dat(:)     = x_eve_syn(:)
        z_eve_dat(:)     = z_eve_syn(:)
        ispec_eve_dat(:) = ispec_eve_syn(:)
        xi_eve_dat(:)    = xi_eve_syn(:)
        gamma_eve_dat(:) = gamma_eve_syn(:)
        otime_dat(:)     = otime_syn(:)

     endif   ! PERT_SOURCE

     ! write events for data to file
     if (.not. READ_IN ) then
        open(19,file=trim(out_dir2)//'events_dat_xz.dat',status='unknown')
        do ievent = 1,nevent
           write(19,'(3f20.10,1i12)') x_eve_dat(ievent), z_eve_dat(ievent), otime_dat(ievent), ievent
        enddo
        close(19)
        if (ISURFACE == 1) then
           open(19,file=trim(out_dir2)//'events_dat_lonlat.dat',status='unknown')
           do ievent = 1,nevent
              write(19,'(2f20.10,1i12)') x_eve_lon_dat(ievent), z_eve_lat_dat(ievent), ievent
           enddo
           close(19)
        endif
     endif

     !stop 'testing'

     !--------------------------------------
     ! receivers

     ! allocation and initialization
     !allocate(rglob_temp(MAX_SR))
     !allocate(x_rec0(MAX_SR),z_rec0(MAX_SR),x_rec_lon0(MAX_SR),z_rec_lat0(MAX_SR),ifilter_rec(MAX_SR))
     x_rec0(:) = 0.0      ; z_rec0(:) = 0.0
     x_rec_lon0(:) = 0.0  ; z_rec_lat0(:) = 0.0
     !rglob_temp(:) = 0

     ! get the lat-lon of the TARGET RECEIVERS

     if (ISURFACE == 1) then   ! enter target lat-lon

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
           read(88,*) (x_rec_lon0(i),z_rec_lat0(i),i = 1,nrec)
           close(88)

        else if (IREC_SPACE == 4) then ! 'regular' mesh of receivers

           ! calculate mesh spacing
           dx = LENGTH/NMESH_REC
           dz = dx
           j = 0
           do xmesh = STATION_GRID_BUFFER+1.0,LENGTH,dx
              do zmesh = STATION_GRID_BUFFER+1.0,HEIGHT,dz
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

     else if (ISURFACE == 0) then   ! enter target x-z

        if (IREC_SPACE == 1) then  ! individual receivers
           stop ' no individual receiver listed'

        else if (IREC_SPACE == 2) then  ! actual station locations

           ! read in (target) receivers
           recfile = trim(in_dir)//'STATIONS_socal1D_rand15'  ! 200 km width
           !recfile = trim(in_dir)//'STATIONS_socal1D_rand30'   ! 400 km width
           open(88,file=trim(recfile),status='unknown')
           read(88,*) nrec
           read(88,*) (x_rec0(i),z_rec0(i),i = 1,nrec)
           close(88)

        else if (IREC_SPACE == 4) then  ! 'regular' line of receivers (BODY waves)

           ! calculate mesh spacing
           dx = LENGTH/NMESH_REC
           j = 0
           do xmesh = dx,LENGTH-dx,dx
              j = j+1
              x_rec0(j) = xmesh
              z_rec0(j) = HEIGHT
           enddo
           nrec = j

           ! receivers at depth (GJI2005 paper)
           !x_rec0(j+1) = 150.0d3; z_rec0(j+1) = 0.75*HEIGHT
           !x_rec0(j+2) = 150.0d3; z_rec0(j+2) = 0.50*HEIGHT
           !x_rec0(j+3) = 150.0d3; z_rec0(j+3) = 0.25*HEIGHT
           !nrec = nrec+3

        else
           stop 'invalid entry for IREC_SPACE'

        endif  ! IREC_SPACE

     endif

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

     if (ISURFACE == 1) then
        ! convert from lat-lon to mesh coordinates (in meters)
        call mesh_geo(nrec,x_rec_lon0,z_rec_lat0,x_rec0,z_rec0,UTM_PROJECTION_ZONE,ILONLAT2MESH)
     endif

     ! filter target points (utm-mesh) -- update nrec
     call station_filter(nrec, x_rec0, z_rec0, ifilter_rec, STATION_GRID_BUFFER)
     !call station_filter(nrec, x_rec, z_rec, dmin_trsh_r, ncoast, coast_x, coast_z)
     !call station_filter_2(nrec, x_rec, z_rec, -1)  ! -1 for left, +1 for right

     if (nrec < 1) stop 'Must have at least one receiver'

     ! allocate vectors
     if (ISURFACE == 1) allocate(x_rec_lon(nrec),z_rec_lat(nrec))
     allocate(x_rec(nrec),z_rec(nrec))
     allocate(rglob(nrec))
     allocate(ispec_rec(nrec),xi_rec(nrec),gamma_rec(nrec))

     ! assign allowable target receivers to (x_rec, z_rec)
     do i = 1,nrec
        x_rec(i) = x_rec0(ifilter_rec(i))
        z_rec(i) = z_rec0(ifilter_rec(i))
     enddo

     ! determine the (eta, xi) corresponding to the target points
     ! this UPDATES x_rec, z_rec; rglob is the index of the closest gridpoint
     call locate_targets(nrec, x_rec, z_rec, rglob, ispec_rec, xi_rec, gamma_rec)

     print *
     print *, 'irec, x, z, iglob, ispec, xi, gamma'
     do i = 1,nrec
        write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_rec(i), z_rec(i), rglob(i), ispec_rec(i), xi_rec(i), gamma_rec(i)
     enddo

     if (ISURFACE == 1) then
        ! convert from mesh to lat-lon
        call mesh_geo(nrec, x_rec_lon, z_rec_lat, x_rec, z_rec, UTM_PROJECTION_ZONE, IMESH2LONLAT)

        ! display target receivers and final receivers -- and the distance between (in meters)
        ! The distance error is due to the UTM conversion, but the lat-lon points are only
        ! used for plotting purposes, so this is fine.
        print *
        print *, ' receivers [x_rec_lon0, z_rec_lat0, x_rec_lon, x_rec_lat, dist (m)] :'
        do i = 1,nrec
           temp1 = x_rec_lon0(ifilter_rec(i))
           temp2 = z_rec_lat0(ifilter_rec(i))
           temp3 = x_rec_lon(i)
           temp4 = z_rec_lat(i)
           temp5 = acos( sin(temp2/DEG)*sin(temp4/DEG) + cos(temp2/DEG)*cos(temp4/DEG)*cos(temp1/DEG - temp3/DEG) )
           write(*,'(i8, 5f17.10)') i, temp1, temp2, temp3, temp4, temp5*6371.0*1000.0
        enddo
     endif

     ! receiver indices for writing the misfit chi_data(rec), summed over all events
     rglob0(1:nrec) = rglob(1:nrec)

     !do i = 1,nrec
     !   print *, sngl(x(rglob(i))/1000.0),sngl(z(rglob(i))/1000.0),sngl(x_rec(i)/1000.0),sngl(z_rec(i)/1000.0)
     !enddo

     !stop 'testing'

     ! deallocate
     !deallocate(x_rec,z_rec)

     ! write receivers to file
     if (ISURFACE == 1) then
        open(18,file=trim(out_dir2)//'recs_lonlat.dat',status='unknown')
        do irec = 1,nrec
           write(18,'(2f20.10,1i12)') x_rec_lon(irec), z_rec_lat(irec), irec
        enddo
        close(18)
     endif
     open(19,file=trim(out_dir2)//'recs_xz.dat',status='unknown')
     do irec = 1,nrec
        write(19,'(2f20.10,1i12)') x_rec(irec), z_rec(irec), irec
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
        call lagrange_poly(xi_eve_syn(ievent),NGLLX,xigll,hxie,hpxie)
        call lagrange_poly(gamma_eve_syn(ievent),NGLLZ,zigll,hgammae,hpgammae)
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
     ! FOR THE OCEAN MICROSEISM SIMULATIONS.

     if (ISURFACE == 1) then

        ! mesh of target points
        x_recf0(:) = 0.0
        z_recf0(:) = 0.0
        dx = LENGTH/9.0
        dz = dx
        i = 0
        do xmesh = 0.0,LENGTH,dx
           do zmesh = 0.0,HEIGHT,dz
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
        do i = 1,nrecf
           x_recf(i) = x_recf0(ifilter_recf(i))
           z_recf(i) = z_recf0(ifilter_recf(i))
           !print *, i, x_recf(i), z_recf(i)
        enddo

        ! fglob is the index of the CLOSEST GRIDPOINT
        call locate_targets(nrecf, x_recf, z_recf, fglob)

        !do i = 1,nrecf
        !  print *, sngl(x(fglob(i))/1000.0),sngl(z(fglob(i))/1000.0),sngl(x_recf(i)/1000.0),sngl(z_recf(i)/1000.0)
        !enddo

        if (0 == 1) then
           ! display target receivers and final fake receivers -- distances in km
           print *
           print *, ' fake receivers [x_rec0, z_rec0, x_rec, x_rec, dist (km)]:'
           do i = 1,nrecf
              temp1 = x_recf0(ifilter_recf(i)) / 1000.0
              temp2 = z_recf0(ifilter_recf(i)) / 1000.0
              temp3 = x_recf(i) / 1000.0
              temp4 = z_recf(i) / 1000.0
              temp5 = sqrt( (temp3-temp1)**2 + (temp4-temp2)**2 )
              write(*,'(i8,5f17.10)') i, temp1, temp2, temp3, temp4, temp5
           enddo
        endif

        ! deallocate: all we need is fglob
        deallocate(x_recf,z_recf)

     endif

     !--------------------------------------
     ! initial model vectors

     !nmod_str = NVAR_STRUCT * NGLOB      ! structure parameters
     nmod_str = NVAR_STRUCT * NLOCAL     ! structure parameters
     nmod_src = NVAR_SOURCE * nevent     ! source parameters
     nmod     = nmod_str + nmod_src       ! total model parameters

     print *, '         NGLOB ', NGLOB
     print *, '        NLOCAL ', NLOCAL
     print *, '          nmod ', nmod
     print *, '      nmod_str ', nmod_str
     print *, '      nmod_src ', nmod_src

     m_inds(:,:) = 0
     m_inds(1,1) = 1                   ; m_inds(1,2) = NLOCAL
     m_inds(2,1) = nmod_str+1          ; m_inds(2,2) = nmod_str+nevent
     m_inds(3,1) = nmod_str+nevent+1   ; m_inds(3,2) = nmod_str+2*nevent
     m_inds(4,1) = nmod_str+2*nevent+1 ; m_inds(4,2) = nmod
!!$     m_inds(1,1) = 1                   ; m_inds(1,2) = NLOCAL
!!$     m_inds(2,1) = NLOCAL+1            ; m_inds(2,2) = nmod_str
!!$     m_inds(3,1) = nmod_str+1          ; m_inds(3,2) = nmod_str+nevent
!!$     m_inds(4,1) = nmod_str+nevent+1   ; m_inds(4,2) = nmod_str+2*nevent
!!$     m_inds(5,1) = nmod_str+2*nevent+1 ; m_inds(5,2) = nmod

     ! write source indexing to file
     open(unit=18,file=trim(out_dir2)//'m_inds.dat',status='unknown')
     do i = 1,NVAR
         write(18,'(2i10)') m_inds(i,1),m_inds(i,2)
     enddo
     close(18)

     ! source indexing for the model vector
     allocate(index_source(NVAR_SOURCE, nevent))
     itemp = 0
     do i = 1,NVAR_SOURCE
        do j = 1,nevent
           itemp = itemp + 1
           index_source(i,j) = itemp
        enddo
     enddo

     ! write source indexing to file
     open(unit=18,file=trim(out_dir2)//'index_source.dat',status='unknown')
     do i = 1,NVAR_SOURCE
        do j = 1,nevent
           write(18,'(4i10)') index_source(i,j),i,j
        enddo
     enddo
     close(18)

     ! allocate model vector
     allocate(m0(nmod),mt(nmod),mdiff(nmod),mprior(nmod),mtarget(nmod))
     allocate(cov_imetric(nmod),cov_model(nmod))
     !allocate(icov_metric(nmod),icov_model(nmod))
     allocate(gradient(nmod),gradient_data(nmod),gradient_model(nmod))
     allocate(m0_vec(nmod),mt_vec(nmod))
     allocate(g0(nmod),gt(nmod),gk(nmod),p0(nmod),pk(nmod))
     !allocate(m0_vec_initial(nmod))

     ! gradient values for each parameter
     allocate(gradient_data_all(nevent,NVAR))
     gradient_data_all(:,:) = 0.0

     !allocate(m_vel(nmod_str),m0_vel(nmod_str),mt_vel(nmod_str))

     ! initialize model vectors
     m0(:)       = 0.0
     mt(:)       = 0.0
     m0_vec(:)   = 0.0
     mt_vec(:)   = 0.0
     mprior(:)   = 0.0
     mtarget(:)  = 0.0

     ! initialize gradient vectors
     g0(:) = 0.0
     gk(:) = 0.0
     gt(:) = 0.0
     p0(:) = 0.0
     pk(:) = 0.0

     !------------------
     ! source-related vectors

     ! scaling for source model coefficients
     allocate(source_gradient(nmod_src))
     allocate(m_src_syn(nmod_src),m_src_dat(nmod_src),m_src_prior(nmod_src))
     !allocate(m_src_syn_vec(nmod_src),m_src_dat_vec(nmod_src))
     !allocate(m_src_syn_vec_initial(nmod_src),m_src_syn_initial(nmod_src))
     !allocate(m_scale_src_all(nmod_src))
     !allocate(m_src(nmod_src),mt_src(nmod_src),m0_src(nmod_src),m_scale_src_all(nmod_src))

     m_src_prior(:) = 0.0
     m_src_dat(:) = 0.0
     m_src_syn(:) = 0.0
     !m_src_dat_vec(:) = 0.0
     !m_src_syn_vec(:) = 0.0
     !m_src_syn_vec_initial(:) = 0.0

     ! fill source model vector for SYNTHETICS AND DATA
     ! NOTE: Only fill the source parameters that are being RUN;
     !       otherwise, the computation of the norm of the model vector
     !       and target model vector will differ based on events
     !       that you are not interested in.
     !       This also avoids the issue of iterating the non-running
     !       sources to new sources.
     do ievent = ievent_min,ievent_max
        !itemp1 = (ievent-1)*3 + 1
        !itemp2 = (ievent-1)*3 + 2
        !itemp3 = (ievent-1)*3 + 3

        itemp1 = index_source(1,ievent)
        itemp2 = index_source(2,ievent)
        itemp3 = index_source(3,ievent)

        m_src_syn(itemp1) = otime_syn(ievent)
        m_src_syn(itemp2) = x_eve_syn(ievent)
        m_src_syn(itemp3) = z_eve_syn(ievent)

        m_src_dat(itemp1) = otime_dat(ievent)
        m_src_dat(itemp2) = x_eve_dat(ievent)
        m_src_dat(itemp3) = z_eve_dat(ievent)

        m_src_prior(itemp1) = otime(ievent)
        m_src_prior(itemp2) = x_eve0_syn(ievent)
        m_src_prior(itemp3) = z_eve0_syn(ievent)
     enddo

     !------------------
     ! PRIOR and INITIAL model vectors

     ! initial source parameters for SYNTHETICS
     !m_src_syn_vec_initial(:) = m_src_syn_vec(:)
     !m_src_syn_initial(:)     = m_src_syn_vec(:) - m_src_syn_vec_initial(:)     ! = 0.

     ! initial source parameters for DATA and SYNTHETICS
     ! perturbations are w.r.t. INITIAL SYNTHETICS
     !m_src_dat(:) = m_src_dat_vec(:) - m_src_syn_vec_initial(:)
     !m_src_syn(:) = m_src_syn_vec(:) - m_src_syn_vec_initial(:)

     ! create m0 -- should be identical to what is "un-done" in the CG algorithm
     temp_local1(:,:,:) = log( beta_syn(:,:,:) / beta0 )
     call local2mvec(temp_local1, nmod_src, m_src_syn, nmod, m0)

     ! create m0_vec
     call local2mvec(beta_syn, nmod_src, m_src_syn, nmod, m0_vec)
     !m0_vec_initial(:) = m0_vec(:)

     ! create mprior
     temp_local1(:,:,:) = 0.0
     call local2mvec(temp_local1, nmod_src, m_src_prior, nmod, mprior)

     ! create mtarget
     temp_local1(:,:,:) = log( beta_dat(:,:,:) / beta0 )
     call local2mvec(temp_local1, nmod_src, m_src_dat, nmod, mtarget)

     ! write out prior, initial, and target models
     if (.not. READ_IN ) then
        open(89,file=trim(out_dir2)//'prior_initial_target_models.dat',status='unknown')
        do i = 1,nmod
           write(89,'(3e24.12)') mprior(i), m0(i), mtarget(i)
        enddo
        close(89)

!!$        open(88,file=trim(out_dir2)//'initial_source.dat',status='unknown')
!!$        do i = 1,nmod_src
!!$           !write(88,'(5e24.12)') m_scale_src_all(i), m_src_syn_vec_initial(i), m_src_syn_initial(i), m_src_dat_vec(i), m_src_dat(i)
!!$           !write(88,'(4e24.12)') m_src_syn_vec_initial(i), m_src_syn_initial(i), m_src_dat_vec(i), m_src_dat(i)
!!$           !write(88,'(4e24.12)') m_src_syn_vec(i), m_src_syn(i), m_src_dat_vec(i), m_src_dat(i)
!!$           write(88,'(3e24.12)') m_src_prior(i), m_src_syn(i), m_src_dat(i)
!!$        enddo
!!$        close(88)
     endif

     !------------------
     ! MODEL COVARIANCE MATRIX = inverse metric tensor
     ! NOTE: This part of the code is crucial for doing joint source-structure inversions,
     !       which require a particular balance between source gradients and structure gradients.
     !       IT IS STILL IN DEVELOPMENT, but will work for the test cases.

     ! (INVERSE) METRIC TENSOR: variances only
     ! NOTE: The inverse metric tensor is a diagonal covariance matrix.

     ! scaling vector for STRUCTURE parameters
     sigma_checker_scale = 2.0   ! to mimic a Gaussian distibution
     if ((IMODEL_SYN == 3) .and. (IMODEL_DAT == 3)) then
        m_scale_str(1) = sigma_beta0
        !m_scale_str(2) = sigma_bulk0
     else
        m_scale_str(1) = (afac/sigma_checker_scale)/100.0    ! dimensionless : beta* = ln(beta/beta0)
        !m_scale_str(2) = (afac/sigma_checker_scale)/100.0    ! dimensionless : bulk* = ln(bulk/bulk0)
     endif

     ! scaling vector for SOURCE parameters
     if ( GJI_PAPER == 1 ) then
        !joint_str_src_grad_ratio = (70000.0 / 5000.0)   ! to match the new version
        !joint_str_src_grad_ratio = 1.0
        m_scale_src(1) = 2*hdur         ! dts is scaled by period (20s)
        m_scale_src(2) = 2*hdur*beta0   ! dxs is scaled by wavelength (70000 m)
        m_scale_src(3) = 2*hdur*beta0   ! dzs is scaled by wavelength (70000 m)
     else
        !joint_str_src_grad_ratio = 1.0
        !m_scale_src(1) = src_pert_time * (20.0/1.0)/(70000.0/5000.0)   ! to match the old version
        m_scale_src(1) = src_pert_time
        m_scale_src(2) = src_pert_dist
        m_scale_src(3) = src_pert_dist
     endif

     sigma_beta = m_scale_str(1)
     !sigma_bulk = m_scale_str(2)
     sigma_ts   = m_scale_src(1)
     sigma_xs   = m_scale_src(2)
     sigma_zs   = m_scale_src(3)

!!$     ! For JOINT inversions, there are several ways to balance the source and structure portions.
!!$     ! A simple way is to incorporate a scaling factor in the sigma-structure terms.
!!$     ! This will control the relative weight of each set of parameters.
!!$     ! NOTE : We approximate our checkerboard-generated structure values by a Gaussian distribution.
!!$     if ( INV_SOURCE == 1 .and. INV_STRUCT_BETA == 1) then
!!$        if ( READ_IN ) then
!!$           joint_str_src_grad_ratio = 1.0
!!$        else
!!$           ! see also scale_struct_gji (F^2)
!!$           joint_str_src_grad_ratio = 10.0
!!$           !joint_str_src_grad_ratio = 6.25
!!$        endif
!!$        joint_total_weight = 0.5         ! split the source and structure parts evenly
!!$
!!$     else   ! basic structure or basic source inversion
!!$        joint_str_src_grad_ratio = 1.0
!!$        joint_total_weight = 1.0
!!$
!!$     endif
!!$     jfac    = 0.5 * (joint_str_src_grad_ratio + 1.0)  ! =1 for non-joint inversions
!!$     jsw_str = joint_str_src_grad_ratio / jfac             ! =1 for non-joint inversions
!!$     jsw_src = 1.0 / jfac                                ! =1 for non-joint inversions

     !------------------------------

     ! initialize to no weights
     covm_weight_parts(:) = 0.0

     ! factors for balancing the model norm term
     ! NOTE: we ignore parts of the model norm that do not participate in the inversion
     if ( INV_STRUCT_BETA == 1 .and. INV_SOURCE_T == 0 .and. INV_SOURCE_X == 0 ) then   ! structure
        fac_str = 1.00
     else if ( INV_STRUCT_BETA == 0 .and. INV_SOURCE_T == 1 .and. INV_SOURCE_X == 0 ) then   ! origin time
        fac_ts = 1.00
     else if ( INV_STRUCT_BETA == 0 .and. INV_SOURCE_T == 0 .and. INV_SOURCE_X == 1 ) then   ! location
        fac_xs = 0.50 ;  fac_ys = 0.50
     else if ( INV_STRUCT_BETA == 1 .and. INV_SOURCE_T == 1 .and. INV_SOURCE_X == 0 ) then   ! structure + origin time
        fac_str = 0.50 ;  fac_ts = 0.50
     else if ( INV_STRUCT_BETA == 1 .and. INV_SOURCE_T == 0 .and. INV_SOURCE_X == 1 ) then   ! structure + location
        fac_str = 0.50 ;  fac_xs = 0.25 ;  fac_ys = 0.25
     else if ( INV_STRUCT_BETA == 0 .and. INV_SOURCE_T == 1 .and. INV_SOURCE_X == 1 ) then   ! source
        fac_ts = 0.50 ;  fac_xs = 0.25 ;  fac_ys = 0.25
     else if ( INV_STRUCT_BETA == 1 .and. INV_SOURCE_T == 1 .and. INV_SOURCE_X == 1 ) then   ! joint
        fac_str = 1.0/2.0  ;  fac_ts = 1.0/6.0 ;  fac_xs = 1.0/6.0  ;  fac_ys = 1.0/6.0
        !fac_str = 0.85  ;  fac_ts = 0.05 ;  fac_xs = 0.05  ;  fac_ys = 0.05
     else
        stop 'you must invert for something'
     endif

     ! NOTE: normalize to ensure that the total norm is 1
     covm_weight_parts(1) = fac_str
     covm_weight_parts(2) = fac_ts
     covm_weight_parts(3) = fac_xs
     covm_weight_parts(4) = fac_ys
     covm_weight_parts(:) = covm_weight_parts(:) / sum(covm_weight_parts)
     !covm_weight_parts(:) = 1.0   ! OLD VERSION

     ! simple way to balance GRADIENT according to norm (NOT norm-squared)
     covm_weight_parts = covm_weight_parts * covm_weight_parts

     ! unbalanced initial gradient values (gradient_norm_all_unbalanced.dat)
     ! ugsq --> unbalanced gradient norm-squared
     covg_weight_parts(:) = 1.0

     !fac_total_g = 1.0
     if ( INV_STRUCT_BETA == 1 .and. INV_SOURCE_T == 1 .and. INV_SOURCE_X == 1 ) then   ! joint

        ! NOTE: we obtain these values from either:
        !         gradient_norm_data_all_unbalanced.dat
        !         gradient_norm_all_unbalanced.dat

        ugsq_str = 1.0 ; ugsq_ts = 1.0 ; ugsq_xs = 1.0 ; ugsq_ys = 1.0

        ! TEST CASES: 1-25 events, 1-5 events, 5-5 event
        !ugsq_str = 0.1296496361d6 ; ugsq_ts = 0.2706999550d4; ugsq_xs = 0.1334372694d4 ; ugsq_ys = 0.1743549898d4  ! Nfac = 2, 25 events

        !ugsq_str = 0.2070999127d6 ; ugsq_ts = 0.3507740673d4 ; ugsq_xs = 0.1940368586d4 ; ugsq_ys = 0.2142118702d4  ! Nfac = 3, 25 events
        !ugsq_str = 0.4247330856d6 ; ugsq_ts = 0.5896527525d4 ; ugsq_xs = 0.2643538184d4 ; ugsq_ys = 0.3342815391d4  ! Nfac = 3, 5 events
        !ugsq_str = 0.1205146534d6 ; ugsq_ts = 0.2926409454d3 ; ugsq_xs = 0.9695606936d3 ; ugsq_ys = 0.7224889563d3  ! Nfac = 3, 1 event

        ! 9200
        ugsq_str = 0.8418200277d6 ; ugsq_ts = 0.4639935107d3 ; ugsq_xs = 0.1635906965d4 ; ugsq_ys = 0.1147402541d4   ! Gaussians, 1 event

        ! 1200 - with new source
        !ugsq_str = 0.1176265806d6 ; ugsq_ts = 0.1823448342d3 ; ugsq_xs = 0.6421163226d2 ; ugsq_ys = 0.2776564145d1   ! Gaussians, 1 event
        ! 1200 - balance full gradient
        !ugsq_str = 0.1176904168d6 ; ugsq_ts = 0.1443574694d3 ; ugsq_xs = 0.8206282043d2 ; ugsq_ys = 0.1632463053d2   ! Gaussians, 1 event

        ! 1300 - with 0.7,0.1,0.1,0.1, old source, MOISMPRIOR
        !ugsq_str = 0.1191920326d6 ; ugsq_ts = 0.4533353832d3 ; ugsq_xs = 0.6191223030d2 ; ugsq_ys = 0.4387672586d2   ! Gaussians, 1 event

        ! 3200,3600
        !ugsq_str = 0.3026988450d6 ; ugsq_ts = 0.8260031407d3 ; ugsq_xs = 0.1187469038d4 ; ugsq_ys = 0.1381078706d4   ! Gaussians, 25 events
        ! 3400,3800
        !ugsq_str = 0.3028636452d6 ; ugsq_ts = 0.8327738622d3  ; ugsq_xs =0.1195006162d4 ; ugsq_ys = 0.1381645351d4   ! Gaussians, 25 events
        ! 3300/1700
        !ugsq_str = 0.3352965529d6 ; ugsq_ts = 0.8068082292d3 ; ugsq_xs = 0.1209790150d4 ; ugsq_ys = 0.1391331978d4  ! Gaussians, 25 events

!!$        ! ad hoc: choose balance among the four parts of the gradient
!!$        fac_str = 1.0/2.0  ;  fac_ts = 1.0/6.0 ;  fac_xs = 1.0/6.0  ;  fac_ys = 1.0/6.0
!!$        !fac_str = 1.0/4.0  ;  fac_ts = 1.0/4.0 ;  fac_xs = 1.0/4.0  ;  fac_ys = 1.0/4.0
!!$        !fac_str = 0.7  ;  fac_ts = 0.1 ;  fac_xs = 0.1  ;  fac_ys = 0.1
!!$        !fac_str = 0.85  ;  fac_ts = 0.05 ;  fac_xs = 0.05  ;  fac_ys = 0.05
!!$
!!$        covg_weight_parts(1) = ugsq_str/fac_str
!!$        covg_weight_parts(2) = ugsq_ts/fac_ts
!!$        covg_weight_parts(3) = ugsq_xs/fac_xs
!!$        covg_weight_parts(4) = ugsq_ys/fac_ys
!!$        covg_weight_parts(:) = covg_weight_parts(:) / sum(covg_weight_parts)

        covg_weight_parts(1) = ugsq_str / covm_weight_parts(1)
        covg_weight_parts(2) = ugsq_ts / covm_weight_parts(2)
        covg_weight_parts(3) = ugsq_xs / covm_weight_parts(3)
        covg_weight_parts(4) = ugsq_ys / covm_weight_parts(4)

        ! COMMENT OUT if you want the gradient-norm parts to sum to ONE
        covg_weight_parts(:) = covg_weight_parts(:) / sum(covg_weight_parts)

     endif

     ! TESTING
     !covm_weight_parts(:) = 1.0
     !covg_weight_parts(:) = 1.0

     ! Re-set any weights that are 0.0 to 1.0; these portions of the covariance matrix
     ! should not play a role, since the corresponding gradients will always be 0.0.
     !if (fac_str < EPS) fac_str = 1.0
     !if (fac_ts < EPS) fac_ts = 1.0
     !if (fac_xs < EPS) fac_xs = 1.0
     !if (fac_ys < EPS) fac_ys = 1.0

     ! Because we do not have perfect coverage, we need to adjust the normalization such that
     ! a perfectly recovered (source or structure) model gives a norm of about 1.0.
     ! The means adjusting the weights of the respective parts, based on the
     ! perfectly recovered model (i.e., no data errors added, no model norm term).
     ! Thus, the norm of the target model will then be somewhat GREATER than 1.0.
     !coverage_str = 0.666 / 0.962
     !coverage_src = 0.946 / 1.018

     ! If the initial and target models are from a Gaussian distribution,
     ! then this factor is not needed.
     coverage_str = 1.0
     coverage_src = 1.0

     cov_model(m_inds(1,1):m_inds(1,2)) = ( sigma_beta )**2 / da_local_vec(:) * AREA * coverage_str
     cov_model(m_inds(2,1):m_inds(2,2)) = sigma_ts**2 * dnevent_run * coverage_src
     cov_model(m_inds(3,1):m_inds(3,2)) = sigma_xs**2 * dnevent_run * coverage_src
     cov_model(m_inds(4,1):m_inds(4,2)) = sigma_zs**2 * dnevent_run * coverage_src

     ! incorporate relative weighting to make the final metric
     ! STRUCTURE: (fac_str / ugsq_str) * ugsq_str --> fac_str
     cov_imetric(m_inds(1,1):m_inds(1,2)) = cov_model(m_inds(1,1):m_inds(1,2)) / covg_weight_parts(1)
     cov_imetric(m_inds(2,1):m_inds(2,2)) = cov_model(m_inds(2,1):m_inds(2,2)) / covg_weight_parts(2)
     cov_imetric(m_inds(3,1):m_inds(3,2)) = cov_model(m_inds(3,1):m_inds(3,2)) / covg_weight_parts(3)
     cov_imetric(m_inds(4,1):m_inds(4,2)) = cov_model(m_inds(4,1):m_inds(4,2)) / covg_weight_parts(4)

     ! METRIC TENSOR: inverse diagponal covariance matrix
     !icov_metric = 1.0 / cov_imetric
     !icov_model  = 1.0 / cov_model       ! unbalanced version

     ! possible stopping criteria based on the target model
     ! NOTE 1: this depends on what you pick as your initial model (prior mean, or a sample)
     ! NOTE 2: THIS IS NOT USED
     !chi_model_stop = 0.5 * model_target_norm
     chi_model_stop = 0.0

     ! possible stopping criteria based on fitting the data
     chi_data_stop = 0.0
     if (ADD_DATA_ERRORS) chi_data_stop = 0.5

     ! write model covariance values to file
     open(unit=19,file=trim(out_dir2)//'sigma_values.dat',status='unknown')
     write(19,'(2e20.10)') sigma_beta, m_scale_str(1)
     !write(19,'(2e20.10)') sigma_bulk, m_scale_str(2)
     write(19,'(2e20.10)') sigma_ts, m_scale_src(1)
     write(19,'(2e20.10)') sigma_xs, m_scale_src(2)
     write(19,'(2e20.10)') sigma_zs, m_scale_src(3)
     close(19)

     open(unit=19,file=trim(out_dir2)//'scaling_values_covm.dat',status='unknown')
     do i=1,NVAR
        write(19,*) 'covm_weight_parts', covm_weight_parts(i)
     enddo
     write(19,*) 'sum_covm_weight_parts', sum(covm_weight_parts)
     do i=1,NVAR
        write(19,*) 'covg_weight_parts', covg_weight_parts(i)
     enddo
     write(19,*) 'sum_covg_weight_parts', sum(covg_weight_parts)

!!$     write(19,*) 'fac_str', fac_str
!!$     write(19,*) 'fac_ts', fac_ts
!!$     write(19,*) 'fac_xs', fac_xs
!!$     write(19,*) 'fac_ys', fac_ys
!!$     write(19,*) 'fac_total', fac_total
     write(19,*) 'ugsq_str', ugsq_str
     write(19,*) 'ugsq_ts', ugsq_ts
     write(19,*) 'ugsq_xs', ugsq_xs
     write(19,*) 'ugsq_ys', ugsq_ys
     write(19,*) 'dnevent_run', dnevent_run
     write(19,*) 'coverage_str', coverage_str
     write(19,*) 'coverage_src', coverage_src
!!$     write(19,*) 'cov_imetric_fac_str', (fac_str / ugsq_str) * fac_total
!!$     write(19,*) 'cov_imetric_fac_ts', (fac_ts / ugsq_ts) * fac_total
!!$     write(19,*) 'cov_imetric_fac_xs', (fac_xs / ugsq_xs) * fac_total
!!$     write(19,*) 'cov_imetric_fac_ys', (fac_ys / ugsq_ys) * fac_total
     close(19)

     open(unit=19,file=trim(out_dir2)//'scaling_values.dat',status='unknown')
     write(19,*) 'GJI_PAPER = ', GJI_PAPER
     write(19,*) 'IRUNZ = ', IRUNZ
     write(19,'(6i10)') 1, NLOCAL, NLOCAL+1, nmod_str, nmod_str+1, nmod
     write(19,'(2i10)') nevent_run, NLOCAL
     write(19,'(5e14.6)') sum(da_local(:,:,:)), sum(da_local_vec(:)), sum(da_global(:)), LENGTH*HEIGHT, AREA
     write(19,'(3e14.6)') minval(da_local_vec(:)), maxval(da_local_vec(:)), sum(da_local_vec(:))/NLOCAL
     write(19,'(3e14.6)') minval(da_global(:)), maxval(da_global(:)), sum(da_global(:))/NGLOB
!!$     write(19,*)
!!$     write(19,*) 'COVERAGE WEIGHTS: '
!!$     write(19,'(2f14.6)') coverage_str, coverage_src
!!$     write(19,*) 'JOINT WEIGHTS (based on unbalanced gradients): '
!!$     write(19,'(4f14.6)') fac_str, fac_ts, fac_xs, fac_ys
!!$     write(19,'(1f14.6)') fac_total
!!$     write(19,'(4e14.6)') ugsq_str, ugsq_ts, ugsq_xs, ugsq_ys
!!$     !write(19,'(1e14.6)') ugsq_src
!!$     write(19,'(4e14.6)') fac_str/ugsq_str, fac_ts/ugsq_ts, fac_xs/ugsq_xs, fac_ys/ugsq_ys
!!$     write(19,'(4e14.6)') (fac_str/ugsq_str)*fac_total, (fac_ts/ugsq_ts)*fac_total, &
!!$                           (fac_xs/ugsq_xs)*fac_total,  (fac_ys/ugsq_ys)*fac_total
!!$     !write(19,'(a20,5e14.6)') 'JOINT WEIGHTS: ', joint_str_src_grad_ratio, joint_total_weight, jfac, jsw_str, jsw_src
!!$     !write(19,'(a20,3e14.6)') 'STRUCTURE: ', jsw_str, joint_total_weight, jsw_str/joint_total_weight
!!$     !write(19,'(a20,3e14.6)') 'SOURCE: ', jsw_src, joint_total_weight, jsw_src/joint_total_weight
     close(19)

     open(unit=19,file=trim(out_dir2)//'cov_model_diagonal.dat',status='unknown')
     do i = 1,nmod
        write(19,'(2e20.10)') cov_model(i), cov_imetric(i)
        !write(19,'(4e20.10)') cov_imetric(i), icov_metric(i), cov_model(i), icov_model(i)
     enddo
     close(19)

     ! write model covariance matrix diagonal to file
     !open(unit=19,file=trim(out_dir2)//'cov_imetric_diagonal.dat',status='unknown')
     !write(19,'(1e20.10)') (cov_imetric(i), i = 1,nmod)
     !close(19)

     ! write inverse model covariance matrix diagonal to file
     !open(unit=19,file=trim(out_dir2)//'icov_metric_diagonal.dat',status='unknown')
     !write(19,'(1e20.10)') (icov_metric(i), i = 1,nmod)
     !close(19)

     !------------------
     ! DATA COVARIANCE MATRIX AND DATA INDEXING

     ! number of measurements
     nmeas     = nevent * nrec * NCOMP
     nmeas_run = nevent_run * nrec * NCOMP
     print *, ' nmeas = nevent * nrec * ncomp = ', nmeas
     print *, ' nmeas_run = nevent_run * nrec * ncomp = ', nmeas_run

     ! data indexing
     allocate(index_data(nevent,nrec,ncomp))
     itemp = 0
     do i = 1,nevent
        do j = 1,nrec
           do k = 1,NCOMP
              itemp = itemp + 1
              index_data(i,j,k) = itemp
           enddo
        enddo
     enddo

     ! write data indexing to file
     open(unit=18,file=trim(out_dir2)//'index_data.dat',status='unknown')
     do i = 1,nevent
        do j = 1,nrec
           do k = 1,NCOMP
              write(18,'(4i10)') index_data(i,j,k),i,j,k
           enddo
        enddo
     enddo
     close(18)

     ! data covariance
     allocate(cov_data(nmeas))
     cov_data(:) = 0.0

     if (IKER == 0) then
        cov_data(:) = SIGMA_WAVEFORM * SIGMA_WAVEFORM * nmeas_run
     else if (IKER == 1) then
        cov_data(:) = SIGMA_DT * SIGMA_DT * nmeas_run
     else if (IKER == 2) then
        cov_data(:) = SIGMA_DLNA * SIGMA_DLNA * nmeas_run
     endif

     if (IKER == 1) then
        open(unit=19,file=trim(out_dir2)//'scaling_values_covd.dat',status='unknown')
        write(19,*) 'ievent_min', ievent_min
        write(19,*) 'ievent_max', ievent_max
        write(19,*) 'nevent_run', nevent_run
        write(19,*) 'nrec', nrec
        write(19,*) 'NCOMP', NCOMP
        write(19,*) 'nmeas_run', nmeas_run
        write(19,*) 'SIGMA_DT', SIGMA_DT
        close(19)
     endif

     ! write data covariance matrix diagonal to file
     open(unit=19,file=trim(out_dir2)//'cov_data_diagonal.dat',status='unknown')
     write(19,'(1e20.10)') (cov_data(i), i = 1,nmeas)
     close(19)

     ! load measurement perturbations
     if ( ADD_DATA_ERRORS ) then

        if (nmeas >= 10000) stop 'measurement error file only contains 10000 lines'

        allocate(measure_pert_vec(nmeas))
        measure_pert_vec(:) = 0.0

        ! By reading in the same set of perturbations, rather than generating
        ! a new distrubution each time, we can be sure that the data errors will
        ! not influence the convergence.
        open(55,file=trim(in_dir)//'sigma_0p1_pert.dat',status='unknown')
        do i = 1,nmeas
           read(55,*) measure_pert_vec(i)
        enddo
        close(55)

        ! check this by writing to file
        open(55,file=trim(out_dir2)//'data_errors_added.dat',status='unknown')
        do i = 1,nmeas
           write(55,'(1e20.12)') measure_pert_vec(i)
        enddo
        close(55)

     endif

     !------------------

     !nrec = 1

     if (ISOURCE_LOG) open(91,file=trim(out_dir2)//'source_vector.log',status='unknown')

     var_red_val = 100.0    ! initialize
     chi_k_val   = 1.0d10     ! initialize to a very high value
     itest = 0
     !============================================
     ! LOOP 2: over the models in the CG optimization
     do istep = 0,2*NITERATION
        !============================================

        imod = (istep - mod(istep,2))/2          ! index into model number

        irun = irun0 + istep
        print *,' ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  == '
        print *,'  istep, imod, irun : ', istep, imod, irun
        if (INV_STRUCT_BETA == 1) print *, '  inverting for structure parameters'
        if (INV_SOURCE_T == 1) print *, '  inverting for source origin times'
        if (INV_SOURCE_X == 1) print *, '  inverting for source locations'
        print *,' ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  == '

        if (READ_IN) then
           out_dir1 = trim(out_dir2)     ! no iteration
        else
           write(out_dir1,'(a,i4.4,a)') trim(out_dir3)//'run_',irun,'/'
           command1 = 'mkdir ' // trim(out_dir1)
           call system(command1)
           !open(19,file='temp.csh',status='unknown')
           !write(19,'(2a)') 'mkdir ',out_dir1
           !close(19)
           !call system('chmod 755 temp.csh ; temp.csh')
        endif

        ! KEY: obtain the structure and source arrays from the model vector:
        !      structure parameters fill the top portion (beta_syn);
        !      source parameters fill the bottom portion (m_src_syn)
        if (itest == 0) then      ! reference model (current model)
           call mvec2local(nmod, nmod_src, m0_vec, beta_syn, m_src_syn)

        else if (itest == 1) then  ! test model
           call mvec2local(nmod, nmod_src, mt_vec, beta_syn, m_src_syn)
        endif

        ! for CG algorithm, update kappa_syn and mu_syn
        if (.not. READ_IN) then
          !kappa_syn = rho_syn * bulk_syn * bulk_syn
          mu_syn    = rho_syn * beta_syn * beta_syn
        endif

        ! read in the sources from another file
        ! NOTE: FOR NOW, WE DO NOT LOAD THE DATA SOURCES -- THEY SHOULD BE IDENTICAL.
        !       (In the future, we might want to modify this to read in ANY data sources.)
        if ( READ_IN .and. (INV_SOURCE_T == 1 .and. INV_SOURCE_X == 1)) then

           write(filename,'(a,i4.4,a)') trim(out_dir2)//'src_syn_m',hmod,'.dat'
           open(unit=18,file=filename,status='unknown')
           m_src_syn(:) = m_src_dat(:)   ! initialize to no perturbation from target sources
           do i = 1,nevent
              !itemp1 = (i-1)*3 + 1     ! origin time
              !itemp2 = (i-1)*3 + 2     ! x position
              !itemp3 = (i-1)*3 + 3     ! z position

              itemp1 = index_source(1,i)
              itemp2 = index_source(2,i)
              itemp3 = index_source(3,i)

              read(18,*) temp1, temp2, &
                   m_src_syn(itemp1), m_src_syn(itemp2), m_src_syn(itemp3), &
                   temp3, temp4, temp5
           enddo
           close(18)

           ! update the entries in the model vector, then re-compute the norm
           !m0(nmod_str+1 : nmod) = m_src_syn_vec(:) - m_src_syn_vec_initial(:)
           m0(nmod_str+1 : nmod) = m_src_syn(:)

        endif

        !--------------------------
        ! compute model norm term in the misfit function
        ! NOTE 1: THREE possible outputs: cov_model, cov_model with covm weights, cov_model with covg weights
        ! NOTE 2: THE FACTOR OF 0.5 IS NOT INCLUDED HERE
        ! NOTE 3: These variables are norm-squared (no square root is taken)

        imnorm = 1   ! =1 for m-mprior ; =0 for g

        ! model vector, cov_model
        filename = trim(out_dir1)//'model_norm_all.dat'
        if (itest == 0) then      ! reference model (current model)
           call compute_norm_sq(filename, imnorm, &
               ievent_min, ievent_max, nevent, index_source, nmod, &
               m0, mprior, cov_model, model_norm_parts, covm_weight_parts)
        else if (itest == 1) then  ! test model
           call compute_norm_sq(filename, imnorm, &
               ievent_min, ievent_max, nevent, index_source, nmod, &
               mt, mprior, cov_model, model_norm_parts, covm_weight_parts)
        endif
        !filename = trim(out_dir1)//'model_norm_all.dat'
        !call write_norm_sq(filename,model_norm_parts)

!!$        ! this output should exactly match gradient_norm_model_all.dat
!!$        call compute_norm_sq(filename, imnorm, &
!!$             ievent_min, ievent_max, nevent, index_source, nmod, &
!!$             m0, mprior, cov_imetric, model_norm_parts)

        ! model_norm is used in the misfit function
        model_norm_struct = model_norm_parts(1)
        model_norm_source = sum( model_norm_parts(2:4) )
        model_norm        = sum( model_norm_parts(1:4) )

        ! target model vector, cov_model
        ! NOTE: THIS SHOULD NOT CHANGE
        filename = trim(out_dir1)//'model_norm_target_all.dat'
        call compute_norm_sq(filename, imnorm, &
             ievent_min, ievent_max, nevent, index_source, nmod, &
             mtarget, mprior, cov_model, model_norm_target_parts, covm_weight_parts)
        !call write_norm_sq(filename,model_norm_target_parts)

        model_norm_target_struct = model_norm_target_parts(1)
        model_norm_target_source = sum( model_norm_target_parts(2:4) )
        model_norm_target        = sum( model_norm_target_parts(1:4) )

        ! compute the difference between the present model and the target model
        mdiff(:) = 0.0
        if (itest == 0) then
           mdiff = mtarget - m0 + mprior
        else if (itest == 1) then
           mdiff = mtarget - mt + mprior
        endif

        ! diff vector, cov_model
        filename = trim(out_dir1)//'model_norm_diff_all.dat'
        call compute_norm_sq(filename, imnorm, &
           ievent_min, ievent_max, nevent, index_source, nmod, &
           mdiff, mprior, cov_model, model_norm_diff_parts, covm_weight_parts)
        !call write_norm_sq(filename,model_norm_diff_parts)

        model_norm_diff_struct = model_norm_diff_parts(1)
        model_norm_diff_source = sum( model_norm_diff_parts(2:4) )
        model_norm_diff        = sum( model_norm_diff_parts(1:4) )

        !stop 'TESTING MODEL NORMS'

        !--------------------------

        ! initialize measurement vector for the current model (or test model)
        !allocate(measure_vec(nmeas,5))   ! waveform, tt-xcor, dlna-xcorr, tt-mtm, dlna-mtm
        allocate(measure_vec(nmeas,5))   ! tt-xcor-pert, tt-xcor, dlna-xcorr-pert, dln-xcorr, waveform
        measure_vec(:,:) = 0.0
        !imeasure = 0                     ! KEY: initialize for each new model

        ! initialize gradient vectors
        gradient(:) = 0.0 ; gradient_data(:) = 0.0 ; gradient_model(:) = 0.0
        source_gradient(:) = 0.0
        btype_kernel_sum = 0.0

        ! initialize misfit
        chi_data(:,:,:,:) = 0.0

        !stop 'testing'

        ! flag to write wave-speed values and source locations to file ONCE per model
        ifirst_event = 1

        !============================================
        ! LOOP 3: over the events
        do ievent = ievent_min, ievent_max
        !do ievent = 5,5    ! 1,5 or 5,5
           !============================================

           print *,'------------------------------------'
           print *,'  EVENT NUMBER ',ievent
           print *,'  istep, imod, irun : ', istep, imod, irun
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

           ! get the lat-lon of the TARGET SOURCES

           if (ISURFACE == 1) then

              if (ISRC_SPACE <= 5) then  ! if you do NOT want a point source from the event list

                 x_src_lon0(:) = 0.0  ; z_src_lat0(:) = 0.0
                 x_src0(:) = 0.0      ; z_src0(:) = 0.0

                 if (ISRC_SPACE == 1) then  ! point source(s)

                    !x_src0(1) = LENGTH/4      ; z_src0(1) = HEIGHT/2
                    !x_src_lon0(1) = -118.5370     ; z_src_lat0(1) = 34.2130   ! Northridge
                    !x_src_lon0(1) = -117.918     ; z_src_lat0(1) = 32.3290   ! 2004/06/15
                    !x_src_lon0(1) = -117.776     ; z_src_lat0(1) = 33.917    ! Yorba Linda

                    !x_src_lon0(1) = -119.0     ; z_src_lat0(1) = 33.0
                    !x_src_lon0(2) = -118.0     ; z_src_lat0(2) = 33.5
                    !x_src_lon0(3) = -119.5     ; z_src_lat0(3) = 34.3

                    x_src_lon0(1) = -117.7616    ; z_src_lat0(1) = 33.9228     ! event kernel for Yorba Linda
                    !x_src_lon0(1) = -116.8480    ; z_src_lat0(1) = 34.3098     ! event kernel for Big Bear

                    !x_src_lon0(1) = x_eve_lon_syn(ievent)
                    !z_src_lat0(1) = z_eve_lat_syn(ievent)
                    nsrc = 1

                 else if (ISRC_SPACE == 2) then  ! finite source segment

                    ! specify the target starting point of the fault, the azimuth, and the length
                    x_src_lon_i = -119.0    ;  z_src_lat_i = 33.0   ;   flen   = 100.0d+03  ! short fault
                    !x_src_lon_i = -118.0   ;  z_src_lat_i = 32.0   ;   flen   = 500.0d+03  ! long fault
                    src_az = -45.0*(PI/180.0)        ! azimuth (clockwise from north)
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
                       print *, x_src0(i)/1000.0, z_src0(i)/1000.0
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
                    dx = dh/2.0
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

                 do i = 1,nsrc
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
                 do i = 1,nsrc
                    x_src(i) = x_src0(ifilter_src(i))
                    z_src(i) = z_src0(ifilter_src(i))
                 enddo

                 ! determine the (eta, xi) corresponding to the target points
                 ! this UPDATES x_src, z_src; sglob is the index of the closest gridpoint
                 call locate_targets(nsrc, x_src, z_src, sglob, ispec_src, xi_src, gamma_src)

                 print *
                 do i = 1,nsrc
                    write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_src(i), z_src(i), sglob(i), ispec_src(i), xi_src(i), gamma_src(i)
                 enddo

                 ! convert from mesh to lat-lon
                 call mesh_geo(nsrc, x_src_lon, z_src_lat, x_src, z_src, UTM_PROJECTION_ZONE, IMESH2LONLAT)

                 ! display target sources and final sources -- and the distance between (in meters)
                 ! The distance error is due to the UTM conversion, but the lat-lon points are only
                 ! used for plotting purposes, so this is fine.
                 print *
                 print *, 'sources [x_src_lon0, z_src_lat0, x_src_lon, x_src_lat, dist (m)]:'
                 do i = 1,nsrc
                    temp1 = x_src_lon0(ifilter_src(i))
                    temp2 = z_src_lat0(ifilter_src(i))
                    temp3 = x_src_lon(i)
                    temp4 = z_src_lat(i)
                    temp5 = acos( sin(temp2/DEG)*sin(temp4/DEG) + cos(temp2/DEG)*cos(temp4/DEG)*cos(temp1/DEG - temp3/DEG) )
                    write(*,'(i8, 5f17.10)') i, temp1, temp2, temp3, temp4, temp5*6371.0*1000.0
                 enddo

                 ! CHT, 05-Dec-2006
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

              else         ! select POINT SOURCE for array of events
                 nsrc = 1
              endif

           else
              if (ISRC_SPACE /= 6) stop ' for body waves, we only select a point source from the specified events'
              nsrc = 1

           endif  ! ISURFACE == 1

           ! allocate vectors (sources for data)
           if (ISURFACE == 1) allocate(x_src_lon_dat(nsrc),z_src_lat_dat(nsrc))
           allocate(x_src_dat(nsrc),z_src_dat(nsrc))
           allocate(sglob_dat(nsrc))
           allocate(ispec_src_dat(nsrc),xi_src_dat(nsrc),gamma_src_dat(nsrc))

           ! allocate 1-D Lagrange interpolators and derivatives
           allocate(hxisd_store(nsrc,NGLLX),hgammasd_store(nsrc,NGLLZ))

           if (ISRC_SPACE <= 5) then  ! if you do NOT want a point source from the event list
              if (ISURFACE == 0) stop ' for body waves, we only select a point source from the specified events'

              ! use same source for data and synthetics
              sglob_dat(1:nsrc)     = sglob(1:nsrc)   ! closest gridpoint
              x_src_dat(1:nsrc)     = x_src(1:nsrc)
              z_src_dat(1:nsrc)     = z_src(1:nsrc)
              x_src_lon_dat(1:nsrc) = x_src_lon(1:nsrc)
              z_src_lat_dat(1:nsrc) = z_src_lat(1:nsrc)

              ispec_src_dat(1:nsrc) = ispec_src(1:nsrc)
              gamma_src_dat(1:nsrc) = gamma_src(1:nsrc)
              hxisd_store(1:nsrc,:) = hxis_store(1:nsrc,:)
              hgammasd_store(1:nsrc,:) = hgammas_store(1:nsrc,:)

              origin_time_dat = origin_time

           else                      ! select POINT SOURCE for array of events

              !-----------------
              ! sources for data (DO NOT allow for iterative perturbations)

              sglob_dat(1) = eglob_dat(ievent)   ! closest gridpoint
              x_src_dat(1) = x_eve_dat(ievent)
              z_src_dat(1) = z_eve_dat(ievent)
              origin_time_dat = otime_dat(ievent)
              if (ISURFACE == 1) then
                 x_src_lon_dat(1) = x_eve_lon_dat(ievent)
                 z_src_lat_dat(1) = z_eve_lat_dat(ievent)
              endif
              ispec_src_dat(1) = ispec_eve_dat(ievent)
              xi_src_dat(1)    = xi_eve_dat(ievent)
              gamma_src_dat(1) = gamma_eve_dat(ievent)
              hxisd_store(1,:) = hxied_store(ievent,:)
              hgammasd_store(1,:) = hgammaed_store(ievent,:)

              !-----------------
              ! sources for synthetics (allows for source perturbations)

              ! source perturbations for this event
              !itemp1 = (ievent-1)*3 + 1
              !itemp2 = (ievent-1)*3 + 2
              !itemp3 = (ievent-1)*3 + 3

              itemp1 = index_source(1,ievent)
              itemp2 = index_source(2,ievent)
              itemp3 = index_source(3,ievent)

              ! allocate vectors (sources for synthetics)
              if (ISURFACE == 1) allocate(x_src_lon(nsrc),z_src_lat(nsrc))
              allocate(x_src(nsrc),z_src(nsrc))
              allocate(sglob(nsrc))
              allocate(ispec_src(nsrc),xi_src(nsrc),gamma_src(nsrc))

              ! KEY: source values for current model and event ievent
              origin_time_syn = m_src_syn(itemp1)     ! new origin time
              x_src(1) = m_src_syn(itemp2)            ! new x position
              z_src(1) = m_src_syn(itemp3)            ! new z position

              ! filter target points (utm-mesh) -- update nsrc
              call station_filter(nsrc, x_src, z_src, ifilter_src, SOURCE_GRID_BUFFER)

              if (nsrc /= 1) stop 'Must be a single point source'

              ! determine the (eta, xi) corresponding to the target points
              ! this UPDATES x_src, z_src; sglob is the index of the closest gridpoint
              call locate_targets(nsrc, x_src, z_src, sglob, ispec_src, xi_src, gamma_src)

              print *
              do i = 1,nsrc
                 write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_src(i), z_src(i), sglob(i), ispec_src(i), xi_src(i), gamma_src(i)
                 write(*,'(i8,2e18.8,2i8,2e18.8)') i, x_src_dat(i), z_src_dat(i), &
                      sglob_dat(i), ispec_src_dat(i), xi_src_dat(i), gamma_src_dat(i)
              enddo

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
              if (ISURFACE == 1) then
                 ! convert from mesh to lat-lon
                 call mesh_geo(nsrc, x_src_lon, z_src_lat, x_src, z_src, UTM_PROJECTION_ZONE, IMESH2LONLAT)

                 print *
                 print *, 'sources [x_src_lon_dat, z_src_lat_dat, x_src_lon, z_src_lat, dist (m)]:'
                 do i = 1,nsrc
                    temp1 = x_src_lon_dat(i)
                    temp2 = z_src_lat_dat(i)
                    temp3 = x_src_lon(i)
                    temp4 = z_src_lat(i)
                    temp5 = acos( sin(temp2/DEG)*sin(temp4/DEG) + cos(temp2/DEG)*cos(temp4/DEG)*cos(temp1/DEG - temp3/DEG) )
                    write(*,'(i8, 5f17.10)') i, temp1, temp2, temp3, temp4, temp5*6371.0*1000.0
                 enddo

              else if (ISURFACE == 0) then
                 print *
                 print *, 'sources [x_src_dat, z_src_dat, x_src, z_src, dist (m)]:'
                 do i = 1,nsrc
                    temp1 = x_src_dat(i)
                    temp2 = z_src_dat(i)
                    temp3 = x_src(i)
                    temp4 = z_src(i)
                    temp5 = sqrt( (temp3-temp1)**2 + (temp4-temp2)**2 )
                    write(*,'(i8, 5e20.10)') i, temp1, temp2, temp3, temp4, temp5
                 enddo
              endif

           endif  ! ISRC_SPACE

           ! source log file
           if (ISOURCE_LOG .and. ISRC_SPACE == 6) then

              ! indices into source vector
              itemp1 = index_source(1,ievent)
              itemp2 = index_source(2,ievent)
              itemp3 = index_source(3,ievent)

              !itemp1 = (ievent-1)*3 + 1     ! origin time
              !itemp2 = (ievent-1)*3 + 2     ! x position
              !itemp3 = (ievent-1)*3 + 3     ! z position

              write(91,*)
              write(91,*) '------------------------'
              write(91,*) 'istep, imod, ievent, irun : '
              write(91,*) istep, imod, ievent, irun
              write(91,*) 'SOURCE MODEL (xs, ys, t0) :'
              write(91,'(a12,3a18)') ' ','m_src_syn','m_src_dat','dat - syn'
              write(91,'(a12,3f18.8)') ' t0, s : ', &
                   m_src_syn(itemp1), &
                   m_src_dat(itemp1), &
                   m_src_dat(itemp1) - m_src_syn(itemp1)
              write(91,'(a12,3f18.8)') ' xs, km : ', &
                   m_src_syn(itemp2)/1000.0, &
                   m_src_dat(itemp2)/1000.0, &
                  (m_src_dat(itemp2) - m_src_syn(itemp2))/1000.0
              write(91,'(a12,3f18.8)') ' zs, km : ', &
                   m_src_syn(itemp3)/1000.0, &
                   m_src_dat(itemp3)/1000.0, &
                  (m_src_dat(itemp3) - m_src_syn(itemp3))/1000.0
           endif

           !--------------------------------------
           ! source time function FOR DATA AND SYNTHETICS

           ! source magnitude (same for data and synthetics)
           if (ISURFACE == 0) then
              ! DEBUG ARRAY SIZE
              f0(1) = FNORM * FOR_X
              f0(2) = FNORM * FOR_Y
              f0(3) = FNORM * FOR_Z
           else if (ISURFACE == 1) then
              f0(1) = FNORM
           else
              stop 'NCOMP must be 1 or 3'
           endif

           ! source time function for DATA (DO NOT allow for iterative origin time perturbations)
           stf_dat(:) = 0.0
           call get_source_time_function(origin_time_dat,stf_dat,ti)

           ! source function for data (includes the source magnitude f0)
           allocate(samp_dat(NSTEP,NCOMP,nsrc))
           samp_dat = 0.0
           do i = 1,nsrc
              do icomp = 1,NCOMP
                 samp_dat(:, icomp, i) = stf_dat(:) * f0(icomp)
              enddo
           enddo

           ! source time function for SYNTHETICS (allows for iterative origin time perturbations)
           stf_syn(:) = 0.0
           call get_source_time_function(origin_time_syn,stf_syn,ti)

           ! source function for synthetics (includes the source magnitude f0)
           allocate(samp(NSTEP,NCOMP,nsrc))
           samp = 0.0
           do i = 1,nsrc
              do icomp = 1,NCOMP
                 samp(:, icomp, i) = stf_syn(:) * f0(icomp)
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
                 write(*,*) '  actual location       : ', sngl(x_src(i)), ', ', sngl(z_src(i))
                 write(*,*) '  closest GLL gridpoint : ', sngl(x(sglob(i))), ', ', sngl(z(sglob(i)))
                 !write(*,*) '  actual location       : ', sngl(x_src_lon(i)), ', ', sngl(z_src_lat(i))
                 !write(*,*) '  closest GLL gridpoint : ', sngl(x_lon(sglob(i))), ', ', sngl(z_lat(sglob(i)))
                 do itime=1,NSTEP
                    write(12,'(1f16.6,1e16.6)') ti(itime), samp_dat(itime,icomp,i)
                 enddo
                 close(12)
              enddo
           enddo

           !--------------------------------------
           ! testing FFT routines

!!$           print *, 'testing the fft subroutine using the source time function'
!!$           call write_spectral_map(samp, nsrc, sglob, trim(out_dir)//'test')
!!$
!!$           print *, 'testing the bandpass filter subroutine using the source time function'
!!$           call filter(ti, samp, nsrc)
!!$
!!$           print *, 'testing the fft subroutine using the bandpass-filtered source time function'
!!$           call write_spectral_map(samp, nsrc, sglob, trim(out_dir)//'test_f')
!!$
!!$           stop

           !--------------------------------------

           ! write source and receivers to file
           ! NOTE THAT THE RECEIVER LABELING IS SIMPLY NUMERICAL ORDER (FOR NOW)
           !write(filename,'(a,i3.3,a)') trim(out_dir)//'sr_e',ievent,'.txt'
           file_src_rec = trim(out_dir)//'sr.txt'
           open(12,file=file_src_rec,status='unknown',iostat=ios)
           if (ios /= 0) stop 'Error opening out_dir/sr.txt'
           write(12,'(a,2e20.10,i10)') ('S ', x_plot(sglob(i)), z_plot(sglob(i)), i, i = 1,nsrc)
           write(12,'(a,2e20.10,i10)') ('R ', x_plot(rglob(i)), z_plot(rglob(i)), i, i = 1,nrec)
           close(12)

           ! the following are not needed, since we have sglob/rglob/fglob, the index vectors
           ! into the sources, receivers, and fake receivers

           deallocate(x_src,z_src)
           deallocate(x_src_dat,z_src_dat)
           if (ISURFACE == 1)  deallocate(x_src_lon,z_src_lat,x_src_lon_dat,z_src_lat_dat)

           !stop 'testing'

           !-----------------------------------------------------
           ! write the current models to file

           ! only write the files for the first event, since they are the same for each event
           if (ifirst_event == 1) then

              ! write velocity structure for data and synthetics to file -- LOCAL LEVEL
              file_structure_dat = 'structure_dat.dat'
              file_structure_syn = 'structure_syn.dat'
              open(unit=18,file=trim(out_dir1)//file_structure_dat,status='unknown')
              open(unit=19,file=trim(out_dir1)//file_structure_syn,status='unknown')
              do ispec = 1,NSPEC
                 do j = 1,NGLLZ
                    do i = 1,NGLLX
                       iglob = ibool(i,j,ispec)

                       ! TARGET MODEL (data)
                       write(18,'(6e20.10)') x_plot(iglob), z_plot(iglob), &
                            kappa_dat(i,j,ispec), mu_dat(i,j,ispec), rho_dat(i,j,ispec), &
                            log( beta_dat(i,j,ispec) / beta0 )

                       ! CURRENT MODEL (synthetics)
                       write(19,'(6e20.10)') x_plot(iglob), z_plot(iglob), &
                            kappa_syn(i,j,ispec), mu_syn(i,j,ispec), rho_syn(i,j,ispec), &
                            log( beta_syn(i,j,ispec) / beta0 )
                    enddo
                 enddo
              enddo
              close(18)
              close(19)

!!$              ! write wave-speed maps to file (data is always the same)
!!$              ! TARGET MODEL (data)
!!$              call local2global(bulk_dat, temp_global1)
!!$              call local2global(beta_dat, temp_global2)
!!$              open(unit=18,file=trim(out_dir1)//file_structure_dat,status='unknown')
!!$              do iglob = 1,NGLOB
!!$                 write(18,'(6e20.10)') x_plot(iglob), z_plot(iglob), temp_global1(iglob), temp_global2(iglob), &
!!$                      log( temp_global1(iglob) / bulk0 )*100.0, log( temp_global2(iglob) / beta0 )*100.0
!!$              enddo
!!$              close(18)
!!$
!!$              ! CURRENT MODEL (synthetics)
!!$              call local2global(bulk_syn, temp_global1)
!!$              call local2global(beta_syn, temp_global2)
!!$              open(unit=19,file=trim(out_dir1)//file_structure_syn,status='unknown')
!!$              do iglob = 1,NGLOB
!!$                 write(19,'(6e20.10)') x_plot(iglob), z_plot(iglob), temp_global1(iglob), temp_global2(iglob), &
!!$                      log( temp_global1(iglob) / bulk0 )*100.0, log( temp_global2(iglob) / beta0 )*100.0
!!$              enddo
!!$              close(19)

              ! write ALL source vectors to file (data is always the same)
              open(unit=19,file=trim(out_dir1)//'src_dat.dat',status='unknown')
              open(unit=20,file=trim(out_dir1)//'src_syn.dat',status='unknown')
              do i = 1,nevent
                 ! indices into source vector
                 itemp1 = index_source(1,i)
                 itemp2 = index_source(2,i)
                 itemp3 = index_source(3,i)

                 !itemp1 = (i-1)*3 + 1     ! origin time
                 !itemp2 = (i-1)*3 + 2     ! x position
                 !itemp3 = (i-1)*3 + 3     ! z position

                 ! longitude and latitude values from INITIAL SYNTHETIC events for plotting only
                 xtemp = x_plot(eglob_syn(i))
                 ztemp = z_plot(eglob_syn(i))

                 ! sources for data
                 write(19,'(8e20.10)') xtemp, ztemp, &
                      m_src_dat(itemp1), m_src_dat(itemp2), m_src_dat(itemp3), &
                      m_src_dat(itemp1) - m_src_dat(itemp1), &
                      m_src_dat(itemp2) - m_src_dat(itemp2), &
                      m_src_dat(itemp3) - m_src_dat(itemp3)

                 ! sources for synthetics
                 write(20,'(8e20.10)') xtemp, ztemp, &
                      m_src_syn(itemp1), m_src_syn(itemp2), m_src_syn(itemp3), &
                      m_src_dat(itemp1) - m_src_syn(itemp1), &
                      m_src_dat(itemp2) - m_src_syn(itemp2), &
                      m_src_dat(itemp3) - m_src_syn(itemp3)
              enddo
              close(19)
              close(20)

              ifirst_event = 0

           endif  ! ifirst_event

           ! source time functions for data and synthetics
           if (WRITE_STF_F) call write_seismogram(samp_dat, nsrc, trim(out_dir)//'stffor_dat')
           if (WRITE_STF_F) call write_seismogram(samp, nsrc, trim(out_dir)//'stffor_syn')

           ! STOPPING CRITERIA : exit program if structure parameters are unrealistic
           ! for either the test model or present model
!!$           if ( maxval( m0*m0 / (SIGMA_FAC**2 * cov_imetric) ) > 1. ) then
!!$              print *, ' maxval( m0*m0 / (SIGMA_FAC**2 * cov_imetric) ) > 1. ) : ', &
!!$                         maxval( m0*m0 / (SIGMA_FAC**2 * cov_imetric) )
!!$              print *, ' index : ', maxloc( m0*m0 / (SIGMA_FAC**2 * cov_imetric) )
!!$              stop ' STOP: model values are too extreme'
!!$           endif
!!$           if ( maxval( mt*mt / (SIGMA_FAC**2 * cov_imetric) ) > 1. ) then
!!$              print *, ' maxval( mt*mt / (SIGMA_FAC**2 * cov_imetric) ) > 1. ) : ', &
!!$                         maxval( mt*mt / (SIGMA_FAC**2 * cov_imetric) )
!!$              print *, ' index : ', maxloc( mt*mt / (SIGMA_FAC**2 * cov_imetric) )
!!$              stop ' STOP: model values are too extreme'
!!$           endif

!!$           vel_min = beta_min / 2.
!!$           vel_max = alpha_max * 2.
!!$           if (itest==1) then
!!$              if ( minval(mt_vec(1:nmod_str)) < vel_min .or. maxval(mt_vec(1:nmod_str)) > vel_max ) then
!!$                 print *, minval(mt_vec(1:nmod_str)), maxval(mt_vec(1:nmod_str))
!!$                 stop ' STOP: test model for structure is too extreme'
!!$              endif
!!$           else
!!$              if ( minval(m0_vec(1:nmod_str)) < vel_min .or. maxval(m0_vec(1:nmod_str)) > vel_max ) then
!!$                 print *, minval(m0_vec(1:nmod_str)), maxval(m0_vec(1:nmod_str))
!!$                 stop ' STOP: current model for structure is too extreme'
!!$              endif
!!$           endif

!stop 'testing -- not running wave simulations'

           !=============================================================================================
           ! ******************************** WAVE SIMULATIONS ******************************************
           !=============================================================================================

           !=========================
           ! DATA (forward wavefield)

           ! compute data for misfit kernels
           allocate(data(NSTEP,NCOMP,nrec))
           data = 0.0

           ! assign structure values for computing the data (always the same)
           kappa = kappa_dat
           mu    = mu_dat
           rho   = rho_dat

           ! DATA (forward wavefield)
           ! nsrc, sglob, samp_dat  : #src, src index, source
           ! nrec, rglob, data :      #rec, rec index, seis
           ! NOTE THAT sglob_dat ALLOWS FOR A PERTURBED SOURCE FROM THE SOURCES USED TO COMPUTE THE SYNTHETICS
           isolver = 1
           idata = 1
           call solver(isolver, idata, &
                nsrc, sglob_dat, ispec_src_dat,  hxisd_store, hgammasd_store, samp_dat, &
                nrec, rglob,     ispec_rec,      hxir_store,  hgammar_store,  data)

           ! write out seismograms at the receivers
           data_tag = 'dat'
           if (WRITE_SEISMO_F) call write_seismogram(data, nrec, trim(out_dir)//data_tag)

           !stop 'testing simulation for data'

           !=========================
           ! SYNTHETICS (forward wavefield)

           !stop 'testing'

           ! forward wavefield
           allocate(syn(NSTEP,NCOMP,nrec))
           syn = 0.0

           ! assign structure values for the synthetics
           ! (bulk and beta may change for each new model)
           !kappa = rho_syn * (alpha_syn * alpha_syn - FOUR_THIRDS * beta_syn * beta_syn )
           kappa = kappa_syn
           mu    = mu_syn
           rho   = rho_syn

           isolver = 1
           idata = 0
           last_frame_name = trim(out_dir)//'last_frame.txt'
           allocate(absorb_field(NSTEP, NCOMP, NGLL, NELE, NABSORB))
           call solver(isolver, idata, &
                nsrc, sglob, ispec_src, hxis_store, hgammas_store, samp, &
                nrec, rglob, ispec_rec, hxir_store, hgammar_store,  syn, trim(last_frame_name), absorb_field)

           ! write out seismograms at the receivers
           syn_tag = 'syn'
           !syn_tag = 'forward'
           if (WRITE_SEISMO_F) call write_seismogram(syn, nrec, trim(out_dir)//syn_tag)

!!$           if (WRITE_SPECTRAL_MAP_F) then
!!$              print *, 'compute and write out forward spectral map '
!!$              call write_spectral_map(syn, nrec, rglob, trim(out_dir)//'spectral_forward',WRITE_SPECTRA_F)
!!$           endif

           !=========================
           ! MEASUREMENTS AND ADJOINT SOURCES

           ! adjoint source function
           allocate(adj_syn(NSTEP,NCOMP,nrec))
           adj_syn = 0.0

           ! initialize time window for adjoint source to be the entire record
           allocate(tstart(nrec),tend(nrec))
           tstart(:) = 0.0
           tend(:)   = NSTEP*DT

           if (ISURFACE == 1) then
              ! compute time windows for measurements based on a HOMOGENEOUS reference model
              ! Window with should be wide enough to capture synthetic and data waveforms,
              ! but narrow enough to exclude suprious boundary reflections.
              print *
              print *, 'cut times for adjoint sources'
              open(19,file=trim(out_dir)//'measurement_time_windows.dat',status='unknown')
              do i = 1,nrec

                 ! d: distance from the source to the receivers (assuming HOMOGENEOUS model)
                 d = sqrt( (x(sglob(1)) - x(rglob(i)))**2 + (z(sglob(1)) - z(rglob(i)))**2  )
                 tcen = tshift + d/beta0

                 ! wider window for data and synthetics that have LARGE time shifts
                 tstart(i) = tcen - HWIN1
                 tend(i)   = tcen + HWIN1

                 write(19,'(i8,3f12.4)') i, tstart(i), tcen, tend(i)
              enddo
              close(19)
           endif

           !stop 'testing'

           if (WRITE_SEISMO_RECONSTRUCT) then
              allocate(data_recon(NSTEP,NCOMP,nrec))
              data_recon = 0.0
           endif

           ! construct the adjoint source (IKER in wave2d_constants.f90)
           !call make_adjoint_source(nrec, syn, tstart, tend, adj_syn, data)
           call mtm_adj(ievent, nrec, syn, tstart, tend, adj_syn, data, data_recon)

           if (WRITE_SEISMO_RECONSTRUCT) call write_seismogram(data_recon, nrec, trim(out_dir)//'dat_recon')

           ! specify which components you want to send back
           if (ISURFACE == 0) then
              if (REV_X == 0) adj_syn(:,1,:) = 0.0
              if (REV_Y == 0) adj_syn(:,2,:) = 0.0
              if (REV_Z == 0) adj_syn(:,3,:) = 0.0
           endif

           ! write out adjoint source time function at the receivers
           stfadj_tag = 'stfadj'
           if (WRITE_STF_A) call write_seismogram(adj_syn, nrec, trim(out_dir)//stfadj_tag)

           ! OUTPUT ASCII FILES --> SAC FILES (make_sac_files.pl)
           ! (1) data, (2) synthetics, (3) adjoint source time function

!!$           if (WRITE_SEISMO_F) then
!!$              filename1 = trim(script_dir)//'make_sac_files.csh'
!!$              filename2 = 'make_sac_files.pl'
!!$              open(19,file=filename1,status='unknown')
!!$              if (IKER <= 4) write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', &
!!$                   trim(out_dir),' ', trim(data_tag)  ,' ','1', tshift
!!$              write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', &
!!$                   trim(out_dir),' ', trim(syn_tag)   ,' ','1', tshift
!!$              if (WRITE_STF_A)   write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', &
!!$                   trim(out_dir),' ', trim(stfadj_tag),' ','1', tshift
!!$              close(19)
!!$              call system('chmod 755 scripts/make_sac_files.csh ; scripts/make_sac_files.csh')
!!$           endif

           ! misfit for each receiver for one event for all components
           print *,'---------------------------'
           print *,' misfit for IKER = ',IKER
           print *,' istep, imod, irun, ievent : ', istep, imod, irun, ievent
           print *,' data_norm2(ievent) = ', sum(chi_data(ievent,:,:,1))
           print *,'---------------------------'

!!$           ! write ALL adjoint sources to file
!!$           do ipick = 0,6
!!$              call mtm_adj(ipick, ievent, nrec, syn, tstart, tend, adj_syn, data)
!!$              write(stfadj_tag,'(a,i1)') 'stfadj-', ipick
!!$              call write_source_function(nrec, ti, adj_syn, rglob, trim(out_dir)//stfadj_tag)
!!$           enddo
!!$
!!$           ! plot ALL adjoint sources
!!$           filename1 = trim(script_dir)//'make_sac_files.csh'
!!$           filename2 = 'make_sac_files.pl'
!!$           open(19,file=filename1,status='unknown')
!!$           write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', trim(out_dir),' ', trim(data_tag)  ,' ',"1",tshift
!!$           write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', trim(out_dir),' ', trim(syn_tag)   ,' ',"1",tshift
!!$           do ipick = 0,6
!!$              write(stfadj_tag,'(a,i1)') 'stfadj-', ipick
!!$              write(19,'(7a,f12.6)') trim(script_dir)//trim(filename2),' ', trim(out_dir),' ', trim(stfadj_tag),' ',"1",tshift
!!$           enddo
!!$           close(19)
!!$           call system('chmod 755 scripts/make_sac_files.csh ; scripts/make_sac_files.csh')

           !stop 'testing adjoint sources'

           !------------
           ! we always evaluate the kernels for the present model
           ! we only evaluate the kernel for the test model if itest==1 or POLY_ORDER==3

           if ( (itest == 0 .or. POLY_ORDER == 3) .and. COMPUTE_KERNELS) then
              print *
              print *, 'compute the kernel via adjoint wavefield interaction'
              print *

              isolver = 3
              allocate(atype_kernel(NGLLX,NGLLZ,NSPEC),btype_kernel(NGLLX,NGLLZ,NSPEC),rtype_kernel(NGLLX,NGLLZ,NSPEC))

              ! initialize source perturbation time series
              allocate(three_source_model(NSTEP,NCOMP,nsrc,NVAR_SOURCE))
              !allocate(three_source_model(NSTEP,NCOMP,nsrc,10))             ! testing
              three_source_model(:,:,:,:) = 0.0

              ! kernels
              ! samp enters as the forward source time function,
              ! and returns as the (time-reversed) adjoint seismograms
              call solver(isolver, idata, &
                   nsrc, sglob, ispec_src, hxis_store, hgammas_store, samp, &
                   nrec, rglob, ispec_rec, hxir_store, hgammar_store, adj_syn, &
                   trim(last_frame_name), absorb_field, &
                   atype_kernel, btype_kernel, rtype_kernel, three_source_model, stf_syn, f0)

              !----------------------

              ! gradient vector for source parameters
              ! NOTE : FILLS ONLY NVAR_SOURCE ENTRIES OF source_gradient FOR EACH ievent
              print *
              print *, ' integrate time series to make the gradient for the source (unscaled)'
              do icomp = 1,NCOMP
                 do ipar = 1,NVAR_SOURCE   ! mislocations AND origin time error
                 !do ipar = 2,3              ! test mislocations only
                 !do ipar = 1,1              ! test origin time error only
                 !do ipar = 1,10             ! testing (see solver.f90)

                    !itemp = (ievent-1)*3 + ipar
                    itemp = index_source(ipar,ievent)

                    write(*,'(3i8,1e20.10)') icomp, ipar, itemp, DT*sum(three_source_model(:,icomp,1,ipar))
                    source_gradient(itemp) = source_gradient(itemp) + DT*sum(three_source_model(:,icomp,1,ipar))

                    ! write out all the time series used for the source inversion
                    write(filename,'(a,i2.2,a,i1.1,a,i2.2)') trim(out_dir)//'pert_src_e',ievent,'_c',icomp,'_p',ipar
                    call write_time_series(three_source_model(:,icomp,1,ipar), NSTEP, filename)
                 enddo
              enddo

              deallocate(three_source_model)
              !----------------------

              ! smooth EACH event kernel for the data subspace method
              if (ISMOOTH_EVENT_KERNEL == 1) then
                 call smooth_function(btype_kernel, GAMMA_SMOOTH_KERNEL, kbeta_smooth)
                 !call smooth_function(atype_kernel, GAMMA_SMOOTH_KERNEL, kbulk_smooth)

              else
                 kbeta_smooth = btype_kernel
                 !kbulk_smooth = atype_kernel
              endif

!!$              ! smooth EACH event kernel for the data subspace method
!!$              if (ISMOOTH_EVENT_KERNEL == 1) then

                 ! write smoothed beta kernel to file
                 open(unit = 13, file = trim(out_dir)//'kernel_basis_smooth', status = 'unknown',iostat=ios)
                 if (ios /= 0) stop 'Error writing smoothed kernel to disk'
                 k = 1
                 ktemp(:) = 0.0
                 do ispec = 1,NSPEC
                    do j = 1,NGLLZ
                       do i = 1,NGLLX
                          iglob = ibool(i,j,ispec)
                          write(13,'(4e16.6)') x_plot(iglob), z_plot(iglob), &
                             kbeta_smooth(i,j,ispec), btype_kernel(i,j,ispec)

                          ktemp(k) = kbeta_smooth(i,j,ispec)   ! used to compute norm of gradient
                          k = k+1
                       enddo
                    enddo
                 enddo
                 close(13)

                 ! store UNBALANCED norm of gradient -- cov_model
                 itemp1 = index_source(1,ievent)
                 itemp2 = index_source(2,ievent)
                 itemp3 = index_source(3,ievent)
                 gradient_data_all(ievent,1) = sum( ktemp**2 * cov_model(m_inds(1,1):m_inds(1,2)) )
                 gradient_data_all(ievent,2) = source_gradient(itemp1)**2 * cov_model(itemp1+nmod_str)
                 gradient_data_all(ievent,3) = source_gradient(itemp2)**2 * cov_model(itemp2+nmod_str)
                 gradient_data_all(ievent,4) = source_gradient(itemp3)**2 * cov_model(itemp3+nmod_str)

                 ! write the gradient parts FOR EACH EVENT to file
                 open(unit=18,file=trim(out_dir)//'gradient_data_unbalanced.dat',status='unknown')
                 write(18,'(4e20.10)') gradient_data_all(ievent,1), gradient_data_all(ievent,2), &
                      gradient_data_all(ievent,3), gradient_data_all(ievent,4)
                 close(18)

!!$              endif

              ! summed kernel (= misfit kernel)
              btype_kernel_sum = btype_kernel_sum + kbeta_smooth

              deallocate(atype_kernel,btype_kernel,rtype_kernel)

              ! write out adjoint seismograms at the original sources
              if (WRITE_SEISMO_A) call write_seismogram(samp, nsrc, trim(out_dir)//'synadj')

           endif

           ! deallocate variables associated with one event

           if (WRITE_SEISMO_RECONSTRUCT) deallocate(data_recon)

           deallocate(data, adj_syn, absorb_field)
           deallocate(syn, tstart, tend)

           deallocate(ispec_src_dat, xi_src_dat, gamma_src_dat)
           deallocate(hxisd_store, hgammasd_store)
           deallocate(samp_dat, sglob_dat)

           deallocate(ispec_src, xi_src, gamma_src)
           deallocate(hxis_store, hgammas_store)
           deallocate(samp, sglob)

        !=====================
        enddo  ! LOOP 3 (ievent, events)
        !=====================

        print *
        print *, ' Done with all the events for irun =', irun
        print *, ' Now we update the model, or compute a test model'

        !  stop 'testing'

        !-----------------------------------------------------
        ! misfit values

        ! sum of all chi values (used in CG inversion) (nsrc = 1)
        data_norm = sum(chi_data(:,:,:,1))
        if (INCLUDE_MODEL_NORM) then
           chi_val = 0.5 * ( data_norm + model_norm )
        else
           chi_val = 0.5 * data_norm
        endif

        ! VARIANCE REDUCTION
        ! NOTE 1: if the new chi (chi_val) is larger than the previous chi (chi_k_val), then VR < 0.
        ! NOTE 2: consider whether chi_data_stop should be removed
        if (istep >= 2) then
           !var_red_val = 100.0 * ( 1.0 - ( (chi_val - chi_data_stop)**2 / (chi_k_val - chi_data_stop)**2 ) )
           !var_red_val = -log( abs(chi_val - chi_data_stop) / abs(chi_k_val - chi_data_stop) )
           var_red_val = log( chi_k_val / chi_val )

           open(19,file=trim(out_dir1)//'variance_reduction.dat',status='unknown')
           write(19,'(1f20.10)') var_red_val, VAR_RED_MIN
           close(19)
        endif

        ! write chi_data values to various files
        call write_chi(out_dir1,nevent,nrec)

        ! write all nevent * nrec * ncomp measurements to file
        open(unit=18,file=trim(out_dir1)//'measure_vec.dat',status='unknown')
        do i = 1,nmeas
           write(18,'(5e20.10)') (measure_vec(i,j),j=1,5)
        enddo
        close(18)
        deallocate(measure_vec)

        if (ISOURCE_LOG) then
           write(91,*) 'SUMMED MISFIT FOR ALL EVENTS'
           write(91,'(a30,1i16)') ' N-data : ', nmeas_run
           write(91,'(a30,1e16.8)') ' data_norm2(m) : ', data_norm
           write(91,'(a30,3e16.8)') ' model_norm2(m) : ', model_norm, model_norm_struct, model_norm_source
           write(91,'(a30,3e16.8)') ' model_norm2_target(m) : ', &
              model_norm_target, model_norm_target_struct, model_norm_target_source
           write(91,'(a30,3e16.8)') ' model_norm2_diff(m) : ', &
              model_norm_diff, model_norm_diff_struct, model_norm_diff_source
           !write(91,'(a30,1e16.8)') ' chi_model_stop(mtarget) : ', chi_model_stop
           write(91,'(a30,1e16.8)') ' chi_data_stop : ', chi_data_stop
           write(91,'(a30,1e16.8)') ' chi(m) : ', chi_val
        endif

        ! STOPPING CRITERIA for CG algorithm
        ! Program will stop after files are written (search stop_cg)
        stop_cg = .false.
        if ( data_norm <= 2.0 * chi_data_stop) then
            stop_cg = .true.
            print *, ' data_norm <= 2(chi_data_stop) : ', data_norm, 2.0*chi_data_stop
        endif

        !if ( itest==0 .and. (chi_val >= chi_k_val) ) then
        if ( itest == 0 .and. var_red_val < VAR_RED_MIN ) then
            stop_cg = .true.
            print *, '   chi_k_val : ', chi_k_val
            print *, '     chi_val : ', chi_val
            print *, ' var_red_val : ', var_red_val
            print *, ' VAR_RED_MIN : ', VAR_RED_MIN
            print *, ' var_red_val < VAR_RED_MIN'
        endif

!!$        if (istep == 0) chi_val_0 = chi_val
!!$        if ( chi_val <= chi_val_0 * CONV_STOP ) then
!!$            print *, ' chi_val : ', chi_val
!!$            print *, ' chi_val_0 * CONV_STOP : ', chi_val_0 * CONV_STOP
!!$            print *, ' chi_val <= chi_val_0 * CONV_STOP'
!!$            stop ' convergence criteria is met.'
!!$        endif


        ! enter this section if you want to compute the gradient and/or iterate
        if ( COMPUTE_KERNELS ) then

           print *, ' Now we compute the gradient of the misfit function (using the misfit kernel)'

           !-----------------------------------------
           ! ASSEMBLE THE MODEL VECTOR

           if (itest == 0) then

              m0(:) = 0.0
              temp_local1(:,:,:) = log( beta_syn(:,:,:) / beta0 )
              call local2mvec(temp_local1, nmod_src, m_src_syn, nmod, m0)

!!$              ! m_src_syn: w.r.t initial source values (m_src_syn_vec_initial)
!!$              m0(:) = 0.0
!!$              temp_local1(:,:,:) = log( beta_syn(:,:,:) / beta0 )
!!$              m_src_syn(:) = m_src_syn_vec(:) - m0_vec_initial(nmod_str+1 : nmod)
!!$              call local2mvec(temp_local1, nmod_src, m_src_syn, nmod, m0)

           endif

           !-----------------------------------------------------
           ! COMPUTE THE GRADIENT OF THE MISFIT function, if the present model is not
           ! a test model or if the CG polynomial is a cubic function

           if ( itest == 0 .or. POLY_ORDER == 3 ) then

              print *, ' Computing the gradient of the misfit function for a given model'
              !gradient(:) = 0.0 ; gradient_data(:) = 0.0 ; gradient_model(:) = 0.0  ! initialized above

              ! source gradient
              open(unit=19,file=trim(out_dir1)//'source_gradient.dat',status='unknown')
              write(19,'(1e20.10)') (source_gradient(i), i = 1,nmod_src)
              close(19)

              if (INV_STRUCT_BETA == 1) then

                 ! norm of parts of unbalanced gradient for all events
                 ! NOTE: within each event directory there is also gradient_data_unbalanced.dat
                 open(unit=18,file=trim(out_dir1)//'gradient_data_all_unbalanced.dat',status='unknown')
                 do i = 1,nevent
                    write(18,'(4e20.10)') gradient_data_all(i,1), gradient_data_all(i,2), &
                       gradient_data_all(i,3), gradient_data_all(i,4)
                 enddo
                 close(18)

                 ! WRITE OUT summed event kernel (= misfit kernel)
                 !call local2global(atype_kernel_sum, temp_global1)
                 call local2global(btype_kernel_sum, temp_global2)
                 open(19,file=trim(out_dir1)//'summed_ker.dat',status='unknown')
                 do iglob = 1,NGLOB
                    write(19,'(3e16.6)') x_plot(iglob), z_plot(iglob), temp_global2(iglob)
                 enddo
                 close(19)

                 ! Smooth the summed event kernel (= misfit kernel) to remove spurious src-rec effects;
                 ! this is done via Gaussian convolution.
                 if ( ISMOOTH_MISFIT_KERNEL == 1) then
                    call smooth_function(btype_kernel_sum, GAMMA_SMOOTH_KERNEL, kbeta_smooth)
                 else
                    kbeta_smooth = btype_kernel_sum
                 endif
                 !kbulk_smooth = 0.0       ! testing for S-wave only

              endif  ! INV_STRUCT_BETA

              !------------------------
              ! assemble the gradient vector

              ! NOTE 1: gradient of data term of misfit function
              ! NOTE 2: structure gradient g_i = K_i * dA_i
              temp_local1 = kbeta_smooth * da_local
              call local2mvec(temp_local1, nmod_src, source_gradient, nmod, gradient_data)

              ! add the gradient term associated with the MODEL NORM in the misfit function
              ! KEY QUESTION: cov_model or cov_imetric? cov_model appears to perform slightly better
              if ( INCLUDE_MODEL_NORM ) then
                 if (itest == 0) then
                    gradient_model(:) = (m0(:) - mprior(:)) / cov_imetric(:)
                    !gradient_model(:) = (m0(:) - mprior(:)) / cov_model(:)
                 else
                    gradient_model(:) = (mt(:) - mprior(:)) / cov_imetric(:)
                    !gradient_model(:) = (mt(:) - mprior(:)) / cov_model(:)
                 endif
              endif

              ! KEY: if NOT inverting for structure or source, then set those parts of the gradient to zero
              if (INV_STRUCT_BETA == 0) then
                  gradient_data(m_inds(1,1) : m_inds(1,2)) = 0.0
                 gradient_model(m_inds(1,1) : m_inds(1,2)) = 0.0
              endif
              if (INV_SOURCE_T == 0) then
                  gradient_data(m_inds(2,1) : m_inds(2,2)) = 0.0
                 gradient_model(m_inds(2,1) : m_inds(2,2)) = 0.0
              endif
              if (INV_SOURCE_X == 0) then   ! both xs and ys
                  gradient_data(m_inds(3,1) : m_inds(4,2)) = 0.0
                 gradient_model(m_inds(3,1) : m_inds(4,2)) = 0.0
              endif

              ! overall gradient
              gradient = gradient_data + gradient_model

              ! NORMS OF THE GRADIENT (C^-1 is used, not C ; mprior is not used)
              imnorm = 0
!!$
!!$              ! normalize -- no cov weights
!!$               call compute_norm_sq(imnorm, ievent_min, ievent_max, nevent, index_source, nmod, &
!!$                      gradient, mprior, cov_model, gradient_norm_parts)
!!$              gnorm = sqrt( sum( gradient_norm_parts ) )
!!$              print *, 'normalizing gradient by ', gnorm
!!$              gradient(:) = gradient(:) / gnorm
!!$              filename1 = trim(out_dir1)//'gradient_norm_all_unit.dat'
!!$              call write_norm_sq(filename1,gradient_norm_parts)
!!$              !stop 'testing gradient normalization'

              ! write gradient vector to file (source and structure parameters)
              !write(19,'(1e20.10)') (gradient(i), i = 1,nmod)
              open(unit=19,file=trim(out_dir1)//'gradient_vec.dat',status='unknown')
              do i = 1,nmod
                 write(19,'(3e20.10)') gradient(i), gradient_data(i), gradient_model(i)
              enddo
              close(19)

              filename1 = trim(out_dir1)//'gradient_norm_all.dat'
              filename2 = trim(out_dir1)//'gradient_norm_data_all.dat'
              filename3 = trim(out_dir1)//'gradient_norm_model_all.dat'

              ! compute the balanced gradient norm
              ! NOTE: these two options should be equivalent (check to be sure)
              if (1 == 1) then    ! cov_imetric
                 call compute_norm_sq(filename1, imnorm, &
                      ievent_min, ievent_max, nevent, index_source, nmod, &
                      gradient, mprior, cov_imetric, gradient_norm_parts)
                 call compute_norm_sq(filename2, imnorm, &
                      ievent_min, ievent_max, nevent, index_source, nmod, &
                      gradient_data, mprior, cov_imetric, gradient_norm_data_parts)
                 call compute_norm_sq(filename3, imnorm, &
                      ievent_min, ievent_max, nevent, index_source, nmod, &
                      gradient_model, mprior, cov_imetric, gradient_norm_model_parts)
                 ! NOTE: gnorm( Cinv(m-mprior) ) should equal mnorm( m-mprior )

              else              ! cov_model plus weights
                 call compute_norm_sq(filename1, imnorm, &
                      ievent_min, ievent_max, nevent, index_source, nmod, &
                      gradient, mprior, cov_model, gradient_norm_parts, covg_weight_parts)
                 call compute_norm_sq(filename2, imnorm, &
                      ievent_min, ievent_max, nevent, index_source, nmod, &
                      gradient_data, mprior, cov_model, gradient_norm_data_parts, covg_weight_parts)
                 call compute_norm_sq(filename3, imnorm, &
                      ievent_min, ievent_max, nevent, index_source, nmod, &
                      gradient_model, mprior, cov_model, gradient_norm_model_parts, covg_weight_parts)
              endif

!!$              ! old scaling used in GJI-2007 paper
!!$              !scale_struct_gji = sqrt( norm_source_2 / (norm_bulk_2 + norm_beta_2) )
!!$              scale_struct_gji = sqrt( gradient_norm_source / gradient_norm_struct )
!!$              open(unit=19,file=trim(out_dir1)//'scaling_gji.dat',status='unknown')
!!$              write(19,'(1e20.10)') scale_struct_gji
!!$              close(19)

!!$              ! write gradient norms to file
!!$              filename1 = trim(out_dir1)//'gradient_norm_all.dat'
!!$              filename2 = trim(out_dir1)//'gradient_norm_data_all.dat'
!!$              filename3 = trim(out_dir1)//'gradient_norm_model_all.dat'
!!$              call write_norm_sq(filename1,gradient_norm_parts)
!!$              call write_norm_sq(filename2,gradient_norm_data_parts)
!!$              call write_norm_sq(filename3,gradient_norm_model_parts)

              ! SAME AS ABOVE, using cov_model WITHOUT weights
              filename1 = trim(out_dir1)//'gradient_norm_all_unbalanced.dat'
              filename2 = trim(out_dir1)//'gradient_norm_data_all_unbalanced.dat'
              filename3 = trim(out_dir1)//'gradient_norm_model_all_unbalanced.dat'
              call compute_norm_sq(filename1, imnorm, &
                 ievent_min, ievent_max, nevent, index_source, nmod, &
                 gradient, mprior, cov_model, gradient_norm_parts)
              call compute_norm_sq(filename2, imnorm, &
                 ievent_min, ievent_max, nevent, index_source, nmod, &
                 gradient_data, mprior, cov_model, gradient_norm_data_parts)
              call compute_norm_sq(filename3, imnorm, &
                 ievent_min, ievent_max, nevent, index_source, nmod, &
                 gradient_model, mprior, cov_model, gradient_norm_model_parts)

!!$              ! write (unbalanced) gradient norms to file
!!$              filename1 = trim(out_dir1)//'gradient_norm_all_unbalanced.dat'
!!$              filename2 = trim(out_dir1)//'gradient_norm_data_all_unbalanced.dat'
!!$              filename3 = trim(out_dir1)//'gradient_norm_model_all_unbalanced.dat'
!!$              call write_norm_sq(filename1,gradient_norm_parts)
!!$              call write_norm_sq(filename2,gradient_norm_data_parts)
!!$              call write_norm_sq(filename3,gradient_norm_model_parts)

              !stop 'TESTING GRADIENT NORMS'

              ! update balance using unweighted gradient norm -- see cov_imetric above
              ! NOTE 1: this will only use gradient for non-test models (itest=0)
              ! NOTE 2: do you want this for the first step only (istep=0) or each step
              if ( (istep == 0) .and. (0 == 1) ) then

                 ! NOTE: choose to balance using the full unbalanced gradient or the data term only
                 gbalance_parts(:) = 0.0
                 gbalance_parts = gradient_norm_data_parts
                 !gbalance_parts = gradient_norm_parts

                 open(unit=19,file=trim(out_dir1)//'cov_imetric_weights.dat',status='unknown')
                 do i = 1,NVAR
                    write(19,'(3e16.6)') gbalance_parts(i), covm_weight_parts(i), &
                        gbalance_parts(i) / covm_weight_parts(i)
                 enddo
                 close(19)

                 cov_imetric(m_inds(1,1):m_inds(1,2)) = cov_model(m_inds(1,1):m_inds(1,2)) &
                     / (gbalance_parts(1) / covm_weight_parts(1))
                 cov_imetric(m_inds(2,1):m_inds(2,2)) = cov_model(m_inds(2,1):m_inds(2,2)) &
                     / (gbalance_parts(2) / covm_weight_parts(2))
                 cov_imetric(m_inds(3,1):m_inds(3,2)) = cov_model(m_inds(3,1):m_inds(3,2)) &
                     / (gbalance_parts(3) / covm_weight_parts(3))
                 cov_imetric(m_inds(4,1):m_inds(4,2)) = cov_model(m_inds(4,1):m_inds(4,2)) &
                     / (gbalance_parts(4) / covm_weight_parts(4))

                 ! write the new cov_imetric to file (over-write previous file)
                 open(unit=19,file=trim(out_dir2)//'cov_model_diagonal.dat',status='unknown')
                 do i = 1,nmod
                    write(19,'(2e20.10)') cov_model(i), cov_imetric(i)
                 enddo
                 close(19)
              endif

!!$              ! scaling for joint inversions
!!$              if ( istep==0 .and. 0==1 ) then
!!$
!!$                 ! scaling for joint source-structure inversions
!!$
!!$                 ! Formulas used for GJI paper
!!$                 !norm_bulk_2   = sum( kbulk_smooth**2 * da_local )
!!$                 !norm_beta_2   = sum( kbeta_smooth**2 * da_local )
!!$                 !norm_source_2 = sum( ( source_gradient(:) * m_scale_src_all(:) )**2 )
!!$
!!$                 ! NORM-SQUARED of the steepest ascent vector : sqrt( ghat C ghat )
!!$                 norm_bulk_2   = sum( kbulk_smooth(:,:,:)**2 * da_local(:,:,:) * (m_scale_str(1))**2 * AREA )
!!$                 norm_beta_2   = sum( kbeta_smooth(:,:,:)**2 * da_local(:,:,:) * (m_scale_str(2))**2 * AREA )
!!$                 norm_source_2 = sum( ( source_gradient(:) * m_scale_src_all(:) )**2 * dnevent_run * NVAR_SOURCE )   ! m_scale_src_all NOT ALLOCATED
!!$
!!$                 ! This shows the balance used in the GJI-2007 paper
!!$                 scale_struct_gji = sqrt( norm_source_2 / (norm_bulk_2 + norm_beta_2) )
!!$
!!$                 !scale_struct = scale_struct_gji
!!$                 scale_struct = 1.0
!!$
!!$                 write(*,*) 'Scaling for joint inversions:'
!!$                 write(*,'(a24,1e14.4)') ' norm_bulk            : ', sqrt( norm_bulk_2 )
!!$                 write(*,'(a24,1e14.4)') ' norm_beta            : ', sqrt( norm_beta_2 )
!!$                 write(*,'(a24,1e14.4)') ' norm_source          : ', sqrt( norm_source_2 )
!!$                 write(*,'(a24,1e14.4)') ' scale_struct_gji (F) : ', scale_struct_gji
!!$                 write(*,'(a24,1e14.4)') ' scale_struct (F)     : ', scale_struct
!!$
!!$                 open(unit=19,file=trim(out_dir1)//'scaling_values.dat',status='unknown')
!!$                 write(19,*) 'GJI_PAPER = ', GJI_PAPER
!!$                 write(19,*) 'IRUNZ = ', IRUNZ
!!$                 write(19,'(6i14)') 1, NLOCAL, NLOCAL+1, nmod_str, nmod_str+1, nmod
!!$                 write(19,'(2i14)') nevent_run, NLOCAL
!!$                 write(19,'(5e14.6)') sum(da_local(:,:,:)), sum(da_local_vec(:)), sum(da_global(:)), LENGTH*HEIGHT, AREA
!!$                 write(19,'(3e14.6)') minval(da_local_vec(:)), maxval(da_local_vec(:)), sum(da_local_vec(:))/NLOCAL
!!$                 write(19,'(3e14.6)') minval(da_global(:)), maxval(da_global(:)), sum(da_global(:))/NGLOB
!!$                 write(19,'(a20,2e14.6)') 'cov_str_fac : ', cov_str_fac, cov_str_fac**2
!!$                 write(19,'(a20,2e14.6)') 'cov_src_fac : ', cov_src_fac, cov_src_fac**2
!!$                 write(19,*) 'SIGMA VALUES :'
!!$                 write(19,'(5e14.6)') m_scale_str(1), m_scale_str(2), m_scale_src(1), m_scale_src(2), m_scale_src(3)
!!$                 write(19,'(5e14.6)') (m_scale_str(1))**2, (m_scale_str(2))**2, &
!!$                      (m_scale_src(1))**2, (m_scale_src(2))**2, (m_scale_src(3))**2
!!$                 write(19,'(5e14.6)') sigma_bulk, sigma_beta, sigma_xs, sigma_zs, sigma_ts
!!$                 write(19,'(5e14.6)') (sigma_bulk)**2, (sigma_beta)**2, (sigma_xs)**2, (sigma_zs)**2, (sigma_ts)**2
!!$                 write(19,'(a20,4e14.6)') 'scale_struct : ', scale_struct, scale_struct**2, 1.0/scale_struct, (1.0/scale_struct)**2
!!$                 write(19,'(a20,4e14.6)') 'scale_struct_gji : ', scale_struct_gji , scale_struct_gji **2, &
!!$                      1.0/scale_struct_gji , (1.0/scale_struct_gji )**2
!!$                 close(19)
!!$              endif

           endif  ! itest == 0 .or. POLY_ORDER == 3

           !====================
           ! CG algorithm: update search direction and model
           !
           ! For details, see Tarantola (2005), Section 6.22
           ! HAT variables: g0, gk, gt
           ! NON-HAT variables: p0, pk, m0, mt

           print *, ' Entering CG algorithm to compute new model or test model'

           if (itest == 0) then      ! if the present kernel is for a REAL model

              chi_k_val = chi_val
              gk(:) = gradient(:)

              !-----------------------------------------

              ! update the current search direction p0 (p0 and pk are NON-hat)
              if (istep == 0) then
                 beta_val = 0.0
                 p0(:) = 0.0      ! initialize
              else
                 beta_val_top = sum((gk(:) - g0(:)) * (cov_imetric(:)*gk(:)) )
                 beta_val_bot = sum(g0(:) * (cov_imetric(:)*g0(:)) )
                 beta_val = beta_val_top / beta_val_bot
              endif
              pk(:) = -cov_imetric(:) * gk(:) + beta_val * p0(:)

              ! write gradient vectors (before re-assigning g0 and gk)
              if (1 == 1) then
                 open(unit=19,file=trim(out_dir1)//'cg_grad_vectors.dat',status='unknown')
                 do i = 1,nmod
                    write(19,'(4e16.6)') g0(i), gk(i), p0(i), pk(i)
                 enddo
                 close(19)
              endif

              ! KEY: test value for line-search to get test model
              ! You must be careful to define the metric for the dot product operation.
              ! Here, we use a diagonal covariance matrix; for a full covariance matrix,
              ! the operation would be : g_k C^kl p_l

              ! The minimum of the test-model parameter should go through 0.5, which
              ! is the misfit if all the data are fit within the uncertainties specified by
              ! the data covariance matrix.
              mu_val = chi_data_stop
              !mu_val = 0.0

              lam_t_val_bot = sum( gk(:) * pk(:) )    ! ghat dot p
              lam_t_val = 2.0*(mu_val - chi_k_val) / lam_t_val_bot

              ! alternative approach (more ad hoc)
              !lam_t_val = -2.0*(1.0 - mu_val)*chi_k_val  / sum( gk(:) * pk(:) )

              !lam_t_val = -2.0*chi_k_val / dot_product(gk, pk)      ! GJI paper
              !istep_switch = 6
              !if (istep < istep_switch)  lam_t_val = -2.0*chi_k_val / dot_product(gk, pk)    ! quadratic extrapolation
              !if (istep >= istep_switch) lam_t_val =    -chi_k_val / dot_product(gk, pk)    ! linear extrapolation

              ! compute test model
              mt(:)  = m0(:) + lam_t_val * pk(:)
              itest  = 1

!!$           ! Re-set the source parameters for the events NOT being run to the initial values.
!!$           ! Otherwise, you will perturb all 25 events, even though they are not being run,
!!$           ! and this will affect the computation of the model norm.
!!$           do i = 1,nevent
!!$              itemp1 = nmod_str + (i-1)*3 + 1
!!$              itemp2 = nmod_str + (i-1)*3 + 2
!!$              itemp3 = nmod_str + (i-1)*3 + 3
!!$              if ( i < ievent_min .or. i > ievent_max ) then
!!$                 mt(itemp1) = 0.0
!!$                 mt(itemp2) = 0.0
!!$                 mt(itemp3) = 0.0
!!$              endif
!!$           enddo

              mt_vec(1 : NLOCAL) = beta0 * exp( mt(1 : NLOCAL) )
              mt_vec(nmod_str+1 : nmod) = mt(nmod_str+1 : nmod)

!!$              ! STRUCTURE PARAMETERS (bulk, then beta)
!!$              if (INV_STRUCT_BETA == 0) then
!!$                 mt_vec(1 : nmod_str) = m0_vec_initial(1 : nmod_str)     ! initial structure
!!$              else
!!$                 mt_vec(1 : NLOCAL) = beta0 * exp( mt(1 : NLOCAL) )
!!$                 !mt_vec(1 : NLOCAL)          = bulk0 * exp( mt(1 : NLOCAL) )
!!$                 !mt_vec(NLOCAL+1 : nmod_str) = beta0 * exp( mt(NLOCAL+1 : nmod_str) )
!!$              endif
!!$
!!$              ! SOURCE PARAMETERS
!!$              if (INV_SOURCE_T == 0 .and. INV_SOURCE_X == 0) then
!!$                 mt_vec(nmod_str+1 : nmod) = m0_vec_initial(nmod_str+1 : nmod)     ! initial source
!!$              else
!!$                 mt_vec(nmod_str+1 : nmod) = mt(nmod_str+1 : nmod) + m0_vec_initial(nmod_str+1 : nmod)
!!$              endif

              print *, 'lam_t_val = ', sngl(lam_t_val)

              ! save vectors for next iteration
              g0(:) = gk(:)    ! gradient
              p0(:) = pk(:)    ! search direction vector

              ! write scalars to file
              open(unit=19,file=trim(out_dir1)//'cg_test_vals.dat',status='unknown')
              write(19,'(5e16.8)') mu_val, lam_t_val_bot, lam_t_val, beta_val, chi_k_val
              close(19)

           else if (itest == 1) then  ! if present kernel is for a test model

              chi_t_val = chi_val

              ! a cubic or quadratic fit requires at least 5 values
              xx1 = 0.0
              xx2 = lam_t_val
              yy1 = chi_k_val
              yy2 = chi_t_val
              g1  = sum(g0(:) * pk(:))
              !g1  = dot_product(g0, pk)              ! GJI paper

              if (POLY_ORDER == 3) then
                 ! use cubic polynomial: six values gives an analytical minimum
                 ! see Matlab scripts cubic_min_4.m and cubic_min.m

                 ! compute gradient of test function (needed for cubic polynomial)
                 gt(:) = gradient(:)

                 g2 = sum(gt(:) * pk(:))
                 !g2 = dot_product(gt,pk)              ! GJI paper

                 ! coefficients of the cubic polynomial
                 Pa = ( -2.0*(yy2-yy1) + (g1+g2)*(xx2-xx1) ) / (xx2-xx1)**3
                 Pb = ( 3.0*(yy2-yy1) - (2.0*g1 + g2)*(xx2-xx1) ) / (xx2-xx1)**2
                 Pc = g1
                 Pd = yy1

                 ! get the analytical minimum
                 qfac = Pb**2.0 - 3.0*Pa*Pc
                 if (Pa /= 0.0 .and. qfac >= 0.0) then
                    xmin = (-Pb + sqrt(qfac)) / (3.0*Pa)
                 else if (Pa == 0.0 .and. Pb /= 0.0) then
                    xmin = -Pc/(2.0*Pb)
                 else
                    print *, 'Pa, Pb, Pc, Pd, qfac: ',Pa,Pb,Pc,Pd,qfac
                    print *, 'xmin: ',xmin
                    stop 'check the cubic input polynomial'
                 endif

                 ! write cubic function parameters to file
                 open(unit=19,file=trim(out_dir1)//'cubic_poly.dat',status='unknown')
                 write(19,'(11e16.8)') xx1,xx2,yy1,yy2,g1,g2,Pa,Pb,Pc,Pd,xmin
                 close(19)

              else
                 ! use quadratic polynomial: five values gives an analytical minimum
                 ! see Matlab script quad_min_4.m

                 ! coefficients of the quadratic polynomial (ax^2 + bx + c)
                 Pa = ((yy2 - yy1) - g1*(xx2 - xx1)) / (xx2**2 - xx1**2)
                 Pb = g1
                 Pc = yy1 - Pa*xx1**2 - Pb*xx1

                 ! get the analytical minimum (the vertex)
                 if (Pa /= 0.0) then
                    xmin = -Pb / (2.0*Pa)
                 else
                    print *, 'Pa, Pb, Pc: ',Pa,Pb,Pc
                    print *, 'xmin: ',xmin
                    stop 'check the quadratic input polynomial'
                 endif

                 ! write quadratic function parameters to file
                 open(unit=19,file=trim(out_dir1)//'quad_poly.dat',status='unknown')
                 write(19,'(9e16.8)') xx1,xx2,yy1,yy2,g1,Pa,Pb,Pc,xmin
                 close(19)

              endif  ! POLY_ORDER == 3

              ! compute updated model
              lam_0_val = xmin
              m0(:) = m0(:) + lam_0_val * pk(:)
              itest = 0

!!$           ! Re-set the source parameters for the events NOT being run to the initial values.
!!$           ! Otherwise, you will perturb all 25 events, even though they are not being run,
!!$           ! and this will affect the computation of the model norm.
!!$           do i = 1,nevent
!!$              itemp1 = nmod_str + (i-1)*3 + 1
!!$              itemp2 = nmod_str + (i-1)*3 + 2
!!$              itemp3 = nmod_str + (i-1)*3 + 3
!!$              if ( i < ievent_min .or. i > ievent_max ) then
!!$                 m0(itemp1) = 0.0
!!$                 m0(itemp2) = 0.0
!!$                 m0(itemp3) = 0.0
!!$              endif
!!$           enddo

              m0_vec(1 : NLOCAL) = beta0 * exp( m0(1 : NLOCAL) )
              m0_vec(nmod_str+1 : nmod) = m0(nmod_str+1 : nmod)

!!$              ! STRUCTURE PARAMETERS
!!$              if (INV_STRUCT_BETA == 0) then
!!$                 m0_vec(1 : nmod_str) = m0_vec_initial(1 : nmod_str)   ! initial structure
!!$              else
!!$                 m0_vec(1 : NLOCAL) = beta0 * exp( m0(1 : NLOCAL) )
!!$                 !m0_vec(1 : NLOCAL)          = bulk0 * exp( m0(1 : NLOCAL) )
!!$                 !m0_vec(NLOCAL+1 : nmod_str) = beta0 * exp( m0(NLOCAL+1 : nmod_str) )
!!$              endif
!!$
!!$              ! SOURCE PARAMETERS
!!$              if (INV_SOURCE_T == 0 .and. INV_SOURCE_X == 0) then
!!$                 m0_vec(nmod_str+1 : nmod) = m0_vec_initial(nmod_str+1 : nmod)      ! initial source
!!$              else
!!$                 m0_vec(nmod_str+1 : nmod) = m0(nmod_str+1 : nmod) + m0_vec_initial(nmod_str+1 : nmod)
!!$              endif

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

           ! FUTURE WORK: if the structure or source parameters are unrealistic,
           !    then we should exit the program before the NEXT wave simulation.
           !    This way, we at least write the unrealistic model to file.

           !====================

           ! stop CG algorithm (search stop_cg)
           if ( stop_cg ) then
              stop ' convergence criteria is met.'
           endif

           !call system("xv model.jpg &")

        endif   ! COMPUTE_KERNELS

        print *, ' Done with LOOP 2 iteration -- on to next model.'

        !==================================
     enddo  ! LOOP 2 (istep, models)
     !==================================

     print *, ' Done with LOOP 2.'

     if (ISOURCE_LOG) close(91)  ! source log file

     ! deallocate event and receiver variables
     if (ISURFACE == 1) deallocate(fglob)

     deallocate(rglob)
     deallocate(ispec_rec,xi_rec,gamma_rec)
     deallocate(x_rec,z_rec)
     deallocate(hxir_store, hgammar_store)
     if (ISURFACE == 1) deallocate(x_rec_lon,z_rec_lat)

     deallocate(eglob_syn)
     deallocate(ispec_eve_syn,xi_eve_syn,gamma_eve_syn)
     deallocate(x_eve_syn,z_eve_syn)
     deallocate(hxie_store,hgammae_store)
     if (ISURFACE == 1) deallocate(x_eve_lon_syn,z_eve_lat_syn)

     deallocate(eglob_dat)
     deallocate(ispec_eve_dat,xi_eve_dat,gamma_eve_dat)
     deallocate(x_eve_dat,z_eve_dat)
     deallocate(hxied_store, hgammaed_store)
     if (ISURFACE == 1) deallocate(x_eve_lon_dat,z_eve_lat_dat)

     deallocate(otime_syn,otime_dat)

     ! conjugate-gradient variables
     deallocate(m0,mt,mdiff,mtarget,mprior)
     deallocate(cov_imetric,cov_model)
     !deallocate(icov_metric,icov_model)
     deallocate(gradient,gradient_data,gradient_model)
     deallocate(gradient_data_all)
     deallocate(g0,gt,gk,p0,pk)
     deallocate(source_gradient)
     deallocate(m0_vec,mt_vec)
     deallocate(m_src_syn,m_src_dat,m_src_prior)
     !deallocate(m_src_syn_vec,m_src_dat_vec)
     !deallocate(m0_vec_initial,m_src_syn_vec_initial,m_src_syn_initial)
     !deallocate(m_scale_src_all)

     deallocate(cov_data,index_data,index_source)
     if (ADD_DATA_ERRORS) deallocate(measure_pert_vec)

     ! print output
     print *
     print *, ' max time-step is ', NSTEP
     print *, '  time-step DT is ', sngl(DT)
     print *, '      max time is ', sngl(dble(NSTEP)*DT), ' s,', sngl(dble(NSTEP)*DT/60.0), ' min'
     print *, ' space-step dh is ', sngl(dh)
     print *
     print *, '         beta0 is ', sngl(beta0)
     print *, '         bulk0 is ', sngl(bulk0)
     print *, '        alpha0 is ', sngl(alpha0)
     print *, '          hdur is ', sngl(hdur)
     print *, '    time shift is ', sngl(tshift)
     print *
     print *, '         NGLOB ', NGLOB
     print *, '        NLOCAL ', NLOCAL
     print *, '          nmod ', nmod
     print *, '      nmod_str ', nmod_str
     print *, '      nmod_src ', nmod_src
     print *

  enddo  ! LOOP 1 (iq, tests)

end program wave2d

