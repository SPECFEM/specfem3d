program wave2d_cmap

  use wave2d_variables
  use wave2d_solver
  use wave2d_sub
  use wave2d_sub2
  implicit none

! wavefields and seismograms
  double precision, dimension(:,:,:), allocatable :: samp, data, syn, adj_syn

! socal coast and shelf data
  integer :: UTM_PROJECTION_ZONE, FINITE_SOURCE
  integer :: nshelf,ncoast,itemp
  double precision :: d,dmin
  double precision, dimension(:), allocatable :: shelf_lat,shelf_lon,shelf_z,shelf_x
  double precision, dimension(:), allocatable :: coast_lat,coast_lon,coast_z,coast_x
  character(len=200) :: shelf_file,coast_file

! model variables
  integer :: Nfac
  double precision :: w_scale, afac

! additional velocity models
  double precision, dimension(NGLOB) :: c_glob_dat, c_glob_syn, m0, m0_vel

! regular mesh for computing phase velocity
!  integer :: ntemp, i_regular
!  double precision, dimension(:), allocatable :: x9,z9,x_lon9,z_lat9,da9
!  double precision, dimension(:), allocatable :: c_glob_dat9, c_glob_syn9, c_glob9
!  double precision :: fac

  double precision :: xtemp,ztemp,dx,dz,xmesh,zmesh
  double precision :: t_target, tcen
  double precision :: junk1,junk2,junk3

  double precision :: temp1,temp2,temp3,temp4

! smoothing
  integer :: igaus
  double precision :: dist2, dtrsh2, xtar, ztar, xcen, zcen
  double precision, dimension(NGLOB) :: k_rough_global, k_smooth_global
  double precision, dimension(NGLOB) :: k_gaus_global, k_gaus_global_ex, k_gaus_int_global, da
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: k_temp
  !double precision, dimension(NLOCAL) :: gk_local, m0_local, mt_local, mt2_local, da_local
  !integer, dimension(NLOCAL) :: ibool2
  character(len=200) :: file_smooth

  character(len=20)  :: data_tag,syn_tag,stf_tag,stfadj_tag
  character(len=200) :: srcfile,recfile,socal_map
  character(len=200) :: socal_dir, model_dir, script_dir, out_dir2, out_dir1, ev_dir
  character(len=200) :: filename, filename1, filename2, file_dat_c, file_syn_c, file_src_rec, file_stf

  integer :: i, j, k, itime, iglob, irec, isrc, ipick, ievent, icomp
  integer :: isolver, irun0, irun, idat, iopt, ispec, istep

  !********* PROGRAM STARTS HERE *********************

  irun0 = 420

  !out_dir2   = "OUTPUT/"
  in_dir     = "INPUT/"
  script_dir = "scripts/"

  socal_dir = '/home/carltape/socal_modes/local_modes/'
  model_dir = 'model_files/'

  UTM_PROJECTION_ZONE = 11

!--------------------------------------

  print *
  print *, ' max time-step is ', NSTEP
  print *, '            DT is ', sngl(DT)
  print *, '      max time is ', sngl(NSTEP*DT), ' s,', sngl(NSTEP*DT/60), ' min'
  print *

!--------------------------------------
! mesher

  ! set up model (gets Jacobian)
  call mesher()

  ! da vector
  da(:) = 0.
  do ispec = 1,NSPEC
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob     = ibool(i,j,ispec)
        da(iglob) = da(iglob) + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
      enddo
    enddo
  enddo

  print *
  write(*,'(a,2f20.10)') '     da(min/max) : ', minval(da), maxval(da)
  write(*,*)             '      sum [ da ] : ', sum(da)
  write(*,*)             ' LENGTH * HEIGHT : ', LENGTH * HEIGHT
  print *
  write(*,*)             '           SIGMA : ', SIGMA
  write(*,*)             '           irun0 : ', irun0
  write(*,*)             '            IKER : ', IKER
  print *
  print *, ' GLL weights:'
  do i=1,NGLLX
     print *, wxgll(i)
  enddo
  do i=1,NGLLZ
     print *, wzgll(i)
  enddo

  ! determine the UTM coordinates of your origin
  print *
  print *, 'UTM check for origin of mesh:'
  print *, LON_MIN,LAT_MIN
  call utm_geo(LON_MIN,LAT_MIN,utm_xmin,utm_zmin,UTM_PROJECTION_ZONE,ILONGLAT2UTM)
  print *, utm_xmin, utm_zmin
  call utm_geo(xtemp,ztemp,utm_xmin,utm_zmin,UTM_PROJECTION_ZONE,IUTM2LONGLAT)
  print *, xtemp,ztemp   ! should match LON_MIN,LAT_MIN
  print *

!!$  ! regular mesh or SPECFEM mesh
!!$  i_regular = 0
!!$
!!$  if (i_regular==1) then
!!$
!!$     ! enter new uniform mesh
!!$     dx = 40.0d+03
!!$     dz = dx
!!$     fac = 10.0d+03
!!$     k = 0
!!$     do xtemp = -fac,LENGTH+fac,dx
!!$        do ztemp = -fac,HEIGHT+fac,dz
!!$           k = k+1
!!$        enddo
!!$     enddo
!!$     ntemp = k
!!$
!!$  else
!!$     ntemp = NGLOB
!!$  endif
!!$
!!$  allocate(x9(ntemp),z9(ntemp),x_lon9(ntemp),z_lat9(ntemp))
!!$  allocate(c_glob9(ntemp),c_glob_syn9(ntemp),c_glob_dat9(ntemp))
!!$  allocate(da9(ntemp))
!!$
!!$  if (i_regular==1) then
!!$
!!$     k = 0
!!$     do xtemp = -fac,LENGTH+fac,dx
!!$        do ztemp = -fac,HEIGHT+fac,dz
!!$     !do i=utm_xmin-fac,utm_xmin+LENGTH+fac,dx
!!$     !   do j=utm_zmin-fac,utm_zmin+HEIGHT+fac,dz
!!$           k = k+1
!!$           x9(k) = xtemp
!!$           z9(k) = ztemp
!!$        enddo
!!$     enddo
!!$
!!$     da9(:) = LENGTH*HEIGHT/ntemp
!!$
!!$  else
!!$     x9(:) = x(:)
!!$     z9(:) = z(:)
!!$     da9(:) = da(:)
!!$  endif

  print *, NGLOB, ' gridpoints for which we compute a phase velocity'
  print *, 'x-range is: ', minval(x), maxval(x)
  print *, 'z-range is: ', minval(z), maxval(z)

  !stop 'testing'

  ! convert global gridpoint mesh coordinates to lat-lon
  x_lon(:) = 0.
  z_lat(:) = 0.
  call mesh_geo(NGLOB,x_lon,z_lat,x,z,UTM_PROJECTION_ZONE,IMESH2LONLAT)
  !write(*,'(2f16.6)') (x_lon(iglob), z_lat(iglob), iglob=1,NGLOB)

  !stop 'testing'

!--------------------------------------
! phase velocity model

  t_target = 2*hdur     ! target period for phase velocity map

     ! write lat-lon gridpoints to file
     filename = trim(socal_dir) // 'socal_input.dat'
     print *, 'Output data file is ',trim(filename)
     open(unit=15,file=filename,status='unknown')
     write(15,*) t_target
     write(15,*) 0
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

  !============================================
  ! smooth the Gaussian according to SIGMA

  k_rough_global(:) = c_glob(:)/c0 - 1.

  ! Gaussian half-width controlling the smoothing (m)
  dtrsh2 = (3.*SIGMA)**2  ! all points outside d^2 are set to zero

  ! EXAMPLE Gaussian smoothing function for one point
  ! find the closest gridpoint to the target point
  xtar = 0.7*LENGTH
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
        k_gaus_global_ex(iglob) = (1./(2*PI*SIGMA**2)) * exp(-dist2 / (2.*SIGMA**2))
  enddo

  ! compute the SMOOTHED kernel
  k_smooth_global(:) = 0.
  k_gaus_int_global(:) = 0.
  do iglob = 1,NGLOB
     ! compute a Gaussian centered at the iglob point
     xcen = x(iglob)
     zcen = z(iglob)
     k_gaus_global(:) = 0.
     do i = 1,NGLOB
        dist2 = (xcen - x(i))**2 + (zcen - z(i))**2
        if (dist2 <= dtrsh2) &
           k_gaus_global(i) = (1./(2.*PI*SIGMA**2)) * exp(-dist2 / (2.*SIGMA**2))
     enddo

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

  do i=1,NGLOB
     print *, k_gaus_int_global(i)
  enddo

  !============================================

  !k_smooth_global(:) = c_glob(:)

  ! c-maps for data
  c_glob_dat(:) = c0*(1. + k_smooth_global(:) )

  ! write data phase velocity map to file
  file_dat_c = 'socal_vel_dat.dat'
  open(unit=18,file=file_dat_c,status='unknown')
  do iglob = 1,NGLOB
     write(18,*) sngl(x_lon(iglob)), sngl(z_lat(iglob)), sngl(c_glob_dat(iglob)), sngl(c_glob_dat(iglob)/c0 - 1.)
  enddo
  close(18)

  ! write syn phase velocity map to file
  ! (reference model is 0% perturbation)
  file_syn_c = 'socal_vel_syn.dat'
  open(unit=18,file=file_syn_c,status='unknown')
  do iglob = 1,NGLOB
     write(18,*) sngl(x_lon(iglob)), sngl(z_lat(iglob)), sngl(c_glob_syn(iglob)), sngl(c_glob_syn(iglob)/c0 - 1.)
  enddo
  close(18)

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

  ! write smooth-related functions to file
  file_smooth = 'fun_smooth.dat'
  open(unit=19,file=file_smooth,status='unknown')
  do iglob = 1,NGLOB
     write(19,'(6e16.6)') x(iglob), z(iglob), &
        k_rough_global(iglob), k_gaus_global_ex(iglob), &
        k_smooth_global(iglob), k_rough_global(iglob) - k_smooth_global(iglob)
  enddo
  close(19)

  ! plot rough function, Gaussian filter, smooth function, and residual
  filename1 = 'get_smooth.csh'
  filename2 = trim(script_dir)//'plot_smoothed_function.pl'
  open(19,file=filename1,status='unknown')
  write(19,'(5a,1e16.6)') trim(filename2),' ', '.',' ',trim(file_smooth), SIGMA
  close(19)
  call system('chmod 755 get_smooth.csh ; get_smooth.csh')

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

  !deallocate(x,z,x_lon,z_lat,c_glob,c_glob_syn,c_glob_dat)

end program wave2d_cmap

