program test_smooth

  use wave2d_variables  ! global variables
  use wave2d_solver     ! mesher

  implicit none

  ! smoothing
  integer :: igaus
  double precision :: d, dmin, dist2, dtrsh2, xtar, ztar, gamma, xcen, zcen
  double precision, dimension(NGLOB) :: k_rough_global, k_smooth_global
  double precision, dimension(NGLOB) :: k_gaus_global, k_gaus_global_ex, k_gaus_int_global
  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: k_temp, k_gaus_local, k_smooth_local, k_rough_local

  ! corners of the grid
  double precision, dimension(4) :: corners_utm_x, corners_utm_z, corners_lon, corners_lat
  double precision :: sq_south, sq_east, sq_north, sq_west, sq_fac, da_mean, sda_mean, m_scale_str, mfac

  ! corners of the elements, in a global vector
  !integer, dimension(NSPEC_CORNER) :: ielement_corner
  !logical, dimension(NGLOB) :: mask_ibool
  !integer :: iglob1, iglob2, iglob3, iglob4

  double precision :: xtemp,ztemp,dx,dz,xmesh,zmesh
  double precision :: junk1,junk2,junk3
  double precision :: temp1,temp2,temp3,temp4,temp5,temp6,temp7

  integer :: i, j, k, iq, iglob, itemp, itype, icorner, ispec
  !integer :: isolver, irun0, irun, idat, iopt, ispec, istep, istep_switch, imod

  !********* PROGRAM STARTS HERE *********************

  ! mesher (gets Jacobian)
  call mesher()

  ! compute da_global, da_local, valence, ielement_corner
  call mesher_additional()

  da_mean = sum(da_global)/NGLOB
  sda_mean = sqrt(da_mean)

  ! global mesh
  open(unit=15,file='global_mesh.dat',status='unknown')
  do iglob = 1,NGLOB
     write(15,'(2i8,3e18.8)') iglob, valence(iglob), x(iglob), z(iglob), da_global(iglob)
  enddo
  close(15)

!!$  ! da vector and valence vector
!!$  da_global(:) = 0.
!!$  valence(:) = 0
!!$  do ispec = 1,NSPEC
!!$    do j = 1,NGLLZ
!!$      do i = 1,NGLLX
!!$        iglob = ibool(i,j,ispec)
!!$        da_local(i,j,ispec) = wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
!!$
!!$        !da_global(iglob) = da_global(iglob) + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
!!$        da_global(iglob) = da_global(iglob) + da_local(i,j,ispec)
!!$        valence(iglob) = valence(iglob) + 1
!!$      enddo
!!$    enddo
!!$  enddo
!!$  da_mean = sum(da_global)/NGLOB
!!$  sda_mean = sqrt(da_mean)
!!$
!!$  ! global mesh
!!$  open(unit=15,file='global_mesh.dat',status='unknown')
!!$  do iglob = 1,NGLOB
!!$    write(15,'(2i8,3e18.8)') iglob, valence(iglob), x(iglob), z(iglob), da_global(iglob)
!!$  enddo
!!$  close(15)
!!$
!!$  ! GLL points defining the corners of elements
!!$  mask_ibool(:) = .false.
!!$  do ispec = 1,NSPEC
!!$     iglob1 = ibool(1,1,ispec)
!!$     iglob2 = ibool(NGLLX,1,ispec)
!!$     iglob3 = ibool(1,NGLLZ,ispec)
!!$     iglob4 = ibool(NGLLX,NGLLZ,ispec)
!!$
!!$      if (.not. mask_ibool(iglob1)) mask_ibool(iglob1) = .true.
!!$      if (.not. mask_ibool(iglob2)) mask_ibool(iglob2) = .true.
!!$      if (.not. mask_ibool(iglob3)) mask_ibool(iglob3) = .true.
!!$      if (.not. mask_ibool(iglob4)) mask_ibool(iglob4) = .true.
!!$  enddo
!!$
!!$  k = 0
!!$  ielement_corner(:) = 0
!!$  do iglob = 1,NGLOB
!!$      if (mask_ibool(iglob)) then
!!$         k = k+1
!!$         ielement_corner(k) = iglob
!!$      endif
!!$  enddo

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
           write(15,'(4i10,3e18.8)') ispec, i, j, valence(iglob), da_global(iglob)/1e6, x(iglob)/1000., z(iglob)/1000.
        enddo
     enddo
     close(15)

     ! GLL points defining the corners of elements
     open(unit=15,file='element_corners.dat',status='unknown')
     do i = 1,NSPEC_CORNER
        iglob = ielement_corner(i)
        write(15,'(1i10,2e18.8)') iglob, x(iglob)/1000., z(iglob)/1000.
     enddo
     close(15)

     stop 'testing'
  endif

  print *
  write(*,'(a,1f20.10)') '     da_min  (m^2) : ', minval(da_global)
  write(*,'(a,1f20.10)') '     da_mean (m^2) : ', da_mean
  write(*,'(a,1f20.10)') '     da_max  (m^2) : ', maxval(da_global)
  write(*,*)             ' sum [ da_global ] : ', sum(da_global)
  write(*,*)             ' sum [ da_local  ] : ', sum(da_local)
  write(*,*)             '   LENGTH * HEIGHT : ', LENGTH * HEIGHT
  print *
  print *, ' GLL weights:'
  do i=1,NGLLX
     print *, wxgll(i)
  enddo
  do i=1,NGLLZ
     print *, wzgll(i)
  enddo

  !-----------------------
  ! apply smoothing

  gamma = 60.0d+03     ! scalelength of smoothing Gaussian

  itype = 3        ! type of integration to test

  ! read in sample function
  k_rough_global(:) = 0.

  ! read in test function
  ! NOTE: this has to have been greated using the same mesh
  open(unit=19,file='/net/denali/scratch1/carltape/OUTPUT/run_0460/fun_smooth.dat',status='unknown')
  do iglob = 1,NGLOB
     read(19,'(8e16.6)') temp1, temp2, temp3, temp4, &
          k_rough_global(iglob), temp5, temp6, temp7
  enddo
  close(19)

  ! construct local version of the unsmoothed kernel
  if (itype == 1) then
     k_rough_local(:,:,:) = 0.
     do ispec = 1,NSPEC
        do j = 1,NGLLZ
           do i = 1,NGLLX
              itemp = ibool(i,j,ispec)
              k_rough_local(i,j,ispec) = k_rough_global(itemp)
           enddo
        enddo
     enddo
  endif

  ! Smoothing function is a Gaussian whose full-width is given by gamma;
  ! all points outside d^2 (dtrsh2) are set to zero.
  dtrsh2 = (1.5*gamma)**2

  ! EXAMPLE Gaussian smoothing function for one point
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
  k_gaus_global_ex(:) = 0.
  do iglob = 1,NGLOB
     dist2 = (xcen - x(iglob))**2 + (zcen - z(iglob))**2
     if (dist2 <= dtrsh2) &
          k_gaus_global_ex(iglob) = (4./(PI*gamma**2)) * exp(-4.*dist2 / (gamma**2))
  enddo

  !------------------------------------
  ! Compute the SMOOTHED kernel by convolving a Gaussian with the UNSMOOTHED kernel.
  ! This involves integrating NGLOB products between a Gaussian and the unsmoothed kernel.

  print *, 'convolving the kernel with a Gaussian...'

  k_smooth_global(:) = 0.

  ! loop over CONVOLUTION POINTS, either GLL points or simply the element corners
  do iglob = 1,NGLOB

     !do icorner = 1,NSPEC_CORNER
     !   iglob = ielement_corner(icorner)

     if (mod(iglob,500) == 0) write(*,*) iglob, ' out of ', NGLOB

     ! compute a Gaussian centered at the iglob point
     ! (part of the Gaussian may be outside the grid)
     xcen = x(iglob)
     zcen = z(iglob)
     if (itype == 1) then

        k_gaus_local(:,:,:) = 0.
        do ispec = 1,NSPEC
           do j = 1,NGLLZ
              do i = 1,NGLLX
                 itemp = ibool(i,j,ispec)
                 dist2 = (xcen - x(itemp))**2 + (zcen - z(itemp))**2
                 if (dist2 <= dtrsh2) then
                    k_gaus_local(i,j,ispec) = (4./(PI*gamma**2)) * exp(-4.*dist2 / (gamma**2))
                    !k_gaus_global(itemp) = k_gaus_local(i,j,ispec)
                 endif
              enddo
           enddo
        enddo
     else

        k_gaus_global(:) = 0.
        do i = 1,NGLOB
           dist2 = (xcen - x(i))**2 + (zcen - z(i))**2
           if (dist2 <= dtrsh2) &
                k_gaus_global(i) = (4./(PI*gamma**2)) * exp(-4.*dist2 / (gamma**2))
        enddo
     endif

     ! There are two primary steps:
     ! (1) Integrate the Gaussian over the grid; this provides the normalization for the Gaussian
     ! and accounts for Gaussians that are partially outside the grid.
     ! (2) Integrate the product of the Gaussian and the rough function.

     if (itype == 1) then            ! local integration with local arrays

        k_gaus_int_global(iglob) = sum( k_gaus_local(:,:,:) * da_local(:,:,:) )
        k_smooth_global(iglob) = sum( k_rough_local(:,:,:) * k_gaus_local(:,:,:) * da_local(:,:,:) ) / k_gaus_int_global(iglob)

     else if (itype == 2) then        ! local integration with global array

        k_temp(:,:,:) = 0.
        do ispec = 1,NSPEC
           do j = 1,NGLLZ
              do i = 1,NGLLX
                 itemp = ibool(i,j,ispec)
                 k_temp(i,j,ispec) = k_gaus_global(itemp) * da_local(i,j,ispec)
              enddo
           enddo
        enddo
        k_gaus_int_global(iglob) = sum( k_temp(:,:,:) )

        k_temp(:,:,:) = 0.
        do ispec = 1,NSPEC
           do j = 1,NGLLZ
              do i = 1,NGLLX
                 itemp = ibool(i,j,ispec)
                 k_temp(i,j,ispec) = k_rough_global(itemp) * k_gaus_global(itemp) * da_local(i,j,ispec)
              enddo
           enddo
        enddo
        k_smooth_global(iglob) = sum( k_temp(:,:,:) ) / k_gaus_int_global(iglob)

     else if (itype == 3) then       ! global integration with global arrays

        k_gaus_int_global(iglob) = sum( k_gaus_global(:) * da_global(:) )
        k_smooth_global(iglob) = sum( k_rough_global(:) * k_gaus_global(:) * da_global(:) ) / k_gaus_int_global(iglob)

     endif

  enddo   ! iglob = 1,NGLOB

  ! write smooth-related functions to file
  open(unit=19,file='fun_smooth.dat',status='unknown')
  do iglob = 1,NGLOB
     !do icorner = 1,NSPEC_CORNER
     !iglob = ielement_corner(icorner)

     write(19,'(9e16.6)') x_lon(iglob), z_lat(iglob), x(iglob), z(iglob), &
          k_rough_global(iglob), k_gaus_global_ex(iglob), &
          k_smooth_global(iglob), k_rough_global(iglob) - k_smooth_global(iglob), &
          k_gaus_int_global(iglob)
  enddo
  close(19)

  ! plot rough function, Gaussian filter, smooth function, and residual
  !filename1 = 'get_smooth.csh'
  !filename2 = trim(script_dir)//'plot_smoothed_function.pl'
  !open(19,file=filename1,status='unknown')
  !write(19,'(5a,1e16.6)') trim(filename2),' ', trim(out_dir1),' ', &
  !   trim(file_smooth), gamma
  !close(19)
  !call system('chmod 755 get_smooth.csh ; get_smooth.csh')

end program test_smooth

