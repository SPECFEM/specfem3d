!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


  subroutine prepare_gravity()

! precomputes gravity factors

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use gravity_perturbation, only: gravity_init

  implicit none

  ! local parameters
  double precision RICB,RCMB,RTOPDDOUBLEPRIME, &
    R80,R220,R400,R600,R670,R771,RMOHO,RMIDDLE_CRUST,ROCEAN
  double precision :: rspl_gravity(NR),gspl(NR),gspl2(NR)
  double precision :: radius,g,dg ! radius_km
  !double precision :: g_cmb_dble,g_icb_dble
  double precision :: rho,drhodr,vp,vs,Qkappa,Qmu
  integer :: nspl_gravity !int_radius
  integer :: i,j,k,iglob,ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing gravity"
    call flush_IMAIN()
  endif

  ! for gravity perturbation calculations
  ! sets up arrays for gravity field
  call gravity_init()

  ! sets up weights needed for integration of gravity
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        wgll_cube(i,j,k) = sngl( wxgll(i)*wygll(j)*wzgll(k) )
      enddo
    enddo
  enddo

  ! store g, rho and dg/dr=dg using normalized radius in lookup table every 100 m
  ! get density and velocity from PREM model using dummy doubling flag
  ! this assumes that the gravity perturbations are small and smooth
  ! and that we can neglect the 3D model and use PREM every 100 m in all cases
  ! this is probably a rather reasonable assumption
  if (GRAVITY) then
    ! allocates gravity arrays
    allocate(minus_deriv_gravity(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2156')
    allocate(minus_g(NGLOB_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2157')
    if (ier /= 0) stop 'error allocating gravity arrays'

    ! sets up spline table
    call make_gravity(nspl_gravity,rspl_gravity,gspl,gspl2, &
                          ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                          R771,RTOPDDOUBLEPRIME,RCMB,RICB)

    ! pre-calculates gravity terms for all global points
    do iglob = 1,NGLOB_AB

      ! normalized radius ( zstore values given in m, negative values for depth)
      radius = ( R_EARTH + zstore(iglob) ) / R_EARTH
      call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g)

      ! use PREM density profile to calculate gravity (fine for other 1D models)
      call model_prem_iso(radius,rho,drhodr,vp,vs,Qkappa,Qmu, &
                        RICB,RCMB,RTOPDDOUBLEPRIME, &
                        R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

      dg = 4.0d0*rho - 2.0d0*g/radius

      ! re-dimensionalize
      g = g * R_EARTH*(PI*GRAV*RHOAV) ! in m / s^2 ( should be around 10 m/s^2)
      dg = dg * R_EARTH*(PI*GRAV*RHOAV) / R_EARTH ! gradient d/dz g , in 1/s^2

      minus_deriv_gravity(iglob) = - dg
      minus_g(iglob) = - g ! in negative z-direction

      ! debug
      !if (iglob == 1 .or. iglob == 1000 .or. iglob == 10000) then
      !  ! re-dimensionalize
      !  radius = radius * R_EARTH ! in m
      !  vp = vp * R_EARTH*dsqrt(PI*GRAV*RHOAV)  ! in m / s
      !  rho = rho  * RHOAV  ! in kg / m^3
      !  print *,'gravity: radius=',radius,'g=',g,'depth=',radius-R_EARTH
      !  print *,'vp=',vp,'rho=',rho,'kappa=',(vp**2) * rho
      !  print *,'minus_g..=',minus_g(iglob)
      !endif
    enddo

  else
    ! allocates dummy gravity arrays
    allocate( minus_deriv_gravity(0), minus_g(0), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2158')
    if (ier /= 0) stop 'error allocating gravity arrays'
  endif

  ! compute the gravity integrals if needed
  if (GRAVITY_INTEGRALS) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...computing gravity integrals'
      call flush_IMAIN()
    endif
    call compute_gravity_integrals()
  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_gravity

!
!-------------------------------------------------------------------------------------------------
!

  ! compute gravity integrals of that slice, and then total integrals for the whole mesh

  subroutine compute_gravity_integrals()

  use constants

  use specfem_par

  implicit none

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  double precision :: jacobianl
  integer :: i,j,k,ispec,ispec_irreg,iglob,ier
  double precision :: xval,yval,zval
  double precision :: xval_squared,yval_squared,zval_squared
  double precision :: x_meshpoint,y_meshpoint,z_meshpoint
  double precision :: distance_squared,distance_cubed, &
                      three_over_distance_squared,one_over_distance_cubed,three_over_distance_fifth_power
  double precision :: common_multiplying_factor,common_mult_times_one_over,common_mult_times_three_over

  integer :: iobservation
  character(len=MAX_STRING_LEN) :: outputname

! read the observation surface
  x_observation(:) = 0.d0
  y_observation(:) = 0.d0
  z_observation(:) = 0.d0

  if (myrank == 0) then
    open(unit=IIN,file=OBSERVATION_GRID_FILE(1:len_trim(OBSERVATION_GRID_FILE)),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error opening file observation_grid_to_use_for_gravity.txt')
    do iobservation = 1,NTOTAL_OBSERVATION
      read(IIN,*) x_observation(iobservation),y_observation(iobservation),z_observation(iobservation)
    enddo
    close(unit=IIN)
  endif

! broadcast the observation surface read
  call bcast_all_dp(x_observation, NTOTAL_OBSERVATION)
  call bcast_all_dp(y_observation, NTOTAL_OBSERVATION)
  call bcast_all_dp(z_observation, NTOTAL_OBSERVATION)

! initialize the gravity arrays
  g_x(:) = 0.d0
  g_y(:) = 0.d0
  g_z(:) = 0.d0

  G_xx(:) = 0.d0
  G_yy(:) = 0.d0
  G_zz(:) = 0.d0
  G_xy(:) = 0.d0
  G_xz(:) = 0.d0
  G_yz(:) = 0.d0

  ! calculates volume of all elements in mesh
  do ispec = 1,NSPEC_AB

    ! print information about number of elements done so far
    if (myrank == 0 .and. (mod(ispec,NSPEC_DISPLAY_INTERVAL) == 0 .or. ispec == 1 .or. ispec == NSPEC_AB)) then
       write(IMAIN,*) 'for gravity integrals, ',ispec,' elements computed out of ',NSPEC_AB
       ! write time stamp file to give information about progression of the calculation of gravity integrals
       write(outputname,"('/timestamp_gravity_calculations_ispec',i7.7,'_out_of_',i7.7)") ispec,NSPEC_AB
       ! timestamp file output
       open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write')
       write(IOUT,*) ispec,' elements done for gravity calculations out of ',NSPEC_AB
       close(unit=IOUT)
    endif

    ispec_irreg = irregular_element_number(ispec)

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ! computes the Jacobian
          if (ispec_irreg /= 0) then
            ! irregular element
            xixl = xixstore(i,j,k,ispec_irreg)
            xiyl = xiystore(i,j,k,ispec_irreg)
            xizl = xizstore(i,j,k,ispec_irreg)
            etaxl = etaxstore(i,j,k,ispec_irreg)
            etayl = etaystore(i,j,k,ispec_irreg)
            etazl = etazstore(i,j,k,ispec_irreg)
            gammaxl = gammaxstore(i,j,k,ispec_irreg)
            gammayl = gammaystore(i,j,k,ispec_irreg)
            gammazl = gammazstore(i,j,k,ispec_irreg)
            ! do this in double precision for accuracy
            jacobianl = 1.d0 / dble(xixl*(etayl*gammazl-etazl*gammayl) &
                          - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                          + xizl*(etaxl*gammayl-etayl*gammaxl))
          else
            !regular element
            jacobianl = 1.d0 / dble(xix_regular*xix_regular*xix_regular)
          endif

          if (CHECK_FOR_NEGATIVE_JACOBIANS .and. jacobianl <= ZERO) &
            stop 'error: negative Jacobian found in integral calculation'

          iglob = ibool(i,j,k,ispec)
          x_meshpoint = xstore(iglob)
          y_meshpoint = ystore(iglob)
          z_meshpoint = zstore(iglob)

          weight = wxgll(i)*wygll(j)*wzgll(k)
          common_multiplying_factor = jacobianl * weight * rhostore(i,j,k,ispec) * GRAV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! beginning of loop on all the data to create
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! loop on all the points in the observation surface
          do iobservation = 1,NTOTAL_OBSERVATION

            xval = x_meshpoint - x_observation(iobservation)
            yval = y_meshpoint - y_observation(iobservation)
            zval = z_meshpoint - z_observation(iobservation)

            xval_squared = xval**2
            yval_squared = yval**2
            zval_squared = zval**2

            distance_squared = xval_squared + yval_squared + zval_squared
            distance_cubed = distance_squared * sqrt(distance_squared)

            three_over_distance_squared = 3.d0 / distance_squared
            one_over_distance_cubed = 1.d0 / distance_cubed
            three_over_distance_fifth_power = three_over_distance_squared * one_over_distance_cubed

            common_mult_times_one_over = common_multiplying_factor * one_over_distance_cubed
            common_mult_times_three_over = common_multiplying_factor * three_over_distance_fifth_power

            g_x(iobservation) = g_x(iobservation) + common_mult_times_one_over * xval
            g_y(iobservation) = g_y(iobservation) + common_mult_times_one_over * yval
            g_z(iobservation) = g_z(iobservation) + common_mult_times_one_over * zval

            G_xx(iobservation) = G_xx(iobservation) &
              + common_mult_times_one_over * (xval_squared * three_over_distance_squared - 1.d0)
            G_yy(iobservation) = G_yy(iobservation) &
              + common_mult_times_one_over * (yval_squared * three_over_distance_squared - 1.d0)
            G_zz(iobservation) = G_zz(iobservation) &
              + common_mult_times_one_over * (zval_squared * three_over_distance_squared - 1.d0)

            G_xy(iobservation) = G_xy(iobservation) + common_mult_times_three_over * xval*yval
            G_xz(iobservation) = G_xz(iobservation) + common_mult_times_three_over * xval*zval
            G_yz(iobservation) = G_yz(iobservation) + common_mult_times_three_over * yval*zval

          enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! end of loop on all the data to create
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        enddo
      enddo
    enddo
  enddo

    ! the result is displayed in Eotvos = 1.e+9 s-2
    G_xx(:) = G_xx(:) * SI_UNITS_TO_EOTVOS
    G_yy(:) = G_yy(:) * SI_UNITS_TO_EOTVOS
    G_zz(:) = G_zz(:) * SI_UNITS_TO_EOTVOS
    G_xy(:) = G_xy(:) * SI_UNITS_TO_EOTVOS
    G_xz(:) = G_xz(:) * SI_UNITS_TO_EOTVOS
    G_yz(:) = G_yz(:) * SI_UNITS_TO_EOTVOS

    ! use an MPI reduction to compute the total value of the integral into a temporary array
    ! and then copy it back into the original array
    call sum_all_1Darray_dp(g_x,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) g_x(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(g_y,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) g_y(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(g_z,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) g_z(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_xx,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_xx(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_yy,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_yy(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_zz,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_zz(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_xy,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_xy(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_xz,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_xz(:) = temporary_array_for_sum(:)

    call sum_all_1Darray_dp(G_yz,temporary_array_for_sum,NTOTAL_OBSERVATION)
    if (myrank == 0) G_yz(:) = temporary_array_for_sum(:)

  !--- print number of points and elements in the mesh for each region
  if (myrank == 0) then

      temporary_array_for_sum(:) = sqrt(g_x(:)**2 + g_y(:)**2 + g_z(:)**2)
      write(IMAIN,*)
      write(IMAIN,*) 'minval of norm of g vector on whole observation surface = ',minval(temporary_array_for_sum),' m.s-2'
      write(IMAIN,*) 'maxval of norm of g vector on whole observation surface = ',maxval(temporary_array_for_sum),' m.s-2'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_xx on whole observation surface = ',minval(G_xx),' Eotvos'
      write(IMAIN,*) 'maxval of G_xx on whole observation surface = ',maxval(G_xx),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_yy on whole observation surface = ',minval(G_yy),' Eotvos'
      write(IMAIN,*) 'maxval of G_yy on whole observation surface = ',maxval(G_yy),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_zz on whole observation surface = ',minval(G_zz),' Eotvos'
      write(IMAIN,*) 'maxval of G_zz on whole observation surface = ',maxval(G_zz),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_xy on whole observation surface = ',minval(G_xy),' Eotvos'
      write(IMAIN,*) 'maxval of G_xy on whole observation surface = ',maxval(G_xy),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_xz on whole observation surface = ',minval(G_xz),' Eotvos'
      write(IMAIN,*) 'maxval of G_xz on whole observation surface = ',maxval(G_xz),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'minval of G_yz on whole observation surface = ',minval(G_yz),' Eotvos'
      write(IMAIN,*) 'maxval of G_yz on whole observation surface = ',maxval(G_yz),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) 'Minval and maxval of trace of G, which in principle should be zero:'
      write(IMAIN,*)
      temporary_array_for_sum(:) = abs(G_xx(:) + G_yy(:) + G_zz(:))
      write(IMAIN,*) 'minval of abs(G_xx + G_yy + G_zz) on whole observation surface = ',minval(temporary_array_for_sum),' Eotvos'
      write(IMAIN,*) 'maxval of abs(G_xx + G_yy + G_zz) on whole observation surface = ',maxval(temporary_array_for_sum),' Eotvos'

      write(IMAIN,*)
      write(IMAIN,*) '-----------------------------'
      write(IMAIN,*)
      write(IMAIN,*) 'displaying the fields computed at observation point = ',iobs_receiver,' out of ',NTOTAL_OBSERVATION
      write(IMAIN,*)
      write(IMAIN,*) 'computed g_x  = ',g_x(iobs_receiver),' m.s-2'
      write(IMAIN,*) 'computed g_y  = ',g_y(iobs_receiver),' m.s-2'
      write(IMAIN,*) 'computed g_z  = ',g_z(iobs_receiver),' m.s-2'
      write(IMAIN,*)
      write(IMAIN,*) 'computed norm of g vector = ',sqrt(g_x(iobs_receiver)**2 + g_y(iobs_receiver)**2 + &
                                                                 g_z(iobs_receiver)**2),' m.s-2'

      write(IMAIN,*)
      write(IMAIN,*) 'computed G_xx = ',G_xx(iobs_receiver),' Eotvos'
      write(IMAIN,*) 'computed G_yy = ',G_yy(iobs_receiver),' Eotvos'
      write(IMAIN,*) 'computed G_zz = ',G_zz(iobs_receiver),' Eotvos'
      write(IMAIN,*)
      write(IMAIN,*) 'G tensor should be traceless, G_xx + G_yy + G_zz = 0.'
      write(IMAIN,*) 'Actual sum obtained = ',G_xx(iobs_receiver) + G_yy(iobs_receiver) + G_zz(iobs_receiver)
      if (max(abs(G_xx(iobs_receiver)),abs(G_yy(iobs_receiver)),abs(G_zz(iobs_receiver))) > TINYVAL) &
           write(IMAIN,*) ' i.e., ',sngl(100.d0*abs(G_xx(iobs_receiver) + G_yy(iobs_receiver) + G_zz(iobs_receiver)) / &
                                     max(abs(G_xx(iobs_receiver)),abs(G_yy(iobs_receiver)),abs(G_zz(iobs_receiver)))), &
                                     '% of max(abs(G_xx),abs(G_yy),abs(G_zz))'
      write(IMAIN,*)
      write(IMAIN,*) 'computed G_xy = ',G_xy(iobs_receiver),' Eotvos'
      write(IMAIN,*) 'computed G_xz = ',G_xz(iobs_receiver),' Eotvos'
      write(IMAIN,*) 'computed G_yz = ',G_yz(iobs_receiver),' Eotvos'

      ! save the results
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_g_x_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) g_x(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_g_y_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) g_y(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_g_z_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) g_z(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_norm_of_g_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) sqrt(g_x(iobservation)**2 + g_y(iobservation)**2 + g_z(iobservation)**2)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_xx_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_xx(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_yy_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_yy(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_zz_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_zz(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_xy_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_xy(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_xz_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_xz(iobservation)
      enddo
      close(unit=IOUT)

      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/results_G_yz_for_GMT.txt',status='unknown',action='write')
      do iobservation = 1,NTOTAL_OBSERVATION
        write(IOUT,*) G_yz(iobservation)
      enddo
      close(unit=IOUT)

  endif

  end subroutine compute_gravity_integrals

