!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine prepare_timerun()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_movie
  
  implicit none

! user info
  if(myrank == 0) then

    write(IMAIN,*)
    if(TOPOGRAPHY) then
      write(IMAIN,*) 'incorporating surface topography'
    else
      write(IMAIN,*) 'no surface topography'
    endif

    write(IMAIN,*)
    if(ATTENUATION) then
      write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'
      if(USE_OLSEN_ATTENUATION) then
        write(IMAIN,*) 'using Olsen''s attenuation'
      else
        write(IMAIN,*) 'not using Olsen''s attenuation'
      endif
    else
      write(IMAIN,*) 'no attenuation'
    endif

    write(IMAIN,*)
    if(ANISOTROPY) then
      write(IMAIN,*) 'incorporating anisotropy'
    else
      write(IMAIN,*) 'no anisotropy'
    endif

    write(IMAIN,*)
    if(OCEANS) then
      write(IMAIN,*) 'incorporating the oceans using equivalent load'
    else
      write(IMAIN,*) 'no oceans'
    endif

  endif

! synchronize all the processes before assembling the mass matrix
! to make sure all the nodes have finished to read their databases
  call sync_all()

! the mass matrix needs to be assembled with MPI here once and for all
  call assemble_MPI_scalar_ext_mesh(NPROC,NGLOB_AB,rmass, &
         buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
         nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,my_neighbours_ext_mesh, &
         request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

  if(myrank == 0) write(IMAIN,*) 'end assembling MPI mass matrix'

! check that mass matrix is positive
  if(minval(rmass(:)) <= 0.) call exit_MPI(myrank,'negative mass matrix term')
  if(OCEANS .and. minval(rmass_ocean_load(:)) <= 0.) &
       call exit_MPI(myrank,'negative ocean load mass matrix term')

! for efficiency, invert final mass matrix once and for all in each slice
  if(OCEANS) rmass_ocean_load(:) = 1. / rmass_ocean_load(:)
  rmass(:) = 1.0 / rmass(:)

! if attenuation is on, shift PREM to right frequency
! rescale mu in PREM to average frequency for attenuation

  if(ATTENUATION) then

! get and store PREM attenuation model
    do iattenuation = 1,NUM_REGIONS_ATTENUATION

      call get_attenuation_model(myrank,iattenuation,tau_mu_dble, &
        tau_sigma_dble,beta_dble,one_minus_sum_beta_dble,factor_scale_dble)

! distinguish between single and double precision for reals
      if(CUSTOM_REAL == SIZE_REAL) then
        tau_mu(iattenuation,:) = sngl(tau_mu_dble(:))
        tau_sigma(iattenuation,:) = sngl(tau_sigma_dble(:))
        beta(iattenuation,:) = sngl(beta_dble(:))
        factor_scale(iattenuation) = sngl(factor_scale_dble)
        one_minus_sum_beta(iattenuation) = sngl(one_minus_sum_beta_dble)
      else
        tau_mu(iattenuation,:) = tau_mu_dble(:)
        tau_sigma(iattenuation,:) = tau_sigma_dble(:)
        beta(iattenuation,:) = beta_dble(:)
        factor_scale(iattenuation) = factor_scale_dble
        one_minus_sum_beta(iattenuation) = one_minus_sum_beta_dble
      endif
    enddo

! rescale shear modulus according to attenuation model
    !pll
    do ispec = 1,NSPEC_AB
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX

! use scaling rule similar to Olsen et al. (2003)          
!! We might need to fix the attenuation part for the anisotropy case
!! At this stage, we turn the ATTENUATION flag off always, and still keep mustore
            if(USE_OLSEN_ATTENUATION) then
              vs_val = mustore(i,j,k,ispec) / rho_vs(i,j,k,ispec)
              call get_attenuation_model_Olsen_sediment( vs_val, iselected )
            else                        
! takes iflag set in (CUBIT) mesh         
              iselected = iflag_attenuation_store(i,j,k,ispec)
            endif
            
! scales only mu             
            scale_factor = factor_scale(iselected)
            mustore(i,j,k,ispec) = mustore(i,j,k,ispec) * scale_factor
            
          enddo
        enddo
      enddo
    enddo

! obsolete, old way...
!pll 
!   do ispec = 1,NSPEC_AB
!    if(not_fully_in_bedrock(ispec)) then
!      do k=1,NGLLZ
!        do j=1,NGLLY
!          do i=1,NGLLX
!
!! distinguish attenuation factors
!   if(flag_sediments(i,j,k,ispec)) then
!
!! use constant attenuation of Q = 90
!! or use scaling rule similar to Olsen et al. (2003)
!! We might need to fix the attenuation part for the anisotropy case
!! At this stage, we turn the ATTENUATION flag off always, and still keep mustore
!     if(USE_OLSEN_ATTENUATION) then
!       vs_val = mustore(i,j,k,ispec) / rho_vs(i,j,k,ispec)
!! use rule Q_mu = constant * v_s
!       Q_mu = OLSEN_ATTENUATION_RATIO * vs_val
!       int_Q_mu = 10 * nint(Q_mu / 10.)
!       if(int_Q_mu < 40) int_Q_mu = 40
!       if(int_Q_mu > 150) int_Q_mu = 150
!
!       if(int_Q_mu == 40) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_40
!       else if(int_Q_mu == 50) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_50
!       else if(int_Q_mu == 60) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_60
!       else if(int_Q_mu == 70) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_70
!       else if(int_Q_mu == 80) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_80
!       else if(int_Q_mu == 90) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_90
!       else if(int_Q_mu == 100) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_100
!       else if(int_Q_mu == 110) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_110
!       else if(int_Q_mu == 120) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_120
!       else if(int_Q_mu == 130) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_130
!       else if(int_Q_mu == 140) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_140
!       else if(int_Q_mu == 150) then
!         iattenuation_sediments = IATTENUATION_SEDIMENTS_150
!       else
!         stop 'incorrect attenuation coefficient'
!       endif
!
!     else
!       iattenuation_sediments = IATTENUATION_SEDIMENTS_90
!     endif
!
!     scale_factor = factor_scale(iattenuation_sediments)
!   else
!     scale_factor = factor_scale(IATTENUATION_BEDROCK)
!   endif
!
!      mustore(i,j,k,ispec) = mustore(i,j,k,ispec) * scale_factor
!
!          enddo
!        enddo
!      enddo
!    endif
!    enddo
    
  endif ! ATTENUATION

! allocate seismogram array
  if (nrec_local > 0) then
    allocate(seismograms_d(NDIM,nrec_local,NSTEP))
    allocate(seismograms_v(NDIM,nrec_local,NSTEP))
    allocate(seismograms_a(NDIM,nrec_local,NSTEP))
! initialize seismograms
    seismograms_d(:,:,:) = 0._CUSTOM_REAL
    seismograms_v(:,:,:) = 0._CUSTOM_REAL
    seismograms_a(:,:,:) = 0._CUSTOM_REAL
    if (SIMULATION_TYPE == 2) then
    ! allocate Frechet derivatives array
      allocate(Mxx_der(nrec_local),Myy_der(nrec_local),Mzz_der(nrec_local),Mxy_der(nrec_local), &
               Mxz_der(nrec_local),Myz_der(nrec_local), sloc_der(NDIM,nrec_local))
      Mxx_der = 0._CUSTOM_REAL
      Myy_der = 0._CUSTOM_REAL
      Mzz_der = 0._CUSTOM_REAL
      Mxy_der = 0._CUSTOM_REAL
      Mxz_der = 0._CUSTOM_REAL
      Myz_der = 0._CUSTOM_REAL
      sloc_der = 0._CUSTOM_REAL
      allocate(seismograms_eps(NDIM,NDIM,nrec_local,NSTEP))
      seismograms_eps(:,:,:,:) = 0._CUSTOM_REAL
    endif
  endif

! initialize arrays to zero
  displ(:,:) = 0._CUSTOM_REAL
  veloc(:,:) = 0._CUSTOM_REAL
  accel(:,:) = 0._CUSTOM_REAL

! put negligible initial value to avoid very slow underflow trapping
  if(FIX_UNDERFLOW_PROBLEM) displ(:,:) = VERYSMALLVAL

!! DK DK array not created yet for CUBIT
! if (SIMULATION_TYPE == 3)  then ! kernel calculation, read in last frame

! open(unit=27,file=trim(prname)//'save_forward_arrays.bin',status='old',action='read',form='unformatted')
! read(27) b_displ
! read(27) b_veloc
! read(27) b_accel

! rho_kl(:,:,:,:) = 0._CUSTOM_REAL
! mu_kl(:,:,:,:) = 0._CUSTOM_REAL
! kappa_kl(:,:,:,:) = 0._CUSTOM_REAL

! endif

! allocate files to save movies and shaking map
  if(MOVIE_SURFACE .or. CREATE_SHAKEMAP) then
    if (USE_HIGHRES_FOR_MOVIES) then
      nmovie_points = NGLLX * NGLLY * NSPEC2D_TOP
    else
      nmovie_points = NGNOD2D * NSPEC2D_TOP
      iorderi(1) = 1
      iorderi(2) = NGLLX
      iorderi(3) = NGLLX
      iorderi(4) = 1
      iorderj(1) = 1
      iorderj(2) = 1
      iorderj(3) = NGLLY
      iorderj(4) = NGLLY
    endif
    allocate(store_val_x(nmovie_points))
    allocate(store_val_y(nmovie_points))
    allocate(store_val_z(nmovie_points))
    allocate(store_val_ux(nmovie_points))
    allocate(store_val_uy(nmovie_points))
    allocate(store_val_uz(nmovie_points))
    allocate(store_val_norm_displ(nmovie_points))
    allocate(store_val_norm_veloc(nmovie_points))
    allocate(store_val_norm_accel(nmovie_points))

    allocate(store_val_x_all(nmovie_points,0:NPROC-1))
    allocate(store_val_y_all(nmovie_points,0:NPROC-1))
    allocate(store_val_z_all(nmovie_points,0:NPROC-1))
    allocate(store_val_ux_all(nmovie_points,0:NPROC-1))
    allocate(store_val_uy_all(nmovie_points,0:NPROC-1))
    allocate(store_val_uz_all(nmovie_points,0:NPROC-1))

! to compute max of norm for shaking map
    store_val_norm_displ(:) = -1.
    store_val_norm_veloc(:) = -1.
    store_val_norm_accel(:) = -1.
  else if (MOVIE_VOLUME) then
    allocate(div(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(curl_x(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(curl_y(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(curl_z(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  endif

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '           time step: ',sngl(DT),' s'
    write(IMAIN,*) 'number of time steps: ',NSTEP
    write(IMAIN,*) 'total simulated time: ',sngl(NSTEP*DT),' seconds'
    write(IMAIN,*)
  endif

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
    deltat = sngl(DT)
  else
    deltat = DT
  endif
  deltatover2 = deltat/2.
  deltatsqover2 = deltat*deltat/2.
  if (SIMULATION_TYPE == 3) then
    if(CUSTOM_REAL == SIZE_REAL) then
      b_deltat = - sngl(DT)
    else
      b_deltat = - DT
    endif
    b_deltatover2 = b_deltat/2.
    b_deltatsqover2 = b_deltat*b_deltat/2.
  endif

! precompute Runge-Kutta coefficients if attenuation
  if(ATTENUATION) then
    tauinv(:,:) = - 1. / tau_sigma(:,:)
    factor_common(:,:) = 2. * beta(:,:) * tauinv(:,:)
    alphaval(:,:) = 1 + deltat*tauinv(:,:) + deltat**2*tauinv(:,:)**2 / 2. + &
      deltat**3*tauinv(:,:)**3 / 6. + deltat**4*tauinv(:,:)**4 / 24.
    betaval(:,:) = deltat / 2. + deltat**2*tauinv(:,:) / 3. + deltat**3*tauinv(:,:)**2 / 8. + deltat**4*tauinv(:,:)**3 / 24.
    gammaval(:,:) = deltat / 2. + deltat**2*tauinv(:,:) / 6. + deltat**3*tauinv(:,:)**2 / 24.
    !if (SIMULATION_TYPE == 3) then
    !  b_alphaval(:,:) = 1 + b_deltat*tauinv(:,:) + b_deltat**2*tauinv(:,:)**2 / 2. + &
    !        b_deltat**3*tauinv(:,:)**3 / 6. + b_deltat**4*tauinv(:,:)**4 / 24.
    !  b_betaval(:,:) = b_deltat / 2. + b_deltat**2*tauinv(:,:) / 3. + &
    !        b_deltat**3*tauinv(:,:)**2 / 8. + b_deltat**4*tauinv(:,:)**3 / 24.
    !  b_gammaval(:,:) = b_deltat / 2. + b_deltat**2*tauinv(:,:) / 6. + &
    !        b_deltat**3*tauinv(:,:)**2 / 24.
    !endif
  endif


  !pll, to put elsewhere
  ! note: currently, they need to be defined here, as they are used in the routine arguments 
  !          for compute_forces_with_Deville()
  allocate(R_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS))
  allocate(R_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS))
  allocate(R_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS))
  allocate(R_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS))
  allocate(R_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS))
  allocate(epsilondev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB))
  allocate(epsilondev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB))
  allocate(epsilondev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB))
  allocate(epsilondev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB))
  allocate(epsilondev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB))

! clear memory variables if attenuation
  if(ATTENUATION) then
  
   ! initialize memory variables for attenuation
    epsilondev_xx(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_yy(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_xy(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_xz(:,:,:,:) = 0._CUSTOM_REAL
    epsilondev_yz(:,:,:,:) = 0._CUSTOM_REAL

    R_xx(:,:,:,:,:) = 0._CUSTOM_REAL
    R_yy(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xy(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xz(:,:,:,:,:) = 0._CUSTOM_REAL
    R_yz(:,:,:,:,:) = 0._CUSTOM_REAL

    if(FIX_UNDERFLOW_PROBLEM) then
      R_xx(:,:,:,:,:) = VERYSMALLVAL
      R_yy(:,:,:,:,:) = VERYSMALLVAL
      R_xy(:,:,:,:,:) = VERYSMALLVAL
      R_xz(:,:,:,:,:) = VERYSMALLVAL
      R_yz(:,:,:,:,:) = VERYSMALLVAL
    endif

!! DK DK array not created yet for CUBIT
!   if (SIMULATION_TYPE == 3) then
!     read(27) b_R_xx
!     read(27) b_R_yy
!     read(27) b_R_xy
!     read(27) b_R_xz
!     read(27) b_R_yz
!     read(27) b_epsilondev_xx
!     read(27) b_epsilondev_yy
!     read(27) b_epsilondev_xy
!     read(27) b_epsilondev_xz
!     read(27) b_epsilondev_yz
!   endif

  endif
  close(27)

! initialize Moho boundary index
! if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
!   ispec2D_moho_top = 0
!   ispec2D_moho_bot = 0
!   k_top = 1
!   k_bot = NGLLZ
! endif



  end subroutine