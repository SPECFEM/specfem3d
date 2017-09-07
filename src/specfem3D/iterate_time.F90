!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine iterate_time()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use gravity_perturbation, only: gravity_timeseries, GRAVITY_SIMULATION

  implicit none

! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
! local parameters
  integer :: ier
! *********************************************************************************

  ! for EXACT_UNDOING_TO_DISK
  integer :: ispec,iglob,i,j,k,counter,record_length
  integer, dimension(:), allocatable :: integer_mask_ibool_exact_undo
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: buffer_for_disk
  character(len=MAX_STRING_LEN) outputname
  ! timing
  double precision, external :: wtime

  !! for FK point for intialization injected wavefield
  real(kind=CUSTOM_REAL) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box

  !----  create a Gnuplot script to display the energy curve in log scale
  if (OUTPUT_ENERGY .and. myrank == 0) then
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set terminal x11'
    write(IOUT_ENERGY,*) '#set terminal wxt'
    write(IOUT_ENERGY,*) '#set terminal postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time step number"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,'(a152)') 'plot "energy.dat" us 1:2 t "Kinetic Energy" w l lc 1, "energy.dat" us 1:3 &
                         &t "Potential Energy" w l lc 2, "energy.dat" us 1:4 t "Total Energy" w l lc 4'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:2 t "Kinetic Energy" w l lc 1'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:3 t "Potential Energy" w l lc 2'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:4 t "Total Energy" w l lc 4'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

 !! CD CD adds this (temporary)
  if (RECIPROCITY_AND_KH_INTEGRAL) open(unit=158,file='KH_integral',status='unknown')

  ! open the file in which we will store the energy curve
  if (OUTPUT_ENERGY .and. myrank == 0) &
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'energy.dat',status='unknown',action='write')


! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation

! ********************************************
!  initial setup for future FK3D calculations
! ********************************************

  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK .and. SIMULATION_TYPE == 1) then

    call nbound(NSPEC_AB,num_abs_boundary_faces,abs_boundary_ispec,ispec_is_elastic,npt)

    !! compute the bottom midle point of the domain


    !! VM VM dealocate in case of severals runs occurs in inverse_problem program
    if (allocated(nbdglb)) deallocate(nbdglb)
    if (allocated(vx_FK))  deallocate(vx_FK)
    if (allocated(vy_FK))  deallocate(vy_FK)
    if (allocated(vz_FK))  deallocate(vz_FK)
    if (allocated(tx_FK))  deallocate(tx_FK)
    if (allocated(ty_FK))  deallocate(ty_FK)
    if (allocated(tz_FK))  deallocate(tz_FK)
    if (allocated(VX_t))   deallocate(VX_t)
    if (allocated(VY_t))   deallocate(VY_t)
    if (allocated(VZ_t))   deallocate(VZ_t)
    if (allocated(TX_t))   deallocate(TX_t)
    if (allocated(TY_t))   deallocate(TY_t)
    if (allocated(TZ_t))   deallocate(TZ_t)

    !! allocate memory for FK solution
    if (npt > 0) then

       allocate(nbdglb(npt))
       allocate(vx_FK(npt),vy_FK(npt),vz_FK(npt),tx_FK(npt),ty_FK(npt),tz_FK(npt))

    else

       allocate(nbdglb(1))
       allocate(vx_FK(1),vy_FK(1),vz_FK(1),tx_FK(1),ty_FK(1),tz_FK(1))

    endif

    call FindBoundaryBox(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)

    if (myrank == 0) then

       call ReadFKModelInput(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box)

    endif

    ! send FK parameters to others MPI slices

    call bcast_all_singlei(kpsv)
    call bcast_all_singlei(nlayer)
    if (myrank > 0) allocate(al_FK(nlayer),be_FK(nlayer),mu_FK(nlayer),h_FK(nlayer))
    call bcast_all_cr(al_FK, nlayer)
    call bcast_all_cr(be_FK, nlayer)
    call bcast_all_cr(mu_FK, nlayer)
    call bcast_all_cr(h_FK, nlayer)
    call bcast_all_singlecr(phi_FK)
    call bcast_all_singlecr(theta_FK)
    call bcast_all_singlecr(ff0)
    call bcast_all_singlecr(xx0)
    call bcast_all_singlecr(yy0)
    call bcast_all_singlecr(zz0)
    call bcast_all_singlecr(tt0)
    call bcast_all_singlecr(Z_REF_for_FK)
    call bcast_all_singlecr(tmax_fk)
    call bcast_all_singlel(stag)


    phi_FK   = phi_FK*acos(-1.d0)/180.d0
    theta_FK = theta_FK*acos(-1.d0)/180.d0

    if (kpsv == 1) then
       p = sin(theta_FK)/al_FK(nlayer)
    else if (kpsv == 2) then
       p = sin(theta_FK)/be_FK(nlayer)
    endif

    tg  = 1.d0/ff0

    !-------------------------------------------------------------
    ! get MPI starting time for FK
    time_start = wtime()

    if (npt > 0) then

       call find_size_of_working_arrays(deltat, tmax_fk, NF_FOR_STORING, NF_FOR_FFT, NPOW_FOR_INTERP, NP_RESAMP, DF_FK)
       if (myrank == 0) write(IMAIN,*) 'ALLOCATE SAVED ARRAYS ', NF_FOR_STORING, NF_FOR_FFT, npt

       !! arrays for storing FK solution --------------------------------------------

       allocate(VX_t(npt,  -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
       if (ier /= 0) stop 'error while allocating VX_t'
       VX_t(:,:)=0._CUSTOM_REAL

       allocate(VY_t(npt,  -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
       if (ier /= 0) stop 'error while allocating VY_t'
       VY_t(:,:)=0._CUSTOM_REAL

       allocate(VZ_t(npt,  -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
       if (ier /= 0) stop 'error while allocating VZ_t'
       VZ_t(:,:)=0._CUSTOM_REAL

       allocate(TX_t(npt,  -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
       if (ier /= 0) stop 'error while allocating TX_t'
       TX_t(:,:)=0._CUSTOM_REAL

       allocate(TY_t(npt,  -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
       if (ier /= 0) stop 'error while allocating TY_t'
       TY_t(:,:)=0._CUSTOM_REAL

       allocate(TZ_t(npt, -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
       if (ier /= 0) stop 'error while allocating TZ_t'
       TZ_t(:,:)=0._CUSTOM_REAL



       call FK3D(myrank, NSPEC_AB, ibool, abs_boundary_ijk, abs_boundary_normal, &
            abs_boundary_ispec, num_abs_boundary_faces, ispec_is_elastic, kpsv, nlayer, nstep, npt, nbdglb, &
            p, phi_FK, xx0, yy0, zz0, tg, &
            tt0, al_FK, be_FK, mu_FK, h_FK, deltat, &
            NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, DF_FK)

    endif

    call synchronize_all()

    !-------------------------------------------------------------
    ! get MPI ending ting time for FK
    time_start = wtime() - time_start
    if (myrank == 0) then
       write(IMAIN, *)
       write(IMAIN, '(a35,1x, f20.2, a7)')  " Elapsed time for FK computation : ",  time_start, " sec. "
       write(IMAIN, *)
       write(IMAIN,*) " ********************************************** "
       call flush_IMAIN()
    endif

    deallocate(al_FK, be_FK, mu_FK, h_FK)

 endif

! *****************************************************
! * end of initial setup for future FK3D calculations *
! *****************************************************
! *********************************************************************************

!
!   s t a r t   t i m e   i t e r a t i o n s
!

  ! synchronize all processes to make sure everybody is ready to start time loop
  call synchronize_all()
  if (myrank == 0) write(IMAIN,*) 'All processes are synchronized before time loop'

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! create an empty file to monitor the start of the simulation
  if (myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop.txt',status='unknown',action='write')
    write(IOUT,*) 'hello, starting time loop'
    close(IOUT)
  endif

  ! initialize variables for writing seismograms
  seismo_offset = it_begin-1
  seismo_current = 0

  ! *********************************************************************************

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  if (EXACT_UNDOING_TO_DISK) then

    if (GPU_MODE) call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK not supported for GPUs')

    if (UNDO_ATTENUATION) &
      call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK needs UNDO_ATTENUATION to be off because it computes the kernel directly instead')

    if (SIMULATION_TYPE == 1 .and. .not. SAVE_FORWARD) &
      call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK requires SAVE_FORWARD if SIMULATION_TYPE == 1')

    if (ANISOTROPIC_KL) call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK requires ANISOTROPIC_KL to be turned off')

!! DK DK determine the largest value of iglob that we need to save to disk,
!! DK DK since we save the upper part of the mesh only in the case of surface-wave kernels
    ! crust_mantle
    allocate(integer_mask_ibool_exact_undo(NGLOB_AB))
    integer_mask_ibool_exact_undo(:) = -1

    counter = 0
    do ispec = 1, NSPEC_AB
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)
!           height = xstore(iglob)
            ! save that element only if it is in the upper part of the mesh
!           if (height >= 3000.d0) then
            if (.true.) then
              ! if this point has not yet been found before
              if (integer_mask_ibool_exact_undo(iglob) == -1) then
                ! create a new unique point
                counter = counter + 1
                integer_mask_ibool_exact_undo(iglob) = counter
              endif
            endif
          enddo
        enddo
      enddo
    enddo

    ! allocate the buffer used to dump a single time step
    allocate(buffer_for_disk(counter))

    ! open the file in which we will dump all the time steps (in a single file)
    write(outputname,"('huge_dumps/proc',i6.6,'_huge_dump_of_all_time_steps.bin')") myrank
    inquire(iolength=record_length) buffer_for_disk
    ! we write to or read from the file depending on the simulation type
    if (SIMULATION_TYPE == 1) then
      open(file=outputname, unit=IFILE_FOR_EXACT_UNDOING, action='write', status='unknown', &
                      form='unformatted', access='direct', recl=record_length)
    else if (SIMULATION_TYPE == 3) then
      open(file=outputname, unit=IFILE_FOR_EXACT_UNDOING, action='read', status='old', &
                      form='unformatted', access='direct', recl=record_length)
    else
      call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK can only be used with SIMULATION_TYPE == 1 or SIMULATION_TYPE == 3')
    endif

  endif ! of if (EXACT_UNDOING_TO_DISK)

  ! get MPI starting
  time_start = wtime()

  it_begin = 1
  it_end = NSTEP
  do it = it_begin,it_end

    ! simulation status output and stability check
    if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
      call check_stability()
    endif

    ! simulation status output and stability check
    if (OUTPUT_ENERGY) then
      if (mod(it,NTSTEP_BETWEEN_OUTPUT_ENERGY) == 0 .or. it == 5 .or. it == NSTEP) call compute_energy()
    endif

    ! updates wavefields using Newmark time scheme
    if (.not. USE_LDDRK) call update_displacement_scheme()

    ! calculates stiffness term
    if (.not. GPU_MODE) then
      ! wavefields on CPU

      ! note: the order of the computations for acoustic and elastic domains is crucial for coupled simulations
      if (SIMULATION_TYPE == 3) then
        ! kernel/adjoint simulations

        ! adjoint wavefields
        if (ELASTIC_SIMULATION .and. ACOUSTIC_SIMULATION) then
          ! coupled acoustic-elastic simulations
          ! 1. elastic domain w/ adjoint wavefields
          call compute_forces_viscoelastic_calling()
          ! 2. acoustic domain w/ adjoint wavefields
          call compute_forces_acoustic_calling()
        else
          ! non-coupled simulations
          ! (purely acoustic or elastic)
          if (ACOUSTIC_SIMULATION) call compute_forces_acoustic_calling()
          if (ELASTIC_SIMULATION) call compute_forces_viscoelastic_calling()
        endif

        ! backward/reconstructed wavefields
        ! acoustic solver
        ! (needs to be done after elastic one)
        if (ACOUSTIC_SIMULATION) call compute_forces_acoustic_backward_calling()
        ! elastic solver
        ! (needs to be done first, before poroelastic one)
        if (ELASTIC_SIMULATION) call compute_forces_viscoelastic_backward_calling()

      else
        ! forward simulations
        do istage = 1, NSTAGE_TIME_SCHEME
          if (USE_LDDRK) call update_displ_lddrk()
          ! 1. acoustic domain
          if (ACOUSTIC_SIMULATION) call compute_forces_acoustic_calling()
          ! 2. elastic domain
          if (ELASTIC_SIMULATION) call compute_forces_viscoelastic_calling()
        enddo
      endif

      ! poroelastic solver
      if (POROELASTIC_SIMULATION) call compute_forces_poroelastic_calling()

    else
      ! wavefields on GPU
      ! acoustic solver
      if (ACOUSTIC_SIMULATION) call compute_forces_acoustic_GPU_calling()
      ! elastic solver
      ! (needs to be done first, before poroelastic one)
      if (ELASTIC_SIMULATION) call compute_forces_viscoelastic_GPU_calling()
    endif

    ! restores last time snapshot saved for backward/reconstruction of wavefields
    ! note: this must be read in after the Newmark time scheme
    if (SIMULATION_TYPE == 3 .and. it == 1) then
      call it_read_forward_arrays()
    endif

    ! calculating gravity field at current timestep
    if (GRAVITY_SIMULATION) call gravity_timeseries()

    ! write the seismograms with time shift (GPU_MODE transfer included)
    call write_seismograms()

    ! adjoint simulations: kernels
    if (SIMULATION_TYPE == 3) then
      call compute_kernels()
    endif

    ! outputs movie files
    call write_movie_output()

    ! first step of noise tomography, i.e., save a surface movie at every time step
    ! modified from the subroutine 'write_movie_surface'
    if (NOISE_TOMOGRAPHY == 1) then
      call noise_save_surface_movie(displ,ibool,noise_surface_movie,it,NSPEC_AB,NGLOB_AB, &
                            num_free_surface_faces,free_surface_ispec,free_surface_ijk,Mesh_pointer,GPU_MODE)
    endif

    ! updates VTK window
    if (VTK_MODE) then
      call it_update_vtkwindow()
    endif

    !! CD CD add this : under validation option
    if (RECIPROCITY_AND_KH_INTEGRAL) then
      if (.not. SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
        call surface_or_volume_integral_on_whole_domain()
        write(158,*) it*DT, integral_boun(1), integral_boun(2), integral_boun(3)
      endif
    endif

  !
  !---- end of time iteration loop
  !
  enddo   ! end of main time loop


  ! close the huge file that contains a dump of all the time steps to disk
  if (EXACT_UNDOING_TO_DISK) close(IFILE_FOR_EXACT_UNDOING)

  call it_print_elapsed_time()

  !! CD CD added this
  if (RECIPROCITY_AND_KH_INTEGRAL) then
    close(158)
    close(237)
    close(238)
  endif

  ! Transfer fields from GPU card to host for further analysis
  if (GPU_MODE) call it_transfer_from_GPU()

!----  close energy file
  if (OUTPUT_ENERGY .and. myrank == 0) close(IOUT_ENERGY)

  end subroutine iterate_time

!
!-------------------------------------------------------------------------------------------------
!

  subroutine surface_or_volume_integral_on_whole_domain()

  use constants

  use specfem_par

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE, INJECTION_TECHNIQUE_TYPE

  implicit none

  ! local parameters
  double precision :: weightpt, jacobianpt

  double precision, dimension(3) :: integr_volloc, integr_bounloc
  double precision, dimension(3,NGLLSQUARE*num_abs_boundary_faces) :: f_integrand_bounloc
  double precision, dimension(3,NGLOB_AB) :: f_integrand_volloc

  real(kind=CUSTOM_REAL) :: xixpt,xiypt,xizpt,etaxpt,etaypt,etazpt,gammaxpt,gammaypt,gammazpt

  integer :: ier,iint,i,j,k,ispec,iglob,igll,iface

!
!-----------------------------------------------------------------------------
!

  ! NOTE : 'f_integrandloc' have to be defined at all 'iglob'
  ! (all the points of the surface, or all the point of the volume)

  if (RECIPROCITY_AND_KH_INTEGRAL) then

    allocate(f_integrand_KH(3,NGLLSQUARE*num_abs_boundary_faces), stat=ier)

    call integrand_for_computing_Kirchoff_Helmholtz_integral()

    do iint = 1,3
      f_integrand_bounloc(iint,:) = f_integrand_KH(iint,:)
    enddo

  endif

  do iint = 1,3
    integr_volloc(iint)  = 0.d0
    integr_bounloc(iint) = 0.d0
  enddo

  if ( (Surf_or_vol_integral == 2) .or. (Surf_or_vol_integral == 3) ) then

    allocate(integral_vol(3), stat=ier)

    ! calculates volume of all elements in mesh
    do ispec = 1, NSPEC_AB

      ! main computation
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX

            weightpt = wxgll(i)*wygll(j)*wzgll(k)

            ! compute the Jacobian
            xixpt    = xix(i,j,k,ispec)
            xiypt    = xiy(i,j,k,ispec)
            xizpt    = xiz(i,j,k,ispec)
            etaxpt   = etax(i,j,k,ispec)
            etaypt   = etay(i,j,k,ispec)
            etazpt   = etaz(i,j,k,ispec)
            gammaxpt = gammax(i,j,k,ispec)
            gammaypt = gammay(i,j,k,ispec)
            gammazpt = gammaz(i,j,k,ispec)

            ! do this in double precision for accuracy
            jacobianpt = 1.d0 / dble(xixpt*(etaypt*gammazpt-etazpt*gammaypt) &
                                    - xiypt*(etaxpt*gammazpt-etazpt*gammaxpt) &
                                    + xizpt*(etaxpt*gammaypt-etaypt*gammaxpt))

            if (CHECK_FOR_NEGATIVE_JACOBIANS .and. jacobianpt <= ZERO) stop &
                                             'error: negative Jacobian found in volume integral calculation'

            iglob = ibool(i,j,k,ispec)

            integr_volloc(:)  = integr_volloc(:) + ( weightpt * jacobianpt * f_integrand_volloc(:,iglob) )

          enddo
        enddo
      enddo
    enddo

    call sum_all_1Darray_dp(integr_volloc, integral_vol, 3)

  endif

!
!---
!

  if ( (Surf_or_vol_integral == 1) .or. (Surf_or_vol_integral == 3) ) then

    allocate(integral_boun(3), stat=ier)

    ! calculates integral on all the surface of the whole domain

    if (num_free_surface_faces > 0) then

      do iface = 1, num_free_surface_faces

        ispec = free_surface_ispec(iface)

        do igll = 1, NGLLSQUARE

          ! gets local indices for GLL point
          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)

          iglob = ibool(i,j,k,ispec)

!!          integr_bounloc(:) = integr_bounloc(:) + (free_surface_jacobian2Dw(igll,iface) * f_integrand_bounloc(:,iglob))
          integr_bounloc(:) = integr_bounloc(:) + (free_surface_jacobian2Dw(igll,iface) * f_integrand_bounloc(:,igll*iface))

        enddo
      enddo
    endif

    if (num_abs_boundary_faces > 0) then

      do iface = 1,num_abs_boundary_faces

        ispec = abs_boundary_ispec(iface)

        do igll = 1,NGLLSQUARE

          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          iglob = ibool(i,j,k,ispec)

!!          integr_bounloc(:) = integr_bounloc(:) + (abs_boundary_jacobian2Dw(igll,iface) * f_integrand_bounloc(:,iglob))
          integr_bounloc(:) = integr_bounloc(:) + (abs_boundary_jacobian2Dw(igll,iface) * f_integrand_bounloc(:,igll*iface))

        enddo
      enddo
    endif

    call sum_all_1Darray_dp(integr_bounloc, integral_boun, 3)

  endif

  if (RECIPROCITY_AND_KH_INTEGRAL) deallocate(f_integrand_KH)

  end subroutine surface_or_volume_integral_on_whole_domain

!
!-------------------------------------------------------------------------------------------------
!
!--- CD CD : subroutine not validated yet
!
!
  subroutine integrand_for_computing_Kirchoff_Helmholtz_integral()

  use constants

  use specfem_par

  use shared_parameters

  implicit none

  ! local parameters
  integer :: itau

!!  integer :: i, iglob, iint

  double precision, dimension(3,NGLLSQUARE*num_abs_boundary_faces) :: convol_displsem_tractaxisem
  double precision, dimension(3,NGLLSQUARE*num_abs_boundary_faces) :: convol_tractsem_displaxisem

  convol_displsem_tractaxisem(:,:) = 0.d0
  convol_tractsem_displaxisem(:,:) = 0.d0

!---
!---

  if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then

    if (old_DSM_coupling_from_Vadim) then

    else
      !! for 2D light version
    endif

  else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then

    call read_axisem_disp_file(num_abs_boundary_faces*NGLLSQUARE)
    call read_specfem_tract_file(NGLLSQUARE*num_abs_boundary_faces)
    call read_specfem_disp_file(NGLLSQUARE*num_abs_boundary_faces)

!!!-- ci-dessous : en cours de modification

    do itau = 1, it
      if ( it > itau .and. (it - itau <= NSTEP) ) then

        convol_displsem_tractaxisem(:,:) = convol_displsem_tractaxisem(:,:) + &
                                           Displ_specfem_time(:,:, itau) * Tract_axisem_time(:,:, it - itau) * dt

        convol_tractsem_displaxisem(:,:) = convol_tractsem_displaxisem(:,:) + &
                                           Tract_specfem_time(:,:, itau) * Displ_axisem_time(:,:, it - itau) * dt

      endif
    enddo

    f_integrand_KH(:,:) = convol_displsem_tractaxisem(:,:) - convol_tractsem_displaxisem(:,:)

  endif

!---
!---

  end subroutine integrand_for_computing_Kirchoff_Helmholtz_integral

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_axisem_disp_file(nb)

  use constants

  use specfem_par, only: it, Displ_axisem_time, RECIPROCITY_AND_KH_INTEGRAL, IIN_displ_axisem

  implicit none

  ! argument
  integer nb

  ! local parameters
  real(kind=CUSTOM_REAL) :: Displ_axisem(3,nb)

  read(IIN_displ_axisem) Displ_axisem

  if (RECIPROCITY_AND_KH_INTEGRAL) Displ_axisem_time(:,:,it) = Displ_axisem(:,:)

  end subroutine read_axisem_disp_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_specfem_disp_file(nb)

  use constants

  use specfem_par, only: it, Displ_specfem_time, RECIPROCITY_AND_KH_INTEGRAL, num_abs_boundary_faces

  implicit none

  ! argument
  integer nb

  ! local parameters
  integer iface, igll
  real(kind=CUSTOM_REAL) :: Displ_specfem(3,nb)

  do iface = 1,num_abs_boundary_faces
    do igll = 1,NGLLSQUARE

      read(238) Displ_specfem(1,igll*iface), Displ_specfem(3,igll*iface), Displ_specfem(3,igll*iface)

    enddo
  enddo

  if (RECIPROCITY_AND_KH_INTEGRAL) Displ_specfem_time(:,:,it) = Displ_specfem(:,:)

  end subroutine read_specfem_disp_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_specfem_tract_file(nb)

  use constants

  use specfem_par, only: it, Tract_specfem_time, RECIPROCITY_AND_KH_INTEGRAL, num_abs_boundary_faces

  implicit none

  ! argument
  integer nb

  ! local parameters
  integer iface, igll
  real(kind=CUSTOM_REAL) :: Tract_specfem(3,nb)

  do iface = 1,num_abs_boundary_faces
    do igll = 1,NGLLSQUARE
      read(237) Tract_specfem(1,igll*iface), Tract_specfem(3,igll*iface), Tract_specfem(3,igll*iface)
    enddo
  enddo

  if (RECIPROCITY_AND_KH_INTEGRAL) Tract_specfem_time(:,:,it) = Tract_specfem(:,:)

  end subroutine read_specfem_tract_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_transfer_from_GPU()

! transfers fields on GPU back onto CPU

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  ! to store forward wave fields
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then

    ! acoustic potentials
    if (ACOUSTIC_SIMULATION) &
      call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                                          potential_dot_acoustic, potential_dot_dot_acoustic, &
                                          Mesh_pointer)

    ! elastic wavefield
    if (ELASTIC_SIMULATION) then
      call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc,accel,Mesh_pointer)

      if (ATTENUATION) &
        call transfer_fields_att_from_device(Mesh_pointer, &
                    R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
                    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                    size(epsilondev_xx))

    endif

  else if (SIMULATION_TYPE == 3) then

    ! to store kernels
    ! acoustic domains
    if (ACOUSTIC_SIMULATION) then
      ! only in case needed...
      !call transfer_b_fields_ac_from_device(NGLOB_AB,b_potential_acoustic, &
      !                      b_potential_dot_acoustic, b_potential_dot_dot_acoustic, Mesh_pointer)

      ! acoustic kernels
      call transfer_kernels_ac_to_host(Mesh_pointer,rho_ac_kl,kappa_ac_kl,NSPEC_AB)
    endif

    ! elastic domains
    if (ELASTIC_SIMULATION) then
      ! only in case needed...
      !call transfer_b_fields_from_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)

      ! elastic kernels
      call transfer_kernels_el_to_host(Mesh_pointer,rho_kl,mu_kl,kappa_kl,cijkl_kl,NSPEC_AB)
    endif

    ! specific noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      call transfer_kernels_noise_to_host(Mesh_pointer,sigma_kl,NSPEC_AB)
    endif

    ! approximative Hessian for preconditioning kernels
    if (APPROXIMATE_HESS_KL) then
      if (ELASTIC_SIMULATION) &
           call transfer_kernels_hess_el_tohost(Mesh_pointer,hess_kl,hess_rho_kl,hess_kappa_kl,hess_mu_kl,NSPEC_AB)
      if (ACOUSTIC_SIMULATION) &
           call transfer_kernels_hess_ac_tohost(Mesh_pointer,hess_ac_kl, hess_rho_ac_kl,hess_kappa_ac_kl,NSPEC_AB)
    endif

  endif

  ! from here on, no gpu data is needed anymore
  ! frees allocated memory on GPU
  call prepare_cleanup_device(Mesh_pointer,ACOUSTIC_SIMULATION,ELASTIC_SIMULATION, &
                              STACEY_ABSORBING_CONDITIONS,NOISE_TOMOGRAPHY,COMPUTE_AND_STORE_STRAIN, &
                              ATTENUATION,ANISOTROPY,APPROXIMATE_OCEAN_LOAD, &
                              APPROXIMATE_HESS_KL)

  end subroutine it_transfer_from_GPU

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_update_vtkwindow()

  implicit none

  stop 'it_update_vtkwindow() not implemented for this code yet'

  end subroutine it_update_vtkwindow

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_read_forward_arrays()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  integer :: ier

! restores last time snapshot saved for backward/reconstruction of wavefields
! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
!          and adjoint sources will become more complicated
!          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields
  if (ADIOS_FOR_FORWARD_ARRAYS) then
    call read_forward_arrays_adios()
  else
    ! reads in wavefields
    open(unit=IIN,file=trim(prname)//'save_forward_arrays.bin',status='old', &
          action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'error: opening save_forward_arrays'
      print *,'path: ',trim(prname)//'save_forward_arrays.bin'
      call exit_mpi(myrank,'error open file save_forward_arrays.bin')
    endif

    if (ACOUSTIC_SIMULATION) then
      read(IIN) b_potential_acoustic
      read(IIN) b_potential_dot_acoustic
      read(IIN) b_potential_dot_dot_acoustic
    endif

    ! elastic wavefields
    if (ELASTIC_SIMULATION) then
      read(IIN) b_displ
      read(IIN) b_veloc
      read(IIN) b_accel
      ! memory variables if attenuation
      if (ATTENUATION) then
        read(IIN) b_R_trace
        read(IIN) b_R_xx
        read(IIN) b_R_yy
        read(IIN) b_R_xy
        read(IIN) b_R_xz
        read(IIN) b_R_yz
        read(IIN) b_epsilondev_trace
        read(IIN) b_epsilondev_xx
        read(IIN) b_epsilondev_yy
        read(IIN) b_epsilondev_xy
        read(IIN) b_epsilondev_xz
        read(IIN) b_epsilondev_yz
      endif
    endif

    ! poroelastic wavefields
    if (POROELASTIC_SIMULATION) then
      read(IIN) b_displs_poroelastic
      read(IIN) b_velocs_poroelastic
      read(IIN) b_accels_poroelastic
      read(IIN) b_displw_poroelastic
      read(IIN) b_velocw_poroelastic
      read(IIN) b_accelw_poroelastic
    endif

    close(IIN)
  endif

  if (GPU_MODE) then
    if (ACOUSTIC_SIMULATION) then
    ! transfers fields onto GPU
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic, &
                                          b_potential_dot_dot_acoustic, &
                                          Mesh_pointer)
    endif
    ! elastic wavefields
    if (ELASTIC_SIMULATION) then
      ! puts elastic wavefield to GPU
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)
      ! memory variables if attenuation
      if (ATTENUATION) then
        call transfer_b_fields_att_to_device(Mesh_pointer, &
                           b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz, &
                           size(b_R_xx), &
                           b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy, &
                           b_epsilondev_xz,b_epsilondev_yz, &
                           size(b_epsilondev_xx))
      endif
    endif
  endif

  end subroutine it_read_forward_arrays

  !!==============================================================================================================================
  !! VM VM READING INPUT FILE FOR FK MODEL ----------------------------------------------------------------------
  subroutine ReadFKModelInput(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box)
    use specfem_par
    use specfem_par_elastic

    real(kind=CUSTOM_REAL), intent(in) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box
    integer                :: ioerr
    character(len=100)     :: keyword, keyvalue, line
    character(len=100)     :: keyword_tmp, incident_wave

    real(kind=CUSTOM_REAL) :: rho_layer, vp_layer, vs_layer, ztop_layer
    real(kind=CUSTOM_REAL) :: Radius_box, wave_length_at_bottom
    real(kind=CUSTOM_REAL), dimension(:), allocatable  :: rho_fk_input, vp_fk_input, vs_fk_input, ztop_fk_input
    integer,  dimension(:), allocatable  :: ilayer_fk_input
    integer  :: ilayer
    logical  :: position_of_wavefront_not_read

    !!--------------------------------------------------------------
    ! # model description :
    ! NLAYER   n # number of layers
    ! LAYER 1  rho, vp ,vs, ztop
    ! LAYER 2  rho, vp, vs, ztop
    ! ...
    ! LAYER n  rho, vp, vs, ztop  # homogenoeus half-space
    ! #
    ! # incident wave description:
    ! INCIDENT_WAVE  "p" or "sv"
    ! BACK_AZITUTH    bazi
    ! INCIDENCE       inc
    ! ORIGIN_WAVEFRONT xx0, yy0, zz0
    ! ORIGIN_TIME      tt0
    ! FREQUENCY_MAX    ff0
    ! TIME_WINDOW      tmax_fk
    !!----------------------------------------------------------------

    !! set default values
    tt0=0.
    tmax_fk=128.
    ff0=0.1
    kpsv=1
    position_of_wavefront_not_read=.true.
    stag=.false.

    !! READING input file
    open(85,file=trim(FKMODEL_FILE))
    do
       read(85, fmt='(a100)',iostat=ioerr) line
       !call remove_begin_blanks(line)
       if (ioerr < 0) exit
       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
       read(line,*) keyword, keyvalue
       select case(trim(keyword))

       case('NLAYER')
          read(line, *) keyword_tmp, nlayer
          allocate(al_FK(nlayer), be_FK(nlayer), mu_FK(nlayer), h_FK(nlayer))
          allocate(rho_fk_input(nlayer), vp_fk_input(nlayer), vs_fk_input(nlayer), ztop_fk_input(nlayer+1), &
               ilayer_fk_input(nlayer+1))
          ilayer_fk_input(:)=-1

       case('LAYER')
          read(line, *) keyword_tmp, ilayer, rho_layer, vp_layer, vs_layer, ztop_layer
          ilayer_fk_input(ilayer)=ilayer
          rho_fk_input(ilayer)=rho_layer
          vp_fk_input(ilayer)=vp_layer
          vs_fk_input(ilayer)=vs_layer
          ztop_fk_input(ilayer)=ztop_layer

       case('INCIDENT_WAVE')
           read(line,*)  keyword_tmp, incident_wave

           select case(trim(incident_wave))
              case ('p', 'P')
                 kpsv=1
              case('sv','SV')
                 kpsv=2
              case default
                 kpsv=1
              end select

       case('BACK_AZIMUTH')
          read(line,*)  keyword_tmp, phi_FK
          phi_FK = - phi_FK - 90.

       case('AZIMUTH')
          read(line,*)  keyword_tmp, phi_FK
          phi_FK = 90. - phi_FK

       case('TAKE_OFF')
          read(line,*)  keyword_tmp, theta_FK

       case('ORIGIN_WAVEFRONT')
           read(line,*)  keyword_tmp, xx0, yy0, zz0
           position_of_wavefront_not_read=.false.

        case('ORIGIN_TIME')
           read(line,*)  keyword_tmp, tt0

       case('FREQUENCY_MAX')
          read(line,*)  keyword_tmp, ff0

       case('TIME_WINDOW')
          read(line,*)  keyword_tmp, tmax_fk

       end select
    !!------------------------------------------------------------------------------------------------------
    enddo

    if (allocated(ilayer_fk_input)) then



       ilayer_fk_input(nlayer+1) = ilayer_fk_input(nlayer)
       ztop_fk_input(nlayer+1)=ztop_fk_input(nlayer)
       Z_ref_for_FK=ztop_fk_input(nlayer)
       do ilayer=1, nlayer
          al_FK(ilayer) = vp_fk_input(ilayer)
          be_FK(ilayer) = vs_fk_input(ilayer)
          mu_FK(ilayer) = rho_fk_input(ilayer) * vs_fk_input(ilayer)**2
          h_FK(ilayer) =  ztop_fk_input(ilayer) - ztop_fk_input(ilayer+1)

          if (ilayer_fk_input(ilayer) == -1) then
             write(*,*) " ERROR READING FK INPUT FILE "
             write(*,*) " MISSING LAYER ", ilayer
             stop
          endif
       enddo

       deallocate(ilayer_fk_input, rho_fk_input, vp_fk_input, vs_fk_input, ztop_fk_input)

    else

       write(*,*) " ERROR READING FK INPUT FILE "
       write(*,*) " NOT BE ABLE TO READ MODEL PROPERTIES "
       stop

    endif

    !! compute position of wave front
    if (position_of_wavefront_not_read) then
       xx0=0.5*(Xmin_box + Xmax_box)
       yy0=0.5*(Ymin_box + Ymax_box)
       Radius_box = sqrt( (Xmin_box - xx0)**2 + (Ymin_box - yy0)**2)

       if (kpsv == 1) then
          wave_length_at_bottom = al_FK(nlayer) / ff0
       else if (kpsv == 2) then
          wave_length_at_bottom = be_FK(nlayer) / ff0
       endif

       zz0 = Zmin_box - Radius_box * sin ( abs (theta_FK) * (acos(-1.d0) / 180.d0)  ) -  &
             wave_length_at_bottom * cos ( abs (theta_FK) * (acos(-1.d0) / 180.d0)  ) -  &
             Z_ref_for_FK

    endif

    write(IMAIN,*) " ********************************************** "
    write(IMAIN,*) "         USING FK INJECTION TECHNIQUE           "
    write(IMAIN,*) " ********************************************** "

    write(IMAIN,*)
    write(IMAIN,*) "         Model : " , nlayer , " layers "
    write(IMAIN,*)

    do ilayer =1, nlayer
       write(IMAIN,'(a7, i3, 3(a6,2x,f8.3), 3x, a9, f18.5 )') &
            'layer ' , ilayer, &
            " rho  =",   mu_FK(ilayer) /  be_FK(ilayer)**2, &
            " vp   =",   al_FK(ilayer), &
            " vs   =",   be_FK(ilayer), &
            " Height =", h_FK(ilayer)
    enddo

    write(IMAIN,*)
    write(IMAIN,*)
    write(IMAIN,*)  " Phi FK :",   phi_FK
    write(IMAIN,*)  " theta FK :", theta_FK
    write(IMAIN,*)

    write(IMAIN,*)
    write(IMAIN,*) " Origin wavefront point FK :", xx0, yy0, zz0
    write(IMAIN,*) " time shift  FK :", tt0
    write(IMAIN,*) " Z reference for FK routine : ", Z_ref_for_FK
    write(IMAIN,*) " Window for FK computing : ", tmax_fk
    write(IMAIN,*) " freqnecy max :", ff0
    write(IMAIN,*) " type of incmoming wave (1=p), (2=sv) :",kpsv
    write(IMAIN,*)
    write(IMAIN,*)
    write(IMAIN,*)

    call flush_IMAIN()
  end subroutine ReadFKModelInput

  !! ----------------------------------
  subroutine FindBoundaryBox(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)
    use specfem_par

    real(kind=CUSTOM_REAL), intent(in out) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box
    real(kind=CUSTOM_REAL)                 :: Xmin_loc, Xmax_loc, Ymin_loc, Ymax_loc, Zmin_loc, Zmax_loc

    Xmin_loc = minval (xstore(:))
    Xmax_loc = maxval (xstore(:))
    Ymin_loc = minval (ystore(:))
    Ymax_loc = maxval (ystore(:))
    Zmin_loc = minval (zstore(:))
    Zmax_loc = maxval (zstore(:))

    call min_all_all_cr(Xmin_loc, Xmin_box)
    call max_all_all_cr(Xmax_loc, Xmax_box)
    call min_all_all_cr(Ymin_loc, Ymin_box)
    call max_all_all_cr(Ymax_loc, Ymax_box)
    call min_all_all_cr(Zmin_loc, Zmin_box)
    call max_all_all_cr(Zmax_loc, Zmax_box)


  end subroutine FindBoundaryBox
  !!==============================================================================================================================
