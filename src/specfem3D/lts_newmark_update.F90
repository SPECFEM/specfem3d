!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

! Local Time Stepping (LTS)
!
! In case you use this local time stepping feature in your study, please reference this work:
!
! Rietmann, M., M. Grote, D. Peter, O. Schenk, 2017
! Newmark local time stepping on high-performance computing architectures,
! Journal of Comp. Physics, 334, p. 308-326.
! https://doi.org/10.1016/j.jcp.2016.11.012
!
! Rietmann, M., B. Ucar, D. Peter, O. Schenk, M. Grote, 2015.
! Load-balanced local time stepping for large-scale wave propagation,
! in: Parallel Distributed Processing Symposium (IPDPS), IEEE International, May 2015.
! https://doi.org/10.1109/IPDPS.2015.10


! Newmark local time step scheme


  subroutine lts_newmark_update_displ()

! updates displacement
!
! note: the LTS newmark update is done at the end in the global step and updates all in one call, i.e. from
!           displ_n -> displ_n+1, veloc_n -> veloc_n+1, given accel_n+1
!
!       this is slightly different to the default Newmark scheme which does two steps, a predictor and corrector step.
!       the predictor is called at the beginning when entering the time loop
!       displ_n -> displ_n+1, veloc_n -> veloc_n1/2, given accel_n.
!       the corrector is called after the compute forces to step veloc_n1/2 -> veloc_n+1, given accel_n+1
!
!       thus, for the default Newmark scheme calling the write_seismograms() after the compute forces
!       will have displ_n+1,veloc_n+1,accel_n+1. note that in the very first iteration, in the predictor step
!       displ_n+1 is updated with veloc_n==accel_n==0.
!       in that case, displ_n+1 is still at zero and the wavefield is for time (it-1) * DT.
!
!       for the LTS scheme, displ_n+1 gets calculated after the global step, based on accel_n+1.
!       thus, even in the first iteration, displ_n+1 corresponds to it * DT, not (it-1) * DT.
!
!       to match the default Newmark scheme, we call here lts_newmark_update_displ() to just copy
!       the updated LTS wavefield, displ_p (coarsest level), into the main wavefield displ (used to output the seismograms)
!       before the global LTS stepping is called.
!       this splits the displ update into the local LTS done after the LTS time stepping, and the copy into the main,
!       similar to the default Newmark scheme.

  !use constants, only: CUSTOM_REAL

  use shared_parameters, only: GPU_MODE
  use specfem_par_elastic, only: displ !,accel
  ! LTS
  use specfem_par_lts, only: num_p_level,displ_p

  implicit none

  if (.not. GPU_MODE) then
    ! on CPU
    ! place final coarse-level step into standard fields
    ! updates displacement
    displ(:,:) = displ_p(:,:,num_p_level)
    ! zeros out acceleration - will be done in global step on specific local levels
    !accel(:,:) = 0.0_CUSTOM_REAL
  else
    ! on GPU
    !#TODO: LTS on GPU newmark update displ
    stop 'LTS on GPU w/ newmark update displ not implemented yet'
  endif

  end subroutine lts_newmark_update_displ

!
!-------------------------------------------------------------------------------------------------
!

  subroutine lts_newmark_update(veloc_p,displ_p,accel,m,ilevel,num_p_level,deltat_lts_local)

! main LTS time update
!
! (requires all p-level accelerations to be computed at this point)

  use constants, only: CUSTOM_REAL,STABILITY_THRESHOLD,NDIM,myrank
  use shared_parameters, only: STACEY_ABSORBING_CONDITIONS,GPU_MODE

  use specfem_par, only: NGLOB_AB

  ! elastic domains
  use specfem_par_elastic, only: veloc

  ! GPU
  !use specfem_par, only: Mesh_pointer

  ! LTS
  use specfem_par_lts, only: iglob_p_refine, p_level_iglob_start, p_level_iglob_end, &
    num_p_level_coarser_to_update, p_level_coarser_to_update, lts_current_m, &
    cmassxyz,rmassxyz,rmassxyz_mod, &
    accel_collected

  ! debug
  use specfem_par, only: it,xstore,ystore,zstore

  implicit none

  integer,intent(in) :: num_p_level

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB,num_p_level),intent(inout) :: displ_p
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB,num_p_level),intent(inout) :: veloc_p
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel

  integer,intent(in) :: ilevel
  real(kind=CUSTOM_REAL),intent(in) :: deltat_lts_local

  ! local parameters
  integer :: m
  integer :: is,ie,is0,ie0
  integer :: iglob_n, iglob

  ! debug
  real(kind=CUSTOM_REAL) :: norm_d,norm_v,norm_a
  integer,dimension(1) :: imax
  logical, parameter :: DEBUG = .false.

  ! we assume contiguous range of p-level nodes
  ! initialization
  is = 0
  ie = 0
  is0 = 0
  ie0 = 0

  ! Local-time stepping update: Newmark scheme

  if (.not. GPU_MODE) then
    ! on CPU

    ! ---- updates to P nodes -----------------------------------

    ! contiguous range of p-level nodes

    ! all up to current level

    ! start index of current p-level
    is = p_level_iglob_start(ilevel)

    ! end index of current p-level
    ie = p_level_iglob_end(ilevel)

    ! gets start index of finest level
    is0 = p_level_iglob_start(1)

    ! end index of next finer p-level
    if (ilevel > 1) then
      ie0 = p_level_iglob_end(ilevel-1)
    endif

    ! checks (done in lts_setup already)
    if (ie < is) stop 'Error lts newmark update: end index of current level invalid'
    if (ie < is0) stop 'Error lts newmark update: start/end index invalid'
    if (ilevel > 1 .and. ie0 < is0) stop 'Error lts newmark update: end index of finer level invalid'

    !
    ! acceleration update
    !
    if (ilevel < num_p_level) then
      ! finer p-levels
      accel(:,is0:ie) = rmassxyz(:,is0:ie) * accel(:,is0:ie)
    else
      ! coarsest p-level
      ! use modified mass matrix here
      ! we note that this is equivalent to: (M+dt/2*C)^{-1} * accel
      accel(:,is0:ie) = rmassxyz_mod(:,is0:ie) * accel(:,is0:ie)
    endif

    !
    ! velocity update
    !
    ! local time step
    if (m == 1 .and. ilevel < num_p_level) then
      ! 1st local time step for finer p-levels
      ! updates finer p-levels up to current p-level
      veloc_p(:,is0:ie,ilevel) = 0.5_CUSTOM_REAL * deltat_lts_local * accel(:,is0:ie)

      if (ilevel > 1) then
        ! only finer p-levels
        veloc_p(:,is0:ie0,ilevel) = veloc_p(:,is0:ie0,ilevel) &
                          + 1.0_CUSTOM_REAL/deltat_lts_local * (displ_p(:,is0:ie0,ilevel-1) - displ_p(:,is0:ie0,ilevel))
        ! only current p-level (excluding coarsest level)
        veloc_p(:,is:ie,ilevel) = veloc_p(:,is:ie,ilevel) &
                          + 1.0_CUSTOM_REAL/deltat_lts_local * displ_p(:,is:ie,ilevel-1)
      endif
    else
      ! after 1st local time step
      if (ilevel < num_p_level) then
        ! finer p-levels
        ! updates finer levels and up to current
        veloc_p(:,is0:ie,ilevel) = veloc_p(:,is0:ie,ilevel) + deltat_lts_local * accel(:,is0:ie)

        ! accumulate contributions from finer levels
        if (ilevel > 1) then
          ! multiple levels
          ! updates all finer nodes
          veloc_p(:,is0:ie0,ilevel) = veloc_p(:,is0:ie0,ilevel) &
                            + 2.0_CUSTOM_REAL/deltat_lts_local * (displ_p(:,is0:ie0,ilevel-1) - displ_p(:,is0:ie0,ilevel))

          ! updates current p-level nodes
          veloc_p(:,is:ie,ilevel) = veloc_p(:,is:ie,ilevel) &
                            + 2.0_CUSTOM_REAL/deltat_lts_local * displ_p(:,is:ie,ilevel-1)
        endif
      else
        ! coarsest p-level

        ! updates on coarse-global level (with Stacey-ABC)
        veloc_p(:,is0:ie,ilevel) = veloc_p(:,is0:ie,ilevel) + deltat_lts_local * accel(:,is0:ie)

        ! these updates are modified by the boundary term (1+Minv*dt/2*C)^{-1}.
        ! This ensure that the method works when fine-elements are on the boundary.
        if (ilevel > 1) then
          ! multiple levels
          ! updates all finer nodes
          if (STACEY_ABSORBING_CONDITIONS) then
            veloc_p(:,is0:ie0,ilevel) = veloc_p(:,is0:ie0,ilevel) &
                              + 1.0_CUSTOM_REAL/(1.0_CUSTOM_REAL+rmassxyz(:,is0:ie0)*cmassxyz(:,is0:ie0)) &
                                * 2.0_CUSTOM_REAL/deltat_lts_local * (displ_p(:,is0:ie0,ilevel-1) - displ_p(:,is0:ie0,ilevel))
          else
            ! note: in case of no absorbing boundaries, cmassxyz is set to zero
            veloc_p(:,is0:ie0,ilevel) = veloc_p(:,is0:ie0,ilevel) &
                              + 2.0_CUSTOM_REAL/deltat_lts_local * (displ_p(:,is0:ie0,ilevel-1) - displ_p(:,is0:ie0,ilevel))
          endif

          ! updates current p-level nodes
          veloc_p(:,is:ie,ilevel) = veloc_p(:,is:ie,ilevel) &
                            + 2.0_CUSTOM_REAL/deltat_lts_local * displ_p(:,is:ie,ilevel-1)
        endif

      endif

    endif

    !
    ! displacement update
    !
    displ_p(:,is0:ie,ilevel) = displ_p(:,is0:ie,ilevel) + deltat_lts_local * veloc_p(:,is0:ie,ilevel)

    ! ---------------- updates to R ----------------------------------------

    ! p-level boundary updates

    ! coarsest level is already finished from above

    if (ilevel < num_p_level) then
      ! finer p-levels
      ! considers contributions from coarser to finer p-levels
      do iglob_n = 1,num_p_level_coarser_to_update(ilevel)
        iglob = p_level_coarser_to_update(iglob_n,ilevel)

        ! checks (done in lts_setup already)
        !if (iglob < ie) call exit_mpi(myrank,"ASSERT: coarser iglob should start in next coarser level")
        !if (iglob > NGLOB_AB) call exit_mpi(myrank,"ASSERT: coarser iglob index is0 out of bounds!")

        !
        ! acceleration update
        !
        if (iglob_p_refine(iglob) > 1) then
          ! node belongs to finer p-levels
          accel(:,iglob) = rmassxyz(:,iglob) * accel(:,iglob)
        else
          ! node belongs to coarsest level
          accel(:,iglob) = rmassxyz_mod(:,iglob) * accel(:,iglob)
        endif

        !
        ! velocity update
        !
        if (m == 1 .and. ilevel < num_p_level) then
          ! 1st local time step for finer p-levels
          veloc_p(:,iglob,ilevel) = 0.5_CUSTOM_REAL * deltat_lts_local * accel(:,iglob)

          if (ilevel > 1) then
            veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) &
                    + 1.0_CUSTOM_REAL/deltat_lts_local * displ_p(:,iglob,ilevel-1)
          endif
        else
          ! all subsequent local time steps after 1st step
          veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + deltat_lts_local * accel(:,iglob)

          ! accumulate contributions from finer levels
          if (ilevel > 1) then
            ! multiple levels
            veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) &
                    + 2.0_CUSTOM_REAL/deltat_lts_local * displ_p(:,iglob,ilevel-1)
          endif
        endif

        !
        ! displacement update
        !
        displ_p(:,iglob,ilevel) = displ_p(:,iglob,ilevel) + deltat_lts_local * veloc_p(:,iglob,ilevel)
      enddo
    endif

    ! collects acceleration a_n+1 = B u_n+1 from diffferent levels
    if (allocated(accel_collected)) then
      ! collects acceleration B u_n+1 from this current level ilevel
      ! this is computed in the first local iteration (m==1) where the initial condition sets u_0 = u_n+1
      if (ilevel < num_p_level) then
        ! only store if accel was computed the very first time this level was called
        !(i.e., lts_current_m(ilevel+1) is still at 1)
        if (m == 1 .and. lts_current_m(ilevel+1) == 1) then
          if (STACEY_ABSORBING_CONDITIONS) then
            ! uses same mass matrix as on coarsest level
            accel_collected(:,is:ie) = rmassxyz_mod(:,is:ie)/rmassxyz(:,is:ie) * accel(:,is:ie)
          else
            accel_collected(:,is:ie) = accel(:,is:ie)
          endif
        endif
      else
        accel_collected(:,is:ie) = accel(:,is:ie)
      endif
    endif

    ! updates global main wavefields
    !
    ! at this point, when calling this update recursively for the last level ilevel==num_p_level, which is the coarsest level,
    ! all p-levels have been updated and we can copy the coarsest-level LTS array into the main wavefield.
    !
    ! note: we will update only the velocity here to match the default Newmark scheme, which uses a predictor and corrector step.
    !       there, displacement will get updated in the predictor step called before the compute forces routines, while velocity gets
    !       updated afterwards.
    !       to match that order, we call lts_newmark_update_displ() explicitly before the main time stepping to update the displ array.
    !
    !this would be both at once:
    ! place final coarse-level step into standard fields
    !if (ilevel == num_p_level) then
    !  displ(:,:) = displ_p(:,:,num_p_level)
    !  veloc(:,:) = veloc_p(:,:,num_p_level)
    !endif
    !
    !now, only veloc to match the default Newmark behavior:
    ! place final coarse-level step into standard fields
    if (ilevel == num_p_level) then
      veloc(:,:) = veloc_p(:,:,num_p_level)

      ! sets accel to collected accel wavefield for (possible) seismogram outputs or shakemaps
      ! note: for the next time loop iteration, accel will be zeroed out again,
      !       thus this change will have no effect on the time marching.
      if (allocated(accel_collected)) then
        accel(:,:) = accel_collected(:,:)
      endif
    endif

  else
    ! on GPU
    !#TODO: LTS on GPU w/ Newmark update
    stop 'LTS on GPU w/ Newmark update not implemented yet'
    ! newmark update
    !call finalize_step_lts(Mesh_pointer,deltat_lts_local,ilevel,m, &
    !                       p_level_iglob_start,p_level_iglob_end,num_p_level_coarser_to_update)
  endif ! GPU_MODE

  ! debug
  if (DEBUG) then
    ! stability check
    ! acceleration
    norm_a = maxval(sqrt(accel(1,:)**2 + accel(2,:)**2 + accel(3,:)**2))
    if (norm_a > STABILITY_THRESHOLD .or. norm_a < 0.0_CUSTOM_REAL .or. norm_a /= norm_a) then
      ! user output
      imax = maxloc(accel(1,:)**2 + accel(2,:)**2 + accel(3,:)**2)
      print *,'Error: simulation blew up in rank',myrank,'at time step: ',it,'local m:',lts_current_m(ilevel),'level:',ilevel
      print *,'  norm of acceleration maximum value:',norm_a
      print *,'  maximum location: iglob = ',imax(1)
      print *,'  x/y/z           :',xstore(imax(1)),ystore(imax(1)),zstore(imax(1))
      print *,'please consider decreasing time step size DT in Par_file...'
      call exit_MPI(myrank,'simulation became unstable and blew up')
    endif

    ! velocity
    norm_v = maxval(sqrt(veloc_p(1,:,ilevel)**2 + veloc_p(2,:,ilevel)**2 + veloc_p(3,:,ilevel)**2))
    if (norm_v > STABILITY_THRESHOLD .or. norm_v < 0.0_CUSTOM_REAL .or. norm_v /= norm_v) then
      ! user output
      imax = maxloc(veloc_p(1,:,ilevel)**2 + veloc_p(2,:,ilevel)**2 + veloc_p(3,:,ilevel)**2)
      print *,'Error: simulation blew up in rank',myrank,'at time step: ',it,'local m:',lts_current_m(ilevel),'level:',ilevel
      print *,'  norm of velocity maximum value:',norm_v
      print *,'  maximum location: iglob = ',imax(1)
      print *,'  x/y/z           :',xstore(imax(1)),ystore(imax(1)),zstore(imax(1))
      print *,'please consider decreasing time step size DT in Par_file...'
      call exit_MPI(myrank,'simulation became unstable and blew up')
    endif

    ! displacement
    norm_d = maxval(sqrt(displ_p(1,:,ilevel)**2 + displ_p(2,:,ilevel)**2 + displ_p(3,:,ilevel)**2))
    if (norm_d > STABILITY_THRESHOLD .or. norm_d < 0.0_CUSTOM_REAL .or. norm_d /= norm_d) then
      ! user output
      imax = maxloc(displ_p(1,:,ilevel)**2 + displ_p(2,:,ilevel)**2 + displ_p(3,:,ilevel)**2)
      print *,'Error: simulation blew up in rank',myrank,'at time step: ',it,'local m:',lts_current_m(ilevel),'level:',ilevel
      print *,'  norm of displacement maximum value:',norm_d
      print *,'  maximum location: iglob = ',imax(1)
      print *,'  x/y/z           :',xstore(imax(1)),ystore(imax(1)),zstore(imax(1))
      print *,'please consider decreasing time step size DT in Par_file...'
      call exit_MPI(myrank,'simulation became unstable and blew up')
    endif
  endif ! DEBUG

  end subroutine lts_newmark_update

!
!-------------------------------------------------------------------------------
!

! unused...

!  subroutine lts_newmark_update_v2(veloc_p,displ_p,accel,rmassx,m,ilevel,num_p_level,deltat_lts_local)
!
!  use specfem_par, only: CUSTOM_REAL, NGLOB_AB, NDIM,GPU_MODE,Mesh_pointer,myrank
!  use specfem_par_lts, only: p_level_iglob_start, p_level_iglob_end, &
!       num_p_level_coarser_to_update, p_level_coarser_to_update
!
!  implicit none
!
!  integer :: num_p_level
!
!  real(kind=CUSTOM_REAL) displ_p(NDIM,NGLOB_AB,num_p_level)
!  real(kind=CUSTOM_REAL) veloc_p(NDIM,NGLOB_AB,num_p_level)
!  real(kind=CUSTOM_REAL) accel(NDIM,NGLOB_AB)
!  real(kind=CUSTOM_REAL) rmassx(NGLOB_AB)
!
!  integer :: m,ilevel,is,ie,ie0, is2, iglob_n, iglob
!  real(kind=CUSTOM_REAL) :: deltat_lts_local
!
!  logical,parameter :: USE_FAST = .true.
!
!  if (GPU_MODE) then
!
!    if (p_level_iglob_start(1) /= 1) call exit_mpi(myrank,"ASSERT FAIL: iglob should start at 1")
!    call finalize_step_lts(Mesh_pointer,deltat_lts_local,ilevel,m, &
!         p_level_iglob_start,p_level_iglob_end,num_p_level_coarser_to_update)
!
!  else ! GPU MODE
!
!    ! ---- updates to P nodes -----------------------------------
!
!    ! contiguous range of p-level nodes
!
!    is = p_level_iglob_start(1)
!    ie = p_level_iglob_end(ilevel)
!
!    if (ilevel > 1) then
!      ie0 = p_level_iglob_end(ilevel-1)
!      is2 = p_level_iglob_start(ilevel)
!    endif
!
!    if (ilevel == num_p_level .and. ie /= NGLOB_AB) call exit_mpi(myrank,"ASSERT: coarsest level must have all iglobs!")
!
!    accel(1,is:ie) = rmassx(is:ie) * accel(1,is:ie)
!    accel(2,is:ie) = rmassx(is:ie) * accel(2,is:ie)
!    accel(3,is:ie) = rmassx(is:ie) * accel(3,is:ie)
!
!    ! call copy_field_from_gpu(Mesh_pointer,accel_2)
!    ! print *, "maxval(rmass*accel)=",maxval(abs(accel)), &
!    ! "maxval(rmass*accel-gpu)",maxval(abs(displ_p(:,:,ilevel))), &
!
!    if (m==1 .and. ilevel < num_p_level) then
!      veloc_p(:,is:ie,ilevel) = deltat_lts_local/2.0 * accel(:,is:ie)
!
!      if (ilevel > 1) then
!        veloc_p(:,is:ie0,ilevel) = veloc_p(:,is:ie0,ilevel) + 1.0/deltat_lts_local * &
!             (displ_p(:,is:ie0,ilevel-1) - displ_p(:,is:ie0,ilevel))
!        veloc_p(:,is2:ie,ilevel) = veloc_p(:,is2:ie,ilevel) + 1.0/deltat_lts_local * &
!             (displ_p(:,is2:ie,ilevel-1))
!      endif
!
!    else ! after 1st step
!      veloc_p(:,is:ie,ilevel) = veloc_p(:,is:ie,ilevel) + deltat_lts_local*accel(:,is:ie)
!
!      ! accumulate contributions from finer levels
!      if (ilevel > 1) then
!        veloc_p(:,is:ie0,ilevel) = veloc_p(:,is:ie0,ilevel) + 2.0/deltat_lts_local * &
!             (displ_p(:,is:ie0,ilevel-1) - displ_p(:,is:ie0,ilevel))
!        veloc_p(:,is2:ie,ilevel) = veloc_p(:,is2:ie,ilevel) + 2.0/deltat_lts_local * &
!             (displ_p(:,is2:ie,ilevel-1))
!      endif
!
!    endif
!
!    displ_p(:,is:ie,ilevel) = displ_p(:,is:ie,ilevel) + deltat_lts_local * veloc_p(:,is:ie,ilevel)
!
!    ! ---------------- updates to R ----------------------------------------
!
!    ! coarsest level is already finished from above
!    if (USE_FAST) then
!      if (ilevel < num_p_level) then
!
!        do iglob_n=1,num_p_level_coarser_to_update(ilevel)
!          iglob = p_level_coarser_to_update(iglob_n,ilevel)
!
!          accel(1,iglob) = rmassx(iglob) * accel(1,iglob)
!          accel(2,iglob) = rmassx(iglob) * accel(2,iglob)
!          accel(3,iglob) = rmassx(iglob) * accel(3,iglob)
!
!          if (m==1 .and. ilevel < num_p_level) then
!
!            veloc_p(:,iglob,ilevel) = deltat_lts_local/2.0 * accel(:,iglob)
!
!            if (ilevel > 1) &
!                 veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + 1.0/deltat_lts_local * &
!                 (displ_p(:,iglob,ilevel-1))
!          else
!            veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + deltat_lts_local*accel(:,iglob)
!
!            ! accumulate contributions from finer levels
!            if (ilevel > 1) &
!                 veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + 2.0/deltat_lts_local * &
!                 (displ_p(:,iglob,ilevel-1))
!          endif
!
!          displ_p(:,iglob,ilevel) = displ_p(:,iglob,ilevel) + deltat_lts_local * veloc_p(:,iglob,ilevel)
!
!        enddo
!      endif
!    else
!      ! coarsest level is already finished from above
!      if (ilevel < num_p_level) then
!
!        do iglob_n=1,num_p_level_coarser_to_update(ilevel)
!          iglob = p_level_coarser_to_update(iglob_n,ilevel)
!
!          accel(1,iglob) = rmassx(iglob) * accel(1,iglob)
!          accel(2,iglob) = rmassx(iglob) * accel(2,iglob)
!          accel(3,iglob) = rmassx(iglob) * accel(3,iglob)
!
!          if (m==1 .and. ilevel < num_p_level) then
!
!            veloc_p(:,iglob,ilevel) = deltat_lts_local/2.0 * accel(:,iglob)
!
!            if (ilevel > 1) &
!                 veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + 1.0/deltat_lts_local * &
!                 (displ_p(:,iglob,ilevel-1) - displ_p(:,iglob,ilevel))
!          else
!            veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + deltat_lts_local*accel(:,iglob)
!
!            ! accumulate contributions from finer levels
!            if (ilevel > 1) &
!                 veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + 2.0/deltat_lts_local * &
!                 (displ_p(:,iglob,ilevel-1) - displ_p(:,iglob,ilevel))
!          endif
!
!          displ_p(:,iglob,ilevel) = displ_p(:,iglob,ilevel) + deltat_lts_local * veloc_p(:,iglob,ilevel)
!
!        enddo
!      endif
!    endif ! fast mode
!
!  endif ! GPU_MODE
!
!
!  end subroutine lts_newmark_update_v2
!

!
!-------------------------------------------------------------------------------
!


! unused...

!  subroutine lts_newmark_update_ref(veloc_p,displ_p,accel,m,ilevel,num_p_level,deltat_lts_local)
!
!  use specfem_par, only: CUSTOM_REAL, NGLOB_AB, NDIM,GPU_MODE,Mesh_pointer,myrank
!  use specfem_par_lts, only: p_level_iglob_start, p_level_iglob_end, &
!       num_p_level_coarser_to_update, p_level_coarser_to_update, &
!       rmassxyz_mod
!
!  implicit none
!
!  integer :: num_p_level
!
!  real(kind=CUSTOM_REAL) displ_p(NDIM,NGLOB_AB,num_p_level)
!  real(kind=CUSTOM_REAL) veloc_p(NDIM,NGLOB_AB,num_p_level)
!  real(kind=CUSTOM_REAL) accel(NDIM,NGLOB_AB)
!
!  integer :: m,ilevel,is,ie,ie0, is2, iglob_n, iglob
!  real(kind=CUSTOM_REAL) :: deltat_lts_local
!
!  logical,parameter :: USE_FAST = .true.
!
!  if (GPU_MODE) then
!
!    if (p_level_iglob_start(1) /= 1) call exit_mpi(myrank,"ASSERT FAIL: iglob should start at 1")
!    call finalize_step_lts(Mesh_pointer,deltat_lts_local,ilevel,m, &
!         p_level_iglob_start,p_level_iglob_end,num_p_level_coarser_to_update)
!
!  else ! GPU MODE
!
!    ! ---- updates to P nodes -----------------------------------
!
!    ! contiguous range of p-level nodes
!
!    is = p_level_iglob_start(1)
!    ie = p_level_iglob_end(ilevel)
!
!    if (ilevel > 1) then
!      ie0 = p_level_iglob_end(ilevel-1)
!      is2 = p_level_iglob_start(ilevel)
!    endif
!
!    if (ilevel == num_p_level .and. ie /= NGLOB_AB) call exit_mpi(myrank,"ASSERT: coarsest level must have all iglobs!")
!
!    accel(1,is:ie) = rmassxyz_mod(1,is:ie) * accel(1,is:ie)
!    accel(2,is:ie) = rmassxyz_mod(2,is:ie) * accel(2,is:ie)
!    accel(3,is:ie) = rmassxyz_mod(3,is:ie) * accel(3,is:ie)
!
!    ! call copy_field_from_gpu(Mesh_pointer,accel_2)
!    ! print *, "maxval(rmass*accel)=",maxval(abs(accel)), &
!    ! "maxval(rmass*accel-gpu)",maxval(abs(displ_p(:,:,ilevel))), &
!
!    if (m==1 .and. ilevel < num_p_level) then
!      veloc_p(:,is:ie,ilevel) = deltat_lts_local/2.0 * accel(:,is:ie)
!
!      if (ilevel > 1) then
!        veloc_p(:,is:ie0,ilevel) = veloc_p(:,is:ie0,ilevel) + 1.0/deltat_lts_local * &
!             (displ_p(:,is:ie0,ilevel-1) - displ_p(:,is:ie0,ilevel))
!        veloc_p(:,is2:ie,ilevel) = veloc_p(:,is2:ie,ilevel) + 1.0/deltat_lts_local * &
!             (displ_p(:,is2:ie,ilevel-1))
!      endif
!
!    else ! after 1st step
!      veloc_p(:,is:ie,ilevel) = veloc_p(:,is:ie,ilevel) + deltat_lts_local*accel(:,is:ie)
!
!      ! accumulate contributions from finer levels
!      if (ilevel > 1) then
!        veloc_p(:,is:ie0,ilevel) = veloc_p(:,is:ie0,ilevel) + 2.0/deltat_lts_local * &
!             (displ_p(:,is:ie0,ilevel-1) - displ_p(:,is:ie0,ilevel))
!        veloc_p(:,is2:ie,ilevel) = veloc_p(:,is2:ie,ilevel) + 2.0/deltat_lts_local * &
!             (displ_p(:,is2:ie,ilevel-1))
!      endif
!
!    endif
!
!    displ_p(:,is:ie,ilevel) = displ_p(:,is:ie,ilevel) + deltat_lts_local * veloc_p(:,is:ie,ilevel)
!
!    ! ---------------- updates to R ----------------------------------------
!
!    ! coarsest level is already finished from above
!    if (USE_FAST) then
!      if (ilevel < num_p_level) then
!
!        do iglob_n=1,num_p_level_coarser_to_update(ilevel)
!          iglob = p_level_coarser_to_update(iglob_n,ilevel)
!
!          accel(1,iglob) = rmassxyz_mod(1,iglob) * accel(1,iglob)
!          accel(2,iglob) = rmassxyz_mod(2,iglob) * accel(2,iglob)
!          accel(3,iglob) = rmassxyz_mod(3,iglob) * accel(3,iglob)
!
!          if (m==1 .and. ilevel < num_p_level) then
!
!            veloc_p(:,iglob,ilevel) = deltat_lts_local/2.0 * accel(:,iglob)
!
!            if (ilevel > 1) &
!                 veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + 1.0/deltat_lts_local * &
!                 (displ_p(:,iglob,ilevel-1))
!          else
!            veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + deltat_lts_local*accel(:,iglob)
!
!            ! accumulate contributions from finer levels
!            if (ilevel > 1) &
!                 veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + 2.0/deltat_lts_local * &
!                 (displ_p(:,iglob,ilevel-1))
!          endif
!
!          displ_p(:,iglob,ilevel) = displ_p(:,iglob,ilevel) + deltat_lts_local * veloc_p(:,iglob,ilevel)
!
!        enddo
!      endif
!    else
!      ! coarsest level is already finished from above
!      if (ilevel < num_p_level) then
!
!        do iglob_n=1,num_p_level_coarser_to_update(ilevel)
!          iglob = p_level_coarser_to_update(iglob_n,ilevel)
!
!          accel(:,iglob) = rmassxyz_mod(:,iglob) * accel(:,iglob)
!
!          if (m==1 .and. ilevel < num_p_level) then
!
!            veloc_p(:,iglob,ilevel) = deltat_lts_local/2.0 * accel(:,iglob)
!
!            if (ilevel > 1) &
!                 veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + 1.0/deltat_lts_local * &
!                 (displ_p(:,iglob,ilevel-1) - displ_p(:,iglob,ilevel))
!          else
!            veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + deltat_lts_local*accel(:,iglob)
!
!            ! accumulate contributions from finer levels
!            if (ilevel > 1) &
!                 veloc_p(:,iglob,ilevel) = veloc_p(:,iglob,ilevel) + 2.0/deltat_lts_local * &
!                 (displ_p(:,iglob,ilevel-1) - displ_p(:,iglob,ilevel))
!          endif
!
!          displ_p(:,iglob,ilevel) = displ_p(:,iglob,ilevel) + deltat_lts_local * veloc_p(:,iglob,ilevel)
!
!        enddo
!      endif
!    endif ! fast mode
!
!  endif ! GPU_MODE
!
!
! end subroutine lts_newmark_update_ref
!

!
!----------------------------------------------------------------------------------
!

! unused...

!  subroutine lts_newmark_update_slow(veloc_p,displ_p,accel,rmassx,m,ilevel,num_p_level,deltat_lts_local)
!
!  use specfem_par, only: CUSTOM_REAL, NGLOB_AB, NDIM
!
!  implicit none
!
!  integer :: num_p_level
!
!  real(kind=CUSTOM_REAL) displ_p(NDIM,NGLOB_AB,num_p_level)
!  real(kind=CUSTOM_REAL) veloc_p(NDIM,NGLOB_AB,num_p_level)
!  real(kind=CUSTOM_REAL) accel(NDIM,NGLOB_AB)
!  real(kind=CUSTOM_REAL) rmassx(NGLOB_AB)
!
!  integer :: m,ilevel
!  real(kind=CUSTOM_REAL) :: deltat_lts_local
!
!  accel(1,:) = rmassx(:) * accel(1,:)
!  accel(2,:) = rmassx(:) * accel(2,:)
!  accel(3,:) = rmassx(:) * accel(3,:)
!
!  ! call copy_field_from_gpu(Mesh_pointer,accel_2)
!  ! print *, "maxval(rmass*accel)=",maxval(abs(accel)), &
!  ! "maxval(rmass*accel-gpu)",maxval(abs(displ_p(:,:,ilevel))), &
!
!  if (m==1 .and. ilevel < num_p_level) then
!    veloc_p(:,:,ilevel) = deltat_lts_local/2.0 * accel(:,:)
!
!    if (ilevel > 1) &
!         veloc_p(:,:,ilevel) = veloc_p(:,:,ilevel) + 1.0/deltat_lts_local * &
!         (displ_p(:,:,ilevel-1) - displ_p(:,:,ilevel))
!  else ! after 1st step
!    veloc_p(:,:,ilevel) = veloc_p(:,:,ilevel) + deltat_lts_local*accel(:,:)
!
!    ! accumulate contributions from finer levels
!    if (ilevel > 1) &
!         veloc_p(:,:,ilevel) = veloc_p(:,:,ilevel) + 2.0/deltat_lts_local * &
!         (displ_p(:,:,ilevel-1) - displ_p(:,:,ilevel))
!  endif
!
!  displ_p(:,:,ilevel) = displ_p(:,:,ilevel) + deltat_lts_local * veloc_p(:,:,ilevel)
!
!  end subroutine lts_newmark_update_slow
!

