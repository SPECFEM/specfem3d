!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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


! time iteration using local time stepping
!
! Authors: Max Rietmann, Daniel Peter

  subroutine lts_iterate_time()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie

  use gravity_perturbation, only: gravity_timeseries, GRAVITY_SIMULATION

  use specfem_par_lts


  implicit none

  integer :: ilevel

  ! checks
  if (.not. LTS_MODE) call exit_MPI(myrank,'Error LTS_MODE flag must be set to .true.')

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'using LTS time iteration loop'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! set initial condition to all levels
  do ilevel = 1,num_p_level
    displ_p(:,:,ilevel) = displ(:,:)
  enddo

  ! time loop
  do it = it_begin,it_end

    ! simulation status output and stability check
    if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
      call check_stability()
    endif

    ! simulation status output and stability check
    if (OUTPUT_ENERGY) then
      if (mod(it,NTSTEP_BETWEEN_OUTPUT_ENERGY) == 0 .or. it == 5 .or. it == NSTEP) call compute_energy()
    endif

    ! calculates stiffness term
    if (SIMULATION_TYPE == 3) then
      ! safety stop
      call exit_MPI(myrank,'Error LTS_MODE for kernel simulations (SIMULATION_TYPE == 3) not implemented yet')
      !call lts_global_step(num_p_level)
    else
      ! LTS Newmark update
      ! (only displ update)
      call lts_newmark_update_displ()

      ! finishes a global time step (with multiple steps in dt/p regions)
      call lts_global_step(num_p_level)
    endif

    ! restores last time snapshot saved for backward/reconstruction of wavefields
    ! note: this must be read in after the Newmark time scheme
    if (SIMULATION_TYPE == 3 .and. it == 1) then
      call read_forward_arrays()
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
    if (MOVIE_SIMULATION) call write_movie_output()

    ! first step of noise tomography, i.e., save a surface movie at every time step
    ! modified from the subroutine 'write_movie_surface'
    if (NOISE_TOMOGRAPHY == 1) then
      call noise_save_surface_movie()
    endif

  !
  !---- end of time iteration loop
  !
  enddo   ! end of main time loop

  end subroutine lts_iterate_time
