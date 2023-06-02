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

! (recursive) time step using local time stepping
!
! Authors: Max Rietmann, Daniel Peter

  recursive subroutine lts_global_step(ilevel)

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use specfem_par_lts

  use fault_solver_dynamic, only: bc_dynflt_set3d_all,SIMULATION_TYPE_DYN
  use fault_solver_kinematic, only: bc_kinflt_set_all,SIMULATION_TYPE_KIN

  ! debugging
  !use fault_solver_dynamic, only: faults

  implicit none

  integer, intent(in) :: ilevel

  ! local parameters
  integer :: iphase,m,p
  integer :: iilevel,p_relative
  real(kind=CUSTOM_REAL) :: deltat_lts_local
  double precision :: dt_local

  !#TODO: LTS w/ fault daniel debug
  !real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB) :: tmp_displ,tmp_veloc
  !real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: tmp_field

  logical :: backward_simulation

  ! debugging
  logical,parameter :: DEBUG_TIME_STEPPING = .false.
  logical,parameter :: DEBUG_VTK_OUTPUT    = .false.
  !integer :: iproc,ier
  !integer :: i,iglob
  !integer :: it_total

  ! safety check
  ! only for forward simulations at the moment
  if (SIMULATION_TYPE /= 1) stop 'LTS only implemented for forward simulation so far.'

  ! forward simulations for now
  backward_simulation = .false.

  ! debug
  !if (myrank == 0) print *,"debug: lts global_step: time step = ",it,"level = ",ilevel,"p = ",p_level(ilevel)

  ! initializes
  ! note: we call this routine lts_global_step recursively.
  !       therefore, we want to re-set the lts_current_m(:) array only for the first time it is called,
  !       which is with ilevel==num_p_level
  if (ilevel == num_p_level) then
    ! keeps track of current local time step
    lts_current_m(:) = 0
    ! local time step increments
    lts_it_local = 0
  endif

  ! set initial condition for non-coarsest levels
  if (ilevel < num_p_level) call lts_set_finer_initial_condition(displ_p,ilevel)

  ! every level but coarsest does 2 or more loops
  do m = 1,p_level_m_loops(ilevel)
    ! keep track of each level for accurate time
    lts_current_m(ilevel) = m

    ! recursive routine call
    ! computes next finer level
    if (ilevel > 1) call lts_global_step(ilevel-1)

    ! once the finer levels finished, we continue here with the current ilevel.
    !
    ! updates local time step increment
    lts_it_local = lts_it_local + 1

    ! sets current p-level elements
    current_lts_elem => p_elem(:,ilevel)

    ! sets current boundary-level elements
    current_lts_boundary_elem => boundary_elem(:,ilevel)

    ! determines local time
    ! (used in attenuation scheme for compute forces)
    !
    ! p (dt/p) refinement number in this specified p-level
    p = p_level(ilevel)
    ! local time step
    deltat_lts_local = real( DT / dble(p), kind=CUSTOM_REAL)

    if (ATTENUATION) then
      ! determines alphaval,betaval,gammaval for runge-kutta scheme based on local time step
      call get_attenuation_memory_values(tau_sigma,deltat_lts_local,alphaval,betaval,gammaval)
    endif

    ! sets accel to zero
    call lts_set_accel_zero(accel,ilevel)

    !#TODO: LTS w/ faults daniel debug:
    !if (FAULT_SIMULATION) then
    !  tmp_displ(:,:) = 0.0_CUSTOM_REAL
    !  tmp_veloc(:,:) = 0.0_CUSTOM_REAL
    !endif

    ! compute a = B*P*d + B*R*d
    !
    ! loops over outer (1) / inner (2) elements
    do iphase = 1,2
      ! compute K*P*d on strictly P elements
      if (.not. GPU_MODE) then
        ! on CPU
        ! sets compute force type
        lts_type_compute_pelem = .true.   ! p-elements only (no boundary elements)

        call compute_forces_viscoelastic(iphase,deltat_lts_local, &
                                         displ_p(:,:,ilevel),veloc_p(:,:,ilevel),accel, &
                                         alphaval,betaval,gammaval, &
                                         R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                         R_trace_lddrk, &
                                         R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                         epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                         epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                                         backward_simulation)
      else
        ! on GPU
        !#TODO: LTS on GPU compute_forces
        stop 'LTS on GPU w/ compute_forces not implemented yet'
        ! if (iphase == 1) then
        !   num_elements = nspec_outer_elastic
        ! else
        !   num_elements = nspec_inner_elastic
        ! endif
        ! iispec = 1
        ! do ispec=1,num_elements
        !   ispec_p = phase_ispec_inner_elastic(ispec,iphase)
        !   if (.not. p_elem(ispec_p,ilevel)) cycle
        !   if (boundary_elem(ispec_p,ilevel) .eqv. .true.) cycle
        !   element_list(iispec) = ispec_p
        !   iispec = iispec + 1
        ! enddo
        !call compute_forces_viscoelastic_cuda_lts(Mesh_pointer,ilevel,iphase)
      endif

      ! compute "assembly" contributions from coarser and finer levels
      ! B*R*d + B*F*d
      if (.not. GPU_MODE) then
        ! on CPU
        call lts_boundary_contribution(displ_p, veloc_p, accel, ilevel, &
                                       p, deltat_lts_local, &
                                       iphase, backward_simulation)
      else
        ! on GPU
        !#TODO: LTS on GPU boundary contribution
        stop 'LTS on GPU w/ boundary contribution not implemented yet'
        !call lts_boundary_contribution_cuda(Mesh_pointer,num_p_level_boundary_nodes(iphase,ilevel), &
        !                                    num_p_level_boundary_ispec(iphase,ilevel), &
        !                                    iphase,ilevel)
      endif

      ! computes additional contributions
      if (iphase == 1) then
        ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
        if (STACEY_ABSORBING_CONDITIONS) then
          ! (optional) it step for storage of boundary contributions
          ! example: level=1,p=2,m=1 -> it_local = 1 - it=1, NSTEP_LOCAL=3 -> it_total = 1
          !          level=1,p=2,m=2 -> it_local = 2                          it_total = 2
          !          level=2,p=1,m=1 -> it_local = 3                          it_total = 3
          !          level=1,p=2,m=1 -> it_local = 1 - it=2, NSTEP_LOCAL=3 -> it_total = 4
          !          level=1,p=2,m=2 -> it_local = 2                          it_total = 5
          !          ..                                                       ..
          !it_total = (it-1)*NSTEP_LOCAL + lts_it_local
          !debug
          !print *,'debug: lts it/it_local/NSTEP_LOCAL/it_total',it,lts_it_local,NSTEP_LOCAL,it_total

          ! absorbing boundary conditions are only done on the coarsest level;
          ! it is run on the full boundary (including other p-levels)
          !
          ! note: in this case, the absorbing boundary contribution is only called once per it,
          !       thus we don't need it_total, but can use it as usual instead for storage.
          !
          !       also, veloc_p(:,:,ilevel==num_p_level) is still at time n and hasn't been updated to n+1.
          !       thus, veloc_p(:,:,num_p_level) is still equal to veloc(:,:).
          !
          !       finally, interpolation of veloc to current local time step (veloc = veloc + 1/2 dt accel) doesn't help here.
          !       this would make the boundary contribution worse with accel being the accel_n wavefield. one would need to
          !       figure out how to properly update accel and veloc_p on refinement-levels (ilevel < num_p_level).
          !
          !       for now, the boundary should be fine if the absorbing boundary surface is all on coarsest mesh elements.
          !       it is therefore a meshing problem to create such a mesh with local refinements only inside the mesh.
          !
          !#TODO: LTS improve the absorbing boundary on local refinement levels.

          if (ilevel == num_p_level) then
            if (.not. GPU_MODE) then
              ! on CPU
              call compute_stacey_viscoelastic_forward(NSPEC_AB,NGLOB_AB,accel, &
                                                       ibool,iphase, &
                                                       abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                       abs_boundary_ijk,abs_boundary_ispec, &
                                                       num_abs_boundary_faces, &
                                                       veloc,rho_vp,rho_vs, &
                                                       ispec_is_elastic, &
                                                       it,b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)
            else
              ! on GPU
              !#TODO: LTS on GPU stacey contribution
              stop 'LTS on GPU w/ Stacey contribution not implemented yet'
              !call compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
              !                                     NSTEP,it, &
              !                                     b_num_abs_boundary_faces,b_reclen_field,b_absorb_field, &
              !                                     Mesh_pointer,1) ! 1 == forward
            endif
          endif
        endif

        ! source contributions
        !
        ! calculate time for this level
        current_lts_time = dble(it-1)*DT - t0
        ! adds local step increments
        do iilevel = ilevel,(num_p_level-1)
          current_lts_time = current_lts_time + dble(lts_current_m(iilevel)-1) * (dble(DT) / dble(p_level(iilevel)))
        enddo

        ! corrects timing for sources in refined LTS levels
        ! note: for better accuracy, shifting the start time helps to improve accuracy;
        !       not sure why this must/should be done for sources in refined p-levels
        !
        ! examples for shifting seismograms
        ! 2-levels:
        ! level 1 - local_p 1/2: shift by +0.02 w/ DT 0.04 ; +0.015  w/ DT 0.03
        !                        -> 1 * 1/2 * 0.02
        ! level 1 - local_p 1/4: shift by +0.05 w/ DT 0.04 ; +0.0375 w/ DT 0.03; +0.025 w/ DT 0.02
        !                        -> 1 * 1/2 * 0.02 + 3 * 1/4 * 0.02
        ! level 1 - local_p 1/8: shift by +0.10 w/ DT 0.04 ; +0.075  w/ DT 0.03; +0.05 w/ DT 0.02
        !                        -> 1 * 1/2 * 0.02 + 7 * 1/4 * 0.02 (+ 1 * 1/4 * 0.02)
        ! 3-levels:
        ! level 1 - local_p 1/8 & 1/4; shift by +0.12 w/ DT 0.04 ; +0.09   w/ DT 0.03 ; +0.06 w/ DT 0.02
        !                              -> 1 * 1/2 * 0.02 + (1 * 1/2 * 0.04 + 3 * 1/4 * 0.04)
        ! level 2 - local_p 1/4;       shift by +0.10 w/ DT 0.04 ; +0.075  w/ DT 0.03 ; +0.05 w/ DT 0.02
        !                              -> 1 * 1/2 * 0.04 + 3 * 1/4 * 0.04
        ! this formula to shift is empirical and has no theoretical basis yet - todo in future to figure out a correct way!
        do iilevel = ilevel,(num_p_level-1)
          ! local time step
          dt_local = dble(DT) / dble(p_level(iilevel))
          current_lts_time = current_lts_time - 0.5d0 * dt_local                         ! corrects m == 1 update?
          ! relative refinements compared to next level:
          ! e.g., 8/4 -> 2, 4/1 -> 4
          p_relative = p_level(iilevel)/p_level(iilevel+1)
          select case(p_relative)
          case (2)
            continue
          case (4)
            current_lts_time = current_lts_time - 3 * 0.25d0 * dt_local                  ! corrects m > 1 updates?
          case (8)
            current_lts_time = current_lts_time - 8 * 0.25d0 * dt_local                  ! corrects m > 1 updates?
          case default
            current_lts_time = current_lts_time - (p_relative - 1) * 0.25d0 * dt_local   ! corrects m > 1 updates?
          end select
        enddo
        !debug
        !print *,'debug: lts global_step: final time it/ilevel/it_local/current_lts_time',it,ilevel,lts_it_local,current_lts_time

        ! adds source term (single-force/moment-tensor solution)
        if (.not. GPU_MODE) then
          ! on CPU
          call compute_add_sources_viscoelastic(accel)
        else
          ! on GPU
          call compute_add_sources_viscoelastic_GPU()
        endif
      endif ! iphase

      ! assemble all the contributions between slices using MPI
      if (iphase == 1) then
        ! sends accel values to corresponding MPI interface neighbors
        call assemble_MPI_vector_async_send_lts(NPROC,NGLOB_AB,accel,ilevel, &
                                                buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                                num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                                my_neighbors_ext_mesh, &
                                                request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
      else
        ! waits for send/receive requests to be completed and assembles values
        ! receives MPI buffers
        call assemble_MPI_vector_async_recv_lts(NPROC,NGLOB_AB,accel,ilevel, &
                                                buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                                max_nibool_interfaces_ext_mesh, &
                                                nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                                request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                                my_neighbors_ext_mesh)
      endif
    enddo ! iphase

    ! debugging
    ! loops over each process
    !do iproc = 0,NPROC-1
    !  if (myrank == iproc) then
    !    ! user output
    !    print *,'debug: lts global_step: rank:',myrank,'ilevel:',ilevel,'m:',m,'before newmark ', &
    !            'max displ:',maxval(displ_p(:,:,ilevel)),'veloc:',maxval(veloc_p(:,:,ilevel)),'accel:',maxval(accel(:,:))
    !  endif
    !  call synchronize_all()
    !enddo

    !#TODO: LTS w/ faults daniel debug: add fault routines for LTS - TODO
    if (FAULT_SIMULATION) then
      ! note: fault rupture routines would need current displ and veloc to update accel accordingly.
      !       however, in LTS the fine-level stepping uses displ_p and veloc_p which are not the same as interpolated veloc.
      !       this can be seen for example by the initial conditions for displ_p and veloc_p, in particular veloc_p is different.
      !
      !       this is a similar problem for absorbing boundary conditions in case those are on refined p-elements.
      !
      ! for now, only update all when on the coarsest level
      if (ilevel == num_p_level) then
        !#TODO: LTS update veloc to current local time step

        ! debug
        !tmp_displ(:,:) = displ(:,:)
        !tmp_veloc(:,:) = veloc(:,:)
        ! displ
        !write(filename,'(a,i6.6,a,i6.6,a,i1.1,a,i2.2)') '../OUTPUT_FILES/proc',myrank,'_displ_it',it,'_level',ilevel,'_m',m
        !call write_VTK_wavefield(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,tmp_displ,filename)
        !print *,'written file: ',trim(filename)//'.vtk'
        ! veloc
        !write(filename,'(a,i6.6,a,i6.6,a,i1.1,a,i2.2)') '../OUTPUT_FILES/proc',myrank,'_veloc_it',it,'_level',ilevel,'_m',m
        !call write_VTK_wavefield(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,tmp_veloc,filename)
        !print *,'written file: ',trim(filename)//'.vtk'
        ! saves points on acoustic-elastic coupling interface
        !write(filename,'(a,i6.6,a,i6.6,a,i1.1,a,i2.2)') '../OUTPUT_FILES/proc',myrank,'_fault1_it',it,'_level',ilevel,'_m',m
        !allocate( tmp_field(3,faults(1)%nglob))
        !do i = 1,faults(1)%nglob
        !  iglob = faults(1)%ibulk1(i)
        !  tmp_field(:,i) = tmp_displ(:,iglob)
        !enddo
        !call write_VTK_data_points_field(NGLOB_AB,xstore,ystore,zstore,faults(1)%ibulk1,faults(1)%nglob,tmp_field,filename)
        !print *,'written file: ',trim(filename)//'.vtk'
        !deallocate(tmp_field)

        if (SIMULATION_TYPE_DYN) call bc_dynflt_set3d_all(accel,veloc,displ,p,ilevel)
        if (SIMULATION_TYPE_KIN) call bc_kinflt_set_all(accel,veloc,displ,p,ilevel)
      endif
    endif

    ! apply M^-1, and step velocity and displacement
    call lts_newmark_update(veloc_p,displ_p,accel,m,ilevel,num_p_level,deltat_lts_local)

    ! debugging
    ! loops over each process
    !do iproc = 0,NPROC-1
    !  if (myrank == iproc) then
    !    ! user output
    !    print *,'debug: lts global_step: rank:',myrank,'ilevel:',ilevel,'m:',m,'after newmark  ', &
    !           'max displ:',maxval(displ_p(:,:,ilevel)),'veloc:',maxval(veloc_p(:,:,ilevel)),'accel:',maxval(accel(:,:))
    !  endif
    !  call synchronize_all()
    !enddo

  enddo ! m steps

  ! debug
  if (DEBUG_VTK_OUTPUT) then
    call lts_debug_vtk_output(ilevel)
  endif

  ! debugging
  if (DEBUG_TIME_STEPPING) then
    ! standard step for testing
    call lts_debug_time_stepping_refstep(ilevel)
    ! compares reference with lts wavefield
    call lts_debug_time_stepping_compare(ilevel)
  endif

  end subroutine lts_global_step

!
!---------------------------------------------------------------------------------------------------
!

  subroutine lts_set_accel_zero(accel,ilevel)

  use constants, only: NDIM,CUSTOM_REAL

  use shared_parameters, only: GPU_MODE

  use specfem_par, only: NGLOB_AB

  ! GPU
  !use specfem_par, only: Mesh_pointer
  !use specfem_par_lts, only: num_p_level

  ! LTS
  use specfem_par_lts, only: p_level_iglob_start, p_level_iglob_end, &
       num_p_level_coarser_to_update, p_level_coarser_to_update

  implicit none

  real(kind=CUSTOM_REAL),intent(inout) :: accel(NDIM,NGLOB_AB)
  integer,intent(in) :: ilevel

  ! local parameters
  integer :: is, ie, iglob_n, iglob

  ! sets accel to zero in appropriate places
  if (.not. GPU_MODE) then
    ! on CPU
    ! start index from finest level
    !is = p_level_iglob_start(1)
    ! start index of current level
    is = p_level_iglob_start(ilevel)
    ! end index of current level
    ie = p_level_iglob_end(ilevel)

    !debug
    !print *,'debug: lts accel zero: ilevel/is/ie',ilevel,is,ie,'num coarser',num_p_level_coarser_to_update(ilevel)

    ! checks
    if (ie < is) stop 'Error lts_set_accel_zero: invalid start/end index'

    ! zeroes out
    accel(:,is:ie) = 0.0_CUSTOM_REAL

    ! cancels out coarser node displacements
    do iglob_n = 1,num_p_level_coarser_to_update(ilevel)
      iglob = p_level_coarser_to_update(iglob_n,ilevel)
      accel(:,iglob) = 0.0_CUSTOM_REAL
    enddo
  else
    ! on GPU
    !#TODO: LTS on GPU compute_forces
    stop 'LTS on GPU w/ lts_set_accel_zero not implemented yet'
    !if (ilevel == num_p_level) then
    !  call lts_zero_accel_cuda(Mesh_pointer)
    !else if (ilevel < num_p_level) then
    !  call lts_set_accel_zero_fast_cuda(Mesh_pointer,ilevel,p_level_iglob_end,num_p_level_coarser_to_update)
    !endif
  endif

  end subroutine lts_set_accel_zero

!
!---------------------------------------------------------------------------------------------------
!

  subroutine lts_set_finer_initial_condition(displ_p,ilevel)

  use constants, only: NDIM,CUSTOM_REAL

  use shared_parameters, only: GPU_MODE

  use specfem_par, only: NGLOB_AB

  ! LTS
  use specfem_par_lts, only: p_level_iglob_start, p_level_iglob_end, &
    num_p_level_coarser_to_update, p_level_coarser_to_update, &
    num_p_level

  ! GPU
  !use specfem_par, only: Mesh_pointer

  implicit none

  real(kind=CUSTOM_REAL),intent(inout) :: displ_p(NDIM,NGLOB_AB,num_p_level)
  integer,intent(in) :: ilevel

  ! local parameters
  integer :: is, ie, iglob_n, iglob

  ! checks that current level is not coarsest one
  if (ilevel >= num_p_level) return

  ! copies over displacement
  ! from next coarser level to this current p-level, as starting displacement
  if (.not. GPU_MODE) then
    ! on CPU
    ! copies all from finest level up to end index of current level
    ! this is not working properly...
    !displ_p(:,:,ilevel) = displ_p(:,:,ilevel+1)
    !
    ! fast/partial memory copy of displacement array
    ! gets start index from finest level up to end index of current level
    is = p_level_iglob_start(1)
    ie = p_level_iglob_end(ilevel)

    ! checks
    if (ie < is) stop 'Error lts_set_finer_initial_condition: invalid start/end index'

    ! copies over displacement from coarser level
    displ_p(:,is:ie,ilevel) = displ_p(:,is:ie,ilevel+1)

    ! cancels out coarser node displacements
    do iglob_n = 1,num_p_level_coarser_to_update(ilevel)
      iglob = p_level_coarser_to_update(iglob_n,ilevel)
      displ_p(:,iglob,ilevel) = 0.0_CUSTOM_REAL
    enddo
  else
    ! on GPU
    !#TODO: LTS on GPU compute_forces
    stop 'LTS on GPU w/ set_finer_initial_condition not implemented yet'
    !call lts_set_finer_initial_condition_cuda(Mesh_pointer,ilevel,p_level_iglob_end,num_p_level_coarser_to_update)
  endif

  end subroutine lts_set_finer_initial_condition

!
!---------------------------------------------------------------------------------------------------
!

! the p-level boundary node code
!
! Handles contributions from nodes that are on the p-level boundary,
! including nodes that are within the current p-level, and those
! which are in both coarser and finer levels.

  subroutine lts_boundary_contribution(displ_p, veloc_p, accel, ilevel, &
                                       p, deltat_lts_local, &
                                       iphase, backward_simulation)

  use constants, only: myrank,CUSTOM_REAL,NDIM

  use specfem_par, only: NGLOB_AB

  use specfem_par_elastic, only: alphaval,betaval,gammaval, &
    R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
    R_trace_lddrk, &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
    epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
    epsilondev_xz,epsilondev_yz,epsilon_trace_over_3

  ! LTS
  use specfem_par_lts, only: p_level_ilevel_map, num_p_level, lts_type_compute_pelem, iglob_p_refine, &
    num_p_level_boundary_nodes,p_level_boundary_node

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB,num_p_level), intent(in) :: displ_p,veloc_p
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(inout) :: accel

  integer, intent(in) :: ilevel

  integer, intent(in) :: p
  real(kind=CUSTOM_REAL), intent(in) :: deltat_lts_local

  integer, intent(in) :: iphase
  logical, intent(in) :: backward_simulation

  ! local parameters
  integer :: iglob
  integer :: p_node, iilevel
  integer :: num_points,ipoin
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ_p_boundary
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: veloc_p_boundary
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel_p_boundary

  ! prepares boundary wavefield arrays

  ! number of global points on this level boundary
  num_points = num_p_level_boundary_nodes(iphase,ilevel)

  ! checks if anything to do
  if (num_points == 0) return

  ! loop only over boundary nodes
  do ipoin = 1,num_points
    iglob = p_level_boundary_node(ipoin,iphase,ilevel)

    ! checks if node belongs to this or coarser p-level
    p_node = iglob_p_refine(iglob)
    if (p_node > p) then
      print *,'Error: boundary node p value is invalid: ',p_node,'should be <= ',p
      print *,'       iglob ',iglob ,'iphase/ilevel',iphase,ilevel !,'i/j/k/ispec',i,j,k,ispec
      print *,'       p_level_ilevel_map = ',p_level_ilevel_map(p_node)
      call exit_mpi(myrank,"Error: Assert(boundary nodes are this level or coarser only)")
    endif

    ! takes corresponding displacement from this and coarser p-levels
    iilevel = p_level_ilevel_map(p_node)

    ! initial displacement
    displ_p_boundary(:,iglob) = displ_p(:,iglob,iilevel)
    ! initial velocity
    veloc_p_boundary(:,iglob) = veloc_p(:,iglob,iilevel)
    ! initial accel - zeroed out
    accel_p_boundary(:,iglob) = 0.0_CUSTOM_REAL
  enddo

  ! sets compute force type
  lts_type_compute_pelem = .false.   ! boundary elements only

  ! computes acceleration for boundary elements
  call compute_forces_viscoelastic(iphase,deltat_lts_local, &
                                   displ_p_boundary,veloc_p_boundary,accel_p_boundary, &
                                   alphaval,betaval,gammaval, &
                                   R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                   R_trace_lddrk, &
                                   R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                   epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                   epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                                   backward_simulation)


  ! fills element values back into "global" array

  ! loop only over boundary nodes
  do ipoin = 1,num_points
    iglob = p_level_boundary_node(ipoin,iphase,ilevel)

    ! adds boundary contribution
    accel(1,iglob) = accel(1,iglob) + accel_p_boundary(1,iglob)
    accel(2,iglob) = accel(2,iglob) + accel_p_boundary(2,iglob)
    accel(3,iglob) = accel(3,iglob) + accel_p_boundary(3,iglob)
  enddo

  end subroutine lts_boundary_contribution

!
!---------------------------------------------------------------------------------------------------
!

  subroutine lts_debug_time_stepping_refstep(ilevel)

! standard step for testing

  use constants
  use specfem_par
  use specfem_par_elastic
  use specfem_par_lts

  implicit none
  integer, intent(in) :: ilevel

  ! local parameters
  integer :: p,m,ier
  real(kind=CUSTOM_REAL) :: deltat_lts_local

  ! allocates reference wavefields
  if (.not. allocated(displ_ref)) then
    ! debug fields
    allocate(displ_ref(NDIM,NGLOB_AB), &
             veloc_ref(NDIM,NGLOB_AB), &
             accel_ref(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating working LTS fields displ_ref,..')
    ! initializes reference wavefields
    displ_ref(:,:) = 0.0_CUSTOM_REAL
    veloc_ref(:,:) = 0.0_CUSTOM_REAL
    accel_ref(:,:) = 0.0_CUSTOM_REAL
  endif

  ! p_elem and boundary_elem arrays for global reference routine
  if (.not. allocated(p_elem_ref)) then
    allocate(p_elem_ref(NSPEC_AB),boundary_elem_ref(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating working LTS fields p_elem_ref,..')
    ! all elements belong to same lts level
    p_elem_ref(:) = .true.
    boundary_elem_ref(:) = .false.
  endif

  ! only steps p-times when coarsest level is reached
  if (ilevel == num_p_level) then
    p = p_level(1)
    deltat_lts_local = real( DT / dble(p), kind=CUSTOM_REAL)
    ! for p smallest steps
    do m = 1,p
      ! Newmark predictor step
      !!displ_ref(:,:) = displ_ref(:,:) + deltat_lts_local * veloc_ref(:,:) + (deltat_lts_local**2)/2.d0 * accel_ref(:,:)
      !!veloc_ref(:,:) = veloc_ref(:,:) + deltat_lts_local/2.d0 * accel_ref(:,:)
      accel_ref(:,:) = 0.0_CUSTOM_REAL

      ! calling routine for viscoelastic forces
      call global_step_ref(displ_ref,veloc_ref,accel_ref,deltat_lts_local,m)

      ! Newmark corrector step
      !!veloc_ref(:,:) = veloc_ref(:,:) + deltat_lts_local/2.d0 * accel_ref(:,:)

      ! note: displ_p corresponds to time n+1 after the global stepping as it is updated based on accel_n+1.
      !       the LTS scheme is slightly different compared to a default Newmark scheme and has only a single step,
      !       not a predictor and corrector step.
      !
      !       with the predictor/corrector Newmark scheme, we use a (non-staggered) form like
      !         u_n+1 = u_n + dt v_n + 1/2 dt**2 a_n
      !         a_n+1 = B u_n+1 + F_n+1
      !         v_n+1 = v_n + 1/2 dt (a_n + a_n+1)
      !
      !       putting u_n+1=.. and v_n+1/2 = v_n + 1/2 dt (a_n) into the predictor step, and
      !       a_n+1=.. and v_n+1 = v_n+1/2 + 1/2 dt (a_n+1) into the corrector.
      !
      !       thus, in the predictor step, veloc is first updated with a term (dt/2 * accel_n), and then with (dt/2 * accel_n+1)
      !       in the corrector step. outputting velocity after this corrector step, means we output a staggered velocity v_1/2,..
      !       since in the very first step, the initial accel_0 == 0.
      !
      !       the LTS scheme updates displ/veloc/accel to n+1 after the global step.
      !       velocity is still at staggered times, v_-1/2,v_1/2,v_3/2,.., displ_p at d_0,d_1,d_2,...
      !       however, displ_p is updated after the global step immediately, still in the same iteration loop,
      !       whereas in the default Newmark (two step) scheme, it us updated when entering a new iteration in the call to
      !       update_displ_Newmark(). therefore, the write_seismograms() will output in the very first loop d_1 equal to zero
      !       since its updated with v_-1/2 and a_0 which all are zero.
      !       here below we update the reference displ_ref to displ_n+1, and veloc to veloc_n/2+1, according to the LTS scheme.

      !! updates displ_ref to compare against displ_p
      !! single step scheme
      !! displ update can be done either before or after veloc update:
      !!displ_ref(:,:) = displ_ref(:,:) + deltat_lts_local * veloc_ref(:,:) + deltat_lts_local**2 * accel_ref(:,:)
      veloc_ref(:,:) = veloc_ref(:,:) + deltat_lts_local * accel_ref(:,:)
      displ_ref(:,:) = displ_ref(:,:) + deltat_lts_local * veloc_ref(:,:)

      !! two-step version?
      !!veloc_ref(:,:) = veloc_ref(:,:) + DT/2.0 * accel_ref(:,:)
      !!displ_ref(:,:) = displ_ref(:,:) + DT/2.0 * veloc_ref(:,:)
    enddo
  endif

  end subroutine lts_debug_time_stepping_refstep

!
!---------------------------------------------------------------------------------------------------
!

  subroutine lts_debug_time_stepping_compare(ilevel)

! compares reference solution against lts step for testing

  use specfem_par
  use specfem_par_elastic
  use specfem_par_lts

  implicit none
  integer, intent(in) :: ilevel

  ! local parameters
  integer :: iproc
  character(len=MAX_STRING_LEN) :: filename
  real(kind=CUSTOM_REAL) :: max_error

  if (ilevel == num_p_level) then
    ! loops over each process
    do iproc = 0,NPROC-1
      if (myrank == iproc) then
        ! maximum difference between reference solution and local-time stepping
        max_error = maxval(sqrt((displ_p(1,:,ilevel)-displ_ref(1,:))**2 &
                              + (displ_p(2,:,ilevel)-displ_ref(2,:))**2 &
                              + (displ_p(3,:,ilevel)-displ_ref(3,:))**2))
        ! user output
        print *,"debug: lts global_step: rank:",myrank,"it",it,"max error = ",max_error, &
                !maxval(sqrt((displ_p(1,:,ilevel)-displ_ref(1,:))**2)), &
                !maxval(sqrt((displ_p(2,:,ilevel)-displ_ref(2,:))**2)), &
                !maxval(sqrt((displ_p(3,:,ilevel)-displ_ref(3,:))**2)), &
                " max lts:",maxval(displ_p(:,:,ilevel)), &
                " max ref:",maxval(displ_ref)

        ! single components
        !print *,"debug: lts global_step: rank:",myrank,"it",it,"max error x/y/z = ", &
        !        maxval(sqrt((displ_p(1,:,ilevel)-displ_ref(1,:))**2)), &
        !        maxval(sqrt((displ_p(2,:,ilevel)-displ_ref(2,:))**2)), &
        !        maxval(sqrt((displ_p(3,:,ilevel)-displ_ref(3,:))**2))

        ! wavefield snapshot
        if (.false.) then
          write(filename,'(a,i6.6,a,i6.6)') '../OUTPUT_FILES/proc',myrank,'_displ_it',it
          call write_VTK_wavefield(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool, &
                                   displ_p(:,:,ilevel),filename)
          print *,'written file: ',trim(filename)//'.vtk'
        endif
      endif
      call synchronize_all()
    enddo

    !debug
    !!to output seismograms w/ reference wavefield instead
    !!displ_p(:,:,ilevel) = displ_ref(:,:)
    !!veloc_p(:,:,ilevel) = veloc_ref(:,:)
    !!accel(:,:) = accel_ref(:,:)
  endif

  end subroutine lts_debug_time_stepping_compare

!
!---------------------------------------------------------------------------------------------------
!

  subroutine lts_debug_vtk_output(ilevel)

! standard step for testing

  use specfem_par
  use specfem_par_elastic
  use specfem_par_lts

  implicit none
  integer, intent(in) :: ilevel

  ! local parameters
  integer :: iproc
  character(len=MAX_STRING_LEN) :: filename

  if (ilevel == num_p_level) then
    ! loops over each process
    do iproc = 0,NPROC-1
      if (myrank == iproc) then
        ! wavefield snapshot
        if (.false. .and. mod(it,5) == 0) then
          ! displ
          write(filename,'(a,i6.6,a,i6.6)') '../OUTPUT_FILES/proc',myrank,'_displ_it',it
          call write_VTK_wavefield(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool, &
                                   displ,filename)
          print *,'written file: ',trim(filename)//'.vtk'

          ! veloc
          write(filename,'(a,i6.6,a,i6.6)') '../OUTPUT_FILES/proc',myrank,'_veloc_it',it
          call write_VTK_wavefield(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool, &
                                   veloc,filename)
          print *,'written file: ',trim(filename)//'.vtk'

          ! accel
          write(filename,'(a,i6.6,a,i6.6)') '../OUTPUT_FILES/proc',myrank,'_accel_it',it
          call write_VTK_wavefield(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool, &
                                   accel,filename)
          print *,'written file: ',trim(filename)//'.vtk'

        endif
      endif
      call synchronize_all()
    enddo
  endif

  end subroutine lts_debug_vtk_output

!
!---------------------------------------------------------------------------------------------------
!

  subroutine global_step_ref(displ_in,veloc_in,accel_in,deltat_lts_local,m_step)

! a non-LTS global step for reference

  use specfem_par
  use specfem_par_elastic

  use specfem_par_lts, only: current_lts_elem,current_lts_boundary_elem,lts_type_compute_pelem, &
    current_lts_time,p_elem_ref,boundary_elem_ref,rmassxyz_mod

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: displ_in,veloc_in
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel_in
  real(kind=CUSTOM_REAL), intent(in) :: deltat_lts_local
  integer, intent(in) :: m_step

  ! local parameters
  integer :: iphase
  logical :: backward_simulation

  ! safety check
  ! only for forward simulations at the moment
  if (SIMULATION_TYPE /= 1) stop 'LTS only implemented for forward simulation so far.'
  backward_simulation = .false.

  ! sets current p-level elements
  current_lts_elem => p_elem_ref(:)

  ! sets current boundary elements
  current_lts_boundary_elem => boundary_elem_ref(:)

  do iphase = 1,2

    ! sets compute force type
    lts_type_compute_pelem = .true.   ! p-elements only (no boundary elements)

    call compute_forces_viscoelastic(iphase,deltat_lts_local, &
                                     displ_in,veloc_in,accel_in, &
                                     alphaval,betaval,gammaval, &
                                     R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                     R_trace_lddrk, &
                                     R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                     epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                     epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                                     backward_simulation)

    ! computes additional contributions
    if (iphase == 1) then
      ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_viscoelastic_forward(NSPEC_AB,NGLOB_AB,accel_in, &
                                                 ibool,iphase, &
                                                 abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                 abs_boundary_ijk,abs_boundary_ispec, &
                                                 num_abs_boundary_faces, &
                                                 veloc_in,rho_vp,rho_vs, &
                                                 ispec_is_elastic, &
                                                 it,b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)
      endif

      ! current time
      ! note: this routine is called p times with m = 1,..,p and the smallest lts time step (DT/p)
      !       to match a global time stepping scheme)
      !
      ! current (local) time (starting at -t0)
      current_lts_time = dble(it-1)*DT - t0
      ! adds time increment for current "local" m step
      current_lts_time = current_lts_time + dble(m_step-1) * deltat_lts_local

      ! adds source contributions
      call compute_add_sources_viscoelastic(accel_in)
    endif

    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
      ! sends accel values to corresponding MPI interface neighbors
      call assemble_MPI_vector_async_send(NPROC,NGLOB_AB,accel_in, &
                                          buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          my_neighbors_ext_mesh, &
                                          request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
    else
      ! waits for send/receive requests to be completed and assembles values
      ! receives MPI buffers
      call assemble_MPI_vector_async_recv(NPROC,NGLOB_AB,accel_in, &
                                          buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                          max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                          my_neighbors_ext_mesh)
    endif

  enddo

  ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  accel_in(:,:) = accel_in(:,:) * rmassxyz_mod(:,:)

  end subroutine global_step_ref
