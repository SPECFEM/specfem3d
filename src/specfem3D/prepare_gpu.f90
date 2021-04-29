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


  subroutine prepare_GPU()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_noise
  use specfem_par_movie

  use fault_solver_common, only: Kelvin_Voigt_eta,USE_KELVIN_VOIGT_DAMPING,Kelvin_Voigt_eta
  use fault_solver_dynamic, only: SIMULATION_TYPE_DYN,fault_transfer_data_GPU,fault_rsf_swf_init_GPU
  use fault_solver_kinematic, only: SIMULATION_TYPE_KIN

  implicit none

  ! local parameters
  real :: free_mb,used_mb,total_mb

  ! checks if anything to do
  if (.not. GPU_MODE) return

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) "preparing fields and constants on GPU devices"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! evaluates memory required
  call prepare_memory_eval_gpu()

  ! prepares general fields on GPU
  ! add GPU support for the C-PML routines
  call prepare_constants_device(Mesh_pointer, &
                                NGLLX, NSPEC_AB, NGLOB_AB, NSPEC_IRREGULAR,irregular_element_number, &
                                xixstore,xiystore,xizstore, &
                                etaxstore,etaystore,etazstore, &
                                gammaxstore,gammaystore,gammazstore, &
                                xix_regular,jacobian_regular,ibool, &
                                num_interfaces_ext_mesh, max_nibool_interfaces_ext_mesh, &
                                nibool_interfaces_ext_mesh, ibool_interfaces_ext_mesh, &
                                hprime_xx,hprimewgll_xx, &
                                wgllwgll_xy, wgllwgll_xz, wgllwgll_yz, &
                                STACEY_ABSORBING_CONDITIONS, &
                                abs_boundary_ispec, abs_boundary_ijk, &
                                abs_boundary_normal, &
                                abs_boundary_jacobian2Dw, &
                                num_abs_boundary_faces, &
                                ispec_is_inner, &
                                NSOURCES, nsources_local, &
                                sourcearrays, islice_selected_source, ispec_selected_source, &
                                ispec_selected_rec, &
                                nrec, nrec_local, &
                                SIMULATION_TYPE, &
                                USE_MESH_COLORING_GPU, &
                                nspec_acoustic,nspec_elastic, &
                                ispec_is_acoustic,ispec_is_elastic, &
                                myrank,SAVE_FORWARD, &
                                hxir_store,hetar_store,hgammar_store, &
                                nu_rec,nu_source, &
                                islice_selected_rec, &
                                nlength_seismogram, &
                                SAVE_SEISMOGRAMS_DISPLACEMENT,SAVE_SEISMOGRAMS_VELOCITY, &
                                SAVE_SEISMOGRAMS_ACCELERATION,SAVE_SEISMOGRAMS_PRESSURE, &
                                NB_RUNS_ACOUSTIC_GPU, &
                                FAULT_SIMULATION)


  ! prepares fields on GPU for acoustic simulations
  if (ACOUSTIC_SIMULATION) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  loading acoustic arrays"
      call flush_IMAIN()
    endif

    call prepare_fields_acoustic_device(Mesh_pointer, &
                                rmass_acoustic,rhostore,kappastore, &
                                num_phase_ispec_acoustic,phase_ispec_inner_acoustic, &
                                NOISE_TOMOGRAPHY,num_free_surface_faces, &
                                free_surface_ispec,free_surface_ijk, &
                                b_reclen_potential,b_absorb_potential, &
                                ELASTIC_SIMULATION, num_coupling_ac_el_faces, &
                                coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                coupling_ac_el_normal,coupling_ac_el_jacobian2Dw, &
                                num_colors_outer_acoustic,num_colors_inner_acoustic, &
                                num_elem_colors_acoustic)

    if (SIMULATION_TYPE == 3) &
      call prepare_fields_acoustic_adj_dev(Mesh_pointer,APPROXIMATE_HESS_KL)

  endif

  ! prepares fields on GPU for elastic simulations
  ! add GPU support for the C-PML routines
  if (ELASTIC_SIMULATION) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  loading elastic arrays"
      call flush_IMAIN()
    endif

    call prepare_fields_elastic_device(Mesh_pointer, &
                                rmassx,rmassy,rmassz, &
                                rho_vp,rho_vs, &
                                kappastore, mustore, &
                                num_phase_ispec_elastic,phase_ispec_inner_elastic, &
                                b_absorb_field,b_reclen_field, &
                                COMPUTE_AND_STORE_STRAIN, &
                                epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                epsilondev_xz,epsilondev_yz, &
                                ATTENUATION, &
                                size(R_xx), &
                                R_xx,R_yy,R_xy,R_xz,R_yz, &
                                factor_common, &
                                R_trace,epsilondev_trace, &
                                factor_common_kappa, &
                                alphaval,betaval,gammaval, &
                                APPROXIMATE_OCEAN_LOAD,rmass_ocean_load, &
                                NOISE_TOMOGRAPHY, &
                                free_surface_normal,free_surface_ispec,free_surface_ijk, &
                                num_free_surface_faces, &
                                ACOUSTIC_SIMULATION, &
                                num_colors_outer_elastic,num_colors_inner_elastic, &
                                num_elem_colors_elastic, &
                                ANISOTROPY, &
                                c11store,c12store,c13store,c14store,c15store,c16store, &
                                c22store,c23store,c24store,c25store,c26store, &
                                c33store,c34store,c35store,c36store, &
                                c44store,c45store,c46store,c55store,c56store,c66store)

    if (SIMULATION_TYPE == 3) &
      call prepare_fields_elastic_adj_dev(Mesh_pointer, &
                                NDIM*NGLOB_AB, &
                                COMPUTE_AND_STORE_STRAIN, &
                                epsilon_trace_over_3, &
                                b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy, &
                                b_epsilondev_xz,b_epsilondev_yz, &
                                b_epsilon_trace_over_3, &
                                ATTENUATION,size(R_xx), &
                                b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz, &
                                b_R_trace,b_epsilondev_trace, &
                                b_alphaval,b_betaval,b_gammaval, &
                                ANISOTROPIC_KL, &
                                APPROXIMATE_HESS_KL)
  endif

  ! prepares fields on GPU for poroelastic simulations
  if (POROELASTIC_SIMULATION) then
    stop 'Poroelastic simulations on GPU not implemented yet'
  endif

  ! prepares needed receiver array for adjoint runs
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  loading adjoint receiver arrays"
      call flush_IMAIN()
    endif
    call prepare_sim2_or_3_const_device(Mesh_pointer,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC, &
                                        hxir_adjstore,hetar_adjstore,hgammar_adjstore, &
                                        nrec,islice_selected_rec,ispec_selected_rec)
  endif

  ! prepares fields on GPU for noise simulations
  if (NOISE_TOMOGRAPHY > 0) then
    ! note: noise tomography is only supported for elastic domains so far.
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  loading noise arrays"
      call flush_IMAIN()
    endif
    ! copies noise  arrays to GPU
    call prepare_fields_noise_device(Mesh_pointer, &
                                NSPEC_AB, NGLOB_AB, &
                                free_surface_ispec, &
                                free_surface_ijk, &
                                num_free_surface_faces, &
                                NOISE_TOMOGRAPHY, &
                                NSTEP,noise_sourcearray, &
                                normal_x_noise,normal_y_noise,normal_z_noise, &
                                mask_noise,free_surface_jacobian2Dw)

  endif ! NOISE_TOMOGRAPHY

  ! prepares gravity arrays
  if (GRAVITY) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  loading gravity"
      call flush_IMAIN()
    endif
    call prepare_fields_gravity_device(Mesh_pointer,GRAVITY, &
                                       minus_deriv_gravity,minus_g,wgll_cube, &
                                       ACOUSTIC_SIMULATION,rhostore)
  endif

  ! prepares fault rupture simulation
  if (FAULT_SIMULATION) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  loading fault simulation"
      call flush_IMAIN()
    endif

    ! dynamic rupture
    if (SIMULATION_TYPE_DYN) then
      ! user output
      if (myrank == 0) then
        write(IMAIN,*) "    dynamic rupture arrays"
        call flush_IMAIN()
      endif
      ! initializes fault data on gpu
      call fault_transfer_data_GPU()

      ! sets up Kelvin-Voigt damping array
      if (USE_KELVIN_VOIGT_DAMPING) then
        call prepare_fault_device(Mesh_pointer,USE_KELVIN_VOIGT_DAMPING,Kelvin_Voigt_eta)
      endif

      ! sets up friction law arrays
      call fault_rsf_swf_init_GPU()
    endif

    ! kinematic rupture
    if (SIMULATION_TYPE_KIN) then
      stop 'Kinematic rupture simulations on GPU not implemented yet'
    endif
  endif

  ! synchronizes processes
  call synchronize_all()

  ! transfer forward and backward fields to device with initial values
  if (myrank == 0) then
    write(IMAIN,*) "  transferring initial wavefield"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! sends initial data to device

  ! puts acoustic initial fields onto GPU
  if (ACOUSTIC_SIMULATION) then
    call transfer_fields_ac_to_device(NGLOB_AB,potential_acoustic, &
                                      potential_dot_acoustic,potential_dot_dot_acoustic,Mesh_pointer)
    if (SIMULATION_TYPE == 3) &
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic,b_potential_dot_dot_acoustic,Mesh_pointer)
  endif

  ! puts elastic initial fields onto GPU
  if (ELASTIC_SIMULATION) then
    ! transfer forward and backward fields to device with initial values
    call transfer_fields_el_to_device(NDIM*NGLOB_AB,displ,veloc,accel,Mesh_pointer)
    if (SIMULATION_TYPE == 3) &
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)
  endif

  ! warnings
  ! obsolete...
  !if (myrank == 0) then
  !  if (SAVE_SEISMOGRAMS_DISPLACEMENT .or. SAVE_SEISMOGRAMS_VELOCITY .or. SAVE_SEISMOGRAMS_ACCELERATION) then
  !    ! warnings
  !    if (.not. ELASTIC_SIMULATION) &
  !      print *, "Warning: Wrong type of seismogram for a pure fluid simulation, use pressure in seismotype"
  !    if (ELASTIC_SIMULATION .and. ACOUSTIC_SIMULATION) &
  !      print *, "Warning: Coupled elastic/fluid simulation has only valid displacement seismograms &
  !              &in elastic domain for GPU simulation"
  !  endif
  !  if (SAVE_SEISMOGRAMS_PRESSURE) then
  !    if (.not. ACOUSTIC_SIMULATION) &
  !      print *, "Warning: Wrong type of seismogram for a pure elastic simulation, use displ veloc or accel in seismotype"
  !    if (ELASTIC_SIMULATION .and. ACOUSTIC_SIMULATION) &
  !      print *, "Warning: Coupled elastic/fluid simulation has only valid pressure seismograms &
  !               &in fluid domain for GPU simulation"
  !  endif
  !endif

  ! synchronizes processes
  call synchronize_all()

  ! outputs GPU usage to files for all processes
  call output_free_device_memory(myrank)

  ! outputs usage for main process
  if (myrank == 0) then
    call get_free_device_memory(free_mb,used_mb,total_mb)
    write(IMAIN,*)
    write(IMAIN,*) "GPU usage: free  =",free_mb," MB",nint(free_mb/total_mb*100.0),"%"
    write(IMAIN,*) "           used  =",used_mb," MB",nint(used_mb/total_mb*100.0),"%"
    write(IMAIN,*) "           total =",total_mb," MB",nint(total_mb/total_mb*100.0),"%"
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  end subroutine prepare_GPU


!
!----------------------------------------------------------------------------------------------
!

  subroutine prepare_memory_eval_gpu()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  use fault_solver_common, only: USE_KELVIN_VOIGT_DAMPING
  use fault_solver_dynamic, only: SIMULATION_TYPE_DYN

  implicit none

  ! local parameters
  double precision :: memory_size,memory_size_glob
  integer,parameter :: NGLL2 = 25
  integer,parameter :: NGLL3 = 125
  integer,parameter :: NGLL3_PADDED = 128

  memory_size = 0.d0

  ! add size of each set of arrays multiplied by the number of such arrays
  ! d_hprime_xx,d_hprimewgll_xx
  memory_size = memory_size + 2.d0 * NGLL2 * dble(CUSTOM_REAL)
  ! padded xix,..gammaz
  memory_size = memory_size + 9.d0 * NGLL3_PADDED * NSPEC_IRREGULAR * dble(CUSTOM_REAL)
  ! ibool + irregular_element_number
  memory_size = memory_size + (NGLL3+1) * NSPEC_AB * dble(SIZE_INTEGER)
  ! d_ibool_interfaces_ext_mesh
  memory_size = memory_size + num_interfaces_ext_mesh * max_nibool_interfaces_ext_mesh * dble(SIZE_INTEGER)
  ! ispec_is_inner
  memory_size = memory_size + NSPEC_AB * dble(SIZE_INTEGER)
  ! d_ispec_is_acoustic
  memory_size = memory_size + NSPEC_AB * dble(SIZE_INTEGER)
  ! d_ispec_is_elastic
  memory_size = memory_size + NSPEC_AB * dble(SIZE_INTEGER)

  if (STACEY_ABSORBING_CONDITIONS) then
    ! d_abs_boundary_ispec
    memory_size = memory_size + num_abs_boundary_faces * dble(SIZE_INTEGER)
    ! d_abs_boundary_ijk
    memory_size = memory_size + num_abs_boundary_faces * 3.d0 * NGLL2 * dble(SIZE_INTEGER)
    ! d_abs_boundary_normal
    memory_size = memory_size + num_abs_boundary_faces * NDIM * NGLL2 * dble(CUSTOM_REAL)
    ! d_abs_boundary_jacobian2Dw
    memory_size = memory_size + num_abs_boundary_faces * NGLL2 * dble(CUSTOM_REAL)
  endif

  ! sources
  ! d_sourcearrays
  memory_size = memory_size + NGLL3 * NSOURCES * NDIM * dble(CUSTOM_REAL)
  ! d_islice_selected_source,d_ispec_selected_source
  memory_size = memory_size + 2.0 * NSOURCES * dble(SIZE_INTEGER)

  ! receivers
  ! d_ispec_selected_rec_loc
  memory_size = memory_size + nrec_local * dble(SIZE_INTEGER)
  ! d_ispec_selected_rec
  memory_size = memory_size + nrec * dble(SIZE_INTEGER)

  ! d_hxir, d_hetar, d_hgammar
  memory_size = memory_size + 3.d0 * NGLLX * nrec_local * dble(CUSTOM_REAL)
  ! d_nu
  memory_size = memory_size + NDIM * NDIM * nrec_local * dble(CUSTOM_REAL)
  ! d_ispec_selected_rec_loc
  memory_size = memory_size + nrec_local * dble(SIZE_INTEGER)

  if (SIMULATION_TYPE == 2) then
    ! d_hxir_adj, d_hetar_adj, d_hgammar_adj
    memory_size = memory_size + 3.d0 * NGLLX * nadj_rec_local * dble(CUSTOM_REAL)
  endif

  ! d_seismograms_d,d_seismograms_v,d_seismograms_a,d_seismograms_p
  if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
    memory_size = memory_size + NDIM * nlength_seismogram * nrec_local * dble(CUSTOM_REAL)
  if (SAVE_SEISMOGRAMS_VELOCITY) &
    memory_size = memory_size + NDIM * nlength_seismogram * nrec_local * dble(CUSTOM_REAL)
  if (SAVE_SEISMOGRAMS_ACCELERATION) &
    memory_size = memory_size + NDIM * nlength_seismogram * nrec_local * dble(CUSTOM_REAL)
  if (SAVE_SEISMOGRAMS_PRESSURE) &
    memory_size = memory_size + nlength_seismogram * nrec_local * dble(CUSTOM_REAL) * NB_RUNS_ACOUSTIC_GPU

  ! acoustic simulations
  if (ACOUSTIC_SIMULATION) then
    ! d_potential_acoustic,d_potential_dot_acoustic,d_potential_dot_dot_acoustic
    memory_size = memory_size + 3.d0 * NGLOB_AB * dble(CUSTOM_REAL) * NB_RUNS_ACOUSTIC_GPU
    ! d_rmass_acoustic
    memory_size = memory_size + NGLOB_AB * dble(CUSTOM_REAL)
    ! padded d_rhostore
    memory_size = memory_size + NGLL3_PADDED * NSPEC_AB * dble(CUSTOM_REAL)
    ! d_kappastore
    memory_size = memory_size + NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
    ! d_phase_ispec_inner_acoustic
    memory_size = memory_size + 2.d0 * num_phase_ispec_acoustic * dble(SIZE_INTEGER)
  endif

  ! elastic simulations
  if (ELASTIC_SIMULATION) then
    ! d_displ,..
    memory_size = memory_size + 3.d0 * NDIM * NGLOB_AB * dble(CUSTOM_REAL)
    ! d_send_accel_buffer
    memory_size = memory_size + 3.d0 * num_interfaces_ext_mesh * max_nibool_interfaces_ext_mesh * dble(CUSTOM_REAL)
    ! d_rmassx,..
    memory_size = memory_size + 3.d0 * NGLOB_AB * dble(CUSTOM_REAL)
    ! d_phase_ispec_inner_elastic
    memory_size = memory_size + 2.d0 * num_phase_ispec_elastic * dble(SIZE_INTEGER)

    if (STACEY_ABSORBING_CONDITIONS .or. PML_CONDITIONS) then
      ! d_rho_vp,..
      memory_size = memory_size + 2.d0 * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
    endif

    ! padded kappav,muv
    memory_size = memory_size + 2.d0 * NGLL3_PADDED * NSPEC_AB * dble(CUSTOM_REAL)

    if (COMPUTE_AND_STORE_STRAIN) then
      ! d_epsilondev_xx,..
      memory_size = memory_size + 6.d0 * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
    endif
    if (ATTENUATION) then
      ! d_R_xx,..
      memory_size = memory_size + 6.d0 * size(R_xx) * dble(CUSTOM_REAL)
      ! d_factor_common, d_factor_common_kappa
      memory_size = memory_size + 2.d0 * N_SLS * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      ! alphaval,..
      memory_size = memory_size + 3.d0 * N_SLS * dble(CUSTOM_REAL)
    endif
    if (ANISOTROPY) then
      ! padded d_c11store,..
      memory_size = memory_size + 21.d0 * NGLL3_PADDED * NSPEC_AB * dble(CUSTOM_REAL)
    endif
    if (APPROXIMATE_OCEAN_LOAD) then
      ! d_rmass_ocean_load
      memory_size = memory_size + NGLOB_AB * dble(CUSTOM_REAL)
      ! d_free_surface_normal
      memory_size = memory_size + 3.d0 * NGLL2 * num_free_surface_faces * dble(CUSTOM_REAL)
    endif
  endif

  ! noise simulations
  if (NOISE_TOMOGRAPHY > 0) then
    ! d_free_surface_ispec
    memory_size = memory_size + num_free_surface_faces * dble(SIZE_INTEGER)
    ! d_free_surface_ijk
    memory_size = memory_size + 3.d0 * NGLL2 * num_free_surface_faces * dble(SIZE_INTEGER)
    ! d_noise_surface_movie
    memory_size = memory_size + 3.d0 * NGLL2 * num_free_surface_faces * dble(CUSTOM_REAL)
    if (NOISE_TOMOGRAPHY == 1) then
      ! d_noise_sourcearray
      memory_size = memory_size + 3.d0 * NGLL3 * NSTEP * dble(CUSTOM_REAL)
    endif
    if (NOISE_TOMOGRAPHY > 1) then
      ! d_normal_x_noise,..
      memory_size = memory_size + 5.d0 * NGLL2 * num_free_surface_faces * dble(CUSTOM_REAL)
    endif
    if (NOISE_TOMOGRAPHY == 3) then
      ! d_sigma_kl
      memory_size = memory_size + NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
    endif
  endif

  if (GRAVITY) then
    ! d_minus_deriv_gravity,d_minus_g
    memory_size = memory_size + 2.d0 * NGLOB_AB * dble(CUSTOM_REAL)
  endif

  if (FAULT_SIMULATION) then
    if (SIMULATION_TYPE_DYN) then
      ! arrays like Flt->B,R,Z,D,V,.. are not considered here yet.
      ! they would need the number of fault points for each fault
      ! memory_size = memory_size + 21.d0 * NGLOB_FLT * dble(CUSTOM_REAL)
      if (USE_KELVIN_VOIGT_DAMPING) then
        ! d_Kelvin_Voigt_eta
        memory_size = memory_size + NSPEC_AB * dble(CUSTOM_REAL)
      endif
    endif
  endif

  ! poor estimate for kernel simulations...
  if (SIMULATION_TYPE == 3) memory_size = 2.d0 * memory_size

  ! maximum of all processes (may vary e.g. due to different nrec_local, num_abs_boundary_faces, ..)
  call max_all_dp(memory_size,memory_size_glob)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  minimum memory requested     : ',memory_size_glob / 1024. / 1024.,'MB per process'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine prepare_memory_eval_gpu
