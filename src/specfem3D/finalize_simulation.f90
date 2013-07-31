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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine finalize_simulation()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use pml_par
  use gravity_perturbation, only : gravity_output, GRAVITY_SIMULATION

  implicit none

  integer :: irec_local

  ! write gravity perturbations
  if (GRAVITY_SIMULATION) call gravity_output()

  ! save last frame

  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
     open(unit=27,file=prname(1:len_trim(prname))//'save_forward_arrays.bin',&
          status='unknown',form='unformatted')

    if( ACOUSTIC_SIMULATION ) then
      write(27) potential_acoustic
      write(27) potential_dot_acoustic
      write(27) potential_dot_dot_acoustic
    endif

    if( ELASTIC_SIMULATION ) then
      write(27) displ
      write(27) veloc
      write(27) accel

      if (ATTENUATION) then
        if(FULL_ATTENUATION_SOLID) write(27) R_trace  !ZN
        write(27) R_xx
        write(27) R_yy
        write(27) R_xy
        write(27) R_xz
        write(27) R_yz
        if(FULL_ATTENUATION_SOLID) write(27) epsilondev_trace !ZN
        write(27) epsilondev_xx
        write(27) epsilondev_yy
        write(27) epsilondev_xy
        write(27) epsilondev_xz
        write(27) epsilondev_yz
      endif
    endif

    if( POROELASTIC_SIMULATION ) then
      write(27) displs_poroelastic
      write(27) velocs_poroelastic
      write(27) accels_poroelastic
      write(27) displw_poroelastic
      write(27) velocw_poroelastic
      write(27) accelw_poroelastic
    endif

    close(27)

! adjoint simulations
  else if (SIMULATION_TYPE == 3) then

    ! adjoint kernels
    call save_adjoint_kernels()

  endif

! closing source time function file
  if(PRINT_SOURCE_TIME_FUNCTION .and. myrank == 0) then
    close(IOSTF)
  endif

! stacey absorbing fields will be reconstructed for adjoint simulations
! using snapshot files of wavefields
  if( STACEY_ABSORBING_CONDITIONS ) then
    ! closes absorbing wavefield saved/to-be-saved by forward simulations
    if( num_abs_boundary_faces > 0 .and. (SIMULATION_TYPE == 3 .or. &
          (SIMULATION_TYPE == 1 .and. SAVE_FORWARD)) ) then

      if( ELASTIC_SIMULATION) call close_file_abs(0) ! close(IOABS)
      if( ACOUSTIC_SIMULATION) call close_file_abs(1) ! close(IOABS_AC)

    endif
  endif

! seismograms and source parameter gradients for (pure type=2) adjoint simulation runs
  if (nrec_local > 0) then
    if (.not. (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3)) then
      ! seismograms
      call write_adj_seismograms2_to_file(myrank,seismograms_eps,number_receiver_global, &
            nrec_local,it,DT,NSTEP,t0,LOCAL_PATH)

      ! source gradients  (for sources in elastic domains)
      do irec_local = 1, nrec_local
        write(outputname,'(a,i5.5)') OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)) // &
            '/src_frechet.',number_receiver_global(irec_local)
        open(unit=27,file=trim(outputname),status='unknown')
        !
        ! r -> z, theta -> -y, phi -> x
        !
        !  Mrr =  Mzz
        !  Mtt =  Myy
        !  Mpp =  Mxx
        !  Mrt = -Myz
        !  Mrp =  Mxz
        !  Mtp = -Mxy
        write(27,*) Mzz_der(irec_local)
        write(27,*) Myy_der(irec_local)
        write(27,*) Mxx_der(irec_local)
        write(27,*) -Myz_der(irec_local)
        write(27,*) Mxz_der(irec_local)
        write(27,*) -Mxy_der(irec_local)
        write(27,*) sloc_der(1,irec_local)
        write(27,*) sloc_der(2,irec_local)
        write(27,*) sloc_der(3,irec_local)
        close(27)
      enddo
    endif
  endif

  ! frees dynamically allocated memory

  if (USE_FORCE_POINT_SOURCE) then
     deallocate(factor_force_source)
     deallocate(comp_dir_vect_source_E)
     deallocate(comp_dir_vect_source_N)
     deallocate(comp_dir_vect_source_Z_UP)
     if (.not. USE_RICKER_TIME_FUNCTION) then
       deallocate(hdur_tiny)
    endif
  endif

  ! mass matrices
  if( ELASTIC_SIMULATION ) then
    deallocate(rmassx)
    deallocate(rmassy)
    deallocate(rmassz)
  endif
  if( ACOUSTIC_SIMULATION ) then
    deallocate(rmass_acoustic)
  endif

  ! C-PML absorbing boundary conditions
  if( PML_CONDITIONS .and. NSPEC_CPML > 0 ) then
     ! outputs informations about C-PML elements in VTK-file format
     call pml_output_VTKs()

     ! deallocates C_PML arrays
     deallocate(CPML_regions)
     deallocate(CPML_to_spec)
     deallocate(is_CPML)
     deallocate(d_store_x)
     deallocate(d_store_y)
     deallocate(d_store_z)
     deallocate(k_store_x)
     deallocate(k_store_y)
     deallocate(k_store_z)
     deallocate(alpha_store)
     deallocate(spec_to_CPML)
     deallocate(CPML_type)

     if( ELASTIC_SIMULATION ) then
       deallocate(PML_dux_dxl)
       deallocate(PML_dux_dyl)
       deallocate(PML_dux_dzl)
       deallocate(PML_duy_dxl)
       deallocate(PML_duy_dyl)
       deallocate(PML_duy_dzl)
       deallocate(PML_duz_dxl)
       deallocate(PML_duz_dyl)
       deallocate(PML_duz_dzl)
       deallocate(PML_dux_dxl_new)
       deallocate(PML_dux_dyl_new)
       deallocate(PML_dux_dzl_new)
       deallocate(PML_duy_dxl_new)
       deallocate(PML_duy_dyl_new)
       deallocate(PML_duy_dzl_new)
       deallocate(PML_duz_dxl_new)
       deallocate(PML_duz_dyl_new)
       deallocate(PML_duz_dzl_new)
       deallocate(rmemory_dux_dxl_x)
       deallocate(rmemory_dux_dyl_x)
       deallocate(rmemory_dux_dzl_x)
       deallocate(rmemory_duy_dxl_x)
       deallocate(rmemory_duy_dyl_x)
       deallocate(rmemory_duz_dxl_x)
       deallocate(rmemory_duz_dzl_x)
       deallocate(rmemory_dux_dxl_y)
       deallocate(rmemory_dux_dyl_y)
       deallocate(rmemory_duy_dxl_y)
       deallocate(rmemory_duy_dyl_y)
       deallocate(rmemory_duy_dzl_y)
       deallocate(rmemory_duz_dyl_y)
       deallocate(rmemory_duz_dzl_y)
       deallocate(rmemory_dux_dxl_z)
       deallocate(rmemory_dux_dzl_z)
       deallocate(rmemory_duy_dyl_z)
       deallocate(rmemory_duy_dzl_z)
       deallocate(rmemory_duz_dxl_z)
       deallocate(rmemory_duz_dyl_z)
       deallocate(rmemory_duz_dzl_z)
       deallocate(rmemory_displ_elastic)
       deallocate(accel_elastic_CPML)
     endif

     if( ACOUSTIC_SIMULATION ) then
       deallocate(PML_dpotential_dxl)
       deallocate(PML_dpotential_dyl)
       deallocate(PML_dpotential_dzl)
       deallocate(PML_dpotential_dxl_new)
       deallocate(PML_dpotential_dyl_new)
       deallocate(PML_dpotential_dzl_new)
       deallocate(rmemory_dpotential_dxl)
       deallocate(rmemory_dpotential_dyl)
       deallocate(rmemory_dpotential_dzl)
       deallocate(rmemory_potential_acoustic)
       deallocate(potential_dot_dot_acoustic_CPML)
     endif

  endif

  deallocate(ibelm_xmin)
  deallocate(ibelm_xmax)
  deallocate(ibelm_ymin)
  deallocate(ibelm_ymax)
  deallocate(ibelm_bottom)
  deallocate(ibelm_top)

! close the main output file
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'End of the simulation'
    write(IMAIN,*)
    close(IMAIN)
  endif

! synchronize all the processes to make sure everybody has finished
  call sync_all()

  end subroutine finalize_simulation
