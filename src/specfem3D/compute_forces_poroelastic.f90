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

! poroelastic solver

subroutine compute_forces_poroelastic()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  integer:: iphase
  logical:: phase_is_inner

! distinguishes two runs: for points on MPI interfaces, and points within the partitions
  do iphase=1,2

    !first for points on MPI interfaces
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

!Note: Contrary to the elastic & acoustic case, the absorbing implementation is within compute_forces
!due to the number of properties which were needed

! solid phase

    call compute_forces_solid( iphase, &
                        NSPEC_AB,NGLOB_AB,displs_poroelastic,accels_poroelastic,&
                        velocs_poroelastic,displw_poroelastic,velocw_poroelastic,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz,&
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wxgll,wygll,wzgll,  &
                        kappaarraystore,rhoarraystore,mustore,etastore,permstore, &
                        phistore,tortstore,jacobian,ibool,&
                        SIMULATION_TYPE,NGLOB_ADJOINT,NSPEC_ADJOINT, &
                        num_phase_ispec_poroelastic,nspec_inner_poroelastic,nspec_outer_poroelastic,&
                        phase_ispec_inner_poroelastic )

! fluid phase

    call compute_forces_fluid( iphase, &
                        NSPEC_AB,NGLOB_AB,displw_poroelastic,accelw_poroelastic,&
                        velocw_poroelastic,displs_poroelastic,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz,&
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wxgll,wygll,wzgll,  &
                        kappaarraystore,rhoarraystore,mustore,etastore,permstore, &
                        phistore,tortstore,jacobian,ibool,&
                        SIMULATION_TYPE,NGLOB_ADJOINT,NSPEC_ADJOINT, &
                        num_phase_ispec_poroelastic,nspec_inner_poroelastic,nspec_outer_poroelastic,&
                        phase_ispec_inner_poroelastic )


    ! adjoint simulations: backward/reconstructed wavefield
    if( SIMULATION_TYPE == 3 ) then
 stop 'adjoint poroelastic simulation not implemented yet'
    endif


! acoustic coupling
    if( ACOUSTIC_SIMULATION ) then
      call compute_coupling_poroelastic_ac(NSPEC_AB,NGLOB_AB, &
                        ibool,accels_poroelastic,accelw_poroelastic, &
                        potential_dot_dot_acoustic, &
                        num_coupling_ac_po_faces, &
                        coupling_ac_po_ispec,coupling_ac_po_ijk, &
                        coupling_ac_po_normal, &
                        coupling_ac_po_jacobian2Dw, &
                        rhoarraystore,phistore,tortstore, &
                        ispec_is_inner,phase_is_inner)

      ! adjoint simulations
      !if( SIMULATION_TYPE == 3 ) &
! 'adjoint acoustic-poroelastic simulation not implemented yet'
!        call compute_coupling_elastic_ac(NSPEC_ADJOINT,NGLOB_ADJOINT, &
!                        ibool,b_accel,b_potential_dot_dot_acoustic, &
!                        num_coupling_ac_el_faces, &
!                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
!                        coupling_ac_el_normal, &
!                        coupling_ac_el_jacobian2Dw, &
!                        ispec_is_inner,phase_is_inner)
    endif

! elastic coupling
!chris: TO DO
    if( ELASTIC_SIMULATION ) &
 stop 'elastic-poroelastic simulation not implemented yet'
!      call compute_coupling_poroelastic_el()

      ! adjoint simulations
!chris: TO DO
!      if( SIMULATION_TYPE == 3 ) &
!        call compute_coupling_elastic_ac(NSPEC_ADJOINT,NGLOB_ADJOINT, &
!                        ibool,b_accel,b_potential_dot_dot_acoustic, &
!                        num_coupling_ac_el_faces, &
!                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
!                        coupling_ac_el_normal, &
!                        coupling_ac_el_jacobian2Dw, &
!                        ispec_is_inner,phase_is_inner)
!    endif

! adds source term (single-force/moment-tensor solution)
    call compute_add_sources_poroelastic( NSPEC_AB,NGLOB_AB, &
                        accels_poroelastic,accelw_poroelastic,&
                        rhoarraystore,phistore,tortstore,&
                        ibool,ispec_is_inner,phase_is_inner, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                        xi_source,eta_source,gamma_source,nu_source, &
                        hdur,hdur_gaussian,tshift_cmt,dt,t0,sourcearrays, &
                        ispec_is_poroelastic,SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                        nrec,islice_selected_rec,ispec_selected_rec, &
                        nadj_rec_local,adj_sourcearrays)

! assemble all the contributions between slices using MPI
    if( phase_is_inner .eqv. .false. ) then
      ! sends accel values to corresponding MPI interface neighbors
      call assemble_MPI_vector_poro_s(NPROC,NGLOB_AB,accels_poroelastic, &
                        accelw_poroelastic,&
                        buffer_send_vector_ext_mesh_s,buffer_recv_vector_ext_mesh_s, &
                        buffer_send_vector_ext_mesh_w,buffer_recv_vector_ext_mesh_w, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_vector_ext_mesh_s,request_recv_vector_ext_mesh_s, &
                        request_send_vector_ext_mesh_w,request_recv_vector_ext_mesh_w)

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
! 'adjoint poroelastic simulation not implemented yet'
!        call assemble_MPI_vector_ext_mesh_s(NPROC,NGLOB_ADJOINT,b_accel, &
!                        b_buffer_send_vector_ext_mesh,b_buffer_recv_vector_ext_mesh, &
!                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
!                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
!                        my_neighbours_ext_mesh, &
!                        b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh)
      endif !adjoint

    else
      ! waits for send/receive requests to be completed and assembles values
! solid phase
      call assemble_MPI_vector_poro_w(NPROC,NGLOB_AB,accels_poroelastic, &
                        accelw_poroelastic,&
                        buffer_recv_vector_ext_mesh_s,buffer_recv_vector_ext_mesh_w, &
                        num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_vector_ext_mesh_s,request_recv_vector_ext_mesh_s, &
                        request_send_vector_ext_mesh_w,request_recv_vector_ext_mesh_w)

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
! 'adjoint poroelastic simulation not implemented yet'
!        call assemble_MPI_vector_ext_mesh_w(NPROC,NGLOB_ADJOINT,b_accel, &
!                        b_buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh,&
!                        max_nibool_interfaces_ext_mesh, &
!                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
!                        b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh)
      endif !adjoint

    endif

    !! DK DK May 2009: removed this because now each slice of a CUBIT + SCOTCH mesh
    !! DK DK May 2009: has a different number of spectral elements and therefore
    !! DK DK May 2009: only the general non-blocking MPI routines assemble_MPI_vector_ext_mesh_s
    !! DK DK May 2009: and assemble_MPI_vector_ext_mesh_w above can be used.
    !! DK DK May 2009: For adjoint runs below (SIMULATION_TYPE == 3) they should be used as well.

  enddo

! solid phase
! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  accels_poroelastic(1,:) = accels_poroelastic(1,:)*rmass_solid_poroelastic(:)
  accels_poroelastic(2,:) = accels_poroelastic(2,:)*rmass_solid_poroelastic(:)
  accels_poroelastic(3,:) = accels_poroelastic(3,:)*rmass_solid_poroelastic(:)

  ! adjoint simulations
  if (SIMULATION_TYPE == 3) then
! 'adjoint poroelastic simulation not implemented yet'
!    b_accels_poroelastic(1,:) = b_accels_poroelastic(1,:)*rmass_solid_poroelastic(:)
!    b_accels_poroelastic(2,:) = b_accels_poroelastic(2,:)*rmass_solid_poroelastic(:)
!    b_accels_poroelastic(3,:) = b_accels_poroelastic(3,:)*rmass_solid_poroelastic(:)
  endif !adjoint

! fluid phase
! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  accelw_poroelastic(1,:) = accelw_poroelastic(1,:)*rmass_fluid_poroelastic(:)
  accelw_poroelastic(2,:) = accelw_poroelastic(2,:)*rmass_fluid_poroelastic(:)
  accelw_poroelastic(3,:) = accelw_poroelastic(3,:)*rmass_fluid_poroelastic(:)

  ! adjoint simulations
  if (SIMULATION_TYPE == 3) then
! 'adjoint poroelastic simulation not implemented yet'
!    b_accelw_poroelastic(1,:) = b_accelw_poroelastic(1,:)*rmass_fluid_poroelastic(:)
!    b_accelw_poroelastic(2,:) = b_accelw_poroelastic(2,:)*rmass_fluid_poroelastic(:)
!    b_accelw_poroelastic(3,:) = b_accelw_poroelastic(3,:)*rmass_fluid_poroelastic(:)
  endif !adjoint

! updates velocities
! Newark finite-difference time scheme with elastic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t) )
!
! where
!   u, v, a are displacement,velocity & acceleration
!   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
!   f denotes a source term (acoustic/elastic)
!   chi_dot_dot is acoustic (fluid) potential ( dotted twice with respect to time)
!
! corrector:
!   updates the velocity term which requires a(t+delta)
! solid phase
  velocs_poroelastic(:,:) = velocs_poroelastic(:,:) + deltatover2*accels_poroelastic(:,:)

! fluid phase
  velocw_poroelastic(:,:) = velocw_poroelastic(:,:) + deltatover2*accelw_poroelastic(:,:)

  ! adjoint simulations
! solid phase
!  if (SIMULATION_TYPE == 3) b_velocs_poroelastic(:,:) = b_velocs_poroelastic(:,:) + &
!                               b_deltatover2*b_accels_poroelastic(:,:)

! fluid phase
!  if (SIMULATION_TYPE == 3) b_velocw_poroelastic(:,:) = b_velocw_poroelastic(:,:) + &
!                               b_deltatover2*b_accelw_poroelastic(:,:)

end subroutine compute_forces_poroelastic

