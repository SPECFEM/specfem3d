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


    if( .NOT. GPU_MODE ) then
! solid phase

    call compute_forces_poro_solid_part( iphase, &
                        NSPEC_AB,NGLOB_AB,displs_poroelastic,accels_poroelastic,&
                        displw_poroelastic,velocw_poroelastic,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz,&
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wxgll,wygll,wzgll,  &
                        kappaarraystore,rhoarraystore,mustore,etastore,permstore, &
                        phistore,tortstore,jacobian,ibool,&
                        epsilonsdev_xx,epsilonsdev_yy,epsilonsdev_xy,&
                        epsilonsdev_xz,epsilonsdev_yz,epsilons_trace_over_3, &
                        SIMULATION_TYPE,NSPEC_ADJOINT, &
                        num_phase_ispec_poroelastic,nspec_inner_poroelastic,nspec_outer_poroelastic,&
                        phase_ispec_inner_poroelastic )

! fluid phase

    call compute_forces_poro_fluid_part( iphase, &
                        NSPEC_AB,NGLOB_AB,displw_poroelastic,accelw_poroelastic,&
                        velocw_poroelastic,displs_poroelastic,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz,&
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wxgll,wygll,wzgll,  &
                        kappaarraystore,rhoarraystore,mustore,etastore,permstore, &
                        phistore,tortstore,jacobian,ibool,&
                        epsilonwdev_xx,epsilonwdev_yy,epsilonwdev_xy,&
                        epsilonwdev_xz,epsilonwdev_yz,epsilonw_trace_over_3, &
                        SIMULATION_TYPE,NSPEC_ADJOINT, &
                        num_phase_ispec_poroelastic,nspec_inner_poroelastic,nspec_outer_poroelastic,&
                        phase_ispec_inner_poroelastic )

    ! adjoint simulations: backward/reconstructed wavefield
    if( SIMULATION_TYPE == 3 ) then
! solid phase

    call compute_forces_poro_solid_part( iphase, &
                        NSPEC_AB,NGLOB_AB,b_displs_poroelastic,b_accels_poroelastic,&
                        b_displw_poroelastic,b_velocw_poroelastic,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz,&
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wxgll,wygll,wzgll,  &
                        kappaarraystore,rhoarraystore,mustore,etastore,permstore, &
                        phistore,tortstore,jacobian,ibool,&
                        b_epsilonsdev_xx,b_epsilonsdev_yy,b_epsilonsdev_xy,&
                        b_epsilonsdev_xz,b_epsilonsdev_yz,b_epsilons_trace_over_3, &
                        SIMULATION_TYPE,NSPEC_ADJOINT, &
                        num_phase_ispec_poroelastic,nspec_inner_poroelastic,nspec_outer_poroelastic,&
                        phase_ispec_inner_poroelastic )

! fluid phase

    call compute_forces_poro_fluid_part( iphase, &
                        NSPEC_AB,NGLOB_AB,b_displw_poroelastic,b_accelw_poroelastic,&
                        b_velocw_poroelastic,b_displs_poroelastic,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz,&
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wxgll,wygll,wzgll,  &
                        kappaarraystore,rhoarraystore,mustore,etastore,permstore, &
                        phistore,tortstore,jacobian,ibool,&
                        b_epsilonwdev_xx,b_epsilonwdev_yy,b_epsilonwdev_xy,&
                        b_epsilonwdev_xz,b_epsilonwdev_yz,b_epsilonw_trace_over_3, &
                        SIMULATION_TYPE,NSPEC_ADJOINT, &
                        num_phase_ispec_poroelastic,nspec_inner_poroelastic,nspec_outer_poroelastic,&
                        phase_ispec_inner_poroelastic )
    endif

    else
      ! on GPU
      call exit_MPI(myrank,'GPU for poroelastic simulation not implemented')
    endif ! GPU_MODE

! adds poroelastic absorbing boundary terms to accelerations (type Stacey conditions)
    if(STACEY_ABSORBING_CONDITIONS) &
      call compute_stacey_poroelastic(NSPEC_AB,NGLOB_AB,accels_poroelastic,accelw_poroelastic, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        velocs_poroelastic,velocw_poroelastic,rho_vpI,rho_vpII,rho_vsI, &
                        rhoarraystore,phistore,tortstore, &
                        ispec_is_poroelastic,SIMULATION_TYPE,SAVE_FORWARD, &
                        NSTEP,it,NGLOB_ADJOINT,b_accels_poroelastic,b_accelw_poroelastic, &
                        b_num_abs_boundary_faces,b_reclen_field_poro,b_absorb_fields, &
                        b_absorb_fieldw)

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
! chris:'adjoint acoustic-poroelastic simulation not implemented yet'
!        call ccmpute_coupling_elastic_ac(NSPEC_ADJOINT,NGLOB_ADJOINT, &
!                        ibool,b_accel,b_potential_dot_dot_acoustic, &
!                        num_coupling_ac_el_faces, &
!                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
!                        coupling_ac_el_normal, &
!                        coupling_ac_el_jacobian2Dw, &
!                        ispec_is_inner,phase_is_inner)
    endif

! elastic coupling
    if( ELASTIC_SIMULATION ) then
      call compute_coupling_poroelastic_el(NSPEC_AB,NGLOB_AB,ibool,&
                        displs_poroelastic,accels_poroelastic,displw_poroelastic,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz,&
                        kappaarraystore,rhoarraystore,mustore, &
                        phistore,tortstore,jacobian,&
                        displ,kappastore, &
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store,&
                        c22store,c23store,c24store,c25store,c26store,c33store,&
                        c34store,c35store,c36store,c44store,c45store,c46store,&
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE,NGLOB_ADJOINT,NSPEC_ADJOINT, &
                        num_coupling_el_po_faces, &
                        coupling_el_po_ispec,coupling_po_el_ispec,&
                        coupling_el_po_ijk,coupling_po_el_ijk, &
                        coupling_el_po_normal, &
                        coupling_el_po_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)

      ! adjoint simulations
! chris:'adjoint elastic-poroelastic simulation not implemented yet'
!      if( SIMULATION_TYPE == 3 ) &
!        call compute_coupling_viscoelastic_ac(NSPEC_ADJOINT,NGLOB_ADJOINT, &
!                        ibool,b_accel,b_potential_dot_dot_acoustic, &
!                        num_coupling_ac_el_faces, &
!                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
!                        coupling_ac_el_normal, &
!                        coupling_ac_el_jacobian2Dw, &
!                        ispec_is_inner,phase_is_inner)
    endif

! adds source term (single-force/moment-tensor solution)
    call compute_add_sources_poroelastic( NSPEC_AB,NGLOB_AB, &
                        accels_poroelastic,accelw_poroelastic,&
                        rhoarraystore,phistore,tortstore,&
                        ibool,ispec_is_inner,phase_is_inner, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                        hdur,hdur_gaussian,tshift_src,dt,t0,sourcearrays, &
                        ispec_is_poroelastic,SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                        nrec,islice_selected_rec,ispec_selected_rec, &
                        nadj_rec_local,adj_sourcearrays,b_accels_poroelastic,b_accelw_poroelastic, &
                        NTSTEP_BETWEEN_READ_ADJSRC)

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
      call assemble_MPI_vector_poro_s(NPROC,NGLOB_ADJOINT,b_accels_poroelastic, &
                        b_accelw_poroelastic,&
                        b_buffer_send_vector_ext_meshs,b_buffer_recv_vector_ext_meshs, &
                        b_buffer_send_vector_ext_meshw,b_buffer_recv_vector_ext_meshw, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        b_request_send_vector_ext_meshs,b_request_recv_vector_ext_meshs, &
                        b_request_send_vector_ext_meshw,b_request_recv_vector_ext_meshw)
      endif !adjoint

    else
      ! waits for send/receive requests to be completed and assembles values
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
      call assemble_MPI_vector_poro_w(NPROC,NGLOB_ADJOINT,b_accels_poroelastic, &
                        b_accelw_poroelastic,&
                        b_buffer_recv_vector_ext_meshs,b_buffer_recv_vector_ext_meshw, &
                        num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        b_request_send_vector_ext_meshs,b_request_recv_vector_ext_meshs, &
                        b_request_send_vector_ext_meshw,b_request_recv_vector_ext_meshw)
      endif !adjoint

    endif
  enddo

! solid phase
! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  accels_poroelastic(1,:) = accels_poroelastic(1,:)*rmass_solid_poroelastic(:)
  accels_poroelastic(2,:) = accels_poroelastic(2,:)*rmass_solid_poroelastic(:)
  accels_poroelastic(3,:) = accels_poroelastic(3,:)*rmass_solid_poroelastic(:)

  ! adjoint simulations
  if (SIMULATION_TYPE == 3) then
    b_accels_poroelastic(1,:) = b_accels_poroelastic(1,:)*rmass_solid_poroelastic(:)
    b_accels_poroelastic(2,:) = b_accels_poroelastic(2,:)*rmass_solid_poroelastic(:)
    b_accels_poroelastic(3,:) = b_accels_poroelastic(3,:)*rmass_solid_poroelastic(:)
  endif !adjoint

! fluid phase
! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  accelw_poroelastic(1,:) = accelw_poroelastic(1,:)*rmass_fluid_poroelastic(:)
  accelw_poroelastic(2,:) = accelw_poroelastic(2,:)*rmass_fluid_poroelastic(:)
  accelw_poroelastic(3,:) = accelw_poroelastic(3,:)*rmass_fluid_poroelastic(:)

  ! adjoint simulations
  if (SIMULATION_TYPE == 3) then
    b_accelw_poroelastic(1,:) = b_accelw_poroelastic(1,:)*rmass_fluid_poroelastic(:)
    b_accelw_poroelastic(2,:) = b_accelw_poroelastic(2,:)*rmass_fluid_poroelastic(:)
    b_accelw_poroelastic(3,:) = b_accelw_poroelastic(3,:)*rmass_fluid_poroelastic(:)
  endif !adjoint

! updates velocities
! Newmark finite-difference time scheme with elastic domains:
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
  if (SIMULATION_TYPE == 3) b_velocs_poroelastic(:,:) = b_velocs_poroelastic(:,:) + &
                               b_deltatover2*b_accels_poroelastic(:,:)

! fluid phase
  if (SIMULATION_TYPE == 3) b_velocw_poroelastic(:,:) = b_velocw_poroelastic(:,:) + &
                               b_deltatover2*b_accelw_poroelastic(:,:)

! elastic coupling
    if( ELASTIC_SIMULATION ) &
      call compute_continuity_disp_po_el(NSPEC_AB,NGLOB_AB,ibool,&
                        accel,veloc,&
                        accels_poroelastic,velocs_poroelastic,&
                        accelw_poroelastic,velocw_poroelastic,&
                        rmass,rmass_solid_poroelastic,&
                        SIMULATION_TYPE,NSPEC_ADJOINT,&
                        num_coupling_el_po_faces,&
                        coupling_el_po_ispec,coupling_el_po_ijk,&
                        deltatover2)


end subroutine compute_forces_poroelastic


!-----------------------------------------------------------------------------------------------------------------

subroutine compute_continuity_disp_po_el(NSPEC_AB,NGLOB_AB,ibool,&
                        accel,veloc,&
                        accels_poroelastic,velocs_poroelastic,&
                        accelw_poroelastic,velocw_poroelastic,&
                        rmass,rmass_solid_poroelastic,&
                        SIMULATION_TYPE,NSPEC_ADJOINT,&
                        num_coupling_el_po_faces,&
                        coupling_el_po_ispec,coupling_el_po_ijk,&
                        deltatover2)
!*******************************************************************************
!         assembling the displacements on the elastic-poro boundaries
!*******************************************************************************

  implicit none
  include 'constants.h'

  integer :: NSPEC_AB,NGLOB_AB

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: veloc,accel
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: velocs_poroelastic,accels_poroelastic
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: velocw_poroelastic,accelw_poroelastic

  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: rmass,rmass_solid_poroelastic

! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  integer :: SIMULATION_TYPE
  integer :: NSPEC_ADJOINT

! elastic-poroelastic coupling surface
  integer :: num_coupling_el_po_faces
  integer :: coupling_el_po_ijk(3,NGLLSQUARE,num_coupling_el_po_faces)
  integer :: coupling_el_po_ispec(num_coupling_el_po_faces)

  real(kind=CUSTOM_REAL),intent(in) :: deltatover2

! local variables
  integer, dimension(NGLOB_AB) :: icount
  integer :: ispec,i,j,k,l,iglob,igll,iface

     icount(:)=ZERO

! loops on all coupling faces
  do iface = 1,num_coupling_el_po_faces

    ! gets corresponding poroelastic spectral element
    ispec = coupling_el_po_ispec(iface)

      ! loops over common GLL points
      do igll = 1, NGLLSQUARE
        i = coupling_el_po_ijk(1,igll,iface)
        j = coupling_el_po_ijk(2,igll,iface)
        k = coupling_el_po_ijk(3,igll,iface)

        ! gets global index of this common GLL point
        iglob = ibool(i,j,k,ispec)
        icount(iglob) = icount(iglob) + 1

        if(icount(iglob) ==1)then

! recovering original velocities and accelerations on boundaries (elastic side)
          veloc(1,iglob) = veloc(1,iglob) - deltatover2*accel(1,iglob)
          veloc(2,iglob) = veloc(2,iglob) - deltatover2*accel(2,iglob)
          veloc(3,iglob) = veloc(3,iglob) - deltatover2*accel(3,iglob)
          accel(1,iglob) = accel(1,iglob) / rmass(iglob)
          accel(2,iglob) = accel(2,iglob) / rmass(iglob)
          accel(3,iglob) = accel(3,iglob) / rmass(iglob)
! recovering original velocities and accelerations on boundaries (poro side)
          velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) - deltatover2*accels_poroelastic(1,iglob)
          velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) - deltatover2*accels_poroelastic(2,iglob)
          velocs_poroelastic(3,iglob) = velocs_poroelastic(3,iglob) - deltatover2*accels_poroelastic(3,iglob)
          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_solid_poroelastic(iglob)
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_solid_poroelastic(iglob)
          accels_poroelastic(3,iglob) = accels_poroelastic(3,iglob) / rmass_solid_poroelastic(iglob)
! assembling accelerations
          accel(1,iglob) = ( accel(1,iglob) + accels_poroelastic(1,iglob) ) / &
                                   ( 1.0/rmass(iglob) +1.0/rmass_solid_poroelastic(iglob) )
          accel(2,iglob) = ( accel(2,iglob) + accels_poroelastic(2,iglob) ) / &
                                   ( 1.0/rmass(iglob) +1.0/rmass_solid_poroelastic(iglob) )
          accel(3,iglob) = ( accel(3,iglob) + accels_poroelastic(3,iglob) ) / &
                                   ( 1.0/rmass(iglob) +1.0/rmass_solid_poroelastic(iglob) )
          accels_poroelastic(1,iglob) = accel(1,iglob)
          accels_poroelastic(2,iglob) = accel(2,iglob)
          accels_poroelastic(3,iglob) = accel(3,iglob)
! updating velocities
          velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) + deltatover2*accels_poroelastic(1,iglob)
          velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) + deltatover2*accels_poroelastic(2,iglob)
          velocs_poroelastic(3,iglob) = velocs_poroelastic(3,iglob) + deltatover2*accels_poroelastic(3,iglob)
          veloc(1,iglob) = veloc(1,iglob) + deltatover2*accel(1,iglob)
          veloc(2,iglob) = veloc(2,iglob) + deltatover2*accel(2,iglob)
          veloc(3,iglob) = veloc(3,iglob) + deltatover2*accel(3,iglob)
! zeros w
          accelw_poroelastic(1,iglob) = 0.d0
          accelw_poroelastic(2,iglob) = 0.d0
          accelw_poroelastic(3,iglob) = 0.d0
          velocw_poroelastic(1,iglob) = 0.d0
          velocw_poroelastic(2,iglob) = 0.d0
          velocw_poroelastic(3,iglob) = 0.d0

   if(SIMULATION_TYPE == 3) then
!chris: to do
   l = NSPEC_ADJOINT ! to avoid compilation warnings
   endif

        endif !if(icount(iglob) ==1)
     enddo ! igll
  enddo ! iface

end subroutine compute_continuity_disp_po_el
