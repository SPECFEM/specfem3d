!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
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

! elastic solver

subroutine compute_forces_elastic()

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

! elastic term
    if(USE_DEVILLE_PRODUCTS) then
      call compute_forces_elastic_Dev(iphase, NSPEC_AB,NGLOB_AB,displ,accel, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        ATTENUATION, &
                        one_minus_sum_beta,factor_common, &
                        alphaval,betaval,gammaval, &
                        NSPEC_ATTENUATION_AB, &
                        R_xx,R_yy,R_xy,R_xz,R_yz, &
                        epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                        epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store,&
                        c22store,c23store,c24store,c25store,c26store,c33store,&
                        c34store,c35store,c36store,c44store,c45store,c46store,&
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE, COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
                        NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ADJOINT,&
                        is_moho_top,is_moho_bot, &
                        dsdx_top,dsdx_bot, &
                        ispec2D_moho_top,ispec2D_moho_bot, &
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic,&
                        phase_ispec_inner_elastic )
    else
      call compute_forces_elastic_noDev( iphase, NSPEC_AB,NGLOB_AB,displ,accel, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        ATTENUATION,&
                        one_minus_sum_beta,factor_common, &
                        alphaval,betaval,gammaval,&
                        NSPEC_ATTENUATION_AB, &
                        R_xx,R_yy,R_xy,R_xz,R_yz, &
                        epsilondev_xx,epsilondev_yy,epsilondev_xy,&
                        epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store,&
                        c22store,c23store,c24store,c25store,c26store,c33store,&
                        c34store,c35store,c36store,c44store,c45store,c46store,&
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE,COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
                        NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ADJOINT,&
                        is_moho_top,is_moho_bot, &
                        dsdx_top,dsdx_bot, &
                        ispec2D_moho_top,ispec2D_moho_bot, &
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic,&
                        phase_ispec_inner_elastic  )
    endif

    ! adjoint simulations: backward/reconstructed wavefield
    if( SIMULATION_TYPE == 3 ) then
      if(USE_DEVILLE_PRODUCTS) then
        call compute_forces_elastic_Dev(iphase, NSPEC_AB,NGLOB_AB, &
                        b_displ,b_accel, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        ATTENUATION, &
                        one_minus_sum_beta,factor_common, &
                        b_alphaval,b_betaval,b_gammaval, &
                        NSPEC_ATTENUATION_AB, &
                        b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz, &
                        b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy, &
                        b_epsilondev_xz,b_epsilondev_yz,b_epsilon_trace_over_3, &
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store,&
                        c22store,c23store,c24store,c25store,c26store,c33store,&
                        c34store,c35store,c36store,c44store,c45store,c46store,&
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE, COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY,&
                        NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ADJOINT,&
                        is_moho_top,is_moho_bot, &
                        b_dsdx_top,b_dsdx_bot, &
                        ispec2D_moho_top,ispec2D_moho_bot, &
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic,&
                        phase_ispec_inner_elastic )
      else
        call compute_forces_elastic_noDev( iphase, NSPEC_AB,NGLOB_AB,&
                        b_displ,b_accel, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        ATTENUATION,&
                        one_minus_sum_beta,factor_common, &
                        b_alphaval,b_betaval,b_gammaval,&
                        NSPEC_ATTENUATION_AB, &
                        b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz, &
                        b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,&
                        b_epsilondev_xz,b_epsilondev_yz,b_epsilon_trace_over_3, &
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store,&
                        c22store,c23store,c24store,c25store,c26store,c33store,&
                        c34store,c35store,c36store,c44store,c45store,c46store,&
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE,COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
                        NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ADJOINT,&
                        is_moho_top,is_moho_bot, &
                        b_dsdx_top,b_dsdx_bot, &
                        ispec2D_moho_top,ispec2D_moho_bot, &
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic,&
                        phase_ispec_inner_elastic  )

      endif
    endif



! adds elastic absorbing boundary term to acceleration (Stacey conditions)
    if(ABSORBING_CONDITIONS) &
      call compute_stacey_elastic(NSPEC_AB,NGLOB_AB,accel, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        veloc,rho_vp,rho_vs, &
                        ispec_is_elastic,SIMULATION_TYPE,myrank,SAVE_FORWARD, &
                        NSTEP,it,NGLOB_ADJOINT,b_accel, &
                        b_num_abs_boundary_faces,b_reclen_field,b_absorb_field )

! acoustic coupling
    if( ACOUSTIC_SIMULATION ) then
      call compute_coupling_elastic_ac(NSPEC_AB,NGLOB_AB, &
                        ibool,accel,potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) &
        call compute_coupling_elastic_ac(NSPEC_ADJOINT,NGLOB_ADJOINT, &
                        ibool,b_accel,b_potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)
    endif


! poroelastic coupling
! not implemented yet
!    if( POROELASTIC_SIMULATION ) &
!      call compute_coupling_elastic_poro()

! adds source term (single-force/moment-tensor solution)
    call compute_add_sources_elastic( NSPEC_AB,NGLOB_AB,accel, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                        xi_source,eta_source,gamma_source,nu_source, &
                        hdur,hdur_gaussian,t_cmt,dt,t0,sourcearrays, &
                        ispec_is_elastic,SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                        nrec,islice_selected_rec,ispec_selected_rec, &
                        nadj_rec_local,adj_sourcearrays,b_accel, &
                        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY )

! assemble all the contributions between slices using MPI
    if( phase_is_inner .eqv. .false. ) then
      ! sends accel values to corresponding MPI interface neighbors
      call assemble_MPI_vector_ext_mesh_s(NPROC,NGLOB_AB,accel, &
                        buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        call assemble_MPI_vector_ext_mesh_s(NPROC,NGLOB_ADJOINT,b_accel, &
                        b_buffer_send_vector_ext_mesh,b_buffer_recv_vector_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh)
      endif !adjoint

    else
      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_vector_ext_mesh_w(NPROC,NGLOB_AB,accel, &
                        buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        call assemble_MPI_vector_ext_mesh_w(NPROC,NGLOB_ADJOINT,b_accel, &
                        b_buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh)
      endif !adjoint

    endif

    !! DK DK May 2009: removed this because now each slice of a CUBIT + SCOTCH mesh
    !! DK DK May 2009: has a different number of spectral elements and therefore
    !! DK DK May 2009: only the general non-blocking MPI routines assemble_MPI_vector_ext_mesh_s
    !! DK DK May 2009: and assemble_MPI_vector_ext_mesh_w above can be used.
    !! DK DK May 2009: For adjoint runs below (SIMULATION_TYPE == 3) they should be used as well.

  enddo

! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  accel(1,:) = accel(1,:)*rmass(:)
  accel(2,:) = accel(2,:)*rmass(:)
  accel(3,:) = accel(3,:)*rmass(:)

  ! adjoint simulations
  if (SIMULATION_TYPE == 3) then
    b_accel(1,:) = b_accel(1,:)*rmass(:)
    b_accel(2,:) = b_accel(2,:)*rmass(:)
    b_accel(3,:) = b_accel(3,:)*rmass(:)
  endif !adjoint


! updates acceleration with ocean load term
  if(OCEANS) then
    call elastic_ocean_load(NSPEC_AB,NGLOB_AB, &
                        ibool,rmass,rmass_ocean_load,accel, &
                        free_surface_normal,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,SIMULATION_TYPE, &
                        NGLOB_ADJOINT,b_accel)
  endif

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
  veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)

  ! adjoint simulations
  if (SIMULATION_TYPE == 3) b_veloc(:,:) = b_veloc(:,:) + b_deltatover2*b_accel(:,:)


end subroutine compute_forces_elastic


!
!-------------------------------------------------------------------------------------------------
!

subroutine elastic_ocean_load(NSPEC_AB,NGLOB_AB, &
                        ibool,rmass,rmass_ocean_load,accel, &
                        free_surface_normal,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,SIMULATION_TYPE, &
                        NGLOB_ADJOINT,b_accel)

! updates acceleration with ocean load term:
! approximates ocean-bottom continuity of pressure & displacement for longer period waves (> ~20s ),
! assuming incompressible fluid column above bathymetry ocean bottom

  implicit none

  include 'constants.h'

  integer :: NSPEC_AB,NGLOB_AB

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB),intent(inout) :: accel
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: rmass,rmass_ocean_load

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! free surface
  integer :: num_free_surface_faces
  real(kind=CUSTOM_REAL) :: free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ispec(num_free_surface_faces)

  ! adjoint simulations
  integer :: SIMULATION_TYPE,NGLOB_ADJOINT
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_ADJOINT):: b_accel

! local parameters
  real(kind=CUSTOM_REAL) :: nx,ny,nz
  real(kind=CUSTOM_REAL) :: additional_term,force_normal_comp
  integer :: i,j,k,ispec,iglob
  integer :: igll,iface
  logical,dimension(NGLOB_AB) :: updated_dof_ocean_load
  ! adjoint locals
  real(kind=CUSTOM_REAL) :: b_additional_term,b_force_normal_comp

  !   initialize the updates
  updated_dof_ocean_load(:) = .false.

  ! for surface elements exactly at the top of the model (ocean bottom)
  do iface = 1,num_free_surface_faces

    ispec = free_surface_ispec(iface)
    do igll = 1, NGLLSQUARE
      i = free_surface_ijk(1,igll,iface)
      j = free_surface_ijk(2,igll,iface)
      k = free_surface_ijk(3,igll,iface)

      ! get global point number
      iglob = ibool(i,j,k,ispec)

      ! only update once
      if(.not. updated_dof_ocean_load(iglob)) then

        ! get normal
        nx = free_surface_normal(1,igll,iface)
        ny = free_surface_normal(2,igll,iface)
        nz = free_surface_normal(3,igll,iface)

        ! make updated component of right-hand side
        ! we divide by rmass() which is 1 / M
        ! we use the total force which includes the Coriolis term above
        force_normal_comp = ( accel(1,iglob)*nx + &
                              accel(2,iglob)*ny + &
                              accel(3,iglob)*nz ) / rmass(iglob)

        additional_term = (rmass_ocean_load(iglob) - rmass(iglob)) * force_normal_comp

        accel(1,iglob) = accel(1,iglob) + additional_term * nx
        accel(2,iglob) = accel(2,iglob) + additional_term * ny
        accel(3,iglob) = accel(3,iglob) + additional_term * nz

        ! adjoint simulations
        if (SIMULATION_TYPE == 3) then
          b_force_normal_comp = ( b_accel(1,iglob)*nx + &
                                  b_accel(2,iglob)*ny + &
                                  b_accel(3,iglob)*nz) / rmass(iglob)
          b_additional_term = (rmass_ocean_load(iglob) - rmass(iglob)) * b_force_normal_comp

          b_accel(1,iglob) = b_accel(1,iglob) + b_additional_term * nx
          b_accel(2,iglob) = b_accel(2,iglob) + b_additional_term * ny
          b_accel(3,iglob) = b_accel(3,iglob) + b_additional_term * nz
        endif !adjoint

        ! done with this point
        updated_dof_ocean_load(iglob) = .true.

      endif

    enddo ! igll
  enddo ! iface

end subroutine elastic_ocean_load

