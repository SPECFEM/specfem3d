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
      call compute_forces_with_Deville(phase_is_inner, NSPEC_AB,NGLOB_AB,displ,accel, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool,ispec_is_inner, &
                        ATTENUATION,USE_OLSEN_ATTENUATION, &
                        one_minus_sum_beta,factor_common,alphaval,betaval,gammaval, &
                        NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
                        epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                        epsilondev_xz,epsilondev_yz,iflag_attenuation_store, &
                        rho_vs,ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store,&
                        c22store,c23store,c24store,c25store,c26store,c33store,&
                        c34store,c35store,c36store,c44store,c45store,c46store,&
                        c55store,c56store,c66store, &
                        ispec_is_elastic )
    else
      call compute_forces_no_Deville( phase_is_inner, NSPEC_AB,NGLOB_AB,displ,accel, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool,ispec_is_inner, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,&
                        one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,&
                        NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
                        epsilondev_xx,epsilondev_yy,epsilondev_xy,&
                        epsilondev_xz,epsilondev_yz,iflag_attenuation_store,&
                        rho_vs,ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store,&
                        c22store,c23store,c24store,c25store,c26store,c33store,&
                        c34store,c35store,c36store,c44store,c45store,c46store,&
                        c55store,c56store,c66store, &
                        ispec_is_elastic )
    endif

! adds elastic absorbing boundary term to acceleration (Stacey conditions)
    if(ABSORBING_CONDITIONS) &
      call elastic_absorbing_boundaries(NSPEC_AB,NGLOB_AB,accel, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        veloc,rho_vp,rho_vs, &
                        ispec_is_elastic )

! acoustic coupling
    if( ACOUSTIC_SIMULATION ) &
      call elastic_coupling_acoustic(NSPEC_AB,NGLOB_AB, &
                        ibool,accel,potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)

! poroelastic coupling
    if( POROELASTIC_SIMULATION ) &
      call elastic_coupling_poroelastic()

! adds source term (single-force/moment-tensor solution)
    call elastic_sources( NSPEC_AB,NGLOB_AB,accel, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                        xi_source,eta_source,gamma_source,nu_source, &
                        hdur,hdur_gaussian,t_cmt,dt,stf,t0,sourcearrays, &
                        ispec_is_elastic  )
    
! assemble all the contributions between slices using MPI
    if( phase_is_inner .eqv. .false. ) then 
      ! sends accel values to corresponding MPI interface neighbors  
      call assemble_MPI_vector_ext_mesh_s(NPROC,NGLOB_AB,accel, &
                        buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
    else
      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_vector_ext_mesh_w(NPROC,NGLOB_AB,accel, &
                        buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
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
  
  !! DK DK array not created yet for CUBIT
  ! if (SIMULATION_TYPE == 3) then
  !   b_accel(1,:) = b_accel(1,:)*rmass(:)
  !   b_accel(2,:) = b_accel(2,:)*rmass(:)
  !   b_accel(3,:) = b_accel(3,:)*rmass(:)
  ! endif


! updates acceleration with ocean load term
  if(OCEANS) then    
    call elastic_ocean_load(NSPEC_AB,NGLOB_AB, &
                        ibool,rmass,rmass_ocean_load,accel, &
                        free_surface_normal,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces)
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

  !! DK DK array not created yet for CUBIT
  ! if (SIMULATION_TYPE == 3) b_veloc(:,:) = b_veloc(:,:) + b_deltatover2*b_accel(:,:)


end subroutine compute_forces_elastic


!
!-------------------------------------------------------------------------------------------------
!

! absorbing boundary term for elastic media (Stacey conditions)

subroutine elastic_absorbing_boundaries(NSPEC_AB,NGLOB_AB,accel, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        veloc,rho_vp,rho_vs, &
                        ispec_is_elastic)

  implicit none

  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB

! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner
  
! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB) :: ispec_is_elastic

! absorbing boundary surface  
  integer :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces) 
  real(kind=CUSTOM_REAL) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces) 
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces) 


! local parameters
  real(kind=CUSTOM_REAL) vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw !weight,jacobianl
  integer :: ispec,iglob,i,j,k,iface,igll
  !integer :: num_gll !,igll_i,igll_j,ispec2D
  

! absorbs absorbing-boundary surface using Stacey condition (Clayton & Enquist)
  do iface=1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

      if( ispec_is_elastic(ispec) ) then
      
        ! reference gll points on boundary face 
        do igll = 1,NGLLSQUARE

          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          ! gets velocity
          iglob=ibool(i,j,k,ispec)
          vx=veloc(1,iglob)
          vy=veloc(2,iglob)
          vz=veloc(3,iglob)

          ! gets associated normal
          nx = abs_boundary_normal(1,igll,iface)
          ny = abs_boundary_normal(2,igll,iface)
          nz = abs_boundary_normal(3,igll,iface)             

          ! velocity component in normal direction (normal points out of element)
          vn = vx*nx + vy*ny + vz*nz
             
          ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it 
          tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
          ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
          tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

          ! gets associated, weighted jacobian 
          jacobianw = abs_boundary_jacobian2Dw(igll,iface)
          
          ! adds stacey term (weak form)
          accel(1,iglob) = accel(1,iglob) - tx*jacobianw
          accel(2,iglob) = accel(2,iglob) - ty*jacobianw
          accel(3,iglob) = accel(3,iglob) - tz*jacobianw

         enddo
      endif ! ispec_is_elastic
    endif ! ispec_is_inner    
  enddo
  
end subroutine elastic_absorbing_boundaries

!
!-------------------------------------------------------------------------------------------------
!

subroutine elastic_coupling_acoustic(NSPEC_AB,NGLOB_AB, &
                        ibool,accel,potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)

! returns the updated acceleration array: accel                        

  implicit none
  include 'constants.h'

  integer :: NSPEC_AB,NGLOB_AB

! displacement and pressure
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic
  
! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! acoustic-elastic coupling surface
  integer :: num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL) :: coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces) 
  real(kind=CUSTOM_REAL) :: coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces)
  integer :: coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces)
  integer :: coupling_ac_el_ispec(num_coupling_ac_el_faces)   

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! local parameters
  real(kind=CUSTOM_REAL) :: pressure
  real(kind=CUSTOM_REAL) :: nx,ny,nz,jacobianw
  
  integer :: iface,igll,ispec,iglob
  integer :: i,j,k
  
! loops on all coupling faces
  do iface = 1,num_coupling_ac_el_faces

    ! gets corresponding spectral element 
    ! (note: can be either acoustic or elastic element, no need to specify since
    !           no material properties are needed for this coupling term)
    ispec = coupling_ac_el_ispec(iface)

    if( ispec_is_inner(ispec) .eqv. phase_is_inner ) then
    
      ! loops over common GLL points
      do igll = 1, NGLLSQUARE
        i = coupling_ac_el_ijk(1,igll,iface)
        j = coupling_ac_el_ijk(2,igll,iface)
        k = coupling_ac_el_ijk(3,igll,iface)
        
        ! gets global index of this common GLL point
        ! (note: should be the same as for corresponding i',j',k',ispec_elastic or ispec_elastic )
        iglob = ibool(i,j,k,ispec)
        
        ! acoustic pressure on global point
        pressure = - potential_dot_dot_acoustic(iglob)

        ! gets associated normal on GLL point
        ! (note convention: pointing outwards of acoustic element)
        nx = coupling_ac_el_normal(1,igll,iface)
        ny = coupling_ac_el_normal(2,igll,iface)
        nz = coupling_ac_el_normal(3,igll,iface)                   
        
        ! gets associated, weighted 2D jacobian 
        ! (note: should be the same for elastic and acoustic element)
        jacobianw = coupling_ac_el_jacobian2Dw(igll,iface)
        
        ! continuity of displacement and pressure on global point
        !
        ! note: newark time scheme together with definition of scalar potential: 
        !          pressure = - chi_dot_dot
        !          requires that this coupling term uses the *UPDATED* pressure (chi_dot_dot), i.e.
        !          pressure at time step [t + delta_t] 
        !          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
        !          it means you have to calculate and update the acoustic pressure first before
        !          calculating this term...
        accel(1,iglob) = accel(1,iglob) + jacobianw*nx*pressure
        accel(2,iglob) = accel(2,iglob) + jacobianw*ny*pressure
        accel(3,iglob) = accel(3,iglob) + jacobianw*nz*pressure
        
      enddo ! igll

    endif
    
  enddo ! iface

end subroutine elastic_coupling_acoustic

!
!-------------------------------------------------------------------------------------------------
!

subroutine elastic_coupling_poroelastic()
  implicit none
 
end subroutine elastic_coupling_poroelastic

!
!-------------------------------------------------------------------------------------------------
!

subroutine elastic_sources( NSPEC_AB,NGLOB_AB,accel, &
                                  ibool,ispec_is_inner,phase_is_inner, &
                                  NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                                  xi_source,eta_source,gamma_source,nu_source, &
                                  hdur,hdur_gaussian,t_cmt,dt,stf,t0,sourcearrays, &
                                  ispec_is_elastic  )

  implicit none

  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! source
  integer :: NSOURCES,myrank,it
  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
  double precision, dimension(3,3,NSOURCES) :: nu_source
  double precision, dimension(NSOURCES) :: hdur,hdur_gaussian,t_cmt 
  double precision :: dt
  real(kind=CUSTOM_REAL), dimension(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrays 

  double precision, external :: comp_source_time_function 

  logical, dimension(NSPEC_AB) :: ispec_is_elastic
  
! local parameters
  double precision :: t0,f0
  double precision :: stf 
  real(kind=CUSTOM_REAL) stf_used 
  integer :: isource,iglob,i,j,k,ispec
  
  do isource = 1,NSOURCES

    !   add the source (only if this proc carries the source)
    if(myrank == islice_selected_source(isource)) then

      ispec = ispec_selected_source(isource)

      if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
  
        if( ispec_is_elastic(ispec) ) then

          if(USE_FORCE_POINT_SOURCE) then

             ! note: for use_force_point_source xi/eta/gamma are in the range [1,NGLL*]
             iglob = ibool(nint(xi_source(isource)), &
                           nint(eta_source(isource)), &
                           nint(gamma_source(isource)), &
                           ispec_selected_source(isource))
                                                      
             f0 = hdur(isource) !! using hdur as a FREQUENCY just to avoid changing CMTSOLUTION file format
             t0 = 1.2d0/f0
             
             if (it == 1 .and. myrank == 0) then
                print *,'using a source of dominant frequency ',f0
                print *,'lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
                print *,'lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
             endif
             
             ! we use nu_source(:,3) here because we want a source normal to the surface.
             ! This is the expression of a Ricker; should be changed according maybe to the Par_file.
             !accel(:,iglob) = accel(:,iglob) + &
             !     sngl(nu_source(:,3,isource) * 10000000.d0 * &
             !            (1.d0-2.d0*PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)) * &
             !     exp(-PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)))
             accel(:,iglob) = accel(:,iglob) + &
                  sngl( nu_source(:,3,isource) * 1.d10 * &
                       (1.d0-2.d0*PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)) * &
                       exp(-PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)) )
             
          else   
             
             stf = comp_source_time_function(dble(it-1)*DT-t0-t_cmt(isource),hdur_gaussian(isource))

             !     distinguish between single and double precision for reals
             if(CUSTOM_REAL == SIZE_REAL) then
                stf_used = sngl(stf)
             else
                stf_used = stf
             endif

             !     add source array
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      iglob = ibool(i,j,k,ispec)
                      accel(:,iglob) = accel(:,iglob) + sourcearrays(isource,:,i,j,k)*stf_used
                   enddo
                enddo
             enddo

          endif ! USE_FORCE_POINT_SOURCE
        endif ! ispec_is_elastic
      endif ! ispec_is_inner     
    endif ! myrank
  
  enddo ! NSOURCES

end subroutine elastic_sources

!
!-------------------------------------------------------------------------------------------------
!

subroutine elastic_ocean_load(NSPEC_AB,NGLOB_AB, &
                        ibool,rmass,rmass_ocean_load,accel, &
                        free_surface_normal,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces)

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

! local parameters
  real(kind=CUSTOM_REAL) :: nx,ny,nz
  real(kind=CUSTOM_REAL) :: additional_term,force_normal_comp
  integer :: i,j,k,ispec,iglob
  integer :: igll,iface
  logical,dimension(NGLOB_AB) :: updated_dof_ocean_load
  
!   initialize the updates
  updated_dof_ocean_load(:) = .false.

! for surface elements exactly at the top of the model (ocean bottom)
!  do ispec2D = 1,NSPEC2D_TOP
  
  do iface = 1,num_free_surface_faces

!! DK DK array not created yet for CUBIT      ispec = ibelm_top(ispec2D)

! only for DOFs exactly at the top of the model (ocean bottom)
!    k = NGLLZ
!    do j = 1,NGLLY
!      do i = 1,NGLLX
    
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
        !! DK DK array not created yet for CUBIT            nx = normal_top(1,i,j,ispec2D)
        !! DK DK array not created yet for CUBIT            ny = normal_top(2,i,j,ispec2D)
        !! DK DK array not created yet for CUBIT            nz = normal_top(3,i,j,ispec2D)
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

        !if (SIMULATION_TYPE == 3) then
        !! DK DK array not created yet for CUBIT
        !             b_force_normal_comp = (b_accel(1,iglob)*nx + &
        !                   b_accel(2,iglob)*ny + b_accel(3,iglob)*nz) / rmass(iglob)
        !  b_additional_term = (rmass_ocean_load(iglob) - rmass(iglob)) * b_force_normal_comp
        !! DK DK array not created yet for CUBIT
        !             b_accel(1,iglob) = b_accel(1,iglob) + b_additional_term * nx
        !             b_accel(2,iglob) = b_accel(2,iglob) + b_additional_term * ny
        !             b_accel(3,iglob) = b_accel(3,iglob) + b_additional_term * nz
        !endif

        ! done with this point
        updated_dof_ocean_load(iglob) = .true.

      endif

!      enddo ! NGLLX
!    enddo ! NGLLY
!  enddo ! NSPEC2D_TOP

    enddo ! igll
  enddo ! iface  

end subroutine elastic_ocean_load

