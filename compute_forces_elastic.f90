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
  use specfem_par_elastic
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
      call compute_forces_with_Deville(phase_is_inner, NSPEC_AB,NGLOB_AB, &
                    displ,accel, &
                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                    hprime_xx,hprime_xxT, &
                    hprimewgll_xx,hprimewgll_xxT, &
                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                    kappastore,mustore,jacobian,ibool, &
                    ispec_is_inner_ext_mesh, &
                    ATTENUATION,USE_OLSEN_ATTENUATION, &
                    one_minus_sum_beta,factor_common,alphaval,betaval,gammaval, &
                    NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
                    epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                    epsilondev_xz,epsilondev_yz,iflag_attenuation_store, &
                    rho_vs, &
                    ANISOTROPY,NSPEC_ANISO, &
                    c11store,c12store,c13store,c14store,c15store,c16store,&
                    c22store,c23store,c24store,c25store,c26store,c33store,&
                    c34store,c35store,c36store,c44store,c45store,c46store,&
                    c55store,c56store,c66store )


      !call compute_forces_with_Deville( phase_is_inner ,NSPEC_AB,NGLOB_AB,&
      !              ATTENUATION,USE_OLSEN_ATTENUATION,displ,accel,&
      !              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
      !              hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
      !              kappastore,mustore,jacobian,ibool,ispec_is_inner_ext_mesh, &
      !              NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
      !              xi_source,eta_source,gamma_source,nu_source, &
      !              hdur,hdur_gaussian,t_cmt,dt,stf,t0,sourcearrays, & 
      !              one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
      !              epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,iflag_attenuation_store,&
      !              ABSORBING_CONDITIONS, &
      !              absorbing_boundary_normal,absorbing_boundary_jacobian2D, &
      !              absorbing_boundary_ijk,absorbing_boundary_ispec, &
      !              num_absorbing_boundary_faces, &                      
      !              veloc,rho_vp,rho_vs)                                 
    else
      call compute_forces_no_Deville( phase_is_inner, NSPEC_AB,NGLOB_AB, &
                    displ,accel, &
                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                    hprime_xx,hprime_yy,hprime_zz, &
                    hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                    kappastore,mustore,jacobian,ibool, &
                    ispec_is_inner_ext_mesh, &
                    ATTENUATION,USE_OLSEN_ATTENUATION,&
                    one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,&
                    NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
                    epsilondev_xx,epsilondev_yy,epsilondev_xy,&
                    epsilondev_xz,epsilondev_yz,iflag_attenuation_store,&
                    rho_vs, &
                    ANISOTROPY,NSPEC_ANISO, &
                    c11store,c12store,c13store,c14store,c15store,c16store,&
                    c22store,c23store,c24store,c25store,c26store,c33store,&
                    c34store,c35store,c36store,c44store,c45store,c46store,&
                    c55store,c56store,c66store)
    endif

! adds elastic absorbing boundary term to acceleration (Stacey conditions)
    if(ABSORBING_CONDITIONS) then 
      call compute_forces_elastic_absorbing_boundaries(NSPEC_AB,NGLOB_AB,accel, &
                    ibool,ispec_is_inner_ext_mesh,phase_is_inner, &
                    absorbing_boundary_normal,absorbing_boundary_jacobian2D, &
                    absorbing_boundary_ijk,absorbing_boundary_ispec, &
                    num_absorbing_boundary_faces, &
                    veloc,rho_vp,rho_vs)
    endif

! adds source term (single-force/moment-tensor solution)
    call compute_forces_elastic_source_term( NSPEC_AB,NGLOB_AB,accel, &
                    ibool,ispec_is_inner_ext_mesh,phase_is_inner, &
                    NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                    xi_source,eta_source,gamma_source,nu_source, &
                    hdur,hdur_gaussian,t_cmt,dt,stf,t0,sourcearrays )


! assemble all the contributions between slices using MPI
    if( phase_is_inner .eqv. .false. ) then
      call assemble_MPI_vector_ext_mesh_s(NPROC,NGLOB_AB,accel, &
                    buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                    num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                    my_neighbours_ext_mesh, &
                    request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
    else
      call assemble_MPI_vector_ext_mesh_w(NPROC,NGLOB_AB,accel, &
                    buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh,&
                    max_nibool_interfaces_ext_mesh, &
                    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                    request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
    endif
  
  enddo
  
! update acceleration 
! points inside processor's partition only
!  if(USE_DEVILLE_PRODUCTS) then
!    call compute_forces_with_Deville( .true., NSPEC_AB,NGLOB_AB,&
!                    ATTENUATION,USE_OLSEN_ATTENUATION,displ,accel,&
!                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
!                    hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
!                    kappastore,mustore,jacobian,ibool,ispec_is_inner_ext_mesh, &
!                    NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
!                    xi_source,eta_source,gamma_source,nu_source, &
!                    hdur,hdur_gaussian,t_cmt,dt,stf,t0,sourcearrays, & 
!                    one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
!                    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,iflag_attenuation_store,&
!                    ABSORBING_CONDITIONS, &
!                    absorbing_boundary_normal,absorbing_boundary_jacobian2D, &
!                    absorbing_boundary_ijk,absorbing_boundary_ispec, &
!                    num_absorbing_boundary_faces, &
!                    veloc,rho_vp,rho_vs)
!  else
!    call compute_forces_no_Deville(NSPEC_AB,NGLOB_AB,displ,accel,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
!       hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
!       kappastore,mustore,jacobian,ibool,ispec_is_inner_ext_mesh,.true., &
!       NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,xi_source,eta_source,gamma_source,nu_source,hdur,dt)
!  endif
!
!! assemble all the contributions between slices using MPI
!    call assemble_MPI_vector_ext_mesh_w(NPROC,NGLOB_AB,accel, &
!                    buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh,&
!                    max_nibool_interfaces_ext_mesh, &
!                    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
!                    request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
  
!! DK DK May 2009: removed this because now each slice of a CUBIT + SCOTCH mesh
!! DK DK May 2009: has a different number of spectral elements and therefore
!! DK DK May 2009: only the general non-blocking MPI routines assemble_MPI_vector_ext_mesh_s
!! DK DK May 2009: and assemble_MPI_vector_ext_mesh_w above can be used.
!! DK DK May 2009: For adjoint runs below (SIMULATION_TYPE == 3) they should be used as well.
! if (SIMULATION_TYPE == 3) call assemble_MPI_vector(b_accel,iproc_xi,iproc_eta,addressing, &
!         iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta, &
!         buffer_send_faces_vector,buffer_received_faces_vector,npoin2D_xi,npoin2D_eta, &
!         NPROC_XI,NPROC_ETA,NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NPOIN2DMAX_XY)

end subroutine compute_forces_elastic


!
!-------------------------------------------------------------------------------------------------
!

! absorbing boundary term for elastic media (Stacey conditions)

subroutine compute_forces_elastic_absorbing_boundaries(NSPEC_AB,NGLOB_AB,accel, &
                            ibool,ispec_is_inner,phase_is_inner, &
                            absorbing_boundary_normal,absorbing_boundary_jacobian2D, &
                            absorbing_boundary_ijk,absorbing_boundary_ispec, &
                            num_absorbing_boundary_faces, &
                            veloc,rho_vp,rho_vs)

  implicit none

  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB

! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! array with derivatives of Lagrange polynomials and precalculated products
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
!  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner
  
! Stacey conditions
!  integer  :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,nspec2D_top
!  integer  :: NSPEC2DMAX_XMIN_XMAX_ext,NSPEC2DMAX_YMIN_YMAX_ext
!  integer, dimension(nspec2D_xmin) :: ibelm_xmin
!  integer, dimension(nspec2D_xmax) :: ibelm_xmax
!  integer, dimension(nspec2D_ymin) :: ibelm_ymin
!  integer, dimension(nspec2D_ymax) :: ibelm_ymax
!  integer, dimension(nspec2D_bottom) :: ibelm_bottom
!  integer, dimension(nspec2D_top) :: ibelm_top

  ! local indices i,j,k of all GLL points on xmin boundary in the element
!  integer :: ibelm_gll_xmin(3,NGLLY,NGLLZ,nspec2D_xmin),ibelm_gll_xmax(3,NGLLY,NGLLZ,nspec2D_xmax), &
!            ibelm_gll_ymin(3,NGLLX,NGLLZ,nspec2D_ymin),ibelm_gll_ymax(3,NGLLX,NGLLZ,nspec2D_ymax), &
!            ibelm_gll_bottom(3,NGLLY,NGLLY,nspec2D_bottom),ibelm_gll_top(3,NGLLY,NGLLY,nspec2D_top)  
  
!  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX_ext) :: nimin,nimax,nkmin_eta
!  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX_ext) :: njmin,njmax,nkmin_xi

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vp,rho_vs

!  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,nspec2D_xmin) :: jacobian2D_xmin
!  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,nspec2D_xmax) :: jacobian2D_xmax
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec2D_ymin) :: jacobian2D_ymin
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec2D_ymax) :: jacobian2D_ymax
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM) :: jacobian2D_bottom
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_top) :: jacobian2D_top
!  
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nspec2D_xmin) :: normal_xmin
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nspec2D_xmax) :: normal_xmax
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nspec2D_ymin) :: normal_ymin
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nspec2D_ymax) :: normal_ymax
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM) :: normal_bottom
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_top) :: normal_top

  integer :: num_absorbing_boundary_faces
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_absorbing_boundary_faces) :: absorbing_boundary_normal
  real(kind=CUSTOM_REAL), dimension(NGLLSQUARE,num_absorbing_boundary_faces) :: absorbing_boundary_jacobian2D
  integer, dimension(3,NGLLSQUARE,num_absorbing_boundary_faces) :: absorbing_boundary_ijk
  integer, dimension(num_absorbing_boundary_faces) :: absorbing_boundary_ispec


! local parameters
  real(kind=CUSTOM_REAL) vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw !weight,jacobianl
  integer :: ispec,iglob,i,j,k,iface,igll
  !integer :: num_gll !,igll_i,igll_j,ispec2D
  

! absorbs absorbing-boundary surface using Stacey condition (Clayton & Enquist)
  do iface=1,num_absorbing_boundary_faces

    ispec = absorbing_boundary_ispec(iface)

    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

      ! reference gll points on boundary face 
      do igll = 1,NGLLSQUARE

        ! gets local indices for GLL point
        i = absorbing_boundary_ijk(1,igll,iface)
        j = absorbing_boundary_ijk(2,igll,iface)
        k = absorbing_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob=ibool(i,j,k,ispec)
        vx=veloc(1,iglob)
        vy=veloc(2,iglob)
        vz=veloc(3,iglob)

        ! gets associated normal
        nx = absorbing_boundary_normal(1,igll,iface)
        ny = absorbing_boundary_normal(2,igll,iface)
        nz = absorbing_boundary_normal(3,igll,iface)             

        ! velocity component in normal direction (normal points out of element)
        vn = vx*nx + vy*ny + vz*nz
           
        ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it 
        tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
        ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
        tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

        ! gets associated, weighted jacobian 
        jacobianw = absorbing_boundary_jacobian2D(igll,iface)
        
        ! adds stacey term (weak form)
        accel(1,iglob) = accel(1,iglob) - tx*jacobianw
        accel(2,iglob) = accel(2,iglob) - ty*jacobianw
        accel(3,iglob) = accel(3,iglob) - tz*jacobianw

       enddo
       
    endif    
  enddo
!
!! old way: assumes box model with absorbing-boundary faces oriented with x,y,z planes
!!   xmin  
!  do ispec2D=1,nspec2D_xmin
!
!    ispec=ibelm_xmin(ispec2D)
!
!    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
!
!! old regular mesh
!!       ! exclude elements that are not on absorbing edges
!!       if(nkmin_xi(1,ispec2D) == 0 .or. njmin(1,ispec2D) == 0) cycle
!!
!!       i=1
!!        do k=nkmin_xi(1,ispec2D),NGLLZ
!!           do j=njmin(1,ispec2D),njmax(1,ispec2D)
!
!! new way, unregular element orientation
!      ! reference gll points on boundary face 
!      do igll_j = 1,NGLLZ
!        do igll_i = 1,NGLLY
!          ! gets local indices for GLL point
!          i = ibelm_gll_xmin(1,igll_i,igll_j,ispec2D)
!          j = ibelm_gll_xmin(2,igll_i,igll_j,ispec2D)
!          k = ibelm_gll_xmin(3,igll_i,igll_j,ispec2D)
!
!          ! gets velocity
!          iglob=ibool(i,j,k,ispec)
!          vx=veloc(1,iglob)
!          vy=veloc(2,iglob)
!          vz=veloc(3,iglob)
!
!          ! gets associated normal
!          nx = normal_xmin(1,igll_i,igll_j,ispec2D)
!          ny = normal_xmin(2,igll_i,igll_j,ispec2D)
!          nz = normal_xmin(3,igll_i,igll_j,ispec2D)             
!          !   nx =  normal_xmin(1,j,k,ispec2D)
!          !   ny =  normal_xmin(2,j,k,ispec2D)
!          !   nz =  normal_xmin(3,j,k,ispec2D)
!
!          ! velocity component in normal direction (normal points out of element)
!          vn = vx*nx + vy*ny + vz*nz
!             
!          ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it 
!          tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
!          ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
!          tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)
!
!          ! gets associated jacobian and 2D weights
!          jacobianl = jacobian2D_xmin(igll_i,igll_j,ispec2D)
!          weight = jacobianl*wgllwgll_yz(igll_i,igll_j)
!
!          ! adds stacey term (weak form)
!          accel(1,iglob) = accel(1,iglob) - tx*weight
!          accel(2,iglob) = accel(2,iglob) - ty*weight
!          accel(3,iglob) = accel(3,iglob) - tz*weight
!
!          enddo
!       enddo
!    end if    
!  enddo
!
!!   xmax
!  do ispec2D=1,nspec2D_xmax
!    
!    ispec=ibelm_xmax(ispec2D)
!    
!    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
!       
!      ! reference gll points on boundary face 
!      do igll_j = 1,NGLLZ
!        do igll_i = 1,NGLLY
!          ! gets local indices for GLL point
!          i = ibelm_gll_xmax(1,igll_i,igll_j,ispec2D)
!          j = ibelm_gll_xmax(2,igll_i,igll_j,ispec2D)
!          k = ibelm_gll_xmax(3,igll_i,igll_j,ispec2D)
!
!          ! gets velocity
!          iglob=ibool(i,j,k,ispec)
!          vx=veloc(1,iglob)
!          vy=veloc(2,iglob)
!          vz=veloc(3,iglob)
!
!          ! gets associated normal
!          nx = normal_xmax(1,igll_i,igll_j,ispec2D)
!          ny = normal_xmax(2,igll_i,igll_j,ispec2D)
!          nz = normal_xmax(3,igll_i,igll_j,ispec2D)             
!
!          ! velocity component in normal direction (normal points out of element)
!          vn = vx*nx + vy*ny + vz*nz
!             
!          ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it 
!          tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
!          ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
!          tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)
!
!          ! gets associated jacobian and 2D weights
!          jacobianl = jacobian2D_xmax(igll_i,igll_j,ispec2D)
!          weight = jacobianl*wgllwgll_yz(igll_i,igll_j)
!
!          ! adds stacey term (weak form)
!          accel(1,iglob) = accel(1,iglob) - tx*weight
!          accel(2,iglob) = accel(2,iglob) - ty*weight
!          accel(3,iglob) = accel(3,iglob) - tz*weight             
!
!        enddo
!      enddo
!    end if
!  enddo
!
!!   ymin
!  do ispec2D=1,nspec2D_ymin
!    
!    ispec=ibelm_ymin(ispec2D)
!    
!    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
!       
!      ! reference gll points on boundary face 
!      do igll_j = 1,NGLLZ
!        do igll_i = 1,NGLLX
!          ! gets local indices for GLL point
!          i = ibelm_gll_ymin(1,igll_i,igll_j,ispec2D)
!          j = ibelm_gll_ymin(2,igll_i,igll_j,ispec2D)
!          k = ibelm_gll_ymin(3,igll_i,igll_j,ispec2D)
!
!          ! gets velocity
!          iglob=ibool(i,j,k,ispec)
!          vx=veloc(1,iglob)
!          vy=veloc(2,iglob)
!          vz=veloc(3,iglob)
!
!          ! gets associated normal
!          nx = normal_ymin(1,igll_i,igll_j,ispec2D)
!          ny = normal_ymin(2,igll_i,igll_j,ispec2D)
!          nz = normal_ymin(3,igll_i,igll_j,ispec2D)             
!
!          ! velocity component in normal direction (normal points out of element)
!          vn = vx*nx + vy*ny + vz*nz
!             
!          ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it 
!          tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
!          ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
!          tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)
!
!          ! gets associated jacobian and 2D weights
!          jacobianl = jacobian2D_ymin(igll_i,igll_j,ispec2D)
!          weight = jacobianl*wgllwgll_xz(igll_i,igll_j)
!
!          ! adds stacey term (weak form)
!          accel(1,iglob) = accel(1,iglob) - tx*weight
!          accel(2,iglob) = accel(2,iglob) - ty*weight
!          accel(3,iglob) = accel(3,iglob) - tz*weight             
!
!        enddo
!      enddo
!       
!    endif
!  enddo
!
!!   ymax
!  do ispec2D=1,nspec2D_ymax
!    
!    ispec=ibelm_ymax(ispec2D)
!
!    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
!
!      ! reference gll points on boundary face 
!      do igll_j = 1,NGLLZ
!        do igll_i = 1,NGLLX
!          ! gets local indices for GLL point
!          i = ibelm_gll_ymax(1,igll_i,igll_j,ispec2D)
!          j = ibelm_gll_ymax(2,igll_i,igll_j,ispec2D)
!          k = ibelm_gll_ymax(3,igll_i,igll_j,ispec2D)
!
!          ! gets velocity
!          iglob=ibool(i,j,k,ispec)
!          vx=veloc(1,iglob)
!          vy=veloc(2,iglob)
!          vz=veloc(3,iglob)
!
!          ! gets associated normal
!          nx = normal_ymax(1,igll_i,igll_j,ispec2D)
!          ny = normal_ymax(2,igll_i,igll_j,ispec2D)
!          nz = normal_ymax(3,igll_i,igll_j,ispec2D)             
!
!          ! velocity component in normal direction (normal points out of element)
!          vn = vx*nx + vy*ny + vz*nz
!             
!          ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it 
!          tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
!          ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
!          tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)
!
!          ! gets associated jacobian and 2D weights
!          jacobianl = jacobian2D_ymax(igll_i,igll_j,ispec2D)
!          weight = jacobianl*wgllwgll_xz(igll_i,igll_j)
!
!          ! adds stacey term (weak form)
!          accel(1,iglob) = accel(1,iglob) - tx*weight
!          accel(2,iglob) = accel(2,iglob) - ty*weight
!          accel(3,iglob) = accel(3,iglob) - tz*weight             
!        enddo
!      enddo
!
!    endif
!  enddo
!
!!   bottom (zmin)
!  do ispec2D=1,NSPEC2D_BOTTOM
!    
!    ispec=ibelm_bottom(ispec2D)
!    
!    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
!
!      ! reference gll points on boundary face 
!      do igll_j = 1,NGLLY
!        do igll_i = 1,NGLLX
!          ! gets local indices for GLL point
!          i = ibelm_gll_bottom(1,igll_i,igll_j,ispec2D)
!          j = ibelm_gll_bottom(2,igll_i,igll_j,ispec2D)
!          k = ibelm_gll_bottom(3,igll_i,igll_j,ispec2D)
!
!          ! gets velocity
!          iglob=ibool(i,j,k,ispec)
!          vx=veloc(1,iglob)
!          vy=veloc(2,iglob)
!          vz=veloc(3,iglob)
!
!          ! gets associated normal
!          nx = normal_bottom(1,igll_i,igll_j,ispec2D)
!          ny = normal_bottom(2,igll_i,igll_j,ispec2D)
!          nz = normal_bottom(3,igll_i,igll_j,ispec2D)             
!
!          ! velocity component in normal direction (normal points out of element)
!          vn = vx*nx + vy*ny + vz*nz
!             
!          ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it 
!          tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
!          ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
!          tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)
!
!          ! gets associated jacobian and 2D weights
!          jacobianl = jacobian2D_bottom(igll_i,igll_j,ispec2D)
!          weight = jacobianl*wgllwgll_xy(igll_i,igll_j)
!
!          ! adds stacey term (weak form)
!          accel(1,iglob) = accel(1,iglob) - tx*weight
!          accel(2,iglob) = accel(2,iglob) - ty*weight
!          accel(3,iglob) = accel(3,iglob) - tz*weight             
!
!        enddo
!      enddo
!      
!    endif
!  enddo
!
!! absorbing at top surface - no free-surface?
!  if( ABSORB_TOP_SURFACE ) then
!    do ispec2D=1,NSPEC2D_TOP
!      
!      ispec=ibelm_top(ispec2D)
!      
!      if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
!
!        ! reference gll points on boundary face 
!        do igll_j = 1,NGLLY
!          do igll_i = 1,NGLLX
!            ! gets local indices for GLL point
!            i = ibelm_gll_top(1,igll_i,igll_j,ispec2D)
!            j = ibelm_gll_top(2,igll_i,igll_j,ispec2D)
!            k = ibelm_gll_top(3,igll_i,igll_j,ispec2D)
!
!            ! gets velocity
!            iglob=ibool(i,j,k,ispec)
!            vx=veloc(1,iglob)
!            vy=veloc(2,iglob)
!            vz=veloc(3,iglob)
!
!            ! gets associated normal
!            nx = normal_top(1,igll_i,igll_j,ispec2D)
!            ny = normal_top(2,igll_i,igll_j,ispec2D)
!            nz = normal_top(3,igll_i,igll_j,ispec2D)             
!
!            ! velocity component in normal direction (normal points out of element)
!            vn = vx*nx + vy*ny + vz*nz
!               
!            ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it 
!            tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
!            ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
!            tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)
!
!            ! gets associated jacobian and 2D weights
!            jacobianl = jacobian2D_top(igll_i,igll_j,ispec2D)
!            weight = jacobianl*wgllwgll_xy(igll_i,igll_j)
!
!            ! adds stacey term (weak form)
!            accel(1,iglob) = accel(1,iglob) - tx*weight
!            accel(2,iglob) = accel(2,iglob) - ty*weight
!            accel(3,iglob) = accel(3,iglob) - tz*weight             
!
!          enddo
!        enddo
!
!      endif
!    enddo
!  endif
  
end subroutine compute_forces_elastic_absorbing_boundaries

!
!-------------------------------------------------------------------------------------------------
!

subroutine compute_forces_elastic_source_term( NSPEC_AB,NGLOB_AB,accel, &
                                  ibool,ispec_is_inner,phase_is_inner, &
                                  NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                                  xi_source,eta_source,gamma_source,nu_source, &
                                  hdur,hdur_gaussian,t_cmt,dt,stf,t0,sourcearrays )

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
  
! local parameters
  double precision :: t0,f0
  double precision :: stf 
  real(kind=CUSTOM_REAL) stf_used 
  integer :: isource,iglob,i,j,k
  
  do isource = 1,NSOURCES

    !   add the source (only if this proc carries the source)
    if(myrank == islice_selected_source(isource)) then

      if (ispec_is_inner(ispec_selected_source(isource)) .eqv. phase_is_inner) then

        if(USE_FORCE_POINT_SOURCE) then
           
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
           !     sngl(nu_source(:,3,isource) * 10000000.d0 * (1.d0-2.d0*PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)) * &
           !     exp(-PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)))
           accel(:,iglob) = accel(:,iglob) + &
                sngl(nu_source(:,3,isource) * 1.d10 * (1.d0-2.d0*PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)) * &
                exp(-PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)))
           
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
                    iglob = ibool(i,j,k,ispec_selected_source(isource))
                    accel(:,iglob) = accel(:,iglob) + sourcearrays(isource,:,i,j,k)*stf_used
                 enddo
              enddo
           enddo

        endif ! USE_FORCE_POINT_SOURCE
      endif ! ispec_is_inner     
    endif ! myrank
  
  enddo ! NSOURCES

end subroutine compute_forces_elastic_source_term

