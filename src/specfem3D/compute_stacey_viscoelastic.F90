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

!=============================================================================
!
! local absorbing boundary condition
!
!=============================================================================
!
! Stacey absorbing boundary
!
! reference: Stacey, 1988, IMPROVED TRANSPARENT BOUNDARY FORMULATIONS FOR THE ELASTIC-WAVE EQUATION,
!            BSSA, 78 (6), p. 2089-2097.
!
!
! Hagstrom-Warburton absorbing boundary
!
! reference: Hagstrom, T. and T. Warburton, 2004,
!            A new auxiliary variable formulation of high-order local radiation boundary conditions:
!            corner compatibility conditions and extensions to first order systems,
!            Wave Motion, 39 (4), p. 327–338.
!
!            Rabinovich et al. 2011,
!            A finite element scheme with a high order absorbing boundary condition for elastodynamics,
!            Comput. Methods Appl. Mech. Eng., 200 (23-24), p. 2048-2066.

module stacey_par

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLSQUARE,NDIM,NDIM2D,MIDX,MIDY,MIDZ, &
                       USE_STACEY_P3
  use shared_parameters, only: NGNOD,NGNOD2D

  implicit none

  private

  ! Hagstrom-Warburton absorbing boundary condition
  !
  ! note: The H-W absorbing condition is a high-order local boundary condition.
  !       At the moment, the implementation of the H-W boundary is experimental and leads to instabilities.
  !       It is turned off by default for now, but left here for further experiments and improvements.

  !       The order of the H-W boundary can be set by parameter HW_ORDER_P.
  !       There are different schemes (Pade, Chebyshev, Newman) presented in the papers how to determine
  !       the coefficients of the boundary terms.
  !       The implementation here uses a Newmark scheme to solve the boundary's extra one-way wave equations.
  !       This local Newmark scheme can further subdivide the time step DT of the (main) simulation
  !       to try increasing the stability.
  !
  !       Instabilities seem to arise first from edge and corner nodes. To apply further compatibility conditions
  !       to these nodes, the flag HW_APPLY_EDGE_CORNER_COMPATIBILITY can be used.
  !
  !       The H-W boundary requires a reference frequency at which it centers the absorption corrections.
  !       This reference frequency can be set by the parameter `f0_FOR_PML` in the Par_file.
  !
  ! for Hagstrom-Warburton (HW) absorbing boundaries
  logical, parameter :: USE_HW_ABC = .false.

  ! parameter for H-W ABC order P
  integer, parameter :: HW_ORDER_P = 2
  ! coefficient schemes
  integer, parameter :: HW_COEFF_PADE = 1
  integer, parameter :: HW_COEFF_CHEBYSHEV = 2
  integer, parameter :: HW_COEFF_NEWMAN = 3
  integer, parameter :: HW_COEFF_SCHEME = HW_COEFF_PADE
  ! number of substeps for auxiliary state updates
  integer, parameter :: HW_ABC_SUBSTEPS = 1
  ! apply edge/corner compatibility conditions
  logical, parameter :: HW_APPLY_EDGE_CORNER_COMPATIBILITY = .true.

  ! free parameter coefficients
  real(kind=CUSTOM_REAL), allocatable, dimension(:) :: hw_a_coeff, hw_sigma
  ! auxiliary variables at each GLL point of each boundary face
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: hw_phi, hw_phi_dot, hw_phi_dotdot
  ! corner-node weighting: 1 / (number of boundary faces sharing each GLL node)
  ! prevents double-counting the phi1 traction at edge/corner nodes (Kucherov & Givoli 2010)
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:) :: hw_corner_weight
  ! for orthogonal basis
  ! 2D shape functions derivatives
  double precision, dimension(:,:,:,:), allocatable :: dershape2D_x,dershape2D_y,dershape2D_z
  ! reference angular frequency
  real(kind=CUSTOM_REAL) :: hw_omega_c

  ! Edge and corner topology for H-W compatibility conditions (Fix C, Fix D)
  type :: edge_info_type
    integer :: id_perp      ! 1 or 2 (which tangential dir is perpendicular to edge)
    integer :: i_edge       ! 1 or NGLL (which GLL row is the edge)
    integer :: sigma        ! +1 or -1 (sign of adjacent face normal)
    integer :: adj_face     ! index of adjacent absorbing face (-1 if none)
    !real(kind=CUSTOM_REAL) :: jacobian_perp  ! Jacobian in perpendicular direction
    ! Component wave speeds: which component uses Vp vs Vs
    ! 1=normal, 2=t1, 3=t2; 1 if Vp, 0 if Vs
    integer :: comp_vp(3)
  end type edge_info_type

  type :: corner_info_type
    integer :: i_edge_t1    ! 1 or NGLL
    integer :: i_edge_t2    ! 1 or NGLL
    integer :: sigma_t1     ! +1 or -1
    integer :: sigma_t2     ! +1 or -1
    integer :: adj_face_t1  ! adjacent face in t1 direction
    integer :: adj_face_t2  ! adjacent face in t2 direction
    ! Component wave speeds for Strang splitting steps
    ! step 1 (t2 direction): 1=normal, 2=t1, 3=t2; 1 if Vp, 0 if Vs
    integer :: comp_vp_step1(3)
    ! step 2 (t1 direction): 1=normal, 2=t1, 3=t2; 1 if Vp, 0 if Vs
    integer :: comp_vp_step2(3)
  end type corner_info_type

  type(edge_info_type), allocatable :: face_edges(:,:)  ! (4, num_faces)
  type(corner_info_type), allocatable :: face_corners(:,:)  ! (4, num_faces)

  ! public
  public :: compute_stacey_viscoelastic

  ! H-W boundary
  public :: USE_HW_ABC, allocate_hw_abc, deallocate_hw_abc, update_hw_abc_states

contains


  subroutine compute_stacey_viscoelastic(NSPEC_AB,NGLOB_AB,accel, &
                                         ibool, &
                                         abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                         abs_boundary_ijk,abs_boundary_ispec, &
                                         num_abs_boundary_faces, &
                                         displ,veloc,rho_vp,rho_vs, &
                                         ispec_is_elastic, &
                                         b_num_abs_boundary_faces,b_absorb_field)

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: displ,veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! adjoint simulations
  integer,intent(in) :: b_num_abs_boundary_faces
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  ! chooses approximation
  if (USE_HW_ABC) then
    ! uses Hagstrom-Warburton absorbing boundary
    call compute_hw_elastic(NSPEC_AB,NGLOB_AB,accel,ibool, &
                            abs_boundary_jacobian2Dw, &
                            abs_boundary_ijk,abs_boundary_ispec, &
                            num_abs_boundary_faces, &
                            displ,veloc,rho_vp,rho_vs, &
                            ispec_is_elastic, &
                            b_num_abs_boundary_faces,b_absorb_field)

  else if (USE_STACEY_P3) then
    ! uses Stacey P3 approximation (second-order)
    call compute_stacey_elastic_p3(NSPEC_AB,NGLOB_AB,accel, &
                                  ibool, &
                                  abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                  abs_boundary_ijk,abs_boundary_ispec, &
                                  num_abs_boundary_faces, &
                                  displ,veloc,rho_vp,rho_vs, &
                                  ispec_is_elastic, &
                                  b_num_abs_boundary_faces,b_absorb_field)

  else
    ! uses Stacey P1 approximation (default)
    call compute_stacey_elastic_p1(NSPEC_AB,NGLOB_AB,accel, &
                                   ibool, &
                                   abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                   abs_boundary_ijk,abs_boundary_ispec, &
                                   num_abs_boundary_faces, &
                                   veloc,rho_vp,rho_vs, &
                                   ispec_is_elastic, &
                                   b_num_abs_boundary_faces,b_absorb_field)
  endif

  end subroutine compute_stacey_viscoelastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_elastic_p1(NSPEC_AB,NGLOB_AB,accel, &
                                       ibool, &
                                       abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                       abs_boundary_ijk,abs_boundary_ispec, &
                                       num_abs_boundary_faces, &
                                       veloc,rho_vp,rho_vs, &
                                       ispec_is_elastic, &
                                       b_num_abs_boundary_faces,b_absorb_field)

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! adjoint simulations
  integer,intent(in) :: b_num_abs_boundary_faces
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  ! local parameters
  real(kind=CUSTOM_REAL) :: vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll

  ! Stacey P1 approximation (default)
  !
  ! Stacey starts with the equations of motions
  !   (1) U_tt = α^2 U_xx + β^2 U_zz + (α^2 - β^2) W_xz
  !   (2) W_tt = β^2 W_xx + α^2 W_zz + (α^2 - β^2) U_xz
  !
  ! After writing out these 2 equations, the paper mentions "where U and W are the vertical and horizontal displacements,
  ! and α and β are the longitudinal and transverse wave velocities".
  ! Given the above equations, this seems incorrect since U_xx gets multiplied by α^2 and W_xx by β^2.
  ! For example, a P-wave traveling along horizontal x-direction would have a non-zero β^2 W_xx,
  ! i.e. traveling with speed β.
  !
  ! The equations point to the following correct definition:
  !   U is horizontal displacement along x-direction, W is vertical displacement along z-direction.
  !
  ! P1 uses equations
  !   (3) U_z = - 1/β U_t
  !   (4) W_z = - 1/α W_t
  ! where U_z denotes the spatial derivate d/dz of U, U_t the temporal derivate d/dt of U.
  !
  ! The boundary is a horizontal internal boundary along which z is maximal and has z as outward normal.
  ! Thus, W is the normal displacement to the boundary and U the tangential displacement component to the boundary.
  !
  ! stress-strain relation for normal stress τ_zz = ρ(α^2 - 2β^2) U_x + ρα^2 W_z  (vertical normal stress)
  !                                          τ_xx = ρ(α^2 - 2β^2) W_z + ρα^2 U_x  (horizontal normal stress)
  !                         and shear stress τ_xz = ρβ^2 (U_z + W_x)
  !
  ! inserting (4) into normal stress relation, the normal stress becomes
  !   τ_zz = − ρα W_t + .. (ignoring additional term ρ(α^2 - 2β^2) U_x)
  ! inserting (3) into shear stress relation, the shear stress becomes
  !   τ_xz = − ρβ U_t + .. (ignoring additional term ρβ^2 W_x)
  !
  ! the first terms τ_zz = − ρα W_t and τ_xz = − ρβ U_t are the P1 approximations.
  !
  ! Normal and shear traction on the horizontal boundary, i.e.,
  !   T_n = τ . n == T_z == τ_zz (normal) and
  !   T_t = τ . t == T_x == τ_xz (tangential),
  ! use velocities like
  !   T_normal      = - rho Vp v_n                (normal)
  !   T_tangential1 = - rho Vs v_tangential       (tangential)
  !
  ! note: traction == impedance x velocity and impedance = density x velocity
  !
  ! units:  [Pa] =  [kg/m^3] [m/s] [m/s] = [kg / m / s^2 ]

! openmp solver
!$OMP PARALLEL if (num_abs_boundary_faces > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(iface,ispec,igll,i,j,k,iglob,vx,vy,vz,vn,nx,ny,nz,tx,ty,tz,jacobianw)
!$OMP DO
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    ! only for elastic domains
    if (.not. ispec_is_elastic(ispec)) cycle

    ! reference GLL points on boundary face
    do igll = 1,NGLLSQUARE
      ! gets local indices for GLL point
      i = abs_boundary_ijk(1,igll,iface)
      j = abs_boundary_ijk(2,igll,iface)
      k = abs_boundary_ijk(3,igll,iface)

      ! gets velocity
      iglob = ibool(i,j,k,ispec)

      vx = veloc(1,iglob)
      vy = veloc(2,iglob)
      vz = veloc(3,iglob)

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
!$OMP ATOMIC
      accel(1,iglob) = accel(1,iglob) - tx*jacobianw
!$OMP ATOMIC
      accel(2,iglob) = accel(2,iglob) - ty*jacobianw
!$OMP ATOMIC
      accel(3,iglob) = accel(3,iglob) - tz*jacobianw

      ! for kernel simulations
      if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
        b_absorb_field(1,igll,iface) = tx*jacobianw
        b_absorb_field(2,igll,iface) = ty*jacobianw
        b_absorb_field(3,igll,iface) = tz*jacobianw
      endif

    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_stacey_elastic_p1

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_elastic_p3(NSPEC_AB,NGLOB_AB,accel, &
                                       ibool, &
                                       abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                       abs_boundary_ijk,abs_boundary_ispec, &
                                       num_abs_boundary_faces, &
                                       displ,veloc,rho_vp,rho_vs, &
                                       ispec_is_elastic, &
                                       b_num_abs_boundary_faces,b_absorb_field)

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE, &
    xstore,ystore,zstore,hprime_xx,rhostore

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: displ,veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! adjoint simulations
  integer,intent(in) :: b_num_abs_boundary_faces
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  ! local parameters
  real(kind=CUSTOM_REAL) :: vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll

  ! for improved Stacey condition
  ! 2D local arrays for surface geometry and displacements
  real(kind=CUSTOM_REAL) :: x_2D(NGLLX, NGLLY), y_2D(NGLLX, NGLLY), z_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: nx_2D(NGLLX, NGLLY), ny_2D(NGLLX, NGLLY), nz_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: t1x_2D(NGLLX, NGLLY), t1y_2D(NGLLX, NGLLY), t1z_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: t2x_2D(NGLLX, NGLLY), t2y_2D(NGLLX, NGLLY), t2z_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: un_2D(NGLLX, NGLLY), ut1_2D(NGLLX, NGLLY), ut2_2D(NGLLX, NGLLY)  ! projected displacement

  ! Derivative arrays
  real(kind=CUSTOM_REAL) :: dx_ds1(NGLLX, NGLLY), dx_ds2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dy_ds1(NGLLX, NGLLY), dy_ds2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dz_ds1(NGLLX, NGLLY), dz_ds2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dun_ds1(NGLLX, NGLLY), dun_ds2(NGLLX, NGLLY)   ! ds1, ds2 derivative of displacement
  real(kind=CUSTOM_REAL) :: dut1_ds1(NGLLX, NGLLY), dut1_ds2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dut2_ds1(NGLLX, NGLLY), dut2_ds2(NGLLX, NGLLY)

  ! for Stacey condition using P3 approximation
  real(kind=CUSTOM_REAL) :: a11, a12, a22, det_a, inv_a11, inv_a12, inv_a22
  real(kind=CUSTOM_REAL) :: dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2
  real(kind=CUSTOM_REAL) :: dt1_un, dt2_un, dt1_ut1, dt2_ut2

  real(kind=CUSTOM_REAL) :: rhol, csl, cpl, fac
  real(kind=CUSTOM_REAL) :: t1x, t1y, t1z, t2x, t2y, t2z, t1_norm, hp1, hp2
  real(kind=CUSTOM_REAL) :: traction_n, traction_t1, traction_t2

  logical :: mask_vary(3)
  integer :: id1, id2, a, b, l
  integer :: face_iglob(NGLLX, NGLLY)
  !debug
  !real(kind=CUSTOM_REAL) :: nt1, nt2, t1t2, rh

  ! Stacey P3 approximation (second-order)
  !
  ! uses equations
  !   (7) U_z = - 1/β U_t  + (β-α)/β W_x
  !   (8) W_z = - 1/α W_t  + (β-α)/α U_x
  !
  ! The boundary is a horizontal internal boundary along which z is maximal and has z as outward normal.
  ! W is the normal displacement to the boundary and U the tangential displacement component to the boundary.
  !
  ! stress-strain relation for normal stress τ_zz = ρ(α^2 - 2β^2) U_x + ρα^2 W_z
  !                         and shear stress τ_xz = ρβ^2 (U_z + W_x)
  !
  ! inserting (8) into normal stress relation, the normal stress becomes
  !   τ_zz = − ρα W_t − ρβ(2β−α) U_x
  ! inserting (7) into shear stress relation, the shear stress becomes
  !   τ_xz = − ρβ U_t + ρβ(2β−α) W_x
  !
  ! the first terms τ_zz = − ρα W_t and τ_xz = − ρβ U_t are the P1 approximations.
  !
  ! general expression in 3D for normal traction becomes
  !   τ_nn = −ρα ∂t un − ρβ(2β−α) (∂ut1/∂t1 + ∂ut2/∂t2)
  ! and shear tractions
  !   τ_t1 = −ρβ ∂t ut1 + ρβ(2β−α) ∂un/∂t1
  !   τ_t2 = −ρβ ∂t ut1 + ρβ(2β−α) ∂un/∂t2
  ! where un is the normal displacement, ut1 and ut2 displacement in tangential directions t1 and t2.
  !
  ! the additional P3 terms use spatial derivatives of displacement like
  !   T_normal      = - rho Vs (2 Vs - Vp)(\partial_tangential1 u_tangential1 + \partial_tangential2 u_tangential2)
  !   T_tangential1 = + rho Vs (2 Vs - Vp) \partial_tangential1 u_n
  !   T_tangential2 = + rho Vs (2 Vs - Vp) \partial_tangential2 u_n
  !
  ! units [Pa] = [kg/m^3] [m/s] [m/s] [ m/m ] = [kg / m / s^2]

! openmp solver
!$OMP PARALLEL if (num_abs_boundary_faces > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(iface,ispec,igll,i,j,k,iglob,vx,vy,vz,vn,nx,ny,nz,tx,ty,tz,jacobianw, &
!$OMP         x_2D,y_2D,z_2D,un_2D,ut1_2D,ut2_2D,nx_2D,ny_2D,nz_2D, &
!$OMP         t1x_2D,t1y_2D,t1z_2D,t2x_2D,t2y_2D,t2z_2D, &
!$OMP         dx_ds1,dx_ds2,dy_ds1,dy_ds2,dz_ds1,dz_ds2, &
!$OMP         dun_ds1,dun_ds2,dut1_ds1,dut1_ds2,dut2_ds1,dut2_ds2, &
!$OMP         a11,a12,a22,det_a,inv_a11,inv_a12,inv_a22, &
!$OMP         dot_g1_t1,dot_g2_t1,dot_g1_t2,dot_g2_t2, &
!$OMP         dt1_un,dt2_un,dt1_ut1,dt2_ut2, &
!$OMP         rhol,csl,cpl,fac,t1x,t1y,t1z,t2x,t2y,t2z,t1_norm,hp1,hp2, &
!$OMP         traction_n,traction_t1,traction_t2,mask_vary,id1,id2,a,b,l,face_iglob)
!$OMP DO
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    ! only for elastic domains
    if (.not. ispec_is_elastic(ispec)) cycle

    ! prepare spatial derivative arrays
    ! find the varying local coordinates for the face to map 1D index igll to 2D grid (a,b)
    mask_vary(:) = .false.
    do igll = 2, NGLLSQUARE
      if (abs_boundary_ijk(1, igll, iface) /= abs_boundary_ijk(1, 1, iface)) mask_vary(1) = .true.
      if (abs_boundary_ijk(2, igll, iface) /= abs_boundary_ijk(2, 1, iface)) mask_vary(2) = .true.
      if (abs_boundary_ijk(3, igll, iface) /= abs_boundary_ijk(3, 1, iface)) mask_vary(3) = .true.
    enddo

    id1 = 0; id2 = 0
    if (.not. mask_vary(1)) then
      id1 = 2; id2 = 3
    else if (.not. mask_vary(2)) then
      id1 = 1; id2 = 3
    else
      id1 = 1; id2 = 2
    endif

    ! gather coordinates and fields onto the 2D grid face
    do igll = 1, NGLLSQUARE
      i = abs_boundary_ijk(1,igll,iface)
      j = abs_boundary_ijk(2,igll,iface)
      k = abs_boundary_ijk(3,igll,iface)
      iglob = ibool(i,j,k,ispec)

      a = abs_boundary_ijk(id1,igll,iface)
      b = abs_boundary_ijk(id2,igll,iface)

      face_iglob(a,b) = iglob

      x_2D(a,b) = xstore(iglob)
      y_2D(a,b) = ystore(iglob)
      z_2D(a,b) = zstore(iglob)

      nx_2D(a,b) = abs_boundary_normal(1,igll,iface)
      ny_2D(a,b) = abs_boundary_normal(2,igll,iface)
      nz_2D(a,b) = abs_boundary_normal(3,igll,iface)
    enddo

    ! projections at each point
    ! (assumes NGLLX == NGLLY == NGLLZ)
    do b = 1,NGLLY
      do a = 1,NGLLX
        iglob = face_iglob(a,b)

        ! pre-compute the orthonormal tangential basis
        ! note: use the face mid-point to evaluate the orthonormal basis
        !       and assign it to all GLL nodes for this face to avoid problems if the boundary is curved
        !nx = nx_2D(MIDX,MIDY)
        !ny = ny_2D(MIDX,MIDY)
        !nz = nz_2D(MIDX,MIDY)
        ! normal
        nx = nx_2D(a,b)
        ny = ny_2D(a,b)
        nz = nz_2D(a,b)

        ! construct right-handed orthonormal tangential basis t1, t2
        ! choose a non-collinear vector to n to compute t1
        ! (Hughes-Moeller method)
        if (abs(nx) < 0.9_CUSTOM_REAL) then
          t1x = 0.0_CUSTOM_REAL
          t1y = -nz
          t1z = ny
        else
          t1x = -ny
          t1y = nx
          t1z = 0.0_CUSTOM_REAL
        endif
        t1_norm = sqrt(t1x*t1x + t1y*t1y + t1z*t1z)

        ! avoid divison by zero
        if (abs(t1_norm) < 1.d-24) t1_norm = 1._CUSTOM_REAL

        ! normalizes t1
        t1x = t1x / t1_norm; t1y = t1y / t1_norm; t1z = t1z / t1_norm

        ! t2 = n x t1
        t2x = ny * t1z - nz * t1y
        t2y = nz * t1x - nx * t1z
        t2z = nx * t1y - ny * t1x

        ! store basis for reuse in igll loop
        t1x_2D(a,b) = t1x; t1y_2D(a,b) = t1y; t1z_2D(a,b) = t1z
        t2x_2D(a,b) = t2x; t2y_2D(a,b) = t2y; t2z_2D(a,b) = t2z

        ! project displacement
        un_2D(a,b)  = displ(1,iglob) * nx  + displ(2,iglob) * ny  + displ(3,iglob) * nz
        ut1_2D(a,b) = displ(1,iglob) * t1x + displ(2,iglob) * t1y + displ(3,iglob) * t1z
        ut2_2D(a,b) = displ(1,iglob) * t2x + displ(2,iglob) * t2y + displ(3,iglob) * t2z
      enddo
    enddo

    ! compute reference derivatives using 1D GLL derivative matrix
    ! (assumes NGLLX == NGLLY == NGLLZ)
    do b = 1, NGLLY
      do a = 1, NGLLX
        dx_ds1(a,b) = 0.0_CUSTOM_REAL; dx_ds2(a,b) = 0.0_CUSTOM_REAL
        dy_ds1(a,b) = 0.0_CUSTOM_REAL; dy_ds2(a,b) = 0.0_CUSTOM_REAL
        dz_ds1(a,b) = 0.0_CUSTOM_REAL; dz_ds2(a,b) = 0.0_CUSTOM_REAL
        dun_ds1(a,b) = 0.0_CUSTOM_REAL; dun_ds2(a,b) = 0.0_CUSTOM_REAL
        dut1_ds1(a,b) = 0.0_CUSTOM_REAL; dut1_ds2(a,b) = 0.0_CUSTOM_REAL
        dut2_ds1(a,b) = 0.0_CUSTOM_REAL; dut2_ds2(a,b) = 0.0_CUSTOM_REAL

        do l = 1, NGLLX
          hp1 = hprime_xx(a,l)
          dx_ds1(a,b) = dx_ds1(a,b) + x_2D(l,b) * hp1
          dy_ds1(a,b) = dy_ds1(a,b) + y_2D(l,b) * hp1
          dz_ds1(a,b) = dz_ds1(a,b) + z_2D(l,b) * hp1
          dun_ds1(a,b) = dun_ds1(a,b) + un_2D(l,b) * hp1
          dut1_ds1(a,b) = dut1_ds1(a,b) + ut1_2D(l,b) * hp1
          dut2_ds1(a,b) = dut2_ds1(a,b) + ut2_2D(l,b) * hp1

          hp2 = hprime_xx(b,l)
          dx_ds2(a,b) = dx_ds2(a,b) + x_2D(a,l) * hp2
          dy_ds2(a,b) = dy_ds2(a,b) + y_2D(a,l) * hp2
          dz_ds2(a,b) = dz_ds2(a,b) + z_2D(a,l) * hp2
          dun_ds2(a,b) = dun_ds2(a,b) + un_2D(a,l) * hp2
          dut1_ds2(a,b) = dut1_ds2(a,b) + ut1_2D(a,l) * hp2
          dut2_ds2(a,b) = dut2_ds2(a,b) + ut2_2D(a,l) * hp2
        enddo
      enddo
    enddo

    ! reference GLL points on boundary face
    do igll = 1,NGLLSQUARE
      ! gets local indices for GLL point
      i = abs_boundary_ijk(1,igll,iface)
      j = abs_boundary_ijk(2,igll,iface)
      k = abs_boundary_ijk(3,igll,iface)

      ! gets velocity
      iglob = ibool(i,j,k,ispec)

      vx = veloc(1,iglob)
      vy = veloc(2,iglob)
      vz = veloc(3,iglob)

      ! local face indexing
      a = abs_boundary_ijk(id1,igll,iface)
      b = abs_boundary_ijk(id2,igll,iface)

      ! gets associated normal
      nx = nx_2D(a,b)
      ny = ny_2D(a,b)
      nz = nz_2D(a,b)

      ! velocity component in normal direction (normal points out of element)
      vn = vx*nx + vy*ny + vz*nz

      ! P1 contribution
      ! Stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it
      !   τ_zz = − ρα W_t     ! normal traction
      !   τ_xz = − ρβ U_t     ! shear traction
      tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
      ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
      tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

      ! additional P3 second-order contribution
      ! metric tensor components a_{alpha beta} = g_alpha . g_beta
      a11 = dx_ds1(a,b)**2 + dy_ds1(a,b)**2 + dz_ds1(a,b)**2
      a12 = dx_ds1(a,b) * dx_ds2(a,b) + dy_ds1(a,b) * dy_ds2(a,b) + dz_ds1(a,b) * dz_ds2(a,b)
      a22 = dx_ds2(a,b)**2 + dy_ds2(a,b)**2 + dz_ds2(a,b)**2

      ! surface Jacobian (squared)
      det_a = a11 * a22 - a12**2

      ! avoid divison by zero
      if (abs(det_a) < 1.d-24) det_a = 1._CUSTOM_REAL

      ! inverse metric tensor
      inv_a11 = a22 / det_a
      inv_a12 = -a12 / det_a
      inv_a22 = a11 / det_a

      ! look up pre-computed basis
      t1x = t1x_2D(a,b); t1y = t1y_2D(a,b); t1z = t1z_2D(a,b)
      t2x = t2x_2D(a,b); t2y = t2y_2D(a,b); t2z = t2z_2D(a,b)

      ! dot-products
      dot_g1_t1 = dx_ds1(a,b) * t1x + dy_ds1(a,b) * t1y + dz_ds1(a,b) * t1z
      dot_g2_t1 = dx_ds2(a,b) * t1x + dy_ds2(a,b) * t1y + dz_ds2(a,b) * t1z

      dot_g1_t2 = dx_ds1(a,b) * t2x + dy_ds1(a,b) * t2y + dz_ds1(a,b) * t2z
      dot_g2_t2 = dx_ds2(a,b) * t2x + dy_ds2(a,b) * t2y + dz_ds2(a,b) * t2z

      dt1_un  = (inv_a11*dun_ds1(a,b) + inv_a12*dun_ds2(a,b)) * dot_g1_t1 &
              + (inv_a12*dun_ds1(a,b) + inv_a22*dun_ds2(a,b)) * dot_g2_t1
      dt2_un  = (inv_a11*dun_ds1(a,b) + inv_a12*dun_ds2(a,b)) * dot_g1_t2 &
              + (inv_a12*dun_ds1(a,b) + inv_a22*dun_ds2(a,b)) * dot_g2_t2

      dt1_ut1 = (inv_a11*dut1_ds1(a,b) + inv_a12*dut1_ds2(a,b)) * dot_g1_t1 &
              + (inv_a12*dut1_ds1(a,b) + inv_a22*dut1_ds2(a,b)) * dot_g2_t1
      dt2_ut2 = (inv_a11*dut2_ds1(a,b) + inv_a12*dut2_ds2(a,b)) * dot_g1_t2 &
              + (inv_a12*dut2_ds1(a,b) + inv_a22*dut2_ds2(a,b)) * dot_g2_t2

      rhol = rhostore(i,j,k,ispec)
      csl = rho_vs(i,j,k,ispec) / rhol
      cpl = rho_vp(i,j,k,ispec) / rhol

      ! Stacey P3
      ! normal stress
      !   τ_zz = − ρα W_t − ρβ(2β−α) U_x
      ! shear stress
      !   τ_xz = − ρβ U_t + ρβ(2β−α) W_x
      !
      ! the first terms τ_zz = − ρα W_t and τ_xz = − ρβ U_t are the P1 approximations.
      ! here we add the additional second terms.
      !
      ! general expression in 3D for normal traction becomes
      !   τ_nn = −ρα ∂t un − ρβ(2β−α) (∂ut1/∂t1 + ∂ut2/∂t2)
      ! and shear tractions
      !   τ_t1 = −ρβ ∂t ut1 + ρβ(2β−α) ∂un/∂t1
      !   τ_t2 = −ρβ ∂t ut2 + ρβ(2β−α) ∂un/∂t2
      ! where un is the normal displacement, ut1 and ut2 displacement in tangential directions t1 and t2.
      !
      ! note: accel gets subtracted by tx, ty, tz, so we add a + sign to normal and a - sign to tangential contribution.
      !       T_normal     = - rho * vs * (2 * vs - vp) (d/dt1 u_t1 + d/dt2 u_t2)
      !       T_tangential = + rho * vs * (2 * vs - vp) d/dx_tangential u_normal
      fac = rho_vs(i,j,k,ispec) * (2.0_CUSTOM_REAL * csl - cpl)

      traction_n = fac * (dt1_ut1 + dt2_ut2)     ! normal traction
      traction_t1 = - fac * dt1_un               ! shear traction
      traction_t2 = - fac * dt2_un

      ! total contribution
      tx = tx + traction_n * nx + traction_t1 * t1x + traction_t2 * t2x
      ty = ty + traction_n * ny + traction_t1 * t1y + traction_t2 * t2y
      tz = tz + traction_n * nz + traction_t1 * t1z + traction_t2 * t2z

      ! gets associated, weighted jacobian
      jacobianw = abs_boundary_jacobian2Dw(igll,iface)

      ! adds stacey term (weak form)
!$OMP ATOMIC
      accel(1,iglob) = accel(1,iglob) - tx*jacobianw
!$OMP ATOMIC
      accel(2,iglob) = accel(2,iglob) - ty*jacobianw
!$OMP ATOMIC
      accel(3,iglob) = accel(3,iglob) - tz*jacobianw

      ! for kernel simulations
      if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
        b_absorb_field(1,igll,iface) = tx*jacobianw
        b_absorb_field(2,igll,iface) = ty*jacobianw
        b_absorb_field(3,igll,iface) = tz*jacobianw
      endif

    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_stacey_elastic_p3


!-------------------------------------------------------------------------------------------------
!
! Hagstrom-Warburton absorbing boundary condition
!
!-------------------------------------------------------------------------------------------------

  subroutine allocate_hw_abc()

  use constants, only: myrank,IMAIN,PI,TWO_PI,NDIM2D

  use shared_parameters, only: f0_FOR_PML

  use specfem_par, only: num_abs_boundary_faces, abs_boundary_ijk, abs_boundary_ispec, abs_boundary_normal, &
                         ibool, NGLOB_AB, xigll, yigll, zigll

  use specfem_par_elastic, only: ispec_is_elastic

  implicit none
  ! local parameters
  integer :: ier, iface, igll, i, j, k, ispec, iglob, id1, id2, id3, a, b, p
  logical :: mask_vary(3)
  integer, allocatable :: node_count(:)
  real(kind=CUSTOM_REAL) :: dot_n,norm_n
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: node_normal
  ! 2D shape functions
  double precision, dimension(:,:,:), allocatable :: shape2D_x,shape2D_y,shape2D_z
  ! corner/edge detection
  integer, allocatable :: face_id1(:), face_id2(:), face_id3(:), face_ispec(:)
  real(kind=CUSTOM_REAL), allocatable :: face_nx_arr(:), face_ny_arr(:), face_nz_arr(:)
  real(kind=CUSTOM_REAL) :: nx, ny, nz
  integer, dimension(NGLLX) :: edge1_nodes, edge2_nodes, edge3_nodes, edge4_nodes
  integer, dimension(:), allocatable :: node_face_count
  integer, dimension(:,:), allocatable :: node_face_list
  integer, dimension(:), allocatable :: node_face_idx
  integer :: max_faces_per_node

  ! vector dot product threshold
  ! if dot product == 1 both normal vectors align, if dot product == 0 vectors are orthogonal (we use 0.1 as threshold)
  real(kind=CUSTOM_REAL), parameter :: THRESHOLD_DOT_PRODUCT = 0.1_CUSTOM_REAL

  ! assemble resulting normal vectors
  logical, parameter :: USE_SHARED_NORMALS = .false.

  ! checks if anything to do
  if (num_abs_boundary_faces <= 0) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  H-W boundary setup:'
    write(IMAIN,*) '    order P: ',HW_ORDER_P
    call flush_IMAIN()
  endif

  ! free parameter coefficients
  ! note: using a setting of all a_j == 1 is equivalent to a Pade approximation,
  !       other choices could be Chebyshev or Newman coefficients.
  !       See:
  !       Hagstrom et al. 2007, LOCAL HIGH-ORDER ABSORBING BOUNDARY CONDITIONS FOR TIME-DEPENDENT WAVES IN GUIDES,
  !       J. Comp. Acoustic, 15 (1), p. 1-22
  allocate(hw_a_coeff(HW_ORDER_P), &
           hw_sigma(HW_ORDER_P),stat=ier)
  if (ier /= 0) stop 'Error allocating H-W ABC coefficients'

  ! note: Hagstrom et al. (2007), LOCAL HIGH-ORDER ABSORBING BOUNDARY CONDITIONS FOR TIME-DEPENDENT WAVES IN GUIDES
  !       for P-order boundary, free parameters a_j with j = 0,1,..,P, and a_0 is always equal to 1;
  !       auxiliary variables phi_j with j = 1,..,P
  !       they state that "for P = 10 and for long times the Newman scheme is the best",
  !       and "for large P and suﬃciently short times, the Pade scheme is much better" probably referring to P > 20.
  !
  !       The relation between coefficient a_p and pole sigma_p follows the universal identity a_p = sigma_p^2 .
  !       * Pade scheme:
  !         this leads to a trivial solution since all a_p == 1 and thus sigma_p = 1.
  !       * Chebyshev: the coefficients are
  !           a_p = 1/2 * (1 + cos((2 p - 1) * PI / (2 P)) )
  !         and the corresponding poles become
  !           sigma_p = cos((2 p - 1) * PI / (4 P))
  !       * Newman: it follows that for coefficients
  !           a_p = exp( -p / sqrt(P))
  !         the poles must be at
  !           sigma_p = exp( -p / (2 * sqrt(P) ))
  !
  ! initialize coefficients
  if (HW_COEFF_SCHEME == HW_COEFF_PADE) then
    ! Pade scheme
    if (myrank == 0) then
      write(IMAIN,*) '    using Padé coefficients'
      call flush_IMAIN()
    endif
    ! a_j coefficients
    hw_a_coeff(:) = 1.0_CUSTOM_REAL
    ! normalized pole positions
    hw_sigma(:) = 1.0_CUSTOM_REAL
  else if (HW_COEFF_SCHEME == HW_COEFF_CHEBYSHEV) then
    ! Chebyshev
    if (myrank == 0) then
      write(IMAIN,*) '    using Chebyshev coefficients'
      call flush_IMAIN()
    endif
    do p = 1,HW_ORDER_P
      ! Hagstrom et al. 2007, eq. (41)
      ! for example: P==1: a_1 = 0.5
      !              P==2: a_1 = 0.853, a_2 = 0.146
      hw_a_coeff(p) = 0.5_CUSTOM_REAL * (1.0_CUSTOM_REAL + cos((2 * p - 1) * PI / (2.0_CUSTOM_REAL * HW_ORDER_P)))
      ! normalized pole positions
      hw_sigma(p) = cos((2 * p - 1) * PI / (4.0_CUSTOM_REAL * HW_ORDER_P))
    enddo
  else if (HW_COEFF_SCHEME == HW_COEFF_NEWMAN) then
    ! Newman
    if (myrank == 0) then
      write(IMAIN,*) '    using Newman coefficients'
      call flush_IMAIN()
    endif
    do p = 1,HW_ORDER_P
      ! for example: P==1: a_1 = 0.367
      !              P==2: a_1 = 0.493, a_2 = 0.243
      hw_a_coeff(p) = exp(- p / sqrt( real(HW_ORDER_P,kind=CUSTOM_REAL)))
      ! normalized pole positions
      ! this is probably not correct yet, would need to re-derive from Hagstrom-Warburton?
      ! for lack of documentation, using the same as a_j...
      hw_sigma(p) = exp(- p / ( 2.0_CUSTOM_REAL * sqrt( real(HW_ORDER_P,kind=CUSTOM_REAL))))
    enddo
  else
    stop 'Error unknown H-W ABC coefficient scheme'
  endif

  ! reference angular frequency
  ! note: when using the dominant source frequency, it centers the higher-order corrections
  !       on the most energetic portion of the wavefield
  if (myrank == 0) then
    write(IMAIN,*) '    reference frequency (from f0_FOR_PML): ',sngl(f0_FOR_PML),'(Hz)'
    call flush_IMAIN()
  endif
  hw_omega_c = TWO_PI * f0_FOR_PML

  ! allocate auxiliary state variables
  allocate(hw_phi(3, NGLLX, NGLLY, num_abs_boundary_faces, HW_ORDER_P), &
           hw_phi_dot(3, NGLLX, NGLLY, num_abs_boundary_faces, HW_ORDER_P), &
           hw_phi_dotdot(3, NGLLX, NGLLY, num_abs_boundary_faces, HW_ORDER_P), stat=ier)
  if (ier /= 0) stop 'Error allocating H-W ABC state variables'
  hw_phi(:,:,:,:,:) = 0.0_CUSTOM_REAL
  hw_phi_dot(:,:,:,:,:) = 0.0_CUSTOM_REAL
  hw_phi_dotdot(:,:,:,:,:) = 0.0_CUSTOM_REAL

  ! H-W boundary has stability problems at corners and edges
  !
  ! allocate corner weights: 1 / (number of elastic boundary faces sharing each GLL node)
  ! this prevents double-counting the phi1 traction at edge and corner nodes
  allocate(hw_corner_weight(NGLLX, NGLLY, num_abs_boundary_faces), stat=ier)
  if (ier /= 0) stop 'Error allocating H-W ABC corner weights'
  hw_corner_weight(:,:,:) = 1.0_CUSTOM_REAL

  ! to check normal on each node to see if this is indeed a corner
  allocate(node_normal(NDIM,NGLOB_AB), stat=ier)
  if (ier /= 0) stop 'Error allocating H-W ABC node_normal array'
  node_normal(:,:) = 0.0_CUSTOM_REAL

  ! count how many elastic boundary faces share each global node
  allocate(node_count(NGLOB_AB), stat=ier)
  if (ier /= 0) stop 'Error allocating H-W ABC node_count'
  node_count(:) = 0

  ! count boundary nodes such that corner nodes have counts > 1
  do iface = 1, num_abs_boundary_faces
    ispec = abs_boundary_ispec(iface)

    ! only for elastic domains
    if (.not. ispec_is_elastic(ispec)) cycle

    do igll = 1, NGLLSQUARE
      i = abs_boundary_ijk(1,igll,iface)
      j = abs_boundary_ijk(2,igll,iface)
      k = abs_boundary_ijk(3,igll,iface)

      iglob = ibool(i,j,k,ispec)

      ! increase node count
      if (node_count(iglob) == 0) then
        ! set normal
        node_normal(1,iglob) = abs_boundary_normal(1,igll,iface)
        node_normal(2,iglob) = abs_boundary_normal(2,igll,iface)
        node_normal(3,iglob) = abs_boundary_normal(3,igll,iface)
        ! increase count
        node_count(iglob) = node_count(iglob) + 1
      else
        ! check normal and if different increase count for corner points
        ! vector dot product
        dot_n = node_normal(1,iglob) * abs_boundary_normal(1,igll,iface) &
                + node_normal(2,iglob) * abs_boundary_normal(2,igll,iface) &
                + node_normal(3,iglob) * abs_boundary_normal(3,igll,iface)
        if (abs(dot_n) < THRESHOLD_DOT_PRODUCT) then
          ! corner
          node_count(iglob) = node_count(iglob) + 1
        endif
      endif
    enddo
  enddo

  ! store corner weights for each face GLL point
  do iface = 1, num_abs_boundary_faces
    ispec = abs_boundary_ispec(iface)

    ! only for elastic domains
    if (.not. ispec_is_elastic(ispec)) cycle

    ! determine face-local 2D index mapping (id1, id2)
    mask_vary(:) = .false.
    do igll = 2, NGLLSQUARE
      if (abs_boundary_ijk(1,igll,iface) /= abs_boundary_ijk(1,1,iface)) mask_vary(1) = .true.
      if (abs_boundary_ijk(2,igll,iface) /= abs_boundary_ijk(2,1,iface)) mask_vary(2) = .true.
      if (abs_boundary_ijk(3,igll,iface) /= abs_boundary_ijk(3,1,iface)) mask_vary(3) = .true.
    enddo
    id1 = 0; id2 = 0
    if (.not. mask_vary(1)) then
      id1 = 2; id2 = 3
    else if (.not. mask_vary(2)) then
      id1 = 1; id2 = 3
    else
      id1 = 1; id2 = 2
    endif

    do igll = 1, NGLLSQUARE
      i = abs_boundary_ijk(1,igll,iface)
      j = abs_boundary_ijk(2,igll,iface)
      k = abs_boundary_ijk(3,igll,iface)

      iglob = ibool(i,j,k,ispec)

      ! local face indexing
      a = abs_boundary_ijk(id1,igll,iface)
      b = abs_boundary_ijk(id2,igll,iface)

      ! corner weights
      if (node_count(iglob) > 1) then
        ! sets weigths to zero to fall back to basic Stacey formulations on corners and edges to be stable
        hw_corner_weight(a, b, iface) = 0.0_CUSTOM_REAL
        ! or damp by weight = 1 / count
        !hw_corner_weight(a, b, iface) = 1.0_CUSTOM_REAL / node_count(iglob)
      endif
    enddo
  enddo

  ! determines a unique normal on shared nodes
  ! by summing normals and normalizing resulting vector
  if (USE_SHARED_NORMALS) then
    node_normal(:,:) = 0.0_CUSTOM_REAL
    node_count(:) = 0

    ! sum normals on shared points
    do iface = 1, num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
      ! only for elastic domains
      if (.not. ispec_is_elastic(ispec)) cycle
      do igll = 1, NGLLSQUARE
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)
        ! accumulate normals and increment count
        node_normal(1,iglob) = node_normal(1,iglob) + abs_boundary_normal(1,igll,iface)
        node_normal(2,iglob) = node_normal(2,iglob) + abs_boundary_normal(2,igll,iface)
        node_normal(3,iglob) = node_normal(3,iglob) + abs_boundary_normal(3,igll,iface)
        node_count(iglob) = node_count(iglob) + 1
      enddo
    enddo

    ! normalize the assembled normals
    do iglob = 1, NGLOB_AB
      if (node_count(iglob) > 0) then
        norm_n = sqrt(node_normal(1,iglob)**2 + node_normal(2,iglob)**2 + node_normal(3,iglob)**2)
        if (norm_n > 1.d-24) then
          node_normal(1,iglob) = node_normal(1,iglob) / norm_n
          node_normal(2,iglob) = node_normal(2,iglob) / norm_n
          node_normal(3,iglob) = node_normal(3,iglob) / norm_n
        endif
      endif
    enddo

    ! update global normal vectors with assembled normals
    do iface = 1, num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
      ! only for elastic domains
      if (.not. ispec_is_elastic(ispec)) cycle
      do igll = 1, NGLLSQUARE
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)
        ! overwrite with assembled normal vector
        abs_boundary_normal(1,igll,iface) = node_normal(1,iglob)
        abs_boundary_normal(2,igll,iface) = node_normal(2,iglob)
        abs_boundary_normal(3,igll,iface) = node_normal(3,iglob)
      enddo
    enddo
  endif

  deallocate(node_count,node_normal)

  ! Build edge and corner topology for H-W compatibility conditions (Fix C, Fix D)
  if (HW_APPLY_EDGE_CORNER_COMPATIBILITY) then
    allocate(face_edges(4, num_abs_boundary_faces), stat=ier)
    if (ier /= 0) stop 'Error allocating face_edges'
    allocate(face_corners(4, num_abs_boundary_faces), stat=ier)
    if (ier /= 0) stop 'Error allocating face_corners'

    ! First pass: store face geometry info
    allocate(face_id1(num_abs_boundary_faces), face_id2(num_abs_boundary_faces), face_id3(num_abs_boundary_faces))
    allocate(face_ispec(num_abs_boundary_faces))
    allocate(face_nx_arr(num_abs_boundary_faces), face_ny_arr(num_abs_boundary_faces), face_nz_arr(num_abs_boundary_faces))

    do iface = 1, num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
      face_ispec(iface) = ispec

      ! only elastic domains
      if (.not. ispec_is_elastic(ispec)) cycle

      ! Determine face-local 2D index mapping (id1, id2, id3)
      mask_vary(:) = .false.
      do igll = 2, NGLLSQUARE
        if (abs_boundary_ijk(1,igll,iface) /= abs_boundary_ijk(1,1,iface)) mask_vary(1) = .true.
        if (abs_boundary_ijk(2,igll,iface) /= abs_boundary_ijk(2,1,iface)) mask_vary(2) = .true.
        if (abs_boundary_ijk(3,igll,iface) /= abs_boundary_ijk(3,1,iface)) mask_vary(3) = .true.
      enddo

      id1 = 0; id2 = 0; id3 = 0
      if (.not. mask_vary(1)) then
        id1 = 2; id2 = 3; id3 = 1
      else if (.not. mask_vary(2)) then
        id1 = 1; id2 = 3; id3 = 2
      else
        id1 = 1; id2 = 2; id3 = 3
      endif

      face_id1(iface) = id1
      face_id2(iface) = id2
      face_id3(iface) = id3

      ! Get face normal
      face_nx_arr(iface) = abs_boundary_normal(1, 1, iface)
      face_ny_arr(iface) = abs_boundary_normal(2, 1, iface)
      face_nz_arr(iface) = abs_boundary_normal(3, 1, iface)
    enddo

    ! Build global node to faces map for edge/corner detection
    allocate(node_face_count(NGLOB_AB))
    node_face_count = 0

    ! First count
    do iface = 1, num_abs_boundary_faces
      ispec = face_ispec(iface)
      ! only elastic domains
      if (.not. ispec_is_elastic(ispec)) cycle

      do igll = 1, NGLLSQUARE
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)
        node_face_count(iglob) = node_face_count(iglob) + 1
      enddo
    enddo

    ! Find max faces per node
    max_faces_per_node = maxval(node_face_count)

    allocate(node_face_list(max_faces_per_node, NGLOB_AB), &
             node_face_idx(NGLOB_AB))
    node_face_list(:,:) = 0
    node_face_idx(:) = 0

    do iface = 1, num_abs_boundary_faces
      ispec = face_ispec(iface)
      ! only elastic domains
      if (.not. ispec_is_elastic(ispec)) cycle

      do igll = 1, NGLLSQUARE
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)
        node_face_idx(iglob) = node_face_idx(iglob) + 1
        node_face_list(node_face_idx(iglob), iglob) = iface
      enddo
    enddo

    ! Now for each face, find adjacent faces for each edge
    do iface = 1, num_abs_boundary_faces
      ispec = face_ispec(iface)
      ! only elastic domains
      if (.not. ispec_is_elastic(ispec)) cycle

      id1 = face_id1(iface)
      id2 = face_id2(iface)
      id3 = face_id3(iface)

      nx = face_nx_arr(iface)
      ny = face_ny_arr(iface)
      nz = face_nz_arr(iface)

      ! For each of the 4 edges of this face, find adjacent face
      ! Edge 1: low id1 (a=1), perpendicular dir = id1
      ! Edge 2: high id1 (a=NGLLX), perpendicular dir = id1
      ! Edge 3: low id2 (b=1), perpendicular dir = id2
      ! Edge 4: high id2 (b=NGLLY), perpendicular dir = id2

      ! Get edge nodes for this face
      edge1_nodes(:) = 0; edge2_nodes(:) = 0; edge3_nodes(:) = 0; edge4_nodes(:) = 0
      do igll = 1, NGLLSQUARE
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)

        a = abs_boundary_ijk(id1,igll,iface)
        b = abs_boundary_ijk(id2,igll,iface)

        ! left/right edge
        if (a == 1) edge1_nodes(b) = iglob
        if (a == NGLLX) edge2_nodes(b) = iglob
        ! bottom/top edge
        if (b == 1) edge3_nodes(a) = iglob
        if (b == NGLLY) edge4_nodes(a) = iglob
      enddo

      ! For each edge, find adjacent face by checking shared nodes
      ! Edge 1: low id1 (a=1)
      call find_adjacent_face(edge1_nodes, iface, 1, face_nx_arr, face_ny_arr, face_nz_arr, &
                              face_id1, face_id2, face_id3, face_ispec, &
                              node_face_list, node_face_count, face_edges)

      ! Edge 2: high id1 (a=NGLLX)
      call find_adjacent_face(edge2_nodes, iface, 2, face_nx_arr, face_ny_arr, face_nz_arr, &
                              face_id1, face_id2, face_id3, face_ispec, &
                              node_face_list, node_face_count, face_edges)

      ! Edge 3: low id2 (b=1)
      call find_adjacent_face(edge3_nodes, iface, 3, face_nx_arr, face_ny_arr, face_nz_arr, &
                              face_id1, face_id2, face_id3, face_ispec, &
                              node_face_list, node_face_count, face_edges)

      ! Edge 4: high id2 (b=NGLLY)
      call find_adjacent_face(edge4_nodes, iface, 4, face_nx_arr, face_ny_arr, face_nz_arr, &
                              face_id1, face_id2, face_id3, face_ispec, &
                              node_face_list, node_face_count, face_edges)
    enddo

    ! Build corners (4 per face)
    do iface = 1, num_abs_boundary_faces
      ispec = face_ispec(iface)
      if (.not. ispec_is_elastic(ispec)) cycle

      id1 = face_id1(iface)
      id2 = face_id2(iface)

      ! Corner 1: (a=1, b=1) - edges 1 & 3
      face_corners(1, iface)%i_edge_t1 = 1
      face_corners(1, iface)%i_edge_t2 = 1
      face_corners(1, iface)%adj_face_t1 = face_edges(1, iface)%adj_face
      face_corners(1, iface)%adj_face_t2 = face_edges(3, iface)%adj_face
      face_corners(1, iface)%sigma_t1 = face_edges(1, iface)%sigma
      face_corners(1, iface)%sigma_t2 = face_edges(3, iface)%sigma

      ! Corner 2: (a=NGLLX, b=1) - edges 2 & 3
      face_corners(2, iface)%i_edge_t1 = NGLLX
      face_corners(2, iface)%i_edge_t2 = 1
      face_corners(2, iface)%adj_face_t1 = face_edges(2, iface)%adj_face
      face_corners(2, iface)%adj_face_t2 = face_edges(3, iface)%adj_face
      face_corners(2, iface)%sigma_t1 = face_edges(2, iface)%sigma
      face_corners(2, iface)%sigma_t2 = face_edges(3, iface)%sigma

      ! Corner 3: (a=1, b=NGLLY) - edges 1 & 4
      face_corners(3, iface)%i_edge_t1 = 1
      face_corners(3, iface)%i_edge_t2 = NGLLY
      face_corners(3, iface)%adj_face_t1 = face_edges(1, iface)%adj_face
      face_corners(3, iface)%adj_face_t2 = face_edges(4, iface)%adj_face
      face_corners(3, iface)%sigma_t1 = face_edges(1, iface)%sigma
      face_corners(3, iface)%sigma_t2 = face_edges(4, iface)%sigma

      ! Corner 4: (a=NGLLX, b=NGLLY) - edges 2 & 4
      face_corners(4, iface)%i_edge_t1 = NGLLX
      face_corners(4, iface)%i_edge_t2 = NGLLY
      face_corners(4, iface)%adj_face_t1 = face_edges(2, iface)%adj_face
      face_corners(4, iface)%adj_face_t2 = face_edges(4, iface)%adj_face
      face_corners(4, iface)%sigma_t1 = face_edges(2, iface)%sigma
      face_corners(4, iface)%sigma_t2 = face_edges(4, iface)%sigma

      ! Set component Vp/Vs for Strang splitting
      ! Step 1 (t2 direction): normal=Vs, t1=Vs, t2=Vp
      ! Step 2 (t1 direction): normal=Vs, t1=Vp, t2=Vs
      face_corners(1, iface)%comp_vp_step1 = (/ 0, 0, 1 /)
      face_corners(2, iface)%comp_vp_step1 = (/ 0, 0, 1 /)
      face_corners(3, iface)%comp_vp_step1 = (/ 0, 0, 1 /)
      face_corners(4, iface)%comp_vp_step1 = (/ 0, 0, 1 /)

      face_corners(1, iface)%comp_vp_step2 = (/ 0, 1, 0 /)
      face_corners(2, iface)%comp_vp_step2 = (/ 0, 1, 0 /)
      face_corners(3, iface)%comp_vp_step2 = (/ 0, 1, 0 /)
      face_corners(4, iface)%comp_vp_step2 = (/ 0, 1, 0 /)
    enddo

    deallocate(face_id1, face_id2, face_id3, face_ispec)
    deallocate(face_nx_arr, face_ny_arr, face_nz_arr)
    deallocate(node_face_count, node_face_list, node_face_idx)
  endif   ! HW_APPLY_EDGE_CORNER_COMPATIBILITY

  ! For orthogonal basis
  ! 2D shape functions and their derivatives
  allocate(shape2D_x(NGNOD2D,NGLLY,NGLLZ), &
           shape2D_y(NGNOD2D,NGLLX,NGLLZ), &
           shape2D_z(NGNOD2D,NGLLX,NGLLY),stat=ier)
  if (ier /= 0) stop 'error allocating array shape2D_x etc.'
  shape2D_x(:,:,:) = 0.d0; shape2D_y(:,:,:) = 0.d0; shape2D_z(:,:,:) = 0.d0

  allocate(dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ), &
           dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ), &
           dershape2D_z(NDIM2D,NGNOD2D,NGLLX,NGLLY),stat=ier)
  if (ier /= 0) stop 'error allocating array dershape2D_x etc.'
  dershape2D_x(:,:,:,:) = 0.d0; dershape2D_y(:,:,:,:) = 0.d0; dershape2D_z(:,:,:,:) = 0.d0

  ! get the 2-D shape functions
  call get_shape2D(shape2D_x,dershape2D_x,yigll,zigll,NGLLY,NGLLZ,NGNOD,NGNOD2D)
  call get_shape2D(shape2D_y,dershape2D_y,xigll,zigll,NGLLX,NGLLZ,NGNOD,NGNOD2D)
  call get_shape2D(shape2D_z,dershape2D_z,xigll,yigll,NGLLX,NGLLY,NGNOD,NGNOD2D)

  ! free temporary arrays
  deallocate(shape2D_x,shape2D_y,shape2D_z)

  end subroutine allocate_hw_abc

!
!-------------------------------------------------------------------------------------------------
!

  subroutine deallocate_hw_abc()

  implicit none

  ! free memory
  if (allocated(hw_a_coeff)) deallocate(hw_a_coeff,hw_sigma)
  if (allocated(hw_phi)) deallocate(hw_phi,hw_phi_dot,hw_phi_dotdot)
  if (allocated(hw_corner_weight)) deallocate(hw_corner_weight)
  if (allocated(dershape2D_x)) deallocate(dershape2D_x,dershape2D_y,dershape2D_z)
  if (allocated(face_edges)) deallocate(face_edges)
  if (allocated(face_corners)) deallocate(face_corners)

  end subroutine deallocate_hw_abc


!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_adjacent_face(edge_nodes, iface, edge_num, &
                                face_nx_arr, face_ny_arr, face_nz_arr, &
                                face_id1, face_id2, face_id3, face_ispec, &
                                node_face_list, node_face_count, face_edges)

  ! Find adjacent face sharing this edge by matching global nodes
  ! edge_num: 1=low id1, 2=high id1, 3=low id2, 4=high id2

  use specfem_par_elastic, only: ispec_is_elastic

  implicit none

  integer, intent(in) :: edge_nodes(NGLLX)
  integer, intent(in) :: iface, edge_num
  real(kind=CUSTOM_REAL), intent(in) :: face_nx_arr(:), face_ny_arr(:), face_nz_arr(:)
  integer, intent(in) :: face_id1(:), face_id2(:), face_id3(:), face_ispec(:)
  integer, intent(in) :: node_face_list(:,:), node_face_count(:)
  type(edge_info_type), intent(inout) :: face_edges(:,:)

  ! local parameters
  integer :: jface, i, k, kk, node_k, matches, max_matches, best_face
  integer :: id1, id2, adj_id3
  real(kind=CUSTOM_REAL) :: nx, ny, nz, adj_nx, adj_ny, adj_nz

  max_matches = 0
  best_face = -1

  ! Look at faces sharing the first node of the edge
  if (node_face_count(edge_nodes(1)) > 1) then
    do i = 1, node_face_count(edge_nodes(1))
      jface = node_face_list(i, edge_nodes(1))
      if (jface == iface) cycle
      if (.not. ispec_is_elastic(face_ispec(jface))) cycle

      ! Count shared nodes
      matches = 0
      do k = 1, NGLLX
        node_k = edge_nodes(k)
        ! Check if this node belongs to jface
        do kk = 1, node_face_count(node_k)
          if (node_face_list(kk, node_k) == jface) then
            matches = matches + 1
            exit
          endif
        enddo
      enddo

      if (matches > max_matches) then
        max_matches = matches
        best_face = jface
      endif
    enddo
  endif

  if (best_face > 0 .and. max_matches >= NGLLX - 1) then
    face_edges(edge_num, iface)%adj_face = best_face

    ! Determine sigma: sign of adjacent face normal in perpendicular direction
    ! Get adjacent face normal
    adj_nx = face_nx_arr(best_face)
    adj_ny = face_ny_arr(best_face)
    adj_nz = face_nz_arr(best_face)

    ! Get this face's normal
    nx = face_nx_arr(iface)
    ny = face_ny_arr(iface)
    nz = face_nz_arr(iface)

    ! The perpendicular direction index is id_perp
    ! The adjacent face's id3 (constant index) should match this face's id_perp
    adj_id3 = face_id3(best_face)

    ! Check if the adjacent face's normal aligns with positive or negative
    ! direction of this face's perpendicular tangential direction
    ! The perpendicular tangential vector is along id_perp
    ! For edge_num 1 (low id1) or 3 (low id2), sigma = -1 if adjacent normal
    ! points in positive id_perp direction
    ! For edge_num 2 (high id1) or 4 (high id2), sigma = +1 if adjacent normal
    ! points in positive id_perp direction

    ! For edge 1 (low side) and 3 (low side): if adjacent normal points
    ! opposite to this face's normal in the perp direction, sigma = -1
    ! For edge 2 (high side) and 4 (high side): if adjacent normal points
    ! same as this face's normal in the perp direction, sigma = +1

    if (edge_num == 1 .or. edge_num == 3) then
      ! Low side: adjacent face normal should point in negative id_perp direction
      face_edges(edge_num, iface)%sigma = -1
    else
      ! High side: adjacent face normal should point in positive id_perp direction
      face_edges(edge_num, iface)%sigma = 1
    endif

    ! Compute Jacobian in perpendicular direction
    ! This is the surface metric component sqrt(a_perp_perp) at the edge
    ! For simplicity, we'll compute it during the update
    !face_edges(edge_num, iface)%jacobian_perp = 1.0_CUSTOM_REAL

    ! Set component Vp/Vs assignment
    ! Component parallel to adjacent face normal -> Vp, others -> Vs
    ! adjacent face normal = (adj_nx, adj_ny, adj_nz)
    ! this face components: 1=normal (nx,ny,nz), 2=t1 (along id1), 3=t2 (along id2)
    ! t1 is along id1 direction, t2 along id2 direction
    id1 = face_id1(iface)
    id2 = face_id2(iface)

    ! initializes vp components
    face_edges(edge_num, iface)%comp_vp = (/ 0, 0, 0 /)

    ! normal component (always Vs for edge condition since it's tangential to adjacent face)
    face_edges(edge_num, iface)%comp_vp(1) = 0
    ! t1 component: Vp if id1 == id_perp of adjacent face (i.e., t1 is normal to adjacent face)
    face_edges(edge_num, iface)%comp_vp(2) = 0
    ! t2 component: Vp if id2 == id_perp of adjacent face
    face_edges(edge_num, iface)%comp_vp(3) = 0

    ! The adjacent face's constant index is adj_id3
    ! If adj_id3 == id1, then t1 is normal to adjacent face -> Vp
    ! If adj_id3 == id2, then t2 is normal to adjacent face -> Vp
    if (adj_id3 == id1) then
      face_edges(edge_num, iface)%comp_vp(2) = 1
    else if (adj_id3 == id2) then
      face_edges(edge_num, iface)%comp_vp(3) = 1
    endif

    ! Store id_perp for the edge (which tangential direction is perpendicular to edge)
    ! Edges 1,2 (constant a) have perp = id1 (1)
    ! Edges 3,4 (constant b) have perp = id2 (2)
    ! index where edge is
    if (edge_num == 1) then
      face_edges(edge_num, iface)%i_edge = 1
      face_edges(edge_num, iface)%id_perp = 1    ! edge at constant a, perp is a (id1)
    else if (edge_num == 2) then
      face_edges(edge_num, iface)%i_edge = NGLLX
      face_edges(edge_num, iface)%id_perp = 1
    else if (edge_num == 3) then
      face_edges(edge_num, iface)%i_edge = 1
      face_edges(edge_num, iface)%id_perp = 2    ! edge at constant b, perp is b (id2)
    else
      face_edges(edge_num, iface)%i_edge = NGLLY
      face_edges(edge_num, iface)%id_perp = 2
    endif

  else
    ! has no adjacent face
    face_edges(edge_num, iface)%adj_face = -1
    face_edges(edge_num, iface)%sigma = 0
  endif

  end subroutine find_adjacent_face


!
!-------------------------------------------------------------------------------------------------
!

  subroutine hw_apply_edge_compatibility_newmark(phi_n, phi_t1, phi_t2, &
                                                 phi_n_dot, phi_t1_dot, phi_t2_dot, &
                                                 phi_n_ddot, phi_t1_ddot, phi_t2_ddot, &
                                                 phi_ini_n, phi_ini_t1, phi_ini_t2, &
                                                 phi_ini_n_dot, phi_ini_t1_dot, phi_ini_t2_dot, &
                                                 phi_ini_n_ddot, phi_ini_t1_ddot, phi_ini_t2_ddot, &
                                                 Vp, Vs, a_j, deltat, &
                                                 hprime_perp, jacobian_perp, &
                                                 NGLL, i_edge, sigma, id_perp, &
                                                 comp_vp)

  ! Apply edge compatibility condition using Newmark scheme
  ! PDE: ∂ₜφ + V*a_j*∂_⊥φ = 0  (where V = Vp or Vs depending on component)
  ! Newmark: Predictor -> Spatial derivative -> Solve for accel -> Corrector
  !
  ! id_perp: which tangential direction (1=i, 2=j, 3=k) is perpendicular to the edge (1=id1/a, 2=id2/b)
  !   - If id_perp == 1: edge is at constant a (i_edge is a-index), derivative in a direction
  !   - If id_perp == 2: edge is at constant b (i_edge is b-index), derivative in b direction

  implicit none

  integer, intent(in) :: NGLL
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(inout) :: phi_n, phi_t1, phi_t2
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(inout) :: phi_n_dot, phi_t1_dot, phi_t2_dot
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(inout) :: phi_n_ddot, phi_t1_ddot, phi_t2_ddot

  ! initial states
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(in) :: phi_ini_n, phi_ini_t1, phi_ini_t2
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(in) :: phi_ini_n_dot, phi_ini_t1_dot, phi_ini_t2_dot
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(in) :: phi_ini_n_ddot, phi_ini_t1_ddot, phi_ini_t2_ddot

  real(kind=CUSTOM_REAL), intent(in) :: Vp, Vs, a_j, deltat
  real(kind=CUSTOM_REAL), intent(in) :: hprime_perp(NGLL,NGLL), jacobian_perp
  integer, intent(in) :: i_edge, sigma, id_perp
  integer, intent(in) :: comp_vp(3)

  ! local parameters
  integer :: k, m
  real(kind=CUSTOM_REAL) :: dphi_n_dot, dphi_t1_dot, dphi_t2_dot
  real(kind=CUSTOM_REAL) :: V_n, V_t1, V_t2
  real(kind=CUSTOM_REAL) :: phi_n_pred, phi_t1_pred, phi_t2_pred
  real(kind=CUSTOM_REAL) :: phi_n_dot_pred, phi_t1_dot_pred, phi_t2_dot_pred
  real(kind=CUSTOM_REAL) :: phi_n_ddot_new, phi_t1_ddot_new, phi_t2_ddot_new
  real(kind=CUSTOM_REAL) :: deltatover2, deltatsqover2

  deltatover2 = 0.5_CUSTOM_REAL * deltat
  deltatsqover2 = 0.5_CUSTOM_REAL * deltat * deltat

  ! Determine wave speeds for each component
  V_n = real(1 - comp_vp(1), kind=CUSTOM_REAL) * Vs + real(comp_vp(1), kind=CUSTOM_REAL) * Vp
  V_t1 = real(1 - comp_vp(2), kind=CUSTOM_REAL) * Vs + real(comp_vp(2), kind=CUSTOM_REAL) * Vp
  V_t2 = real(1 - comp_vp(3), kind=CUSTOM_REAL) * Vs + real(comp_vp(3), kind=CUSTOM_REAL) * Vp

  ! Loop over nodes along the edge
  do k = 1, NGLL
    ! Predictor step
    if (id_perp == 1) then
      ! Edge at constant a (i_edge is a-index), perpendicular is a direction
      phi_n_pred = phi_ini_n(i_edge, k) + deltat * phi_ini_n_dot(i_edge, k) + deltatsqover2 * phi_ini_n_ddot(i_edge, k)
      phi_t1_pred = phi_ini_t1(i_edge, k) + deltat * phi_ini_t1_dot(i_edge, k) + deltatsqover2 * phi_ini_t1_ddot(i_edge, k)
      phi_t2_pred = phi_ini_t2(i_edge, k) + deltat * phi_ini_t2_dot(i_edge, k) + deltatsqover2 * phi_ini_t2_ddot(i_edge, k)

      phi_n_dot_pred = phi_ini_n_dot(i_edge, k) + deltatover2 * phi_ini_n_ddot(i_edge, k)
      phi_t1_dot_pred = phi_ini_t1_dot(i_edge, k) + deltatover2 * phi_ini_t1_ddot(i_edge, k)
      phi_t2_dot_pred = phi_ini_t2_dot(i_edge, k) + deltatover2 * phi_ini_t2_ddot(i_edge, k)

      ! Spatial derivative in b direction (perp) of velocity (Newmark scheme)
      dphi_n_dot = 0.0_CUSTOM_REAL
      dphi_t1_dot = 0.0_CUSTOM_REAL
      dphi_t2_dot = 0.0_CUSTOM_REAL
      do m = 1, NGLL  ! perpendicular derivative at (a=i_edge, b=k)
        ! Predicted velocity at interior point (m, k): phi_dot_pred = phi_dot + dt/2 * phi_ddot
        dphi_n_dot = dphi_n_dot + hprime_perp(i_edge, m) * (phi_ini_n_dot(m, k) + deltatover2 * phi_ini_n_ddot(m, k))
        dphi_t1_dot = dphi_t1_dot + hprime_perp(i_edge, m) * (phi_ini_t1_dot(m, k) + deltatover2 * phi_ini_t1_ddot(m, k))
        dphi_t2_dot = dphi_t2_dot + hprime_perp(i_edge, m) * (phi_ini_t2_dot(m, k) + deltatover2 * phi_ini_t2_ddot(m, k))
      enddo
    else
      ! id_perp == 2: Edge at constant b (i_edge is b-index), perpendicular is b direction
      phi_n_pred = phi_ini_n(k, i_edge) + deltat * phi_ini_n_dot(k, i_edge) + deltatsqover2 * phi_ini_n_ddot(k, i_edge)
      phi_t1_pred = phi_ini_t1(k, i_edge) + deltat * phi_ini_t1_dot(k, i_edge) + deltatsqover2 * phi_ini_t1_ddot(k, i_edge)
      phi_t2_pred = phi_ini_t2(k, i_edge) + deltat * phi_ini_t2_dot(k, i_edge) + deltatsqover2 * phi_ini_t2_ddot(k, i_edge)

      phi_n_dot_pred = phi_ini_n_dot(k, i_edge) + deltatover2 * phi_ini_n_ddot(k, i_edge)
      phi_t1_dot_pred = phi_ini_t1_dot(k, i_edge) + deltatover2 * phi_ini_t1_ddot(k, i_edge)
      phi_t2_dot_pred = phi_ini_t2_dot(k, i_edge) + deltatover2 * phi_ini_t2_ddot(k, i_edge)

      ! Spatial derivative in a direction (perp) of velocity (Newmark scheme)
      dphi_n_dot = 0.0_CUSTOM_REAL
      dphi_t1_dot = 0.0_CUSTOM_REAL
      dphi_t2_dot = 0.0_CUSTOM_REAL
      do m = 1, NGLL
        ! Predicted velocity at interior point (m, k): phi_dot_pred = phi_dot + dt/2 * phi_ddot
        dphi_n_dot = dphi_n_dot + hprime_perp(i_edge, m) * (phi_ini_n_dot(k, m) + deltatover2 * phi_ini_n_ddot(k, m))
        dphi_t1_dot = dphi_t1_dot + hprime_perp(i_edge, m) * (phi_ini_t1_dot(k, m) + deltatover2 * phi_ini_t1_ddot(k, m))
        dphi_t2_dot = dphi_t2_dot + hprime_perp(i_edge, m) * (phi_ini_t2_dot(k, m) + deltatover2 * phi_ini_t2_ddot(k, m))
      enddo
    endif

    dphi_n_dot = dphi_n_dot / jacobian_perp
    dphi_t1_dot = dphi_t1_dot / jacobian_perp
    dphi_t2_dot = dphi_t2_dot / jacobian_perp

    ! PDE: ∂ₜₜ φ + V*a_j*σ*∂_⊥ (∂ₜ φ) = 0
    phi_n_ddot_new = -V_n * a_j * real(sigma, kind=CUSTOM_REAL) * dphi_n_dot
    phi_t1_ddot_new = -V_t1 * a_j * real(sigma, kind=CUSTOM_REAL) * dphi_t1_dot
    phi_t2_ddot_new = -V_t2 * a_j * real(sigma, kind=CUSTOM_REAL) * dphi_t2_dot

    ! Corrector
    if (id_perp == 1) then
      ! Edge at constant a (i_edge is a-index), perpendicular is a direction
      phi_n(i_edge, k) = phi_n_pred
      phi_t1(i_edge, k) = phi_t1_pred
      phi_t2(i_edge, k) = phi_t2_pred

      phi_n_dot(i_edge, k) = phi_n_dot_pred + deltatover2 * phi_n_ddot_new
      phi_t1_dot(i_edge, k) = phi_t1_dot_pred + deltatover2 * phi_t1_ddot_new
      phi_t2_dot(i_edge, k) = phi_t2_dot_pred + deltatover2 * phi_t2_ddot_new

      phi_n_ddot(i_edge, k) = phi_n_ddot_new
      phi_t1_ddot(i_edge, k) = phi_t1_ddot_new
      phi_t2_ddot(i_edge, k) = phi_t2_ddot_new
    else
      phi_n(k, i_edge) = phi_n_pred
      phi_t1(k, i_edge) = phi_t1_pred
      phi_t2(k, i_edge) = phi_t2_pred

      phi_n_dot(k, i_edge) = phi_n_dot_pred + deltatover2 * phi_n_ddot_new
      phi_t1_dot(k, i_edge) = phi_t1_dot_pred + deltatover2 * phi_t1_ddot_new
      phi_t2_dot(k, i_edge) = phi_t2_dot_pred + deltatover2 * phi_t2_ddot_new

      phi_n_ddot(k, i_edge) = phi_n_ddot_new
      phi_t1_ddot(k, i_edge) = phi_t1_ddot_new
      phi_t2_ddot(k, i_edge) = phi_t2_ddot_new
    endif
  enddo

  end subroutine hw_apply_edge_compatibility_newmark


!
!-------------------------------------------------------------------------------------------------
!

  subroutine hw_apply_corner_compatibility_newmark(phi_n, phi_t1, phi_t2, &
                                                   phi_n_dot, phi_t1_dot, phi_t2_dot, &
                                                   phi_n_ddot, phi_t1_ddot, phi_t2_ddot, &
                                                   phi_ini_n, phi_ini_t1, phi_ini_t2, &
                                                   phi_ini_n_dot, phi_ini_t1_dot, phi_ini_t2_dot, &
                                                   phi_ini_n_ddot, phi_ini_t1_ddot, phi_ini_t2_ddot, &
                                                   Vp, Vs, a_j, deltat, &
                                                   hprime_t1, hprime_t2, jacobian_t1, jacobian_t2, &
                                                   NGLL, i_edge_t1, i_edge_t2, &
                                                   sigma_t1, sigma_t2, &
                                                   comp_vp_step1, comp_vp_step2)

  ! Apply corner compatibility using Strang splitting + Newmark
  ! PDE: (∂ₜ + V*a_j*∂_t1)(∂ₜ + V*a_j*∂_t2)φ = 0
  ! Step 1: (∂ₜ + V*a_j*∂_t2)φ = 0
  ! Step 2: (∂ₜ + V*a_j*∂_t1)φ = 0

  implicit none

  integer, intent(in) :: NGLL
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(inout) :: phi_n, phi_t1, phi_t2
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(inout) :: phi_n_dot, phi_t1_dot, phi_t2_dot
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(inout) :: phi_n_ddot, phi_t1_ddot, phi_t2_ddot

  ! initial states
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(in) :: phi_ini_n, phi_ini_t1, phi_ini_t2
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(in) :: phi_ini_n_dot, phi_ini_t1_dot, phi_ini_t2_dot
  real(kind=CUSTOM_REAL), dimension(NGLL,NGLL), intent(in) :: phi_ini_n_ddot, phi_ini_t1_ddot, phi_ini_t2_ddot

  real(kind=CUSTOM_REAL), intent(in) :: Vp, Vs, a_j, deltat
  real(kind=CUSTOM_REAL), intent(in) :: hprime_t1(NGLL,NGLL), hprime_t2(NGLL,NGLL)
  real(kind=CUSTOM_REAL), intent(in) :: jacobian_t1, jacobian_t2
  integer, intent(in) :: i_edge_t1, i_edge_t2
  integer, intent(in) :: sigma_t1, sigma_t2
  integer, intent(in) :: comp_vp_step1(3), comp_vp_step2(3)

  ! local parameters
  real(kind=CUSTOM_REAL) :: dphi_n_t2, dphi_t1_t2, dphi_t2_t2
  real(kind=CUSTOM_REAL) :: dphi_n_t1, dphi_t1_t1, dphi_t2_t1
  real(kind=CUSTOM_REAL) :: V_n_s1, V_t1_s1, V_t2_s1
  real(kind=CUSTOM_REAL) :: V_n_s2, V_t1_s2, V_t2_s2
  real(kind=CUSTOM_REAL) :: phi_n_star, phi_t1_star, phi_t2_star
  real(kind=CUSTOM_REAL) :: phi_n_dot_star, phi_t1_dot_star, phi_t2_dot_star
  real(kind=CUSTOM_REAL) :: phi_n_ddot_star, phi_t1_ddot_star, phi_t2_ddot_star
  real(kind=CUSTOM_REAL) :: deltatover2, deltatsqover2
  integer :: m

  deltatover2 = 0.5_CUSTOM_REAL * deltat
  deltatsqover2 = 0.5_CUSTOM_REAL * deltat * deltat

  ! === STEP 1: t2 direction (∂ₜ + V*a_j*∂_t2)φ = 0 ===
  ! Wave speeds for step 1: normal=Vs, t1=Vs, t2=Vp
  V_n_s1 = real(1 - comp_vp_step1(1), kind=CUSTOM_REAL) * Vs + real(comp_vp_step1(1), kind=CUSTOM_REAL) * Vp
  V_t1_s1 = real(1 - comp_vp_step1(2), kind=CUSTOM_REAL) * Vs + real(comp_vp_step1(2), kind=CUSTOM_REAL) * Vp
  V_t2_s1 = real(1 - comp_vp_step1(3), kind=CUSTOM_REAL) * Vs + real(comp_vp_step1(3), kind=CUSTOM_REAL) * Vp

  ! Predictor for step 1
  phi_n_star = phi_ini_n(i_edge_t1, i_edge_t2) + deltat * phi_ini_n_dot(i_edge_t1, i_edge_t2) &
               + deltatsqover2 * phi_ini_n_ddot(i_edge_t1, i_edge_t2)
  phi_t1_star = phi_ini_t1(i_edge_t1, i_edge_t2) + deltat * phi_ini_t1_dot(i_edge_t1, i_edge_t2) &
                + deltatsqover2 * phi_ini_t1_ddot(i_edge_t1, i_edge_t2)
  phi_t2_star = phi_ini_t2(i_edge_t1, i_edge_t2) + deltat * phi_ini_t2_dot(i_edge_t1, i_edge_t2) &
                + deltatsqover2 * phi_ini_t2_ddot(i_edge_t1, i_edge_t2)

  phi_n_dot_star = phi_ini_n_dot(i_edge_t1, i_edge_t2) + deltatover2 * phi_ini_n_ddot(i_edge_t1, i_edge_t2)
  phi_t1_dot_star = phi_ini_t1_dot(i_edge_t1, i_edge_t2) + deltatover2 * phi_ini_t1_ddot(i_edge_t1, i_edge_t2)
  phi_t2_dot_star = phi_ini_t2_dot(i_edge_t1, i_edge_t2) + deltatover2 * phi_ini_t2_ddot(i_edge_t1, i_edge_t2)

  ! Spatial derivative in t2 direction at predicted state
  dphi_n_t2 = 0.0_CUSTOM_REAL
  dphi_t1_t2 = 0.0_CUSTOM_REAL
  dphi_t2_t2 = 0.0_CUSTOM_REAL
  do m = 1, NGLL
    ! Predicted velocity at interior point (m, k): phi_dot_pred = phi_dot + dt/2 * phi_ddot
    dphi_n_t2 = dphi_n_t2 + hprime_t2(i_edge_t2, m) &
                * (phi_ini_n_dot(i_edge_t1, m) + deltatover2 * phi_ini_n_ddot(i_edge_t1, m))
    dphi_t1_t2 = dphi_t1_t2 + hprime_t2(i_edge_t2, m) &
                 * (phi_ini_t1_dot(i_edge_t1, m) + deltatover2 * phi_ini_t1_ddot(i_edge_t1, m))
    dphi_t2_t2 = dphi_t2_t2 + hprime_t2(i_edge_t2, m) &
                 * (phi_ini_t2_dot(i_edge_t1, m) + deltatover2 * phi_ini_t2_ddot(i_edge_t1, m))
  enddo
  dphi_n_t2 = dphi_n_t2 / jacobian_t2
  dphi_t1_t2 = dphi_t1_t2 / jacobian_t2
  dphi_t2_t2 = dphi_t2_t2 / jacobian_t2

  ! Acceleration for step 1
  phi_n_ddot_star = -V_n_s1 * a_j * real(sigma_t2, kind=CUSTOM_REAL) * dphi_n_t2
  phi_t1_ddot_star = -V_t1_s1 * a_j * real(sigma_t2, kind=CUSTOM_REAL) * dphi_t1_t2
  phi_t2_ddot_star = -V_t2_s1 * a_j * real(sigma_t2, kind=CUSTOM_REAL) * dphi_t2_t2

  ! Corrector for step 1
  phi_n_star = phi_n_star
  phi_t1_star = phi_t1_star
  phi_t2_star = phi_t2_star

  phi_n_dot_star = phi_n_dot_star + deltatover2 * phi_n_ddot_star
  phi_t1_dot_star = phi_t1_dot_star + deltatover2 * phi_t1_ddot_star
  phi_t2_dot_star = phi_t2_dot_star + deltatover2 * phi_t2_ddot_star

  ! === STEP 2: t1 direction (∂ₜ + V*a_j*∂_t1)φ = 0 ===
  ! Wave speeds for step 2: normal=Vs, t1=Vp, t2=Vs
  V_n_s2 = real(1 - comp_vp_step2(1), kind=CUSTOM_REAL) * Vs + real(comp_vp_step2(1), kind=CUSTOM_REAL) * Vp
  V_t1_s2 = real(1 - comp_vp_step2(2), kind=CUSTOM_REAL) * Vs + real(comp_vp_step2(2), kind=CUSTOM_REAL) * Vp
  V_t2_s2 = real(1 - comp_vp_step2(3), kind=CUSTOM_REAL) * Vs + real(comp_vp_step2(3), kind=CUSTOM_REAL) * Vp

  ! Predictor for step 2 (using step 1 corrected values)
  phi_n(i_edge_t1, i_edge_t2) = phi_n_star
  phi_t1(i_edge_t1, i_edge_t2) = phi_t1_star
  phi_t2(i_edge_t1, i_edge_t2) = phi_t2_star

  phi_n_dot(i_edge_t1, i_edge_t2) = phi_n_dot_star
  phi_t1_dot(i_edge_t1, i_edge_t2) = phi_t1_dot_star
  phi_t2_dot(i_edge_t1, i_edge_t2) = phi_t2_dot_star

  ! Spatial derivative in t1 direction
  dphi_n_t1 = 0.0_CUSTOM_REAL
  dphi_t1_t1 = 0.0_CUSTOM_REAL
  dphi_t2_t1 = 0.0_CUSTOM_REAL
  do m = 1, NGLL
    ! uses velocity (phi_*_dot) values corrected from step 1
    dphi_n_t1 = dphi_n_t1 + hprime_t1(i_edge_t1, m) * phi_n_dot(m, i_edge_t2)
    dphi_t1_t1 = dphi_t1_t1 + hprime_t1(i_edge_t1, m) * phi_t1_dot(m, i_edge_t2)
    dphi_t2_t1 = dphi_t2_t1 + hprime_t1(i_edge_t1, m) * phi_t2_dot(m, i_edge_t2)
  enddo
  dphi_n_t1 = dphi_n_t1 / jacobian_t1
  dphi_t1_t1 = dphi_t1_t1 / jacobian_t1
  dphi_t2_t1 = dphi_t2_t1 / jacobian_t1

  ! Acceleration for step 2
  phi_n_ddot(i_edge_t1, i_edge_t2) = -V_n_s2 * a_j * real(sigma_t1, kind=CUSTOM_REAL) * dphi_n_t1
  phi_t1_ddot(i_edge_t1, i_edge_t2) = -V_t1_s2 * a_j * real(sigma_t1, kind=CUSTOM_REAL) * dphi_t1_t1
  phi_t2_ddot(i_edge_t1, i_edge_t2) = -V_t2_s2 * a_j * real(sigma_t1, kind=CUSTOM_REAL) * dphi_t2_t1

  ! Final corrector for step 2
  phi_n_dot(i_edge_t1, i_edge_t2) = phi_n_dot(i_edge_t1, i_edge_t2) + deltatover2 * phi_n_ddot(i_edge_t1, i_edge_t2)
  phi_t1_dot(i_edge_t1, i_edge_t2) = phi_t1_dot(i_edge_t1, i_edge_t2) + deltatover2 * phi_t1_ddot(i_edge_t1, i_edge_t2)
  phi_t2_dot(i_edge_t1, i_edge_t2) = phi_t2_dot(i_edge_t1, i_edge_t2) + deltatover2 * phi_t2_ddot(i_edge_t1, i_edge_t2)

  end subroutine hw_apply_corner_compatibility_newmark


!
!-------------------------------------------------------------------------------------------------
!

  subroutine determine_face_ortho_basis(iface,ispec,nx_2D,ny_2D,nz_2D,t1x_2D,t1y_2D,t1z_2D,t2x_2D,t2y_2D,t2z_2D)

! determines an orthonormal basis for a given GLL point on a face.
! this computes first the outward pointing normal (n) and the tangent vector along xi-direction (t_xi), and
! then constructs the second tangent vector based on cross-products of n and t_xi.
!
! see also subroutine `compute_jacobian_2D()` in file get_jacobian_boundaries.f90 how the normal gets computed.

  use specfem_par, only: ibool,xstore,ystore,zstore,abs_boundary_ijk

  implicit none

  integer, intent(in) :: iface,ispec
  real(kind=CUSTOM_REAL), intent(out) :: nx_2D(NGLLX,NGLLY),ny_2D(NGLLX,NGLLY),nz_2D(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL), intent(out) :: t1x_2D(NGLLX,NGLLY),t1y_2D(NGLLX,NGLLY),t1z_2D(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL), intent(out) :: t2x_2D(NGLLX,NGLLY),t2y_2D(NGLLX,NGLLY),t2z_2D(NGLLX,NGLLY)

  ! local parameters
  ! face corners
  double precision :: xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)
  double precision :: dershape2D(NDIM2D,NGNOD2D,NGLLX,NGLLY)
  double precision :: xxi,xeta,yxi,yeta,zxi,zeta
  double precision :: unx,uny,unz,jacobian
  double precision :: nx,ny,nz,t1x,t1y,t1z,t2x,t2y,t2z,t1_norm
  double precision :: face_n(NDIM),v_tmp(NDIM),tmp
  logical :: mask_vary(3)
  integer :: id1, id2, id3, igll, a, b, ia, face_id, iglob

  ! check that the anchor number is correct
  if (NGNOD2D /= 4 .and. NGNOD2D /= 9) stop 'Surface elements should have 4 or 9 control nodes'

  ! assumes NGLLX == NGLLY == NGLLZ
  if (NGLLX /= NGLLY .or. NGLLX /= NGLLZ) stop 'Surface elements should have NGLLX == NGLLY == NGLLZ'

  ! find the varying local coordinates for the face to map 1D index igll to 2D grid (a,b)
  mask_vary(:) = .false.
  do igll = 2, NGLLSQUARE
    if (abs_boundary_ijk(1, igll, iface) /= abs_boundary_ijk(1, 1, iface)) mask_vary(1) = .true.  ! i-index changes
    if (abs_boundary_ijk(2, igll, iface) /= abs_boundary_ijk(2, 1, iface)) mask_vary(2) = .true.  ! j-index changes
    if (abs_boundary_ijk(3, igll, iface) /= abs_boundary_ijk(3, 1, iface)) mask_vary(3) = .true.  ! k-index changes
  enddo

  ! id1 - first changing index (1==i or 2==j or 3==k)
  ! id2 - second changing index
  ! id3 - constant index
  id1 = 0; id2 = 0; id3 = 0
  if (.not. mask_vary(1)) then
    id1 = 2; id2 = 3; id3 = 1
  else if (.not. mask_vary(2)) then
    id1 = 1; id2 = 3; id3 = 2
  else
    id1 = 1; id2 = 2; id3 = 3
  endif

  ! check constant index to determine face id (1==xmin,2==xmax,3==ymin,4==ymax,5==zbottom,6==ztop)
  face_id = 0
  if (id3 == 1) then
    ! i-index is constant
    if (abs_boundary_ijk(id3, 1, iface) == 1) then
      face_id = 1  ! xmin
    else
      face_id = 2  ! xmax
    endif
    dershape2D(:,:,:,:) = dershape2D_x(:,:,:,:)
  else if (id3 == 2) then
    ! j-index is constant
    if (abs_boundary_ijk(id3, 1, iface) == 1) then
      face_id = 3  ! ymin
    else
      face_id = 4  ! ymax
    endif
    dershape2D(:,:,:,:) = dershape2D_y(:,:,:,:)
  else
    ! k-index is constant
    if (abs_boundary_ijk(id3, 1, iface) == 1) then
      face_id = 5  ! zmin
    else
      face_id = 6  ! zmax
    endif
    dershape2D(:,:,:,:) = dershape2D_z(:,:,:,:)
  endif

  ! set corner points on reference face
  select case (face_id)
  case (1)
    ! xmin
    xelm(1) = xstore( ibool(1,1,1,ispec) )
    yelm(1) = ystore( ibool(1,1,1,ispec) )
    zelm(1) = zstore( ibool(1,1,1,ispec) )
    xelm(2) = xstore( ibool(1,NGLLY,1,ispec) )
    yelm(2) = ystore( ibool(1,NGLLY,1,ispec) )
    zelm(2) = zstore( ibool(1,NGLLY,1,ispec) )
    xelm(3) = xstore( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(3) = ystore( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(3) = zstore( ibool(1,NGLLY,NGLLZ,ispec) )
    xelm(4) = xstore( ibool(1,1,NGLLZ,ispec) )
    yelm(4) = ystore( ibool(1,1,NGLLZ,ispec) )
    zelm(4) = zstore( ibool(1,1,NGLLZ,ispec) )
    if (NGNOD2D == 9) then
      xelm(5) = xstore( ibool(1,MIDY,1,ispec) )
      yelm(5) = ystore( ibool(1,MIDY,1,ispec) )
      zelm(5) = zstore( ibool(1,MIDY,1,ispec) )
      xelm(6) = xstore( ibool(1,NGLLY,MIDZ,ispec) )
      yelm(6) = ystore( ibool(1,NGLLY,MIDZ,ispec) )
      zelm(6) = zstore( ibool(1,NGLLY,MIDZ,ispec) )
      xelm(7) = xstore( ibool(1,MIDY,NGLLZ,ispec) )
      yelm(7) = ystore( ibool(1,MIDY,NGLLZ,ispec) )
      zelm(7) = zstore( ibool(1,MIDY,NGLLZ,ispec) )
      xelm(8) = xstore( ibool(1,1,MIDZ,ispec) )
      yelm(8) = ystore( ibool(1,1,MIDZ,ispec) )
      zelm(8) = zstore( ibool(1,1,MIDZ,ispec) )
      xelm(9) = xstore( ibool(1,MIDY,MIDZ,ispec) )
      yelm(9) = ystore( ibool(1,MIDY,MIDZ,ispec) )
      zelm(9) = zstore( ibool(1,MIDY,MIDZ,ispec) )
    endif
  case (2)
    ! xmax
    xelm(1) = xstore( ibool(NGLLX,1,1,ispec) )
    yelm(1) = ystore( ibool(NGLLX,1,1,ispec) )
    zelm(1) = zstore( ibool(NGLLX,1,1,ispec) )
    xelm(2) = xstore( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(2) = ystore( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(2) = zstore( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(3) = xstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3) = ystore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3) = zstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4) = xstore( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(4) = ystore( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(4) = zstore( ibool(NGLLX,1,NGLLZ,ispec) )
    if (NGNOD2D == 9) then
      xelm(5) = xstore( ibool(NGLLX,MIDY,1,ispec) )
      yelm(5) = ystore( ibool(NGLLX,MIDY,1,ispec) )
      zelm(5) = zstore( ibool(NGLLX,MIDY,1,ispec) )
      xelm(6) = xstore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      yelm(6) = ystore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      zelm(6) = zstore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      xelm(7) = xstore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      yelm(7) = ystore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      zelm(7) = zstore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      xelm(8) = xstore( ibool(NGLLX,1,MIDZ,ispec) )
      yelm(8) = ystore( ibool(NGLLX,1,MIDZ,ispec) )
      zelm(8) = zstore( ibool(NGLLX,1,MIDZ,ispec) )
      xelm(9) = xstore( ibool(NGLLX,MIDY,MIDZ,ispec) )
      yelm(9) = ystore( ibool(NGLLX,MIDY,MIDZ,ispec) )
      zelm(9) = zstore( ibool(NGLLX,MIDY,MIDZ,ispec) )
    endif
  case (3)
    ! ymin
    xelm(1) = xstore( ibool(1,1,1,ispec) )
    yelm(1) = ystore( ibool(1,1,1,ispec) )
    zelm(1) = zstore( ibool(1,1,1,ispec) )
    xelm(2) = xstore( ibool(NGLLX,1,1,ispec) )
    yelm(2) = ystore( ibool(NGLLX,1,1,ispec) )
    zelm(2) = zstore( ibool(NGLLX,1,1,ispec) )
    xelm(3) = xstore( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(3) = ystore( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(3) = zstore( ibool(NGLLX,1,NGLLZ,ispec) )
    xelm(4) = xstore( ibool(1,1,NGLLZ,ispec) )
    yelm(4) = ystore( ibool(1,1,NGLLZ,ispec) )
    zelm(4) = zstore( ibool(1,1,NGLLZ,ispec) )
    if (NGNOD2D == 9) then
      xelm(5) = xstore( ibool(MIDX,1,1,ispec) )
      yelm(5) = ystore( ibool(MIDX,1,1,ispec) )
      zelm(5) = zstore( ibool(MIDX,1,1,ispec) )
      xelm(6) = xstore( ibool(NGLLX,1,MIDZ,ispec) )
      yelm(6) = ystore( ibool(NGLLX,1,MIDZ,ispec) )
      zelm(6) = zstore( ibool(NGLLX,1,MIDZ,ispec) )
      xelm(7) = xstore( ibool(MIDX,1,NGLLZ,ispec) )
      yelm(7) = ystore( ibool(MIDX,1,NGLLZ,ispec) )
      zelm(7) = zstore( ibool(MIDX,1,NGLLZ,ispec) )
      xelm(8) = xstore( ibool(1,1,MIDZ,ispec) )
      yelm(8) = ystore( ibool(1,1,MIDZ,ispec) )
      zelm(8) = zstore( ibool(1,1,MIDZ,ispec) )
      xelm(9) = xstore( ibool(MIDX,1,MIDZ,ispec) )
      yelm(9) = ystore( ibool(MIDX,1,MIDZ,ispec) )
      zelm(9) = zstore( ibool(MIDX,1,MIDZ,ispec) )
    endif
  case (4)
    ! ymax
    xelm(1) = xstore( ibool(1,NGLLY,1,ispec) )
    yelm(1) = ystore( ibool(1,NGLLY,1,ispec) )
    zelm(1) = zstore( ibool(1,NGLLY,1,ispec) )
    xelm(2) = xstore( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(2) = ystore( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(2) = zstore( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(3) = xstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3) = ystore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3) = zstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4) = xstore( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(4) = ystore( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(4) = zstore( ibool(1,NGLLY,NGLLZ,ispec) )
    if (NGNOD2D == 9) then
      xelm(5) = xstore( ibool(MIDX,NGLLY,1,ispec) )
      yelm(5) = ystore( ibool(MIDX,NGLLY,1,ispec) )
      zelm(5) = zstore( ibool(MIDX,NGLLY,1,ispec) )
      xelm(6) = xstore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      yelm(6) = ystore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      zelm(6) = zstore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      xelm(7) = xstore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      yelm(7) = ystore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      zelm(7) = zstore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      xelm(8) = xstore( ibool(1,NGLLY,MIDZ,ispec) )
      yelm(8) = ystore( ibool(1,NGLLY,MIDZ,ispec) )
      zelm(8) = zstore( ibool(1,NGLLY,MIDZ,ispec) )
      xelm(9) = xstore( ibool(MIDX,NGLLY,MIDZ,ispec) )
      yelm(9) = ystore( ibool(MIDX,NGLLY,MIDZ,ispec) )
      zelm(9) = zstore( ibool(MIDX,NGLLY,MIDZ,ispec) )
    endif
  case (5)
    ! zmin
    xelm(1) = xstore( ibool(1,1,1,ispec) )
    yelm(1) = ystore( ibool(1,1,1,ispec) )
    zelm(1) = zstore( ibool(1,1,1,ispec) )
    xelm(2) = xstore( ibool(NGLLX,1,1,ispec) )
    yelm(2) = ystore( ibool(NGLLX,1,1,ispec) )
    zelm(2) = zstore( ibool(NGLLX,1,1,ispec) )
    xelm(3) = xstore( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(3) = ystore( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(3) = zstore( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(4) = xstore( ibool(1,NGLLY,1,ispec) )
    yelm(4) = ystore( ibool(1,NGLLY,1,ispec) )
    zelm(4) = zstore( ibool(1,NGLLY,1,ispec) )
    if (NGNOD2D == 9) then
      xelm(5) = xstore( ibool(MIDX,1,1,ispec) )
      yelm(5) = ystore( ibool(MIDX,1,1,ispec) )
      zelm(5) = zstore( ibool(MIDX,1,1,ispec) )
      xelm(6) = xstore( ibool(NGLLX,MIDY,1,ispec) )
      yelm(6) = ystore( ibool(NGLLX,MIDY,1,ispec) )
      zelm(6) = zstore( ibool(NGLLX,MIDY,1,ispec) )
      xelm(7) = xstore( ibool(MIDX,NGLLY,1,ispec) )
      yelm(7) = ystore( ibool(MIDX,NGLLY,1,ispec) )
      zelm(7) = zstore( ibool(MIDX,NGLLY,1,ispec) )
      xelm(8) = xstore( ibool(1,MIDY,1,ispec) )
      yelm(8) = ystore( ibool(1,MIDY,1,ispec) )
      zelm(8) = zstore( ibool(1,MIDY,1,ispec) )
      xelm(9) = xstore( ibool(MIDX,MIDY,1,ispec) )
      yelm(9) = ystore( ibool(MIDX,MIDY,1,ispec) )
      zelm(9) = zstore( ibool(MIDX,MIDY,1,ispec) )
    endif
  case (6)
    ! zmax
    xelm(1) = xstore( ibool(1,1,NGLLZ,ispec) )
    yelm(1) = ystore( ibool(1,1,NGLLZ,ispec) )
    zelm(1) = zstore( ibool(1,1,NGLLZ,ispec) )
    xelm(2) = xstore( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(2) = ystore( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(2) = zstore( ibool(NGLLX,1,NGLLZ,ispec) )
    xelm(3) = xstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3) = ystore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3) = zstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4) = xstore( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(4) = ystore( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(4) = zstore( ibool(1,NGLLY,NGLLZ,ispec) )
    if (NGNOD2D == 9) then
      xelm(5) = xstore( ibool(MIDX,1,NGLLZ,ispec) )
      yelm(5) = ystore( ibool(MIDX,1,NGLLZ,ispec) )
      zelm(5) = zstore( ibool(MIDX,1,NGLLZ,ispec) )
      xelm(6) = xstore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      yelm(6) = ystore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      zelm(6) = zstore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      xelm(7) = xstore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      yelm(7) = ystore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      zelm(7) = zstore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      xelm(8) = xstore( ibool(1,MIDY,NGLLZ,ispec) )
      yelm(8) = ystore( ibool(1,MIDY,NGLLZ,ispec) )
      zelm(8) = zstore( ibool(1,MIDY,NGLLZ,ispec) )
      xelm(9) = xstore( ibool(MIDX,MIDY,NGLLZ,ispec) )
      yelm(9) = ystore( ibool(MIDX,MIDY,NGLLZ,ispec) )
      zelm(9) = zstore( ibool(MIDX,MIDY,NGLLZ,ispec) )
    endif
  case default
    stop 'Error face id for 2D basis'
  end select

  ! determines initial orientation given by three corners on the face
  ! cross-product of vectors from corner 1 to corner 2 and from corner 1 to corner 3
  face_n(1) =   (yelm(2)-yelm(1)) * (zelm(3)-zelm(1)) - (zelm(2)-zelm(1)) * (yelm(3)-yelm(1))
  face_n(2) = - (xelm(2)-xelm(1)) * (zelm(3)-zelm(1)) + (zelm(2)-zelm(1)) * (xelm(3)-xelm(1))
  face_n(3) =   (xelm(2)-xelm(1)) * (yelm(3)-yelm(1)) - (yelm(2)-yelm(1)) * (xelm(3)-xelm(1))

  tmp = dsqrt( face_n(1)*face_n(1) + face_n(2)*face_n(2) + face_n(3)*face_n(3) )
  if (tmp <= 0.d0) stop 'Error get element face normal'

  face_n(:) = face_n(:) / tmp

  ! checks that this normal direction is outwards of element:
  ! takes additional corner out of face plane and determines scalar product (dot product) to normal
  iglob = 0
  select case (face_id)
  case (1) ! opposite to xmin face
    iglob = ibool(NGLLX,1,1,ispec)
  case (2) ! opposite to xmax face
    iglob = ibool(1,1,1,ispec)
  case (3) ! opposite to ymin face
    iglob = ibool(1,NGLLY,1,ispec)
  case (4) ! opposite to ymax face
    iglob = ibool(1,1,1,ispec)
  case (5) ! opposite to bottom
    iglob = ibool(1,1,NGLLZ,ispec)
  case (6) ! opposite to top
    iglob = ibool(1,1,1,ispec)
  case default
    stop 'Error get element face face_id value'
  end select
  ! vector from corner 1 to this opposite one
  v_tmp(1) = xstore(iglob) - xelm(1)
  v_tmp(2) = ystore(iglob) - yelm(1)
  v_tmp(3) = zstore(iglob) - zelm(1)

  ! scalar product (dot product)
  tmp = v_tmp(1)*face_n(1) + v_tmp(2)*face_n(2) + v_tmp(3)*face_n(3)

  ! makes sure normal points outwards, that is points away from this additional corner
  ! and scalar product (dot product) is negative
  if (tmp > 0.d0) then
    face_n(:) = - face_n(:)
  endif

  ! determine orthonormal basis
  ! (assumes NGLLX == NGLLY == NGLLZ)
  do b = 1,NGLLY
    do a = 1,NGLLX
      xxi = 0.d0
      xeta = 0.d0
      yxi = 0.d0
      yeta = 0.d0
      zxi = 0.d0
      zeta = 0.d0

      do ia = 1,NGNOD2D
        xxi = xxi + dershape2D(1,ia,a,b) * xelm(ia)
        xeta = xeta + dershape2D(2,ia,a,b) * xelm(ia)
        yxi = yxi + dershape2D(1,ia,a,b) * yelm(ia)
        yeta = yeta + dershape2D(2,ia,a,b) * yelm(ia)
        zxi = zxi + dershape2D(1,ia,a,b) * zelm(ia)
        zeta = zeta + dershape2D(2,ia,a,b) * zelm(ia)
      enddo

      ! calculate the unnormalized normal to the boundary
      unx = yxi * zeta - yeta * zxi
      uny = zxi * xeta - zeta * xxi
      unz = xxi * yeta - xeta * yxi

      jacobian = dsqrt(unx*unx + uny*uny + unz*unz)
      if (jacobian <= 0.d0) stop '2D Jacobian undefined in routine determine_face_basis'

      ! normalize normal vector
      nx = unx / jacobian
      ny = uny / jacobian
      nz = unz / jacobian

      ! double-check with face normal to have normal pointing outwards
      ! determines orientation of normal and flips direction such that normal points outwards
      tmp = face_n(1) * nx + face_n(2) * ny + face_n(3) * nz
      if (tmp < 0.d0) then
        ! swap direction
        nx = - nx; ny = - ny; nz = - nz
      endif

      ! first tangent: normalize the xi-direction tangent
      t1_norm = dsqrt(xxi*xxi + yxi*yxi + zxi*zxi)
      if (t1_norm <= 0.d0) stop '2D vector xi undefined in routine determine_face_basis'

      t1x = xxi / t1_norm
      t1y = yxi / t1_norm
      t1z = zxi / t1_norm

      ! second tangent: n x t1 (orthogonal to both, unit length by construction)
      ! note: for distorted, bended elements, the vector in eta-direction might not be orthogonal to xi-direction vector.
      !       thus, the cross-product is used between normal and first tangential vector.
      t2x = ny * t1z - nz * t1y
      t2y = nz * t1x - nx * t1z
      t2z = nx * t1y - ny * t1x

      ! fill return
      nx_2D(a,b) = real(nx, kind=CUSTOM_REAL)
      ny_2D(a,b) = real(ny, kind=CUSTOM_REAL)
      nz_2D(a,b) = real(nz, kind=CUSTOM_REAL)

      t1x_2D(a,b) = real(t1x, kind=CUSTOM_REAL)
      t1y_2D(a,b) = real(t1y, kind=CUSTOM_REAL)
      t1z_2D(a,b) = real(t1z, kind=CUSTOM_REAL)

      t2x_2D(a,b) = real(t2x, kind=CUSTOM_REAL)
      t2y_2D(a,b) = real(t2y, kind=CUSTOM_REAL)
      t2z_2D(a,b) = real(t2z, kind=CUSTOM_REAL)
    enddo
  enddo

  end subroutine

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_tangential_derivs(F, dt1_F, dt2_F, &
                                       inv_a11, inv_a12, inv_a22, &
                                       dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)

  use specfem_par, only: hprime_xx

  implicit none
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY), intent(in) :: F
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY), intent(out) :: dt1_F, dt2_F
  ! metric and coordinate derivatives
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY), intent(in) :: inv_a11, inv_a12, inv_a22
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY), intent(in) :: dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2

  ! local parameters
  real(kind=CUSTOM_REAL) :: dF_ds1, dF_ds2, hp1, hp2
  integer :: a, b, l

  ! derivatives in tangential directions
  do b = 1, NGLLY
    do a = 1, NGLLX
      dF_ds1 = 0.0_CUSTOM_REAL
      dF_ds2 = 0.0_CUSTOM_REAL

      ! assumes NGLLX == NGLLY
      do l = 1, NGLLX
        hp1 = hprime_xx(a, l)
        dF_ds1 = dF_ds1 + F(l, b) * hp1
        hp2 = hprime_xx(b, l)
        dF_ds2 = dF_ds2 + F(a, l) * hp2
      enddo

      ! derivative of field F along tangential direction t1 (d/dt1)
      dt1_F(a, b) = (inv_a11(a, b) * dF_ds1 + inv_a12(a, b) * dF_ds2) * dot_g1_t1(a, b) &
                  + (inv_a12(a, b) * dF_ds1 + inv_a22(a, b) * dF_ds2) * dot_g2_t1(a, b)

      ! derivative of field F along tangential direction t2 (d/dt2)
      dt2_F(a, b) = (inv_a11(a, b) * dF_ds1 + inv_a12(a, b) * dF_ds2) * dot_g1_t2(a, b) &
                  + (inv_a12(a, b) * dF_ds1 + inv_a22(a, b) * dF_ds2) * dot_g2_t2(a, b)
    enddo
  enddo

  end subroutine compute_tangential_derivs

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_surface_laplacian_weak(F, laplace_F, jacobian_2D, inv_a11, inv_a12, inv_a22, wgll_2D)

  use specfem_par, only: hprime_xx

  implicit none
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY), intent(in) :: F
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY), intent(out) :: laplace_F
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY), intent(in) :: jacobian_2D
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY), intent(in) :: inv_a11, inv_a12, inv_a22
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY), intent(in) :: wgll_2D

  ! local parameters
  real(kind=CUSTOM_REAL) :: dF_ds1(NGLLX, NGLLY), dF_ds2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: flux1(NGLLX, NGLLY), flux2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: stiff_F, hp1, hp2
  integer :: a, b, l

  ! reference derivatives and contravariant metric fluxes
  do b = 1, NGLLY
    do a = 1, NGLLX
      dF_ds1(a,b) = 0.0_CUSTOM_REAL
      dF_ds2(a,b) = 0.0_CUSTOM_REAL

      ! assumes NGLLX == NGLLY
      do l = 1, NGLLX
        hp1 = hprime_xx(a, l)
        dF_ds1(a,b) = dF_ds1(a,b) + F(l, b) * hp1
        hp2 = hprime_xx(b, l)
        dF_ds2(a,b) = dF_ds2(a,b) + F(a, l) * hp2
      enddo

      flux1(a,b) = jacobian_2D(a,b) * (inv_a11(a,b) * dF_ds1(a,b) + inv_a12(a,b) * dF_ds2(a,b))
      flux2(a,b) = jacobian_2D(a,b) * (inv_a12(a,b) * dF_ds1(a,b) + inv_a22(a,b) * dF_ds2(a,b))
    enddo
  enddo

  ! weak surface Laplacian: M^-1 * (-K F)
  do b = 1, NGLLY
    do a = 1, NGLLX
      stiff_F = 0.0_CUSTOM_REAL

      do l = 1, NGLLX
        hp1 = hprime_xx(l, a)
        stiff_F = stiff_F + wgll_2D(l,b) * flux1(l,b) * hp1
        hp2 = hprime_xx(l, b)
        stiff_F = stiff_F + wgll_2D(a,l) * flux2(a,l) * hp2
      enddo

      laplace_F(a,b) = - stiff_F / (jacobian_2D(a,b) * wgll_2D(a,b))
    enddo
  enddo

  end subroutine compute_surface_laplacian_weak

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_hw_elastic(NSPEC_AB,NGLOB_AB,accel,ibool, &
                                abs_boundary_jacobian2Dw, &
                                abs_boundary_ijk,abs_boundary_ispec, &
                                num_abs_boundary_faces, &
                                displ,veloc,rho_vp,rho_vs, &
                                ispec_is_elastic, &
                                b_num_abs_boundary_faces,b_absorb_field)

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE,hprime_xx,xstore,ystore,zstore,rhostore,DT

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: displ,veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! adjoint simulations
  integer,intent(in) :: b_num_abs_boundary_faces
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  ! local parameters
  real(kind=CUSTOM_REAL) :: vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll

  ! 2D local arrays for surface geometry and displacements
  real(kind=CUSTOM_REAL) :: x_2D(NGLLX, NGLLY), y_2D(NGLLX, NGLLY), z_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: nx_2D(NGLLX, NGLLY), ny_2D(NGLLX, NGLLY), nz_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: t1x_2D(NGLLX, NGLLY), t1y_2D(NGLLX, NGLLY), t1z_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: t2x_2D(NGLLX, NGLLY), t2y_2D(NGLLX, NGLLY), t2z_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: u_n(NGLLX, NGLLY), u_t1(NGLLX, NGLLY), u_t2(NGLLX, NGLLY)  ! projected displacement
  real(kind=CUSTOM_REAL) :: phi1_n(NGLLX, NGLLY), phi1_t1(NGLLX, NGLLY), phi1_t2(NGLLX, NGLLY)  ! projected aux phi_1

  ! Derivative
  real(kind=CUSTOM_REAL) :: dx_ds1, dx_ds2
  real(kind=CUSTOM_REAL) :: dy_ds1, dy_ds2
  real(kind=CUSTOM_REAL) :: dz_ds1, dz_ds2

  ! Metric tensor components
  real(kind=CUSTOM_REAL) :: a11, a12, a22
  real(kind=CUSTOM_REAL) :: inv_a11(NGLLX, NGLLY), inv_a12(NGLLX, NGLLY), inv_a22(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dot_g1_t1(NGLLX, NGLLY), dot_g2_t1(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dot_g1_t2(NGLLX, NGLLY), dot_g2_t2(NGLLX, NGLLY)

  ! Tangential derivatives
  real(kind=CUSTOM_REAL) :: dt1_u_n(NGLLX, NGLLY), dt2_u_n(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dt1_u_t1(NGLLX, NGLLY), dt2_u_t1(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dt1_u_t2(NGLLX, NGLLY), dt2_u_t2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dt1_phi1_n(NGLLX, NGLLY), dt2_phi1_n(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dt1_phi1_t1(NGLLX, NGLLY), dt2_phi1_t1(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dt1_phi1_t2(NGLLX, NGLLY), dt2_phi1_t2(NGLLX, NGLLY)

  real(kind=CUSTOM_REAL) :: rhol, csl, cpl, fac, mu
  real(kind=CUSTOM_REAL) :: deltat, deltatsqover2
  real(kind=CUSTOM_REAL) :: t1x, t1y, t1z, t2x, t2y, t2z, t1_norm, hp1, hp2, det_a
  real(kind=CUSTOM_REAL) :: traction_n, traction_t1, traction_t2
  real(kind=CUSTOM_REAL) :: div_u_t, div_phi1_t

  logical :: mask_vary(3)
  integer :: id1, id2, a, b, l
  integer :: face_iglob(NGLLX, NGLLY)

  ! Newmark predictor factors for phi1
  deltat = real(DT, kind=CUSTOM_REAL)
  deltatsqover2 = 0.5_CUSTOM_REAL * deltat * deltat

  ! loop over all boundary faces
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    ! only for elastic domains
    if (.not. ispec_is_elastic(ispec)) cycle

    ! prepare spatial derivative arrays
    ! find the varying local coordinates for the face to map 1D index igll to 2D grid (a,b)
    mask_vary(:) = .false.
    do igll = 2, NGLLSQUARE
      if (abs_boundary_ijk(1, igll, iface) /= abs_boundary_ijk(1, 1, iface)) mask_vary(1) = .true.
      if (abs_boundary_ijk(2, igll, iface) /= abs_boundary_ijk(2, 1, iface)) mask_vary(2) = .true.
      if (abs_boundary_ijk(3, igll, iface) /= abs_boundary_ijk(3, 1, iface)) mask_vary(3) = .true.
    enddo

    id1 = 0; id2 = 0
    if (.not. mask_vary(1)) then
      id1 = 2; id2 = 3
    else if (.not. mask_vary(2)) then
      id1 = 1; id2 = 3
    else
      id1 = 1; id2 = 2
    endif

    ! gather coordinates and fields onto the 2D grid face
    do igll = 1, NGLLSQUARE
      i = abs_boundary_ijk(1,igll,iface)
      j = abs_boundary_ijk(2,igll,iface)
      k = abs_boundary_ijk(3,igll,iface)
      iglob = ibool(i,j,k,ispec)

      a = abs_boundary_ijk(id1,igll,iface)
      b = abs_boundary_ijk(id2,igll,iface)

      face_iglob(a,b) = iglob

      x_2D(a,b) = xstore(iglob)
      y_2D(a,b) = ystore(iglob)
      z_2D(a,b) = zstore(iglob)

      !nx_2D(a,b) = abs_boundary_normal(1,igll,iface)
      !ny_2D(a,b) = abs_boundary_normal(2,igll,iface)
      !nz_2D(a,b) = abs_boundary_normal(3,igll,iface)
    enddo

    ! orthonormal basis at each point
    call determine_face_ortho_basis(iface,ispec,nx_2D,ny_2D,nz_2D,t1x_2D,t1y_2D,t1z_2D,t2x_2D,t2y_2D,t2z_2D)

    ! projections at each point
    ! (assumes NGLLX == NGLLY == NGLLZ)
    do b = 1,NGLLY
      do a = 1,NGLLX
        iglob = face_iglob(a,b)

        ! pre-compute the orthonormal tangential basis
        if (.false.) then
          ! note: use the face mid-point to evaluate the orthonormal basis
          !       and assign it to all GLL nodes for this face to avoid problems if the boundary is curved
          !nx = nx_2D(MIDX,MIDY)
          !ny = ny_2D(MIDX,MIDY)
          !nz = nz_2D(MIDX,MIDY)
          ! normal
          nx = nx_2D(a,b)
          ny = ny_2D(a,b)
          nz = nz_2D(a,b)

          ! construct right-handed orthonormal tangential basis t1, t2
          ! choose a non-collinear vector to n to compute t1
          ! (Hughes-Moeller method)
          if (abs(nx) < 0.9_CUSTOM_REAL) then
            t1x = 0.0_CUSTOM_REAL
            t1y = -nz
            t1z = ny
          else
            t1x = -ny
            t1y = nx
            t1z = 0.0_CUSTOM_REAL
          endif
          t1_norm = sqrt(t1x*t1x + t1y*t1y + t1z*t1z)

          ! avoid divison by zero
          if (abs(t1_norm) < 1.d-24) t1_norm = 1._CUSTOM_REAL

          ! normalizes t1
          t1x = t1x / t1_norm; t1y = t1y / t1_norm; t1z = t1z / t1_norm

          ! t2 = n x t1
          t2x = ny * t1z - nz * t1y
          t2y = nz * t1x - nx * t1z
          t2z = nx * t1y - ny * t1x

          ! store basis for reuse in igll loop
          t1x_2D(a,b) = t1x; t1y_2D(a,b) = t1y; t1z_2D(a,b) = t1z
          t2x_2D(a,b) = t2x; t2y_2D(a,b) = t2y; t2z_2D(a,b) = t2z
        else
          ! gets pre-calculated basis
          ! normal
          nx = nx_2D(a,b)
          ny = ny_2D(a,b)
          nz = nz_2D(a,b)

          ! tangential t1 & t2
          t1x = t1x_2D(a,b)
          t1y = t1y_2D(a,b)
          t1z = t1z_2D(a,b)

          t2x = t2x_2D(a,b)
          t2y = t2y_2D(a,b)
          t2z = t2z_2D(a,b)
        endif

        ! project displacement
        u_n(a,b)  = displ(1,iglob) * nx  + displ(2,iglob) * ny  + displ(3,iglob) * nz
        u_t1(a,b) = displ(1,iglob) * t1x + displ(2,iglob) * t1y + displ(3,iglob) * t1z
        u_t2(a,b) = displ(1,iglob) * t2x + displ(2,iglob) * t2y + displ(3,iglob) * t2z

        ! auxiliary state for first-order state (P==1), predicted to current time t+deltat
        ! (Newmark predictor: phi_pred = phi + dt*phi_dot + dt^2/2 * phi_dotdot)
        phi1_n(a,b)  = hw_phi(1,a,b,iface,1) + deltat * hw_phi_dot(1,a,b,iface,1) + deltatsqover2 * hw_phi_dotdot(1,a,b,iface,1)
        phi1_t1(a,b) = hw_phi(2,a,b,iface,1) + deltat * hw_phi_dot(2,a,b,iface,1) + deltatsqover2 * hw_phi_dotdot(2,a,b,iface,1)
        phi1_t2(a,b) = hw_phi(3,a,b,iface,1) + deltat * hw_phi_dot(3,a,b,iface,1) + deltatsqover2 * hw_phi_dotdot(3,a,b,iface,1)
      enddo
    enddo

    ! compute reference derivatives using 1D GLL derivative matrix
    ! (assumes NGLLX == NGLLY == NGLLZ)
    do b = 1, NGLLY
      do a = 1, NGLLX
        dx_ds1 = 0.0_CUSTOM_REAL; dx_ds2 = 0.0_CUSTOM_REAL
        dy_ds1 = 0.0_CUSTOM_REAL; dy_ds2 = 0.0_CUSTOM_REAL
        dz_ds1 = 0.0_CUSTOM_REAL; dz_ds2 = 0.0_CUSTOM_REAL

        do l = 1, NGLLX
          hp1 = hprime_xx(a,l)
          dx_ds1 = dx_ds1 + x_2D(l,b) * hp1
          dy_ds1 = dy_ds1 + y_2D(l,b) * hp1
          dz_ds1 = dz_ds1 + z_2D(l,b) * hp1

          hp2 = hprime_xx(b,l)
          dx_ds2 = dx_ds2 + x_2D(a,l) * hp2
          dy_ds2 = dy_ds2 + y_2D(a,l) * hp2
          dz_ds2 = dz_ds2 + z_2D(a,l) * hp2
        enddo

        ! metric tensor components a_{alpha beta} = g_alpha . g_beta
        a11 = dx_ds1**2 + dy_ds1**2 + dz_ds1**2
        a12 = dx_ds1 * dx_ds2 + dy_ds1 * dy_ds2 + dz_ds1 * dz_ds2
        a22 = dx_ds2**2 + dy_ds2**2 + dz_ds2**2

        ! surface Jacobian (squared)
        det_a = a11 * a22 - a12**2

        ! avoid divison by zero
        if (abs(det_a) < 1.d-24) det_a = 1._CUSTOM_REAL

        ! inverse metric tensor
        inv_a11(a,b) = a22 / det_a
        inv_a12(a,b) = -a12 / det_a
        inv_a22(a,b) = a11 / det_a

        ! look up pre-computed basis
        t1x = t1x_2D(a,b); t1y = t1y_2D(a,b); t1z = t1z_2D(a,b)
        t2x = t2x_2D(a,b); t2y = t2y_2D(a,b); t2z = t2z_2D(a,b)

        ! dot-products
        dot_g1_t1(a,b) = dx_ds1 * t1x + dy_ds1 * t1y + dz_ds1 * t1z
        dot_g2_t1(a,b) = dx_ds2 * t1x + dy_ds2 * t1y + dz_ds2 * t1z

        dot_g1_t2(a,b) = dx_ds1 * t2x + dy_ds1 * t2y + dz_ds1 * t2z
        dot_g2_t2(a,b) = dx_ds2 * t2x + dy_ds2 * t2y + dz_ds2 * t2z
      enddo
    enddo

    ! tangential derivatives
    call compute_tangential_derivs(u_n, dt1_u_n, dt2_u_n, &
                                   inv_a11, inv_a12, inv_a22, &
                                   dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)
    call compute_tangential_derivs(u_t1, dt1_u_t1, dt2_u_t1, &
                                   inv_a11, inv_a12, inv_a22, &
                                   dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)
    call compute_tangential_derivs(u_t2, dt1_u_t2, dt2_u_t2, &
                                   inv_a11, inv_a12, inv_a22, &
                                   dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)

    call compute_tangential_derivs(phi1_n, dt1_phi1_n, dt2_phi1_n, &
                                   inv_a11, inv_a12, inv_a22, &
                                   dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)
    call compute_tangential_derivs(phi1_t1, dt1_phi1_t1, dt2_phi1_t1, &
                                   inv_a11, inv_a12, inv_a22, &
                                   dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)
    call compute_tangential_derivs(phi1_t2, dt1_phi1_t2, dt2_phi1_t2, &
                                   inv_a11, inv_a12, inv_a22, &
                                   dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)

    ! reference GLL points on boundary face
    do igll = 1,NGLLSQUARE
      ! gets local indices for GLL point
      i = abs_boundary_ijk(1,igll,iface)
      j = abs_boundary_ijk(2,igll,iface)
      k = abs_boundary_ijk(3,igll,iface)

      ! gets velocity
      iglob = ibool(i,j,k,ispec)

      vx = veloc(1,iglob)
      vy = veloc(2,iglob)
      vz = veloc(3,iglob)

      ! local face indexing
      a = abs_boundary_ijk(id1,igll,iface)
      b = abs_boundary_ijk(id2,igll,iface)

      ! looks up pre-computed basis
      ! normal
      nx = nx_2D(a,b)
      ny = ny_2D(a,b)
      nz = nz_2D(a,b)
      ! tangential vectors
      t1x = t1x_2D(a,b); t1y = t1y_2D(a,b); t1z = t1z_2D(a,b)
      t2x = t2x_2D(a,b); t2y = t2y_2D(a,b); t2z = t2z_2D(a,b)

      ! velocity component in normal direction (normal points out of element)
      vn = vx*nx + vy*ny + vz*nz

      ! first-order contribution
      ! (dash-pot equivalent to P1 from Stacey)
      ! velocity vector component * vp * rho in normal direction + vs * rho component tangential to it
      tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
      ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
      tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

      ! H-W higher-order contributions
      ! factors
      rhol = rhostore(i,j,k,ispec)
      csl = rho_vs(i,j,k,ispec) / rhol
      cpl = rho_vp(i,j,k,ispec) / rhol
      mu = rhol * csl**2

      fac = rho_vs(i,j,k,ispec) * (2.0_CUSTOM_REAL * csl - cpl)

      ! surface divergence
      div_u_t = dt1_u_t1(a,b) + dt2_u_t2(a,b)
      div_phi1_t = dt1_phi1_t1(a,b) + dt2_phi1_t2(a,b)

      ! note: as we subtract the total traction from accel, we switch signs of +/- terms.
      !       the term for normal traction `hw_a_coeff(1) * fac * div_u_t` and `hw_a_coeff(1) * fac * dt1_u_n(a,b)` for
      !       the shear tractions are equivalent to the P3 Stacey approximation.
      !       new are the auxiliary state term corrections `**_phi1_**` for the Hagstrom-Warburton boundary.

      ! w/out corner weights
      ! normal traction
      !   traction_n = hw_a_coeff(1) * fac * div_u_t - mu * div_phi1_t
      ! shear tractions
      !   traction_t1 = - hw_a_coeff(1) * fac * dt1_u_n(a,b) - mu * dt1_phi1_n(a,b)
      !   traction_t2 = - hw_a_coeff(1) * fac * dt2_u_n(a,b) - mu * dt2_phi1_n(a,b)
      !
      ! w/ corner weights
      !   the first correction (`hw_a_coeff(1) * fac * ...` terms) is applied at full weight.
      !   the H-W state phi1 correction (`mu * ...` terms) is weighted by hw_corner_weight so that
      !   corner/edge nodes shared by multiple faces do not accumulate the phi1 contribution
      !   multiple times (Kucherov & Givoli 2010, partition-of-unity at corner nodes).
      !
      ! correction (from physical displacement u, no corner weighting needed)
      traction_n  =   hw_a_coeff(1) * fac * div_u_t
      traction_t1 = - hw_a_coeff(1) * fac * dt1_u_n(a,b)
      traction_t2 = - hw_a_coeff(1) * fac * dt2_u_n(a,b)

      ! H-W auxiliary state correction
      if (HW_APPLY_EDGE_CORNER_COMPATIBILITY) then
        ! with  edge/corner compatibility conditions, no corner weighting needed
        traction_n  = traction_n  - mu * div_phi1_t
        traction_t1 = traction_t1 - mu * dt1_phi1_n(a,b)
        traction_t2 = traction_t2 - mu * dt2_phi1_n(a,b)
      else
        ! uses corner weighting
        traction_n  = traction_n  - mu * div_phi1_t * hw_corner_weight(a,b,iface)
        traction_t1 = traction_t1 - mu * dt1_phi1_n(a,b) * hw_corner_weight(a,b,iface)
        traction_t2 = traction_t2 - mu * dt2_phi1_n(a,b) * hw_corner_weight(a,b,iface)
      endif

      ! total contribution
      tx = tx + traction_n * nx + traction_t1 * t1x + traction_t2 * t2x
      ty = ty + traction_n * ny + traction_t1 * t1y + traction_t2 * t2y
      tz = tz + traction_n * nz + traction_t1 * t1z + traction_t2 * t2z

      ! gets associated, weighted jacobian
      jacobianw = abs_boundary_jacobian2Dw(igll,iface)

      ! adds boundary term (weak form)
      accel(1,iglob) = accel(1,iglob) - tx*jacobianw
      accel(2,iglob) = accel(2,iglob) - ty*jacobianw
      accel(3,iglob) = accel(3,iglob) - tz*jacobianw

      ! for kernel simulations
      if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
        b_absorb_field(1,igll,iface) = tx*jacobianw
        b_absorb_field(2,igll,iface) = ty*jacobianw
        b_absorb_field(3,igll,iface) = tz*jacobianw
      endif
    enddo
  enddo

  end subroutine compute_hw_elastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_hw_abc_states(displ, veloc, accel)

  use specfem_par, only: DT, NGLOB_AB, ibool, hprime_xx, xstore, ystore, zstore, rhostore, &
    wxgll, wygll, wzgll, &
    num_abs_boundary_faces, abs_boundary_ijk, abs_boundary_ispec ! abs_boundary_normal

  use specfem_par_elastic, only: ispec_is_elastic, rho_vp, rho_vs

  implicit none

  ! inputs
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: displ, veloc, accel

  ! local parameters
  real(kind=CUSTOM_REAL) :: nx,ny,nz,t1_norm,hp1,hp2,det_a
  integer :: ispec,iglob,i,j,k,iface,igll,p,prev_p,isubstep

  ! 2D local arrays for surface geometry and displacements
  real(kind=CUSTOM_REAL) :: x_2D(NGLLX, NGLLY), y_2D(NGLLX, NGLLY), z_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: nx_2D(NGLLX, NGLLY), ny_2D(NGLLX, NGLLY), nz_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: t1x_2D(NGLLX, NGLLY), t1y_2D(NGLLX, NGLLY), t1z_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: t2x_2D(NGLLX, NGLLY), t2y_2D(NGLLX, NGLLY), t2z_2D(NGLLX, NGLLY)

  ! indices 2D
  integer :: i_2D(NGLLX, NGLLY), j_2D(NGLLX, NGLLY), k_2D(NGLLX, NGLLY)

  ! projected variables
  real(kind=CUSTOM_REAL) :: u_n(NGLLX, NGLLY), u_t1(NGLLX, NGLLY), u_t2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: v_n(NGLLX, NGLLY), v_t1(NGLLX, NGLLY), v_t2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: a_n(NGLLX, NGLLY), a_t1(NGLLX, NGLLY), a_t2(NGLLX, NGLLY)

  ! predictor arrays
  real(kind=CUSTOM_REAL) :: hw_phi_pred(NDIM, NGLLX, NGLLY, HW_ORDER_P)
  real(kind=CUSTOM_REAL) :: hw_phi_dot_pred(NDIM, NGLLX, NGLLY, HW_ORDER_P)
  real(kind=CUSTOM_REAL) :: hw_phi_dotdot_new(NDIM, NGLLX, NGLLY, HW_ORDER_P)

  ! previous order variables (p-1)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: phi_prev_n,phi_prev_t1,phi_prev_t2
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: phi_dot_prev_n,phi_dot_prev_t1,phi_dot_prev_t2
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: phi_dotdot_prev_n,phi_dotdot_prev_t1,phi_dotdot_prev_t2

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,HW_ORDER_P) :: phi_ini_n, phi_ini_t1, phi_ini_t2
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,HW_ORDER_P) :: phi_ini_n_dot, phi_ini_t1_dot, phi_ini_t2_dot
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,HW_ORDER_P) :: phi_ini_n_ddot, phi_ini_t1_ddot, phi_ini_t2_ddot

  ! derivative arrays
  real(kind=CUSTOM_REAL) :: dx_ds1, dx_ds2
  real(kind=CUSTOM_REAL) :: dy_ds1, dy_ds2
  real(kind=CUSTOM_REAL) :: dz_ds1, dz_ds2

  ! Metric tensor components
  real(kind=CUSTOM_REAL) :: a11(NGLLX, NGLLY), a12(NGLLX, NGLLY), a22(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: inv_a11(NGLLX, NGLLY), inv_a12(NGLLX, NGLLY), inv_a22(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: jacobian_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: wgll_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dot_g1_t1(NGLLX, NGLLY), dot_g2_t1(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dot_g1_t2(NGLLX, NGLLY), dot_g2_t2(NGLLX, NGLLY)

  ! Tangential derivatives for source terms and Laplacians
  real(kind=CUSTOM_REAL) :: dt1_phi_dot_prev_t1(NGLLX, NGLLY), dt2_phi_dot_prev_t1(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dt1_phi_dot_prev_t2(NGLLX, NGLLY), dt2_phi_dot_prev_t2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dt1_phi_dot_prev_n(NGLLX, NGLLY), dt2_phi_dot_prev_n(NGLLX, NGLLY)

  real(kind=CUSTOM_REAL) :: dt1_phi_t1(NGLLX, NGLLY), dt2_phi_t1(NGLLX, NGLLY)

  real(kind=CUSTOM_REAL) :: dt1_phi_t2(NGLLX, NGLLY), dt2_phi_t2(NGLLX, NGLLY)

  real(kind=CUSTOM_REAL) :: div_phi_t(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dt1_div_phi_t(NGLLX, NGLLY), dt2_div_phi_t(NGLLX, NGLLY)

  real(kind=CUSTOM_REAL) :: div_phi_dot_prev_t(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: laplace_phi_n(NGLLX, NGLLY), laplace_phi_t1(NGLLX, NGLLY), laplace_phi_t2(NGLLX, NGLLY)

  real(kind=CUSTOM_REAL) :: S_n(NGLLX, NGLLY), S_t1(NGLLX, NGLLY), S_t2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: spatial_term_n(NGLLX, NGLLY), spatial_term_t1(NGLLX, NGLLY), spatial_term_t2(NGLLX, NGLLY)

  real(kind=CUSTOM_REAL) :: deltat, deltat_hw, deltatover2_hw, deltatsqover2_hw
  real(kind=CUSTOM_REAL) :: rhol, csl, cpl, coeff_S_n, coeff_S_t, damp_p
  real(kind=CUSTOM_REAL) :: t1x, t1y, t1z, t2x, t2y, t2z

  logical :: mask_vary(3)
  integer :: id1, id2, a, b, l
  integer :: face_iglob(NGLLX, NGLLY)

  ! corner/edge
  integer :: icorner, iedge, i_edge_local, id_perp_edge
  integer :: i_edge_t1_local, i_edge_t2_local
  integer :: iglob_edge, iglob_corner, mid_a, mid_b
  real(kind=CUSTOM_REAL) :: jacobian_perp_edge
  real(kind=CUSTOM_REAL) :: cpl_edge, csl_edge
  real(kind=CUSTOM_REAL) :: jacobian_t1_corner, jacobian_t2_corner
  real(kind=CUSTOM_REAL) :: cpl_corner, csl_corner

  ! uses Newmark time-scheme to update auxiliary state variables
  deltat = real(DT, kind=CUSTOM_REAL)

  if (HW_ABC_SUBSTEPS < 1) stop 'Error HW_ABC_SUBSTEPS must be >= 1'

  deltat_hw = deltat / real(HW_ABC_SUBSTEPS, kind=CUSTOM_REAL)
  deltatover2_hw = 0.5_CUSTOM_REAL * deltat_hw
  deltatsqover2_hw = 0.5_CUSTOM_REAL * deltat_hw * deltat_hw

  ! loops over all faces
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    ! only for elastic domains
    if (.not. ispec_is_elastic(ispec)) cycle

    ! prepare spatial derivative arrays
    ! find the varying local coordinates for the face to map 1D index igll to 2D grid (a,b)
    mask_vary(:) = .false.
    do igll = 2, NGLLSQUARE
      if (abs_boundary_ijk(1, igll, iface) /= abs_boundary_ijk(1, 1, iface)) mask_vary(1) = .true.
      if (abs_boundary_ijk(2, igll, iface) /= abs_boundary_ijk(2, 1, iface)) mask_vary(2) = .true.
      if (abs_boundary_ijk(3, igll, iface) /= abs_boundary_ijk(3, 1, iface)) mask_vary(3) = .true.
    enddo

    id1 = 0; id2 = 0
    if (.not. mask_vary(1)) then
      id1 = 2; id2 = 3
    else if (.not. mask_vary(2)) then
      id1 = 1; id2 = 3
    else
      id1 = 1; id2 = 2
    endif

    ! local 2D GLL quadrature weights for this face orientation
    do b = 1, NGLLY
      do a = 1, NGLLX
        if (id1 == 2 .and. id2 == 3) then
          wgll_2D(a,b) = wygll(a) * wzgll(b)
        else if (id1 == 1 .and. id2 == 3) then
          wgll_2D(a,b) = wxgll(a) * wzgll(b)
        else
          wgll_2D(a,b) = wxgll(a) * wygll(b)
        endif
      enddo
    enddo

    ! gather coordinates and fields onto the 2D grid face
    do igll = 1, NGLLSQUARE
      i = abs_boundary_ijk(1,igll,iface)
      j = abs_boundary_ijk(2,igll,iface)
      k = abs_boundary_ijk(3,igll,iface)
      iglob = ibool(i,j,k,ispec)

      a = abs_boundary_ijk(id1,igll,iface)
      b = abs_boundary_ijk(id2,igll,iface)

      face_iglob(a,b) = iglob

      x_2D(a,b) = xstore(iglob)
      y_2D(a,b) = ystore(iglob)
      z_2D(a,b) = zstore(iglob)

      !nx_2D(a,b) = abs_boundary_normal(1,igll,iface)
      !ny_2D(a,b) = abs_boundary_normal(2,igll,iface)
      !nz_2D(a,b) = abs_boundary_normal(3,igll,iface)

      i_2D(a,b) = i
      j_2D(a,b) = j
      k_2D(a,b) = k
    enddo

    ! orthonormal basis at each point
    call determine_face_ortho_basis(iface,ispec,nx_2D,ny_2D,nz_2D,t1x_2D,t1y_2D,t1z_2D,t2x_2D,t2y_2D,t2z_2D)

    ! projections at each point
    ! (assumes NGLLX == NGLLY == NGLLZ)
    do b = 1,NGLLY
      do a = 1,NGLLX
        iglob = face_iglob(a,b)

        ! pre-compute the orthonormal tangential basis
        if (.false.) then
          ! note: use the face mid-point to evaluate the orthonormal basis
          !       and assign it to all GLL nodes for this face to avoid problems if the boundary is curved
          !nx = nx_2D(MIDX,MIDY)
          !ny = ny_2D(MIDX,MIDY)
          !nz = nz_2D(MIDX,MIDY)
          ! normal
          nx = nx_2D(a,b)
          ny = ny_2D(a,b)
          nz = nz_2D(a,b)

          ! construct right-handed orthonormal tangential basis t1, t2
          ! choose a non-collinear vector to n to compute t1
          ! (Hughes-Moeller method)
          if (abs(nx) < 0.9_CUSTOM_REAL) then
            t1x = 0.0_CUSTOM_REAL
            t1y = -nz
            t1z = ny
          else
            t1x = -ny
            t1y = nx
            t1z = 0.0_CUSTOM_REAL
          endif
          t1_norm = sqrt(t1x*t1x + t1y*t1y + t1z*t1z)

          ! avoid divison by zero
          if (abs(t1_norm) < 1.d-24) t1_norm = 1._CUSTOM_REAL

          ! normalizes t1
          t1x = t1x / t1_norm; t1y = t1y / t1_norm; t1z = t1z / t1_norm

          ! t2 = n x t1
          t2x = ny * t1z - nz * t1y
          t2y = nz * t1x - nx * t1z
          t2z = nx * t1y - ny * t1x

          ! store basis for reuse in igll loop
          t1x_2D(a,b) = t1x; t1y_2D(a,b) = t1y; t1z_2D(a,b) = t1z
          t2x_2D(a,b) = t2x; t2y_2D(a,b) = t2y; t2z_2D(a,b) = t2z
        else
          ! gets pre-calculated basis
          ! normal
          nx = nx_2D(a,b)
          ny = ny_2D(a,b)
          nz = nz_2D(a,b)

          ! tangential t1 & t2
          t1x = t1x_2D(a,b)
          t1y = t1y_2D(a,b)
          t1z = t1z_2D(a,b)

          t2x = t2x_2D(a,b)
          t2y = t2y_2D(a,b)
          t2z = t2z_2D(a,b)
        endif

        ! project displacement/velocity/acceleration
        u_n(a,b)  = displ(1,iglob) * nx  + displ(2,iglob) * ny  + displ(3,iglob) * nz
        u_t1(a,b) = displ(1,iglob) * t1x + displ(2,iglob) * t1y + displ(3,iglob) * t1z
        u_t2(a,b) = displ(1,iglob) * t2x + displ(2,iglob) * t2y + displ(3,iglob) * t2z

        v_n(a,b)  = veloc(1,iglob) * nx  + veloc(2,iglob) * ny  + veloc(3,iglob) * nz
        v_t1(a,b) = veloc(1,iglob) * t1x + veloc(2,iglob) * t1y + veloc(3,iglob) * t1z
        v_t2(a,b) = veloc(1,iglob) * t2x + veloc(2,iglob) * t2y + veloc(3,iglob) * t2z

        a_n(a,b)  = accel(1,iglob) * nx  + accel(2,iglob) * ny  + accel(3,iglob) * nz
        a_t1(a,b) = accel(1,iglob) * t1x + accel(2,iglob) * t1y + accel(3,iglob) * t1z
        a_t2(a,b) = accel(1,iglob) * t2x + accel(2,iglob) * t2y + accel(3,iglob) * t2z
      enddo
    enddo

    ! compute reference derivatives using 1D GLL derivative matrix
    ! (assumes NGLLX == NGLLY == NGLLZ)
    do b = 1, NGLLY
      do a = 1, NGLLX
        dx_ds1 = 0.0_CUSTOM_REAL; dx_ds2 = 0.0_CUSTOM_REAL
        dy_ds1 = 0.0_CUSTOM_REAL; dy_ds2 = 0.0_CUSTOM_REAL
        dz_ds1 = 0.0_CUSTOM_REAL; dz_ds2 = 0.0_CUSTOM_REAL

        do l = 1, NGLLX
          hp1 = hprime_xx(a,l)
          dx_ds1 = dx_ds1 + x_2D(l,b) * hp1
          dy_ds1 = dy_ds1 + y_2D(l,b) * hp1
          dz_ds1 = dz_ds1 + z_2D(l,b) * hp1

          hp2 = hprime_xx(b,l)
          dx_ds2 = dx_ds2 + x_2D(a,l) * hp2
          dy_ds2 = dy_ds2 + y_2D(a,l) * hp2
          dz_ds2 = dz_ds2 + z_2D(a,l) * hp2
        enddo

        ! metric tensor components a_{alpha beta} = g_alpha . g_beta
        a11(a,b) = dx_ds1**2 + dy_ds1**2 + dz_ds1**2
        a12(a,b) = dx_ds1 * dx_ds2 + dy_ds1 * dy_ds2 + dz_ds1 * dz_ds2
        a22(a,b) = dx_ds2**2 + dy_ds2**2 + dz_ds2**2

        ! surface Jacobian (squared)
        det_a = a11(a,b) * a22(a,b) - a12(a,b)**2

        ! avoid divison by zero
        if (abs(det_a) < 1.d-24) det_a = 1._CUSTOM_REAL

        ! inverse metric tensor
        inv_a11(a,b) = a22(a,b) / det_a
        inv_a12(a,b) = -a12(a,b) / det_a
        inv_a22(a,b) = a11(a,b) / det_a
        jacobian_2D(a,b) = sqrt(det_a)

        ! look up pre-computed basis and normal
        nx = nx_2D(a,b)
        ny = ny_2D(a,b)
        nz = nz_2D(a,b)
        t1x = t1x_2D(a,b); t1y = t1y_2D(a,b); t1z = t1z_2D(a,b)
        t2x = t2x_2D(a,b); t2y = t2y_2D(a,b); t2z = t2z_2D(a,b)

        ! dot-products
        dot_g1_t1(a,b) = dx_ds1 * t1x + dy_ds1 * t1y + dz_ds1 * t1z
        dot_g2_t1(a,b) = dx_ds2 * t1x + dy_ds2 * t1y + dz_ds2 * t1z

        dot_g1_t2(a,b) = dx_ds1 * t2x + dy_ds1 * t2y + dz_ds1 * t2z
        dot_g2_t2(a,b) = dx_ds2 * t2x + dy_ds2 * t2y + dz_ds2 * t2z
      enddo
    enddo

    do isubstep = 1, HW_ABC_SUBSTEPS

      if (HW_APPLY_EDGE_CORNER_COMPATIBILITY) then
        ! store initial states
        do p = 1, HW_ORDER_P
          phi_ini_n(:,:,p) = hw_phi(1,:,:,iface,p)
          phi_ini_t1(:,:,p) = hw_phi(2,:,:,iface,p)
          phi_ini_t2(:,:,p) = hw_phi(3,:,:,iface,p)

          phi_ini_n_dot(:,:,p) = hw_phi_dot(1,:,:,iface,p)
          phi_ini_t1_dot(:,:,p) = hw_phi_dot(2,:,:,iface,p)
          phi_ini_t2_dot(:,:,p) = hw_phi_dot(3,:,:,iface,p)

          phi_ini_n_ddot(:,:,p) = hw_phi_dotdot(1,:,:,iface,p)
          phi_ini_t1_ddot(:,:,p) = hw_phi_dotdot(2,:,:,iface,p)
          phi_ini_t2_ddot(:,:,p) = hw_phi_dotdot(3,:,:,iface,p)
        enddo
      endif

      ! Newmark time-scheme
      ! (Predictor): Advance all variables for p = 1 ... HW_ORDER_P
      do p = 1, HW_ORDER_P
        do b = 1, NGLLY
          do a = 1, NGLLX
            ! displacement
            hw_phi_pred(1, a, b, p) = hw_phi(1, a, b, iface, p) &
                                      + deltat_hw * hw_phi_dot(1, a, b, iface, p) &
                                      + deltatsqover2_hw * hw_phi_dotdot(1, a, b, iface, p)
            hw_phi_pred(2, a, b, p) = hw_phi(2, a, b, iface, p) &
                                      + deltat_hw * hw_phi_dot(2, a, b, iface, p) &
                                      + deltatsqover2_hw * hw_phi_dotdot(2, a, b, iface, p)
            hw_phi_pred(3, a, b, p) = hw_phi(3, a, b, iface, p) &
                                      + deltat_hw * hw_phi_dot(3, a, b, iface, p) &
                                      + deltatsqover2_hw * hw_phi_dotdot(3, a, b, iface, p)

            ! velocity
            hw_phi_dot_pred(1, a, b, p) = hw_phi_dot(1, a, b, iface, p) + deltatover2_hw * hw_phi_dotdot(1, a, b, iface, p)
            hw_phi_dot_pred(2, a, b, p) = hw_phi_dot(2, a, b, iface, p) + deltatover2_hw * hw_phi_dotdot(2, a, b, iface, p)
            hw_phi_dot_pred(3, a, b, p) = hw_phi_dot(3, a, b, iface, p) + deltatover2_hw * hw_phi_dotdot(3, a, b, iface, p)

            ! acceleration
            hw_phi_dotdot_new(1, a, b, p) = 0.0_CUSTOM_REAL
            hw_phi_dotdot_new(2, a, b, p) = 0.0_CUSTOM_REAL
            hw_phi_dotdot_new(3, a, b, p) = 0.0_CUSTOM_REAL
          enddo
        enddo
      enddo   ! HW_ORDER_P

      ! auxiliary state updates (Recursive Loop)
      do p = 1, HW_ORDER_P
        if (p == 1) then
          ! first-order state (P==1)
          do b = 1, NGLLY
            do a = 1, NGLLX
              ! displacement
              phi_prev_n(a, b) = u_n(a, b)
              phi_prev_t1(a, b) = u_t1(a, b)
              phi_prev_t2(a, b) = u_t2(a, b)

              ! velocity
              phi_dot_prev_n(a, b) = v_n(a, b)
              phi_dot_prev_t1(a, b) = v_t1(a, b)
              phi_dot_prev_t2(a, b) = v_t2(a, b)

              ! acceleration
              phi_dotdot_prev_n(a, b) = a_n(a, b)
              phi_dotdot_prev_t1(a, b) = a_t1(a, b)
              phi_dotdot_prev_t2(a, b) = a_t2(a, b)
            enddo
          enddo
        else
          prev_p = p - 1  ! previous order
          do b = 1, NGLLY
            do a = 1, NGLLX
              phi_prev_n(a, b) = hw_phi_pred(1, a, b, prev_p)
              phi_prev_t1(a, b) = hw_phi_pred(2, a, b, prev_p)
              phi_prev_t2(a, b) = hw_phi_pred(3, a, b, prev_p)

              phi_dot_prev_n(a, b) = hw_phi_dot_pred(1, a, b, prev_p) + deltatover2_hw * hw_phi_dotdot_new(1, a, b, prev_p)
              phi_dot_prev_t1(a, b) = hw_phi_dot_pred(2, a, b, prev_p) + deltatover2_hw * hw_phi_dotdot_new(2, a, b, prev_p)
              phi_dot_prev_t2(a, b) = hw_phi_dot_pred(3, a, b, prev_p) + deltatover2_hw * hw_phi_dotdot_new(3, a, b, prev_p)

              phi_dotdot_prev_n(a, b) = hw_phi_dotdot_new(1, a, b, prev_p)
              phi_dotdot_prev_t1(a, b) = hw_phi_dotdot_new(2, a, b, prev_p)
              phi_dotdot_prev_t2(a, b) = hw_phi_dotdot_new(3, a, b, prev_p)
            enddo
          enddo
        endif

        ! spatial derivatives of previous order velocities (for coupling/source terms)
        call compute_tangential_derivs(phi_dot_prev_n, dt1_phi_dot_prev_n, dt2_phi_dot_prev_n, &
                                       inv_a11, inv_a12, inv_a22, &
                                       dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)
        call compute_tangential_derivs(phi_dot_prev_t1, dt1_phi_dot_prev_t1, dt2_phi_dot_prev_t1, &
                                       inv_a11, inv_a12, inv_a22, &
                                       dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)
        call compute_tangential_derivs(phi_dot_prev_t2, dt1_phi_dot_prev_t2, dt2_phi_dot_prev_t2, &
                                       inv_a11, inv_a12, inv_a22, &
                                       dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)

        ! spatial derivatives of current order predicted variables (for grad(div))
        call compute_tangential_derivs(hw_phi_pred(2, :, :, p), dt1_phi_t1, dt2_phi_t1, &
                                       inv_a11, inv_a12, inv_a22, &
                                       dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)

        call compute_tangential_derivs(hw_phi_pred(3, :, :, p), dt1_phi_t2, dt2_phi_t2, &
                                       inv_a11, inv_a12, inv_a22, &
                                       dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)

        call compute_surface_laplacian_weak(hw_phi_pred(1, :, :, p), laplace_phi_n, &
                                            jacobian_2D, inv_a11, inv_a12, inv_a22, wgll_2D)
        call compute_surface_laplacian_weak(hw_phi_pred(2, :, :, p), laplace_phi_t1, &
                                            jacobian_2D, inv_a11, inv_a12, inv_a22, wgll_2D)
        call compute_surface_laplacian_weak(hw_phi_pred(3, :, :, p), laplace_phi_t2, &
                                            jacobian_2D, inv_a11, inv_a12, inv_a22, wgll_2D)

        do b = 1, NGLLY
          do a = 1, NGLLX
            i = i_2D(a,b)
            j = j_2D(a,b)
            k = k_2D(a,b)

            rhol = rhostore(i,j,k,ispec)
            csl = rho_vs(i,j,k,ispec) / rhol        ! Vs == sqrt( mu/rho )
            cpl = rho_vp(i,j,k,ispec) / rhol        ! Vp == sqrt( (lambda + 2 mu)/rho )

            coeff_S_n = (cpl**2 - csl**2) / cpl     ! coeff_S_n == (lambda + mu)/rho / vp
            coeff_S_t = (cpl**2 - csl**2) / csl     ! coeff_S_t == (lambda + mu)/rho / vs

            div_phi_dot_prev_t(a,b) = dt1_phi_dot_prev_t1(a,b) + dt2_phi_dot_prev_t2(a,b)

            ! previous-order damping coefficients
            if (p == 1) then
              ! For the source at p=1, the previous-order damping uses the same pole as p=1 by convention
              damp_p  = hw_sigma(1) * hw_omega_c            ! units: 1/s
            else
              damp_p  = hw_sigma(p-1) * hw_omega_c          ! units: 1/s
            endif

            S_n(a,b) = hw_a_coeff(p) * phi_dotdot_prev_n(a,b) - damp_p * phi_dot_prev_n(a,b) &
                       + coeff_S_n * div_phi_dot_prev_t(a,b)

            S_t1(a,b) = hw_a_coeff(p) * phi_dotdot_prev_t1(a,b) - damp_p * phi_dot_prev_t1(a,b) &
                        + coeff_S_t * dt1_phi_dot_prev_n(a,b)

            S_t2(a,b) = hw_a_coeff(p) * phi_dotdot_prev_t2(a,b) - damp_p * phi_dot_prev_t2(a,b) &
                        + coeff_S_t * dt2_phi_dot_prev_n(a,b)
          enddo
        enddo

        ! gradient of divergence for t
        div_phi_t(:,:) = dt1_phi_t1(:,:) + dt2_phi_t2(:,:)

        call compute_tangential_derivs(div_phi_t, dt1_div_phi_t, dt2_div_phi_t, &
                                       inv_a11, inv_a12, inv_a22, &
                                       dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2)

        do b = 1, NGLLY
          do a = 1, NGLLX
            i = i_2D(a,b)
            j = j_2D(a,b)
            k = k_2D(a,b)

            rhol = rhostore(i,j,k,ispec)
            csl = rho_vs(i,j,k,ispec) / rhol
            cpl = rho_vp(i,j,k,ispec) / rhol

            spatial_term_n(a,b) = csl**2 * laplace_phi_n(a,b)
            spatial_term_t1(a,b) = csl**2 * laplace_phi_t1(a,b) + (cpl**2 - csl**2) * dt1_div_phi_t(a,b)
            spatial_term_t2(a,b) = csl**2 * laplace_phi_t2(a,b) + (cpl**2 - csl**2) * dt2_div_phi_t(a,b)

            ! damping coefficient
            damp_p = hw_sigma(p) * hw_omega_c   ! same denominator for normal and tangential when sigma_p is scalar

            ! semi-implicit ADE inversion
            ! normal direction
            hw_phi_dotdot_new(1, a, b, p) = (S_n(a,b) - damp_p * hw_phi_dot_pred(1, a, b, p) + spatial_term_n(a,b)) &
                                            / (hw_a_coeff(p) + deltatover2_hw * damp_p)
            ! tangential directions
            hw_phi_dotdot_new(2, a, b, p) = (S_t1(a,b) - damp_p * hw_phi_dot_pred(2, a, b, p) + spatial_term_t1(a,b)) &
                                            / (hw_a_coeff(p) + deltatover2_hw * damp_p)
            hw_phi_dotdot_new(3, a, b, p) = (S_t2(a,b) - damp_p * hw_phi_dot_pred(3, a, b, p) + spatial_term_t2(a,b)) &
                                            / (hw_a_coeff(p) + deltatover2_hw * damp_p)
          enddo
        enddo
      enddo  ! HW_ORDER_P

      ! Newmark (Corrector)
      do p = 1, HW_ORDER_P
        do b = 1, NGLLY
          do a = 1, NGLLX
            hw_phi(1, a, b, iface, p) = hw_phi_pred(1, a, b, p)
            hw_phi(2, a, b, iface, p) = hw_phi_pred(2, a, b, p)
            hw_phi(3, a, b, iface, p) = hw_phi_pred(3, a, b, p)

            hw_phi_dot(1, a, b, iface, p) = hw_phi_dot_pred(1, a, b, p) + deltatover2_hw * hw_phi_dotdot_new(1, a, b, p)
            hw_phi_dot(2, a, b, iface, p) = hw_phi_dot_pred(2, a, b, p) + deltatover2_hw * hw_phi_dotdot_new(2, a, b, p)
            hw_phi_dot(3, a, b, iface, p) = hw_phi_dot_pred(3, a, b, p) + deltatover2_hw * hw_phi_dotdot_new(3, a, b, p)

            hw_phi_dotdot(1, a, b, iface, p) = hw_phi_dotdot_new(1, a, b, p)
            hw_phi_dotdot(2, a, b, iface, p) = hw_phi_dotdot_new(2, a, b, p)
            hw_phi_dotdot(3, a, b, iface, p) = hw_phi_dotdot_new(3, a, b, p)
          enddo
        enddo
      enddo  ! HW_ORDER_P

      ! === EDGE COMPATIBILITY CONDITIONS (Fix C) ===
      if (HW_APPLY_EDGE_CORNER_COMPATIBILITY) then
        ! Apply after face-interior update, before next substep
        do p = 1, HW_ORDER_P
          do iedge = 1, 4
            if (face_edges(iedge, iface)%adj_face > 0) then
              i_edge_local = face_edges(iedge, iface)%i_edge
              id_perp_edge = face_edges(iedge, iface)%id_perp

              ! Get cpl, csl at the edge using the midpoint of the edge
              ! For edge at constant a (id_perp=1): use midpoint in b direction
              ! For edge at constant b (id_perp=2): use midpoint in a direction
              mid_a = MIDX
              mid_b = MIDY
              if (id_perp_edge == 1) then
                ! Edge at constant a (i_edge is a-index), perp is a direction (id1)
                ! Use midpoint along edge for material properties
                i = i_2D(i_edge_local, mid_b)
                j = j_2D(i_edge_local, mid_b)
                k = k_2D(i_edge_local, mid_b)
                iglob = face_iglob(i_edge_local, mid_b)  ! approximate
              else
                ! Edge at constant b (i_edge is b-index), perp is b direction (id2)
                i = i_2D(mid_a, i_edge_local)
                j = j_2D(mid_a, i_edge_local)
                k = k_2D(mid_a, i_edge_local)
                iglob = face_iglob(mid_a,i_edge_local)  ! approximate
              endif
              ! Get actual global index
              iglob_edge = ibool(i, j, k, ispec)
              rhol = rhostore(i, j, k, ispec)
              csl_edge = rho_vs(i, j, k, ispec) / rhol
              cpl_edge = rho_vp(i, j, k, ispec) / rhol

              ! jacobian_perp: sqrt of metric component in perpendicular direction
              if (id_perp_edge == 1) then
                ! Perpendicular is a direction (id1) -> sqrt(a11)
                jacobian_perp_edge = sqrt(a11(i_edge_local, mid_b))
              else
                ! Perpendicular is b direction (id2) -> sqrt(a22)
                jacobian_perp_edge = sqrt(a22(mid_a, i_edge_local))
              endif

              call hw_apply_edge_compatibility_newmark(hw_phi(1,:,:,iface,p), hw_phi(2,:,:,iface,p), hw_phi(3,:,:,iface,p), &
                                    hw_phi_dot(1,:,:,iface,p), hw_phi_dot(2,:,:,iface,p), hw_phi_dot(3,:,:,iface,p), &
                                    hw_phi_dotdot(1,:,:,iface,p), hw_phi_dotdot(2,:,:,iface,p), hw_phi_dotdot(3,:,:,iface,p), &
                                    phi_ini_n(:,:,p), phi_ini_t1(:,:,p), phi_ini_t2(:,:,p), &
                                    phi_ini_n_dot(:,:,p), phi_ini_t1_dot(:,:,p), phi_ini_t2_dot(:,:,p), &
                                    phi_ini_n_ddot(:,:,p), phi_ini_t1_ddot(:,:,p), phi_ini_t2_ddot(:,:,p), &
                                    cpl_edge, csl_edge, hw_a_coeff(p), deltat_hw, &
                                    hprime_xx, jacobian_perp_edge, &
                                    NGLLX, i_edge_local, face_edges(iedge, iface)%sigma, id_perp_edge, &
                                    face_edges(iedge, iface)%comp_vp)
            endif
          enddo
        enddo
      endif

      ! === CORNER COMPATIBILITY CONDITIONS (Fix D) ===
      if (HW_APPLY_EDGE_CORNER_COMPATIBILITY) then
        ! Apply after edge conditions
        do p = 1, HW_ORDER_P
          do icorner = 1, 4
            if (face_corners(icorner, iface)%adj_face_t1 > 0 .and. &
                face_corners(icorner, iface)%adj_face_t2 > 0) then
              i_edge_t1_local = face_corners(icorner, iface)%i_edge_t1
              i_edge_t2_local = face_corners(icorner, iface)%i_edge_t2

              ! Get cpl, csl at the corner
              i = i_2D(i_edge_t1_local, i_edge_t2_local)
              j = j_2D(i_edge_t1_local, i_edge_t2_local)
              k = k_2D(i_edge_t1_local, i_edge_t2_local)

              iglob_corner = ibool(i, j, k, ispec)
              rhol = rhostore(i, j, k, ispec)
              csl_corner = rho_vs(i, j, k, ispec) / rhol
              cpl_corner = rho_vp(i, j, k, ispec) / rhol

              ! jacobian_t1 = sqrt(a11) at corner (derivative in t1/id1/a direction)
              jacobian_t1_corner = sqrt(a11(i_edge_t1_local, i_edge_t2_local))
              ! jacobian_t2 = sqrt(a22) at corner (derivative in t2/id2/b direction)
              jacobian_t2_corner = sqrt(a22(i_edge_t1_local, i_edge_t2_local))

              call hw_apply_corner_compatibility_newmark(hw_phi(1,:,:,iface,p), hw_phi(2,:,:,iface,p), hw_phi(3,:,:,iface,p), &
                                      hw_phi_dot(1,:,:,iface,p), hw_phi_dot(2,:,:,iface,p), hw_phi_dot(3,:,:,iface,p), &
                                      hw_phi_dotdot(1,:,:,iface,p), hw_phi_dotdot(2,:,:,iface,p), hw_phi_dotdot(3,:,:,iface,p), &
                                      phi_ini_n(:,:,p), phi_ini_t1(:,:,p), phi_ini_t2(:,:,p), &
                                      phi_ini_n_dot(:,:,p), phi_ini_t1_dot(:,:,p), phi_ini_t2_dot(:,:,p), &
                                      phi_ini_n_ddot(:,:,p), phi_ini_t1_ddot(:,:,p), phi_ini_t2_ddot(:,:,p), &
                                      cpl_corner, csl_corner, hw_a_coeff(p), deltat_hw, &
                                      hprime_xx, hprime_xx, jacobian_t1_corner, jacobian_t2_corner, &
                                      NGLLX, i_edge_t1_local, i_edge_t2_local, &
                                      face_corners(icorner, iface)%sigma_t1, face_corners(icorner, iface)%sigma_t2, &
                                      face_corners(icorner, iface)%comp_vp_step1, face_corners(icorner, iface)%comp_vp_step2)
            endif
          enddo
        enddo
      endif

    enddo  ! HW_ABC_SUBSTEPS
  enddo

  end subroutine update_hw_abc_states

end module stacey_par


!
!-------------------------------------------------------------------------------------------------
!

! for elastic solver

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_forward(NSPEC_AB,NGLOB_AB,accel, &
                                                 ibool,iphase, &
                                                 abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                 abs_boundary_ijk,abs_boundary_ispec, &
                                                 num_abs_boundary_faces,veloc,rho_vp,rho_vs, &
                                                 ispec_is_elastic, &
                                                 it,b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLSQUARE,NDIM,IOABS,SAVE_RUN_BOUN_FOR_KH_INTEGRAL

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE

  ! Kirchoff-Helmholtz integrals
  use specfem_par_elastic, only: displ

  ! wavefield injection
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE
  ! boundary injection wavefield parts for saving together with b_absorb_field
  use specfem_par_coupling, only: b_boundary_injection_field

  ! Stacey approximations
  use stacey_par, only: compute_stacey_viscoelastic

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! communication overlap
  integer,intent(in) :: iphase

  ! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! adjoint simulations
  integer,intent(in) :: it
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  ! local parameters
  integer :: ispec,iglob,i,j,k,iface,igll

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! injecting boundary wavefield
  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE == 1) then
    ! adds boundary contribution from injected wavefield
    call compute_coupled_injection_contribution_el(NGLOB_AB,accel,iphase,it)
  endif

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Engquist)
  call compute_stacey_viscoelastic(NSPEC_AB,NGLOB_AB,accel, &
                                   ibool, &
                                   abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                   abs_boundary_ijk,abs_boundary_ispec, &
                                   num_abs_boundary_faces, &
                                   displ,veloc,rho_vp,rho_vs, &
                                   ispec_is_elastic, &
                                   b_num_abs_boundary_faces,b_absorb_field)

  ! for kernel simulations: stores absorbed wavefield part
  if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
    ! adds boundary injection contribution to stacey contribution before saving to disk
    ! this avoids storing the boundary injection arrays as a separate file
    ! for kernel simulations to reconstruct forward wavefields.
    if (COUPLE_WITH_INJECTION_TECHNIQUE) then
      b_absorb_field(:,:,:) = b_absorb_field(:,:,:) + b_boundary_injection_field(:,:,:)
    endif

    ! writes out absorbing boundary value
    call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
  endif

  ! Kirchoff-Helmholtz integrals
  if (SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
      if (ispec_is_elastic(ispec)) then
        ! reference GLL points on boundary face
        do igll = 1,NGLLSQUARE
          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)
          iglob = ibool(i,j,k,ispec)
          write(237) b_absorb_field(1,igll,iface), b_absorb_field(2,igll,iface), b_absorb_field(3,igll,iface)
          write(238) displ(1,iglob), displ(2,iglob), displ(3,iglob)
        enddo
      endif
    enddo
  endif

  ! wavefield injection
  ! not used and implemented yet..
  !if (COUPLE_WITH_INJECTION_TECHNIQUE) then !! To verify for NOBU version
  !  if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
  !     if (iphase == 1) it_dsm = it_dsm + 1
  !     !! TODO: maybe call integrand_for_computing_Kirchoff_Helmholtz_integral here
  !  endif
  !endif

  end subroutine compute_stacey_viscoelastic_forward

!
!=====================================================================
!

! for elastic solver

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_backward(NSPEC_AB, &
                                                  ibool,iphase, &
                                                  abs_boundary_ijk,abs_boundary_ispec, &
                                                  num_abs_boundary_faces, &
                                                  ispec_is_elastic,SIMULATION_TYPE, &
                                                  NSTEP,it,NGLOB_ADJOINT,b_accel, &
                                                  b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants
  use specfem_par, only: myrank

  implicit none

  integer,intent(in) :: NSPEC_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! communication overlap
  integer,intent(in) :: iphase

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! adjoint simulations
  integer,intent(in) :: SIMULATION_TYPE
  integer,intent(in) :: NSTEP,it,NGLOB_ADJOINT
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_ADJOINT),intent(inout) :: b_accel

  ! local parameters
  integer :: ispec,iglob,i,j,k,iface,igll

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling routine compute_stacey_viscoelastic_backward() with wrong SIMULATION_TYPE')

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! adjoint simulations:
  ! reads in absorbing boundary array (when first phase is running)
  ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
  call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Engquist)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_elastic(ispec)) then
      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

        ! adjoint simulations
        b_accel(:,iglob) = b_accel(:,iglob) - b_absorb_field(:,igll,iface)
      enddo
    endif ! ispec_is_elastic
  enddo

  end subroutine compute_stacey_viscoelastic_backward

!
!=====================================================================
!

  subroutine compute_stacey_viscoelastic_backward_undoatt(NSPEC_AB,NGLOB_AB,b_accel,b_veloc, &
                                                          ibool,iphase, &
                                                          abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                          abs_boundary_ijk,abs_boundary_ispec, &
                                                          num_abs_boundary_faces, &
                                                          rho_vp,rho_vs,ispec_is_elastic)

  use constants
  use specfem_par, only: myrank,SIMULATION_TYPE

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: b_accel
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: b_veloc

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! communication overlap
  integer,intent(in) :: iphase

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rho_vp,rho_vs
  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! local parameters
  real(kind=CUSTOM_REAL) :: vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling routine compute_stacey_viscoelastic_backward() with wrong SIMULATION_TYPE')

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Engquist)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_elastic(ispec)) then
      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

        vx = b_veloc(1,iglob)
        vy = b_veloc(2,iglob)
        vz = b_veloc(3,iglob)

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
        b_accel(1,iglob) = b_accel(1,iglob) - tx*jacobianw
        b_accel(2,iglob) = b_accel(2,iglob) - ty*jacobianw
        b_accel(3,iglob) = b_accel(3,iglob) - tz*jacobianw
      enddo
    endif ! ispec_is_elastic
  enddo

  end subroutine compute_stacey_viscoelastic_backward_undoatt

!
!=====================================================================
!

! for elastic solver on GPU

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
                                             NSTEP,it, &
                                             b_num_abs_boundary_faces,b_reclen_field,b_absorb_field,Mesh_pointer, &
                                             FORWARD_OR_ADJOINT)

  use constants

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,UNDO_ATTENUATION_AND_OR_PML

  ! wavefield injection
  ! use specfem_par, only: NGLOB_AB
  ! use specfem_par_elastic, only: accel
  ! boundary injection wavefield parts for saving together with b_absorb_field
  use specfem_par_coupling, only: b_boundary_injection_field

  implicit none

  ! communication overlap
  integer,intent(in) :: iphase

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces

  ! adjoint simulations
  integer,intent(in) :: NSTEP,it
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  ! GPU_MODE variables
  integer(kind=8),intent(in) :: Mesh_pointer
  integer, intent(in) :: FORWARD_OR_ADJOINT

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! wavefield injection
  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE == 1) then
    ! note: wavefield injection arrays are only available on CPU
    !
    !       as a quick work-around, we transfer here the velocity and acceleration arrays between GPU-CPU and vice versa,
    !       until the full injection contribution will be implemented on the GPU side as CUDA kernels.
    !       this is slowing down the simulation a bit.
    !
    ! transfers acceleration to the CPU
    ! call transfer_accel_from_device(NDIM*NGLOB_AB,accel, Mesh_pointer)

    ! ! adds boundary contribution from injected wavefield
    ! call compute_coupled_injection_contribution_el(NGLOB_AB,accel,iphase,it)

    ! ! transfers updated acceleration field back to the GPU
    ! call transfer_accel_to_device(NDIM*NGLOB_AB,accel, Mesh_pointer)

    ! now handled on GPU
    call compute_coupled_injection_contribution_el_GPU(iphase,Mesh_pointer)
  endif

  if (UNDO_ATTENUATION_AND_OR_PML) then
    ! no need to store boundaries on disk
    ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
    call compute_stacey_viscoelastic_undoatt_cuda(Mesh_pointer,iphase,FORWARD_OR_ADJOINT)
  else
    ! adjoint simulations:
    if (SIMULATION_TYPE == 3) then
      ! reads in absorbing boundary array (when first phase is running)
      ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
      call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)
    endif !adjoint

    call compute_stacey_viscoelastic_cuda(Mesh_pointer,iphase,b_absorb_field,FORWARD_OR_ADJOINT)

    ! for kernel simulations: stores absorbed wavefield part
    if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
      ! adds boundary injection contribution to stacey contribution before saving to disk
      ! this avoids storing the boundary injection arrays as a separate file
      ! for kernel simulations to reconstruct forward wavefields.
      if (COUPLE_WITH_INJECTION_TECHNIQUE) then
        b_absorb_field(:,:,:) = b_absorb_field(:,:,:) + b_boundary_injection_field(:,:,:)
      endif

      ! writes out absorbing boundary value
      call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
    endif
  endif

  ! wavefield injection
  ! not used and implemented yet..
  !if (COUPLE_WITH_INJECTION_TECHNIQUE) then !! To verify for NOBU version
  !  if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
  !     if (iphase == 1) it_dsm = it_dsm + 1
  !     !! TODO: maybe call integrand_for_computing_Kirchoff_Helmholtz_integral here
  !  endif
  !endif

  end subroutine compute_stacey_viscoelastic_GPU


!=============================================================================
!
! For coupling with external code
!
!=============================================================================

  subroutine compute_coupled_injection_contribution_el(NGLOB_AB,accel,iphase,it)

  use constants

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE

  use specfem_par, only: abs_boundary_normal,abs_boundary_jacobian2Dw, &
    abs_boundary_ijk,abs_boundary_ispec, &
    num_abs_boundary_faces

  use specfem_par, only: ibool
  use specfem_par_elastic, only: ispec_is_elastic,rho_vp,rho_vs

  ! boundary coupling
  !use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,INJECTION_TECHNIQUE_TYPE,RECIPROCITY_AND_KH_INTEGRAL
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,INJECTION_TECHNIQUE_TYPE
  ! use specfem_par_coupling, only: it_dsm, &
  !   Veloc_dsm_boundary, Tract_dsm_boundary, Veloc_axisem, Tract_axisem, Tract_axisem_time, &
  !   Veloc_specfem, Tract_specfem
  use specfem_par_coupling, only: it_dsm, &
    Veloc_dsm_boundary, Tract_dsm_boundary, Veloc_axisem, Tract_axisem, &
    Veloc_specfem, Tract_specfem
  ! FK3D calculation
  use specfem_par_coupling, only: ipt_table, NP_RESAMP, Veloc_FK, Tract_FK
  ! boundary injection wavefield parts for saving together with b_absorb_field
  use specfem_par_coupling, only: b_boundary_injection_field

  implicit none

  integer,intent(in) :: NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel

  ! communication overlap
  integer,intent(in) :: iphase

  ! adjoint simulations
  integer,intent(in) :: it

  ! local parameters
  real(kind=CUSTOM_REAL) :: vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll
  ! for the FK3D calculation
  ! FK surface
  integer :: ipt, ii, kk, iim1, iip1, iip2
  real(kind=CUSTOM_REAL) :: cs1,cs2,cs3,cs4,w
  real(kind=CUSTOM_REAL) :: vx_FK,vy_FK,vz_FK,tx_FK,ty_FK,tz_FK

!! comment from Vadim Monteiller, Feb 2017:

! txxbd est calcule dans add_to_compute_stacey_viscoelastic_1.F90 qui lui-meme appelle une subroutine qui se trouve
! dans add_to_compute_stacey_viscoelastic_11.F90. En fait je ne peux pas directement stocker txxbd en memoire
! (sinon il faut 500Go de ram, mais j'ai pas ca). Donc tous les 100 pas de temps je fais une fft
! et je prends la partie du sismo qui m'interesse. C'est a ce moment qu txxbd est rempli. C'est le call suivant qui fait ca:
!
! call store_next_FK_solution( VX_f, VY_f, VZ_f, TX_f, TY_f, TZ_f, &
!                    WKS_CMPLX_FOR_FFT, WKS_REAL_FOR_FFT, NF_FOR_STORING, &
!                    NF_FOR_FFT, NTIME_BETWEEN_FFT, NPOW_FOR_FFT, &
!                    vxbd, vybd, vzbd, txxbd, tyybd, tzzbd, npt, it, deltat)

! je prends la solution en frequence : VX_f, VY_f, VZ_f, TX_f, TY_f, TZ_f
! et je la sors en temps : vxbd, vybd, vzbd, txxbd, tyybd, tzzbd, pour les 100 pas de temps suivants.

! la subroutine store_next_FK_solution est dans add_to_compute_stacey_viscoelastic_11.F90 ligne 882

! on stocke directement la traction (on fait toujours ca avec DSM et AxiSEM aussi),
! ca evite de stocker 6 composantes de stress, surtout qu'on a des problemes de memoire.

  ! safety checks
  if (.not. COUPLE_WITH_INJECTION_TECHNIQUE) return

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! only for forward wavefield
  if (SIMULATION_TYPE /= 1) return

  ! injecting boundary
  !
  ! note: Following Tong et al. (2014, GJI, 197 (1); https://doi.org/10.1093/gji/ggt508) or
  !       Tong et al. (2014, GRL, 41; https://doi.org/10.1002/2014GL061644), coupling the wavefield
  !       by injecting the displacement and stresses at the mesh boundary can be merged together
  !       with the Clayton-Engquist absorbing boundary by:
  !            ( T_total - T_FK ) * n = - rho alpha [ n * d\dt(u_total - u_FK)]n - rho beta [ t * d/dt(u_total - u_FK)]t
  !       see equation (1) in Tong et al.'s GJI paper; equivalently using velocities v:
  !            ( T_total - T_FK ) * n = - rho alpha [ n * (v_total - v_FK)]n - rho beta [ t * (v_total - v_FK)]t
  !        in this case, the velocity rather than the displacement is needed and the Clayton-Engquist formulation modified.
  !        this allows to absorb the outgoing wavefield components, while still injecting an external wavefield.
  !
  !       -> option A: (original implementation) use the above modified Clayton-Engquist formulation to inject the wavefield
  !                    and absorb outgoing components
  !
  !       Changing the Clayton-Engquist boundary expression however would require to inject the wavefield with a
  !       different formulation. thus, instead of using the modified equation (1) from above, an equivalent way
  !       would be to inject the displacement and stress from the wavefield before the absorbing boundary gets computed,
  !       separating absorbing boundary implementations from the wavefield injection.
  !
  !       -> option B: 1. take only T_FK as a right-hand-side term to update acceleration
  !                    2. apply the absorbing boundary stresses based on the injected velocity:
  !                         T_abs_inj * n =  - rho alpha [ n * (- v_FK)]n - rho beta [ t * (- v_FK)]t
  !                    3. add the final stress as an additional right-hand-side term due to the injection:
  !                         T_final * n = (- T_FK + T_abs_inj) * n
  !                    this might help to separate absorbing boundary and wavefield injection, and could also be slightly
  !                    faster in case this coupling is not used as we avoid additional if-statements within the loops.
  !
  !       we now implement option B, as this also allows an easier integration for GPU simulations.
  !       thus, coupling with wavefield injection now is supported for both GPU and kernel simulations.
  !
  !       however, coupling boundary wavefields still requires to have Stacey absorbing boundaries set as well.
  !       this is due to the fact that we only call this routine within the Stacey routine,
  !       and also because we store the coupling contribution together with the Stacey ones for reconstructing the wavefield
  !       in kernels simulations.

  ! gets velocity & stress for boundary points
  select case(INJECTION_TECHNIQUE_TYPE)
  ! this is now read in by routine fetch_injection_wavefield() ...
  ! NQDU comment: don't need
  !
  ! case (INJECTION_TECHNIQUE_IS_DSM)
  !   ! DSM coupling
  !   if (old_DSM_coupling_from_Vadim) then
  !     if (mod(it_dsm,Ntime_step_dsm+1) == 0 .or. it == 1) then
  !       call read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)
  !     endif
  !   else
  !     !! MODIFS DE NOBU 2D
  !   endif
  ! case (INJECTION_TECHNIQUE_IS_AXISEM)
  !   ! AxiSEM coupling
  !   call read_axisem_file(Veloc_axisem,Tract_axisem,num_abs_boundary_faces*NGLLSQUARE)
  !   !! CD CD add this
  !   if (RECIPROCITY_AND_KH_INTEGRAL) Tract_axisem_time(:,:,it) = Tract_axisem(:,:)
  ! case (INJECTION_TECHNIQUE_IS_SPECFEM)
  !   ! SPECFEM coupling
  !   call read_specfem_file(Veloc_specfem,Tract_specfem,num_abs_boundary_faces*NGLLSQUARE,it)

  case (INJECTION_TECHNIQUE_IS_FK)
    ! FK coupling
    !! find indices
    ! example:
    !   np_resamp = 1 and it = 1,2,3,4,5,6, ..
    !   --> ii = 1,2,3,4,5,6,..,NSTEP
    !   np_resamp = 2 and it = 1,2,3,4,5,6, ..
    !   --> ii = 1,1,2,2,3,3,..,NSTEP/2
    ii = floor( real(it + NP_RESAMP - 1) / real( NP_RESAMP))
    ! example:
    !       kk = 1,2,1,2,1,2,,..
    kk = it - (ii-1) * NP_RESAMP
    ! example:
    !       w = 0,1/2,0,1/2,..
    w = dble(kk-1) / dble(NP_RESAMP)

    ! Cubic spline values
    cs4 = w*w*w/6.d0
    cs1 = 1.d0/6.d0 + w*(w-1.d0)/2.d0 - cs4
    cs3 = w + cs1 - 2.d0*cs4
    cs2 = 1.d0 - cs1 - cs3 - cs4

    ! interpolation indices
    iim1 = ii-1        ! 0,..
    iip1 = ii+1        ! 2,..
    iip2 = ii+2        ! 3,..
  end select

  ! wavefield injection
  ! option B: with injection of velocity/stress before the absorbing boundary condition
  !           as separate terms to the right-hand-side
  do iface = 1,num_abs_boundary_faces
    ispec = abs_boundary_ispec(iface)
    if (ispec_is_elastic(ispec)) then
      ! GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! global index
        iglob = ibool(i,j,k,ispec)

        ! gets injected velocity & stress
        ! velocity and stresses would be subtracted from total, therefore we add a minus sign when getting the values
        select case(INJECTION_TECHNIQUE_TYPE)
        case (INJECTION_TECHNIQUE_IS_DSM)                 !! To verify for NOBU version
          ! DSM coupling
          ! velocity
          it_dsm = it_dsm - 1

          vx = - Veloc_dsm_boundary(1,it_dsm,igll,iface)
          vy = - Veloc_dsm_boundary(2,it_dsm,igll,iface)
          vz = - Veloc_dsm_boundary(3,it_dsm,igll,iface)
          ! stress
          tx = - Tract_dsm_boundary(1,it_dsm,igll,iface)
          ty = - Tract_dsm_boundary(2,it_dsm,igll,iface)
          tz = - Tract_dsm_boundary(3,it_dsm,igll,iface)

          it_dsm = it_dsm + 1

        case (INJECTION_TECHNIQUE_IS_AXISEM)
          ! AxiSEM coupling
          ipt = igll + NGLLSQUARE*(iface - 1)
          ! velocity
          vx = - Veloc_axisem(1,ipt)
          vy = - Veloc_axisem(2,ipt)
          vz = - Veloc_axisem(3,ipt)
          ! stress
          tx = - Tract_axisem(1,ipt)
          ty = - Tract_axisem(2,ipt)
          tz = - Tract_axisem(3,ipt)

        case (INJECTION_TECHNIQUE_IS_SPECFEM)
          ! SPECFEM coupling
          ! indexing (same as Axisem)
          ipt = igll + NGLLSQUARE*(iface - 1)
          ! velocity
          vx = - Veloc_specfem(1,ipt)
          vy = - Veloc_specfem(2,ipt)
          vz = - Veloc_specfem(3,ipt)
          ! stress
          tx = - Tract_specfem(1,ipt)
          ty = - Tract_specfem(2,ipt)
          tz = - Tract_specfem(3,ipt)

        case (INJECTION_TECHNIQUE_IS_FK)
          ! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
          ! point index using table lookup
          ipt = ipt_table(igll,iface)

          ! interpolates velocity/stress
          vx_FK = cs1 * Veloc_FK(1,ipt,iim1) + cs2 * Veloc_FK(1,ipt,ii) + cs3 * Veloc_FK(1,ipt,iip1) + cs4 * Veloc_FK(1,ipt,iip2)
          vy_FK = cs1 * Veloc_FK(2,ipt,iim1) + cs2 * Veloc_FK(2,ipt,ii) + cs3 * Veloc_FK(2,ipt,iip1) + cs4 * Veloc_FK(2,ipt,iip2)
          vz_FK = cs1 * Veloc_FK(3,ipt,iim1) + cs2 * Veloc_FK(3,ipt,ii) + cs3 * Veloc_FK(3,ipt,iip1) + cs4 * Veloc_FK(3,ipt,iip2)

          tx_FK = cs1 * Tract_FK(1,ipt,iim1) + cs2 * Tract_FK(1,ipt,ii) + cs3 * Tract_FK(1,ipt,iip1) + cs4 * Tract_FK(1,ipt,iip2)
          ty_FK = cs1 * Tract_FK(2,ipt,iim1) + cs2 * Tract_FK(2,ipt,ii) + cs3 * Tract_FK(2,ipt,iip1) + cs4 * Tract_FK(2,ipt,iip2)
          tz_FK = cs1 * Tract_FK(3,ipt,iim1) + cs2 * Tract_FK(3,ipt,ii) + cs3 * Tract_FK(3,ipt,iip1) + cs4 * Tract_FK(3,ipt,iip2)

          ! velocity
          vx = - vx_FK
          vy = - vy_FK
          vz = - vz_FK
          ! stress
          tx = - tx_FK
          ty = - ty_FK
          tz = - tz_FK
        end select

        ! computes absorbing boundary for injected velocity
        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! gets associated normal
        nx = abs_boundary_normal(1,igll,iface)
        ny = abs_boundary_normal(2,igll,iface)
        nz = abs_boundary_normal(3,igll,iface)

        ! velocity component in normal direction (normal points out of element)
        vn = vx*nx + vy*ny + vz*nz

        ! adds stacey term to injected stresses:
        ! velocity vector component * vp * rho in normal direction + vs * rho component tangential to it
        tx = tx + rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
        ty = ty + rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
        tz = tz + rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

        ! adds final stress term for injected wavefield (weak form)
        accel(1,iglob) = accel(1,iglob) - tx*jacobianw
        accel(2,iglob) = accel(2,iglob) - ty*jacobianw
        accel(3,iglob) = accel(3,iglob) - tz*jacobianw

        ! for kernel simulations: stores contribution to buffer array and add it to stacey buffer before saving to disk
        if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
          b_boundary_injection_field(1,igll,iface) = tx*jacobianw
          b_boundary_injection_field(2,igll,iface) = ty*jacobianw
          b_boundary_injection_field(3,igll,iface) = tz*jacobianw
        endif
      enddo
    endif ! elastic
  enddo

  end subroutine compute_coupled_injection_contribution_el


!=============================================================================
!
! For coupling with external code, GPU version
!
!=============================================================================

  subroutine compute_coupled_injection_contribution_el_GPU(iphase,Mesh_pointer)

  use constants
  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE
  use specfem_par, only: num_abs_boundary_faces
  ! boundary coupling
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE
  ! FK3D calculation
  use specfem_par_coupling, only: b_boundary_injection_field

  implicit none
  ! communication overlap
  integer,intent(in) :: iphase
  ! GPU_MODE variables
  integer(kind=8),intent(in) :: Mesh_pointer

  ! safety checks
  if (.not. COUPLE_WITH_INJECTION_TECHNIQUE) return

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! only for forward wavefield
  if (SIMULATION_TYPE /= 1) return

  ! compute contribution in device
  call compute_coupled_injection_contribution_el_device(Mesh_pointer,b_boundary_injection_field, &
                                                        SAVE_STACEY)

  end subroutine compute_coupled_injection_contribution_el_GPU
