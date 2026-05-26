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

! for elastic solver

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_forward(NSPEC_AB,NGLOB_AB,accel, &
                                                 ibool,iphase, &
                                                 abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                 abs_boundary_ijk,abs_boundary_ispec, &
                                                 num_abs_boundary_faces,veloc,rho_vp,rho_vs, &
                                                 ispec_is_elastic, &
                                                 it,b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE, &
    xstore,ystore,zstore,hprime_xx,rhostore

  ! Kirchoff-Helmholtz integrals
  use specfem_par_elastic, only: displ

  ! wavefield injection
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE
  ! boundary injection wavefield parts for saving together with b_absorb_field
  use specfem_par_coupling, only: b_boundary_injection_field

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
  real(kind=CUSTOM_REAL) :: vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll

  ! for improved Stacey condition using P3 approximation
  ! 2D local arrays for surface geometry and displacements
  real(kind=CUSTOM_REAL) :: x_2D(NGLLX, NGLLY), y_2D(NGLLX, NGLLY), z_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: un_2D(NGLLX, NGLLY), ut1_2D(NGLLX, NGLLY), ut2_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: nx_2D(NGLLX, NGLLY), ny_2D(NGLLX, NGLLY), nz_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: t1x_2D(NGLLX, NGLLY), t1y_2D(NGLLX, NGLLY), t1z_2D(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: t2x_2D(NGLLX, NGLLY), t2y_2D(NGLLX, NGLLY), t2z_2D(NGLLX, NGLLY)

  ! Derivative arrays
  real(kind=CUSTOM_REAL) :: dx_ds1(NGLLX, NGLLY), dx_ds2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dy_ds1(NGLLX, NGLLY), dy_ds2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dz_ds1(NGLLX, NGLLY), dz_ds2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dun_ds1(NGLLX, NGLLY), dun_ds2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dut1_ds1(NGLLX, NGLLY), dut1_ds2(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: dut2_ds1(NGLLX, NGLLY), dut2_ds2(NGLLX, NGLLY)

  real(kind=CUSTOM_REAL) :: a11, a12, a22, det_a, inv_a11, inv_a12, inv_a22
  real(kind=CUSTOM_REAL) :: dot_g1_t1, dot_g2_t1, dot_g1_t2, dot_g2_t2
  real(kind=CUSTOM_REAL) :: dt1_un, dt2_un, dt1_ut1, dt2_ut2

  real(kind=CUSTOM_REAL) :: rhol, csl, cpl
  real(kind=CUSTOM_REAL) :: t1x, t1y, t1z, t2x, t2y, t2z, t1_norm, hp1, hp2
  real(kind=CUSTOM_REAL) :: t2_n2, t2_t1, t2_t2

  logical :: mask_vary(3)
  integer :: id1, id2, a, b, l
  integer :: face_iglob(NGLLX, NGLLY)
  !debug
  !real(kind=CUSTOM_REAL) :: nt1, nt2, t1t2, rh

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
  if (.not. USE_SECOND_ORDER_STACEY) then
    ! uses Stacey P1 approximation
    !
    ! uses velocities like
    !   T_n = - rho Vp v_n                          (normal)
    !   T_tangential1 = - rho Vs v_tangential       (tangential)
    !
    ! units:  [Pa] =  [kg/m^3] [m/s] [m/s] = [kg / m / s^2 ]

! openmp solver
!$OMP PARALLEL if (num_abs_boundary_faces > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(iface,ispec,igll,i,j,k,iglob,vx,vy,vz,vn,nx,ny,nz,tx,ty,tz,jacobianw)
!$OMP DO
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
      endif ! ispec_is_elastic
    enddo
!$OMP ENDDO
!$OMP END PARALLEL

  else
    ! uses Stacey P3 approximation (second-order)
    !
    ! uses spatial derivatives of displacement like
    !   T_n = - rho Vs (2 Vs - Vp)( \partial_tangential1 u_tangential1 + \partial_tangential2 u_tangential2 )
    !
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
!$OMP         rhol,csl,cpl,t1x,t1y,t1z,t2x,t2y,t2z,t1_norm,hp1,hp2, &
!$OMP         t2_n2,t2_t1,t2_t2,mask_vary,id1,id2,a,b,l,face_iglob)
!$OMP DO
    do iface = 1,num_abs_boundary_faces

      ispec = abs_boundary_ispec(iface)

      if (ispec_is_elastic(ispec)) then

        ! prepare P3 approximation spatial derivative arrays
        ! Find the varying local coordinates for the face to map 1D index igll to 2D grid (a,b)
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

        ! Gather coordinates and fields onto the 2D grid face
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

        ! Pre-compute the orthonormal tangential basis and projections at each point
        ! (assumes NGLLX == NGLLY == NGLLZ)
        do b = 1,NGLLY
          do a = 1,NGLLX
            iglob = face_iglob(a,b)

            ! normal
            nx = nx_2D(a,b)
            ny = ny_2D(a,b)
            nz = nz_2D(a,b)

            ! Construct right-handed orthonormal tangential basis t1, t2
            ! Choose a non-collinear vector to n to compute t1
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

            ! debug
            !if (myrank == 0) &
            !  print '("debug:",i4," iface ",i8," iglob ",i8,"t1",3(1x,f12.6)," t2",3(1x,f12.6))', &
            !        myrank,iface,iglob,t1x,t1y,t1z,t2x,t2y,t2z
            ! check orthonormal basis orientation
            !nt1 = nx*t1x + ny*t1y + nz*t1z     ! dot-product n · t1
            !nt2 = nx*t2x + ny*t2y + nz*t2z     !             n · t2
            !t1t2 = t1x*t2x + t1y*t2y + t1z*t2z !             t1 · t2
            !rh = nx*(t1y*t2z - t1z*t2y) + ny*(t1z*t2x - t1x*t2z) + nz*(t1x*t2y - t1y*t2x)  ! right-handedness: n · (t1 x t2)
            !if (abs(rh - 1.d0) > 1.d-20 .or. abs(nt1) > 1.d-20 .or. abs(nt2) > 1.d-20 .or. abs(t1t2) > 1.d-20) &
            !  print '("debug:",i4," iface ",i8," iglob ",i8,3(1x,f12.6),1x,"rh",f12.6)',  &
            !          myrank,iface,iglob,nt1,nt2,t1t2,rh

            ! Store basis for reuse in igll loop
            t1x_2D(a,b) = t1x; t1y_2D(a,b) = t1y; t1z_2D(a,b) = t1z
            t2x_2D(a,b) = t2x; t2y_2D(a,b) = t2y; t2z_2D(a,b) = t2z

            ! Project displacement
            un_2D(a,b)  = displ(1,iglob) * nx  + displ(2,iglob) * ny  + displ(3,iglob) * nz
            ut1_2D(a,b) = displ(1,iglob) * t1x + displ(2,iglob) * t1y + displ(3,iglob) * t1z
            ut2_2D(a,b) = displ(1,iglob) * t2x + displ(2,iglob) * t2y + displ(3,iglob) * t2z
          enddo
        enddo

        ! Compute reference derivatives using 1D GLL derivative matrix
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

          ! gets associated normal
          nx = abs_boundary_normal(1,igll,iface)
          ny = abs_boundary_normal(2,igll,iface)
          nz = abs_boundary_normal(3,igll,iface)

          ! velocity component in normal direction (normal points out of element)
          vn = vx*nx + vy*ny + vz*nz

          ! P1 contribution
          ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it
          tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
          ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
          tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

          ! additional P3 second-order contribution
          a = abs_boundary_ijk(id1,igll,iface)
          b = abs_boundary_ijk(id2,igll,iface)

          ! Metric tensor components a_{alpha beta} = g_alpha . g_beta
          a11 = dx_ds1(a,b)**2 + dy_ds1(a,b)**2 + dz_ds1(a,b)**2
          a12 = dx_ds1(a,b) * dx_ds2(a,b) + dy_ds1(a,b) * dy_ds2(a,b) + dz_ds1(a,b) * dz_ds2(a,b)
          a22 = dx_ds2(a,b)**2 + dy_ds2(a,b)**2 + dz_ds2(a,b)**2

          det_a = a11*a22 - a12**2

          ! avoid divison by zero
          if (abs(det_a) < 1.d-24) det_a = 1._CUSTOM_REAL

          ! Inverse metric tensor
          inv_a11 = a22 / det_a
          inv_a12 = -a12 / det_a
          inv_a22 = a11 / det_a

          ! look up pre-computed basis and normal
          nx = nx_2D(a,b)
          ny = ny_2D(a,b)
          nz = nz_2D(a,b)
          t1x = t1x_2D(a,b); t1y = t1y_2D(a,b); t1z = t1z_2D(a,b)
          t2x = t2x_2D(a,b); t2y = t2y_2D(a,b); t2z = t2z_2D(a,b)

          dot_g1_t1 = dx_ds1(a,b)*t1x + dy_ds1(a,b)*t1y + dz_ds1(a,b)*t1z
          dot_g2_t1 = dx_ds2(a,b)*t1x + dy_ds2(a,b)*t1y + dz_ds2(a,b)*t1z

          dot_g1_t2 = dx_ds1(a,b)*t2x + dy_ds1(a,b)*t2y + dz_ds1(a,b)*t2z
          dot_g2_t2 = dx_ds2(a,b)*t2x + dy_ds2(a,b)*t2y + dz_ds2(a,b)*t2z

          dt1_un  = (inv_a11*dun_ds1(a,b) + inv_a12*dun_ds2(a,b)) * dot_g1_t1 &
                  + (inv_a12*dun_ds1(a,b) + inv_a22*dun_ds2(a,b)) * dot_g2_t1
          dt2_un  = (inv_a11*dun_ds1(a,b) + inv_a12*dun_ds2(a,b)) * dot_g1_t2 &
                  + (inv_a12*dun_ds1(a,b) + inv_a22*dun_ds2(a,b)) * dot_g2_t2

          dt1_ut1 = (inv_a11*dut1_ds1(a,b) + inv_a12*dut1_ds2(a,b)) * dot_g1_t1 &
                  + (inv_a12*dut1_ds1(a,b) + inv_a22*dut1_ds2(a,b)) * dot_g2_t1
          dt2_ut2 = (inv_a11*dut2_ds1(a,b) + inv_a12*dut2_ds2(a,b)) * dot_g1_t2 &
                  + (inv_a12*dut2_ds1(a,b) + inv_a22*dut2_ds2(a,b)) * dot_g2_t2

          ! T_normal2     = - rho * vs * (2 * vs - vp) (du_t1/dt1 + du_t2/dt2)
          ! T_tangential2 = + rho * vs * (2 * vs - vp) du_normal / dx_tangential
          rhol = rhostore(i,j,k,ispec)
          csl = rho_vs(i,j,k,ispec) / rhol
          cpl = rho_vp(i,j,k,ispec) / rhol

          ! note: accel gets subtracted by tx, ty, tz, so we add -T^(2) to tx, ty, tz
          t2_n2 = rho_vs(i,j,k,ispec) * (2.0_CUSTOM_REAL * csl - cpl) * (dt1_ut1 + dt2_ut2)
          t2_t1 = -rho_vs(i,j,k,ispec) * (2.0_CUSTOM_REAL * csl - cpl) * dt1_un
          t2_t2 = -rho_vs(i,j,k,ispec) * (2.0_CUSTOM_REAL * csl - cpl) * dt2_un

          tx = tx + t2_n2 * nx + t2_t1 * t1x + t2_t2 * t2x
          ty = ty + t2_n2 * ny + t2_t1 * t1y + t2_t2 * t2y
          tz = tz + t2_n2 * nz + t2_t1 * t1z + t2_t2 * t2z

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
      endif ! ispec_is_elastic
    enddo
!$OMP ENDDO
!$OMP END PARALLEL

  endif  ! USE_SECOND_ORDER_STACEY

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
  !! CD CD added this
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

