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

! defines mass matrices

  subroutine define_mass_matrices_elastic(nglob,nspec,nspec_irregular,ibool,rhostore, &
                                          jacobianstore,irregular_element_number,jacobian_regular, &
                                          wxgll,wygll,wzgll,ispec_is_elastic, &
                                          rmass)

! returns precomputed elastic mass matrix

  use constants, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL

  implicit none

  integer, intent(in) :: nglob
  integer, intent(in) :: nspec
  integer, intent(in) :: nspec_irregular

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: rhostore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_irregular), intent(in) :: jacobianstore
  integer, dimension(nspec), intent(in) :: irregular_element_number
  real(kind=CUSTOM_REAL), intent(in) :: jacobian_regular

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX), intent(in) :: wxgll
  double precision, dimension(NGLLY), intent(in) :: wygll
  double precision, dimension(NGLLZ), intent(in) :: wzgll

  logical, dimension(nspec), intent(in) :: ispec_is_elastic

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: rmass

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl
  integer :: ispec,ispec_irreg,i,j,k,iglob

  ! elastic mass matrix
  rmass(:) = 0._CUSTOM_REAL

  do ispec = 1,nspec
    ! elastic domain
    if (ispec_is_elastic(ispec)) then

      ispec_irreg = irregular_element_number(ispec)
      if (ispec_irreg == 0) jacobianl = jacobian_regular

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)

            weight = wxgll(i)*wygll(j)*wzgll(k)
            if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

            rmass(iglob) = rmass(iglob) + &
                           real( dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)),kind=CUSTOM_REAL)
          enddo
        enddo
      enddo
    endif
  enddo

  end subroutine define_mass_matrices_elastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_mass_matrices_acoustic(nglob,nspec,nspec_irregular,ibool,kappastore, &
                                           jacobianstore,irregular_element_number,jacobian_regular, &
                                           wxgll,wygll,wzgll,ispec_is_acoustic, &
                                           rmass_acoustic)

! returns precomputed acoustic mass matrix

  use constants, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL

  implicit none

  integer, intent(in) :: nglob
  integer, intent(in) :: nspec
  integer, intent(in) :: nspec_irregular

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: kappastore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_irregular), intent(in) :: jacobianstore
  integer, dimension(nspec), intent(in) :: irregular_element_number
  real(kind=CUSTOM_REAL), intent(in) :: jacobian_regular

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX), intent(in) :: wxgll
  double precision, dimension(NGLLY), intent(in) :: wygll
  double precision, dimension(NGLLZ), intent(in) :: wzgll

  logical, dimension(nspec), intent(in) :: ispec_is_acoustic

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: rmass_acoustic

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl
  integer :: ispec,ispec_irreg,i,j,k,iglob

  ! acoustic mass matrix
  rmass_acoustic(:) = 0._CUSTOM_REAL

  do ispec = 1,nspec
    ! acoustic domain
    if (ispec_is_acoustic(ispec)) then

      ispec_irreg = irregular_element_number(ispec)
      if (ispec_irreg == 0) jacobianl = jacobian_regular

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)

            weight = wxgll(i)*wygll(j)*wzgll(k)
            if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

            ! distinguish between single and double precision for reals
            rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                                    real( dble(jacobianl) * weight / dble(kappastore(i,j,k,ispec)),kind=CUSTOM_REAL)
          enddo
        enddo
      enddo
    endif
  enddo

  end subroutine define_mass_matrices_acoustic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_mass_matrices_poroelastic(nglob,nspec,nspec_irregular,ibool,rhoarraystore,phistore,tortstore, &
                                              jacobianstore,irregular_element_number,jacobian_regular, &
                                              wxgll,wygll,wzgll,ispec_is_poroelastic, &
                                              rmass_solid_poroelastic,rmass_fluid_poroelastic)

! returns precomputed poroelastic mass matrices

  use constants, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL

  implicit none

  integer, intent(in) :: nglob
  integer, intent(in) :: nspec
  integer, intent(in) :: nspec_irregular

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: rhoarraystore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: phistore,tortstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_irregular), intent(in) :: jacobianstore
  integer, dimension(nspec), intent(in) :: irregular_element_number
  real(kind=CUSTOM_REAL), intent(in) :: jacobian_regular

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX), intent(in) :: wxgll
  double precision, dimension(NGLLY), intent(in) :: wygll
  double precision, dimension(NGLLZ), intent(in) :: wzgll

  logical, dimension(nspec), intent(in) :: ispec_is_poroelastic

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: rmass_solid_poroelastic,rmass_fluid_poroelastic

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl
  integer :: ispec,ispec_irreg,i,j,k,iglob
  real(kind=CUSTOM_REAL) :: rho_s, rho_f, rho_bar, phi, tort

  ! poroelastic mass matrices
  rmass_solid_poroelastic(:) = 0._CUSTOM_REAL
  rmass_fluid_poroelastic(:) = 0._CUSTOM_REAL

  do ispec = 1,nspec
    ! poroelastic domain
    if (ispec_is_poroelastic(ispec)) then

      ispec_irreg = irregular_element_number(ispec)
      if (ispec_irreg == 0) jacobianl = jacobian_regular

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)

            weight = wxgll(i)*wygll(j)*wzgll(k)
            if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

            rho_s = rhoarraystore(1,i,j,k,ispec)
            rho_f = rhoarraystore(2,i,j,k,ispec)
            phi = phistore(i,j,k,ispec)
            tort = tortstore(i,j,k,ispec)

            rho_bar = (1._CUSTOM_REAL-phi)*rho_s + phi*rho_f

            ! for the solid mass matrix
            rmass_solid_poroelastic(iglob) = rmass_solid_poroelastic(iglob) + &
                  real( dble(jacobianl) * weight * dble(rho_bar - phi*rho_f/tort), &
                       kind=CUSTOM_REAL)

            ! for the fluid mass matrix
            rmass_fluid_poroelastic(iglob) = rmass_fluid_poroelastic(iglob) + &
                  real( dble(jacobianl) * weight * dble(rho_bar*rho_f*tort - phi*rho_f*rho_f) / dble(rho_bar*phi), &
                       kind=CUSTOM_REAL)
          enddo
        enddo
      enddo
    endif
  enddo

  end subroutine define_mass_matrices_poroelastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_mass_matrices_ocean_load(nglob,nspec,ibool,xstore_unique,ystore_unique,zstore_unique, &
                                             num_free_surface_faces,free_surface_ispec,free_surface_ijk,free_surface_jacobian2Dw, &
                                             NX_TOPO,NY_TOPO,itopo_bathy, &
                                             ispec_is_elastic,rmass_ocean_load)

! returns precomputed mass matrix for ocean load

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLSQUARE, &
    MINIMUM_THICKNESS_3D_OCEANS,RHO_APPROXIMATE_OCEAN_LOAD

  use shared_parameters, only: TOPOGRAPHY

  implicit none

  integer,intent(in) :: nglob
  integer,intent(in) :: nspec

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  real(kind=CUSTOM_REAL),dimension(nglob) :: xstore_unique,ystore_unique,zstore_unique

  integer, intent(in) :: num_free_surface_faces
  integer, dimension(num_free_surface_faces), intent(in) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces), intent(in) :: free_surface_ijk
  real(kind=CUSTOM_REAL), dimension(NGLLSQUARE,num_free_surface_faces), intent(in) :: free_surface_jacobian2Dw

  integer, intent(in) :: NX_TOPO,NY_TOPO
  integer, dimension(NX_TOPO,NY_TOPO), intent(in) :: itopo_bathy

  logical, dimension(nspec), intent(in) :: ispec_is_elastic

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: rmass_ocean_load

  ! local parameters
  double precision :: weight
  double precision :: elevation
  double precision :: height_oceans
  integer :: ix_oceans,iy_oceans,iz_oceans,ispec_oceans,ispec2D,igll,iglob

  real(kind=CUSTOM_REAL) :: xloc,yloc,loc_elevation

  ! contribution of the oceans for surface elements exactly at ocean bottom
  rmass_ocean_load(:) = 0._CUSTOM_REAL

  do ispec2D = 1,num_free_surface_faces

    ! gets corresponding ispec
    ispec_oceans = free_surface_ispec(ispec2D)

    ! only adds contribution if top boundary is elastic, no need to add this approximate calculation
    ! if top is already acoustic/poroelastic
    if (ispec_is_elastic(ispec_oceans)) then

      do igll = 1,NGLLSQUARE
        ix_oceans = free_surface_ijk(1,igll,ispec2D)
        iy_oceans = free_surface_ijk(2,igll,ispec2D)
        iz_oceans = free_surface_ijk(3,igll,ispec2D)

        iglob = ibool(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

        ! compute local height of oceans
        if (TOPOGRAPHY) then

          ! takes elevation from topography file
          xloc = xstore_unique(iglob)
          yloc = ystore_unique(iglob)

          call get_topo_bathy_elevation(xloc,yloc,loc_elevation, &
                                        itopo_bathy,NX_TOPO,NY_TOPO)

          elevation = dble(loc_elevation)

        else

          ! takes elevation from z-coordinate of mesh point
          elevation = zstore_unique(iglob)

        endif

        ! suppress positive elevation, which means no oceans
        if (elevation >= - MINIMUM_THICKNESS_3D_OCEANS) then
          height_oceans = 0.d0
        else
          height_oceans = dabs(elevation)
        endif

        ! take into account inertia of water column
        weight = dble(free_surface_jacobian2Dw(igll,ispec2D)) * dble(RHO_APPROXIMATE_OCEAN_LOAD) * height_oceans

        ! distinguish between single and double precision for reals
        rmass_ocean_load(iglob) = rmass_ocean_load(iglob) + real(weight,kind=CUSTOM_REAL)

      enddo ! igll
    endif ! ispec_is_elastic
  enddo ! num_free_surface_faces

  end subroutine define_mass_matrices_ocean_load

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_mass_matrices_Stacey_elastic(nglob,nspec,DT,ibool,rho_vp,rho_vs, &
                                                 num_abs_boundary_faces,abs_boundary_ispec,abs_boundary_ijk, &
                                                 abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                 ispec_is_elastic, &
                                                 rmassx, rmassy, rmassz)

! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
!
! thus, the mass matrix must be replaced by three mass matrices including the "C" damping matrix

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,NGLLSQUARE

  implicit none

  integer,intent(in) :: nglob
  integer,intent(in) :: nspec

  double precision, intent(in) :: DT

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: rho_vp,rho_vs

  integer, intent(in) :: num_abs_boundary_faces
  integer, dimension(num_abs_boundary_faces), intent(in) :: abs_boundary_ispec
  integer, dimension(3,NGLLSQUARE,num_abs_boundary_faces), intent(in) :: abs_boundary_ijk
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_abs_boundary_faces), intent(in) :: abs_boundary_normal
  real(kind=CUSTOM_REAL), dimension(NGLLSQUARE,num_abs_boundary_faces), intent(in) :: abs_boundary_jacobian2Dw

  logical, dimension(nspec), intent(in) :: ispec_is_elastic

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: rmassx,rmassy,rmassz

  ! local parameters
  real(kind=CUSTOM_REAL) :: jacobianw
  real(kind=CUSTOM_REAL) :: deltat,deltatover2
  real(kind=CUSTOM_REAL) :: tx,ty,tz
  real(kind=CUSTOM_REAL) :: nx,ny,nz,vn
  integer :: ispec,iglob,i,j,k,iface,igll

  ! elastic domains

  ! use the non-dimensional time step to make the mass matrix correction
  deltat = real(DT,kind=CUSTOM_REAL)
  deltatover2 = real(0.5d0*DT,kind=CUSTOM_REAL)

  ! adds contributions to mass matrix to stabilize Stacey conditions
  rmassx(:) = 0._CUSTOM_REAL
  rmassy(:) = 0._CUSTOM_REAL
  rmassz(:) = 0._CUSTOM_REAL

  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    ! elastic elements
    if (ispec_is_elastic(ispec)) then
      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets associated normal
        nx = abs_boundary_normal(1,igll,iface)
        ny = abs_boundary_normal(2,igll,iface)
        nz = abs_boundary_normal(3,igll,iface)

        vn = deltatover2*(nx+ny+nz)

        ! C*deltat/2 contributions
        tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(deltatover2-vn*nx)
        ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(deltatover2-vn*ny)
        tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(deltatover2-vn*nz)

        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

        ! assembles mass matrix on global points
        rmassx(iglob) = rmassx(iglob) + tx*jacobianw
        rmassy(iglob) = rmassy(iglob) + ty*jacobianw
        rmassz(iglob) = rmassz(iglob) + tz*jacobianw
      enddo
    endif ! elastic
  enddo

  end subroutine define_mass_matrices_Stacey_elastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_mass_matrices_Stacey_acoustic(nglob,nspec,DT,ibool,rho_vp, &
                                                  num_abs_boundary_faces,abs_boundary_ispec,abs_boundary_ijk, &
                                                  abs_boundary_jacobian2Dw, &
                                                  ispec_is_acoustic, &
                                                  rmassz_acoustic)

! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
!
! thus, the acoustic mass matrix must be replaced by a mass matrix including the "C" damping matrix

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLSQUARE

  implicit none

  integer,intent(in) :: nglob
  integer,intent(in) :: nspec

  double precision, intent(in) :: DT

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: rho_vp

  integer, intent(in) :: num_abs_boundary_faces
  integer, dimension(num_abs_boundary_faces), intent(in) :: abs_boundary_ispec
  integer, dimension(3,NGLLSQUARE,num_abs_boundary_faces), intent(in) :: abs_boundary_ijk
  real(kind=CUSTOM_REAL), dimension(NGLLSQUARE,num_abs_boundary_faces), intent(in) :: abs_boundary_jacobian2Dw

  logical, dimension(nspec), intent(in) :: ispec_is_acoustic

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: rmassz_acoustic

  ! local parameters
  real(kind=CUSTOM_REAL) :: jacobianw
  real(kind=CUSTOM_REAL) :: deltat,deltatover2
  real(kind=CUSTOM_REAL) :: sn
  integer :: ispec,iglob,i,j,k,iface,igll

  ! acoustic domains

  ! use the non-dimensional time step to make the mass matrix correction
  deltat = real(DT,kind=CUSTOM_REAL)
  deltatover2 = real(0.5d0*DT,kind=CUSTOM_REAL)

  ! adds contributions to mass matrix to stabilize Stacey conditions
  rmassz_acoustic(:) = 0._CUSTOM_REAL

  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    ! acoustic element
    if (ispec_is_acoustic(ispec)) then

      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! C * DT/2 contribution
        sn = deltatover2/rho_vp(i,j,k,ispec)

        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! gets global index
        iglob = ibool(i,j,k,ispec)

        rmassz_acoustic(iglob) = rmassz_acoustic(iglob) + jacobianw*sn
      enddo
    endif ! acoustic
  enddo

  end subroutine define_mass_matrices_Stacey_acoustic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_mass_matrices_pml_elastic(nglob,nspec,nspec_irregular,DT,ibool,rhostore, &
                                              jacobianstore,irregular_element_number,jacobian_regular, &
                                              wxgll,wygll,wzgll,ispec_is_elastic, &
                                              nspec_cpml,is_CPML,CPML_regions,CPML_to_spec, &
                                              d_store_x,d_store_y,d_store_z, &
                                              K_store_x,K_store_y,K_store_z, &
                                              rmass)

! returns precomputed elastic mass matrix w/ PML elements

  use constants, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL

  ! PML
  use constants, only: CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: nglob
  integer, intent(in) :: nspec
  integer, intent(in) :: nspec_irregular

  double precision, intent(in) :: DT

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: rhostore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_irregular), intent(in) :: jacobianstore
  integer, dimension(nspec), intent(in) :: irregular_element_number
  real(kind=CUSTOM_REAL), intent(in) :: jacobian_regular

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX), intent(in) :: wxgll
  double precision, dimension(NGLLY), intent(in) :: wygll
  double precision, dimension(NGLLZ), intent(in) :: wzgll

  logical, dimension(nspec), intent(in) :: ispec_is_elastic

  ! PML
  integer, intent(in) :: nspec_cpml
  logical, dimension(nspec), intent(in) :: is_CPML
  integer, dimension(nspec_cpml), intent(in) :: CPML_regions,CPML_to_spec
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_cpml), intent(in) :: d_store_x,d_store_y,d_store_z
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_cpml), intent(in) :: K_store_x,K_store_y,K_store_z

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: rmass

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl,deltat,mat_coef
  integer :: ispec,ispec_irreg,iglob,i,j,k,ispec_CPML

  ! use the non-dimensional time step to make the mass matrix correction
  deltat = real(DT,kind=CUSTOM_REAL)

  ! elastic mass matrix
  rmass(:) = 0._CUSTOM_REAL

  ! loops over physical mesh elements
  do ispec = 1,nspec
    if (.not. is_CPML(ispec) .and. ispec_is_elastic(ispec)) then

      ispec_irreg = irregular_element_number(ispec)
      if (ispec_irreg == 0) jacobianl = jacobian_regular

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ! defines the material coefficient associated to the domain
            mat_coef = rhostore(i,j,k,ispec)
            iglob = ibool(i,j,k,ispec)

            weight = wxgll(i)*wygll(j)*wzgll(k)
            if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

            rmass(iglob) = rmass(iglob) + real( dble(jacobianl) * weight * dble(mat_coef),kind=CUSTOM_REAL)
          enddo
        enddo
      enddo
    endif
  enddo

  ! loops over C-PML elements
  do ispec_CPML = 1,nspec_cpml
    ispec = CPML_to_spec(ispec_CPML)

    if (is_CPML(ispec) .and. ispec_is_elastic(ispec)) then

      ispec_irreg = irregular_element_number(ispec)
      if (ispec_irreg == 0) jacobianl = jacobian_regular

      if (CPML_regions(ispec_CPML) == CPML_X_ONLY) then
        ! X_surface C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = rhostore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)

              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass(iglob) = rmass(iglob) + &
                real( dble(jacobianl) * weight * dble(mat_coef) * &
                     (dble(K_store_x(i,j,k,ispec_CPML)) + dble(d_store_x(i,j,k,ispec_CPML)) * deltat / 2.d0),kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_Y_ONLY) then
        ! Y_surface C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = rhostore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)

              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass(iglob) = rmass(iglob) + &
                real( dble(jacobianl) * weight * dble(mat_coef) * &
                     (dble(K_store_y(i,j,k,ispec_CPML)) + dble(d_store_y(i,j,k,ispec_CPML)) * deltat / 2.d0),kind=CUSTOM_REAL)
            enddo
          enddo
        enddo


      else if (CPML_regions(ispec_CPML) == CPML_Z_ONLY) then
        ! Z_surface C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = rhostore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)

              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass(iglob) = rmass(iglob) + &
                real( dble(jacobianl) * weight * dble(mat_coef) * &
                     (dble(K_store_z(i,j,k,ispec_CPML)) + dble(d_store_z(i,j,k,ispec_CPML)) * deltat / 2.d0),kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_XY_ONLY) then
        ! XY_edge C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = rhostore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)

              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass(iglob) = rmass(iglob) + &
                real( dble(jacobianl) * weight * dble(mat_coef) * &
                     (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) + &
                      (dble(d_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) + &
                       dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML))) * deltat / 2.d0),kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_XZ_ONLY) then
        ! XZ_edge C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = rhostore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)

              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass(iglob) = rmass(iglob) + &
                real( dble(jacobianl) * weight * dble(mat_coef) * &
                     (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                      (dble(d_store_x(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                       dble(d_store_z(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML))) * deltat / 2.d0),kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_YZ_ONLY) then
        ! YZ_edge C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = rhostore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)

              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass(iglob) = rmass(iglob) + &
                real( dble(jacobianl) * weight * dble(mat_coef) * &
                     (dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                      (dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                       dble(d_store_z(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML))) * deltat / 2.d0) ,kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_XYZ) then
        ! XYZ_corner C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = rhostore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)

              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass(iglob) = rmass(iglob) + &
                real( dble(jacobianl) * weight * dble(mat_coef) * &
                     (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) * &
                      dble(K_store_z(i,j,k,ispec_CPML)) + (dble(d_store_x(i,j,k,ispec_CPML)) * &
                      dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                      dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML)) * &
                      dble(K_store_z(i,j,k,ispec_CPML)) + dble(d_store_z(i,j,k,ispec_CPML)) * &
                      dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML))) * deltat / 2.d0) ,kind=CUSTOM_REAL)
            enddo
          enddo
        enddo
      else
        stop 'error in PML mesh file'
      endif
    endif
  enddo ! do ispec_CPML=1,nspec_cpml

  end subroutine define_mass_matrices_pml_elastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_mass_matrices_pml_acoustic(nglob,nspec,nspec_irregular,DT,ibool,kappastore, &
                                               jacobianstore,irregular_element_number,jacobian_regular, &
                                               wxgll,wygll,wzgll,ispec_is_acoustic, &
                                               nspec_cpml,is_CPML,CPML_regions,CPML_to_spec, &
                                               d_store_x,d_store_y,d_store_z, &
                                               K_store_x,K_store_y,K_store_z, &
                                               rmass_acoustic)

! returns precomputed acoustic mass matrix w/ PML elements

  use constants, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL

  ! PML
  use constants, only: CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
  CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: nglob
  integer, intent(in) :: nspec
  integer, intent(in) :: nspec_irregular

  double precision, intent(in) :: DT

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: kappastore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_irregular), intent(in) :: jacobianstore
  integer, dimension(nspec), intent(in) :: irregular_element_number
  real(kind=CUSTOM_REAL), intent(in) :: jacobian_regular

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX), intent(in) :: wxgll
  double precision, dimension(NGLLY), intent(in) :: wygll
  double precision, dimension(NGLLZ), intent(in) :: wzgll

  logical, dimension(nspec), intent(in) :: ispec_is_acoustic

  ! PML
  integer, intent(in) :: nspec_cpml
  logical, dimension(nspec), intent(in) :: is_CPML
  integer, dimension(nspec_cpml), intent(in) :: CPML_regions,CPML_to_spec
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_cpml), intent(in) :: d_store_x,d_store_y,d_store_z
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_cpml), intent(in) :: K_store_x,K_store_y,K_store_z

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: rmass_acoustic

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl,deltat,mat_coef
  integer :: ispec,ispec_irreg,iglob,i,j,k,ispec_CPML

  ! use the non-dimensional time step to make the mass matrix correction
  deltat = real(DT,kind=CUSTOM_REAL)

  ! acoustic mass matrix
  rmass_acoustic(:) = 0._CUSTOM_REAL

  ! loops over physical mesh elements
  do ispec = 1,nspec
    if (.not. is_CPML(ispec) .and. ispec_is_acoustic(ispec)) then

      ispec_irreg = irregular_element_number(ispec)
      if (ispec_irreg == 0) jacobianl = jacobian_regular

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ! defines the material coefficient associated to the domain
            mat_coef = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)
            iglob = ibool(i,j,k,ispec)
            weight = wxgll(i)*wygll(j)*wzgll(k)
            if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

            rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                         real( dble(jacobianl) * weight * dble(mat_coef) ,kind=CUSTOM_REAL)
          enddo
        enddo
      enddo
    endif
  enddo

  ! loops over C-PML elements
  do ispec_CPML=1,nspec_cpml
    ispec = CPML_to_spec(ispec_CPML)

    if (is_CPML(ispec) .and. ispec_is_acoustic(ispec)) then

      ispec_irreg = irregular_element_number(ispec)
      if (ispec_irreg == 0) jacobianl = jacobian_regular

      if (CPML_regions(ispec_CPML) == CPML_X_ONLY) then
        ! X_surface C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)
              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                        real( dble(jacobianl) * weight * dble(mat_coef) * &
                            (dble(K_store_x(i,j,k,ispec_CPML)) + dble(d_store_x(i,j,k,ispec_CPML)) &
                            * deltat / 2.d0) ,kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_Y_ONLY) then
        ! Y_surface C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)
              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                        real( dble(jacobianl) * weight * dble(mat_coef) * &
                            (dble(K_store_y(i,j,k,ispec_CPML)) + dble(d_store_y(i,j,k,ispec_CPML)) &
                            * deltat / 2.d0) ,kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_Z_ONLY) then
        ! Z_surface C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)
              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                        real( dble(jacobianl) * weight * dble(mat_coef) * &
                            (dble(K_store_z(i,j,k,ispec_CPML)) + dble(d_store_z(i,j,k,ispec_CPML)) &
                            * deltat / 2.d0) ,kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_XY_ONLY) then
        ! XY_edge C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)
              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                        real( dble(jacobianl) * weight * dble(mat_coef) * &
                            (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) + &
                            (dble(d_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) + &
                            dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML))) &
                            * deltat / 2.d0) ,kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_XZ_ONLY) then
        ! XZ_edge C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)
              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                        real( dble(jacobianl) * weight * dble(mat_coef) * &
                            (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                            (dble(d_store_x(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                            dble(d_store_z(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML))) &
                            * deltat / 2.d0) ,kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_YZ_ONLY) then
        ! YZ_edge C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)
              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                        real( dble(jacobianl) * weight * dble(mat_coef) * &
                            (dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                            (dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                            dble(d_store_z(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML))) &
                            * deltat / 2.d0) ,kind=CUSTOM_REAL)
            enddo
          enddo
        enddo

      else if (CPML_regions(ispec_CPML) == CPML_XYZ) then
        ! XYZ_corner C-PML
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! defines the material coefficient associated to the domain
              mat_coef = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)
              iglob = ibool(i,j,k,ispec)
              weight = wxgll(i)*wygll(j)*wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

              rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                        real( dble(jacobianl) * weight * dble(mat_coef) * &
                            (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) * &
                            dble(K_store_z(i,j,k,ispec_CPML)) + (dble(d_store_x(i,j,k,ispec_CPML)) * &
                            dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                            dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML)) * &
                            dble(K_store_z(i,j,k,ispec_CPML)) + dble(d_store_z(i,j,k,ispec_CPML)) * &
                            dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML))) * &
                            deltat / 2.d0) ,kind=CUSTOM_REAL)
            enddo
          enddo
        enddo
      else
        stop 'error in PML mesh file'
      endif
    endif
  enddo ! do ispec_CPML=1,nspec_cpml

  end subroutine define_mass_matrices_pml_acoustic
