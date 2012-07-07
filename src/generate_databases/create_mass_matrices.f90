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

  subroutine create_mass_matrices(nglob,nspec,ibool)

! returns precomputed mass matrix in rmass array

  use create_regions_mesh_ext_par
  implicit none

! number of spectral elements in each block
  integer :: nspec
  integer :: nglob

! arrays with the mesh global indices
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl
  integer :: ispec,i,j,k,iglob,ier
  real(kind=CUSTOM_REAL) :: rho_s, rho_f, rho_bar, phi, tort

! allocates memory
  allocate(rmass(nglob),stat=ier); if(ier /= 0) stop 'error in allocate rmass'
  allocate(rmass_acoustic(nglob),stat=ier); if(ier /= 0) stop 'error in allocate rmass_acoustic'
  allocate(rmass_solid_poroelastic(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(rmass_fluid_poroelastic(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'

! creates mass matrix
  rmass(:) = 0._CUSTOM_REAL
  rmass_acoustic(:) = 0._CUSTOM_REAL
  rmass_solid_poroelastic(:) = 0._CUSTOM_REAL
  rmass_fluid_poroelastic(:) = 0._CUSTOM_REAL

  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)

          weight = wxgll(i)*wygll(j)*wzgll(k)
          jacobianl = jacobianstore(i,j,k,ispec)

! acoustic mass matrix
          if( ispec_is_acoustic(ispec) ) then
            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then
              rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                    sngl( dble(jacobianl) * weight / dble(kappastore(i,j,k,ispec)) )
            else
               rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                    jacobianl * weight / kappastore(i,j,k,ispec)
            endif
          endif

! elastic mass matrix
          if( ispec_is_elastic(ispec) ) then
            if(CUSTOM_REAL == SIZE_REAL) then
              rmass(iglob) = rmass(iglob) + &
                    sngl( dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)) )
            else
               rmass(iglob) = rmass(iglob) + &
                    jacobianl * weight * rhostore(i,j,k,ispec)
            endif
          endif

! poroelastic mass matrices
          if( ispec_is_poroelastic(ispec) ) then

            rho_s = rhoarraystore(1,i,j,k,ispec)
            rho_f = rhoarraystore(2,i,j,k,ispec)
            phi = phistore(i,j,k,ispec)
            tort = tortstore(i,j,k,ispec)
            rho_bar = (1._CUSTOM_REAL-phi)*rho_s + phi*rho_f

            if(CUSTOM_REAL == SIZE_REAL) then
              ! for the solid mass matrix
              rmass_solid_poroelastic(iglob) = rmass_solid_poroelastic(iglob) + &
                  sngl( dble(jacobianl) * weight * dble(rho_bar - phi*rho_f/tort) )

              ! for the fluid mass matrix
              rmass_fluid_poroelastic(iglob) = rmass_fluid_poroelastic(iglob) + &
                  sngl( dble(jacobianl) * weight * dble(rho_bar*rho_f*tort - &
                                              phi*rho_f*rho_f)/dble(rho_bar*phi) )
            else
              rmass_solid_poroelastic(iglob) = rmass_solid_poroelastic(iglob) + &
                  jacobianl * weight * (rho_bar - phi*rho_f/tort)

              rmass_fluid_poroelastic(iglob) = rmass_fluid_poroelastic(iglob) + &
                  jacobianl * weight * (rho_bar*rho_f*tort - &
                                              phi*rho_f*rho_f) / (rho_bar*phi)
            endif

          endif

        enddo
      enddo
    enddo
  enddo ! nspec

  end subroutine create_mass_matrices

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_ocean_load(nglob,nspec,ibool,OCEANS,TOPOGRAPHY, &
                                            UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                                            NX_TOPO,NY_TOPO,itopo_bathy)

! returns precomputed mass matrix in rmass array

  use create_regions_mesh_ext_par
  implicit none

  ! number of spectral elements in each block
  integer :: nspec
  integer :: nglob

  ! arrays with the mesh global indices
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  logical :: OCEANS,TOPOGRAPHY

  ! use integer array to store topography values
  integer :: UTM_PROJECTION_ZONE
  logical :: SUPPRESS_UTM_PROJECTION

  integer :: NX_TOPO,NY_TOPO
  integer, dimension(NX_TOPO,NY_TOPO) :: itopo_bathy

  ! local parameters
  double precision :: weight
  double precision :: elevation
  double precision :: height_oceans
  integer :: ix_oceans,iy_oceans,iz_oceans,ispec_oceans,ispec2D,igll,iglobnum
  integer :: ier

  real(kind=CUSTOM_REAL) :: xloc,yloc,loc_elevation

  ! creates ocean load mass matrix
  if(OCEANS) then

    ! adding ocean load mass matrix at ocean bottom
    NGLOB_OCEAN = nglob
    allocate(rmass_ocean_load(NGLOB_OCEAN),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rmass_ocean_load'

    ! create ocean load mass matrix for degrees of freedom at ocean bottom
    rmass_ocean_load(:) = 0._CUSTOM_REAL

    ! add contribution of the oceans for surface elements exactly at ocean bottom
    do ispec2D = 1,num_free_surface_faces

      ! gets corresponding ispec
      ispec_oceans = free_surface_ispec(ispec2D)

      ! only adds contribution if top boundary is elastic, no need to add this approximate calculation
      ! if top is already acoustic/poroelastic
      if( ispec_is_elastic(ispec_oceans) ) then

        do igll=1,NGLLSQUARE
          ix_oceans = free_surface_ijk(1,igll,ispec2D)
          iy_oceans = free_surface_ijk(2,igll,ispec2D)
          iz_oceans = free_surface_ijk(3,igll,ispec2D)

          iglobnum=ibool(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

          ! compute local height of oceans
          if( TOPOGRAPHY ) then

            ! takes elevation from topography file
            xloc = xstore_dummy(iglobnum)
            yloc = ystore_dummy(iglobnum)

            call get_topo_bathy_elevation(xloc,yloc,loc_elevation, &
                                        itopo_bathy,NX_TOPO,NY_TOPO, &
                                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION)

            elevation = dble(loc_elevation)

          else

            ! takes elevation from z-coordinate of mesh point
            elevation = zstore_dummy(iglobnum)

          endif

          ! suppress positive elevation, which means no oceans
          if(elevation >= - MINIMUM_THICKNESS_3D_OCEANS) then
            height_oceans = 0.d0
          else
            height_oceans = dabs(elevation)
          endif

          ! take into account inertia of water column
          weight = dble( free_surface_jacobian2Dw(igll,ispec2D)) &
                   * dble(RHO_OCEANS) * height_oceans

          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + sngl(weight)
          else
            rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + weight
          endif

        enddo ! igll
      endif ! ispec_is_elastic
    enddo ! num_free_surface_faces

    ! adds regular mass matrix to ocean load contribution
    rmass_ocean_load(:) = rmass_ocean_load(:) + rmass(:)

  else

    ! allocate dummy array if no oceans
    NGLOB_OCEAN = 1
    allocate(rmass_ocean_load(NGLOB_OCEAN),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array rmass_ocean_load'

  endif

  end subroutine create_mass_matrices_ocean_load
