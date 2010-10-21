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

! allocates memory
  allocate(rmass(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(rmass_acoustic(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'
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
            
            stop 'poroelastic mass matrices not implemented yet'
            
            !rho_solid = density(1,kmato(ispec))
            !rho_fluid = density(2,kmato(ispec))
            !phi = porosity(kmato(ispec))
            !tort = tortuosity(kmato(ispec))
            !rho_bar = (1._CUSTOM_REAL-phil)*rhol_s + phil*rhol_f          
            !
            !if(CUSTOM_REAL == SIZE_REAL) then            
            !  ! for the solid mass matrix
            !  rmass_solid_poroelastic(iglob) = rmass_solid_poroelastic(iglob) + &
            !      sngl( dble(jacobianl) * weight * dble(rho_bar - phi*rho_fluid/tort) )
            !  
            !  ! for the fluid mass matrix
            !  rmass_fluid_poroelastic(iglob) = rmass_fluid_poroelastic(iglob) + &
            !      sngl( dble(jacobianl) * weight * dble(rho_bar*rho_fluid*tort - &
            !                                  phi*rho_fluid*rho_fluid)/dble(rho_bar*phi) )            
            !else
            !  rmass_solid_poroelastic(iglob) = rmass_solid_poroelastic(iglob) + &
            !      jacobianl * weight * (rho_bar - phi*rho_fluid/tort)
            !  
            !  rmass_fluid_poroelastic(iglob) = rmass_fluid_poroelastic(iglob) + &
            !      jacobianl * weight * (rho_bar*rho_fluid*tort - &
            !                                  phi*rho_fluid*rho_fluid) / (rho_bar*phi) 
            !endif
          endif
          
        enddo
      enddo
    enddo
  enddo ! nspec  

  end subroutine create_mass_matrices

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_ocean_load(nglob,nspec,ibool,OCEANS,&
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,NX_TOPO,NY_TOPO, &
                        ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO, &
                        itopo_bathy)

! returns precomputed mass matrix in rmass array
  
  use create_regions_mesh_ext_par 
  implicit none

! number of spectral elements in each block
  integer :: nspec
  integer :: nglob
  
! arrays with the mesh global indices
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  logical :: OCEANS

! use integer array to store topography values
  integer :: UTM_PROJECTION_ZONE
  logical :: SUPPRESS_UTM_PROJECTION
  integer :: NX_TOPO,NY_TOPO
  double precision :: ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO
  integer, dimension(NX_TOPO,NY_TOPO) :: itopo_bathy

  
! local parameters
  double precision :: weight
  double precision :: xval,yval,long,lat,elevation
  double precision :: height_oceans
  double precision :: long_corner,lat_corner,ratio_xi,ratio_eta
  integer :: ix_oceans,iy_oceans,iz_oceans,ispec_oceans,ispec2D,igll,iglobnum
  integer :: icornerlong,icornerlat

! creates ocean load mass matrix
  if(OCEANS) then

    ! adding ocean load mass matrix at ocean bottom
    NGLOB_OCEAN = nglob
    allocate(rmass_ocean_load(NGLOB_OCEAN))

    ! create ocean load mass matrix for degrees of freedom at ocean bottom
    rmass_ocean_load(:) = 0._CUSTOM_REAL

    ! add contribution of the oceans for surface elements exactly at ocean bottom
    do ispec2D = 1,num_free_surface_faces

      ispec_oceans = free_surface_ispec(ispec2D)

      ! only adds contribution if top boundary is elastic, no need to add this approximate calculation
      ! if top is already acoustic/poroelastic
      if( ispec_is_elastic(ispec_oceans) ) then

        do igll=1,NGLLSQUARE
          ix_oceans = free_surface_ijk(1,igll,ispec2D)
          iy_oceans = free_surface_ijk(1,igll,ispec2D)
          iz_oceans = free_surface_ijk(1,igll,ispec2D)
        
          iglobnum=ibool(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

          ! compute local height of oceans

          ! get coordinates of current point
          xval = xstore_dummy(iglobnum)
          yval = ystore_dummy(iglobnum)

          ! project x and y in UTM back to long/lat since topo file is in long/lat
          call utm_geo(long,lat,xval,yval,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)

          ! get coordinate of corner in bathy/topo model
          icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
          icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1

          ! avoid edge effects and extend with identical point if outside model
          if(icornerlong < 1) icornerlong = 1
          if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
          if(icornerlat < 1) icornerlat = 1
          if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1

          ! compute coordinates of corner
          long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
          lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO

          ! compute ratio for interpolation
          ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
          ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO

          ! avoid edge effects
          if(ratio_xi < 0.) ratio_xi = 0.
          if(ratio_xi > 1.) ratio_xi = 1.
          if(ratio_eta < 0.) ratio_eta = 0.
          if(ratio_eta > 1.) ratio_eta = 1.

          ! interpolate elevation at current point
          elevation = &
                itopo_bathy(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
                itopo_bathy(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
                itopo_bathy(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
                itopo_bathy(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

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

    ! add regular mass matrix to ocean load contribution
    rmass_ocean_load(:) = rmass_ocean_load(:) + rmass(:)

  else

    ! allocate dummy array if no oceans
    NGLOB_OCEAN = 1
    allocate(rmass_ocean_load(NGLOB_OCEAN))

  endif

  end subroutine create_mass_matrices_ocean_load


