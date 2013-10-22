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

  subroutine create_mass_matrices(nglob,nspec,ibool,PML_CONDITIONS,STACEY_ABSORBING_CONDITIONS)

! returns precomputed mass matrix in rmass array

  use create_regions_mesh_ext_par

  implicit none

! number of spectral elements in each block
  integer :: nspec
  integer :: nglob

! arrays with the mesh global indices
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! C-PML flag
  logical :: PML_CONDITIONS

! Stacey boundary condition flag
  logical :: STACEY_ABSORBING_CONDITIONS

! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl
  integer :: ispec,i,j,k,iglob,ier
  real(kind=CUSTOM_REAL) :: rho_s, rho_f, rho_bar, phi, tort

  ! elastic domains
  if( ELASTIC_SIMULATION ) then
     ! allocates memory
     allocate(rmass(nglob),stat=ier); if(ier /= 0) stop 'error allocating array rmass'
     rmass(:) = 0._CUSTOM_REAL

     ! returns elastic mass matrix
     if( PML_CONDITIONS ) then
        call create_mass_matrices_pml_elastic(nspec,ibool)
     else
        do ispec=1,nspec
           if( ispec_is_elastic(ispec) ) then
              do k=1,NGLLZ
                 do j=1,NGLLY
                    do i=1,NGLLX
                       iglob = ibool(i,j,k,ispec)

                       weight = wxgll(i)*wygll(j)*wzgll(k)
                       jacobianl = jacobianstore(i,j,k,ispec)

                       if(CUSTOM_REAL == SIZE_REAL) then
                          rmass(iglob) = rmass(iglob) + &
                               sngl( dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)) )
                       else
                          rmass(iglob) = rmass(iglob) + &
                               jacobianl * weight * rhostore(i,j,k,ispec)
                       endif
                    enddo
                 enddo
              enddo
           endif
        enddo
     endif
  endif

  ! acoustic domains
  if( ACOUSTIC_SIMULATION) then
     ! allocates memory
     allocate(rmass_acoustic(nglob),stat=ier); if(ier /= 0) stop 'error allocating array rmass_acoustic'
     rmass_acoustic(:) = 0._CUSTOM_REAL

     ! returns acoustic mass matrix
     if( PML_CONDITIONS ) then
        call create_mass_matrices_pml_acoustic(nspec,ibool)
     else
        do ispec=1,nspec
           if( ispec_is_acoustic(ispec) ) then
              do k=1,NGLLZ
                 do j=1,NGLLY
                    do i=1,NGLLX
                       iglob = ibool(i,j,k,ispec)

                       weight = wxgll(i)*wygll(j)*wzgll(k)
                       jacobianl = jacobianstore(i,j,k,ispec)

                       ! distinguish between single and double precision for reals
                       if(CUSTOM_REAL == SIZE_REAL) then
                          rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                               sngl( dble(jacobianl) * weight / dble(kappastore(i,j,k,ispec)) )
                       else
                          rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                               jacobianl * weight / kappastore(i,j,k,ispec)
                       endif
                    enddo
                 enddo
              enddo
           endif
        enddo
     endif
  endif

  ! poroelastic domains
  if( POROELASTIC_SIMULATION) then
    ! allocates memory
    allocate(rmass_solid_poroelastic(nglob),stat=ier); if(ier /= 0) stop 'error in allocate rmass_solid_poroelastic'
    allocate(rmass_fluid_poroelastic(nglob),stat=ier); if(ier /= 0) stop 'error in allocate rmass_fluid_poroelastic'
    rmass_solid_poroelastic(:) = 0._CUSTOM_REAL
    rmass_fluid_poroelastic(:) = 0._CUSTOM_REAL

    ! poroelastic mass matrices
    do ispec=1,nspec
      if( ispec_is_poroelastic(ispec) ) then
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              iglob = ibool(i,j,k,ispec)

              weight = wxgll(i)*wygll(j)*wzgll(k)
              jacobianl = jacobianstore(i,j,k,ispec)

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
            enddo
          enddo
        enddo
      endif
    enddo ! nspec
  endif

  ! Stacey absorbing conditions (adds C*deltat/2 contribution to the mass matrices on Stacey edges)
  if(STACEY_ABSORBING_CONDITIONS)call create_mass_matrices_Stacey(nglob,nspec,ibool)

  ! ocean load mass matrix
  call create_mass_matrices_ocean_load(nglob,nspec,ibool)

  end subroutine create_mass_matrices

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_ocean_load(nglob,nspec,ibool)

! returns precomputed mass matrix in rmass array

  use generate_databases_par, only: &
    APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
    NX_TOPO,NY_TOPO,itopo_bathy,myrank

  use create_regions_mesh_ext_par

  implicit none

  ! number of spectral elements in each block
  integer :: nspec
  integer :: nglob

  ! arrays with the mesh global indices
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! local parameters
  double precision :: weight
  double precision :: elevation
  double precision :: height_oceans
  integer :: ix_oceans,iy_oceans,iz_oceans,ispec_oceans,ispec2D,igll,iglobnum
  integer :: ier

  real(kind=CUSTOM_REAL) :: xloc,yloc,loc_elevation

  ! creates ocean load mass matrix
  if(APPROXIMATE_OCEAN_LOAD) then

    if( myrank == 0) then
      write(IMAIN,*) '  ...creating ocean load mass matrix '
    endif

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
                   * dble(RHO_APPROXIMATE_OCEAN_LOAD) * height_oceans

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

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_Stacey(nglob,nspec,ibool)

! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix

  use generate_databases_par, only: DT
  use create_regions_mesh_ext_par

  implicit none

  integer :: nglob
  integer :: nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! local parameters
  real(kind=CUSTOM_REAL) :: jacobianw
  real(kind=CUSTOM_REAL) :: deltat,deltatover2
  real(kind=CUSTOM_REAL) :: tx,ty,tz,sn
  real(kind=CUSTOM_REAL) :: nx,ny,nz,vn
  integer :: ispec,iglob,i,j,k,iface,igll,ier

  ! checks if anything to do
  if( num_abs_boundary_faces > 0 ) then
    nglob_xy = nglob
  else
    nglob_xy = 1
  endif

  ! elastic domains
  if( ELASTIC_SIMULATION ) then
    allocate( rmassx(nglob_xy), rmassy(nglob_xy), rmassz(nglob_xy), stat=ier)
    if(ier /= 0) stop 'error in allocate 21'
    rmassx(:) = 0._CUSTOM_REAL
    rmassy(:) = 0._CUSTOM_REAL
    rmassz(:) = 0._CUSTOM_REAL
  endif

  ! acoustic domains
  if( ACOUSTIC_SIMULATION ) then
    allocate( rmassz_acoustic(nglob_xy), stat=ier)
    if(ier /= 0) stop 'error in allocate 22'
    rmassz_acoustic(:) = 0._CUSTOM_REAL
  endif

  ! use the non-dimensional time step to make the mass matrix correction
  if(CUSTOM_REAL == SIZE_REAL) then
    deltat = sngl(DT)
    deltatover2 = sngl(0.5d0*deltat)
  else
    deltat = DT
    deltatover2 = 0.5d0*deltat
  endif

  ! adds contributions to mass matrix to stabilize Stacey conditions
  do iface=1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    ! elastic elements
    if( ispec_is_elastic(ispec)) then
      ! reference gll points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

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

        ! assembles mass matrix on global points
        rmassx(iglob) = rmassx(iglob) + tx*jacobianw
        rmassy(iglob) = rmassy(iglob) + ty*jacobianw
        rmassz(iglob) = rmassz(iglob) + tz*jacobianw
      enddo
    endif ! elastic

    ! acoustic element
    if( ispec_is_acoustic(ispec) ) then

      ! reference gll points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets global index
        iglob=ibool(i,j,k,ispec)

        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! C * DT/2 contribution
        sn = deltatover2/rho_vp(i,j,k,ispec)

        rmassz_acoustic(iglob) = rmassz_acoustic(iglob) + jacobianw*sn
      enddo
    endif ! acoustic
  enddo

  end subroutine create_mass_matrices_Stacey

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_pml_elastic(nspec,ibool)

    use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,SIZE_REAL,CUSTOM_REAL,DT,&
                                      is_CPML,CPML_regions,d_store_x,d_store_y,d_store_z, &
                                      K_store_x,K_store_y,K_store_z,nspec_cpml,CPML_to_spec,&
                                      CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                                      CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ
    use create_regions_mesh_ext_par,only : rmass,rhostore,jacobianstore,wxgll,wygll,wzgll,ispec_is_elastic

    implicit none

    integer, intent(in) :: nspec
    integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool

    ! local parameters
    double precision :: weight
    real(kind=CUSTOM_REAL) :: jacobianl,deltat,mat_coef
    integer :: ispec,iglob,i,j,k,ispec_CPML

    ! use the non-dimensional time step to make the mass matrix correction
    if(CUSTOM_REAL == SIZE_REAL) then
       deltat = sngl(DT)
    else
       deltat = DT
    endif

    ! loops over physical mesh elements
    do ispec=1,nspec
       if( .not. is_CPML(ispec) .and. ispec_is_elastic(ispec) ) then
          do k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX
                   ! defines the material coefficient associated to the domain
                   mat_coef = rhostore(i,j,k,ispec)
                   iglob = ibool(i,j,k,ispec)
                   weight = wxgll(i)*wygll(j)*wzgll(k)
                   jacobianl = jacobianstore(i,j,k,ispec)

                   if(CUSTOM_REAL == SIZE_REAL) then
                      rmass(iglob) = rmass(iglob) + &
                           sngl( dble(jacobianl) * weight * dble(mat_coef) )
                   else
                      rmass(iglob) = rmass(iglob) + &
                           jacobianl * weight * mat_coef
                   endif
                enddo
             enddo
          enddo
       endif
    enddo

    ! loops over C-PML elements
    do ispec_CPML=1,nspec_cpml
       ispec = CPML_to_spec(ispec_CPML)

       if( is_CPML(ispec) .and. ispec_is_elastic(ispec) ) then
          ! X_surface C-PML
          if( CPML_regions(ispec_CPML) == CPML_X_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = rhostore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass(iglob) = rmass(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_x(i,j,k,ispec_CPML)) + dble(d_store_x(i,j,k,ispec_CPML)) * deltat / 2.d0) )
                      else
                         rmass(iglob) = rmass(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_x(i,j,k,ispec_CPML) + d_store_x(i,j,k,ispec_CPML) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! Y_surface C-PML
          else if( CPML_regions(ispec_CPML) == CPML_Y_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = rhostore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass(iglob) = rmass(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_y(i,j,k,ispec_CPML)) + dble(d_store_y(i,j,k,ispec_CPML)) * deltat / 2.d0) )
                      else
                         rmass(iglob) = rmass(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_y(i,j,k,ispec_CPML) + d_store_y(i,j,k,ispec_CPML) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! Z_surface C-PML
          else if( CPML_regions(ispec_CPML) == CPML_Z_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = rhostore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass(iglob) = rmass(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_z(i,j,k,ispec_CPML)) + dble(d_store_z(i,j,k,ispec_CPML)) * deltat / 2.d0) )
                      else
                         rmass(iglob) = rmass(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_z(i,j,k,ispec_CPML) + d_store_z(i,j,k,ispec_CPML) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! XY_edge C-PML
          else if( CPML_regions(ispec_CPML) == CPML_XY_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = rhostore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass(iglob) = rmass(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) + &
                              (dble(d_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) + &
                              dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML))) * deltat / 2.d0) )
                      else
                         rmass(iglob) = rmass(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_x(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML) + &
                              (d_store_x(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML) + &
                              d_store_y(i,j,k,ispec_CPML) * K_store_x(i,j,k,ispec_CPML)) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! XZ_edge C-PML
          else if( CPML_regions(ispec_CPML) == CPML_XZ_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = rhostore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass(iglob) = rmass(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                              (dble(d_store_x(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                              dble(d_store_z(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML))) * deltat / 2.d0) )
                      else
                         rmass(iglob) = rmass(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_x(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              (d_store_x(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              d_store_z(i,j,k,ispec_CPML) * K_store_x(i,j,k,ispec_CPML)) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! YZ_edge C-PML
          else if( CPML_regions(ispec_CPML) == CPML_YZ_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = rhostore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass(iglob) = rmass(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                              (dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                              dble(d_store_z(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML))) * deltat / 2.d0) )
                      else
                         rmass(iglob) = rmass(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_y(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              (d_store_y(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              d_store_z(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML)) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! XYZ_corner C-PML
          else if( CPML_regions(ispec_CPML) == CPML_XYZ ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = rhostore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass(iglob) = rmass(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) * &
                              dble(K_store_z(i,j,k,ispec_CPML)) + (dble(d_store_x(i,j,k,ispec_CPML)) * &
                              dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                              dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML)) * &
                              dble(K_store_z(i,j,k,ispec_CPML)) + dble(d_store_z(i,j,k,ispec_CPML)) * &
                              dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML))) * &
                              deltat / 2.d0) )
                      else
                         rmass(iglob) = rmass(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_x(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              (d_store_x(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              d_store_y(i,j,k,ispec_CPML) * K_store_x(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              d_store_z(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML)) * &
                              deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo
          else
             stop 'error in PML mesh file'
          endif
       endif
    enddo ! do ispec_CPML=1,nspec_cpml

  end subroutine create_mass_matrices_pml_elastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_pml_acoustic(nspec,ibool)

    use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,SIZE_REAL,CUSTOM_REAL,DT,&
                                      is_CPML,CPML_regions,d_store_x,d_store_y,d_store_z, &
                                      K_store_x,K_store_y,K_store_z,nspec_cpml,CPML_to_spec,&
                                      CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                                      CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ
    use create_regions_mesh_ext_par,only : rmass_acoustic,kappastore,jacobianstore,wxgll,wygll,wzgll,ispec_is_acoustic

    implicit none

    integer, intent(in) :: nspec
    integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool

    ! local parameters
    double precision :: weight
    real(kind=CUSTOM_REAL) :: jacobianl,deltat,mat_coef
    integer :: ispec,iglob,i,j,k,ispec_CPML

    ! use the non-dimensional time step to make the mass matrix correction
    if(CUSTOM_REAL == SIZE_REAL) then
       deltat = sngl(DT)
    else
       deltat = DT
    endif

    ! loops over physical mesh elements
    do ispec=1,nspec
       if( .not. is_CPML(ispec) .and. ispec_is_acoustic(ispec) ) then
          do k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX
                   ! defines the material coefficient associated to the domain
                   mat_coef = 1.d0 / kappastore(i,j,k,ispec)
                   iglob = ibool(i,j,k,ispec)
                   weight = wxgll(i)*wygll(j)*wzgll(k)
                   jacobianl = jacobianstore(i,j,k,ispec)

                   if(CUSTOM_REAL == SIZE_REAL) then
                      rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                           sngl( dble(jacobianl) * weight * dble(mat_coef) )
                   else
                      rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                           jacobianl * weight * mat_coef
                   endif
                enddo
             enddo
          enddo
       endif
    enddo

    ! loops over C-PML elements
    do ispec_CPML=1,nspec_cpml
       ispec = CPML_to_spec(ispec_CPML)

       if( is_CPML(ispec) .and. ispec_is_acoustic(ispec) ) then
          ! X_surface C-PML
          if( CPML_regions(ispec_CPML) == CPML_X_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = 1.d0 / kappastore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_x(i,j,k,ispec_CPML)) + dble(d_store_x(i,j,k,ispec_CPML)) * deltat / 2.d0) )
                      else
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_x(i,j,k,ispec_CPML) + d_store_x(i,j,k,ispec_CPML) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! Y_surface C-PML
          else if( CPML_regions(ispec_CPML) == CPML_Y_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = 1.d0 / kappastore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_y(i,j,k,ispec_CPML)) + dble(d_store_y(i,j,k,ispec_CPML)) * deltat / 2.d0) )
                      else
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_y(i,j,k,ispec_CPML) + d_store_y(i,j,k,ispec_CPML) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! Z_surface C-PML
          else if( CPML_regions(ispec_CPML) == CPML_Z_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = 1.d0 / kappastore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_z(i,j,k,ispec_CPML)) + dble(d_store_z(i,j,k,ispec_CPML)) * deltat / 2.d0) )
                      else
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_z(i,j,k,ispec_CPML) + d_store_z(i,j,k,ispec_CPML) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! XY_edge C-PML
          else if( CPML_regions(ispec_CPML) == CPML_XY_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = 1.d0 / kappastore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) + &
                              (dble(d_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) + &
                              dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML))) * deltat / 2.d0) )
                      else
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_x(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML) + &
                              (d_store_x(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML) + &
                              d_store_y(i,j,k,ispec_CPML) * K_store_x(i,j,k,ispec_CPML)) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! XZ_edge C-PML
          else if( CPML_regions(ispec_CPML) == CPML_XZ_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = 1.d0 / kappastore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                              (dble(d_store_x(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                              dble(d_store_z(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML))) * deltat / 2.d0) )
                      else
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_x(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              (d_store_x(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              d_store_z(i,j,k,ispec_CPML) * K_store_x(i,j,k,ispec_CPML)) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! YZ_edge C-PML
          else if( CPML_regions(ispec_CPML) == CPML_YZ_ONLY ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = 1.d0 / kappastore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                              (dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                              dble(d_store_z(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML))) * deltat / 2.d0) )
                      else
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_y(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              (d_store_y(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              d_store_z(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML)) * deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo

             ! XYZ_corner C-PML
          else if( CPML_regions(ispec_CPML) == CPML_XYZ ) then
             do k=1,NGLLZ
                do j=1,NGLLY
                   do i=1,NGLLX
                      ! defines the material coefficient associated to the domain
                      mat_coef = 1.d0 / kappastore(i,j,k,ispec)
                      iglob = ibool(i,j,k,ispec)
                      weight = wxgll(i)*wygll(j)*wzgll(k)
                      jacobianl = jacobianstore(i,j,k,ispec)

                      if(CUSTOM_REAL == SIZE_REAL) then
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              sngl( dble(jacobianl) * weight * dble(mat_coef) * &
                              (dble(K_store_x(i,j,k,ispec_CPML)) * dble(K_store_y(i,j,k,ispec_CPML)) * &
                              dble(K_store_z(i,j,k,ispec_CPML)) + (dble(d_store_x(i,j,k,ispec_CPML)) * &
                              dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML)) + &
                              dble(d_store_y(i,j,k,ispec_CPML)) * dble(K_store_x(i,j,k,ispec_CPML)) * &
                              dble(K_store_z(i,j,k,ispec_CPML)) + dble(d_store_z(i,j,k,ispec_CPML)) * &
                              dble(K_store_y(i,j,k,ispec_CPML)) * dble(K_store_z(i,j,k,ispec_CPML))) * &
                              deltat / 2.d0) )
                      else
                         rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                              jacobianl * weight * mat_coef * &
                              (K_store_x(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              (d_store_x(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              d_store_y(i,j,k,ispec_CPML) * K_store_x(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML) + &
                              d_store_z(i,j,k,ispec_CPML) * K_store_y(i,j,k,ispec_CPML) * K_store_z(i,j,k,ispec_CPML)) * &
                              deltat / 2.d0)
                      endif
                   enddo
                enddo
             enddo
          else
             stop 'error in PML mesh file'
          endif
       endif
    enddo ! do ispec_CPML=1,nspec_cpml

  end subroutine create_mass_matrices_pml_acoustic
