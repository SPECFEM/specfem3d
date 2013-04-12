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
!
! United States and French Government Sponsorship Acknowledged.

subroutine pml_set_local_dampingcoeff(myrank,xstore,ystore,zstore)

  ! calculates damping profiles and auxiliary coefficients on C-PML points

  use generate_databases_par, only: ibool,NGLOB_AB,d_store_x,d_store_y,d_store_z, &
                                    K_store_x,K_store_y,K_store_z,alpha_store,CPML_to_spec, &
                                    CPML_width_x,CPML_width_y,CPML_width_z,NPOWER,K_MAX_PML, &
                                    CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,nspec_cpml,PML_INSTEAD_OF_FREE_SURFACE, &
                                    IMAIN,FOUR_THIRDS,CPML_REGIONS,f0_FOR_PML,PI, &
                                    CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  use create_regions_mesh_ext_par, only: kappastore,mustore,rhostore,rho_vp,ispec_is_acoustic,ispec_is_elastic, &
                                         ELASTIC_SIMULATION, ACOUSTIC_SIMULATION

  implicit none

  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: xstore,ystore,zstore

  ! local parameters
  integer :: i,j,k,ispec,iglob,ispec_CPML,ier

  ! JC JC: Remove the parameter definition here and make the calculation of ALPHA_MAX_PML automatic
  !        by recovering the value of hdur in FORCESOLUTION/CMTSOLUTION
  real(kind=CUSTOM_REAL) :: ALPHA_MAX_PML

  real(kind=CUSTOM_REAL) :: pml_damping_profile_l,dist,vp
  real(kind=CUSTOM_REAL) :: xoriginleft,xoriginright,yoriginfront,yoriginback,zoriginbottom,zorigintop
  real(kind=CUSTOM_REAL) :: abscissa_in_PML_x,abscissa_in_PML_y,abscissa_in_PML_z
  real(kind=CUSTOM_REAL) :: d_x,d_y,d_z,k_x,k_y,k_z,alpha_x,alpha_y,alpha_z
  real(kind=CUSTOM_REAL) :: x_min,x_min_all,y_min,y_min_all,z_min,z_min_all,&
                            x_max,x_max_all,y_max,y_max_all,z_max,z_max_all,&
                            x_origin,y_origin,z_origin
  real(kind=CUSTOM_REAL) :: CPML_width_x_left, CPML_width_x_right,&
                            CPML_width_y_front,CPML_width_y_back,&
                            CPML_width_z_top,CPML_width_z_bottom,&
                            CPML_x_left, CPML_x_right,&
                            CPML_y_front,CPML_y_back,&
                            CPML_z_top,CPML_z_bottom,&
                            CPML_width_x_left_max_all, CPML_width_x_right_max_all,&
                            CPML_width_y_front_max_all,CPML_width_y_back_max_all,&
                            CPML_width_z_top_max_all,CPML_width_z_bottom_max_all

  ! stores damping profiles
  allocate(d_store_x(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if(ier /= 0) stop 'error allocating array d_store_x'
  allocate(d_store_y(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if(ier /= 0) stop 'error allocating array d_store_y'
  allocate(d_store_z(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if(ier /= 0) stop 'error allocating array d_store_z'

  ! stores auxiliary coefficients
  allocate(K_store_x(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if(ier /= 0) stop 'error allocating array K_store_x'
  allocate(K_store_y(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if(ier /= 0) stop 'error allocating array K_store_y'
  allocate(K_store_z(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if(ier /= 0) stop 'error allocating array K_store_z'
  allocate(alpha_store(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if(ier /= 0) stop 'error allocating array alpha_store'

  d_store_x = 0._CUSTOM_REAL
  d_store_y = 0._CUSTOM_REAL
  d_store_z = 0._CUSTOM_REAL

  K_store_x = 0._CUSTOM_REAL
  K_store_y = 0._CUSTOM_REAL
  K_store_z = 0._CUSTOM_REAL

  alpha_store = 0._CUSTOM_REAL

  if(ELASTIC_SIMULATION .and. .not. ACOUSTIC_SIMULATION)then
    ALPHA_MAX_PML = PI*f0_FOR_PML*1.0  ! from Festa and Vilotte (2005)
  elseif(ACOUSTIC_SIMULATION .and. .not. ELASTIC_SIMULATION)then
    ALPHA_MAX_PML = PI*f0_FOR_PML*2.0
  elseif(ELASTIC_SIMULATION .and. ACOUSTIC_SIMULATION)then
    ALPHA_MAX_PML = PI*f0_FOR_PML*2.0
  else
    stop 'Currently PML is not implemented for POROELASTIC_SIMULATION'
  endif

  x_min = minval(xstore(:))
  x_max = maxval(xstore(:))
  y_min = minval(ystore(:))
  y_max = maxval(ystore(:))
  z_min = minval(zstore(:))
  z_max = maxval(zstore(:))

  x_min_all = 10.e30
  x_max_all = -10.e30
  y_min_all = 10.e30
  y_max_all = -10.e30
  z_min_all = 10.e30
  z_max_all = -10.e30

  call min_all_all_cr(x_min,x_min_all)
  call min_all_all_cr(y_min,y_min_all)
  call min_all_all_cr(z_min,z_min_all)

  call max_all_all_cr(x_max,x_max_all)
  call max_all_all_cr(y_max,y_max_all)
  call max_all_all_cr(z_max,z_max_all)

  x_origin = (x_min_all + x_max_all)/2.0
  y_origin = (y_min_all + y_max_all)/2.0
  z_origin = (z_max_all + z_min_all)/2.0

!Calculation of CPML_width_x,CPML_width_y,CPML_width_Z
!we assume CPML_width_x,CPML_width_y,CPML_width_Z are constants inside PML layer

  CPML_width_x_left = 0.0
  CPML_width_x_right = 0.0
  CPML_width_y_front = 0.0
  CPML_width_y_back = 0.0
  CPML_width_z_top = 0.0
  CPML_width_z_bottom = 0.0

  CPML_x_right = x_max_all
  CPML_x_left = x_min_all
  CPML_y_front = y_max_all
  CPML_y_back = y_min_all
  CPML_z_top = z_max_all
  CPML_z_bottom = z_min_all

  do ispec_CPML=1,nspec_cpml
     ispec = CPML_to_spec(ispec_CPML)
     do k=1,NGLLZ
        do j=1,NGLLY
           do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            if( CPML_regions(ispec_CPML) == CPML_X_ONLY ) then
              if(xstore(iglob) - x_origin > 0.d0)then
                if(xstore(iglob) - x_origin <= CPML_x_right - x_origin )then
                   CPML_x_right = xstore(iglob)
                endif
              else
                if(abs(xstore(iglob) - x_origin) <= abs(CPML_x_left-x_origin))then
                   CPML_x_left = xstore(iglob)
                endif
              endif
             endif

            if( CPML_regions(ispec_CPML) == CPML_Y_ONLY ) then
              if(ystore(iglob) - y_origin > 0.d0)then
                if(ystore(iglob) - y_origin <= CPML_y_front - y_origin )then
                   CPML_y_front = ystore(iglob)
                endif
              else
                if(abs(ystore(iglob) - y_origin) <= abs(CPML_y_back-y_origin))then
                   CPML_y_back = ystore(iglob)
                endif
              endif
             endif

            if( CPML_regions(ispec_CPML) == CPML_Z_ONLY ) then
              if(zstore(iglob) - z_origin > 0.d0)then
                if(zstore(iglob) - z_origin <= CPML_z_top - z_origin )then
                   CPML_z_top = zstore(iglob)
                endif
              else
                if(abs(zstore(iglob) - z_origin) <= abs(CPML_z_bottom-z_origin))then
                   CPML_z_bottom = zstore(iglob)
                endif
              endif
             endif

            enddo
         enddo
      enddo
  enddo

  CPML_width_x_right = x_max_all - CPML_x_right
  CPML_width_x_left = CPML_x_left - x_min_all
  CPML_width_y_front = y_max_all - CPML_y_front
  CPML_width_y_back = CPML_y_back - y_min_all
  CPML_width_z_top = z_max_all - CPML_z_top
  CPML_width_z_bottom = CPML_z_bottom - z_min_all

  call max_all_all_cr(CPML_width_x_left,CPML_width_x_left_max_all)
  call max_all_all_cr(CPML_width_x_right,CPML_width_x_right_max_all)
  call max_all_all_cr(CPML_width_y_front,CPML_width_y_front_max_all)
  call max_all_all_cr(CPML_width_y_back,CPML_width_y_back_max_all)
  call max_all_all_cr(CPML_width_z_top,CPML_width_z_top_max_all)
  call max_all_all_cr(CPML_width_z_bottom,CPML_width_z_bottom_max_all)

  xoriginleft   = x_min_all + CPML_width_x_left_max_all
  xoriginright  = x_max_all - CPML_width_x_right_max_all
  yoriginback   = y_min_all + CPML_width_y_back_max_all
  yoriginfront  = y_max_all - CPML_width_y_front_max_all
  zoriginbottom = z_min_all + CPML_width_z_bottom_max_all

  CPML_width_x = max(CPML_width_x_left_max_all,CPML_width_x_right_max_all)
  CPML_width_y = max(CPML_width_y_front_max_all,CPML_width_y_back_max_all)
  CPML_width_z = max(CPML_width_z_bottom_max_all,CPML_width_z_top_max_all )

  if( PML_INSTEAD_OF_FREE_SURFACE ) then
     zorigintop = z_max_all - CPML_width_z_top_max_all
  endif

  ! user output
  if( myrank == 0 ) then
     write(IMAIN,*)
     write(IMAIN,*) 'Boundary values of X-/Y-/Z-regions'
     write(IMAIN,*) minval(xstore(:)), maxval(xstore(:))
     write(IMAIN,*) minval(ystore(:)), maxval(ystore(:))
     write(IMAIN,*) minval(zstore(:)), maxval(zstore(:))
     write(IMAIN,*)
     write(IMAIN,*) 'Origins of right/left X-surface C-PML',xoriginright,xoriginleft
     write(IMAIN,*) 'Origins of front/back Y-surface C-PML',yoriginfront,yoriginback
     write(IMAIN,*) 'Origin of bottom Z-surface C-PML',zoriginbottom
     if( PML_INSTEAD_OF_FREE_SURFACE ) then
        write(IMAIN,*) 'Origin of top Z-surface C-PML',zorigintop
     end if
     write(IMAIN,*)
     write(IMAIN,*) 'CPML_width_x: ',CPML_width_x
     write(IMAIN,*) 'CPML_width_y: ',CPML_width_y
     write(IMAIN,*) 'CPML_width_z: ',CPML_width_z
     write(IMAIN,*)
  endif
  call sync_all()

  ! loops over all C-PML elements
  do ispec_CPML=1,nspec_cpml
     ispec = CPML_to_spec(ispec_CPML)

     do k=1,NGLLZ
        do j=1,NGLLY
           do i=1,NGLLX
              ! calculates P-velocity
              if( ispec_is_acoustic(ispec) ) then
                 vp = sqrt( kappastore(i,j,k,ispec)/rhostore(i,j,k,ispec) )
              else if( ispec_is_elastic(ispec) ) then
                 vp = (FOUR_THIRDS * mustore(i,j,k,ispec) + kappastore(i,j,k,ispec)) / rho_vp(i,j,k,ispec)
              else
                 print*,'element index',ispec
                 print*,'C-PML element index ',ispec_CPML
                 call exit_mpi(myrank,'C-PML error: element has an unvalid P-velocity')
              endif

              iglob = ibool(i,j,k,ispec)

              if( CPML_regions(ispec_CPML) == CPML_X_ONLY ) then
                 !------------------------------------------------------------------------------
                 !---------------------------- X-surface C-PML ---------------------------------
                 !------------------------------------------------------------------------------

                 if( xstore(iglob) - x_origin > 0.d0 ) then

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xstore(iglob) - xoriginright

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML element's GLL point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                 elseif( xstore(iglob) - x_origin < 0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xoriginleft - xstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML grid point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                 else
                    stop "there is error in mesh of CPML-layer x"
                 endif

                 !! DK DK define an alias for y and z variable names (which are the same)
                 !  stores damping profiles and auxiliary coefficients at the C-PML element's GLL points
                 K_store_x(i,j,k,ispec_CPML) = K_x
                 d_store_x(i,j,k,ispec_CPML) = d_x

                 K_store_y(i,j,k,ispec_CPML) = 1.d0
                 d_store_y(i,j,k,ispec_CPML) = 0.d0

                 K_store_z(i,j,k,ispec_CPML) = 1.d0
                 d_store_z(i,j,k,ispec_CPML) = 0.d0

                 alpha_store(i,j,k,ispec_CPML) = alpha_x

              elseif( CPML_regions(ispec_CPML) == CPML_Y_ONLY ) then
                 !------------------------------------------------------------------------------
                 !---------------------------- Y-surface C-PML ---------------------------------
                 !------------------------------------------------------------------------------

                 if( ystore(iglob) - y_origin > 0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = ystore(iglob) - yoriginfront

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                 elseif( ystore(iglob) - y_origin < 0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = yoriginback - ystore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                 else
                    stop "there is error in mesh of  CPML-layer y"

                 endif

                 !! DK DK define an alias for y and z variable names (which are the same)
                 !  stores damping profiles and auxiliary coefficients at the C-PML element's GLL points
                 K_store_x(i,j,k,ispec_CPML) = 1.d0
                 d_store_x(i,j,k,ispec_CPML) = 0.d0

                 K_store_y(i,j,k,ispec_CPML) = K_y
                 d_store_y(i,j,k,ispec_CPML) = d_y

                 K_store_z(i,j,k,ispec_CPML) = 1.d0
                 d_store_z(i,j,k,ispec_CPML) = 0.d0

                 alpha_store(i,j,k,ispec_CPML) = alpha_y

              elseif( CPML_regions(ispec_CPML) == CPML_Z_ONLY ) then
                 !------------------------------------------------------------------------------
                 !---------------------------- Z-surface C-PML ---------------------------------
                 !------------------------------------------------------------------------------

                 if( zstore(iglob) - z_origin > 0.d0 ) then
                    if( PML_INSTEAD_OF_FREE_SURFACE ) then
                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_z = zstore(iglob) - zorigintop

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_z / CPML_width_z

                       ! gets damping profile at the C-PML element's GLL point
                       d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                       alpha_z = ALPHA_MAX_PML / 2.d0
                       K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       if( d_z < 0.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                       if( K_z < 1.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif
                    endif
                 elseif( zstore(iglob) - z_origin < 0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_z = zoriginbottom - zstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_z / CPML_width_z

                    ! gets damping profile at the C-PML element's GLL point
                    d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                    alpha_z = ALPHA_MAX_PML / 2.d0
                    K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_z < 0.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                    if( K_z < 1.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                 else
                    stop "there is error in mesh of CPML-layer z"
                 endif

                 !! DK DK define an alias for y and z variable names (which are the same)
                 !  stores damping profiles and auxiliary coefficients at the C-PML element's GLL points
                 K_store_x(i,j,k,ispec_CPML) = 1.d0
                 d_store_x(i,j,k,ispec_CPML) = 0.d0

                 K_store_y(i,j,k,ispec_CPML) = 1.d0
                 d_store_y(i,j,k,ispec_CPML) = 0.d0

                 K_store_z(i,j,k,ispec_CPML) = K_z
                 d_store_z(i,j,k,ispec_CPML) = d_z

                 alpha_store(i,j,k,ispec_CPML) = alpha_z

              elseif( CPML_regions(ispec_CPML) == CPML_XY_ONLY ) then
                 !------------------------------------------------------------------------------
                 !---------------------------- XY-edge C-PML -----------------------------------
                 !------------------------------------------------------------------------------

                 if( xstore(iglob) - x_origin>0.d0 .and. ystore(iglob) - y_origin>0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xstore(iglob) - xoriginright

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML element's GLL point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = ystore(iglob) - yoriginfront

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                 elseif( xstore(iglob) - x_origin>0.d0 .and. ystore(iglob) - y_origin<0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xstore(iglob) - xoriginright

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML element's GLL point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = yoriginback - ystore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                 elseif( xstore(iglob) - x_origin<0.d0 .and. ystore(iglob) - y_origin>0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xoriginleft - xstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML element's GLL point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = ystore(iglob) - yoriginfront

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                 elseif( xstore(iglob)  - x_origin <0.d0 .and. ystore(iglob)  - y_origin <0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xoriginleft - xstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML element's GLL point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = yoriginback - ystore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                 else
                    stop "there is error in mesh of CPML-layer xy"
                 endif

                 !! DK DK define an alias for y and z variable names (which are the same)
                 !  stores damping profiles and auxiliary coefficients at the C-PML element's GLL points
                 K_store_x(i,j,k,ispec_CPML) = K_x
                 d_store_x(i,j,k,ispec_CPML) = d_x

                 K_store_y(i,j,k,ispec_CPML) = K_y
                 d_store_y(i,j,k,ispec_CPML) = d_y


                 K_store_z(i,j,k,ispec_CPML) = 1.d0
                 d_store_z(i,j,k,ispec_CPML) = 0.d0

                 alpha_store(i,j,k,ispec_CPML) = ALPHA_MAX_PML / 2.d0

              elseif( CPML_regions(ispec_CPML) == CPML_XZ_ONLY ) then
                 !------------------------------------------------------------------------------
                 !---------------------------- XZ-edge C-PML -----------------------------------
                 !------------------------------------------------------------------------------

                 if( xstore(iglob) - x_origin>0.d0 .and. zstore(iglob) - z_origin>0.d0 ) then
                    if( PML_INSTEAD_OF_FREE_SURFACE ) then
                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_x = xstore(iglob) - xoriginright

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_x / CPML_width_x

                       ! gets damping profile at the C-PML element's GLL point
                       d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                       alpha_x = ALPHA_MAX_PML / 2.d0
                       K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_z = zstore(iglob) - zorigintop

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_z / CPML_width_z

                       ! gets damping profile at the C-PML element's GLL point
                       d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                       alpha_z = ALPHA_MAX_PML / 2.d0
                       K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       if( d_x < 0.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( K_x < 1.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( d_z < 0.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                       if( K_z < 1.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                    endif

                 elseif( xstore(iglob) - x_origin>0.d0 .and. zstore(iglob) - z_origin<0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xstore(iglob) - xoriginright

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML element's GLL point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_z = zoriginbottom - zstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_z / CPML_width_z

                    ! gets damping profile at the C-PML element's GLL point
                    d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                    alpha_z = ALPHA_MAX_PML / 2.d0
                    K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( d_z < 0.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                    if( K_z < 1.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                 elseif( xstore(iglob) - x_origin<0.d0 .and. zstore(iglob) - z_origin>0.d0 ) then
                    if( PML_INSTEAD_OF_FREE_SURFACE ) then
                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_x = xoriginleft - xstore(iglob)

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_x / CPML_width_x

                       ! gets damping profile at the C-PML element's GLL point
                       d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                       alpha_x = ALPHA_MAX_PML / 2.d0
                       K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_z = zstore(iglob) - zorigintop

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_z / CPML_width_z

                       ! gets damping profile at the C-PML element's GLL point
                       d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                       alpha_z = ALPHA_MAX_PML / 2.d0
                       K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       if( d_x < 0.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( K_x < 1.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( d_z < 0.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                       if( K_z < 1.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                    endif

                 elseif( xstore(iglob) - x_origin<0.d0 .and. zstore(iglob) - z_origin<0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xoriginleft - xstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML element's GLL point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_z = zoriginbottom - zstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_z / CPML_width_z

                    ! gets damping profile at the C-PML element's GLL point
                    d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                    alpha_z = ALPHA_MAX_PML / 2.d0
                    K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( d_z < 0.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                    if( K_z < 1.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                 else
                    stop "there is error in mesh of CPML-layer xz"
                 endif

                 !! DK DK define an alias for y and z variable names (which are the same)
                 !  stores damping profiles and auxiliary coefficients at the C-PML element's GLL points
                 K_store_x(i,j,k,ispec_CPML) = K_x
                 d_store_x(i,j,k,ispec_CPML) = d_x

                 K_store_y(i,j,k,ispec_CPML) = 1.d0
                 d_store_y(i,j,k,ispec_CPML) = 0.d0

                 K_store_z(i,j,k,ispec_CPML) = K_z
                 d_store_z(i,j,k,ispec_CPML) = d_z

                 alpha_store(i,j,k,ispec_CPML) = ALPHA_MAX_PML / 2.d0

              elseif( CPML_regions(ispec_CPML) == CPML_YZ_ONLY ) then
                 !------------------------------------------------------------------------------
                 !---------------------------- YZ-edge C-PML -----------------------------------
                 !------------------------------------------------------------------------------

                 if( ystore(iglob) - y_origin>0.d0 .and. zstore(iglob) - z_origin>0.d0 ) then
                    if( PML_INSTEAD_OF_FREE_SURFACE ) then
                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_y = ystore(iglob) - yoriginfront

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_y / CPML_width_y

                       ! gets damping profile at the C-PML element's GLL point
                       d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                       alpha_y = ALPHA_MAX_PML / 2.d0
                       K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_z = zstore(iglob) - zorigintop

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_z / CPML_width_z

                       ! gets damping profile at the C-PML element's GLL point
                       d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                       alpha_z = ALPHA_MAX_PML / 2.d0
                       K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       if( d_y < 0.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( K_y < 1.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif
 
                       if( d_z < 0.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                       if( K_z < 1.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif
                    endif

                 elseif( ystore(iglob) - y_origin>0.d0 .and. zstore(iglob) - z_origin<0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = ystore(iglob) - yoriginfront

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_z = zoriginbottom - zstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_z / CPML_width_z

                    ! gets damping profile at the C-PML element's GLL point
                    d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                    alpha_z = ALPHA_MAX_PML / 2.d0
                    K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( d_z < 0.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                    if( K_z < 1.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                 elseif( ystore(iglob) - y_origin<0.d0 .and. zstore(iglob) - z_origin>0.d0 ) then
                    if( PML_INSTEAD_OF_FREE_SURFACE ) then
                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_y = yoriginback - ystore(iglob)

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_y / CPML_width_y

                       ! gets damping profile at the C-PML element's GLL point
                       d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                       alpha_y = ALPHA_MAX_PML / 2.d0
                       K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_z = zstore(iglob) - zorigintop

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_z / CPML_width_z

                       ! gets damping profile at the C-PML element's GLL point
                       d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                       alpha_z = ALPHA_MAX_PML / 2.d0
                       K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       if( d_y < 0.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( K_y < 1.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( d_z < 0.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                       if( K_z < 1.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                    endif

                 elseif( ystore(iglob) - y_origin<0.d0 .and. zstore(iglob) - z_origin<0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = yoriginback - ystore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_z = zoriginbottom - zstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_z / CPML_width_z

                    ! gets damping profile at the C-PML element's GLL point
                    d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                    alpha_z = ALPHA_MAX_PML / 2.d0
                    K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( d_z < 0.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                    if( K_z < 1.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                 else
                    stop "there is error in mesh of CPML-layer yz"
                 endif

                 !! DK DK define an alias for y and z variable names (which are the same)
                 K_store_x(i,j,k,ispec_CPML) = 1.d0
                 d_store_x(i,j,k,ispec_CPML) = 0.d0

                 K_store_y(i,j,k,ispec_CPML) = K_y
                 d_store_y(i,j,k,ispec_CPML) = d_y

                 K_store_z(i,j,k,ispec_CPML) = K_z
                 d_store_z(i,j,k,ispec_CPML) = d_z

                 alpha_store(i,j,k,ispec_CPML) = ALPHA_MAX_PML / 2.d0

              elseif( CPML_regions(ispec_CPML) == CPML_XYZ ) then
                 !------------------------------------------------------------------------------
                 !---------------------------- XYZ-corner C-PML --------------------------------
                 !------------------------------------------------------------------------------

                 if( xstore(iglob) - x_origin>0.d0 .and. ystore(iglob) - y_origin>0.d0 .and. zstore(iglob) - z_origin>0.d0 ) then
                    if( PML_INSTEAD_OF_FREE_SURFACE ) then
                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_x = xstore(iglob) - xoriginright

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_x / CPML_width_x

                       ! gets damping profile at the C-PML grid point
                       d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                       alpha_x = ALPHA_MAX_PML / 2.d0
                       K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_y = ystore(iglob) - yoriginfront

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_y / CPML_width_y

                       ! gets damping profile at the C-PML element's GLL point
                       d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                       alpha_y = ALPHA_MAX_PML / 2.d0
                       K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_z = zstore(iglob) - zorigintop

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_z / CPML_width_z

                       ! gets damping profile at the C-PML element's GLL point
                       d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                       alpha_z = ALPHA_MAX_PML / 2.d0
                       K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       if( d_x < 0.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( K_x < 1.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( d_y < 0.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( K_y < 1.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( d_z < 0.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                       if( K_z < 1.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif
                    endif

                 elseif( xstore(iglob) - x_origin>0.d0 .and. ystore(iglob) - y_origin>0.d0 .and. zstore(iglob) - z_origin<0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xstore(iglob) - xoriginright

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML grid point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = ystore(iglob) - yoriginfront

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_z = zoriginbottom - zstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_z / CPML_width_z

                    ! gets damping profile at the C-PML element's GLL point
                    d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                    alpha_z = ALPHA_MAX_PML / 2.d0
                    K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( d_z < 0.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                    if( K_z < 1.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif


                 elseif( xstore(iglob) - x_origin>0.d0 .and. ystore(iglob) - y_origin<0.d0 .and. zstore(iglob) - z_origin>0.d0 ) then
                    if( PML_INSTEAD_OF_FREE_SURFACE ) then
                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_x = xstore(iglob) - xoriginright

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_x / CPML_width_x

                       ! gets damping profile at the C-PML element's GLL point
                       d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                       alpha_x = ALPHA_MAX_PML / 2.d0
                       K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_y = yoriginback - ystore(iglob)

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_y / CPML_width_y

                       ! gets damping profile at the C-PML element's GLL point
                       d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                       alpha_y = ALPHA_MAX_PML / 2.d0
                       K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER


                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_z = zstore(iglob) - zorigintop

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_z / CPML_width_z

                       ! gets damping profile at the C-PML element's GLL point
                       d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                       alpha_z = ALPHA_MAX_PML / 2.d0
                       K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       if( d_x < 0.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( K_x < 1.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( d_y < 0.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( K_y < 1.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( d_z < 0.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                       if( K_z < 1.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                    endif

                 elseif( xstore(iglob) - x_origin>0.d0 .and. ystore(iglob) - y_origin<0.d0 .and. zstore(iglob) - z_origin < 0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xstore(iglob) - xoriginright

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML element's GLL point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER


                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = yoriginback - ystore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER


                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_z = zoriginbottom - zstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_z / CPML_width_z

                    ! gets damping profile at the C-PML element's GLL point
                    d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                    alpha_z = ALPHA_MAX_PML / 2.d0
                    K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( d_z < 0.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                    if( K_z < 1.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                 elseif( xstore(iglob) - x_origin<0.d0 .and. ystore(iglob) - y_origin>0.d0 .and. zstore(iglob) - z_origin>0.d0 ) then
                    if( PML_INSTEAD_OF_FREE_SURFACE ) then
                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_x = xoriginleft - xstore(iglob)

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_x / CPML_width_x

                       ! gets damping profile at the C-PML grid point
                       d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                       alpha_x = ALPHA_MAX_PML / 2.d0
                       K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_y = ystore(iglob) - yoriginfront

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_y / CPML_width_y

                       ! gets damping profile at the C-PML element's GLL point
                       d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                       alpha_y = ALPHA_MAX_PML / 2.d0
                       K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_z = zstore(iglob) - zorigintop

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_z / CPML_width_z

                       ! gets damping profile at the C-PML element's GLL point
                       d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                       alpha_z = ALPHA_MAX_PML / 2.d0
                       K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       if( d_x < 0.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( K_x < 1.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( d_y < 0.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( K_y < 1.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( d_z < 0.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                       if( K_z < 1.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                    endif

                 elseif( xstore(iglob) - x_origin<0.d0 .and. ystore(iglob) - y_origin>0.d0 .and. zstore(iglob) - z_origin<0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xoriginleft - xstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML grid point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = ystore(iglob) - yoriginfront

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_z = zoriginbottom - zstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_z / CPML_width_z

                    ! gets damping profile at the C-PML element's GLL point
                    d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                    alpha_z = ALPHA_MAX_PML / 2.d0
                    K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( d_z < 0.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                    if( K_z < 1.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                 elseif( xstore(iglob) - x_origin<0.d0 .and. ystore(iglob) - y_origin<0.d0 .and. zstore(iglob) - z_origin>0.d0 ) then
                    if( PML_INSTEAD_OF_FREE_SURFACE ) then
                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_x = xoriginleft - xstore(iglob)

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_x / CPML_width_x

                       ! gets damping profile at the C-PML element's GLL point
                       d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                       alpha_x = ALPHA_MAX_PML / 2.d0
                       K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_y = yoriginback - ystore(iglob)

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_y / CPML_width_y

                       ! gets damping profile at the C-PML element's GLL point
                       d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                       alpha_y = ALPHA_MAX_PML / 2.d0
                       K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       ! gets abscissa of current grid point along the damping profile
                       abscissa_in_PML_z = zstore(iglob) - zorigintop

                       ! determines distance to C-PML/mesh interface
                       dist = abscissa_in_PML_z / CPML_width_z

                       ! gets damping profile at the C-PML element's GLL point
                       d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                       alpha_z = ALPHA_MAX_PML / 2.d0
                       K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                       if( d_x < 0.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( K_x < 1.d0 ) then
                          d_x = 0.d0
                          K_x = 1.d0
                       endif

                       if( d_y < 0.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( K_y < 1.d0 ) then
                          d_y = 0.d0
                          K_y = 1.d0
                       endif

                       if( d_z < 0.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                       if( K_z < 1.d0 ) then
                          d_z = 0.d0
                          K_z = 1.d0
                       endif

                    endif

                 elseif( xstore(iglob) - x_origin<0.d0 .and. ystore(iglob) - y_origin<0.d0 .and. zstore(iglob) - z_origin < 0.d0 ) then
                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_x = xoriginleft - xstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_x / CPML_width_x

                    ! gets damping profile at the C-PML element's GLL point
                    d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                    alpha_x = ALPHA_MAX_PML / 2.d0
                    K_x = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_y = yoriginback - ystore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_y / CPML_width_y

                    ! gets damping profile at the C-PML element's GLL point
                    d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                    alpha_y = ALPHA_MAX_PML / 2.d0
                    K_y = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    ! gets abscissa of current grid point along the damping profile
                    abscissa_in_PML_z = zoriginbottom - zstore(iglob)

                    ! determines distance to C-PML/mesh interface
                    dist = abscissa_in_PML_z / CPML_width_z

                    ! gets damping profile at the C-PML element's GLL point
                    d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                    alpha_z = ALPHA_MAX_PML / 2.d0
                    K_z = 1.d0 + (K_MAX_PML - 1.d0) * dist**NPOWER

                    if( d_x < 0.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( K_x < 1.d0 ) then
                       d_x = 0.d0
                       K_x = 1.d0
                    endif

                    if( d_y < 0.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( K_y < 1.d0 ) then
                       d_y = 0.d0
                       K_y = 1.d0
                    endif

                    if( d_z < 0.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                    if( K_z < 1.d0 ) then
                       d_z = 0.d0
                       K_z = 1.d0
                    endif

                 else
                    stop "there is error in mesh of CPML-layer xyz"
                 endif

                 !! DK DK define an alias for y and z variable names (which are the same)
                 K_store_x(i,j,k,ispec_CPML) = K_x
                 d_store_x(i,j,k,ispec_CPML) = d_x

                 K_store_y(i,j,k,ispec_CPML) = K_y
                 d_store_y(i,j,k,ispec_CPML) = d_y

                 K_store_z(i,j,k,ispec_CPML) = K_z
                 d_store_z(i,j,k,ispec_CPML) = d_z

                 alpha_store(i,j,k,ispec_CPML) = ALPHA_MAX_PML / 2.d0

              endif
           enddo
        enddo
     enddo
  enddo !ispec_CPML

end subroutine pml_set_local_dampingcoeff

!
!-------------------------------------------------------------------------------------------------
!

function pml_damping_profile_l(myrank,iglob,dist,vp,delta)

  ! defines d, the damping profile at the C-PML element's GLL point for a given:
  !   dist:  distance to C-PML/mesh interface
  !   vp:    P-velocity
  !   delta: thickness of the C-PML layer

  use generate_databases_par, only: CUSTOM_REAL,NPOWER,CPML_Rcoef,damping_factor

  implicit none

  integer, intent(in) :: myrank,iglob

  real(kind=CUSTOM_REAL), intent(in) :: dist,vp,delta

  real(kind=CUSTOM_REAL) :: pml_damping_profile_l

  ! gets damping profile
  if( NPOWER >= 1 ) then
     ! INRIA research report section 6.1:  http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
     pml_damping_profile_l = - ((NPOWER + 1) * vp * log(CPML_Rcoef) / (2.d0 * delta) * damping_factor) * dist**NPOWER
  else
     call exit_mpi(myrank,'C-PML error: NPOWER must be greater than or equal to 1')
  endif

  ! checks coordinates of C-PML points and thickness of C-PML layer
  if( delta < dist ) then
     print*,'C-PML point ',iglob
     print*,'distance to C-PML/mesh interface ',dist
     print*,'C-PML thickness ',delta
     call exit_mpi(myrank,'C-PML error: distance to C-PML/mesh interface is bigger than thickness of C-PML layer')
  endif

end function pml_damping_profile_l

