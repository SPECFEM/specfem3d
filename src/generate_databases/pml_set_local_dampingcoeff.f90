!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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
                                    K_store_x,K_store_y,K_store_z,alpha_store_x,alpha_store_y,alpha_store_z,CPML_to_spec, &
                                    CPML_width_x,CPML_width_y,CPML_width_z,NPOWER,&
                                    CUSTOM_REAL,SIZE_REAL,NGLLX,NGLLY,NGLLZ,nspec_cpml,PML_INSTEAD_OF_FREE_SURFACE, &
                                    IMAIN,CPML_REGIONS,f0_FOR_PML,PI, &
                                    CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ,&
                                    SIMULATION_TYPE,SAVE_FORWARD,nspec => NSPEC_AB,is_CPML,&
                                    mask_ibool_interior_domain,nglob_interface_PML_acoustic,points_interface_PML_acoustic,&
                                    nglob_interface_PML_elastic,points_interface_PML_elastic,&
                                    ZERO,ONE,TWO,HUGEVAL

  use create_regions_mesh_ext_par, only: rhostore,rho_vp,ispec_is_acoustic,ispec_is_elastic, &
                                         ELASTIC_SIMULATION, ACOUSTIC_SIMULATION

  implicit none

  integer, intent(in) :: myrank
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: xstore,ystore,zstore

  ! local parameters
  integer :: i,j,k,ispec,iglob,ispec_CPML,ier

  real(kind=CUSTOM_REAL) :: ALPHA_MAX_PML_x,ALPHA_MAX_PML_y,ALPHA_MAX_PML_z
  real(kind=CUSTOM_REAL), parameter :: K_MAX_PML = ONE,K_MIN_PML= ONE
  real(kind=CUSTOM_REAL), parameter :: const_for_tune_pml_damping_profile = 0.009_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: pml_damping_profile_l,dist,vp
  real(kind=CUSTOM_REAL) :: xoriginleft,xoriginright,yoriginfront,yoriginback,zoriginbottom,zorigintop
  real(kind=CUSTOM_REAL) :: abscissa_in_PML_x,abscissa_in_PML_y,abscissa_in_PML_z
  real(kind=CUSTOM_REAL) :: d_x,d_y,d_z,k_x,k_y,k_z,alpha_x,alpha_y,alpha_z,beta_x,beta_y,beta_z, &
                            meanval_1
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
                            CPML_width_z_top_max_all,CPML_width_z_bottom_max_all,&
                            vp_max,vp_max_all


  ! checks number of PML elements
  if (count(is_CPML(:)) /= NSPEC_CPML) then
    print *,'Error in slice ',myrank,': number of PML elements ',NSPEC_CPML,' but only ',count(is_CPML(:)),' flags set'
    stop 'Error C-PML array has invalid number of PML flags set'
  endif

  ! checks if C-PML flags assigned correctly
  do ispec_CPML = 1,NSPEC_CPML
    ispec = CPML_to_spec(ispec_CPML)
    if (.not. is_CPML(ispec)) stop 'Error found C-PML element with invalid PML flag'
  enddo

  ! stores damping profiles
  allocate(d_store_x(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) stop 'error allocating array d_store_x'
  allocate(d_store_y(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) stop 'error allocating array d_store_y'
  allocate(d_store_z(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) stop 'error allocating array d_store_z'

  ! stores auxiliary coefficients
  allocate(K_store_x(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) stop 'error allocating array K_store_x'
  allocate(K_store_y(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) stop 'error allocating array K_store_y'
  allocate(K_store_z(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) stop 'error allocating array K_store_z'
  allocate(alpha_store_x(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) stop 'error allocating array alpha_store_x'
  allocate(alpha_store_y(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) stop 'error allocating array alpha_store_y'
  allocate(alpha_store_z(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) stop 'error allocating array alpha_store_z'

  K_store_x(:,:,:,:) = ZERO
  K_store_y(:,:,:,:) = ZERO
  K_store_z(:,:,:,:) = ZERO

  d_store_x(:,:,:,:) = ZERO
  d_store_y(:,:,:,:) = ZERO
  d_store_z(:,:,:,:) = ZERO

  alpha_store_x(:,:,:,:) = ZERO
  alpha_store_y(:,:,:,:) = ZERO
  alpha_store_z(:,:,:,:) = ZERO

! from Festa and Vilotte (2005)
! Slight different from them is that we purposely set ALPHA_MAX_PML_x < ALPHA_MAX_PML_y < ALPHA_MAX_PML_z
! to facilitate the remove of singularities in PML damping parameters.
  ALPHA_MAX_PML_x = PI*f0_FOR_PML * 0.9_CUSTOM_REAL
  ALPHA_MAX_PML_y = PI*f0_FOR_PML * ONE
  ALPHA_MAX_PML_z = PI*f0_FOR_PML * 1.1_CUSTOM_REAL

! Assuming the computational domain is convex and can be approximatly seen as a box
! Calculation of origin of whole computational domain
  x_min = minval(xstore(:))
  x_max = maxval(xstore(:))
  y_min = minval(ystore(:))
  y_max = maxval(ystore(:))
  z_min = minval(zstore(:))
  z_max = maxval(zstore(:))

  if (CUSTOM_REAL == SIZE_REAL) then
    x_min_all = HUGEVAL
    y_min_all = HUGEVAL
    z_min_all = HUGEVAL

    x_max_all = - HUGEVAL
    y_max_all = - HUGEVAL
    z_max_all = - HUGEVAL
  else
    x_min_all = HUGEVAL
    y_min_all = HUGEVAL
    z_min_all = HUGEVAL

    x_max_all = - HUGEVAL
    y_max_all = - HUGEVAL
    z_max_all = - HUGEVAL
  endif

  call min_all_all_cr(x_min,x_min_all)
  call min_all_all_cr(y_min,y_min_all)
  call min_all_all_cr(z_min,z_min_all)

  call max_all_all_cr(x_max,x_max_all)
  call max_all_all_cr(y_max,y_max_all)
  call max_all_all_cr(z_max,z_max_all)

  x_origin = (x_min_all + x_max_all) / TWO
  y_origin = (y_min_all + y_max_all) / TWO
  z_origin = (z_max_all + z_min_all) / TWO

! Assuming CPML_width_x,CPML_width_y,CPML_width_Z are constants inside PML layer
! Calculation of width of PML along x, y and z direction, such as CPML_width_x,CPML_width_y,CPML_width_Z
  CPML_width_x_left   = ZERO
  CPML_width_x_right  = ZERO
  CPML_width_y_front  = ZERO
  CPML_width_y_back   = ZERO
  CPML_width_z_top    = ZERO
  CPML_width_z_bottom = ZERO

  CPML_x_right  = x_max_all
  CPML_x_left   = x_min_all
  CPML_y_front  = y_max_all
  CPML_y_back   = y_min_all
  CPML_z_top    = z_max_all
  CPML_z_bottom = z_min_all

  do ispec_CPML = 1,nspec_cpml

    ispec = CPML_to_spec(ispec_CPML)

    do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX

      iglob = ibool(i,j,k,ispec)

      if (CPML_regions(ispec_CPML) == CPML_X_ONLY  .or. CPML_regions(ispec_CPML) == CPML_XY_ONLY .or. &
          CPML_regions(ispec_CPML) == CPML_XZ_ONLY .or. CPML_regions(ispec_CPML) == CPML_XYZ) then
        if (xstore(iglob) - x_origin > ZERO) then
          if (xstore(iglob) - x_origin <= CPML_x_right - x_origin) then
            CPML_x_right = xstore(iglob)
          endif
        else
          if (abs(xstore(iglob) - x_origin) <= abs(CPML_x_left-x_origin)) then
            CPML_x_left = xstore(iglob)
          endif
        endif
      endif

      if (CPML_regions(ispec_CPML) == CPML_Y_ONLY  .or. CPML_regions(ispec_CPML) == CPML_XY_ONLY .or. &
          CPML_regions(ispec_CPML) == CPML_YZ_ONLY .or. CPML_regions(ispec_CPML) == CPML_XYZ) then
        if (ystore(iglob) - y_origin > ZERO) then
          if (ystore(iglob) - y_origin <= CPML_y_front - y_origin) then
            CPML_y_front = ystore(iglob)
          endif
        else
          if (abs(ystore(iglob) - y_origin) <= abs(CPML_y_back-y_origin)) then
            CPML_y_back = ystore(iglob)
          endif
        endif
      endif

      if (CPML_regions(ispec_CPML) == CPML_Z_ONLY .or. CPML_regions(ispec_CPML) == CPML_YZ_ONLY .or. &
        CPML_regions(ispec_CPML) == CPML_XZ_ONLY .or.  CPML_regions(ispec_CPML) == CPML_XYZ) then
        if (zstore(iglob) - z_origin > ZERO) then
          if (zstore(iglob) - z_origin <= CPML_z_top - z_origin) then
            CPML_z_top = zstore(iglob)
          endif
        else
          if (abs(zstore(iglob) - z_origin) <= abs(CPML_z_bottom-z_origin)) then
            CPML_z_bottom = zstore(iglob)
          endif
        endif
      endif

    enddo; enddo; enddo

  enddo

  CPML_width_x_right  = x_max_all - CPML_x_right
  CPML_width_x_left   = CPML_x_left - x_min_all
  CPML_width_y_front  = y_max_all - CPML_y_front
  CPML_width_y_back   = CPML_y_back - y_min_all
  CPML_width_z_top    = z_max_all - CPML_z_top
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
  CPML_width_z = max(CPML_width_z_bottom_max_all,CPML_width_z_top_max_all)

  if (PML_INSTEAD_OF_FREE_SURFACE) then
    zorigintop = z_max_all - CPML_width_z_top_max_all
  endif

! Calculation of maximum p velocity inside PML
  vp_max = ZERO
  do ispec_CPML = 1,nspec_cpml
    ispec = CPML_to_spec(ispec_CPML)
    do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
      vp = rho_vp(i,j,k,ispec) / rhostore(i,j,k,ispec)
      if (vp >= vp_max) then
        vp_max = vp
      endif
    enddo; enddo; enddo
  enddo

  call max_all_all_cr(vp_max,vp_max_all)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Boundary values of X-/Y-/Z-regions'
    write(IMAIN,*) minval(xstore(:)), maxval(xstore(:))
    write(IMAIN,*) minval(ystore(:)), maxval(ystore(:))
    write(IMAIN,*) minval(zstore(:)), maxval(zstore(:))
    write(IMAIN,*)
    write(IMAIN,*) 'Origins of right/left X-surface C-PML',xoriginright,xoriginleft
    write(IMAIN,*) 'Origins of front/back Y-surface C-PML',yoriginfront,yoriginback
    write(IMAIN,*) 'Origin of bottom Z-surface C-PML',zoriginbottom
    if (PML_INSTEAD_OF_FREE_SURFACE) then
      write(IMAIN,*) 'Origin of top Z-surface C-PML',zorigintop
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'CPML_width_x: ',CPML_width_x
    write(IMAIN,*) 'CPML_width_y: ',CPML_width_y
    write(IMAIN,*) 'CPML_width_z: ',CPML_width_z
    write(IMAIN,*)
  endif

  call synchronize_all()

  ! loops over all C-PML elements
  do ispec_CPML = 1,nspec_cpml

    ispec = CPML_to_spec(ispec_CPML)

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          ! calculates P-velocity
          if (ispec_is_acoustic(ispec)) then
! vp = rho_vp(i,j,k,ispec)/rhostore(i,j,k,ispec)
! For convenience only, when computing the damping profile inside PML,
! we set the required variable "vp" to be constant and equal to "vp_max_all"
            vp = vp_max_all
          else if (ispec_is_elastic(ispec)) then
! vp = rho_vp(i,j,k,ispec)/rhostore(i,j,k,ispec)
! For convenience only, when computing the damping profile inside PML,
! we set the required variable "vp" to be constant and equal to "vp_max_all"
            vp = vp_max_all
          else
            print *,'element index',ispec
            print *,'C-PML element index ',ispec_CPML
            call exit_mpi(myrank,'C-PML error: element has an unvalid P-velocity')
          endif

          iglob = ibool(i,j,k,ispec)

          if (CPML_regions(ispec_CPML) == CPML_X_ONLY) then
            !------------------------------------------------------------------------------
            !---------------------------- X-surface C-PML ---------------------------------
            !------------------------------------------------------------------------------
            if (xstore(iglob) - x_origin > ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xstore(iglob) - xoriginright

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML element's GLL point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! avoid d_x to be less than zero
              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

            else if (xstore(iglob) - x_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xoriginleft - xstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML grid point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

            else
              stop "there is error in mesh of CPML-layer x"
            endif

            !  stores damping profiles and auxiliary coefficients at the C-PML element's GLL points
            if (alpha_x < ZERO) stop "there is error in mesh of CPML-layer x"

            K_store_x(i,j,k,ispec_CPML) = K_x
            d_store_x(i,j,k,ispec_CPML) = d_x
            alpha_store_x(i,j,k,ispec_CPML) = alpha_x

            K_store_y(i,j,k,ispec_CPML) = ONE
            d_store_y(i,j,k,ispec_CPML) = ZERO
            alpha_store_y(i,j,k,ispec_CPML) = ZERO

            K_store_z(i,j,k,ispec_CPML) = ONE
            d_store_z(i,j,k,ispec_CPML) = ZERO
            alpha_store_z(i,j,k,ispec_CPML) = ZERO

          else if (CPML_regions(ispec_CPML) == CPML_Y_ONLY) then
            !------------------------------------------------------------------------------
            !---------------------------- Y-surface C-PML ---------------------------------
            !------------------------------------------------------------------------------
            if (ystore(iglob) - y_origin > ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = ystore(iglob) - yoriginfront

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

            else if (ystore(iglob) - y_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = yoriginback - ystore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

            else
              stop "there is error in mesh of  CPML-layer y"
            endif

            !  stores damping profiles and auxiliary coefficients at the C-PML element's GLL points
            if (alpha_y < ZERO) stop "there is error in mesh of  CPML-layer y"

            K_store_x(i,j,k,ispec_CPML) = ONE
            d_store_x(i,j,k,ispec_CPML) = ZERO
            alpha_store_x(i,j,k,ispec_CPML) = ZERO

            K_store_y(i,j,k,ispec_CPML) = K_y
            d_store_y(i,j,k,ispec_CPML) = d_y
            alpha_store_y(i,j,k,ispec_CPML) = alpha_y

            K_store_z(i,j,k,ispec_CPML) = ONE
            d_store_z(i,j,k,ispec_CPML) = ZERO
            alpha_store_z(i,j,k,ispec_CPML) = ZERO

          else if (CPML_regions(ispec_CPML) == CPML_Z_ONLY) then
            !------------------------------------------------------------------------------
            !---------------------------- Z-surface C-PML ---------------------------------
            !------------------------------------------------------------------------------
            if (zstore(iglob) - z_origin > ZERO) then

              if (PML_INSTEAD_OF_FREE_SURFACE) then
                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_z = zstore(iglob) - zorigintop

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_z / CPML_width_z

                ! gets damping profile at the C-PML element's GLL point
                K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

                if (K_z < ONE .or. d_z < ZERO) then
                  K_z = ONE; d_z = ZERO
                endif

              endif

            else if (zstore(iglob) - z_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_z = zoriginbottom - zstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_z / CPML_width_z

              ! gets damping profile at the C-PML element's GLL point
              K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
              alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

              if (K_z < ONE .or. d_z < ZERO) then
                K_z = ONE; d_z = ZERO
              endif

            else
              stop "there is error in mesh of CPML-layer z"
            endif

            !  stores damping profiles and auxiliary coefficients at the C-PML element's GLL points
            if (alpha_z < ZERO) stop "there is error in mesh of CPML-layer z"

            K_store_x(i,j,k,ispec_CPML) = ONE
            d_store_x(i,j,k,ispec_CPML) = ZERO
            alpha_store_x(i,j,k,ispec_CPML) = ZERO

            K_store_y(i,j,k,ispec_CPML) = ONE
            d_store_y(i,j,k,ispec_CPML) = ZERO
            alpha_store_y(i,j,k,ispec_CPML) = ZERO

            K_store_z(i,j,k,ispec_CPML) = K_z
            d_store_z(i,j,k,ispec_CPML) = d_z
            alpha_store_z(i,j,k,ispec_CPML) = alpha_z

          else if (CPML_regions(ispec_CPML) == CPML_XY_ONLY) then
            !------------------------------------------------------------------------------
            !---------------------------- XY-surface C-PML --------------------------------
            !------------------------------------------------------------------------------
            if (xstore(iglob) - x_origin > ZERO .and. ystore(iglob) - y_origin > ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xstore(iglob) - xoriginright

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML element's GLL point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = ystore(iglob) - yoriginfront

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

            else if (xstore(iglob) - x_origin > ZERO .and. ystore(iglob) - y_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xstore(iglob) - xoriginright

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML element's GLL point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = yoriginback - ystore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

            else if (xstore(iglob) - x_origin < ZERO .and. ystore(iglob) - y_origin > ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xoriginleft - xstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML element's GLL point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = ystore(iglob) - yoriginfront

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

            else if (xstore(iglob) - x_origin < ZERO .and. ystore(iglob) - y_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xoriginleft - xstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML element's GLL point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = yoriginback - ystore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

            else
              stop "there is error in mesh of CPML-layer xy"
            endif

            !  stores damping profiles and auxiliary coefficients at the C-PML element's GLL points
            if (alpha_x < ZERO .or. alpha_y < ZERO) stop "there is error in mesh of CPML-layer xy"

            K_store_x(i,j,k,ispec_CPML) = K_x
            d_store_x(i,j,k,ispec_CPML) = d_x
            alpha_store_x(i,j,k,ispec_CPML) = alpha_x

            K_store_y(i,j,k,ispec_CPML) = K_y
            d_store_y(i,j,k,ispec_CPML) = d_y
            alpha_store_y(i,j,k,ispec_CPML) = alpha_y

            K_store_z(i,j,k,ispec_CPML) = ONE
            d_store_z(i,j,k,ispec_CPML) = ZERO
            alpha_store_z(i,j,k,ispec_CPML) = ZERO

          else if (CPML_regions(ispec_CPML) == CPML_XZ_ONLY) then
            !------------------------------------------------------------------------------
            !---------------------------- XZ-surface C-PML --------------------------------
            !------------------------------------------------------------------------------
            if (xstore(iglob) - x_origin > ZERO .and. zstore(iglob) - z_origin > ZERO) then
              if (PML_INSTEAD_OF_FREE_SURFACE) then
                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_x = xstore(iglob) - xoriginright

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_x / CPML_width_x

                ! gets damping profile at the C-PML element's GLL point
                K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_z = zstore(iglob) - zorigintop

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_z / CPML_width_z

                ! gets damping profile at the C-PML element's GLL point
                K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

                if (K_x < ONE .or. d_x < ZERO) then
                  K_x = ONE; d_x = ZERO
                endif

                if (K_z < ONE .or. d_z < ZERO) then
                  K_z = ONE; d_z = ZERO
                endif

              endif

            else if (xstore(iglob) - x_origin > ZERO .and. zstore(iglob) - z_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xstore(iglob) - xoriginright

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML element's GLL point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_z = zoriginbottom - zstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_z / CPML_width_z

              ! gets damping profile at the C-PML element's GLL point
              K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
              alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

              if (K_z < ONE .or. d_z < ZERO) then
                K_z = ONE; d_z = ZERO
              endif

            else if (xstore(iglob) - x_origin < ZERO .and. zstore(iglob) - z_origin > ZERO) then
              if (PML_INSTEAD_OF_FREE_SURFACE) then
                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_x = xoriginleft - xstore(iglob)

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_x / CPML_width_x

                ! gets damping profile at the C-PML element's GLL point
                K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_z = zstore(iglob) - zorigintop

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_z / CPML_width_z

                ! gets damping profile at the C-PML element's GLL point
                K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

                if (K_x < ONE .or. d_x < ZERO) then
                  K_x = ONE; d_x = ZERO
                endif

                if (K_z < ONE .or. d_z < ZERO) then
                  K_z = ONE; d_z = ZERO
                endif
              endif

            else if (xstore(iglob) - x_origin < ZERO .and. zstore(iglob) - z_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xoriginleft - xstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML element's GLL point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_z = zoriginbottom - zstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_z / CPML_width_z

              ! gets damping profile at the C-PML element's GLL point
              K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
              alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

              if (K_z < ONE .or. d_z < ZERO) then
                K_z = ONE; d_z = ZERO
              endif
            else
              stop "there is error in mesh of CPML-layer xz"
            endif

            !  stores damping profiles and auxiliary coefficients at the C-PML element's GLL points
            if (alpha_x < ZERO .or. alpha_z < ZERO) stop "there is error in mesh of CPML-layer xz"

            K_store_x(i,j,k,ispec_CPML) = K_x
            d_store_x(i,j,k,ispec_CPML) = d_x
            alpha_store_x(i,j,k,ispec_CPML) = alpha_x

            K_store_y(i,j,k,ispec_CPML) = ONE
            d_store_y(i,j,k,ispec_CPML) = ZERO
            alpha_store_y(i,j,k,ispec_CPML) = ZERO

            K_store_z(i,j,k,ispec_CPML) = K_z
            d_store_z(i,j,k,ispec_CPML) = d_z
            alpha_store_z(i,j,k,ispec_CPML) = alpha_z

          else if (CPML_regions(ispec_CPML) == CPML_YZ_ONLY) then
            !------------------------------------------------------------------------------
            !---------------------------- YZ-surface C-PML --------------------------------
            !------------------------------------------------------------------------------
            if (ystore(iglob) - y_origin > ZERO .and. zstore(iglob) - z_origin > ZERO) then
              if (PML_INSTEAD_OF_FREE_SURFACE) then
                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_y = ystore(iglob) - yoriginfront

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_y / CPML_width_y

                ! gets damping profile at the C-PML element's GLL point
                K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_z = zstore(iglob) - zorigintop

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_z / CPML_width_z

                ! gets damping profile at the C-PML element's GLL point
                K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

                if (d_y < ZERO  .or. K_y < ONE) then
                  K_y = ONE; d_y = ZERO
                endif

                if (d_z < ZERO  .or. K_z < ONE) then
                  K_z = ONE; d_z = ZERO
                endif
              endif

            else if (ystore(iglob) - y_origin > ZERO .and. zstore(iglob) - z_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = ystore(iglob) - yoriginfront

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_z = zoriginbottom - zstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_z / CPML_width_z

              ! gets damping profile at the C-PML element's GLL point
              K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
              alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

              if (K_z < ONE .or. d_z < ZERO) then
                K_z = ONE; d_z = ZERO
              endif

            else if (ystore(iglob) - y_origin < ZERO .and. zstore(iglob) - z_origin > ZERO) then
              if (PML_INSTEAD_OF_FREE_SURFACE) then
                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_y = yoriginback - ystore(iglob)

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_y / CPML_width_y

                ! gets damping profile at the C-PML element's GLL point
                K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_z = zstore(iglob) - zorigintop

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_z / CPML_width_z

                ! gets damping profile at the C-PML element's GLL point
                K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

                if (K_y < ONE .or. d_y < ZERO) then
                  K_y = ONE; d_y = ZERO
                endif

                if (K_z < ONE .or. d_z < ZERO) then
                  K_z = ONE; d_z = ZERO
                endif
              endif

            else if (ystore(iglob) - y_origin < ZERO .and. zstore(iglob) - z_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = yoriginback - ystore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_z = zoriginbottom - zstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_z / CPML_width_z

              ! gets damping profile at the C-PML element's GLL point
              K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
              alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

              if (K_z < ONE .or. d_z < ZERO) then
                K_z = ONE; d_z = ZERO
              endif

            else
              stop "there is error in mesh of CPML-layer yz"
            endif

            !! DK DK define an alias for y and z variable names (which are the same)
            if (alpha_y < ZERO .or. alpha_z < ZERO) stop "there is error in mesh of CPML-layer yz"

            K_store_x(i,j,k,ispec_CPML) = ONE
            d_store_x(i,j,k,ispec_CPML) = ZERO
            alpha_store_x(i,j,k,ispec_CPML) = ZERO

            K_store_y(i,j,k,ispec_CPML) = K_y
            d_store_y(i,j,k,ispec_CPML) = d_y
            alpha_store_y(i,j,k,ispec_CPML) = alpha_y

            K_store_z(i,j,k,ispec_CPML) = K_z
            d_store_z(i,j,k,ispec_CPML) = d_z
            alpha_store_z(i,j,k,ispec_CPML) = alpha_z

          else if (CPML_regions(ispec_CPML) == CPML_XYZ) then
            !------------------------------------------------------------------------------
            !---------------------------- XYZ-surface C-PML -------------------------------
            !------------------------------------------------------------------------------
            if (xstore(iglob) - x_origin > ZERO .and. ystore(iglob) - y_origin>ZERO .and. &
                zstore(iglob) - z_origin > ZERO) then
              if (PML_INSTEAD_OF_FREE_SURFACE) then
                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_x = xstore(iglob) - xoriginright

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_x / CPML_width_x

                ! gets damping profile at the C-PML grid point
                K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_y = ystore(iglob) - yoriginfront

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_y / CPML_width_y

                ! gets damping profile at the C-PML element's GLL point
                K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_z = zstore(iglob) - zorigintop

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_z / CPML_width_z

                ! gets damping profile at the C-PML element's GLL point
                K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

                if (K_x < ONE .or. d_x < ZERO) then
                  K_x = ONE; d_x = ZERO
                endif

                if (K_y < ONE .or. d_y < ZERO) then
                  K_y = ONE; d_y = ZERO
                endif

                if (K_z < ONE .or. d_z < ZERO) then
                  K_z = ONE; d_z = ZERO
                endif

              endif

            else if (xstore(iglob) - x_origin > ZERO .and. ystore(iglob) - y_origin > ZERO .and. &
                     zstore(iglob) - z_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xstore(iglob) - xoriginright

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML grid point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = ystore(iglob) - yoriginfront

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_z = zoriginbottom - zstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_z / CPML_width_z

              ! gets damping profile at the C-PML element's GLL point
              K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
              alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

              if (K_z < ONE .or. d_z < ZERO) then
                K_z = ONE; d_z = ZERO
              endif

            else if (xstore(iglob) - x_origin > ZERO .and. ystore(iglob) - y_origin < ZERO .and. &
                     zstore(iglob) - z_origin > ZERO) then
              if (PML_INSTEAD_OF_FREE_SURFACE) then
                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_x = xstore(iglob) - xoriginright

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_x / CPML_width_x

                ! gets damping profile at the C-PML element's GLL point
                K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_y = yoriginback - ystore(iglob)

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_y / CPML_width_y

                ! gets damping profile at the C-PML element's GLL point
                K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_z = zstore(iglob) - zorigintop

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_z / CPML_width_z

                ! gets damping profile at the C-PML element's GLL point
                K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

                if (K_x < ONE .or. d_x < ZERO) then
                  K_x = ONE; d_x = ZERO
                endif

                if (K_y < ONE .or. d_y < ZERO) then
                  K_y = ONE; d_y = ZERO
                endif

                if (K_z < ONE .or. d_z < ZERO) then
                  K_z = ONE; d_z = ZERO
                endif
              endif

            else if (xstore(iglob) - x_origin > ZERO .and. ystore(iglob) - y_origin < ZERO .and. &
                     zstore(iglob) - z_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xstore(iglob) - xoriginright

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML element's GLL point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = yoriginback - ystore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_z = zoriginbottom - zstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_z / CPML_width_z

              ! gets damping profile at the C-PML element's GLL point
              K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
              alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

              if (K_z < ONE .or. d_z < ZERO) then
                K_z = ONE; d_z = ZERO
              endif

            else if (xstore(iglob) - x_origin < ZERO .and. ystore(iglob) - y_origin > ZERO .and. &
                     zstore(iglob) - z_origin > ZERO) then
              if (PML_INSTEAD_OF_FREE_SURFACE) then
                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_x = xoriginleft - xstore(iglob)

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_x / CPML_width_x

                ! gets damping profile at the C-PML grid point
                K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_y = ystore(iglob) - yoriginfront

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_y / CPML_width_y

                ! gets damping profile at the C-PML element's GLL point
                K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_z = zstore(iglob) - zorigintop

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_z / CPML_width_z

                ! gets damping profile at the C-PML element's GLL point
                K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

                if (K_x < ONE .or. d_x < ZERO) then
                  K_x = ONE; d_x = ZERO
                endif

                if (K_y < ONE .or. d_y < ZERO) then
                  K_y = ONE; d_y = ZERO
                endif

                if (K_z < ONE .or. d_z < ZERO) then
                  K_z = ONE; d_z = ZERO
                endif
              endif

            else if (xstore(iglob) - x_origin < ZERO .and. ystore(iglob) - y_origin > ZERO .and. &
                     zstore(iglob) - z_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xoriginleft - xstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML grid point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = ystore(iglob) - yoriginfront

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_z = zoriginbottom - zstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_z / CPML_width_z

              ! gets damping profile at the C-PML element's GLL point
              K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
              alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

              if (K_z < ONE .or. d_z < ZERO) then
                K_z = ONE; d_z = ZERO
              endif

            else if (xstore(iglob) - x_origin < ZERO .and. ystore(iglob) - y_origin < ZERO .and. &
                     zstore(iglob) - z_origin > ZERO) then
              if (PML_INSTEAD_OF_FREE_SURFACE) then
                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_x = xoriginleft - xstore(iglob)

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_x / CPML_width_x

                ! gets damping profile at the C-PML element's GLL point
                K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
                alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_y = yoriginback - ystore(iglob)

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_y / CPML_width_y

                ! gets damping profile at the C-PML element's GLL point
                K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
                alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

                ! gets abscissa of current grid point along the damping profile
                abscissa_in_PML_z = zstore(iglob) - zorigintop

                ! determines distance to C-PML/mesh interface
                dist = abscissa_in_PML_z / CPML_width_z

                ! gets damping profile at the C-PML element's GLL point
                K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
                d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
                alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

                if (K_x < ONE .or. d_x < ZERO) then
                  K_x = ONE; d_x = ZERO
                endif

                if (K_y < ONE .or. d_y < ZERO) then
                  K_y = ONE; d_y = ZERO
                endif

                if (K_z < ONE .or. d_z < ZERO) then
                  K_z = ONE; d_z = ZERO
                endif
              endif

            else if (xstore(iglob) - x_origin < ZERO .and. ystore(iglob) - y_origin < ZERO .and. &
                     zstore(iglob) - z_origin < ZERO) then
              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_x = xoriginleft - xstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_x / CPML_width_x

              ! gets damping profile at the C-PML element's GLL point
              K_x = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_x = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
              alpha_x = ALPHA_MAX_PML_x * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_y = yoriginback - ystore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_y / CPML_width_y

              ! gets damping profile at the C-PML element's GLL point
              K_y = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_y = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
              alpha_y = ALPHA_MAX_PML_y * (ONE - dist)

              ! gets abscissa of current grid point along the damping profile
              abscissa_in_PML_z = zoriginbottom - zstore(iglob)

              ! determines distance to C-PML/mesh interface
              dist = abscissa_in_PML_z / CPML_width_z

              ! gets damping profile at the C-PML element's GLL point
              K_z = K_MIN_PML + (K_MAX_PML - ONE) * dist
              d_z = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
              alpha_z = ALPHA_MAX_PML_z * (ONE - dist)

              if (K_x < ONE .or. d_x < ZERO) then
                K_x = ONE; d_x = ZERO
              endif

              if (K_y < ONE .or. d_y < ZERO) then
                K_y = ONE; d_y = ZERO
              endif

              if (K_z < ONE .or. d_z < ZERO) then
                K_z = ONE; d_z = ZERO
              endif

            else
              stop "there is error in mesh of CPML-layer xyz"
            endif

            !! DK DK define an alias for y and z variable names (which are the same)
            if (alpha_x < ZERO .or. alpha_y < ZERO .or. alpha_z < ZERO) stop "there is error in mesh of CPML-layer xyz"

            K_store_x(i,j,k,ispec_CPML) = K_x
            d_store_x(i,j,k,ispec_CPML) = d_x
            alpha_store_x(i,j,k,ispec_CPML) = alpha_x

            K_store_y(i,j,k,ispec_CPML) = K_y
            d_store_y(i,j,k,ispec_CPML) = d_y
            alpha_store_y(i,j,k,ispec_CPML) = alpha_y

            K_store_z(i,j,k,ispec_CPML) = K_z
            d_store_z(i,j,k,ispec_CPML) = d_z
            alpha_store_z(i,j,k,ispec_CPML) = alpha_z


          endif
        enddo
      enddo
    enddo
  enddo !ispec_CPML

  do ispec_CPML = 1,nspec_cpml
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          if (CPML_regions(ispec_CPML) == CPML_XY_ONLY) then

            K_x = K_store_x(i,j,k,ispec_CPML)
            d_x = d_store_x(i,j,k,ispec_CPML)
            alpha_x = alpha_store_x(i,j,k,ispec_CPML)

            K_y = K_store_y(i,j,k,ispec_CPML)
            d_y = d_store_y(i,j,k,ispec_CPML)
            alpha_y = alpha_store_y(i,j,k,ispec_CPML)

            call seperate_two_changeable_value(alpha_x,alpha_y,ALPHA_MAX_PML_y)

            if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
              stop 'there is error in seperation of alpha_x, alpha_y'
            endif

            beta_x = alpha_x + d_x / K_x
            beta_y = alpha_y + d_y / K_y

            call seperate_two_value_with_one_changeable(beta_x,alpha_y,ALPHA_MAX_PML_y)
            call seperate_two_value_with_one_changeable(beta_y,alpha_x,ALPHA_MAX_PML_y)

            if (abs(beta_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                abs(beta_y - alpha_x) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
              stop 'there is error in seperation of beta_x and alpha_y, beta_y and alpha_x'
            endif

            d_x = (beta_x - alpha_x) * K_x
            d_y = (beta_y - alpha_y) * K_y

            d_store_x(i,j,k,ispec_CPML) = d_x
            alpha_store_x(i,j,k,ispec_CPML) = alpha_x

            d_store_y(i,j,k,ispec_CPML) = d_y
            alpha_store_y(i,j,k,ispec_CPML) = alpha_y

          else if (CPML_regions(ispec_CPML) == CPML_XZ_ONLY) then

            K_x = K_store_x(i,j,k,ispec_CPML)
            d_x = d_store_x(i,j,k,ispec_CPML)
            alpha_x = alpha_store_x(i,j,k,ispec_CPML)

            K_z = K_store_z(i,j,k,ispec_CPML)
            d_z = d_store_z(i,j,k,ispec_CPML)
            alpha_z = alpha_store_z(i,j,k,ispec_CPML)

            call seperate_two_changeable_value(alpha_x,alpha_z,ALPHA_MAX_PML_y)

            if (abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
              stop 'there is error in seperation of alpha_x, alpha_z'
            endif

            beta_x = alpha_x + d_x / K_x
            beta_z = alpha_z + d_z / K_z

            call seperate_two_value_with_one_changeable(beta_x,alpha_z,ALPHA_MAX_PML_y)
            call seperate_two_value_with_one_changeable(beta_z,alpha_x,ALPHA_MAX_PML_y)

            if (abs(beta_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                abs(beta_z - alpha_x) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
              stop 'there is error in seperation of beta_x and alpha_z, beta_z and alpha_z'
            endif

            d_x = (beta_x - alpha_x) * K_x
            d_z = (beta_z - alpha_z) * K_z

            d_store_x(i,j,k,ispec_CPML) = d_x
            alpha_store_x(i,j,k,ispec_CPML) = alpha_x

            d_store_z(i,j,k,ispec_CPML) = d_z
            alpha_store_z(i,j,k,ispec_CPML) = alpha_z

          else if (CPML_regions(ispec_CPML) == CPML_YZ_ONLY) then

            K_y = K_store_y(i,j,k,ispec_CPML)
            d_y = d_store_y(i,j,k,ispec_CPML)
            alpha_y = alpha_store_y(i,j,k,ispec_CPML)

            K_z = K_store_z(i,j,k,ispec_CPML)
            d_z = d_store_z(i,j,k,ispec_CPML)
            alpha_z = alpha_store_z(i,j,k,ispec_CPML)

            call seperate_two_changeable_value(alpha_y,alpha_z,ALPHA_MAX_PML_y)

            if (abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
              stop 'there is error in seperation of alpha_y, alpha_z'
            endif

            beta_y = alpha_y + d_y / K_y
            beta_z = alpha_z + d_z / K_z

            call seperate_two_value_with_one_changeable(beta_y,alpha_z,ALPHA_MAX_PML_y)
            call seperate_two_value_with_one_changeable(beta_z,alpha_y,ALPHA_MAX_PML_y)

            if (abs(beta_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                abs(beta_z - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
              stop 'there is error in seperation of beta_y and alpha_z, beta_z and alpha_y'
            endif

            d_y = (beta_y - alpha_y) * K_y
            d_z = (beta_z - alpha_z) * K_z

            d_store_y(i,j,k,ispec_CPML) = d_y
            alpha_store_y(i,j,k,ispec_CPML) = alpha_y

            d_store_z(i,j,k,ispec_CPML) = d_z
            alpha_store_z(i,j,k,ispec_CPML) = alpha_z

          else if(CPML_regions(ispec_CPML) == CPML_XYZ) then

            K_x = K_store_x(i,j,k,ispec_CPML)
            d_x = d_store_x(i,j,k,ispec_CPML)
            alpha_x = alpha_store_x(i,j,k,ispec_CPML)

            K_y = K_store_y(i,j,k,ispec_CPML)
            d_y = d_store_y(i,j,k,ispec_CPML)
            alpha_y = alpha_store_y(i,j,k,ispec_CPML)

            K_z = K_store_z(i,j,k,ispec_CPML)
            d_z = d_store_z(i,j,k,ispec_CPML)
            alpha_z = alpha_store_z(i,j,k,ispec_CPML)

            if (alpha_x == alpha_y .and. alpha_y == alpha_z) then
              meanval_1 = (alpha_x + alpha_y + alpha_z) / 3.0_CUSTOM_REAL + 0.011_CUSTOM_REAL * ALPHA_MAX_PML_y
              alpha_x = meanval_1 - 0.01_CUSTOM_REAL * ALPHA_MAX_PML_y
              alpha_y = meanval_1
              alpha_z = meanval_1 + 0.01_CUSTOM_REAL * ALPHA_MAX_PML_y

              if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
                stop 'error in seperation of alpha_x, alpha_y, alpha_z in case alpha_x == alpha_y .and. alpha_y == alpha_z'
              endif

            else if (alpha_x == alpha_y .and. alpha_y /= alpha_z) then
              call seperate_three_changeable_value_with_two_equal(alpha_x,alpha_y,alpha_z,ALPHA_MAX_PML_y)

              if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
                stop 'error in seperation of alpha_x, alpha_y, alpha_z in case alpha_x == alpha_y .and. alpha_y /= alpha_z'
              endif

            else if (alpha_y == alpha_z .and. alpha_x /= alpha_y) then
              call seperate_three_changeable_value_with_two_equal(alpha_y,alpha_z,alpha_x,ALPHA_MAX_PML_y)

              if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
                stop 'error in seperation of alpha_x, alpha_y, alpha_z in case alpha_y == alpha_z .and. alpha_x /= alpha_y'
              endif

            else if (alpha_x == alpha_z .and. alpha_y /= alpha_z) then
              call seperate_three_changeable_value_with_two_equal(alpha_x,alpha_z,alpha_y,ALPHA_MAX_PML_y)

              if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or.&
                  abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
                stop 'error in seperation of alpha_x, alpha_y, alpha_z in case alpha_x == alpha_z .and. alpha_y /= alpha_z'
              endif

            else if (alpha_x > alpha_y .and. alpha_y > alpha_z) then
              call seperate_three_sequential_changeable_value(alpha_x,alpha_y,alpha_z,ALPHA_MAX_PML_y)

              if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
                stop 'error in seperation of alpha_x, alpha_y, alpha_z in case alpha_x > alpha_y .and. alpha_y > alpha_z'
              endif

            else if (alpha_x > alpha_z .and. alpha_z > alpha_y) then
              call seperate_three_sequential_changeable_value(alpha_x,alpha_z,alpha_y,ALPHA_MAX_PML_y)

              if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
                stop 'error in seperation of alpha_x, alpha_y, alpha_z in case alpha_x > alpha_z .and. alpha_z > alpha_y'
              endif

            else if (alpha_y > alpha_z .and. alpha_z > alpha_x) then
              call seperate_three_sequential_changeable_value(alpha_y,alpha_z,alpha_x,ALPHA_MAX_PML_y)

              if(abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                 abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                 abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
                stop 'error in seperation of alpha_x, alpha_y, alpha_z in case alpha_y > alpha_z .and. alpha_z > alpha_x'
              endif

            else if (alpha_y > alpha_x .and. alpha_x > alpha_z) then
              call seperate_three_sequential_changeable_value(alpha_y,alpha_x,alpha_z,ALPHA_MAX_PML_y)

              if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
                stop 'error in seperation of alpha_x, alpha_y, alpha_z in case alpha_y > alpha_x .and. alpha_x > alpha_z'
              endif

            else if (alpha_z > alpha_x .and. alpha_x > alpha_y) then
              call seperate_three_sequential_changeable_value(alpha_z,alpha_x,alpha_y,ALPHA_MAX_PML_y)

              if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                  abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
                stop 'error in seperation of alpha_x, alpha_y, alpha_z in case alpha_z > alpha_x .and. alpha_x > alpha_y'
              endif

            else if (alpha_z > alpha_y .and. alpha_y > alpha_x) then
              call seperate_three_sequential_changeable_value(alpha_z,alpha_y,alpha_x,ALPHA_MAX_PML_y)

              if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y.or.&
                  abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y.or.&
                  abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
                stop 'error in seperation of alpha_x, alpha_y, alpha_z in case alpha_z > alpha_y .and. alpha_y > alpha_x'
              endif

            else
              stop 'there is error in alpha division'
            endif

            if (abs(alpha_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                abs(alpha_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                abs(alpha_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
              stop 'there is error in seperation of alpha_x, alpha_y, alpha_z '
            endif

            beta_x = alpha_x + d_x / K_x
            beta_y = alpha_y + d_y / K_y
            beta_z = alpha_z + d_z / K_z

            call seperate_betai_alphaj_alphak(beta_x,alpha_y,alpha_z,ALPHA_MAX_PML_y)
            call seperate_betai_alphaj_alphak(beta_y,alpha_z,alpha_x,ALPHA_MAX_PML_y)
            call seperate_betai_alphaj_alphak(beta_z,alpha_x,alpha_y,ALPHA_MAX_PML_y)

            if (abs(beta_x - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                abs(beta_x - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
              stop 'there is error in seperation of beta_x, alpha_y,alpha_z '
            endif

            if (abs(beta_y - alpha_x) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y .or. &
                abs(beta_y - alpha_z) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
              stop 'there is error in seperation of beta_y, alpha_x,alpha_z '
            endif

            if (abs(beta_z - alpha_x) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y.or. &
                abs(beta_z - alpha_y) < const_for_tune_pml_damping_profile * ALPHA_MAX_PML_y) then
              stop 'there is error in seperation of beta_z, alpha_x,alpha_y '
            endif

            d_x = (beta_x - alpha_x) * K_x
            d_y = (beta_y - alpha_y) * K_y
            d_z = (beta_z - alpha_z) * K_z

            d_store_x(i,j,k,ispec_CPML) = d_x
            alpha_store_x(i,j,k,ispec_CPML) = alpha_x
            d_store_y(i,j,k,ispec_CPML) = d_y
            alpha_store_y(i,j,k,ispec_CPML) = alpha_y
            d_store_z(i,j,k,ispec_CPML) = d_z
            alpha_store_z(i,j,k,ispec_CPML) = alpha_z

          endif

        enddo
      enddo
    enddo

  enddo


! --------------------------------------------------------------------------------------------
! for adjoint tomography
! create the array store the points on interface between PML and interior computational domain
! --------------------------------------------------------------------------------------------

  if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then

    !mask all points belong interior computational domain
    allocate(mask_ibool_interior_domain(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'error allocating array mask_ibool_interior_domain'
    mask_ibool_interior_domain = .false.
    do ispec = 1,nspec
      if (.not. is_CPML(ispec)) then
        do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
          mask_ibool_interior_domain(ibool(i,j,k,ispec)) = .true.
        enddo; enddo; enddo
      endif
    enddo

    !------------------------------------------------------
    !  begin of acoustic domain
    !------------------------------------------------------
    nglob_interface_PML_acoustic = 0

    if (ACOUSTIC_SIMULATION) then

      do ispec = 1,nspec
        if (ispec_is_acoustic(ispec) .and. is_CPML(ispec)) then
          do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
            if (mask_ibool_interior_domain(ibool(i,j,k,ispec))) then
              nglob_interface_PML_acoustic = nglob_interface_PML_acoustic + 1
            endif
          enddo; enddo; enddo
        endif
      enddo

      if (nglob_interface_PML_acoustic > 0) then
        allocate(points_interface_PML_acoustic(nglob_interface_PML_acoustic),stat=ier)
        if (ier /= 0) stop 'error allocating array points_interface_PML_acoustic'
        points_interface_PML_acoustic = 0
        nglob_interface_PML_acoustic = 0
        do ispec=1,nspec
          if (ispec_is_acoustic(ispec) .and. is_CPML(ispec)) then
            do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
              if (mask_ibool_interior_domain(ibool(i,j,k,ispec))) then
                nglob_interface_PML_acoustic = nglob_interface_PML_acoustic + 1
                points_interface_PML_acoustic(nglob_interface_PML_acoustic) = ibool(i,j,k,ispec)
              endif
            enddo; enddo; enddo
          endif
        enddo
      endif

    endif
    !------------------------------------------------------
    ! end of acoustic domains
    !------------------------------------------------------

    !------------------------------------------------------
    ! begin of elastic domains
    !------------------------------------------------------
    nglob_interface_PML_elastic = 0

    if (ELASTIC_SIMULATION) then

      do ispec=1,nspec
        if (ispec_is_elastic(ispec) .and. is_CPML(ispec)) then
          do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
            if (mask_ibool_interior_domain(ibool(i,j,k,ispec))) then
              nglob_interface_PML_elastic = nglob_interface_PML_elastic + 1
            endif
          enddo; enddo; enddo
        endif
      enddo

      if (nglob_interface_PML_elastic > 0) then
        allocate(points_interface_PML_elastic(nglob_interface_PML_elastic),stat=ier)
        if (ier /= 0) stop 'error allocating array points_interface_PML_elastic'
        points_interface_PML_elastic = 0
        nglob_interface_PML_elastic = 0
        do ispec = 1,nspec
          if (ispec_is_elastic(ispec) .and. is_CPML(ispec)) then
            do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
              if (mask_ibool_interior_domain(ibool(i,j,k,ispec))) then
                nglob_interface_PML_elastic = nglob_interface_PML_elastic + 1
                points_interface_PML_elastic(nglob_interface_PML_elastic) = ibool(i,j,k,ispec)
              endif
            enddo; enddo; enddo
          endif
        enddo
      endif

    endif
    !------------------------------------------------------
    ! end of elastic domains
    !------------------------------------------------------
  endif
end subroutine pml_set_local_dampingcoeff

!
!-------------------------------------------------------------------------------------------------
!

function pml_damping_profile_l(myrank,iglob,dist,vp,delta)

  ! defines d, the damping profile at the C-PML element's GLL point for a given:
  !   dist:  distance to C-PML/mesh interface
  !   vp:    P-velocity
  !   delta: thickness of the C-PML layer

  use generate_databases_par, only: CUSTOM_REAL,NPOWER,CPML_Rcoef,TWO

  implicit none

  integer, intent(in) :: myrank,iglob

  real(kind=CUSTOM_REAL), intent(in) :: dist,vp,delta

  real(kind=CUSTOM_REAL) :: pml_damping_profile_l

  ! gets damping profile
  if (NPOWER >= 1) then
     ! In INRIA research report section 6.1:  http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
     ! pml_damping_profile_l = - ((NPOWER + 1) * vp * log(CPML_Rcoef) / (TWO * delta)) * dist**(NPOWER)
     ! due to tests it is more accurate to use following definition in case NPOWER = 1 defined in constants.h.in
    pml_damping_profile_l = - ((NPOWER + 1) * vp * log(CPML_Rcoef) / (TWO * delta)) &
                            * dist**(1.2_CUSTOM_REAL * NPOWER)
  else
    call exit_mpi(myrank,'C-PML error: NPOWER must be greater than or equal to 1')
  endif

  ! checks coordinates of C-PML points and thickness of C-PML layer
  ! the distance is a relative distance here (a ratio), thus we compare to 1 instead of comparing to the thickness of the PML
  if (dist > 1._CUSTOM_REAL ) then
    print *,'C-PML point ',iglob
    print *,'distance to C-PML/mesh interface ',dist
    print *,'C-PML thickness ',delta
    call exit_mpi(myrank,'C-PML error: distance to C-PML/mesh interface is bigger than thickness of C-PML layer')
  endif

end function pml_damping_profile_l

subroutine seperate_two_changeable_value(value_a,value_b,ALPHA_MAX_PML_y)

  use generate_databases_par, only: CUSTOM_REAL,TWO
  implicit none
  real(kind=CUSTOM_REAL), intent(in) :: ALPHA_MAX_PML_y
  real(kind=CUSTOM_REAL) :: value_a,value_b
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation = 0.009_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_two = 0.01_CUSTOM_REAL
  if (abs(value_a - value_b) < const_for_sepertation * ALPHA_MAX_PML_y) then
     if(value_a >= value_b) then
        value_a = (value_a + value_b) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
     else
        value_b = (value_a + value_b) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
     endif
  endif

end subroutine seperate_two_changeable_value

subroutine seperate_three_changeable_value_with_two_equal(value_a,value_b,value_c,ALPHA_MAX_PML_y)

  use generate_databases_par, only: CUSTOM_REAL,TWO
  implicit none
  real(kind=CUSTOM_REAL), intent(in) :: ALPHA_MAX_PML_y
  real(kind=CUSTOM_REAL) :: value_a,value_b,value_c
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation = 0.009_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_two = 0.01_CUSTOM_REAL

  if (value_a /= value_b) then
    stop 'value_a should be equal to value_b'
  endif

  value_a = value_b + const_for_sepertation_two * ALPHA_MAX_PML_y

  if (abs(value_b - value_c) < const_for_sepertation * ALPHA_MAX_PML_y) then
    if (value_b >= value_c) then
      value_b = (value_b + value_c) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
      value_a = value_b + const_for_sepertation_two * ALPHA_MAX_PML_y
    else
      value_c = (value_b + value_c) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
      value_a = value_c + const_for_sepertation_two * ALPHA_MAX_PML_y
    endif
  endif

end subroutine seperate_three_changeable_value_with_two_equal

subroutine seperate_two_value_with_one_changeable(value_a,value_b,ALPHA_MAX_PML_y)

  use generate_databases_par, only: CUSTOM_REAL,TWO
  implicit none
  real(kind=CUSTOM_REAL), intent(in) :: ALPHA_MAX_PML_y,value_b
  real(kind=CUSTOM_REAL) :: value_a
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation = 0.009_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_two = 0.01_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_three = 0.015_CUSTOM_REAL

  if (abs(value_a - value_b) < const_for_sepertation * ALPHA_MAX_PML_y) then
    if (value_a >= value_b) then
      value_a = (value_a + value_b) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
    else
      value_a = (value_a + value_b) / TWO + const_for_sepertation_three * ALPHA_MAX_PML_y
    endif
  endif

end subroutine seperate_two_value_with_one_changeable

subroutine seperate_three_sequential_changeable_value(value_max,value_mean,value_min,ALPHA_MAX_PML_y)

  use generate_databases_par, only: CUSTOM_REAL,TWO
  implicit none
  real(kind=CUSTOM_REAL), intent(in) :: ALPHA_MAX_PML_y
  real(kind=CUSTOM_REAL) :: value_max,value_mean,value_min
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation = 0.009_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_two = 0.01_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_three = 0.015_CUSTOM_REAL

  if (.not. (value_max > value_mean .and. value_mean > value_min)) then
    stop 'It should be value_max > value_mean and value_mean > value_min'
  endif

  if (abs(value_max - value_mean) >= const_for_sepertation * ALPHA_MAX_PML_y .and. &
      abs(value_mean - value_min) >= const_for_sepertation * ALPHA_MAX_PML_y) then
    ! To facilitate error check, we left the above if sentence, but leave no operation inside purposely
  else if (abs(value_max - value_mean) < const_for_sepertation * ALPHA_MAX_PML_y .and. &
           abs(value_mean - value_min) >= const_for_sepertation * ALPHA_MAX_PML_y) then
    value_max = (value_max + value_mean) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
  else if (abs(value_max - value_mean) >= const_for_sepertation * ALPHA_MAX_PML_y .and. &
           abs(value_mean - value_min) < const_for_sepertation * ALPHA_MAX_PML_y) then
    value_mean = (value_mean + value_min) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
    if (abs(value_max - value_mean) < const_for_sepertation * ALPHA_MAX_PML_y) then
      if (max(value_max,value_mean) == value_max) then
        value_max = (value_max + value_mean) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
      else if (max(value_max,value_mean) == value_mean) then
           value_mean = (value_mean + value_max) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
      endif
    endif
  else if (abs(value_max - value_mean) < const_for_sepertation * ALPHA_MAX_PML_y .and. &
           abs(value_mean - value_min) < const_for_sepertation * ALPHA_MAX_PML_y) then
    value_mean = (value_mean + value_min) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y

    if (abs(value_max - value_min) < const_for_sepertation * ALPHA_MAX_PML_y) then
      if (max(value_max,value_min) == value_max) then
        value_max = (value_max + value_min) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
      else if (max(value_max,value_min) == value_min) then
        value_max = (value_max + value_min) / TWO + const_for_sepertation_three * ALPHA_MAX_PML_y
      endif
    endif

    if (abs(value_max - value_mean) < const_for_sepertation * ALPHA_MAX_PML_y) then
      if (max(value_max,value_mean) == value_max) then
        value_max = (value_max + value_mean) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
      else if (max(value_max,value_mean) == value_mean) then
        value_max = (value_max + value_mean) / TWO + const_for_sepertation_three * ALPHA_MAX_PML_y
      endif
    endif

  endif
end subroutine seperate_three_sequential_changeable_value

subroutine seperate_three_sequential_value_with_only_meanvalue_changeable(value_max,value_mean,value_min,ALPHA_MAX_PML_y)

  use generate_databases_par, only: CUSTOM_REAL,TWO
  implicit none
  real(kind=CUSTOM_REAL), intent(in) :: ALPHA_MAX_PML_y,value_max,value_min
  real(kind=CUSTOM_REAL) :: value_mean
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation = 0.009_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_two = 0.01_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_three = 0.015_CUSTOM_REAL

  if (.not. (value_max > value_mean .and. value_mean > value_min)) then
    stop 'It should be value_max > value_mean and value_mean > value_min'
  endif

  if (abs(value_max - value_min) < const_for_sepertation * ALPHA_MAX_PML_y ) then
    stop 'It should be abs(value_max - value_mean) >= 0.009_CUSTOM_REAL * ALPHA_MAX_PML_y'
  endif

  if (abs(value_max - value_mean) >= const_for_sepertation * ALPHA_MAX_PML_y .and. &
      abs(value_mean - value_min) >= const_for_sepertation * ALPHA_MAX_PML_y) then
    ! To facilitate error check, we left the above if sentence, but leave no operation inside purposely
  else if (abs(value_max - value_mean) < const_for_sepertation * ALPHA_MAX_PML_y .and. &
           abs(value_mean - value_min) >= const_for_sepertation * ALPHA_MAX_PML_y) then
    value_mean = (value_max + value_mean) / TWO + const_for_sepertation_three * ALPHA_MAX_PML_y
  else if (abs(value_max - value_mean) >= const_for_sepertation * ALPHA_MAX_PML_y .and. &
           abs(value_mean - value_min) < const_for_sepertation * ALPHA_MAX_PML_y) then
    value_mean = (value_min + value_mean) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
    if (abs(value_max - value_mean) < const_for_sepertation * ALPHA_MAX_PML_y) then
      if (value_mean >= value_max) then
        value_mean = (value_mean + value_max) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
      else
        value_mean = (value_mean + value_max) / TWO + const_for_sepertation_three * ALPHA_MAX_PML_y
      endif
    endif
  else if (abs(value_max - value_mean) < const_for_sepertation * ALPHA_MAX_PML_y .and. &
           abs(value_mean - value_min) < const_for_sepertation * ALPHA_MAX_PML_y) then
    value_mean = (value_min + value_mean) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
    if (abs(value_max - value_mean) < const_for_sepertation * ALPHA_MAX_PML_y) then
      if (value_mean >= value_max) then
        value_mean = (value_mean + value_max) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
      else
        value_mean = (value_mean + value_max) / TWO + const_for_sepertation_three * ALPHA_MAX_PML_y
      endif
    endif
  endif
end subroutine seperate_three_sequential_value_with_only_meanvalue_changeable

subroutine seperate_three_sequential_value_with_only_minvalue_changeable(value_max,value_mean,value_min,ALPHA_MAX_PML_y)

  use generate_databases_par, only: CUSTOM_REAL,TWO
  implicit none
  real(kind=CUSTOM_REAL), intent(in) :: ALPHA_MAX_PML_y,value_max,value_mean
  real(kind=CUSTOM_REAL) :: value_min
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation = 0.009_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_two = 0.01_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_three = 0.015_CUSTOM_REAL

  if (.not. (value_max > value_mean .and. value_mean > value_min)) then
    stop 'It should be value_max > value_mean and value_mean > value_min'
  endif

  if (abs(value_max - value_mean) < const_for_sepertation * ALPHA_MAX_PML_y) then
    stop 'It should be abs(value_max - value_mean) >= 0.009_CUSTOM_REAL * ALPHA_MAX_PML_y'
  endif

  if (abs(value_mean - value_min) < const_for_sepertation * ALPHA_MAX_PML_y) then
    value_min = (value_mean + value_min) / TWO + const_for_sepertation_three * ALPHA_MAX_PML_y
  endif

  if (abs(value_max - value_min) < const_for_sepertation * ALPHA_MAX_PML_y) then
    if (value_min >= value_max) then
      value_min = (value_min + value_max) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
    else
      value_min = (value_min + value_max) / TWO + const_for_sepertation_three * ALPHA_MAX_PML_y
    endif
  endif
end subroutine seperate_three_sequential_value_with_only_minvalue_changeable

subroutine seperate_betai_alphaj_alphak(beta_x,alpha_y,alpha_z,ALPHA_MAX_PML_y)

  use generate_databases_par, only: CUSTOM_REAL,TWO
  implicit none
  real(kind=CUSTOM_REAL) :: beta_x,alpha_y,alpha_z,ALPHA_MAX_PML_y
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation = 0.009_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_two = 0.01_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: const_for_sepertation_three = 0.015_CUSTOM_REAL

  if (beta_x == alpha_y) then
    beta_x = beta_x + const_for_sepertation_two * ALPHA_MAX_PML_y
    call seperate_two_value_with_one_changeable(beta_x,alpha_z,ALPHA_MAX_PML_y)
  else if (beta_x == alpha_z) then
    beta_x = beta_x + const_for_sepertation_two * ALPHA_MAX_PML_y
    call seperate_two_value_with_one_changeable(beta_x,alpha_y,ALPHA_MAX_PML_y)
  else if (beta_x > alpha_y .and. alpha_y > alpha_z) then
    if (abs(beta_x - alpha_y) < const_for_sepertation * ALPHA_MAX_PML_y) then
      beta_x = (beta_x + alpha_y) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
    endif
  else if (beta_x > alpha_z .and. alpha_z > alpha_y) then
    if (abs(beta_x - alpha_z) < const_for_sepertation * ALPHA_MAX_PML_y) then
      beta_x = (beta_x + alpha_z) / TWO + const_for_sepertation_two * ALPHA_MAX_PML_y
    endif
  else if (alpha_z > beta_x .and. beta_x > alpha_y) then
    call seperate_three_sequential_value_with_only_meanvalue_changeable &
                  (alpha_z,beta_x,alpha_y,ALPHA_MAX_PML_y)
  else if (alpha_z > alpha_y .and. alpha_y > beta_x) then
    call seperate_three_sequential_value_with_only_minvalue_changeable &
                  (alpha_z,alpha_y,beta_x,ALPHA_MAX_PML_y)
  else if (alpha_y > beta_x .and. beta_x > alpha_z) then
    call seperate_three_sequential_value_with_only_meanvalue_changeable &
                  (alpha_y,beta_x,alpha_z,ALPHA_MAX_PML_y)
  else if (alpha_y > alpha_z .and. alpha_z > beta_x) then
    call seperate_three_sequential_value_with_only_minvalue_changeable &
                  (alpha_y,alpha_z,beta_x,ALPHA_MAX_PML_y)
  endif
end subroutine seperate_betai_alphaj_alphak


