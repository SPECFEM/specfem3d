!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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


  subroutine pml_set_local_dampingcoeff(xstore,ystore,zstore)

  ! calculates damping profiles and auxiliary coefficients on C-PML points

  use constants, only: myrank,ZERO,ONE,TWO,HUGEVAL

  use shared_parameters, only: ACOUSTIC_SIMULATION, ELASTIC_SIMULATION

  use generate_databases_par, only: ibool,NGLOB_AB,d_store_x,d_store_y,d_store_z, &
                                    K_store_x,K_store_y,K_store_z,alpha_store_x,alpha_store_y,alpha_store_z,CPML_to_spec, &
                                    CPML_width_x,CPML_width_y,CPML_width_z,min_distance_between_CPML_parameter, &
                                    CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,nspec_cpml,PML_INSTEAD_OF_FREE_SURFACE, &
                                    IMAIN,CPML_REGIONS,f0_FOR_PML,PI, &
                                    CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ, &
                                    SIMULATION_TYPE,SAVE_FORWARD,nspec => NSPEC_AB,is_CPML, &
                                    mask_ibool_interior_domain,nglob_interface_PML_acoustic,points_interface_PML_acoustic, &
                                    nglob_interface_PML_elastic,points_interface_PML_elastic

  use create_regions_mesh_ext_par, only: rhostore,rho_vp,ispec_is_acoustic,ispec_is_elastic

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: xstore,ystore,zstore

  ! local parameters
  integer :: i,j,k,ispec,iglob,ispec_CPML,ier

  real(kind=CUSTOM_REAL) :: ALPHA_MAX_PML_x,ALPHA_MAX_PML_y,ALPHA_MAX_PML_z
  real(kind=CUSTOM_REAL), parameter :: K_MAX_PML = ONE,K_MIN_PML= ONE
  real(kind=CUSTOM_REAL) :: pml_damping_profile_l,dist,vp
  real(kind=CUSTOM_REAL) :: xoriginleft,xoriginright,yoriginfront,yoriginback,zoriginbottom,zorigintop
  real(kind=CUSTOM_REAL) :: abscissa_in_PML_x,abscissa_in_PML_y,abscissa_in_PML_z
  real(kind=CUSTOM_REAL) :: d_x,d_y,d_z,k_x,k_y,k_z,alpha_x,alpha_y,alpha_z,beta_x,beta_y,beta_z
  real(kind=CUSTOM_REAL) :: x_min,x_min_all,y_min,y_min_all,z_min,z_min_all, &
                            x_max,x_max_all,y_max,y_max_all,z_max,z_max_all, &
                            x_origin,y_origin,z_origin
  real(kind=CUSTOM_REAL) :: CPML_width_x_left, CPML_width_x_right, &
                            CPML_width_y_front,CPML_width_y_back, &
                            CPML_width_z_top,CPML_width_z_bottom, &
                            CPML_x_left, CPML_x_right, &
                            CPML_y_front,CPML_y_back, &
                            CPML_z_top,CPML_z_bottom, &
                            CPML_width_x_left_max_all, CPML_width_x_right_max_all, &
                            CPML_width_y_front_max_all,CPML_width_y_back_max_all, &
                            CPML_width_z_top_max_all,CPML_width_z_bottom_max_all, &
                            vp_max,vp_max_all

! for robust parameter separation of PML damping parameter
  real(kind=CUSTOM_REAL) :: distance_min,distance_min_glob, &
                            const_for_separation_two,const_for_separation_four,maxtemp,mintemp, &
                            min_distance_between_CPML_parameter_glob
  real(kind=CUSTOM_REAL) :: x1,x2,y1,y2,z1,z2
  integer :: iglob1,iglob2
! for robust parameter separation of PML damping parameter

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
  if (nspec_cpml > 0) then
    allocate(d_store_x(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 841')
    if (ier /= 0) stop 'error allocating array d_store_x'
    allocate(d_store_y(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 842')
    if (ier /= 0) stop 'error allocating array d_store_y'
    allocate(d_store_z(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 843')
    if (ier /= 0) stop 'error allocating array d_store_z'

    ! stores auxiliary coefficients
    allocate(K_store_x(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 844')
    if (ier /= 0) stop 'error allocating array K_store_x'
    allocate(K_store_y(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 845')
    if (ier /= 0) stop 'error allocating array K_store_y'
    allocate(K_store_z(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 846')
    if (ier /= 0) stop 'error allocating array K_store_z'
    allocate(alpha_store_x(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 847')
    if (ier /= 0) stop 'error allocating array alpha_store_x'
    allocate(alpha_store_y(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 848')
    if (ier /= 0) stop 'error allocating array alpha_store_y'
    allocate(alpha_store_z(NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 849')
    if (ier /= 0) stop 'error allocating array alpha_store_z'
  else
    ! dummy allocations
    allocate(d_store_x(1,1,1,1),d_store_y(1,1,1,1),d_store_z(1,1,1,1), &
             K_store_x(1,1,1,1),K_store_y(1,1,1,1),K_store_z(1,1,1,1), &
             alpha_store_x(1,1,1,1),alpha_store_y(1,1,1,1),alpha_store_z(1,1,1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy arrays'
  endif
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

  x_min_all = HUGEVAL
  y_min_all = HUGEVAL
  z_min_all = HUGEVAL

  x_max_all = - HUGEVAL
  y_max_all = - HUGEVAL
  z_max_all = - HUGEVAL

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

      if (CPML_regions(ispec_CPML) == CPML_X_ONLY .or. CPML_regions(ispec_CPML) == CPML_XY_ONLY .or. &
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

      if (CPML_regions(ispec_CPML) == CPML_Y_ONLY .or. CPML_regions(ispec_CPML) == CPML_XY_ONLY .or. &
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
        CPML_regions(ispec_CPML) == CPML_XZ_ONLY .or. CPML_regions(ispec_CPML) == CPML_XYZ) then
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
    write(IMAIN,*) 'Boundary values of X-/Y-/Z-regions:'
    write(IMAIN,*) '  X: ',x_min_all, x_max_all
    write(IMAIN,*) '  Y: ',y_min_all, y_max_all
    write(IMAIN,*) '  Z: ',z_min_all, z_max_all
    write(IMAIN,*)
    write(IMAIN,*) '  Origins of left/right X-surface C-PML',xoriginleft,xoriginright
    write(IMAIN,*) '  Origins of back/front Y-surface C-PML',yoriginback,yoriginfront
    write(IMAIN,*) '  Origin of bottom Z-surface C-PML     ',zoriginbottom
    if (PML_INSTEAD_OF_FREE_SURFACE) then
      write(IMAIN,*) '  Origin of top Z-surface C-PML        ',zorigintop
    endif
    write(IMAIN,*)
    write(IMAIN,*) '  CPML_width_x: ',CPML_width_x
    write(IMAIN,*) '  CPML_width_y: ',CPML_width_y
    write(IMAIN,*) '  CPML_width_z: ',CPML_width_z
    write(IMAIN,*)
    write(IMAIN,*) '  maximum Vp in C-PML',vp_max_all
    write(IMAIN,*)
    call flush_IMAIN()
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

                if (d_y < ZERO .or. K_y < ONE) then
                  K_y = ONE; d_y = ZERO
                endif

                if (d_z < ZERO .or. K_z < ONE) then
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
            if (xstore(iglob) - x_origin > ZERO .and. ystore(iglob) - y_origin > ZERO .and. &
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

! for robust parameter separation of PML damping parameter
  distance_min = HUGEVAL
  distance_min_glob = HUGEVAL
  min_distance_between_CPML_parameter_glob = HUGEVAL
  do ispec_CPML = 1,nspec_cpml
    ispec = CPML_to_spec(ispec_CPML)
  ! loops over all GLL points
  ! (combines directions to speed up calculations)
    do k=1,NGLLZ-1
      do j=1,NGLLY-1
        do i=1,NGLLX-1
          ! reference point
          iglob1 = ibool(i,j,k,ispec)
          x1 = xstore(iglob1)
          y1 = ystore(iglob1)
          z1 = zstore(iglob1)

          ! along X
          iglob2 = ibool(i+1,j,k,ispec)
          x2 = xstore(iglob2)
          y2 = ystore(iglob2)
          z2 = zstore(iglob2)

          dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

          if (dist < distance_min) distance_min = dist

          ! along Y
          iglob2 = ibool(i,j+1,k,ispec)
          x2 = xstore(iglob2)
          y2 = ystore(iglob2)
          z2 = zstore(iglob2)

          dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

          if (dist < distance_min) distance_min = dist

          ! along Z
          iglob2 = ibool(i,j,k+1,ispec)
          x2 = xstore(iglob2)
          y2 = ystore(iglob2)
          z2 = zstore(iglob2)

          dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)

          if (dist < distance_min) distance_min = dist

        enddo
      enddo
    enddo
  enddo
  distance_min = sqrt(distance_min)
  call min_all_all_cr(distance_min,distance_min_glob)
  if (myrank == 0) then
    if (distance_min_glob <= 0.0_CUSTOM_REAL) call exit_mpi(myrank,"error: GLL points minimum distance")
  endif
  min_distance_between_CPML_parameter = ALPHA_MAX_PML_x * distance_min_glob / &
                                        max(CPML_width_x,CPML_width_y,CPML_width_z) / 8._CUSTOM_REAL
  call min_all_all_cr(min_distance_between_CPML_parameter,min_distance_between_CPML_parameter_glob)
  min_distance_between_CPML_parameter = min_distance_between_CPML_parameter_glob
  const_for_separation_two = min_distance_between_CPML_parameter * 2._CUSTOM_REAL
  const_for_separation_four = min_distance_between_CPML_parameter * 4._CUSTOM_REAL
  min_distance_between_CPML_parameter = min_distance_between_CPML_parameter

! for robust parameter separation of PML damping parameter

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

            if (abs(alpha_x - alpha_y) < min_distance_between_CPML_parameter) then
              call seperate_two_changeable_value(alpha_x,alpha_y,const_for_separation_two)
            endif

            beta_x = alpha_x + d_x / K_x
            beta_y = alpha_y + d_y / K_y

            if (abs(beta_x - alpha_y) < min_distance_between_CPML_parameter) then
              call seperate_two_value_with_one_changeable(beta_x,alpha_y,const_for_separation_two, &
                                                          const_for_separation_four)
            endif

            if (abs(beta_y - alpha_x) < min_distance_between_CPML_parameter) then
              call seperate_two_value_with_one_changeable(beta_y,alpha_x,const_for_separation_two, &
                                                          const_for_separation_four)
            endif

            if (abs(alpha_x - alpha_y) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of alpha_x, alpha_y'
            endif

            if (abs(beta_x - alpha_y) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of beta_x and alpha_y'
            endif

            if (abs(beta_y - alpha_x) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of beta_y and alpha_x'
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

            if (abs(alpha_x - alpha_z) < min_distance_between_CPML_parameter) then
              call seperate_two_changeable_value(alpha_x,alpha_z,const_for_separation_two)
            endif

            beta_x = alpha_x + d_x / K_x
            beta_z = alpha_z + d_z / K_z

            if (abs(beta_x - alpha_z) < min_distance_between_CPML_parameter) then
              call seperate_two_value_with_one_changeable(beta_x,alpha_z,const_for_separation_two, &
                                                          const_for_separation_four)
            endif

            if (abs(beta_z - alpha_x) < min_distance_between_CPML_parameter) then
              call seperate_two_value_with_one_changeable(beta_z,alpha_x,const_for_separation_two, &
                                                          const_for_separation_four)
            endif

            if (abs(alpha_x - alpha_z) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of alpha_x, alpha_z'
            endif

            if (abs(beta_x - alpha_z) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of beta_x and alpha_z'
            endif

            if (abs(beta_z - alpha_x) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of beta_z and alpha_z'
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

            if (abs(alpha_y - alpha_z) < min_distance_between_CPML_parameter) then
              call seperate_two_changeable_value(alpha_y,alpha_z,const_for_separation_two)
            endif

            beta_y = alpha_y + d_y / K_y
            beta_z = alpha_z + d_z / K_z

            if (abs(beta_y - alpha_z) < min_distance_between_CPML_parameter) then
              call seperate_two_value_with_one_changeable(beta_y,alpha_z,const_for_separation_two, &
                                                          const_for_separation_four)
            endif

            if (abs(beta_z - alpha_y) < min_distance_between_CPML_parameter) then
              call seperate_two_value_with_one_changeable(beta_z,alpha_y,const_for_separation_two, &
                                                          const_for_separation_four)
            endif

            if (abs(alpha_y - alpha_z) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of alpha_y, alpha_z'
            endif

            if (abs(beta_y - alpha_z) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of beta_y and alpha_z'
            endif

            if (abs(beta_z - alpha_y) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of beta_z and alpha_y'
            endif

            d_y = (beta_y - alpha_y) * K_y
            d_z = (beta_z - alpha_z) * K_z

            d_store_y(i,j,k,ispec_CPML) = d_y
            alpha_store_y(i,j,k,ispec_CPML) = alpha_y

            d_store_z(i,j,k,ispec_CPML) = d_z
            alpha_store_z(i,j,k,ispec_CPML) = alpha_z

          else if (CPML_regions(ispec_CPML) == CPML_XYZ) then

            K_x = K_store_x(i,j,k,ispec_CPML)
            d_x = d_store_x(i,j,k,ispec_CPML)
            alpha_x = alpha_store_x(i,j,k,ispec_CPML)

            K_y = K_store_y(i,j,k,ispec_CPML)
            d_y = d_store_y(i,j,k,ispec_CPML)
            alpha_y = alpha_store_y(i,j,k,ispec_CPML)

            K_z = K_store_z(i,j,k,ispec_CPML)
            d_z = d_store_z(i,j,k,ispec_CPML)
            alpha_z = alpha_store_z(i,j,k,ispec_CPML)

            if (abs(alpha_x - alpha_y) < min_distance_between_CPML_parameter) then
              if (alpha_x > alpha_y) then
                alpha_x = alpha_y + const_for_separation_two
              else
                alpha_y = alpha_x + const_for_separation_two
              endif
              maxtemp = max(alpha_x, alpha_y)
              mintemp = min(alpha_x, alpha_y)
              if (alpha_z > maxtemp) then
                 if (abs(alpha_z - maxtemp) < min_distance_between_CPML_parameter) then
                   alpha_z = maxtemp + const_for_separation_two
                 endif
              else if (alpha_z < mintemp) then
                 if (abs(alpha_z - mintemp) < min_distance_between_CPML_parameter) then
                   if (alpha_x > alpha_y) then
                     alpha_x = alpha_z + const_for_separation_four
                     alpha_y = alpha_z + const_for_separation_two
                   else
                     alpha_y = alpha_z + const_for_separation_four
                     alpha_x = alpha_z + const_for_separation_two
                   endif
                 endif
              else
                 if (alpha_x > alpha_y) then
                   alpha_x = alpha_y + const_for_separation_four
                   alpha_z = alpha_y + const_for_separation_two
                 else
                   alpha_y = alpha_x + const_for_separation_four
                   alpha_z = alpha_x + const_for_separation_two
                 endif
              endif
            endif

            if (abs(alpha_x - alpha_z) < min_distance_between_CPML_parameter) then
              if (alpha_x > alpha_z) then
                alpha_x = alpha_z + const_for_separation_two
              else
                alpha_z = alpha_x + const_for_separation_two
              endif
              maxtemp = max(alpha_x, alpha_z)
              mintemp = min(alpha_x, alpha_z)
              if (alpha_y > maxtemp) then
                 if (abs(alpha_y - maxtemp) < min_distance_between_CPML_parameter) then
                   alpha_y = maxtemp + const_for_separation_two
                 endif
              else if (alpha_y < mintemp) then
                 if (abs(alpha_y - mintemp) < min_distance_between_CPML_parameter) then
                   if (alpha_x > alpha_z) then
                     alpha_x = alpha_y + const_for_separation_four
                     alpha_z = alpha_y + const_for_separation_two
                   else
                     alpha_z = alpha_y + const_for_separation_four
                     alpha_x = alpha_y + const_for_separation_two
                   endif
                 endif
              else
                 if (alpha_x > alpha_z) then
                   alpha_x = alpha_z + const_for_separation_four
                   alpha_y = alpha_z + const_for_separation_two
                 else
                   alpha_z = alpha_x + const_for_separation_four
                   alpha_y = alpha_x + const_for_separation_two
                 endif
              endif
            endif

            if (abs(alpha_y - alpha_z) < min_distance_between_CPML_parameter) then
              if (alpha_y > alpha_z) then
                alpha_y = alpha_z + const_for_separation_two
              else
                alpha_z = alpha_y + const_for_separation_two
              endif
              maxtemp = max(alpha_y, alpha_z)
              mintemp = min(alpha_y, alpha_z)
              if (alpha_x > maxtemp) then
                 if (abs(alpha_x - maxtemp) < min_distance_between_CPML_parameter) then
                   alpha_x = maxtemp + const_for_separation_two
                 endif
              else if (alpha_x < mintemp) then
                 if (abs(alpha_x - mintemp) < min_distance_between_CPML_parameter) then
                   if (alpha_y > alpha_z) then
                     alpha_y = alpha_x + const_for_separation_four
                     alpha_z = alpha_x + const_for_separation_two
                   else
                     alpha_z = alpha_x + const_for_separation_four
                     alpha_y = alpha_x + const_for_separation_two
                   endif
                 endif
              else
                 if (alpha_y > alpha_z) then
                   alpha_y = alpha_z + const_for_separation_four
                   alpha_x = alpha_z + const_for_separation_two
                 else
                   alpha_z = alpha_y + const_for_separation_four
                   alpha_x = alpha_y + const_for_separation_two
                 endif
              endif
            endif

            if (abs(alpha_x - alpha_y) < min_distance_between_CPML_parameter .or. &
                abs(alpha_y - alpha_z) < min_distance_between_CPML_parameter .or. &
                abs(alpha_x - alpha_z) < min_distance_between_CPML_parameter) then
              stop 'error in separation of alpha_x, alpha_y, alpha_z'
            endif

            beta_x = alpha_x + d_x / K_x
            maxtemp = max(alpha_y, alpha_z)
            mintemp = min(alpha_y, alpha_z)
            if (beta_x > maxtemp) then
               if (abs(beta_x- maxtemp) < min_distance_between_CPML_parameter) then
                 beta_x = maxtemp + const_for_separation_two
               endif
            else if (beta_x < mintemp) then
               if (abs(beta_x - mintemp) < min_distance_between_CPML_parameter) then
                 if (alpha_y > alpha_z) then
                   beta_x = alpha_y + const_for_separation_two
                 else
                   beta_x = alpha_z + const_for_separation_two
                 endif
               endif
            else
               if (abs(beta_x - maxtemp) < min_distance_between_CPML_parameter) then
                  beta_x = maxtemp + const_for_separation_two
               endif

               if (abs(beta_x - mintemp) < min_distance_between_CPML_parameter) then
                  beta_x = mintemp + const_for_separation_two
                  if (abs(beta_x - maxtemp) < min_distance_between_CPML_parameter) then
                     beta_x = maxtemp + const_for_separation_two
                  endif
               endif
            endif

            beta_y = alpha_y + d_y / K_y

            maxtemp = max(alpha_x, alpha_z)
            mintemp = min(alpha_x, alpha_z)
            if (beta_y > maxtemp) then
               if (abs(beta_y - maxtemp) < min_distance_between_CPML_parameter) then
                 beta_y = maxtemp + const_for_separation_two
               endif
            else if (beta_y < mintemp) then
               if (abs(beta_y - mintemp) < min_distance_between_CPML_parameter) then
                 if (alpha_x > alpha_z) then
                   beta_y = alpha_x + const_for_separation_two
                 else
                   beta_y = alpha_z + const_for_separation_two
                 endif
               endif
            else
               if (abs(beta_y - maxtemp) < min_distance_between_CPML_parameter) then
                  beta_y = maxtemp + const_for_separation_two
               endif

               if (abs(beta_y - mintemp) < min_distance_between_CPML_parameter) then
                  beta_y = mintemp + const_for_separation_two
                  if (abs(beta_y - maxtemp) < min_distance_between_CPML_parameter) then
                     beta_y = maxtemp + const_for_separation_two
                  endif
               endif
            endif

            beta_z = alpha_z + d_z / K_z
            maxtemp = max(alpha_x, alpha_y)
            mintemp = min(alpha_x, alpha_y)
            if (beta_z > maxtemp) then
               if (abs(beta_z - maxtemp) < min_distance_between_CPML_parameter) then
                 beta_z = maxtemp + const_for_separation_two
               endif
            else if (beta_z < mintemp) then
               if (abs(beta_z - mintemp) < min_distance_between_CPML_parameter) then
                 if (alpha_x > alpha_y) then
                   beta_z = alpha_x + const_for_separation_two
                 else
                   beta_z = alpha_y + const_for_separation_two
                 endif
               endif
            else
               if (abs(beta_z - maxtemp) < min_distance_between_CPML_parameter) then
                  beta_z = maxtemp + const_for_separation_two
               endif

               if (abs(beta_z - mintemp) < min_distance_between_CPML_parameter) then
                  beta_z = mintemp + const_for_separation_two
                  if (abs(beta_z - maxtemp) < min_distance_between_CPML_parameter) then
                     beta_z = maxtemp + const_for_separation_two
                  endif
               endif
            endif

            if (abs(beta_x - alpha_y) < min_distance_between_CPML_parameter .or. &
                abs(beta_x - alpha_z) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of beta_x,alpha_y,alpha_z '
            endif

            if (abs(beta_y - alpha_x) < min_distance_between_CPML_parameter .or. &
                abs(beta_y - alpha_z) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of beta_y, alpha_x,alpha_z '
            endif

            if (abs(beta_z - alpha_x) < min_distance_between_CPML_parameter .or. &
                abs(beta_z - alpha_y) < min_distance_between_CPML_parameter) then
              stop 'there is an error in the separation of beta_z,alpha_x,alpha_y '
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
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 850')
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
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 851')
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
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 852')
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

  use generate_databases_par, only: CUSTOM_REAL,SIZE_REAL,NPOWER,CPML_Rcoef,TWO

  implicit none

  integer, intent(in) :: myrank,iglob

  real(kind=CUSTOM_REAL), intent(in) :: dist,vp,delta

  real(kind=CUSTOM_REAL) :: pml_damping_profile_l

  ! gets damping profile
  if (NPOWER >= 1) then
     ! In INRIA research report section 6.1:  http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
     ! pml_damping_profile_l = - ((NPOWER + 1) * vp * log(CPML_Rcoef) / (TWO * delta)) * dist**(NPOWER)
     ! based on tests it is more accurate to use the following definition when NPOWER = 1 as defined in constants.h.in
    if (CUSTOM_REAL == SIZE_REAL) then
      pml_damping_profile_l = - sngl(((NPOWER + 1.d0) * dble(vp) * log(CPML_Rcoef) / &
                                          (TWO * dble(delta))) * dble(dist)**(1.2d0 * NPOWER))
    else
      pml_damping_profile_l = - ((NPOWER + 1.d0) * vp * log(CPML_Rcoef) / (TWO * delta)) * dist**(1.2d0 * NPOWER)
    endif
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

!
!-------------------------------------------------------------------------------------------------
!

  subroutine seperate_two_changeable_value(value_a,value_b,const_for_separation_two)

  use generate_databases_par, only: CUSTOM_REAL

  implicit none
  real(kind=CUSTOM_REAL), intent(in)  :: const_for_separation_two
  real(kind=CUSTOM_REAL) :: value_a,value_b
  if (value_a >= value_b) then
     value_a = value_b + const_for_separation_two
  else
     value_b = value_a + const_for_separation_two
  endif

  end subroutine seperate_two_changeable_value

!
!-------------------------------------------------------------------------------------------------
!

  subroutine seperate_two_value_with_one_changeable(value_a,value_b,const_for_separation_two, &
                                                    const_for_separation_four)

  use generate_databases_par, only: CUSTOM_REAL

  implicit none
  real(kind=CUSTOM_REAL), intent(in) :: value_b, const_for_separation_two, const_for_separation_four
  real(kind=CUSTOM_REAL) :: value_a

  if (value_a >= value_b) then
    value_a = value_b + const_for_separation_two
  else
    value_a = value_b + const_for_separation_four
  endif

  end subroutine seperate_two_value_with_one_changeable
