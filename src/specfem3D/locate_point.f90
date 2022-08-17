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


!--------------------------------------------------------------------------------------------------------------------
!  locates target point in mesh.
!--------------------------------------------------------------------------------------------------------------------

  subroutine locate_point_in_mesh(x_target, y_target, z_target, &
                                  POINT_CAN_BE_BURIED, elemsize_max_glob, &
                                  ispec_selected, xi_found, eta_found, gamma_found, &
                                  x_found, y_found, z_found, domain, nu_point, final_distance_squared)

  use constants
  use specfem_par, only: ibool,myrank,NSPEC_AB,NGLOB_AB, &
                         xstore,ystore,zstore, &
                         ispec_is_surface_external_mesh,iglob_is_surface_external_mesh, &
                         xyz_midpoints

  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_poroelastic, only: ispec_is_poroelastic

  use kdtree_search, only: kdtree_search_num_nodes,kdtree_search_index, &
    kdtree_count_nearest_n_neighbors,kdtree_get_nearest_n_neighbors,kdtree_find_nearest_neighbor

  implicit none

  double precision,                      intent(in)     :: x_target, y_target, z_target
  logical,                               intent(in)     :: POINT_CAN_BE_BURIED
  real(kind=CUSTOM_REAL),                intent(in)     :: elemsize_max_glob

  double precision,                      intent(inout)    :: x_found,  y_found,  z_found
  double precision,                      intent(inout)    :: xi_found, eta_found, gamma_found
  integer,                               intent(inout)    :: ispec_selected, domain
  double precision, dimension(NDIM,NDIM),intent(inout)    :: nu_point
  double precision,                      intent(inout)    :: final_distance_squared

  ! locals
  integer :: ispec, iglob, i, j, k
  ! location search
  double precision :: maximal_elem_size_squared, dist_squared
  double precision :: distmin_squared
  double precision :: final_distance_squared_this_element
  double precision :: x,y,z
  double precision :: xi,eta,gamma

  integer :: ix_initial_guess, iy_initial_guess, iz_initial_guess
  integer :: imin,imax,jmin,jmax,kmin,kmax

!! DK DK dec 2017: also loop on all the elements in contact with the initial guess element to improve accuracy of estimate
  logical, dimension(:),allocatable :: flag_topological    ! making array allocatable, otherwise will crash for large meshes
  integer :: number_of_mesh_elements_for_the_initial_guess
  logical :: use_adjacent_elements_search

! dynamic array
!  integer, dimension(:), allocatable :: array_of_all_elements_of_ispec_selected
!  logical, parameter :: USE_SINGLE_PASS = .false.
!
! static array, max size is 50 (typically found around 18 to 24 elements)
  integer, dimension(50) :: array_of_all_elements_of_ispec_selected   ! 15% faster
  logical, parameter :: USE_SINGLE_PASS = .true.

  ! kdtree search
  ! tree nodes search
  double precision,dimension(3) :: xyz_target
  double precision :: r_search,dist_nearest
  ! search elements
  integer :: ielem,num_elem_local,ier,ispec_nearest
  logical :: use_brute_force_search
  logical :: is_better_location

  !debug
  !print *,'locate point in mesh: ',x_target, y_target, z_target

  ! sets maximal element size for search
  ! use 10 times the distance as a criterion for source detection
  maximal_elem_size_squared = (10. * elemsize_max_glob)**2

  ! INITIALIZE LOCATION --------
  ispec_selected   = 1   !! first element by default
  ix_initial_guess = 1
  iy_initial_guess = 1
  iz_initial_guess = 1

  ! puts initial xi/eta/gamma outside of element
  xi_found = 2.d0
  eta_found = 2.d0
  gamma_found = 2.d0

  ! set distance to huge initial value
  distmin_squared = HUGEVAL

  ! point search type
  if (DO_BRUTE_FORCE_POINT_SEARCH .or. (POINT_CAN_BE_BURIED .eqv. .false.)) then
    use_brute_force_search = .true.
  else
    use_brute_force_search = .false.
  endif

  if (use_brute_force_search) then
    ! brute-force search always loops over whole mesh slice
    num_elem_local = NSPEC_AB
  else
    ! kdtree search
    kdtree_search_num_nodes = 0

    ! searches target location
    xyz_target(1) = x_target
    xyz_target(2) = y_target
    xyz_target(3) = z_target

    !debug
    !print *,'find nearest neighbor',xyz_target

    ! finds closest element in slice
    call kdtree_find_nearest_neighbor(xyz_target,ispec_nearest,dist_nearest)

    !debug
    !print *,'closest element',ispec_nearest,'distance = ',dist_nearest

    ! single element for first detection loop
    ! in case point location is outside of element, it will find a closest GLL point at the element boundary
    ! and search again in adjacent elements
    num_elem_local = 1
  endif

!! DK DK dec 2017 changed this
  ! ! set which GLL points will be considered during research, depending on the possibility to bury the point or not
  ! if (.not. POINT_CAN_BE_BURIED) then
  imin  = 1
  imax  = NGLLX
  jmin  = 1
  jmax  = NGLLY
  kmin  = 1
  kmax  = NGLLZ
  ! else
  !   imin  = 2
  !   imax  = NGLLX-1
  !   jmin  = 2
  !   jmax  = NGLLY-1
  !   kmin  = 2
  !   kmax  = NGLLZ-1
  ! endif

  ! find the element candidate that may contain the target point
  do ielem = 1, num_elem_local

    ! search element
    if (use_brute_force_search) then
      ! brute-force search
      ispec = ielem
    else
      ! kd-tree search element
      ispec = ispec_nearest
    endif

    ! skip buried elements in case station must lie on surface
    if (.not. POINT_CAN_BE_BURIED .and. .not. ispec_is_surface_external_mesh(ispec)) cycle

    ! distance to element midpoint
    dist_squared = (x_target - xyz_midpoints(1,ispec))*(x_target - xyz_midpoints(1,ispec)) &
                 + (y_target - xyz_midpoints(2,ispec))*(y_target - xyz_midpoints(2,ispec)) &
                 + (z_target - xyz_midpoints(3,ispec))*(z_target - xyz_midpoints(3,ispec))

    ! we compare squared distances instead of distances to speed up the comparison by avoiding computing a square root
    if (dist_squared > maximal_elem_size_squared) cycle ! exclude elements that are too far from target

    ! find closest GLL point form target
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          iglob = ibool(i,j,k,ispec)

          ! skip inner GLL points in case station must lie on surface
          if (.not. POINT_CAN_BE_BURIED .and. .not. iglob_is_surface_external_mesh(iglob)) cycle

          x = dble(xstore(iglob))
          y = dble(ystore(iglob))
          z = dble(zstore(iglob))

          ! distance to GLL point
          dist_squared = (x_target - x)*(x_target - x) &
                       + (y_target - y)*(y_target - y) &
                       + (z_target - z)*(z_target - z)

          if (dist_squared < distmin_squared) then
            distmin_squared = dist_squared
            ispec_selected  = ispec
            ix_initial_guess = i
            iy_initial_guess = j
            iz_initial_guess = k

            x_found = x
            y_found = y
            z_found = z
          endif

        enddo
      enddo
    enddo

  enddo

  ! get the rotation matrix that will be used to rotate
  ! initializes rotation matrix
  ! East
  nu_point(1,1) = 1.d0
  nu_point(1,2) = 0.d0
  nu_point(1,3) = 0.d0
  ! North
  nu_point(2,1) = 0.d0
  nu_point(2,2) = 1.d0
  nu_point(2,3) = 0.d0
  ! Vertical
  nu_point(3,1) = 0.d0
  nu_point(3,2) = 0.d0
  nu_point(3,3) = 1.d0
  ! -- source force vector/receiver seismogram -- if the point is on the surface
  call define_station_rotation_matrix(POINT_CAN_BE_BURIED, &
                                      ix_initial_guess,iy_initial_guess,iz_initial_guess, &
                                      ispec_selected,nu_point)

  ! checks if initial guess lies on an element boundary
  if ((ix_initial_guess == 1 .or. ix_initial_guess == NGLLX) .or. &
      (iy_initial_guess == 1 .or. iy_initial_guess == NGLLY) .or. &
      (iz_initial_guess == 1 .or. iz_initial_guess == NGLLZ)) then
    ! best guess close to a GLL point on an element boundary
    ! checks if at surface
    if ((ix_initial_guess > 1 .and. ix_initial_guess < NGLLX) .and. &
        (iy_initial_guess > 1 .and. iy_initial_guess < NGLLY) .and. &
        (iz_initial_guess == NGLLZ)) then
      ! checks if point on surface, but within element
      iglob = ibool(ix_initial_guess,iy_initial_guess,iz_initial_guess,ispec_selected)
      if (iglob_is_surface_external_mesh(iglob)) then
        ! skip adjacent element
        use_adjacent_elements_search = .false.
      else
        ! not on a surface, check with top element
        use_adjacent_elements_search = .true.
      endif
    else
      ! include searching adjacent elements for better position
      use_adjacent_elements_search = .true.
    endif
  else
    ! best guess lies within the selected element, no need for checking adjacent elements
    use_adjacent_elements_search = .false.
  endif

  ! debug
  !print *,'locate ',use_adjacent_elements_search,ix_initial_guess,iy_initial_guess,iz_initial_guess

  if (use_adjacent_elements_search) then
    !! DK DK dec 2017: also loop on all the elements in contact with the initial guess element to improve accuracy of estimate

    if (use_brute_force_search) then
      ! brute-force search always loops over whole mesh slice
      num_elem_local = NSPEC_AB
    else
      ! kdtree search
      ! for adjacent search
      ! searches sphere around target location
      r_search = 0.8 * elemsize_max_glob

      ! gets number of tree points within search radius
      ! (within search sphere)
      call kdtree_count_nearest_n_neighbors(xyz_target,r_search,num_elem_local)

      ! sets n-search number of nodes
      kdtree_search_num_nodes = num_elem_local

      ! debug
      !print *,'  total number of search elements: ',num_elem_local,POINT_CAN_BE_BURIED,r_search,elemsize_max_glob

      ! allocates search array
      if (kdtree_search_num_nodes > 0) then
        allocate(kdtree_search_index(kdtree_search_num_nodes),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2418')
        if (ier /= 0) stop 'Error allocating array kdtree_search_index'

        ! finds closest n points in mesh
        ! (within search sphere)
        call kdtree_get_nearest_n_neighbors(xyz_target,r_search,num_elem_local)
      else
        ! no near target locations have been found
        ! it is outside of this slice

        !print *,'kdtree search: could not find close enough elements for adjacency'
        !print *,'  num_elem_local: ',num_elem_local,'r_search = ',r_search,'elemsize = ',elemsize_max_glob

        ! starts with dummy element
        allocate(kdtree_search_index(1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 2419')
        if (ier /= 0) stop 'Error allocating array kdtree_search_index'
        kdtree_search_index(1) = 1 ! first element

        ! dummy search in this slice
        num_elem_local = 1
        kdtree_search_num_nodes = num_elem_local
      endif
    endif

    ! flagging corners
    allocate(flag_topological(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating flag_topological array'
    flag_topological(:) = .false.

    ! mark the eight corners of the initial guess element
    flag_topological(ibool(1,1,1,ispec_selected)) = .true.
    flag_topological(ibool(NGLLX,1,1,ispec_selected)) = .true.
    flag_topological(ibool(NGLLX,NGLLY,1,ispec_selected)) = .true.
    flag_topological(ibool(1,NGLLY,1,ispec_selected)) = .true.

    flag_topological(ibool(1,1,NGLLZ,ispec_selected)) = .true.
    flag_topological(ibool(NGLLX,1,NGLLZ,ispec_selected)) = .true.
    flag_topological(ibool(NGLLX,NGLLY,NGLLZ,ispec_selected)) = .true.
    flag_topological(ibool(1,NGLLY,NGLLZ,ispec_selected)) = .true.

    if (USE_SINGLE_PASS) then
      ! assumes a maximum of a 100 adjacent elements
      !
      ! first store the initial guess itself
      number_of_mesh_elements_for_the_initial_guess = 1
      array_of_all_elements_of_ispec_selected(number_of_mesh_elements_for_the_initial_guess) = ispec_selected

      ! then store all the others
      do ielem = 1,num_elem_local
        if (use_brute_force_search) then
          ! brute-force search
          ispec = ielem
        else
          ! kd-tree search element
          ispec = kdtree_search_index(ielem)
        endif

        ! omit selected element since already included
        if (ispec == ispec_selected) cycle

        ! loop on the eight corners only, no need to loop on the rest since we just want to detect adjacency
        do k = 1,NGLLZ,NGLLZ-1
          do j = 1,NGLLY,NGLLY-1
            do i = 1,NGLLX,NGLLX-1
              if (flag_topological(ibool(i,j,k,ispec))) then
                ! this element is in contact with the initial guess
                number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
                ! check
                if (number_of_mesh_elements_for_the_initial_guess > 100) stop 'Error must increase array size in locate_point.f90'

                array_of_all_elements_of_ispec_selected(number_of_mesh_elements_for_the_initial_guess) = ispec

                ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
                goto 707

              endif
            enddo
          enddo
        enddo

707     continue

      enddo ! ielem

    else
      ! note: this involves 2-passes over all elements and dynamic memory allocation
      !       which can be very slow for large meshes

      ! loop on all the elements to count how many are shared with the initial guess
      number_of_mesh_elements_for_the_initial_guess = 1
      do ielem = 1,num_elem_local
        if (use_brute_force_search) then
          ! brute-force search
          ispec = ielem
        else
          ! kd-tree search element
          ispec = kdtree_search_index(ielem)
        endif

        if (ispec == ispec_selected) cycle

        ! loop on the eight corners only, no need to loop on the rest since we just want to detect adjacency
        do k = 1,NGLLZ,NGLLZ-1
          do j = 1,NGLLY,NGLLY-1
            do i = 1,NGLLX,NGLLX-1
              if (flag_topological(ibool(i,j,k,ispec))) then
                ! this element is in contact with the initial guess
                number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1

                ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
                goto 700

              endif
            enddo
          enddo
        enddo

700     continue

      enddo ! ielem

      ! debug
      !print *,'number_of_mesh_elements_for_the_initial_guess = ',number_of_mesh_elements_for_the_initial_guess

      ! now that we know the number of elements, we can allocate the list of elements and create it
      !allocate(array_of_all_elements_of_ispec_selected(number_of_mesh_elements_for_the_initial_guess))

      ! first store the initial guess itself
      number_of_mesh_elements_for_the_initial_guess = 1
      array_of_all_elements_of_ispec_selected(number_of_mesh_elements_for_the_initial_guess) = ispec_selected

      ! then store all the others
      do ielem = 1,num_elem_local
        if (use_brute_force_search) then
          ! brute-force search
          ispec = ielem
        else
          ! kd-tree search element
          ispec = kdtree_search_index(ielem)
        endif

        if (ispec == ispec_selected) cycle

        ! loop on the eight corners only, no need to loop on the rest since we just want to detect adjacency
        do k = 1,NGLLZ,NGLLZ-1
          do j = 1,NGLLY,NGLLY-1
            do i = 1,NGLLX,NGLLX-1
              if (flag_topological(ibool(i,j,k,ispec))) then
                ! this element is in contact with the initial guess
                number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
                array_of_all_elements_of_ispec_selected(number_of_mesh_elements_for_the_initial_guess) = ispec

                ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
                goto 800

              endif
            enddo
          enddo
        enddo

800     continue

      enddo ! ielem
    endif ! USE_SINGLE_PASS

    ! debug
    !print *,'adjacent elements ',number_of_mesh_elements_for_the_initial_guess

    ! frees search results
    if (.not. use_brute_force_search) then
      if (kdtree_search_num_nodes > 0) then
        deallocate(kdtree_search_index)
        kdtree_search_num_nodes = 0
      endif
    endif

    deallocate(flag_topological)

  else
    ! no need for adjacent element, only loop within initial guess
    number_of_mesh_elements_for_the_initial_guess = 1
    array_of_all_elements_of_ispec_selected(1) = ispec_selected
  endif

!! DK DK dec 2017
  final_distance_squared = HUGEVAL

  do i = 1,number_of_mesh_elements_for_the_initial_guess

!! DK DK dec 2017 set initial guess in the middle of the element, since we computed the true one only for the true initial guess
!! DK DK dec 2017 the nonlinear process below will converge anyway
    if (i > 1) then
      ix_initial_guess = MIDX
      iy_initial_guess = MIDY
      iz_initial_guess = MIDZ
    endif

    ispec = array_of_all_elements_of_ispec_selected(i)

    ! gets xi/eta/gamma and corresponding x/y/z coordinates
    call find_local_coordinates(x_target,y_target,z_target,xi,eta,gamma,x,y,z, &
                                ispec,ix_initial_guess,iy_initial_guess,iz_initial_guess)

    ! compute final distance squared between asked and found
    final_distance_squared_this_element = (x_target-x)*(x_target-x) + (y_target-y)*(y_target-y) + (z_target-z)*(z_target-z)

    ! flag to determine if taking new position
    is_better_location = .false.

    ! if we found a close point
    if (final_distance_squared_this_element < 1.d-10) then
      ! checks if position within element
      if (abs(xi) <= 1.d0 .and. abs(eta) <= 1.d0 .and. abs(gamma) <= 1.d0) then
        ! takes position if old position has xi/eta/gamma outside of element
        if (abs(xi_found) > 1.d0 .or. abs(eta_found) > 1.d0 .or. abs(gamma_found) > 1.d0) then
          ! old position is out of element, new inside an element; let's take the new one
          is_better_location = .true.
        endif
        ! takes position if old position is in acoustic element and new one in elastic one
        ! prefers having station in elastic elements over acoustic at interfaces
        if (ispec_is_acoustic(ispec_selected) .and. ispec_is_elastic(ispec)) then
          is_better_location = .true.
        endif
      endif
    endif

    ! if we have found an element that gives a shorter distance
    if (final_distance_squared_this_element < final_distance_squared) then
      ! check if location is better in case almost similar distances
      ! since we allow xi/eta/gamma to go a bit outside of the element range [-1.01,1.01],
      ! we keep the element where xi/eta/gamma is still within element
      !
      ! note: this can avoid problems when points on interfaces, say between acoustic/elastic, have been choosen.
      !       the point search takes whatever element comes last with a best position, even if the point
      !       has been located slightly away from the element. thus, we try to check if the previous location
      !       was accurate enough with coordinates still within an element.
      is_better_location = .true.
      if (abs(final_distance_squared_this_element - final_distance_squared) < 1.d-15) then
        if (abs(xi) > 1.d0 .or. abs(eta) > 1.d0 .or. abs(gamma) > 1.d0) then
          if (abs(xi_found) > 1.d0 .or. abs(eta_found) > 1.d0 .or. abs(gamma_found) > 1.d0) then
            ! old position is also out of element, let's take the new one
            is_better_location = .true.
          else
            ! old position is still within element and new distance has only little improved, so keep old one
            is_better_location = .false.
          endif
        endif
      endif
    endif

    ! store information about the point found
    if (is_better_location) then
      ! note: xi/eta/gamma will be in range [-1,1]
      ispec_selected = ispec

      xi_found = xi
      eta_found = eta
      gamma_found = gamma

      x_found = x
      y_found = y
      z_found = z

      !   store final distance squared between asked and found
      final_distance_squared = final_distance_squared_this_element
    endif

!! DK DK dec 2017
  enddo

  ! sets whether acoustic (1) or elastic (2)
  if (ispec_is_acoustic( ispec_selected )) then
    domain = IDOMAIN_ACOUSTIC
  else if (ispec_is_elastic( ispec_selected )) then
    domain = IDOMAIN_ELASTIC
  else if (ispec_is_poroelastic( ispec_selected )) then
    domain = IDOMAIN_POROELASTIC
  else
    domain = 0
  endif

!  deallocate(array_of_all_elements_of_ispec_selected)

  end subroutine locate_point_in_mesh

!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_local_coordinates(x_target,y_target,z_target,xi,eta,gamma,x,y,z, &
                                    ispec,ix_initial_guess,iy_initial_guess,iz_initial_guess)

  use constants, only: HUGEVAL,NUM_ITER
  use specfem_par, only: xstore,ystore,zstore,ibool,NGNOD,anchor_iax,anchor_iay,anchor_iaz, &
                         xigll,yigll,zigll

  implicit none

  double precision,intent(in) :: x_target,y_target,z_target
  double precision,intent(out) :: xi,eta,gamma
  double precision,intent(out) :: x,y,z

  integer,intent(in) :: ispec,ix_initial_guess,iy_initial_guess,iz_initial_guess

  ! local parameters
  integer :: ia,iter_loop
  integer :: iglob

  ! coordinates of the control points of the surface element
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  double precision :: dx,dy,dz,dx_min,dy_min,dz_min,d_sq,d_min_sq
  double precision :: dxi,deta,dgamma
  double precision :: xix,xiy,xiz
  double precision :: etax,etay,etaz
  double precision :: gammax,gammay,gammaz

  ! initializes
  x = 0.d0
  y = 0.d0
  z = 0.d0

  ! define coordinates of the control points of the element
  do ia = 1,NGNOD
    iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
    xelm(ia) = dble(xstore(iglob))
    yelm(ia) = dble(ystore(iglob))
    zelm(ia) = dble(zstore(iglob))
  enddo

  ! use initial guess in xi and eta
  xi = xigll(ix_initial_guess)
  eta = yigll(iy_initial_guess)
  gamma = zigll(iz_initial_guess)

  d_min_sq = HUGEVAL
  dx_min = HUGEVAL
  dy_min = HUGEVAL
  dz_min = HUGEVAL

  ! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

    ! recompute Jacobian for the new point
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                            xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,NGNOD)

    ! compute distance to target location
    dx = - (x - x_target)
    dy = - (y - y_target)
    dz = - (z - z_target)

    !debug
    !print *,'  iter ',iter_loop,'dx',sngl(dx),sngl(dx_min),'dy',sngl(dy),sngl(dy_min),'dz',sngl(dz),sngl(dz_min),d_min_sq

    ! distance squared
    d_sq = dx*dx + dy*dy + dz*dz

    ! compute increments
    if (d_sq < d_min_sq) then
      d_min_sq = d_sq
      dx_min = dx
      dy_min = dy
      dz_min = dz

      dxi = xix*dx + xiy*dy + xiz*dz
      deta = etax*dx + etay*dy + etaz*dz
      dgamma = gammax*dx + gammay*dy + gammaz*dz
    else
      ! new position is worse than old one, no change necessary
      dxi = 0.d0
      deta = 0.d0
      dgamma = 0.d0
    endif

    ! decreases step length if step is large
    if ((dxi*dxi + deta*deta + dgamma*dgamma) > 1.0d0) then
      dxi = dxi * 0.33333333333d0
      deta = deta * 0.33333333333d0
      dgamma = dgamma * 0.33333333333d0
    endif
    ! alternative: impose limit on increments (seems to result in slightly less accurate locations)
    !if (abs(dxi) > 0.3d0 ) dxi = sign(1.0d0,dxi)*0.3d0
    !if (abs(deta) > 0.3d0 ) deta = sign(1.0d0,deta)*0.3d0
    !if (abs(dgamma) > 0.3d0 ) dgamma = sign(1.0d0,dgamma)*0.3d0

    !debug
    !print *,'  dxi/..',(dxi**2 + deta**2 + dgamma**2),dxi,deta,dgamma

    ! update values
    xi = xi + dxi
    eta = eta + deta
    gamma = gamma + dgamma

    ! impose that we stay in that element
    ! (useful if user gives a receiver outside the mesh for instance)
    ! we can go slightly outside the [1,1] segment since with finite elements
    ! the polynomial solution is defined everywhere
    ! can be useful for convergence of iterative scheme with distorted elements
    !
    ! note: we are more conservative here than in the 3D_GLOBE version, since the mesh can have voids or complex shapes,
    !       and we want the point location not to be too far off the mesh volume
    if (xi > 1.01d0) xi     =  1.01d0
    if (xi < -1.01d0) xi     = -1.01d0
    if (eta > 1.01d0) eta    =  1.01d0
    if (eta < -1.01d0) eta    = -1.01d0
    if (gamma > 1.01d0) gamma  =  1.01d0
    if (gamma < -1.01d0) gamma  = -1.01d0

  ! end of non linear iterations
  enddo

  ! compute final coordinates of point found
  call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,NGNOD)

  end subroutine find_local_coordinates


!--------------------------------------------------------------------------------------------------------------------
!  Define the rotation matrix in the selected point
!--------------------------------------------------------------------------------------------------------------------

  subroutine define_station_rotation_matrix(POINT_CAN_BE_BURIED,ix_initial_guess,iy_initial_guess,iz_initial_guess, &
                                            ispec_selected,nu_point)

  use constants
  use specfem_par, only: ibool,xstore,ystore,zstore,iglob_is_surface_external_mesh

  logical,                                intent(in)  :: POINT_CAN_BE_BURIED
  integer,                                intent(in)  :: ix_initial_guess,iy_initial_guess,iz_initial_guess,ispec_selected
  double precision, dimension(NDIM,NDIM), intent(inout) :: nu_point

  ! local parameters
  ! for surface locating and normal computing with external mesh
  real(kind=CUSTOM_REAL), dimension(NDIM) :: u_vector,v_vector,w_vector
  integer                                 :: pt0_ix,pt0_iy,pt0_iz,pt1_ix,pt1_iy,pt1_iz,pt2_ix,pt2_iy,pt2_iz

  ! Rotation matrix is diagonal if the point is inside the mesh, or if decided by user
  if (POINT_CAN_BE_BURIED .or. (.not. EXTERNAL_MESH_RECEIVERS_NORMAL)) then
    ! East
    nu_point(1,1) = 1.d0
    nu_point(1,2) = 0.d0
    nu_point(1,3) = 0.d0
    ! North
    nu_point(2,1) = 0.d0
    nu_point(2,2) = 1.d0
    nu_point(2,3) = 0.d0
    ! Vertical
    nu_point(3,1) = 0.d0
    nu_point(3,2) = 0.d0
    nu_point(3,3) = 1.d0

  else
    ! note: the following computes the normal at the corner of the element.
    !       this normal is still fine for linear elements, however, for quadratic elements (NGNOD == 27)
    !       there might be slight differences to the actual normal at the target point.

    ! get normal to the face of the hexahedra if receiver is on the surface
    pt0_ix = -1
    pt0_iy = -1
    pt0_iz = -1
    pt1_ix = -1
    pt1_iy = -1
    pt1_iz = -1
    pt2_ix = -1
    pt2_iy = -1
    pt2_iz = -1
    ! we get two vectors of the face (three points) to compute the normal
    if (ix_initial_guess == 1 .and. iglob_is_surface_external_mesh(ibool(1,2,2,ispec_selected))) then
      pt0_ix = 1
      pt0_iy = NGLLY
      pt0_iz = 1
      pt1_ix = 1
      pt1_iy = 1
      pt1_iz = 1
      pt2_ix = 1
      pt2_iy = NGLLY
      pt2_iz = NGLLZ
    endif
    if (ix_initial_guess == NGLLX .and. iglob_is_surface_external_mesh(ibool(NGLLX,2,2,ispec_selected))) then
      pt0_ix = NGLLX
      pt0_iy = 1
      pt0_iz = 1
      pt1_ix = NGLLX
      pt1_iy = NGLLY
      pt1_iz = 1
      pt2_ix = NGLLX
      pt2_iy = 1
      pt2_iz = NGLLZ
    endif

    if (iy_initial_guess == 1 .and. iglob_is_surface_external_mesh(ibool(2,1,2,ispec_selected))) then
      pt0_ix = 1
      pt0_iy = 1
      pt0_iz = 1
      pt1_ix = NGLLX
      pt1_iy = 1
      pt1_iz = 1
      pt2_ix = 1
      pt2_iy = 1
      pt2_iz = NGLLZ
    endif
    if (iy_initial_guess == NGLLY .and. iglob_is_surface_external_mesh(ibool(2,NGLLY,2,ispec_selected))) then
      pt0_ix = NGLLX
      pt0_iy = NGLLY
      pt0_iz = 1
      pt1_ix = 1
      pt1_iy = NGLLY
      pt1_iz = 1
      pt2_ix = NGLLX
      pt2_iy = NGLLY
      pt2_iz = NGLLZ
    endif

    if (iz_initial_guess == 1 .and. iglob_is_surface_external_mesh(ibool(2,2,1,ispec_selected))) then
      pt0_ix = NGLLX
      pt0_iy = 1
      pt0_iz = 1
      pt1_ix = 1
      pt1_iy = 1
      pt1_iz = 1
      pt2_ix = NGLLX
      pt2_iy = NGLLY
      pt2_iz = 1
    endif
    if (iz_initial_guess == NGLLZ .and. iglob_is_surface_external_mesh(ibool(2,2,NGLLZ,ispec_selected))) then
      pt0_ix = 1
      pt0_iy = 1
      pt0_iz = NGLLZ
      pt1_ix = NGLLX
      pt1_iy = 1
      pt1_iz = NGLLZ
      pt2_ix = 1
      pt2_iy = NGLLY
      pt2_iz = NGLLZ
    endif

    if (pt0_ix < 0 .or. pt0_iy < 0 .or. pt0_iz < 0 .or. &
        pt1_ix < 0 .or. pt1_iy < 0 .or. pt1_iz < 0 .or. &
        pt2_ix < 0 .or. pt2_iy < 0 .or. pt2_iz < 0) then
      stop 'error in computing normal for receivers.'
    endif

    u_vector(1) = xstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected)) &
                - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))
    u_vector(2) = ystore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected)) &
                - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))
    u_vector(3) = zstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected)) &
                - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))

    v_vector(1) = xstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected)) &
                - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))
    v_vector(2) = ystore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected)) &
                - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))
    v_vector(3) = zstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected)) &
                - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))

    ! cross product
    w_vector(1) = u_vector(2)*v_vector(3) - u_vector(3)*v_vector(2)
    w_vector(2) = u_vector(3)*v_vector(1) - u_vector(1)*v_vector(3)
    w_vector(3) = u_vector(1)*v_vector(2) - u_vector(2)*v_vector(1)

    ! normalize vector w
    w_vector(:) = w_vector(:)/sqrt(w_vector(1)**2+w_vector(2)**2+w_vector(3)**2)

    ! build the two other vectors for a direct base: we normalize u, and v=w^u
    u_vector(:) = u_vector(:)/sqrt(u_vector(1)**2+u_vector(2)**2+u_vector(3)**2)
    v_vector(1) = w_vector(2)*u_vector(3) - w_vector(3)*u_vector(2)
    v_vector(2) = w_vector(3)*u_vector(1) - w_vector(1)*u_vector(3)
    v_vector(3) = w_vector(1)*u_vector(2) - w_vector(2)*u_vector(1)

    ! build rotation matrice nu
    ! East (u)
    nu_point(1,1) = u_vector(1)
    nu_point(1,2) = v_vector(1)
    nu_point(1,3) = w_vector(1)
    ! North (v)
    nu_point(2,1) = u_vector(2)
    nu_point(2,2) = v_vector(2)
    nu_point(2,3) = w_vector(2)
    ! Vertical (w)
    nu_point(3,1) = u_vector(3)
    nu_point(3,2) = v_vector(3)
    nu_point(3,3) = w_vector(3)

  endif

  end subroutine define_station_rotation_matrix
