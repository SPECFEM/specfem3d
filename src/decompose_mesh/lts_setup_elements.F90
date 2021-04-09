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


  subroutine lts_setup_elements()

! puts elements into separate refinement levels for local time stepping

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,COURANT_SUGGESTED, &
    LTS_SINGLE_P_LEVEL,LTS_TWO_P_LEVEL,LTS_OVERLAP_REGION, &
    LTS_DECREASE_DT,LTS_STABILITY_MARGIN_DT,LTS_SAFETY_MARGIN

  use decompose_mesh_par, only: elmnts,nnodes,nodes_coords,nspec,DT,NGNOD,LOCAL_PATH, &
    ispec_p_refine, num_p_level, p_level, num_ispec_level, LTS_MODE

  implicit none

  ! local parameters
  ! ratio of minimum distance of GLL points, depending on NGLL
  double precision,dimension(15) :: percent_GLL

  ! velocity, element size
  real(kind=CUSTOM_REAL) :: vp,vs
  real(kind=CUSTOM_REAL) :: distance_min
  real(kind=CUSTOM_REAL) :: elemsize

  ! time steps
  double precision :: dtmin,dtmax
  double precision :: time_step
  double precision :: deltat_lts

  real(kind=CUSTOM_REAL),dimension(:),allocatable :: ispec_time_step

  ! p refinements ( local time step = dt/p )
  integer :: p
  integer :: p_min,p_max
  integer :: p_fine,p_coarse

  integer, dimension(:),allocatable :: p_level_relative

  integer :: nglob
  integer,dimension(:),allocatable :: iglob_p_refine

  integer,dimension(:),allocatable :: tmp_i
  integer :: i,j,ispec,iglob,ilevel
  integer :: ier

  logical :: found

  ! element flags
  logical,dimension(:,:),allocatable :: p_elem
  logical,dimension(:,:),allocatable :: boundary_elem

  integer,dimension(:),allocatable :: num_ispec_level_b

  double precision :: total_work,total_work_lts,total_work_lts_b
  double precision,dimension(:),allocatable :: xstore_dummy,ystore_dummy,zstore_dummy

  ! vtk output
  logical,parameter :: DEBUG_VTK_OUTPUT = .false.
  character(len=MAX_STRING_LEN) :: filename

  ! initializes
  num_p_level = 1

  ! user output
  if (LTS_MODE) then
    print *
    print *,'local time stepping: turned ON'
  else
    print *
    print *,'local time stepping: turned OFF'
    print *
    ! nothing to do
    return
  endif

  ! at this point, we don't have all GLL points yet, just anchor points (NGNOD per element)
  nglob = nnodes

  ! user output
  print *,'  setting up elements for local time stepping'
  print *,'  number of elements:',nspec
  print *,'  number of nodes   :',nglob

  ! define percentage of smallest distance between GLL points for NGLL points
  ! percentages were computed by calling the GLL points routine for each degree
  percent_GLL(1) = 0.d0
  percent_GLL(2) = 100.d0
  percent_GLL(3) = 50.d0
  percent_GLL(4) = 27.639320225002102d0
  percent_GLL(5) = 17.267316464601141d0
  percent_GLL(6) = 11.747233803526763d0
  percent_GLL(7) = 8.4888051860716516d0
  percent_GLL(8) = 6.4129925745196719d0
  percent_GLL(9) = 5.0121002294269914d0
  percent_GLL(10) = 4.0233045916770571d0
  percent_GLL(11) = 3.2999284795970416d0
  percent_GLL(12) = 2.7550363888558858d0
  percent_GLL(13) = 2.3345076678918053d0
  percent_GLL(14) = 2.0032477366369594d0
  percent_GLL(15) = 1.7377036748080721d0
  ! convert to real percentage
  percent_GLL(:) = percent_GLL(:) / 100.d0
  if (NGLLX < 2 .or. NGLLX > 15) stop 'error value of NGLLX not supported yet'

  ! allocates memory for p-refinement per element
  allocate(ispec_p_refine(nspec),stat=ier)
  if (ier /= 0) stop 'error allocating array ispec_p_refine'

  ! temporary arrays
  allocate(ispec_time_step(nspec),stat=ier)
  if (ier /= 0) stop 'error allocating array ispec_time_step'
  allocate(iglob_p_refine(nglob),stat=ier)
  if (ier /= 0) stop 'error allocating array iglob_p_refine'

  ! estimated global, stable time step based on CFL condition
  dtmin = 1.e30
  dtmax = 0.0
  ispec_time_step(:) = 0.0
  do ispec = 1, nspec
    ! determines velocities within this element
    call lts_get_elem_vpvs(ispec,vp,vs)

    ! determines size of element
    call lts_get_elem_size(ispec,elemsize)

    ! minimum distance based on GLL point distribution
    distance_min = elemsize * percent_GLL(NGLLX)

    ! estimated time step based on CFL condition
    time_step = COURANT_SUGGESTED * distance_min / max( vp,vs )

    ! stores value
    ispec_time_step(ispec) = time_step

    ! determines coarsest time step
    if (dtmax < time_step) then
      dtmax = time_step
    endif
    ! determines finest time step
    if (dtmin > time_step) then
      dtmin = time_step
    endif
  enddo
  print *,'  estimated time step min   = ',sngl(dtmin),' seconds'
  print *,'  estimated time step max   = ',sngl(dtmax),' seconds'
  print *,'  estimated time step ratio = ',sngl(dtmax/dtmin)
  print *

  ! users sets DT = .. in Par_file, we use this as a threshold limit for the minimum time step size
  ! (useful in case CFL is underestimating size)
  if (dtmin > DT) then
    print *,'  DT time step size set by Par_file: DT = ',sngl(DT),' limits coarsest time step size!'
    if (dtmin / DT >= 2.0) then
      print *,'  please set DT to higher value: ',sngl(dtmin),' otherwise LTS will not be optimal!!!'
      print *
    endif
  endif

  print *,'  suggested minimum DT time step = ',sngl(dtmin)
  print *

  ! we will use now the smallest time step estimate and put elements
  ! into bins with increasing time step size (power of 2 of smallest time step)

  ! setup p' refinements for each global point: local time step == dtmin * p'
  iglob_p_refine(:) = 1000000000
  do ispec = 1,nspec

    ! gets estimated time step based on CFL condition
    time_step = ispec_time_step(ispec)

    ! uses a 5% safety margin for binning
    time_step = ( 1.d0 - LTS_SAFETY_MARGIN ) * time_step

    ! local time step level compared to global time step size
    p = floor( time_step / dtmin )
    if (p < 1) p = 1

    ! debug
    !print *,'ispec',ispec,'time step=',time_step,dtmin,'     p refinement = ',p

    ! debug
    ! uniform distribution for testing
    if (LTS_SINGLE_P_LEVEL) then
      p = 1
    endif

    ! two-level distribution for testing
    if (LTS_TWO_P_LEVEL) then
      if (p > 1) p = 2
    endif

    ! enforces p' to be power of 2, find next smallest power of 2:
    ! example: 2 -> 2, 3 -> 2, 4 -> 4, 5 -> 4, .., 7 -> 4, 8 -> 8, 9 -> 8, ..
    if (p > 1) then
      i = 0
      do while ( 2**i <= p )
        i = i + 1
      enddo
      p = 2**(i-1)
    endif

    ! sets p' refinements on global points
    ! note: only corner points are defined, no GLL points yet (points can be shared between elements)
    do i = 1,NGNOD
      iglob = elmnts(i,ispec)
      ! note: here local time step == dtmin * p' ,thus sets p' only if smaller;
      !       we will switch to dt/p afterwards
      if (p < iglob_p_refine(iglob)) iglob_p_refine(iglob) = p
    enddo
  enddo
  p_min = minval(iglob_p_refine(:))
  p_max = maxval(iglob_p_refine(:))

  ! check
  if (p_min >= 1000000000) then
    print *,'error p value minimum: min/max = ',p_min,p_max
    stop 'error minimum p refinement must be 1, please check element distribution...'
  endif
  if (p_max >= 1000000000) then
    print *,'error p value maximum: min/max = ',p_min,p_max
    stop 'error maximum p refinement must be at least 1, please check element distribution...'
  endif

  ! coarsest time step for LTS (largest multiple p of smallest time step)
  deltat_lts = dtmin * dble(p_max)

  ! takes away percent of optimal time step for LTS stability
  ! LTS must decrease time step for stability depending on p depth
  ! to avoid this decrease, one can try to overlap the fine region by one/two elements
  if (LTS_DECREASE_DT) then
    ! sets new stable coarse time step for LTS
    deltat_lts = ( 1.d0 - LTS_STABILITY_MARGIN_DT ) * deltat_lts
    ! user output
    print *, '  decreases cfl time step by ',LTS_STABILITY_MARGIN_DT*100.0,'percent'
  endif

  ! user output
  print *, '  suggested global coarsest time step       = ',deltat_lts
  print *

  ! re-orders p' -> p = 1 / p'
  ! note: 1 == coarsest, p_max == smallest such that local time step == dt / p
  do iglob = 1,nglob
    p = iglob_p_refine(iglob)
    iglob_p_refine(iglob) = p_max / p
  enddo

  ! re-gets min/max p value on global nodes over all processes
  p_min = minval(iglob_p_refine(:))
  p_max = maxval(iglob_p_refine(:))

  ! from here on, a higher p-value indicates a smaller local time step dt/p

  ! sets p refinement for element
  ! (adds also elements touching global points with higher refinement)
  ispec_p_refine(:) = 0
  do ispec = 1,nspec
    do i = 1,NGNOD
      iglob = elmnts(i,ispec)
      p = iglob_p_refine(iglob)
      ! note: now local time step == dt/p
      if (p > ispec_p_refine(ispec)) ispec_p_refine(ispec) = p
    enddo
  enddo
  ! checks that all elements have a valid p-value assigned
  if (minval(ispec_p_refine(:)) == 0) stop 'error ispec_p_refine has zero entry'

  ! user output
  print *,'  p refinement of nodes   : min/max = ',p_min,p_max
  print *,'  p refinement of elements: min/max = ',minval(ispec_p_refine(:)),maxval(ispec_p_refine(:))
  print *

  ! adds coarser elements which overlap to finer p-level
  if (LTS_OVERLAP_REGION) then
    call lts_add_overlap_elements(nglob,nspec,iglob_p_refine,ispec_p_refine)

    p_min = minval(iglob_p_refine(:))
    p_max = maxval(iglob_p_refine(:))

    ! user output
    print *,'with overlap:'
    print *,'  p refinement of nodes   : min/max = ',p_min,p_max
    print *,'  p refinement of elements: min/max = ',minval(ispec_p_refine(:)),maxval(ispec_p_refine(:))
    print *
  endif

  ! counts number of needed different levels of refinement
  allocate(tmp_i(nspec),stat=ier)
  if (ier /= 0) stop 'error allocating array tmp_i'

  tmp_i(:) = 0
  num_p_level = 0
  do ispec = 1,nspec
    p = ispec_p_refine(ispec)
    ! search in table
    found = .false.
    do i = 1,num_p_level
      if (p == tmp_i(i)) then
        found = .true.
        exit
      endif
    enddo
    ! inserts new p level
    if (.not. found ) then
      ! increases levels
      num_p_level = num_p_level + 1
      ! adds new entry into sorted list
      tmp_i(num_p_level) = p
      ! sorts list in decreasing order (like in: 8 4 3 1)
      do i = 2, num_p_level
        j = i - 1
        p = tmp_i(i)
        do while (j > 1 .and. tmp_i(j) < p)
          tmp_i(j+1) = tmp_i(j)
          j = j - 1
        enddo
        if (j == 1 .and. tmp_i(j) < p) then
          tmp_i(j+1) = tmp_i(j)
          j = j - 1
        endif
        tmp_i(j+1) = p
      enddo
    endif
  enddo
  if (num_p_level < 1) stop 'error no p level found'

  ! number of local time step refinements for each LTS domain
  allocate(p_level(num_p_level),stat=ier)
  if (ier /= 0) stop 'error allocating array p_level'
  p_level(:) = tmp_i(1:num_p_level)

  ! free temporary array
  deallocate(tmp_i)

  ! checks array
  if (minval(p_level(:)) > 1) stop 'error coarsest p level starts with p > 1'
  if (maxval(p_level(:)) /= p_max) stop 'error smallest p level has wrong maximum value'
  if (.not. p_max == p_level(1) ) stop 'error first p_level does not start with maximum p value'

  ! sets up relative depth (in power of 2 ) of fine region to enclosing coarse one
  ! like from 16 8 2 1 -> 2 2 1 1 which works with the p_level_grid stepping
  allocate(p_level_relative(num_p_level),stat=ier)
  if (ier /= 0) stop 'error allocating array p_level_relative'

  p_level_relative(:) = 1
  do i = 1,num_p_level-1
    p_fine = p_level(i)
    p_coarse = p_level(i+1)
    ! relative depth of fine region to coarse one
    ! (how many steps to reach coarser one)
    if (mod(p_fine,2*p_coarse) == 0) then
      p = p_fine / (2*p_coarse)
    else
      do while (mod(p_fine,2*p_coarse) /= 0 )
        p_fine = p_fine + 1
        p = p_fine / ( 2*p_coarse)
      enddo
    endif
    p_level_relative(i) = p
    print *,'  level',i,':  fine/coarse p refinement = ',p_fine,p_coarse,'relative = ',p
  enddo
  ! finest level must increase time steps by factor 2 with current p_level_grid stepping scheme
  if (num_p_level > 1) then
    p_level_relative(1) = 2 * p_level_relative(1)
  endif

  ! user output
  print *,'  number of p-levels        : ',num_p_level
  print *,'  p-level array             : ',p_level(:)
  print *,'  p-level relative array    : ',p_level_relative(:)

  ! determines elements which belong to p-levels
  allocate(p_elem(nspec,num_p_level),stat=ier)
  if (ier /= 0) stop 'error allocating array p_elem'
  allocate(boundary_elem(nspec,num_p_level),stat=ier)
  if (ier /= 0) stop 'error allocating array boundary_elem'

  p_elem(:,:) = .false.
  boundary_elem(:,:) = .false.

  ! loops over refinement levels
  do ilevel = 1, num_p_level
    ! p refinement number in this specified p-level
    p = p_level(ilevel)

    ! sets element flags in this p-level
    do ispec = 1,nspec
      do i = 1,NGNOD
        iglob = elmnts(i,ispec)
        if (iglob_p_refine(iglob) == p) p_elem(ispec,ilevel) = .true.
      enddo
    enddo
    !if (itime == 1) then
    !  print *,'p elements:'
    !  print *,p_elem(:)
    !endif

    ! sets any boundary element flags with surrounding levels
    do ispec = 1,nspec
      ! elements with p-level global points
      ! note: also elements with only one node are p-elements
      if (p_elem(ispec,ilevel) .eqv. .true. ) then
        do i = 1,NGNOD
          iglob = elmnts(i,ispec)
          if (iglob_p_refine(iglob) /= p) boundary_elem(ispec,ilevel) = .true.
        enddo
      endif
    enddo
  enddo

  ! counts number of elements for each level
  allocate(num_ispec_level(num_p_level))
  num_ispec_level(:) = 0
  do ilevel = 1,num_p_level
    p = p_level(ilevel)
    num_ispec_level(ilevel) = count(ispec_p_refine(:) == p)
  enddo
  if (minval(num_ispec_level(:)) == 0) stop 'error p level without element found'

  ! user output
  print *,'  p-level number of elements: ',num_ispec_level(:)

  ! theoretical speed-up values (valid only for all elements belonging to same domain type)
  ! "pure" speed-up, without taking account of additional coarse/fine boundary contributions
  total_work_lts = 0.0
  do ilevel = 1,num_p_level
    p = p_level(ilevel)
    total_work_lts = total_work_lts + num_ispec_level(ilevel) * p
  enddo
  ! work without lts
  total_work = sum( num_ispec_level(:)) * maxval(p_level(:))

  ! counts boundary elements
  allocate(num_ispec_level_b(num_p_level))
  num_ispec_level_b(:) = 0
  do ilevel = 1,num_p_level
    do ispec = 1,nspec
      if (boundary_elem(ispec,ilevel)) then
        num_ispec_level_b(ilevel) = num_ispec_level_b(ilevel) + 1
      endif
    enddo
  enddo
  ! work with additional coarse/fine boundary contributions
  total_work_lts_b = total_work_lts
  do ilevel = 1,num_p_level
    p = p_level(ilevel)
    total_work_lts_b = total_work_lts_b + num_ispec_level_b(ilevel) * p
  enddo

  ! user output
  print *,'  p-level boundary elements : ',num_ispec_level_b(:)
  print *
  print *,'  theoretical speed-up value: ',sngl(total_work / total_work_lts),'(without boundary contributions)'
  print *,'  theoretical speed-up value: ',sngl(total_work / total_work_lts_b),'(with boundary contributions)'
  print *

  ! for vtk output
  if (DEBUG_VTK_OUTPUT) then
    ! p-refinements
    filename = trim(LOCAL_PATH)//'/lts_ispec_p_refine'
    allocate(xstore_dummy(nglob),ystore_dummy(nglob),zstore_dummy(nglob))
    xstore_dummy(:) = nodes_coords(1,:)
    ystore_dummy(:) = nodes_coords(2,:)
    zstore_dummy(:) = nodes_coords(3,:)
    call write_VTK_data_ngnod_elem_i(nspec,nglob,NGNOD,xstore_dummy,ystore_dummy,zstore_dummy, &
                                     elmnts,ispec_p_refine,filename)
    print *,'  written file: ',trim(filename)//'.vtk'
    print *
    deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
  endif

  ! free memory
  deallocate(p_level_relative)
  deallocate(p_elem,boundary_elem)
  deallocate(iglob_p_refine)
  deallocate(num_ispec_level_b)

  end subroutine lts_setup_elements

!-------------------------------

  subroutine lts_get_elem_vpvs(ispec,vp,vs)

  ! calculates the min/max velocity value of the specified element (ispec) for acoustic / elastic domains
  use constants, only: CUSTOM_REAL

  use decompose_mesh_par, only: LTS_MODE,mat,mat_prop

  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: vp,vs
  integer,intent(in) :: ispec

  ! local parameters
  integer :: num_mat,idomain_id

  ! get material_id for this element
  num_mat = mat(1,ispec)

  ! check
  if (LTS_MODE) then
    if (num_mat < 0) stop 'error element with undefined material property not supported yet in LTS'
  endif

  ! gets domain id ( 1 == acoustic / 2 == elastic / 3 == poroelastic
  idomain_id = mat_prop(7,num_mat)

  ! check
  if (LTS_MODE) then
    if (idomain_id < 1 .or. idomain_id > 2) stop 'error element domain id not supported yet in LTS'
  endif

  ! material is elastic or acoustic
  vp = mat_prop(2,num_mat)
  vs = mat_prop(3,num_mat)

  end subroutine lts_get_elem_vpvs

!-------------------------------

  subroutine lts_get_elem_size(ispec,elemsize)

  ! calculates the min/max size of the specified element (ispec)
  use constants, only: CUSTOM_REAL,HUGEVAL

  use decompose_mesh_par, only: NGNOD,elmnts,nodes_coords

  implicit none

  integer, intent(in) :: ispec
  real(kind=CUSTOM_REAL),intent(out) :: elemsize

  ! local parameters
  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max
  real(kind=CUSTOM_REAL) :: dx2,x0,y0,z0
  integer :: i,j,iglob_a,iglob_b

  ! initializes
  elemsize_min = HUGEVAL
  elemsize_max = -HUGEVAL

  ! loops over nodes of element
  do i = 1,NGNOD
    iglob_a = elmnts(i,ispec)

    x0 = nodes_coords(1,iglob_a)
    y0 = nodes_coords(2,iglob_a)
    z0 = nodes_coords(3,iglob_a)

    ! loops over all other corners
    do j = 1,NGNOD
      ! coordinates
      iglob_b = elmnts(j,ispec)

      ! skip reference node
      ! distances between points
      if (iglob_a /= iglob_b) then
        dx2 = (x0 - nodes_coords(1,iglob_b))**2 &
            + (y0 - nodes_coords(2,iglob_b))**2 &
            + (z0 - nodes_coords(3,iglob_b))**2
        if (dx2 < elemsize_min) elemsize_min = dx2
        if (dx2 > elemsize_max) elemsize_max = dx2
      endif
    enddo
  enddo

  ! returns smallest distance estimate
  elemsize = sqrt(elemsize_min)

  end subroutine lts_get_elem_size

