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


  subroutine setup_color_perm(myrank,nspec,nglob,ibool,ANISOTROPY,SAVE_MESH_FILES)

! sets up mesh coloring and permutes elements

  use create_regions_mesh_ext_par
  implicit none

  integer :: myrank,nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  logical :: ANISOTROPY,SAVE_MESH_FILES

  ! local parameters
  integer, dimension(:), allocatable :: perm
  integer :: ier

  ! user output
  if(myrank == 0) then
    write(IMAIN,*) '     use coloring = ',USE_MESH_COLORING_GPU
  endif

  ! initializes
  num_colors_outer_acoustic = 0
  num_colors_inner_acoustic = 0
  num_colors_outer_elastic = 0
  num_colors_inner_elastic = 0

  ! mesh coloring
  if( USE_MESH_COLORING_GPU ) then

    ! creates coloring of elements
    allocate(perm(nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating temporary perm array'
    perm(:) = 0

    ! acoustic domains
    if( ACOUSTIC_SIMULATION ) then
      if( myrank == 0) then
        write(IMAIN,*) '     acoustic domains: '
      endif
      call setup_color(myrank,nspec,nglob,ibool,perm, &
                          ispec_is_acoustic,1, &
                          num_phase_ispec_acoustic,phase_ispec_inner_acoustic, &
                          SAVE_MESH_FILES)
    else
      ! allocates dummy arrays
      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
      if( ier /= 0 ) stop 'error allocating num_elem_colors_acoustic array'
    endif

    ! elastic domains
    if( ELASTIC_SIMULATION ) then
      if( myrank == 0) then
        write(IMAIN,*) '     elastic domains: '
      endif
      call setup_color(myrank,nspec,nglob,ibool,perm, &
                          ispec_is_elastic,2, &
                          num_phase_ispec_elastic,phase_ispec_inner_elastic, &
                          SAVE_MESH_FILES)
    else
      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if( ier /= 0 ) stop 'error allocating num_elem_colors_elastic array'
    endif

    ! checks: after all domains are done
    if(minval(perm) /= 1) &
      call exit_MPI(myrank, 'minval(perm) should be 1')
    if(maxval(perm) /= max(num_phase_ispec_acoustic,num_phase_ispec_elastic)) &
      call exit_MPI(myrank, 'maxval(perm) should be max(num_phase_ispec_..)')

    ! sorts array according to permutation
    call sync_all()
    if(myrank == 0) then
      write(IMAIN,*) '     mesh permutation:'
    endif
    call setup_permutation(myrank,nspec,nglob,ibool,ANISOTROPY,perm, &
                              SAVE_MESH_FILES)

    deallocate(perm)

  else

    ! allocates dummy arrays
    allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
    if( ier /= 0 ) stop 'error allocating num_elem_colors_acoustic array'
    allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
    if( ier /= 0 ) stop 'error allocating num_elem_colors_elastic array'

  endif ! USE_MESH_COLORING_GPU

  end subroutine setup_color_perm

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_color(myrank,nspec,nglob,ibool,perm, &
                            ispec_is_d,idomain, &
                            num_phase_ispec_d,phase_ispec_inner_d, &
                            SAVE_MESH_FILES)

! sets up mesh coloring

  use create_regions_mesh_ext_par
  implicit none

  integer :: myrank,nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  integer, dimension(nspec) :: perm

  ! wrapper array for ispec is in domain:
  ! idomain acoustic == 1 or elastic == 2
  integer :: idomain
  logical, dimension(nspec) :: ispec_is_d
  integer :: num_phase_ispec_d
  integer, dimension(num_phase_ispec_d,2) :: phase_ispec_inner_d

  logical :: SAVE_MESH_FILES

  ! local parameters
  ! added for color permutation
  integer :: nb_colors_outer_elements,nb_colors_inner_elements
  integer, dimension(:), allocatable :: num_of_elems_in_this_color
  integer, dimension(:), allocatable :: color
  integer, dimension(:), allocatable :: first_elem_number_in_this_color
  logical, dimension(:), allocatable :: is_on_a_slice_edge

  integer :: nspec_outer,nspec_inner,nspec_domain
  integer :: nspec_outer_min_global,nspec_outer_max_global
  integer :: nb_colors,nb_colors_min,nb_colors_max

  integer :: icolor,ispec,ispec_counter
  integer :: ispec_inner,ispec_outer
  integer :: ier

  character(len=2),dimension(2) :: str_domain = (/ "ac", "el" /)
  character(len=256) :: filename

  logical, parameter :: DEBUG = .false.

  !!!! David Michea: detection of the edges, coloring and permutation separately

  ! implement mesh coloring for GPUs if needed, to create subsets of disconnected elements
  ! to remove dependencies and the need for atomic operations in the sum of
  ! elemental contributions in the solver

  ! allocates temporary array with colors
  allocate(color(nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating temporary color array'
  allocate(first_elem_number_in_this_color(MAX_NUMBER_OF_COLORS + 1),stat=ier)
  if( ier /= 0 ) stop 'error allocating first_elem_number_in_this_color array'

  ! flags for elements on outer rims
  ! opposite to what is stored in ispec_is_inner
  allocate(is_on_a_slice_edge(nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating is_on_a_slice_edge array'
  do ispec = 1,nspec
    is_on_a_slice_edge(ispec) = .not. ispec_is_inner(ispec)
  enddo

  ! fast element coloring scheme
  call get_perm_color_faster(is_on_a_slice_edge,ispec_is_d, &
                            ibool,perm,color, &
                            nspec,nglob, &
                            nb_colors_outer_elements,nb_colors_inner_elements, &
                            nspec_outer,nspec_inner,nspec_domain, &
                            first_elem_number_in_this_color, &
                            myrank)

  ! for the last color, the next color is fictitious and its first (fictitious) element number is nspec + 1
  first_elem_number_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements + 1) &
    = nspec_domain + 1

  allocate(num_of_elems_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements),stat=ier)
  if( ier /= 0 ) stop 'error allocating num_of_elems_in_this_color array'

  num_of_elems_in_this_color(:) = 0
  do icolor = 1, nb_colors_outer_elements + nb_colors_inner_elements
    num_of_elems_in_this_color(icolor) = first_elem_number_in_this_color(icolor+1) - first_elem_number_in_this_color(icolor)
  enddo

  ! check that the sum of all the numbers of elements found in each color is equal
  ! to the total number of elements in the mesh
  if(sum(num_of_elems_in_this_color) /= nspec_domain) then
    print *,'error number of elements in this color:',idomain
    print *,'rank: ',myrank,' nspec = ',nspec_domain
    print *,'  total number of elements in all the colors of the mesh = ', &
      sum(num_of_elems_in_this_color)
    call exit_MPI(myrank, 'incorrect total number of elements in all the colors of the mesh')
  endif

  ! check that the sum of all the numbers of elements found in each color for the outer elements is equal
  ! to the total number of outer elements found in the mesh
  if(sum(num_of_elems_in_this_color(1:nb_colors_outer_elements)) /= nspec_outer) then
    print *,'error number of outer elements in this color:',idomain
    print *,'rank: ',myrank,' nspec_outer = ',nspec_outer
    print*,'nb_colors_outer_elements = ',nb_colors_outer_elements
    print *,'total number of elements in all the colors of the mesh for outer elements = ', &
      sum(num_of_elems_in_this_color(1:nb_colors_outer_elements))
    call exit_MPI(myrank, 'incorrect total number of elements in all the colors of the mesh for outer elements')
  endif

  ! debug: file output
  if( SAVE_MESH_FILES .and. DEBUG ) then
    filename = prname(1:len_trim(prname))//'color_'//str_domain(idomain)
    call write_VTK_data_elem_i(nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                              color,filename)
  endif

  deallocate(first_elem_number_in_this_color)
  deallocate(is_on_a_slice_edge)
  deallocate(color)

  ! debug: no mesh coloring, only creates dummy coloring arrays
  if( DEBUG ) then
    nb_colors_outer_elements = 0
    nb_colors_inner_elements = 0
    ispec_counter = 0

    ! first generate all the outer elements
    do ispec = 1,nspec
      if( ispec_is_d(ispec) ) then
        if( ispec_is_inner(ispec) .eqv. .false. ) then
          ispec_counter = ispec_counter + 1
          perm(ispec) = ispec_counter
        endif
      endif
    enddo

    ! store total number of outer elements
    nspec_outer = ispec_counter

    ! only single color
    if(nspec_outer > 0 ) nb_colors_outer_elements = 1

    ! then generate all the inner elements
    do ispec = 1,nspec
      if( ispec_is_d(ispec) ) then
        if( ispec_is_inner(ispec) .eqv. .true. ) then
          ispec_counter = ispec_counter + 1
          perm(ispec) = ispec_counter - nspec_outer ! starts again at 1
        endif
      endif
    enddo
    nspec_inner = ispec_counter - nspec_outer

    ! only single color
    if(nspec_inner > 0 ) nb_colors_inner_elements = 1

    allocate(num_of_elems_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements),stat=ier)
    if( ier /= 0 ) stop 'error allocating num_of_elems_in_this_color array'

    if( nspec_outer > 0 ) num_of_elems_in_this_color(1) = nspec_outer
    if( nspec_inner > 0 ) num_of_elems_in_this_color(2) = nspec_inner
  endif ! debug

  ! debug: saves mesh coloring numbers into files
  if( DEBUG ) then
    ! debug file output
    filename = prname(1:len_trim(prname))//'num_of_elems_in_this_color_'//str_domain(idomain)//'.dat'
    open(unit=99,file=trim(filename),status='unknown',iostat=ier)
    if( ier /= 0 ) stop 'error opening num_of_elems_in_this_color file'
    ! number of colors for outer elements
    write(99,*) nb_colors_outer_elements
    ! number of colors for inner elements
    write(99,*) nb_colors_inner_elements
    ! number of elements in each color
    ! outer elements
    do icolor = 1, nb_colors_outer_elements + nb_colors_inner_elements
      write(99,*) num_of_elems_in_this_color(icolor)
    enddo
    close(99)
  endif

  ! sets up domain coloring arrays
  select case(idomain)
  case( 1 )
    ! acoustic domains
    num_colors_outer_acoustic = nb_colors_outer_elements
    num_colors_inner_acoustic = nb_colors_inner_elements

    allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
    if( ier /= 0 ) stop 'error allocating num_elem_colors_acoustic array'

    num_elem_colors_acoustic(:) = num_of_elems_in_this_color(:)

  case( 2 )
    ! elastic domains
    num_colors_outer_elastic = nb_colors_outer_elements
    num_colors_inner_elastic = nb_colors_inner_elements

    allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
    if( ier /= 0 ) stop 'error allocating num_elem_colors_elastic array'

    num_elem_colors_elastic(:) = num_of_elems_in_this_color(:)

  case default
    stop 'error idomain not recognized'
  end select

  ! sets up elements for loops in simulations
  ispec_inner = 0
  ispec_outer = 0
  do ispec = 1, nspec
    ! only elements in this domain
    if( ispec_is_d(ispec) ) then

      ! sets phase_ispec arrays with ordering of elements
      if( ispec_is_inner(ispec) .eqv. .false. ) then
        ! outer elements
        ispec_outer = perm(ispec)

        ! checks
        if( ispec_outer < 1 .or. ispec_outer > num_phase_ispec_d ) then
          print*,'error outer permutation:',idomain
          print*,'rank:',myrank,'  ispec_inner = ',ispec_outer
          print*,'num_phase_ispec_d = ',num_phase_ispec_d
          call exit_MPI(myrank,'error outer acoustic permutation')
        endif

        phase_ispec_inner_d(ispec_outer,1) = ispec

      else
        ! inner elements
        ispec_inner = perm(ispec)

        ! checks
        if( ispec_inner < 1 .or. ispec_inner > num_phase_ispec_d ) then
          print*,'error inner permutation:',idomain
          print*,'rank:',myrank,'  ispec_inner = ',ispec_inner
          print*,'num_phase_ispec_d = ',num_phase_ispec_d
          call exit_MPI(myrank,'error inner acoustic permutation')
        endif

        phase_ispec_inner_d(ispec_inner,2) = ispec

      endif
    endif
  enddo

  ! total number of colors
  nb_colors = nb_colors_inner_elements + nb_colors_outer_elements
  call min_all_i(nb_colors,nb_colors_min)
  call max_all_i(nb_colors,nb_colors_max)

  ! user output
  call min_all_i(nspec_outer,nspec_outer_min_global)
  call max_all_i(nspec_outer,nspec_outer_max_global)
  call min_all_i(nspec_outer,nspec_outer_min_global)
  call max_all_i(nspec_outer,nspec_outer_max_global)
  if(myrank == 0) then
    write(IMAIN,*) '       colors min = ',nb_colors_min
    write(IMAIN,*) '       colors max = ',nb_colors_max
    write(IMAIN,*) '       outer elements: min = ',nspec_outer_min_global
    write(IMAIN,*) '       outer elements: max = ',nspec_outer_max_global
  endif

  ! debug: outputs permutation array as vtk file
  if( DEBUG ) then
    filename = prname(1:len_trim(prname))//'perm_'//str_domain(idomain)
    call write_VTK_data_elem_i(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        perm,filename)
  endif

  deallocate(num_of_elems_in_this_color)

  end subroutine setup_color

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_permutation(myrank,nspec,nglob,ibool,ANISOTROPY,perm, &
                                  SAVE_MESH_FILES)

  use create_regions_mesh_ext_par
  implicit none

  integer :: myrank,nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  logical :: ANISOTROPY

  integer, dimension(nspec),intent(inout) :: perm

  logical :: SAVE_MESH_FILES

  ! local parameters
  ! added for sorting
  integer, dimension(:,:,:,:), allocatable :: temp_array_int
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: temp_array_real
  logical, dimension(:), allocatable :: temp_array_logical_1D

  integer, dimension(:), allocatable :: temp_perm_global
  logical, dimension(:), allocatable :: mask_global

  integer :: icolor,icounter,ispec,ielem,ier,i
  integer :: iface,old_ispec,new_ispec

  character(len=256) :: filename

  logical,parameter :: DEBUG = .false.

  ! sorts array according to permutation
  allocate(temp_perm_global(nspec),stat=ier)
  if( ier /= 0 ) stop 'error temp_perm_global array'

  ! global ordering
  temp_perm_global(:) = 0
  icounter = 0

  ! fills global permutation array
  ! starts with elastic elements
  if( ELASTIC_SIMULATION ) then
    ! first outer elements coloring
    ! phase element counter
    ielem = 0
    do icolor = 1,num_colors_outer_elastic
      ! loops through elements
      do i = 1,num_elem_colors_elastic(icolor)
        ielem = ielem + 1
        ispec = phase_ispec_inner_elastic(ielem,1) ! 1 <-- first phase, outer elements
        ! reorders elements
        icounter = icounter + 1
        temp_perm_global(ispec) = icounter
        ! resets to new order
        phase_ispec_inner_elastic(ielem,1) = icounter
      enddo
    enddo
    ! inner elements coloring
    ielem = 0
    do icolor = num_colors_outer_elastic+1,num_colors_outer_elastic+num_colors_inner_elastic
      ! loops through elements
      do i = 1,num_elem_colors_elastic(icolor)
        ielem = ielem + 1
        ispec = phase_ispec_inner_elastic(ielem,2) ! 2 <-- second phase, inner elements
        ! reorders elements
        icounter = icounter + 1
        temp_perm_global(ispec) = icounter
        ! resets to new order
        phase_ispec_inner_elastic(ielem,2) = icounter
      enddo
    enddo
  endif

  ! continues with acoustic elements
  if( ACOUSTIC_SIMULATION ) then
    ! first outer elements coloring
    ! phase element counter
    ielem = 0
    do icolor = 1,num_colors_outer_acoustic
      ! loops through elements
      do i = 1,num_elem_colors_acoustic(icolor)
        ielem = ielem + 1
        ispec = phase_ispec_inner_acoustic(ielem,1) ! 1 <-- first phase, outer elements
        ! reorders elements
        icounter = icounter + 1
        temp_perm_global(ispec) = icounter
        ! resets to new order
        phase_ispec_inner_acoustic(ielem,1) = icounter
      enddo
    enddo
    ! inner elements coloring
    ielem = 0
    do icolor = num_colors_outer_acoustic+1,num_colors_outer_acoustic+num_colors_inner_acoustic
      ! loops through elements
      do i = 1,num_elem_colors_acoustic(icolor)
        ielem = ielem + 1
        ispec = phase_ispec_inner_acoustic(ielem,2) ! 2 <-- second phase, inner elements
        ! reorders elements
        icounter = icounter + 1
        temp_perm_global(ispec) = icounter
        ! resets to new order
        phase_ispec_inner_acoustic(ielem,2) = icounter
      enddo
    enddo
  endif

  ! checks
  if( icounter /= nspec ) then
    print*,'error temp perm: ',icounter,nspec
    stop 'error temporary global permutation incomplete'
  endif

  ! checks perm entries
  if(minval(temp_perm_global) /= 1) call exit_MPI(myrank, 'minval(temp_perm_global) should be 1')
  if(maxval(temp_perm_global) /= nspec) call exit_MPI(myrank, 'maxval(temp_perm_global) should be nspec')

  ! checks if every element was uniquely set
  allocate(mask_global(nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating temporary mask_global'
  mask_global(:) = .false.
  icounter = 0 ! counts permutations
  do ispec = 1, nspec
    new_ispec = temp_perm_global(ispec)
    ! checks bounds
    if( new_ispec < 1 .or. new_ispec > nspec ) call exit_MPI(myrank,'error temp_perm_global ispec bounds')
    ! checks if already set
    if( mask_global(new_ispec) ) then
      print*,'error temp_perm_global:',ispec,new_ispec,'element already set'
      call exit_MPI(myrank,'error global permutation')
    else
      mask_global(new_ispec) = .true.
    endif
    ! counts permutations
    if( new_ispec /= ispec ) icounter = icounter + 1
  enddo

  ! checks number of set elements
  if( count(mask_global(:)) /= nspec ) then
    print*,'error temp_perm_global:',count(mask_global(:)),nspec,'permutation incomplete'
    call exit_MPI(myrank,'error global permutation incomplete')
  endif
  deallocate(mask_global)

  ! user output
  if(myrank == 0) then
    write(IMAIN,*) '       number of permutations = ',icounter
  endif

  ! outputs permutation array as vtk file
  if( SAVE_MESH_FILES .and. DEBUG ) then
    filename = prname(1:len_trim(prname))//'perm_global'
    call write_VTK_data_elem_i(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        temp_perm_global,filename)
  endif

  ! store as new permutation
  perm(:) = temp_perm_global(:)
  deallocate(temp_perm_global)

  ! permutes all required mesh arrays according to new ordering

  ! permutation of ibool
  allocate(temp_array_int(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating temporary temp_array_int'
  call permute_elements_integer(ibool,temp_array_int,perm,nspec)
  deallocate(temp_array_int)

  ! element domain flags
  allocate(temp_array_logical_1D(nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating temporary temp_array_logical_1D'
  call permute_elements_logical1D(ispec_is_acoustic,temp_array_logical_1D,perm,nspec)
  call permute_elements_logical1D(ispec_is_elastic,temp_array_logical_1D,perm,nspec)
  call permute_elements_logical1D(ispec_is_poroelastic,temp_array_logical_1D,perm,nspec)
  call permute_elements_logical1D(ispec_is_inner,temp_array_logical_1D,perm,nspec)
  deallocate(temp_array_logical_1D)

  ! mesh arrays
  allocate(temp_array_real(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating temporary temp_array_real'
  call permute_elements_real(xixstore,temp_array_real,perm,nspec)
  call permute_elements_real(xiystore,temp_array_real,perm,nspec)
  call permute_elements_real(xizstore,temp_array_real,perm,nspec)
  call permute_elements_real(etaxstore,temp_array_real,perm,nspec)
  call permute_elements_real(etaystore,temp_array_real,perm,nspec)
  call permute_elements_real(etazstore,temp_array_real,perm,nspec)
  call permute_elements_real(gammaxstore,temp_array_real,perm,nspec)
  call permute_elements_real(gammaystore,temp_array_real,perm,nspec)
  call permute_elements_real(gammazstore,temp_array_real,perm,nspec)
  call permute_elements_real(jacobianstore,temp_array_real,perm,nspec)

  ! material parameters
  call permute_elements_real(kappastore,temp_array_real,perm,nspec)
  call permute_elements_real(mustore,temp_array_real,perm,nspec)

  ! acoustic arrays
  if( ACOUSTIC_SIMULATION ) then
    call permute_elements_real(rhostore,temp_array_real,perm,nspec)
  endif

  ! elastic arrays
  if( ELASTIC_SIMULATION ) then
    call permute_elements_real(rho_vp,temp_array_real,perm,nspec)
    call permute_elements_real(rho_vs,temp_array_real,perm,nspec)
    if( ANISOTROPY ) then
      call permute_elements_real(c11store,temp_array_real,perm,nspec)
      call permute_elements_real(c12store,temp_array_real,perm,nspec)
      call permute_elements_real(c13store,temp_array_real,perm,nspec)
      call permute_elements_real(c14store,temp_array_real,perm,nspec)
      call permute_elements_real(c15store,temp_array_real,perm,nspec)
      call permute_elements_real(c16store,temp_array_real,perm,nspec)
      call permute_elements_real(c22store,temp_array_real,perm,nspec)
      call permute_elements_real(c23store,temp_array_real,perm,nspec)
      call permute_elements_real(c24store,temp_array_real,perm,nspec)
      call permute_elements_real(c25store,temp_array_real,perm,nspec)
      call permute_elements_real(c33store,temp_array_real,perm,nspec)
      call permute_elements_real(c34store,temp_array_real,perm,nspec)
      call permute_elements_real(c35store,temp_array_real,perm,nspec)
      call permute_elements_real(c36store,temp_array_real,perm,nspec)
      call permute_elements_real(c44store,temp_array_real,perm,nspec)
      call permute_elements_real(c45store,temp_array_real,perm,nspec)
      call permute_elements_real(c46store,temp_array_real,perm,nspec)
      call permute_elements_real(c55store,temp_array_real,perm,nspec)
      call permute_elements_real(c56store,temp_array_real,perm,nspec)
      call permute_elements_real(c66store,temp_array_real,perm,nspec)
    endif
  endif
  deallocate(temp_array_real)

  ! poroelastic arrays
  if( POROELASTIC_SIMULATION ) then
    stop 'mesh permutation for poroelastic simulations not supported yet'
  endif

  ! boundary surface
  if( num_abs_boundary_faces > 0 ) then
    do iface = 1,num_abs_boundary_faces
      old_ispec = abs_boundary_ispec(iface)
      new_ispec = perm(old_ispec)
      abs_boundary_ispec(iface) = new_ispec
    enddo
  endif

  ! free surface
  if( num_free_surface_faces > 0 ) then
    do iface = 1,num_free_surface_faces
      old_ispec = free_surface_ispec(iface)
      new_ispec = perm(old_ispec)
      free_surface_ispec(iface) = new_ispec
    enddo
  endif

  ! coupling surface
  if( num_coupling_ac_el_faces > 0 ) then
    do iface = 1,num_coupling_ac_el_faces
      old_ispec = coupling_ac_el_ispec(iface)
      new_ispec = perm(old_ispec)
      coupling_ac_el_ispec(iface) = new_ispec
    enddo
  endif
  if( num_coupling_ac_po_faces > 0 ) then
    do iface = 1,num_coupling_ac_po_faces
      old_ispec = coupling_ac_po_ispec(iface)
      new_ispec = perm(old_ispec)
      coupling_ac_po_ispec(iface) = new_ispec
    enddo
  endif
  if( num_coupling_el_po_faces > 0 ) then
    do iface = 1,num_coupling_el_po_faces
      ! elastic-poroelastic
      old_ispec = coupling_el_po_ispec(iface)
      new_ispec = perm(old_ispec)
      coupling_el_po_ispec(iface) = new_ispec
      ! poroelastic-elastic
      old_ispec = coupling_po_el_ispec(iface)
      new_ispec = perm(old_ispec)
      coupling_po_el_ispec(iface) = new_ispec
    enddo
  endif

  ! moho surface
  if( NSPEC2D_MOHO > 0 ) then
    allocate(temp_array_logical_1D(nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating temporary temp_array_logical_1D'
    call permute_elements_logical1D(is_moho_top,temp_array_logical_1D,perm,nspec)
    call permute_elements_logical1D(is_moho_bot,temp_array_logical_1D,perm,nspec)
    deallocate(temp_array_logical_1D)
    do iface = 1,NSPEC2D_MOHO
      ! top
      old_ispec = ibelm_moho_top(iface)
      new_ispec = perm(old_ispec)
      ibelm_moho_top(iface) = new_ispec
      ! bottom
      old_ispec = ibelm_moho_bot(iface)
      new_ispec = perm(old_ispec)
      ibelm_moho_bot(iface) = new_ispec
    enddo
  endif

  end subroutine setup_permutation

