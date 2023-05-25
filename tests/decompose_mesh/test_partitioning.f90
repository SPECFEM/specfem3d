
!! DK DK put this because I do not know how to fix the rules.mk dependencies
  include "../../src/shared/serial.f90"

program test_partitioning

  use constants, only: NDIM
  use decompose_mesh_par
  use fault_scotch, only: ANY_FAULT

  implicit none
  integer :: i,ispec,inode,ier
  integer :: itype
  integer :: ix,imod
  double precision :: width,xloc

  print *,'program: test_partitioning'

  ! ficticious setup
  nspec = 12
  NGNOD = 8

  ! allocates mesh elements
  allocate(elmnts(NGNOD,nspec),stat=ier)
  if (ier /= 0) stop 'error allocating array elmnts'
  elmnts(:,:) = 0

  print *,'nspec  = ',nspec
  print *,'NGNOD  = ',NGNOD

  ! element ids
  inode = 0
  do ispec = 1,nspec
    do i = 1,NGNOD
      if (ispec == 1) then
        ! adds new point
        inode = inode + 1
        elmnts(i,ispec) = inode
      else
        ! connects left side with previous element right side
        select case(i)
        case(1)
          elmnts(i,ispec) = elmnts(5,ispec-1)
        case (2)
          elmnts(i,ispec) = elmnts(6,ispec-1)
        case (3)
          elmnts(i,ispec) = elmnts(7,ispec-1)
        case (4)
          elmnts(i,ispec) = elmnts(8,ispec-1)
        case default
          ! adds new point
          inode = inode + 1
          elmnts(i,ispec) = inode
        end select
      endif
    enddo
  enddo

  nnodes = inode
  print *,'nnodes = ',nnodes

  ! ficticious coordinates
  allocate(nodes_coords(NDIM,nnodes),stat=ier)
  if (ier /= 0) stop 'error allocating array nodes_coords'
  nodes_coords(:,:) = 0.0
  width = 2000.0
  xloc = 0.0
  do inode = 1,nnodes
    imod = mod(inode,4)
    ix = int((inode-1)/4)
    ! x-coord
    if (imod == 1) then
      if (ix > 0) then
        if (ix <= 4) then
          xloc = xloc + width ! adds default element size
        else
          xloc = xloc + width / 4.0 ! smaller elements
        endif
      endif
    endif
    nodes_coords(1,inode) = xloc

    ! y-coord
    select case (imod)
    case (0,1)
      nodes_coords(2,inode) = 0.0
    case(2,3)
      nodes_coords(2,inode) = width
    end select
    ! z-coord
    select case (imod)
    case (1,2)
      nodes_coords(3,inode) = 0.0
    case(0,3)
      nodes_coords(3,inode) = width
    end select
  enddo
  print *,'nodes coords:'
  do ispec = 1,nspec
    do i = 1,NGNOD
      inode = elmnts(i,ispec)
      print *,'  xyz',ispec,i,nodes_coords(1,inode),nodes_coords(2,inode),nodes_coords(3,inode)
    enddo
  enddo
  print *

  ! initializes
  nsize = 0
  sup_neighbor = 0

  ! valence routine
  call check_valence()

  ! checks
  print *,'nsize = ',nsize
  if (nsize /= 2) then
    print *,'error valence: nsize =',nsize,'should be 2'
    stop 1
  else
    print *,'  result is correct'
  endif

  ! assigns elastic material
  allocate(mat(2,nspec),stat=ier)
  if (ier /= 0) stop 'error allocating array mat'
  mat(:,:) = 0
  do ispec = 1, nspec
    mat(1,ispec) = 1
  enddo

  count_def_mat = 1
  count_undef_mat = 0
  allocate(mat_prop(17,count_def_mat),stat=ier)
  if (ier /= 0) stop 'error allocating array mat_prop'
  allocate(undef_mat_prop(6,count_undef_mat),stat=ier)
  if (ier /= 0) stop 'error allocating array undef_mat_prop'
  mat_prop(:,:) = 0.d0
  undef_mat_prop(:,:) = ''

  num_mat = 1
  mat_prop(1,num_mat) = 2500.d0   ! rho
  mat_prop(2,num_mat) = 3200.d0   ! vp
  mat_prop(3,num_mat) = 1100.d0   ! vs
  mat_prop(4,num_mat) = 9000.d0   ! qkappa
  mat_prop(5,num_mat) = 200.d0    ! qmu
  mat_prop(6,num_mat) = 0         ! aniso_flag
  mat_prop(7,num_mat) = 1         ! idomain_id (1 == acoustic)

  ANY_FAULT = .false.
  nparts = 4
  num_p_level = 1

  ! without LTS
  LTS_MODE = .false.

  print *,'LTS mode: ',LTS_MODE
  print *

  ! partitioning routine
  print *
  print *,'calling partitioning routine...'
  print *
  do itype = 1,2
    ! partitioner
    select case (itype)
    case (1)
      PARTITIONING_TYPE = 1 ! SCOTCH
    case (2)
      PARTITIONING_TYPE = 4 ! ROWS_PART
    end select

    print *,'partitioning_type: ',PARTITIONING_TYPE
    print *

    call decompose_mesh()

    ! check
    print *,'partitioning results:'
    do i = 0,nparts-1
      print *,'partition ',i
      do ispec = 1,nspec
        if (part(ispec) == i) print *,'  contains element ',ispec
        ! checks if element belongs to no partition
        if (part(ispec) < 0) then
          print *,'error partitioning: element ',ispec,'has invalid partition number: ',part(ispec)
          stop 1
        endif
      enddo
    enddo
    print *

    ! re-set values
    elmnts(:,:) = elmnts(:,:) + 1
    ! frees array
    deallocate(part)

  enddo

  ! with LTS
  LTS_MODE = .true.
  DT = 0.1

  print *,'LTS mode: ',LTS_MODE
  print *
  print *,'calling lts_setup_elements() routine...'
  print *
  call lts_setup_elements()

  ! checks
  print *,'num_p_level = ',num_p_level
  if (num_p_level /= 2) then
    print *,'error partitioning: num_p_level =',num_p_level,'should be 2'
    stop 1
  else
    print *,'  result is correct'
  endif

  ! partitioning routine
  print *
  print *,'calling partitioning routine...'
  print *

  PARTITIONING_TYPE = 1 ! SCOTCH
  print *,'partitioning_type: ',PARTITIONING_TYPE
  print *

  call decompose_mesh()

  ! check
  print *,'partitioning results:'
  do i = 0,nparts-1
    print *,'partition ',i
    do ispec = 1,nspec
      if (part(ispec) == i) print *,'  contains element ',ispec
      ! checks if element belongs to no partition
      if (part(ispec) < 0) then
        print *,'error partitioning: element ',ispec,'has invalid partition number: ',part(ispec)
        stop 1
      endif
    enddo
  enddo
  print *

  deallocate(part)
  deallocate(nodes_coords)
  deallocate(elmnts)
  deallocate(mat,mat_prop,undef_mat_prop)

  ! done
  print *,'test_partitioning done successfully'

end program test_partitioning

