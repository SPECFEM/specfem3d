
!! DK DK put this because I do not know how to fix the rules.mk dependencies
  include "../../src/shared/serial.f90"

program test_partitioning

  use decompose_mesh

  implicit none
  integer :: i

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
      ! connects last element with first one
      if (ispec == nspec .and. i == 5) then
        elmnts(i,ispec) = elmnts(1,1)
      else if (ispec == nspec .and. i == 6) then
        elmnts(i,ispec) = elmnts(2,1)
      else if (ispec == nspec .and. i == 7) then
        elmnts(i,ispec) = elmnts(3,1)
      else if (ispec == nspec .and. i == 8) then
        elmnts(i,ispec) = elmnts(4,1)
      else
        inode = inode + 1
        elmnts(i,ispec) = inode
      endif
    enddo
  enddo

  nnodes = inode
  print *,'nnodes = ',nnodes

  ! initializes
  nsize = 0
  sup_neighbour = 0

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
  allocate(mat_prop(16,count_def_mat),stat=ier)
  if (ier /= 0) stop 'error allocating array mat_prop'
  allocate(undef_mat_prop(6,count_undef_mat),stat=ier)
  if (ier /= 0) stop 'error allocating array undef_mat_prop'
  mat_prop(:,:) = 0.d0
  undef_mat_prop(:,:) = ''

  num_mat = 1
  mat_prop(1,num_mat) = 2500.d0   ! rho
  mat_prop(2,num_mat) = 3200.d0   ! vp
  mat_prop(3,num_mat) = 1100.d0   ! vs
  mat_prop(4,num_mat) = 200.d0    ! qmu
  mat_prop(5,num_mat) = 0         ! aniso_flag
  mat_prop(6,num_mat) = 1         ! idomain_id
  mat_prop(7,num_mat) = 9000.d0   ! qkappa

  ANY_FAULT = .false.
  nparts = 4

  ! partitioning routine
  print *,'calling partitioning routine...'

  call scotch_partitioning()

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

  ! this crashes, not sure why...
  !deallocate(elmnts)
  !deallocate(mat,mat_prop,undef_mat_prop)

  ! done
  print *,'test_partitioning done successfully'

end program test_partitioning

