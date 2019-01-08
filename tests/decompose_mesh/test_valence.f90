
!! DK DK put this because I do not know how to fix the rules.mk dependencies
  include "../../src/shared/serial.f90"

program test_valence

  use decompose_mesh_par

  implicit none
  integer :: i,ispec,inode,ier

  print *,'program: test_valence'

  ! ficticious setup 1
  print *,'test 1:'
  print *,'-------'
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
      inode = inode + 1
      elmnts(i,ispec) = inode
    enddo
  enddo

  nnodes = inode
  print *,'nnodes = ',nnodes

  ! initializes
  nsize = -1

  ! valence routine
  call check_valence()

  ! checks
  print *,'nsize = ',nsize
  if (nsize /= 1) then
    print *,'error valence: nsize =',nsize,'should be 1'
    stop 1
  else
    print *,'  result is correct'
  endif

  deallocate(elmnts)

  ! ficticious setup 2
  print *,'test 2:'
  print *,'-------'
  nspec = 12
  NGNOD = 27

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
      ! connects arbitrary nodes
      if (ispec == 3 .and. i == 3) then
        elmnts(i,ispec) = 2
      else if (ispec == 4 .and. i == 9) then
        elmnts(i,ispec) = 2
      else if (ispec == 5 .and. i == 12) then
        elmnts(i,ispec) = 2
      else if (ispec == 9 .and. i == 17) then
        elmnts(i,ispec) = 2
      else if (ispec == 12 .and. i == 27) then
        elmnts(i,ispec) = 2
      else
        inode = inode + 1
        elmnts(i,ispec) = inode
      endif
    enddo
  enddo

  nnodes = inode
  print *,'nnodes = ',nnodes

  ! initializes
  nsize = -1

  ! valence routine
  call check_valence()

  ! checks
  print *,'nsize = ',nsize
  if (nsize /= 6) then
    print *,'error valence: nsize =',nsize,'should be 6'
    stop 1
  else
    print *,'  result is correct'
  endif

  deallocate(elmnts)

  ! done
  print *,'test_valence done successfully'

end program test_valence

