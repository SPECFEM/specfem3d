!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 1
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology October 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_MPI_cutplanes_eta(myrank,prname,nspec,iMPIcut_eta,ibool, &
                        xstore,ystore,zstore,mask_ibool,npointot, &
                        NSPEC2D_A_XI,NSPEC2D_B_XI)

! this routine detects cut planes along eta
! In principle the left cut plane of the first slice
! and the right cut plane of the last slice are not used
! in the solver except if we want to have periodic conditions

  implicit none

  include "constants.h"

  integer nspec,myrank
  integer NSPEC2D_A_XI,NSPEC2D_B_XI

  logical iMPIcut_eta(2,nspec)

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! logical mask used to create arrays iboolleft_eta and iboolright_eta
  integer npointot
  logical mask_ibool(npointot)

! global element numbering
  integer ispec

! MPI cut-plane element numbering
  integer ispecc1,ispecc2,npoin2D_eta,ix,iy,iz
  integer nspec2Dtheor1,nspec2Dtheor2

! processor identification
  character(len=150) prname

! theoretical number of surface elements in the buffers
! cut planes along eta=constant correspond to XI faces
      nspec2Dtheor1 = NSPEC2D_A_XI
      nspec2Dtheor2 = NSPEC2D_B_XI

! write the MPI buffers for the left and right edges of the slice
! and the position of the points to check that the buffers are fine

!
! determine if the element falls on the left MPI cut plane
!

! global point number and coordinates left MPI cut-plane
  open(unit=10,file=prname(1:len_trim(prname))//'iboolleft_eta.txt',status='unknown')

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

! nb of global points shared with the other slice
  npoin2D_eta = 0

! nb of elements in this cut-plane
  ispecc1=0

  do ispec=1,nspec
  if(iMPIcut_eta(1,ispec)) then

    ispecc1=ispecc1+1

! loop on all the points in that 2-D element, including edges
  iy = 1
  do ix=1,NGLLX
      do iz=1,NGLLZ

! select point, if not already selected
  if(.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
      mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
      npoin2D_eta = npoin2D_eta + 1

      write(10,*) ibool(ix,iy,iz,ispec),xstore(ix,iy,iz,ispec), &
              ystore(ix,iy,iz,ispec),zstore(ix,iy,iz,ispec)
  endif

      enddo
  enddo

  endif
  enddo

! put flag to indicate end of the list of points
  write(10,*) '0 0  0.  0.  0.'

! write total number of points
  write(10,*) npoin2D_eta

  close(10)

! compare number of surface elements detected to analytical value
  if(ispecc1 /= nspec2Dtheor1 .and. ispecc1 /= nspec2Dtheor2) &
    call exit_MPI(myrank,'error MPI cut-planes detection in eta=left')

!
! determine if the element falls on the right MPI cut plane
!

! global point number and coordinates right MPI cut-plane
  open(unit=10,file=prname(1:len_trim(prname))//'iboolright_eta.txt',status='unknown')

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

! nb of global points shared with the other slice
  npoin2D_eta = 0

! nb of elements in this cut-plane
  ispecc2=0

  do ispec=1,nspec
  if(iMPIcut_eta(2,ispec)) then

    ispecc2=ispecc2+1

! loop on all the points in that 2-D element, including edges
  iy = NGLLY
  do ix=1,NGLLX
      do iz=1,NGLLZ

! select point, if not already selected
  if(.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
      mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
      npoin2D_eta = npoin2D_eta + 1

      write(10,*) ibool(ix,iy,iz,ispec),xstore(ix,iy,iz,ispec), &
              ystore(ix,iy,iz,ispec),zstore(ix,iy,iz,ispec)
  endif

      enddo
  enddo

  endif
  enddo

! put flag to indicate end of the list of points
  write(10,*) '0 0  0.  0.  0.'

! write total number of points
  write(10,*) npoin2D_eta

  close(10)

! compare number of surface elements detected to analytical value
  if(ispecc2 /= nspec2Dtheor1 .and. ispecc2 /= nspec2Dtheor2) &
    call exit_MPI(myrank,'error MPI cut-planes detection in eta=right')

  end subroutine get_MPI_cutplanes_eta

