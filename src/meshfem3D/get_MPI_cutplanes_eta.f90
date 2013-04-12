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

  subroutine get_MPI_cutplanes_eta(myrank,prname,nspec,iMPIcut_eta,ibool, &
                        xstore,ystore,zstore,mask_ibool,npointot, &
                        NSPEC2D_A_XI,NSPEC2D_B_XI)

! this routine detects cut planes along eta
! In principle the left cut plane of the first slice
! and the right cut plane of the last slice are not used
! in the solver except if we want to have periodic conditions

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  integer nspec,myrank
  integer NSPEC2D_A_XI,NSPEC2D_B_XI

  logical iMPIcut_eta(2,nspec)

  integer ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)

  double precision xstore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision ystore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision zstore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)

! logical mask used to create arrays iboolleft_eta and iboolright_eta
  integer npointot
  logical mask_ibool(npointot)

! global element numbering
  integer ispec

! MPI cut-plane element numbering
  integer ispecc1,ispecc2,npoin2D_eta,ix,iy,iz
  integer nspec2Dtheor1,nspec2Dtheor2

! processor identification
  character(len=256) prname

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
  do ix=1,NGLLX_M
      do iz=1,NGLLZ_M

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
  iy = NGLLY_M
  do ix=1,NGLLX_M
      do iz=1,NGLLZ_M

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

