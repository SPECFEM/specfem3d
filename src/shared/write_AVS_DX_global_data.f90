!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

! create AVS or DX 3D data for the slice, to be recombined in postprocessing
  subroutine write_AVS_DX_global_data(myrank,prname,nspec,ibool,idoubling, &
                 xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot)

  implicit none

  include "constants.h"

  integer nspec,myrank
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! logical mask used to output global points only once
  integer npointot
  logical mask_ibool(npointot)

! numbering of global AVS or DX points
  integer num_ibool_AVS_DX(npointot)

  integer ispec
  integer iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer npoin,numpoin

! processor identification
  character(len=256) prname

! ------------------------------------

! writing points
  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpoints.txt',status='unknown')

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.


! mark global AVS or DX points
  do ispec=1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)
    mask_ibool(iglob1) = .true.
    mask_ibool(iglob2) = .true.
    mask_ibool(iglob3) = .true.
    mask_ibool(iglob4) = .true.
    mask_ibool(iglob5) = .true.
    mask_ibool(iglob6) = .true.
    mask_ibool(iglob7) = .true.
    mask_ibool(iglob8) = .true.
  enddo

! count global number of AVS or DX points
  npoin = count(mask_ibool(:))

! number of points in AVS or DX file
  write(10,*) npoin

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

! output global AVS or DX points
  numpoin = 0
  do ispec=1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)
    if(.not. mask_ibool(iglob1)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob1) = numpoin
      write(10,*) numpoin,sngl(xstore(1,1,1,ispec)), &
              sngl(ystore(1,1,1,ispec)),sngl(zstore(1,1,1,ispec))
    endif
    if(.not. mask_ibool(iglob2)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob2) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,1,1,ispec)), &
              sngl(ystore(NGLLX,1,1,ispec)),sngl(zstore(NGLLX,1,1,ispec))
    endif
    if(.not. mask_ibool(iglob3)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob3) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,NGLLY,1,ispec)), &
              sngl(ystore(NGLLX,NGLLY,1,ispec)),sngl(zstore(NGLLX,NGLLY,1,ispec))
    endif
    if(.not. mask_ibool(iglob4)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob4) = numpoin
      write(10,*) numpoin,sngl(xstore(1,NGLLY,1,ispec)), &
              sngl(ystore(1,NGLLY,1,ispec)),sngl(zstore(1,NGLLY,1,ispec))
    endif
    if(.not. mask_ibool(iglob5)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob5) = numpoin
      write(10,*) numpoin,sngl(xstore(1,1,NGLLZ,ispec)), &
              sngl(ystore(1,1,NGLLZ,ispec)),sngl(zstore(1,1,NGLLZ,ispec))
    endif
    if(.not. mask_ibool(iglob6)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob6) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,1,NGLLZ,ispec)), &
              sngl(ystore(NGLLX,1,NGLLZ,ispec)),sngl(zstore(NGLLX,1,NGLLZ,ispec))
    endif
    if(.not. mask_ibool(iglob7)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob7) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec)), &
              sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec)),sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
    endif
    if(.not. mask_ibool(iglob8)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob8) = numpoin
      write(10,*) numpoin,sngl(xstore(1,NGLLY,NGLLZ,ispec)), &
              sngl(ystore(1,NGLLY,NGLLZ,ispec)),sngl(zstore(1,NGLLY,NGLLZ,ispec))
    endif
    mask_ibool(iglob1) = .true.
    mask_ibool(iglob2) = .true.
    mask_ibool(iglob3) = .true.
    mask_ibool(iglob4) = .true.
    mask_ibool(iglob5) = .true.
    mask_ibool(iglob6) = .true.
    mask_ibool(iglob7) = .true.
    mask_ibool(iglob8) = .true.
  enddo

! check that number of global points output is okay
  if(numpoin /= npoin) &
    call exit_MPI(myrank,'incorrect number of global points in AVS or DX file creation')

  close(10)

! writing elements
  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelements.txt',status='unknown')

! number of elements in AVS or DX file
  write(10,*) nspec

! output global AVS or DX elements
  do ispec=1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)
    write(10,*) ispec,idoubling(ispec),num_ibool_AVS_DX(iglob1), &
                  num_ibool_AVS_DX(iglob2),num_ibool_AVS_DX(iglob3), &
                  num_ibool_AVS_DX(iglob4),num_ibool_AVS_DX(iglob5), &
                  num_ibool_AVS_DX(iglob6),num_ibool_AVS_DX(iglob7), &
                  num_ibool_AVS_DX(iglob8)
  enddo

  close(10)

  end subroutine write_AVS_DX_global_data

