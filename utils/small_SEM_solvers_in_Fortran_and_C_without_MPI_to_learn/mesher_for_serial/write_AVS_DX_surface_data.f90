!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

! create AVS or DX 2D data for the surface of the model
! to be recombined in postprocessing
  subroutine write_AVS_DX_surface_data(myrank,prname,nspec,iboun, &
     ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot)

  implicit none

  include "constants.h"

  integer nspec,myrank
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

  logical iboun(6,nspec)

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! logical mask used to output global points only once
  integer npointot
  logical mask_ibool(npointot)

! numbering of global AVS or DX points
  integer num_ibool_AVS_DX(npointot)

  integer ispec
  integer, dimension(8) :: iglobval
  integer npoin,numpoin,nspecface,ispecface

! processor identification
  character(len=150) prname

! writing points
  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpointssurface.txt',status='unknown')

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  nspecface = 0

! mark global AVS or DX points
  do ispec=1,nspec
! only if at the surface (top plane)
  if(iboun(6,ispec)) then

    iglobval(5)=ibool(1,1,NGLLZ,ispec)
    iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
    iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

! element is at the surface
    nspecface = nspecface + 1
    mask_ibool(iglobval(5)) = .true.
    mask_ibool(iglobval(6)) = .true.
    mask_ibool(iglobval(7)) = .true.
    mask_ibool(iglobval(8)) = .true.

  endif
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
! only if at the surface
  if(iboun(6,ispec)) then

    iglobval(5)=ibool(1,1,NGLLZ,ispec)
    iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
    iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

! top face
  if(iboun(6,ispec)) then

    if(.not. mask_ibool(iglobval(5))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(5)) = numpoin
      write(10,*) numpoin,sngl(xstore(1,1,NGLLZ,ispec)), &
              sngl(ystore(1,1,NGLLZ,ispec)),sngl(zstore(1,1,NGLLZ,ispec))
    endif

    if(.not. mask_ibool(iglobval(6))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(6)) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,1,NGLLZ,ispec)), &
              sngl(ystore(NGLLX,1,NGLLZ,ispec)),sngl(zstore(NGLLX,1,NGLLZ,ispec))
    endif

    if(.not. mask_ibool(iglobval(7))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(7)) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec)), &
              sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec)),sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
    endif

    if(.not. mask_ibool(iglobval(8))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(8)) = numpoin
      write(10,*) numpoin,sngl(xstore(1,NGLLY,NGLLZ,ispec)), &
              sngl(ystore(1,NGLLY,NGLLZ,ispec)),sngl(zstore(1,NGLLY,NGLLZ,ispec))
    endif

    mask_ibool(iglobval(5)) = .true.
    mask_ibool(iglobval(6)) = .true.
    mask_ibool(iglobval(7)) = .true.
    mask_ibool(iglobval(8)) = .true.

  endif

  endif
  enddo

! check that number of global points output is okay
  if(numpoin /= npoin) &
    call exit_MPI(myrank,'incorrect number of global points in AVS or DX file creation')

  close(10)

! output global AVS or DX elements

! writing elements
  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelementssurface.txt',status='unknown')

! number of elements in AVS or DX file
  write(10,*) nspecface

  ispecface = 0
  do ispec=1,nspec
! only if at the surface
  if(iboun(6,ispec)) then

    iglobval(5)=ibool(1,1,NGLLZ,ispec)
    iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
    iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

! top face
    ispecface = ispecface + 1
    write(10,*) ispecface,idoubling(ispec),num_ibool_AVS_DX(iglobval(5)), &
                  num_ibool_AVS_DX(iglobval(6)),num_ibool_AVS_DX(iglobval(7)), &
                  num_ibool_AVS_DX(iglobval(8))

  endif
  enddo

! check that number of surface elements output is okay
  if(ispecface /= nspecface) &
    call exit_MPI(myrank,'incorrect number of surface elements in AVS or DX file creation')

  close(10)

  end subroutine write_AVS_DX_surface_data

