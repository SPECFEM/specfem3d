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

  subroutine get_MPI_cutplanes_eta(myrank,prname,nspec,iMPIcut_eta,ibool, &
                        xstore,ystore,zstore,mask_ibool,npointot, &
                        NSPEC2D_XI_FACE,iregion,NGLOB2DMAX_XY,mask_ibool2,npoin2D_eta)

! this routine detects cut planes along eta
! In principle the left cut plane of the first slice
! and the right cut plane of the last slice are not used
! in the solver except if we want to have periodic conditions

  implicit none

  include "constants.h"

  logical mask_ibool2(npointot)

  integer nspec,myrank,nglob,ipoin2D,NGLOB2DMAX_XY,iregion
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_XI_FACE


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
  integer nspec2Dtheor

! processor identification
  character(len=150) prname

! arrays for sorting routine
  integer, dimension(:), allocatable :: ind,ninseg,iglob,locval,iwork
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: ibool_selected
  double precision, dimension(:), allocatable :: xstore_selected,ystore_selected,zstore_selected


! allocate arrays for message buffers with maximum size
! define maximum size for message buffers
  if (USE_MESH_COLORING_INNER_OUTER) then
    allocate(ibool_selected(NGLOB2DMAX_XY))
    allocate(xstore_selected(NGLOB2DMAX_XY))
    allocate(ystore_selected(NGLOB2DMAX_XY))
    allocate(zstore_selected(NGLOB2DMAX_XY))
    allocate(ind(NGLOB2DMAX_XY))
    allocate(ninseg(NGLOB2DMAX_XY))
    allocate(iglob(NGLOB2DMAX_XY))
    allocate(locval(NGLOB2DMAX_XY))
    allocate(ifseg(NGLOB2DMAX_XY))
    allocate(iwork(NGLOB2DMAX_XY))
    allocate(work(NGLOB2DMAX_XY))
  endif

! theoretical number of surface elements in the buffers
! cut planes along eta=constant correspond to XI faces
      nspec2Dtheor = NSPEC2D_XI_FACE(iregion,1)

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
                if (USE_MESH_COLORING_INNER_OUTER) then
                  ibool_selected(npoin2D_eta) = ibool(ix,iy,iz,ispec)
                  xstore_selected(npoin2D_eta) = xstore(ix,iy,iz,ispec)
                  ystore_selected(npoin2D_eta) = ystore(ix,iy,iz,ispec)
                  zstore_selected(npoin2D_eta) = zstore(ix,iy,iz,ispec)
                else
                  write(10,*) ibool(ix,iy,iz,ispec)
                endif
            endif
          enddo
      enddo
    endif
  enddo

  if (USE_MESH_COLORING_INNER_OUTER) then
    call sort_array_coordinates(npoin2D_eta,xstore_selected,ystore_selected,zstore_selected, &
            ibool_selected,iglob,locval,ifseg,nglob,ind,ninseg,iwork,work)

    do ipoin2D=1,npoin2D_eta
        write(10,*) ibool_selected(ipoin2D)
        mask_ibool2(ibool_selected(ipoin2D)) = .true.
    enddo
  endif

! put flag to indicate end of the list of points
  write(10,*) 0

! write total number of points
  write(10,*) npoin2D_eta

  close(10)

! compare number of surface elements detected to analytical value
  if(ispecc1 /= nspec2Dtheor) call exit_MPI(myrank,'error MPI cut-planes detection in eta=left')

!
! determine if the element falls on the right MPI cut plane
!
      nspec2Dtheor = NSPEC2D_XI_FACE(iregion,2)

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
              if (USE_MESH_COLORING_INNER_OUTER) then
                ibool_selected(npoin2D_eta) = ibool(ix,iy,iz,ispec)
                xstore_selected(npoin2D_eta) = xstore(ix,iy,iz,ispec)
                ystore_selected(npoin2D_eta) = ystore(ix,iy,iz,ispec)
                zstore_selected(npoin2D_eta) = zstore(ix,iy,iz,ispec)
              else
                write(10,*) ibool(ix,iy,iz,ispec)
              endif
          endif
        enddo
      enddo
    endif
  enddo

  if (USE_MESH_COLORING_INNER_OUTER) then
    call sort_array_coordinates(npoin2D_eta,xstore_selected,ystore_selected,zstore_selected, &
            ibool_selected,iglob,locval,ifseg,nglob,ind,ninseg,iwork,work)

    do ipoin2D=1,npoin2D_eta
        write(10,*) ibool_selected(ipoin2D)
        mask_ibool2(ibool_selected(ipoin2D)) = .true.
    enddo
  endif

! put flag to indicate end of the list of points
  write(10,*) 0

! write total number of points
  write(10,*) npoin2D_eta

  close(10)

! compare number of surface elements detected to analytical value
  if(ispecc2 /= nspec2Dtheor) call exit_MPI(myrank,'error MPI cut-planes detection in eta=right')

  if (USE_MESH_COLORING_INNER_OUTER) then
    deallocate(ibool_selected)
    deallocate(xstore_selected)
    deallocate(ystore_selected)
    deallocate(zstore_selected)
    deallocate(ind)
    deallocate(ninseg)
    deallocate(iglob)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(iwork)
    deallocate(work)
  endif

  end subroutine get_MPI_cutplanes_eta

