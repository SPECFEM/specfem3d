!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! output some statistics about the mesh

  subroutine save_mesh_statistics(nspec,nglob,NEX_XI,NEX_ETA,NPROC,UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX)

  implicit none

  include "constants.h"

  integer nspec,nglob
  integer NEX_XI,NEX_ETA,NPROC

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX

  write(IMAIN,*)
  write(IMAIN,*) 'mesh statistics:'
  write(IMAIN,*) '---------------'
  write(IMAIN,*)
  write(IMAIN,*) 'number of processors = ',NPROC
  write(IMAIN,*)
  write(IMAIN,*) 'total number of elements per AB slice = ',nspec
  write(IMAIN,*) 'total number of points per AB slice = ',nglob
  write(IMAIN,*)
  write(IMAIN,*) 'total for full mesh:'
  write(IMAIN,*) '-------------------'
  write(IMAIN,*)
  write(IMAIN,*) 'exact total number of spectral elements in entire mesh = ',NPROC*nspec
  write(IMAIN,*) 'approximate total number of points in entire mesh = ',dble(NPROC)*dble(nglob)
! there are 3 DOFs in solid regions
  write(IMAIN,*) 'approximate total number of degrees of freedom in entire mesh = ',3.d0*dble(NPROC)*dble(nglob)
  write(IMAIN,*)
  write(IMAIN,*) 'resolution of the mesh at the surface:'
  write(IMAIN,*) '-------------------------------------'
  write(IMAIN,*)
  write(IMAIN,*) 'spectral elements along X = ',NEX_XI
  write(IMAIN,*) 'spectral elements along Y = ',NEX_ETA
  write(IMAIN,*) 'GLL points along X = ',NEX_XI*(NGLLX-1) + 1
  write(IMAIN,*) 'GLL points along Y = ',NEX_ETA*(NGLLY-1) + 1
  write(IMAIN,*) 'average distance between points along X in m = ',sngl(UTM_X_MAX-UTM_X_MIN)/real(NEX_XI*(NGLLX-1))
  write(IMAIN,*) 'average distance between points along Y in m = ',sngl(UTM_Y_MAX-UTM_Y_MIN)/real(NEX_ETA*(NGLLY-1))
  write(IMAIN,*)

  end subroutine save_mesh_statistics

