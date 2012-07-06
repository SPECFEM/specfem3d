!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
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

  subroutine compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM, &
      NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB,USE_REGULAR_MESH)

  implicit none

  include "constants.h"

! parameters read from parameter file
  integer NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA
  integer NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM

! parameters to be computed based upon parameters above read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER

  integer NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB

  integer NEX_DOUBLING_SEDIM_XI,NEX_DOUBLING_SEDIM_ETA
  integer NEX_DOUBLING_SEDIM_PER_PROC_XI,NEX_DOUBLING_SEDIM_PER_PROC_ETA
  integer NSPEC2D_DOUBLING_A_XI,NSPEC2D_DOUBLING_A_ETA
  integer NSPEC2D_DOUBLING_B_XI,NSPEC2D_DOUBLING_B_ETA
  integer NSPEC_DOUBLING_AB
  integer NUM_DOUBLING_BRICKS
  integer NUM2D_DOUBLING_BRICKS_XI,NUM2D_DOUBLING_BRICKS_ETA
  integer nglob_no_doubling_volume,nglob_no_doubling_surface
  integer nblocks_xi,nblocks_eta
  integer nglob_surface_typeA,nglob_surface_typeB
  integer NSPEC1D_RADIAL_BEDROCK,NPOIN1D_RADIAL_BEDROCK

  integer NSPEC_NO_DOUBLING,NSPEC2D_NO_DOUBLING_XI,NSPEC2D_NO_DOUBLING_ETA

  logical USE_REGULAR_MESH

!
!--- case of a regular mesh
!
  if(USE_REGULAR_MESH) then

! total number of spectral elements along radius
  NER = NER_SEDIM

! number of elements horizontally in each slice (i.e. per processor)
! these two values MUST be equal in all cases
  NEX_PER_PROC_XI = NEX_XI / NPROC_XI
  NEX_PER_PROC_ETA = NEX_ETA / NPROC_ETA

! total number of processors in each of the six chunks
  NPROC = NPROC_XI * NPROC_ETA

! exact number of spectral elements without doubling layers
  NSPEC_NO_DOUBLING = NEX_XI*NEX_ETA*NER_SEDIM

! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%

! exact number of surface elements for a chunk without doubling layers

  NSPEC2D_NO_DOUBLING_XI = NEX_PER_PROC_XI*NER_SEDIM

  NSPEC2D_NO_DOUBLING_ETA = NEX_PER_PROC_ETA*NER_SEDIM

! exact number of spectral elements
  NSPEC_AB = NSPEC_NO_DOUBLING / NPROC

! exact number of surface elements for faces A and B along XI and ETA
  NSPEC2D_A_XI = NSPEC2D_NO_DOUBLING_XI
  NSPEC2D_B_XI = NSPEC2D_NO_DOUBLING_XI
  NSPEC2D_A_ETA = NSPEC2D_NO_DOUBLING_ETA
  NSPEC2D_B_ETA = NSPEC2D_NO_DOUBLING_ETA

! exact number of surface elements on the bottom and top boundaries
! and theoretical number of spectral elements in radial direction

  NSPEC2D_TOP = NEX_XI*NEX_ETA / NPROC
  NSPEC2D_BOTTOM = NSPEC2D_TOP

  NSPEC1D_RADIAL_BEDROCK = NER

! face with max number of elements is type B here
! maximum number of surface elements on vertical boundaries of the slices
  NSPEC2DMAX_XMIN_XMAX = NSPEC2D_B_ETA
  NSPEC2DMAX_YMIN_YMAX = NSPEC2D_B_XI

! theoretical number of Gauss-Lobatto points in radial direction
  NPOIN1D_RADIAL_BEDROCK = NSPEC1D_RADIAL_BEDROCK*(NGLLZ-1)+1

! 2-D addressing and buffers for summation between slices
! we add one to number of points because of the flag after the last point
  NPOIN2DMAX_XMIN_XMAX = NSPEC2DMAX_XMIN_XMAX*NGLLY*NGLLZ + 1
  NPOIN2DMAX_YMIN_YMAX = NSPEC2DMAX_YMIN_YMAX*NGLLX*NGLLZ + 1

! exact number of global points
  NGLOB_AB = (NEX_PER_PROC_XI*(NGLLX-1)+1) * (NEX_PER_PROC_ETA*(NGLLY-1)+1) * (NER*(NGLLZ-1)+1)

!
!--- case of a non-regular mesh with mesh doublings
!
  else

! total number of spectral elements along radius
  NER = NER_BOTTOM_MOHO + NER_MOHO_16 + NER_16_BASEMENT + NER_BASEMENT_SEDIM + NER_SEDIM

! number of elements horizontally in each slice (i.e. per processor)
! these two values MUST be equal in all cases
  NEX_PER_PROC_XI = NEX_XI / NPROC_XI
  NEX_PER_PROC_ETA = NEX_ETA / NPROC_ETA

! total number of processors in each of the six chunks
  NPROC = NPROC_XI * NPROC_ETA

! number of spectral elements at the bottom of the doubling below the moho
  NEX_DOUBLING_SEDIM_XI=NEX_XI/2
  NEX_DOUBLING_SEDIM_ETA=NEX_ETA/2
  NEX_DOUBLING_SEDIM_PER_PROC_XI=NEX_PER_PROC_XI/2
  NEX_DOUBLING_SEDIM_PER_PROC_ETA=NEX_PER_PROC_ETA/2

! exact number of spectral elements without doubling layers
  NSPEC_NO_DOUBLING = &
     (NEX_DOUBLING_SEDIM_XI*NEX_DOUBLING_SEDIM_ETA*(NER_BASEMENT_SEDIM/2-3) &
    +(NEX_XI/4)*(NEX_ETA/4)*(NER_16_BASEMENT/2-3) &
    +(NEX_XI/4)*(NEX_ETA/4)*(NER_MOHO_16/2) &
    +(NEX_XI/4)*(NEX_ETA/4)*(NER_BOTTOM_MOHO/4)) + NEX_XI*NEX_ETA*NER_SEDIM

! exact number of spectral elements in the doubling regions

! number of elementary bricks in the two regions with doubling
  NUM_DOUBLING_BRICKS = ((NEX_XI/4)*(NEX_ETA/4) &
        +NEX_DOUBLING_SEDIM_XI*NEX_DOUBLING_SEDIM_ETA)/4

! for type AB, each doubling brick contains 40 elements on 3 levels
  NSPEC_DOUBLING_AB=40*NUM_DOUBLING_BRICKS

! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%

! exact number of surface elements for a chunk without doubling layers

  NSPEC2D_NO_DOUBLING_XI = &
      NEX_DOUBLING_SEDIM_PER_PROC_XI*(NER_BASEMENT_SEDIM/2-3) &
     +(NEX_PER_PROC_XI/4)*(NER_16_BASEMENT/2-3) &
     +(NEX_PER_PROC_XI/4)*(NER_MOHO_16/2) &
     +(NEX_PER_PROC_XI/4)*(NER_BOTTOM_MOHO/4) + NEX_PER_PROC_XI*NER_SEDIM

  NSPEC2D_NO_DOUBLING_ETA = &
       NEX_DOUBLING_SEDIM_PER_PROC_ETA*(NER_BASEMENT_SEDIM/2-3) &
     +(NEX_PER_PROC_ETA/4)*(NER_16_BASEMENT/2-3) &
     +(NEX_PER_PROC_ETA/4)*(NER_MOHO_16/2) &
     +(NEX_PER_PROC_ETA/4)*(NER_BOTTOM_MOHO/4) + NEX_PER_PROC_ETA*NER_SEDIM

! exact number of surface elements in the doubling regions

! number of elementary bricks in the two regions with doubling
  NUM2D_DOUBLING_BRICKS_XI = ((NEX_PER_PROC_XI/4) &
        +NEX_DOUBLING_SEDIM_PER_PROC_XI)/2

  NUM2D_DOUBLING_BRICKS_ETA = ((NEX_PER_PROC_ETA/4) &
        +NEX_DOUBLING_SEDIM_PER_PROC_ETA)/2

! for type A, each doubling brick contains 10 elements on 3 levels
  NSPEC2D_DOUBLING_A_XI=10*NUM2D_DOUBLING_BRICKS_XI
  NSPEC2D_DOUBLING_A_ETA=10*NUM2D_DOUBLING_BRICKS_ETA

! for type B, each doubling brick contains 12 elements on 3 levels
  NSPEC2D_DOUBLING_B_XI=12*NUM2D_DOUBLING_BRICKS_XI
  NSPEC2D_DOUBLING_B_ETA=12*NUM2D_DOUBLING_BRICKS_ETA

! exact number of spectral elements
  NSPEC_AB = (NSPEC_NO_DOUBLING + NSPEC_DOUBLING_AB) / NPROC

! exact number of surface elements for faces A and B
! along XI and ETA for doubling region
  NSPEC2D_A_XI = NSPEC2D_NO_DOUBLING_XI + NSPEC2D_DOUBLING_A_XI
  NSPEC2D_B_XI = NSPEC2D_NO_DOUBLING_XI + NSPEC2D_DOUBLING_B_XI
  NSPEC2D_A_ETA = NSPEC2D_NO_DOUBLING_ETA + NSPEC2D_DOUBLING_A_ETA
  NSPEC2D_B_ETA = NSPEC2D_NO_DOUBLING_ETA + NSPEC2D_DOUBLING_B_ETA

! exact number of surface elements on the bottom and top boundaries
! and theoretical number of spectral elements in radial direction

  NSPEC2D_TOP = NEX_XI*NEX_ETA / NPROC
  NSPEC2D_BOTTOM = (NEX_XI/4)*(NEX_ETA/4) / NPROC

  NSPEC1D_RADIAL_BEDROCK = (NER_BASEMENT_SEDIM+NER_16_BASEMENT+NER_MOHO_16)/2 + NER_BOTTOM_MOHO/4

! face with max number of elements is type B here
! maximum number of surface elements on vertical boundaries of the slices
  NSPEC2DMAX_XMIN_XMAX = NSPEC2D_B_ETA
  NSPEC2DMAX_YMIN_YMAX = NSPEC2D_B_XI

! theoretical number of Gauss-Lobatto points in radial direction
  NPOIN1D_RADIAL_BEDROCK = NSPEC1D_RADIAL_BEDROCK*(NGLLZ-1)+1

! 2-D addressing and buffers for summation between slices
! we add one to number of points because of the flag after the last point
  NPOIN2DMAX_XMIN_XMAX = NSPEC2DMAX_XMIN_XMAX*NGLLY*NGLLZ + 1
  NPOIN2DMAX_YMIN_YMAX = NSPEC2DMAX_YMIN_YMAX*NGLLX*NGLLZ + 1

! exact number of global points

! case of the doubling regions
! formulas computed using Mathematica for the basic 3D doubling brick

! compute number of points in blocks with no doubling
! exclude the three surfaces in contact with the doubling regions
  nglob_no_doubling_volume = (4*(NGLLX-1)+1)*(4*(NGLLX-1)+1)*((NER_BASEMENT_SEDIM/2-3 )*(NGLLX-1)-1) &
    +(2*(NGLLX-1)+1)*(2*(NGLLX-1)+1)*(((NER_16_BASEMENT/2+NER_MOHO_16/2+NER_BOTTOM_MOHO/4)-3)*(NGLLX-1)+0)

! number of basic blocks in each slice
  nblocks_xi = NEX_PER_PROC_XI / 8
  nblocks_eta = NEX_PER_PROC_ETA / 8

  NGLOB_AB = nblocks_xi*nblocks_eta*(200*NGLLX**3 - 484*NGLLX**2 + 392*NGLLX - 106 + nglob_no_doubling_volume)

! same thing for 2D surfaces for the three types of faces
  nglob_no_doubling_surface = (4*(NGLLX-1)+1)*((NER_BASEMENT_SEDIM/2-3)*(NGLLX-1)-1) &
    +(2*(NGLLX-1)+1)*(((NER_16_BASEMENT/2+NER_MOHO_16/2+NER_BOTTOM_MOHO/4)-3)*(NGLLX-1)+0)

  nglob_surface_typeA = 30*NGLLX**2 - 45 * NGLLX + 17
  nglob_surface_typeB = 36*NGLLX**2 - 57 * NGLLX + 23

! final number of points in volume obtained by removing planes counted twice
  NGLOB_AB = NGLOB_AB &
     - (nblocks_xi-1)*nblocks_eta*(nglob_surface_typeA + nglob_no_doubling_surface) &
     - (nblocks_eta-1)*nblocks_xi*(nglob_surface_typeB + nglob_no_doubling_surface) &
     + (nblocks_eta-1)*(nblocks_xi-1)*NPOIN1D_RADIAL_BEDROCK

! add number of points in the sediments
  NGLOB_AB = NGLOB_AB + (NEX_PER_PROC_XI*(NGLLX-1)+1) &
    *(NEX_PER_PROC_ETA*(NGLLY-1)+1)*(NER_SEDIM*(NGLLZ-1)+0)

  endif ! end of section for non-regular mesh with doublings

  end subroutine compute_parameters

