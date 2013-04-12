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

  subroutine compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
                               NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                               NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
                               NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
                               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB,&
                               USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  ! parameters read from parameter file
  integer NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA
  integer NER_REGULAR1,NER_REGULAR2,NER_REGULAR3
  integer NGLOB_DOUBLING,NGLOB_DOUBLING_1,NGLOB_DOUBLING_2,NGLOB_NO_DOUBLING
  integer NDOUBLINGS
  integer ner_doublings(2)

  ! parameters to be computed based upon parameters above read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER
  integer NEX_DOUBLING_ABOVE_PER_PROC_XI,NEX_DOUBLING_ABOVE_PER_PROC_ETA
  integer NEX_DOUBLING_ABOVE_XI,NEX_DOUBLING_ABOVE_ETA

  integer NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
       NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
       NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB

  integer NSPEC2D_DOUBLING_A_XI,NSPEC2D_DOUBLING_A_ETA
  integer NSPEC2D_DOUBLING_B_XI,NSPEC2D_DOUBLING_B_ETA
  integer NSPEC_DOUBLING_AB
  integer NUM_DOUBLING_BRICKS
  integer NUM2D_DOUBLING_BRICKS_XI,NUM2D_DOUBLING_BRICKS_ETA
  integer NSPEC1D_RADIAL_BEDROCK,NPOIN1D_RADIAL_BEDROCK

  integer NSPEC_NO_DOUBLING,NSPEC2D_NO_DOUBLING_XI,NSPEC2D_NO_DOUBLING_ETA

  logical USE_REGULAR_MESH

  ! number of elements horizontally in each slice (i.e. per processor)
  ! these two values MUST be equal in all cases
  NEX_PER_PROC_XI = NEX_XI / NPROC_XI
  NEX_PER_PROC_ETA = NEX_ETA / NPROC_ETA

  ! total number of processors in each of the six chunks
  NPROC = NPROC_XI * NPROC_ETA

  !
  !--- case of a regular mesh
  !
  if(USE_REGULAR_MESH) then

     ! exact number of spectral elements without doubling layers
     NSPEC_NO_DOUBLING = NEX_XI*NEX_ETA*NER

     ! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%

     ! exact number of surface elements for a chunk without doubling layers

     NSPEC2D_NO_DOUBLING_XI = NEX_PER_PROC_XI*NER

     NSPEC2D_NO_DOUBLING_ETA = NEX_PER_PROC_ETA*NER

     ! exact number of spectral elements
     !  NSPEC_AB = NSPEC_NO_DOUBLING / NPROC
     NSPEC_AB = NER * NEX_PER_PROC_XI * NEX_PER_PROC_ETA

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
     NPOIN1D_RADIAL_BEDROCK = NSPEC1D_RADIAL_BEDROCK*(NGLLZ_M-1)+1

     ! 2-D addressing and buffers for summation between slices
     ! we add one to number of points because of the flag after the last point
     NPOIN2DMAX_XMIN_XMAX = NSPEC2DMAX_XMIN_XMAX*NGLLY_M*NGLLZ_M + 1
     NPOIN2DMAX_YMIN_YMAX = NSPEC2DMAX_YMIN_YMAX*NGLLX_M*NGLLZ_M + 1

     ! exact number of global points
     NGLOB_AB = (NEX_PER_PROC_XI*(NGLLX_M-1)+1) * (NEX_PER_PROC_ETA*(NGLLY_M-1)+1) * (NER*(NGLLZ_M-1)+1)

     !
     !--- case of a non-regular mesh with mesh doublings
     !
  else if(NDOUBLINGS == 1) then

     ! number of spectral elements at the bottom of the doubling below the moho
     NEX_DOUBLING_ABOVE_XI=NEX_XI/2
     NEX_DOUBLING_ABOVE_ETA=NEX_ETA/2
     NEX_DOUBLING_ABOVE_PER_PROC_XI=NEX_PER_PROC_XI/2
     NEX_DOUBLING_ABOVE_PER_PROC_ETA=NEX_PER_PROC_ETA/2

     ! exact number of spectral elements without doubling layers
     NER_REGULAR2 = NER - ner_doublings(1)
     NER_REGULAR1 = ner_doublings(1) - 2

     NSPEC_NO_DOUBLING = &
          (NEX_XI/2)*(NEX_ETA/2)*NER_REGULAR1 &
          + NEX_XI*NEX_ETA*NER_REGULAR2

     ! exact number of spectral elements in the doubling regions

     ! number of elementary bricks in the two regions with doubling
     NUM_DOUBLING_BRICKS = (NEX_DOUBLING_ABOVE_XI*NEX_DOUBLING_ABOVE_ETA)/4

     ! for type AB, each doubling brick contains 32 elements on 2 levels
     NSPEC_DOUBLING_AB=32*NUM_DOUBLING_BRICKS

     ! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%

     ! exact number of surface elements for a chunk without doubling layers

     NSPEC2D_NO_DOUBLING_XI = &
          (NEX_PER_PROC_XI/2)*NER_REGULAR1 &
          + NEX_PER_PROC_XI*NER_REGULAR2

     NSPEC2D_NO_DOUBLING_ETA = &
          (NEX_PER_PROC_ETA/2)*NER_REGULAR1 &
          +NEX_PER_PROC_ETA*NER_REGULAR2

     ! exact number of surface elements in the doubling regions

     ! number of elementary bricks in the two regions with doubling
     NUM2D_DOUBLING_BRICKS_XI = (NEX_PER_PROC_XI/2 &
          + NEX_DOUBLING_ABOVE_PER_PROC_XI)/2

     NUM2D_DOUBLING_BRICKS_ETA = (NEX_PER_PROC_ETA/2 &
          + NEX_DOUBLING_ABOVE_PER_PROC_ETA)/2

     ! for type A, each doubling brick contains 10 elements on 3 levels
     NSPEC2D_DOUBLING_A_XI=10*NUM2D_DOUBLING_BRICKS_XI
     NSPEC2D_DOUBLING_A_ETA=10*NUM2D_DOUBLING_BRICKS_ETA

     ! for type B, each doubling brick contains 12 elements on 3 levels
     NSPEC2D_DOUBLING_B_XI=10*NUM2D_DOUBLING_BRICKS_XI
     NSPEC2D_DOUBLING_B_ETA=10*NUM2D_DOUBLING_BRICKS_ETA

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
     NSPEC2D_BOTTOM = (NEX_XI/2)*(NEX_ETA/2) / NPROC

     ! face with max number of elements is type B here
     ! maximum number of surface elements on vertical boundaries of the slices
     NSPEC2DMAX_XMIN_XMAX = NSPEC2D_B_ETA
     NSPEC2DMAX_YMIN_YMAX = NSPEC2D_B_XI

     ! theoretical number of Gauss-Lobatto points in radial direction
     !  NPOIN1D_RADIAL_BEDROCK = NSPEC1D_RADIAL_BEDROCK*(NGLLZ_M-1)+1

     ! 2-D addressing and buffers for summation between slices
     ! we add one to number of points because of the flag after the last point
     NPOIN2DMAX_XMIN_XMAX = NSPEC2DMAX_XMIN_XMAX*NGLLY_M*NGLLZ_M + 1
     NPOIN2DMAX_YMIN_YMAX = NSPEC2DMAX_YMIN_YMAX*NGLLX_M*NGLLZ_M + 1

     ! exact number of global points

     ! case of the doubling regions

     NGLOB_DOUBLING_1 = (NEX_PER_PROC_XI/4)*(NEX_PER_PROC_ETA/4)*(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
          - ((NEX_PER_PROC_XI/4-1)*(NEX_PER_PROC_ETA/4) &
          + (NEX_PER_PROC_ETA/4-1)*(NEX_PER_PROC_XI/4)) * (8*NGLLX_M**2-11*NGLLX_M+4) &
          + (NEX_PER_PROC_XI/4-1)*(NEX_PER_PROC_ETA/4-1)*(NGLLX_M+1)

     ! NGLOB_DOUBLING_1 = (NEX_PER_PROC_XI/2)*(NEX_PER_PROC_ETA/2)*(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
     !  - ((NEX_PER_PROC_XI/2-1)*(NEX_PER_PROC_ETA/2) &
     !  + (NEX_PER_PROC_ETA/2-1)*(NEX_PER_PROC_XI/2)) * (8*NGLLX_M**2-11*NGLLX_M+4)&
     !  + (NEX_PER_PROC_XI/2-1)*(NEX_PER_PROC_ETA/2-1)*NGLLX_M

     ! NGLOB_DOUBLING_2 =(NEX_PER_PROC_XI/8)*(NEX_PER_PROC_ETA/8)*(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
     !  - ((NEX_PER_PROC_XI/8-1)*(NEX_PER_PROC_ETA/8) &
     !  + (NEX_PER_PROC_ETA/8-1)*(NEX_PER_PROC_XI/8)) * (8*NGLLX_M**2-11*NGLLX_M+4)&
     !  + (NEX_PER_PROC_XI/8-1)*(NEX_PER_PROC_ETA/8-1)*NGLLX_M

     NGLOB_DOUBLING = NGLOB_DOUBLING_1 !+ NGLOB_DOUBLING_2

     ! compute number of points in blocks with no doubling
     ! exclude the four surfaces in contact with the doubling regions

     NGLOB_NO_DOUBLING = (NEX_PER_PROC_XI/2 + 1)*(NEX_PER_PROC_ETA/2 + 1)*NER_REGULAR1 &
          + (NEX_PER_PROC_XI + 1)*(NEX_PER_PROC_ETA + 1)*NER_REGULAR2

     NGLOB_AB = NGLOB_DOUBLING + NGLOB_NO_DOUBLING

  else if (NDOUBLINGS == 2) then

     ! number of spectral elements at the bottom of the doubling below the moho
     NEX_DOUBLING_ABOVE_XI=NEX_XI/2
     NEX_DOUBLING_ABOVE_ETA=NEX_ETA/2
     NEX_DOUBLING_ABOVE_PER_PROC_XI=NEX_PER_PROC_XI/2
     NEX_DOUBLING_ABOVE_PER_PROC_ETA=NEX_PER_PROC_ETA/2

     ! exact number of spectral elements without doubling layers
     NER_REGULAR3 = NER - ner_doublings(1)
     NER_REGULAR2 = (ner_doublings(1) - 2) - ner_doublings(2)
     NER_REGULAR1 = ner_doublings(2) - 2

     NSPEC_NO_DOUBLING = &
          (NEX_XI/4)*(NEX_ETA/4)*NER_REGULAR1 &
          + (NEX_XI/2)*(NEX_ETA/2)*NER_REGULAR2 &
          + NEX_XI*NEX_ETA*NER_REGULAR3

     ! exact number of spectral elements in the doubling regions

     ! number of elementary bricks in the two regions with doubling
     NUM_DOUBLING_BRICKS = (NEX_DOUBLING_ABOVE_XI*NEX_DOUBLING_ABOVE_ETA)/4&
          + (NEX_DOUBLING_ABOVE_XI*NEX_DOUBLING_ABOVE_ETA)/16

     ! for type AB, each doubling brick contains 32 elements on 2 levels
     NSPEC_DOUBLING_AB=32*NUM_DOUBLING_BRICKS


     ! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%

     ! exact number of surface elements for a chunk without doubling layers

     NSPEC2D_NO_DOUBLING_XI = &
          (NEX_PER_PROC_XI/4)*NER_REGULAR1 &
          + (NEX_PER_PROC_XI/2)*NER_REGULAR2 &
          + NEX_PER_PROC_XI*NER_REGULAR3

     NSPEC2D_NO_DOUBLING_ETA = &
          (NEX_PER_PROC_ETA/4)*NER_REGULAR1 &
          + (NEX_PER_PROC_ETA/2)*NER_REGULAR2 &
          + NEX_PER_PROC_ETA*NER_REGULAR3

     ! exact number of surface elements in the doubling regions

     ! number of elementary bricks in the two regions with doubling
     NUM2D_DOUBLING_BRICKS_XI = (NEX_PER_PROC_XI/2 &
          + NEX_DOUBLING_ABOVE_PER_PROC_XI)/2

     NUM2D_DOUBLING_BRICKS_ETA = (NEX_PER_PROC_ETA/2 &
          + NEX_DOUBLING_ABOVE_PER_PROC_ETA)/2

     ! for type A, each doubling brick contains 10 elements on 3 levels
     NSPEC2D_DOUBLING_A_XI=7*NUM2D_DOUBLING_BRICKS_XI
     NSPEC2D_DOUBLING_A_ETA=7*NUM2D_DOUBLING_BRICKS_ETA

     ! for type B, each doubling brick contains 12 elements on 3 levels
     NSPEC2D_DOUBLING_B_XI=10*NUM2D_DOUBLING_BRICKS_XI
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

     !  NSPEC1D_RADIAL_BEDROCK = (NER_BASEMENT_SEDIM+NER_16_BASEMENT+NER_MOHO_16)/2 + NER_BOTTOM_MOHO/4

     ! face with max number of elements is type B here
     ! maximum number of surface elements on vertical boundaries of the slices
     NSPEC2DMAX_XMIN_XMAX = NSPEC2D_B_ETA
     NSPEC2DMAX_YMIN_YMAX = NSPEC2D_B_XI

     !debug
     !print*,'nspec minmax:',NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_yMAX

     ! theoretical number of Gauss-Lobatto points in radial direction
     !  NPOIN1D_RADIAL_BEDROCK = NSPEC1D_RADIAL_BEDROCK*(NGLLZ_M-1)+1

     ! 2-D addressing and buffers for summation between slices
     ! we add one to number of points because of the flag after the last point
     NPOIN2DMAX_XMIN_XMAX = NSPEC2DMAX_XMIN_XMAX*NGLLY_M*NGLLZ_M + 1
     NPOIN2DMAX_YMIN_YMAX = NSPEC2DMAX_YMIN_YMAX*NGLLX_M*NGLLZ_M + 1


     ! exact number of global points

     ! case of the doubling regions

     NGLOB_DOUBLING_1 = (NEX_PER_PROC_XI/8)*(NEX_PER_PROC_ETA/8)*(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
          - ((NEX_PER_PROC_XI/8-1)*(NEX_PER_PROC_ETA/8) &
          + (NEX_PER_PROC_ETA/8-1)*(NEX_PER_PROC_XI/8)) * (8*NGLLX_M**2-11*NGLLX_M+4)&
          + (NEX_PER_PROC_XI/8-1)*(NEX_PER_PROC_ETA/8-1)*(NGLLX_M+1)

     NGLOB_DOUBLING_2 = (NEX_PER_PROC_XI/4)*(NEX_PER_PROC_ETA/4)*(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
          - ((NEX_PER_PROC_XI/4-1)*(NEX_PER_PROC_ETA/4) &
          + (NEX_PER_PROC_ETA/4-1)*(NEX_PER_PROC_XI/4)) * (8*NGLLX_M**2-11*NGLLX_M+4)&
          + (NEX_PER_PROC_XI/4-1)*(NEX_PER_PROC_ETA/4-1)*(NGLLX_M+1)

     NGLOB_DOUBLING = NGLOB_DOUBLING_1 + NGLOB_DOUBLING_2

     ! compute number of points in blocks with no doubling
     ! exclude the four surfaces in contact with the doubling regions

     NGLOB_NO_DOUBLING = (NEX_PER_PROC_XI/4 + 1)*(NEX_PER_PROC_ETA/4 + 1)*NER_REGULAR1 &
          + (NEX_PER_PROC_XI/2 + 1)*(NEX_PER_PROC_ETA/2 + 1)*(NER_REGULAR2 - 1) &
          + (NEX_PER_PROC_XI + 1)*(NEX_PER_PROC_ETA + 1)*NER_REGULAR3

     NGLOB_AB = NGLOB_DOUBLING + NGLOB_NO_DOUBLING


  endif ! end of section for non-regular mesh with doublings

end subroutine compute_parameters

