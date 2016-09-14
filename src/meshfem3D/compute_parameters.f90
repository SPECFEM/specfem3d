!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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
                                NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB, &
                                USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

  use constants

  implicit none

  include "constants_meshfem3D.h"

  ! parameters read from parameter file
  integer,intent(in) :: NER
  integer,intent(in) :: NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA
  integer,intent(in) :: NDOUBLINGS

  logical,intent(in):: USE_REGULAR_MESH
  integer,dimension(NDOUBLINGS),intent(in) :: ner_doublings

  ! parameters to be computed based upon parameters above read from file
  integer,intent(out) :: NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer,intent(out) :: NSPEC_AB,NGLOB_AB

  integer,intent(out) :: NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_A_ETA,NSPEC2D_B_ETA

  integer,intent(out) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP
  integer,intent(out) :: NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX

  ! local parameters
  integer :: NGLOB_DOUBLING,NGLOB_DOUBLING_SINGLE
  integer :: NER_REGULAR_SINGLE
  integer :: NGLOB_NO_DOUBLING
  integer :: NSPEC_NO_DOUBLING,NSPEC2D_NO_DOUBLING_XI,NSPEC2D_NO_DOUBLING_ETA

  !integer :: NGLOB_DOUBLING_1,NGLOB_DOUBLING_2
  !integer :: NEX_DOUBLING_ABOVE_PER_PROC_XI,NEX_DOUBLING_ABOVE_PER_PROC_ETA
  !integer :: NEX_DOUBLING_ABOVE_XI,NEX_DOUBLING_ABOVE_ETA
  !integer :: NER_REGULAR1,NER_REGULAR2,NER_REGULAR3

  integer :: NUM_DOUBLING_BRICKS
  integer :: NUM2D_DOUBLING_BRICKS_XI,NUM2D_DOUBLING_BRICKS_ETA

  integer :: NSPEC2D_DOUBLING_A_XI,NSPEC2D_DOUBLING_A_ETA
  integer :: NSPEC2D_DOUBLING_B_XI,NSPEC2D_DOUBLING_B_ETA
  integer :: NSPEC_DOUBLING_AB

  ! doubling regions
  integer :: isubregion,max_subregion
  integer :: doubling_factor,max_doubling_factor
  integer :: idoubling_bottom,idoubling_top
  integer :: nz_top,nz_bottom,nz_doubling
  integer :: one_exclude

  ! number of elements horizontally in each slice (i.e. per processor)
  NEX_PER_PROC_XI = NEX_XI / NPROC_XI
  NEX_PER_PROC_ETA = NEX_ETA / NPROC_ETA

  ! total number of processors
  NPROC = NPROC_XI * NPROC_ETA

  ! initializes outputs
  NSPEC_AB = 0
  NGLOB_AB = 0

  NSPEC2D_A_XI = 0
  NSPEC2D_B_XI = 0
  NSPEC2D_A_ETA = 0
  NSPEC2D_B_ETA = 0

  NSPEC2DMAX_XMIN_XMAX = 0
  NSPEC2DMAX_YMIN_YMAX = 0
  NSPEC2D_BOTTOM = 0
  NSPEC2D_TOP = 0
  NPOIN2DMAX_XMIN_XMAX = 0
  NPOIN2DMAX_YMIN_YMAX = 0

  !
  !--- case of a regular mesh
  !
  if (USE_REGULAR_MESH) then

    ! exact number of spectral elements without doubling layers
    NSPEC_NO_DOUBLING = NEX_XI * NEX_ETA * NER

    ! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%

    ! exact number of surface elements for a chunk without doubling layers

    NSPEC2D_NO_DOUBLING_XI = NEX_PER_PROC_XI*NER

    NSPEC2D_NO_DOUBLING_ETA = NEX_PER_PROC_ETA*NER

    ! exact number of spectral elements per slice
    NSPEC_AB = NER * NEX_PER_PROC_XI * NEX_PER_PROC_ETA

    ! exact number of surface elements for faces A and B along XI and ETA
    NSPEC2D_A_XI = NSPEC2D_NO_DOUBLING_XI
    NSPEC2D_B_XI = NSPEC2D_NO_DOUBLING_XI
    NSPEC2D_A_ETA = NSPEC2D_NO_DOUBLING_ETA
    NSPEC2D_B_ETA = NSPEC2D_NO_DOUBLING_ETA

    ! exact number of surface elements on the bottom and top boundaries
    ! and theoretical number of spectral elements in radial direction

    NSPEC2D_TOP = NEX_XI * NEX_ETA / NPROC
    NSPEC2D_BOTTOM = NSPEC2D_TOP

    !NSPEC1D_RADIAL_BEDROCK = NER

    ! face with max number of elements is type B here
    ! maximum number of surface elements on vertical boundaries of the slices
    NSPEC2DMAX_XMIN_XMAX = NSPEC2D_B_ETA
    NSPEC2DMAX_YMIN_YMAX = NSPEC2D_B_XI

    ! theoretical number of Gauss-Lobatto points in radial direction
    !NPOIN1D_RADIAL_BEDROCK = NSPEC1D_RADIAL_BEDROCK*(NGLLZ_M-1)+1

    ! 2-D addressing and buffers for summation between slices
    ! we add one to number of points because of the flag after the last point
    NPOIN2DMAX_XMIN_XMAX = NSPEC2DMAX_XMIN_XMAX*NGLLY_M*NGLLZ_M + 1
    NPOIN2DMAX_YMIN_YMAX = NSPEC2DMAX_YMIN_YMAX*NGLLX_M*NGLLZ_M + 1

    ! exact number of global points per slice
    NGLOB_AB = (NEX_PER_PROC_XI*(NGLLX_M-1)+1) * (NEX_PER_PROC_ETA*(NGLLY_M-1)+1) * (NER*(NGLLZ_M-1)+1)

  else

    !
    !--- case of a non-regular mesh with mesh doublings
    !

    ! initializes element counters
    NSPEC_NO_DOUBLING = 0
    NSPEC_DOUBLING_AB = 0

    ! initializes surface element counters
    NSPEC2D_NO_DOUBLING_XI = 0
    NSPEC2D_NO_DOUBLING_ETA = 0

    NSPEC2D_DOUBLING_A_XI = 0
    NSPEC2D_DOUBLING_A_ETA = 0

    NSPEC2D_DOUBLING_B_XI = 0
    NSPEC2D_DOUBLING_B_ETA = 0

    ! initializes global point counters
    NGLOB_NO_DOUBLING = 0
    NGLOB_DOUBLING = 0

    ! isubregion:
    !   example NDOUBLINGS = 1: isubregion = 1,2,3
    !   example NDOUBLINGS = 2: isubregion = 1,2,3,4,5
    ! the values goes from 1 to 2 * NDOUBLINGS + 1
    !
    ! maximum number of subregions: isubregion is in range [1,max_subregion]
    max_subregion = 2 * NDOUBLINGS + 1

    ! loops over all subregions
    do isubregion = 1,max_subregion
      ! factor for doubling for subregion:
      !   for NDOUBLINGS = 1, decreases 2, 2, 1 for subregion 1,2,3
      !   for NDOUBLINGS = 2, decreases 4, 4, 2, 2, 1 for subregion 1,2,3,4,5
      !   for NDOUBLINGS = 3, decreases 8,.., 4,.., 2,.., 1 for subregion 1,..,3,..,5,..7
      doubling_factor = 2**( int((max_subregion - isubregion + 1) / 2) )

      ! debug
      !print *,'subregion',isubregion,'doubling_factor = ',doubling_factor

      ! doubling layer index:
      !   example NDOUBLINGS = 2 -> for isubregion = 2: take only 1. doubling layer ner_doubling(1)
      !   example NDOUBLINGS = 3 -> for isubregion = 2: take 2. doubling layer ner_doubling(2) (counting from bottom up)
      !                                            = 4: take 1. doubling layer ner_doubling(1)
      if (NDOUBLINGS == 1) then
        ! single doubling layer
        idoubling_bottom = 1
        idoubling_top = 1
      else
        ! bottom doubling zone (current doubling region index)
        idoubling_bottom = NDOUBLINGS - int(isubregion / 2) + 1
        ! top doubling zone
        idoubling_top = idoubling_bottom - 1
      endif

      ! debug
      !print *,'idoubling_bottom = ',idoubling_bottom,'idoubling_top = ',idoubling_top

      ! determines increments
      if (modulo(isubregion,2) == 1) then
        ! regular layer
        ! bottom and top element layer number
        nz_bottom = 0
        if (isubregion > 1) nz_bottom = ner_doublings(idoubling_bottom)

        nz_top = NER - 1
        if (isubregion < max_subregion) nz_top = ner_doublings(idoubling_top) - 2 - 1

        ! number of regular element layers in vertical direction
        NER_REGULAR_SINGLE = nz_top - nz_bottom + 1

        ! elements in no doubling regions
        NSPEC_NO_DOUBLING = NSPEC_NO_DOUBLING &
          + (NEX_XI/doubling_factor) * (NEX_ETA/doubling_factor) * NER_REGULAR_SINGLE

        ! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%
        ! exact number of surface elements for a chunk without doubling layers
        NSPEC2D_NO_DOUBLING_XI = NSPEC2D_NO_DOUBLING_XI  &
          + (NEX_PER_PROC_XI/doubling_factor) * NER_REGULAR_SINGLE

        NSPEC2D_NO_DOUBLING_ETA = NSPEC2D_NO_DOUBLING_ETA  &
          + (NEX_PER_PROC_ETA/doubling_factor) * NER_REGULAR_SINGLE

        ! exact number of global points
        ! compute number of points in blocks with no doubling
        ! exclude the surfaces in contact with the doubling regions
        one_exclude = 0
        if (isubregion > 1 .and. isubregion < max_subregion) one_exclude = 1

        NGLOB_NO_DOUBLING = NGLOB_NO_DOUBLING &
          + (NEX_PER_PROC_XI/doubling_factor + 1) * (NEX_PER_PROC_ETA/doubling_factor + 1) * (NER_REGULAR_SINGLE - one_exclude)

      else
        ! doubling layer
        ! element layer number of doubling zone
        nz_doubling = ner_doublings(idoubling_bottom) - 2

        ! number of regular element layers in vertical direction
        NUM_DOUBLING_BRICKS = NEX_XI/(doubling_factor * 2) * NEX_ETA/(doubling_factor * 2)

        ! number of elementary bricks in the regions with doubling
        ! for type AB, each doubling brick contains 32 elements on 2 levels
        NSPEC_DOUBLING_AB =  NSPEC_DOUBLING_AB + 32 * NUM_DOUBLING_BRICKS

        ! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%
        ! exact number of surface elements in the doubling regions
        ! number of elementary bricks in regions with doubling
        !
        ! note: NSPEC2D_** values are only used to estimate the maximum values for NSPEC2DMAX_**
        !       thus, an overestimation is still fine
        NUM2D_DOUBLING_BRICKS_XI = 2 * NEX_PER_PROC_XI/(doubling_factor * 2)

        NUM2D_DOUBLING_BRICKS_ETA = 2 * NEX_PER_PROC_ETA/(doubling_factor * 2)

        ! old way:
        ! daniel: not sure where these factors came from, and why they change from one doubling to two doubling layers.
        !         i will take a factor 10 to estimate surfaces for a doubling brick
        !
        ! NDOUBLINGS == 2
        !NSPEC2D_DOUBLING_A_XI = 10 * NUM2D_DOUBLING_BRICKS_XI
        !NSPEC2D_DOUBLING_A_ETA = 10 * NUM2D_DOUBLING_BRICKS_ETA
        !
        !NSPEC2D_DOUBLING_B_XI = 10 * NUM2D_DOUBLING_BRICKS_XI
        !NSPEC2D_DOUBLING_B_ETA = 10 * NUM2D_DOUBLING_BRICKS_ETA

        ! NDOUBLINGS == 3
        !NSPEC2D_DOUBLING_A_XI = 7 * NUM2D_DOUBLING_BRICKS_XI
        !NSPEC2D_DOUBLING_A_ETA = 7 * NUM2D_DOUBLING_BRICKS_ETA
        !
        !NSPEC2D_DOUBLING_B_XI = 10 * NUM2D_DOUBLING_BRICKS_XI
        !NSPEC2D_DOUBLING_B_ETA = 12 * NUM2D_DOUBLING_BRICKS_ETA

        NSPEC2D_DOUBLING_A_XI = NSPEC2D_DOUBLING_A_XI + 10 * NUM2D_DOUBLING_BRICKS_XI
        NSPEC2D_DOUBLING_A_ETA = NSPEC2D_DOUBLING_A_ETA + 10 * NUM2D_DOUBLING_BRICKS_ETA

        NSPEC2D_DOUBLING_B_XI = NSPEC2D_DOUBLING_B_XI + 10 * NUM2D_DOUBLING_BRICKS_XI
        NSPEC2D_DOUBLING_B_ETA = NSPEC2D_DOUBLING_B_ETA  + 10 * NUM2D_DOUBLING_BRICKS_ETA

        ! exact number of global points
        ! case of the doubling regions
        ! current, single doubling layer points
        NGLOB_DOUBLING_SINGLE = &
          (NEX_PER_PROC_XI/(doubling_factor * 2))*(NEX_PER_PROC_ETA/(doubling_factor * 2)) &
              *(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
          - ((NEX_PER_PROC_XI/(doubling_factor * 2)-1)*(NEX_PER_PROC_ETA/(doubling_factor * 2)) &
          + (NEX_PER_PROC_ETA/(doubling_factor * 2)-1)*(NEX_PER_PROC_XI/(doubling_factor * 2)))&
              * (8*NGLLX_M**2 - 11*NGLLX_M + 4) &
          + (NEX_PER_PROC_XI/(doubling_factor * 2)-1)*(NEX_PER_PROC_ETA/(doubling_factor * 2)-1)*(NGLLX_M+1)

        ! total doubling layer points
        NGLOB_DOUBLING = NGLOB_DOUBLING + NGLOB_DOUBLING_SINGLE

      endif
    enddo ! subregions

    ! exact number of spectral elements per slice
    NSPEC_AB = (NSPEC_NO_DOUBLING + NSPEC_DOUBLING_AB) / NPROC

    ! (exact) number of surface elements for faces A and B
    ! along XI and ETA for all regions
    !
    ! note: NSPEC2D_** values are only used to estimate the maximum values for NSPEC2DMAX_**
    NSPEC2D_A_XI = NSPEC2D_NO_DOUBLING_XI + NSPEC2D_DOUBLING_A_XI
    NSPEC2D_B_XI = NSPEC2D_NO_DOUBLING_XI + NSPEC2D_DOUBLING_B_XI

    NSPEC2D_A_ETA = NSPEC2D_NO_DOUBLING_ETA + NSPEC2D_DOUBLING_A_ETA
    NSPEC2D_B_ETA = NSPEC2D_NO_DOUBLING_ETA + NSPEC2D_DOUBLING_B_ETA

    ! maximum doubling factor
    max_doubling_factor = 2**( int(max_subregion / 2) )

    ! exact number of surface elements on the bottom and top boundaries
    NSPEC2D_TOP = NEX_XI * NEX_ETA / NPROC
    NSPEC2D_BOTTOM = (NEX_XI/max_doubling_factor) * (NEX_ETA/max_doubling_factor) / NPROC

    ! face with max number of elements (used to be type B)
    ! maximum number of surface elements on vertical boundaries of the slices
    NSPEC2DMAX_XMIN_XMAX = NSPEC2D_B_ETA
    NSPEC2DMAX_YMIN_YMAX = NSPEC2D_B_XI

    ! 2-D addressing and buffers for summation between slices
    ! we add one to number of points because of the flag after the last point
    NPOIN2DMAX_XMIN_XMAX = NSPEC2DMAX_XMIN_XMAX * NGLLY_M * NGLLZ_M + 1
    NPOIN2DMAX_YMIN_YMAX = NSPEC2DMAX_YMIN_YMAX * NGLLX_M * NGLLZ_M + 1

    ! exact number of global points
    NGLOB_AB = NGLOB_DOUBLING + NGLOB_NO_DOUBLING

! > original hard-coded way
!    if (NDOUBLINGS == 1) then
!
!      ! number of spectral elements at the bottom of the doubling
!      NEX_DOUBLING_ABOVE_XI = NEX_XI/2
!      NEX_DOUBLING_ABOVE_ETA = NEX_ETA/2
!      NEX_DOUBLING_ABOVE_PER_PROC_XI = NEX_PER_PROC_XI/2
!      NEX_DOUBLING_ABOVE_PER_PROC_ETA = NEX_PER_PROC_ETA/2
!
!      ! exact number of spectral elements without doubling layers
!      NER_REGULAR2 = NER - ner_doublings(1)
!      NER_REGULAR1 = ner_doublings(1) - 2
!
!      NSPEC_NO_DOUBLING = &
!          (NEX_XI/2)*(NEX_ETA/2)*NER_REGULAR1 + NEX_XI*NEX_ETA*NER_REGULAR2
!
!      ! exact number of spectral elements in the doubling regions
!
!      ! number of elementary bricks in the two regions with doubling
!      NUM_DOUBLING_BRICKS = (NEX_DOUBLING_ABOVE_XI*NEX_DOUBLING_ABOVE_ETA)/4
!
!      ! for type AB, each doubling brick contains 32 elements on 2 levels
!      NSPEC_DOUBLING_AB=32*NUM_DOUBLING_BRICKS
!
!      ! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%
!
!      ! exact number of surface elements for a chunk without doubling layers
!
!      NSPEC2D_NO_DOUBLING_XI = &
!          (NEX_PER_PROC_XI/2)*NER_REGULAR1 + NEX_PER_PROC_XI*NER_REGULAR2
!
!      NSPEC2D_NO_DOUBLING_ETA = &
!          (NEX_PER_PROC_ETA/2)*NER_REGULAR1 + NEX_PER_PROC_ETA*NER_REGULAR2
!
!      ! exact number of surface elements in the doubling regions
!
!      ! number of elementary bricks in the two regions with doubling
!      NUM2D_DOUBLING_BRICKS_XI = (NEX_PER_PROC_XI/2 &
!          + NEX_DOUBLING_ABOVE_PER_PROC_XI)/2
!
!      NUM2D_DOUBLING_BRICKS_ETA = (NEX_PER_PROC_ETA/2 &
!          + NEX_DOUBLING_ABOVE_PER_PROC_ETA)/2
!
!      ! for type A, each doubling brick contains 10 elements on 3 levels
!      NSPEC2D_DOUBLING_A_XI=10*NUM2D_DOUBLING_BRICKS_XI
!      NSPEC2D_DOUBLING_A_ETA=10*NUM2D_DOUBLING_BRICKS_ETA
!
!      ! for type B, each doubling brick contains 12 elements on 3 levels
!      NSPEC2D_DOUBLING_B_XI=10*NUM2D_DOUBLING_BRICKS_XI
!      NSPEC2D_DOUBLING_B_ETA=10*NUM2D_DOUBLING_BRICKS_ETA
!
!      ! exact number of spectral elements
!      NSPEC_AB = (NSPEC_NO_DOUBLING + NSPEC_DOUBLING_AB) / NPROC
!
!      ! exact number of surface elements for faces A and B
!      ! along XI and ETA for doubling region
!      NSPEC2D_A_XI = NSPEC2D_NO_DOUBLING_XI + NSPEC2D_DOUBLING_A_XI
!      NSPEC2D_B_XI = NSPEC2D_NO_DOUBLING_XI + NSPEC2D_DOUBLING_B_XI
!      NSPEC2D_A_ETA = NSPEC2D_NO_DOUBLING_ETA + NSPEC2D_DOUBLING_A_ETA
!      NSPEC2D_B_ETA = NSPEC2D_NO_DOUBLING_ETA + NSPEC2D_DOUBLING_B_ETA
!
!      ! exact number of surface elements on the bottom and top boundaries
!      ! and theoretical number of spectral elements in radial direction
!
!      NSPEC2D_TOP = NEX_XI*NEX_ETA / NPROC
!      NSPEC2D_BOTTOM = (NEX_XI/2)*(NEX_ETA/2) / NPROC
!
!      ! face with max number of elements is type B here
!      ! maximum number of surface elements on vertical boundaries of the slices
!      NSPEC2DMAX_XMIN_XMAX = NSPEC2D_B_ETA
!      NSPEC2DMAX_YMIN_YMAX = NSPEC2D_B_XI
!
!      ! theoretical number of Gauss-Lobatto points in radial direction
!      !  NPOIN1D_RADIAL_BEDROCK = NSPEC1D_RADIAL_BEDROCK*(NGLLZ_M-1)+1
!
!      ! 2-D addressing and buffers for summation between slices
!      ! we add one to number of points because of the flag after the last point
!      NPOIN2DMAX_XMIN_XMAX = NSPEC2DMAX_XMIN_XMAX*NGLLY_M*NGLLZ_M + 1
!      NPOIN2DMAX_YMIN_YMAX = NSPEC2DMAX_YMIN_YMAX*NGLLX_M*NGLLZ_M + 1
!
!      ! exact number of global points
!
!      ! case of the doubling regions
!
!      NGLOB_DOUBLING_1 = (NEX_PER_PROC_XI/4)*(NEX_PER_PROC_ETA/4)*(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
!          - ((NEX_PER_PROC_XI/4-1)*(NEX_PER_PROC_ETA/4) &
!          + (NEX_PER_PROC_ETA/4-1)*(NEX_PER_PROC_XI/4)) * (8*NGLLX_M**2-11*NGLLX_M+4) &
!          + (NEX_PER_PROC_XI/4-1)*(NEX_PER_PROC_ETA/4-1)*(NGLLX_M+1)
!
!      ! NGLOB_DOUBLING_1 = (NEX_PER_PROC_XI/2)*(NEX_PER_PROC_ETA/2)*(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
!      !  - ((NEX_PER_PROC_XI/2-1)*(NEX_PER_PROC_ETA/2) &
!      !  + (NEX_PER_PROC_ETA/2-1)*(NEX_PER_PROC_XI/2)) * (8*NGLLX_M**2-11*NGLLX_M+4)&
!      !  + (NEX_PER_PROC_XI/2-1)*(NEX_PER_PROC_ETA/2-1)*NGLLX_M
!
!      ! NGLOB_DOUBLING_2 =(NEX_PER_PROC_XI/8)*(NEX_PER_PROC_ETA/8)*(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
!      !  - ((NEX_PER_PROC_XI/8-1)*(NEX_PER_PROC_ETA/8) &
!      !  + (NEX_PER_PROC_ETA/8-1)*(NEX_PER_PROC_XI/8)) * (8*NGLLX_M**2-11*NGLLX_M+4)&
!      !  + (NEX_PER_PROC_XI/8-1)*(NEX_PER_PROC_ETA/8-1)*NGLLX_M
!
!      NGLOB_DOUBLING = NGLOB_DOUBLING_1 !+ NGLOB_DOUBLING_2
!
!      ! compute number of points in blocks with no doubling
!      ! exclude the four surfaces in contact with the doubling regions
!
!      NGLOB_NO_DOUBLING = (NEX_PER_PROC_XI/2 + 1)*(NEX_PER_PROC_ETA/2 + 1)*NER_REGULAR1 &
!          + (NEX_PER_PROC_XI + 1)*(NEX_PER_PROC_ETA + 1)*NER_REGULAR2
!
!      NGLOB_AB = NGLOB_DOUBLING + NGLOB_NO_DOUBLING
!
!    else if (NDOUBLINGS == 2) then
!
!      ! number of spectral elements at the bottom of the doubling below the moho
!      NEX_DOUBLING_ABOVE_XI=NEX_XI/2
!      NEX_DOUBLING_ABOVE_ETA=NEX_ETA/2
!      NEX_DOUBLING_ABOVE_PER_PROC_XI=NEX_PER_PROC_XI/2
!      NEX_DOUBLING_ABOVE_PER_PROC_ETA=NEX_PER_PROC_ETA/2
!
!      ! exact number of spectral elements without doubling layers
!      NER_REGULAR3 = NER - ner_doublings(1)
!      NER_REGULAR2 = (ner_doublings(1) - 2) - ner_doublings(2)
!      NER_REGULAR1 = ner_doublings(2) - 2
!
!      NSPEC_NO_DOUBLING = &
!          (NEX_XI/4)*(NEX_ETA/4)*NER_REGULAR1 &
!          + (NEX_XI/2)*(NEX_ETA/2)*NER_REGULAR2 &
!          + NEX_XI*NEX_ETA*NER_REGULAR3
!
!      ! exact number of spectral elements in the doubling regions
!
!      ! number of elementary bricks in the two regions with doubling
!      NUM_DOUBLING_BRICKS = (NEX_DOUBLING_ABOVE_XI*NEX_DOUBLING_ABOVE_ETA)/4&
!          + (NEX_DOUBLING_ABOVE_XI*NEX_DOUBLING_ABOVE_ETA)/16
!
!      ! for type AB, each doubling brick contains 32 elements on 2 levels
!      NSPEC_DOUBLING_AB=32*NUM_DOUBLING_BRICKS
!
!
!      ! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%
!
!      ! exact number of surface elements for a chunk without doubling layers
!
!      NSPEC2D_NO_DOUBLING_XI = &
!          (NEX_PER_PROC_XI/4)*NER_REGULAR1 &
!          + (NEX_PER_PROC_XI/2)*NER_REGULAR2 &
!          + NEX_PER_PROC_XI*NER_REGULAR3
!
!      NSPEC2D_NO_DOUBLING_ETA = &
!          (NEX_PER_PROC_ETA/4)*NER_REGULAR1 &
!          + (NEX_PER_PROC_ETA/2)*NER_REGULAR2 &
!          + NEX_PER_PROC_ETA*NER_REGULAR3
!
!      ! exact number of surface elements in the doubling regions
!
!      ! number of elementary bricks in the two regions with doubling
!      NUM2D_DOUBLING_BRICKS_XI = (NEX_PER_PROC_XI/2 &
!          + NEX_DOUBLING_ABOVE_PER_PROC_XI)/2
!
!      NUM2D_DOUBLING_BRICKS_ETA = (NEX_PER_PROC_ETA/2 &
!          + NEX_DOUBLING_ABOVE_PER_PROC_ETA)/2
!
!      ! for type A, each doubling brick contains 10 elements on 3 levels
!      NSPEC2D_DOUBLING_A_XI=7*NUM2D_DOUBLING_BRICKS_XI
!      NSPEC2D_DOUBLING_A_ETA=7*NUM2D_DOUBLING_BRICKS_ETA
!
!      ! for type B, each doubling brick contains 12 elements on 3 levels
!      NSPEC2D_DOUBLING_B_XI=10*NUM2D_DOUBLING_BRICKS_XI
!      NSPEC2D_DOUBLING_B_ETA=12*NUM2D_DOUBLING_BRICKS_ETA
!
!      ! exact number of spectral elements
!      NSPEC_AB = (NSPEC_NO_DOUBLING + NSPEC_DOUBLING_AB) / NPROC
!
!      ! exact number of surface elements for faces A and B
!      ! along XI and ETA for doubling region
!      NSPEC2D_A_XI = NSPEC2D_NO_DOUBLING_XI + NSPEC2D_DOUBLING_A_XI
!      NSPEC2D_B_XI = NSPEC2D_NO_DOUBLING_XI + NSPEC2D_DOUBLING_B_XI
!      NSPEC2D_A_ETA = NSPEC2D_NO_DOUBLING_ETA + NSPEC2D_DOUBLING_A_ETA
!      NSPEC2D_B_ETA = NSPEC2D_NO_DOUBLING_ETA + NSPEC2D_DOUBLING_B_ETA
!
!      ! exact number of surface elements on the bottom and top boundaries
!      ! and theoretical number of spectral elements in radial direction
!
!      NSPEC2D_TOP = NEX_XI*NEX_ETA / NPROC
!      NSPEC2D_BOTTOM = (NEX_XI/4)*(NEX_ETA/4) / NPROC
!
!      !  NSPEC1D_RADIAL_BEDROCK = (NER_BASEMENT_SEDIM+NER_16_BASEMENT+NER_MOHO_16)/2 + NER_BOTTOM_MOHO/4
!
!      ! face with max number of elements is type B here
!      ! maximum number of surface elements on vertical boundaries of the slices
!      NSPEC2DMAX_XMIN_XMAX = NSPEC2D_B_ETA
!      NSPEC2DMAX_YMIN_YMAX = NSPEC2D_B_XI
!
!      !debug
!      !print *,'nspec minmax:',NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_yMAX
!
!      ! theoretical number of Gauss-Lobatto points in radial direction
!      !  NPOIN1D_RADIAL_BEDROCK = NSPEC1D_RADIAL_BEDROCK*(NGLLZ_M-1)+1
!
!      ! 2-D addressing and buffers for summation between slices
!      ! we add one to number of points because of the flag after the last point
!      NPOIN2DMAX_XMIN_XMAX = NSPEC2DMAX_XMIN_XMAX*NGLLY_M*NGLLZ_M + 1
!      NPOIN2DMAX_YMIN_YMAX = NSPEC2DMAX_YMIN_YMAX*NGLLX_M*NGLLZ_M + 1
!
!
!      ! exact number of global points
!
!      ! case of the doubling regions
!
!      NGLOB_DOUBLING_1 = (NEX_PER_PROC_XI/8)*(NEX_PER_PROC_ETA/8)*(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
!          - ((NEX_PER_PROC_XI/8-1)*(NEX_PER_PROC_ETA/8) &
!          + (NEX_PER_PROC_ETA/8-1)*(NEX_PER_PROC_XI/8)) * (8*NGLLX_M**2-11*NGLLX_M+4)&
!          + (NEX_PER_PROC_XI/8-1)*(NEX_PER_PROC_ETA/8-1)*(NGLLX_M+1)
!
!      NGLOB_DOUBLING_2 = (NEX_PER_PROC_XI/4)*(NEX_PER_PROC_ETA/4)*(32*NGLLX_M**3 - 70*NGLLX_M**2 + 52*NGLLX_M - 13) &
!          - ((NEX_PER_PROC_XI/4-1)*(NEX_PER_PROC_ETA/4) &
!          + (NEX_PER_PROC_ETA/4-1)*(NEX_PER_PROC_XI/4)) * (8*NGLLX_M**2-11*NGLLX_M+4)&
!          + (NEX_PER_PROC_XI/4-1)*(NEX_PER_PROC_ETA/4-1)*(NGLLX_M+1)
!
!      NGLOB_DOUBLING = NGLOB_DOUBLING_1 + NGLOB_DOUBLING_2
!
!      ! compute number of points in blocks with no doubling
!      ! exclude the four surfaces in contact with the doubling regions
!
!      NGLOB_NO_DOUBLING = (NEX_PER_PROC_XI/4 + 1)*(NEX_PER_PROC_ETA/4 + 1)*NER_REGULAR1 &
!          + (NEX_PER_PROC_XI/2 + 1)*(NEX_PER_PROC_ETA/2 + 1)*(NER_REGULAR2 - 1) &
!          + (NEX_PER_PROC_XI + 1)*(NEX_PER_PROC_ETA + 1)*NER_REGULAR3
!
!      NGLOB_AB = NGLOB_DOUBLING + NGLOB_NO_DOUBLING
!
!    else
!      stop 'Error NDOUBLINGS value not supported yet'
!    endif
! < original hardcoded way

  endif ! end of section for non-regular mesh with doublings

end subroutine compute_parameters

