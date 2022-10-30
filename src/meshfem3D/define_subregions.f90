!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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


  subroutine define_model_regions(NEX_PER_PROC_XI,NEX_PER_PROC_ETA,iproc_xi,iproc_eta, &
                                  isubregion,nbsubregions,subregions, &
                                  iaddx,iaddy,iaddz,ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
                                  num_material)

  use constants

  implicit none

  integer,intent(in) :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer,intent(in) :: iproc_xi,iproc_eta

  !  definition of the different regions of the model in the mesh (nx,ny,nz)
  !  #1 #2 : nx_begining,nx_end
  !  #3 #4 : ny_begining,ny_end
  !  #5 #6 : nz_begining,nz_end
  !     #7 : material number
  integer,intent(in) :: isubregion,nbsubregions
  integer,intent(in) :: subregions(nbsubregions,7)

  integer,intent(out) :: ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir
  integer,intent(out) :: iax,iay,iar

  ! topology of the elements
  integer,intent(inout) :: iaddx(NGNOD_EIGHT_CORNERS)
  integer,intent(inout) :: iaddy(NGNOD_EIGHT_CORNERS)
  integer,intent(inout) :: iaddz(NGNOD_EIGHT_CORNERS)

  integer,intent(out) :: num_material

  call usual_hex_nodes(NGNOD_EIGHT_CORNERS,iaddx,iaddy,iaddz)

  ix1=2*(subregions(isubregion,1) - iproc_xi*NEX_PER_PROC_XI - 1)
  if (ix1 < 0) ix1 = 0

  ix2=2*(subregions(isubregion,2) - iproc_xi*NEX_PER_PROC_XI - 1)
  if (ix2 > 2*(NEX_PER_PROC_XI - 1)) ix2 = 2*(NEX_PER_PROC_XI - 1)
  dix=2

  iy1=2*(subregions(isubregion,3) - iproc_eta*NEX_PER_PROC_ETA - 1)
  if (iy1 < 0) iy1 = 0

  iy2=2*(subregions(isubregion,4) - iproc_eta*NEX_PER_PROC_ETA - 1)
  if (iy2 > 2*(NEX_PER_PROC_ETA - 1)) iy2 = 2*(NEX_PER_PROC_ETA - 1)
  diy=2

  ir1=2*(subregions(isubregion,5) - 1)
  ir2=2*(subregions(isubregion,6) - 1)
  dir=2

  iax=1
  iay=1
  iar=1

  num_material = subregions(isubregion,7)

  end subroutine define_model_regions

!
!------------------------------------------------------------------------------------
!

  subroutine define_mesh_regions(USE_REGULAR_MESH,isubregion,NER, &
                                 NEX_PER_PROC_XI,NEX_PER_PROC_ETA,iproc_xi,iproc_eta, &
                                 NDOUBLINGS,ner_doublings, &
                                 iaddx,iaddy,iaddz,ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar)

  use constants

  implicit none

  logical,intent(in) :: USE_REGULAR_MESH

  integer,intent(in) :: isubregion
  integer,intent(in) :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER
  integer,intent(in) :: iproc_xi,iproc_eta

  integer,intent(out) :: ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir
  integer,intent(out) :: iax,iay,iar

  integer,intent(in) :: NDOUBLINGS
  integer,dimension(NDOUBLINGS),intent(in) :: ner_doublings

  ! topology of the elements
  integer,intent(inout) :: iaddx(NGNOD_EIGHT_CORNERS)
  integer,intent(inout) :: iaddy(NGNOD_EIGHT_CORNERS)
  integer,intent(inout) :: iaddz(NGNOD_EIGHT_CORNERS)

  ! local parameters
  integer :: doubling_factor
  integer :: idoubling_bottom,idoubling_top
  integer :: nz_top,nz_bottom,nz_doubling
  integer :: max_subregion

  ! hex node addressing uses units of 2 for corners
  integer,parameter :: IUNIT = 2

  ! to avoid compiler warnings
  integer :: idummy

  idummy = iproc_xi
  idummy = iproc_eta

! **************

  ! sets up node addressing
  ! note: usual hex nodes addressing uses units of 2 for corner nodes
  call usual_hex_nodes(NGNOD_EIGHT_CORNERS,iaddx,iaddy,iaddz)

!
!--- case of a regular mesh
!
  if (USE_REGULAR_MESH) then
    ! loop over elements in xi-direction
    ix1 = 0
    ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
    dix = IUNIT

    ! loop over elements in eta-direction
    iy1 = 0
    iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
    diy = IUNIT

    ! loop over elements in vertical direction
    ! (takes full region)
    ir1 = 0
    ir2 = IUNIT*(NER - 1)
    dir = IUNIT

    ! element increments
    iax = 1
    iay = 1
    iar = 1

  else

    ! isubregion:
    !   example NDOUBLINGS = 1: isubregion = 1,2,3
    !   example NDOUBLINGS = 2: isubregion = 1,2,3,4,5
    ! the values goes from 1 to 2 * NDOUBLINGS + 1
    !
    ! maximum number of subregions: isubregion is in range [1,max_subregion]
    max_subregion = 2 * NDOUBLINGS + 1

    ! unit increment: element addressing uses unit length 2
    !   regular element      -> di = 2 * 1
    !   double-size element  -> di = 2 * 2 = 4 = IUNIT * doubling_factor
    !   4-fold sized element -> di = 2 * 4 = 8
    !
    ! example NDOUBLINGS = 1:
    !   -> increment for 1. subregion = 2 (double)  = 2 ** 1 = 2 ** [( 3 - 1 + 1) / 2]
    !                    2. subregion = 2 doubling  = 2 ** 1 = 2 ** [( 3 - 2 + 1) / 2]
    !                    3. subregion = 1 (regular) = 2 ** 0 = 2 ** [( 3 - 3 + 1) / 2]
    !
    ! example for regular regions only:
    ! NDOUBLINGS = 2:
    !   -> increment for 1. subregion = 4 (4-fold)  = 2 ** 2 = 2 ** [( 5 - 1) / 2]
    !                    3. subregion = 2 (double)  = 2 ** 1 = 2 ** [( 5 - 3) / 2]
    !                    5. subregion = 1 (regular) = 2 ** 0 = 2 ** [( 5 - 5) / 2]
    ! NDOUBLINGS = 3:
    !   -> increment for 1. subregion = 8 (8-fold)  = 2 ** [( 7 - 1) / 2]
    !                    3. subregion = 4 (4-fold)  = 2 ** [( 7 - 3) / 2]
    !                    5. subregion = 2 (double)  = 2 ** [( 7 - 5) / 2]
    !                    7. subregion = 1 (regular) = 2 ** [( 7 - 7) / 2]
    !
    ! factor for doubling: 1, 2, 4, 8, ..
    doubling_factor = 2**( int((max_subregion - isubregion + 1) / 2) )

    ! debug
    !if (myrank == 0) print *,'subregion',isubregion,'doubling_factor = ',doubling_factor

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
    !if (myrank == 0) print *,'  idoubling_bottom = ',idoubling_bottom,'idoubling_top = ',idoubling_top

    ! determines increments
    if (modulo(isubregion,2) == 1) then
      ! regular layer

      ! checks indices
      if (isubregion > 1) then
        if (idoubling_bottom < 1 .or. idoubling_bottom > NDOUBLINGS) &
          call exit_MPI(myrank,'Error invalid idoubling_bottom index in regular layer')
      endif
      if (isubregion < max_subregion) then
        if (idoubling_top < 1 .or. idoubling_top > NDOUBLINGS) &
          call exit_MPI(myrank,'Error invalid idoubling_top index in regular layer')
      endif

      ! bottom and top element layer number
      nz_bottom = 0
      if (isubregion > 1) nz_bottom = ner_doublings(idoubling_bottom)

      nz_top = NER - 1
      if (isubregion < max_subregion) nz_top = ner_doublings(idoubling_top) - 2 - 1

      ! debug
      !if (myrank == 0) print *,'  nz_bottom = ',nz_bottom,'nz_top = ',nz_top

      ! loop over elements in xi-direction
      ix1 = 0
      ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
      ! steps over 4 unit lengths (instead of just 2) to double element size
      dix = IUNIT * doubling_factor

      ! loop over elements in eta-direction
      iy1 = 0
      iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
      ! steps over 4 unit lengths (instead of just 2) to double element size
      diy = IUNIT * doubling_factor

      ! loop over elements in vertical direction
      ! only bottom up to one before doubling layer
      ir1 = IUNIT * nz_bottom
      ir2 = IUNIT * nz_top
      dir = IUNIT

      ! element increments
      iax = doubling_factor
      iay = doubling_factor
      iar = 1
    else
      ! doubling layer

      ! checks index
      if (idoubling_bottom < 1 .or. idoubling_bottom > NDOUBLINGS) &
        call exit_MPI(myrank,'Error invalid idoubling_bottom index in doubling layer')

      ! element layer number of doubling zone
      nz_doubling = ner_doublings(idoubling_bottom) - 2

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  doubling layer positioned at layer: ',ner_doublings(idoubling_bottom)
      endif

      ! debug
      !if (myrank == 0) print *,'  nz_doubling = ',nz_doubling

      ! loop over elements in xi-direction
      ix1 = 0
      ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
      dix = IUNIT * doubling_factor * 2

      ! loop over elements in eta-direction
      iy1 = 0
      iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
      diy = IUNIT * doubling_factor * 2

      ! loop over elements in vertical direction
      ir1 = IUNIT * nz_doubling
      ir2 = IUNIT * nz_doubling
      dir = IUNIT

      ! element increments
      iax = doubling_factor * 2
      iay = doubling_factor * 2
      iar = 2
    endif

    ! debug
    !print *,'ir1/ir2/dir = ',ir1,ir2,dir,' ix1/ix2/dix = ',ix1,ix2,dix, &
    !        ' iy1/iy2/diy = ',iy1,iy2,diy,' iax/iay/iar = ',iax,iay,iar


! > original hard-coded way
!    if (NDOUBLINGS == 1) then
!
!      select case (isubregion)
!      case (1)
!        ! regular region with double-sized elements
!        ! loop over elements in xi-direction
!        ix1 = 0
!        ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
!        ! steps over 4 unit lengths (instead of just 2) to double element size
!        dix = IUNIT * 2
!
!        ! loop over elements in eta-direction
!        iy1 = 0
!        iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
!        ! steps over 4 unit lengths (instead of just 2) to double element size
!        diy = IUNIT * 2
!
!        ! loop over elements in vertical direction
!        ! only bottom up to one before doubling layer
!        ir1 = 0
!        ir2 = IUNIT*(ner_doublings(1) - 2 -1)
!        dir = IUNIT
!
!        ! element increments
!        iax = 2
!        iay = 2
!        iar = 1
!
!      case (2)
!        ! doubling layer
!        ! loop over elements in xi-direction
!        ix1 = 0
!        ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
!        dix = IUNIT * 2 * 2
!
!        ! loop over elements in eta-direction
!        iy1 = 0
!        iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
!        diy = IUNIT * 2 * 2
!
!        ! loop over elements in vertical direction
!        ir1 = IUNIT*(ner_doublings(1) - 2)
!        ir2 = IUNIT*(ner_doublings(1) - 2)
!        dir = IUNIT
!
!        ! element increments
!        iax = 4
!        iay = 4
!        iar = 2
!
!      case (3)
!        ! regular region with regular-sized elements (to reach NEX count at surface)
!        ! loop over elements in xi-direction
!        ix1 = 0
!        ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
!        dix = IUNIT
!
!        ! loop over elements in eta-direction
!        iy1 = 0
!        iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
!        diy = IUNIT
!
!        ! loop over elements in vertical direction
!        ir1 = IUNIT*(ner_doublings(1))
!        ir2 = IUNIT*(NER - 1)
!        dir = IUNIT
!
!        ! element increments
!        iax = 1
!        iay = 1
!        iar = 1
!
!      case default
!         stop 'Wrong number of meshregions'
!
!      end select
!
!    else if (NDOUBLINGS == 2) then
!
!      select case (isubregion)
!      case (1)
!        ! regular region, bottom with 4-fold sized elements (4 * 2 unit lengths = 8 di increments)
!        ! loop over elements in xi-direction
!        ix1 = 0
!        ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
!        dix = IUNIT * 4
!
!        ! loop over elements in eta-direction
!        iy1 = 0
!        iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
!        diy = IUNIT * 4
!
!        ! loop over elements in vertical direction
!        ir1 = 0
!        ir2 = IUNIT*(ner_doublings(2) - 2 -1)
!        dir = IUNIT
!
!        ! element increments
!        iax = 4
!        iay = 4
!        iar = 1
!
!      case (2)
!        ! 1. doubling layer
!        ! loop over elements in xi-direction
!        ix1 = 0
!        ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
!        dix = IUNIT * 4 * 2
!
!        ! loop over elements in eta-direction
!        iy1 = 0
!        iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
!        diy = IUNIT * 4 * 2
!
!        ! loop over elements in vertical direction
!        ! note: ner_doublings is sorted in decreasing order, e.g. (14,10),
!        !       i.e. ner_doublings(2) is smaller. we count from bottom up, it's thus the first doubling we reach
!        ir1 = IUNIT*(ner_doublings(2) - 2)
!        ir2 = IUNIT*(ner_doublings(2) - 2)
!        dir = IUNIT
!
!        ! element increments
!        iax = 8
!        iay = 8
!        iar = 2
!
!      case (3)
!        ! regular region, middle with double-sized elements
!        ! loop over elements in xi-direction
!        ix1 = 0
!        ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
!        dix = IUNIT * 2
!
!        ! loop over elements in eta-direction
!        iy1 = 0
!        iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
!        diy = IUNIT * 2
!
!        ! loop over elements in vertical direction
!        ir1 = IUNIT*ner_doublings(2)
!        ir2 = IUNIT*(ner_doublings(1) -2 -1)
!        dir = IUNIT
!
!        ! element increments
!        iax = 2
!        iay = 2
!        iar = 1
!
!      case (4)
!        ! 2. doubling layer
!        ! loop over elements in xi-direction
!        ix1 = 0
!        ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
!        dix = IUNIT * 2 * 2
!
!        ! loop over elements in eta-direction
!        iy1 = 0
!        iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
!        diy = IUNIT * 2 * 2
!
!        ! loop over elements in vertical direction
!        ir1 = IUNIT*(ner_doublings(1) - 2)
!        ir2 = IUNIT*(ner_doublings(1) - 2)
!        dir = IUNIT
!
!        ! element increments
!        iax = 4
!        iay = 4
!        iar = 2
!
!      case (5)
!        ! regular region, top with regular sized elements (to reach NEX count)
!        ! loop over elements in xi-direction
!        ix1 = 0
!        ix2 = IUNIT*(NEX_PER_PROC_XI - 1)
!        dix = IUNIT
!
!        ! loop over elements in eta-direction
!        iy1 = 0
!        iy2 = IUNIT*(NEX_PER_PROC_ETA - 1)
!        diy = IUNIT
!
!        ! loop over elements in vertical direction
!        ir1 = IUNIT*(ner_doublings(1))
!        ir2 = IUNIT*(NER - 1)
!        dir = IUNIT
!
!        ! element increments
!        iax = 1
!        iay = 1
!        iar = 1
!
!      case default
!        stop 'Wrong number of meshregions'
!      end select
!
!    else
!      stop 'Wrong number of doublings'
!    endif
! < original hard-coded way

  endif

  end subroutine define_mesh_regions

