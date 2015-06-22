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

subroutine get_jacobian_discontinuities(myrank,ispec,ix_elem,iy_elem,rmin,rmax,r1,r2,r3,r4,r5,r6,r7,r8, &
                     xstore,ystore,zstore,dershape2D_bottom, &
                     ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                     normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                     ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                     NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,r_moho,r_400,r_670, &
                     is_superbrick,USE_ONE_LAYER_SB,ispec_superbrick,nex_eta_moho,HONOR_1D_SPHERICAL_MOHO)

  implicit none

  include 'constants.h'

  ! input
  integer myrank, ispec, ix_elem, iy_elem
  double precision rmin,rmax
  double precision xstore(NGLLX,NGLLY,NGLLZ)
  double precision ystore(NGLLX,NGLLY,NGLLZ)
  double precision zstore(NGLLX,NGLLY,NGLLZ)
  double precision dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY)
  integer NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670, nex_eta_moho, ispec_superbrick
  double precision r_moho, r_400, r_670
  logical :: is_superbrick, USE_ONE_LAYER_SB,HONOR_1D_SPHERICAL_MOHO

  ! output
  integer ispec2D_moho_top, ispec2D_moho_bot, ispec2D_400_top, ispec2D_400_bot, ispec2D_670_top, ispec2D_670_bot
  integer,dimension(NSPEC2D_MOHO) :: ibelm_moho_top, ibelm_moho_bot
  integer,dimension(NSPEC2D_400) :: ibelm_400_top, ibelm_400_bot
  integer,dimension(NSPEC2D_670) :: ibelm_670_top, ibelm_670_bot
  real(kind=CUSTOM_REAL) :: normal_moho(NDIM,NGLLX,NGLLY,NSPEC2D_MOHO), jacobian2D_moho(NGLLX,NGLLY,NSPEC2D_MOHO)
  real(kind=CUSTOM_REAL) :: normal_400(NDIM,NGLLX,NGLLY,NSPEC2D_400), jacobian2D_400(NGLLX,NGLLY,NSPEC2D_400)
  real(kind=CUSTOM_REAL) :: normal_670(NDIM,NGLLX,NGLLY,NSPEC2D_670), jacobian2D_670(NGLLX,NGLLY,NSPEC2D_670)

  ! local variables
  double precision, dimension(NGNOD2D) :: xelm2, yelm2, zelm2
  double precision :: r1, r2, r3, r4, r5, r6, r7, r8
  double precision :: target_moho_high, target_moho_low, target_400_high, target_400_low, target_670_high, target_670_low
  integer :: nele_sub_block, ispec_list(16), map_irem_ix_12(8), map_irem_ix_34(8), map_irem_iy_odd(8), map_irem_iy_even(8)
  integer :: map_isub_ix(4), map_isub_iy(4), map_ix(NSPEC_DOUBLING_SUPERBRICK),  map_iy(NSPEC_DOUBLING_SUPERBRICK)
  integer :: i, ispec_superbrick_current, isub_block, irem_block, irem_ix, irem_iy, ix,iy,ix_top,iy_top, ispec2D_moho_bot_map

  ! ======================


  ! find the coordinates of 9 nodes for the bottom surface element to compute the jacobian if needed
  xelm2(1)=xstore(1,1,1)
  yelm2(1)=ystore(1,1,1)
  zelm2(1)=zstore(1,1,1)
  xelm2(2)=xstore(NGLLX,1,1)
  yelm2(2)=ystore(NGLLX,1,1)
  zelm2(2)=zstore(NGLLX,1,1)
  xelm2(3)=xstore(NGLLX,NGLLY,1)
  yelm2(3)=ystore(NGLLX,NGLLY,1)
  zelm2(3)=zstore(NGLLX,NGLLY,1)
  xelm2(4)=xstore(1,NGLLY,1)
  yelm2(4)=ystore(1,NGLLY,1)
  zelm2(4)=zstore(1,NGLLY,1)
  xelm2(5)=xstore((NGLLX+1)/2,1,1)
  yelm2(5)=ystore((NGLLX+1)/2,1,1)
  zelm2(5)=zstore((NGLLX+1)/2,1,1)
  xelm2(6)=xstore(NGLLX,(NGLLY+1)/2,1)
  yelm2(6)=ystore(NGLLX,(NGLLY+1)/2,1)
  zelm2(6)=zstore(NGLLX,(NGLLY+1)/2,1)
  xelm2(7)=xstore((NGLLX+1)/2,NGLLY,1)
  yelm2(7)=ystore((NGLLX+1)/2,NGLLY,1)
  zelm2(7)=zstore((NGLLX+1)/2,NGLLY,1)
  xelm2(8)=xstore(1,(NGLLY+1)/2,1)
  yelm2(8)=ystore(1,(NGLLY+1)/2,1)
  zelm2(8)=zstore(1,(NGLLY+1)/2,1)
  xelm2(9)=xstore((NGLLX+1)/2,(NGLLY+1)/2,1)
  yelm2(9)=ystore((NGLLX+1)/2,(NGLLY+1)/2,1)
  zelm2(9)=zstore((NGLLX+1)/2,(NGLLY+1)/2,1)

! radii to determine if an element is on the discontinuity or not
  target_moho_high = r_moho * (ONE + SMALLVAL)
  target_moho_low = r_moho * (ONE - SMALLVAL)
  target_400_high = r_400 * (ONE + SMALLVAL)
  target_400_low = r_400 * (ONE - SMALLVAL)
  target_670_high = r_670 * (ONE + SMALLVAL)
  target_670_low = r_670 * (ONE - SMALLVAL)

! setup the mapping array for superbrick case (only invoked for Moho bottom)
  if (is_superbrick) then
    map_irem_ix_12=(/2,2,0,1,0,1,0,0/)
    map_irem_ix_34=(/1,1,0,2,0,2,0,0/)
    map_irem_iy_odd=(/1,2,0,1,0,2,0,0/)
    map_irem_iy_even=(/2,1,0,2,0,1,0,0/)
    if (USE_ONE_LAYER_SB) then
      nele_sub_block = 7
      ispec_list=(/1,2,4,6,8,9,11,13,15,16,18,20,22,23,25,27/)
   else
      nele_sub_block = 8
      ispec_list=(/1,2,4,6,9,10,12,14,17,18,20,22,25,26,28,30/)
    endif
    map_isub_ix=(/2,2,1,1/)
    map_isub_iy=(/2,1,2,1/)

    map_ix(1:NSPEC_DOUBLING_SUPERBRICK) = 0
    map_iy(1:NSPEC_DOUBLING_SUPERBRICK) = 0

    do i = 1, 16
      ispec_superbrick_current=ispec_list(i)
      isub_block = ispec_superbrick_current/nele_sub_block + 1
      irem_block = mod(ispec_superbrick_current,nele_sub_block)

      if (isub_block > 2) then
        irem_ix = map_irem_ix_34(irem_block)
      else
        irem_ix = map_irem_ix_12(irem_block)
      endif
      if (mod(isub_block,2) == 0) then
        irem_iy = map_irem_iy_even(irem_block)
      else
        irem_iy = map_irem_iy_odd(irem_block)
      endif
      map_ix(ispec_list(i)) = (map_isub_ix(isub_block) - 1) * 2 + irem_ix
      map_iy(ispec_list(i)) = (map_isub_iy(isub_block) - 1) * 2 + irem_iy
!      if (ispec_superbrick == 1 .and. myrank == 0) &
!                 write(*,'(10i4)') i, ispec_list(i), map_ix(ispec_list(i)), map_iy(ispec_list(i))
    enddo
  endif

! determine if the elements are on the discontinuity, and calculate the boundary jaocobian if needed
  if (.not. is_superbrick) then

! Moho top
    if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO .and. &
               abs(rmin-r_moho)/r_moho < SMALLVAL .and. r1 < target_moho_high .and. r2 < target_moho_high &
               .and. r3 < target_moho_high .and. r4 < target_moho_high) then
        ispec2D_moho_top = ispec2D_moho_top + 1
        ibelm_moho_top(ispec2D_moho_top) = ispec
        call compute_jacobian_2D(myrank,ispec2D_moho_top,xelm2,yelm2,zelm2,dershape2D_bottom, &
                   jacobian2D_moho,normal_moho,NGLLX,NGLLY,NSPEC2D_MOHO)
! 400 top
  else if (abs(rmin-r_400)/r_400 < SMALLVAL .and. r1 < target_400_high .and. r2 < target_400_high &
             .and. r3 < target_400_high .and. r4 < target_400_high) then
    ispec2D_400_top = ispec2D_400_top + 1
    ibelm_400_top(ispec2D_400_top) = ispec
    call compute_jacobian_2D(myrank,ispec2D_400_top,xelm2,yelm2,zelm2,dershape2D_bottom, &
               jacobian2D_400,normal_400,NGLLX,NGLLY,NSPEC2D_400)

! 400 bot
  else if (abs(rmax-r_400)/r_400 < SMALLVAL .and. r5 > target_400_low .and. r6 > target_400_low &
             .and. r7 > target_400_low .and. r8 > target_400_low) then
    ispec2D_400_bot = ispec2D_400_bot + 1
    ibelm_400_bot(ispec2D_400_bot) = ispec

! 670 top
  else if (abs(rmin-r_670)/r_670 < SMALLVAL .and. r1 < target_670_high .and. r2 < target_670_high &
             .and. r3 < target_670_high .and. r4 < target_670_high) then
    ispec2D_670_top = ispec2D_670_top + 1
    ibelm_670_top(ispec2D_670_top) = ispec
    call compute_jacobian_2D(myrank,ispec2D_670_top,xelm2,yelm2,zelm2,dershape2D_bottom, &
               jacobian2D_670,normal_670,NGLLX,NGLLY,NSPEC2D_670)
! 670 bot
  else if (abs(rmax-r_670)/r_670 < SMALLVAL .and. r5 > target_670_low .and. r6 > target_670_low &
             .and. r7 > target_670_low .and. r8 > target_670_low) then
    ispec2D_670_bot = ispec2D_670_bot + 1
    ibelm_670_bot(ispec2D_670_bot) = ispec
  endif

  else ! superbrick case
    ! Moho bot (special care should be taken to deal with mapping 2D element indices)
    if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO .and. &
               abs(rmax-r_moho)/r_moho < SMALLVAL .and. r5 > target_moho_low .and. r6 > target_moho_low &
               .and. r7 > target_moho_low .and. r8 > target_moho_low) then
      ispec2D_moho_bot = ispec2D_moho_bot + 1
      ix=map_ix(ispec_superbrick)
      iy=map_iy(ispec_superbrick)
      if (ix == 0 .or. iy == 0) call exit_mpi(myrank, 'Check (ix,iy) on the Moho bot is 0')
      ix_top = (ix_elem - 1)  + ix
      iy_top = (iy_elem - 1)  + iy
      ispec2D_moho_bot_map = (ix_top - 1) * nex_eta_moho + iy_top
!      if (myrank == 0) write(*,'(10i6)') ix_elem, iy_elem, ispec_superbrick, ix, iy, ix_top, iy_top, ispec2D_moho_bot_map
      ibelm_moho_bot(ispec2D_moho_bot_map) = ispec
    endif
  endif


end subroutine get_jacobian_discontinuities

