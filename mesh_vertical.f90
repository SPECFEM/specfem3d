!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 2
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine mesh_vertical(myrank,rn,NER,NER_BOTTOM_MOHO,NER_MOHO_16, &
                           NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM, &
!! DK DK UGLY modif z_top by Emmanuel Chaljub here
!! DK DK UGLY modif Manu removed                           z_top, &
                           Z_DEPTH_BLOCK,Z_BASEMENT_SURFACE,Z_DEPTH_MOHO,MOHO_MAP_LUPEI,MODEL)

! create the vertical mesh, honoring the major discontinuities in the basin

  implicit none

  include "constants.h"

  integer myrank
  integer NER,NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM
  logical MOHO_MAP_LUPEI
  double precision Z_DEPTH_BLOCK,Z_BASEMENT_SURFACE,Z_DEPTH_MOHO
  double precision rn(0:2*NER)

  character(len=150) MODEL

!! DK DK UGLY modif z_top by Emmanuel Chaljub here
!! DK DK UGLY modif Manu removed  double precision z_top

  integer npr,ir,NER_Lacq

  npr = -1

  if(MODEL /= 'Lacq_gas_field_France') then

!
!--- bottom of the mesh (Z_DEPTH_BLOCK) to Moho
!
  do ir=0,2*NER_BOTTOM_MOHO-1
    npr=npr+1
    rn(npr)=(Z_DEPTH_MOHO-Z_DEPTH_BLOCK)*dble(ir)/dble(2*NER_BOTTOM_MOHO)
  enddo

! do not use d16km when Moho map is honored
  if(MOHO_MAP_LUPEI) then

!
!--- Moho to modified basement surface
!
    do ir=0,2*(NER_MOHO_16+NER_16_BASEMENT)-1
      npr=npr+1
      rn(npr)=(Z_DEPTH_MOHO-Z_DEPTH_BLOCK) + (Z_BASEMENT_SURFACE-Z_DEPTH_MOHO)*dble(ir)/dble(2*(NER_MOHO_16+NER_16_BASEMENT))
    enddo

  else
!
!--- Moho to d16km
!
    do ir=0,2*NER_MOHO_16-1
      npr=npr+1
      rn(npr)=(Z_DEPTH_MOHO-Z_DEPTH_BLOCK) + (DEPTH_16km_SOCAL-Z_DEPTH_MOHO)*dble(ir)/dble(2*NER_MOHO_16)
    enddo
!
!--- d16km to modified basement surface
!
    do ir=0,2*NER_16_BASEMENT-1
      npr=npr+1
      rn(npr)=(DEPTH_16km_SOCAL-Z_DEPTH_BLOCK) + (Z_BASEMENT_SURFACE-DEPTH_16km_SOCAL)*dble(ir)/dble(2*NER_16_BASEMENT)
    enddo

  endif

!
!--- modified basement surface to surface of model (topography/bathymetry)
!
! also create last point exactly at the surface
! other regions above stop one point below
  do ir=0,2*(NER_BASEMENT_SEDIM+NER_SEDIM) - 0
    npr=npr+1
    rn(npr)=(Z_BASEMENT_SURFACE-Z_DEPTH_BLOCK) + &
!! DK DK UGLY modif z_top by Emmanuel Chaljub here
!! DK DK UGLY suppressed Manu's modif and put old code back because better mesh
!! DK DK UGLY investigate this in detail one day
 (Z_SURFACE-Z_BASEMENT_SURFACE)*dble(ir)/dble(2*(NER_BASEMENT_SEDIM+NER_SEDIM))
!! DK DK UGLY modif Manu removed     (z_top-Z_BASEMENT_SURFACE)*dble(ir)/dble(2*(NER_BASEMENT_SEDIM+NER_SEDIM))
  enddo

! for Lacq (France) gas field, case of Lacq_gas_field_France
  else if(MODEL == 'Lacq_gas_field_France') then

!
!--- bottom of the mesh (Z_DEPTH_BLOCK) to main mesh interface (called Moho for compatibility)
!
  do ir=0,2*NER_BOTTOM_MOHO-1
    npr=npr+1
    rn(npr)=(Z_DEPTH_MOHO-Z_DEPTH_BLOCK)*dble(ir)/dble(2*NER_BOTTOM_MOHO)
  enddo

!
!--- main mesh interface (called Moho for compatibility) to surface of model (topography)
!
! also create last point exactly at the surface
! other regions above stop one point below
  NER_Lacq = NER_BASEMENT_SEDIM + NER_SEDIM + NER_16_BASEMENT + NER_MOHO_16
  do ir=0,2*NER_Lacq - 0
    npr=npr+1
    rn(npr)=(Z_DEPTH_MOHO-Z_DEPTH_BLOCK) + (Z_SURFACE-Z_DEPTH_MOHO)*dble(ir)/dble(2*NER_Lacq)
  enddo

  endif  ! end of case of Lacq_gas_field_France

! normalize depths
!! DK DK UGLY modif z_top by Emmanuel Chaljub here
!! DK DK UGLY suppressed Manu's modif and put old code back because better mesh
!! DK DK UGLY investigate this in detail one day
!! DK DK UGLY modif Manu removed  rn(:) = rn(:) / (z_top-Z_DEPTH_BLOCK)
!! DK DK UGLY modif Manu removed
  rn(:) = rn(:) / (Z_SURFACE-Z_DEPTH_BLOCK)

! check that the mesh that has been generated is correct
  if(npr /= 2*NER) call exit_MPI(myrank,'incorrect intervals for basin')

! check that vertical spacing makes sense
  do ir=0,2*NER-1
    if(rn(ir+1) < rn(ir)) call exit_MPI(myrank,'incorrect vertical spacing for basin')
  enddo

  end subroutine mesh_vertical

