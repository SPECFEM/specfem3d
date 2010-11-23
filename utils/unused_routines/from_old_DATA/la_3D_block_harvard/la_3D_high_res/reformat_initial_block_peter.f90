
  program read_gocad_block_extract

!! DK DK UGLY Peter's block is inverted in vertical direction
!! DK DK UGLY this code puts it back into the right format

! when compiling this code on a PC, use "pgf90 -byteswapio " to byte-swap
! the original binary file, which was created on an SGI

  implicit none

  include "../../../constants.h"

!
! new hi-res Voxet Peter July 29, 2002
!

! AXIS_O 371052.25 3774000 400
! AXIS_U 46000 0 0
! AXIS_V 0 -48750 0
! AXIS_W 0 0 -9900
! AXIS_MIN 0 0 0
! AXIS_MAX 1 1 1
! AXIS_N 185 196 100
!

  real vp_block_gocad_sngl(0:NX_GOCAD_HR-1,0:NY_GOCAD_HR-1,0:NZ_GOCAD_HR-1)
  real vp_block_gocad_sngl_clean(0:NX_GOCAD_HR-1,0:NY_GOCAD_HR-1,0:NZ_GOCAD_HR-1)

! use integer array to store topography values
  integer itopo_bathy_basin(NX_TOPO,NY_TOPO)
  integer iclosestlong,iclosestlat
  double precision elevation,max_error
  double precision lat,long

  integer ix,iy,iz,iz_found,ipoin,ispec,npoin,nspec
  integer icount_undefined,imaterial

  double precision vpmin,vpmax
  double precision xcoord,ycoord,zcoord
  double precision zsedim_found
  integer irecord

  print *
  print *,'reading velocity block from Gocad voxet'
  print *

  print *
  print *,'number of points in block NX_GOCAD_HR,NY_GOCAD_HR,NZ_GOCAD_HR = ',NX_GOCAD_HR,NY_GOCAD_HR,NZ_GOCAD_HR
  print *,'total points in block NX_GOCAD_HR*NY_GOCAD_HR*NZ_GOCAD_HR = ',NX_GOCAD_HR*NY_GOCAD_HR*NZ_GOCAD_HR
  print *,'size of block in bytes = ',4*NX_GOCAD_HR*NY_GOCAD_HR*NZ_GOCAD_HR
  print *

!
!--- read Vp
!
! Gocad stores array in C-style, we read in Fortran so we need to transpose
  open(unit=27,file='LA_HR_VINT@@',status='old',access='direct',recl=4)
  irecord = 0
  do iz = 0,NZ_GOCAD_HR-1
  do iy = 0,NY_GOCAD_HR-1
  do ix = 0,NX_GOCAD_HR-1
    irecord = irecord + 1
    read(27,rec=irecord) vp_block_gocad_sngl(ix,iy,iz)

! use only one convention for threshold: vp > 6500. means fictitious
   if(vp_block_gocad_sngl(ix,iy,iz) < 0.1) vp_block_gocad_sngl(ix,iy,iz) = 6501.
   if(vp_block_gocad_sngl(ix,iy,iz) > 6499.) vp_block_gocad_sngl(ix,iy,iz) = 6501.

! invert Y and Z axes back to normal
    vp_block_gocad_sngl_clean(ix,NY_GOCAD_HR-1-iy,NZ_GOCAD_HR-1-iz) = vp_block_gocad_sngl(ix,iy,iz)

  enddo
  enddo
  enddo
  close(27)

! write block under different name
  open(unit=27,file='LA_HR_VINT_clean@@',status='unknown',access='direct',recl=4)
  irecord = 0
  do iz = 0,NZ_GOCAD_HR-1
  do iy = 0,NY_GOCAD_HR-1
  do ix = 0,NX_GOCAD_HR-1
    irecord = irecord + 1
    write(27,rec=irecord) vp_block_gocad_sngl_clean(ix,iy,iz)
  enddo
  enddo
  enddo
  close(27)

  end program read_gocad_block_extract


